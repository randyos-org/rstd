//! Triangular inverse, the products `U Uᴴ` and `Lᴴ L`, and triangular
//! solves: OpenBLAS's `interface/lapack/trtri.c`, `trti2.c`, `lauum.c`,
//! `lauu2.c`, `trtrs.c` and their `z` forms, over `lapack/trtri/`,
//! `lapack/trti2/`, `lapack/lauum/`, `lapack/lauu2/` and `lapack/trtrs/`
//! (OpenBLAS v0.3.34-317). The semantics are Reference-LAPACK's `xTRTRI`,
//! `xTRTI2`, `xLAUUM`, `xLAUU2` and `xTRTRS` (`lapack-netlib/SRC/`).
//!
//! The blocked forms call the level-3 drivers directly, as upstream calls
//! its kernels, so Accelerate never runs them; upstream's `lauum_U_single.c`
//! and `lauum_L_single.c` fuse their `herk` and `trmm` per column block, the
//! same arithmetic. The unblocked `trti2` and `lauu2` are loops, upstream's
//! `trmv`, `scal`, `dot` and `gemv` in the reference's order.
//!
//! Upstream takes `trtri_[UL]_parallel.c` and `lauum_[UL]_parallel.c`
//! (other blocks, and for `trtri` another order) when it has more than one
//! thread; here the size alone picks them (`chol.threads`), and the drivers'
//! workers keep the bits. `trtrs_parallel.c` is the single form's `trsm` on
//! threads, as here.
//!
//! Where they differ, the reference wins: `trtrs` checks its arguments in
//! the reference's order (upstream reports a bad `diag` before a bad
//! `uplo`), and a zero on the diagonal is found by `==`, where upstream asks
//! `amin` and `iamin`, whose NaN handling varies by kernel.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const driver = @import("../blas/blas3/driver.zig");
const triangular = @import("../blas/blas2/triangular.zig");
const Vtable = @import("../blas/Vtable.zig");
const chol = @import("cholesky.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Diag = types.Diag;
const Trans = types.Trans;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = chol.scalar;
const add = chol.add;
const mul = chol.mul;
const conj = chol.conj;

/// Bad arguments of `trtri`, `trti2`, `lauum` and `lauu2`.
pub const ArgError = chol.SquareError;
/// What `trti2`, `lauu2` and `lauum` return.
pub const LauumError = ArgError || Allocator.Error || Io.Cancelable;
/// What `trtri` returns: `Singular` for a zero on the diagonal.
pub const TrtriError = LauumError || error{Singular};
/// What `trtrs` returns: `Singular` for a zero on the diagonal.
pub const TrtrsError = error{ InvalidUplo, InvalidTrans, InvalidN, InvalidNrhs, InvalidLda, InvalidLdb, Singular } || Allocator.Error || Io.Cancelable;

// ---- trtri, trti2 ----------------------------------------------------------------

/// `?TRTRI`: `A^-1` of the `uplo` triangular `n x n` `A`, in place. A unit
/// `diag` is taken as ones and never read or written; the other triangle is
/// never read.
///
/// For a non-unit `diag`, a zero on the diagonal returns `error.Singular`
/// before `A` is touched, the first one's index (1-based) in `pivot`, the
/// reference's `INFO`. `pivot` is written only then.
pub fn trtri(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, diag: Diag, n: isize, a: [*]T, lda: isize, pivot: ?*usize) TrtriError!void {
    comptime chol.floatType(T, "trtri");
    if (!try chol.checkSquare(T, uplo, n, lda)) return;
    const un: usize = @intCast(n);
    const la: usize = @intCast(lda);
    if (diag == .n) if (firstZero(T, un, a, la)) |k| {
        if (pivot) |p| p.* = k;
        return error.Singular;
    };
    return inverse(allocator, io, T, uplo == .u, diag == .u, un, a, la, chol.threads(.trtri, T, un));
}

/// `?TRTI2`: `trtri` unblocked (`trti2_U.c` and the rest, as loops). No
/// singularity check, as in the reference: a zero on the diagonal divides.
pub fn trti2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, diag: Diag, n: isize, a: [*]T, lda: isize) LauumError!void {
    comptime chol.floatType(T, "trti2");
    _ = allocator;
    _ = io;
    if (!try chol.checkSquare(T, uplo, n, lda)) return;
    unblockedInverse(T, uplo == .u, diag == .u, @intCast(n), a, @intCast(lda));
}

/// The first zero on the diagonal of the `n x n` `a` (1-based), or null.
pub fn firstZero(comptime T: type, n: usize, a: [*]const T, lda: usize) ?usize {
    for (0..n) |i| if (chol.isZero(T, a[i * (lda + 1)])) return i + 1;
    return null;
}

/// `trtri`'s algorithm: with `want` (`chol.threads`), the parallel one on
/// up to that many workers, else the single one.
pub fn inverse(allocator: Allocator, io: Io, comptime T: type, upper: bool, unit: bool, n: usize, a: [*]T, lda: usize, want: ?u32) driver.Error!void {
    const vt = chol.table(T);
    const blk = driver.Blocking.of(T, vt);
    if (want) |w| return inverseParallel(allocator, io, T, vt, upper, unit, chol.dtbOf(vt), blk.q, n, a, lda, w);
    return inverseSingle(allocator, io, T, vt, upper, unit, blk.q, n, a, lda);
}

/// `trtri_U_single.c` and `trtri_L_single.c`: `trti2` up to `Q` (lower: below
/// it), else blocks of `Q`, left to right (lower: right to left), each with
/// its off-diagonal block multiplied by the inverse already made and solved
/// by its own diagonal block (`TRMM_LN*`, `TRSM_RN*` with -1), then
/// inverted.
fn inverseSingle(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, unit: bool, q: usize, n: usize, a: [*]T, lda: usize) driver.Error!void {
    if (if (upper) n <= q else n < q) return unblockedInverse(T, upper, unit, n, a, lda);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    // upper from the first block, lower from the last
    var j: usize = if (upper) 0 else (n - 1) / q * q;
    while (true) {
        const jb = @min(n - j, q);
        const d = a + j * (lda + 1);
        if (upper) {
            if (j > 0) {
                const b = a + j * lda;
                try driver.trmm(allocator, io, T, vt, true, true, false, false, unit, j, jb, one, a, lda, b, lda, 1);
                try driver.trsm(allocator, io, T, vt, false, true, false, false, unit, j, jb, minus, d, lda, b, lda, 1);
            }
        } else {
            const m = n - j - jb;
            if (m > 0) {
                const b = d + jb;
                try driver.trmm(allocator, io, T, vt, true, false, false, false, unit, m, jb, one, b + jb * lda, lda, b, lda, 1);
                try driver.trsm(allocator, io, T, vt, false, false, false, false, unit, m, jb, minus, d, lda, b, lda, 1);
            }
        }
        unblockedInverse(T, upper, unit, jb, d, lda);
        if (upper) {
            j += jb;
            if (j >= n) return;
        } else {
            if (j == 0) return;
            j -= q;
        }
    }
}

/// `trtri_U_parallel.c` and `trtri_L_parallel.c`: `trti2` up to
/// `dtb_entries`, else blocks of `Q`, or a quarter of `n` below `4 * Q`.
/// Each block's off-diagonal part is solved by its diagonal block (with
/// -1), that block inverted recursively, then the rest updated and
/// multiplied by it (`GEMM_NN`, `TRMM_LN*`), on up to `want` workers.
fn inverseParallel(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, unit: bool, dtb: usize, q: usize, n: usize, a: [*]T, lda: usize, want: u32) driver.Error!void {
    if (n <= dtb) return unblockedInverse(T, upper, unit, n, a, lda);
    const blocking = if (n < 4 * q) (n + 3) / 4 else q;
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    var i: usize = if (upper) 0 else (n - 1) / blocking * blocking;
    while (true) {
        const bk = @min(n - i, blocking);
        const d = a + i * (lda + 1);
        if (upper) {
            const rest = n - i - bk;
            if (i > 0) try driver.trsm(allocator, io, T, vt, false, true, false, false, unit, i, bk, minus, d, lda, a + i * lda, lda, want);
            try inverseParallel(allocator, io, T, vt, upper, unit, dtb, q, bk, d, lda, want);
            if (rest > 0) {
                const b = d + bk * lda;
                if (i > 0) try driver.gemm(allocator, io, T, vt, false, false, .n, i, rest, bk, one, a + i * lda, lda, b, lda, one, a + (i + bk) * lda, lda, want);
                try driver.trmm(allocator, io, T, vt, true, true, false, false, unit, bk, rest, one, d, lda, b, lda, want);
            }
        } else {
            const m = n - bk - i;
            if (m > 0) try driver.trsm(allocator, io, T, vt, false, false, false, false, unit, m, bk, minus, d, lda, d + bk, lda, want);
            try inverseParallel(allocator, io, T, vt, upper, unit, dtb, q, bk, d, lda, want);
            if (i > 0) {
                if (m > 0) try driver.gemm(allocator, io, T, vt, false, false, .n, m, i, bk, one, d + bk, lda, a + i, lda, one, a + i + bk, lda, want);
                try driver.trmm(allocator, io, T, vt, true, false, false, false, unit, bk, i, one, d, lda, a + i, lda, want);
            }
        }
        if (upper) {
            i += bk;
            if (i >= n) return;
        } else {
            if (i == 0) return;
            i -= blocking;
        }
    }
}

/// `trti2_U.c`, `trti2_L.c`, `ztrti2_U.c`, `ztrti2_L.c`: column by column
/// (upper left to right, lower right to left), the diagonal element
/// inverted, then the column below (upper: above) it multiplied by the
/// inverse already made (`TRMV_N*`) and by `-1 / A(j, j)` (`SCAL_K`).
fn unblockedInverse(comptime T: type, upper: bool, unit: bool, n: usize, a: [*]T, lda: usize) void {
    for (0..n) |jj| {
        const j = if (upper) jj else n - 1 - jj;
        const dj = &a[j * (lda + 1)];
        const ajj: T = if (unit) scalar(T, 1) else recip(T, dj.*);
        if (!unit) dj.* = ajj;
        const x = a + j * lda;
        const s = chol.neg(T, ajj);
        if (upper) {
            // `x := U(:j, :j) x`, a column at a time, as `xTRMV`
            for (0..j) |k| {
                const t = x[k];
                const col = a + k * lda;
                for (0..k) |i| x[i] = add(T, x[i], mul(T, t, col[i]));
                if (!unit) x[k] = mul(T, t, col[k]);
            }
            for (0..j) |i| x[i] = mul(T, s, x[i]);
        } else {
            // `x := L(j+1:, j+1:) x`, from the last column back
            var k = n;
            while (k > j + 1) {
                k -= 1;
                const t = x[k];
                const col = a + k * lda;
                for (k + 1..n) |i| x[i] = add(T, x[i], mul(T, t, col[i]));
                if (!unit) x[k] = mul(T, t, col[k]);
            }
            for (j + 1..n) |i| x[i] = mul(T, s, x[i]);
        }
    }
}

/// `1 / v`: complex as `ztrti2_U.c` has it, by the ratio of the smaller part
/// to the larger.
fn recip(comptime T: type, v: T) T {
    if (comptime !isComplex(T)) return 1 / v;
    if (@abs(v.re) >= @abs(v.im)) {
        const ratio = v.im / v.re;
        const den = 1 / (v.re * (1 + ratio * ratio));
        return .init(den, -ratio * den);
    }
    const ratio = v.re / v.im;
    const den = 1 / (v.im * (1 + ratio * ratio));
    return .init(ratio * den, -den);
}

// ---- lauum, lauu2 ----------------------------------------------------------------

/// `?LAUUM`: `U Uᴴ` (`u`) or `Lᴴ L` (`l`) of the triangular `n x n` factor
/// in `A`'s `uplo` triangle, into that triangle. The other triangle is never
/// read. Complex: the diagonal comes out real, as `xLAUUM`'s.
pub fn lauum(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize) LauumError!void {
    comptime chol.floatType(T, "lauum");
    if (!try chol.checkSquare(T, uplo, n, lda)) return;
    const un: usize = @intCast(n);
    return product(allocator, io, T, uplo == .u, un, a, @intCast(lda), chol.threads(.lauum, T, un));
}

/// `?LAUU2`: `lauum` unblocked (`lauu2_U.c` and the rest, as loops).
pub fn lauu2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize) LauumError!void {
    comptime chol.floatType(T, "lauu2");
    _ = allocator;
    _ = io;
    if (!try chol.checkSquare(T, uplo, n, lda)) return;
    unblockedProduct(T, uplo == .u, @intCast(n), a, @intCast(lda));
}

/// `lauum`'s algorithm: with `want` (`chol.threads`), the parallel one on
/// up to that many workers, else the single one.
pub fn product(allocator: Allocator, io: Io, comptime T: type, upper: bool, n: usize, a: [*]T, lda: usize, want: ?u32) driver.Error!void {
    const vt = chol.table(T);
    const blk = driver.Blocking.of(T, vt);
    if (want) |w| return productParallel(allocator, io, T, vt, upper, blk, n, a, lda, w);
    return productSingle(allocator, io, T, vt, upper, chol.dtbOf(vt), blk.q, n, a, lda);
}

/// `lauum_U_single.c` and `lauum_L_single.c`: `lauu2` up to `dtb_entries`,
/// else blocks of `Q`, or a quarter of `n` up to `4 * Q`, each first adding
/// its off-diagonal part's product into the part already made, then
/// multiplied by its diagonal block, then that block recursively.
fn productSingle(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, dtb: usize, q: usize, n: usize, a: [*]T, lda: usize) driver.Error!void {
    if (n <= dtb) return unblockedProduct(T, upper, n, a, lda);
    const blocking = if (n <= 4 * q) (n + 3) / 4 else q;
    var j: usize = 0;
    while (j < n) : (j += blocking) {
        const bk = @min(n - j, blocking);
        try offDiagonal(allocator, io, T, vt, upper, j, bk, a, lda, 1);
        try productSingle(allocator, io, T, vt, upper, dtb, q, bk, a + j * (lda + 1), lda);
    }
}

/// `lauum_U_parallel.c` and `lauum_L_parallel.c`: the single form up to
/// `2 * unroll_n`, else blocks of half `n` rounded up to `unroll_n`, at most
/// `Q`, each block as there, recursively, on up to `want` workers. Upstream
/// drops to the single form for one thread too; here only the size does.
fn productParallel(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, blk: driver.Blocking, n: usize, a: [*]T, lda: usize, want: u32) driver.Error!void {
    if (n <= 2 * blk.un) return productSingle(allocator, io, T, vt, upper, chol.dtbOf(vt), blk.q, n, a, lda);
    const blocking = @min((n / 2 + blk.un - 1) / blk.un * blk.un, blk.q);
    var i: usize = 0;
    while (i < n) : (i += blocking) {
        const bk = @min(n - i, blocking);
        try offDiagonal(allocator, io, T, vt, upper, i, bk, a, lda, want);
        try productParallel(allocator, io, T, vt, upper, blk, bk, a + i * (lda + 1), lda, want);
    }
}

/// Block `j` (`bk` wide) of `lauum`: with `C` its columns above the
/// diagonal (`u`) or rows left of it, `A(:j, :j) += C Cᴴ` (`HERK_UN`; `l`:
/// `Cᴴ C`, `HERK_LC`), then `C := C U(j, j)ᴴ` (`TRMM_RCUN`; `l`:
/// `L(j, j)ᴴ C`, `TRMM_LCLN`).
fn offDiagonal(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, j: usize, bk: usize, a: [*]T, lda: usize, want: u32) driver.Error!void {
    if (j == 0) return;
    const complex = comptime isComplex(T);
    const one = scalar(T, 1);
    const d = a + j * (lda + 1);
    const c = if (upper) a + j * lda else a + j;
    try driver.syrk(allocator, io, T, vt, complex, upper, !upper, j, bk, one, c, lda, null, 0, one, a, lda, want);
    if (upper) {
        try driver.trmm(allocator, io, T, vt, false, true, true, complex, false, j, bk, one, d, lda, c, lda, want);
    } else {
        try driver.trmm(allocator, io, T, vt, true, false, true, complex, false, bk, j, one, d, lda, c, lda, want);
    }
}

/// `lauu2_U.c`, `lauu2_L.c`, `zlauu2_U.c`, `zlauu2_L.c`: for each `i`, column
/// `i` to the diagonal (`l`: row `i`) times the real `A(i, i)` (`SCAL_K`),
/// then, but for the last, `A(i, i)` plus the rest of row (`l`: column)
/// `i`'s squared magnitudes (`DOTU_K`, `DOTC_K`; the diagonal comes out
/// real) and the column (`l`: row) above (left of) it plus the block right
/// of (below) it times that row conjugated (`GEMV_N`, `GEMV_O`; `GEMV_T`,
/// `GEMV_U`).
fn unblockedProduct(comptime T: type, upper: bool, n: usize, a: [*]T, lda: usize) void {
    const F = Float(T);
    for (0..n) |i| {
        const d = &a[i * (lda + 1)];
        const aii = chol.re(T, d.*);
        // column `i` (`l`: row `i`) up to and with the diagonal
        const s: usize = if (upper) 1 else lda;
        const v = if (upper) a + i * lda else a + i;
        for (0..i + 1) |k| v[k * s] = chol.mulReal(T, v[k * s], aii);
        if (i + 1 == n) continue;
        // the rest of row `i` (`l`: column `i`)
        const rs: usize = if (upper) lda else 1;
        const r = if (upper) a + i + (i + 1) * lda else a + i + 1 + i * lda;
        const len = n - i - 1;
        var sum: F = 0;
        for (0..len) |k| sum += chol.abs2(T, r[k * rs]);
        d.* = scalar(T, chol.re(T, d.*) + sum);
        if (upper) {
            // `A(:i, i) += A(:i, i+1:) conj(A(i, i+1:))`, a column at a time
            for (0..len) |k| {
                const t = conj(T, r[k * lda]);
                const col = a + (i + 1 + k) * lda;
                for (0..i) |row| v[row] = add(T, v[row], mul(T, t, col[row]));
            }
        } else {
            // `A(i, :i) += A(i+1:, :i)ᵀ conj(A(i+1:, i))`, a dot each
            for (0..i) |c| {
                const col = a + i + 1 + c * lda;
                var t = scalar(T, 0);
                for (0..len) |k| t = add(T, t, mul(T, col[k], conj(T, r[k])));
                v[c * lda] = add(T, v[c * lda], t);
            }
        }
    }
}

// ---- trtrs -----------------------------------------------------------------------

/// `?TRTRS`: solves `op(A) X = B` for the `n x nrhs` `X`, in place of `B`,
/// `A` triangular (`uplo`, `diag`). `trans` is `n`, `t`, `r` (`conj(A)`) or
/// `c` (`Aᴴ`); real `r` and `c` are `n` and `t`, as upstream takes them,
/// where the reference takes no `r`. One right-hand side is `trsv`, as
/// upstream's `trtrs_single.c` has it.
///
/// For a non-unit `diag`, a zero on the diagonal returns `error.Singular`
/// before `B` is touched (even for `nrhs = 0`, as the reference), its index
/// (1-based) in `pivot`.
pub fn trtrs(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, nrhs: isize, a: [*]const T, lda: isize, b: [*]T, ldb: isize, pivot: ?*usize) TrtrsError!void {
    comptime chol.floatType(T, "trtrs");
    if (uplo.isRev()) return error.InvalidUplo;
    if (trans.conjX()) return error.InvalidTrans;
    if (n < 0) return error.InvalidN;
    if (nrhs < 0) return error.InvalidNrhs;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (ldb < @max(1, n)) return error.InvalidLdb;
    if (n == 0) return;
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(n, nrhs, ldb, @sizeOf(T))) return error.InvalidLdb;
    const m: usize = @intCast(n);
    const la: usize = @intCast(lda);
    if (diag == .n) if (firstZero(T, m, a, la)) |k| {
        if (pivot) |p| p.* = k;
        return error.Singular;
    };
    if (nrhs == 0) return;
    const vt = chol.table(T);
    if (nrhs == 1) return triangular.trsvTable(allocator, io, T, vt, uplo, trans, diag, m, a, la, b);
    const want = chol.threads(.trtrs, T, m) orelse 1;
    try driver.trsm(allocator, io, T, vt, true, uplo == .u, trans.transposed(), trans.conjA(), diag == .u, m, @intCast(nrhs), scalar(T, 1), a, la, b, @intCast(ldb), want);
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const blas = @import("../blas/__root__.zig");
const check = @import("../blas/_test/blas1/check.zig");
const guard = @import("../blas/_test/blas2/guard.zig");

const W = chol.W;
const test_types = chol.test_types;
const poisoned = chol.poisoned;
const isNanT = chol.isNanT;

const nan = std.math.nan(f64);

/// A test triangle in `upper`'s side of `n x n`, stored at `lda = n + 2`:
/// the diagonal in [2.5, 3.5), the rest within `1 / n`, so it and its
/// inverse are well conditioned; `w` keeps it whole in `W(T)`, a unit
/// diagonal as ones. The other triangle, the padding and a unit diagonal
/// hold NaN.
fn Tri(comptime T: type) type {
    return struct {
        n: usize,
        lda: usize,
        upper: bool,
        unit: bool,
        a: []T,
        w: []W(T),

        const Self = @This();

        fn init(al: Allocator, n: usize, upper: bool, unit: bool, seed: u64) !Self {
            const lda = n + 2;
            const a = try al.alloc(T, @max(1, lda * n));
            errdefer al.free(a);
            const w = try al.alloc(W(T), n * n);
            errdefer al.free(w);
            const v = try al.alloc(T, n * n);
            defer al.free(v);
            check.fill(T, v, seed);
            const scale: Float(T) = @floatCast(1 / @as(f64, @floatFromInt(@max(n, 1))));
            for (0..n) |j| for (0..n) |i| {
                const x = v[i + j * n];
                const e = if (if (upper) i > j else i < j)
                    scalar(T, 0)
                else if (i == j)
                    scalar(T, if (unit) 1 else 3 + chol.re(T, x) / 2)
                else
                    chol.mulReal(T, x, scale);
                w[i + j * n] = check.wide(T, e);
            };
            var self: Self = .{ .n = n, .lda = lda, .upper = upper, .unit = unit, .a = a, .w = w };
            self.store();
            return self;
        }

        fn deinit(self: *Self, al: Allocator) void {
            al.free(self.a);
            al.free(self.w);
        }

        /// Whether the routines read and write `(i, j)` of `a`.
        fn stored(self: *const Self, i: usize, j: usize) bool {
            if (i >= self.n) return false;
            if (i == j) return !self.unit;
            return if (self.upper) i < j else i > j;
        }

        fn store(self: *Self) void {
            @memset(self.a, poisoned(T));
            for (0..self.n) |j| for (0..self.n) |i| if (self.stored(i, j)) {
                self.a[i + j * self.lda] = chol.toT(T, self.w[i + j * self.n]);
            };
        }

        /// The computed matrix's `(i, j)`: ones on a unit diagonal, zeros
        /// outside the triangle.
        fn got(self: *const Self, i: usize, j: usize) W(T) {
            if (i == j and self.unit) return check.wide(T, scalar(T, 1));
            if (!self.stored(i, j)) return check.wide(T, scalar(T, 0));
            return check.wide(T, self.a[i + j * self.lda]);
        }

        fn expectUntouched(self: *const Self, what: []const u8) !void {
            for (0..self.n) |j| for (0..self.lda) |i| if (!self.stored(i, j) and !isNanT(T, self.a[i + j * self.lda])) {
                check.say("{s} {s} {s} n={d}: wrote ({d},{d}), outside its triangle\n", .{ what, @typeName(T), if (self.upper) "u" else "l", self.n, i, j });
                return error.TestUnexpectedResult;
            };
        }

        /// `A X = I` within `8 (n + 2)` ulps of `|A| |X|`, `X` the computed
        /// inverse, nothing outside the triangle written.
        fn expectInverse(self: *const Self) !void {
            const n = self.n;
            const eps: f128 = std.math.floatEps(Float(T));
            const ulps: f128 = @floatFromInt(8 * (n + 2));
            for (0..n) |j| for (0..n) |i| {
                if (if (self.upper) i > j else i < j) continue;
                var s: W(T) = check.wide(T, scalar(T, if (i == j) -1 else 0));
                var mag: f128 = 0;
                for (@min(i, j)..@max(i, j) + 1) |k| {
                    const x = self.w[i + k * n];
                    const y = self.got(k, j);
                    s = check.addW(s, check.mulW(x, y));
                    mag += check.abs1(x) * check.abs1(y);
                }
                const d = check.abs1(s);
                if (!(d <= ulps * eps * mag)) {
                    check.say("trtri {s} {s}{s} n={d}: (A X - I)({d},{d}) = {e}, allowed {e}\n", .{ @typeName(T), if (self.upper) "u" else "l", if (self.unit) "u" else "n", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
                    return error.TestUnexpectedResult;
                }
            };
            try self.expectUntouched("trtri");
        }

        /// The stored triangle against `U Uᴴ` or `Lᴴ L` of `w` within
        /// `2 (n + 2)` ulps of the product in magnitude, a complex diagonal
        /// exactly real, nothing outside the triangle written.
        fn expectProduct(self: *const Self) !void {
            const n = self.n;
            const eps: f128 = std.math.floatEps(Float(T));
            const ulps: f128 = @floatFromInt(2 * (n + 2));
            for (0..n) |j| for (0..n) |i| {
                if (!self.stored(i, j)) continue;
                var s: W(T) = check.wide(T, scalar(T, 0));
                var mag: f128 = 0;
                for (@max(i, j)..n) |k| {
                    // `U(i, k) conj(U(j, k))`, or `conj(L(k, i)) L(k, j)`
                    const x = if (self.upper) self.w[i + k * n] else check.conjW(self.w[k + i * n]);
                    const y = if (self.upper) check.conjW(self.w[j + k * n]) else self.w[k + j * n];
                    s = check.addW(s, check.mulW(x, y));
                    mag += check.abs1(x) * check.abs1(y);
                }
                const g = self.a[i + j * self.lda];
                const d = check.abs1(check.subW(check.wide(T, g), s));
                if (!(d <= ulps * eps * mag)) {
                    check.say("lauum {s} {s} n={d}: ({d},{d}) off by {e}, allowed {e}\n", .{ @typeName(T), if (self.upper) "u" else "l", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
                    return error.TestUnexpectedResult;
                }
                if (comptime isComplex(T)) if (i == j and g.im != 0) {
                    check.say("lauum {s} n={d}: ({d},{d}) has imaginary part {e}\n", .{ @typeName(T), n, i, i, g.im });
                    return error.TestUnexpectedResult;
                };
            };
            try self.expectUntouched("lauum");
        }
    };
}

/// Sizes around every block edge: on `chol.shrink`'s blocking (`Q` 3,
/// `dtb_entries` 4) and on a core's own.
const shrunk_sizes = [_]usize{ 1, 2, 3, 4, 5, 8, 9, 12, 13, 14, 23 };
const own_sizes = [_]usize{ 1, 33, 70, 130 };

/// `f` over every test type, triangle, blocking and size, with the single
/// algorithm (`want` null) and the parallel one on one worker.
fn eachTri(vt: *const Vtable, comptime unit_too: bool, comptime f: anytype) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]*const Vtable{ chol.shrink(vt), vt }, 0..) |v, vi| {
        const saved = blas.testUse(v);
        defer blas.testRestore(saved);
        const sizes: []const usize = if (vi == 0) &shrunk_sizes else &own_sizes;
        for (sizes) |n| for ([_]bool{ true, false }) |upper| for ([_]bool{ false, true }) |unit| {
            if (unit and !unit_too) continue;
            if (check.rig() and n > 70) continue;
            for ([_]?u32{ null, 1 }) |want| {
                var t = try Tri(T).init(al, n, upper, unit, n * 13 + @as(u64, @intFromBool(upper)) * 2 + @intFromBool(unit));
                defer t.deinit(al);
                f(T, &t, want) catch |e| {
                    check.say("{s}, blocking {s}, {s} algorithm\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own", if (want == null) "single" else "parallel" });
                    return e;
                };
            }
        };
    };
}

fn inverseCase(comptime T: type, t: *Tri(T), want: ?u32) !void {
    var poison: guard.Poison = .{ .child = testing.allocator };
    try inverse(poison.allocator(), check.io, T, t.upper, t.unit, t.n, t.a.ptr, t.lda, want);
    try t.expectInverse();
}

fn inverseCases(vt: *const Vtable) anyerror!void {
    return eachTri(vt, true, inverseCase);
}

test "trtri: A X = I, both triangles and diagonals, both algorithms, on every linked core" {
    try check.eachCore(inverseCases);
}

fn productCase(comptime T: type, t: *Tri(T), want: ?u32) !void {
    var poison: guard.Poison = .{ .child = testing.allocator };
    try product(poison.allocator(), check.io, T, t.upper, t.n, t.a.ptr, t.lda, want);
    try t.expectProduct();
}

fn productCases(vt: *const Vtable) anyerror!void {
    return eachTri(vt, false, productCase);
}

test "lauum: U Uᴴ and Lᴴ L, both algorithms, on every linked core" {
    try check.eachCore(productCases);
}

test "trti2 and lauu2 are the unblocked forms, and every size picks an algorithm" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| for ([_]bool{ true, false }) |upper| {
        var t = try Tri(T).init(al, 9, upper, false, 3);
        defer t.deinit(al);
        try trti2(al, check.io, T, if (upper) .u else .l, .n, 9, t.a.ptr, @intCast(t.lda));
        try t.expectInverse();
        t.store();
        try lauu2(al, check.io, T, if (upper) .u else .l, 9, t.a.ptr, @intCast(t.lda));
        try t.expectProduct();
        // the public routines, past upstream's thresholds too
        for ([_]usize{ 40, 200 }) |n| {
            if (check.rig() and n > 70) continue;
            var p = try Tri(T).init(al, n, upper, false, n);
            defer p.deinit(al);
            try trtri(al, check.io, T, if (upper) .u else .l, .n, @intCast(n), p.a.ptr, @intCast(p.lda), null);
            try p.expectInverse();
            p.store();
            try lauum(al, check.io, T, if (upper) .u else .l, @intCast(n), p.a.ptr, @intCast(p.lda));
            try p.expectProduct();
        }
    };
}

test "trtri: a zero on the diagonal, the first one's index, nothing written; a unit diagonal is not read" {
    const al = testing.failing_allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]Uplo{ .u, .l }) |uplo| {
        var a: [16]T = @splat(scalar(T, 0));
        for (0..4) |i| a[i * 5] = scalar(T, 2);
        a[5] = scalar(T, 0);
        a[15] = scalar(T, 0);
        // a part only zero is not zero
        if (comptime isComplex(T)) a[0] = .init(0, 1);
        const before = a;
        var pivot: usize = 99;
        try testing.expectError(error.Singular, trtri(al, check.io, T, uplo, .n, 4, &a, 4, &pivot));
        try testing.expectEqual(@as(usize, 2), pivot);
        try testing.expectEqualSlices(T, &before, &a);
        // `u`: the zeros are ones; `pivot` is written only on failure
        pivot = 99;
        try trtri(testing.allocator, check.io, T, uplo, .u, 4, &a, 4, &pivot);
        try testing.expectEqual(@as(usize, 99), pivot);
        for (0..4) |i| try testing.expect(chol.sameBits(T, before[i * 5], a[i * 5]));
    };
}

test "trtri, trti2, lauum, lauu2: the first bad argument, nothing read for n = 0" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    try testing.expectError(error.InvalidUplo, trtri(fa, check.io, f64, .v, .n, -1, &a, 0, null));
    try testing.expectError(error.InvalidN, trtri(fa, check.io, f64, .u, .n, -1, &a, 0, null));
    try testing.expectError(error.InvalidLda, trtri(fa, check.io, f64, .l, .u, 2, &a, 1, null));
    try testing.expectError(error.InvalidLda, trtri(fa, check.io, f64, .l, .u, 2, &a, std.math.maxInt(isize) / 2, null));
    try testing.expectError(error.InvalidUplo, trti2(fa, check.io, f64, .m, .n, -1, &a, 0));
    try testing.expectError(error.InvalidN, trti2(fa, check.io, f64, .u, .u, -1, &a, 0));
    try testing.expectError(error.InvalidLda, trti2(fa, check.io, f64, .u, .n, 2, &a, 1));
    try testing.expectError(error.InvalidUplo, lauum(fa, check.io, f64, .v, -1, &a, 0));
    try testing.expectError(error.InvalidN, lauum(fa, check.io, f64, .l, -1, &a, 0));
    try testing.expectError(error.InvalidLda, lauum(fa, check.io, f64, .l, 2, &a, 1));
    try testing.expectError(error.InvalidUplo, lauu2(fa, check.io, f64, .m, -1, &a, 0));
    try testing.expectError(error.InvalidLda, lauu2(fa, check.io, f64, .u, 0, &a, 0));
    try trtri(fa, check.io, f64, .u, .n, 0, &a, 1, null);
    try trti2(fa, check.io, f64, .l, .n, 0, &a, 1);
    try lauum(fa, check.io, f64, .u, 0, &a, 1);
    try lauu2(fa, check.io, f64, .l, 0, &a, 1);
    for (a) |v| try testing.expect(std.math.isNan(v));
}

/// `trtrs` on `t` for `nrhs` right-hand sides at `ldb = n + 1` (the padding
/// row NaN): `op(A) X - B` within `8 (n + 2)` ulps of `|op(A)| |X|`, the
/// padding and `A` untouched.
fn solveCase(comptime T: type, al: Allocator, t: *const Tri(T), trans: Trans, nrhs: usize, seed: u64) !void {
    const n = t.n;
    const ldb = n + 1;
    const b = try al.alloc(T, ldb * nrhs);
    defer al.free(b);
    const b0 = try al.alloc(T, n * nrhs);
    defer al.free(b0);
    check.fill(T, b0, seed);
    for (0..nrhs) |j| {
        @memcpy(b[j * ldb ..][0..n], b0[j * n ..][0..n]);
        b[j * ldb + n] = poisoned(T);
    }
    var poison: guard.Poison = .{ .child = al };
    const uplo: Uplo = if (t.upper) .u else .l;
    const diag: Diag = if (t.unit) .u else .n;
    try trtrs(poison.allocator(), check.io, T, uplo, trans, diag, @intCast(n), @intCast(nrhs), t.a.ptr, @intCast(t.lda), b.ptr, @intCast(ldb), null);
    const tr = trans.transposed();
    const cj = isComplex(T) and trans.conjA();
    const eps: f128 = std.math.floatEps(Float(T));
    const ulps: f128 = @floatFromInt(8 * (n + 2));
    for (0..nrhs) |j| {
        try testing.expect(isNanT(T, b[j * ldb + n]));
        for (0..n) |i| {
            var s: W(T) = check.wide(T, b0[i + j * n]);
            var mag: f128 = 0;
            for (0..n) |k| {
                const v = if (tr) t.w[k + i * n] else t.w[i + k * n];
                const x = check.wide(T, b[k + j * ldb]);
                s = check.subW(s, check.mulW(if (cj) check.conjW(v) else v, x));
                mag += check.abs1(v) * check.abs1(x);
            }
            const d = check.abs1(s);
            if (!(d <= ulps * eps * mag)) {
                check.say("trtrs {s} {s}{s} {s} n={d} nrhs={d}: residual ({d},{d}) {e}, allowed {e}\n", .{ @typeName(T), @tagName(uplo), @tagName(diag), @tagName(trans), n, nrhs, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
                return error.TestUnexpectedResult;
            }
        }
    }
    try t.expectUntouched("trtrs");
}

fn solveCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const saved = blas.testUse(vt);
        defer blas.testRestore(saved);
        for ([_]usize{ 1, 2, 7, 33, 70, 130 }) |n| for ([_]bool{ true, false }) |upper| for ([_]bool{ false, true }) |unit| {
            if (check.rig() and n > 70) continue;
            var t = try Tri(T).init(al, n, upper, unit, n * 7 + @intFromBool(unit));
            defer t.deinit(al);
            for ([_]Trans{ .n, .t, .r, .c }) |trans| for ([_]usize{ 1, 3 }) |nrhs| {
                try solveCase(T, al, &t, trans, nrhs, n + nrhs);
            };
        };
    }
}

test "trtrs: op(A) X = B, every op, both triangles and diagonals, one and several right-hand sides, on every linked core" {
    try check.eachCore(solveCases);
}

test "trtrs: a zero on the diagonal before B is touched, even with no right-hand side; the first bad argument" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = .{ 1, 2, 3, 0 };
    var b: [4]f64 = @splat(nan);
    var pivot: usize = 99;
    try testing.expectError(error.Singular, trtrs(fa, check.io, f64, .u, .n, .n, 2, 2, &a, 2, &b, 2, &pivot));
    try testing.expectEqual(@as(usize, 2), pivot);
    pivot = 99;
    try testing.expectError(error.Singular, trtrs(fa, check.io, f64, .l, .t, .n, 2, 0, &a, 2, &b, 2, &pivot));
    try testing.expectEqual(@as(usize, 2), pivot);
    // a unit diagonal is not read
    try trtrs(fa, check.io, f64, .u, .n, .u, 2, 0, &a, 2, &b, 2, &pivot);
    for (b) |v| try testing.expect(std.math.isNan(v));
    // the reference's order: `uplo` before `diag` and `trans`
    try testing.expectError(error.InvalidUplo, trtrs(fa, check.io, f64, .v, .o, .n, -1, -1, &a, 0, &b, 0, null));
    try testing.expectError(error.InvalidTrans, trtrs(fa, check.io, f64, .u, .o, .n, -1, -1, &a, 0, &b, 0, null));
    try testing.expectError(error.InvalidN, trtrs(fa, check.io, f64, .u, .n, .n, -1, -1, &a, 0, &b, 0, null));
    try testing.expectError(error.InvalidNrhs, trtrs(fa, check.io, f64, .u, .n, .n, 2, -1, &a, 0, &b, 0, null));
    try testing.expectError(error.InvalidLda, trtrs(fa, check.io, f64, .u, .n, .n, 2, 1, &a, 1, &b, 0, null));
    try testing.expectError(error.InvalidLdb, trtrs(fa, check.io, f64, .u, .n, .n, 2, 1, &a, 2, &b, 1, null));
    try testing.expectError(error.InvalidLdb, trtrs(fa, check.io, f64, .u, .n, .n, 2, 2, &a, 2, &b, std.math.maxInt(isize) / 2, null));
    try trtrs(fa, check.io, f64, .u, .n, .n, 0, 2, &a, 1, &b, 1, null);
}

/// Inputs of a threaded-bits case: a triangle, copied into each run's `out`.
fn Threaded(comptime T: type) type {
    return struct {
        t: *const Tri(T),

        fn inv(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            @memcpy(out, c.t.a);
            try inverse(testing.allocator, io, T, c.t.upper, c.t.unit, c.t.n, out.ptr, c.t.lda, want);
        }

        fn prod(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            @memcpy(out, c.t.a);
            try product(testing.allocator, io, T, c.t.upper, c.t.n, out.ptr, c.t.lda, want);
        }

        fn solve(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            // `trtrs` threads by `blas.config`'s workers
            const saved = blas.config.cpu_count;
            defer blas.config.cpu_count = saved;
            blas.config.cpu_count = want;
            const n = c.t.n;
            check.fill(T, out, 9);
            try trtrs(testing.allocator, io, T, if (c.t.upper) .u else .l, .c, .n, @intCast(n), @intCast(out.len / n), c.t.a.ptr, @intCast(c.t.lda), out.ptr, @intCast(n), null);
        }
    };
}

test "trtri's and lauum's parallel algorithms and trtrs give one worker's bits at any worker count" {
    _ = blas.coreName();
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        const n: usize = if (check.rig()) 70 else 130;
        var t = try Tri(T).init(al, n, upper, false, 21);
        defer t.deinit(al);
        const c: Threaded(T) = .{ .t = &t };
        try chol.expectSameBits(T, c, Threaded(T).inv, t.a.len, "trtri");
        try chol.expectSameBits(T, c, Threaded(T).prod, t.a.len, "lauum");
        try chol.expectSameBits(T, c, Threaded(T).solve, n * 37, "trtrs");
    };
}
