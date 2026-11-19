//! Triangular inverse, `U Uᴴ` and `Lᴴ L`, and triangular solves: OpenBLAS
//! v0.3.34-317's `interface/lapack/trtri.c`, `trti2.c`, `lauum.c`,
//! `lauu2.c`, `trtrs.c` and their `z` forms over `lapack/trtri/` and the
//! rest, with Reference-LAPACK's semantics (`lapack-netlib/SRC/`).
//! Departures and threading: design.md's LAPACK section.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const driver = @import("../blas/blas3/driver.zig");
const triangular = @import("../blas/blas2/triangular.zig");
const Vtable = @import("../blas/Vtable.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Diag = types.Diag;
const Trans = types.Trans;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const add = ops.add;
const mul = ops.mul;
const conj = ops.conj;

/// What `trti2` and `lauu2` return: a reversed `uplo`, a slice too short or
/// a leading dimension too small.
pub const ArgError = interface.TriangleError;
/// What `lauum` returns: `ArgError`, and the allocator's and cancellation.
pub const LauumError = ArgError || Allocator.Error || Io.Cancelable;
/// What `trtri` returns: `Singular` for a zero on the diagonal.
pub const TrtriError = LauumError || interface.SingularError;
/// What `trtrs` returns: `Singular` for a zero on the diagonal, a `trans`
/// it does not take, and `TrtriError`'s others.
pub const TrtrsError = TrtriError || interface.TransError;

// ---- trtri, trti2 ----------------------------------------------------------------

/// `?TRTRI`: `A^-1` of the `uplo` triangular `n x n` `A`, in place. A unit
/// `diag` is taken as ones and never read or written; the other triangle is
/// never read.
///
/// For a non-unit `diag`, a zero on the diagonal returns `error.Singular`
/// before `A` is touched.
pub fn trtri(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, diag: Diag, n: usize, a: []T, lda: usize) TrtriError!void {
    comptime ops.requireFloat(T, "trtri");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (n == 0) return;
    if (diag == .n and firstZero(T, n, a.ptr, lda) != null) return error.Singular;
    return inverse(allocator, io, T, uplo == .u, diag == .u, n, a.ptr, lda, interface.threads(.trtri, T, n));
}

/// `?TRTI2`: `trtri` unblocked (`trti2_U.c` and the rest, as loops). No
/// singularity check, as in the reference: a zero on the diagonal divides.
pub fn trti2(comptime T: type, uplo: Uplo, diag: Diag, n: usize, a: []T, lda: usize) ArgError!void {
    comptime ops.requireFloat(T, "trti2");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (n == 0) return;
    unblockedInverse(T, uplo == .u, diag == .u, n, a.ptr, lda);
}

/// The first zero on the diagonal of the `n x n` `a` (0-based), or null.
pub fn firstZero(comptime T: type, n: usize, a: [*]const T, lda: usize) ?usize {
    for (0..n) |i| if (ops.isZero(T, a[i * (lda + 1)])) return i;
    return null;
}

/// `trtri`'s algorithm: with `want` (`interface.threads`), the parallel one on
/// up to that many workers, else the single one.
pub fn inverse(allocator: Allocator, io: Io, comptime T: type, upper: bool, unit: bool, n: usize, a: [*]T, lda: usize, want: ?u32) driver.Error!void {
    const vt = interface.table(T);
    const blk = driver.Blocking.of(T, vt);
    if (want) |w| return inverseParallel(allocator, io, T, vt, upper, unit, interface.dtbOf(vt), blk.q, n, a, lda, w);
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
pub fn unblockedInverse(comptime T: type, upper: bool, unit: bool, n: usize, a: [*]T, lda: usize) void {
    for (0..n) |jj| {
        const j = if (upper) jj else n - 1 - jj;
        const dj = &a[j * (lda + 1)];
        const ajj: T = if (unit) scalar(T, 1) else recip(T, dj.*);
        if (!unit) dj.* = ajj;
        const x = a + j * lda;
        const s = ops.neg(T, ajj);
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
pub fn lauum(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) LauumError!void {
    comptime ops.requireFloat(T, "lauum");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (n == 0) return;
    return product(allocator, io, T, uplo == .u, n, a.ptr, lda, interface.threads(.lauum, T, n));
}

/// `?LAUU2`: `lauum` unblocked (`lauu2_U.c` and the rest, as loops).
pub fn lauu2(comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) ArgError!void {
    comptime ops.requireFloat(T, "lauu2");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (n == 0) return;
    unblockedProduct(T, uplo == .u, n, a.ptr, lda);
}

/// `lauum`'s algorithm: with `want` (`interface.threads`), the parallel one on
/// up to that many workers, else the single one.
pub fn product(allocator: Allocator, io: Io, comptime T: type, upper: bool, n: usize, a: [*]T, lda: usize, want: ?u32) driver.Error!void {
    const vt = interface.table(T);
    const blk = driver.Blocking.of(T, vt);
    if (want) |w| return productParallel(allocator, io, T, vt, upper, blk, n, a, lda, w);
    return productSingle(allocator, io, T, vt, upper, interface.dtbOf(vt), blk.q, n, a, lda);
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
    if (n <= 2 * blk.un) return productSingle(allocator, io, T, vt, upper, interface.dtbOf(vt), blk.q, n, a, lda);
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
pub fn unblockedProduct(comptime T: type, upper: bool, n: usize, a: [*]T, lda: usize) void {
    const F = Float(T);
    for (0..n) |i| {
        const d = &a[i * (lda + 1)];
        const aii = ops.re(T, d.*);
        // column `i` (`l`: row `i`) up to and with the diagonal
        const s: usize = if (upper) 1 else lda;
        const v = if (upper) a + i * lda else a + i;
        for (0..i + 1) |k| v[k * s] = ops.mulReal(T, v[k * s], aii);
        if (i + 1 == n) continue;
        // the rest of row `i` (`l`: column `i`)
        const rs: usize = if (upper) lda else 1;
        const r = if (upper) a + i + (i + 1) * lda else a + i + 1 + i * lda;
        const len = n - i - 1;
        var sum: F = 0;
        for (0..len) |k| sum += ops.abs2(T, r[k * rs]);
        d.* = scalar(T, ops.re(T, d.*) + sum);
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
/// before `B` is touched, even for `nrhs = 0`, as the reference.
pub fn trtrs(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: usize, nrhs: usize, a: []const T, lda: usize, b: []T, ldb: usize) TrtrsError!void {
    comptime ops.requireFloat(T, "trtrs");
    try interface.takes(trans);
    try interface.holdsTriangle(uplo, n, a.len, lda);
    try util.holds(b.len, n, nrhs, ldb);
    if (n == 0) return;
    if (diag == .n and firstZero(T, n, a.ptr, lda) != null) return error.Singular;
    if (nrhs == 0) return;
    const vt = interface.table(T);
    if (nrhs == 1) return triangular.trsvTable(allocator, io, T, vt, uplo, trans, diag, n, a.ptr, lda, b.ptr);
    const want = interface.threads(.trtrs, T, n) orelse 1;
    try driver.trsm(allocator, io, T, vt, true, uplo == .u, trans.transposed(), trans.conjA(), diag == .u, n, nrhs, scalar(T, 1), a.ptr, lda, b.ptr, ldb, want);
}
