//! ReLAPACK's recursive `sygst` and `hegst` (see `relapack.zig` for the source).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const driver = @import("../blas/blas3/driver.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const relapack = @import("relapack.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Side = types.Side;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const Error = driver.Error;
const crossover = relapack.crossover;
const ctrans = relapack.ctrans;
const split = relapack.split;
const symm = relapack.symm;
const syrk = relapack.syrk;
const Tview = relapack.Tview;
const trsm = relapack.trsm;

/// What `sygst` and `hegst` return: bad arguments in the reference's order,
/// `error.OutOfMemory` or `error.Canceled`.
/// The generalized eigenproblem `sygst` and `hegst` bring to standard form:
/// LAPACK's `itype` 1, 2 and 3.
pub const Problem = enum {
    /// `A x = λ B x`
    ax_lbx,
    /// `A B x = λ x`
    abx_lx,
    /// `B A x = λ x`
    bax_lx,
};

/// What `sygst` and `hegst` return: a reversed `uplo`, a slice too short or
/// a leading dimension too small, and the allocator's and cancellation.
pub const SygstError = interface.TriangleError || Allocator.Error || Io.Cancelable;

/// `?SYGST` by ReLAPACK's recursion (`relapack/src/dsygst.c`, `ssygst.c`;
/// OpenBLAS v0.3.34-317): the symmetric-definite `problem` to standard form,
/// `B` holding `potrf`'s factor in the same `uplo` triangle:
/// `A := inv(Uᵀ) A inv(U)` or `inv(L) A inv(Lᵀ)` for `ax_lbx`, `U A Uᵀ` or
/// `Lᵀ A L` for the other two, into `A`'s
/// `uplo` triangle. Real `T`; the other triangles are never read.
///
/// The base case is Reference-LAPACK's `xSYGS2` (`gs2`). The workspace
/// ReLAPACK takes from `malloc` (`n1 (n - n1)`, `XSYGST_ALLOW_MALLOC`) comes
/// from `allocator`; it holds `-1/2 B21 A11` (or its like) once, where
/// otherwise `symm` runs twice. Every level uses it, as `zhegst.c`'s `>=`
/// has it; `dsygst.c`, `ssygst.c` and `chegst.c` test `>`, which leaves the
/// top level's exactly sized buffer unused.
pub fn sygst(allocator: Allocator, io: Io, comptime T: type, problem: Problem, uplo: Uplo, n: usize, a: []T, lda: usize, b: []const T, ldb: usize) SygstError!void {
    comptime ops.requireFloat(T, "sygst");
    if (comptime isComplex(T)) @compileError("sygst: " ++ @typeName(T) ++ " is complex; hegst takes it");
    return gstImpl(allocator, io, T, problem, uplo, n, a, lda, b, ldb);
}

/// `?HEGST` (`relapack/src/zhegst.c`, `chegst.c`): `sygst` for the
/// Hermitian-definite problem, `ᴴ` for `ᵀ`; `A`'s diagonal comes out real.
/// Complex only.
pub fn hegst(allocator: Allocator, io: Io, comptime T: type, problem: Problem, uplo: Uplo, n: usize, a: []T, lda: usize, b: []const T, ldb: usize) SygstError!void {
    comptime ops.requireComplex(T, "hegst");
    return gstImpl(allocator, io, T, problem, uplo, n, a, lda, b, ldb);
}

fn gstImpl(al: Allocator, io: Io, comptime T: type, problem: Problem, uplo: Uplo, n: usize, a: []T, lda: usize, b: []const T, ldb: usize) SygstError!void {
    try interface.holdsTriangle(uplo, n, a.len, lda);
    try util.holds(b.len, n, n, ldb);
    if (n == 0) return;
    const n1 = split(T, n);
    const work = try al.alloc(T, n1 * (n - n1));
    defer al.free(work);
    try gstRec(al, io, T, problem == .ax_lbx, uplo, n, a.ptr, lda, b.ptr, ldb, work);
}

/// `RELAPACK_?sygst_rec`; `inv` for `ax_lbx`.
fn gstRec(al: Allocator, io: Io, comptime T: type, inv: bool, uplo: Uplo, n: usize, a: [*]T, lda: usize, b: [*]const T, ldb: usize, work: []T) Error!void {
    if (n <= crossover) return gs2(T, inv, uplo == .u, n, a, lda, b, ldb);
    const zero = scalar(T, 0);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const half = scalar(T, 0.5);
    const mhalf = scalar(T, -0.5);
    const ct = ctrans(T);
    const n1 = split(T, n);
    const n2 = n - n1;
    const a_tr = a + lda * n1;
    const a_bl = a + n1;
    const a_br = a + lda * n1 + n1;
    const b_tr = b + ldb * n1;
    const b_bl = b + n1;
    const b_br = b + ldb * n1 + n1;
    try gstRec(al, io, T, inv, uplo, n1, a, lda, b, ldb, work);
    // `T`, computed once into the workspace, or `symm` twice
    const use = work.len >= n2 * n1;
    const lower = uplo == .l;
    // the off-diagonal block, its rows and columns
    const x = if (lower) a_bl else a_tr;
    const xm = if (lower) n2 else n1;
    const xn = if (lower) n1 else n2;
    if (inv) {
        if (lower) {
            // A_BL = A_BL / B_TL'
            try trsm(al, io, T, true, .r, .l, ct, .n, n2, n1, one, b, ldb, a_bl, lda);
        } else {
            // A_TR = B_TL' \ A_TR
            try trsm(al, io, T, true, .l, .u, ct, .n, n1, n2, one, b, ldb, a_tr, lda);
        }
        // X -= 1/2 B_BL A_TL (`u`: A_TL B_TR), twice around the rank-2 update
        const side: Side = if (lower) .r else .l;
        const s = if (lower) b_bl else b_tr;
        if (use) {
            try symm(al, io, T, true, side, uplo, xm, xn, mhalf, a, lda, s, ldb, zero, work.ptr, xm);
            addBlock(T, xm, xn, work.ptr, xm, x, lda);
        } else try symm(al, io, T, true, side, uplo, xm, xn, mhalf, a, lda, s, ldb, one, x, lda);
        // A_BR -= A_BL B_BL' + B_BL A_BL' (`u`: A_TR' B_TR + B_TR' A_TR)
        try syrk(al, io, T, true, uplo, if (lower) .n else ct, n2, n1, minus, x, lda, s, ldb, one, a_br, lda);
        if (use) {
            addBlock(T, xm, xn, work.ptr, xm, x, lda);
        } else try symm(al, io, T, true, side, uplo, xm, xn, mhalf, a, lda, s, ldb, one, x, lda);
        if (lower) {
            // A_BL = B_BR \ A_BL
            try trsm(al, io, T, true, .l, .l, .n, .n, n2, n1, one, b_br, ldb, a_bl, lda);
        } else {
            // A_TR = A_TR / B_BR
            try trsm(al, io, T, true, .r, .u, .n, .n, n1, n2, one, b_br, ldb, a_tr, lda);
        }
    } else {
        if (lower) {
            // A_BL = A_BL B_TL
            try trsm(al, io, T, false, .r, .l, .n, .n, n2, n1, one, b, ldb, a_bl, lda);
        } else {
            // A_TR = B_TL A_TR
            try trsm(al, io, T, false, .l, .u, .n, .n, n1, n2, one, b, ldb, a_tr, lda);
        }
        // X += 1/2 A_BR B_BL (`u`: B_TR A_BR), twice around the rank-2 update
        const side: Side = if (lower) .l else .r;
        const s = if (lower) b_bl else b_tr;
        if (use) {
            try symm(al, io, T, true, side, uplo, xm, xn, half, a_br, lda, s, ldb, zero, work.ptr, xm);
            addBlock(T, xm, xn, work.ptr, xm, x, lda);
        } else try symm(al, io, T, true, side, uplo, xm, xn, half, a_br, lda, s, ldb, one, x, lda);
        // A_TL += A_BL' B_BL + B_BL' A_BL (`u`: A_TR B_TR' + B_TR A_TR')
        try syrk(al, io, T, true, uplo, if (lower) ct else .n, n1, n2, one, x, lda, s, ldb, one, a, lda);
        if (use) {
            addBlock(T, xm, xn, work.ptr, xm, x, lda);
        } else try symm(al, io, T, true, side, uplo, xm, xn, half, a_br, lda, s, ldb, one, x, lda);
        if (lower) {
            // A_BL = B_BR' A_BL
            try trsm(al, io, T, false, .l, .l, ct, .n, n2, n1, one, b_br, ldb, a_bl, lda);
        } else {
            // A_TR = A_TR B_BR'
            try trsm(al, io, T, false, .r, .u, ct, .n, n1, n2, one, b_br, ldb, a_tr, lda);
        }
    }
    try gstRec(al, io, T, inv, uplo, n2, a_br, lda, b_br, ldb, work);
}

/// `Y += X` on `m x n` (`?AXPY` with 1, a column at a time).
fn addBlock(comptime T: type, m: usize, n: usize, x: [*]const T, ldx: usize, y: [*]T, ldy: usize) void {
    for (0..n) |j| for (0..m) |i| {
        y[i + j * ldy] = ops.add(T, y[i + j * ldy], x[i + j * ldx]);
    };
}

/// Reference-LAPACK's `xSYGS2` and `xHEGS2` (`lapack-netlib/SRC/dsygs2.f`,
/// `zhegs2.f`; OpenBLAS v0.3.34-317), the unblocked `sygst`, with its
/// `syr2`/`her2`, `trsv`, `trmv`, `axpy` and `scal` as loops in the
/// reference BLAS's order. `upper` runs the lower algorithm on the
/// transposes, the conjugates of the reference's upper vectors (it conjugates
/// a row to get a column), so the same arithmetic up to their loop orders.
fn gs2(comptime T: type, inv: bool, upper: bool, n: usize, a_: [*]T, lda: usize, b_: [*]const T, ldb: usize) void {
    const herm = comptime isComplex(T);
    const a: Tview(T) = .{ .p = a_, .ld = lda, .t = upper };
    const b: Tview(T) = .{ .p = @constCast(b_), .ld = ldb, .t = upper };
    const cj = struct {
        fn f(v: T) T {
            return if (comptime herm) v.conjugate() else v;
        }
    }.f;
    const mul = ops.mul;
    const add = ops.add;
    for (0..n) |k| {
        const akk0 = ops.re(T, a.at(k, k).*);
        const bkk = ops.re(T, b.at(k, k).*);
        if (inv) {
            // inv(L) A inv(Lᴴ): the column below `k`, then the rest
            const akk = akk0 / (bkk * bkk);
            a.at(k, k).* = scalar(T, akk);
            if (k + 1 == n) continue;
            const r = 1 / bkk;
            for (k + 1..n) |i| a.at(i, k).* = ops.mulReal(T, a.at(i, k).*, r);
            const c = -0.5 * akk;
            for (k + 1..n) |i| a.at(i, k).* = add(T, a.at(i, k).*, ops.mulReal(T, b.at(i, k).*, c));
            rank2(T, a, k + 1, n, k, b, k, scalar(T, -1));
            for (k + 1..n) |i| a.at(i, k).* = add(T, a.at(i, k).*, ops.mulReal(T, b.at(i, k).*, c));
            // x := inv(B22) x, `xTRSV` `N`
            for (k + 1..n) |j| {
                const xj = a.at(j, k);
                if (ops.isZero(T, xj.*)) continue;
                xj.* = ops.div(T, xj.*, b.at(j, j).*);
                const t = xj.*;
                for (j + 1..n) |i| a.at(i, k).* = ops.sub(T, a.at(i, k).*, mul(T, t, b.at(i, j).*));
            }
        } else {
            // Lᴴ A L: row `k` left of the diagonal, as a conjugated column
            // x := B11ᴴ x, `xTRMV` `C`
            for (0..k) |j| {
                var t = cj(a.at(k, j).*);
                t = mul(T, t, cj(b.at(j, j).*));
                for (j + 1..k) |i| t = add(T, t, mul(T, cj(b.at(i, j).*), cj(a.at(k, i).*)));
                a.at(k, j).* = cj(t);
            }
            const c = 0.5 * akk0;
            for (0..k) |j| a.at(k, j).* = cj(add(T, cj(a.at(k, j).*), ops.mulReal(T, cj(b.at(k, j).*), c)));
            rank2Row(T, a, k, b, scalar(T, 1));
            for (0..k) |j| a.at(k, j).* = cj(ops.mulReal(T, add(T, cj(a.at(k, j).*), ops.mulReal(T, cj(b.at(k, j).*), c)), bkk));
            a.at(k, k).* = scalar(T, akk0 * (bkk * bkk));
        }
    }
}

/// `xSYR2` or `xHER2` (lower) with `alpha` on `A(lo:hi, lo:hi)`, `x` column
/// `xc` of `A` and `y` column `yc` of `B`, both from row `lo`.
fn rank2(comptime T: type, a: Tview(T), lo: usize, hi: usize, xc: usize, b: Tview(T), yc: usize, alpha: T) void {
    const herm = comptime isComplex(T);
    for (lo..hi) |j| {
        const xj = a.at(j, xc).*;
        const yj = b.at(j, yc).*;
        if (ops.isZero(T, xj) and ops.isZero(T, yj)) {
            if (comptime herm) a.at(j, j).* = scalar(T, a.at(j, j).re);
            continue;
        }
        const t1 = ops.mul(T, alpha, ops.conj(T, yj));
        const t2 = ops.conj(T, ops.mul(T, alpha, xj));
        if (comptime herm) {
            const d = ops.add(T, ops.mul(T, xj, t1), ops.mul(T, yj, t2));
            a.at(j, j).* = scalar(T, a.at(j, j).re + d.re);
        }
        for (if (herm) j + 1 else j..hi) |i| {
            a.at(i, j).* = ops.add(T, a.at(i, j).*, ops.add(T, ops.mul(T, a.at(i, xc).*, t1), ops.mul(T, b.at(i, yc).*, t2)));
        }
    }
}

/// `rank2` with `x` and `y` the conjugated rows `k` of `A` and `B` left of
/// the diagonal, on `A(:k, :k)`.
fn rank2Row(comptime T: type, a: Tview(T), k: usize, b: Tview(T), alpha: T) void {
    const herm = comptime isComplex(T);
    const cj = ops.conj;
    for (0..k) |j| {
        const xj = cj(T, a.at(k, j).*);
        const yj = cj(T, b.at(k, j).*);
        if (ops.isZero(T, xj) and ops.isZero(T, yj)) {
            if (comptime herm) a.at(j, j).* = scalar(T, a.at(j, j).re);
            continue;
        }
        const t1 = ops.mul(T, alpha, cj(T, yj));
        const t2 = cj(T, ops.mul(T, alpha, xj));
        if (comptime herm) {
            const d = ops.add(T, ops.mul(T, xj, t1), ops.mul(T, yj, t2));
            a.at(j, j).* = scalar(T, a.at(j, j).re + d.re);
        }
        for (if (herm) j + 1 else j..k) |i| {
            a.at(i, j).* = ops.add(T, a.at(i, j).*, ops.add(T, ops.mul(T, cj(T, a.at(k, i).*), t1), ops.mul(T, cj(T, b.at(k, i).*), t2)));
        }
    }
}
