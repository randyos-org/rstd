//! LU with complete pivoting, `P A Q = L U`, and its solve: faer v0.24.4
//! (`faer/src/linalg/lu/full_pivoting/factor.rs`, `lu_in_place_unblocked`
//! and `lu_in_place`; `solve.rs`; MIT). LAPACK has none but `?GETC2`,
//! which perturbs small pivots; this does not.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const threading = @import("../blas/threading.zig");
const laswp = @import("../blas/ext/laswp/kernel.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const lu = @import("lu.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;
const Float = types.Float;
const Mat = ops.Mat;

/// What `getrfComplete` returns: a slice too short or a leading dimension
/// too small, and cancellation. A rank below `min(m, n)` is no error.
pub const GetrfCompleteError = interface.ShapeError || Io.Cancelable;

/// `P A Q = L U` for the `m x n` `A` in `a`, `L` (unit diagonal implied)
/// below the diagonal and `U` on and above it, each pivot the largest
/// `|re| + |im|` of what is left. Row `i` was interchanged with row
/// `ipiv[i]`, column `i` with column `jpiv[i]` (0-based, `min(m, n)` each).
/// Returns the rank `r`: it stops where what is left is all below the
/// smallest normal number, its pivots then the identity and the trailing
/// block as the updates left it. The trailing updates run on up to `ger`'s
/// workers, by columns, with the same bits for any count.
pub fn getrfComplete(io: Io, comptime T: type, m: usize, n: usize, a: []T, lda: usize, ipiv: []usize, jpiv: []usize) GetrfCompleteError!usize {
    comptime ops.requireFloat(T, "getrfComplete");
    try util.holds(a.len, m, n, lda);
    const k = @min(m, n);
    if (ipiv.len < k or jpiv.len < k) return error.SliceTooShort;
    if (k == 0) return 0;
    return factor(io, T, .{ .ptr = a.ptr, .rows = m, .cols = n, .ld = lda }, ipiv[0..k], jpiv[0..k]);
}

/// The largest `|re| + |im|` found and where; the first of equals, in
/// column-major order (faer's strict `>`).
fn Best(comptime F: type) type {
    return struct {
        row: usize = 0,
        col: usize = 0,
        val: F = 0,

        fn better(x: @This(), y: @This()) @This() {
            return if (y.val > x.val) y else x;
        }
    };
}

/// `lu_in_place_unblocked`, column-major (`transpose` false).
fn factor(io: Io, comptime T: type, a: Mat(T), ipiv: []usize, jpiv: []usize) Io.Cancelable!usize {
    const F = Float(T);
    const kn = ipiv.len;
    var best: Best(F) = .{};
    for (0..a.cols) |j| best = best.better(scan(T, a, 0, j));
    for (0..kn) |k| {
        // faer's `max_score < min_positive()`
        if (best.val < std.math.floatMin(F)) {
            for (k..kn) |i| {
                ipiv[i] = i;
                jpiv[i] = i;
            }
            return k;
        }
        ipiv[k] = best.row;
        jpiv[k] = best.col;
        if (best.row != k) for (0..a.cols) |j| std.mem.swap(T, a.at(k, j), a.at(best.row, j));
        if (best.col != k) ops.swapCols(T, a, k, best.col);
        const inv = ops.recip(T, a.get(k, k));
        for (k + 1..a.rows) |i| a.at(i, k).* = ops.mul(T, a.get(i, k), inv);
        if (k + 1 == kn) break;
        best = try trailing(io, T, a, k);
    }
    return kn;
}

/// The best of column `j` from row `from` on.
fn scan(comptime T: type, a: Mat(T), from: usize, j: usize) Best(Float(T)) {
    var b: Best(Float(T)) = .{ .col = j };
    for (from..a.rows) |i| {
        const v = ops.abs1(T, a.get(i, j));
        if (v > b.val) b = .{ .row = i, .col = j, .val = v };
    }
    return b;
}

/// `A11 -= A10(:, k) A01(k, :)` below and right of `(k, k)`, and the best
/// of the result: faer's `rank_one_update_and_best_in_matrix`, its
/// threads by columns, folded in column order.
fn trailing(io: Io, comptime T: type, a: Mat(T), k: usize) Io.Cancelable!Best(Float(T)) {
    const Part = struct {
        a: Mat(T),
        k: usize,

        fn columns(s: @This(), c: threading.Chunk) Best(Float(T)) {
            var b: Best(Float(T)) = .{ .col = s.k + 1 + c.start };
            const l = s.a.col(s.k);
            for (s.k + 1 + c.start..s.k + 1 + c.end) |j| {
                const r = s.a.get(s.k, j);
                const col = s.a.col(j);
                for (s.k + 1..s.a.rows) |i| {
                    const x = ops.sub(T, col[i], ops.mul(T, l[i], r));
                    col[i] = x;
                    const v = ops.abs1(T, x);
                    if (v > b.val) b = .{ .row = i, .col = j, .val = v };
                }
            }
            return b;
        }
    };
    const cols = a.cols - k - 1;
    const work = threading.product(&.{ a.rows - k - 1, cols });
    var p = threading.plan(.ger, T, work, cols, .even, 1);
    p.workers = threading.share(p.workers, threading.matrix_grain, work *| @sizeOf(T));
    return threading.reduce(io, &p, Best(Float(T)), Part{ .a = a, .k = k }, Part.columns, Best(Float(T)).better);
}

/// What `getrsComplete` returns: a `trans` it does not take, a pivot
/// outside `0..n - 1` (checked before `B` is touched), a slice too short
/// or a leading dimension too small, and the allocator's and cancellation.
pub const GetrsCompleteError = interface.TransError || interface.PivotError || interface.ShapeError || Allocator.Error || Io.Cancelable;

/// Solves `op(A) X = B` for the `n x nrhs` `X`, in place of `B`, given
/// `getrfComplete`'s `P A Q = L U` of the square `A` (faer's
/// `solve_in_place_with_conj`, `solve_transpose_in_place_with_conj`):
/// `X = Q U⁻¹ L⁻¹ P B`, transposed `X = Pᵀ op(L)⁻¹ op(U)⁻¹ Qᵀ B`. `trans`
/// as `getrs`'s. No singularity check: past a rank below `n`, `U`'s
/// diagonal is tiny or zero.
pub fn getrsComplete(allocator: Allocator, io: Io, comptime T: type, trans: Trans, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const usize, jpiv: []const usize, b: []T, ldb: usize) GetrsCompleteError!void {
    comptime ops.requireFloat(T, "getrsComplete");
    try interface.takes(trans);
    try util.holds(a.len, n, n, lda);
    try util.holds(b.len, n, nrhs, ldb);
    try interface.pivotsIn(ipiv, n);
    try interface.pivotsIn(jpiv, n);
    if (n == 0 or nrhs == 0) return;
    const want = interface.threads(.getrs, T, n) orelse 1;
    // `lu.solve` takes `P` (`ipiv`) on the side `getrs` does
    const t = trans.transposed();
    if (t) laswp.swapRows(T, nrhs, b.ptr, ldb, 0, n, jpiv.ptr, false);
    try lu.solve(allocator, io, T, trans, n, nrhs, a.ptr, lda, ipiv.ptr, b.ptr, ldb, want);
    if (!t) laswp.swapRows(T, nrhs, b.ptr, ldb, 0, n, jpiv.ptr, true);
}
