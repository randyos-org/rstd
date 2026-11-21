//! The sparse tests' helpers: random matrices beside their dense twins, and
//! dense references in plain loops.

const std = @import("std");
const matrix = @import("../../../../linalg/matrix.zig");
const types = @import("../../blas/types.zig");
const ops = @import("../../lapack/ops.zig");
const check = @import("../../blas/_test/check.zig");
const sparse = @import("../__root__.zig");

const Allocator = std.mem.Allocator;
const Trans = types.Trans;

/// The element types every test covers.
pub const all_types = check.sd_types;

/// A matrix as CSC and dense.
pub fn Pair(comptime T: type) type {
    return struct {
        s: sparse.CompressedColumns(T),
        d: matrix.Mat(T),

        pub fn deinit(p: *@This(), al: Allocator) void {
            p.s.deinit(al);
            p.d.deinit(al);
        }
    };
}

/// `rows` by `cols`, about `per_mille` of 1000 positions set, in random
/// order with repeats, each repeat summed as `toCompressedColumns` sums it.
pub fn random(comptime T: type, al: Allocator, rows: usize, cols: usize, per_mille: u32, seed: u64) !Pair(T) {
    var prng: std.Random.DefaultPrng = .init(seed);
    const r = prng.random();
    var t = try sparse.Triplets(T).init(rows, cols);
    defer t.deinit(al);
    const d: matrix.Mat(T) = try .zeros(al, rows, cols);
    errdefer d.deinit(al);
    const count = rows * cols * per_mille / 1000;
    var v: [1]T = undefined;
    for (0..count) |k| {
        const i = r.uintLessThan(usize, rows);
        const j = r.uintLessThan(usize, cols);
        check.fill(T, &v, seed *% 7919 +% k);
        try t.append(al, i, j, v[0]);
        const at = &d.data[i + j * rows];
        at.* = ops.add(T, at.*, v[0]);
    }
    return .{ .s = try t.toCompressedColumns(al), .d = d };
}

/// Whether `a` holds the CSC invariants: column starts rising to `nnz`,
/// rows ascending and in range in each column.
pub fn valid(comptime T: type, a: sparse.CompressedColumns(T)) bool {
    if (a.col_start.len != @as(usize, a.cols) + 1 or a.col_start[0] != 0 or a.col_start[a.cols] != a.nonzeros() or a.row.len != a.value.len) return false;
    for (0..a.cols) |j| {
        if (a.col_start[j] > a.col_start[j + 1]) return false;
        const rows = a.row[a.col_start[j]..a.col_start[j + 1]];
        for (rows, 0..) |i, k| if (i >= a.rows or (k > 0 and rows[k - 1] >= i)) return false;
    }
    return true;
}

/// `a` as a fresh dense matrix, through `toDense`.
pub fn dense(comptime T: type, al: Allocator, a: sparse.CompressedColumns(T)) !matrix.Mat(T) {
    const m: matrix.Mat(T) = try .init(al, a.rows, a.cols);
    errdefer m.deinit(al);
    try a.toDense(m.data, @max(a.rows, 1));
    return m;
}

/// `got`'s dense form against `want`, exactly.
pub fn expectDense(comptime T: type, al: Allocator, want: matrix.Mat(T), got: sparse.CompressedColumns(T)) !void {
    try std.testing.expect(valid(T, got));
    const g = try dense(T, al, got);
    defer g.deinit(al);
    try std.testing.expectEqual(want.rows, g.rows);
    try std.testing.expectEqual(want.cols, g.cols);
    try std.testing.expectEqualSlices(u8, std.mem.sliceAsBytes(want.data), std.mem.sliceAsBytes(g.data));
}

fn opA(comptime T: type, t: Trans, v: T) T {
    return if (t.conjA()) ops.conj(T, v) else v;
}

/// `y := alpha op(A) x + beta y` over dense `a`, in plain loops.
pub fn mulVec(comptime T: type, a: matrix.Mat(T), t: Trans, alpha: T, x: []const T, beta: T, y: []T) void {
    const tr = t.transposed();
    const ny = if (tr) a.cols else a.rows;
    const nx = if (tr) a.rows else a.cols;
    for (0..ny) |i| {
        var s = ops.zero(T);
        for (0..nx) |j| {
            const e = if (tr) a.at(j, i) else a.at(i, j);
            const xj = if (t.conjX()) ops.conj(T, x[j]) else x[j];
            s = ops.add(T, s, ops.mul(T, opA(T, t, e), xj));
        }
        const as = ops.mul(T, alpha, s);
        y[i] = if (ops.isZero(T, beta)) as else ops.add(T, as, ops.mul(T, beta, y[i]));
    }
}

/// `A B`, dense, in plain loops.
pub fn mulDense(comptime T: type, al: Allocator, a: matrix.Mat(T), b: matrix.Mat(T)) !matrix.Mat(T) {
    const c: matrix.Mat(T) = try .zeros(al, a.rows, b.cols);
    for (0..b.cols) |j| for (0..a.rows) |i| {
        var s = ops.zero(T);
        for (0..a.cols) |k| s = ops.add(T, s, ops.mul(T, a.at(i, k), b.at(k, j)));
        c.data[i + j * a.rows] = s;
    };
    return c;
}

/// The tolerance for sums of `terms` products of values in [-2, 2), scaled
/// by `scale`.
pub fn tol(comptime T: type, terms: usize, scale: f64) types.Float(T) {
    return check.tolerance(types.Float(T), terms) * @as(types.Float(T), @floatCast(scale));
}

/// `n` random values from `seed`.
pub fn vec(comptime T: type, al: Allocator, n: usize, seed: u64) ![]T {
    const v = try al.alloc(T, n);
    check.fill(T, v, seed);
    return v;
}
