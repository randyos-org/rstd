//! `solveTri` against dense products: `op(A) x` solved back to `x`.

const std = @import("std");
const matrix = @import("../../matrix.zig");
const types = @import("../../blas/types.zig");
const ops = @import("../../lapack/ops.zig");
const check = @import("../../blas/_test/check.zig");
const sparse = @import("../__root__.zig");
const util = @import("util.zig");

const testing = std.testing;
const Trans = types.Trans;
const Uplo = types.Uplo;
const Diag = types.Diag;

/// A random `n` by `n` matrix, entries below 1/2 in magnitude off the
/// diagonal and at least 4 on it, so its triangles are well conditioned,
/// unit diagonal or not; the other triangle filled too.
fn diagonallyHeavy(comptime T: type, al: std.mem.Allocator, n: usize, seed: u64) !util.Pair(T) {
    var p = try util.random(T, al, n, n, 300, seed);
    errdefer p.deinit(al);
    // exact: a power of two
    const quarter = check.val(T, 0.25, 0);
    p.s.scal(quarter);
    for (p.d.data) |*v| v.* = ops.mul(T, quarter, v.*);
    var t = try sparse.Triplets(T).init(n, n);
    defer t.deinit(al);
    for (0..n) |j| for (p.s.col_start[j]..p.s.col_start[j + 1]) |k| {
        try t.append(al, p.s.row[k], j, p.s.value[k]);
    };
    for (0..n) |i| {
        const d = check.val(T, if (i % 2 == 0) 4 else -5, 1);
        try t.append(al, i, i, d);
        p.d.data[i + i * n] = ops.add(T, p.d.data[i + i * n], d);
    }
    p.s.deinit(al);
    p.s = try t.toCsc(al);
    return p;
}

/// The triangle `uplo` of dense `a`, its diagonal ones for a unit `diag`.
fn triangle(comptime T: type, al: std.mem.Allocator, a: matrix.Mat(T), uplo: Uplo, diag: Diag) !matrix.Mat(T) {
    const n = a.rows;
    const t: matrix.Mat(T) = try .zeros(al, n, n);
    for (0..n) |j| for (0..n) |i| {
        const in = if (uplo == .l) i >= j else i <= j;
        if (in) t.data[i + j * n] = if (i == j and diag == .u) ops.one(T) else a.at(i, j);
    };
    return t;
}

test "solveTri: lower and upper, every trans, unit or not" {
    const al = testing.allocator;
    inline for (util.all_types) |T| {
        for ([_]usize{ 1, 2, 9, 40 }, 0..) |n, k| {
            var p = try diagonallyHeavy(T, al, n, k + 90);
            defer p.deinit(al);
            const x0 = try util.vec(T, al, n + 1, k);
            defer al.free(x0);
            const x = try al.alloc(T, n + 1);
            defer al.free(x);
            for ([_]Uplo{ .l, .u }) |uplo| for ([_]Diag{ .n, .u }) |diag| {
                const tri = try triangle(T, al, p.d, uplo, diag);
                defer tri.deinit(al);
                for ([_]Trans{ .n, .t, .r, .c }) |t| {
                    // b := op(tri) x0, then solved in place
                    @memcpy(x, x0);
                    util.mulVec(T, tri, t, ops.one(T), x0, ops.zero(T), x);
                    try p.s.solveTri(uplo, t, diag, x);
                    check.expectCloseSlices(T, x0[0..n], x[0..n], util.tol(T, n, 4)) catch |e| {
                        check.say("{s} n {d} uplo {s} diag {s} trans {s}\n", .{ @typeName(T), n, @tagName(uplo), @tagName(diag), @tagName(t) });
                        return e;
                    };
                    try testing.expectEqual(x0[n], x[n]);
                }
            };
        }
    }
}

test "solveTri: singular, argument and dimension checks" {
    const al = testing.allocator;
    var t = try sparse.Triplets(f64).init(3, 3);
    defer t.deinit(al);
    try t.append(al, 0, 0, 2);
    try t.append(al, 1, 0, 1);
    try t.append(al, 2, 2, 0);
    try t.append(al, 0, 2, 7);
    var a = try t.toCsc(al);
    defer a.deinit(al);
    var x = [_]f64{ 4, 5, 6 };
    // (1, 1) missing, (2, 2) a stored zero: x untouched
    try testing.expectError(error.Singular, a.solveTri(.l, .n, .n, &x));
    try testing.expectError(error.Singular, a.solveTri(.u, .t, .n, &x));
    try testing.expectEqualSlices(f64, &.{ 4, 5, 6 }, &x);
    // a unit diagonal reads neither
    try a.solveTri(.l, .n, .u, &x);
    try testing.expectEqualSlices(f64, &.{ 4, 1, 6 }, &x);
    try a.solveTri(.u, .n, .u, &x);
    try testing.expectEqualSlices(f64, &.{ -38, 1, 6 }, &x);
    try testing.expectError(error.InvalidUplo, a.solveTri(.m, .n, .u, &x));
    try testing.expectError(error.InvalidTrans, a.solveTri(.l, .o, .u, &x));
    try testing.expectError(error.DimensionMismatch, a.solveTri(.l, .n, .u, x[0..2]));
    var r = try sparse.Triplets(f64).init(3, 2);
    defer r.deinit(al);
    var b = try r.toCsc(al);
    defer b.deinit(al);
    try testing.expectError(error.DimensionMismatch, b.solveTri(.l, .n, .u, &x));
    var z = try sparse.Triplets(types.Complex32).init(0, 0);
    defer z.deinit(al);
    var e = try z.toCsc(al);
    defer e.deinit(al);
    var none: [0]types.Complex32 = .{};
    try e.solveTri(.u, .c, .n, &none);
}
