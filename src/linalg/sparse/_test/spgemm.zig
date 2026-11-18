//! `spgemm`, its symbolic and numeric passes, against dense loops.

const std = @import("std");
const matrix = @import("../../matrix.zig");
const types = @import("../../blas/types.zig");
const ops = @import("../../lapack/ops.zig");
const check = @import("../../blas/_test/check.zig");
const sparse = @import("../__root__.zig");
const util = @import("util.zig");

const testing = std.testing;

test "spgemm against dense loops" {
    const al = testing.allocator;
    inline for (util.all_types) |T| {
        for ([_][3]usize{ .{ 1, 1, 1 }, .{ 9, 6, 11 }, .{ 20, 30, 7 }, .{ 0, 5, 3 }, .{ 4, 0, 6 }, .{ 5, 7, 0 } }, 0..) |s, k| {
            var a = try util.random(T, al, s[0], s[1], 250, k + 60);
            defer a.deinit(al);
            var b = try util.random(T, al, s[1], s[2], 250, k + 70);
            defer b.deinit(al);
            var c = try a.s.spgemm(al, b.s);
            defer c.deinit(al);
            try testing.expect(util.valid(T, c));
            try testing.expectEqual(a.s.rows, c.rows);
            try testing.expectEqual(b.s.cols, c.cols);
            const want = try util.mulDense(T, al, a.d, b.d);
            defer want.deinit(al);
            const got = try c.toDense(al);
            defer got.deinit(al);
            try check.expectCloseSlices(T, want.data, got.data, util.tol(T, s[1], 1));
            // the pattern is the Boolean product's: no entry without a path
            for (0..c.cols) |j| for (c.col_start[j]..c.col_start[j + 1]) |q| {
                const i = c.row[q];
                var path = false;
                for (b.s.row[b.s.col_start[j]..b.s.col_start[j + 1]]) |kk| {
                    for (a.s.row[a.s.col_start[kk]..a.s.col_start[kk + 1]]) |r| path = path or r == i;
                }
                try testing.expect(path);
            };
        }
    }
}

test "spgemm: the numeric pass reuses the symbolic pattern" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| {
        var a = try util.random(T, al, 15, 12, 300, 1);
        defer a.deinit(al);
        var b = try util.random(T, al, 12, 10, 300, 2);
        defer b.deinit(al);
        var c = try a.s.spgemmSymbolic(al, b.s);
        defer c.deinit(al);
        for (c.value) |v| try testing.expect(ops.isZero(T, v));
        // new values, same patterns
        for (a.s.value, 0..) |*v, k| v.* = check.val(T, @floatFromInt(k % 5), 1);
        for (b.s.value, 0..) |*v, k| v.* = check.val(T, @floatFromInt(k % 3), -1);
        const ad = try a.s.toDense(al);
        defer ad.deinit(al);
        const bd = try b.s.toDense(al);
        defer bd.deinit(al);
        try a.s.spgemmNumeric(al, b.s, c);
        const want = try util.mulDense(T, al, ad, bd);
        defer want.deinit(al);
        try util.expectDense(T, al, want, c);
        // again: the slots start from zero each time
        try a.s.spgemmNumeric(al, b.s, c);
        try util.expectDense(T, al, want, c);
    }
}

test "spgemm: dimension and pattern checks" {
    const al = testing.allocator;
    var a = try util.random(f64, al, 6, 5, 400, 3);
    defer a.deinit(al);
    var b = try util.random(f64, al, 5, 4, 400, 4);
    defer b.deinit(al);
    try testing.expectError(error.DimensionMismatch, b.s.spgemm(al, a.s));
    try testing.expectError(error.DimensionMismatch, b.s.spgemmSymbolic(al, b.s));
    var c = try a.s.spgemmSymbolic(al, b.s);
    defer c.deinit(al);
    try testing.expectError(error.DimensionMismatch, a.s.spgemmNumeric(al, a.s, c));
    // a pattern with an entry `c` lacks
    var t = try sparse.Triplets(f64).init(6, 5);
    defer t.deinit(al);
    for (0..6) |i| for (0..5) |j| try t.append(al, i, j, 1);
    var full = try t.toCsc(al);
    defer full.deinit(al);
    var bt = try sparse.Triplets(f64).init(5, 4);
    defer bt.deinit(al);
    for (0..4) |j| try bt.append(al, j, j, 1);
    var diag = try bt.toCsc(al);
    defer diag.deinit(al);
    var cd = try a.s.spgemmSymbolic(al, diag);
    defer cd.deinit(al);
    if (cd.nonzeros() < 24) try testing.expectError(error.PatternMismatch, full.spgemmNumeric(al, diag, cd));
}
