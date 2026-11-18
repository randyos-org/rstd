//! Building and taking apart: `Triplets`, dense conversions, `scal`,
//! columns as sparse vectors, extraction, transposes and the CSR twin.

const std = @import("std");
const matrix = @import("../../matrix.zig");
const types = @import("../../blas/types.zig");
const ops = @import("../../lapack/ops.zig");
const check = @import("../../blas/_test/check.zig");
const sparse = @import("../__root__.zig");
const util = @import("util.zig");

const testing = std.testing;
const Mat = matrix.Mat;

test "Triplets: duplicates sum in insertion order, rows ascend" {
    const al = testing.allocator;
    var t = try sparse.Triplets(f64).init(4, 3);
    defer t.deinit(al);
    try t.append(al, 3, 1, 5);
    try t.append(al, 0, 1, 1e16);
    try t.append(al, 2, 0, 7);
    try t.append(al, 0, 1, 1);
    try t.append(al, 1, 1, 6);
    try t.append(al, 0, 1, -1e16);
    try t.append(al, 2, 2, 3);
    try t.append(al, 2, 2, -3);
    var a = try t.toCsc(al);
    defer a.deinit(al);
    try testing.expect(util.valid(f64, a));
    try testing.expectEqualSlices(u32, &.{ 0, 1, 4, 5 }, a.col_start);
    try testing.expectEqualSlices(u32, &.{ 2, 0, 1, 3, 2 }, a.row);
    // (1e16 + 1) - 1e16 rounds to 0, as inserted; a zero sum stays stored
    try testing.expectEqualSlices(f64, &.{ 7, 0, 6, 5, 0 }, a.value);
}

test "Triplets: index checks and empty shapes" {
    const al = testing.allocator;
    var t = try sparse.Triplets(f64).init(3, 2);
    defer t.deinit(al);
    try t.append(al, 2, 1, 4);
    try testing.expectError(error.IndexOutOfRange, t.append(al, 3, 0, 1));
    try testing.expectError(error.IndexOutOfRange, t.append(al, 0, 2, 1));
    try testing.expectEqual(@as(usize, 1), t.x.items.len);
    if (comptime @bitSizeOf(usize) > 32) {
        try testing.expectError(error.TooLarge, sparse.Triplets(f64).init(@as(usize, std.math.maxInt(u32)) + 1, 1));
    }
    inline for (.{ .{ 0, 0 }, .{ 0, 3 }, .{ 3, 0 }, .{ 2, 2 } }) |s| {
        var e = try sparse.Triplets(types.Complex64).init(s[0], s[1]);
        defer e.deinit(al);
        var a = try e.toCsc(al);
        defer a.deinit(al);
        try testing.expect(util.valid(types.Complex64, a));
        try testing.expectEqual(@as(usize, 0), a.nonzeros());
        const d = try a.toDense(al);
        defer d.deinit(al);
        try testing.expectEqual(@as(usize, s[0] * s[1]), d.data.len);
    }
}

test "random matrices match their dense twins" {
    const al = testing.allocator;
    inline for (util.all_types) |T| {
        for ([_][2]usize{ .{ 1, 1 }, .{ 7, 3 }, .{ 20, 31 }, .{ 64, 9 } }, 0..) |s, n| {
            var p = try util.random(T, al, s[0], s[1], 300, n);
            defer p.deinit(al);
            try util.expectDense(T, al, p.d, p.s);
        }
    }
}

test "fromDense drops magnitudes at most drop_below; toDense restores" {
    const al = testing.allocator;
    inline for (util.all_types) |T| {
        const F = types.Float(T);
        const vals = [_][2]f64{ .{ 0, 0 }, .{ 1, 0 }, .{ -0.5, 0 }, .{ 0, 0.25 }, .{ 3, 4 }, .{ 0.2, 0.2 }, .{ 0.1, 0 }, .{ -2, 0 } };
        var data: [8]T = undefined;
        for (&data, vals) |*d, v| d.* = check.val(T, v[0], if (comptime types.isComplex(T)) v[1] else 0);
        // 3 by 2 at ld 4: row 3 is not the matrix
        const view = try matrix.ConstView(T).init(&data, 3, 2, 4);
        var a = try sparse.Csc(T).fromDense(al, view, 0);
        defer a.deinit(al);
        try testing.expect(util.valid(T, a));
        try testing.expectEqualSlices(u32, &.{ 0, 2, 5 }, a.col_start);
        try testing.expectEqualSlices(u32, &.{ 1, 2, 0, 1, 2 }, a.row);
        const d = try a.toDense(al);
        defer d.deinit(al);
        for (0..3) |i| for (0..2) |j| try testing.expectEqual(view.at(i, j), d.at(i, j));
        // |0.1| and real 0.2 at most 0.25, |0.2 + 0.2i| past it
        var b = try sparse.Csc(T).fromDense(al, view, @as(F, 0.25));
        defer b.deinit(al);
        try testing.expectEqualSlices(u32, if (comptime types.isComplex(T)) &.{ 1, 2, 0, 1 } else &.{ 1, 2, 0 }, b.row);
        // a negative drop_below keeps every entry, zeros too
        var c = try sparse.Csc(T).fromDense(al, view, -1);
        defer c.deinit(al);
        try testing.expectEqual(@as(usize, 6), c.nonzeros());
    }
    if (comptime @bitSizeOf(usize) > 32) {
        const big: matrix.ConstView(f64) = .{ .data = &.{}, .rows = @as(usize, std.math.maxInt(u32)) + 1, .cols = 0, .ld = @as(usize, std.math.maxInt(u32)) + 1 };
        try testing.expectError(error.TooLarge, sparse.Csc(f64).fromDense(al, big, 0));
    }
    var e = try sparse.Csc(f32).fromDense(al, .{ .data = &.{}, .rows = 0, .cols = 4, .ld = 1 }, 0);
    defer e.deinit(al);
    try testing.expect(util.valid(f32, e));
}

test "scal, and columns as sparse vectors: axpy, dot, dotc" {
    const al = testing.allocator;
    inline for (util.all_types) |T| {
        var p = try util.random(T, al, 40, 6, 250, 11);
        defer p.deinit(al);
        const alpha = check.val(T, 0.75, -1.25);
        const y0 = try util.vec(T, al, 41, 5);
        defer al.free(y0);
        for (0..6) |j| {
            const col = p.s.col(j);
            try testing.expectEqual(@as(usize, 40), col.len);
            // axpy against the dense column, over 40 of y's 41 elements
            const y = try al.dupe(T, y0);
            defer al.free(y);
            try col.axpy(alpha, y);
            for (0..40) |i| try check.expectClose(T, ops.add(T, y0[i], ops.mul(T, alpha, p.d.at(i, j))), y[i], util.tol(T, 2, 1));
            try testing.expectEqual(y0[40], y[40]);
            var d = ops.zero(T);
            var dc = ops.zero(T);
            for (0..40) |i| {
                d = ops.add(T, d, ops.mul(T, p.d.at(i, j), y0[i]));
                dc = ops.add(T, dc, ops.mulConj(T, p.d.at(i, j), y0[i]));
            }
            try check.expectClose(T, d, try col.dot(y0), util.tol(T, 40, 1));
            try check.expectClose(T, dc, try col.dotc(y0), util.tol(T, 40, 1));
            try testing.expectError(error.DimensionMismatch, col.dot(y0[0..39]));
            try testing.expectError(error.DimensionMismatch, col.axpy(alpha, y[0..39]));
        }
        p.s.scal(alpha);
        for (p.d.data) |*v| v.* = ops.mul(T, alpha, v.*);
        try util.expectDense(T, al, p.d, p.s);
    }
}

test "sub and selectCols copy the parts asked for" {
    const al = testing.allocator;
    inline for (util.all_types) |T| {
        var p = try util.random(T, al, 13, 11, 400, 3);
        defer p.deinit(al);
        for ([_][4]usize{ .{ 0, 0, 13, 11 }, .{ 2, 3, 5, 6 }, .{ 12, 10, 1, 1 }, .{ 4, 4, 0, 3 }, .{ 13, 11, 0, 0 } }) |s| {
            var b = try p.s.sub(al, s[0], s[1], s[2], s[3]);
            defer b.deinit(al);
            const want: Mat(T) = try .zeros(al, s[2], s[3]);
            defer want.deinit(al);
            for (0..s[2]) |i| for (0..s[3]) |j| {
                want.data[i + j * s[2]] = p.d.at(s[0] + i, s[1] + j);
            };
            try util.expectDense(T, al, want, b);
        }
        try testing.expectError(error.DimensionMismatch, p.s.sub(al, 1, 0, 13, 1));
        try testing.expectError(error.DimensionMismatch, p.s.sub(al, 0, 11, 1, 1));
        try testing.expectError(error.DimensionMismatch, p.s.sub(al, std.math.maxInt(usize), 0, 2, 1));
        const pick = [_]u32{ 10, 0, 4, 4, 7 };
        var c = try p.s.selectCols(al, &pick);
        defer c.deinit(al);
        const want: Mat(T) = try .init(al, 13, pick.len);
        defer want.deinit(al);
        for (pick, 0..) |q, j| @memcpy(want.data[j * 13 ..][0..13], p.d.data[q * 13 ..][0..13]);
        try util.expectDense(T, al, want, c);
        try testing.expectError(error.IndexOutOfRange, p.s.selectCols(al, &.{ 1, 11 }));
        var none = try p.s.selectCols(al, &.{});
        defer none.deinit(al);
        try testing.expectEqual(@as(usize, 0), none.cols);
    }
}

test "transpose, ctranspose and the CSR twin" {
    const al = testing.allocator;
    inline for (util.all_types) |T| {
        var p = try util.random(T, al, 9, 14, 300, 8);
        defer p.deinit(al);
        var t = try p.s.transpose(al);
        defer t.deinit(al);
        var h = try p.s.ctranspose(al);
        defer h.deinit(al);
        const want_t: Mat(T) = try .init(al, 14, 9);
        defer want_t.deinit(al);
        const want_h: Mat(T) = try .init(al, 14, 9);
        defer want_h.deinit(al);
        for (0..9) |i| for (0..14) |j| {
            want_t.data[j + i * 14] = p.d.at(i, j);
            // a zero's conjugate is a negative zero; nothing stored there
            const v = p.d.at(i, j);
            want_h.data[j + i * 14] = if (ops.isZero(T, v)) v else ops.conj(T, v);
        };
        try util.expectDense(T, al, want_t, t);
        try util.expectDense(T, al, want_h, h);
        // the CSR of A: row i as a sparse vector is A's row
        var r = try p.s.toCsr(al);
        defer r.deinit(al);
        try testing.expectEqual(@as(usize, 9), r.rows);
        for (0..9) |i| {
            const row = r.row(i);
            var k: usize = 0;
            for (0..14) |j| {
                if (k < row.index.len and row.index[k] == j) {
                    try testing.expectEqual(p.d.at(i, j), row.value[k]);
                    k += 1;
                } else try testing.expect(ops.isZero(T, p.d.at(i, j)));
            }
        }
        var back = try r.toCsc(al);
        defer back.deinit(al);
        try util.expectDense(T, al, p.d, back);
        // a CSC of A is a CSR of Aᵀ, the same arrays
        const v = p.s.asCsrOfTranspose();
        try testing.expectEqual(p.s.row.ptr, v.col.ptr);
        try testing.expectEqual(p.s.cols, v.rows);
        try testing.expectEqual(p.s.row.ptr, v.asCscOfTranspose().row.ptr);
    }
}
