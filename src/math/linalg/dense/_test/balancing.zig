//! Tests: `balance`.

const std = @import("std");
const matrix = @import("../../matrix.zig");
const vector = @import("../../vector.zig");
const ops = @import("../../../kernels/lapack/ops.zig");
const balancing = @import("../balancing.zig");

const testing = std.testing;
const Mat = matrix.Mat;
const Real = vector.Real;
const C64 = std.math.Complex(f64);
const test_io = std.Io.failing;
const balance = balancing.balance;

/// `B(p[i], p[j]) = A(i, j) s[j] / s[i]` exactly, `s` powers of 2, `p` a
/// permutation, `T(i, p[i]) = s[i]`.
fn expectBalanced(comptime T: type, a: Mat(T), f: balancing.Balance(T)) !void {
    const al = testing.allocator;
    const n = a.rows;
    const seen = try al.alloc(bool, n);
    defer al.free(seen);
    @memset(seen, false);
    for (f.p, f.s) |c, x| {
        seen[c] = true;
        try testing.expectEqual(@as(Real(T), 1), std.math.frexp(x).significand * 2);
    }
    for (seen) |s| try testing.expect(s);
    for (0..n) |j| for (0..n) |i| {
        const want = ops.divReal(T, ops.mulReal(T, a.at(i, j), f.s[j]), f.s[i]);
        try testing.expectEqual(want, f.b.at(f.p[i], f.p[j]));
    };
    const t = try f.t(al);
    defer t.deinit(al);
    for (0..n) |i| try testing.expectEqual(ops.scalar(T, f.s[i]), t.at(i, f.p[i]));
}

fn badlyScaled(comptime T: type, n: usize) !Mat(T) {
    const a = try Mat(T).init(testing.allocator, n, n);
    for (0..n) |j| for (0..n) |i| {
        const e: Real(T) = @floatFromInt(@as(i64, @intCast(j)) - @as(i64, @intCast(i)));
        const v = std.math.pow(Real(T), 16, 2 * e) * (1 + @as(Real(T), @floatFromInt((i * 7 + j * 3) % 5)) / 8);
        a.data[i + j * n] = ops.scalar(T, v);
    };
    return a;
}

test "balance: B = T⁻¹ A T exactly, scalings powers of 2, a badly scaled matrix evened out" {
    const al = testing.allocator;
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ 1, 4, 9 }) |n| for ([_]balancing.BalanceMode{ .both, .noperm }) |mode| {
        const a = try badlyScaled(T, n);
        defer a.deinit(al);
        const f = try balance(al, test_io, a, mode);
        defer f.deinit(al);
        try expectBalanced(T, a, f);
        if (mode == .noperm) for (f.p, 0..) |c, i| try testing.expectEqual(i, c);
        // the corners' 16^(2 (n - 1)) spread shrinks
        if (n > 1) {
            const before = ops.abs1(T, a.at(0, n - 1)) / ops.abs1(T, a.at(n - 1, 0));
            const after = ops.abs1(T, f.b.at(f.p[0], f.p[n - 1])) / ops.abs1(T, f.b.at(f.p[n - 1], f.p[0]));
            try testing.expect(after < before);
        }
    };
}

test "balance: an isolated eigenvalue is permuted out; errors" {
    const al = testing.allocator;
    // row 1 is zero off the diagonal: its eigenvalue, 5, is isolated
    const a = try Mat(f64).init(al, 3, 3);
    defer a.deinit(al);
    @memcpy(a.data, &[_]f64{ 1, 0, 3, 2, 5, 4, 6, 0, 7 });
    const f = try balance(al, test_io, a, .both);
    defer f.deinit(al);
    try expectBalanced(f64, a, f);
    try testing.expectEqual(@as(usize, 2), f.p[1]);
    const wide = try Mat(f64).zeros(al, 2, 3);
    defer wide.deinit(al);
    try testing.expectError(error.NotSquare, balance(al, test_io, wide, .both));
    // a NaN in the block left to scale
    a.data[0] = std.math.nan(f64);
    try testing.expectError(error.NanInA, balance(al, test_io, a, .both));
}
