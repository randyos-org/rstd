//! Tests: `lasq1` (dqds) against `bdsqr`'s QR iteration.

const std = @import("std");
const check = @import("../../blas/_test/check.zig");
const lasq = @import("../lasq.zig");
const bdsqr = @import("../bdsqr.zig");

const testing = std.testing;

/// `d` and `e` (`n - 1`) by dqds and by the QR iteration, each value within
/// `tol` of the other relatively (absolutely below `floor`).
fn agree(comptime F: type, d0: []const F, e0: []const F, tol: F) !void {
    const al = testing.allocator;
    const n = d0.len;
    const d = try al.dupe(F, d0);
    defer al.free(d);
    const e = try al.dupe(F, e0);
    defer al.free(e);
    const work = try al.alloc(F, 4 * n);
    defer al.free(work);
    try testing.expectEqual(lasq.Outcome.converged, lasq.lasq1(F, d, e, work));
    const qd = try al.dupe(F, d0);
    defer al.free(qd);
    const qe = try al.alloc(F, n);
    defer al.free(qe);
    @memcpy(qe[0 .. n - 1], e0);
    qe[n - 1] = 0;
    try bdsqr.qrAlgorithm(F, qd, qe, null, null);
    const floor = std.math.floatMin(F) / std.math.floatEps(F);
    for (d, qd, 0..) |x, y, i| {
        try testing.expect(x >= 0);
        if (i > 0) try testing.expect(x <= d[i - 1]);
        const err = @abs(x - y);
        if (!(err <= tol * @max(y, floor))) {
            check.say("lasq1 {s} n={d} i={d}: dqds {e}, qr {e}\n", .{ @typeName(F), n, i, x, y });
            return error.TestUnexpectedResult;
        }
    }
}

test "lasq1: random bidiagonals agree with the QR iteration" {
    const al = testing.allocator;
    inline for (.{ f32, f64 }) |F| for ([_]usize{ 3, 4, 5, 10, 25, 60 }) |n| {
        const d = try al.alloc(F, n);
        defer al.free(d);
        const e = try al.alloc(F, n - 1);
        defer al.free(e);
        check.fill(F, d, n * 3 + 1);
        check.fill(F, e, n * 3 + 2);
        try agree(F, d, e, 40 * @as(F, @floatFromInt(n)) * std.math.floatEps(F));
    };
}

test "lasq1: a graded bidiagonal keeps its small values relatively" {
    const n = 14;
    var d: [n]f64 = undefined;
    var e: [n - 1]f64 = undefined;
    for (0..n) |i| d[i] = std.math.pow(f64, 10, -@as(f64, @floatFromInt(i)));
    for (0..n - 1) |i| e[i] = 0.5 * std.math.pow(f64, 10, -@as(f64, @floatFromInt(i)) - 0.5);
    try agree(f64, &d, &e, 100 * n * std.math.floatEps(f64));
    // and reversed
    std.mem.reverse(f64, &d);
    std.mem.reverse(f64, &e);
    try agree(f64, &d, &e, 100 * n * std.math.floatEps(f64));
}

test "lasq1: diagonal, a zero on the diagonal, n of 1 and 2" {
    var work: [24]f64 = undefined;
    // diagonal: the sorted absolute values
    var d = [_]f64{ -3, 1, 2, -0.5 };
    var e = [_]f64{ 0, 0, 0 };
    try testing.expectEqual(lasq.Outcome.converged, lasq.lasq1(f64, &d, &e, &work));
    try testing.expectEqualSlices(f64, &.{ 3, 2, 1, 0.5 }, &d);
    // a zero diagonal entry: a zero singular value
    var dz = [_]f64{ 1, 0, 2, 1.5, 0.75 };
    var ez = [_]f64{ 0.5, 0.25, 1, 0.125 };
    try testing.expectEqual(lasq.Outcome.converged, lasq.lasq1(f64, &dz, &ez, &work));
    try testing.expect(dz[4] <= 1e-15 * dz[0]);
    try agree(f64, &.{ 1, 0, 2, 1.5, 0.75 }, &.{ 0.5, 0.25, 1, 0.125 }, 1e-13);
    var d1 = [_]f64{-2};
    try testing.expectEqual(lasq.Outcome.converged, lasq.lasq1(f64, &d1, &.{}, &work));
    try testing.expectEqual(2, d1[0]);
    var d2 = [_]f64{ 3, -4 };
    var e2 = [_]f64{0};
    try testing.expectEqual(lasq.Outcome.converged, lasq.lasq1(f64, &d2, &e2, &work));
    try testing.expectEqualSlices(f64, &.{ 4, 3 }, &d2);
    // a NaN: nothing usable
    var dn = [_]f64{ 1, std.math.nan(f64), 2 };
    var en = [_]f64{ 1, 1 };
    try testing.expectEqual(lasq.Outcome.failed, lasq.lasq1(f64, &dn, &en, &work));
}
