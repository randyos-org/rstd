//! Tests: `gerqf`, `orgrq` and `ormrq`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const test_io = @import("../../blas/_test/io.zig");
const blas = @import("../../blas/__root__.zig");
const lcheck = @import("check.zig");
const ops = @import("../ops.zig");
const rq = @import("../rq.zig");

const testing = std.testing;
const Allocator = std.mem.Allocator;
const Mat = ops.Mat;
const Float = types.Float;
const Side = types.Side;

/// `drqt01.f`'s ratios: `‖R - A Qᴴ‖ / (n ‖A‖ eps)` and `‖I - Q Qᴴ‖ /
/// (n eps)` below the threshold, `Q` (`n x n`) from `orgrq`; then
/// `drqt03.f`'s: `ormrq` against `Q` for every side and `op`.
fn rqCase(comptime T: type, al: Allocator, m: usize, n: usize, seed: u64) !void {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const k = @min(m, n);
    const a0 = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, a0);
    check.fill(T, a0.slice(), seed);
    const a = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, a);
    a.copyFrom(a0);
    const tau = try al.alloc(T, @max(k, 1));
    defer al.free(tau);
    try rq.gerqf(al, check.io, T, m, n, a.slice(), a.ld, tau);
    var anorm: F = 0;
    for (a0.slice()) |v| anorm = @max(anorm, ops.abs(T, v));
    // the reflections, A's last k rows, into Q's last k
    const q = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, q);
    q.fill(ops.zero(T));
    q.sub(n - k, 0, k, n).copyFrom(a.sub(m - k, 0, k, n));
    try rq.orgrq(al, check.io, T, n, n, k, q.slice(), n, tau);
    const nf: F = @floatFromInt(@max(n, 1));
    // R = A Qᴴ, zero below the trapezoid ending at (m - 1, n - 1)
    const r = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, r);
    lcheck.mul(T, false, a0, true, q, r);
    const want = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, want);
    for (0..n) |j| for (0..m) |i| {
        want.at(i, j).* = if (j + m >= i + n) a.get(i, j) else ops.zero(T);
    };
    try lcheck.expectNear(T, "gerqf R", r, want, lcheck.thresh * nf * eps * @max(anorm, 1));
    const g = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, g);
    lcheck.mul(T, false, q, true, q, g);
    const id = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, id);
    id.identity();
    try lcheck.expectNear(T, "orgrq Q Qᴴ", g, id, lcheck.thresh * nf * eps);
    if (k == 0) return;
    const ops_list: []const types.Trans = if (comptime types.isComplex(T)) &.{ .n, .c } else &.{ .n, .t };
    const refl = a.sub(m - k, 0, k, n);
    for ([_]Side{ .l, .r }) |side| for (ops_list) |tr| {
        const cr = if (side == .l) n else 5;
        const cc = if (side == .l) 5 else n;
        const c = try ops.alloc(T, al, cr, cc);
        defer ops.free(T, al, c);
        check.fill(T, c.slice(), seed + 9);
        const want_c = try ops.alloc(T, al, cr, cc);
        defer ops.free(T, al, want_c);
        if (side == .l) lcheck.mul(T, tr != .n, q, false, c, want_c) else lcheck.mul(T, false, c, tr != .n, q, want_c);
        try rq.ormrq(al, check.io, T, side, tr, cr, cc, k, refl.slice(), refl.ld, tau, c.slice(), c.ld);
        try lcheck.expectNear(T, "ormrq", c, want_c, lcheck.thresh * nf * eps * 4);
    };
}

test "lapack gerqf: R = A Qᴴ and Q Qᴴ = I (orgrq), ormrq on both sides, square, wide, tall, blocked" {
    const al = testing.allocator;
    inline for (lcheck.test_types) |T| {
        for ([_][2]usize{ .{ 1, 1 }, .{ 1, 4 }, .{ 3, 5 }, .{ 5, 3 }, .{ 10, 10 }, .{ 40, 70 }, .{ 70, 40 } }) |mn| {
            try rqCase(T, al, mn[0], mn[1], mn[0] * 89 + mn[1]);
        }
    }
}

test "lapack gerqf: a zero matrix gives tau = 0; empty; bad arguments" {
    const al = testing.allocator;
    inline for (lcheck.test_types) |T| {
        var a = [_]T{ ops.zero(T), ops.zero(T), ops.zero(T), ops.zero(T) };
        var tau = [_]T{ ops.one(T), ops.one(T) };
        try rq.gerqf(al, check.io, T, 2, 2, &a, 2, &tau);
        for (tau) |t| try testing.expect(ops.isZero(T, t));
        for (a) |v| try testing.expect(ops.isZero(T, v));
        try rq.gerqf(al, check.io, T, 0, 2, &a, 1, &tau);
    }
    var a: [4]f64 = undefined;
    var tau: [2]f64 = undefined;
    try testing.expectError(error.LeadingDimensionTooSmall, rq.gerqf(al, check.io, f64, 2, 2, &a, 1, &tau));
    try testing.expectError(error.SliceTooShort, rq.gerqf(al, check.io, f64, 2, 3, &a, 2, &tau));
    try testing.expectError(error.SliceTooShort, rq.gerqf(al, check.io, f64, 2, 2, &a, 2, tau[0..1]));
    try testing.expectError(error.TooManyRows, rq.orgrq(al, check.io, f64, 2, 1, 1, &a, 2, &tau));
    try testing.expectError(error.TooManyReflections, rq.orgrq(al, check.io, f64, 2, 2, 3, &a, 2, &tau));
    try testing.expectError(error.UnsupportedTranspose, rq.ormrq(al, check.io, f64, .l, .r, 2, 2, 1, &a, 1, &tau, &a, 2));
    try testing.expectError(error.TooManyReflections, rq.ormrq(al, check.io, f64, .r, .n, 2, 1, 2, &a, 2, &tau, &a, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, rq.ormrq(al, check.io, f64, .l, .n, 2, 2, 1, &a, 1, &tau, &a, 1));
}

test "lapack gerqf: one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        const m = 90;
        const n = 300;
        const a0 = try al.alloc(T, m * n);
        defer al.free(a0);
        check.fill(T, a0, 7);
        var outs: [3][]T = undefined;
        for (&outs) |*o| o.* = try al.alloc(T, m * n + m);
        defer for (outs) |o| al.free(o);
        for (0..test_io.worker_runs.len) |c| {
            const io = runs.io(c);
            @memcpy(outs[c][0 .. m * n], a0);
            try rq.gerqf(al, io, T, m, n, outs[c][0 .. m * n], m, outs[c][m * n ..]);
        }
        try test_io.expectSameBytes(outs);
    }
}
