//! Tests: `getrfComplete` and `getrsComplete`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const test_io = @import("../../blas/_test/io.zig");
const threading = @import("../../blas/threading.zig");
const blas = @import("../../blas/__root__.zig");
const ops = @import("../ops.zig");
const complete = @import("../getrf_complete.zig");

const testing = std.testing;
const Float = types.Float;
const Trans = types.Trans;

/// `‖P A Q - L U‖` over the entries (faer's `reconstruct.rs`, then the
/// comparison), `f` holding the factors of the `m x n` `a`.
fn residual(comptime T: type, m: usize, n: usize, a: []const T, f: []const T, ipiv: []const usize, jpiv: []const usize) !Float(T) {
    const k = @min(m, n);
    // `P A Q`: the interchanges in order
    const pa = try testing.allocator.dupe(T, a);
    defer testing.allocator.free(pa);
    for (0..k) |i| {
        if (ipiv[i] != i) for (0..n) |j| std.mem.swap(T, &pa[i + j * m], &pa[ipiv[i] + j * m]);
        if (jpiv[i] != i) for (0..m) |r| std.mem.swap(T, &pa[r + i * m], &pa[r + jpiv[i] * m]);
    }
    var worst: Float(T) = 0;
    for (0..n) |j| for (0..m) |i| {
        var s = ops.zero(T);
        for (0..@min(i + 1, j + 1, k)) |l| {
            const lv = if (l == i) ops.one(T) else f[i + l * m];
            s = ops.add(T, s, ops.mul(T, lv, f[l + j * m]));
        }
        worst = @max(worst, ops.abs(T, ops.sub(T, s, pa[i + j * m])));
    };
    return worst;
}

test "getrfComplete: P A Q = L U, L bounded, every type and shape" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| for ([_][2]usize{ .{ 1, 1 }, .{ 1, 5 }, .{ 6, 1 }, .{ 7, 7 }, .{ 20, 9 }, .{ 9, 20 }, .{ 33, 33 } }) |s| {
        const m = s[0];
        const n = s[1];
        const k = @min(m, n);
        const a = try al.alloc(T, m * n);
        defer al.free(a);
        check.fill(T, a, m * 13 + n);
        const f = try al.dupe(T, a);
        defer al.free(f);
        var ipiv: [33]usize = undefined;
        var jpiv: [33]usize = undefined;
        try testing.expectEqual(k, try complete.getrfComplete(check.io, T, m, n, f, m, &ipiv, &jpiv));
        var big: Float(T) = 0;
        for (a) |v| big = @max(big, ops.abs(T, v));
        const tol = 50 * @as(Float(T), @floatFromInt(m + n)) * std.math.floatEps(Float(T)) * big;
        const res = try residual(T, m, n, a, f, ipiv[0..k], jpiv[0..k]);
        if (!(res <= tol)) {
            check.say("getrfComplete {s} {d}x{d}: residual {e}, tol {e}\n", .{ @typeName(T), m, n, res, tol });
            return error.TestUnexpectedResult;
        }
        // each pivot is the largest `|re| + |im|` of its column: `|L| <= 1`,
        // complex `<= sqrt(2)`
        const bound: Float(T) = if (comptime types.isComplex(T)) std.math.sqrt2 else 1;
        for (0..k) |j| for (j + 1..m) |i| {
            try testing.expect(ops.abs(T, f[i + j * m]) <= bound * (1 + 4 * std.math.floatEps(Float(T))));
        };
    };
}

test "getrfComplete: an exact rank below min(m, n), and a zero matrix" {
    inline for (.{ f64, types.Complex64 }) |T| {
        // rank 3: three random columns, the rest zero
        const m = 8;
        const n = 6;
        var a: [m * n]T = @splat(ops.zero(T));
        check.fill(T, a[0 .. 3 * m], 5);
        var f = a;
        var ipiv: [n]usize = undefined;
        var jpiv: [n]usize = undefined;
        try testing.expectEqual(@as(usize, 3), try complete.getrfComplete(check.io, T, m, n, &f, m, &ipiv, &jpiv));
        for (3..n) |i| {
            try testing.expectEqual(i, ipiv[i]);
            try testing.expectEqual(i, jpiv[i]);
        }
        try testing.expect(try residual(T, m, n, &a, &f, &ipiv, &jpiv) <= 1e-13);
        var z: [12]T = @splat(ops.zero(T));
        try testing.expectEqual(@as(usize, 0), try complete.getrfComplete(check.io, T, 3, 4, &z, 3, ipiv[0..3], jpiv[0..3]));
        try testing.expectEqualSlices(usize, &.{ 0, 1, 2 }, ipiv[0..3]);
        try testing.expectEqualSlices(usize, &.{ 0, 1, 2 }, jpiv[0..3]);
    }
    var a = [_]f64{ 1, 2 };
    var piv: [1]usize = undefined;
    try testing.expectError(error.SliceTooShort, complete.getrfComplete(check.io, f64, 2, 2, &a, 2, &piv, &piv));
    try testing.expectError(error.LeadingDimensionTooSmall, complete.getrfComplete(check.io, f64, 2, 1, &a, 1, &piv, &piv));
    var none: [0]usize = .{};
    try testing.expectError(error.SliceTooShort, complete.getrfComplete(check.io, f64, 1, 2, &a, 1, &none, &piv));
}

test "getrfComplete: one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    const m = 90;
    const n = 70;
    const a0 = try al.alloc(types.Complex64, m * n);
    defer al.free(a0);
    check.fill(types.Complex64, a0, 31);
    var outs: [test_io.worker_runs.len][]types.Complex64 = undefined;
    for (&outs) |*o| o.* = try al.alloc(types.Complex64, m * n);
    defer for (outs) |o| al.free(o);
    var pivs: [test_io.worker_runs.len][2 * n]usize = undefined;
    for (0..test_io.worker_runs.len) |c| {
        @memcpy(outs[c], a0);
        // many chunks for every run; the fold does not depend on them
        threading.test_serial_max = 64;
        _ = try complete.getrfComplete(runs.io(c), types.Complex64, m, n, outs[c], m, pivs[c][0..n], pivs[c][n..]);
    }
    try test_io.expectSameBytes(outs);
    var ps: [test_io.worker_runs.len][]const usize = undefined;
    for (&ps, &pivs) |*q, *p| q.* = p;
    try test_io.expectSameBytes(ps);
}

test "getrsComplete: op(A) X = B from getrfComplete, every trans" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| for ([_]usize{ 1, 5, 24 }) |n| {
        const nrhs = 3;
        const a = try al.alloc(T, n * n);
        defer al.free(a);
        check.fill(T, a, n * 7);
        for (0..n) |i| a[i + i * n] = ops.add(T, a[i + i * n], ops.scalar(T, 2));
        const f = try al.dupe(T, a);
        defer al.free(f);
        var ipiv: [24]usize = undefined;
        var jpiv: [24]usize = undefined;
        try testing.expectEqual(n, try complete.getrfComplete(check.io, T, n, n, f, n, &ipiv, &jpiv));
        const b = try al.alloc(T, n * nrhs);
        defer al.free(b);
        check.fill(T, b, n + 1);
        for ([_]Trans{ .n, .t, .r, .c }) |trans| {
            const x = try al.dupe(T, b);
            defer al.free(x);
            try complete.getrsComplete(al, check.io, T, trans, n, nrhs, f, n, ipiv[0..n], jpiv[0..n], x, n);
            var worst: Float(T) = 0;
            var xbig: Float(T) = 0;
            for (x) |v| xbig = @max(xbig, ops.abs(T, v));
            for (0..nrhs) |c| for (0..n) |i| {
                var s = ops.zero(T);
                for (0..n) |l| {
                    var e = if (trans.transposed()) a[l + i * n] else a[i + l * n];
                    if (trans.conjA()) e = ops.conj(T, e);
                    s = ops.add(T, s, ops.mul(T, e, x[l + c * n]));
                }
                worst = @max(worst, ops.abs(T, ops.sub(T, s, b[i + c * n])));
            };
            const tol = 200 * @as(Float(T), @floatFromInt(n * n)) * std.math.floatEps(Float(T)) * @max(xbig, 1);
            if (!(worst <= tol)) {
                check.say("getrsComplete {s} {t} n={d}: residual {e}, tol {e}\n", .{ @typeName(T), trans, n, worst, tol });
                return error.TestUnexpectedResult;
            }
        }
    };
    var a = [_]f64{ 1, 0, 0, 1 };
    var b = [_]f64{ 3, 4 };
    const bad = [_]usize{ 0, 2 };
    const ok = [_]usize{ 0, 1 };
    try testing.expectError(error.PivotOutOfRange, complete.getrsComplete(al, check.io, f64, .n, 2, 1, &a, 2, &ok, &bad, &b, 2));
    try testing.expectEqualSlices(f64, &.{ 3, 4 }, &b);
    try testing.expectError(error.SliceTooShort, complete.getrsComplete(al, check.io, f64, .n, 2, 2, &a, 2, &ok, &ok, &b, 2));
}
