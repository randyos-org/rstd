//! Tests of `stebr.zig` and its parts (`laed4.zig`, `laed8_br.zig`,
//! `laed9_br.zig`, `auxiliary.zig`'s `lamrg`): eigenvalues against a wider
//! type's QR iteration and against the 1-2-1 matrix's closed form.

const std = @import("std");
const check = @import("../../blas/_test/check.zig");
const test_io = @import("../../blas/_test/io.zig");
const threading = @import("../../blas/threading.zig");
const blas = @import("../../blas/__root__.zig");
const auxiliary = @import("../auxiliary.zig");
const tevd = @import("../tridiag_evd.zig");
const stebr_dc = @import("../stebr.zig");
const laed4_solver = @import("../laed4.zig");
const evd = @import("../evd.zig");
const lcheck = @import("check.zig");
const ops = @import("../ops.zig");

const testing = std.testing;
const Allocator = std.mem.Allocator;
const stebr = stebr_dc.stebr;

const Kind = enum { random, clustered, graded, toeplitz, wilkinson, glued, split };

/// The race's matrices (issue #508) and the classic hard ones, from `seed`.
fn build(comptime F: type, kind: Kind, d: []F, e: []F, seed: u64) void {
    const n = d.len;
    var prng = std.Random.DefaultPrng.init(seed);
    const r = prng.random();
    const eps = std.math.floatEps(F);
    for (d, 0..) |*x, i| {
        const u = r.float(F) * 2 - 1;
        const v = r.float(F) * 2 - 1;
        const off = i + 1 < n;
        switch (kind) {
            .random, .split => {
                x.* = u;
                if (off) e[i] = if (kind == .split and (i % 7 == 6 or (i >= n / 3 and i < n / 3 + 5))) 0 else v;
            },
            // the race's: eigenvalues near 1
            .clustered => {
                x.* = 1 + 1e-10 * u;
                if (off) e[i] = 1e-6 * v;
            },
            .graded => {
                x.* = std.math.pow(F, 10, -@as(F, @floatFromInt(i % 12))) * u;
                if (off) e[i] = std.math.pow(F, 10, -@as(F, @floatFromInt((i % 12) + 1))) * v;
            },
            .toeplitz => {
                x.* = 2;
                if (off) e[i] = -1;
            },
            // W+: |m - i| down the diagonal, a close pair at the top
            .wilkinson => {
                x.* = @abs(@as(F, @floatFromInt(n / 2)) - @as(F, @floatFromInt(i)));
                if (off) e[i] = 1;
            },
            // W+ of order 21 glued by 1000 eps, past the split test
            .glued => {
                x.* = @abs(10 - @as(F, @floatFromInt(i % 21)));
                if (off) e[i] = if (i % 21 == 20) 1000 * eps else 1;
            },
        }
    }
}

/// The largest error over `eps ‖T‖max`.
fn errorIn(comptime F: type, comptime W: type, got: []const F, want: []const W, tmax: W) f64 {
    var m: W = 0;
    for (got, want) |g, w| m = @max(m, @abs(@as(W, g) - w));
    return @floatCast(m / (std.math.floatEps(F) * @max(tmax, std.math.floatMin(W))));
}

fn absMax(comptime F: type, d: []const F, e: []const F) F {
    return @max(ops.absMax(F, d), ops.absMax(F, e));
}

/// Eigenvalues by `stebr` in `F` and by QR iteration in `W` on the same
/// (`F`-rounded) tridiagonal; the error in units of `eps ‖T‖max`, and the
/// `F` QR path's for comparison.
fn compare(comptime F: type, comptime W: type, al: Allocator, kind: Kind, n: usize, seed: u64) !struct { stebr: f64, qr: f64 } {
    const d = try al.alloc(F, n);
    defer al.free(d);
    const e = try al.alloc(F, n);
    defer al.free(e);
    e[n - 1] = 0;
    build(F, kind, d, e, seed);
    const dw = try al.alloc(W, n);
    defer al.free(dw);
    const ew = try al.alloc(W, n);
    defer al.free(ew);
    for (d, e, dw, ew) |x, y, *xw, *yw| {
        xw.* = x;
        yw.* = y;
    }
    const tmax: W = absMax(W, dw, ew);
    if (kind == .toeplitz) {
        // 2 - 2 cos(k π / (n + 1))
        for (dw, 1..) |*x, k| x.* = 2 - 2 * @cos(@as(W, @floatFromInt(k)) * std.math.pi / @as(W, @floatFromInt(n + 1)));
    } else try tevd.qrAlgorithm(W, dw, ew[0 .. n - 1], null);
    const d2 = try al.dupe(F, d);
    defer al.free(d2);
    const e2 = try al.dupe(F, e);
    defer al.free(e2);
    try tevd.qrAlgorithm(F, d2, e2[0 .. n - 1], null);
    try stebr(al, F, d, e[0 .. n - 1]);
    for (1..n) |i| if (!(d[i - 1] <= d[i])) {
        check.say("stebr {s} {t} n={d}: out of order at {d}\n", .{ @typeName(F), kind, n, i });
        return error.TestUnexpectedResult;
    };
    return .{ .stebr = errorIn(F, W, d, dw, tmax), .qr = errorIn(F, W, d2, dw, tmax) };
}

/// The bound: `8 sqrt(n)` units of `eps ‖T‖max`, past QR iteration's own.
fn expectBound(comptime F: type, kind: Kind, n: usize, r: anytype) !void {
    const bound = 8 * @sqrt(@as(f64, @floatFromInt(n))) + 8;
    if (!(r.stebr <= bound)) {
        check.say("stebr {s} {t} n={d}: error {d:.1} eps |T|, QR {d:.1}, bound {d:.1}\n", .{ @typeName(F), kind, n, r.stebr, r.qr, bound });
        return error.TestUnexpectedResult;
    }
}

test "stebr: eigenvalues within 8 sqrt(n) eps |T| of a wider QR iteration, every kind, f32 and f64, n 1 to 300" {
    const al = testing.allocator;
    inline for (.{ .{ f32, f64 }, .{ f64, f128 } }) |p| {
        for ([_]usize{ 1, 2, 3, 25, 26, 27, 51, 100, 300 }) |n| for (std.enums.values(Kind)) |kind| {
            // f128 by software: its reference stops at 100
            if (p[1] == f128 and n > 100 and kind != .toeplitz) continue;
            const r = try compare(p[0], p[1], al, kind, n, n * 31 + @backingInt(kind));
            try expectBound(p[0], kind, n, r);
        };
    }
}

test "stebr: past a thousand, against f64 QR iteration (f32) and the 1-2-1 closed form" {
    const al = testing.allocator;
    for ([_]usize{ 1000, 2049 }) |n| {
        for (std.enums.values(Kind)) |kind| try expectBound(f32, kind, n, try compare(f32, f64, al, kind, n, n + @backingInt(kind)));
        try expectBound(f64, .toeplitz, n, try compare(f64, f128, al, .toeplitz, n, n));
    }
}

test "stebr: zero, diagonal and a split into ones; the input untouched below 2" {
    const al = testing.allocator;
    var d: [40]f64 = @splat(0);
    var e: [39]f64 = @splat(0);
    try stebr(al, f64, &d, &e);
    for (d) |x| try testing.expectEqual(@as(f64, 0), x);
    for (&d, 0..) |*x, i| x.* = @floatFromInt((i * 17) % 40);
    try stebr(al, f64, &d, &e);
    for (d, 0..) |x, i| try testing.expectEqual(@as(f64, @floatFromInt(i)), x);
    var one = [_]f64{-3};
    try stebr(testing.failing_allocator, f64, &one, e[0..0]);
    try testing.expectEqual(@as(f64, -3), one[0]);
}

test "stebr: a NaN ends (NoConvergence or NaN out), as QR iteration" {
    const al = testing.allocator;
    var d: [60]f64 = undefined;
    var e: [60]f64 = undefined;
    build(f64, .random, &d, &e, 3);
    d[30] = std.math.nan(f64);
    stebr(al, f64, &d, e[0..59]) catch |x| try testing.expectEqual(error.NoConvergence, x);
}

test "stebr: OutOfMemory" {
    var d: [60]f64 = undefined;
    var e: [60]f64 = undefined;
    build(f64, .random, &d, &e, 4);
    try testing.expectError(error.OutOfMemory, stebr(testing.failing_allocator, f64, &d, e[0..59]));
}

test "laed4: every root of random secular equations against f128's, n 1 to 40" {
    inline for (.{ f32, f64 }) |F| {
        for ([_]usize{ 1, 2, 3, 4, 7, 40 }) |n| for (0..4) |seed| {
            var prng = std.Random.DefaultPrng.init(n * 4 + seed);
            const r = prng.random();
            var d: [40]F = undefined;
            var z: [40]F = undefined;
            var d128: [40]f128 = undefined;
            var z128: [40]f128 = undefined;
            var acc: F = 0;
            for (0..n) |i| {
                acc += 0.01 + r.float(F);
                d[i] = acc;
                z[i] = r.float(F) * 2 - 1;
            }
            var zz: F = 0;
            for (z[0..n]) |x| zz += x * x;
            for (z[0..n]) |*x| x.* /= @sqrt(zz);
            for (0..n) |i| {
                d128[i] = d[i];
                z128[i] = z[i];
            }
            const rho: F = 0.5 + r.float(F);
            var delta: [40]F = undefined;
            var delta128: [40]f128 = undefined;
            for (0..n) |i| {
                const lam = try laed4_solver.laed4(F, i, d[0..n], z[0..n], delta[0..n], rho);
                const ref = try laed4_solver.laed4(f128, i, d128[0..n], z128[0..n], delta128[0..n], rho);
                const scale: f128 = @max(@abs(ref), d128[n - 1]);
                try testing.expect(@abs(@as(f128, lam) - ref) <= 16 * std.math.floatEps(F) * scale);
                // between its poles; rounded onto one at worst
                try testing.expect(lam >= d[i] and (i == n - 1 or lam <= d[i + 1]));
                // n <= 2: delta is the unit eigenvector
                if (n <= 2) {
                    var s: F = 0;
                    for (delta[0..n]) |x| s += x * x;
                    try testing.expectApproxEqAbs(@as(F, 1), s, 8 * std.math.floatEps(F));
                }
            }
        };
    }
}

test "lamrg: an ascending and a descending list merged, ties first-list first" {
    const a = [_]f64{ 1, 3, 5, 6, 4, 3, 2 };
    var idx: [7]usize = undefined;
    auxiliary.lamrg(f64, 3, 4, &a, true, false, &idx);
    try testing.expectEqualSlices(usize, &.{ 0, 6, 1, 5, 4, 2, 3 }, &idx);
    auxiliary.lamrg(f64, 0, 3, a[0..3], true, true, idx[0..3]);
    try testing.expectEqualSlices(usize, &.{ 0, 1, 2 }, idx[0..3]);
    auxiliary.lamrg(f64, 3, 0, a[4..7], false, true, idx[0..3]);
    try testing.expectEqualSlices(usize, &.{ 2, 1, 0 }, idx[0..3]);
}

test "stebr: syevd's values alone give one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    const n = 300;
    const a0 = try ops.alloc(f64, al, n, n);
    defer ops.free(f64, al, a0);
    lcheck.hermitian(f64, a0, 11);
    var outs: [3][]f64 = undefined;
    for (&outs) |*o| o.* = try al.alloc(f64, n);
    defer for (outs) |o| al.free(o);
    const a = try al.alloc(f64, n * n);
    defer al.free(a);
    for (0..test_io.worker_runs.len) |c| {
        @memcpy(a, a0.ptr[0 .. n * n]);
        const io = runs.io(c);
        // the same chunks for every run, one worker's too: the reduction's
        // `symv` folds partials by chunk, whose count follows the size alone
        threading.test_serial_max = 4096;
        try evd.syevd(al, io, f64, .n, .l, n, a, n, outs[c]);
    }
    try test_io.expectSameBytes(outs);
}
