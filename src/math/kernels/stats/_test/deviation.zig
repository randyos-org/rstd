//! `deviationSums` against a direct sum in `f128`, in every `Summation`
//! mode, its NaN and infinities, and its bits at any worker count.

const std = @import("std");
const types = @import("../../blas/types.zig");
const strided = @import("../../blas/strided.zig");
const threading = @import("../../blas/threading.zig");
const fast = @import("../../blas/blas1/fast.zig");
const check = @import("../../blas/_test/check.zig");
const test_io = @import("../../blas/_test/io.zig");
const api = @import("../../blas/blas1/api.zig");
const stats = @import("../__root__.zig");

const testing = std.testing;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const ConstStrided = strided.ConstStrided;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const Summation = api.Summation;

/// `Σ(x − m)` (`re`, `im`) and `Σ|x − m|²` in `f128`, and the sums of the
/// terms' magnitudes.
const Want = struct { re: f128 = 0, im: f128 = 0, sq: f128 = 0, mag: f128 = 0 };

fn want(comptime T: type, m: T, x: ConstStrided(T)) Want {
    var r: Want = .{};
    for (0..x.n) |i| {
        const a = x.get(i);
        if (comptime isComplex(T)) {
            const dr = @as(f128, a.re) - @as(f128, m.re);
            const di = @as(f128, a.im) - @as(f128, m.im);
            r.re += dr;
            r.im += di;
            r.sq += dr * dr + di * di;
            r.mag += @abs(dr) + @abs(di);
        } else {
            const d = @as(f128, a) - @as(f128, m);
            r.re += d;
            r.sq += d * d;
            r.mag += @abs(d);
        }
    }
    return r;
}

/// The error bound per mode, in epsilons of the terms' magnitudes: the
/// final rounding to `Float(T)` with a chain of `n` additions in the
/// accumulating type, compensated about one rounding, `.fast` as
/// `fast.zig`'s tests bound it.
fn bound(comptime T: type, s: Summation, n: usize, mag: f128) f128 {
    const F = Float(T);
    const len = n * @as(usize, if (isComplex(T)) 2 else 1);
    const e: f128 = std.math.floatEps(F);
    const wide: f128 = if (F == f32) std.math.floatEps(f64) else e;
    const k: f128 = @floatFromInt(len + 2);
    return mag * switch (s) {
        .plain => 2 * e + k * wide,
        .compensated => 4 * e,
        .fast => if (comptime fast.serves(T)) @as(f128, @floatFromInt(len / 4 + 2 * fast.W + 10)) * e else 2 * e + k * wide,
    };
}

fn parts(comptime T: type, r: stats.DeviationSums(T)) [3]f128 {
    if (comptime isComplex(T)) return .{ r.sum.re, r.sum.im, r.sum_sq };
    return .{ r.sum, 0, r.sum_sq };
}

const all_types = .{ f32, f64, Complex32, Complex64 };
const modes = [_]Summation{ .plain, .compensated, .fast };

test "deviationSums: within each mode's bound of an f128 sum" {
    const al = testing.allocator;
    const n_max = 9000;
    const sizes = [_]usize{ 0, 1, 2, 7, 63, 64, 65, 1000, 4097, n_max };
    const strides = [_]isize{ 1, 2, -3 };
    inline for (all_types) |T| {
        const buf = try al.alloc(T, 3 * n_max);
        defer al.free(buf);
        check.fill(T, buf, 11);
        // an offset the deviations remove, as a mean would
        const m = check.val(T, 0.125, -0.25);
        for (sizes) |n| for (strides) |st| {
            const x = try ConstStrided(T).init(buf[0..if (n == 0) 0 else (n - 1) * @abs(st) + 1], n, st);
            const w = want(T, m, x);
            for (modes) |s| {
                const got = parts(T, try stats.deviationSums(check.io, T, s, m, x));
                const b = bound(T, s, n, w.mag);
                const bq = bound(T, s, n, w.sq);
                if (!(@abs(got[0] - w.re) <= b and @abs(got[1] - w.im) <= b and @abs(got[2] - w.sq) <= bq)) {
                    std.debug.print("deviationSums {s} {t} n={d} stride={d}: {any}, want {d} {d} {d}\n", .{ @typeName(T), s, n, st, got, @as(f64, @floatCast(w.re)), @as(f64, @floatCast(w.im)), @as(f64, @floatCast(w.sq)) });
                    return error.TestUnexpectedResult;
                }
            }
        };
    }
}

test "deviationSums: a NaN or an infinity makes both sums NaN or infinite" {
    inline for (all_types) |T| {
        for (modes) |s| {
            var x: [100]T = undefined;
            check.fill(T, &x, 3);
            x[57] = check.val(T, std.math.nan(f64), 0);
            const r = parts(T, try stats.deviationSums(check.io, T, s, check.val(T, 0, 0), .of(&x)));
            try testing.expect(std.math.isNan(r[0]) and std.math.isNan(r[2]));
            x[57] = check.val(T, std.math.inf(f64), 0);
            const q = parts(T, try stats.deviationSums(check.io, T, s, check.val(T, 0, 0), .of(&x)));
            // infinite, or NaN where a compensation meets it
            try testing.expect(!std.math.isFinite(q[0]) and !std.math.isFinite(q[2]));
        }
    }
}

test "deviationSums: one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    const n = (1 << 16) + 37;
    inline for (all_types) |T| {
        const x = try al.alloc(T, 3 * n);
        defer al.free(x);
        check.fillMixed(T, x, 5);
        const m = check.val(T, 0.5, 0.25);
        const views = [_]ConstStrided(T){ .of(x[0..n]), try .init(x, n, -3) };
        for (views) |v| for (modes) |s| {
            var outs: [test_io.worker_runs.len][3]f128 = undefined;
            for (&outs, 0..) |*o, c| {
                const io = runs.io(c);
                // the same chunks at every count: only the workers differ
                threading.test_serial_max = 4096;
                o.* = parts(T, try stats.deviationSums(io, T, s, m, v));
            }
            for (outs[1..]) |o| for (0..3) |k| {
                try testing.expectEqual(@as(u128, @bitCast(outs[0][k])), @as(u128, @bitCast(o[k])));
            };
        };
    }
}
