//! The `compensated` section (T138; only when named, any host): `sum`,
//! `asum` and `dot` summed plain against compensated (`blas.Summation`),
//! `f32` and `f64`, at n = 100 to 10^7, on rstd's `Io` (`--threads`): the
//! price of the compensated sums, for choosing the defaults.

const std = @import("std");
const Io = std.Io;
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const libs = @import("libs.zig");

const blas = rstd.linalg.blas;
const Summation = blas.Summation;
const fail = libs.fail;

const timing: bench.Timing = .{ .min_ns = 20 * std.time.ns_per_ms, .batches = 7 };
const sizes = [_]usize{ 100, 1000, 10_000, 100_000, 1_000_000, 10_000_000 };
const Op = enum { sum, asum, dot };

pub fn run(io: Io, gpa: std.mem.Allocator) !void {
    _ = io;
    try bb.clock(libs.rstd_io, "compensated start");
    std.debug.print(
        \\
        \\ns per element, best of {d} batches of at least {d} ms; rstd's `Io` ({s}); ratios are compensated over plain
        \\
        \\| type | n | sum | sum comp. | ratio | asum | asum comp. | ratio | dot | dot comp. | ratio |
        \\| ---- | -: | --: | --------: | ----: | ---: | ---------: | ----: | --: | --------: | ----: |
        \\
    , .{ timing.batches, @divTrunc(timing.min_ns, std.time.ns_per_ms), libs.io_note });
    inline for (.{ f32, f64 }) |T| {
        for (sizes) |n| {
            if (n > libs.n_max) continue;
            try race(T, gpa, n);
        }
    }
}

fn race(comptime T: type, gpa: std.mem.Allocator, n: usize) !void {
    const xy = try gpa.alloc(T, 2 * n);
    defer gpa.free(xy);
    var rng = std.Random.DefaultPrng.init(n);
    for (xy) |*v| v.* = rng.random().float(T) * 2 - 1;
    var ctx: [6]Call(T) = undefined;
    for (&ctx, 0..) |*c, i| c.* = .{ .x = xy[0..n], .y = xy[n..], .op = @fromBackingInt(@as(u2, @intCast(i / 2))), .s = if (i % 2 == 0) .plain else .compensated };
    var calls: [6]bench.Call = undefined;
    for (&calls, &ctx) |*c, *x| c.* = .of(x);
    var secs: [6]f64 = undefined;
    try bench.race(libs.rstd_io, timing, &calls, &secs);
    if (libs.lib_err) |err| return err;
    var ns: [6]f64 = undefined;
    for (&ns, secs) |*o, s| o.* = s * 1e9 / @as(f64, @floatFromInt(n));
    std.debug.print("| {s} | {d} | {d:.3} | {d:.3} | {d:.2} | {d:.3} | {d:.3} | {d:.2} | {d:.3} | {d:.3} | {d:.2} |\n", .{ @typeName(T), n, ns[0], ns[1], ns[1] / ns[0], ns[2], ns[3], ns[3] / ns[2], ns[4], ns[5], ns[5] / ns[4] });
}

fn Call(comptime T: type) type {
    return struct {
        x: []const T,
        y: []const T,
        op: Op,
        s: Summation,

        pub fn run(o: *anyopaque) void {
            const c: *@This() = @ptrCast(@alignCast(o));
            const io = libs.rstd_io;
            const r = switch (c.op) {
                .sum => blas.sumWith(io, T, .of(c.x), c.s),
                .asum => blas.asumWith(io, T, .of(c.x), c.s),
                .dot => blas.dotWith(io, T, .of(c.x), .of(c.y), c.s),
            };
            // kept, so the call is not dropped
            std.mem.doNotOptimizeAway(r catch |e| return fail(e));
        }
    };
}
