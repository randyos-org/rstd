//! The `potrf` section (T193; only when named, any host): rstd's own upper
//! `potrf`, Accelerate's routing off, every type at n = 32 to 1024: `potf2`
//! (unblocked), `potrf` on one worker (the ceiling at 1), and `potrf` on
//! rstd's `Io` (`--threads`), with the gain of the workers. What
//! `lapack/interface.zig`'s `team_grain` was measured with.

const std = @import("std");
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const libs = @import("libs.zig");

const blas = rstd.math.kernels.blas;
const lapack = rstd.math.kernels.lapack;
const fail = libs.fail;

const timing: bench.Timing = .{ .min_ns = 20 * std.time.ns_per_ms, .batches = 9 };
const sizes = [_]usize{ 32, 64, 128, 256, 512, 1024 };
const types = .{ f32, f64, std.math.Complex(f32), std.math.Complex(f64) };
const Way = enum { unblocked, one, all };
const ways = std.enums.values(Way);

pub fn run() !void {
    const io = libs.rstd_io;
    _ = blas.useAccelerate(false);
    try bb.clock(io, "potrf start");
    std.debug.print(
        \\
        \\rstd's own upper `potrf`, GFLOP/s by n³/3 (four times as many complex), best of {d} batches of at least {d} ms; `A` is restored before each call (counted); core {s}, {d} cpus, ceiling {?d}
        \\
        \\| type | n | potf2 | one worker | workers | workers/one |
        \\| ---- | -: | ----: | ---------: | ------: | ----------: |
        \\
    , .{ timing.batches, @divTrunc(timing.min_ns, std.time.ns_per_ms), blas.coreName(), blas.config.cpu_count, blas.config.threadCeiling() });
    inline for (types) |T| for (sizes) |n| {
        if (n > libs.n_max) continue;
        try race(T, n);
    };
}

fn Call(comptime T: type) type {
    return struct {
        a: []T,
        a0: []const T,
        n: usize,
        way: Way,

        pub fn run(p: *anyopaque) void {
            const c: *@This() = @ptrCast(@alignCast(p));
            @memcpy(c.a, c.a0);
            const saved = blas.config.num_threads;
            defer blas.config.num_threads = saved;
            switch (c.way) {
                .unblocked => return lapack.potf2(T, .u, c.n, c.a, c.n) catch |e| fail(e),
                .one => blas.config.num_threads = 1,
                .all => {},
            }
            lapack.potrf(libs.rstd_alloc, libs.rstd_io, T, .u, c.n, c.a, c.n) catch |e| fail(e);
        }
    };
}

fn race(comptime T: type, n: usize) !void {
    const complex = comptime @typeInfo(T) == .@"struct";
    const F = if (complex) @FieldType(T, "re") else T;
    const gpa = libs.rstd_alloc;
    const a0 = try gpa.alloc(T, n * n);
    defer gpa.free(a0);
    const a = try gpa.alloc(T, n * n);
    defer gpa.free(a);
    var rng = std.Random.DefaultPrng.init(n);
    bb.fillRandom(F, @as([*]F, @ptrCast(a0.ptr))[0 .. n * n * (if (complex) 2 else 1)], rng.random());
    // Hermitian, diagonally dominant
    for (0..n) |j| {
        for (0..j) |i| a0[j + i * n] = if (complex) a0[i + j * n].conjugate() else a0[i + j * n];
        a0[j + j * n] = if (complex) .init(@floatFromInt(2 * n), 0) else @floatFromInt(n);
    }
    var ctx: [ways.len]Call(T) = undefined;
    var calls: [ways.len]bench.Call = undefined;
    for (&ctx, &calls, ways) |*c, *k, w| {
        c.* = .{ .a = a, .a0 = a0, .n = n, .way = w };
        k.* = .of(c);
    }
    var secs: [ways.len]f64 = undefined;
    try bench.race(libs.rstd_io, timing, &calls, &secs);
    if (libs.lib_err) |err| return err;
    const nf: f64 = @floatFromInt(n);
    const flops = nf * nf * nf / 3 * @as(f64, if (complex) 4 else 1) / 1e9;
    std.debug.print("| {s} | {d} | {d:.1} | {d:.1} | {d:.1} | {d:.2} |\n", .{ if (complex) (if (F == f32) "c32" else "c64") else @typeName(T), n, flops / secs[0], flops / secs[1], flops / secs[2], secs[1] / secs[2] });
}
