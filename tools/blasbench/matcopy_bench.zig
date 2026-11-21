//! The `matcopy` section (T185; only when named, any host): the real
//! out-of-place copy of an n by n matrix, `f32` and `f64`, n = 64 to 4096,
//! as a pure copy (`alpha` 1) and scaled. Transposed: the plain loops
//! (`ext/matcopy.zig`'s before the tiles, written out here) against the
//! tiled kernel in one call (`ext/transpose.zig`) and the public `omatcopy`
//! on rstd's `Io` (`--threads`), which bands it. Not transposed: the kernel
//! in one call against the public `omatcopy`.

const std = @import("std");
const Io = std.Io;
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const libs = @import("libs.zig");

const blas = rstd.math.kernels.blas;
const fail = libs.fail;

const timing: bench.Timing = .{ .min_ns = 20 * std.time.ns_per_ms, .batches = 9 };
const sizes = [_]usize{ 64, 128, 256, 512, 1024, 2048, 4096 };
const Way = enum { t_plain, t_tiled, t_public, n_kernel, n_public };
const ways = std.enums.values(Way);
/// The float, the transposed kernel's slot and the plain one's.
const rows = .{ .{ f32, .somatcopy_k_ct, .somatcopy_k_cn }, .{ f64, .domatcopy_k_ct, .domatcopy_k_cn } };

pub fn run(io: Io, gpa: std.mem.Allocator) !void {
    _ = io;
    try bb.clock(libs.rstd_io, "matcopy start");
    std.debug.print(
        \\
        \\ns per element, best of {d} batches of at least {d} ms; `public` is `omatcopy` on rstd's `Io` (core {s}, {d} cpus, ceiling {?d}), the rest one call on one thread
        \\
        \\| type | alpha | n | t plain | t tiled | plain/tiled | t public | tiled/public | n kernel | n public | kernel/public |
        \\| ---- | ----: | -: | ------: | ------: | ----------: | -------: | -----------: | -------: | -------: | ------------: |
        \\
    , .{ timing.batches, @divTrunc(timing.min_ns, std.time.ns_per_ms), blas.coreName(), blas.config.cpu_count, blas.config.threadCeiling() });
    inline for (rows) |r| {
        const F = r[0];
        for ([_]F{ 1, 1.5 }) |alpha| for (sizes) |n| {
            if (n > libs.n_max) continue;
            try race(F, r[1], r[2], gpa, n, alpha);
        };
    }
}

fn race(comptime F: type, comptime ct: blas.Slot, comptime cn: blas.Slot, gpa: std.mem.Allocator, n: usize, alpha: F) !void {
    const ab = try gpa.alloc(F, 2 * n * n);
    defer gpa.free(ab);
    var rng = std.Random.DefaultPrng.init(n);
    bb.fillRandom(F, ab, rng.random());
    var ctx: [ways.len]Call(F, ct, cn) = undefined;
    for (&ctx, ways) |*c, w| c.* = .{ .a = ab[0 .. n * n], .b = ab[n * n ..], .n = n, .alpha = alpha, .way = w };
    var calls: [ways.len]bench.Call = undefined;
    for (&calls, &ctx) |*c, *x| c.* = .of(x);
    var secs: [ways.len]f64 = undefined;
    try bench.race(libs.rstd_io, timing, &calls, &secs);
    if (libs.lib_err) |err| return err;
    var ns: [ways.len]f64 = undefined;
    for (&ns, secs) |*o, s| o.* = s * 1e9 / @as(f64, @floatFromInt(n * n));
    std.debug.print("| {s} | {d} | {d} | {d:.3} | {d:.3} | {d:.2} | {d:.3} | {d:.2} | {d:.3} | {d:.3} | {d:.2} |\n", .{ @typeName(F), alpha, n, ns[0], ns[1], ns[0] / ns[1], ns[2], ns[1] / ns[2], ns[3], ns[4], ns[3] / ns[4] });
}

/// The transposed kernel's loops before the tiles: each element scaled as
/// `ext/matcopy.zig` scales it.
noinline fn plain(comptime F: type, n: usize, alpha: F, a: [*]const F, b: [*]F) void {
    for (0..n) |j| for (0..n) |i| {
        const x = a[i + j * n];
        b[j + i * n] = if (alpha == 0) 0 else if (alpha == 1) x else alpha * x;
    };
}

fn Call(comptime F: type, comptime ct: blas.Slot, comptime cn: blas.Slot) type {
    return struct {
        a: []F,
        b: []F,
        n: usize,
        alpha: F,
        way: Way,

        pub fn run(o: *anyopaque) void {
            const c: *@This() = @ptrCast(@alignCast(o));
            const io = libs.rstd_io;
            const n = c.n;
            switch (c.way) {
                .t_plain => plain(F, n, c.alpha, c.a.ptr, c.b.ptr),
                .t_tiled => _ = blas.kernel(ct)(@intCast(n), @intCast(n), c.alpha, c.a.ptr, @intCast(n), c.b.ptr, @intCast(n)),
                .t_public => blas.omatcopy(io, F, .t, n, n, c.alpha, c.a, n, c.b, n) catch |e| return fail(e),
                .n_kernel => _ = blas.kernel(cn)(@intCast(n), @intCast(n), c.alpha, c.a.ptr, @intCast(n), c.b.ptr, @intCast(n)),
                .n_public => blas.omatcopy(io, F, .n, n, n, c.alpha, c.a, n, c.b, n) catch |e| return fail(e),
            }
        }
    };
}
