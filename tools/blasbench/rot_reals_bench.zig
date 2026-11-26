//! The `rotreals` section (T183; only when named, any host): complex `rot`
//! by a real rotation through its own slot (`csrot_k`, `zdrot_k`) against
//! the real slot (`srot_k`, `drot_k`) on the same `2n` reals, unit strides,
//! one call each, at n = 1000 to 10^6: whether `blas1/rot_reals.zig`'s
//! route (OpenBLAS's, from n = 4096 on x86_64) pays on this core.

const std = @import("std");
const Io = std.Io;
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const libs = @import("libs.zig");

const blas = rstd.math.kernels.blas;

const timing: bench.Timing = .{ .min_ns = 20 * std.time.ns_per_ms, .batches = 15 };
const sizes = [_]usize{ 1000, 2000, 4096, 10_000, 100_000, 1_000_000 };
const Slot = blas.Slot;
/// The float, its name, the complex slot and the real one.
const rows = .{ .{ f32, "c32", .csrot_k, .srot_k }, .{ f64, "c64", .zdrot_k, .drot_k } };

pub fn run(io: Io, gpa: std.mem.Allocator) !void {
    _ = io;
    try bb.clock(libs.rstd_io, "rotreals start");
    std.debug.print(
        \\
        \\core {s}; ns per complex element, best of {d} batches of at least {d} ms, one call per kernel, unit strides; speedup is complex over real
        \\
        \\| type | n | complex slot | real slot, 2n | speedup |
        \\| ---- | -: | -----------: | ------------: | ------: |
        \\
    , .{ blas.coreName(), timing.batches, @divTrunc(timing.min_ns, std.time.ns_per_ms) });
    inline for (rows) |r| {
        for (sizes) |n| {
            if (n > libs.n_max) continue;
            try race(r[0], r[1], r[2], r[3], gpa, n);
        }
    }
}

fn race(comptime F: type, name: []const u8, comptime complex: Slot, comptime real: Slot, gpa: std.mem.Allocator, n: usize) !void {
    const xy = try gpa.alloc(F, 4 * n);
    defer gpa.free(xy);
    var rng = std.Random.DefaultPrng.init(n);
    bb.fillRandom(F, xy, rng.random());
    var ctx = [2]Call(F, complex, real){ .{ .x = xy[0 .. 2 * n], .y = xy[2 * n ..], .reals = false }, .{ .x = xy[0 .. 2 * n], .y = xy[2 * n ..], .reals = true } };
    const calls = [2]bench.Call{ .of(&ctx[0]), .of(&ctx[1]) };
    var secs: [2]f64 = undefined;
    try bench.race(libs.rstd_io, timing, &calls, &secs);
    const per = 1e9 / @as(f64, @floatFromInt(n));
    std.debug.print("| {s} | {d} | {d:.3} | {d:.3} | {d:.2} |\n", .{ name, n, secs[0] * per, secs[1] * per, secs[0] / secs[1] });
}

/// A unit rotation, so repeated calls keep the data's size.
fn Call(comptime F: type, comptime complex: Slot, comptime real: Slot) type {
    return struct {
        x: []F,
        y: []F,
        reals: bool,

        pub fn run(o: *anyopaque) void {
            const c: *@This() = @ptrCast(@alignCast(o));
            if (c.reals) {
                _ = blas.kernel(real)(@intCast(c.x.len), c.x.ptr, 1, c.y.ptr, 1, 0.6, 0.8);
            } else {
                _ = blas.kernel(complex)(@intCast(c.x.len / 2), c.x.ptr, 1, c.y.ptr, 1, 0.6, 0.8);
            }
        }
    };
}
