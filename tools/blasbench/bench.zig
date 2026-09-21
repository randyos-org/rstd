//! Races BLAS kernel implementations against each other over a size sweep.
//!
//! Every implementation is called through the reference's own kernel ABI
//! (`daxpy_k` and friends), so OpenBLAS's assembled `.S` can be linked in and
//! compared directly with no wrapper between them.
//!
//! Correctness is checked once per implementation, before timing. Timing then
//! runs a batch of calls inside one clock reading, because a single call on a
//! small vector is shorter than the clock costs to read; the batch size grows
//! until it is long enough to measure, and the best of several batches is
//! reported, the run least disturbed by the rest of the machine.

const std = @import("std");
const Io = std.Io;

/// The reference's kernel signature: `n`, two unused, `alpha`, then the
/// vectors and their strides, then two more the axpy kernels ignore.
pub const Kernel = *const fn (
    n: isize,
    dummy0: isize,
    dummy1: isize,
    alpha: f64,
    x: [*]f64,
    inc_x: isize,
    y: [*]f64,
    inc_y: isize,
    dummy: [*]f64,
    dummy2: isize,
) callconv(.c) i32;

pub const Impl = struct {
    name: []const u8,
    run: Kernel,
};

/// What one call costs, for turning seconds into a rate.
pub const Work = struct {
    /// Floating-point operations for a length of `n`.
    flops: *const fn (n: usize) f64,
    /// Bytes moved for a length of `n`; axpy reads x and y and writes y.
    bytes: *const fn (n: usize) f64,
};

/// Small enough that repeated application cannot make `y` grow out of range.
const alpha: f64 = 1e-9;

fn nanos(io: Io) i96 {
    return Io.Clock.awake.now(io).toNanoseconds();
}

/// Seconds for one call, from the fastest batch.
fn timeOne(io: Io, k: Kernel, n: usize, x: []f64, y: []f64, min_ns: i96) f64 {
    _ = k(@intCast(n), 0, 0, alpha, x.ptr, 1, y.ptr, 1, x.ptr, 0); // warm up

    var reps: usize = 1;
    while (true) {
        const t0 = nanos(io);
        for (0..reps) |_| _ = k(@intCast(n), 0, 0, alpha, x.ptr, 1, y.ptr, 1, x.ptr, 0);
        const dt = nanos(io) - t0;
        if (dt >= min_ns) break;
        reps *= 2;
    }

    var best: f64 = std.math.inf(f64);
    for (0..5) |_| {
        const t0 = nanos(io);
        for (0..reps) |_| _ = k(@intCast(n), 0, 0, alpha, x.ptr, 1, y.ptr, 1, x.ptr, 0);
        const dt = nanos(io) - t0;
        const s = @as(f64, @floatFromInt(dt)) / 1e9 / @as(f64, @floatFromInt(reps));
        if (s < best) best = s;
    }
    return best;
}

/// The largest relative difference from `y0 + alpha * x` after one call.
fn check(k: Kernel, n: usize, x: []f64, y: []f64, y0: []const f64) f64 {
    @memcpy(y, y0);
    _ = k(@intCast(n), 0, 0, alpha, x.ptr, 1, y.ptr, 1, x.ptr, 0);
    var worst: f64 = 0;
    for (0..n) |i| {
        const want = y0[i] + alpha * x[i];
        const r = @abs(y[i] - want) / @max(@abs(want), 1e-300);
        if (r > worst) worst = r;
    }
    for (n..y.len) |i| if (y[i] != y0[i]) {
        worst = std.math.inf(f64); // wrote past `n`
    };
    return worst;
}

/// Runs every implementation at every size, checking results then timing.
pub fn sweep(
    io: Io,
    gpa: std.mem.Allocator,
    impls: []const Impl,
    sizes: []const usize,
    work: Work,
    min_ns: i96,
) !void {
    const max_n = sizes[sizes.len - 1];
    const x = try gpa.alignedAlloc(f64, .@"64", max_n);
    const y = try gpa.alignedAlloc(f64, .@"64", max_n);
    const y0 = try gpa.alloc(f64, max_n);

    var rng = std.Random.DefaultPrng.init(12345);
    const r = rng.random();
    for (x, y0) |*a, *b| {
        a.* = r.float(f64) * 2 - 1;
        b.* = r.float(f64) * 2 - 1;
    }

    std.debug.print("GB/s (GFLOP/s), best of 5 batches\n{s:>10}", .{"n"});
    for (impls) |im| std.debug.print("{s:>18}", .{im.name});
    std.debug.print("\n", .{});

    for (sizes) |n| {
        std.debug.print("{d:>10}", .{n});
        for (impls) |im| {
            const err = check(im.run, n, x[0..n], y[0..n], y0[0..n]);
            if (err > 1e-12) {
                std.debug.print("{s:>18}", .{"WRONG"});
                continue;
            }
            @memcpy(y, y0);
            const s = timeOne(io, im.run, n, x[0..n], y[0..n], min_ns);
            var buf: [32]u8 = undefined;
            const t = std.fmt.bufPrint(&buf, "{d:.0} ({d:.1})", .{
                work.bytes(n) / s / 1e9,
                work.flops(n) / s / 1e9,
            }) catch "?";
            std.debug.print("{s:>18}", .{t});
        }
        std.debug.print("\n", .{});
    }
}

pub const axpy_work: Work = .{
    .flops = struct {
        fn f(n: usize) f64 {
            return @as(f64, @floatFromInt(n)) * 2; // one multiply, one add
        }
    }.f,
    .bytes = struct {
        fn f(n: usize) f64 {
            return @as(f64, @floatFromInt(n)) * 8 * 3; // read x, read y, write y
        }
    }.f,
};
