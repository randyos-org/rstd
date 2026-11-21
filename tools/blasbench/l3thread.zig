//! `l3thread` (only when named): rstd's public `gemm`, square `n`, serial and
//! at 2, 4, 8 and 16 workers, medians in GFLOPS; the bench behind
//! `src/math/kernels/blas/threading.md`'s level-3 timings. Moved from
//! `blas3/threaded.zig`'s test, which forced each count through the driver:
//! here a count is rstd's ceiling (`blas.config.num_threads`), so its size
//! rule and `cpu_count` may use fewer.

const std = @import("std");
const Io = std.Io;
const rstd = @import("rstd");
const bb = @import("blas_bench.zig");

const blas = rstd.linalg.blas;

const counts = [_]u32{ 2, 4, 8, 16 };

pub fn run(gpa: std.mem.Allocator) !void {
    // its own pool: the bench's `Io` holds `--threads`, 1 by default
    var threaded: Io.Threaded = .init(gpa, .{});
    defer threaded.deinit();
    const io = threaded.io();
    try bb.clock(io, "l3thread start");
    std.debug.print("core {s}, cpu count {d}\n", .{ blas.coreName(), blas.config.cpu_count });
    const saved = blas.config.num_threads;
    defer blas.config.num_threads = saved;
    // scratch from a buffer touched once, so no call pays for fresh pages
    const scratch = try gpa.alloc(u8, 256 << 20);
    defer gpa.free(scratch);
    @memset(scratch, 0);
    var fba: std.heap.FixedBufferAllocator = .init(scratch);
    var prng: std.Random.DefaultPrng = .init(0x13);
    inline for (.{ f64, f32 }) |T| for ([_]usize{ 256, 512, 1000, 2000 }) |n| {
        const a = try gpa.alloc(T, n * n);
        defer gpa.free(a);
        const b = try gpa.alloc(T, n * n);
        defer gpa.free(b);
        const c = try gpa.alloc(T, n * n);
        defer gpa.free(c);
        bb.fillRandom(T, a, prng.random());
        bb.fillRandom(T, b, prng.random());
        @memset(c, 0);
        const reps: usize = if (n >= 2000) 9 else if (n >= 1000) 15 else 25;
        // serial, then each count, interleaved; the first round warms up
        var times: [1 + counts.len][25]f64 = undefined;
        for (0..reps + 1) |round| for (0..times.len) |x| {
            fba.reset();
            blas.config.num_threads = if (x == 0) 1 else @intCast(counts[x - 1]);
            const t0 = Io.Clock.awake.now(io).toNanoseconds();
            try blas.gemm(fba.allocator(), io, T, .n, .n, n, n, n, 1, a, n, b, n, 0, c, n);
            if (round > 0) times[x][round - 1] = @floatFromInt(Io.Clock.awake.now(io).toNanoseconds() - t0);
        };
        const flops = 2.0 * @as(f64, @floatFromInt(n * n * n));
        std.debug.print("{s} n={d}: serial {d:.1}", .{ @typeName(T), n, flops / median(times[0][0..reps]) });
        for (counts, 0..) |w, i| std.debug.print(" | w={d} {d:.1}", .{ w, flops / median(times[1 + i][0..reps]) });
        std.debug.print(" GFLOPS\n", .{});
    };
}

fn median(x: []f64) f64 {
    std.mem.sort(f64, x, {}, std.sort.asc(f64));
    return x[x.len / 2];
}
