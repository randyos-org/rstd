//! The `l2thread` section (T265, issue #845; only when named, any host): the
//! level-2 routines `accroute` races (`accroute/level2.zig`), rstd's own
//! kernels on one worker against `--threads`' workers in one process, the two
//! calls' batches in turn, over orders around their thread rules' `serial_max`
//! (`threading/rules.zig`): where a second worker starts to pay. Leaves
//! Accelerate's routing off.

const std = @import("std");
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const libs = @import("libs.zig");
const aw = @import("accroute/work.zig");
const a2 = @import("accroute/level2.zig");

const Io = std.Io;
const blas = rstd.math.kernels.blas;
const Prec = bb.Prec;
const precs = bb.precs;
const fail = libs.fail;
const Work = aw.Work;

const orders = [_]usize{ 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048 };

pub fn run(io: Io, gpa: std.mem.Allocator) !void {
    _ = io;
    // initializes first, which reads the ceiling and the CPUs
    _ = blas.useAccelerate(false);
    const many = blas.config.threadCeiling() orelse blas.config.cpu_count;
    if (many < 2) return error.NeedsThreads;
    try bb.clock(libs.rstd_io, "l2thread start");
    std.debug.print(
        \\
        \\ns a call, best of {d} batches of at least {d} ms, rstd's own kernels on 1 worker and on {d} in turn; `many/one` is the {d}-worker rate over one's: over 1 where more workers pay
        \\
        \\| routine | p | n | one | many | many/one |
        \\| ------- | - | -: | --: | ---: | -------: |
        \\
    , .{ bb.l1_timing.batches, @divTrunc(bb.l1_timing.min_ns, std.time.ns_per_ms), many, many });
    inline for (a2.cases) |C| {
        if (libs.wanted(C.name)) inline for (precs) |p| {
            if (comptime std.mem.findScalar(u8, C.letters, p.name[0]) != null) {
                for (orders) |n| {
                    if (n > libs.n_max) continue;
                    try race(C, p, gpa, n, @intCast(many));
                }
            }
        };
    }
}

/// One call on `workers` workers (the thread ceiling).
fn Call(comptime C: type, comptime p: Prec) type {
    return struct {
        w: *Work(p),
        workers: i32,
        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            blas.config.num_threads = s.workers;
            s.w.restore(C.restore);
            C.call(p, s.w) catch |e| fail(e);
        }
    };
}

fn race(comptime C: type, comptime p: Prec, gpa: std.mem.Allocator, n: usize, many: i32) !void {
    var w = try Work(p).init(gpa, C.shape, .{ .n = n });
    defer w.deinit(gpa);
    w.fill(C.fill);
    w.reset();
    const before = blas.config.num_threads;
    defer blas.config.num_threads = before;
    var one: Call(C, p) = .{ .w = &w, .workers = 1 };
    var lots: Call(C, p) = .{ .w = &w, .workers = many };
    var secs: [2]f64 = undefined;
    try bench.race(libs.rstd_io, bb.l1_timing, &.{ .of(&one), .of(&lots) }, &secs);
    if (libs.lib_err) |err| return err;
    std.debug.print("| {s} | {s} | {d} | {d:.0} | {d:.0} | {d:.2} |\n", .{ C.name, p.name, n, secs[0] * 1e9, secs[1] * 1e9, secs[0] / secs[1] });
}
