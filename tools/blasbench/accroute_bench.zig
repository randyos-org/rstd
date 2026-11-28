//! The `accroute` section (T257, issue #797; only when named, macOS): rstd's
//! public routines with Accelerate's routing on and off (`useAccelerate`) in
//! one process, the two calls' batches in turn, for the slots the sweep of
//! issue #772 did not time. `routed/own` is the routed call's rate over rstd's
//! own: under 0.97, `apple_accel/speed.zig` keeps rstd's kernels for the slot.

const std = @import("std");
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const libs = @import("libs.zig");
const rc = @import("ext/race.zig");
const aw = @import("accroute/work.zig");
const a1 = @import("accroute/level1.zig");
const a2 = @import("accroute/level2.zig");
const a3 = @import("accroute/level3.zig");

const Io = std.Io;
const blas = rstd.math.kernels.blas;
const Prec = bb.Prec;
const precs = bb.precs;
const fail = libs.fail;
const Size = aw.Size;
const Work = aw.Work;
const sizes = aw.sizes;

const cases = a1.cases ++ a2.cases ++ a3.cases;

pub fn run(io: Io, gpa: std.mem.Allocator) !void {
    _ = io;
    if (!blas.useAccelerate(true)) return error.AccelerateNotRouted;
    try bb.clock(libs.rstd_io, "accroute start");
    std.debug.print(
        \\
        \\ns a call, best of the batches (level 1 and 2: {d} of at least {d} ms; level 3: {d} of {d} ms), routed and own in turn; rstd threads up to {?d} (null: the CPUs), Accelerate single-threaded in rstd's process; `routed/own` is the routed call's rate over rstd's own kernels'; `diff` the routed result's largest difference from rstd's own, in eps of the type relative to its largest entry; unit strides but `inc`, upper, no transpose, non-unit; `x` restored before each triangular call (counted)
        \\
        \\| routine | p | n | inc | routed | own | routed/own | diff |
        \\| ------- | - | -: | --: | -----: | --: | ---------: | ---: |
        \\
    , .{ bb.l1_timing.batches, @divTrunc(bb.l1_timing.min_ns, std.time.ns_per_ms), bb.l3_timing.batches, @divTrunc(bb.l3_timing.min_ns, std.time.ns_per_ms), blas.config.threadCeiling() });
    inline for (cases) |C| {
        if (libs.wanted(C.name)) inline for (precs) |p| {
            if (comptime std.mem.findScalar(u8, C.letters, p.name[0]) != null) {
                for (sizes(C.shape, p)) |sz| {
                    if (sz.n > libs.n_max) continue;
                    try race(C, p, gpa, sz);
                }
            }
        };
    }
    _ = blas.useAccelerate(true);
}

/// One call with routing set as `routed` says.
fn Call(comptime C: type, comptime p: Prec) type {
    return struct {
        w: *Work(p),
        routed: bool,
        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            _ = blas.useAccelerate(s.routed);
            s.w.restore(C.restore);
            C.call(p, s.w) catch |e| fail(e);
        }
    };
}

fn race(comptime C: type, comptime p: Prec, gpa: std.mem.Allocator, sz: Size) !void {
    var w = try Work(p).init(gpa, C.shape, sz);
    defer w.deinit(gpa);
    w.fill(C.fill);
    var on: Call(C, p) = .{ .w = &w, .routed = true };
    var off: Call(C, p) = .{ .w = &w, .routed = false };
    // each from fresh operands once, to compare
    w.reset();
    Call(C, p).run(&off);
    const mine = try gpa.dupe(p.F, w.out(C.out));
    defer gpa.free(mine);
    w.reset();
    Call(C, p).run(&on);
    if (libs.lib_err) |err| return err;
    const diff = rc.diffNorm(p.F, mine, w.out(C.out));
    w.reset();
    const calls = [_]bench.Call{ .of(&on), .of(&off) };
    const t = if (C.shape == .mm) bb.l3_timing else bb.l1_timing;
    var secs: [2]f64 = undefined;
    try bench.race(libs.rstd_io, t, &calls, &secs);
    if (libs.lib_err) |err| return err;
    std.debug.print("| {s} | {s} | {d} | {d} | {d:.0} | {d:.0} | {d:.2} | {d:.1} |\n", .{ C.name, p.name, sz.n, sz.inc, secs[0] * 1e9, secs[1] * 1e9, secs[1] / secs[0], diff });
}
