//! The `fft` section (N20, issue #155; only when named, any host): `fft.Plan`
//! and `RealPlan` made with vDSP off (rstd's own transforms) against made
//! with it on (vDSP's on macOS), one transform a call, the two calls' batches
//! in turn, at every length vDSP takes up to 2^20 (`f 2^k`, `f` in 1, 3, 5,
//! 15), `f32` and `f64`. `vdsp/own` is vDSP's rate over rstd's own; plans take
//! vDSP where two runs had it at 0.97 or more (`taken`; `vdsp.faster`, #856).

const std = @import("std");
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const libs = @import("libs.zig");
const rc = @import("ext/race.zig");

const Io = std.Io;
const fft = rstd.math.kernels.fft;
const Complex = std.math.Complex;
const fail = libs.fail;

const timing: bench.Timing = .{ .min_ns = 10 * std.time.ns_per_ms, .batches = 15 };
/// `f 2^k` up to 2^20, `k` from 3, by `f`.
const sizes = blk: {
    var a: [80]usize = undefined;
    var c = 0;
    for ([_]usize{ 1, 3, 5, 15 }) |f| for (3..21) |k| if (f << k <= 1 << 20) {
        a[c] = f << k;
        c += 1;
    };
    const out = a[0..c].*;
    break :blk out;
};

pub fn run(io: Io, gpa: std.mem.Allocator) !void {
    _ = io;
    try bb.clock(libs.rstd_io, "fft start");
    std.debug.print(
        \\
        \\ns a transform, best of {d} batches of at least {d} ms, vDSP's plan and rstd's own in turn, one thread; `vdsp/own` is the vDSP plan's rate over rstd's own; `diff` the largest difference of the two outputs in eps of the type relative to the largest entry; `kind` rstd's own plan's; `taken` whether a plan made by default takes vDSP
        \\
        \\| transform | p | n | kind | taken | vdsp | own | vdsp/own | diff |
        \\| --------- | - | -: | ---- | ----- | ---: | --: | -------: | ---: |
        \\
    , .{ timing.batches, @divTrunc(timing.min_ns, std.time.ns_per_ms) });
    inline for (.{ f32, f64 }) |T| {
        for (sizes) |n| try raceComplex(T, gpa, n);
        for (sizes) |n| try raceReal(T, gpa, n);
    }
    _ = fft.useAccelerate(true);
}

/// Whether a plan of `P` made by default at `n` takes vDSP.
fn taken(comptime P: type, gpa: std.mem.Allocator, n: usize) ![]const u8 {
    var p = try P.init(gpa, n);
    defer p.deinit(gpa);
    return if (p.routed()) "yes" else "no";
}

/// Both plans of one kind at `n`, vDSP's at any of its lengths; null where
/// vDSP did not take the length.
fn plans(comptime P: type, gpa: std.mem.Allocator, n: usize) !?[2]P {
    _ = fft.useAccelerate(false);
    var own = try P.init(gpa, n);
    errdefer own.deinit(gpa);
    if (!fft.useAccelerate(true)) {
        own.deinit(gpa);
        return null;
    }
    fft.useAccelerateEverywhere(true);
    defer fft.useAccelerateEverywhere(false);
    var routed = try P.init(gpa, n);
    if (!routed.routed()) {
        routed.deinit(gpa);
        own.deinit(gpa);
        return null;
    }
    return .{ routed, own };
}

fn CallC(comptime T: type) type {
    return struct {
        plan: *const fft.Plan(T),
        data: []Complex(T),
        scratch: []Complex(T),
        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.plan.transform(.failing, s.data, s.scratch, .forward, .{}, .ortho) catch |e| fail(e);
        }
    };
}

fn raceComplex(comptime T: type, gpa: std.mem.Allocator, n: usize) !void {
    var ps = (try plans(fft.Plan(T), gpa, n)) orelse return;
    defer for (&ps) |*p| p.deinit(gpa);
    const x = try gpa.alloc(Complex(T), n);
    defer gpa.free(x);
    var rng = std.Random.DefaultPrng.init(n);
    for (x) |*e| e.* = .init(rng.random().float(T) - 0.5, rng.random().float(T) - 0.5);
    var bufs: [2][]Complex(T) = undefined;
    var scr: [2][]Complex(T) = undefined;
    for (0..2) |i| {
        bufs[i] = try gpa.dupe(Complex(T), x);
        scr[i] = try gpa.alloc(Complex(T), ps[i].scratchLen());
    }
    defer for (0..2) |i| {
        gpa.free(bufs[i]);
        gpa.free(scr[i]);
    };
    var calls: [2]CallC(T) = undefined;
    for (0..2) |i| calls[i] = .{ .plan = &ps[i], .data = bufs[i], .scratch = scr[i] };
    for (&calls) |*c| CallC(T).run(c);
    const diff = rc.diffNorm(T, floats(T, bufs[1]), floats(T, bufs[0]));
    // ortho keeps the size, so repeated calls stay finite
    var secs: [2]f64 = undefined;
    try bench.race(libs.rstd_io, timing, &.{ .of(&calls[0]), .of(&calls[1]) }, &secs);
    if (libs.lib_err) |e| return e;
    std.debug.print("| complex | {s} | {d} | {s} | {s} | {d:.0} | {d:.0} | {d:.2} | {d:.1} |\n", .{ @typeName(T), n, @tagName(ps[1].kind), try taken(fft.Plan(T), gpa, n), secs[0] * 1e9, secs[1] * 1e9, secs[1] / secs[0], diff });
}

fn CallR(comptime T: type) type {
    return struct {
        plan: *const fft.RealPlan(T),
        in: []T,
        out: []Complex(T),
        scratch: []T,
        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.plan.forward(.failing, s.out, s.in, s.scratch, .{}, .backward) catch |e| fail(e);
        }
    };
}

fn raceReal(comptime T: type, gpa: std.mem.Allocator, n: usize) !void {
    var ps = (try plans(fft.RealPlan(T), gpa, n)) orelse return;
    defer for (&ps) |*p| p.deinit(gpa);
    const x = try gpa.alloc(T, n);
    defer gpa.free(x);
    var rng = std.Random.DefaultPrng.init(n);
    for (x) |*e| e.* = rng.random().float(T) - 0.5;
    var outs: [2][]Complex(T) = undefined;
    var scr: [2][]T = undefined;
    for (0..2) |i| {
        outs[i] = try gpa.alloc(Complex(T), n / 2 + 1);
        scr[i] = try gpa.alloc(T, ps[i].scratchLen());
    }
    defer for (0..2) |i| {
        gpa.free(outs[i]);
        gpa.free(scr[i]);
    };
    var calls: [2]CallR(T) = undefined;
    for (0..2) |i| calls[i] = .{ .plan = &ps[i], .in = x, .out = outs[i], .scratch = scr[i] };
    for (&calls) |*c| CallR(T).run(c);
    const diff = rc.diffNorm(T, floats(T, outs[1]), floats(T, outs[0]));
    var secs: [2]f64 = undefined;
    try bench.race(libs.rstd_io, timing, &.{ .of(&calls[0]), .of(&calls[1]) }, &secs);
    if (libs.lib_err) |e| return e;
    std.debug.print("| real | {s} | {d} | {s} | {s} | {d:.0} | {d:.0} | {d:.2} | {d:.1} |\n", .{ @typeName(T), n, @tagName(ps[1].kind), try taken(fft.RealPlan(T), gpa, n), secs[0] * 1e9, secs[1] * 1e9, secs[1] / secs[0], diff });
}

/// A complex slice as its floats.
fn floats(comptime T: type, v: []const Complex(T)) []const T {
    return @as([*]const T, @ptrCast(v.ptr))[0 .. 2 * v.len];
}
