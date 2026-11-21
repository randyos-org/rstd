//! `blas1/interface.zig`'s tests at the float types without vtable slots.

const std = @import("std");
const types = @import("../../types.zig");
const util = @import("../../util.zig");
const check = @import("../check.zig");
const reference = @import("../reference.zig");
const blas2_interface = @import("../../blas2/interface.zig");
const interface = @import("../../blas1/interface.zig");

const isComplex = types.isComplex;
const Float = types.Float;
const zero = util.zero;
const copy = interface.copy;
const swap = interface.swap;
const scal = interface.scal;
const dot = interface.dot;
const dotc = interface.dotc;
const nrm2 = interface.nrm2;
const asum = interface.asum;
const sum = interface.sum;
const index = interface.index;
const extreme = interface.extreme;
const rot = interface.rot;
const rotm = interface.rotm;
const axpy = interface.axpy;
const axpby = interface.axpby;
const rotg = interface.rotg;

/// Level 1 over `T`, which has no vtable slots, against an `f128` reference:
/// `x` at stride 2 and `y` at stride -1.
fn anyFloat(comptime T: type) !void {
    const W = reference.Wide(T);
    const n = 37;
    const ulps = 2 * (n + 2);
    var x: [2 * n]T = undefined;
    var y: [n]T = undefined;
    check.fill(T, &x, 11);
    check.fill(T, &y, 12);
    var wx: [n]W = undefined;
    var wy: [n]W = undefined;
    for (0..n) |i| {
        wx[i] = reference.wide(T, x[2 * i]);
        wy[i] = reference.wide(T, y[n - 1 - i]);
    }
    const alpha: T = if (comptime isComplex(T)) .init(1.25, -0.5) else 1.25;
    const beta: T = if (comptime isComplex(T)) .init(-0.75, 0.25) else -0.75;
    const wa = reference.wide(T, alpha);
    const wb = reference.wide(T, beta);

    // reductions
    var d = zero(W);
    var dc = zero(W);
    var dscale: f128 = 0;
    var ssq: f128 = 0;
    var s1: f128 = 0;
    var best: usize = 0;
    for (0..n) |i| {
        d = reference.addW(d, reference.mulW(wx[i], wy[i]));
        dc = reference.addW(dc, reference.mulW(reference.conjW(wx[i]), wy[i]));
        dscale += reference.abs1(wx[i]) * reference.abs1(wy[i]);
        ssq += if (comptime isComplex(T)) wx[i].re * wx[i].re + wx[i].im * wx[i].im else wx[i] * wx[i];
        s1 += reference.abs1(wx[i]);
        if (reference.abs1(wx[i]) > reference.abs1(wx[best])) best = i;
    }
    try reference.expectNear(T, d, try dot(check.io, T, n, &x, 2, &y, -1), ulps, dscale);
    if (comptime isComplex(T)) try reference.expectNear(T, dc, try dotc(check.io, T, n, &x, 2, &y, -1), ulps, dscale);
    const norm = reference.sqrtW(ssq);
    try reference.expectNear(Float(T), norm, try nrm2(check.io, T, n, &x, 2), ulps + reference.rootUlps(Float(T)), norm);
    try reference.expectNear(Float(T), s1, try asum(check.io, T, n, &x, 2), ulps, s1);
    try reference.expectNear(Float(T), reference.abs1(wx[best]), try extreme(check.io, T, .largest_magnitude, n, &x, 2), 4, s1);
    // the first of near-equal keys may differ; its key must be the largest
    const at = try index(check.io, T, .largest_magnitude, n, &x, 2);
    try reference.expectNear(Float(T), reference.abs1(wx[best]), @floatCast(reference.abs1(wx[@intCast(at - 1)])), 4, reference.abs1(wx[best]));

    // updates, element by element
    var y2 = y;
    try axpy(check.io, T, n, alpha, &x, 2, &y2, -1);
    for (0..n) |i| {
        const want = reference.addW(reference.mulW(wa, wx[i]), wy[i]);
        try reference.expectNear(T, want, y2[n - 1 - i], 8, reference.abs1(wa) * reference.abs1(wx[i]) + reference.abs1(wy[i]));
    }
    y2 = y;
    try axpby(check.io, T, n, alpha, &x, 2, beta, &y2, -1);
    for (0..n) |i| {
        const want = reference.addW(reference.mulW(wa, wx[i]), reference.mulW(wb, wy[i]));
        try reference.expectNear(T, want, y2[n - 1 - i], 8, reference.abs1(wa) * reference.abs1(wx[i]) + reference.abs1(wb) * reference.abs1(wy[i]));
    }
    var x2 = x;
    try scal(check.io, T, n, alpha, &x2, 2);
    for (0..n) |i| try reference.expectNear(T, reference.mulW(wa, wx[i]), x2[2 * i], 8, reference.abs1(wa) * reference.abs1(wx[i]));
    const c: Float(T) = 0.6;
    const sn: Float(T) = 0.8;
    x2 = x;
    y2 = y;
    try rot(check.io, T, n, &x2, 2, &y2, -1, c, sn);
    const wc = reference.wide(Float(T), c);
    const ws = reference.wide(Float(T), sn);
    for (0..n) |i| {
        const mag = reference.abs1(wx[i]) + reference.abs1(wy[i]);
        const rx = if (comptime isComplex(T)) W.init(wc * wx[i].re + ws * wy[i].re, wc * wx[i].im + ws * wy[i].im) else wc * wx[i] + ws * wy[i];
        const ry = if (comptime isComplex(T)) W.init(wc * wy[i].re - ws * wx[i].re, wc * wy[i].im - ws * wx[i].im) else wc * wy[i] - ws * wx[i];
        try reference.expectNear(T, rx, x2[2 * i], 8, mag);
        try reference.expectNear(T, ry, y2[n - 1 - i], 8, mag);
    }
    x2 = x;
    y2 = y;
    try swap(check.io, T, n, &x2, 2, &y2, -1);
    try copy(check.io, T, n, &x2, 2, &y2, -1);
    for (0..n) |i| {
        try reference.expectNear(T, wy[i], x2[2 * i], 0, 0);
        try reference.expectNear(T, wy[i], y2[n - 1 - i], 0, 0);
    }
    if (comptime isComplex(T)) return;

    var ps = [5]T{ -1, 2, 3, -1, 0.5 };
    x2 = x;
    y2 = y;
    try rotm(check.io, T, n, &x2, 2, &y2, -1, &ps);
    for (0..n) |i| {
        const mag = 2 * @abs(wx[i]) + @abs(wy[i]);
        try reference.expectNear(T, 2 * wx[i] - wy[i], x2[2 * i], 8, mag);
        try reference.expectNear(T, 3 * wx[i] + 0.5 * wy[i], y2[n - 1 - i], 8, mag);
    }
    try reference.expectNear(T, sumOf(T, &wx), try sum(check.io, T, n, &x, 2), ulps, s1);
    var ga: T = 3;
    var gb: T = 4;
    var gc: T = undefined;
    var gs: T = undefined;
    rotg(T, &ga, &gb, &gc, &gs);
    try reference.expectNear(T, 5, ga, 4, 5);
    try reference.expectNear(T, 0.6, gc, 4, 1);
    try reference.expectNear(T, 0.8, gs, 4, 1);
}

fn sumOf(comptime T: type, w: []const reference.Wide(T)) f128 {
    var acc: f128 = 0;
    for (w) |v| acc += v;
    return acc;
}

test "interface: level 1 at any float type runs the portable kernels built for it" {
    inline for (check.any_floats) |F| {
        inline for (.{ F, std.math.Complex(F) }) |T| anyFloat(T) catch |e| {
            check.say("{s}\n", .{@typeName(T)});
            return e;
        };
    }
}

test "interface: f16 reductions accumulate wide, and counts past its range stay finite" {
    const H = f16;
    const C = std.math.Complex(H);
    const al = std.testing.allocator;
    // past 2048 an f16 sum of ones stalls; past 65504 it overflows
    const n = 5000;
    const ones = try al.alloc(H, n);
    defer al.free(ones);
    @memset(ones, 1);
    try check.expectClose(H, 5000, try asum(check.io, H, n, ones.ptr, 1), 0);
    try check.expectClose(H, 5000, try sum(check.io, H, n, ones.ptr, 1), 0);
    try check.expectClose(H, 5000, try dot(check.io, H, n, ones.ptr, 1, ones.ptr, 1), 0);
    const big = try al.alloc(H, 4000);
    defer al.free(big);
    for (big, 0..) |*v, i| v.* = if (i < 2000) 40 else -40;
    try check.expectClose(H, 0, try sum(check.io, H, 4000, big.ptr, 1), 0);
    const cones = try al.alloc(C, 3000);
    defer al.free(cones);
    @memset(cones, .init(1, 0));
    try check.expectClose(C, .init(3000, 0), try dot(check.io, C, 3000, cones.ptr, 1, cones.ptr, 1), 0);
    try check.expectClose(C, .init(3000, 0), try dotc(check.io, C, 3000, cones.ptr, 1, cones.ptr, 1), 0);
    // gemv's transposed sums, real and complex
    var y = [_]H{0};
    try blas2_interface.gemv(al, check.io, H, .t, n, 1, 1, ones.ptr, n, ones.ptr, 1, 0, &y, 1);
    try check.expectClose(H, 5000, y[0], 0);
    var cy = [_]C{.init(0, 0)};
    try blas2_interface.gemv(al, check.io, C, .c, 3000, 1, .init(1, 0), cones.ptr, 3000, cones.ptr, 1, .init(0, 0), &cy, 1);
    try check.expectClose(C, .init(3000, 0), cy[0], 0);
    // a zero stride: sqrt(n) |x[0]| with n past f16's range
    var x = [_]H{0.01};
    const want: f64 = @sqrt(70000.0) * @as(f64, @as(H, 0.01));
    try check.expectClose(f64, want, try nrm2(check.io, H, 70000, &x, 0), 4e-3);
    var cx = [_]C{.init(0.01, 0)};
    try check.expectClose(f64, want, try nrm2(check.io, C, 70000, &cx, 0), 4e-3);
    var yz = [_]H{0};
    try axpy(check.io, H, 70000, 1, &x, 0, &yz, 0);
    try check.expectClose(f64, 700, yz[0], 1);
}
