//! The plans of `../__root__.zig` against a direct DFT, their round trips,
//! batches and errors; a stub's tests skip.

const std = @import("std");
const fft = @import("../__root__.zig");

const Complex = std.math.Complex;
const Plan = fft.Plan;
const RealPlan = fft.RealPlan;
const StaticPlan = fft.StaticPlan;
const TrigTransformPlan = fft.TrigTransformPlan;
const TrigTransformKind = fft.TrigTransformKind;
const ModifiedCosineTransform = fft.ModifiedCosineTransform;
const BatchLayouts = fft.BatchLayouts;

const testing = std.testing;

/// A stub's tests skip until its body lands.
fn skip(err: anyerror) anyerror {
    return if (err == error.NotImplemented) error.SkipZigTest else err;
}

/// A test signal with no symmetry to hide behind.
fn signal(comptime T: type, out: []Complex(T)) void {
    for (out, 0..) |*x, i| {
        const t: f64 = @floatFromInt(i);
        x.* = .init(@floatCast(@sin(1 + 0.7 * t)), @floatCast(@cos(0.3 * t * t)));
    }
}

/// The forward DFT by its definition, in f64.
fn directDft(comptime T: type, out: []Complex(f64), in: []const Complex(T)) void {
    const n = in.len;
    for (out, 0..) |*o, k| {
        var acc: Complex(f64) = .init(0, 0);
        for (in, 0..) |x, j| {
            const angle = -2 * std.math.pi * @as(f64, @floatFromInt((j * k) % n)) / @as(f64, @floatFromInt(n));
            const w: Complex(f64) = .init(@cos(angle), @sin(angle));
            acc = acc.add(w.mul(.init(@floatCast(x.re), @floatCast(x.im))));
        }
        o.* = acc;
    }
}

fn expectNear(comptime T: type, want: []const Complex(f64), got: []const Complex(T), tolerance: f64) !void {
    for (want, got) |w, g| {
        try testing.expectApproxEqAbs(w.re, g.re, tolerance);
        try testing.expectApproxEqAbs(w.im, g.im, tolerance);
    }
}

test "init rejects n = 0" {
    try testing.expectError(error.ZeroLength, Plan(f32).init(testing.allocator, 0));
    try testing.expectError(error.ZeroLength, RealPlan(f64).init(testing.allocator, 0));
    try testing.expectError(error.ZeroLength, TrigTransformPlan(f64).init(testing.allocator, 0, .dct2));
    try testing.expectError(error.ZeroLength, ModifiedCosineTransform(f32).init(testing.allocator, 0));
}

test "Plan matches a direct DFT, n = 1 to 64" {
    var in: [64]Complex(f32) = undefined;
    var out: [64]Complex(f32) = undefined;
    var want: [64]Complex(f64) = undefined;
    for (1..65) |n| {
        var plan = Plan(f32).init(testing.allocator, n) catch |e| return skip(e);
        defer plan.deinit(testing.allocator);
        const scratch = try testing.allocator.alloc(Complex(f32), plan.scratchLen());
        defer testing.allocator.free(scratch);
        signal(f32, in[0..n]);
        directDft(f32, want[0..n], in[0..n]);
        try plan.transformInto(testing.io, out[0..n], in[0..n], scratch, .forward, .{}, .backward);
        try expectNear(f32, want[0..n], out[0..n], 1e-5 * @as(f64, @floatFromInt(n)));
    }
}

test "Plan round trips, and ortho keeps the energy" {
    const sizes = [_]usize{ 1, 2, 7, 12, 97, 128, 360, 1009 };
    for (sizes) |n| {
        var plan = Plan(f64).init(testing.allocator, n) catch |e| return skip(e);
        defer plan.deinit(testing.allocator);
        const scratch = try testing.allocator.alloc(Complex(f64), plan.scratchLen());
        defer testing.allocator.free(scratch);
        const x = try testing.allocator.alloc(Complex(f64), n);
        defer testing.allocator.free(x);
        const y = try testing.allocator.alloc(Complex(f64), n);
        defer testing.allocator.free(y);
        signal(f64, x);
        @memcpy(y, x);

        try plan.transform(testing.io, y, scratch, .forward, .{}, .backward);
        try plan.transform(testing.io, y, scratch, .backward, .{}, .backward);
        for (x, y) |a, b| {
            try testing.expectApproxEqAbs(a.re, b.re, 1e-12);
            try testing.expectApproxEqAbs(a.im, b.im, 1e-12);
        }

        // Parseval
        try plan.transform(testing.io, y, scratch, .forward, .{}, .ortho);
        var ex: f64 = 0;
        var ey: f64 = 0;
        for (x, y) |a, b| {
            ex += a.re * a.re + a.im * a.im;
            ey += b.re * b.re + b.im * b.im;
        }
        try testing.expectApproxEqRel(ex, ey, 1e-12);
    }
}

test "Plan batch: interleaved channels" {
    const n = 16;
    var plan = Plan(f64).init(testing.allocator, n) catch |e| return skip(e);
    defer plan.deinit(testing.allocator);
    const scratch = try testing.allocator.alloc(Complex(f64), 2 * plan.scratchLen());
    defer testing.allocator.free(scratch);

    var both: [2 * n]Complex(f64) = undefined;
    signal(f64, &both);
    var want: [2][n]Complex(f64) = undefined;
    for (&want, 0..) |*w, c| {
        var channel: [n]Complex(f64) = undefined;
        for (&channel, 0..) |*x, i| x.* = both[2 * i + c];
        directDft(f64, w, &channel);
    }
    try plan.transform(testing.io, &both, scratch, .forward, .{ .count = 2, .stride = 2, .distance = 1 }, .backward);
    for (want, 0..) |w, c| {
        var got: [n]Complex(f64) = undefined;
        for (&got, 0..) |*x, i| x.* = both[2 * i + c];
        try expectNear(f64, &w, &got, 1e-12);
    }
}

test "RealPlan matches the complex DFT's half" {
    var re: [33]f64 = undefined;
    var full: [33]Complex(f64) = undefined;
    var want: [33]Complex(f64) = undefined;
    var half: [17]Complex(f64) = undefined;
    var back: [33]f64 = undefined;
    for (1..34) |n| {
        var plan = RealPlan(f64).init(testing.allocator, n) catch |e| return skip(e);
        defer plan.deinit(testing.allocator);
        const scratch = try testing.allocator.alloc(f64, plan.scratchLen());
        defer testing.allocator.free(scratch);
        signal(f64, full[0..n]);
        for (re[0..n], full[0..n]) |*r, *c| {
            r.* = c.re;
            c.im = 0;
        }
        directDft(f64, want[0..n], full[0..n]);
        const m = n / 2 + 1;
        try plan.forward(testing.io, half[0..m], re[0..n], scratch, .{}, .backward);
        try expectNear(f64, want[0..m], half[0..m], 1e-12);
        try plan.backward(testing.io, back[0..n], half[0..m], scratch, .{}, .backward);
        for (re[0..n], back[0..n]) |a, b| try testing.expectApproxEqAbs(a, b, 1e-12);
    }
}

test "RealPlan unordered: circular convolution" {
    const n = 64;
    var plan = RealPlan(f32).init(testing.allocator, n) catch |e| return skip(e);
    defer plan.deinit(testing.allocator);
    const scratch = try testing.allocator.alloc(f32, plan.scratchLen());
    defer testing.allocator.free(scratch);

    var a: [n]f32 = undefined;
    var b: [n]f32 = undefined;
    for (&a, &b, 0..) |*x, *y, i| {
        const t: f32 = @floatFromInt(i);
        x.* = @sin(0.4 * t);
        y.* = if (i < 5) 1 else 0; // a short box filter
    }
    var want: [n]f32 = @splat(0);
    for (&want, 0..) |*w, k| {
        for (0..n) |j| w.* += a[j] * b[(k + n - j) % n];
    }

    var fa: [n]f32 = undefined;
    var fb: [n]f32 = undefined;
    var acc: [n]f32 = @splat(0);
    var got: [n]f32 = undefined;
    plan.forwardUnordered(&fa, &a, scratch) catch |e| return skip(e);
    try plan.forwardUnordered(&fb, &b, scratch);
    try plan.convolveAccumulate(&acc, &fa, &fb, 1.0 / @as(f32, n));
    try plan.backwardUnordered(&got, &acc, scratch);
    for (want, got) |w, g| try testing.expectApproxEqAbs(w, g, 1e-4);
}

test "StaticPlan matches Plan" {
    const n = 12;
    const static: StaticPlan(n, f64) = .{};
    var plan = Plan(f64).init(testing.allocator, n) catch |e| return skip(e);
    defer plan.deinit(testing.allocator);
    const scratch = try testing.allocator.alloc(Complex(f64), @max(plan.scratchLen(), static.scratchLen()));
    defer testing.allocator.free(scratch);

    var in: [n]Complex(f64) = undefined;
    signal(f64, &in);
    var want: [n]Complex(f64) = undefined;
    var got: [n]Complex(f64) = undefined;
    try plan.transformInto(testing.io, &want, &in, scratch, .forward, .{}, .backward);
    static.transformInto(testing.io, &got, &in, scratch, .forward, .{}, .backward) catch |e| return skip(e);
    try testing.expectEqualSlices(Complex(f64), &want, &got); // dispatch by size alone: same bits
}

test "TrigTransformPlan: ortho DCT-II then DCT-III, and DCT-IV twice, give back the input" {
    const n = 10;
    const pairs = [_][2]TrigTransformKind{ .{ .dct2, .dct3 }, .{ .dst2, .dst3 }, .{ .dct4, .dct4 }, .{ .dst4, .dst4 } };
    for (pairs) |pair| {
        var there = TrigTransformPlan(f64).init(testing.allocator, n, pair[0]) catch |e| return skip(e);
        defer there.deinit(testing.allocator);
        var back = try TrigTransformPlan(f64).init(testing.allocator, n, pair[1]);
        defer back.deinit(testing.allocator);
        const scratch = try testing.allocator.alloc(f64, @max(there.scratchLen(), back.scratchLen()));
        defer testing.allocator.free(scratch);

        var x: [n]f64 = undefined;
        for (&x, 0..) |*v, i| v.* = @sin(1 + 0.7 * @as(f64, @floatFromInt(i)));
        var y = x;
        try there.transform(testing.io, &y, scratch, .{}, .ortho);
        try back.transform(testing.io, &y, scratch, .{}, .ortho);
        for (x, y) |a, b| try testing.expectApproxEqAbs(a, b, 1e-12);
    }
}

test "ModifiedCosineTransform: TDAC with a sine window" {
    const n = 16;
    var mdct = ModifiedCosineTransform(f64).init(testing.allocator, n) catch |e| return skip(e);
    defer mdct.deinit(testing.allocator);
    const scratch = try testing.allocator.alloc(f64, mdct.scratchLen());
    defer testing.allocator.free(scratch);

    var window: [2 * n]f64 = undefined;
    for (&window, 0..) |*w, i| w.* = @sin(std.math.pi * (@as(f64, @floatFromInt(i)) + 0.5) / (2 * n));
    var x: [4 * n]f64 = undefined;
    for (&x, 0..) |*v, i| v.* = @cos(0.37 * @as(f64, @floatFromInt(i)));

    // frames at hops of n; the middle 2 n samples see two frames each
    var y: [4 * n]f64 = @splat(0);
    for (0..3) |f| {
        var frame: [2 * n]f64 = undefined;
        for (&frame, window, x[f * n ..][0 .. 2 * n]) |*v, w, s| v.* = w * s;
        var coefficients: [n]f64 = undefined;
        try mdct.forward(&coefficients, &frame, scratch);
        try mdct.backward(&frame, &coefficients, scratch);
        for (frame, window, y[f * n ..][0 .. 2 * n]) |v, w, *o| o.* += w * v;
    }
    for (x[n .. 3 * n], y[n .. 3 * n]) |a, b| try testing.expectApproxEqAbs(a, b, 1e-12);
}

test "Plan matches a direct DFT at every float type, mixed radix, the generic pass and Bluestein's" {
    const al = testing.allocator;
    inline for (.{ f16, f32, f64, f80, f128 }) |T| {
        const eps: f64 = @floatCast(std.math.floatEps(T));
        const sizes: []const usize = if (T == f16) &.{ 1, 2, 3, 5, 8, 12, 13, 32 } else &.{ 1, 2, 3, 5, 7, 8, 11, 12, 13, 16, 30, 97, 128, 360, 1009 };
        for (sizes) |n| {
            var p = try Plan(T).init(al, n);
            defer p.deinit(al);
            const scratch = try al.alloc(Complex(T), p.scratchLen());
            defer al.free(scratch);
            const in = try al.alloc(Complex(T), n);
            defer al.free(in);
            const out = try al.alloc(Complex(T), n);
            defer al.free(out);
            const want = try al.alloc(Complex(f64), n);
            defer al.free(want);
            signal(T, in);
            directDft(T, want, in);
            try p.transformInto(testing.io, out, in, scratch, .forward, .{}, .backward);
            const nf: f64 = @floatFromInt(n);
            // the f64 reference bounds the wider types
            try expectNear(T, want, out, @max(8 * eps * nf, 1e-13 * nf));
        }
    }
}

test "Plan: 97 takes the generic pass, 1009 Bluestein's" {
    var a = try Plan(f64).init(testing.allocator, 97);
    defer a.deinit(testing.allocator);
    try testing.expect(a.kind == .pack);
    var b = try Plan(f64).init(testing.allocator, 1009);
    defer b.deinit(testing.allocator);
    try testing.expect(b.kind == .blue);
}

test "Plan: short scratch or data, and a count mismatch" {
    const al = testing.allocator;
    var p = try Plan(f64).init(al, 8);
    defer p.deinit(al);
    const scratch = try al.alloc(Complex(f64), p.scratchLen());
    defer al.free(scratch);
    var data: [16]Complex(f64) = @splat(.init(0, 0));
    try testing.expectError(error.SliceTooShort, p.transform(testing.io, &data, scratch[1..], .forward, .{}, .backward));
    try testing.expectError(error.SliceTooShort, p.transform(testing.io, data[0..7], scratch, .forward, .{}, .backward));
    try testing.expectError(error.SliceTooShort, p.transform(testing.io, &data, scratch, .forward, .{ .count = 3 }, .backward));
    try p.transform(testing.io, &data, scratch, .forward, .{ .count = 2 }, .backward);
    try p.transform(testing.io, &data, scratch, .forward, .{ .count = 0 }, .backward);
    try testing.expectError(error.CountMismatch, p.transformInto(testing.io, &data, &data, scratch, .forward, .{ .in = .{ .count = 2 } }, .backward));
}

test "RealPlan matches the direct DFT's half at every float type, radix 2 to 5, the general passes and Bluestein's" {
    const al = testing.allocator;
    inline for (.{ f16, f32, f64, f80, f128 }) |T| {
        const eps: f64 = @floatCast(std.math.floatEps(T));
        const sizes: []const usize = if (T == f16) &.{ 1, 2, 3, 4, 5, 6, 7, 9, 16, 30 } else &.{ 1, 2, 3, 4, 5, 6, 7, 9, 12, 13, 16, 30, 49, 97, 100, 360, 1009 };
        for (sizes) |n| {
            var p = try RealPlan(T).init(al, n);
            defer p.deinit(al);
            const scratch = try al.alloc(T, p.scratchLen());
            defer al.free(scratch);
            const full = try al.alloc(Complex(T), n);
            defer al.free(full);
            const re = try al.alloc(T, n);
            defer al.free(re);
            const want = try al.alloc(Complex(f64), n);
            defer al.free(want);
            const m = n / 2 + 1;
            const half = try al.alloc(Complex(T), m);
            defer al.free(half);
            const back = try al.alloc(T, n);
            defer al.free(back);
            signal(T, full);
            for (re, full) |*r, *c| {
                r.* = c.re;
                c.im = 0;
            }
            directDft(T, want, full);
            try p.forward(testing.io, half, re, scratch, .{}, .backward);
            const nf: f64 = @floatFromInt(n);
            const tol = @max(8 * eps * nf, 1e-13 * nf);
            try expectNear(T, want[0..m], half, tol);
            try p.backward(testing.io, back, half, scratch, .{}, .backward);
            for (re, back) |a, b| try testing.expectApproxEqAbs(@as(f64, @floatCast(a)), @as(f64, @floatCast(b)), tol);
        }
    }
}

test "RealPlan: 97 takes the general passes, 1009 Bluestein's; strided channels and ortho" {
    const al = testing.allocator;
    var a = try RealPlan(f64).init(al, 97);
    defer a.deinit(al);
    try testing.expect(a.kind == .pack);
    var b = try RealPlan(f64).init(al, 1009);
    defer b.deinit(al);
    try testing.expect(b.kind == .blue);

    // two interleaved real channels into two interleaved half spectra
    const n = 30;
    const m = n / 2 + 1;
    var p = try RealPlan(f64).init(al, n);
    defer p.deinit(al);
    const scratch = try al.alloc(f64, p.scratchLen());
    defer al.free(scratch);
    var x: [2 * n]f64 = undefined;
    for (&x, 0..) |*v, i| v.* = @sin(0.3 * @as(f64, @floatFromInt(i * i)));
    var spec: [2 * m]Complex(f64) = undefined;
    const lay: BatchLayouts = .{ .in = .{ .count = 2, .stride = 2, .distance = 1 }, .out = .{ .count = 2, .stride = 2, .distance = 1 } };
    try p.forward(testing.io, &spec, &x, scratch, lay, .ortho);
    for (0..2) |c| {
        var one: [n]f64 = undefined;
        for (&one, 0..) |*v, i| v.* = x[2 * i + c];
        var want: [m]Complex(f64) = undefined;
        try p.forward(testing.io, &want, &one, scratch, .{}, .ortho);
        for (0..m) |j| {
            try testing.expectEqual(want[j].re, spec[2 * j + c].re);
            try testing.expectEqual(want[j].im, spec[2 * j + c].im);
        }
    }
    var y: [2 * n]f64 = undefined;
    try p.backward(testing.io, &y, &spec, scratch, lay, .ortho);
    for (x, y) |u, v| try testing.expectApproxEqAbs(u, v, 1e-13);
    try testing.expectError(error.SliceTooShort, p.forward(testing.io, spec[0..m], &x, scratch, lay, .backward));
    try testing.expectError(error.CountMismatch, p.backward(testing.io, &y, &spec, scratch, .{ .in = .{ .count = 2 } }, .backward));
}

test "a batch runs as many at once as the scratch holds, with one slot's bits" {
    const al = testing.allocator;
    const n = 24;
    const count = 7;
    var p = try Plan(f64).init(al, n);
    defer p.deinit(al);
    var r = try RealPlan(f64).init(al, n);
    defer r.deinit(al);
    const slots = 4;
    const scratch = try al.alloc(Complex(f64), slots * p.scratchLen());
    defer al.free(scratch);
    const rscratch = try al.alloc(f64, slots * r.scratchLen());
    defer al.free(rscratch);
    var x: [n * count]Complex(f64) = undefined;
    signal(f64, &x);
    var one = x;
    var many = x;
    try p.transform(testing.io, &one, scratch[0..p.scratchLen()], .forward, .{ .count = count }, .backward);
    try p.transform(testing.io, &many, scratch, .forward, .{ .count = count }, .backward);
    try testing.expectEqualSlices(Complex(f64), &one, &many);
    var re: [n * count]f64 = undefined;
    for (&re, x) |*v, c| v.* = c.re;
    var h1: [(n / 2 + 1) * count]Complex(f64) = undefined;
    var h2: [(n / 2 + 1) * count]Complex(f64) = undefined;
    const lay: BatchLayouts = .{ .in = .{ .count = count }, .out = .{ .count = count } };
    try r.forward(testing.io, &h1, &re, rscratch[0..r.scratchLen()], lay, .backward);
    try r.forward(testing.io, &h2, &re, rscratch, lay, .backward);
    try testing.expectEqualSlices(Complex(f64), &h1, &h2);
}
