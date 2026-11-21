//! `Plan(f32)` and `RealPlan(f32)` at pffft's sizes (`pffft/`) against
//! PocketFFT's f32 results, the routing, and the unordered forms' errors.

const std = @import("std");
const fft = @import("../__root__.zig");
const cfftp = @import("../cfftp/__root__.zig");
const rfftp = @import("../rfftp/__root__.zig");
const pffft = @import("../pffft/__root__.zig");

const testing = std.testing;
const Complex = std.math.Complex;

fn signal(out: []f32, seed: f32) void {
    for (out, 0..) |*x, i| {
        const t: f32 = @floatFromInt(i);
        x.* = @sin(seed + 0.7 * t) + 0.5 * @cos(0.3 * t * t);
    }
}

test "Plan(f32) at pffft's sizes matches PocketFFT's, both ways" {
    const al = testing.allocator;
    const sizes = [_]usize{ 16, 32, 48, 80, 144, 240, 1024, 4800 };
    for (sizes) |n| {
        var p = try fft.Plan(f32).init(al, n);
        defer p.deinit(al);
        try testing.expect(p.kind == .simd);
        var q = try cfftp.Cfftp(f32).init(al, n);
        defer q.deinit(al);
        const scratch = try al.alloc(Complex(f32), p.scratchLen());
        defer al.free(scratch);
        const x = try al.alloc(Complex(f32), 3 * n);
        defer al.free(x);
        const in = x[0..n];
        const got = x[n .. 2 * n];
        const want = x[2 * n ..];
        signal(std.mem.bytesAsSlice(f32, std.mem.sliceAsBytes(in)), 1);
        const nf: f32 = @floatFromInt(n);
        const tol = 2e-6 * nf;
        inline for (.{ fft.Direction.forward, .backward }) |dir| {
            @memcpy(want, in);
            q.exec(want, scratch[0..n], 1, dir == .forward);
            try p.transformInto(testing.io, got, in, scratch, dir, .{}, .{ .norm = .forward });
            for (want, got) |w, g| {
                const f: f32 = if (dir == .forward) 1 / nf else 1;
                try testing.expectApproxEqAbs(w.re * f, g.re, tol);
                try testing.expectApproxEqAbs(w.im * f, g.im, tol);
            }
        }
        // a round trip
        try p.transformInto(testing.io, got, in, scratch, .forward, .{}, .{});
        try p.transform(testing.io, got, scratch, .backward, .{}, .{});
        for (in, got) |a, b| {
            try testing.expectApproxEqAbs(a.re, b.re, 1e-5);
            try testing.expectApproxEqAbs(a.im, b.im, 1e-5);
        }
    }
}

test "RealPlan(f32) at pffft's sizes matches PocketFFT's, both ways" {
    const al = testing.allocator;
    const sizes = [_]usize{ 32, 64, 96, 160, 480, 1024, 2400 };
    for (sizes) |n| {
        var p = try fft.RealPlan(f32).init(al, n);
        defer p.deinit(al);
        try testing.expect(p.kind == .simd);
        var q = try rfftp.Rfftp(f32).init(al, n);
        defer q.deinit(al);
        const m = n / 2 + 1;
        const scratch = try al.alloc(f32, p.scratchLen());
        defer al.free(scratch);
        const re = try al.alloc(f32, 3 * n);
        defer al.free(re);
        const half = try al.alloc(Complex(f32), m);
        defer al.free(half);
        const x = re[0..n];
        const hc = re[n .. 2 * n];
        const back = re[2 * n ..];
        signal(x, 2);
        @memcpy(hc, x);
        q.exec(hc, scratch[0..n], 1, true);
        try p.forward(testing.io, half, x, scratch, .{}, .{});
        const tol = 2e-6 * @as(f32, @floatFromInt(n));
        try testing.expectApproxEqAbs(hc[0], half[0].re, tol);
        try testing.expectEqual(@as(f32, 0), half[0].im);
        for (1..n / 2) |j| {
            try testing.expectApproxEqAbs(hc[2 * j - 1], half[j].re, tol);
            try testing.expectApproxEqAbs(hc[2 * j], half[j].im, tol);
        }
        try testing.expectApproxEqAbs(hc[n - 1], half[n / 2].re, tol);
        try testing.expectEqual(@as(f32, 0), half[n / 2].im);

        // backward: PocketFFT's from the same halfcomplex values
        try p.backward(testing.io, back, half, scratch, .{}, .{ .norm = .forward });
        q.exec(hc, scratch[0..n], 1, false);
        for (hc, back) |w, g| try testing.expectApproxEqAbs(w, g, tol);
        try p.backward(testing.io, back, half, scratch, .{}, .{});
        for (x, back) |a, b| try testing.expectApproxEqAbs(a, b, 1e-5);
    }
}

test "pffft's routing: f32 at its sizes only" {
    const al = testing.allocator;
    const Case = struct { n: usize, simd: bool };
    for ([_]Case{ .{ .n = 16, .simd = true }, .{ .n = 8, .simd = false }, .{ .n = 24, .simd = false }, .{ .n = 112, .simd = false }, .{ .n = 4800, .simd = true } }) |c| {
        var p = try fft.Plan(f32).init(al, c.n);
        defer p.deinit(al);
        try testing.expectEqual(c.simd, p.kind == .simd);
    }
    for ([_]Case{ .{ .n = 32, .simd = true }, .{ .n = 16, .simd = false }, .{ .n = 48, .simd = false }, .{ .n = 224, .simd = false }, .{ .n = 7200, .simd = true } }) |c| {
        var p = try fft.RealPlan(f32).init(al, c.n);
        defer p.deinit(al);
        try testing.expectEqual(c.simd, p.kind == .simd);
    }
    var d = try fft.Plan(f64).init(al, 64);
    defer d.deinit(al);
    try testing.expect(d.kind == .pack);
    var r = try fft.RealPlan(f64).init(al, 64);
    defer r.deinit(al);
    try testing.expect(r.kind == .pack);
}

test "RealPlan unordered: errors for other types, sizes and short slices" {
    const al = testing.allocator;
    var buf: [128]f32 = @splat(0);
    var d = try fft.RealPlan(f64).init(al, 64);
    defer d.deinit(al);
    var dbuf: [128]f64 = @splat(0);
    try testing.expectError(error.UnsupportedType, d.forwardUnordered(dbuf[0..64], dbuf[0..64], &dbuf));
    try testing.expectError(error.UnsupportedType, d.convolveAccumulate(dbuf[0..64], dbuf[0..64], dbuf[0..64], 1));
    try testing.expectError(error.UnsupportedType, d.backwardUnordered(dbuf[0..64], dbuf[0..64], &dbuf));

    var s = try fft.RealPlan(f32).init(al, 48);
    defer s.deinit(al);
    try testing.expectError(error.UnsupportedSize, s.forwardUnordered(buf[0..48], buf[0..48], &buf));
    try testing.expectError(error.UnsupportedSize, s.convolveAccumulate(buf[0..48], buf[0..48], buf[0..48], 1));
    try testing.expectError(error.UnsupportedSize, s.backwardUnordered(buf[0..48], buf[0..48], &buf));

    var p = try fft.RealPlan(f32).init(al, 32);
    defer p.deinit(al);
    try testing.expectError(error.SliceTooShort, p.forwardUnordered(buf[0..31], buf[32..64], buf[64..128]));
    try testing.expectError(error.SliceTooShort, p.forwardUnordered(buf[0..32], buf[32..64], buf[64..127]));
    try testing.expectError(error.SliceTooShort, p.convolveAccumulate(buf[0..32], buf[32..63], buf[64..96], 1));
    try testing.expectError(error.SliceTooShort, p.backwardUnordered(buf[0..32], buf[32..63], buf[64..128]));
}

test "RealPlan unordered: an in-place round trip, and convolution at a 3 and 5 size" {
    const al = testing.allocator;
    const n = 480;
    var p = try fft.RealPlan(f32).init(al, n);
    defer p.deinit(al);
    const scratch = try al.alloc(f32, p.scratchLen());
    defer al.free(scratch);
    var a: [n]f32 = undefined;
    var b: [n]f32 = undefined;
    signal(&a, 3);
    for (&b, 0..) |*y, i| y.* = if (i < 7) @as(f32, @floatFromInt(i + 1)) else 0;
    var want: [n]f32 = @splat(0);
    for (&want, 0..) |*w, k| {
        for (0..n) |j| w.* += a[j] * b[(k + n - j) % n];
    }
    var x = a;
    try p.forwardUnordered(&x, &x, scratch);
    try p.backwardUnordered(&x, &x, scratch);
    for (a, x) |u, v| try testing.expectApproxEqAbs(u * n, v, 1e-3);

    var fa: [n]f32 = undefined;
    var fb: [n]f32 = undefined;
    try p.forwardUnordered(&fa, &a, scratch);
    try p.forwardUnordered(&fb, &b, scratch);
    var acc: [n]f32 = @splat(0);
    try p.convolveAccumulate(&acc, &fa, &fb, 0.5 / @as(f32, n));
    try p.convolveAccumulate(&acc, &fa, &fb, 0.5 / @as(f32, n)); // accumulates
    try p.backwardUnordered(&acc, &acc, scratch);
    for (want, acc) |w, g| try testing.expectApproxEqAbs(w, g, 1e-3);
}

test "pffft's 8-lane steps match its 4-lane ones, ordered, unordered and convolving" {
    const al = testing.allocator;
    const Case = struct { n: usize, kind: pffft.Kind };
    const cases = [_]Case{
        .{ .n = 128, .kind = .real },    .{ .n = 384, .kind = .real },     .{ .n = 640, .kind = .real },
        .{ .n = 3840, .kind = .real },   .{ .n = 64, .kind = .complex },   .{ .n = 192, .kind = .complex },
        .{ .n = 320, .kind = .complex }, .{ .n = 4800, .kind = .complex },
    };
    for (cases) |c| {
        var s4 = try pffft.Setup.initLanes(al, c.n, c.kind, 4);
        defer s4.deinit(al);
        var s8 = try pffft.Setup.initLanes(al, c.n, c.kind, 8);
        defer s8.deinit(al);
        const len = if (c.kind == .real) c.n else 2 * c.n;
        const buf = try al.alloc(f32, 7 * len);
        defer al.free(buf);
        const x = buf[0..len];
        const y4 = buf[len .. 2 * len];
        const y8 = buf[2 * len .. 3 * len];
        const w = buf[3 * len .. 4 * len];
        const z = buf[4 * len .. 5 * len];
        const acc4 = buf[5 * len .. 6 * len];
        const acc8 = buf[6 * len ..];
        signal(x, 2);
        const nf: f32 = @floatFromInt(c.n);
        const tol = 4e-6 * nf;
        // ordered, both ways; backward is unscaled
        s4.transform(x.ptr, y4.ptr, w.ptr, true, true);
        s8.transform(x.ptr, y8.ptr, w.ptr, true, true);
        for (y4, y8) |a, b| try testing.expectApproxEqAbs(a, b, tol);
        s8.transform(y8.ptr, z.ptr, w.ptr, false, true);
        for (x, z) |a, b| try testing.expectApproxEqAbs(nf * a, b, tol);
        // unordered: a round trip, then a circular convolution
        s8.transform(x.ptr, y8.ptr, w.ptr, true, false);
        s8.transform(y8.ptr, z.ptr, w.ptr, false, false);
        for (x, z) |a, b| try testing.expectApproxEqAbs(nf * a, b, tol);
        s4.transform(x.ptr, y4.ptr, w.ptr, true, false);
        @memset(acc4, 0);
        @memset(acc8, 0);
        s4.convolveAccumulate(y4.ptr, y4.ptr, acc4.ptr, 1 / nf);
        s8.convolveAccumulate(y8.ptr, y8.ptr, acc8.ptr, 1 / nf);
        s4.transform(acc4.ptr, y4.ptr, w.ptr, false, false);
        s8.transform(acc8.ptr, y8.ptr, w.ptr, false, false);
        for (y4, y8) |a, b| try testing.expectApproxEqAbs(a, b, 4 * tol);
    }
}
