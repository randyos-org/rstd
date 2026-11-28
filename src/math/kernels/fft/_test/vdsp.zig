//! Tests: `Plan` and `RealPlan` on vDSP (`vdsp.zig`) against rstd's own
//! transforms, the lengths vDSP takes, and the switch.

const std = @import("std");
const sysinfo = @import("builtin");
const fft = @import("../__root__.zig");
const vdsp = @import("../vdsp.zig");
const pffft = @import("../pffft/__root__.zig");

const testing = std.testing;
const Complex = std.math.Complex;
const Io = std.Io;

/// A batch of 3 transforms, every other element, 1 apart.
const layout: fft.BatchLayout = .{ .count = 3, .stride = 2 };

fn fill(comptime T: type, x: []T, seed: u64) void {
    var rng = std.Random.DefaultPrng.init(seed);
    for (x) |*e| e.* = rng.random().float(T) * 2 - 1;
}

/// `got` within `n` eps of `want`, relative to the largest entry.
fn near(comptime T: type, want: []const T, got: []const T, n: usize) !void {
    var big: T = 0;
    for (want) |w| big = @max(big, @abs(w));
    const tol = 8 * std.math.floatEps(T) * @as(T, @floatFromInt(n)) * @max(big, 1);
    for (want, got) |w, g| try testing.expectApproxEqAbs(w, g, tol);
}

fn floats(comptime T: type, v: []Complex(T)) []T {
    return @as([*]T, @ptrCast(v.ptr))[0 .. 2 * v.len];
}

/// The two plans of `P` at `n`: vDSP's (where it takes `n`, faster or not)
/// and rstd's own.
fn both(comptime P: type, n: usize) ![2]P {
    const al = testing.allocator;
    _ = fft.useAccelerate(false);
    var own = try P.init(al, n);
    errdefer own.deinit(al);
    _ = fft.useAccelerate(true);
    fft.useAccelerateEverywhere(true);
    defer fft.useAccelerateEverywhere(false);
    const routed = try P.init(al, n);
    try testing.expect(!own.routed());
    return .{ routed, own };
}

test "vDSP's lengths: f 2^k, f in 1, 3, 5, 15, k at least 3 complex, 4 real" {
    for ([_]usize{ 8, 16, 24, 40, 48, 120, 1 << 20, 15 << 10 }) |n| try testing.expect(vdsp.validLength(n, false));
    for ([_]usize{ 1, 2, 4, 6, 12, 7 << 4, 9 << 4, 25 << 4, 100 }) |n| try testing.expect(!vdsp.validLength(n, false));
    for ([_]usize{ 16, 48, 80, 240 }) |n| try testing.expect(vdsp.validLength(n, true));
    for ([_]usize{ 8, 24, 40, 120 }) |n| try testing.expect(!vdsp.validLength(n, true));
}

test "Plan on vDSP matches rstd's own: every type, normalization and direction, a strided batch" {
    const al = testing.allocator;
    inline for (.{ f32, f64 }) |T| for ([_]usize{ 8, 16, 24, 40, 48, 120, 1024, 3 << 10, 15 << 6 }) |n| {
        var ps = try both(fft.Plan(T), n);
        defer for (&ps) |*p| p.deinit(al);
        // vDSP comes before pffft
        try testing.expectEqual(vdsp.available, ps[0].routed());
        const len = layout.count * n * layout.stride;
        const in = try al.alloc(Complex(T), len);
        defer al.free(in);
        fill(T, floats(T, in), n);
        var out: [2][]Complex(T) = undefined;
        var scratch: [2][]Complex(T) = undefined;
        for (0..2) |i| {
            out[i] = try al.alloc(Complex(T), len);
            @memset(out[i], .init(0, 0));
            scratch[i] = try al.alloc(Complex(T), ps[i].scratchLen());
        }
        defer for (0..2) |i| {
            al.free(out[i]);
            al.free(scratch[i]);
        };
        inline for (.{ fft.Direction.forward, .backward }) |dir| for ([_]fft.Normalization{ .backward, .ortho, .forward }) |norm| {
            for (0..2) |i| try ps[i].transformInto(testing.io, out[i], in, scratch[i], dir, .{ .in = layout, .out = layout }, norm);
            try near(T, floats(T, out[1]), floats(T, out[0]), n);
        };
    };
}

test "RealPlan on vDSP matches rstd's own: every type and normalization, both ways" {
    const al = testing.allocator;
    inline for (.{ f32, f64 }) |T| for ([_]usize{ 16, 48, 80, 240, 1024, 3 << 10, 15 << 6 }) |n| {
        var ps = try both(fft.RealPlan(T), n);
        defer for (&ps) |*p| p.deinit(al);
        // `f32` at pffft's sizes keeps pffft, which alone has the unordered forms
        const pffft_size = T == f32 and pffft.validSize(n, .real);
        try testing.expectEqual(vdsp.available and !pffft_size, ps[0].routed());
        const m = n / 2 + 1;
        const x = try al.alloc(T, n);
        defer al.free(x);
        fill(T, x, n + 1);
        var half: [2][]Complex(T) = undefined;
        var back: [2][]T = undefined;
        var scratch: [2][]T = undefined;
        for (0..2) |i| {
            half[i] = try al.alloc(Complex(T), m);
            back[i] = try al.alloc(T, n);
            scratch[i] = try al.alloc(T, ps[i].scratchLen());
        }
        defer for (0..2) |i| {
            al.free(half[i]);
            al.free(back[i]);
            al.free(scratch[i]);
        };
        for ([_]fft.Normalization{ .backward, .ortho, .forward }) |norm| {
            for (0..2) |i| try ps[i].forward(testing.io, half[i], x, scratch[i], .{}, norm);
            try near(T, floats(T, half[1]), floats(T, half[0]), n);
            for (0..2) |i| try ps[i].backward(testing.io, back[i], half[1], scratch[i], .{}, norm);
            try near(T, back[1], back[0], n);
        }
    };
}

test "a routed batch: the same bits on Io.failing and on testing.io with many slots" {
    const al = testing.allocator;
    _ = fft.useAccelerate(true);
    var p = try fft.Plan(f64).init(al, 256);
    defer p.deinit(al);
    const count = 8;
    const in = try al.alloc(Complex(f64), count * 256);
    defer al.free(in);
    fill(f64, floats(f64, in), 9);
    const a = try al.alloc(Complex(f64), in.len);
    defer al.free(a);
    const b = try al.alloc(Complex(f64), in.len);
    defer al.free(b);
    const s1 = try al.alloc(Complex(f64), p.scratchLen());
    defer al.free(s1);
    const s8 = try al.alloc(Complex(f64), count * p.scratchLen());
    defer al.free(s8);
    const l: fft.BatchLayout = .{ .count = count };
    try p.transformInto(Io.failing, a, in, s1, .forward, .{ .in = l, .out = l }, .backward);
    try p.transformInto(testing.io, b, in, s8, .forward, .{ .in = l, .out = l }, .backward);
    try testing.expectEqualSlices(u8, std.mem.sliceAsBytes(a), std.mem.sliceAsBytes(b));
}

test "useAccelerate(false) keeps rstd's own transforms, and other lengths never route" {
    const al = testing.allocator;
    try testing.expect(!fft.useAccelerate(false));
    var p = try fft.Plan(f64).init(al, 64);
    defer p.deinit(al);
    try testing.expect(!p.routed());
    try testing.expectEqual(vdsp.available, fft.useAccelerate(true));
    var q = try fft.Plan(f64).init(al, 7 << 4);
    defer q.deinit(al);
    try testing.expect(!q.routed());
    var h = try fft.Plan(f16).init(al, 64);
    defer h.deinit(al);
    try testing.expect(!h.routed());
}

test "plans take vDSP where it ran faster: per architecture" {
    const al = testing.allocator;
    const Case = struct { real: bool, double: bool, n: usize, arm: bool, x86: bool };
    const cases = [_]Case{
        .{ .real = false, .double = false, .n = 1024, .arm = false, .x86 = false },
        .{ .real = false, .double = false, .n = 3840, .arm = true, .x86 = false },
        .{ .real = false, .double = false, .n = 3072, .arm = false, .x86 = false },
        .{ .real = false, .double = false, .n = 120, .arm = true, .x86 = true },
        .{ .real = false, .double = false, .n = 3 << 15, .arm = false, .x86 = false },
        .{ .real = false, .double = false, .n = 1 << 17, .arm = true, .x86 = false },
        .{ .real = false, .double = true, .n = 4096, .arm = true, .x86 = true },
        .{ .real = false, .double = true, .n = 3 << 10, .arm = false, .x86 = true },
        .{ .real = true, .double = false, .n = 16, .arm = false, .x86 = true },
        .{ .real = true, .double = false, .n = 240, .arm = true, .x86 = true },
        .{ .real = true, .double = true, .n = 1024, .arm = true, .x86 = true },
        .{ .real = true, .double = true, .n = 1 << 15, .arm = false, .x86 = true },
        .{ .real = true, .double = true, .n = 15 << 4, .arm = false, .x86 = true },
    };
    _ = fft.useAccelerate(true);
    inline for (.{ false, true }) |real| inline for (.{ f32, f64 }) |T| for (cases) |c| {
        if (c.real != real or c.double != (T == f64)) continue;
        const want = if (sysinfo.cpu.arch == .aarch64) c.arm else c.x86;
        try testing.expectEqual(want, vdsp.faster(T, real, c.n));
        var p = try (if (real) fft.RealPlan(T) else fft.Plan(T)).init(al, c.n);
        defer p.deinit(al);
        try testing.expectEqual(vdsp.available and want, p.routed());
    };
}

test "default f32 plans keep pffft's accuracy around 65536: within 8 eps of an f64 transform" {
    // vDSP's arm64 transform past 65536 at 3, 5, 15 times 2^k was 18 to 24 (#870)
    const al = testing.allocator;
    _ = fft.useAccelerate(true);
    for ([_]usize{ 15 << 12, 1 << 17, 3 << 15, 5 << 14 }) |n| {
        const x = try al.alloc(Complex(f32), n);
        defer al.free(x);
        fill(f32, floats(f32, x), n);
        const want = try al.alloc(Complex(f64), n);
        defer al.free(want);
        for (x, want) |a, *b| b.* = .init(a.re, a.im);
        var pd = try fft.Plan(f64).init(al, n);
        defer pd.deinit(al);
        const sd = try al.alloc(Complex(f64), pd.scratchLen());
        defer al.free(sd);
        try pd.transform(testing.io, want, sd, .forward, .{}, .backward);
        var p = try fft.Plan(f32).init(al, n);
        defer p.deinit(al);
        const s = try al.alloc(Complex(f32), p.scratchLen());
        defer al.free(s);
        try p.transform(testing.io, x, s, .forward, .{}, .backward);
        var big: f64 = 0;
        var worst: f64 = 0;
        for (x, want) |a, b| {
            big = @max(big, @abs(b.re), @abs(b.im));
            worst = @max(worst, @abs(a.re - b.re), @abs(a.im - b.im));
        }
        try testing.expect(worst <= 8 * std.math.floatEps(f32) * big);
    }
}
