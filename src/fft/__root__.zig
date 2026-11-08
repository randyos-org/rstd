//! Fast Fourier transforms (rstd's design.md, "FFT"): PocketFFT's algorithms
//! in Zig, pffft's SIMD. Stubs for now: bodies return `error.NotImplemented`,
//! their tests skip. A plan depends only on `n`, the element type and the
//! core; immutable, it serves every thread. Transforms allocate nothing.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;
const Io = std.Io;

/// `forward` is `exp(-2πi jk/n)`, as NumPy and FFTW.
pub const Direction = enum { forward, backward };

/// Which direction divides by `n`; `ortho` divides both by `√n`.
pub const Norm = enum { backward, ortho, forward };

/// NumPy's default: a round trip gives back its input.
pub const Options = struct { norm: Norm = .backward };

/// A batch of transforms in a slice: transform `k`'s element `i` is at
/// `k * distance + i * stride`. `distance` null is `len * stride`, `len`
/// being that slice's elements per transform (`n`, or `n / 2 + 1` for a half
/// spectrum). Transforms must not overlap.
pub const Layout = struct { count: usize = 1, stride: usize = 1, distance: ?usize = null };

/// An out-of-place batch: both sides' counts must match.
pub const Layouts = struct { in: Layout = .{}, out: Layout = .{} };

pub const InitError = error{ OutOfMemory, ZeroLength, NotImplemented };
pub const TransformError = error{ SliceTooShort, Canceled, NotImplemented };
pub const TransformIntoError = TransformError || error{CountMismatch};

/// One transform at a time, no `io`: audio frames and the unordered forms.
pub const FrameError = error{ SliceTooShort, NotImplemented };

/// Complex, any `n`: PocketFFT's mixed radix (2, 3, 4, 5, 7, 11), Bluestein's
/// for large prime factors.
pub fn Plan(comptime T: type) type {
    comptime checkFloat(T);
    return struct {
        n: usize,
        scratch_len: usize,

        const Self = @This();
        const C = Complex(T);

        pub fn init(allocator: Allocator, n: usize) InitError!Self {
            _ = allocator;
            if (n == 0) return error.ZeroLength;
            return error.NotImplemented;
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            _ = allocator;
            self.* = undefined;
        }

        /// Scratch per transform, any layout. A batch runs as many transforms
        /// at once as the scratch holds, so the bits never depend on it.
        pub fn scratchLen(self: Self) usize {
            return self.scratch_len;
        }

        pub fn transform(self: *const Self, io: Io, data: []C, scratch: []C, direction: Direction, layout: Layout, options: Options) TransformError!void {
            _ = self;
            _ = io;
            _ = data;
            _ = scratch;
            _ = direction;
            _ = layout;
            _ = options;
            return error.NotImplemented;
        }

        /// `out` and `in` must not overlap.
        pub fn transformInto(self: *const Self, io: Io, out: []C, in: []const C, scratch: []C, direction: Direction, layouts: Layouts, options: Options) TransformIntoError!void {
            _ = self;
            _ = io;
            _ = out;
            _ = in;
            _ = scratch;
            _ = direction;
            _ = options;
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            return error.NotImplemented;
        }
    };
}

/// Real input: `n` reals to the `n / 2 + 1` values of the half spectrum
/// (NumPy's `rfft`), and back.
pub fn RealPlan(comptime T: type) type {
    comptime checkFloat(T);
    return struct {
        n: usize,
        scratch_len: usize,

        const Self = @This();
        const C = Complex(T);

        pub fn init(allocator: Allocator, n: usize) InitError!Self {
            _ = allocator;
            if (n == 0) return error.ZeroLength;
            return error.NotImplemented;
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            _ = allocator;
            self.* = undefined;
        }

        /// Scratch per transform, as `Plan`'s.
        pub fn scratchLen(self: Self) usize {
            return self.scratch_len;
        }

        pub fn forward(self: *const Self, io: Io, out: []C, in: []const T, scratch: []T, layouts: Layouts, options: Options) TransformIntoError!void {
            _ = self;
            _ = io;
            _ = out;
            _ = in;
            _ = scratch;
            _ = options;
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            return error.NotImplemented;
        }

        /// Reads only the half spectrum; the imaginary parts at 0 (and at
        /// `n / 2` for even `n`) are ignored, as NumPy's `irfft`.
        pub fn backward(self: *const Self, io: Io, out: []T, in: []const C, scratch: []T, layouts: Layouts, options: Options) TransformIntoError!void {
            _ = self;
            _ = io;
            _ = out;
            _ = in;
            _ = scratch;
            _ = options;
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            return error.NotImplemented;
        }

        /// The spectrum as `n` reals in the plan's own order (pffft's), for
        /// fast convolution: no reordering either way. Unnormalized.
        pub fn forwardUnordered(self: *const Self, out: []T, in: []const T, scratch: []T) FrameError!void {
            _ = self;
            _ = out;
            _ = in;
            _ = scratch;
            return error.NotImplemented;
        }

        /// `acc += scale * a * b`, element-wise in the unordered spectrum.
        /// `scale` 1 / n makes a round trip a circular convolution.
        pub fn convolveAccumulate(self: *const Self, acc: []T, a: []const T, b: []const T, scale: T) FrameError!void {
            _ = self;
            _ = acc;
            _ = a;
            _ = b;
            _ = scale;
            return error.NotImplemented;
        }

        /// Unnormalized.
        pub fn backwardUnordered(self: *const Self, out: []T, in: []const T, scratch: []T) FrameError!void {
            _ = self;
            _ = out;
            _ = in;
            _ = scratch;
            return error.NotImplemented;
        }
    };
}

/// `Plan` for a comptime `n`: factors and twiddles at comptime, unrolled, no
/// allocator; the same calls on a value `.{}`.
pub fn StaticPlan(comptime n: usize, comptime T: type) type {
    comptime checkFloat(T);
    if (n == 0) @compileError("rstd.fft.StaticPlan: n is 0");
    return struct {
        const Self = @This();
        const C = Complex(T);

        /// Stub: PocketFFT's ping-pong buffer; Bluestein's sizes need more.
        pub fn scratchLen(_: Self) usize {
            return n;
        }

        pub fn transform(self: Self, io: Io, data: []C, scratch: []C, direction: Direction, layout: Layout, options: Options) TransformError!void {
            _ = self;
            _ = io;
            _ = data;
            _ = scratch;
            _ = direction;
            _ = layout;
            _ = options;
            return error.NotImplemented;
        }

        pub fn transformInto(self: Self, io: Io, out: []C, in: []const C, scratch: []C, direction: Direction, layouts: Layouts, options: Options) TransformIntoError!void {
            _ = self;
            _ = io;
            _ = out;
            _ = in;
            _ = scratch;
            _ = direction;
            _ = options;
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            return error.NotImplemented;
        }
    };
}

/// `RealPlan` for a comptime `n`, as `StaticPlan`.
pub fn StaticRealPlan(comptime n: usize, comptime T: type) type {
    comptime checkFloat(T);
    if (n == 0) @compileError("rstd.fft.StaticRealPlan: n is 0");
    return struct {
        const Self = @This();
        const C = Complex(T);

        /// Stub, as `StaticPlan`'s.
        pub fn scratchLen(_: Self) usize {
            return n;
        }

        pub fn forward(self: Self, io: Io, out: []C, in: []const T, scratch: []T, layouts: Layouts, options: Options) TransformIntoError!void {
            _ = self;
            _ = io;
            _ = out;
            _ = in;
            _ = scratch;
            _ = options;
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            return error.NotImplemented;
        }

        pub fn backward(self: Self, io: Io, out: []T, in: []const C, scratch: []T, layouts: Layouts, options: Options) TransformIntoError!void {
            _ = self;
            _ = io;
            _ = out;
            _ = in;
            _ = scratch;
            _ = options;
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            return error.NotImplemented;
        }
    };
}

/// DCT and DST types 1 to 4, as SciPy's `dct` and `dst`: `.backward` leaves
/// the forward unscaled, `.ortho` makes the matrix orthogonal. A type 2's
/// inverse is the type 3, and the reverse.
pub const TrigKind = enum { dct1, dct2, dct3, dct4, dst1, dst2, dst3, dst4 };

pub fn TrigPlan(comptime T: type) type {
    comptime checkFloat(T);
    return struct {
        n: usize,
        kind: TrigKind,
        scratch_len: usize,

        const Self = @This();

        pub fn init(allocator: Allocator, n: usize, kind: TrigKind) InitError!Self {
            _ = allocator;
            _ = kind;
            if (n == 0) return error.ZeroLength;
            return error.NotImplemented;
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            _ = allocator;
            self.* = undefined;
        }

        pub fn scratchLen(self: Self) usize {
            return self.scratch_len;
        }

        pub fn transform(self: *const Self, io: Io, data: []T, scratch: []T, layout: Layout, options: Options) TransformError!void {
            _ = self;
            _ = io;
            _ = data;
            _ = scratch;
            _ = layout;
            _ = options;
            return error.NotImplemented;
        }
    };
}

/// The MDCT (audio codecs) on a DCT-IV: `2 n` samples to `n` coefficients
/// and back. The caller windows; with a Princen-Bradley window (sine, KBD),
/// windowed overlap-add at hop `n` gives back the input (TDAC).
pub fn Mdct(comptime T: type) type {
    comptime checkFloat(T);
    return struct {
        n: usize,
        scratch_len: usize,

        const Self = @This();

        pub fn init(allocator: Allocator, n: usize) InitError!Self {
            _ = allocator;
            if (n == 0) return error.ZeroLength;
            return error.NotImplemented;
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            _ = allocator;
            self.* = undefined;
        }

        pub fn scratchLen(self: Self) usize {
            return self.scratch_len;
        }

        /// `out[0..n]` from `in[0 .. 2 n]`.
        pub fn forward(self: *const Self, out: []T, in: []const T, scratch: []T) FrameError!void {
            _ = self;
            _ = out;
            _ = in;
            _ = scratch;
            return error.NotImplemented;
        }

        /// `out[0 .. 2 n]` from `in[0..n]`, still aliased; overlap-add removes it.
        pub fn backward(self: *const Self, out: []T, in: []const T, scratch: []T) FrameError!void {
            _ = self;
            _ = out;
            _ = in;
            _ = scratch;
            return error.NotImplemented;
        }
    };
}

fn checkFloat(comptime T: type) void {
    if (@typeInfo(T) != .float) @compileError("rstd.fft: " ++ @typeName(T) ++ " is not a float type");
}

// ---- tests ----

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
            acc = acc.add(w.mul(.init(x.re, x.im)));
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
    try testing.expectError(error.ZeroLength, TrigPlan(f64).init(testing.allocator, 0, .dct2));
    try testing.expectError(error.ZeroLength, Mdct(f32).init(testing.allocator, 0));
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
        try plan.transformInto(testing.io, out[0..n], in[0..n], scratch, .forward, .{}, .{});
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

        try plan.transform(testing.io, y, scratch, .forward, .{}, .{});
        try plan.transform(testing.io, y, scratch, .backward, .{}, .{});
        for (x, y) |a, b| {
            try testing.expectApproxEqAbs(a.re, b.re, 1e-12);
            try testing.expectApproxEqAbs(a.im, b.im, 1e-12);
        }

        // Parseval
        try plan.transform(testing.io, y, scratch, .forward, .{}, .{ .norm = .ortho });
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
    try plan.transform(testing.io, &both, scratch, .forward, .{ .count = 2, .stride = 2, .distance = 1 }, .{});
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
        try plan.forward(testing.io, half[0..m], re[0..n], scratch, .{}, .{});
        try expectNear(f64, want[0..m], half[0..m], 1e-12);
        try plan.backward(testing.io, back[0..n], half[0..m], scratch, .{}, .{});
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
    try plan.forwardUnordered(&fa, &a, scratch);
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
    try plan.transformInto(testing.io, &want, &in, scratch, .forward, .{}, .{});
    try static.transformInto(testing.io, &got, &in, scratch, .forward, .{}, .{});
    try testing.expectEqualSlices(Complex(f64), &want, &got); // dispatch by size alone: same bits
}

test "TrigPlan: ortho DCT-II then DCT-III, and DCT-IV twice, give back the input" {
    const n = 10;
    const pairs = [_][2]TrigKind{ .{ .dct2, .dct3 }, .{ .dst2, .dst3 }, .{ .dct4, .dct4 }, .{ .dst4, .dst4 } };
    for (pairs) |pair| {
        var there = TrigPlan(f64).init(testing.allocator, n, pair[0]) catch |e| return skip(e);
        defer there.deinit(testing.allocator);
        var back = try TrigPlan(f64).init(testing.allocator, n, pair[1]);
        defer back.deinit(testing.allocator);
        const scratch = try testing.allocator.alloc(f64, @max(there.scratchLen(), back.scratchLen()));
        defer testing.allocator.free(scratch);

        var x: [n]f64 = undefined;
        for (&x, 0..) |*v, i| v.* = @sin(1 + 0.7 * @as(f64, @floatFromInt(i)));
        var y = x;
        try there.transform(testing.io, &y, scratch, .{}, .{ .norm = .ortho });
        try back.transform(testing.io, &y, scratch, .{}, .{ .norm = .ortho });
        for (x, y) |a, b| try testing.expectApproxEqAbs(a, b, 1e-12);
    }
}

test "Mdct: TDAC with a sine window" {
    const n = 16;
    var mdct = Mdct(f64).init(testing.allocator, n) catch |e| return skip(e);
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
