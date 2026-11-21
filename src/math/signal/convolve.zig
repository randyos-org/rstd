//! Uniformly partitioned convolution (UPOLS: Wefers, "Partitioned
//! convolution algorithms for real-time auralization", 2015, section 5.2):
//! a long filter cut into `block`-sample parts, each convolved by
//! overlap-save in `fft.RealPlan`'s unordered spectra (pffft's) through a
//! frequency-domain delay line. Latency is one block; nothing allocates per
//! block. f32 at pffft's sizes (`2 block = 32 2^a 3^b 5^c`).

const std = @import("std");
const fft = @import("../kernels/fft/__root__.zig");

const Allocator = std.mem.Allocator;

pub const InitError = Allocator.Error || error{
    /// `block` is 0, or the filter is empty.
    ZeroLength,
    /// `2 block` is not one of pffft's real sizes.
    UnsupportedSize,
    /// Not f32: pffft's f64 SIMD is later.
    UnsupportedType,
};

pub const Error = error{SliceTooShort};

pub fn Convolver(comptime T: type) type {
    return struct {
        plan: fft.RealPlan(T),
        block: usize,
        /// Partitions: the filter's spectra, and the delay line's.
        parts: usize,
        /// The delay line's newest spectrum.
        head: usize,
        /// The filter's spectra, then the delay line's, `parts` each; the
        /// input window (the last two blocks), the sum and the plan's scratch.
        buf: []align(16) T,

        const Self = @This();

        /// The filter's spectra, from `filter`; `block` samples per `process`.
        pub fn init(allocator: Allocator, filter: []const T, block: usize) InitError!Self {
            if (block == 0 or filter.len == 0) return error.ZeroLength;
            const n = 2 * block;
            var plan = fft.RealPlan(T).init(allocator, n) catch |e| switch (e) {
                error.OutOfMemory => |x| return x,
                error.ZeroLength, error.NotImplemented => unreachable,
            };
            errdefer plan.deinit(allocator);
            const parts = std.math.divCeil(usize, filter.len, block) catch unreachable;
            const buf = try allocator.alignedAlloc(T, .@"16", 2 * parts * n + 2 * n + plan.scratchLen());
            errdefer allocator.free(buf);
            @memset(buf, 0);
            var self: Self = .{ .plan = plan, .block = block, .parts = parts, .head = 0, .buf = buf };
            // each part zero-padded to `n`, in the sum's slot until transformed
            const pad = self.sum();
            for (0..parts) |p| {
                const piece = filter[p * block .. @min(filter.len, (p + 1) * block)];
                @memcpy(pad[0..piece.len], piece);
                @memset(pad[piece.len..], 0);
                self.plan.forwardUnordered(self.filterSpectrum(p), pad, self.scratch()) catch |e| return switch (e) {
                    error.SliceTooShort => unreachable,
                    else => |x| x,
                };
            }
            return self;
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            self.plan.deinit(allocator);
            allocator.free(self.buf);
            self.* = undefined;
        }

        /// Silence: the delay line and input window cleared, the filter kept.
        pub fn reset(self: *Self) void {
            const n = 2 * self.block;
            @memset(self.buf[self.parts * n ..][0 .. self.parts * n + n], 0);
            self.head = 0;
        }

        /// The next `block` samples of `in` convolved into `out[0..block]`,
        /// which may be `in`: output lags input by one block.
        pub fn process(self: *Self, out: []T, in: []const T) Error!void {
            const b = self.block;
            const n = 2 * b;
            if (in.len < b or out.len < b) return error.SliceTooShort;
            const window = self.inputWindow();
            @memmove(window[0..b], window[b..n]);
            @memcpy(window[b..n], in[0..b]);
            self.head = (self.head + 1) % self.parts;
            self.plan.forwardUnordered(self.lineSpectrum(self.head), window, self.scratch()) catch unreachable;
            // part p meets the input from p blocks ago, in a fixed order
            const acc = self.sum();
            @memset(acc, 0);
            const scale = 1 / @as(T, @floatFromInt(n));
            for (0..self.parts) |p| {
                const age = (self.head + self.parts - p) % self.parts;
                self.plan.convolveAccumulate(acc, self.lineSpectrum(age), self.filterSpectrum(p), scale) catch unreachable;
            }
            self.plan.backwardUnordered(acc, acc, self.scratch()) catch unreachable;
            // overlap-save: the second half is free of the circular wrap
            @memcpy(out[0..b], acc[b..n]);
        }

        fn filterSpectrum(self: Self, p: usize) []T {
            const n = 2 * self.block;
            return self.buf[p * n ..][0..n];
        }

        fn lineSpectrum(self: Self, k: usize) []T {
            const n = 2 * self.block;
            return self.buf[(self.parts + k) * n ..][0..n];
        }

        fn inputWindow(self: Self) []T {
            const n = 2 * self.block;
            return self.buf[2 * self.parts * n ..][0..n];
        }

        fn sum(self: Self) []T {
            const n = 2 * self.block;
            return self.buf[(2 * self.parts + 1) * n ..][0..n];
        }

        fn scratch(self: Self) []T {
            const n = 2 * self.block;
            return self.buf[(2 * self.parts + 2) * n ..];
        }
    };
}

// ---- tests ----

const testing = std.testing;

/// `y[t] = sum h[k] x[t - k]`, in f64.
fn direct(y: []f64, h: []const f32, x: []const f32) void {
    for (y, 0..) |*o, t| {
        var s: f64 = 0;
        for (h, 0..) |hk, k| {
            if (k > t) break;
            s += @as(f64, hk) * x[t - k];
        }
        o.* = s;
    }
}

test "partitioned convolution matches direct convolution" {
    var prng: std.Random.DefaultPrng = .init(31);
    const r = prng.random();
    // a partial last part, and one shorter than a block
    for ([_]usize{ 100, 16, 5 }) |taps| {
        const block = 16;
        const blocks = 12;
        const h = try testing.allocator.alloc(f32, taps);
        defer testing.allocator.free(h);
        for (h) |*v| v.* = r.float(f32) * 2 - 1;
        var x: [blocks * block]f32 = undefined;
        for (&x) |*v| v.* = r.float(f32) * 2 - 1;
        var want: [x.len]f64 = undefined;
        direct(&want, h, &x);

        var c = try Convolver(f32).init(testing.allocator, h, block);
        defer c.deinit(testing.allocator);
        var got: [x.len]f32 = undefined;
        for (0..blocks) |k| try c.process(got[k * block ..][0..block], x[k * block ..][0..block]);
        for (got, want) |g, w| try testing.expectApproxEqAbs(w, g, 1e-5 * @as(f64, @floatFromInt(taps)));

        // after reset, the same input gives the same output, in place
        c.reset();
        var again = x;
        for (0..blocks) |k| try c.process(again[k * block ..][0..block], again[k * block ..][0..block]);
        try testing.expectEqualSlices(f32, &got, &again);
    }
}

test "partitioned convolution: sizes, types and short slices" {
    const h = [_]f32{ 1, 0.5 };
    try testing.expectError(error.ZeroLength, Convolver(f32).init(testing.allocator, &h, 0));
    try testing.expectError(error.ZeroLength, Convolver(f32).init(testing.allocator, &.{}, 16));
    // 2 * 17 is not 32 2^a 3^b 5^c
    try testing.expectError(error.UnsupportedSize, Convolver(f32).init(testing.allocator, &h, 17));
    try testing.expectError(error.UnsupportedType, Convolver(f64).init(testing.allocator, &.{ 1, 0.5 }, 16));

    var c = try Convolver(f32).init(testing.allocator, &h, 48);
    defer c.deinit(testing.allocator);
    var buf: [47]f32 = @splat(0);
    try testing.expectError(error.SliceTooShort, c.process(&buf, &buf));
}
