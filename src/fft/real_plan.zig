//! `fft.RealPlan`: `n` reals to the `n / 2 + 1` values of the half spectrum
//! (NumPy's `rfft`) and back, PocketFFT's `pocketfft_r`
//! (`pocketfft_hdronly.h`, cpp branch, commit 4d865cba, BSD-3, LICENSE):
//! `rfftp/` unless its cost model prefers Bluestein's on complex values,
//! through FFTPACK's halfcomplex order either way. The unordered forms are
//! pffft's, still stubs (N10).

const std = @import("std");
const fft = @import("__root__.zig");
const cmplx = @import("cmplx.zig");
const util = @import("util.zig");
const complex_plan = @import("plan.zig");
const rfftp = @import("rfftp/__root__.zig");
const bluestein = @import("bluestein.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const C = cmplx.C;
const Layouts = fft.Layouts;
const Options = fft.Options;

pub fn RealPlan(comptime T: type) type {
    if (@typeInfo(T) != .float) @compileError("rstd.fft: " ++ @typeName(T) ++ " is not a float type");
    return struct {
        n: usize,
        scratch_len: usize,
        kind: Kind,

        const Self = @This();
        const Kind = union(enum) { pack: rfftp.Rfftp(T), blue: bluestein.Bluestein(T) };

        pub fn init(allocator: Allocator, n: usize) fft.InitError!Self {
            if (n == 0) return error.ZeroLength;
            // `pocketfft_r`'s choice
            const lpf = if (n < 50) 0 else util.largestPrimeFactor(n);
            if (lpf * lpf > n) blue: {
                const good = util.goodSizeComplex(2 * n - 1) orelse break :blue;
                const comp1 = 0.5 * util.costGuess(n);
                // fudge factor that appears to give good overall performance
                const comp2 = 2 * util.costGuess(good) * 1.5;
                if (comp2 >= comp1) break :blue;
                const b = bluestein.Bluestein(T).init(allocator, n) catch |e| switch (e) {
                    error.OutOfMemory => |x| return x,
                    error.TooLarge => break :blue,
                };
                // the gathered reals, then Bluestein's complex scratch as reals
                return .{ .n = n, .scratch_len = n + 2 * b.realScratchLen(), .kind = .{ .blue = b } };
            }
            const p = try rfftp.Rfftp(T).init(allocator, n);
            // the gathered reals and the ping-pong buffer
            return .{ .n = n, .scratch_len = 2 * n, .kind = .{ .pack = p } };
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            switch (self.kind) {
                inline else => |*k| k.deinit(allocator),
            }
            self.* = undefined;
        }

        /// Scratch per transform, as `Plan`'s.
        pub fn scratchLen(self: Self) usize {
            return self.scratch_len;
        }

        pub fn forward(self: *const Self, io: Io, out: []C(T), in: []const T, scratch: []T, layouts: Layouts, options: Options) fft.TransformIntoError!void {
            _ = io;
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            const n = self.n;
            const m = n / 2 + 1;
            if (scratch.len < self.scratch_len) return error.SliceTooShort;
            try complex_plan.covers(in.len, n, layouts.in);
            try complex_plan.covers(out.len, m, layouts.out);
            const din = layouts.in.distance orelse (n * layouts.in.stride);
            const dout = layouts.out.distance orelse (m * layouts.out.stride);
            const buf = scratch[0..n];
            const f = complex_plan.factor(T, n, true, options.norm);
            for (0..layouts.in.count) |k| {
                for (buf, 0..) |*b, i| b.* = in[k * din + i * layouts.in.stride];
                self.run(buf, scratch[n..self.scratch_len], f, true);
                // the halfcomplex order to the half spectrum
                const o = out[k * dout ..];
                const s = layouts.out.stride;
                o[0] = .{ .re = buf[0], .im = 0 };
                for (1..m) |j| o[j * s] = .{ .re = buf[2 * j - 1], .im = if (2 * j < n) buf[2 * j] else 0 };
            }
        }

        /// Reads only the half spectrum; the imaginary parts at 0 (and at
        /// `n / 2` for even `n`) are ignored, as NumPy's `irfft`.
        pub fn backward(self: *const Self, io: Io, out: []T, in: []const C(T), scratch: []T, layouts: Layouts, options: Options) fft.TransformIntoError!void {
            _ = io;
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            const n = self.n;
            const m = n / 2 + 1;
            if (scratch.len < self.scratch_len) return error.SliceTooShort;
            try complex_plan.covers(in.len, m, layouts.in);
            try complex_plan.covers(out.len, n, layouts.out);
            const din = layouts.in.distance orelse (m * layouts.in.stride);
            const dout = layouts.out.distance orelse (n * layouts.out.stride);
            const buf = scratch[0..n];
            const f = complex_plan.factor(T, n, false, options.norm);
            for (0..layouts.in.count) |k| {
                const h = in[k * din ..];
                const s = layouts.in.stride;
                buf[0] = h[0].re;
                for (1..m) |j| {
                    buf[2 * j - 1] = h[j * s].re;
                    if (2 * j < n) buf[2 * j] = h[j * s].im;
                }
                self.run(buf, scratch[n..self.scratch_len], f, false);
                for (buf, 0..) |b, i| out[k * dout + i * layouts.out.stride] = b;
            }
        }

        /// One transform of `buf` in the halfcomplex order, `rest` its space.
        fn run(self: *const Self, buf: []T, rest: []T, f: T, fwd: bool) void {
            switch (self.kind) {
                .pack => |*p| p.exec(buf, rest[0..self.n], f, fwd),
                .blue => |*b| {
                    // the reals' storage, as complex values
                    const bytes = std.mem.sliceAsBytes(rest);
                    const cs: []C(T) = @alignCast(std.mem.bytesAsSlice(C(T), bytes[0 .. bytes.len / @sizeOf(C(T)) * @sizeOf(C(T))]));
                    b.execReal(buf, cs, f, fwd);
                },
            }
        }

        /// The spectrum as `n` reals in the plan's own order (pffft's), for
        /// fast convolution: no reordering either way. Unnormalized.
        pub fn forwardUnordered(self: *const Self, out: []T, in: []const T, scratch: []T) fft.FrameError!void {
            _ = self;
            _ = out;
            _ = in;
            _ = scratch;
            return error.NotImplemented;
        }

        /// `acc += scale * a * b`, element-wise in the unordered spectrum.
        /// `scale` 1 / n makes a round trip a circular convolution.
        pub fn convolveAccumulate(self: *const Self, acc: []T, a: []const T, b: []const T, scale: T) fft.FrameError!void {
            _ = self;
            _ = acc;
            _ = a;
            _ = b;
            _ = scale;
            return error.NotImplemented;
        }

        /// Unnormalized.
        pub fn backwardUnordered(self: *const Self, out: []T, in: []const T, scratch: []T) fft.FrameError!void {
            _ = self;
            _ = out;
            _ = in;
            _ = scratch;
            return error.NotImplemented;
        }
    };
}
