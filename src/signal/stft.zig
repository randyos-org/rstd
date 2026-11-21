//! The short-time Fourier transform over `fft.RealPlan`, and its inverse by
//! weighted overlap-add (MATLAB's `stft` and `istft` with `'wola'`): frame
//! `k` starts at `k * hop`, is windowed, zero-padded to `nfft` and
//! transformed to its `nfft / 2 + 1` half spectrum. Nothing allocates per
//! frame; the plan, the window and the scratch are made up front.

const std = @import("std");
const fft = @import("../math/kernels/fft/__root__.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Complex = std.math.Complex;

pub const InitError = Allocator.Error || error{
    /// `nfft` or `hop` is 0, or the window is empty.
    ZeroLength,
    /// The window is longer than `nfft`.
    WindowTooLong,
};

pub const Error = error{ SliceTooShort, Canceled };

pub fn Stft(comptime T: type) type {
    return struct {
        plan: fft.RealPlan(T),
        window: []T,
        hop: usize,
        nfft: usize,

        const Self = @This();

        /// `window` (copied, at most `nfft` points), `hop` samples between
        /// frames.
        pub fn init(allocator: Allocator, window: []const T, hop: usize, nfft: usize) InitError!Self {
            if (window.len == 0 or hop == 0 or nfft == 0) return error.ZeroLength;
            if (window.len > nfft) return error.WindowTooLong;
            var plan = fft.RealPlan(T).init(allocator, nfft) catch |e| switch (e) {
                error.OutOfMemory => |x| return x,
                error.ZeroLength, error.NotImplemented => unreachable,
            };
            errdefer plan.deinit(allocator);
            return .{ .plan = plan, .window = try allocator.dupe(T, window), .hop = hop, .nfft = nfft };
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            self.plan.deinit(allocator);
            allocator.free(self.window);
            self.* = undefined;
        }

        /// Values per frame's spectrum: `nfft / 2 + 1`.
        pub fn bins(self: Self) usize {
            return self.nfft / 2 + 1;
        }

        /// Frames of a signal of `len` samples: those that fit whole, as
        /// MATLAB's; none when the window is longer.
        pub fn frames(self: Self, len: usize) usize {
            if (len < self.window.len) return 0;
            return (len - self.window.len) / self.hop + 1;
        }

        /// The frame buffer and the plan's.
        pub fn scratchLen(self: Self) usize {
            return self.nfft + self.plan.scratchLen();
        }

        /// `out[k * bins() ..][0..bins()]` is frame `k`'s spectrum, for every
        /// frame of `in` (`frames(in.len)`).
        pub fn forward(self: *const Self, io: Io, out: []Complex(T), in: []const T, scratch: []T) Error!void {
            const nf = self.frames(in.len);
            const m = self.bins();
            if (out.len < nf * m or scratch.len < self.scratchLen()) return error.SliceTooShort;
            const buf = scratch[0..self.nfft];
            const rest = scratch[self.nfft..];
            for (0..nf) |k| {
                const at = in[k * self.hop ..][0..self.window.len];
                for (buf[0..self.window.len], at, self.window) |*b, x, w| b.* = x * w;
                @memset(buf[self.window.len..], 0);
                self.plan.forward(io, out[k * m ..][0..m], buf, rest, .{}, .backward) catch |e| return narrow(e);
            }
        }

        /// The signal from `frames` spectra in `in` (laid out as `forward`'s)
        /// into `out[0 .. (frames - 1) hop + window.len]`: each frame back,
        /// windowed again and overlap-added, every sample divided by the sum
        /// of the squared windows over it (0 where that sum is).
        pub fn inverse(self: *const Self, io: Io, out: []T, in: []const Complex(T), frame_count: usize, scratch: []T) Error!void {
            const m = self.bins();
            const wl = self.window.len;
            if (frame_count == 0) return;
            const len = (frame_count - 1) * self.hop + wl;
            if (out.len < len or in.len < frame_count * m or scratch.len < self.scratchLen()) return error.SliceTooShort;
            const buf = scratch[0..self.nfft];
            const rest = scratch[self.nfft..];
            @memset(out[0..len], 0);
            for (0..frame_count) |k| {
                self.plan.backward(io, buf, in[k * m ..][0..m], rest, .{}, .backward) catch |e| return narrow(e);
                for (out[k * self.hop ..][0..wl], buf[0..wl], self.window) |*o, x, w| o.* += x * w;
            }
            for (out[0..len], 0..) |*o, i| {
                // frames k with k hop <= i < k hop + wl
                const last = @min(i / self.hop, frame_count - 1);
                const first = if (i + 1 > wl) (i + 1 - wl + self.hop - 1) / self.hop else 0;
                var s: T = 0;
                var k = first;
                while (k <= last) : (k += 1) {
                    const w = self.window[i - k * self.hop];
                    s += w * w;
                }
                o.* = if (s > 0) o.* / s else 0;
            }
        }

        fn narrow(e: fft.TransformIntoError) Error {
            return switch (e) {
                error.Canceled => error.Canceled,
                // sized above; the plan is real
                error.SliceTooShort, error.CountMismatch, error.NotImplemented => unreachable,
            };
        }
    };
}

// ---- tests ----

const testing = std.testing;
const window_ = @import("window.zig");

test "Stft: each frame's spectrum is the windowed frame's, and the inverse gives the signal back" {
    const al = testing.allocator;
    const n = 300;
    var x: [n]f64 = undefined;
    for (&x, 0..) |*v, i| v.* = @sin(0.37 * @as(f64, @floatFromInt(i))) + 0.25 * @cos(0.011 * @as(f64, @floatFromInt(i * i)));
    var win: [64]f64 = undefined;
    try window_.window(f64, &win, .hann, .periodic);
    for ([_][2]usize{ .{ 16, 64 }, .{ 32, 128 }, .{ 48, 64 } }) |hn| {
        var s = try Stft(f64).init(al, &win, hn[0], hn[1]);
        defer s.deinit(al);
        const nf = s.frames(n);
        const m = s.bins();
        const spec = try al.alloc(Complex(f64), nf * m);
        defer al.free(spec);
        const scratch = try al.alloc(f64, s.scratchLen());
        defer al.free(scratch);
        try s.forward(testing.io, spec, &x, scratch);
        // frame 2 by itself
        var frame: [128]f64 = @splat(0);
        for (frame[0..64], x[2 * hn[0] ..][0..64], win) |*f, v, w| f.* = v * w;
        var want: [65]Complex(f64) = undefined;
        try s.plan.forward(testing.io, want[0..m], frame[0..hn[1]], scratch[hn[1]..], .{}, .backward);
        for (want[0..m], spec[2 * m ..][0..m]) |a, b| {
            try testing.expectEqual(a.re, b.re);
            try testing.expectEqual(a.im, b.im);
        }
        var y: [n]f64 = undefined;
        try s.inverse(testing.io, &y, spec, nf, scratch);
        // wherever a frame's window is nonzero over the sample
        const len = (nf - 1) * hn[0] + 64;
        for (x[1..len], y[1..len]) |a, b| try testing.expectApproxEqAbs(a, b, 1e-12);
    }
}

test "Stft: init and slice errors" {
    const al = testing.allocator;
    const w = [_]f32{ 1, 1, 1, 1 };
    try testing.expectError(error.ZeroLength, Stft(f32).init(al, &w, 0, 8));
    try testing.expectError(error.WindowTooLong, Stft(f32).init(al, &w, 2, 3));
    var s = try Stft(f32).init(al, &w, 2, 8);
    defer s.deinit(al);
    try testing.expectEqual(@as(usize, 0), s.frames(3));
    try testing.expectEqual(@as(usize, 3), s.frames(8));
    var x: [8]f32 = @splat(1);
    var spec: [10]Complex(f32) = undefined;
    const scratch = try al.alloc(f32, s.scratchLen());
    defer al.free(scratch);
    try testing.expectError(error.SliceTooShort, s.forward(testing.io, &spec, &x, scratch));
    try testing.expectError(error.SliceTooShort, s.forward(testing.io, &spec, &x, scratch[1..]));
}
