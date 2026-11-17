//! `fft.Plan`: the complex 1-D FFT of any length, PocketFFT's `pocketfft_c`
//! (`pocketfft_hdronly.h`, cpp branch, commit 4d865cba, BSD-3, LICENSE):
//! mixed radix (`cfftp/`) unless its cost model prefers Bluestein's
//! (`bluestein.zig`) for a large prime factor. Batches run one after
//! another (N9 makes them concurrent).

const std = @import("std");
const fft = @import("__root__.zig");
const cmplx = @import("cmplx.zig");
const util = @import("util.zig");
const cfftp = @import("cfftp/__root__.zig");
const bluestein = @import("bluestein.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const C = cmplx.C;
const Direction = fft.Direction;
const Layout = fft.Layout;
const Layouts = fft.Layouts;
const Options = fft.Options;

pub fn Plan(comptime T: type) type {
    if (@typeInfo(T) != .float) @compileError("rstd.fft: " ++ @typeName(T) ++ " is not a float type");
    return struct {
        n: usize,
        scratch_len: usize,
        kind: Kind,

        const Self = @This();
        const Kind = union(enum) { pack: cfftp.Cfftp(T), blue: bluestein.Bluestein(T) };

        pub fn init(allocator: Allocator, n: usize) fft.InitError!Self {
            if (n == 0) return error.ZeroLength;
            // `pocketfft_c`'s choice
            const lpf = if (n < 50) 0 else util.largestPrimeFactor(n);
            if (lpf * lpf > n) blue: {
                const good = util.goodSizeComplex(2 * n - 1) orelse break :blue;
                const comp1 = util.costGuess(n);
                // fudge factor that appears to give good overall performance
                const comp2 = 2 * util.costGuess(good) * 1.5;
                if (comp2 >= comp1) break :blue;
                const b = bluestein.Bluestein(T).init(allocator, n) catch |e| switch (e) {
                    error.OutOfMemory => |x| return x,
                    error.TooLarge => break :blue,
                };
                // one gathered transform, then Bluestein's own
                return .{ .n = n, .scratch_len = n + b.scratchLen(), .kind = .{ .blue = b } };
            }
            const p = try cfftp.Cfftp(T).init(allocator, n);
            // one gathered transform and the ping-pong buffer
            return .{ .n = n, .scratch_len = 2 * n, .kind = .{ .pack = p } };
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            switch (self.kind) {
                inline else => |*k| k.deinit(allocator),
            }
            self.* = undefined;
        }

        /// Scratch per transform, any layout. A batch runs as many transforms
        /// at once as the scratch holds, so the bits never depend on it.
        pub fn scratchLen(self: Self) usize {
            return self.scratch_len;
        }

        pub fn transform(self: *const Self, io: Io, data: []C(T), scratch: []C(T), direction: Direction, layout: Layout, options: Options) fft.TransformError!void {
            _ = io;
            const n = self.n;
            if (scratch.len < self.scratch_len) return error.SliceTooShort;
            try covers(data.len, n, layout);
            const dist = layout.distance orelse (n * layout.stride);
            const buf = scratch[0..n];
            const rest = scratch[n..self.scratch_len];
            for (0..layout.count) |k| {
                const base = k * dist;
                for (buf, 0..) |*b, i| b.* = data[base + i * layout.stride];
                self.run(buf, rest, direction, options);
                for (buf, 0..) |b, i| data[base + i * layout.stride] = b;
            }
        }

        /// `out` and `in` must not overlap.
        pub fn transformInto(self: *const Self, io: Io, out: []C(T), in: []const C(T), scratch: []C(T), direction: Direction, layouts: Layouts, options: Options) fft.TransformIntoError!void {
            _ = io;
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            const n = self.n;
            if (scratch.len < self.scratch_len) return error.SliceTooShort;
            try covers(in.len, n, layouts.in);
            try covers(out.len, n, layouts.out);
            const din = layouts.in.distance orelse (n * layouts.in.stride);
            const dout = layouts.out.distance orelse (n * layouts.out.stride);
            const buf = scratch[0..n];
            const rest = scratch[n..self.scratch_len];
            for (0..layouts.in.count) |k| {
                for (buf, 0..) |*b, i| b.* = in[k * din + i * layouts.in.stride];
                self.run(buf, rest, direction, options);
                for (buf, 0..) |b, i| out[k * dout + i * layouts.out.stride] = b;
            }
        }

        /// One contiguous transform in place, `rest` its working space.
        fn run(self: *const Self, c: []C(T), rest: []C(T), direction: Direction, options: Options) void {
            const fwd = direction == .forward;
            const f = factor(T, self.n, fwd, options.norm);
            switch (self.kind) {
                .pack => |*p| p.exec(c, rest, f, fwd),
                .blue => |*b| b.exec(c, rest, f, fwd),
            }
        }
    };
}

/// The scale by `norm`: `1 / n` on the side it names, `1 / √n` both ways
/// for `ortho`.
pub fn factor(comptime T: type, n: usize, fwd: bool, norm: fft.Norm) T {
    const nf: T = @floatFromInt(n);
    return switch (norm) {
        .backward => if (fwd) 1 else 1 / nf,
        .forward => if (fwd) 1 / nf else 1,
        .ortho => 1 / @sqrt(nf),
    };
}

/// A layout's last element, `(count - 1) distance + (len - 1) stride`, is in
/// a slice of `have`.
pub fn covers(have: usize, len: usize, layout: Layout) error{SliceTooShort}!void {
    if (layout.count == 0) return;
    const dist = layout.distance orelse (std.math.mul(usize, len, layout.stride) catch return error.SliceTooShort);
    const a = std.math.mul(usize, layout.count - 1, dist) catch return error.SliceTooShort;
    const b = std.math.mul(usize, len - 1, layout.stride) catch return error.SliceTooShort;
    const last = std.math.add(usize, a, b) catch return error.SliceTooShort;
    if (last >= have) return error.SliceTooShort;
}
