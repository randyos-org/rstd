//! `fft.Plan`: the complex 1-D FFT of any length, PocketFFT's `pocketfft_c`
//! (`pocketfft_hdronly.h`, cpp branch, commit 4d865cba, BSD-3, LICENSE):
//! mixed radix (`cfftp/`) unless its cost model prefers Bluestein's
//! (`bluestein.zig`) for a large prime factor. For f32 at pffft's sizes
//! (`n = 16 2^a 3^b 5^c`, at most `2^26`), pffft's SIMD instead
//! (`pffft/`). A batch runs through `io` (`batch.zig`).

const std = @import("std");
const fft = @import("__root__.zig");
const cmplx = @import("cmplx.zig");
const util = @import("util.zig");
const cfftp = @import("cfftp/__root__.zig");
const bluestein = @import("bluestein.zig");
const pffft = @import("pffft/__root__.zig");
const batch_ = @import("batch.zig");

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
        const Kind = union(enum) { pack: cfftp.Cfftp(T), blue: bluestein.Bluestein(T), simd: pffft.Simd(T) };

        pub fn init(allocator: Allocator, n: usize) fft.InitError!Self {
            if (n == 0) return error.ZeroLength;
            if (T == f32 and pffft.validSize(n, .complex)) {
                const s = try pffft.Setup.init(allocator, n, .complex);
                // one gathered transform and pffft's work buffer
                return .{ .n = n, .scratch_len = 2 * n, .kind = .{ .simd = s } };
            }
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
            if (scratch.len < self.scratch_len) return error.SliceTooShort;
            try covers(data.len, self.n, layout);
            // each transform is gathered whole before it is written back
            try self.batch(io, data, data, scratch, direction, layout, layout, options);
        }

        /// `out` and `in` must not overlap.
        pub fn transformInto(self: *const Self, io: Io, out: []C(T), in: []const C(T), scratch: []C(T), direction: Direction, layouts: Layouts, options: Options) fft.TransformIntoError!void {
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            if (scratch.len < self.scratch_len) return error.SliceTooShort;
            try covers(in.len, self.n, layouts.in);
            try covers(out.len, self.n, layouts.out);
            try self.batch(io, out, in, scratch, direction, layouts.in, layouts.out, options);
        }

        /// The checked batch, through `io` (`batch.zig`).
        fn batch(self: *const Self, io: Io, out: []C(T), in: []const C(T), scratch: []C(T), direction: Direction, lin: Layout, lout: Layout, options: Options) Io.Cancelable!void {
            const n = self.n;
            const Ctx = struct {
                plan: *const Self,
                out: []C(T),
                in: []const C(T),
                scratch: []C(T),
                direction: Direction,
                lin: Layout,
                lout: Layout,
                din: usize,
                dout: usize,
                options: Options,

                fn one(c: *const @This(), k: usize, slot: usize) void {
                    const sl = c.plan.scratch_len;
                    const s = c.scratch[slot * sl ..][0..sl];
                    const buf = s[0..c.plan.n];
                    for (buf, 0..) |*b, i| b.* = c.in[k * c.din + i * c.lin.stride];
                    c.plan.run(buf, s[c.plan.n..], c.direction, c.options);
                    for (buf, 0..) |b, i| c.out[k * c.dout + i * c.lout.stride] = b;
                }
            };
            const ctx: Ctx = .{
                .plan = self,
                .out = out,
                .in = in,
                .scratch = scratch,
                .direction = direction,
                .lin = lin,
                .lout = lout,
                .din = lin.distance orelse (n * lin.stride),
                .dout = lout.distance orelse (n * lout.stride),
                .options = options,
            };
            try batch_.run(io, lin.count, scratch.len / self.scratch_len, &ctx, Ctx.one);
        }

        /// One contiguous transform in place, `rest` its working space.
        fn run(self: *const Self, c: []C(T), rest: []C(T), direction: Direction, options: Options) void {
            const fwd = direction == .forward;
            const f = factor(T, self.n, fwd, options.norm);
            switch (self.kind) {
                .pack => |*p| p.exec(c, rest, f, fwd),
                .blue => |*b| b.exec(c, rest, f, fwd),
                .simd => |*s| if (T == f32) {
                    const fl: [*]f32 = @ptrCast(c.ptr);
                    s.transform(fl, fl, @ptrCast(rest.ptr), fwd, true);
                    if (f != 1) for (c) |*x| {
                        x.* = cmplx.scale(T, x.*, f);
                    };
                } else unreachable,
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
