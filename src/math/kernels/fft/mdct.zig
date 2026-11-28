//! `fft.ModifiedCosineTransform`: the MDCT, `X(k) = Σ x(j) cos(π/n (j + 1/2
//! + n/2)(k + 1/2))` over `2 n` samples, and its inverse with `2 / n`, by the
//! DCT-IV of the samples folded to `n` (`(-c' - d, a - b')` of the quarters
//! `a b c d`, `'` reversing; Princen and Bradley's TDAC, as Wikipedia's
//! "Modified discrete cosine transform" derives it) on `TrigTransformPlan`.
//! An even `n` only, for whole quarters.

const std = @import("std");
const fft = @import("__root__.zig");
const trig_plan = @import("trig_plan.zig");

const Allocator = std.mem.Allocator;

/// `fft.InitError`, and `OddLength` for an odd `n`.
pub const InitError = fft.InitError || error{OddLength};

pub fn ModifiedCosineTransform(comptime T: type) type {
    if (@typeInfo(T) != .float) @compileError("rstd.math.kernels.fft: " ++ @typeName(T) ++ " is not a float type");
    return struct {
        n: usize,
        scratch_len: usize,
        dct4: trig_plan.TrigTransformPlan(T),

        const Self = @This();

        pub fn init(allocator: Allocator, n: usize) InitError!Self {
            if (n == 0) return error.ZeroLength;
            if (n % 2 == 1) return error.OddLength;
            const p = trig_plan.TrigTransformPlan(T).init(allocator, n, .dct4) catch |e| switch (e) {
                error.LengthTooShort => unreachable,
                else => |x| return x,
            };
            return .{ .n = n, .scratch_len = p.scratchLen(), .dct4 = p };
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            self.dct4.deinit(allocator);
            self.* = undefined;
        }

        pub fn scratchLen(self: Self) usize {
            return self.scratch_len;
        }

        /// `out[0..n]` from `in[0 .. 2 n]`; the slices must not overlap.
        pub fn forward(self: *const Self, out: []T, in: []const T, scratch: []T) fft.FrameError!void {
            const n = self.n;
            const h = n / 2;
            if (out.len < n or in.len < 2 * n or scratch.len < self.scratch_len) return error.SliceTooShort;
            for (0..h) |j| {
                out[j] = -in[n + h - 1 - j] - in[n + h + j];
                out[h + j] = in[j] - in[n - 1 - j];
            }
            // the DCT-IV's `2 Σ`, halved
            self.dct4.exec(out[0..n], scratch[n..], 0.5, false);
        }

        /// `out[0 .. 2 n]` from `in[0..n]`, scaled by `2 / n` and still aliased:
        /// windowed overlap-add at hop `n` with a Princen-Bradley window
        /// (sine, KBD) removes the aliasing and gives back the input. The
        /// slices must not overlap.
        pub fn backward(self: *const Self, out: []T, in: []const T, scratch: []T) fft.FrameError!void {
            const n = self.n;
            const h = n / 2;
            if (out.len < 2 * n or in.len < n or scratch.len < self.scratch_len) return error.SliceTooShort;
            const v = scratch[0..n];
            @memcpy(v, in[0..n]);
            const nf: T = @floatFromInt(n);
            self.dct4.exec(v, scratch[n..], 1 / nf, false);
            for (0..h) |j| out[j] = v[h + j];
            for (h..3 * h) |j| out[j] = -v[3 * h - 1 - j];
            for (3 * h..2 * n) |j| out[j] = -v[j - 3 * h];
        }
    };
}
