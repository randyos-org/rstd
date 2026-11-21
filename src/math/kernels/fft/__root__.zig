//! `rstd.math.kernels.fft` (`../../design.md`, #611): PocketFFT's algorithms
//! in Zig, pffft's SIMD (`pffft/`, f32 at its sizes). `Plan` and `RealPlan`
//! are done (`plan.zig`, `real_plan.zig`); the rest are stubs whose
//! bodies return `error.NotImplemented`, their tests skip. A plan depends only
//! on `n`, the element type and the core; immutable, it serves every thread.
//! Transforms allocate nothing.

const std = @import("std");
const complex_plan = @import("plan.zig");
const real_plan = @import("real_plan.zig");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;
const Io = std.Io;

/// `forward` is `exp(-2πi jk/n)`, as NumPy and FFTW.
pub const Direction = enum { forward, backward };

/// Which direction divides by `n`; `ortho` divides both by `√n`.
pub const Normalization = enum { backward, ortho, forward };

/// A batch of transforms in a slice: transform `k`'s element `i` is at
/// `k * distance + i * stride`. `distance` null is `len * stride`, `len`
/// being that slice's elements per transform (`n`, or `n / 2 + 1` for a half
/// spectrum). Transforms must not overlap.
pub const BatchLayout = struct { count: usize = 1, stride: usize = 1, distance: ?usize = null };

/// An out-of-place batch: both sides' counts must match.
pub const BatchLayouts = struct { in: BatchLayout = .{}, out: BatchLayout = .{} };

pub const InitError = error{ OutOfMemory, ZeroLength, NotImplemented };
pub const TransformError = error{ SliceTooShort, Canceled, NotImplemented };
pub const TransformIntoError = TransformError || error{CountMismatch};

/// One transform at a time, no `io`: audio frames.
pub const FrameError = error{ SliceTooShort, NotImplemented };

/// `RealPlan`'s unordered forms, pffft's: f32 only (`UnsupportedType`;
/// f64 SIMD is later), at pffft's real sizes, `n = 32 2^a 3^b 5^c` up to
/// `2^26` (`UnsupportedSize`).
pub const UnorderedError = error{ SliceTooShort, UnsupportedSize, UnsupportedType };

/// Complex, any `n`: PocketFFT's mixed radix (2, 3, 4, 5, 7, 8, 11 and a
/// generic pass), Bluestein's for large prime factors; f32 at pffft's
/// sizes, `16 2^a 3^b 5^c`, pffft's SIMD (`plan.zig`).
pub const Plan = complex_plan.Plan;

/// Real input: `n` reals to the `n / 2 + 1` values of the half spectrum
/// (NumPy's `rfft`), and back: PocketFFT's `rfftp`, or Bluestein's; f32 at
/// pffft's sizes, `32 2^a 3^b 5^c`, pffft's SIMD, with its unordered forms
/// for fast convolution (`real_plan.zig`).
pub const RealPlan = real_plan.RealPlan;

/// `Plan` for a comptime `n`: factors and twiddles at comptime, unrolled, no
/// allocator; the same calls on a value `.{}`.
pub fn StaticPlan(comptime n: usize, comptime T: type) type {
    comptime checkFloat(T);
    if (n == 0) @compileError("rstd.math.kernels.fft.StaticPlan: n is 0");
    return struct {
        const Self = @This();
        const C = Complex(T);

        /// Stub: PocketFFT's ping-pong buffer; Bluestein's sizes need more.
        pub fn scratchLen(_: Self) usize {
            return n;
        }

        pub fn transform(self: Self, io: Io, data: []C, scratch: []C, direction: Direction, layout: BatchLayout, normalization: Normalization) TransformError!void {
            _ = self;
            _ = io;
            _ = data;
            _ = scratch;
            _ = direction;
            _ = layout;
            _ = normalization;
            return error.NotImplemented;
        }

        pub fn transformInto(self: Self, io: Io, out: []C, in: []const C, scratch: []C, direction: Direction, layouts: BatchLayouts, normalization: Normalization) TransformIntoError!void {
            _ = self;
            _ = io;
            _ = out;
            _ = in;
            _ = scratch;
            _ = direction;
            _ = normalization;
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            return error.NotImplemented;
        }
    };
}

/// `RealPlan` for a comptime `n`, as `StaticPlan`.
pub fn StaticRealPlan(comptime n: usize, comptime T: type) type {
    comptime checkFloat(T);
    if (n == 0) @compileError("rstd.math.kernels.fft.StaticRealPlan: n is 0");
    return struct {
        const Self = @This();
        const C = Complex(T);

        /// Stub, as `StaticPlan`'s.
        pub fn scratchLen(_: Self) usize {
            return n;
        }

        pub fn forward(self: Self, io: Io, out: []C, in: []const T, scratch: []T, layouts: BatchLayouts, normalization: Normalization) TransformIntoError!void {
            _ = self;
            _ = io;
            _ = out;
            _ = in;
            _ = scratch;
            _ = normalization;
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            return error.NotImplemented;
        }

        pub fn backward(self: Self, io: Io, out: []T, in: []const C, scratch: []T, layouts: BatchLayouts, normalization: Normalization) TransformIntoError!void {
            _ = self;
            _ = io;
            _ = out;
            _ = in;
            _ = scratch;
            _ = normalization;
            if (layouts.in.count != layouts.out.count) return error.CountMismatch;
            return error.NotImplemented;
        }
    };
}

/// DCT and DST types 1 to 4, as SciPy's `dct` and `dst`: `.backward` leaves
/// the forward unscaled, `.ortho` makes the matrix orthogonal. A type 2's
/// inverse is the type 3, and the reverse.
pub const TrigTransformKind = enum { dct1, dct2, dct3, dct4, dst1, dst2, dst3, dst4 };

pub fn TrigTransformPlan(comptime T: type) type {
    comptime checkFloat(T);
    return struct {
        n: usize,
        kind: TrigTransformKind,
        scratch_len: usize,

        const Self = @This();

        pub fn init(allocator: Allocator, n: usize, kind: TrigTransformKind) InitError!Self {
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

        pub fn transform(self: *const Self, io: Io, data: []T, scratch: []T, layout: BatchLayout, normalization: Normalization) TransformError!void {
            _ = self;
            _ = io;
            _ = data;
            _ = scratch;
            _ = layout;
            _ = normalization;
            return error.NotImplemented;
        }
    };
}

/// The MDCT (audio codecs) on a DCT-IV: `2 n` samples to `n` coefficients
/// and back. The caller windows; with a Princen-Bradley window (sine, KBD),
/// windowed overlap-add at hop `n` gives back the input (TDAC).
pub fn ModifiedCosineTransform(comptime T: type) type {
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
    if (@typeInfo(T) != .float) @compileError("rstd.math.kernels.fft: " ++ @typeName(T) ++ " is not a float type");
}

test {
    _ = @import("util.zig");
    _ = @import("pffft/__root__.zig");
    _ = @import("_test/pffft.zig");
    _ = @import("_test/plan.zig");
}
