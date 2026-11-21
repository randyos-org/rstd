//! `rstd.math.fft` (`design.md`, #611): transforms of whole vectors, as
//! NumPy's `fft`, `ifft`, `rfft` and `irfft`. Each plans per call over
//! `kernels.fft` and allocates its result; the `*With` siblings take
//! `FftOptions`.

const std = @import("std");
const kernels_fft = @import("kernels/fft/__root__.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Complex = std.math.Complex;
const Direction = kernels_fft.Direction;

/// The `*With` transforms' options; each default is the plain function's.
pub const FftOptions = struct {
    /// NumPy's default: the inverse divides by `n`, so a round trip gives
    /// back its input.
    normalization: kernels_fft.Normalization = .backward,
};

pub const FftError = error{ OutOfMemory, ZeroLength, Canceled };
pub const IrfftError = FftError || error{LengthMismatch};

/// The DFT of `x`, `exp(-2πi jk/n)`.
pub fn fft(allocator: Allocator, io: Io, comptime T: type, x: []const Complex(T)) FftError![]Complex(T) {
    return complexTransform(allocator, io, T, x, .forward, .{});
}

/// The inverse DFT of `x`, divided by `n`.
pub fn ifft(allocator: Allocator, io: Io, comptime T: type, x: []const Complex(T)) FftError![]Complex(T) {
    return complexTransform(allocator, io, T, x, .backward, .{});
}

/// `fft`, scaled as `options` asks.
pub fn fftWith(allocator: Allocator, io: Io, comptime T: type, x: []const Complex(T), options: FftOptions) FftError![]Complex(T) {
    return complexTransform(allocator, io, T, x, .forward, options);
}

/// `ifft`, scaled as `options` asks.
pub fn ifftWith(allocator: Allocator, io: Io, comptime T: type, x: []const Complex(T), options: FftOptions) FftError![]Complex(T) {
    return complexTransform(allocator, io, T, x, .backward, options);
}

/// The `n / 2 + 1` values of the half spectrum of the `n` reals in `x`.
pub fn rfft(allocator: Allocator, io: Io, comptime T: type, x: []const T) FftError![]Complex(T) {
    return rfftWith(allocator, io, T, x, .{});
}

/// `n` reals from the half spectrum `x` (`n / 2 + 1` values, else
/// `error.LengthMismatch`), divided by `n`. The imaginary parts at 0, and at
/// `n / 2` for even `n`, are ignored, as NumPy's.
pub fn irfft(allocator: Allocator, io: Io, comptime T: type, n: usize, x: []const Complex(T)) IrfftError![]T {
    return irfftWith(allocator, io, T, n, x, .{});
}

/// `rfft`, scaled as `options` asks.
pub fn rfftWith(allocator: Allocator, io: Io, comptime T: type, x: []const T, options: FftOptions) FftError![]Complex(T) {
    var plan = kernels_fft.RealPlan(T).init(allocator, x.len) catch |e| return initError(e);
    defer plan.deinit(allocator);
    const out = try allocator.alloc(Complex(T), x.len / 2 + 1);
    errdefer allocator.free(out);
    const scratch = try allocator.alloc(T, plan.scratchLen());
    defer allocator.free(scratch);
    plan.forward(io, out, x, scratch, .{}, options.normalization) catch |e| return transformError(e);
    return out;
}

/// `irfft`, scaled as `options` asks.
pub fn irfftWith(allocator: Allocator, io: Io, comptime T: type, n: usize, x: []const Complex(T), options: FftOptions) IrfftError![]T {
    var plan = kernels_fft.RealPlan(T).init(allocator, n) catch |e| return initError(e);
    defer plan.deinit(allocator);
    if (x.len != n / 2 + 1) return error.LengthMismatch;
    const out = try allocator.alloc(T, n);
    errdefer allocator.free(out);
    const scratch = try allocator.alloc(T, plan.scratchLen());
    defer allocator.free(scratch);
    plan.backward(io, out, x, scratch, .{}, options.normalization) catch |e| return transformError(e);
    return out;
}

fn complexTransform(allocator: Allocator, io: Io, comptime T: type, x: []const Complex(T), direction: Direction, options: FftOptions) FftError![]Complex(T) {
    var plan = kernels_fft.Plan(T).init(allocator, x.len) catch |e| return initError(e);
    defer plan.deinit(allocator);
    const out = try allocator.alloc(Complex(T), x.len);
    errdefer allocator.free(out);
    const scratch = try allocator.alloc(Complex(T), plan.scratchLen());
    defer allocator.free(scratch);
    plan.transformInto(io, out, x, scratch, direction, .{}, options.normalization) catch |e| return transformError(e);
    return out;
}

/// `Plan` and `RealPlan` have no stub path, so never `NotImplemented`.
fn initError(e: kernels_fft.InitError) FftError {
    return switch (e) {
        error.OutOfMemory => error.OutOfMemory,
        error.ZeroLength => error.ZeroLength,
        error.NotImplemented => unreachable,
    };
}

/// Buffers sized here and one transform per call: only `Canceled` is left.
fn transformError(e: kernels_fft.TransformIntoError) FftError {
    return switch (e) {
        error.Canceled => error.Canceled,
        error.SliceTooShort, error.CountMismatch, error.NotImplemented => unreachable,
    };
}

test {
    _ = @import("_test/fft.zig");
}
