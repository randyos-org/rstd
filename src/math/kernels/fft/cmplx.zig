//! PocketFFT's `cmplx` operations (`pocketfft_hdronly.h`, cpp branch,
//! commit 4d865cba, BSD-3, LICENSE), in its operation order so the bits
//! follow from the source.

const std = @import("std");

pub fn C(comptime T: type) type {
    return std.math.Complex(T);
}

pub inline fn add(comptime T: type, a: C(T), b: C(T)) C(T) {
    return .{ .re = a.re + b.re, .im = a.im + b.im };
}

pub inline fn sub(comptime T: type, a: C(T), b: C(T)) C(T) {
    return .{ .re = a.re - b.re, .im = a.im - b.im };
}

pub inline fn scale(comptime T: type, a: C(T), s: T) C(T) {
    return .{ .re = a.re * s, .im = a.im * s };
}

/// `special_mul`: `a conj(w)` forward, `a w` backward.
pub inline fn twiddle(comptime T: type, comptime fwd: bool, a: C(T), w: C(T)) C(T) {
    return if (fwd)
        .{ .re = a.re * w.re + a.im * w.im, .im = a.im * w.re - a.re * w.im }
    else
        .{ .re = a.re * w.re - a.im * w.im, .im = a.re * w.im + a.im * w.re };
}

/// `ROTX90`: `-i a` forward, `i a` backward.
pub inline fn rot90(comptime T: type, comptime fwd: bool, a: C(T)) C(T) {
    return if (fwd) .{ .re = a.im, .im = -a.re } else .{ .re = -a.im, .im = a.re };
}

const hsqt2 = 0.707106781186547524400844362104849;

/// `ROTX45`.
pub inline fn rot45(comptime T: type, comptime fwd: bool, a: C(T)) C(T) {
    const h: T = hsqt2;
    return if (fwd) .{ .re = h * (a.re + a.im), .im = h * (a.im - a.re) } else .{ .re = h * (a.re - a.im), .im = h * (a.im + a.re) };
}

/// `ROTX135`.
pub inline fn rot135(comptime T: type, comptime fwd: bool, a: C(T)) C(T) {
    const h: T = hsqt2;
    return if (fwd) .{ .re = h * (a.im - a.re), .im = h * (-a.re - a.im) } else .{ .re = h * (-a.re - a.im), .im = h * (a.re - a.im) };
}
