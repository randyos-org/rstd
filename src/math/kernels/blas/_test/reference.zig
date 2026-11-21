//! The wide reference the tests compute results in, for any float type: `T`
//! at `f128`'s precision, its arithmetic, and the checks' error bounds in ulps.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../types.zig");
const check = @import("check.zig");

const isComplex = types.isComplex;
const FloatOf = types.Float;

// ---- error bounds of the wide checks, in ulps of the result --------------------

/// An element-wise result (`scal`, `axpy`, `axpby`, `rot`, `rotm`).
pub const elem_ulps: f128 = 4;
/// An element of a rank update: `A + alpha x y`, and the second product of a
/// rank-2 one.
pub const update_ulps: f128 = 8;

/// A dot of `n` terms.
pub fn dotUlps(n: usize) f128 {
    return @floatFromInt(4 * (n + 2));
}

/// A sum of `n` terms (`asum`, `sum`, `nrm2`'s squares, a level-3 element).
pub fn sumUlps(n: usize) f128 {
    return @floatFromInt(2 * (n + 2));
}

/// An element of a triangular solve of `n` unknowns, which carries the
/// errors of the ones before it.
pub fn solveUlps(n: usize) f128 {
    return @floatFromInt(16 * (n + 2));
}

/// How far, relative, a found extreme may be from its key: 2 eps.
pub fn keyTol(comptime F: type) f128 {
    return 2 * std.math.floatEps(F);
}

/// The least divisor magnitude a scale of the wide checks takes, so a zero
/// pivot gives a finite one.
pub const tiny_mag: f128 = 1e-300;

// ---- a wide reference, for any float type ----------------------------------

/// `T` at the precision of the `f128` reference.
pub fn Wide(comptime T: type) type {
    return if (isComplex(T)) std.math.Complex(f128) else f128;
}

pub fn wide(comptime T: type, v: T) Wide(T) {
    return if (comptime isComplex(T)) .init(v.re, v.im) else v;
}

/// `@sqrt` of a `Wide` value. sparc64's compiler-rt `_Qp_sqrt` takes the root
/// in `f64` (issue #184); two Newton steps there restore `f128`'s precision.
pub fn sqrtW(x: f128) f128 {
    var y = @sqrt(x);
    if (comptime sysinfo.target.cpu.arch != .sparc64) return y;
    if (!(x > 0) or !std.math.isFinite(x)) return y;
    // the guess is taken in f64: keep `x` in its range
    if (x < 0x1p-1000) return sqrtW(x * 0x1p1000) * 0x1p-500;
    if (x > 0x1p1000) return sqrtW(x * 0x1p-1000) * 0x1p500;
    inline for (0..2) |_| y = (y + x / y) / 2;
    return y;
}

/// Ulps of `F` the library's own square roots may add: on sparc64 its `f128`
/// roots keep `f64`'s precision (issue #184).
pub fn rootUlps(comptime F: type) f128 {
    return if (sysinfo.target.cpu.arch == .sparc64 and F == f128) 0x1p60 else 0;
}

/// `|re| + |im|` of a `Wide` value, `|v|` for a real one.
pub fn abs1(v: anytype) f128 {
    return if (comptime isComplex(@TypeOf(v))) @abs(v.re) + @abs(v.im) else @abs(v);
}

/// `abs1` of the finite parts: a tolerance's scale, where the results a
/// special part reaches are NaN or infinite and need none.
pub fn mag(v: anytype) f128 {
    if (comptime isComplex(@TypeOf(v))) return mag(v.re) + mag(v.im);
    return if (std.math.isFinite(v)) @abs(v) else 0;
}

pub fn mulW(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.mul(b) else a * b;
}

pub fn addW(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.add(b) else a + b;
}

pub fn subW(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.sub(b) else a - b;
}

pub fn divW(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.div(b) else a / b;
}

/// `f` times `a`, part by part: a real factor, with no cross terms.
pub fn scaleW(f: f128, a: anytype) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) .init(f * a.re, f * a.im) else f * a;
}

pub fn conjW(a: anytype) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.conjugate() else a;
}

/// `got` against the reference `want`, within `ulps` units of `T`'s
/// precision of `scale`: a backward-error bound, `scale` being the sum of the
/// magnitudes that went into `want`.
pub fn expectNear(comptime T: type, want: Wide(T), got: T, ulps: f128, scale: f128) !void {
    const F = FloatOf(T);
    const tol = ulps * std.math.floatEps(F) * scale + std.math.floatMin(F);
    try check.expectClose(Wide(T), want, wide(T, got), tol);
}
