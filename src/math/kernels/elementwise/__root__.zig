//! `rstd.math.kernels.elementwise` (`../../design.md`, #614): functions over
//! arrays, MKL's VM and Accelerate's vForce. One input: `y[i] = f(x[i])`,
//! `y` may be `x`; two: `z[i] = f(x[i], y[i])`, `z` may be either. The
//! transcendental functions are Zig's builtins, with the accuracy of Zig's
//! `std.math` (musl's algorithms) on every lane; the arithmetic, `sqrt`,
//! `abs` and the roundings are exact (correctly rounded). Real floats, and
//! complex for the arithmetic and `abs`; the rest is #626.

const std = @import("std");
const types = @import("../blas/types.zig");
const strided = @import("../blas/strided.zig");
const map = @import("map.zig");

const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const ConstStrided = strided.ConstStrided;
const Strided = strided.Strided;

/// A length mismatch between the views, or cancellation.
pub const Error = map.Error;

/// `x + y`, correctly rounded.
pub fn add(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void {
    return map.binary(io, T, .add, x, y, z);
}

/// `x − y`, correctly rounded.
pub fn sub(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void {
    return map.binary(io, T, .sub, x, y, z);
}

/// `x · y`, correctly rounded for reals; complex as `std.math.Complex.mul`.
pub fn mul(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void {
    return map.binary(io, T, .mul, x, y, z);
}

/// `x / y`, correctly rounded for reals; complex as `std.math.Complex.div`.
pub fn div(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void {
    return map.binary(io, T, .div, x, y, z);
}

/// `√x`, correctly rounded.
pub fn sqrt(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .sqrt, x, y);
}

/// `eˣ`, as `@exp`.
pub fn exp(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .exp, x, y);
}

/// `2ˣ`, as `@exp2`.
pub fn exp2(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .exp2, x, y);
}

/// `ln x`, as `@log`.
pub fn log(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .log, x, y);
}

/// `log₂ x`, as `@log2`.
pub fn log2(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .log2, x, y);
}

/// `log₁₀ x`, as `@log10`.
pub fn log10(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .log10, x, y);
}

/// `sin x`, as `@sin`.
pub fn sin(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .sin, x, y);
}

/// `cos x`, as `@cos`.
pub fn cos(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .cos, x, y);
}

/// `tan x`, as `@tan`.
pub fn tan(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .tan, x, y);
}

/// `|x|`, exact; complex: the modulus by `std.math.hypot`, without overflow.
pub fn abs(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(Float(T))) Error!void {
    return map.unary(io, T, .abs, x, y);
}

/// `⌊x⌋`, exact.
pub fn floor(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .floor, x, y);
}

/// `⌈x⌉`, exact.
pub fn ceil(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .ceil, x, y);
}

/// `x` rounded toward zero, exact.
pub fn trunc(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .trunc, x, y);
}

/// `x` rounded to nearest, ties away from zero (Fortran's `ANINT`), exact.
pub fn round(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void {
    return real(io, T, .round, x, y);
}

fn real(io: Io, comptime T: type, comptime op: map.Unary, x: ConstStrided(T), y: Strided(T)) Error!void {
    if (comptime isComplex(T)) @compileError("rstd.math.kernels.elementwise." ++ @tagName(op) ++ ": complex is #626");
    return map.unary(io, T, op, x, y);
}

test {
    _ = @import("_test/elementwise.zig");
}
