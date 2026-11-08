const std = @import("std");
const linalg = @import("../../../linalg/__root__.zig");
const vec = linalg.vector;
const V3f32 = vec.V3f32;

pub const radians = std.math.degreesToRadians;
pub const degrees = std.math.radiansToDegrees;
pub const sin = std.math.sin;
pub const cos = std.math.cos;
pub const tan = std.math.tan;
pub const asin = std.math.asin;
pub const acos = std.math.acos;
pub const atan = std.math.atan;
pub const pow = std.math.pow;
pub const exp = std.math.exp;
pub const log = std.math.log;
pub const exp2 = std.math.exp2;
pub const log2 = std.math.log2;
pub const sqrt = std.math.sqrt;
pub const sign = std.math.sign;
pub const floor = std.math.floor;
pub const ceil = std.math.ceil;
pub const mod = std.math.mod;
pub const clamp = std.math.clamp;

/// Scalar `s` as `T`: splatted when `T` is a vector.
inline fn like(comptime T: type, s: anytype) T {
    return if (@typeInfo(T) == .vector) @splat(s) else s;
}

pub inline fn inversesqrt(x: anytype) @TypeOf(x) {
    return like(@TypeOf(x), 1) / @sqrt(x);
}

pub fn stepScalar(edge: f32, x: f32) f32 {
    return if (x < edge) 0 else 1;
}

pub inline fn fract(x: anytype) @TypeOf(x) {
    return x - @floor(x);
}

/// Median of three scalars; the standard way to recover a signed distance
/// from a 3-channel MSDF texture sample
pub inline fn median3(a: anytype, b: @TypeOf(a), c: @TypeOf(a)) @TypeOf(a) {
    return @max(@min(a, b), @min(@max(a, b), c));
}

pub const mix = vec.mix;
pub const dot = vec.dot;
pub const dot2 = vec.dot2;
pub const length = vec.length;
pub const distance = vec.distance;
pub const cross = vec.cross;
pub const sum = vec.sum;

pub inline fn normalize(v: V3f32) V3f32 {
    return vec.norm(v) catch @splat(0);
}

pub fn smoothstepScalar(edge0: f32, edge1: f32, x: f32) f32 {
    const t = std.math.clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

pub inline fn reflect(i: anytype, n: anytype) @TypeOf(i) {
    const T = @TypeOf(i);
    return i - like(T, 2 * dot(i, n)) * n;
}

pub fn refract(i: anytype, n: anytype, eta: f32) @TypeOf(i) {
    const T = @TypeOf(i);
    const dotni = dot(n, i);
    const k = 1.0 - eta * eta * (1.0 - dotni * dotni);
    return if (k < 0.0) like(T, 0) else like(T, eta) * i - like(T, eta * dotni + @sqrt(k)) * n;
}

pub inline fn faceforward(n: anytype, i: anytype, nref: anytype) @TypeOf(n) {
    return if (dot(nref, i) < 0.0) n else -n;
}

// matrixCompMult()
// lessThan()
// lessThanEqual()
// greaterThan()
// greaterThanEqual()
// equal()
// notEqual()
// any()
// all()
// not()
