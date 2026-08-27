const std = @import("std");
const types = @import("types.zig");
const V2f32 = types.V2f32;
const V3f32 = types.V3f32;

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

pub inline fn inversesqrt(x: anytype) @TypeOf(x) {
    return 1 / @sqrt(x);
}

pub fn mix(a: anytype, b: anytype, t: f32) @TypeOf(a) {
    const tvec: @TypeOf(a) = @splat(t);
    return a + (b - a) * tvec;
}

pub fn step_scalar(edge: f32, x: f32) f32 {
    return if (x < edge) 0 else 1;
}

pub inline fn fract(x: anytype) @TypeOf(x) {
    return x - @floor(x);
}

pub inline fn dot(a: anytype, b: anytype) f32 {
    return sum(a * b);
}

pub inline fn sum(v: anytype) @TypeOf(v[0]) {
    return @reduce(.Add, v);
}

pub inline fn dot2(v: anytype) f32 {
    return sum(v * v);
}

pub inline fn length(v: anytype) f32 {
    return @sqrt(dot2(v));
}

pub fn smoothstep_scalar(edge0: f32, edge1: f32, x: f32) f32 {
    const t = std.math.clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

pub inline fn distance(p0: anytype, p1: anytype) f32 {
    return length(p1 - p0);
}

pub inline fn cross(a: V3f32, b: V3f32) V3f32 {
    return .{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

pub inline fn normalize(v: V3f32) V3f32 {
    return v / length(v);
}

pub inline fn reflect(i: anytype, n: anytype) @TypeOf(i) {
    return i - 2 * dot(i, n) * n;
}

pub fn refract(i: anytype, n: anytype, eta: f32) @TypeOf(i) {
    const dotni = dot(n, i);
    const k = 1.0 - eta * eta * (1.0 - dotni * dotni);
    return if (k < 0.0) 0 * i else eta * i - (eta * dotni + @sqrt(k)) * n;
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
