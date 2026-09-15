//! Generic vector-math routines shared by every consumer that works with
//! `@Vector` types, regardless of component count or scalar type - CPU-side
//! geometry math builds its domain-specific API on top of these instead of
//! keeping its own copy of the same handful of formulas.

const std = @import("std");

pub const matlab = @import("matlab/__root__.zig");

pub const V2f32 = @Vector(2, f32);
pub const V3f32 = @Vector(3, f32);
pub const V4f32 = @Vector(4, f32);

pub const V2f64 = @Vector(2, f64);
pub const V3f64 = @Vector(3, f64);
pub const V4f64 = @Vector(4, f64);

pub inline fn splatVec(N: usize, v: anytype) @Vector(N, @TypeOf(v)) {
    return splatVecType(@TypeOf(v), N, v);
}

pub inline fn splatVecF32(N: usize, v: f32) @Vector(N, f32) {
    return splatVecType(f32, N, v);
}

pub inline fn splatV2f32(v: f32) V2f32 {
    return splatVecF32(2, v);
}

pub inline fn splatV3f32(v: f32) V3f32 {
    return splatVecF32(3, v);
}

pub inline fn splatV4f32(v: f32) V4f32 {
    return splatVecF32(4, v);
}

pub inline fn splatV2f64(v: f64) V2f64 {
    return splatVecF64(2, v);
}

pub inline fn splatV3f64(v: f64) V3f64 {
    return splatVecF64(3, v);
}

pub inline fn splatV4f64(v: f64) V4f64 {
    return splatVecF64(4, v);
}

pub inline fn splatVecF64(N: usize, v: f64) @Vector(N, f64) {
    return splatVecType(f64, N, v);
}

pub inline fn splatVecType(comptime T: type, N: usize, v: T) @Vector(N, T) {
    return @as(@Vector(N, T), @splat(v));
}

pub inline fn splatAs(typeSrcVec: anytype, v: anytype) @TypeOf(typeSrcVec) {
    return @as(@TypeOf(typeSrcVec), @splat(v));
}

pub inline fn sum(v: anytype) @TypeOf(v[0]) {
    return @reduce(.Add, v);
}

pub inline fn zeroVecType(N: usize, T: type) @Vector(N, T) {
    return @as(@Vector(N, T), @splat(@as(T, 0)));
}

pub inline fn zeroVecF32(N: usize) @Vector(N, f32) {
    return zeroVecType(N, f32);
}

pub inline fn zeroVecF64(N: usize) @Vector(N, f64) {
    return zeroVecType(N, f64);
}

/// Fused multiply-accumulate chain rather than an elementwise multiply
/// followed by a horizontal add: one fewer independently-rounded step,
/// which matters for the curve/distance math this feeds (see msdf/math.zig).
pub inline fn dot(a: anytype, b: @TypeOf(a)) @TypeOf(a[0]) {
    const T = @TypeOf(a[0]);
    const len = @typeInfo(@TypeOf(a)).vector.len;
    var acc: T = a[0] * b[0];
    inline for (1..len) |i| acc = @mulAdd(T, a[i], b[i], acc);
    return acc;
}

pub inline fn dot2(v: anytype) @TypeOf(v[0]) {
    return dot(v, v);
}

pub inline fn length(v: anytype) @TypeOf(v[0]) {
    return @sqrt(dot2(v));
}

pub inline fn distance(a: anytype, b: @TypeOf(a)) @TypeOf(a[0]) {
    return length(b - a);
}

/// Normalizes `v`, returning `zero_fallback` instead of NaN/inf if `v` is
/// the zero vector.
pub inline fn norm(v: anytype, zero_fallback: @TypeOf(v)) @TypeOf(v) {
    const len = length(v);
    if (len == 0) return zero_fallback;
    return v / splatAs(v, len);
}

/// 2D perp-dot product ("2D cross product"): the scalar a 3D cross product's
/// z-component would be if both inputs were embedded in the xy-plane.
/// Positive when `b` is counter-clockwise from `a`.
pub inline fn cross2(a: anytype, b: @TypeOf(a)) @TypeOf(a[0]) {
    return a[0] * b[1] - a[1] * b[0];
}

/// Standard 3D cross product.
pub inline fn cross3(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return .{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

/// Linear interpolation, generic over `a`/`b` being a scalar or a vector,
/// and `t` being a scalar or a per-component vector weight.
pub fn mix(a: anytype, b: anytype, t: anytype) @TypeOf(a, b) {
    const BaseType = @TypeOf(a, b);
    const WeightType = @TypeOf(t);
    const weight_info = @typeInfo(WeightType);
    switch (@typeInfo(BaseType)) {
        .float, .comptime_float => {
            if (weight_info != .float and weight_info != .comptime_float)
                @compileError("Invalid weight type, float base types require float weight types");
            return @mulAdd(@TypeOf(a, b), 1 - t, a, t * b);
        },
        .vector => {
            const weight: BaseType = if (weight_info == .vector) t else @splat(t);
            const one: BaseType = @splat(1);
            return @mulAdd(@TypeOf(a, b), one - weight, a, weight * b);
        },
        else => @compileError("Invalid base type, only floats and their vectors are supported"),
    }
}

test "dot sanity" {
    const ret = dot(V3f32{ 1, 2, 3 }, V3f32{ 4, 5, 6 });
    try std.testing.expectApproxEqAbs(32, ret, 0.001);
}

test "dot/length/distance agree with the obvious 3-4-5 triangle" {
    const a: V2f32 = .{ 0, 0 };
    const b: V2f32 = .{ 3, 4 };
    try std.testing.expectEqual(@as(f32, 25), dot2(b - a));
    try std.testing.expectEqual(@as(f32, 5), length(b - a));
    try std.testing.expectEqual(@as(f32, 5), distance(a, b));
}

test "cross2 is positive for a counter-clockwise turn" {
    const a: V2f64 = .{ 1, 0 };
    const b: V2f64 = .{ 0, 1 };
    try std.testing.expectEqual(@as(f64, 1), cross2(a, b));
}

test "cross3 matches the standard basis identity x cross y = z" {
    const x: V3f32 = .{ 1, 0, 0 };
    const y: V3f32 = .{ 0, 1, 0 };
    try std.testing.expectEqual(V3f32{ 0, 0, 1 }, cross3(x, y));
}

test "cross sanity" {
    const ret = cross3(V3f32{ 1, 2, 3 }, V3f32{ 5, 11, 34 });
    try std.testing.expectApproxEqAbs(35, ret[0], 0.0001);
    try std.testing.expectApproxEqAbs(-19, ret[1], 0.0001);
    try std.testing.expectApproxEqAbs(1, ret[2], 0.0001);
}

test "normalizeSafe falls back on the zero vector instead of producing NaN" {
    const zero: V2f64 = .{ 0, 0 };
    const fallback: V2f64 = .{ 0, 1 };
    try std.testing.expectEqual(fallback, norm(zero, fallback));

    const v: V2f64 = .{ 0, 5 };
    try std.testing.expectEqual(V2f64{ 0, 1 }, norm(v, fallback));
}
