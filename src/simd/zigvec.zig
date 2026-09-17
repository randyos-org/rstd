//! Generic vector-math routines shared by every consumer that works with
//! `@Vector` types, regardless of component count or scalar type - CPU-side
//! geometry math builds its domain-specific API on top of these instead of
//! keeping its own copy of the same handful of formulas.

const std = @import("std");
const simd = std.simd;

pub const f32_simd_lanes = std.simd.suggestVectorLength(f32) orelse 1;
pub const SimdVecF32 = @Vector(f32_simd_lanes, f32);

pub const f64_simd_lanes = std.simd.suggestVectorLength(f64) orelse 1;
pub const SimdVecF64 = @Vector(f64_simd_lanes, f64);

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

pub inline fn zeroVecType(N: usize, T: type) @Vector(N, T) {
    return @as(@Vector(N, T), @splat(@as(T, 0)));
}

pub inline fn zeroVecF32(N: usize) @Vector(N, f32) {
    return zeroVecType(N, f32);
}

pub inline fn zeroVecF64(N: usize) @Vector(N, f64) {
    return zeroVecType(N, f64);
}

pub inline fn zeroAs(typeSrcVec: anytype) @TypeOf(typeSrcVec) {
    return splatAs(typeSrcVec, 0);
}
