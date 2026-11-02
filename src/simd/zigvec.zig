//! The target's preferred SIMD widths for `f32` and `f64`, and their vector
//! types. Fixed-size vectors (`V2f32` and the rest) are `linalg.vector`'s.

const rsimd = @import("__root__.zig");

pub const f32_simd_lanes = rsimd.vectorLen(f32);
pub const SimdVecF32 = @Vector(f32_simd_lanes, f32);

pub const f64_simd_lanes = rsimd.vectorLen(f64);
pub const SimdVecF64 = @Vector(f64_simd_lanes, f64);
