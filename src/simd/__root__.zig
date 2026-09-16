//! Generic vector-math routines shared by every consumer that works with
//! `@Vector` types, regardless of component count or scalar type - CPU-side
//! geometry math builds its domain-specific API on top of these instead of
//! keeping its own copy of the same handful of formulas.

const std = @import("std");

pub const zigvec = @import("zigvec.zig");
pub const reduce = @import("reduce.zig");

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
