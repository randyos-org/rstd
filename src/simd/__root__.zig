//! Generic vector-math routines shared by every consumer that works with
//! `@Vector` types, regardless of component count or scalar type - CPU-side
//! geometry math builds its domain-specific API on top of these instead of
//! keeping its own copy of the same handful of formulas.

const std = @import("std");
const sysinfo = @import("builtin");
const blas_fma = @import("../linalg/blas/fma.zig");

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
            return mulAdd(@TypeOf(a, b), 1 - t, a, t * b);
        },
        .vector => {
            const weight: BaseType = if (weight_info == .vector) t else @splat(t);
            const one: BaseType = @splat(1);
            return mulAdd(@TypeOf(a, b), one - weight, a, weight * b);
        },
        else => @compileError("Invalid base type, only floats and their vectors are supported"),
    }
}

/// Multiply-add, fused or not as the build's BLAS policy says (by default:
/// fused when the target has an FMA unit, plain `a*b+c` otherwise, so `@mulAdd`
/// never degrades into a libm `fma()` call). The decision lives in
/// `linalg/blas/fma.zig`, which mirrors how the reference's per-core compile
/// flags decide it.
pub inline fn mulAdd(comptime T: type, a: T, b: T, c: T) T {
    return blas_fma.mulAdd(T, a, b, c);
}

/// Preferred SIMD vector length (element count, not bits) for `T` on the
/// compile target -- e.g. 2 for `f64` on plain SSE2, 8 once AVX-512 is in
/// play. Thin wrapper over `std.simd.suggestVectorLength`, falling back to
/// scalar (1 lane) where it has no suggestion, so every caller gets a
/// usable `@Vector(vectorLen(T), T)` without its own null-handling.
///
/// This is a *compile-time* choice tied to whatever CPU model the build was
/// configured for (Zig's default x86_64 "baseline" model only guarantees
/// SSE2, so this returns 2 for `f64` unless `-Dcpu=` says otherwise) -- not
/// runtime dispatch. A binary built this way hits an illegal instruction on
/// hardware narrower than what it assumes, exactly like a microarch-locked
/// OpenBLAS `TARGET=SKYLAKEX` build. Runtime "detect the CPU and pick a
/// kernel" (`DYNAMIC_ARCH`) is a separate, deliberately deferred problem.
pub fn vectorLen(comptime T: type) comptime_int {
    return comptime std.simd.suggestVectorLength(T) orelse 1;
}

const zigstd = @import("../zigstd.zig");

test "mulAdd matches plain arithmetic regardless of which path it takes" {
    try zigstd.testExpectEqual(@as(f64, 14.0), mulAdd(f64, 2.0, 3.0, 8.0));
    try zigstd.testExpectEqual(@as(f32, 14.0), mulAdd(f32, 2.0, 3.0, 8.0));
}

test {
    zigstd.testRefAllDecls(@This());
}
