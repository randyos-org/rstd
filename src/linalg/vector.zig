//! `rstd.linalg`'s fixed-size vector math (design.md, "Planned: the public
//! `rstd.linalg` layer"): free functions over a `Vec(N, T)` or an `[N]T`, `T`
//! a float or, for arrays, its `std.math.Complex` (interleaved), inferring `N`
//! and `T`. All inline; shaderlib re-exports them under GLSL names.

const std = @import("std");
const sysinfo = @import("builtin");
const fma = @import("blas/fma.zig");

/// SPIR-V targets (rgl's shaders -- see rgl/shaders/) don't support f64,
/// which `std.math.hypot`'s overflow/underflow-safe fallback in `length`
/// below needs internally. Shader inputs (pixel-space positions,
/// normalized directions, SDF distances) never approach the extremes that
/// fallback exists for, so `length` skips it entirely there rather than
/// failing to compile.
const on_spirv = switch (sysinfo.target.cpu.arch) {
    .spirv32, .spirv64 => true,
    else => false,
};

/// `@Vector(N, T)`, real `T` only; a complex fixed-size vector is `[N]Complex(T)`.
pub fn Vec(comptime N: comptime_int, comptime T: type) type {
    return @Vector(N, T);
}

pub const V2f32 = Vec(2, f32);
pub const V3f32 = Vec(3, f32);
pub const V4f32 = Vec(4, f32);
pub const V2f64 = Vec(2, f64);
pub const V3f64 = Vec(3, f64);
pub const V4f64 = Vec(4, f64);

/// The element type of a `Vec` or an array.
pub fn Elem(comptime V: type) type {
    return switch (@typeInfo(V)) {
        .vector => |info| info.child,
        .array => |info| info.child,
        else => @compileError("expects a Vec or an array, found " ++ @typeName(V)),
    };
}

/// The length of a `Vec` or an array.
pub fn count(comptime V: type) comptime_int {
    return switch (@typeInfo(V)) {
        .vector => |info| info.len,
        .array => |info| info.len,
        else => @compileError("expects a Vec or an array, found " ++ @typeName(V)),
    };
}

/// Whether `T` is a `std.math.Complex`.
pub fn isComplex(comptime T: type) bool {
    return @typeInfo(T) == .@"struct" and @hasField(T, "re") and
        T == std.math.Complex(@FieldType(T, "re"));
}

/// The real type under element type `T`.
pub fn Real(comptime T: type) type {
    return if (isComplex(T)) @FieldType(T, "re") else T;
}

/// `Σ a b` over the real parts in storage order (`re`, then `im`, per complex
/// element), each product after the first fused into the sum as `mulAdd` says.
inline fn partDot(a: anytype, b: @TypeOf(a)) Real(Elem(@TypeOf(a))) {
    const R = Real(Elem(@TypeOf(a)));
    if (comptime !isComplex(Elem(@TypeOf(a)))) {
        var acc: R = a[0] * b[0];
        inline for (1..count(@TypeOf(a))) |i| acc = mulAdd(R, a[i], b[i], acc);
        return acc;
    }
    var acc: R = a[0].re * b[0].re;
    acc = mulAdd(R, a[0].im, b[0].im, acc);
    inline for (1..count(@TypeOf(a))) |i| {
        acc = mulAdd(R, a[i].re, b[i].re, acc);
        acc = mulAdd(R, a[i].im, b[i].im, acc);
    }
    return acc;
}

/// Every element `x`.
pub inline fn splat(comptime V: type, x: Elem(V)) V {
    return @splat(x);
}

/// Every element 0.
pub inline fn zero(comptime V: type) V {
    const E = Elem(V);
    return splat(V, if (comptime isComplex(E)) E.init(0, 0) else 0);
}

/// `Σ conj(a[i]) b[i]`: `aᴴ b` for complex, as MATLAB's `dot`. In element
/// order, each product after the first fused into the sum as `mulAdd` says;
/// a complex real part adds `re re` then `im im`, its imaginary part `re im`
/// then `-im re`.
pub inline fn dot(a: anytype, b: @TypeOf(a)) Elem(@TypeOf(a)) {
    const T = Elem(@TypeOf(a));
    if (comptime !isComplex(T)) return partDot(a, b);
    const R = Real(T);
    var im: R = a[0].re * b[0].im;
    im = mulAdd(R, -a[0].im, b[0].re, im);
    inline for (1..count(@TypeOf(a))) |i| {
        im = mulAdd(R, a[i].re, b[i].im, im);
        im = mulAdd(R, -a[i].im, b[i].re, im);
    }
    return .{ .re = partDot(a, b), .im = im };
}

/// `dot(v, v)`'s real part, `Σ |v[i]|²`, as the real type.
pub inline fn dot2(v: anytype) Real(Elem(@TypeOf(v))) {
    return partDot(v, v);
}

/// The 2-norm (MATLAB's `norm`, GLSL's `length`). Tries `sqrt(dot2(v))`
/// first and falls back to a `std.math.hypot` fold on overflow or underflow
/// (except on SPIR-V, which can't use that fallback).
pub inline fn length(v: anytype) Real(Elem(@TypeOf(v))) {
    const sq = dot2(v);
    const looks_fine = on_spirv or (std.math.isFinite(sq) and (sq != 0 or isZero(v)));
    if (looks_fine) return @sqrt(sq);

    if (comptime !isComplex(Elem(@TypeOf(v)))) {
        var acc = @abs(v[0]);
        inline for (1..count(@TypeOf(v))) |i| acc = std.math.hypot(acc, v[i]);
        return acc;
    }
    var acc = std.math.hypot(@abs(v[0].re), v[0].im);
    inline for (1..count(@TypeOf(v))) |i| acc = std.math.hypot(std.math.hypot(acc, v[i].re), v[i].im);
    return acc;
}

inline fn isZero(v: anytype) bool {
    const complex = comptime isComplex(Elem(@TypeOf(v)));
    inline for (0..count(@TypeOf(v))) |i| {
        if (if (complex) v[i].re != 0 or v[i].im != 0 else v[i] != 0) return false;
    }
    return true;
}

/// `length(b - a)`.
pub inline fn distance(a: anytype, b: @TypeOf(a)) Real(Elem(@TypeOf(a))) {
    return length(diff(b, a));
}

/// `a - b`, for arrays too.
inline fn diff(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    if (@typeInfo(@TypeOf(a)) == .vector) return a - b;
    var r: @TypeOf(a) = undefined;
    inline for (0..count(@TypeOf(a))) |i| {
        r[i] = if (comptime isComplex(Elem(@TypeOf(a)))) a[i].sub(b[i]) else a[i] - b[i];
    }
    return r;
}

/// Normalizes `v` to `v / length(v)` (GLSL's `normalize`; MATLAB's `norm` is
/// `length`). A fallback is the caller's: `norm(v) catch fallback`.
pub inline fn norm(v: anytype) error{ZeroVector}!@TypeOf(v) {
    const V = @TypeOf(v);
    const len = length(v);
    if (len == 0) return error.ZeroVector;
    if (@typeInfo(V) == .vector) return v / splat(V, len);
    var r: V = undefined;
    inline for (0..count(V)) |i| {
        r[i] = if (comptime isComplex(Elem(V))) .{ .re = v[i].re / len, .im = v[i].im / len } else v[i] / len;
    }
    return r;
}

/// 2D perp-dot product ("2D cross product"): the scalar a 3D cross product's
/// z-component would be if both inputs were embedded in the xy-plane.
/// Positive when `b` is counter-clockwise from `a`.
pub inline fn cross2(a: anytype, b: @TypeOf(a)) Elem(@TypeOf(a)) {
    return a[0] * b[1] - a[1] * b[0];
}

/// Standard 3D cross product.
pub inline fn cross(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
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

/// Multiply-add, fused or not as the root's `std_float_options` says
/// (`zigcfg.FloatOptions`; by default where the target fuses in hardware),
/// as the fork's `@Matrix` products are (`linalg/blas/fma.zig`'s `mulAdd`).
pub inline fn mulAdd(comptime T: type, a: T, b: T, c: T) T {
    return fma.mulAdd(T, a, b, c);
}

/// `Σ v[i]`, in element order.
pub inline fn sum(v: anytype) Elem(@TypeOf(v)) {
    var acc = v[0];
    inline for (1..count(@TypeOf(v))) |i| {
        acc = if (comptime isComplex(Elem(@TypeOf(v)))) acc.add(v[i]) else acc + v[i];
    }
    return acc;
}

/// `Π v[i]`, in element order.
pub inline fn prod(v: anytype) Elem(@TypeOf(v)) {
    var acc = v[0];
    inline for (1..count(@TypeOf(v))) |i| {
        acc = if (comptime isComplex(Elem(@TypeOf(v)))) acc.mul(v[i]) else acc * v[i];
    }
    return acc;
}

/// The largest element, or the first NaN (the BLAS's `max`: `@reduce(.Max)`
/// skips NaNs); ties keep the first. Under SPIR-V an `FMax` chain, whose NaN
/// result is undefined.
pub inline fn max(v: anytype) Elem(@TypeOf(v)) {
    return extreme(true, v);
}

/// The smallest element, or the first NaN, as `max`.
pub inline fn min(v: anytype) Elem(@TypeOf(v)) {
    return extreme(false, v);
}

inline fn extreme(comptime largest: bool, v: anytype) Elem(@TypeOf(v)) {
    if (comptime isComplex(Elem(@TypeOf(v)))) @compileError("max and min order real elements only");
    var best = v[0];
    if (on_spirv) {
        inline for (1..count(@TypeOf(v))) |i| best = if (largest) @max(best, v[i]) else @min(best, v[i]);
        return best;
    }
    if (std.math.isNan(best)) return best;
    inline for (1..count(@TypeOf(v))) |i| {
        const x = v[i];
        if (std.math.isNan(x)) return x;
        if (if (largest) x > best else x < best) best = x;
    }
    return best;
}

const zigstd = @import("../zigstd.zig");

test "dot sanity" {
    const ret = dot(V3f32{ 1, 2, 3 }, V3f32{ 4, 5, 6 });
    try zigstd.testExpectApproxEqAbs(32, ret, 0.001);
}

test "dot/length/distance agree with the obvious 3-4-5 triangle" {
    const a: V2f32 = .{ 0, 0 };
    const b: V2f32 = .{ 3, 4 };
    try zigstd.testExpectEqual(@as(f32, 25), dot2(b - a));
    try zigstd.testExpectEqual(@as(f32, 5), length(b - a));
    try zigstd.testExpectEqual(@as(f32, 5), distance(a, b));
}

test "length doesn't overflow to inf on huge components" {
    // sqrt(dot2(v)) would square 1e300 to 1e600, which overflows f64
    // (max ~1.8e308) to +inf; the hypot-fold never squares anything that
    // large without going through hypot's own overflow-safe rescale first.
    const v: V2f64 = .{ 1e300, 1e300 };
    const want = @sqrt(2.0) * 1e300;
    try zigstd.testExpectApproxEqAbs(want, length(v), want * 1e-12);
    try zigstd.testExpectTrue(std.math.isFinite(length(v)));
}

test "length doesn't underflow to 0 on tiny components" {
    // sqrt(dot2(v)) would square 1e-200 to 1e-400, which underflows f64
    // (smallest subnormal ~4.9e-324) to exactly 0, losing the vector
    // entirely; the hypot-fold stays nonzero.
    const v: V2f64 = .{ 1e-200, 1e-200 };
    const want = @sqrt(2.0) * 1e-200;
    try zigstd.testExpectApproxEqAbs(want, length(v), want * 1e-12);
    try zigstd.testExpectTrue(length(v) > 0);
}

test "length stays safe across a 3-component fold with mixed extreme magnitudes" {
    // Exercises the fold across more than one hypot call (unlike the 2D
    // cases above, which only ever call hypot once): a huge component, a
    // tiny one, and an ordinary one all combined in the same vector.
    const v: V3f64 = .{ 1e300, 1e-200, 3.0 };
    const want = 1e300; // the other two components are negligible next to it
    try zigstd.testExpectApproxEqAbs(want, length(v), want * 1e-9);
    try zigstd.testExpectTrue(std.math.isFinite(length(v)));
}

test "cross2 is positive for a counter-clockwise turn" {
    const a: V2f64 = .{ 1, 0 };
    const b: V2f64 = .{ 0, 1 };
    try zigstd.testExpectEqual(@as(f64, 1), cross2(a, b));
}

test "cross matches the standard basis identity x cross y = z" {
    const x: V3f32 = .{ 1, 0, 0 };
    const y: V3f32 = .{ 0, 1, 0 };
    try zigstd.testExpectEqual(V3f32{ 0, 0, 1 }, cross(x, y));
}

test "cross sanity" {
    const ret = cross(V3f32{ 1, 2, 3 }, V3f32{ 5, 11, 34 });
    try zigstd.testExpectApproxEqAbs(35, ret[0], 0.0001);
    try zigstd.testExpectApproxEqAbs(-19, ret[1], 0.0001);
    try zigstd.testExpectApproxEqAbs(1, ret[2], 0.0001);
}

test "norm fails on the zero vector instead of producing NaN" {
    const zero_v: V2f64 = .{ 0, 0 };
    try zigstd.testExpectError(error.ZeroVector, norm(zero_v));
    try zigstd.testExpectEqual(V2f64{ 0, 1 }, try norm(V2f64{ 0, 5 }));
}

test "mulAdd matches plain arithmetic regardless of which path it takes" {
    try zigstd.testExpectEqual(@as(f64, 14.0), mulAdd(f64, 2.0, 3.0, 8.0));
    try zigstd.testExpectEqual(@as(f32, 14.0), mulAdd(f32, 2.0, 3.0, 8.0));
}

/// The formulas as of 4a8caada (`linalg.vector`, `rsimd.mix`), which `N` 2 to
/// 4 must keep bit for bit (font and SDF output).
const Before = struct {
    fn dot(a: anytype, b: @TypeOf(a)) @TypeOf(a[0]) {
        const T = @TypeOf(a[0]);
        var acc: T = a[0] * b[0];
        inline for (1..@typeInfo(@TypeOf(a)).vector.len) |i| acc = fma.mulAdd(T, a[i], b[i], acc);
        return acc;
    }

    fn length(v: anytype) @TypeOf(v[0]) {
        const T = @TypeOf(v[0]);
        const sq = Before.dot(v, v);
        const zero_v: @TypeOf(v) = @splat(0);
        if (std.math.isFinite(sq) and (sq != 0 or !@reduce(.Or, v != zero_v))) return @sqrt(sq);
        var acc: T = @abs(v[0]);
        inline for (1..@typeInfo(@TypeOf(v)).vector.len) |i| acc = std.math.hypot(acc, v[i]);
        return acc;
    }

    fn norm(v: anytype, zero_fallback: @TypeOf(v)) @TypeOf(v) {
        const len = Before.length(v);
        if (len == 0) return zero_fallback;
        return v / @as(@TypeOf(v), @splat(len));
    }

    fn mix(a: anytype, b: @TypeOf(a), t: @TypeOf(a[0])) @TypeOf(a) {
        const V = @TypeOf(a);
        const weight: V = @splat(t);
        const one: V = @splat(1);
        return fma.mulAdd(V, one - weight, a, weight * b);
    }
};

fn expectSameBits(want: anytype, got: @TypeOf(want)) !void {
    const T = @TypeOf(want);
    const U = switch (@typeInfo(T)) {
        .vector => |info| @Vector(info.len, @Int(.unsigned, @bitSizeOf(info.child))),
        else => @Int(.unsigned, @bitSizeOf(T)),
    };
    try zigstd.testExpectEqual(@as(U, @bitCast(want)), @as(U, @bitCast(got)));
}

fn expectAsBefore(comptime T: type, comptime n: comptime_int, random: std.Random) !void {
    const scales = [_]f64{ 1, 1e-3, 1e3, 1e30, 1e-30, 1e200, 1e-200 };
    const V = Vec(n, T);
    for (0..200) |trial| {
        const s: T = @floatCast(scales[trial % scales.len]);
        var a: V = undefined;
        var b: V = undefined;
        inline for (0..n) |i| {
            a[i] = (random.float(T) * 2 - 1) * s;
            b[i] = (random.float(T) * 2 - 1) * s;
        }
        const t = random.float(T);
        try expectSameBits(Before.dot(a, b), dot(a, b));
        try expectSameBits(Before.dot(a, a), dot2(a));
        try expectSameBits(Before.length(a), length(a));
        try expectSameBits(Before.length(b - a), distance(a, b));
        try expectSameBits(Before.norm(a, zero(V)), norm(a) catch zero(V));
        try expectSameBits(Before.mix(a, b, t), mix(a, b, t));
        // an array holds the same numbers to the same bits
        const arr_a: [n]T = a;
        const arr_b: [n]T = b;
        try expectSameBits(dot(a, b), dot(arr_a, arr_b));
        try expectSameBits(length(a), length(arr_a));
        try expectSameBits(norm(a) catch zero(V), @as(V, norm(arr_a) catch zero([n]T)));
    }
}

test "N 2 to 4 keep the formulas' bits as before" {
    var prng: std.Random.DefaultPrng = .init(0x5eed);
    inline for (.{ f32, f64 }) |T| inline for (2..5) |n| try expectAsBefore(T, n, prng.random());
}

test "complex arrays: dot conjugates its first operand" {
    const C = std.math.Complex(f64);
    const a = [2]C{ .init(1, 2), .init(3, 4) };
    const b = [2]C{ .init(5, 6), .init(7, 8) };
    try zigstd.testExpectEqual(C.init(70, -8), dot(a, b));
    try zigstd.testExpectEqual(@as(f64, 30), dot2(a));
    try zigstd.testExpectEqual(@as(f64, 5), length([2]C{ .init(3, 4), .init(0, 0) }));
    try zigstd.testExpectEqual(@as(f64, 5), distance([1]C{.init(1, 1)}, [1]C{.init(4, 5)}));
    const u = try norm([2]C{ .init(0, 3), .init(4, 0) });
    try zigstd.testExpectEqual(C.init(0, 0.6), u[0]);
    try zigstd.testExpectEqual(C.init(0.8, 0), u[1]);
    try zigstd.testExpectError(error.ZeroVector, norm(zero([3]C)));
    try zigstd.testExpectEqual(C.init(4, 6), sum(a));
    try zigstd.testExpectEqual(C.init(-5, 10), prod(a));
    try zigstd.testExpectEqual(C.init(2, -1), splat([4]C, .init(2, -1))[3]);
}

test "sum, prod, max and min in element order" {
    try zigstd.testExpectEqual(@as(f32, 10), sum(V4f32{ 1, 2, 3, 4 }));
    try zigstd.testExpectEqual(@as(f32, 24), prod([4]f32{ 1, 2, 3, 4 }));
    try zigstd.testExpectEqual(@as(f64, 7), max(V3f64{ -1, 7, 3 }));
    try zigstd.testExpectEqual(@as(f64, -1), min(V3f64{ -1, 7, 3 }));
}

test "max and min: the first NaN wins, ties keep the first" {
    const nan1: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_0001));
    const nan2: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_0002));
    try expectSameBits(nan1, max(V4f64{ 1, nan1, 9, nan2 }));
    try expectSameBits(nan1, min([4]f64{ 1, nan1, -9, nan2 }));
    try expectSameBits(@as(f64, -0.0), max(V2f64{ -0.0, 0.0 }));
    try expectSameBits(@as(f64, 0.0), min(V2f64{ 0.0, -0.0 }));
}

test {
    zigstd.testRefAllDecls(@This());
}
