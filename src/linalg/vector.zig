//! Generic vector-math routines shared by every consumer that works with
//! `@Vector` types, regardless of component count or scalar type - CPU-side
//! geometry math builds its domain-specific API on top of these instead of
//! keeping its own copy of the same handful of formulas.

const std = @import("std");
const zigstd = @import("../zigstd.zig");
const rsimd = @import("../simd/__root__.zig");
const zigvec = rsimd.zigvec;

const splatAs = zigvec.splatAs;
const zeroAs = zigvec.zeroAs;
const V2f32 = zigvec.V2f32;
const V3f32 = zigvec.V3f32;
const V4f32 = zigvec.V4f32;
const V2f64 = zigvec.V2f64;
const V3f64 = zigvec.V3f64;
const V4f64 = zigvec.V4f64;

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

/// Euclidean norm. Tries the cheap `sqrt(dot2(v))` path unconditionally
/// first and falls back to std.math.hypot if overflow/underflow.
pub inline fn length(v: anytype) @TypeOf(v[0]) {
    const T = @TypeOf(v[0]);
    const sq = dot2(v);
    const looks_fine = std.math.isFinite(sq) and
        (sq != 0 or !@reduce(.Or, v != zeroAs(v)));
    if (looks_fine) return @sqrt(sq);

    const len = @typeInfo(@TypeOf(v)).vector.len;
    var acc: T = @abs(v[0]);
    inline for (1..len) |i| acc = std.math.hypot(acc, v[i]);
    return acc;
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

test "cross3 matches the standard basis identity x cross y = z" {
    const x: V3f32 = .{ 1, 0, 0 };
    const y: V3f32 = .{ 0, 1, 0 };
    try zigstd.testExpectEqual(V3f32{ 0, 0, 1 }, cross3(x, y));
}

test "cross sanity" {
    const ret = cross3(V3f32{ 1, 2, 3 }, V3f32{ 5, 11, 34 });
    try zigstd.testExpectApproxEqAbs(35, ret[0], 0.0001);
    try zigstd.testExpectApproxEqAbs(-19, ret[1], 0.0001);
    try zigstd.testExpectApproxEqAbs(1, ret[2], 0.0001);
}

test "normalizeSafe falls back on the zero vector instead of producing NaN" {
    const zero: V2f64 = .{ 0, 0 };
    const fallback: V2f64 = .{ 0, 1 };
    try zigstd.testExpectEqual(fallback, norm(zero, fallback));

    const v: V2f64 = .{ 0, 5 };
    try zigstd.testExpectEqual(V2f64{ 0, 1 }, norm(v, fallback));
}
