//! Tests of `../vector.zig`.

const std = @import("std");
const fma = @import("../../kernels/blas/fma.zig");
const vroot = @import("../vector.zig");

const Elem = vroot.Elem;
const V2f32 = vroot.V2f32;
const V2f64 = vroot.V2f64;
const V3f32 = vroot.V3f32;
const V3f64 = vroot.V3f64;
const V4f32 = vroot.V4f32;
const V4f64 = vroot.V4f64;
const Vec = vroot.Vec;
const count = vroot.count;
const cross = vroot.cross;
const cross2 = vroot.cross2;
const distance = vroot.distance;
const dot = vroot.dot;
const dot2 = vroot.dot2;
const interleave = vroot.interleave;
const length = vroot.length;
const max = vroot.max;
const min = vroot.min;
const mix = vroot.mix;
const mulAdd = vroot.mulAdd;
const norm = vroot.norm;
const prod = vroot.prod;
const reverse = vroot.reverse;
const reversePairs = vroot.reversePairs;
const splat = vroot.splat;
const sum = vroot.sum;
const swapHalves = vroot.swapHalves;
const uninterleave = vroot.uninterleave;
const zero = vroot.zero;

const zigstd = @import("../../../zigstd.zig");

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

const expectSameBits = @import("../../kernels/blas/_test/bitwise.zig").expectSameBits;

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

test "pffft's lane shuffles: SSE's lane orders at 4, and inverse pairs at 2 and 8" {
    const a: V4f32 = .{ 0, 1, 2, 3 };
    const b: V4f32 = .{ 4, 5, 6, 7 };
    // _mm_unpacklo_ps / _mm_unpackhi_ps
    const il = interleave(a, b);
    try zigstd.testExpectEqual(V4f32{ 0, 4, 1, 5 }, il[0]);
    try zigstd.testExpectEqual(V4f32{ 2, 6, 3, 7 }, il[1]);
    // _MM_SHUFFLE(2,0,2,0) and (3,1,3,1)
    const un = uninterleave(a, b);
    try zigstd.testExpectEqual(V4f32{ 0, 2, 4, 6 }, un[0]);
    try zigstd.testExpectEqual(V4f32{ 1, 3, 5, 7 }, un[1]);
    // _mm_shuffle_ps(b, a, _MM_SHUFFLE(3,2,1,0))
    try zigstd.testExpectEqual(V4f32{ 4, 5, 2, 3 }, swapHalves(a, b));
    try zigstd.testExpectEqual(V4f32{ 3, 2, 1, 0 }, reverse(a));
    try zigstd.testExpectEqual(V4f32{ 2, 3, 0, 1 }, reversePairs(a));
    inline for (.{ Vec(2, f64), Vec(8, f32) }) |V| {
        const x: V = std.simd.iota(Elem(V), count(V));
        const y = x + @as(V, @splat(100));
        const back = uninterleave(interleave(x, y)[0], interleave(x, y)[1]);
        try zigstd.testExpectEqual(x, back[0]);
        try zigstd.testExpectEqual(y, back[1]);
        try zigstd.testExpectEqual(x, reverse(reverse(x)));
        try zigstd.testExpectEqual(x, reversePairs(reversePairs(x)));
    }
}

test "max and min: the first NaN wins, ties keep the first" {
    const nan1: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_0001));
    const nan2: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_0002));
    try expectSameBits(nan1, max(V4f64{ 1, nan1, 9, nan2 }));
    try expectSameBits(nan1, min([4]f64{ 1, nan1, -9, nan2 }));
    try expectSameBits(@as(f64, -0.0), max(V2f64{ -0.0, 0.0 }));
    try expectSameBits(@as(f64, 0.0), min(V2f64{ 0.0, -0.0 }));
}
