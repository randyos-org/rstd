//! `SUM`: the sum of `x` (`re + im` for complex). The portable kernels
//! (`*_generic`) and tests of every core's kernels against a plain loop.

const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

const BlasLong = types.BlasLong;
const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

/// `kernel/arm/sum.c` and `zsum.c`: the sum of `x`, of `re + im` for complex
/// (unlike `asum`, no absolute value); an element is `w` floats. Up to 32
/// bits the sum is `f64`, as in `asum`.
fn Sum(comptime F: type, comptime w: comptime_int) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong) callconv(types.kernelCc(F)) F {
            if (n <= 0 or inc_x <= 0) return 0;
            var sum: if (@bitSizeOf(F) <= 32) f64 else F = 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                inline for (0..w) |c| sum += util.at(F, x, i * inc_x * w + c).*;
            }
            return @floatCast(sum);
        }
    };
}

pub const ssum_generic = &Sum(f32, 1).k;
pub const dsum_generic = &Sum(f64, 1).k;
pub const qsum_generic = &Sum(XDouble, 1).k;
pub const csum_generic = &Sum(f32, 2).k;
pub const zsum_generic = &Sum(f64, 2).k;
pub const xsum_generic = &Sum(XDouble, 2).k;

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    const w = if (types.isComplex(T)) 2 else 1;
    return struct {
        pub const sum_k = &Sum(F, w).k;
    };
}

/// Points every sum slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.ssum_k = ssum_generic;
    vt.dsum_k = dsum_generic;
    vt.qsum_k = qsum_generic;
    vt.csum_k = csum_generic;
    vt.zsum_k = zsum_generic;
    vt.xsum_k = xsum_generic;
}

// ---- tests -----------------------------------------------------------------

const std = @import("std");
const check = @import("../../_test/check.zig");
const bitwise = @import("../../_test/bitwise.zig");
const reference = @import("../../_test/reference.zig");

/// `kern` against a plain loop in `f64`; `T` is the element type.
fn checkSum(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    for (check.sizes) |n| {
        for ([_]usize{ 1, 2, 3 }) |inc| {
            var x: [check.capacity]T = undefined;
            check.fill(T, &x, n + 1);
            const W = reference.Wide(F);
            var acc: W = 0;
            for (0..n) |i| {
                const v = x[i * inc];
                acc += if (comptime complex) @as(W, v.re) + @as(W, v.im) else @as(W, v);
            }
            const want = acc;
            const got = kern(@intCast(n), @as([*]F, @ptrCast(&x)), @intCast(inc));
            check.expectCloseRef(F, want, got, check.tolerance(F, 2 * n)) catch |e| return check.case(e, "sum", T, n, .{ inc, 0 });
        }
    }
}

fn checkAll(vt: *const Vtable) !void {
    var failed: ?anyerror = null;
    check.keep(&failed, checkSum(f32, vt.ssum_k.?));
    check.keep(&failed, checkSum(f64, vt.dsum_k.?));
    check.keep(&failed, checkSum(std.math.Complex(f32), vt.csum_k.?));
    check.keep(&failed, checkSum(std.math.Complex(f64), vt.zsum_k.?));
    if (failed) |e| return e;
}

test "portable sum kernels match a plain loop" {
    const vt = check.portableVtable();
    try checkAll(vt);
    try checkSum(XDouble, vt.qsum_k.?);
    try checkSum(std.math.Complex(XDouble), vt.xsum_k.?);
}

test "portable sum kernels at the extended type match a plain loop" {
    try checkSum(check.Ext, &Sum(check.Ext, 1).k);
    try checkSum(std.math.Complex(check.Ext), &Sum(check.Ext, 2).k);
}

test "every linked core's sum kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}

/// `vector.zig`'s `sum` order written out plainly over the `w * n` floats
/// (`w` 2 for complex): unit stride, float `i` of each full `2 * W` into
/// lane `i % W` of sum `i / W % 2`; strided, element `i` of each full group
/// into lane `i % 4` (real), or into sum `i % 2` by parts (complex). Then the sums lane by lane, the lanes and
/// the rest in order, in `f64` for `f32`.
fn vectorOrder(comptime F: type, comptime w: comptime_int, n: usize, x: []const F, inc: usize) F {
    const D = if (F == f32) f64 else F;
    const W = blas.vector.Real(F, 4).W;
    if (inc == 1) {
        const len = w * n;
        var acc: [2][W]D = @splat(@splat(0));
        var i: usize = 0;
        while (i + 2 * W <= len) : (i += 2 * W) for (0..2) |q| for (0..W) |k| {
            acc[q][k] += x[i + q * W + k];
        };
        var s: D = acc[0][0] + acc[1][0];
        for (1..W) |k| s += acc[0][k] + acc[1][k];
        while (i < len) : (i += 1) s += x[i];
        return @floatCast(s);
    }
    if (w == 1) {
        var acc: [4]D = @splat(0);
        var i: usize = 0;
        while (i + 4 <= n) : (i += 4) for (0..4) |k| {
            acc[k] += x[(i + k) * inc];
        };
        var s: D = acc[0];
        for (1..4) |k| s += acc[k];
        while (i < n) : (i += 1) s += x[i * inc];
        return @floatCast(s);
    }
    var acc: [2][2]D = @splat(@splat(0));
    for (0..n) |i| for (0..2) |c| {
        acc[i % 2][c] += x[2 * i * inc + c];
    };
    return @floatCast((acc[0][0] + acc[1][0]) + (acc[0][1] + acc[1][1]));
}

fn checkVectorSum(comptime F: type, comptime w: comptime_int) !void {
    const K = if (w == 2) blas.vector.Complex(F, 4) else blas.vector.Real(F, 4);
    const nan = std.math.nan(F);
    const inf = std.math.inf(F);
    const tiny = std.math.floatTrueMin(F);
    const specials = [_][2]F{ .{ nan, 1 }, .{ inf, 1 }, .{ inf, -inf }, .{ tiny, -tiny }, .{ tiny, 3 * tiny }, .{ std.math.floatMax(F), std.math.floatMax(F) } };
    var x: [w * check.capacity]F = undefined;
    for ([_]usize{ 0, 1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33, 65, 100, 257 }) |n| {
        for (0..specials.len + 2) |p| for ([_]usize{ 1, 2, 3 }) |inc| {
            var prng = std.Random.DefaultPrng.init(n * 31 + p);
            const r = prng.random();
            // magnitudes over many binades, so the order shows in the bits;
            // past `specials`, all denormals, then none
            for (&x) |*v| v.* = if (p == specials.len) tiny * @as(F, @floatFromInt(r.intRangeAtMost(i32, -99, 99))) else (r.float(F) - 0.5) * std.math.ldexp(@as(F, 1), r.intRangeAtMost(i32, -20, 20));
            if (p < specials.len and n > 0) for (specials[p]) |v| {
                x[w * inc * r.uintLessThan(usize, n) + r.uintLessThan(usize, w)] = v;
            };
            const want = vectorOrder(F, w, n, &x, inc);
            const got = K.sum(@intCast(n), &x, @intCast(inc));
            if (!(std.math.isNan(want) and std.math.isNan(got))) try bitwise.expectBits(F, want, got, if (w == 2) "vector sum, complex" else "vector sum", n, inc);
            // a negative or zero stride sums nothing, as the portable kernels
            for ([_]BlasLong{ 0, -1, -2 }) |neg| try std.testing.expectEqual(@as(F, 0), K.sum(@intCast(n), &x, neg));
        };
    }
}

test "the @Vector sum kernels add in their stated order, bit for bit, special values and denormals too" {
    // not where lanes flush, which no core's kernels use, nor on the x87,
    // whose bits hang on spills
    inline for (.{ f32, f64 }) |F| {
        if (comptime util.vectorIeee(F) and !util.onX87(f64)) inline for (.{ 1, 2 }) |w| try checkVectorSum(F, w);
    }
}

test "the @Vector f32 sum rounds once: 4096 thirds sum exactly in f64" {
    if (comptime !util.vectorIeee(f32)) return;
    var x: [4096]f32 = @splat(1.0 / 3.0);
    const want: f32 = @floatCast(@as(f64, x.len) * x[0]);
    try std.testing.expectEqual(want, blas.vector.Real(f32, 4).sum(x.len, &x, 1));
    try std.testing.expectEqual(want, blas.vector.Complex(f32, 4).sum(x.len / 2, &x, 1));
}
