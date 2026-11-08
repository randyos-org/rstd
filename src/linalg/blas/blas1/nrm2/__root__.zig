//! `NRM2`: the Euclidean norm of `x`. The portable kernels (`*_generic`)
//! and tests of every core's kernels against a plain loop.

const types = @import("../../types.zig");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

const BlasLong = types.BlasLong;
const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

/// The Euclidean norm of `n` elements of `w` floats. Up to 32 bits, the
/// plain sum of squares in `f64`, where every square is exact and none
/// overflows or underflows; at 64 bits, Blue's algorithm (`util.blue`). Wider,
/// x86_64 `nrm2.S`'s and `znrm2.S`'s plain sum of squares (`util.extendedSum`),
/// unscaled: it overflows where they do. Only a zero stride is refused: a
/// negative one walks backwards.
fn Nrm2(comptime F: type, comptime w: comptime_int) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong) callconv(types.kernelCc(F)) F {
            if (n <= 0 or inc_x == 0) return 0;
            if (comptime types.isExtended(F)) return @sqrt(util.extendedSum(F, w, .square, n, x, inc_x));
            if (w == 1 and n == 1) return @abs(x[0]);
            if (comptime @bitSizeOf(F) > 32) return util.blue(F, w, n, x, inc_x);
            var s: f64 = 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                inline for (0..w) |c| {
                    const v: f64 = util.at(F, x, i * inc_x * w + c).*;
                    s += v * v;
                }
            }
            return @floatCast(@sqrt(s));
        }
    };
}

pub const snrm2_generic = &Nrm2(f32, 1).k;
pub const dnrm2_generic = &Nrm2(f64, 1).k;
pub const qnrm2_generic = &Nrm2(XDouble, 1).k;
pub const cnrm2_generic = &Nrm2(f32, 2).k;
pub const znrm2_generic = &Nrm2(f64, 2).k;
pub const xnrm2_generic = &Nrm2(XDouble, 2).k;

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    const w = if (types.isComplex(T)) 2 else 1;
    return struct {
        pub const nrm2_k = &Nrm2(F, w).k;
    };
}

/// Points every nrm2 slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.snrm2_k = snrm2_generic;
    vt.dnrm2_k = dnrm2_generic;
    vt.qnrm2_k = qnrm2_generic;
    vt.cnrm2_k = cnrm2_generic;
    vt.znrm2_k = znrm2_generic;
    vt.xnrm2_k = xnrm2_generic;
}

// ---- tests -----------------------------------------------------------------

const zigstd = @import("../../../../zigstd.zig");
const std = @import("std");
const check = @import("../../_test/check.zig");
const bitwise = @import("../../_test/bitwise.zig");
const pattern = @import("../../_test/pattern.zig");
const reference = @import("../../_test/reference.zig");

/// `kern` against a plain loop in `f64`; `T` is the element type. A negative
/// stride starts from the last element, as the interface passes it.
fn checkNrm2(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    for (check.sizes) |n| {
        for ([_]isize{ 1, 2, 3, -1, -2 }) |s| {
            const inc = @abs(s);
            var x: [check.capacity]T = undefined;
            check.fill(T, &x, n + 1);
            const W = reference.Wide(F);
            var acc: W = 0;
            for (0..n) |i| {
                const v = x[i * inc];
                acc += if (comptime complex) @as(W, v.re) * v.re + @as(W, v.im) * v.im else @as(W, v) * v;
            }
            const want = reference.sqrtW(acc);
            const start = if (s < 0) (n - 1) * inc else 0;
            const got = kern(@intCast(n), @as([*]F, @ptrCast(&x[start])), s);
            check.expectCloseRef(F, want, got, check.tolerance(F, 2 * n) * @max(1, want)) catch |e| return check.case(e, if (s < 0) "nrm2 (stride negated)" else "nrm2", T, n, .{ inc, 0 });
        }
    }
}

fn checkAll(vt: *const Vtable) !void {
    try checkNrm2(f32, vt.snrm2_k.?);
    try checkNrm2(f64, vt.dnrm2_k.?);
    try checkNrm2(std.math.Complex(f32), vt.cnrm2_k.?);
    try checkNrm2(std.math.Complex(f64), vt.znrm2_k.?);
}

/// Magnitudes by element, as powers of two: a cycle through `exps`.
const Range = struct { name: []const u8, exps: []const i32 };

/// The exponents of the ranges `accuracy` takes `F` through: denormal, past
/// the square root of the smallest and largest normal, near the largest, and
/// mixed, which take each branch of Blue's algorithm.
fn ranges(comptime F: type) [9]Range {
    const lo = std.math.floatExponentMin(F);
    const hi = std.math.floatExponentMax(F);
    const dn = lo - 4;
    const tiny = @divFloor(lo, 2) - 4;
    const big = @divFloor(hi, 2) + 4;
    return .{
        .{ .name = "medium", .exps = &.{0} },
        .{ .name = "denormal", .exps = &.{dn} },
        .{ .name = "tiny", .exps = &.{tiny} },
        .{ .name = "big", .exps = &.{big} },
        .{ .name = "near the largest", .exps = &.{hi - 3} },
        .{ .name = "denormal and medium", .exps = &.{ dn, 0 } },
        .{ .name = "medium and big", .exps = &.{ 0, big } },
        .{ .name = "denormal and big", .exps = &.{ dn, big } },
        .{ .name = "tiny, medium and big", .exps = &.{ tiny, 0, big } },
    };
}

/// `kern` against an `f128` reference over `ranges`, at sizes up to 2^14
/// and strides negative too: within one ulp for `s` and `c`, whose squares
/// sum exactly in `f64`, and `n w / 2 + 2` for `d` and `z` (Blue's, or the
/// `@Vector` plain sum where it is safe); infinite past the range.
fn accuracy(comptime T: type, kern: anytype, al: std.mem.Allocator) !void {
    const F = check.FloatOf(T);
    const w = if (comptime check.isComplex(T)) 2 else 1;
    const long: usize = if (check.rig()) 1 << 12 else 1 << 14;
    const buf = try al.alloc(F, w * 3 * long);
    defer al.free(buf);
    for (ranges(F)) |r| for ([_]usize{ 1, 2, 3, 7, 16, 33, 100, long }) |n| for ([_]isize{ 1, 3, -2 }) |s| {
        if (n == long and s != 1) continue;
        const inc = @abs(s);
        check.fill(F, buf[0 .. w * n * inc], n + r.exps.len);
        var acc: f128 = 0;
        for (0..n) |i| for (0..w) |c| {
            const at = (i * inc) * w + c;
            buf[at] = std.math.ldexp(buf[at], r.exps[i % r.exps.len]);
            acc += @as(f128, buf[at]) * buf[at];
        };
        const want = reference.sqrtW(acc);
        const start = if (s < 0) (n - 1) * inc * w else 0;
        const got = kern(@intCast(n), buf.ptr + start, s);
        const ulps: f128 = if (F == f32) 1 else @floatFromInt(n * w / 2 + 2);
        pattern.component(F, want, got, ulps, want) catch |e| {
            check.say("{s}: ", .{r.name});
            return check.case(e, "nrm2", T, n, .{ inc, 0 });
        };
    };
}

fn checkAccuracy(vt: *const Vtable) !void {
    const al = std.testing.allocator;
    try accuracy(f32, vt.snrm2_k.?, al);
    try accuracy(f64, vt.dnrm2_k.?, al);
    try accuracy(std.math.Complex(f32), vt.cnrm2_k.?, al);
    try accuracy(std.math.Complex(f64), vt.znrm2_k.?, al);
}

test "portable nrm2: within an ulp for s and c, of n/2 for d and z, denormal to near the largest" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAccuracy(&vt);
}

test "every linked core's nrm2: within an ulp for s and c, of n/2 for d and z, denormal to near the largest" {
    try check.eachCore(checkAccuracy);
}

test "portable nrm2 kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkNrm2(XDouble, vt.qnrm2_k.?);
    try checkNrm2(std.math.Complex(XDouble), vt.xnrm2_k.?);
}

test "portable nrm2 kernels at the extended type match a plain loop" {
    try checkNrm2(check.Ext, &Nrm2(check.Ext, 1).k);
    try checkNrm2(std.math.Complex(check.Ext), &Nrm2(check.Ext, 2).k);
}

test "every linked core's nrm2 kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}

test "portable nrm2: infinite for infinities, NaN for a NaN among them" {
    const inf = std.math.inf(f64);
    var x = [_]f64{ 1, -inf, 2, inf, 3 };
    try std.testing.expect(std.math.isPositiveInf(dnrm2_generic(5, &x, 1)));
    var z = [_]f64{ inf, 1, 0, -inf };
    try std.testing.expect(std.math.isPositiveInf(znrm2_generic(2, &z, 1)));
    x[2] = std.math.nan(f64);
    try std.testing.expect(std.math.isNan(dnrm2_generic(5, &x, 1)));
}

test "nrm2 past 64 bits is the plain sum of squares, overflow included; at 64 bits or less, it does not overflow" {
    // by width, so a 64-bit `c_longdouble` of its own takes Blue's algorithm
    inline for (.{ f32, f64, c_longdouble, f80, f128 }) |F| {
        const big = if (comptime types.isExtended(F)) std.math.ldexp(@as(F, 1), 16000) else std.math.floatMax(F) / 2;
        var x = [_]F{ big, big };
        const got = Nrm2(F, 1).k(2, &x, 1);
        if (comptime types.isExtended(F)) {
            try std.testing.expect(std.math.isInf(got));
        } else {
            try std.testing.expectApproxEqRel(big * @sqrt(@as(F, 2)), got, 1e-6);
        }
    }
}

test "the extended nrm2 sums fold in nrm2.S's and znrm2.S's order, bit for bit" {
    inline for (.{ f80, f128 }) |F| {
        var x: [check.capacity * 2]F = undefined;
        for (check.sizes) |n| for ([_]usize{ 1, 3 }) |inc| {
            check.fill(F, &x, n + 11);
            const len: BlasLong = @intCast(n);
            const li: BlasLong = @intCast(inc);
            inline for (.{ 1, 2 }) |w| {
                try bitwise.expectBits(F, check.plainExtendedSum(F, w, .square, n, &x, inc), util.extendedSum(F, w, .square, len, &x, li), "nrm2", n, inc);
            }
        };
    }
}

test "the extended sum of squares at f80 carries the full significand" {
    // 1 + 2^-63 squared is 1 + 2^-62 in f80 and 1 in f64
    var x = [_]f80{ 1 + std.math.ldexp(@as(f80, 1), -63), 0, 0, 0, 0 };
    try zigstd.testExpectEqual(1 + std.math.ldexp(@as(f80, 1), -62), util.extendedSum(f80, 1, .square, 5, &x, 1));
}

test "the extended nrm2 sums give the values traced by hand from nrm2.S and znrm2.S" {
    // as asum's (`asum/`), with squares: 2^32 squared is 2^64, and 2^64 + 1
    // rounds to 2^64 in f80
    const F = f80;
    const t = std.math.ldexp(@as(F, 1), 64);
    const h = std.math.ldexp(@as(F, 1), 32);
    var nx = [_]F{ 1, 1, h, 0, 0, 0, 0, 0, 1 };
    try zigstd.testExpectEqual(t, util.extendedSum(F, 1, .square, 8, &nx, 1));
    try zigstd.testExpectEqual(t + 4, util.extendedSum(F, 1, .square, 9, &nx, 1));
    var znx = [_]F{ 1, 1, h, 0, 0, 0, 0, 0, 0, 1 };
    try zigstd.testExpectEqual(t + 2, util.extendedSum(F, 2, .square, 5, &znx, 1));
}

test "portable nrm2 kernels neither overflow nor underflow" {
    inline for (.{ .{ f32, snrm2_generic }, .{ f64, dnrm2_generic } }) |p| {
        const T = p[0];
        const big: T = std.math.floatMax(T) / 2;
        const tiny: T = std.math.floatMin(T);
        var x = [_]T{ big, big, big };
        try check.expectClose(f64, @as(f64, big) * @sqrt(3.0), p[1](3, &x, 1), @as(f64, big) * 1e-6);
        x = .{ tiny, tiny, tiny };
        try check.expectClose(f64, @as(f64, tiny) * @sqrt(3.0), p[1](3, &x, 1), @as(f64, tiny) * 1e-6);
        x = .{ 0, 0, 0 };
        try std.testing.expectEqual(@as(T, 0), p[1](3, &x, 1));
    }
}
