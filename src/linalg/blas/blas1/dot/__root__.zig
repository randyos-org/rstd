//! `DOT`: `x . y` (`dotu`/`dotc` for complex, `dsdot` summing `f32` in `f64`).
//! The vtable wrappers, the portable kernels (`*_generic`), and tests of every
//! core's kernels against a plain loop.
//!
//! Argument names are copied from the OpenBLAS kernel sources (see the
//! `CNAME(...)` definitions under `kernel/<arch>/*.c`). Several are named
//! `dummy*` in the reference: kernel-ABI padding slots this shape does not use.

const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");

const BlasLong = types.BlasLong;
const HF16 = types.HF16;
const BF16 = types.BF16;
const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

pub fn sbdot_k(n: BlasLong, x: [*]BF16, inc_x: BlasLong, y: [*]BF16, inc_y: BlasLong) f32 {
    return blas.kernel(.sbdot_k)(n, x, inc_x, y, inc_y);
}
pub fn dsbdot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) f64 {
    return blas.kernel(.dsbdot_k)(n, x, inc_x, y, inc_y);
}
pub fn sdot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) f32 {
    return blas.kernel(.sdot_k)(n, x, inc_x, y, inc_y);
}
pub fn ddot_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong) f64 {
    return blas.kernel(.ddot_k)(n, x, inc_x, y, inc_y);
}
pub fn dsdot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) f64 {
    return blas.kernel(.dsdot_k)(n, x, inc_x, y, inc_y);
}
pub fn qdot_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong) XDouble {
    return blas.kernel(.qdot_k)(n, x, inc_x, y, inc_y);
}
pub fn cdotu_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) Complex32 {
    return blas.kernel(.cdotu_k)(n, x, inc_x, y, inc_y).toComplex();
}
pub fn cdotc_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) Complex32 {
    return blas.kernel(.cdotc_k)(n, x, inc_x, y, inc_y).toComplex();
}
pub fn zdotu_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong) Complex64 {
    return blas.kernel(.zdotu_k)(n, x, inc_x, y, inc_y).toComplex();
}
pub fn zdotc_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong) Complex64 {
    return blas.kernel(.zdotc_k)(n, x, inc_x, y, inc_y).toComplex();
}
pub fn xdotu_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong) ComplexXD {
    return blas.kernel(.xdotu_k)(n, x, inc_x, y, inc_y).toComplex();
}
pub fn xdotc_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong) ComplexXD {
    return blas.kernel(.xdotc_k)(n, x, inc_x, y, inc_y).toComplex();
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

/// `x . y` after `kernel/arm/dot.c`, summed sequentially in `Acc` (`f64` for
/// `sdot` too, so it is at least as accurate as the reference's `float` or
/// `double` products) and returned as `R` (`dsdot` returns the `f64` sum).
/// Not bit-identical to `generic/dot.c`'s unrolled sums; see the readme.
/// Wider than `f64`, it is x86_64 `dot.S`'s x87 sum instead (`extendedDot`).
fn Dot(comptime F: type, comptime Acc: type, comptime R: type) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) callconv(types.kernelCc(F)) R {
            if (comptime types.isExtended(F)) return extendedDot(F, n, x, inc_x, y, inc_y);
            var sum: Acc = 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) sum += @as(Acc, util.at(F, x, i * inc_x).*) * @as(Acc, util.at(F, y, i * inc_y).*);
            return @floatCast(sum);
        }
    };
}

/// `kernel/arm/zdot.c`: `x . y`, or `conj(x) . y` (`dotc`), summed in
/// `types.Acc(F)`. Returned as `CComplex` over a vtable type, else `Complex`.
fn Zdot(comptime F: type, comptime conj: bool) type {
    const A = types.Acc(F);
    const R = if (types.hasSlots(F)) types.CComplex(F) else std.math.Complex(F);
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) callconv(types.kernelCc(F)) R {
            if (comptime types.isExtended(F)) {
                const z = extendedZdot(F, conj, n, x, inc_x, y, inc_y);
                return .{ .re = z[0], .im = z[1] };
            }
            var re: A = 0;
            var im: A = 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const xr: A = util.at(F, x, i * inc_x * 2).*;
                const xi: A = util.at(F, x, i * inc_x * 2 + 1).*;
                const yr: A = util.at(F, y, i * inc_y * 2).*;
                const yi: A = util.at(F, y, i * inc_y * 2 + 1).*;
                if (conj) {
                    re += xr * yr + xi * yi;
                    im -= xi * yr - xr * yi;
                } else {
                    re += xr * yr - xi * yi;
                    im += xi * yr + xr * yi;
                }
            }
            return .{ .re = @floatCast(re), .im = @floatCast(im) };
        }
    };
}

/// x86_64 `dot.S`'s x87 sum, which every float wider than `f64` takes: four
/// accumulators by `i % 4`, the tail into the first, folded
/// `(a3 + a1) + (a2 + a0)`, nothing fused. At `f80` it gives the assembly's
/// bits; `kernel/x86/`'s is in the same order.
fn extendedDot(comptime F: type, n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) F {
    var a: [4]F = @splat(0);
    var i: BlasLong = 0;
    const n4 = n & -4;
    while (i < n4) : (i += 4) {
        inline for (0..4) |c| a[c] += util.at(F, x, (i + c) * inc_x).* * util.at(F, y, (i + c) * inc_y).*;
    }
    while (i < n) : (i += 1) a[0] += util.at(F, x, i * inc_x).* * util.at(F, y, i * inc_y).*;
    return (a[3] + a[1]) + (a[2] + a[0]);
}

/// `zdot.S`'s x87 sum: each of an element's four products in its own
/// accumulator, combined at the end; `conj` conjugates `x`. Real and
/// imaginary parts.
fn extendedZdot(comptime F: type, comptime conj: bool, n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) [2]F {
    var a: [4]F = @splat(0);
    var i: BlasLong = 0;
    while (i < n) : (i += 1) {
        const xr = util.at(F, x, i * inc_x * 2).*;
        const xi = util.at(F, x, i * inc_x * 2 + 1).*;
        const yr = util.at(F, y, i * inc_y * 2).*;
        const yi = util.at(F, y, i * inc_y * 2 + 1).*;
        a[0] += xr * yr;
        a[1] += xr * yi;
        a[2] += xi * yr;
        a[3] += xi * yi;
    }
    return if (conj) .{ a[0] + a[3], a[1] - a[2] } else .{ a[0] - a[3], a[1] + a[2] };
}

pub const sdot_generic = &Dot(f32, f64, f32).k;
pub const ddot_generic = &Dot(f64, f64, f64).k;
pub const dsdot_generic = &Dot(f32, f64, f64).k;
pub const qdot_generic = &Dot(XDouble, XDouble, XDouble).k;
pub const cdotu_generic = &Zdot(f32, false).k;
pub const cdotc_generic = &Zdot(f32, true).k;
pub const zdotu_generic = &Zdot(f64, false).k;
pub const zdotc_generic = &Zdot(f64, true).k;
pub const xdotu_generic = &Zdot(XDouble, false).k;
pub const xdotc_generic = &Zdot(XDouble, true).k;

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    return struct {
        pub const dot_k = &Dot(F, types.Acc(F), F).k;
        pub const dotu_k = &Zdot(F, false).k;
        pub const dotc_k = &Zdot(F, true).k;
    };
}

/// Points every dot slot at the portable kernels, except `sbdot` and
/// `dsbdot`: `sbdot` is `ext/bf16`'s, and `dsbdot` aliases `dsdot` there
/// (`kernel/setparam-ref.c`), as every bfloat16-named level-1 slot does.
pub fn installGeneric(vt: *Vtable) void {
    vt.sdot_k = sdot_generic;
    vt.ddot_k = ddot_generic;
    vt.dsdot_k = dsdot_generic;
    vt.qdot_k = qdot_generic;
    vt.cdotu_k = cdotu_generic;
    vt.cdotc_k = cdotc_generic;
    vt.zdotu_k = zdotu_generic;
    vt.zdotc_k = zdotc_generic;
    vt.xdotu_k = xdotu_generic;
    vt.xdotc_k = xdotc_generic;
}

// ---- tests -----------------------------------------------------------------

const sysinfo = @import("builtin");
const zigstd = @import("../../../../zigstd.zig");
const kernel = @import("kernel.zig");
const fma = @import("../../fma.zig");
const check = @import("../../_test/blas1/check.zig");

/// `kern` against a plain loop in `f64`; `T` is the element type, `conj`
/// conjugates `x` (`dotc`).
fn checkDot(comptime T: type, comptime conj: bool, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    for (check.sizes) |n| {
        var x: [check.capacity]T = undefined;
        var y: [check.capacity]T = undefined;
        check.fill(T, &x, n + 1);
        check.fill(T, &y, n + 7);
        for (check.zero_strides) |inc| {
            const W = check.Wide(F);
            var re: W = 0;
            var im: W = 0;
            for (0..n) |i| {
                const a = x[check.pos(n, inc[0], i)];
                const b = y[check.pos(n, inc[1], i)];
                if (comptime complex) {
                    const ar: W = a.re;
                    const ai: W = if (conj) -a.im else a.im;
                    re += ar * b.re - ai * b.im;
                    im += ai * b.re + ar * b.im;
                } else {
                    re += @as(W, a) * @as(W, b);
                }
            }
            const px: [*]F = @ptrCast(&x[check.pos(n, inc[0], 0)]);
            const got = kern(@intCast(n), px, inc[0], @as([*]F, @ptrCast(&y[check.pos(n, inc[1], 0)])), inc[1]);
            const R = @TypeOf(got);
            const tol = check.tolerance(if (complex) F else R, 4 * n);
            const ok = if (comptime complex)
                complexRef(F, re, im, got, tol)
            else
                check.expectCloseRef(R, re, got, tol);
            ok catch |e| return check.case(e, if (conj) "dotc" else "dot", T, n, inc);
        }
    }
}

fn complexRef(comptime F: type, re: check.Wide(F), im: check.Wide(F), got: anytype, tol: F) !void {
    try check.expectCloseRef(F, re, got.re, tol);
    try check.expectCloseRef(F, im, got.im, tol);
}

fn checkAll(vt: *const Vtable) !void {
    var failed: ?anyerror = null;
    check.keep(&failed, checkDot(f32, false, vt.sdot_k.?));
    check.keep(&failed, checkDot(f64, false, vt.ddot_k.?));
    check.keep(&failed, checkDot(f32, false, vt.dsdot_k.?));
    check.keep(&failed, checkDot(std.math.Complex(f32), false, vt.cdotu_k.?));
    check.keep(&failed, checkDot(std.math.Complex(f32), true, vt.cdotc_k.?));
    check.keep(&failed, checkDot(std.math.Complex(f64), false, vt.zdotu_k.?));
    check.keep(&failed, checkDot(std.math.Complex(f64), true, vt.zdotc_k.?));
    if (failed) |e| return e;
}

test "portable dot kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkDot(XDouble, false, vt.qdot_k.?);
    try checkDot(std.math.Complex(XDouble), false, vt.xdotu_k.?);
    try checkDot(std.math.Complex(XDouble), true, vt.xdotc_k.?);
}

test "portable dot kernels at the extended type match a plain loop" {
    const E = check.Ext;
    try checkDot(E, false, &Dot(E, E, E).k);
    try checkDot(std.math.Complex(E), false, &Zdot(E, false).k);
    try checkDot(std.math.Complex(E), true, &Zdot(E, true).k);
}

test "every linked core's dot kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}

/// `sdot` sums in `f64`: within half an ulp of the exact result, and the
/// `f64` sum's error, which an `f32` sum, in any lanes, exceeds by far;
/// `dsdot` is that sum unrounded.
fn checkSdotWide(vt: *const Vtable) !void {
    const al = std.testing.allocator;
    const long: usize = if (check.rig()) 1 << 12 else 1 << 14;
    const x = try al.alloc(f32, 3 * long);
    defer al.free(x);
    const y = try al.alloc(f32, 3 * long);
    defer al.free(y);
    for ([_]bool{ false, true }) |positive| for ([_]usize{ 1000, long }) |n| for ([_]usize{ 1, 3 }) |inc| {
        check.fill(f32, x, n);
        check.fill(f32, y, n + 1);
        if (positive) for (x, y) |*a, *b| {
            a.* = 1 + @abs(a.*) / 2;
            b.* = 1 + @abs(b.*) / 2;
        };
        var want: f128 = 0;
        var mag: f128 = 0;
        for (0..n) |i| {
            const p = @as(f128, x[i * inc]) * y[i * inc];
            want += p;
            mag += @abs(p);
        }
        const got = vt.sdot_k.?(@intCast(n), x.ptr, @intCast(inc), y.ptr, @intCast(inc));
        // `dsdot` is the same sum, unrounded, so a threaded `sdot`'s chunks
        // are its kernel's (`threading.md`)
        const wide = vt.dsdot_k.?(@intCast(n), x.ptr, @intCast(inc), y.ptr, @intCast(inc));
        if (@as(u32, @bitCast(got)) != @as(u32, @bitCast(@as(f32, @floatCast(wide))))) {
            check.say("sdot {e}, dsdot {e}: ", .{ got, wide });
            return check.case(error.TestUnexpectedResult, "dsdot", f32, n, .{ inc, inc });
        }
        if (@abs(wide - want) > @as(f128, @floatFromInt(n)) * 0x1p-52 * mag) {
            check.say("dsdot {e}, want {e}: ", .{ wide, @as(f64, @floatCast(want)) });
            return check.case(error.TestUnexpectedResult, "dsdot", f32, n, .{ inc, inc });
        }
        const ulp = std.math.ldexp(@as(f128, 1), std.math.ilogb(@as(f32, @floatCast(want))) - 23);
        const tol = ulp / 2 + @as(f128, @floatFromInt(n)) * 0x1p-52 * mag;
        if (@abs(got - want) > tol) {
            check.say("got {e}, want {e} within {e}: ", .{ @as(f64, got), @as(f64, @floatCast(want)), @as(f64, @floatCast(tol)) });
            return check.case(error.TestUnexpectedResult, "sdot", f32, n, .{ inc, inc });
        }
    };
}

test "portable sdot: within half an ulp, summed in f64, dsdot its sum" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkSdotWide(&vt);
}

test "every linked core's sdot: within half an ulp, summed in f64, dsdot its sum" {
    try check.eachCore(checkSdotWide);
}

test "the C drivers around the reference's fallback blocks match a plain loop" {
    try checkDot(f64, false, kernel.entry(.ddot_k, 16, null));
    inline for (.{ .{ .cdotu_k, false, 16 }, .{ .cdotc_k, true, 16 }, .{ .zdotu_k, false, 8 }, .{ .zdotc_k, true, 8 } }) |r| {
        const F = kernel.form(r[0]).E;
        try checkDot(std.math.Complex(F), r[1], kernel.entry(r[0], r[2], null));
    }
}

/// The reduction of `ddot_microk_skylakex-2.c` (`sdot_microk_skylakex-2.c` for
/// `f32`) lane by lane: `fm` is the contraction, `wide` the 512-bit pass.
fn skylakexModel(comptime T: type, comptime wide: bool, comptime fuse: bool, n: usize, x: []const T, y: []const T) T {
    const l256 = 32 / @sizeOf(T);
    const l512 = 64 / @sizeOf(T);
    var acc: [4][l256]T = @splat(@splat(0));
    var i: usize = 0;
    if (wide) {
        var w: [4][l512]T = @splat(@splat(0));
        const step = 4 * l512;
        while (i < n & ~@as(usize, step - 1)) : (i += step) {
            for (0..4) |g| for (0..l512) |l| {
                w[g][l] = fma.madd(fuse, T, x[i + g * l512 + l], y[i + g * l512 + l], w[g][l]);
            };
        }
        for (0..4) |g| for (0..l256) |l| {
            acc[g][l] = w[g][l] + w[g][l + l256];
        };
    }
    while (i < n) : (i += 4 * l256) {
        for (0..4) |g| for (0..l256) |l| {
            acc[g][l] = fma.madd(fuse, T, x[i + g * l256 + l], y[i + g * l256 + l], acc[g][l]);
        };
    }
    var h: [l256 / 2]T = undefined;
    for (&h, 0..) |*e, l| {
        const a = ((acc[0][l] + acc[1][l]) + acc[2][l]) + acc[3][l];
        const b = ((acc[0][l + l256 / 2] + acc[1][l + l256 / 2]) + acc[2][l + l256 / 2]) + acc[3][l + l256 / 2];
        e.* = a + b;
    }
    return if (T == f64) h[0] + h[1] else (h[0] + h[1]) + (h[2] + h[3]);
}

test "x86_64's vector blocks sum in the reference's order, bit for bit" {
    if (comptime sysinfo.target.cpu.arch != .x86_64) return error.SkipZigTest;
    const vb = @import("x86_64/kernel.zig");
    inline for (.{.{ .ddot_k, f64, 16 }}) |r| {
        const T = r[1];
        var x: [1024]T = undefined;
        var y: [1024]T = undefined;
        check.fill(T, &x, 3);
        check.fill(T, &y, 5);
        inline for (.{ false, true }) |wide| inline for (.{ false, true }) |fuse| {
            const B = vb.VectorBlockWith(r[0], wide, fuse);
            var n: usize = r[2];
            while (n <= 1024) : (n += r[2]) {
                const want = skylakexModel(T, wide, fuse, n, &x, &y);
                const got = B.run(@intCast(n), &x, &y);
                const U = if (T == f64) u64 else u32;
                zigstd.testExpectEqual(@as(U, @bitCast(want)), @as(U, @bitCast(got))) catch |e| {
                    check.say("{s} wide={} fuse={} n={d}: ", .{ @tagName(r[0]), wide, fuse, n });
                    return e;
                };
            }
        };
    }
}

/// The extended sums written out plainly: accumulator by position, then folded.
const PlainExtended = struct {
    fn dot(comptime F: type, n: usize, x: []const F, ix: usize, y: []const F, iy: usize) F {
        var a: [4]F = @splat(0);
        for (0..n) |i| a[if (i < n / 4 * 4) i % 4 else 0] += x[i * ix] * y[i * iy];
        return (a[3] + a[1]) + (a[2] + a[0]);
    }

    fn zdot(comptime F: type, conj: bool, n: usize, x: []const F, ix: usize, y: []const F, iy: usize) [2]F {
        var a: [4]F = @splat(0);
        for (0..n) |i| for (0..2) |p| for (0..2) |q| {
            a[2 * p + q] += x[2 * i * ix + p] * y[2 * i * iy + q];
        };
        return if (conj) .{ a[0] + a[3], a[1] - a[2] } else .{ a[0] - a[3], a[1] + a[2] };
    }
};

test "the extended dot sums fold in dot.S's and zdot.S's order, bit for bit" {
    inline for (.{ f80, f128 }) |F| {
        var x: [check.capacity * 2]F = undefined;
        var y: [check.capacity * 2]F = undefined;
        for (check.sizes) |n| for ([_]usize{ 1, 3 }) |inc| {
            check.fill(F, &x, n + 11);
            check.fill(F, &y, n + 12);
            const len: BlasLong = @intCast(n);
            const li: BlasLong = @intCast(inc);
            try check.expectBits(F, PlainExtended.dot(F, n, &x, inc, &y, 1), extendedDot(F, len, &x, li, &y, 1), "dot", n, inc);
            inline for (.{ false, true }) |conj| {
                const want = PlainExtended.zdot(F, conj, n, &x, inc, &y, 2);
                const got = extendedZdot(F, conj, len, &x, li, &y, 2);
                try check.expectBits(F, want[0], got[0], "zdot re", n, inc);
                try check.expectBits(F, want[1], got[1], "zdot im", n, inc);
            }
        };
    }
}

test "the extended dot at f80 carries the full significand" {
    // 1 + 2^-63 squared is 1 + 2^-62 in f80 and 1 in f64
    var x = [_]f80{ 1 + std.math.ldexp(@as(f80, 1), -63), 0, 0, 0, 0 };
    try zigstd.testExpectEqual(1 + std.math.ldexp(@as(f80, 1), -62), extendedDot(f80, 5, &x, 1, &x, 1));
}

test "the extended dot sums give the values traced by hand from dot.S and zdot.S" {
    // 2^64 + 1 rounds to 2^64 in f80, so each value below depends on which
    // accumulator a term lands in and on the fold order. Traced once from
    // kernel/x86_64/: dot `faddp %st, %st(k)` per i % 4, tail into st(0),
    // folded (a3 + a1) + (a2 + a0); zdot a0..a3 = xr*yr, xr*yi, xi*yr, xi*yi,
    // re = a0 - a3 (`fsubp` is fsubr in AT&T), im = a1 + a2, conj a0 + a3 and
    // a1 - a2.
    const F = f80;
    const t = std.math.ldexp(@as(F, 1), 64);
    var ones: [8]F = @splat(1);
    // dot: a0 = 1 + 1 + 1 (tail), a1 = t, a2 = 1, a3 = -t: (0) + (1 + 3)
    var dx = [_]F{ 1, t, 1, -t, 1, 1 };
    try zigstd.testExpectEqual(@as(F, 4), extendedDot(F, 6, &dx, 1, &ones, 1));
    // zdot of (t, t), (1, -1) with (1, 1), (1, 1): a0 = t, a1 = t, a2 = t - 1, a3 = t - 1
    var zx = [_]F{ t, t, 1, -1 };
    const u = extendedZdot(F, false, 2, &zx, 1, &ones, 1);
    try zigstd.testExpectEqual(@as(F, 1), u[0]);
    try zigstd.testExpectEqual(2 * t, u[1]); // 2t - 1 rounds to even
    const c = extendedZdot(F, true, 2, &zx, 1, &ones, 1);
    try zigstd.testExpectEqual(2 * t, c[0]);
    try zigstd.testExpectEqual(@as(F, 1), c[1]);
}
