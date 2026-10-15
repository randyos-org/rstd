//! `DOT`: `x . y` (`dotu`/`dotc` for complex, `dsdot` summing `f32` in `f64`).
//! The vtable wrappers, the portable kernels (`*_generic`), and tests of every
//! core's kernels against a plain loop.
//!
//! Argument names are copied from the OpenBLAS kernel sources (see the
//! `CNAME(...)` definitions under `kernel/<arch>/*.c`). Several are named
//! `dummy*` in the reference: kernel-ABI padding slots this shape does not use.

const blas = @import("../../__root__.zig");
const sysinfo = @import("builtin");
const zigstd = @import("../../../../zigstd.zig");
const kernel = @import("kernel.zig");
const uk = @import("../uk_block.zig");
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
const check = @import("../check.zig");
const x87 = @import("../x87.zig");

/// `x . y` after `kernel/arm/dot.c`, summed sequentially in `Acc` (`f64` for
/// `sdot` too, so it is at least as accurate as the reference's `float` or
/// `double` products) and returned as `R` (`dsdot` returns the `f64` sum).
/// Not bit-identical to `generic/dot.c`'s unrolled sums; see the readme.
/// Wider than `f64`, it is x86_64 `dot.S`'s x87 sum instead (`x87.zig`).
fn Dot(comptime F: type, comptime Acc: type, comptime R: type) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) callconv(types.kernelCc(F)) R {
            if (comptime types.isExtended(F)) return x87.dot(F, n, x, inc_x, y, inc_y);
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
                const z = x87.zdot(F, conj, n, x, inc_x, y, inc_y);
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

/// `kern` against a plain loop in `f64`; `T` is the element type, `conj`
/// conjugates `x` (`dotc`). `narrow` allows a `dsdot` the error of summing
/// `f32` products: the reference's blocks do, 32 at a time.
fn checkDot(comptime T: type, comptime conj: bool, kern: anytype) !void {
    try checkDotTol(T, conj, kern, false);
}

fn checkDotTol(comptime T: type, comptime conj: bool, kern: anytype, comptime narrow: bool) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    for (check.sizes) |n| {
        var x: [check.capacity]T = undefined;
        var y: [check.capacity]T = undefined;
        check.fill(T, &x, n + 1);
        check.fill(T, &y, n + 7);
        for (check.strides) |inc| {
            const W = check.Wide(F);
            var re: W = 0;
            var im: W = 0;
            for (0..n) |i| {
                const a = x[i * inc[0]];
                const b = y[i * inc[1]];
                if (comptime complex) {
                    const ar: W = a.re;
                    const ai: W = if (conj) -a.im else a.im;
                    re += ar * b.re - ai * b.im;
                    im += ai * b.re + ar * b.im;
                } else {
                    re += @as(W, a) * @as(W, b);
                }
            }
            const got = kern(@intCast(n), @as([*]F, @ptrCast(&x)), @intCast(inc[0]), @as([*]F, @ptrCast(&y)), @intCast(inc[1]));
            const R = @TypeOf(got);
            const tol = check.tolerance(if (complex or narrow) F else R, 4 * n);
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
    try checkDot(f32, false, vt.sdot_k.?);
    try checkDot(f64, false, vt.ddot_k.?);
    try checkDotTol(f32, false, vt.dsdot_k.?, true);
    try checkDot(std.math.Complex(f32), false, vt.cdotu_k.?);
    try checkDot(std.math.Complex(f32), true, vt.cdotc_k.?);
    try checkDot(std.math.Complex(f64), false, vt.zdotu_k.?);
    try checkDot(std.math.Complex(f64), true, vt.zdotc_k.?);
}

test "portable dot kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkDot(f32, false, vt.dsdot_k.?);
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

test "the C drivers around the reference's fallback blocks match a plain loop" {
    inline for (.{ .{ "sdot", 32 }, .{ "ddot", 16 } }) |r| {
        try checkDot(uk.Float(r[0]), false, kernel.entry(r[0], r[1], null));
    }
    try checkDotTol(f32, false, kernel.entry("dsdot", 32, null), true);
    inline for (.{ .{ "cdotu", false, 16 }, .{ "cdotc", true, 16 }, .{ "zdotu", false, 8 }, .{ "zdotc", true, 8 } }) |r| {
        const F = uk.Float(r[0]);
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
                w[g][l] = kernel.fmadd(fuse, T, x[i + g * l512 + l], y[i + g * l512 + l], w[g][l]);
            };
        }
        for (0..4) |g| for (0..l256) |l| {
            acc[g][l] = w[g][l] + w[g][l + l256];
        };
    }
    while (i < n) : (i += 4 * l256) {
        for (0..4) |g| for (0..l256) |l| {
            acc[g][l] = kernel.fmadd(fuse, T, x[i + g * l256 + l], y[i + g * l256 + l], acc[g][l]);
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
    if (comptime sysinfo.cpu.arch != .x86_64) return error.SkipZigTest;
    const vb = @import("x86_64/kernel.zig");
    inline for (.{ .{ "ddot", f64, 16 }, .{ "sdot", f32, 32 } }) |r| {
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
                    check.say("{s} wide={} fuse={} n={d}: ", .{ r[0], wide, fuse, n });
                    return e;
                };
            }
        };
    }
}
