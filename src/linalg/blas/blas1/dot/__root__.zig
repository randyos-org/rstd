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
const BlasInt = types.BlasInt;
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

/// `kernel/arm/dot.c`: `x . y`, summed in `Acc` (`double` in the reference,
/// for `sdot` too) and returned as `R` (`dsdot` returns the `f64` sum).
fn Dot(comptime F: type, comptime Acc: type, comptime R: type) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) callconv(.c) R {
            var sum: Acc = 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) sum += @as(Acc, util.at(F, x, i * inc_x).*) * @as(Acc, util.at(F, y, i * inc_y).*);
            return @floatCast(sum);
        }
    };
}

/// `kernel/arm/zdot.c`: `x . y`, or `conj(x) . y` (`dotc`), summed in `F`.
fn Zdot(comptime F: type, comptime conj: bool) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) callconv(.c) types.CComplex(F) {
            var re: F = 0;
            var im: F = 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const xr = util.at(F, x, i * inc_x * 2).*;
                const xi = util.at(F, x, i * inc_x * 2 + 1).*;
                const yr = util.at(F, y, i * inc_y * 2).*;
                const yi = util.at(F, y, i * inc_y * 2 + 1).*;
                if (conj) {
                    re += xr * yr + xi * yi;
                    im -= xi * yr - xr * yi;
                } else {
                    re += xr * yr - xi * yi;
                    im += xi * yr + xr * yi;
                }
            }
            return .{ .re = re, .im = im };
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

/// Points every dot slot but bfloat16's at the portable kernels.
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
/// conjugates `x` (`dotc`).
fn checkDot(comptime T: type, comptime conj: bool, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    for (check.sizes) |n| {
        var x: [check.capacity]T = undefined;
        var y: [check.capacity]T = undefined;
        check.fill(T, &x, n + 1);
        check.fill(T, &y, n + 7);
        for (check.strides) |inc| {
            var re: f64 = 0;
            var im: f64 = 0;
            for (0..n) |i| {
                const a = x[i * inc[0]];
                const b = y[i * inc[1]];
                if (comptime complex) {
                    const ar: f64 = a.re;
                    const ai: f64 = if (conj) -a.im else a.im;
                    re += ar * b.re - ai * b.im;
                    im += ai * b.re + ar * b.im;
                } else {
                    re += @as(f64, a) * @as(f64, b);
                }
            }
            const got = kern(@intCast(n), @as([*]F, @ptrCast(&x)), @intCast(inc[0]), @as([*]F, @ptrCast(&y)), @intCast(inc[1]));
            const R = @TypeOf(got);
            const tol = check.tolerance(if (complex) F else R, 4 * n);
            const ok = if (comptime complex)
                check.expectClose(std.math.Complex(F), .init(@floatCast(re), @floatCast(im)), got.toComplex(), tol)
            else
                check.expectClose(R, @floatCast(re), got, tol);
            ok catch |e| return check.case(e, if (conj) "dotc" else "dot", T, n, inc);
        }
    }
}

fn checkAll(vt: *const Vtable) !void {
    try checkDot(f32, false, vt.sdot_k.?);
    try checkDot(f64, false, vt.ddot_k.?);
    try checkDot(f32, false, vt.dsdot_k.?);
    try checkDot(std.math.Complex(f32), false, vt.cdotu_k.?);
    try checkDot(std.math.Complex(f32), true, vt.cdotc_k.?);
    try checkDot(std.math.Complex(f64), false, vt.zdotu_k.?);
    try checkDot(std.math.Complex(f64), true, vt.zdotc_k.?);
}

test "portable dot kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkDot(XDouble, false, vt.qdot_k.?);
    try checkDot(std.math.Complex(XDouble), false, vt.xdotu_k.?);
    try checkDot(std.math.Complex(XDouble), true, vt.xdotc_k.?);
}

test "every linked core's dot kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
