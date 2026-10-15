//! `ROT`: the plane rotation `x := c*x + s*y`, `y := c*y - s*x`. The vtable
//! wrappers, the portable kernels (`*_generic`), and tests of every core's
//! kernels against a plain loop.
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

pub fn sbrot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, c: f32, s: f32) i32 {
    return blas.kernel(.sbrot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn srot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, c: f32, s: f32) i32 {
    return blas.kernel(.srot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn drot_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, c: f64, s: f64) i32 {
    return blas.kernel(.drot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn qrot_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, c: XDouble, s: XDouble) i32 {
    return blas.kernel(.qrot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn csrot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, c: f32, s: f32) i32 {
    return blas.kernel(.csrot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn zdrot_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, c: f64, s: f64) i32 {
    return blas.kernel(.zdrot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn xqrot_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, c: XDouble, s: XDouble) i32 {
    return blas.kernel(.xqrot_k)(n, x, inc_x, y, inc_y, c, s);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../check.zig");

/// `kernel/arm/rot.c` and `zrot.c`: the plane rotation `x := c*x + s*y`,
/// `y := c*y - s*x`; an element is `w` floats, each rotated alike (`csrot`
/// and `zdrot` take a real `c` and `s`).
fn Rot(comptime F: type, comptime w: comptime_int) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, c: F, s: F) callconv(.c) i32 {
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                inline for (0..w) |j| {
                    const px = util.at(F, x, i * inc_x * w + j);
                    const py = util.at(F, y, i * inc_y * w + j);
                    const t = c * px.* + s * py.*;
                    py.* = c * py.* - s * px.*;
                    px.* = t;
                }
            }
            return 0;
        }
    };
}

pub const srot_generic = &Rot(f32, 1).k;
pub const drot_generic = &Rot(f64, 1).k;
pub const qrot_generic = &Rot(XDouble, 1).k;
pub const csrot_generic = &Rot(f32, 2).k;
pub const zdrot_generic = &Rot(f64, 2).k;
pub const xqrot_generic = &Rot(XDouble, 2).k;

/// Points every rot slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.srot_k = srot_generic;
    vt.drot_k = drot_generic;
    vt.qrot_k = qrot_generic;
    vt.csrot_k = csrot_generic;
    vt.zdrot_k = zdrot_generic;
    vt.xqrot_k = xqrot_generic;
}

// ---- tests -----------------------------------------------------------------

/// `kern` against a plain loop; `T` is the element type, real or complex.
fn checkRot(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    const c: F = 0.8;
    const s: F = -0.6;
    for (check.sizes) |n| {
        for (check.strides) |inc| {
            var x: [check.capacity]T = undefined;
            var y: [check.capacity]T = undefined;
            check.fill(T, &x, n + 1);
            check.fill(T, &y, n + 7);
            var want_x = x;
            var want_y = y;
            for (0..n) |i| {
                const a = x[i * inc[0]];
                const b = y[i * inc[1]];
                if (comptime check.isComplex(T)) {
                    want_x[i * inc[0]] = .init(c * a.re + s * b.re, c * a.im + s * b.im);
                    want_y[i * inc[1]] = .init(c * b.re - s * a.re, c * b.im - s * a.im);
                } else {
                    want_x[i * inc[0]] = c * a + s * b;
                    want_y[i * inc[1]] = c * b - s * a;
                }
            }
            _ = kern(@intCast(n), @as([*]F, @ptrCast(&x)), @intCast(inc[0]), @as([*]F, @ptrCast(&y)), @intCast(inc[1]), c, s);
            const tol = check.tolerance(F, 2);
            check.expectCloseSlices(T, &want_x, &x, tol) catch |e| return check.case(e, "rot x", T, n, inc);
            check.expectCloseSlices(T, &want_y, &y, tol) catch |e| return check.case(e, "rot y", T, n, inc);
        }
    }
}

fn checkAll(vt: *const Vtable) !void {
    try checkRot(f32, vt.srot_k.?);
    try checkRot(f64, vt.drot_k.?);
    try checkRot(std.math.Complex(f32), vt.csrot_k.?);
    try checkRot(std.math.Complex(f64), vt.zdrot_k.?);
}

test "portable rot kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkRot(XDouble, vt.qrot_k.?);
    try checkRot(std.math.Complex(XDouble), vt.xqrot_k.?);
}

test "every linked core's rot kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
