//! `AXPBY`: `y := alpha*x + beta*y`. The vtable wrappers, the portable
//! kernels (`*_generic`), and tests of every core's kernels against a plain
//! loop.
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

pub fn saxpby_k(n: BlasLong, alpha: f32, x: [*]f32, inc_x: BlasLong, beta: f32, y: [*]f32, inc_y: BlasLong) i32 {
    return blas.kernel(.saxpby_k)(n, alpha, x, inc_x, beta, y, inc_y);
}
pub fn daxpby_k(n: BlasLong, alpha: f64, x: [*]f64, inc_x: BlasLong, beta: f64, y: [*]f64, inc_y: BlasLong) i32 {
    return blas.kernel(.daxpby_k)(n, alpha, x, inc_x, beta, y, inc_y);
}
pub fn caxpby_k(n: BlasLong, alpha_r: f32, alpha_i: f32, x: [*]f32, inc_x: BlasLong, beta_r: f32, beta_i: f32, y: [*]f32, inc_y: BlasLong) i32 {
    return blas.kernel(.caxpby_k)(n, alpha_r, alpha_i, x, inc_x, beta_r, beta_i, y, inc_y);
}
pub fn zaxpby_k(n: BlasLong, alpha_r: f64, alpha_i: f64, x: [*]f64, inc_x: BlasLong, beta_r: f64, beta_i: f64, y: [*]f64, inc_y: BlasLong) i32 {
    return blas.kernel(.zaxpby_k)(n, alpha_r, alpha_i, x, inc_x, beta_r, beta_i, y, inc_y);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../check.zig");

/// `kernel/arm/axpby.c` and `zaxpby.c`: `y := alpha*x + beta*y`. A zero
/// alpha does not read `x`, and a zero beta does not read `y`.
fn Axpby(comptime F: type) type {
    return struct {
        fn real(n: BlasLong, alpha: F, x: [*]F, inc_x: BlasLong, beta: F, y: [*]F, inc_y: BlasLong) callconv(.c) i32 {
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const py = util.at(F, y, i * inc_y);
                const ax: F = if (alpha == 0) 0 else alpha * util.at(F, x, i * inc_x).*;
                py.* = if (beta == 0) ax else ax + beta * py.*;
            }
            return 0;
        }
        fn complex(n: BlasLong, alpha_r: F, alpha_i: F, x: [*]F, inc_x: BlasLong, beta_r: F, beta_i: F, y: [*]F, inc_y: BlasLong) callconv(.c) i32 {
            const no_alpha = alpha_r == 0 and alpha_i == 0;
            const no_beta = beta_r == 0 and beta_i == 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const yr = util.at(F, y, i * inc_y * 2);
                const yi = util.at(F, y, i * inc_y * 2 + 1);
                var re: F = 0;
                var im: F = 0;
                if (!no_alpha) {
                    const xr = util.at(F, x, i * inc_x * 2).*;
                    const xi = util.at(F, x, i * inc_x * 2 + 1).*;
                    re = alpha_r * xr - alpha_i * xi;
                    im = alpha_r * xi + alpha_i * xr;
                }
                if (!no_beta) {
                    const r = beta_r * yr.* - beta_i * yi.*;
                    im += beta_r * yi.* + beta_i * yr.*;
                    re += r;
                }
                yr.* = re;
                yi.* = im;
            }
            return 0;
        }
    };
}

pub const saxpby_generic = &Axpby(f32).real;
pub const daxpby_generic = &Axpby(f64).real;
pub const caxpby_generic = &Axpby(f32).complex;
pub const zaxpby_generic = &Axpby(f64).complex;

/// Points every axpby slot at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.saxpby_k = saxpby_generic;
    vt.daxpby_k = daxpby_generic;
    vt.caxpby_k = caxpby_generic;
    vt.zaxpby_k = zaxpby_generic;
}

// ---- tests -----------------------------------------------------------------

/// `kern` against a plain loop for each of alpha and beta zero or not; `T` is
/// the element type, real or complex.
fn checkAxpby(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    const values = if (complex) [_]T{ T.init(1.5, -0.75), T.init(0, 0) } else [_]T{ 1.5, 0 };
    const betas = if (complex) [_]T{ T.init(-0.5, 0.25), T.init(0, 0) } else [_]T{ -0.5, 0 };
    for (values) |alpha| for (betas) |beta| {
        for (check.sizes) |n| {
            for (check.strides) |inc| {
                var x: [check.capacity]T = undefined;
                var y: [check.capacity]T = undefined;
                check.fill(T, &x, n + 1);
                check.fill(T, &y, n + 7);
                var want = y;
                for (0..n) |i| {
                    const a = x[i * inc[0]];
                    const b = y[i * inc[1]];
                    want[i * inc[1]] = if (complex)
                        a.mul(alpha).add(b.mul(beta))
                    else
                        alpha * a + beta * b;
                }
                const px: [*]F = @ptrCast(&x);
                const py: [*]F = @ptrCast(&y);
                if (comptime complex) {
                    _ = kern(@intCast(n), alpha.re, alpha.im, px, @intCast(inc[0]), beta.re, beta.im, py, @intCast(inc[1]));
                } else {
                    _ = kern(@intCast(n), alpha, px, @intCast(inc[0]), beta, py, @intCast(inc[1]));
                }
                check.expectCloseSlices(T, &want, &y, check.tolerance(F, 4)) catch |e| {
                    check.say("alpha {any}, beta {any}: ", .{ alpha, beta });
                    return check.case(e, "axpby", T, n, inc);
                };
            }
        }
    };
}

fn checkAll(vt: *const Vtable) !void {
    try checkAxpby(f32, vt.saxpby_k.?);
    try checkAxpby(f64, vt.daxpby_k.?);
    try checkAxpby(std.math.Complex(f32), vt.caxpby_k.?);
    try checkAxpby(std.math.Complex(f64), vt.zaxpby_k.?);
}

test "portable axpby kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
}

test "every linked core's axpby kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
