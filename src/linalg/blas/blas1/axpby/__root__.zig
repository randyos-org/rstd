//! `AXPBY`: `y := alpha*x + beta*y`. The portable kernels (`*_generic`)
//! and tests of every core's kernels against a plain loop.

const types = @import("../../types.zig");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

const BlasLong = types.BlasLong;

// ---- portable kernels ------------------------------------------------------

/// `kernel/arm/axpby.c` and `zaxpby.c`: `y := alpha*x + beta*y`. A zero
/// alpha does not read `x`, and a zero beta does not read `y`.
fn Axpby(comptime F: type) type {
    return struct {
        fn real(n: BlasLong, alpha: F, x: [*]F, inc_x: BlasLong, beta: F, y: [*]F, inc_y: BlasLong) callconv(types.kernelCc(F)) i32 {
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const py = util.at(F, y, i * inc_y);
                const ax: F = if (alpha == 0) 0 else alpha * util.at(F, x, i * inc_x).*;
                py.* = if (beta == 0) ax else ax + beta * py.*;
            }
            return 0;
        }
        fn complex(n: BlasLong, alpha_r: F, alpha_i: F, x: [*]F, inc_x: BlasLong, beta_r: F, beta_i: F, y: [*]F, inc_y: BlasLong) callconv(types.kernelCc(F)) i32 {
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

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    return struct {
        pub const axpby_k = if (types.isComplex(T)) &Axpby(F).complex else &Axpby(F).real;
    };
}

/// Points every axpby slot at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.saxpby_k = saxpby_generic;
    vt.daxpby_k = daxpby_generic;
    vt.caxpby_k = caxpby_generic;
    vt.zaxpby_k = zaxpby_generic;
}

// ---- tests -----------------------------------------------------------------

const std = @import("std");
const check = @import("../../_test/check.zig");

/// `kern` against a plain loop for each of alpha and beta zero or not; `T` is
/// the element type, real or complex.
fn checkAxpby(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    const values = if (complex) [_]T{ T.init(1.5, -0.75), T.init(0, 0) } else [_]T{ 1.5, 0 };
    const betas = if (complex) [_]T{ T.init(-0.5, 0.25), T.init(0, 0) } else [_]T{ -0.5, 0 };
    for (values) |alpha| for (betas) |beta| {
        for (check.sizes) |n| {
            for (check.zero_strides) |inc| {
                var x: [check.capacity]T = undefined;
                var y: [check.capacity]T = undefined;
                check.fill(T, &x, n + 1);
                check.fill(T, &y, n + 7);
                var want = y;
                for (0..n) |i| {
                    const iy = check.pos(n, inc[1], i);
                    const a = x[check.pos(n, inc[0], i)];
                    const b = want[iy];
                    want[iy] = if (complex)
                        a.mul(alpha).add(b.mul(beta))
                    else
                        alpha * a + beta * b;
                }
                const px: [*]F = @ptrCast(&x[check.pos(n, inc[0], 0)]);
                const py: [*]F = @ptrCast(&y[check.pos(n, inc[1], 0)]);
                if (comptime complex) {
                    _ = kern(@intCast(n), alpha.re, alpha.im, px, inc[0], beta.re, beta.im, py, inc[1]);
                } else {
                    _ = kern(@intCast(n), alpha, px, inc[0], beta, py, inc[1]);
                }
                check.expectStrided(T, &want, &y, check.tolerance(F, 4), n, inc) catch |e| {
                    check.say("alpha {any}, beta {any}: ", .{ alpha, beta });
                    return check.case(e, "axpby", T, n, inc);
                };
            }
        }
    };
}

fn checkAll(vt: *const Vtable) !void {
    var failed: ?anyerror = null;
    check.keep(&failed, checkAxpby(f32, vt.saxpby_k.?));
    check.keep(&failed, checkAxpby(f64, vt.daxpby_k.?));
    check.keep(&failed, checkAxpby(std.math.Complex(f32), vt.caxpby_k.?));
    check.keep(&failed, checkAxpby(std.math.Complex(f64), vt.zaxpby_k.?));
    if (failed) |e| return e;
}

test "portable axpby kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
}

test "portable axpby kernels at the extended type match a plain loop" {
    try checkAxpby(check.Ext, &Axpby(check.Ext).real);
    try checkAxpby(std.math.Complex(check.Ext), &Axpby(check.Ext).complex);
}

test "every linked core's axpby kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
