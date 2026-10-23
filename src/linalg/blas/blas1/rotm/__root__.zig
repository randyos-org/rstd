//! `ROTM`: the modified Givens rotation `param` names. The vtable wrappers,
//! the portable kernels (`*_generic`), and tests of every core's kernels
//! against a plain loop.
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

pub fn sbrotm_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, param: [*]f32) i32 {
    return blas.kernel(.sbrotm_k)(n, x, inc_x, y, inc_y, param);
}
pub fn srotm_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, param: [*]f32) i32 {
    return blas.kernel(.srotm_k)(n, x, inc_x, y, inc_y, param);
}
pub fn drotm_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, param: [*]f64) i32 {
    return blas.kernel(.drotm_k)(n, x, inc_x, y, inc_y, param);
}
pub fn qrotm_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, param: [*]XDouble) i32 {
    return blas.kernel(.qrotm_k)(n, x, inc_x, y, inc_y, param);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../check.zig");

/// The matrix `param` names: `param[0]` is the flag, -2 the identity, -1 the
/// full `[h11 h12; h21 h22]` (`param[1..5]` = h11, h21, h12, h22), 0
/// `[1 h12; h21 1]` and 1 `[h11 1; -1 h22]`.
fn matrix(comptime F: type, param: [*]const F) ?[4]F {
    const flag = param[0];
    if (flag == -2) return null;
    if (flag < 0) return .{ param[1], param[3], param[2], param[4] };
    if (flag == 0) return .{ 1, param[3], param[2], 1 };
    return .{ param[1], 1, -1, param[4] };
}

/// `kernel/generic/rotm.c`: `x := h11*x + h12*y`, `y := h21*x + h22*y`. A
/// negative stride starts from the far end, as reference BLAS does.
fn Rotm(comptime F: type) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, param: [*]F) callconv(types.kernelCc(F)) i32 {
            if (n <= 0) return 0;
            const h = matrix(F, param) orelse return 0;
            var ix: BlasLong = if (inc_x < 0) (1 - n) * inc_x else 0;
            var iy: BlasLong = if (inc_y < 0) (1 - n) * inc_y else 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const px = util.at(F, x, ix);
                const py = util.at(F, y, iy);
                const a = px.*;
                const b = py.*;
                px.* = a * h[0] + b * h[1];
                py.* = a * h[2] + b * h[3];
                ix += inc_x;
                iy += inc_y;
            }
            return 0;
        }
    };
}

pub const srotm_generic = &Rotm(f32).k;
pub const drotm_generic = &Rotm(f64).k;
pub const qrotm_generic = &Rotm(XDouble).k;

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    return struct {
        pub const rotm_k = &Rotm(F).k;
    };
}

/// Points every rotm slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.srotm_k = srotm_generic;
    vt.drotm_k = drotm_generic;
    vt.qrotm_k = qrotm_generic;
}

// ---- tests -----------------------------------------------------------------

/// `kern` against a plain loop, for every flag.
fn checkRotm(comptime F: type, kern: anytype) !void {
    for ([_]F{ -2, -1, 0, 1 }) |flag| {
        var param = [5]F{ flag, 0.75, -0.5, 0.25, 1.25 };
        const h = matrix(F, &param);
        for (check.sizes) |n| {
            for (check.zero_strides) |inc| {
                var x: [check.capacity]F = undefined;
                var y: [check.capacity]F = undefined;
                check.fill(F, &x, n + 1);
                check.fill(F, &y, n + 7);
                var want_x = x;
                var want_y = y;
                if (h) |m| for (0..n) |i| {
                    const ix = check.pos(n, inc[0], i);
                    const iy = check.pos(n, inc[1], i);
                    const a = want_x[ix];
                    const b = want_y[iy];
                    want_x[ix] = a * m[0] + b * m[1];
                    want_y[iy] = a * m[2] + b * m[3];
                };
                // the near end, as the interface passes it
                _ = kern(@intCast(n), &x, inc[0], &y, inc[1], &param);
                const tol = check.tolerance(F, 2);
                check.expectStrided(F, &want_x, &x, tol, n, inc) catch |e| return check.case(e, "rotm x", F, n, inc);
                check.expectStrided(F, &want_y, &y, tol, n, inc) catch |e| return check.case(e, "rotm y", F, n, inc);
            }
        }
    }
}

fn checkAll(vt: *const Vtable) !void {
    var failed: ?anyerror = null;
    check.keep(&failed, checkRotm(f32, vt.srotm_k.?));
    check.keep(&failed, checkRotm(f64, vt.drotm_k.?));
    if (failed) |e| return e;
}

test "portable rotm kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkRotm(XDouble, vt.qrotm_k.?);
}

test "portable rotm kernels at the extended type match a plain loop" {
    try checkRotm(check.Ext, &Rotm(check.Ext).k);
}

test "every linked core's rotm kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
