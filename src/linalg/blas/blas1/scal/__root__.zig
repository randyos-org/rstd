//! `SCAL`: `x := alpha*x`. The vtable wrappers, the portable kernels
//! (`*_generic`), and tests of every core's kernels against a plain loop.
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

pub fn bscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da: BF16, x: [*]BF16, inc_x: BlasLong, y: [*]BF16, inc_y: BlasLong, dummy: [*]BF16, dummy2: BlasLong) i32 {
    return blas.kernel(.bscal_k)(n, dummy0, dummy1, da, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn sbscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas.kernel(.sbscal_k)(n, dummy0, dummy1, da, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn sscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas.kernel(.sscal_k)(n, dummy0, dummy1, da, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn dscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, dummy: [*]f64, dummy2: BlasLong) i32 {
    return blas.kernel(.dscal_k)(n, dummy0, dummy1, da, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn qscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, dummy: [*]XDouble, dummy2: BlasLong) i32 {
    return blas.kernel(.qscal_k)(n, dummy0, dummy1, da, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn cscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da_r: f32, da_i: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas.kernel(.cscal_k)(n, dummy0, dummy1, da_r, da_i, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn zscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da_r: f64, da_i: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, dummy: [*]f64, dummy2: BlasLong) i32 {
    return blas.kernel(.zscal_k)(n, dummy0, dummy1, da_r, da_i, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn xscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da_r: XDouble, da_i: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, dummy: [*]XDouble, dummy2: BlasLong) i32 {
    return blas.kernel(.xscal_k)(n, dummy0, dummy1, da_r, da_i, x, inc_x, y, inc_y, dummy, dummy2);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

/// `kernel/arm/scal.c` and `zscal.c`: `x := alpha*x`. `flag` is the last
/// argument (`dummy2`): with 0, a zero alpha zeroes `x`; with 1, which the
/// interface passes, a zero alpha turns NaN and infinity into NaN instead.
fn Scal(comptime F: type) type {
    return struct {
        fn real(n: BlasLong, _: BlasLong, _: BlasLong, da: F, x: [*]F, inc_x: BlasLong, _: [*]F, _: BlasLong, _: [*]F, flag: BlasLong) callconv(types.kernelCc(F)) i32 {
            if (n <= 0 or inc_x <= 0) return 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const v = util.at(F, x, i * inc_x);
                v.* = if (da != 0) da * v.* else if (flag == 0 or std.math.isFinite(v.*)) 0 else std.math.nan(F);
            }
            return 0;
        }
        fn complex(n: BlasLong, _: BlasLong, _: BlasLong, da_r: F, da_i: F, x: [*]F, inc_x: BlasLong, _: [*]F, _: BlasLong, _: [*]F, flag: BlasLong) callconv(types.kernelCc(F)) i32 {
            if (n <= 0 or inc_x <= 0) return 0;
            const zero = flag == 0 and da_r == 0 and da_i == 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const re = util.at(F, x, i * inc_x * 2);
                const im = util.at(F, x, i * inc_x * 2 + 1);
                if (zero) {
                    re.* = 0;
                    im.* = 0;
                } else {
                    const t = da_r * re.* - da_i * im.*;
                    im.* = da_r * im.* + da_i * re.*;
                    re.* = t;
                }
            }
            return 0;
        }
    };
}

pub const sscal_generic = &Scal(f32).real;
pub const dscal_generic = &Scal(f64).real;
pub const qscal_generic = &Scal(XDouble).real;
pub const cscal_generic = &Scal(f32).complex;
pub const zscal_generic = &Scal(f64).complex;
pub const xscal_generic = &Scal(XDouble).complex;

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    return struct {
        pub const scal_k = if (types.isComplex(T)) &Scal(F).complex else &Scal(F).real;
    };
}

/// Points every scal slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.sscal_k = sscal_generic;
    vt.dscal_k = dscal_generic;
    vt.qscal_k = qscal_generic;
    vt.cscal_k = cscal_generic;
    vt.zscal_k = zscal_generic;
    vt.xscal_k = xscal_generic;
}

// ---- tests -----------------------------------------------------------------

const check = @import("../../_test/blas1/check.zig");
const zigstd = @import("../../../../zigstd.zig");

/// `kern` against a plain loop, as the interface calls it (flag 1), for a
/// general and a zero alpha; `T` is the element type, real or complex.
fn checkScal(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    const alphas = if (comptime complex) [_]T{ T.init(1.5, -0.75), T.init(0, 0) } else [_]T{ 1.5, 0 };
    for (alphas) |alpha| {
        for (check.sizes) |n| {
            for ([_]usize{ 1, 2, 3 }) |inc| {
                var x: [check.capacity]T = undefined;
                check.fill(T, &x, n + 1);
                var want = x;
                for (0..n) |i| {
                    const v = &want[i * inc];
                    if (comptime complex) {
                        v.* = .init(alpha.re * v.re - alpha.im * v.im, alpha.re * v.im + alpha.im * v.re);
                    } else {
                        v.* = alpha * v.*;
                    }
                }
                const px: [*]F = @ptrCast(&x);
                if (comptime complex) {
                    _ = kern(@intCast(n), 0, 0, alpha.re, alpha.im, px, @intCast(inc), px, 0, px, 1);
                } else {
                    _ = kern(@intCast(n), 0, 0, alpha, px, @intCast(inc), px, 0, px, 1);
                }
                check.expectCloseSlices(T, &want, &x, check.tolerance(F, 1)) catch |e| return check.case(e, "scal", T, n, .{ inc, 0 });
            }
        }
    }
}

/// A zero alpha on NaN and infinity: flag 0, as the level-2 drivers pass it,
/// stores zeros; flag 1, the interface's, multiplies, giving NaN. Both, so a
/// kernel that reads the flag from the wrong place fails one.
fn checkFlag(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    const n = 9;
    for ([_]BlasLong{ 0, 1 }) |flag| for ([_]usize{ 1, 2 }) |inc| {
        var x: [2 * n]T = undefined;
        for (&x, 0..) |*v, i| {
            const f: F = switch (i % 4) {
                0 => 1,
                1 => std.math.nan(F),
                2 => std.math.inf(F),
                else => -2,
            };
            v.* = if (comptime complex) T.init(f, f) else f;
        }
        const px: [*]F = @ptrCast(&x);
        if (comptime complex) {
            _ = kern(n, 0, 0, 0, 0, px, @intCast(inc), px, 0, px, flag);
        } else {
            _ = kern(n, 0, 0, 0, px, @intCast(inc), px, 0, px, flag);
        }
        for (0..n) |i| {
            const v = x[i * inc];
            const re = if (comptime complex) v.re else v;
            const special = (i * inc) % 4 == 1 or (i * inc) % 4 == 2;
            const ok = if (flag == 1 and special) std.math.isNan(re) else re == 0;
            if (!ok) {
                check.say("flag {d}, element {d} is {e}: ", .{ flag, i, re });
                return check.case(error.TestUnexpectedResult, "scal, zero alpha", T, n, .{ inc, 0 });
            }
        }
    };
}

fn checkAll(vt: *const Vtable) !void {
    try checkScal(f32, vt.sscal_k.?);
    try checkScal(f64, vt.dscal_k.?);
    try checkScal(std.math.Complex(f32), vt.cscal_k.?);
    try checkScal(std.math.Complex(f64), vt.zscal_k.?);
    try checkFlag(f32, vt.sscal_k.?);
    try checkFlag(f64, vt.dscal_k.?);
    try checkFlag(std.math.Complex(f32), vt.cscal_k.?);
    try checkFlag(std.math.Complex(f64), vt.zscal_k.?);
}

test "portable scal kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkScal(XDouble, vt.qscal_k.?);
    try checkScal(std.math.Complex(XDouble), vt.xscal_k.?);
}

test "portable scal kernels at the extended type match a plain loop" {
    try checkScal(check.Ext, &Scal(check.Ext).real);
    try checkScal(std.math.Complex(check.Ext), &Scal(check.Ext).complex);
}

test "portable scal: a zero alpha with flag 1 turns NaN and infinity into NaN" {
    var x = [_]f64{ 1, std.math.nan(f64), std.math.inf(f64), -2 };
    _ = dscal_generic(4, 0, 0, 0, &x, 1, &x, 0, &x, 1);
    try check.expectCloseSlices(f64, &.{ 0, 0, 0, 0 }, &.{ x[0], 0, 0, x[3] }, 0);
    try zigstd.testExpectTrue(std.math.isNan(x[1]) and std.math.isNan(x[2]));
}

test "every linked core's scal kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
