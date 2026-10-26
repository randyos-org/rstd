//! `SUM`: the sum of `x` (`re + im` for complex). The vtable wrappers, the portable
//! kernels (`*_generic`), and tests of every core's kernels against a plain loop.
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

pub fn sbsum_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.sbsum_k)(n, x, inc_x);
}
pub fn ssum_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.ssum_k)(n, x, inc_x);
}
pub fn dsum_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas.kernel(.dsum_k)(n, x, inc_x);
}
pub fn qsum_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas.kernel(.qsum_k)(n, x, inc_x);
}
pub fn csum_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.csum_k)(n, x, inc_x);
}
pub fn zsum_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas.kernel(.zsum_k)(n, x, inc_x);
}
pub fn xsum_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas.kernel(.xsum_k)(n, x, inc_x);
}

// ---- portable kernels ------------------------------------------------------

const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

/// `kernel/arm/sum.c` and `zsum.c`: the sum of `x`, of `re + im` for complex
/// (unlike `asum`, no absolute value); an element is `w` floats.
fn Sum(comptime F: type, comptime w: comptime_int) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong) callconv(types.kernelCc(F)) F {
            if (n <= 0 or inc_x <= 0) return 0;
            var sum: types.Acc(F) = 0;
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
const check = @import("../../_test/blas1/check.zig");

/// `kern` against a plain loop in `f64`; `T` is the element type.
fn checkSum(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    for (check.sizes) |n| {
        for ([_]usize{ 1, 2, 3 }) |inc| {
            var x: [check.capacity]T = undefined;
            check.fill(T, &x, n + 1);
            const W = check.Wide(F);
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
    try checkSum(f32, vt.ssum_k.?);
    try checkSum(f64, vt.dsum_k.?);
    try checkSum(std.math.Complex(f32), vt.csum_k.?);
    try checkSum(std.math.Complex(f64), vt.zsum_k.?);
}

test "portable sum kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
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
