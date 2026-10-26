//! `COPY`: `y := x`. The vtable wrappers, the portable kernels (`*_generic`),
//! and tests of every core's kernels against a plain loop.
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

pub fn sbcopy_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) i32 {
    return blas.kernel(.sbcopy_k)(n, x, inc_x, y, inc_y);
}
pub fn scopy_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) i32 {
    return blas.kernel(.scopy_k)(n, x, inc_x, y, inc_y);
}
pub fn dcopy_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong) i32 {
    return blas.kernel(.dcopy_k)(n, x, inc_x, y, inc_y);
}
pub fn qcopy_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong) i32 {
    return blas.kernel(.qcopy_k)(n, x, inc_x, y, inc_y);
}
pub fn ccopy_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) i32 {
    return blas.kernel(.ccopy_k)(n, x, inc_x, y, inc_y);
}
pub fn zcopy_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong) i32 {
    return blas.kernel(.zcopy_k)(n, x, inc_x, y, inc_y);
}
pub fn xcopy_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong) i32 {
    return blas.kernel(.xcopy_k)(n, x, inc_x, y, inc_y);
}

// ---- portable kernels ------------------------------------------------------

const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

/// `kernel/arm/copy.c` and `zcopy.c`: `y := x`; an element is `w` floats.
fn Copy(comptime F: type, comptime w: comptime_int) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) callconv(types.kernelCc(F)) i32 {
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                inline for (0..w) |c| util.at(F, y, i * inc_y * w + c).* = util.at(F, x, i * inc_x * w + c).*;
            }
            return 0;
        }
    };
}

pub const scopy_generic = &Copy(f32, 1).k;
pub const dcopy_generic = &Copy(f64, 1).k;
pub const qcopy_generic = &Copy(XDouble, 1).k;
pub const ccopy_generic = &Copy(f32, 2).k;
pub const zcopy_generic = &Copy(f64, 2).k;
pub const xcopy_generic = &Copy(XDouble, 2).k;

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    const w = if (types.isComplex(T)) 2 else 1;
    return struct {
        pub const copy_k = &Copy(F, w).k;
    };
}

/// Points every copy slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.scopy_k = scopy_generic;
    vt.dcopy_k = dcopy_generic;
    vt.qcopy_k = qcopy_generic;
    vt.ccopy_k = ccopy_generic;
    vt.zcopy_k = zcopy_generic;
    vt.xcopy_k = xcopy_generic;
}

// ---- tests -----------------------------------------------------------------

const std = @import("std");
const check = @import("../../_test/blas1/check.zig");

/// `kern` against a plain loop; `T` is the element type, real or complex.
fn checkCopy(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    for (check.sizes) |n| {
        var x: [check.capacity]T = undefined;
        var y: [check.capacity]T = undefined;
        check.fill(T, &x, n + 1);
        check.fill(T, &y, n + 7);
        for (check.zero_strides) |inc| {
            var want = y;
            var got = y;
            for (0..n) |i| want[check.pos(n, inc[1], i)] = x[check.pos(n, inc[0], i)];
            const px: [*]F = @ptrCast(&x[check.pos(n, inc[0], 0)]);
            _ = kern(@intCast(n), px, inc[0], @ptrCast(&got[check.pos(n, inc[1], 0)]), inc[1]);
            check.expectCloseSlices(T, &want, &got, 0) catch |e| return check.case(e, "copy", T, n, inc);
        }
    }
}

fn checkAll(vt: *const Vtable) !void {
    var failed: ?anyerror = null;
    check.keep(&failed, checkCopy(f32, vt.scopy_k.?));
    check.keep(&failed, checkCopy(f64, vt.dcopy_k.?));
    check.keep(&failed, checkCopy(std.math.Complex(f32), vt.ccopy_k.?));
    check.keep(&failed, checkCopy(std.math.Complex(f64), vt.zcopy_k.?));
    if (failed) |e| return e;
}

test "portable copy kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkCopy(XDouble, vt.qcopy_k.?);
    try checkCopy(std.math.Complex(XDouble), vt.xcopy_k.?);
}

test "portable copy kernels at the extended type match a plain loop" {
    try checkCopy(check.Ext, &Copy(check.Ext, 1).k);
    try checkCopy(std.math.Complex(check.Ext), &Copy(check.Ext, 2).k);
}

test "every linked core's copy kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
