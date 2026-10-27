//! `MIN`: the smallest `x`. The vtable wrappers, the portable kernels
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

pub fn sbmin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.sbmin_k)(n, x, inc_x);
}
pub fn smin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.smin_k)(n, x, inc_x);
}
pub fn dmin_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas.kernel(.dmin_k)(n, x, inc_x);
}
pub fn qmin_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas.kernel(.qmin_k)(n, x, inc_x);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");

const find: extremum.Find = .{ .abs = false, .max = false };

pub const smin_generic = &extremum.Extremum(f32, 1, find).value;
pub const dmin_generic = &extremum.Extremum(f64, 1, find).value;
pub const qmin_generic = &extremum.Extremum(XDouble, 1, find).value;

/// Points every min slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.smin_k = smin_generic;
    vt.dmin_k = dmin_generic;
    vt.qmin_k = qmin_generic;
}

// ---- tests -----------------------------------------------------------------

const check = @import("../../_test/blas1/check.zig");

fn checkAll(vt: *const Vtable) !void {
    try check.checkExtreme(f32, find, "smin", vt.smin_k.?);
    try check.checkExtreme(f64, find, "dmin", vt.dmin_k.?);
}

test "portable min kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try check.checkExtreme(XDouble, find, "qmin", vt.qmin_k.?);
}

test "every linked core's min kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
