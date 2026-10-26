//! `MAX`: the largest `x`. The vtable wrappers, the portable kernels
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

pub fn sbmax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.sbmax_k)(n, x, inc_x);
}
pub fn smax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.smax_k)(n, x, inc_x);
}
pub fn dmax_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas.kernel(.dmax_k)(n, x, inc_x);
}
pub fn qmax_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas.kernel(.qmax_k)(n, x, inc_x);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");

const find: extremum.Find = .{ .abs = false, .max = true };

pub const smax_generic = &extremum.Extremum(f32, 1, find).value;
pub const dmax_generic = &extremum.Extremum(f64, 1, find).value;
pub const qmax_generic = &extremum.Extremum(XDouble, 1, find).value;

/// Points every max slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.smax_k = smax_generic;
    vt.dmax_k = dmax_generic;
    vt.qmax_k = qmax_generic;
}

// ---- tests -----------------------------------------------------------------

const check = @import("../../_test/blas1/check.zig");

fn checkAll(vt: *const Vtable) !void {
    try check.checkExtreme(f32, find, "smax", vt.smax_k.?);
    try check.checkExtreme(f64, find, "dmax", vt.dmax_k.?);
}

test "portable max kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try check.checkExtreme(XDouble, find, "qmax", vt.qmax_k.?);
}

test "every linked core's max kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
