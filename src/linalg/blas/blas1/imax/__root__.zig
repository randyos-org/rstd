//! `IMAX`: the 1-based index of the largest `x`. The vtable wrappers, the portable kernels
//! (`*_generic`), and tests of every core's kernels against a plain loop.
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

pub fn isbmax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas.kernel(.isbmax_k)(n, x, inc_x);
}
pub fn ismax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas.kernel(.ismax_k)(n, x, inc_x);
}
pub fn idmax_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) BlasLong {
    return blas.kernel(.idmax_k)(n, x, inc_x);
}
pub fn iqmax_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) BlasLong {
    return blas.kernel(.iqmax_k)(n, x, inc_x);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");
const check = @import("../check.zig");

const find: extremum.Find = .{ .abs = false, .max = true };

pub const ismax_generic = &extremum.Extremum(f32, 1, find).index;
pub const idmax_generic = &extremum.Extremum(f64, 1, find).index;
pub const iqmax_generic = &extremum.Extremum(XDouble, 1, find).index;

/// Points every imax slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.ismax_k = ismax_generic;
    vt.idmax_k = idmax_generic;
    vt.iqmax_k = iqmax_generic;
}

// ---- tests -----------------------------------------------------------------

fn checkAll(vt: *const Vtable) !void {
    try extremum.checkKernel(f32, find, "ismax", vt.ismax_k.?);
    try extremum.checkKernel(f64, find, "idmax", vt.idmax_k.?);
}

test "portable imax kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try extremum.checkKernel(XDouble, find, "iqmax", vt.iqmax_k.?);
}

test "every linked core's imax kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
