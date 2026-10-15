//! `IMIN`: the 1-based index of the smallest `x`. The vtable wrappers, the portable kernels
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

pub fn isbmin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas.kernel(.isbmin_k)(n, x, inc_x);
}
pub fn ismin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas.kernel(.ismin_k)(n, x, inc_x);
}
pub fn idmin_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) BlasLong {
    return blas.kernel(.idmin_k)(n, x, inc_x);
}
pub fn iqmin_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) BlasLong {
    return blas.kernel(.iqmin_k)(n, x, inc_x);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");
const check = @import("../check.zig");

const find: extremum.Find = .{ .abs = false, .max = false };

pub const ismin_generic = &extremum.Extremum(f32, 1, find).index;
pub const idmin_generic = &extremum.Extremum(f64, 1, find).index;
pub const iqmin_generic = &extremum.Extremum(XDouble, 1, find).index;

/// Points every imin slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.ismin_k = ismin_generic;
    vt.idmin_k = idmin_generic;
    vt.iqmin_k = iqmin_generic;
}

// ---- tests -----------------------------------------------------------------

fn checkAll(vt: *const Vtable) !void {
    try extremum.checkKernel(f32, find, "ismin", vt.ismin_k.?);
    try extremum.checkKernel(f64, find, "idmin", vt.idmin_k.?);
}

test "portable imin kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try extremum.checkKernel(XDouble, find, "iqmin", vt.iqmin_k.?);
}

test "every linked core's imin kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
