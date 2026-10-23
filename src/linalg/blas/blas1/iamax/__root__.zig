//! `IAMAX`: the 1-based index of the largest `|x|`. The vtable wrappers, the portable kernels
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

pub fn isbamax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas.kernel(.isbamax_k)(n, x, inc_x);
}
pub fn isamax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas.kernel(.isamax_k)(n, x, inc_x);
}
pub fn idamax_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) BlasLong {
    return blas.kernel(.idamax_k)(n, x, inc_x);
}
pub fn iqamax_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) BlasLong {
    return blas.kernel(.iqamax_k)(n, x, inc_x);
}
pub fn icamax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas.kernel(.icamax_k)(n, x, inc_x);
}
pub fn izamax_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) BlasLong {
    return blas.kernel(.izamax_k)(n, x, inc_x);
}
pub fn ixamax_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) BlasLong {
    return blas.kernel(.ixamax_k)(n, x, inc_x);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");
const check = @import("../check.zig");

const find: extremum.Find = .{ .abs = true, .max = true };

pub const isamax_generic = &extremum.Extremum(f32, 1, find).index;
pub const idamax_generic = &extremum.Extremum(f64, 1, find).index;
pub const iqamax_generic = &extremum.Extremum(XDouble, 1, find).index;
pub const icamax_generic = &extremum.Extremum(f32, 2, find).index;
pub const izamax_generic = &extremum.Extremum(f64, 2, find).index;
pub const ixamax_generic = &extremum.Extremum(XDouble, 2, find).index;

/// Points every iamax slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.isamax_k = isamax_generic;
    vt.idamax_k = idamax_generic;
    vt.iqamax_k = iqamax_generic;
    vt.icamax_k = icamax_generic;
    vt.izamax_k = izamax_generic;
    vt.ixamax_k = ixamax_generic;
}

// ---- tests -----------------------------------------------------------------

fn checkAll(vt: *const Vtable) !void {
    try check.checkExtreme(f32, find, "isamax", vt.isamax_k.?);
    try check.checkExtreme(f64, find, "idamax", vt.idamax_k.?);
    try check.checkExtreme(std.math.Complex(f32), find, "icamax", vt.icamax_k.?);
    try check.checkExtreme(std.math.Complex(f64), find, "izamax", vt.izamax_k.?);
}

test "portable iamax kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try check.checkExtreme(XDouble, find, "iqamax", vt.iqamax_k.?);
    try check.checkExtreme(std.math.Complex(XDouble), find, "ixamax", vt.ixamax_k.?);
}

test "every linked core's iamax kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
