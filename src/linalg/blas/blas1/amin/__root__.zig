//! `AMIN`: the smallest `|x|`. The vtable wrappers, the portable kernels
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

pub fn sbamin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.sbamin_k)(n, x, inc_x);
}
pub fn samin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.samin_k)(n, x, inc_x);
}
pub fn damin_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas.kernel(.damin_k)(n, x, inc_x);
}
pub fn qamin_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas.kernel(.qamin_k)(n, x, inc_x);
}
pub fn camin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.camin_k)(n, x, inc_x);
}
pub fn zamin_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas.kernel(.zamin_k)(n, x, inc_x);
}
pub fn xamin_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas.kernel(.xamin_k)(n, x, inc_x);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");
const check = @import("../check.zig");

const find: extremum.Find = .{ .abs = true, .max = false };

pub const samin_generic = &extremum.Extremum(f32, 1, find).value;
pub const damin_generic = &extremum.Extremum(f64, 1, find).value;
pub const qamin_generic = &extremum.Extremum(XDouble, 1, find).value;
pub const camin_generic = &extremum.Extremum(f32, 2, find).value;
pub const zamin_generic = &extremum.Extremum(f64, 2, find).value;
pub const xamin_generic = &extremum.Extremum(XDouble, 2, find).value;

/// Points every amin slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.samin_k = samin_generic;
    vt.damin_k = damin_generic;
    vt.qamin_k = qamin_generic;
    vt.camin_k = camin_generic;
    vt.zamin_k = zamin_generic;
    vt.xamin_k = xamin_generic;
}

// ---- tests -----------------------------------------------------------------

fn checkAll(vt: *const Vtable) !void {
    try check.checkExtreme(f32, find, "samin", vt.samin_k.?);
    try check.checkExtreme(f64, find, "damin", vt.damin_k.?);
    try check.checkExtreme(std.math.Complex(f32), find, "camin", vt.camin_k.?);
    try check.checkExtreme(std.math.Complex(f64), find, "zamin", vt.zamin_k.?);
}

test "portable amin kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try check.checkExtreme(XDouble, find, "qamin", vt.qamin_k.?);
    try check.checkExtreme(std.math.Complex(XDouble), find, "xamin", vt.xamin_k.?);
}

test "every linked core's amin kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
