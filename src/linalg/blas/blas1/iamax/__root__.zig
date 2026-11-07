//! `IAMAX`: the 1-based index of the largest `|x|`. The portable kernels
//! (`*_generic`) and tests of every core's kernels against a plain loop.

const types = @import("../../types.zig");

const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");

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

const std = @import("std");
const check = @import("../../_test/blas1/check.zig");

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
