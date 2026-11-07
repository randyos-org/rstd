//! `AMAX`: the largest `|x|`. The portable kernels (`*_generic`)
//! and tests of every core's kernels against a plain loop.

const types = @import("../../types.zig");

const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");

const find: extremum.Find = .{ .abs = true, .max = true };

pub const samax_generic = &extremum.Extremum(f32, 1, find).value;
pub const damax_generic = &extremum.Extremum(f64, 1, find).value;
pub const qamax_generic = &extremum.Extremum(XDouble, 1, find).value;
pub const camax_generic = &extremum.Extremum(f32, 2, find).value;
pub const zamax_generic = &extremum.Extremum(f64, 2, find).value;
pub const xamax_generic = &extremum.Extremum(XDouble, 2, find).value;

/// Points every amax slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.samax_k = samax_generic;
    vt.damax_k = damax_generic;
    vt.qamax_k = qamax_generic;
    vt.camax_k = camax_generic;
    vt.zamax_k = zamax_generic;
    vt.xamax_k = xamax_generic;
}

// ---- tests -----------------------------------------------------------------

const std = @import("std");
const check = @import("../../_test/blas1/check.zig");

fn checkAll(vt: *const Vtable) !void {
    try check.checkExtreme(f32, find, "samax", vt.samax_k.?);
    try check.checkExtreme(f64, find, "damax", vt.damax_k.?);
    try check.checkExtreme(std.math.Complex(f32), find, "camax", vt.camax_k.?);
    try check.checkExtreme(std.math.Complex(f64), find, "zamax", vt.zamax_k.?);
}

test "portable amax kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try check.checkExtreme(XDouble, find, "qamax", vt.qamax_k.?);
    try check.checkExtreme(std.math.Complex(XDouble), find, "xamax", vt.xamax_k.?);
}

test "every linked core's amax kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
