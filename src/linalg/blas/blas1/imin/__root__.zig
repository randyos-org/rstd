//! `IMIN`: the 1-based index of the smallest `x`. The portable kernels
//! (`*_generic`) and tests of every core's kernels against a plain loop.

const types = @import("../../types.zig");

const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");

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

const check = @import("../../_test/blas1/check.zig");

fn checkAll(vt: *const Vtable) !void {
    try check.checkExtreme(f32, find, "ismin", vt.ismin_k.?);
    try check.checkExtreme(f64, find, "idmin", vt.idmin_k.?);
}

test "portable imin kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try check.checkExtreme(XDouble, find, "iqmin", vt.iqmin_k.?);
}

test "every linked core's imin kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
