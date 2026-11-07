//! `IMAX`: the 1-based index of the largest `x`. The portable kernels
//! (`*_generic`) and tests of every core's kernels against a plain loop.

const types = @import("../../types.zig");

const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");

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

const check = @import("../../_test/blas1/check.zig");

fn checkAll(vt: *const Vtable) !void {
    try check.checkExtreme(f32, find, "ismax", vt.ismax_k.?);
    try check.checkExtreme(f64, find, "idmax", vt.idmax_k.?);
}

test "portable imax kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try check.checkExtreme(XDouble, find, "iqmax", vt.iqmax_k.?);
}

test "every linked core's imax kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
