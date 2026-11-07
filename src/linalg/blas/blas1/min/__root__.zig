//! `MIN`: the smallest `x`. The portable kernels (`*_generic`)
//! and tests of every core's kernels against a plain loop.

const types = @import("../../types.zig");

const XDouble = types.XDouble;

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
