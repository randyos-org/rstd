//! `AMIN`: the smallest `|x|`. The portable kernels (`*_generic`)
//! and tests of every core's kernels against a plain loop.

const types = @import("../../types.zig");

const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");

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

const std = @import("std");
const check = @import("../../_test/blas1/check.zig");

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
