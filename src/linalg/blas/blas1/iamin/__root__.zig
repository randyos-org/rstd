//! `IAMIN`: the 1-based index of the smallest `|x|`. The portable kernels
//! (`*_generic`) and tests of every core's kernels against a plain loop.

const types = @import("../../types.zig");

const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");

const find: extremum.Find = .{ .abs = true, .max = false };

pub const isamin_generic = &extremum.Extremum(f32, 1, find).index;
pub const idamin_generic = &extremum.Extremum(f64, 1, find).index;
pub const iqamin_generic = &extremum.Extremum(XDouble, 1, find).index;
pub const icamin_generic = &extremum.Extremum(f32, 2, find).index;
pub const izamin_generic = &extremum.Extremum(f64, 2, find).index;
pub const ixamin_generic = &extremum.Extremum(XDouble, 2, find).index;

/// Points every iamin slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.isamin_k = isamin_generic;
    vt.idamin_k = idamin_generic;
    vt.iqamin_k = iqamin_generic;
    vt.icamin_k = icamin_generic;
    vt.izamin_k = izamin_generic;
    vt.ixamin_k = ixamin_generic;
}

// ---- tests -----------------------------------------------------------------

const std = @import("std");
const check = @import("../../_test/blas1/check.zig");

fn checkAll(vt: *const Vtable) !void {
    try check.checkExtreme(f32, find, "isamin", vt.isamin_k.?);
    try check.checkExtreme(f64, find, "idamin", vt.idamin_k.?);
    try check.checkExtreme(std.math.Complex(f32), find, "icamin", vt.icamin_k.?);
    try check.checkExtreme(std.math.Complex(f64), find, "izamin", vt.izamin_k.?);
}

test "portable iamin kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try check.checkExtreme(XDouble, find, "iqamin", vt.iqamin_k.?);
    try check.checkExtreme(std.math.Complex(XDouble), find, "ixamin", vt.ixamin_k.?);
}

test "every linked core's iamin kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
