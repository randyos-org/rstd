//! `IAMIN`: the 1-based index of the smallest `|x|`. The vtable wrappers, the portable kernels
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

pub fn isbamin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas.kernel(.isbamin_k)(n, x, inc_x);
}
pub fn isamin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas.kernel(.isamin_k)(n, x, inc_x);
}
pub fn idamin_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) BlasLong {
    return blas.kernel(.idamin_k)(n, x, inc_x);
}
pub fn iqamin_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) BlasLong {
    return blas.kernel(.iqamin_k)(n, x, inc_x);
}
pub fn icamin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas.kernel(.icamin_k)(n, x, inc_x);
}
pub fn izamin_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) BlasLong {
    return blas.kernel(.izamin_k)(n, x, inc_x);
}
pub fn ixamin_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) BlasLong {
    return blas.kernel(.ixamin_k)(n, x, inc_x);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const extremum = @import("../extremum.zig");
const check = @import("../check.zig");

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

fn checkAll(vt: *const Vtable) !void {
    try extremum.checkKernel(f32, find, "isamin", vt.isamin_k.?);
    try extremum.checkKernel(f64, find, "idamin", vt.idamin_k.?);
    try extremum.checkKernel(std.math.Complex(f32), find, "icamin", vt.icamin_k.?);
    try extremum.checkKernel(std.math.Complex(f64), find, "izamin", vt.izamin_k.?);
}

test "portable iamin kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try extremum.checkKernel(XDouble, find, "iqamin", vt.iqamin_k.?);
    try extremum.checkKernel(std.math.Complex(XDouble), find, "ixamin", vt.ixamin_k.?);
}

test "every linked core's iamin kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
