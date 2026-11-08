//! `HEMV` kernels: `y += alpha*S*x` for a Hermitian `S` stored as one
//! triangle (`V`, `M`: its transpose). The portable kernels and tests of every
//! core's kernels against plain arithmetic.

const types = @import("../../types.zig");
const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const sym = @import("../symmetric.zig");
const portable = @import("../../portable.zig");

const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

/// The variants by op: lower triangle, and the transposed (`HEMVREV`)
/// matrix.
const variants = .{ .{ .hemv_U, false, sym.Kind.hermitian }, .{ .hemv_L, true, sym.Kind.hermitian }, .{ .hemv_V, false, sym.Kind.hermitian_rev }, .{ .hemv_M, true, sym.Kind.hermitian_rev } };

/// The portable kernels over `T`, any float's `Complex`, by the slot they
/// would fill (`portable.kernel`): the direct loops, which need no buffer.
pub fn Portable(comptime T: type) type {
    return struct {
        pub const hemv_U = &sym.Direct(T, false, .hermitian).cplx;
        pub const hemv_L = &sym.Direct(T, true, .hermitian).cplx;
        pub const hemv_V = &sym.Direct(T, false, .hermitian_rev).cplx;
        pub const hemv_M = &sym.Direct(T, true, .hermitian_rev).cplx;
    };
}

/// Points every hemv slot at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    inline for (.{ .{ "c", f32 }, .{ "z", f64 }, .{ "x", XDouble } }) |p| {
        inline for (variants) |v| portable.slotPtr(vt, .from(p[0]), v[0]).* = &sym.Auto(p[0], std.math.Complex(p[1]), v[1], v[2]).cplx;
    }
}

// ---- tests -----------------------------------------------------------------

fn checkAll(vt: *const Vtable) !void {
    guard.use(vt);
    defer guard.restore();
    try checkSlots(vt);
}

fn checkSlots(vt: *const Vtable) !void {
    inline for (.{ .{ "c", f32 }, .{ "z", f64 } }) |p| {
        inline for (variants) |v| try sym.checkBoth(p[0], std.math.Complex(p[1]), v[1], v[2], portable.slotPtr(vt, .from(p[0]), v[0]).*.?, portable.slotName(.from(p[0]), v[0]));
    }
}

const check = @import("../../_test/check.zig");
const guard = @import("../../_test/guard.zig");

test "portable hemv kernels match plain arithmetic" {
    const vt = guard.portableVtable();
    guard.use(vt);
    defer guard.restore();
    try checkSlots(vt);
    inline for (variants) |v| try sym.checkBoth("x", std.math.Complex(XDouble), v[1], v[2], portable.slotPtr(vt, .x, v[0]).*.?, portable.slotName(.x, v[0]));
}

test "every linked core's hemv kernels match plain arithmetic, where this machine runs them" {
    try check.eachCore(checkAll);
}
