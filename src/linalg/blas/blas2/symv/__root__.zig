//! `SYMV` kernels: `y += alpha*S*x` for a symmetric `S` stored as one
//! triangle. The portable kernels and tests of every core's kernels against
//! plain arithmetic.

const types = @import("../../types.zig");
const std = @import("std");
const sysinfo = @import("builtin");
const x86_64 = @import("x86_64/kernel.zig");
const ported = @import("driver.zig");
const Vtable = @import("../../Vtable.zig");
const sym = @import("../symmetric.zig");
const portable = @import("../../portable.zig");

const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`): the direct loops, which need no
/// buffer.
pub fn Portable(comptime T: type) type {
    const c = types.isComplex(T);
    return struct {
        pub const symv_U = if (c) &sym.Direct(T, false, .symmetric).cplx else &sym.Direct(T, false, .symmetric).real;
        pub const symv_L = if (c) &sym.Direct(T, true, .symmetric).cplx else &sym.Direct(T, true, .symmetric).real;
    };
}

/// Points every symv slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    inline for (.{ .{ "s", f32 }, .{ "d", f64 }, .{ "q", XDouble } }) |p| {
        portable.slotPtr(vt, .from(p[0]), .symv_U).* = &sym.Auto(p[0], p[1], false, .symmetric).real;
        portable.slotPtr(vt, .from(p[0]), .symv_L).* = &sym.Auto(p[0], p[1], true, .symmetric).real;
    }
    inline for (.{ .{ "c", f32 }, .{ "z", f64 }, .{ "x", XDouble } }) |p| {
        portable.slotPtr(vt, .from(p[0]), .symv_U).* = &sym.Auto(p[0], std.math.Complex(p[1]), false, .symmetric).cplx;
        portable.slotPtr(vt, .from(p[0]), .symv_L).* = &sym.Auto(p[0], std.math.Complex(p[1]), true, .symmetric).cplx;
    }
}

// ---- tests -----------------------------------------------------------------

fn checkAll(vt: *const Vtable) !void {
    guard.use(vt);
    defer guard.restore();
    try checkSlots(vt);
}

fn checkSlots(vt: *const Vtable) !void {
    inline for (check.sd_letters) |p| {
        try sym.checkBoth(p[0], p[1], false, .symmetric, portable.slotPtr(vt, .from(p[0]), .symv_U).*.?, p[0] ++ "symv_U");
        try sym.checkBoth(p[0], p[1], true, .symmetric, portable.slotPtr(vt, .from(p[0]), .symv_L).*.?, p[0] ++ "symv_L");
    }
}

const check = @import("../../_test/check.zig");
const guard = @import("../../_test/guard.zig");

test "portable symv kernels match plain arithmetic" {
    const vt = check.portableVtable();
    guard.use(vt);
    defer guard.restore();
    try checkSlots(vt);
    inline for (.{ .{ "q", XDouble }, .{ "x", std.math.Complex(XDouble) } }) |p| {
        try sym.checkBoth(p[0], p[1], false, .symmetric, portable.slotPtr(vt, .from(p[0]), .symv_U).*.?, p[0] ++ "symv_U");
        try sym.checkBoth(p[0], p[1], true, .symmetric, portable.slotPtr(vt, .from(p[0]), .symv_L).*.?, p[0] ++ "symv_L");
    }
}

test "every linked core's symv kernels match plain arithmetic, where this machine runs them" {
    try check.eachCore(checkAll);
}

test "x86_64's symv ports, every core's block shape, match plain arithmetic as @Vector code on this machine" {
    if (comptime sysinfo.target.cpu.arch != .x86_64) return error.SkipZigTest;
    const P = x86_64;
    inline for (comptime std.meta.tags(P.Micro)) |micro| {
        var vt = check.portableVtable().*;
        inline for (.{ .ssymv_U, .ssymv_L, .dsymv_U, .dsymv_L }) |s| @field(vt, @tagName(s)) = P.Ported(micro).entry(s);
        try checkAll(&vt);
    }
}

test "arm64's asimd symv port matches plain arithmetic, as @Vector code on this machine" {
    if (comptime sysinfo.target.cpu.arch != .aarch64) return error.SkipZigTest;
    var vt = check.portableVtable().*;
    inline for (.{ .ssymv_U, .ssymv_L, .dsymv_U, .dsymv_L }) |s| @field(vt, @tagName(s)) = ported.arm64(s);
    try checkAll(&vt);
}
