//! `SYMV` kernels: `y += alpha*S*x` for a symmetric `S` stored as one
//! triangle. The vtable wrappers, the portable kernels, and tests of every
//! core's kernels against plain arithmetic.
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

pub fn sbsymv_L(m: BlasLong, offset: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.sbsymv_L)(m, offset, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn sbsymv_U(m: BlasLong, offset: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.sbsymv_U)(m, offset, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn ssymv_L(m: BlasLong, offset: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.ssymv_L)(m, offset, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn ssymv_U(m: BlasLong, offset: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.ssymv_U)(m, offset, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn dsymv_L(m: BlasLong, offset: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.dsymv_L)(m, offset, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn dsymv_U(m: BlasLong, offset: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.dsymv_U)(m, offset, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn qsymv_L(m: BlasLong, offset: BlasLong, alpha: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.qsymv_L)(m, offset, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn qsymv_U(m: BlasLong, offset: BlasLong, alpha: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.qsymv_U)(m, offset, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn csymv_L(m: BlasLong, offset: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.csymv_L)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn csymv_U(m: BlasLong, offset: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.csymv_U)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zsymv_L(m: BlasLong, offset: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zsymv_L)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zsymv_U(m: BlasLong, offset: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zsymv_U)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xsymv_L(m: BlasLong, offset: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xsymv_L)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xsymv_U(m: BlasLong, offset: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xsymv_U)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const sym = @import("../symmetric.zig");
const guard = @import("../guard.zig");

/// Points every symv slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    inline for (.{ .{ "s", f32 }, .{ "d", f64 }, .{ "q", XDouble } }) |p| {
        @field(vt, p[0] ++ "symv_U") = &sym.Auto(p[0], p[1], false, .symmetric).real;
        @field(vt, p[0] ++ "symv_L") = &sym.Auto(p[0], p[1], true, .symmetric).real;
    }
    inline for (.{ .{ "c", f32 }, .{ "z", f64 }, .{ "x", XDouble } }) |p| {
        @field(vt, p[0] ++ "symv_U") = &sym.Auto(p[0], std.math.Complex(p[1]), false, .symmetric).cplx;
        @field(vt, p[0] ++ "symv_L") = &sym.Auto(p[0], std.math.Complex(p[1]), true, .symmetric).cplx;
    }
}

// ---- tests -----------------------------------------------------------------

fn checkAll(vt: *const Vtable) !void {
    guard.use(vt);
    defer guard.restore();
    try checkSlots(vt);
}

fn checkSlots(vt: *const Vtable) !void {
    inline for (.{ .{ "s", f32 }, .{ "d", f64 }, .{ "c", std.math.Complex(f32) }, .{ "z", std.math.Complex(f64) } }) |p| {
        try sym.checkBoth(p[0], p[1], false, .symmetric, @field(vt, p[0] ++ "symv_U").?, p[0] ++ "symv_U");
        try sym.checkBoth(p[0], p[1], true, .symmetric, @field(vt, p[0] ++ "symv_L").?, p[0] ++ "symv_L");
    }
}

test "portable symv kernels match plain arithmetic" {
    const vt = guard.portableVtable();
    guard.use(vt);
    defer guard.restore();
    try checkSlots(vt);
    inline for (.{ .{ "q", XDouble }, .{ "x", std.math.Complex(XDouble) } }) |p| {
        try sym.checkBoth(p[0], p[1], false, .symmetric, @field(vt, p[0] ++ "symv_U").?, p[0] ++ "symv_U");
        try sym.checkBoth(p[0], p[1], true, .symmetric, @field(vt, p[0] ++ "symv_L").?, p[0] ++ "symv_L");
    }
}

test "every linked core's symv kernels match plain arithmetic, where this machine runs them" {
    try @import("../../blas1/check.zig").eachCore(checkAll);
}
