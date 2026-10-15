//! `HEMV` kernels: `y += alpha*S*x` for a Hermitian `S` stored as one
//! triangle (`V`, `M`: its transpose). The vtable wrappers, the portable
//! kernels, and tests of every core's kernels against plain arithmetic.
//!
//! Argument names are copied from the OpenBLAS kernel sources (see the
//! `CNAME(...)` definitions under `kernel/<arch>/*.c`). Several are named
//! `dummy*` in the reference: kernel-ABI padding slots this shape does not use.

const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");

const BlasLong = types.BlasLong;
const BlasInt = types.BlasInt;
const HF16 = types.HF16;
const BF16 = types.BF16;
const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

pub fn chemv_L(m: BlasLong, offset: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.chemv_L)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn chemv_U(m: BlasLong, offset: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.chemv_U)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn chemv_M(m: BlasLong, offset: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.chemv_M)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn chemv_V(m: BlasLong, offset: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.chemv_V)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zhemv_L(m: BlasLong, offset: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zhemv_L)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zhemv_U(m: BlasLong, offset: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zhemv_U)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zhemv_M(m: BlasLong, offset: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zhemv_M)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zhemv_V(m: BlasLong, offset: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zhemv_V)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xhemv_L(m: BlasLong, offset: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xhemv_L)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xhemv_U(m: BlasLong, offset: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xhemv_U)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xhemv_M(m: BlasLong, offset: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xhemv_M)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xhemv_V(m: BlasLong, offset: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xhemv_V)(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const sym = @import("../symmetric.zig");

/// The variants by suffix: lower triangle, and the transposed (`HEMVREV`)
/// matrix.
const variants = .{ .{ "U", false, sym.Kind.hermitian }, .{ "L", true, sym.Kind.hermitian }, .{ "V", false, sym.Kind.hermitian_rev }, .{ "M", true, sym.Kind.hermitian_rev } };

/// Points every hemv slot at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    inline for (.{ .{ "c", f32 }, .{ "z", f64 }, .{ "x", XDouble } }) |p| {
        inline for (variants) |v| @field(vt, p[0] ++ "hemv_" ++ v[0]) = &sym.Kernel(std.math.Complex(p[1]), v[1], v[2]).cplx;
    }
}

// ---- tests -----------------------------------------------------------------

fn checkAll(vt: *const Vtable) !void {
    inline for (.{ .{ "c", f32 }, .{ "z", f64 } }) |p| {
        inline for (variants) |v| try sym.checkKernel(std.math.Complex(p[1]), v[1], v[2], @field(vt, p[0] ++ "hemv_" ++ v[0]).?, p[0] ++ "hemv_" ++ v[0]);
    }
}

test "portable hemv kernels match plain arithmetic" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    inline for (variants) |v| try sym.checkKernel(std.math.Complex(XDouble), v[1], v[2], @field(vt, "xhemv_" ++ v[0]).?, "xhemv_" ++ v[0]);
}

test "every linked core's hemv kernels match plain arithmetic, where this machine runs them" {
    try @import("../../blas1/check.zig").eachCore(checkAll);
}
