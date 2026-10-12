//! Vtable wrappers for the `symv` kernels.
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
