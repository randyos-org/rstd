//! Vtable wrappers for the `hemv` kernels.
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
    return blas._vtable.chemv_L.?(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn chemv_U(m: BlasLong, offset: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.chemv_U.?(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn chemv_M(m: BlasLong, offset: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.chemv_M.?(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn chemv_V(m: BlasLong, offset: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.chemv_V.?(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zhemv_L(m: BlasLong, offset: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.zhemv_L.?(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zhemv_U(m: BlasLong, offset: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.zhemv_U.?(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zhemv_M(m: BlasLong, offset: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.zhemv_M.?(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zhemv_V(m: BlasLong, offset: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.zhemv_V.?(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xhemv_L(m: BlasLong, offset: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.xhemv_L.?(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xhemv_U(m: BlasLong, offset: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.xhemv_U.?(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xhemv_M(m: BlasLong, offset: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.xhemv_M.?(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xhemv_V(m: BlasLong, offset: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.xhemv_V.?(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
