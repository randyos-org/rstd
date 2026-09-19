//! Vtable wrappers for the `gemv` kernels.
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

pub fn shgemv_n(m: BlasLong, n: BlasLong, alpha: f32, a: [*]HF16, lda: BlasLong, x: [*]HF16, inc_x: BlasLong, beta: f32, y: [*]f32, inc_y: BlasLong) i32 {
    return blas._vtable.shgemv_n.?(m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);
}
pub fn shgemv_t(m: BlasLong, n: BlasLong, alpha: f32, a: [*]HF16, lda: BlasLong, x: [*]HF16, inc_x: BlasLong, beta: f32, y: [*]f32, inc_y: BlasLong) i32 {
    return blas._vtable.shgemv_t.?(m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);
}
pub fn bgemv_n(m: BlasLong, n: BlasLong, alpha: BF16, a: [*]BF16, lda: BlasLong, x: [*]BF16, inc_x: BlasLong, beta: BF16, y: [*]BF16, inc_y: BlasLong) i32 {
    return blas._vtable.bgemv_n.?(m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);
}
pub fn bgemv_t(m: BlasLong, n: BlasLong, alpha: BF16, a: [*]BF16, lda: BlasLong, x: [*]BF16, inc_x: BlasLong, beta: BF16, y: [*]BF16, inc_y: BlasLong) i32 {
    return blas._vtable.bgemv_t.?(m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);
}
pub fn sbgemv_n(m: BlasLong, n: BlasLong, alpha: f32, a: [*]BF16, lda: BlasLong, x: [*]BF16, inc_x: BlasLong, beta: f32, y: [*]f32, inc_y: BlasLong) i32 {
    return blas._vtable.sbgemv_n.?(m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);
}
pub fn sbgemv_t(m: BlasLong, n: BlasLong, alpha: f32, a: [*]BF16, lda: BlasLong, x: [*]BF16, inc_x: BlasLong, beta: f32, y: [*]f32, inc_y: BlasLong) i32 {
    return blas._vtable.sbgemv_t.?(m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);
}
pub fn sgemv_n(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.sgemv_n.?(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn sgemv_t(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.sgemv_t.?(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn dgemv_n(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.dgemv_n.?(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn dgemv_t(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.dgemv_t.?(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn qgemv_n(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.qgemv_n.?(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn qgemv_t(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.qgemv_t.?(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_n(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.cgemv_n.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_t(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.cgemv_t.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_r(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.cgemv_r.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_c(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.cgemv_c.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_o(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.cgemv_o.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_u(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.cgemv_u.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_s(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.cgemv_s.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_d(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas._vtable.cgemv_d.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_n(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.zgemv_n.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_t(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.zgemv_t.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_r(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.zgemv_r.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_c(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.zgemv_c.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_o(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.zgemv_o.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_u(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.zgemv_u.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_s(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.zgemv_s.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_d(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas._vtable.zgemv_d.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_n(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.xgemv_n.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_t(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.xgemv_t.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_r(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.xgemv_r.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_c(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.xgemv_c.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_o(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.xgemv_o.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_u(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.xgemv_u.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_s(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.xgemv_s.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_d(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas._vtable.xgemv_d.?(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
