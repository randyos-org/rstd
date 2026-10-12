//! Vtable wrappers for the `ger` kernels.
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

pub fn sbger_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, a: [*]f32, lda: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.sbger_k)(m, n, dummy1, alpha, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn sger_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, a: [*]f32, lda: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.sger_k)(m, n, dummy1, alpha, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn dger_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, a: [*]f64, lda: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.dger_k)(m, n, dummy1, alpha, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn qger_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, a: [*]XDouble, lda: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.qger_k)(m, n, dummy1, alpha, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn cgeru_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, a: [*]f32, lda: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgeru_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn cgerc_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, a: [*]f32, lda: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgerc_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn cgerv_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, a: [*]f32, lda: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgerv_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn cgerd_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, a: [*]f32, lda: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgerd_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn zgeru_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, a: [*]f64, lda: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgeru_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn zgerc_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, a: [*]f64, lda: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgerc_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn zgerv_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, a: [*]f64, lda: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgerv_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn zgerd_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, a: [*]f64, lda: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgerd_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn xgeru_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, a: [*]XDouble, lda: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgeru_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn xgerc_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, a: [*]XDouble, lda: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgerc_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn xgerv_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, a: [*]XDouble, lda: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgerv_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn xgerd_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, a: [*]XDouble, lda: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgerd_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
