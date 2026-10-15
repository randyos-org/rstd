//! Vtable wrappers for the `gemm3m` kernels.
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

pub fn cgemm3m_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm3m_kernel)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn cgemm3m_incopyb(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.cgemm3m_incopyb)(m, n, a, lda, b);
}
pub fn cgemm3m_incopyr(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.cgemm3m_incopyr)(m, n, a, lda, b);
}
pub fn cgemm3m_incopyi(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.cgemm3m_incopyi)(m, n, a, lda, b);
}
pub fn cgemm3m_itcopyb(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.cgemm3m_itcopyb)(m, n, a, lda, b);
}
pub fn cgemm3m_itcopyr(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.cgemm3m_itcopyr)(m, n, a, lda, b);
}
pub fn cgemm3m_itcopyi(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.cgemm3m_itcopyi)(m, n, a, lda, b);
}
pub fn cgemm3m_oncopyb(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.cgemm3m_oncopyb)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn cgemm3m_oncopyr(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.cgemm3m_oncopyr)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn cgemm3m_oncopyi(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.cgemm3m_oncopyi)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn cgemm3m_otcopyb(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.cgemm3m_otcopyb)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn cgemm3m_otcopyr(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.cgemm3m_otcopyr)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn cgemm3m_otcopyi(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.cgemm3m_otcopyi)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn zgemm3m_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm3m_kernel)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn zgemm3m_incopyb(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zgemm3m_incopyb)(m, n, a, lda, b);
}
pub fn zgemm3m_incopyr(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zgemm3m_incopyr)(m, n, a, lda, b);
}
pub fn zgemm3m_incopyi(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zgemm3m_incopyi)(m, n, a, lda, b);
}
pub fn zgemm3m_itcopyb(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zgemm3m_itcopyb)(m, n, a, lda, b);
}
pub fn zgemm3m_itcopyr(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zgemm3m_itcopyr)(m, n, a, lda, b);
}
pub fn zgemm3m_itcopyi(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zgemm3m_itcopyi)(m, n, a, lda, b);
}
pub fn zgemm3m_oncopyb(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zgemm3m_oncopyb)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn zgemm3m_oncopyr(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zgemm3m_oncopyr)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn zgemm3m_oncopyi(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zgemm3m_oncopyi)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn zgemm3m_otcopyb(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zgemm3m_otcopyb)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn zgemm3m_otcopyr(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zgemm3m_otcopyr)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn zgemm3m_otcopyi(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zgemm3m_otcopyi)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn xgemm3m_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong) i32 {
    return blas.kernel(.xgemm3m_kernel)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn xgemm3m_incopyb(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm3m_incopyb)(m, n, a, lda, b);
}
pub fn xgemm3m_incopyr(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm3m_incopyr)(m, n, a, lda, b);
}
pub fn xgemm3m_incopyi(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm3m_incopyi)(m, n, a, lda, b);
}
pub fn xgemm3m_itcopyb(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm3m_itcopyb)(m, n, a, lda, b);
}
pub fn xgemm3m_itcopyr(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm3m_itcopyr)(m, n, a, lda, b);
}
pub fn xgemm3m_itcopyi(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm3m_itcopyi)(m, n, a, lda, b);
}
pub fn xgemm3m_oncopyb(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm3m_oncopyb)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn xgemm3m_oncopyr(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm3m_oncopyr)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn xgemm3m_oncopyi(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm3m_oncopyi)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn xgemm3m_otcopyb(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm3m_otcopyb)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn xgemm3m_otcopyr(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm3m_otcopyr)(m, n, a, lda, alpha_r, alpha_i, b);
}
pub fn xgemm3m_otcopyi(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm3m_otcopyi)(m, n, a, lda, alpha_r, alpha_i, b);
}
