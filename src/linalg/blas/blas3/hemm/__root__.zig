//! Vtable wrappers for the `hemm` kernels.
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

pub fn chemm_iutcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.chemm_iutcopy)(m, n, a, lda, posX, posY, b);
}
pub fn chemm_iltcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.chemm_iltcopy)(m, n, a, lda, posX, posY, b);
}
pub fn chemm_outcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.chemm_outcopy)(m, n, a, lda, posX, posY, b);
}
pub fn chemm_oltcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.chemm_oltcopy)(m, n, a, lda, posX, posY, b);
}
pub fn zhemm_iutcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zhemm_iutcopy)(m, n, a, lda, posX, posY, b);
}
pub fn zhemm_iltcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zhemm_iltcopy)(m, n, a, lda, posX, posY, b);
}
pub fn zhemm_outcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zhemm_outcopy)(m, n, a, lda, posX, posY, b);
}
pub fn zhemm_oltcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zhemm_oltcopy)(m, n, a, lda, posX, posY, b);
}
pub fn xhemm_iutcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm_iutcopy)(m, n, a, lda, posX, posY, b);
}
pub fn xhemm_iltcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm_iltcopy)(m, n, a, lda, posX, posY, b);
}
pub fn xhemm_outcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm_outcopy)(m, n, a, lda, posX, posY, b);
}
pub fn xhemm_oltcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm_oltcopy)(m, n, a, lda, posX, posY, b);
}
