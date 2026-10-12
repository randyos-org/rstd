//! Vtable wrappers for the `hemm3m` kernels.
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

pub fn chemm3m_iucopyb(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.chemm3m_iucopyb)(m, n, a, lda, posX, posY, b);
}
pub fn chemm3m_ilcopyb(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.chemm3m_ilcopyb)(m, n, a, lda, posX, posY, b);
}
pub fn chemm3m_iucopyr(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.chemm3m_iucopyr)(m, n, a, lda, posX, posY, b);
}
pub fn chemm3m_ilcopyr(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.chemm3m_ilcopyr)(m, n, a, lda, posX, posY, b);
}
pub fn chemm3m_iucopyi(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.chemm3m_iucopyi)(m, n, a, lda, posX, posY, b);
}
pub fn chemm3m_ilcopyi(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.chemm3m_ilcopyi)(m, n, a, lda, posX, posY, b);
}
pub fn chemm3m_oucopyb(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.chemm3m_oucopyb)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn chemm3m_olcopyb(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.chemm3m_olcopyb)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn chemm3m_oucopyr(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.chemm3m_oucopyr)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn chemm3m_olcopyr(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.chemm3m_olcopyr)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn chemm3m_oucopyi(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.chemm3m_oucopyi)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn chemm3m_olcopyi(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.chemm3m_olcopyi)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn zhemm3m_iucopyb(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zhemm3m_iucopyb)(m, n, a, lda, posX, posY, b);
}
pub fn zhemm3m_ilcopyb(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zhemm3m_ilcopyb)(m, n, a, lda, posX, posY, b);
}
pub fn zhemm3m_iucopyr(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zhemm3m_iucopyr)(m, n, a, lda, posX, posY, b);
}
pub fn zhemm3m_ilcopyr(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zhemm3m_ilcopyr)(m, n, a, lda, posX, posY, b);
}
pub fn zhemm3m_iucopyi(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zhemm3m_iucopyi)(m, n, a, lda, posX, posY, b);
}
pub fn zhemm3m_ilcopyi(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zhemm3m_ilcopyi)(m, n, a, lda, posX, posY, b);
}
pub fn zhemm3m_oucopyb(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zhemm3m_oucopyb)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn zhemm3m_olcopyb(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zhemm3m_olcopyb)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn zhemm3m_oucopyr(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zhemm3m_oucopyr)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn zhemm3m_olcopyr(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zhemm3m_olcopyr)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn zhemm3m_oucopyi(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zhemm3m_oucopyi)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn zhemm3m_olcopyi(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zhemm3m_olcopyi)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn xhemm3m_iucopyb(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm3m_iucopyb)(m, n, a, lda, posX, posY, b);
}
pub fn xhemm3m_ilcopyb(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm3m_ilcopyb)(m, n, a, lda, posX, posY, b);
}
pub fn xhemm3m_iucopyr(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm3m_iucopyr)(m, n, a, lda, posX, posY, b);
}
pub fn xhemm3m_ilcopyr(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm3m_ilcopyr)(m, n, a, lda, posX, posY, b);
}
pub fn xhemm3m_iucopyi(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm3m_iucopyi)(m, n, a, lda, posX, posY, b);
}
pub fn xhemm3m_ilcopyi(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm3m_ilcopyi)(m, n, a, lda, posX, posY, b);
}
pub fn xhemm3m_oucopyb(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm3m_oucopyb)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn xhemm3m_olcopyb(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm3m_olcopyb)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn xhemm3m_oucopyr(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm3m_oucopyr)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn xhemm3m_olcopyr(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm3m_olcopyr)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn xhemm3m_oucopyi(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm3m_oucopyi)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn xhemm3m_olcopyi(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xhemm3m_olcopyi)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
