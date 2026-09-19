//! Vtable wrappers for the `symm` kernels.
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

pub fn sbsymm_iutcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbsymm_iutcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbsymm_iltcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbsymm_iltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbsymm_outcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbsymm_outcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbsymm_oltcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbsymm_oltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ssymm_direct_alpha_betaLU(M: BlasLong, N: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, B: [*]f32, strideB: BlasLong, beta: f32, R: [*]f32, strideR: BlasLong) void {
    blas._vtable.ssymm_direct_alpha_betaLU.?(M, N, alpha, A, strideA, B, strideB, beta, R, strideR);
}
pub fn ssymm_direct_alpha_betaLL(M: BlasLong, N: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, B: [*]f32, strideB: BlasLong, beta: f32, R: [*]f32, strideR: BlasLong) void {
    blas._vtable.ssymm_direct_alpha_betaLL.?(M, N, alpha, A, strideA, B, strideB, beta, R, strideR);
}
pub fn ssymm_iutcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ssymm_iutcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ssymm_iltcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ssymm_iltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ssymm_outcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ssymm_outcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ssymm_oltcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ssymm_oltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dsymm_iutcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dsymm_iutcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dsymm_iltcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dsymm_iltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dsymm_outcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dsymm_outcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dsymm_oltcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dsymm_oltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qsymm_iutcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qsymm_iutcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qsymm_iltcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qsymm_iltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qsymm_outcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qsymm_outcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qsymm_oltcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qsymm_oltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn csymm_iutcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.csymm_iutcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn csymm_iltcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.csymm_iltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn csymm_outcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.csymm_outcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn csymm_oltcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.csymm_oltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn zsymm_iutcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.zsymm_iutcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn zsymm_iltcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.zsymm_iltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn zsymm_outcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.zsymm_outcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn zsymm_oltcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.zsymm_oltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xsymm_iutcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xsymm_iutcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xsymm_iltcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xsymm_iltcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xsymm_outcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xsymm_outcopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xsymm_oltcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xsymm_oltcopy.?(m, n, a, lda, posX, posY, b);
}
