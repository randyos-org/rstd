//! Vtable wrappers for the `neg` kernels.
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

pub fn sbneg_tcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbneg_tcopy)(m, n, a, lda, b);
}
pub fn sneg_tcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sneg_tcopy)(m, n, a, lda, b);
}
pub fn dneg_tcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dneg_tcopy)(m, n, a, lda, b);
}
pub fn qneg_tcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qneg_tcopy)(m, n, a, lda, b);
}
pub fn cneg_tcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.cneg_tcopy)(m, n, a, lda, b);
}
pub fn zneg_tcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zneg_tcopy)(m, n, a, lda, b);
}
pub fn xneg_tcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xneg_tcopy)(m, n, a, lda, b);
}
