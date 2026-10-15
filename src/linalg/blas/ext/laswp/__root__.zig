//! Vtable wrappers for the `laswp` kernels.
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

pub fn sblaswp_ncopy(n: BlasLong, k1: BlasLong, k2: BlasLong, a: [*]f32, lda: BlasLong, ipiv: [*]isize, buffer: [*]f32) i32 {
    return blas.kernel(.sblaswp_ncopy)(n, k1, k2, a, lda, ipiv, buffer);
}
pub fn slaswp_ncopy(n: BlasLong, k1: BlasLong, k2: BlasLong, a: [*]f32, lda: BlasLong, ipiv: [*]isize, buffer: [*]f32) i32 {
    return blas.kernel(.slaswp_ncopy)(n, k1, k2, a, lda, ipiv, buffer);
}
pub fn dlaswp_ncopy(n: BlasLong, k1: BlasLong, k2: BlasLong, a: [*]f64, lda: BlasLong, ipiv: [*]isize, buffer: [*]f64) i32 {
    return blas.kernel(.dlaswp_ncopy)(n, k1, k2, a, lda, ipiv, buffer);
}
pub fn qlaswp_ncopy(n: BlasLong, k1: BlasLong, k2: BlasLong, a: [*]XDouble, lda: BlasLong, ipiv: [*]isize, buffer: [*]XDouble) i32 {
    return blas.kernel(.qlaswp_ncopy)(n, k1, k2, a, lda, ipiv, buffer);
}
pub fn claswp_ncopy(n: BlasLong, k1: BlasLong, k2: BlasLong, a: [*]f32, lda: BlasLong, ipiv: [*]isize, buffer: [*]f32) i32 {
    return blas.kernel(.claswp_ncopy)(n, k1, k2, a, lda, ipiv, buffer);
}
pub fn zlaswp_ncopy(n: BlasLong, k1: BlasLong, k2: BlasLong, a: [*]f64, lda: BlasLong, ipiv: [*]isize, buffer: [*]f64) i32 {
    return blas.kernel(.zlaswp_ncopy)(n, k1, k2, a, lda, ipiv, buffer);
}
pub fn xlaswp_ncopy(n: BlasLong, k1: BlasLong, k2: BlasLong, a: [*]XDouble, lda: BlasLong, ipiv: [*]isize, buffer: [*]XDouble) i32 {
    return blas.kernel(.xlaswp_ncopy)(n, k1, k2, a, lda, ipiv, buffer);
}
