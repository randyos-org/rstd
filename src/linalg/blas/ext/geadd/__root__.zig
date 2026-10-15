//! Vtable wrappers for the `geadd` kernels.
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

pub fn sgeadd_k(rows: BlasLong, cols: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, beta: f32, b: [*]f32, ldb: BlasLong, transa: BlasLong, transb: BlasLong) i32 {
    return blas.kernel(.sgeadd_k)(rows, cols, alpha, a, lda, beta, b, ldb, transa, transb);
}
pub fn dgeadd_k(rows: BlasLong, cols: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, beta: f64, b: [*]f64, ldb: BlasLong, transa: BlasLong, transb: BlasLong) i32 {
    return blas.kernel(.dgeadd_k)(rows, cols, alpha, a, lda, beta, b, ldb, transa, transb);
}
pub fn cgeadd_k(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, beta_r: f32, beta_i: f32, b: [*]f32, ldb: BlasLong, transa: BlasLong, transb: BlasLong) i32 {
    return blas.kernel(.cgeadd_k)(rows, cols, alpha_r, alpha_i, a, lda, beta_r, beta_i, b, ldb, transa, transb);
}
pub fn zgeadd_k(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, beta_r: f64, beta_i: f64, b: [*]f64, ldb: BlasLong, transa: BlasLong, transb: BlasLong) i32 {
    return blas.kernel(.zgeadd_k)(rows, cols, alpha_r, alpha_i, a, lda, beta_r, beta_i, b, ldb, transa, transb);
}
