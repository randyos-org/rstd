//! Vtable wrappers for the `omatcopy` kernels, column-major only: the row-major
//! `?omatcopy_k_r*` ones (`ORDER` `R`) are left out. The portable kernels and
//! their tests are in `../matcopy.zig`.
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

pub fn somatcopy_k_cn(rows: BlasLong, cols: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas.kernel(.somatcopy_k_cn)(rows, cols, alpha, a, lda, b, ldb);
}
pub fn somatcopy_k_ct(rows: BlasLong, cols: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas.kernel(.somatcopy_k_ct)(rows, cols, alpha, a, lda, b, ldb);
}
pub fn domatcopy_k_cn(rows: BlasLong, cols: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas.kernel(.domatcopy_k_cn)(rows, cols, alpha, a, lda, b, ldb);
}
pub fn domatcopy_k_ct(rows: BlasLong, cols: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas.kernel(.domatcopy_k_ct)(rows, cols, alpha, a, lda, b, ldb);
}
pub fn comatcopy_k_cn(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas.kernel(.comatcopy_k_cn)(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn comatcopy_k_ct(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas.kernel(.comatcopy_k_ct)(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn comatcopy_k_cnc(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas.kernel(.comatcopy_k_cnc)(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn comatcopy_k_ctc(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas.kernel(.comatcopy_k_ctc)(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn zomatcopy_k_cn(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas.kernel(.zomatcopy_k_cn)(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn zomatcopy_k_ct(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas.kernel(.zomatcopy_k_ct)(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn zomatcopy_k_cnc(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas.kernel(.zomatcopy_k_cnc)(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn zomatcopy_k_ctc(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas.kernel(.zomatcopy_k_ctc)(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
