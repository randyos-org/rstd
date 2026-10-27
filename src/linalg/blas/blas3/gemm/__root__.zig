//! Vtable wrappers for the `gemm` kernels.
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

pub fn shgemm_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f32, ba: [*]HF16, bb: [*]HF16, c: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.shgemm_kernel)(bm, bn, bk, alpha, ba, bb, c, ldc);
}
pub fn shgemm_beta(m: BlasLong, n: BlasLong, dummy1: BlasLong, beta: f32, dummy2: [*]HF16, dummy3: BlasLong, dummy4: [*]HF16, dummy5: BlasLong, c: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.shgemm_beta)(m, n, dummy1, beta, dummy2, dummy3, dummy4, dummy5, c, ldc);
}
pub fn shgemm_incopy(m: BlasLong, n: BlasLong, a: [*]HF16, lda: BlasLong, b: [*]HF16) i32 {
    return blas.kernel(.shgemm_incopy)(m, n, a, lda, b);
}
pub fn shgemm_itcopy(m: BlasLong, n: BlasLong, a: [*]HF16, lda: BlasLong, b: [*]HF16) i32 {
    return blas.kernel(.shgemm_itcopy)(m, n, a, lda, b);
}
pub fn shgemm_oncopy(m: BlasLong, n: BlasLong, a: [*]HF16, lda: BlasLong, b: [*]HF16) i32 {
    return blas.kernel(.shgemm_oncopy)(m, n, a, lda, b);
}
pub fn shgemm_otcopy(m: BlasLong, n: BlasLong, a: [*]HF16, lda: BlasLong, b: [*]HF16) i32 {
    return blas.kernel(.shgemm_otcopy)(m, n, a, lda, b);
}
pub fn bgemm_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: BF16, ba: [*]BF16, bb: [*]BF16, c: [*]BF16, ldc: BlasLong) i32 {
    return blas.kernel(.bgemm_kernel)(bm, bn, bk, alpha, ba, bb, c, ldc);
}
pub fn bgemm_beta(m: BlasLong, n: BlasLong, dummy1: BlasLong, beta: BF16, dummy2: [*]BF16, dummy3: BlasLong, dummy4: [*]BF16, dummy5: BlasLong, c: [*]BF16, ldc: BlasLong) i32 {
    return blas.kernel(.bgemm_beta)(m, n, dummy1, beta, dummy2, dummy3, dummy4, dummy5, c, ldc);
}
pub fn bgemm_incopy(m: BlasLong, n: BlasLong, a: [*]BF16, lda: BlasLong, b: [*]BF16) i32 {
    return blas.kernel(.bgemm_incopy)(m, n, a, lda, b);
}
pub fn bgemm_itcopy(m: BlasLong, n: BlasLong, a: [*]BF16, lda: BlasLong, b: [*]BF16) i32 {
    return blas.kernel(.bgemm_itcopy)(m, n, a, lda, b);
}
pub fn bgemm_oncopy(m: BlasLong, n: BlasLong, a: [*]BF16, lda: BlasLong, b: [*]BF16) i32 {
    return blas.kernel(.bgemm_oncopy)(m, n, a, lda, b);
}
pub fn bgemm_otcopy(m: BlasLong, n: BlasLong, a: [*]BF16, lda: BlasLong, b: [*]BF16) i32 {
    return blas.kernel(.bgemm_otcopy)(m, n, a, lda, b);
}
pub fn sbgemm_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f32, ba: [*]BF16, bb: [*]BF16, c: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sbgemm_kernel)(bm, bn, bk, alpha, ba, bb, c, ldc);
}
pub fn sbgemm_beta(m: BlasLong, n: BlasLong, dummy1: BlasLong, beta: f32, dummy2: [*]BF16, dummy3: BlasLong, dummy4: [*]BF16, dummy5: BlasLong, c: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sbgemm_beta)(m, n, dummy1, beta, dummy2, dummy3, dummy4, dummy5, c, ldc);
}
pub fn sbgemm_incopy(m: BlasLong, n: BlasLong, a: [*]BF16, lda: BlasLong, b: [*]BF16) i32 {
    return blas.kernel(.sbgemm_incopy)(m, n, a, lda, b);
}
pub fn sbgemm_itcopy(m: BlasLong, n: BlasLong, a: [*]BF16, lda: BlasLong, b: [*]BF16) i32 {
    return blas.kernel(.sbgemm_itcopy)(m, n, a, lda, b);
}
pub fn sbgemm_oncopy(m: BlasLong, n: BlasLong, a: [*]BF16, lda: BlasLong, b: [*]BF16) i32 {
    return blas.kernel(.sbgemm_oncopy)(m, n, a, lda, b);
}
pub fn sbgemm_otcopy(m: BlasLong, n: BlasLong, a: [*]BF16, lda: BlasLong, b: [*]BF16) i32 {
    return blas.kernel(.sbgemm_otcopy)(m, n, a, lda, b);
}
pub fn sbgemm_small_matrix_permit(transa: i32, transb: i32, m: BlasLong, n: BlasLong, k: BlasLong, alpha: f32, beta: f32) i32 {
    return blas.kernel(.sbgemm_small_matrix_permit)(transa, transb, m, n, k, alpha, beta);
}
pub fn sbgemm_small_kernel_nn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]BF16, lda: BlasLong, alpha: f32, B: [*]BF16, ldb: BlasLong, beta: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sbgemm_small_kernel_nn)(m, n, k, A, lda, alpha, B, ldb, beta, C, ldc);
}
pub fn sbgemm_small_kernel_nt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]BF16, lda: BlasLong, alpha: f32, B: [*]BF16, ldb: BlasLong, beta: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sbgemm_small_kernel_nt)(m, n, k, A, lda, alpha, B, ldb, beta, C, ldc);
}
pub fn sbgemm_small_kernel_tn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]BF16, lda: BlasLong, alpha: f32, B: [*]BF16, ldb: BlasLong, beta: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sbgemm_small_kernel_tn)(m, n, k, A, lda, alpha, B, ldb, beta, C, ldc);
}
pub fn sbgemm_small_kernel_tt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]BF16, lda: BlasLong, alpha: f32, B: [*]BF16, ldb: BlasLong, beta: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sbgemm_small_kernel_tt)(m, n, k, A, lda, alpha, B, ldb, beta, C, ldc);
}
pub fn sbgemm_small_kernel_b0_nn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]BF16, lda: BlasLong, alpha: f32, B: [*]BF16, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sbgemm_small_kernel_b0_nn)(m, n, k, A, lda, alpha, B, ldb, C, ldc);
}
pub fn sbgemm_small_kernel_b0_nt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]BF16, lda: BlasLong, alpha: f32, B: [*]BF16, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sbgemm_small_kernel_b0_nt)(m, n, k, A, lda, alpha, B, ldb, C, ldc);
}
pub fn sbgemm_small_kernel_b0_tn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]BF16, lda: BlasLong, alpha: f32, B: [*]BF16, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sbgemm_small_kernel_b0_tn)(m, n, k, A, lda, alpha, B, ldb, C, ldc);
}
pub fn sbgemm_small_kernel_b0_tt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]BF16, lda: BlasLong, alpha: f32, B: [*]BF16, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sbgemm_small_kernel_b0_tt)(m, n, k, A, lda, alpha, B, ldb, C, ldc);
}
pub fn sgemm_direct(M: BlasLong, N: BlasLong, K: BlasLong, A: [*]f32, strideA: BlasLong, B: [*]f32, strideB: BlasLong, R: [*]f32, strideR: BlasLong) void {
    blas.kernel(.sgemm_direct)(M, N, K, A, strideA, B, strideB, R, strideR);
}
pub fn sgemm_direct_performant(M: BlasLong, N: BlasLong, K: BlasLong) i32 {
    return blas.kernel(.sgemm_direct_performant)(M, N, K);
}
pub fn sgemm_direct_alpha_beta(M: BlasLong, N: BlasLong, K: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, B: [*]f32, strideB: BlasLong, beta: f32, R: [*]f32, strideR: BlasLong) void {
    blas.kernel(.sgemm_direct_alpha_beta)(M, N, K, alpha, A, strideA, B, strideB, beta, R, strideR);
}
pub fn sme_sgemm_kernel(ta: [*:0]u8, tb: [*:0]u8, m: BlasLong, n: BlasLong, k: BlasLong, alpha: [*]f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong, beta: [*]f32, c: [*]f32, ldc: BlasLong) void {
    blas.kernel(.sme_sgemm_kernel)(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}
pub fn sgemm_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sgemm_kernel)(bm, bn, bk, alpha, ba, bb, c, ldc);
}
pub fn sgemm_beta(m: BlasLong, n: BlasLong, dummy1: BlasLong, beta: f32, dummy2: [*]f32, dummy3: BlasLong, dummy4: [*]f32, dummy5: BlasLong, c: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sgemm_beta)(m, n, dummy1, beta, dummy2, dummy3, dummy4, dummy5, c, ldc);
}
pub fn sgemm_incopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sgemm_incopy)(m, n, a, lda, b);
}
pub fn sgemm_itcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sgemm_itcopy)(m, n, a, lda, b);
}
pub fn sgemm_oncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sgemm_oncopy)(m, n, a, lda, b);
}
pub fn sgemm_otcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sgemm_otcopy)(m, n, a, lda, b);
}
pub fn sgemm_small_matrix_permit(transa: i32, transb: i32, m: BlasLong, n: BlasLong, k: BlasLong, alpha: f32, beta: f32) i32 {
    return blas.kernel(.sgemm_small_matrix_permit)(transa, transb, m, n, k, alpha, beta);
}
pub fn sgemm_small_kernel_nn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha: f32, B: [*]f32, ldb: BlasLong, beta: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sgemm_small_kernel_nn)(m, n, k, A, lda, alpha, B, ldb, beta, C, ldc);
}
pub fn sgemm_small_kernel_nt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha: f32, B: [*]f32, ldb: BlasLong, beta: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sgemm_small_kernel_nt)(m, n, k, A, lda, alpha, B, ldb, beta, C, ldc);
}
pub fn sgemm_small_kernel_tn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha: f32, B: [*]f32, ldb: BlasLong, beta: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sgemm_small_kernel_tn)(m, n, k, A, lda, alpha, B, ldb, beta, C, ldc);
}
pub fn sgemm_small_kernel_tt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha: f32, B: [*]f32, ldb: BlasLong, beta: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sgemm_small_kernel_tt)(m, n, k, A, lda, alpha, B, ldb, beta, C, ldc);
}
pub fn sgemm_small_kernel_b0_nn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sgemm_small_kernel_b0_nn)(m, n, k, A, lda, alpha, B, ldb, C, ldc);
}
pub fn sgemm_small_kernel_b0_nt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sgemm_small_kernel_b0_nt)(m, n, k, A, lda, alpha, B, ldb, C, ldc);
}
pub fn sgemm_small_kernel_b0_tn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sgemm_small_kernel_b0_tn)(m, n, k, A, lda, alpha, B, ldb, C, ldc);
}
pub fn sgemm_small_kernel_b0_tt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.sgemm_small_kernel_b0_tt)(m, n, k, A, lda, alpha, B, ldb, C, ldc);
}
pub fn sme_dgemm_kernel(ta: [*:0]const u8, tb: [*:0]const u8, m: BlasLong, n: BlasLong, k: BlasLong, alpha: [*]const f64, a: [*]const f64, lda: BlasLong, b: [*]const f64, ldb: BlasLong, beta: [*]const f64, c: [*]f64, ldc: BlasLong) void {
    blas.kernel(.sme_dgemm_kernel)(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}
pub fn dgemm_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.dgemm_kernel)(bm, bn, bk, alpha, ba, bb, c, ldc);
}
pub fn dgemm_beta(m: BlasLong, n: BlasLong, dummy1: BlasLong, beta: f64, dummy2: [*]f64, dummy3: BlasLong, dummy4: [*]f64, dummy5: BlasLong, c: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.dgemm_beta)(m, n, dummy1, beta, dummy2, dummy3, dummy4, dummy5, c, ldc);
}
pub fn dgemm_incopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dgemm_incopy)(m, n, a, lda, b);
}
pub fn dgemm_itcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dgemm_itcopy)(m, n, a, lda, b);
}
pub fn dgemm_oncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dgemm_oncopy)(m, n, a, lda, b);
}
pub fn dgemm_otcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dgemm_otcopy)(m, n, a, lda, b);
}
pub fn dgemm_small_matrix_permit(transa: i32, transb: i32, m: BlasLong, n: BlasLong, k: BlasLong, alpha: f64, beta: f64) i32 {
    return blas.kernel(.dgemm_small_matrix_permit)(transa, transb, m, n, k, alpha, beta);
}
pub fn dgemm_small_kernel_nn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha: f64, B: [*]f64, ldb: BlasLong, beta: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.dgemm_small_kernel_nn)(m, n, k, A, lda, alpha, B, ldb, beta, C, ldc);
}
pub fn dgemm_small_kernel_nt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha: f64, B: [*]f64, ldb: BlasLong, beta: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.dgemm_small_kernel_nt)(m, n, k, A, lda, alpha, B, ldb, beta, C, ldc);
}
pub fn dgemm_small_kernel_tn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha: f64, B: [*]f64, ldb: BlasLong, beta: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.dgemm_small_kernel_tn)(m, n, k, A, lda, alpha, B, ldb, beta, C, ldc);
}
pub fn dgemm_small_kernel_tt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha: f64, B: [*]f64, ldb: BlasLong, beta: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.dgemm_small_kernel_tt)(m, n, k, A, lda, alpha, B, ldb, beta, C, ldc);
}
pub fn dgemm_small_kernel_b0_nn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.dgemm_small_kernel_b0_nn)(m, n, k, A, lda, alpha, B, ldb, C, ldc);
}
pub fn dgemm_small_kernel_b0_nt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.dgemm_small_kernel_b0_nt)(m, n, k, A, lda, alpha, B, ldb, C, ldc);
}
pub fn dgemm_small_kernel_b0_tn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.dgemm_small_kernel_b0_tn)(m, n, k, A, lda, alpha, B, ldb, C, ldc);
}
pub fn dgemm_small_kernel_b0_tt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.dgemm_small_kernel_b0_tt)(m, n, k, A, lda, alpha, B, ldb, C, ldc);
}
pub fn qgemm_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong) i32 {
    return blas.kernel(.qgemm_kernel)(bm, bn, bk, alpha, ba, bb, c, ldc);
}
pub fn qgemm_beta(m: BlasLong, n: BlasLong, dummy1: BlasLong, beta: XDouble, dummy2: [*]XDouble, dummy3: BlasLong, dummy4: [*]XDouble, dummy5: BlasLong, c: [*]XDouble, ldc: BlasLong) i32 {
    return blas.kernel(.qgemm_beta)(m, n, dummy1, beta, dummy2, dummy3, dummy4, dummy5, c, ldc);
}
pub fn qgemm_incopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qgemm_incopy)(m, n, a, lda, b);
}
pub fn qgemm_itcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qgemm_itcopy)(m, n, a, lda, b);
}
pub fn qgemm_oncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qgemm_oncopy)(m, n, a, lda, b);
}
pub fn qgemm_otcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qgemm_otcopy)(m, n, a, lda, b);
}
pub fn sme_cgemm_kernel(ta: [*:0]const u8, tb: [*:0]const u8, m: BlasLong, n: BlasLong, k: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]const f32, lda: BlasLong, b: [*]const f32, ldb: BlasLong, beta_r: f32, beta_i: f32, c: [*]f32, ldc: BlasLong) void {
    blas.kernel(.sme_cgemm_kernel)(ta, tb, m, n, k, alpha_r, alpha_i, a, lda, b, ldb, beta_r, beta_i, c, ldc);
}
pub fn cgemm_kernel_n(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_kernel_n)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn cgemm_kernel_l(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_kernel_l)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn cgemm_kernel_r(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_kernel_r)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn cgemm_kernel_b(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_kernel_b)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn cgemm_beta(m: BlasLong, n: BlasLong, dummy1: BlasLong, beta_r: f32, beta_i: f32, dummy2: [*]f32, dummy3: BlasLong, dummy4: [*]f32, dummy5: BlasLong, c: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_beta)(m, n, dummy1, beta_r, beta_i, dummy2, dummy3, dummy4, dummy5, c, ldc);
}
pub fn cgemm_incopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.cgemm_incopy)(m, n, a, lda, b);
}
pub fn cgemm_itcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.cgemm_itcopy)(m, n, a, lda, b);
}
pub fn cgemm_oncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.cgemm_oncopy)(m, n, a, lda, b);
}
pub fn cgemm_otcopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.cgemm_otcopy)(m, n, a, lda, b);
}
pub fn cgemm_small_matrix_permit(transa: i32, transb: i32, m: BlasLong, n: BlasLong, k: BlasLong, alpha0: f32, alpha1: f32, beta0: f32, beta1: f32) i32 {
    return blas.kernel(.cgemm_small_matrix_permit)(transa, transb, m, n, k, alpha0, alpha1, beta0, beta1);
}
pub fn cgemm_small_kernel_nn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_nn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_nt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_nt)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_nr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_nr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_nc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_nc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_tn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_tn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_tt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_tt)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_tr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_tr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_tc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_tc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_rn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_rn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_rt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_rt)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_rr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_rr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_rc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_rc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_cn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_cn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_ct(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_ct)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_cr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_cr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_cc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, beta0: f32, beta1: f32, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_cc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn cgemm_small_kernel_b0_nn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_nn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_nt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_nt)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_nr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_nr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_nc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_nc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_tn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_tn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_tt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_tt)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_tr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_tr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_tc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_tc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_rn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_rn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_rt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_rt)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_rr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_rr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_rc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_rc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_cn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_cn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_ct(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_ct)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_cr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_cr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn cgemm_small_kernel_b0_cc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f32, lda: BlasLong, alpha0: f32, alpha1: f32, B: [*]f32, ldb: BlasLong, C: [*]f32, ldc: BlasLong) i32 {
    return blas.kernel(.cgemm_small_kernel_b0_cc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn sme_zgemm_kernel(ta: [*:0]const u8, tb: [*:0]const u8, m: BlasLong, n: BlasLong, k: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]const f64, lda: BlasLong, b: [*]const f64, ldb: BlasLong, beta_r: f64, beta_i: f64, c: [*]f64, ldc: BlasLong) void {
    blas.kernel(.sme_zgemm_kernel)(ta, tb, m, n, k, alpha_r, alpha_i, a, lda, b, ldb, beta_r, beta_i, c, ldc);
}
pub fn zgemm_kernel_n(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_kernel_n)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn zgemm_kernel_l(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_kernel_l)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn zgemm_kernel_r(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_kernel_r)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn zgemm_kernel_b(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_kernel_b)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn zgemm_beta(m: BlasLong, n: BlasLong, dummy1: BlasLong, beta_r: f64, beta_i: f64, dummy2: [*]f64, dummy3: BlasLong, dummy4: [*]f64, dummy5: BlasLong, c: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_beta)(m, n, dummy1, beta_r, beta_i, dummy2, dummy3, dummy4, dummy5, c, ldc);
}
pub fn zgemm_incopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zgemm_incopy)(m, n, a, lda, b);
}
pub fn zgemm_itcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zgemm_itcopy)(m, n, a, lda, b);
}
pub fn zgemm_oncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zgemm_oncopy)(m, n, a, lda, b);
}
pub fn zgemm_otcopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zgemm_otcopy)(m, n, a, lda, b);
}
pub fn zgemm_small_matrix_permit(transa: i32, transb: i32, m: BlasLong, n: BlasLong, k: BlasLong, alpha0: f64, alpha1: f64, beta0: f64, beta1: f64) i32 {
    return blas.kernel(.zgemm_small_matrix_permit)(transa, transb, m, n, k, alpha0, alpha1, beta0, beta1);
}
pub fn zgemm_small_kernel_nn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_nn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_nt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_nt)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_nr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_nr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_nc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_nc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_tn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_tn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_tt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_tt)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_tr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_tr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_tc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_tc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_rn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_rn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_rt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_rt)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_rr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_rr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_rc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_rc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_cn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_cn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_ct(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_ct)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_cr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_cr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_cc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, beta0: f64, beta1: f64, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_cc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, beta0, beta1, C, ldc);
}
pub fn zgemm_small_kernel_b0_nn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_nn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_nt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_nt)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_nr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_nr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_nc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_nc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_tn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_tn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_tt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_tt)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_tr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_tr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_tc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_tc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_rn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_rn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_rt(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_rt)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_rr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_rr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_rc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_rc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_cn(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_cn)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_ct(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_ct)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_cr(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_cr)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn zgemm_small_kernel_b0_cc(m: BlasLong, n: BlasLong, k: BlasLong, A: [*]f64, lda: BlasLong, alpha0: f64, alpha1: f64, B: [*]f64, ldb: BlasLong, C: [*]f64, ldc: BlasLong) i32 {
    return blas.kernel(.zgemm_small_kernel_b0_cc)(m, n, k, A, lda, alpha0, alpha1, B, ldb, C, ldc);
}
pub fn xgemm_kernel_n(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong) i32 {
    return blas.kernel(.xgemm_kernel_n)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn xgemm_kernel_l(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong) i32 {
    return blas.kernel(.xgemm_kernel_l)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn xgemm_kernel_r(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong) i32 {
    return blas.kernel(.xgemm_kernel_r)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn xgemm_kernel_b(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong) i32 {
    return blas.kernel(.xgemm_kernel_b)(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc);
}
pub fn xgemm_beta(m: BlasLong, n: BlasLong, dummy1: BlasLong, beta_r: XDouble, beta_i: XDouble, dummy2: [*]XDouble, dummy3: BlasLong, dummy4: [*]XDouble, dummy5: BlasLong, c: [*]XDouble, ldc: BlasLong) i32 {
    return blas.kernel(.xgemm_beta)(m, n, dummy1, beta_r, beta_i, dummy2, dummy3, dummy4, dummy5, c, ldc);
}
pub fn xgemm_incopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm_incopy)(m, n, a, lda, b);
}
pub fn xgemm_itcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm_itcopy)(m, n, a, lda, b);
}
pub fn xgemm_oncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm_oncopy)(m, n, a, lda, b);
}
pub fn xgemm_otcopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xgemm_otcopy)(m, n, a, lda, b);
}

pub const params = @import("params.zig");
pub const kernel = @import("kernel.zig");

/// Points the level-3 slots at the portable kernels (`kernel.zig`).
pub const installGeneric = kernel.installGeneric;

const check = @import("../../_test/blas3/check.zig");
const Vtable = @import("../../Vtable.zig");

fn cases(vt: *const Vtable) !void {
    try check.routineAll(vt, .gemm);
}

const std = @import("std");
const sysinfo = @import("builtin");
const driver = @import("../driver.zig");
const portable = @import("../../portable.zig");
const blas1check = @import("../../_test/blas1/check.zig");

test "portable gemm kernels match the wide reference through the driver" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try cases(&vt);
}

test "every linked core's gemm kernels match the wide reference, where this machine runs them" {
    try blas1check.eachCore(cases);
}

test "every routine over tiles other targets choose: NEON's, 8 x 3 complex" {
    // natively only: the rigs run their own tiles, and the freestanding ones
    // have no room for four more sets of kernels
    if (comptime sysinfo.target.os.tag == .freestanding) return error.SkipZigTest;
    if (blas1check.rig()) return error.SkipZigTest;
    var vt: Vtable = .{};
    portable.install(&vt);
    kernel.installTile(&vt, "s", f32, .{ .vl = 4, .mv = 4, .nr = 6 });
    kernel.installTile(&vt, "d", f64, .{ .vl = 2, .mv = 4, .nr = 6 });
    kernel.installTile(&vt, "c", std.math.Complex(f32), .{ .vl = 2, .mv = 4, .nr = 3 });
    kernel.installTile(&vt, "z", std.math.Complex(f64), .{ .vl = 1, .mv = 4, .nr = 3 });
    try check.coreCases(&vt);
}

test "the unrolls are the installed kernels' own after a core's blocking" {
    const K = kernel.Kernels(f32, kernel.geometry(f32));
    var vt: Vtable = .{};
    portable.install(&vt);
    // a core whose object is not linked keeps the portable kernels
    params.installCore(&vt, .x86_64, .skylakex, .{});
    portable.finish(&vt);
    try std.testing.expectEqual(K.unroll_m, vt.sgemm_unroll_m);
    try std.testing.expectEqual(K.unroll_n, vt.sgemm_unroll_n);
    try std.testing.expectEqual(@max(K.unroll_m, K.unroll_n), vt.sgemm_unroll_mn);
}

test "negative or zero blocking falls back to the generic sizes" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    vt.dgemm_p = -7;
    vt.dgemm_q = 0;
    vt.dgemm_r = -100;
    const blk = driver.Blocking.of(f64, &vt);
    try std.testing.expect(blk.p > 0 and blk.q == 128 and blk.r == 4096);
    try check.routineCases(f64, &vt, .gemm);
}

test "a packing buffer on the stack is poisoned in tests: a copy that writes nothing shows" {
    const Nop = struct {
        fn copy(_: BlasLong, _: BlasLong, _: [*]f64, _: BlasLong, _: [*]f64) callconv(.c) i32 {
            return 0;
        }
    };
    var vt: Vtable = .{};
    installGeneric(&vt);
    vt.dgemm_itcopy = &Nop.copy;
    var a: [4]f64 = @splat(1);
    var c: [4]f64 = @splat(0);
    try driver.gemm(std.testing.failing_allocator, f64, &vt, false, false, .n, 2, 2, 2, 1, &a, 2, &a, 2, 0, &c, 2);
    for (c) |v| try std.testing.expect(std.math.isNan(v));
}

test {
    _ = kernel;
    _ = params;
}
