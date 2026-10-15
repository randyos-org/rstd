//! Vtable wrappers for the `syr2k` kernels.
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

pub fn ssyr2k_direct_alpha_betaUN(N: BlasLong, K: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, B: [*]f32, strideB: BlasLong, beta: f32, R: [*]f32, strideR: BlasLong) void {
    blas.kernel(.ssyr2k_direct_alpha_betaUN)(N, K, alpha, A, strideA, B, strideB, beta, R, strideR);
}
pub fn ssyr2k_direct_alpha_betaUT(N: BlasLong, K: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, B: [*]f32, strideB: BlasLong, beta: f32, R: [*]f32, strideR: BlasLong) void {
    blas.kernel(.ssyr2k_direct_alpha_betaUT)(N, K, alpha, A, strideA, B, strideB, beta, R, strideR);
}
pub fn ssyr2k_direct_alpha_betaLN(N: BlasLong, K: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, B: [*]f32, strideB: BlasLong, beta: f32, R: [*]f32, strideR: BlasLong) void {
    blas.kernel(.ssyr2k_direct_alpha_betaLN)(N, K, alpha, A, strideA, B, strideB, beta, R, strideR);
}
pub fn ssyr2k_direct_alpha_betaLT(N: BlasLong, K: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, B: [*]f32, strideB: BlasLong, beta: f32, R: [*]f32, strideR: BlasLong) void {
    blas.kernel(.ssyr2k_direct_alpha_betaLT)(N, K, alpha, A, strideA, B, strideB, beta, R, strideR);
}
