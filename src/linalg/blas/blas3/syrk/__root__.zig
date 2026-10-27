//! Vtable wrappers for the `syrk` kernels.
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

pub fn ssyrk_direct_alpha_betaUN(N: BlasLong, K: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, beta: f32, C: [*]f32, strideC: BlasLong) void {
    blas.kernel(.ssyrk_direct_alpha_betaUN)(N, K, alpha, A, strideA, beta, C, strideC);
}
pub fn ssyrk_direct_alpha_betaUT(N: BlasLong, K: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, beta: f32, C: [*]f32, strideC: BlasLong) void {
    blas.kernel(.ssyrk_direct_alpha_betaUT)(N, K, alpha, A, strideA, beta, C, strideC);
}
pub fn ssyrk_direct_alpha_betaLN(N: BlasLong, K: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, beta: f32, C: [*]f32, strideC: BlasLong) void {
    blas.kernel(.ssyrk_direct_alpha_betaLN)(N, K, alpha, A, strideA, beta, C, strideC);
}
pub fn ssyrk_direct_alpha_betaLT(N: BlasLong, K: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, beta: f32, C: [*]f32, strideC: BlasLong) void {
    blas.kernel(.ssyrk_direct_alpha_betaLT)(N, K, alpha, A, strideA, beta, C, strideC);
}

const check = @import("../../_test/blas3/check.zig");
const Vtable = @import("../../Vtable.zig");

fn cases(vt: *const Vtable) !void {
    try check.routineAll(vt, .syrk);
    try check.routineAll(vt, .herk);
}

const blas1check = @import("../../_test/blas1/check.zig");
const portable = @import("../../portable.zig");

test "portable syrk and herk match the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's syrk and herk match the wide reference, where this machine runs them" {
    try blas1check.eachCore(cases);
}
