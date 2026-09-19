//! Vtable wrappers for the `omatcopy` kernels.
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

pub fn somatcopy_k_cn(rows: BlasLong, cols: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas._vtable.somatcopy_k_cn.?(rows, cols, alpha, a, lda, b, ldb);
}
pub fn somatcopy_k_ct(rows: BlasLong, cols: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas._vtable.somatcopy_k_ct.?(rows, cols, alpha, a, lda, b, ldb);
}
pub fn somatcopy_k_rn(rows: BlasLong, cols: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas._vtable.somatcopy_k_rn.?(rows, cols, alpha, a, lda, b, ldb);
}
pub fn somatcopy_k_rt(rows: BlasLong, cols: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas._vtable.somatcopy_k_rt.?(rows, cols, alpha, a, lda, b, ldb);
}
pub fn domatcopy_k_cn(rows: BlasLong, cols: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas._vtable.domatcopy_k_cn.?(rows, cols, alpha, a, lda, b, ldb);
}
pub fn domatcopy_k_ct(rows: BlasLong, cols: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas._vtable.domatcopy_k_ct.?(rows, cols, alpha, a, lda, b, ldb);
}
pub fn domatcopy_k_rn(rows: BlasLong, cols: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas._vtable.domatcopy_k_rn.?(rows, cols, alpha, a, lda, b, ldb);
}
pub fn domatcopy_k_rt(rows: BlasLong, cols: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas._vtable.domatcopy_k_rt.?(rows, cols, alpha, a, lda, b, ldb);
}
pub fn comatcopy_k_cn(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas._vtable.comatcopy_k_cn.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn comatcopy_k_ct(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas._vtable.comatcopy_k_ct.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn comatcopy_k_rn(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas._vtable.comatcopy_k_rn.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn comatcopy_k_rt(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas._vtable.comatcopy_k_rt.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn comatcopy_k_cnc(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas._vtable.comatcopy_k_cnc.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn comatcopy_k_ctc(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas._vtable.comatcopy_k_ctc.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn comatcopy_k_rnc(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas._vtable.comatcopy_k_rnc.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn comatcopy_k_rtc(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, b: [*]f32, ldb: BlasLong) i32 {
    return blas._vtable.comatcopy_k_rtc.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn zomatcopy_k_cn(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas._vtable.zomatcopy_k_cn.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn zomatcopy_k_ct(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas._vtable.zomatcopy_k_ct.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn zomatcopy_k_rn(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas._vtable.zomatcopy_k_rn.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn zomatcopy_k_rt(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas._vtable.zomatcopy_k_rt.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn zomatcopy_k_cnc(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas._vtable.zomatcopy_k_cnc.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn zomatcopy_k_ctc(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas._vtable.zomatcopy_k_ctc.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn zomatcopy_k_rnc(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas._vtable.zomatcopy_k_rnc.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
pub fn zomatcopy_k_rtc(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, b: [*]f64, ldb: BlasLong) i32 {
    return blas._vtable.zomatcopy_k_rtc.?(rows, cols, alpha_r, alpha_i, a, lda, b, ldb);
}
