//! Vtable wrappers for the `imatcopy` kernels.
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

pub fn simatcopy_k_cn(rows: BlasLong, cols: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong) i32 {
    return blas._vtable.simatcopy_k_cn.?(rows, cols, alpha, a, lda);
}
pub fn simatcopy_k_ct(rows: BlasLong, cols: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong) i32 {
    return blas._vtable.simatcopy_k_ct.?(rows, cols, alpha, a, lda);
}
pub fn simatcopy_k_rn(rows: BlasLong, cols: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong) i32 {
    return blas._vtable.simatcopy_k_rn.?(rows, cols, alpha, a, lda);
}
pub fn simatcopy_k_rt(rows: BlasLong, cols: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong) i32 {
    return blas._vtable.simatcopy_k_rt.?(rows, cols, alpha, a, lda);
}
pub fn dimatcopy_k_cn(rows: BlasLong, cols: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong) i32 {
    return blas._vtable.dimatcopy_k_cn.?(rows, cols, alpha, a, lda);
}
pub fn dimatcopy_k_ct(rows: BlasLong, cols: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong) i32 {
    return blas._vtable.dimatcopy_k_ct.?(rows, cols, alpha, a, lda);
}
pub fn dimatcopy_k_rn(rows: BlasLong, cols: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong) i32 {
    return blas._vtable.dimatcopy_k_rn.?(rows, cols, alpha, a, lda);
}
pub fn dimatcopy_k_rt(rows: BlasLong, cols: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong) i32 {
    return blas._vtable.dimatcopy_k_rt.?(rows, cols, alpha, a, lda);
}
pub fn cimatcopy_k_cn(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong) i32 {
    return blas._vtable.cimatcopy_k_cn.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn cimatcopy_k_ct(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong) i32 {
    return blas._vtable.cimatcopy_k_ct.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn cimatcopy_k_rn(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong) i32 {
    return blas._vtable.cimatcopy_k_rn.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn cimatcopy_k_rt(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong) i32 {
    return blas._vtable.cimatcopy_k_rt.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn cimatcopy_k_cnc(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong) i32 {
    return blas._vtable.cimatcopy_k_cnc.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn cimatcopy_k_ctc(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong) i32 {
    return blas._vtable.cimatcopy_k_ctc.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn cimatcopy_k_rnc(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong) i32 {
    return blas._vtable.cimatcopy_k_rnc.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn cimatcopy_k_rtc(rows: BlasLong, cols: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong) i32 {
    return blas._vtable.cimatcopy_k_rtc.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn zimatcopy_k_cn(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong) i32 {
    return blas._vtable.zimatcopy_k_cn.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn zimatcopy_k_ct(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong) i32 {
    return blas._vtable.zimatcopy_k_ct.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn zimatcopy_k_rn(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong) i32 {
    return blas._vtable.zimatcopy_k_rn.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn zimatcopy_k_rt(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong) i32 {
    return blas._vtable.zimatcopy_k_rt.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn zimatcopy_k_cnc(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong) i32 {
    return blas._vtable.zimatcopy_k_cnc.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn zimatcopy_k_ctc(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong) i32 {
    return blas._vtable.zimatcopy_k_ctc.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn zimatcopy_k_rnc(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong) i32 {
    return blas._vtable.zimatcopy_k_rnc.?(rows, cols, alpha_r, alpha_i, a, lda);
}
pub fn zimatcopy_k_rtc(rows: BlasLong, cols: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong) i32 {
    return blas._vtable.zimatcopy_k_rtc.?(rows, cols, alpha_r, alpha_i, a, lda);
}
