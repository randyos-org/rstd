//! Vtable wrappers for the `trmm` kernels.
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

pub fn sbtrmm_kernel_RN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.sbtrmm_kernel_RN.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn sbtrmm_kernel_RT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.sbtrmm_kernel_RT.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn sbtrmm_kernel_LN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.sbtrmm_kernel_LN.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn sbtrmm_kernel_LT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.sbtrmm_kernel_LT.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn sbtrmm_iunucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_iunucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_iunncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_iunncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_iutucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_iutucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_iutncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_iutncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_ilnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_ilnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_iltucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_iltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_iltncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_iltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_ounucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_ounucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_ounncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_ounncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_outucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_outucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_outncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_outncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_olnucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_olnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_olnncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_olnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_oltucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_oltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn sbtrmm_oltncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.sbtrmm_oltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_direct_LNUN(M: BlasLong, N: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, B: [*]f32, strideB: BlasLong) void {
    blas._vtable.strmm_direct_LNUN.?(M, N, alpha, A, strideA, B, strideB);
}
pub fn strmm_direct_LNLN(M: BlasLong, N: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, B: [*]f32, strideB: BlasLong) void {
    blas._vtable.strmm_direct_LNLN.?(M, N, alpha, A, strideA, B, strideB);
}
pub fn strmm_direct_LTUN(M: BlasLong, N: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, B: [*]f32, strideB: BlasLong) void {
    blas._vtable.strmm_direct_LTUN.?(M, N, alpha, A, strideA, B, strideB);
}
pub fn strmm_direct_LTLN(M: BlasLong, N: BlasLong, alpha: f32, A: [*]f32, strideA: BlasLong, B: [*]f32, strideB: BlasLong) void {
    blas._vtable.strmm_direct_LTLN.?(M, N, alpha, A, strideA, B, strideB);
}
pub fn strmm_kernel_RN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.strmm_kernel_RN.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn strmm_kernel_RT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.strmm_kernel_RT.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn strmm_kernel_LN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.strmm_kernel_LN.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn strmm_kernel_LT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.strmm_kernel_LT.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn strmm_iunucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_iunucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_iunncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_iunncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_iutucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_iutucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_iutncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_iutncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_ilnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_ilnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_iltucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_iltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_iltncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_iltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_ounucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_ounucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_ounncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_ounncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_outucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_outucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_outncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_outncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_olnucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_olnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_olnncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_olnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_oltucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_oltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn strmm_oltncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.strmm_oltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_kernel_RN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.dtrmm_kernel_RN.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn dtrmm_kernel_RT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.dtrmm_kernel_RT.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn dtrmm_kernel_LN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.dtrmm_kernel_LN.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn dtrmm_kernel_LT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.dtrmm_kernel_LT.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn dtrmm_iunucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_iunucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_iunncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_iunncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_iutucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_iutucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_iutncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_iutncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_ilnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_ilnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_iltucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_iltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_iltncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_iltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_ounucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_ounucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_ounncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_ounncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_outucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_outucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_outncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_outncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_olnucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_olnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_olnncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_olnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_oltucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_oltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn dtrmm_oltncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.dtrmm_oltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_kernel_RN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.qtrmm_kernel_RN.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn qtrmm_kernel_RT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.qtrmm_kernel_RT.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn qtrmm_kernel_LN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.qtrmm_kernel_LN.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn qtrmm_kernel_LT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.qtrmm_kernel_LT.?(bm, bn, bk, alpha, ba, bb, c, ldc, offset);
}
pub fn qtrmm_iunucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_iunucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_iunncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_iunncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_iutucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_iutucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_iutncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_iutncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_ilnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_ilnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_iltucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_iltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_iltncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_iltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_ounucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_ounucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_ounncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_ounncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_outucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_outucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_outncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_outncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_olnucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_olnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_olnncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_olnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_oltucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_oltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn qtrmm_oltncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.qtrmm_oltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_kernel_RN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ctrmm_kernel_RN.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ctrmm_kernel_RT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ctrmm_kernel_RT.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ctrmm_kernel_RR(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ctrmm_kernel_RR.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ctrmm_kernel_RC(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ctrmm_kernel_RC.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ctrmm_kernel_LN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ctrmm_kernel_LN.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ctrmm_kernel_LT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ctrmm_kernel_LT.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ctrmm_kernel_LR(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ctrmm_kernel_LR.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ctrmm_kernel_LC(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f32, alpha_i: f32, ba: [*]f32, bb: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ctrmm_kernel_LC.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ctrmm_iunucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_iunucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_iunncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_iunncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_iutucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_iutucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_iutncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_iutncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_ilnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_ilnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_iltucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_iltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_iltncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_iltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_ounucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_ounucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_ounncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_ounncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_outucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_outucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_outncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_outncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_olnucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_olnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_olnncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_olnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_oltucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_oltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ctrmm_oltncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas._vtable.ctrmm_oltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_kernel_RN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ztrmm_kernel_RN.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ztrmm_kernel_RT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ztrmm_kernel_RT.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ztrmm_kernel_RR(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ztrmm_kernel_RR.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ztrmm_kernel_RC(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ztrmm_kernel_RC.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ztrmm_kernel_LN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ztrmm_kernel_LN.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ztrmm_kernel_LT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ztrmm_kernel_LT.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ztrmm_kernel_LR(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ztrmm_kernel_LR.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ztrmm_kernel_LC(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: f64, alpha_i: f64, ba: [*]f64, bb: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.ztrmm_kernel_LC.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn ztrmm_iunucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_iunucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_iunncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_iunncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_iutucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_iutucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_iutncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_iutncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_ilnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_ilnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_iltucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_iltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_iltncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_iltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_ounucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_ounucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_ounncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_ounncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_outucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_outucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_outncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_outncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_olnucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_olnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_olnncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_olnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_oltucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_oltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn ztrmm_oltncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas._vtable.ztrmm_oltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_kernel_RN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.xtrmm_kernel_RN.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn xtrmm_kernel_RT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.xtrmm_kernel_RT.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn xtrmm_kernel_RR(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.xtrmm_kernel_RR.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn xtrmm_kernel_RC(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.xtrmm_kernel_RC.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn xtrmm_kernel_LN(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.xtrmm_kernel_LN.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn xtrmm_kernel_LT(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.xtrmm_kernel_LT.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn xtrmm_kernel_LR(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.xtrmm_kernel_LR.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn xtrmm_kernel_LC(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: XDouble, alpha_i: XDouble, ba: [*]XDouble, bb: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas._vtable.xtrmm_kernel_LC.?(bm, bn, bk, alpha_r, alpha_i, ba, bb, c, ldc, offset);
}
pub fn xtrmm_iunucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_iunucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_iunncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_iunncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_iutucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_iutucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_iutncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_iutncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_ilnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_ilnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_iltucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_iltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_iltncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_iltncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_ounucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_ounucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_ounncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_ounncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_outucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_outucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_outncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_outncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_olnucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_olnucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_olnncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_olnncopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_oltucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_oltucopy.?(m, n, a, lda, posX, posY, b);
}
pub fn xtrmm_oltncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas._vtable.xtrmm_oltncopy.?(m, n, a, lda, posX, posY, b);
}
