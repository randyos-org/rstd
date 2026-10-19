//! Vtable wrappers for the `trsm` kernels.
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

pub fn sbtrsm_kernel_LN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.sbtrsm_kernel_LN)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn sbtrsm_kernel_LT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.sbtrsm_kernel_LT)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn sbtrsm_kernel_RN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.sbtrsm_kernel_RN)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn sbtrsm_kernel_RT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.sbtrsm_kernel_RT)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn sbtrsm_iunucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_iunucopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_iunncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_iunncopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_iutucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_iutucopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_iutncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_iutncopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_ilnucopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_ilnncopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_iltucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_iltucopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_iltncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_iltncopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_ounucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_ounucopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_ounncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_ounncopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_outucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_outucopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_outncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_outncopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_olnucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_olnucopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_olnncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_olnncopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_oltucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_oltucopy)(m, n, a, lda, offset, b);
}
pub fn sbtrsm_oltncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.sbtrsm_oltncopy)(m, n, a, lda, offset, b);
}
pub fn strsm_kernel_LN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.strsm_kernel_LN)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn strsm_kernel_LT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.strsm_kernel_LT)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn strsm_kernel_RN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.strsm_kernel_RN)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn strsm_kernel_RT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.strsm_kernel_RT)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn strsm_iunucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_iunucopy)(m, n, a, lda, offset, b);
}
pub fn strsm_iunncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_iunncopy)(m, n, a, lda, offset, b);
}
pub fn strsm_iutucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_iutucopy)(m, n, a, lda, offset, b);
}
pub fn strsm_iutncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_iutncopy)(m, n, a, lda, offset, b);
}
pub fn strsm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_ilnucopy)(m, n, a, lda, offset, b);
}
pub fn strsm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_ilnncopy)(m, n, a, lda, offset, b);
}
pub fn strsm_iltucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_iltucopy)(m, n, a, lda, offset, b);
}
pub fn strsm_iltncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_iltncopy)(m, n, a, lda, offset, b);
}
pub fn strsm_ounucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_ounucopy)(m, n, a, lda, offset, b);
}
pub fn strsm_ounncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_ounncopy)(m, n, a, lda, offset, b);
}
pub fn strsm_outucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_outucopy)(m, n, a, lda, offset, b);
}
pub fn strsm_outncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_outncopy)(m, n, a, lda, offset, b);
}
pub fn strsm_olnucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_olnucopy)(m, n, a, lda, offset, b);
}
pub fn strsm_olnncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_olnncopy)(m, n, a, lda, offset, b);
}
pub fn strsm_oltucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_oltucopy)(m, n, a, lda, offset, b);
}
pub fn strsm_oltncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.strsm_oltncopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_kernel_LN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.dtrsm_kernel_LN)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn dtrsm_kernel_LT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.dtrsm_kernel_LT)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn dtrsm_kernel_RN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.dtrsm_kernel_RN)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn dtrsm_kernel_RT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.dtrsm_kernel_RT)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn dtrsm_iunucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_iunucopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_iunncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_iunncopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_iutucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_iutucopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_iutncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_iutncopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_ilnucopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_ilnncopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_iltucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_iltucopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_iltncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_iltncopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_ounucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_ounucopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_ounncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_ounncopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_outucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_outucopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_outncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_outncopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_olnucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_olnucopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_olnncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_olnncopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_oltucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_oltucopy)(m, n, a, lda, offset, b);
}
pub fn dtrsm_oltncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.dtrsm_oltncopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_kernel_LN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: XDouble, a: [*]XDouble, b: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.qtrsm_kernel_LN)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn qtrsm_kernel_LT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: XDouble, a: [*]XDouble, b: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.qtrsm_kernel_LT)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn qtrsm_kernel_RN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: XDouble, a: [*]XDouble, b: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.qtrsm_kernel_RN)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn qtrsm_kernel_RT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: XDouble, a: [*]XDouble, b: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.qtrsm_kernel_RT)(m, n, k, dummy1, a, b, c, ldc, offset);
}
pub fn qtrsm_iunucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_iunucopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_iunncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_iunncopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_iutucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_iutucopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_iutncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_iutncopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_ilnucopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_ilnncopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_iltucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_iltucopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_iltncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_iltncopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_ounucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_ounucopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_ounncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_ounncopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_outucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_outucopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_outncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_outncopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_olnucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_olnucopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_olnncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_olnncopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_oltucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_oltucopy)(m, n, a, lda, offset, b);
}
pub fn qtrsm_oltncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.qtrsm_oltncopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_kernel_LN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, dummy2: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ctrsm_kernel_LN)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ctrsm_kernel_LT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, dummy2: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ctrsm_kernel_LT)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ctrsm_kernel_LR(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, dummy2: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ctrsm_kernel_LR)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ctrsm_kernel_LC(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, dummy2: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ctrsm_kernel_LC)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ctrsm_kernel_RN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, dummy2: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ctrsm_kernel_RN)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ctrsm_kernel_RT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, dummy2: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ctrsm_kernel_RT)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ctrsm_kernel_RR(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, dummy2: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ctrsm_kernel_RR)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ctrsm_kernel_RC(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f32, dummy2: f32, a: [*]f32, b: [*]f32, c: [*]f32, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ctrsm_kernel_RC)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ctrsm_iunucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_iunucopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_iunncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_iunncopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_iutucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_iutucopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_iutncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_iutncopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_ilnucopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_ilnncopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_iltucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_iltucopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_iltncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_iltncopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_ounucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_ounucopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_ounncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_ounncopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_outucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_outucopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_outncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_outncopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_olnucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_olnucopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_olnncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_olnncopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_oltucopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_oltucopy)(m, n, a, lda, offset, b);
}
pub fn ctrsm_oltncopy(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, offset: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.ctrsm_oltncopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_kernel_LN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f64, dummy2: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ztrsm_kernel_LN)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ztrsm_kernel_LT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f64, dummy2: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ztrsm_kernel_LT)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ztrsm_kernel_LR(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f64, dummy2: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ztrsm_kernel_LR)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ztrsm_kernel_LC(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f64, dummy2: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ztrsm_kernel_LC)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ztrsm_kernel_RN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f64, dummy2: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ztrsm_kernel_RN)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ztrsm_kernel_RT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f64, dummy2: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ztrsm_kernel_RT)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ztrsm_kernel_RR(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f64, dummy2: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ztrsm_kernel_RR)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ztrsm_kernel_RC(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: f64, dummy2: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.ztrsm_kernel_RC)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn ztrsm_iunucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_iunucopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_iunncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_iunncopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_iutucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_iutucopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_iutncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_iutncopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_ilnucopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_ilnncopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_iltucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_iltucopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_iltncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_iltncopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_ounucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_ounucopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_ounncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_ounncopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_outucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_outucopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_outncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_outncopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_olnucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_olnucopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_olnncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_olnncopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_oltucopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_oltucopy)(m, n, a, lda, offset, b);
}
pub fn ztrsm_oltncopy(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, offset: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.ztrsm_oltncopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_kernel_LN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: XDouble, dummy2: XDouble, a: [*]XDouble, b: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.xtrsm_kernel_LN)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn xtrsm_kernel_LT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: XDouble, dummy2: XDouble, a: [*]XDouble, b: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.xtrsm_kernel_LT)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn xtrsm_kernel_LR(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: XDouble, dummy2: XDouble, a: [*]XDouble, b: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.xtrsm_kernel_LR)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn xtrsm_kernel_LC(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: XDouble, dummy2: XDouble, a: [*]XDouble, b: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.xtrsm_kernel_LC)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn xtrsm_kernel_RN(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: XDouble, dummy2: XDouble, a: [*]XDouble, b: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.xtrsm_kernel_RN)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn xtrsm_kernel_RT(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: XDouble, dummy2: XDouble, a: [*]XDouble, b: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.xtrsm_kernel_RT)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn xtrsm_kernel_RR(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: XDouble, dummy2: XDouble, a: [*]XDouble, b: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.xtrsm_kernel_RR)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn xtrsm_kernel_RC(m: BlasLong, n: BlasLong, k: BlasLong, dummy1: XDouble, dummy2: XDouble, a: [*]XDouble, b: [*]XDouble, c: [*]XDouble, ldc: BlasLong, offset: BlasLong) i32 {
    return blas.kernel(.xtrsm_kernel_RC)(m, n, k, dummy1, dummy2, a, b, c, ldc, offset);
}
pub fn xtrsm_iunucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_iunucopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_iunncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_iunncopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_iutucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_iutucopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_iutncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_iutncopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_ilnucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_ilnucopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_ilnncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_ilnncopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_iltucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_iltucopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_iltncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_iltncopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_ounucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_ounucopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_ounncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_ounncopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_outucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_outucopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_outncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_outncopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_olnucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_olnucopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_olnncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_olnncopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_oltucopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_oltucopy)(m, n, a, lda, offset, b);
}
pub fn xtrsm_oltncopy(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, offset: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xtrsm_oltncopy)(m, n, a, lda, offset, b);
}

const check = @import("../check.zig");
const blas1check = @import("../../blas1/check.zig");
const Vtable = @import("../../Vtable.zig");

fn cases(vt: *const Vtable) !void {
    try check.routineAll(vt, .trsm);
}

test "portable trsm matches the wide reference through the driver" {
    var vt: Vtable = .{};
    @import("../../portable.zig").install(&vt);
    try cases(&vt);
}

test "every linked core's trsm matches the wide reference, where this machine runs it" {
    try blas1check.eachCore(cases);
}
