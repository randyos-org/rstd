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

const std = @import("std");
const check = @import("../../_test/blas3/check.zig");
const Vtable = @import("../../Vtable.zig");
const driver = @import("../driver.zig");
const iface = @import("../interface.zig");
const triangular = @import("../../blas2/triangular.zig");
const corespec = @import("../../cpu/_corespec.zig");
const registry = @import("../../asm_registry.zig");
const util = @import("../../util.zig");

fn cases(vt: *const Vtable) !void {
    try check.routineAll(vt, .trsm);
}

const blas1check = @import("../../_test/blas1/check.zig");
const portable = @import("../../portable.zig");

test "portable trsm matches the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's trsm matches the wide reference, where this machine runs it" {
    try blas1check.eachCore(cases);
}

const testing = std.testing;
const precisions = .{ f32, f64, Complex32, Complex64 };

/// A small integer, complex or not, as `[re, im]`.
const Int = [2]i64;

fn imul(x: Int, y: Int) Int {
    return .{ x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0] };
}

/// `t` times `v`, exact for a small `v`.
fn scaled(comptime T: type, v: Int, t: types.Float(T)) T {
    const re = @as(types.Float(T), @floatFromInt(v[0])) * t;
    const im = @as(types.Float(T), @floatFromInt(v[1])) * t;
    return if (comptime types.isComplex(T)) .init(re, im) else re;
}

/// `trsm` where `1 / A(i, i)` overflows: `A` is small integers times `t`, a
/// denormal whose reciprocal is infinite, the solution `X` small integers,
/// and `B` is `op(A) X` or `X op(A)`, so every step is exact and the wide
/// reference is `X` itself. Every side, triangle and transpose.
fn tinyDiagonal(comptime T: type, vt: *const Vtable) !void {
    const complex = comptime types.isComplex(T);
    const F = types.Float(T);
    const t = std.math.ldexp(@as(F, 1), std.math.floatExponentMin(F) - 4);
    try testing.expect(std.math.isInf(1 / t));
    const m = 9;
    const n = 7;
    const nan = std.math.nan(F);
    const ts: []const check.Trans = if (complex) &.{ .n, .t, .r, .c } else &.{ .n, .t };
    for ([_]bool{ true, false }) |left| for ([_]bool{ true, false }) |upper| for (ts) |tr| {
        const ka: usize = if (left) m else n;
        // `A` as integers, then stored scaled, NaN where it is not stored
        var ai: [m * m]Int = undefined;
        var a: [m * m]T = undefined;
        for (0..ka) |c| for (0..ka) |r| {
            const v: Int = if (r == c)
                (if (r % 2 == 0) (if (complex) .{ 1, 1 } else .{ 1, 0 }) else (if (complex) .{ 2, -1 } else .{ -2, 0 }))
            else
                .{ @as(i64, @intCast((r + 2 * c) % 3)) - 1, if (complex) @intCast((r + c) % 2) else 0 };
            const kept = if (upper) r <= c else r >= c;
            ai[r + c * ka] = if (kept) v else .{ 0, 0 };
            a[r + c * ka] = if (kept) scaled(T, v, t) else if (complex) .init(nan, nan) else nan;
        };
        const opA = struct {
            fn f(p: []const Int, k: usize, trc: check.Trans, i: usize, l: usize) Int {
                const v = if (trc.transposed()) p[l + i * k] else p[i + l * k];
                return if (trc.conjugated()) .{ v[0], -v[1] } else v;
            }
        }.f;
        var x: [m * n]Int = undefined;
        for (0..n) |j| for (0..m) |i| {
            x[i + j * m] = .{ @as(i64, @intCast((3 * i + j) % 5)) - 2, if (complex) @as(i64, @intCast((i + 2 * j) % 3)) - 1 else 0 };
        };
        var b: [m * n]T = undefined;
        for (0..n) |j| for (0..m) |i| {
            var sum: Int = .{ 0, 0 };
            for (0..ka) |l| {
                const p = if (left) imul(opA(&ai, ka, tr, i, l), x[l + j * m]) else imul(x[i + l * m], opA(&ai, ka, tr, l, j));
                sum = .{ sum[0] + p[0], sum[1] + p[1] };
            }
            b[i + j * m] = scaled(T, sum, t);
        };
        const one: T = if (complex) .init(1, 0) else 1;
        try driver.trsm(testing.allocator, check.io, T, vt, left, upper, tr.transposed(), tr.conjugated(), false, m, n, one, &a, ka, &b, m);
        for (0..n) |j| for (0..m) |i| {
            const want = scaled(T, x[i + j * m], 1);
            const got = b[i + j * m];
            const ok = if (complex) got.re == want.re and got.im == want.im else got == want;
            if (!ok) {
                std.debug.print("{s} left={} upper={} trans={t} at ({d},{d}): got {any}, want {any}\n", .{ @typeName(T), left, upper, tr, i, j, got, want });
                return error.TestUnexpectedResult;
            }
        };
    };
}

fn tinyCases(vt: *const Vtable) !void {
    inline for (precisions) |T| try tinyDiagonal(T, vt);
}

test "trsm divides by a diagonal whose reciprocal overflows, and stays finite and exact" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try tinyCases(&vt);
}

test "every linked core's trsm divides by a tiny diagonal, where this machine runs it" {
    try blas1check.eachCore(tinyCases);
}

/// `trsv` and a one-column `trsm` on the same system, every option, through
/// the interfaces: the same bits. `A` is diagonal, its other stored elements
/// zero, so only the divisions round (the updates' orders differ); every
/// third diagonal element is so small its reciprocal overflows.
fn sameAsTrsv(comptime T: type) !void {
    const complex = comptime types.isComplex(T);
    const F = types.Float(T);
    if (!sameUnit(F)) return;
    const n = 37;
    const t = std.math.ldexp(@as(F, 1), std.math.floatExponentMin(F) - 4);
    const zero: T = if (complex) .init(0, 0) else 0;
    const one: T = if (complex) .init(1, 0) else 1;
    const nan: T = if (complex) .init(std.math.nan(F), std.math.nan(F)) else std.math.nan(F);
    var d: [n]T = undefined;
    var rhs: [n]T = undefined;
    blas1check.fill(T, &d, 5);
    blas1check.fill(T, &rhs, 6);
    for (0..n) |i| {
        // a third tiny, all at least 1 in magnitude before scaling
        const s: F = if (i % 3 == 0) t else 1;
        if (complex) {
            d[i] = .init((d[i].re + std.math.copysign(@as(F, 1), d[i].re)) * s, d[i].im * s);
            rhs[i] = .init(rhs[i].re * s, rhs[i].im * s);
        } else {
            d[i] = (d[i] + std.math.copysign(@as(F, 1), d[i])) * s;
            rhs[i] *= s;
        }
    }
    const ts: []const types.Trans = if (complex) &.{ .n, .t, .r, .c } else &.{ .n, .t };
    var a: [n * n]T = undefined;
    for ([_]bool{ true, false }) |upper| for (ts) |tr| for ([_]bool{ false, true }) |unit| {
        for (0..n) |c| for (0..n) |r| {
            a[r + c * n] = if (r == c) d[r] else if (if (upper) r < c else r > c) zero else nan;
        };
        var x = rhs;
        var y = rhs;
        const uplo: types.Uplo = if (upper) .u else .l;
        const diag: types.Diag = if (unit) .u else .n;
        try triangular.trsv(testing.allocator, check.io, T, uplo, tr, diag, n, &a, n, &x, 1);
        try iface.trsm(testing.allocator, check.io, T, .l, uplo, tr, diag, n, 1, one, &a, n, &y, n);
        for (x) |u| {
            const fin = if (complex) std.math.isFinite(u.re) and std.math.isFinite(u.im) else std.math.isFinite(u);
            try testing.expect(fin);
        }
        try testing.expectEqualSlices(u8, std.mem.asBytes(&x), std.mem.asBytes(&y));
    };
}

/// Whether the selected core's `trsm` computes `F` where `trsv` does. An
/// x87 build's `trsv` divides on the x87, and the objects of its cores with
/// SSE2 in SSE2 (`build/blaskern.zig`); every feature an x86 spec lists
/// implies SSE2. An unlinked core runs the target's own object.
fn sameUnit(comptime F: type) bool {
    if (comptime !util.onX87(F)) return true;
    const core = blas.coreType() orelse return true;
    if (!registry.isLinked(core)) return true;
    const spec = corespec.of(core);
    return spec.model == null and spec.features.len == 0;
}

test "trsv and trsm give the same bits on the same one-column system" {
    inline for (precisions) |T| try sameAsTrsv(T);
}

/// `trsm` with a unit diagonal, which skips the division, against the same
/// solve dividing by a stored diagonal of ones: the same bits, every side,
/// triangle and transpose. The unit run's diagonal is NaN, so never read.
fn unitCases(vt: *const Vtable) !void {
    inline for (precisions) |T| {
        const complex = comptime types.isComplex(T);
        const F = types.Float(T);
        const m = 23;
        const n = 17;
        const nan: T = if (complex) .init(std.math.nan(F), std.math.nan(F)) else std.math.nan(F);
        const one: T = if (complex) .init(1, 0) else 1;
        var a: [m * m]T = undefined;
        var rhs: [m * n]T = undefined;
        blas1check.fill(T, &a, 7);
        blas1check.fill(T, &rhs, 8);
        const ts: []const check.Trans = if (complex) &.{ .n, .t, .r, .c } else &.{ .n, .t };
        for ([_]bool{ true, false }) |left| for ([_]bool{ true, false }) |upper| for (ts) |tr| {
            const ka: usize = if (left) m else n;
            var au = a;
            var an = a;
            for (0..ka) |i| {
                au[i + i * ka] = nan;
                an[i + i * ka] = one;
            }
            var bu = rhs;
            var bn = rhs;
            try driver.trsm(testing.allocator, check.io, T, vt, left, upper, tr.transposed(), tr.conjugated(), true, m, n, one, &au, ka, &bu, m);
            try driver.trsm(testing.allocator, check.io, T, vt, left, upper, tr.transposed(), tr.conjugated(), false, m, n, one, &an, ka, &bn, m);
            testing.expectEqualSlices(u8, std.mem.asBytes(&bn), std.mem.asBytes(&bu)) catch |e| {
                std.debug.print("{s} left={} upper={} trans={t}\n", .{ @typeName(T), left, upper, tr });
                return e;
            };
        };
    }
}

test "trsm on a unit diagonal gives the bits of dividing by ones" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try unitCases(&vt);
}

test "every linked core's trsm on a unit diagonal gives the bits of dividing by ones" {
    try blas1check.eachCore(unitCases);
}
