//! Zig translation of OpenBLAS's `gotoblas_t` (the per-core kernel table),
//! mechanically derived from `common_param.h` in
//! https://github.com/OpenMathLib/OpenBLAS at commit
//! 9724481b599bb75fd260a9e7e991b445ef02a0bd (v0.3.34-127; `common_param.h` is
//! unchanged at v0.3.34-317).
//! Re-derive from that source rather than hand-editing if it goes stale.
//!
//! Every field is present as an optional pointer, null where it does not
//! apply, in place of the reference's `BUILD_*` and `ARCH_*` guards. Left out:
//! the row-major and CBLAS-only kernels (`?gerv_k`, `?gerd_k`, `*_direct*`,
//! `?omatcopy_k_r*`, `?imatcopy_k_r*`) and `sbger_k`, `sbsymv_L`, `sbsymv_U`,
//! which nothing calls (`kernel/setparam-ref.c:112-113`); they are in the
//! git-ignored `unused/`. Imports stay in `src/linalg/blas/`, for the
//! per-core objects; its test is in `__root__.zig`.

const types = @import("types.zig");

const BlasLong = types.BlasLong;
const HF16 = types.HF16;
const BF16 = types.BF16;
const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;
const CComplex32 = types.CComplex32;
const CComplex64 = types.CComplex64;
const CComplexXD = types.CComplexXD;

dtb_entries: i32 = 0,
/// `SYMV_P`, the block of `kernel/generic/symv_k.c` and `zhemv_k.c`. A compile-time
/// constant of each core's `param.h` upstream; 0 here is the generic target's 16.
symv_p: i32 = 0,
switch_ratio: i32 = 0,
divide_rate: i32 = 0,
divide_limit: i32 = 0,
preferred_size: i32 = 0,
offsetA: i32 = 0,
offsetB: i32 = 0,
@"align": i32 = 0,

//#if BUILD_HFLOAT16 == 1
shgemm_p: i32 = 0,
shgemm_q: i32 = 0,
shgemm_r: i32 = 0,
shgemm_unroll_m: i32 = 0,
shgemm_unroll_n: i32 = 0,
shgemm_unroll_mn: i32 = 0,

shgemm_kernel: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]HF16, [*]HF16, [*]f32, BlasLong) callconv(.c) i32 = null,
shgemm_beta: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]HF16, BlasLong, [*]HF16, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,

shgemm_incopy: ?*const fn (BlasLong, BlasLong, [*]HF16, BlasLong, [*]HF16) callconv(.c) i32 = null,
shgemm_itcopy: ?*const fn (BlasLong, BlasLong, [*]HF16, BlasLong, [*]HF16) callconv(.c) i32 = null,
shgemm_oncopy: ?*const fn (BlasLong, BlasLong, [*]HF16, BlasLong, [*]HF16) callconv(.c) i32 = null,
shgemm_otcopy: ?*const fn (BlasLong, BlasLong, [*]HF16, BlasLong, [*]HF16) callconv(.c) i32 = null,

shgemv_n: ?*const fn (BlasLong, BlasLong, f32, [*]HF16, BlasLong, [*]HF16, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
shgemv_t: ?*const fn (BlasLong, BlasLong, f32, [*]HF16, BlasLong, [*]HF16, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
//#endif //BUILD_HFLOAT16 == 1

//#if BUILD_BFLOAT16 == 1
bgemm_p: i32 = 0,
bgemm_q: i32 = 0,
bgemm_r: i32 = 0,
bgemm_unroll_m: i32 = 0,
bgemm_unroll_n: i32 = 0,
bgemm_unroll_mn: i32 = 0,
bgemm_align_k: i32 = 0,

sbgemm_p: i32 = 0,
sbgemm_q: i32 = 0,
sbgemm_r: i32 = 0,
sbgemm_unroll_m: i32 = 0,
sbgemm_unroll_n: i32 = 0,
sbgemm_unroll_mn: i32 = 0,
sbgemm_align_k: i32 = 0,
/// 0 default, 1 for device support amx.
need_amxtile_permission: i32 = 0,

sbstobf16_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]BF16, BlasLong) callconv(.c) void = null,
sbdtobf16_k: ?*const fn (BlasLong, [*]f64, BlasLong, [*]BF16, BlasLong) callconv(.c) void = null,
sbf16tos_k: ?*const fn (BlasLong, [*]BF16, BlasLong, [*]f32, BlasLong) callconv(.c) void = null,
dbf16tod_k: ?*const fn (BlasLong, [*]BF16, BlasLong, [*]f64, BlasLong) callconv(.c) void = null,

sbamax_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
sbamin_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
sbmax_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
sbmin_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
isbamax_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) BlasLong = null,
isbamin_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) BlasLong = null,
isbmax_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) BlasLong = null,
isbmin_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) BlasLong = null,

sbnrm2_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
sbasum_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
sbsum_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
sbcopy_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
sbdot_k: ?*const fn (BlasLong, [*]BF16, BlasLong, [*]BF16, BlasLong) callconv(.c) f32 = null,
dsbdot_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) f64 = null,

sbrot_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, f32, f32) callconv(.c) i32 = null,
sbrotm_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,

bscal_k: ?*const fn (BlasLong, BlasLong, BlasLong, BF16, [*]BF16, BlasLong, [*]BF16, BlasLong, [*]BF16, BlasLong) callconv(.c) i32 = null,
sbaxpy_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
sbscal_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
sbswap_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,

bgemv_n: ?*const fn (BlasLong, BlasLong, BF16, [*]BF16, BlasLong, [*]BF16, BlasLong, BF16, [*]BF16, BlasLong) callconv(.c) i32 = null,
bgemv_t: ?*const fn (BlasLong, BlasLong, BF16, [*]BF16, BlasLong, [*]BF16, BlasLong, BF16, [*]BF16, BlasLong) callconv(.c) i32 = null,

sbgemv_n: ?*const fn (BlasLong, BlasLong, f32, [*]BF16, BlasLong, [*]BF16, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
sbgemv_t: ?*const fn (BlasLong, BlasLong, f32, [*]BF16, BlasLong, [*]BF16, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,

bgemm_kernel: ?*const fn (BlasLong, BlasLong, BlasLong, BF16, [*]BF16, [*]BF16, [*]BF16, BlasLong) callconv(.c) i32 = null,
bgemm_beta: ?*const fn (BlasLong, BlasLong, BlasLong, BF16, [*]BF16, BlasLong, [*]BF16, BlasLong, [*]BF16, BlasLong) callconv(.c) i32 = null,

bgemm_incopy: ?*const fn (BlasLong, BlasLong, [*]BF16, BlasLong, [*]BF16) callconv(.c) i32 = null,
bgemm_itcopy: ?*const fn (BlasLong, BlasLong, [*]BF16, BlasLong, [*]BF16) callconv(.c) i32 = null,
bgemm_oncopy: ?*const fn (BlasLong, BlasLong, [*]BF16, BlasLong, [*]BF16) callconv(.c) i32 = null,
bgemm_otcopy: ?*const fn (BlasLong, BlasLong, [*]BF16, BlasLong, [*]BF16) callconv(.c) i32 = null,

sbgemm_kernel: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]BF16, [*]BF16, [*]f32, BlasLong) callconv(.c) i32 = null,
sbgemm_beta: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]BF16, BlasLong, [*]BF16, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,

sbgemm_incopy: ?*const fn (BlasLong, BlasLong, [*]BF16, BlasLong, [*]BF16) callconv(.c) i32 = null,
sbgemm_itcopy: ?*const fn (BlasLong, BlasLong, [*]BF16, BlasLong, [*]BF16) callconv(.c) i32 = null,
sbgemm_oncopy: ?*const fn (BlasLong, BlasLong, [*]BF16, BlasLong, [*]BF16) callconv(.c) i32 = null,
sbgemm_otcopy: ?*const fn (BlasLong, BlasLong, [*]BF16, BlasLong, [*]BF16) callconv(.c) i32 = null,

sbtrsm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
sbtrsm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
sbtrsm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
sbtrsm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,

sbtrsm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_outucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_outncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrsm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,

sbtrmm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
sbtrmm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
sbtrmm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
sbtrmm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,

sbtrmm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_outucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_outncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbtrmm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,

sbsymm_iutcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbsymm_iltcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbsymm_outcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
sbsymm_oltcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,

sbneg_tcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
sblaswp_ncopy: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, [*]isize, [*]f32) callconv(.c) i32 = null,

//#ifdef SMALL_MATRIX_OPT // (bfloat16)
sbgemm_small_matrix_permit: ?*const fn (i32, i32, BlasLong, BlasLong, BlasLong, f32, f32) callconv(.c) i32 = null,

sbgemm_small_kernel_nn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]BF16, BlasLong, f32, [*]BF16, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
sbgemm_small_kernel_nt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]BF16, BlasLong, f32, [*]BF16, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
sbgemm_small_kernel_tn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]BF16, BlasLong, f32, [*]BF16, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
sbgemm_small_kernel_tt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]BF16, BlasLong, f32, [*]BF16, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,

sbgemm_small_kernel_b0_nn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]BF16, BlasLong, f32, [*]BF16, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
sbgemm_small_kernel_b0_nt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]BF16, BlasLong, f32, [*]BF16, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
sbgemm_small_kernel_b0_tn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]BF16, BlasLong, f32, [*]BF16, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
sbgemm_small_kernel_b0_tt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]BF16, BlasLong, f32, [*]BF16, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
//#endif //def SMALL_MATRIX_OPT
//#endif //BUILD_BFLOAT16 == 1

// classic single/double/complex/complex16 GEMM tiling params
//#if (BUILD_SINGLE == 1) || (BUILD_DOUBLE == 1) || (BUILD_COMPLEX == 1) || (BUILD_COMPLEX16 == 1)
sgemm_p: i32 = 0,
sgemm_q: i32 = 0,
sgemm_r: i32 = 0,
sgemm_unroll_m: i32 = 0,
sgemm_unroll_n: i32 = 0,
sgemm_unroll_mn: i32 = 0,
//#endif //(BUILD_SINGLE == 1) || (BUILD_DOUBLE == 1) || (BUILD_COMPLEX == 1) || (BUILD_COMPLEX16 == 1)

exclusive_cache: i32 = 0,

// single precision (s)
//#if (BUILD_SINGLE == 1) || (BUILD_COMPLEX == 1)
samax_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
samin_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
smax_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
smin_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
//#endif //(BUILD_SINGLE == 1) || (BUILD_COMPLEX == 1)

//#if (BUILD_SINGLE==1) || (BUILD_DOUBLE ==1) || (BUILD_COMPLEX==1)
isamax_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) BlasLong = null,
//#endif //(BUILD_SINGLE==1) || (BUILD_DOUBLE ==1) || (BUILD_COMPLEX==1)
//#if (BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
isamin_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) BlasLong = null,
ismax_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) BlasLong = null,
ismin_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) BlasLong = null,
snrm2_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
sasum_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
//#endif //(BUILD_SINGLE==1) || (BUILD_COMPLEX==1)
//#if (BUILD_SINGLE==1)
ssum_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
//#endif //(BUILD_SINGLE==1)
//#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1
scopy_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
sdot_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
//double (*dsdot_k) (BLASLONG, float *, BLASLONG, float *, BLASLONG);

srot_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, f32, f32) callconv(.c) i32 = null,
srotm_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
//#endif //(BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1
//#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
saxpy_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
//#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1) || (BUILD_COMPLEX16==1)
sscal_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1) || (BUILD_COMPLEX16==1)
//#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
sswap_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,

sgemv_n: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
sgemv_t: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
//#endif //(BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
//#if (BUILD_SINGLE==1)
sger_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
ssymv_L: ?*const fn (BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
ssymv_U: ?*const fn (BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
//#endif //(BUILD_SINGLE==1)

//#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
//#ifdef ARCH_ARM64
sme_sgemm_kernel: ?*const fn ([*:0]u8, [*:0]u8, BlasLong, BlasLong, BlasLong, [*]f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, [*]f32, BlasLong) callconv(.c) void = null,
//#endif //def ARCH_ARM64

sgemm_kernel: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong) callconv(.c) i32 = null,
sgemm_beta: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,

sgemm_incopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
sgemm_itcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
sgemm_oncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
sgemm_otcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
//#endif //(BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
//#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
//#ifdef SMALL_MATRIX_OPT //(single)
sgemm_small_matrix_permit: ?*const fn (i32, i32, BlasLong, BlasLong, BlasLong, f32, f32) callconv(.c) i32 = null,

sgemm_small_kernel_nn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, [*]f32, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
sgemm_small_kernel_nt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, [*]f32, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
sgemm_small_kernel_tn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, [*]f32, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
sgemm_small_kernel_tt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, [*]f32, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,

sgemm_small_kernel_b0_nn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
sgemm_small_kernel_b0_nt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
sgemm_small_kernel_b0_tn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
sgemm_small_kernel_b0_tt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
//#endif //def SMALL_MATRIX_OPT

strsm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
strsm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
strsm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
strsm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,

strsm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_outucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_outncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strsm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
//#endif //(BUILD_SINGLE==1) || (BUILD_DOUBLE==1) || (BUILD_COMPLEX==1)
//#if (BUILD_SINGLE==1)
strmm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
strmm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
strmm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
strmm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,

strmm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_outucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_outncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
strmm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,

ssymm_iutcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ssymm_iltcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ssymm_outcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ssymm_oltcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,

sneg_tcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
slaswp_ncopy: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, [*]isize, [*]f32) callconv(.c) i32 = null,
//#endif //(BUILD_SINGLE==1)

// double precision (d)
//#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
dgemm_p: i32 = 0,
dgemm_q: i32 = 0,
dgemm_r: i32 = 0,
dgemm_unroll_m: i32 = 0,
dgemm_unroll_n: i32 = 0,
dgemm_unroll_mn: i32 = 0,
//#endif //(BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)

//#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
damax_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
damin_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
dmax_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
dmin_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
idamax_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) BlasLong = null,
idamin_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) BlasLong = null,
idmax_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) BlasLong = null,
idmin_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) BlasLong = null,

dnrm2_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
dasum_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
//#endif //(BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
//#if (BUILD_DOUBLE==1)
dsum_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
//#endif //(BUILD_DOUBLE==1)
//#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
dcopy_k: ?*const fn (BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
ddot_k: ?*const fn (BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
//#endif //(BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
//#if (BUILD_SINGLE==1) || (BUILD_DOUBLE==1)
dsdot_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) f64 = null,
//#endif //(BUILD_SINGLE==1) || (BUILD_DOUBLE==1)
//#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
drot_k: ?*const fn (BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, f64, f64) callconv(.c) i32 = null,
drotm_k: ?*const fn (BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
daxpy_k: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
dscal_k: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
dswap_k: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
dgemv_n: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
dgemv_t: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
//#endif //(BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
//#if (BUILD_DOUBLE==1)
dger_k: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,

dsymv_L: ?*const fn (BlasLong, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
dsymv_U: ?*const fn (BlasLong, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
//#endif //(BUILD_DOUBLE==1)
//#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
//#ifdef ARCH_ARM64
/// ARCH_ARM64 only
sme_dgemm_kernel: ?*const fn ([*:0]const u8, [*:0]const u8, BlasLong, BlasLong, BlasLong, [*]const f64, [*]const f64, BlasLong, [*]const f64, BlasLong, [*]const f64, [*]f64, BlasLong) callconv(.c) void = null,
//#endif //def ARCH_ARM64
dgemm_kernel: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, [*]f64, [*]f64, BlasLong) callconv(.c) i32 = null,
dgemm_beta: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,

dgemm_incopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
dgemm_itcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
dgemm_oncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
dgemm_otcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
//#endif //(BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
//#if (BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
//#ifdef SMALL_MATRIX_OPT //(double)
dgemm_small_matrix_permit: ?*const fn (i32, i32, BlasLong, BlasLong, BlasLong, f64, f64) callconv(.c) i32 = null,

dgemm_small_kernel_nn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, [*]f64, BlasLong, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
dgemm_small_kernel_nt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, [*]f64, BlasLong, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
dgemm_small_kernel_tn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, [*]f64, BlasLong, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
dgemm_small_kernel_tt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, [*]f64, BlasLong, f64, [*]f64, BlasLong) callconv(.c) i32 = null,

dgemm_small_kernel_b0_nn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
dgemm_small_kernel_b0_nt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
dgemm_small_kernel_b0_tn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
dgemm_small_kernel_b0_tt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
//#endif //def SMALL_MATRIX_OPT
//#endif //(BUILD_DOUBLE==1) || (BUILD_COMPLEX16==1)
//#if (BUILD_DOUBLE==1)
dtrsm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
dtrsm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
dtrsm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
dtrsm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,

dtrsm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_outucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_outncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrsm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,

dtrmm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
dtrmm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
dtrmm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
dtrmm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,

dtrmm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_outucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_outncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dtrmm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,

dsymm_iutcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dsymm_iltcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dsymm_outcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
dsymm_oltcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,

dneg_tcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
dlaswp_ncopy: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, [*]isize, [*]f64) callconv(.c) i32 = null,
//#endif //(BUILD_DOUBLE==1)

// EXPRECISION (q-family; XDouble is f64 unless types.exprecision or types.quad_precision is on).
// The q/x slots hold the portable Zig kernels on every core, never the x87 assembly:
// `xdotu_k`/`xdotc_k` return an extern struct, through a hidden pointer at f80, where
// the assembly returns `_Complex long double` in st(0)/st(1).
//#ifdef EXPRECISION
qgemm_p: i32 = 0,
qgemm_q: i32 = 0,
qgemm_r: i32 = 0,
qgemm_unroll_m: i32 = 0,
qgemm_unroll_n: i32 = 0,
qgemm_unroll_mn: i32 = 0,

qamax_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
qamin_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
qmax_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
qmin_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
iqamax_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) BlasLong = null,
iqamin_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) BlasLong = null,
iqmax_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) BlasLong = null,
iqmin_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) BlasLong = null,

qnrm2_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
qasum_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
qsum_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
qcopy_k: ?*const fn (BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) i32 = null,
qdot_k: ?*const fn (BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
qrot_k: ?*const fn (BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, XDouble, XDouble) callconv(.c) i32 = null,
qrotm_k: ?*const fn (BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,

qaxpy_k: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) i32 = null,
qscal_k: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) i32 = null,
qswap_k: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) i32 = null,

qgemv_n: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qgemv_t: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qger_k: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,

qsymv_L: ?*const fn (BlasLong, BlasLong, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qsymv_U: ?*const fn (BlasLong, BlasLong, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,

qgemm_kernel: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong) callconv(.c) i32 = null,
qgemm_beta: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) i32 = null,

qgemm_incopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qgemm_itcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qgemm_oncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qgemm_otcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,

qtrsm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
qtrsm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
qtrsm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
qtrsm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,

qtrsm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_outucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_outncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrsm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,

qtrmm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
qtrmm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
qtrmm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
qtrmm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,

qtrmm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_outucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_outncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qtrmm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,

qsymm_iutcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qsymm_iltcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qsymm_outcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qsymm_oltcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,

qneg_tcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
qlaswp_ncopy: ?*const fn (BlasLong, BlasLong, BlasLong, [*]XDouble, BlasLong, [*]isize, [*]XDouble) callconv(.c) i32 = null,
//#endif //def EXPRECISION

// single-precision complex (c)
//#if (BUILD_COMPLEX==1)
cgemm_p: i32 = 0,
cgemm_q: i32 = 0,
cgemm_r: i32 = 0,
cgemm_unroll_m: i32 = 0,
cgemm_unroll_n: i32 = 0,
cgemm_unroll_mn: i32 = 0,

camax_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
camin_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,

icamax_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) BlasLong = null,
icamin_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) BlasLong = null,

cnrm2_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
casum_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,
csum_k: ?*const fn (BlasLong, [*]f32, BlasLong) callconv(.c) f32 = null,

ccopy_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cdotu_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) CComplex32 = null,
cdotc_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) CComplex32 = null,

csrot_k: ?*const fn (BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, f32, f32) callconv(.c) i32 = null,
caxpy_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
caxpyc_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cscal_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cswap_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,

cgemv_n: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemv_t: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemv_r: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemv_c: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemv_o: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemv_u: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemv_s: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemv_d: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgeru_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgerc_k: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,

csymv_L: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
csymv_U: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
chemv_L: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
chemv_U: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
chemv_M: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
chemv_V: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,

//#ifdef ARCH_ARM64
/// ARCH_ARM64 only
sme_cgemm_kernel: ?*const fn ([*:0]const u8, [*:0]const u8, BlasLong, BlasLong, BlasLong, f32, f32, [*]const f32, BlasLong, [*]const f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) void = null,
//#endif //def ARCH_ARM64
cgemm_kernel_n: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_kernel_l: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_kernel_r: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_kernel_b: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_beta: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,

cgemm_incopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemm_itcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemm_oncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemm_otcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,

// #ifdef SMALL_MATRIX_OPT //(single complex)
cgemm_small_matrix_permit: ?*const fn (i32, i32, BlasLong, BlasLong, BlasLong, f32, f32, f32, f32) callconv(.c) i32 = null,

cgemm_small_kernel_nn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_nt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_nr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_nc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,

cgemm_small_kernel_tn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_tt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_tr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_tc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,

cgemm_small_kernel_rn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_rt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_rr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_rc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,

cgemm_small_kernel_cn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_ct: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_cr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_cc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,

cgemm_small_kernel_b0_nn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_b0_nt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_b0_nr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_b0_nc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,

cgemm_small_kernel_b0_tn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_b0_tt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_b0_tr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_b0_tc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,

cgemm_small_kernel_b0_rn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_b0_rt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_b0_rr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_b0_rc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,

cgemm_small_kernel_b0_cn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_b0_ct: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_b0_cr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
cgemm_small_kernel_b0_cc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
//#endif //def SMALL_MATRIX_OPT

ctrsm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrsm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrsm_kernel_LR: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrsm_kernel_LC: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrsm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrsm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrsm_kernel_RR: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrsm_kernel_RC: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,

ctrsm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_outucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_outncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrsm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,

ctrmm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrmm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrmm_kernel_RR: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrmm_kernel_RC: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrmm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrmm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrmm_kernel_LR: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,
ctrmm_kernel_LC: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong, BlasLong) callconv(.c) i32 = null,

ctrmm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_outucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_outncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
ctrmm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,

csymm_iutcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
csymm_iltcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
csymm_outcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
csymm_oltcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,

chemm_iutcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
chemm_iltcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
chemm_outcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
chemm_oltcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,

cgemm3m_p: i32 = 0,
cgemm3m_q: i32 = 0,
cgemm3m_r: i32 = 0,
cgemm3m_unroll_m: i32 = 0,
cgemm3m_unroll_n: i32 = 0,
cgemm3m_unroll_mn: i32 = 0,

cgemm3m_kernel: ?*const fn (BlasLong, BlasLong, BlasLong, f32, f32, [*]f32, [*]f32, [*]f32, BlasLong) callconv(.c) i32 = null,

cgemm3m_incopyb: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemm3m_incopyr: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemm3m_incopyi: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemm3m_itcopyb: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemm3m_itcopyr: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
cgemm3m_itcopyi: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,

cgemm3m_oncopyb: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
cgemm3m_oncopyr: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
cgemm3m_oncopyi: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
cgemm3m_otcopyb: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
cgemm3m_otcopyr: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
cgemm3m_otcopyi: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,

csymm3m_iucopyb: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
csymm3m_ilcopyb: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
csymm3m_iucopyr: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
csymm3m_ilcopyr: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
csymm3m_iucopyi: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
csymm3m_ilcopyi: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,

csymm3m_oucopyb: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
csymm3m_olcopyb: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
csymm3m_oucopyr: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
csymm3m_olcopyr: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
csymm3m_oucopyi: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
csymm3m_olcopyi: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,

chemm3m_iucopyb: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
chemm3m_ilcopyb: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
chemm3m_iucopyr: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
chemm3m_ilcopyr: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
chemm3m_iucopyi: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,
chemm3m_ilcopyi: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, [*]f32) callconv(.c) i32 = null,

chemm3m_oucopyb: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
chemm3m_olcopyb: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
chemm3m_oucopyr: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
chemm3m_olcopyr: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
chemm3m_oucopyi: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,
chemm3m_olcopyi: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, BlasLong, BlasLong, f32, f32, [*]f32) callconv(.c) i32 = null,

cneg_tcopy: ?*const fn (BlasLong, BlasLong, [*]f32, BlasLong, [*]f32) callconv(.c) i32 = null,
claswp_ncopy: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f32, BlasLong, [*]isize, [*]f32) callconv(.c) i32 = null,
//#endif //(BUILD_COMPLEX==1)

// double-precision complex (z)
//#if (BUILD_COMPLEX16 == 1)
zgemm_p: i32 = 0,
zgemm_q: i32 = 0,
zgemm_r: i32 = 0,
zgemm_unroll_m: i32 = 0,
zgemm_unroll_n: i32 = 0,
zgemm_unroll_mn: i32 = 0,

zamax_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
zamin_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
izamax_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) BlasLong = null,
izamin_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) BlasLong = null,

znrm2_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
zasum_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
zsum_k: ?*const fn (BlasLong, [*]f64, BlasLong) callconv(.c) f64 = null,
zcopy_k: ?*const fn (BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zdotu_k: ?*const fn (BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) CComplex64 = null,
zdotc_k: ?*const fn (BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) CComplex64 = null,
zdrot_k: ?*const fn (BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, f64, f64) callconv(.c) i32 = null,

zaxpy_k: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zaxpyc_k: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zscal_k: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zswap_k: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,

zgemv_n: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemv_t: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemv_r: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemv_c: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemv_o: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemv_u: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemv_s: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemv_d: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgeru_k: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgerc_k: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,

zsymv_L: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zsymv_U: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zhemv_L: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zhemv_U: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zhemv_M: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zhemv_V: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,

//#ifdef ARCH_ARM64
/// ARCH_ARM64 only
sme_zgemm_kernel: ?*const fn ([*:0]const u8, [*:0]const u8, BlasLong, BlasLong, BlasLong, f64, f64, [*]const f64, BlasLong, [*]const f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) void = null,
//#endif //def ARCH_ARM64

zgemm_kernel_n: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_kernel_l: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_kernel_r: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_kernel_b: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_beta: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,

zgemm_incopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemm_itcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemm_oncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemm_otcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,

//#ifdef SMALL_MATRIX_OPT //(double complex)
zgemm_small_matrix_permit: ?*const fn (i32, i32, BlasLong, BlasLong, BlasLong, f64, f64, f64, f64) callconv(.c) i32 = null,

zgemm_small_kernel_nn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_nt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_nr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_nc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,

zgemm_small_kernel_tn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_tt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_tr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_tc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,

zgemm_small_kernel_rn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_rt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_rr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_rc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,

zgemm_small_kernel_cn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_ct: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_cr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_cc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,

zgemm_small_kernel_b0_nn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_b0_nt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_b0_nr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_b0_nc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,

zgemm_small_kernel_b0_tn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_b0_tt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_b0_tr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_b0_tc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,

zgemm_small_kernel_b0_rn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_b0_rt: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_b0_rr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_b0_rc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,

zgemm_small_kernel_b0_cn: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_b0_ct: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_b0_cr: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zgemm_small_kernel_b0_cc: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
//#endif //def SMALL_MATRIX_OPT

ztrsm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrsm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrsm_kernel_LR: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrsm_kernel_LC: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrsm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrsm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrsm_kernel_RR: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrsm_kernel_RC: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,

ztrsm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_outucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_outncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrsm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,

ztrmm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrmm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrmm_kernel_RR: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrmm_kernel_RC: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrmm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrmm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrmm_kernel_LR: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,
ztrmm_kernel_LC: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong, BlasLong) callconv(.c) i32 = null,

ztrmm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_outucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_outncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
ztrmm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,

zsymm_iutcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zsymm_iltcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zsymm_outcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zsymm_oltcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,

zhemm_iutcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zhemm_iltcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zhemm_outcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zhemm_oltcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,

zgemm3m_p: i32 = 0,
zgemm3m_q: i32 = 0,
zgemm3m_r: i32 = 0,
zgemm3m_unroll_m: i32 = 0,
zgemm3m_unroll_n: i32 = 0,
zgemm3m_unroll_mn: i32 = 0,

zgemm3m_kernel: ?*const fn (BlasLong, BlasLong, BlasLong, f64, f64, [*]f64, [*]f64, [*]f64, BlasLong) callconv(.c) i32 = null,

zgemm3m_incopyb: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemm3m_incopyr: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemm3m_incopyi: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemm3m_itcopyb: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemm3m_itcopyr: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zgemm3m_itcopyi: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,

zgemm3m_oncopyb: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zgemm3m_oncopyr: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zgemm3m_oncopyi: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zgemm3m_otcopyb: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zgemm3m_otcopyr: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zgemm3m_otcopyi: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,

zsymm3m_iucopyb: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zsymm3m_ilcopyb: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zsymm3m_iucopyr: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zsymm3m_ilcopyr: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zsymm3m_iucopyi: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zsymm3m_ilcopyi: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,

zsymm3m_oucopyb: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zsymm3m_olcopyb: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zsymm3m_oucopyr: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zsymm3m_olcopyr: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zsymm3m_oucopyi: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zsymm3m_olcopyi: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,

zhemm3m_iucopyb: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zhemm3m_ilcopyb: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zhemm3m_iucopyr: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zhemm3m_ilcopyr: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zhemm3m_iucopyi: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,
zhemm3m_ilcopyi: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, [*]f64) callconv(.c) i32 = null,

zhemm3m_oucopyb: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zhemm3m_olcopyb: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zhemm3m_oucopyr: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zhemm3m_olcopyr: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zhemm3m_oucopyi: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,
zhemm3m_olcopyi: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, BlasLong, BlasLong, f64, f64, [*]f64) callconv(.c) i32 = null,

zneg_tcopy: ?*const fn (BlasLong, BlasLong, [*]f64, BlasLong, [*]f64) callconv(.c) i32 = null,
zlaswp_ncopy: ?*const fn (BlasLong, BlasLong, BlasLong, [*]f64, BlasLong, [*]isize, [*]f64) callconv(.c) i32 = null,
//#endif //(BUILD_COMPLEX16 == 1)

// EXPRECISION (x-family; ComplexXD is Complex64 unless types.exprecision or types.quad_precision is on)
//#ifdef EXPRECISION
xgemm_p: i32 = 0,
xgemm_q: i32 = 0,
xgemm_r: i32 = 0,
xgemm_unroll_m: i32 = 0,
xgemm_unroll_n: i32 = 0,
xgemm_unroll_mn: i32 = 0,

xamax_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
xamin_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
ixamax_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) BlasLong = null,
ixamin_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) BlasLong = null,

xnrm2_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
xasum_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
xsum_k: ?*const fn (BlasLong, [*]XDouble, BlasLong) callconv(.c) XDouble = null,
xcopy_k: ?*const fn (BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) i32 = null,
xdotu_k: ?*const fn (BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) CComplexXD = null,
xdotc_k: ?*const fn (BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) CComplexXD = null,
xqrot_k: ?*const fn (BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, XDouble, XDouble) callconv(.c) i32 = null,

xaxpy_k: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) i32 = null,
xaxpyc_k: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) i32 = null,
xscal_k: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) i32 = null,
xswap_k: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) i32 = null,

xgemv_n: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemv_t: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemv_r: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemv_c: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemv_o: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemv_u: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemv_s: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemv_d: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgeru_k: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgerc_k: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,

xsymv_L: ?*const fn (BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xsymv_U: ?*const fn (BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xhemv_L: ?*const fn (BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xhemv_U: ?*const fn (BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xhemv_M: ?*const fn (BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xhemv_V: ?*const fn (BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,

xgemm_kernel_n: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong) callconv(.c) i32 = null,
xgemm_kernel_l: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong) callconv(.c) i32 = null,
xgemm_kernel_r: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong) callconv(.c) i32 = null,
xgemm_kernel_b: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong) callconv(.c) i32 = null,
xgemm_beta: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, BlasLong, [*]XDouble, BlasLong, [*]XDouble, BlasLong) callconv(.c) i32 = null,

xgemm_incopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemm_itcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemm_oncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemm_otcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,

xtrsm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrsm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrsm_kernel_LR: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrsm_kernel_LC: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrsm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrsm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrsm_kernel_RR: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrsm_kernel_RC: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,

xtrsm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_outucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_outncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrsm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,

xtrmm_kernel_RN: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrmm_kernel_RT: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrmm_kernel_RR: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrmm_kernel_RC: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrmm_kernel_LN: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrmm_kernel_LT: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrmm_kernel_LR: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,
xtrmm_kernel_LC: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong, BlasLong) callconv(.c) i32 = null,

xtrmm_iunucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_iunncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_iutucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_iutncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_ilnucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_ilnncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_iltucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_iltncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_ounucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_ounncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_outucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_outncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_olnucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_olnncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_oltucopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xtrmm_oltncopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,

xsymm_iutcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xsymm_iltcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xsymm_outcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xsymm_oltcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,

xhemm_iutcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xhemm_iltcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xhemm_outcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xhemm_oltcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,

xgemm3m_p: i32 = 0,
xgemm3m_q: i32 = 0,
xgemm3m_r: i32 = 0,
xgemm3m_unroll_m: i32 = 0,
xgemm3m_unroll_n: i32 = 0,
xgemm3m_unroll_mn: i32 = 0,

xgemm3m_kernel: ?*const fn (BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble, [*]XDouble, [*]XDouble, BlasLong) callconv(.c) i32 = null,

xgemm3m_incopyb: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemm3m_incopyr: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemm3m_incopyi: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemm3m_itcopyb: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemm3m_itcopyr: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xgemm3m_itcopyi: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,

xgemm3m_oncopyb: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xgemm3m_oncopyr: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xgemm3m_oncopyi: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xgemm3m_otcopyb: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xgemm3m_otcopyr: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xgemm3m_otcopyi: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,

xsymm3m_iucopyb: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xsymm3m_ilcopyb: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xsymm3m_iucopyr: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xsymm3m_ilcopyr: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xsymm3m_iucopyi: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xsymm3m_ilcopyi: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,

xsymm3m_oucopyb: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xsymm3m_olcopyb: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xsymm3m_oucopyr: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xsymm3m_olcopyr: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xsymm3m_oucopyi: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xsymm3m_olcopyi: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,

xhemm3m_iucopyb: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xhemm3m_ilcopyb: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xhemm3m_iucopyr: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xhemm3m_ilcopyr: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xhemm3m_iucopyi: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xhemm3m_ilcopyi: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, [*]XDouble) callconv(.c) i32 = null,

xhemm3m_oucopyb: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xhemm3m_olcopyb: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xhemm3m_oucopyr: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xhemm3m_olcopyr: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xhemm3m_oucopyi: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,
xhemm3m_olcopyi: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, BlasLong, BlasLong, XDouble, XDouble, [*]XDouble) callconv(.c) i32 = null,

xneg_tcopy: ?*const fn (BlasLong, BlasLong, [*]XDouble, BlasLong, [*]XDouble) callconv(.c) i32 = null,
xlaswp_ncopy: ?*const fn (BlasLong, BlasLong, BlasLong, [*]XDouble, BlasLong, [*]isize, [*]XDouble) callconv(.c) i32 = null,
//#endif //def EXPRECISION

init: ?*const fn () callconv(.c) void = null,

snum_opt: i32 = 0,
dnum_opt: i32 = 0,
qnum_opt: i32 = 0,
//#if (BUILD_SINGLE==1)
saxpby_k: ?*const fn (BlasLong, f32, [*]f32, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_SINGLE==1)
//#if (BUILD_DOUBLE==1)
daxpby_k: ?*const fn (BlasLong, f64, [*]f64, BlasLong, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_DOUBLE==1)
//#if (BUILD_COMPLEX==1)
caxpby_k: ?*const fn (BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_COMPLEX==1)
//#if (BUILD_COMPLEX16==1)
zaxpby_k: ?*const fn (BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_COMPLEX16==1)

//#if (BUILD_SINGLE==1)
somatcopy_k_cn: ?*const fn (BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
somatcopy_k_ct: ?*const fn (BlasLong, BlasLong, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_SINGLE==1)

//#if (BUILD_DOUBLE==1)
domatcopy_k_cn: ?*const fn (BlasLong, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
domatcopy_k_ct: ?*const fn (BlasLong, BlasLong, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_DOUBLE==1)

//#if (BUILD_COMPLEX==1)
comatcopy_k_cn: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
comatcopy_k_ct: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,

comatcopy_k_cnc: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
comatcopy_k_ctc: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, [*]f32, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_COMPLEX==1)

//#if (BUILD_COMPLEX16==1)
zomatcopy_k_cn: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zomatcopy_k_ct: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,

zomatcopy_k_cnc: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
zomatcopy_k_ctc: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, [*]f64, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_COMPLEX16==1)

//#if (BUILD_SINGLE==1)
simatcopy_k_cn: ?*const fn (BlasLong, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
simatcopy_k_ct: ?*const fn (BlasLong, BlasLong, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_SINGLE==1)

//#if (BUILD_DOUBLE==1)
dimatcopy_k_cn: ?*const fn (BlasLong, BlasLong, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
dimatcopy_k_ct: ?*const fn (BlasLong, BlasLong, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_DOUBLE==1)

//#if (BUILD_COMPLEX==1)
cimatcopy_k_cn: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cimatcopy_k_ct: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,

cimatcopy_k_cnc: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
cimatcopy_k_ctc: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_COMPLEX==1)

//#if (BUILD_COMPLEX16==1)
zimatcopy_k_cn: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zimatcopy_k_ct: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,

zimatcopy_k_cnc: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
zimatcopy_k_ctc: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_COMPLEX16==1)

//#if (BUILD_SINGLE==1)
sgeadd_k: ?*const fn (BlasLong, BlasLong, f32, [*]f32, BlasLong, f32, [*]f32, BlasLong, BlasLong, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_SINGLE==1)
//#if (BUILD_DOUBLE==1)
dgeadd_k: ?*const fn (BlasLong, BlasLong, f64, [*]f64, BlasLong, f64, [*]f64, BlasLong, BlasLong, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_DOUBLE==1)
//#if (BUILD_COMPLEX==1)
cgeadd_k: ?*const fn (BlasLong, BlasLong, f32, f32, [*]f32, BlasLong, f32, f32, [*]f32, BlasLong, BlasLong, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_COMPLEX==1)
//#if (BUILD_COMPLEX16==1)
zgeadd_k: ?*const fn (BlasLong, BlasLong, f64, f64, [*]f64, BlasLong, f64, f64, [*]f64, BlasLong, BlasLong, BlasLong) callconv(.c) i32 = null,
//#endif //(BUILD_COMPLEX16==1)

// #define FUNC_OFFSET(func) (size_t)(&((gotoblas_t *)NULL)->func)

// #define DTB_ENTRIES  gotoblas -> dtb_entries
// #define GEMM_OFFSET_A gotoblas -> offsetA
// #define GEMM_OFFSET_B gotoblas -> offsetB
// #define GEMM_ALIGN gotoblas -> align

// #define HAVE_EX_L2 gotoblas -> exclusive_cache

// #if (BUILD_HFLOAT16==1)
//     #define SHGEMM_P    gotoblas -> shgemm_p
//     #define SHGEMM_Q    gotoblas -> shgemm_q
//     #define SHGEMM_R    gotoblas -> shgemm_r
//     #define SHGEMM_UNROLL_M gotoblas -> shgemm_unroll_m
//     #define SHGEMM_UNROLL_N gotoblas -> shgemm_unroll_n
//     #define SHGEMM_UNROLL_MN gotoblas -> shgemm_unroll_mn
// #endif

// #if (BUILD_BFLOAT16==1)
//     #define BGEMM_P  gotoblas -> bgemm_p
//     #define BGEMM_Q  gotoblas -> bgemm_q
//     #define BGEMM_R  gotoblas -> bgemm_r
//     #define BGEMM_UNROLL_M gotoblas -> bgemm_unroll_m
//     #define BGEMM_UNROLL_N gotoblas -> bgemm_unroll_n
//     #define BGEMM_UNROLL_MN gotoblas -> bgemm_unroll_mn

//     #define SBGEMM_P  gotoblas -> sbgemm_p
//     #define SBGEMM_Q  gotoblas -> sbgemm_q
//     #define SBGEMM_R  gotoblas -> sbgemm_r
//     #define SBGEMM_UNROLL_M gotoblas -> sbgemm_unroll_m
//     #define SBGEMM_UNROLL_N gotoblas -> sbgemm_unroll_n
//     #define SBGEMM_UNROLL_MN gotoblas -> sbgemm_unroll_mn
// #endif

// #if (BUILD_SINGLE==1)
//     #define SGEMM_P  gotoblas -> sgemm_p
//     #define SGEMM_Q  gotoblas -> sgemm_q
//     #define SGEMM_R  gotoblas -> sgemm_r
//     #define SGEMM_UNROLL_M gotoblas -> sgemm_unroll_m
//     #define SGEMM_UNROLL_N gotoblas -> sgemm_unroll_n
//     #define SGEMM_UNROLL_MN gotoblas -> sgemm_unroll_mn
// #endif

// #if (BUILD_DOUBLE==1)
//     #define DGEMM_P  gotoblas -> dgemm_p
//     #define DGEMM_Q  gotoblas -> dgemm_q
//     #define DGEMM_R  gotoblas -> dgemm_r
//     #define DGEMM_UNROLL_M gotoblas -> dgemm_unroll_m
//     #define DGEMM_UNROLL_N gotoblas -> dgemm_unroll_n
//     #define DGEMM_UNROLL_MN gotoblas -> dgemm_unroll_mn
//     #if (BUILD_SINGLE != 1)
//     #define SGEMM_P  gotoblas -> sgemm_p
//     #define SGEMM_Q  gotoblas -> sgemm_q
//     #define SGEMM_R  1024
//     #define SGEMM_UNROLL_M gotoblas -> sgemm_unroll_m
//     #define SGEMM_UNROLL_N gotoblas -> sgemm_unroll_n
//     #define SGEMM_UNROLL_MN gotoblas -> sgemm_unroll_mn
//     #endif
// #endif

// #define QGEMM_P  gotoblas -> qgemm_p
// #define QGEMM_Q  gotoblas -> qgemm_q
// #define QGEMM_R  gotoblas -> qgemm_r
// #define QGEMM_UNROLL_M gotoblas -> qgemm_unroll_m
// #define QGEMM_UNROLL_N gotoblas -> qgemm_unroll_n
// #define QGEMM_UNROLL_MN gotoblas -> qgemm_unroll_mn

// #if (BUILD_COMPLEX==1)
//     #define CGEMM_P  gotoblas -> cgemm_p
//     #define CGEMM_Q  gotoblas -> cgemm_q
//     #define CGEMM_R  gotoblas -> cgemm_r
//     #define CGEMM_UNROLL_M gotoblas -> cgemm_unroll_m
//     #define CGEMM_UNROLL_N gotoblas -> cgemm_unroll_n
//     #define CGEMM_UNROLL_MN gotoblas -> cgemm_unroll_mn
//     #if (BUILD_SINGLE != 1)
//     #define SGEMM_P  gotoblas -> sgemm_p
//     #define SGEMM_Q  gotoblas -> sgemm_q
//     #define SGEMM_R  1024
//     #define SGEMM_UNROLL_M gotoblas -> sgemm_unroll_m
//     #define SGEMM_UNROLL_N gotoblas -> sgemm_unroll_n
//     #define SGEMM_UNROLL_MN gotoblas -> sgemm_unroll_mn
//     #endif
// #endif

// #if (BUILD_COMPLEX16==1)
//     #define ZGEMM_P  gotoblas -> zgemm_p
//     #define ZGEMM_Q  gotoblas -> zgemm_q
//     #define ZGEMM_R  gotoblas -> zgemm_r
//     #define ZGEMM_UNROLL_M gotoblas -> zgemm_unroll_m
//     #define ZGEMM_UNROLL_N gotoblas -> zgemm_unroll_n
//     #define ZGEMM_UNROLL_MN gotoblas -> zgemm_unroll_mn
//     #if (BUILD_DOUBLE != 1)
//     #define DGEMM_P  gotoblas -> dgemm_p
//     #define DGEMM_Q  gotoblas -> dgemm_q
//     #define DGEMM_R  1024
//     #define DGEMM_UNROLL_M gotoblas -> dgemm_unroll_m
//     #define DGEMM_UNROLL_N gotoblas -> dgemm_unroll_n
//     #define DGEMM_UNROLL_MN gotoblas -> dgemm_unroll_mn
//     #endif
//     #if (BUILD_COMPLEX != 1)
//     #define CGEMM_P  gotoblas -> cgemm_p
//     #define CGEMM_Q  gotoblas -> cgemm_q
//     #define CGEMM_R  gotoblas -> cgemm_r
//     #define CGEMM_UNROLL_M gotoblas -> cgemm_unroll_m
//     #define CGEMM_UNROLL_N gotoblas -> cgemm_unroll_n
//     #define CGEMM_UNROLL_MN gotoblas -> cgemm_unroll_mn
//     #endif
// #endif

// #define XGEMM_P  gotoblas -> xgemm_p
// #define XGEMM_Q  gotoblas -> xgemm_q
// #define XGEMM_R  gotoblas -> xgemm_r
// #define XGEMM_UNROLL_M gotoblas -> xgemm_unroll_m
// #define XGEMM_UNROLL_N gotoblas -> xgemm_unroll_n
// #define XGEMM_UNROLL_MN gotoblas -> xgemm_unroll_mn

// #define CGEMM3M_P  gotoblas -> cgemm3m_p
// #define CGEMM3M_Q  gotoblas -> cgemm3m_q
// #define CGEMM3M_R  gotoblas -> cgemm3m_r
// #define CGEMM3M_UNROLL_M gotoblas -> cgemm3m_unroll_m
// #define CGEMM3M_UNROLL_N gotoblas -> cgemm3m_unroll_n
// #define CGEMM3M_UNROLL_MN gotoblas -> cgemm3m_unroll_mn

// #define ZGEMM3M_P  gotoblas -> zgemm3m_p
// #define ZGEMM3M_Q  gotoblas -> zgemm3m_q
// #define ZGEMM3M_R  gotoblas -> zgemm3m_r
// #define ZGEMM3M_UNROLL_M gotoblas -> zgemm3m_unroll_m
// #define ZGEMM3M_UNROLL_N gotoblas -> zgemm3m_unroll_n
// #define ZGEMM3M_UNROLL_MN gotoblas -> zgemm3m_unroll_mn

// #define XGEMM3M_P  gotoblas -> xgemm3m_p
// #define XGEMM3M_Q  gotoblas -> xgemm3m_q
// #define XGEMM3M_R  gotoblas -> xgemm3m_r
// #define XGEMM3M_UNROLL_M gotoblas -> xgemm3m_unroll_m
// #define XGEMM3M_UNROLL_N gotoblas -> xgemm3m_unroll_n
// #define XGEMM3M_UNROLL_MN gotoblas -> xgemm3m_unroll_mn
