//! Accelerate's routines by precision: which field of the bound table
//! (`binding.Table`) each routine over each element type calls.

const std = @import("std");
const types = @import("../types.zig");
const binding = @import("binding.zig");

const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

/// A field of the bound table, named as its symbol.
pub const Symbol = std.meta.FieldEnum(binding.Table(i64));

/// The routines rstd sends to Accelerate; `dot` is real, `dotu` complex.
pub const Routine = enum {
    // level 1
    copy,
    swap,
    scal,
    rot,
    rotm,
    axpy,
    axpby,
    dot,
    dotu,
    dotc,
    nrm2,
    asum,

    // level 2
    gemv,
    gbmv,
    ger,
    geru,
    gerc,
    symv,
    hemv,
    sbmv,
    hbmv,
    spmv,
    hpmv,
    trmv,
    trsv,
    tbmv,
    tbsv,
    tpmv,
    tpsv,
    syr,
    her,
    spr,
    hpr,
    syr2,
    her2,
    spr2,
    hpr2,

    // level 3
    gemm,
    symm,
    hemm,
    syrk,
    herk,
    syr2k,
    her2k,
    trmm,
    trsm,

    // LAPACK
    potrf,
    getrf,
    getrs,
};

/// The symbols of `r` over `f32`, `f64`, `Complex32` and `Complex64`; null
/// where it has none.
fn symbols(comptime r: Routine) [4]?Symbol {
    return switch (r) {
        .copy => .{ .cblas_scopy, .cblas_dcopy, .cblas_ccopy, .cblas_zcopy },
        .swap => .{ .cblas_sswap, .cblas_dswap, .cblas_cswap, .cblas_zswap },
        .scal => .{ .cblas_sscal, .cblas_dscal, .cblas_cscal, .cblas_zscal },
        .rot => .{ .cblas_srot, .cblas_drot, .cblas_csrot, .cblas_zdrot },
        .rotm => .{ .cblas_srotm, .cblas_drotm, null, null },
        .axpy => .{ .cblas_saxpy, .cblas_daxpy, .cblas_caxpy, .cblas_zaxpy },
        .axpby => .{ .catlas_saxpby, .catlas_daxpby, .catlas_caxpby, .catlas_zaxpby },
        .dot => .{ .cblas_sdot, .cblas_ddot, null, null },
        .dotu => .{ null, null, .cblas_cdotu_sub, .cblas_zdotu_sub },
        .dotc => .{ null, null, .cblas_cdotc_sub, .cblas_zdotc_sub },
        .nrm2 => .{ .cblas_snrm2, .cblas_dnrm2, .cblas_scnrm2, .cblas_dznrm2 },
        .asum => .{ .cblas_sasum, .cblas_dasum, .cblas_scasum, .cblas_dzasum },
        .gemv => .{ .cblas_sgemv, .cblas_dgemv, .cblas_cgemv, .cblas_zgemv },
        .gbmv => .{ .cblas_sgbmv, .cblas_dgbmv, .cblas_cgbmv, .cblas_zgbmv },
        .ger => .{ .cblas_sger, .cblas_dger, null, null },
        .geru => .{ null, null, .cblas_cgeru, .cblas_zgeru },
        .gerc => .{ null, null, .cblas_cgerc, .cblas_zgerc },
        .symv => .{ .cblas_ssymv, .cblas_dsymv, null, null },
        .hemv => .{ null, null, .cblas_chemv, .cblas_zhemv },
        .sbmv => .{ .cblas_ssbmv, .cblas_dsbmv, null, null },
        .hbmv => .{ null, null, .cblas_chbmv, .cblas_zhbmv },
        .spmv => .{ .cblas_sspmv, .cblas_dspmv, null, null },
        .hpmv => .{ null, null, .cblas_chpmv, .cblas_zhpmv },
        .trmv => .{ .cblas_strmv, .cblas_dtrmv, .cblas_ctrmv, .cblas_ztrmv },
        .trsv => .{ .cblas_strsv, .cblas_dtrsv, .cblas_ctrsv, .cblas_ztrsv },
        .tbmv => .{ .cblas_stbmv, .cblas_dtbmv, .cblas_ctbmv, .cblas_ztbmv },
        .tbsv => .{ .cblas_stbsv, .cblas_dtbsv, .cblas_ctbsv, .cblas_ztbsv },
        .tpmv => .{ .cblas_stpmv, .cblas_dtpmv, .cblas_ctpmv, .cblas_ztpmv },
        .tpsv => .{ .cblas_stpsv, .cblas_dtpsv, .cblas_ctpsv, .cblas_ztpsv },
        .syr => .{ .cblas_ssyr, .cblas_dsyr, null, null },
        .her => .{ null, null, .cblas_cher, .cblas_zher },
        .spr => .{ .cblas_sspr, .cblas_dspr, null, null },
        .hpr => .{ null, null, .cblas_chpr, .cblas_zhpr },
        .syr2 => .{ .cblas_ssyr2, .cblas_dsyr2, null, null },
        .her2 => .{ null, null, .cblas_cher2, .cblas_zher2 },
        .spr2 => .{ .cblas_sspr2, .cblas_dspr2, null, null },
        .hpr2 => .{ null, null, .cblas_chpr2, .cblas_zhpr2 },
        .gemm => .{ .cblas_sgemm, .cblas_dgemm, .cblas_cgemm, .cblas_zgemm },
        .symm => .{ .cblas_ssymm, .cblas_dsymm, .cblas_csymm, .cblas_zsymm },
        .hemm => .{ null, null, .cblas_chemm, .cblas_zhemm },
        .syrk => .{ .cblas_ssyrk, .cblas_dsyrk, .cblas_csyrk, .cblas_zsyrk },
        .herk => .{ null, null, .cblas_cherk, .cblas_zherk },
        .syr2k => .{ .cblas_ssyr2k, .cblas_dsyr2k, .cblas_csyr2k, .cblas_zsyr2k },
        .her2k => .{ null, null, .cblas_cher2k, .cblas_zher2k },
        .trmm => .{ .cblas_strmm, .cblas_dtrmm, .cblas_ctrmm, .cblas_ztrmm },
        .trsm => .{ .cblas_strsm, .cblas_dtrsm, .cblas_ctrsm, .cblas_ztrsm },
        .potrf => .{ .spotrf_, .dpotrf_, .cpotrf_, .zpotrf_ },
        .getrf => .{ .sgetrf_, .dgetrf_, .cgetrf_, .zgetrf_ },
        .getrs => .{ .sgetrs_, .dgetrs_, .cgetrs_, .zgetrs_ },
    };
}

/// The symbol of `r` over `T`; a compile error where Accelerate has none.
pub fn routine(comptime r: Routine, comptime T: type) Symbol {
    const i = switch (T) {
        f32 => 0,
        f64 => 1,
        Complex32 => 2,
        Complex64 => 3,
        else => @compileError("Accelerate has no " ++ @typeName(T) ++ " routines"),
    };
    const sym = comptime symbols(r)[i];
    return sym orelse @compileError("Accelerate has no " ++ @tagName(r) ++ " over " ++ @typeName(T));
}
