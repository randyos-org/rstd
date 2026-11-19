//! Accelerate's CBLAS table and its loading (`apple_accel.zig`).

const std = @import("std");
const sysinfo = @import("builtin");
const dll = @import("../../../dll.zig");
const types = @import("../types.zig");

/// Where `load` can find Accelerate.
pub const available = sysinfo.target.os.tag == .macos;
/// x86_64 macOS, whose Accelerate has quirks of its own (design.md).
pub const x86_64 = sysinfo.target.cpu.arch == .x86_64;
/// Whether a real upper `?potrf_` writes the strictly lower triangle: seen on
/// Apple Silicon (macOS 27.2), not on x86_64 (issue #537).
pub const potrf_scratches_lower = sysinfo.target.cpu.arch == .aarch64;

pub const Order = enum(c_int) { row_major = 101, col_major = 102 };
pub const Transpose = enum(c_int) { no_trans = 111, trans = 112, conj_trans = 113 };
pub const Uplo = enum(c_int) { upper = 121, lower = 122 };
pub const Diag = enum(c_int) { non_unit = 131, unit = 132 };
pub const Side = enum(c_int) { left = 141, right = 142 };

/// `BLAS_THREADING` (macOS 15+).
pub const Threading = enum(c_uint) { multi_threaded = 0, single_threaded = 1, _ };

/// One routine's C signature for element `T` and integer `I`. Complex
/// scalars and arrays are `void *` in CBLAS.
fn Sig(comptime I: type, comptime T: type) type {
    const cx = types.isComplex(T);
    return struct {
        const R = types.Float(T);
        /// A scalar argument: by value if real, by pointer if complex.
        const S = if (cx) *const anyopaque else T;
        const X = if (cx) *const anyopaque else [*]const T;
        const Y = if (cx) *anyopaque else [*]T;
        const O = Order;
        const Tr = Transpose;
        const U = Uplo;
        const Dg = Diag;
        const Sd = Side;

        // level 1
        const Dot = *const fn (I, X, I, X, I) callconv(.c) T;
        const DotSub = *const fn (I, X, I, X, I, *anyopaque) callconv(.c) void;
        const Norm = *const fn (I, X, I) callconv(.c) R;
        const Index = *const fn (I, X, I) callconv(.c) I;
        const Swap = *const fn (I, Y, I, Y, I) callconv(.c) void;
        const Copy = *const fn (I, X, I, Y, I) callconv(.c) void;
        const Axpy = *const fn (I, S, X, I, Y, I) callconv(.c) void;
        const Axpby = *const fn (I, S, X, I, S, Y, I) callconv(.c) void;
        const Scal = *const fn (I, S, Y, I) callconv(.c) void;
        const ScalReal = *const fn (I, R, Y, I) callconv(.c) void;
        const Rot = *const fn (I, Y, I, Y, I, R, R) callconv(.c) void;
        const Rotm = *const fn (I, Y, I, Y, I, [*]const R) callconv(.c) void;
        const Rotg = if (cx)
            *const fn (*anyopaque, *const anyopaque, *R, *anyopaque) callconv(.c) void
        else
            *const fn (*R, *R, *R, *R) callconv(.c) void;
        const Rotmg = *const fn (*R, *R, *R, R, [*]R) callconv(.c) void;

        // level 2
        const Gemv = *const fn (O, Tr, I, I, S, X, I, X, I, S, Y, I) callconv(.c) void;
        const Gbmv = *const fn (O, Tr, I, I, I, I, S, X, I, X, I, S, Y, I) callconv(.c) void;
        const Trmv = *const fn (O, U, Tr, Dg, I, X, I, Y, I) callconv(.c) void;
        const Tbmv = *const fn (O, U, Tr, Dg, I, I, X, I, Y, I) callconv(.c) void;
        const Tpmv = *const fn (O, U, Tr, Dg, I, X, Y, I) callconv(.c) void;
        const Symv = *const fn (O, U, I, S, X, I, X, I, S, Y, I) callconv(.c) void;
        const Sbmv = *const fn (O, U, I, I, S, X, I, X, I, S, Y, I) callconv(.c) void;
        const Spmv = *const fn (O, U, I, S, X, X, I, S, Y, I) callconv(.c) void;
        const Ger = *const fn (O, I, I, S, X, I, X, I, Y, I) callconv(.c) void;
        /// `syr` and `her`: `alpha` is real.
        const Syr = *const fn (O, U, I, R, X, I, Y, I) callconv(.c) void;
        const Spr = *const fn (O, U, I, R, X, I, Y) callconv(.c) void;
        const Syr2 = *const fn (O, U, I, S, X, I, X, I, Y, I) callconv(.c) void;
        const Spr2 = *const fn (O, U, I, S, X, I, X, I, Y) callconv(.c) void;

        // level 3
        const Gemm = *const fn (O, Tr, Tr, I, I, I, S, X, I, X, I, S, Y, I) callconv(.c) void;
        const Symm = *const fn (O, Sd, U, I, I, S, X, I, X, I, S, Y, I) callconv(.c) void;
        const Syrk = *const fn (O, U, Tr, I, I, S, X, I, S, Y, I) callconv(.c) void;
        /// `herk`: `alpha` and `beta` are real.
        const Herk = *const fn (O, U, Tr, I, I, R, X, I, R, Y, I) callconv(.c) void;
        const Syr2k = *const fn (O, U, Tr, I, I, S, X, I, X, I, S, Y, I) callconv(.c) void;
        /// `her2k`: `beta` is real.
        const Her2k = *const fn (O, U, Tr, I, I, S, X, I, X, I, R, Y, I) callconv(.c) void;
        const Trmm = *const fn (O, Sd, U, Tr, Dg, I, I, S, X, I, Y, I) callconv(.c) void;

        // LAPACK: Fortran's, every argument by pointer
        const Potrf = *const fn (*const u8, *const I, Y, *const I, *I) callconv(.c) void;
    };
}

/// The CBLAS routines rstd has, and the LAPACK ones it routes (`?potrf_`),
/// with `I` its integer (`i64` for ILP64,
/// `c_int` otherwise). `catlas_?axpby` is Apple's `axpby`. A field is null
/// where this Accelerate has no such symbol.
pub fn Table(comptime I: type) type {
    const s = Sig(I, f32);
    const d = Sig(I, f64);
    const c = Sig(I, types.Complex32);
    const z = Sig(I, types.Complex64);
    return struct {
        // level 1
        cblas_sdot: ?s.Dot,
        cblas_ddot: ?d.Dot,
        cblas_sdsdot: ?*const fn (I, f32, [*]const f32, I, [*]const f32, I) callconv(.c) f32,
        cblas_dsdot: ?*const fn (I, [*]const f32, I, [*]const f32, I) callconv(.c) f64,
        cblas_cdotu_sub: ?c.DotSub,
        cblas_cdotc_sub: ?c.DotSub,
        cblas_zdotu_sub: ?z.DotSub,
        cblas_zdotc_sub: ?z.DotSub,
        cblas_snrm2: ?s.Norm,
        cblas_dnrm2: ?d.Norm,
        cblas_scnrm2: ?c.Norm,
        cblas_dznrm2: ?z.Norm,
        cblas_sasum: ?s.Norm,
        cblas_dasum: ?d.Norm,
        cblas_scasum: ?c.Norm,
        cblas_dzasum: ?z.Norm,
        cblas_isamax: ?s.Index,
        cblas_idamax: ?d.Index,
        cblas_icamax: ?c.Index,
        cblas_izamax: ?z.Index,
        cblas_isamin: ?s.Index,
        cblas_idamin: ?d.Index,
        cblas_icamin: ?c.Index,
        cblas_izamin: ?z.Index,
        cblas_sswap: ?s.Swap,
        cblas_dswap: ?d.Swap,
        cblas_cswap: ?c.Swap,
        cblas_zswap: ?z.Swap,
        cblas_scopy: ?s.Copy,
        cblas_dcopy: ?d.Copy,
        cblas_ccopy: ?c.Copy,
        cblas_zcopy: ?z.Copy,
        cblas_saxpy: ?s.Axpy,
        cblas_daxpy: ?d.Axpy,
        cblas_caxpy: ?c.Axpy,
        cblas_zaxpy: ?z.Axpy,
        catlas_saxpby: ?s.Axpby,
        catlas_daxpby: ?d.Axpby,
        catlas_caxpby: ?c.Axpby,
        catlas_zaxpby: ?z.Axpby,
        cblas_sscal: ?s.Scal,
        cblas_dscal: ?d.Scal,
        cblas_cscal: ?c.Scal,
        cblas_zscal: ?z.Scal,
        cblas_csscal: ?c.ScalReal,
        cblas_zdscal: ?z.ScalReal,
        cblas_srot: ?s.Rot,
        cblas_drot: ?d.Rot,
        cblas_csrot: ?c.Rot,
        cblas_zdrot: ?z.Rot,
        cblas_srotm: ?s.Rotm,
        cblas_drotm: ?d.Rotm,
        cblas_srotg: ?s.Rotg,
        cblas_drotg: ?d.Rotg,
        cblas_crotg: ?c.Rotg,
        cblas_zrotg: ?z.Rotg,
        cblas_srotmg: ?s.Rotmg,
        cblas_drotmg: ?d.Rotmg,

        // level 2
        cblas_sgemv: ?s.Gemv,
        cblas_dgemv: ?d.Gemv,
        cblas_cgemv: ?c.Gemv,
        cblas_zgemv: ?z.Gemv,
        cblas_sgbmv: ?s.Gbmv,
        cblas_dgbmv: ?d.Gbmv,
        cblas_cgbmv: ?c.Gbmv,
        cblas_zgbmv: ?z.Gbmv,
        cblas_strmv: ?s.Trmv,
        cblas_dtrmv: ?d.Trmv,
        cblas_ctrmv: ?c.Trmv,
        cblas_ztrmv: ?z.Trmv,
        cblas_strsv: ?s.Trmv,
        cblas_dtrsv: ?d.Trmv,
        cblas_ctrsv: ?c.Trmv,
        cblas_ztrsv: ?z.Trmv,
        cblas_stbmv: ?s.Tbmv,
        cblas_dtbmv: ?d.Tbmv,
        cblas_ctbmv: ?c.Tbmv,
        cblas_ztbmv: ?z.Tbmv,
        cblas_stbsv: ?s.Tbmv,
        cblas_dtbsv: ?d.Tbmv,
        cblas_ctbsv: ?c.Tbmv,
        cblas_ztbsv: ?z.Tbmv,
        cblas_stpmv: ?s.Tpmv,
        cblas_dtpmv: ?d.Tpmv,
        cblas_ctpmv: ?c.Tpmv,
        cblas_ztpmv: ?z.Tpmv,
        cblas_stpsv: ?s.Tpmv,
        cblas_dtpsv: ?d.Tpmv,
        cblas_ctpsv: ?c.Tpmv,
        cblas_ztpsv: ?z.Tpmv,
        cblas_ssymv: ?s.Symv,
        cblas_dsymv: ?d.Symv,
        cblas_chemv: ?c.Symv,
        cblas_zhemv: ?z.Symv,
        cblas_ssbmv: ?s.Sbmv,
        cblas_dsbmv: ?d.Sbmv,
        cblas_chbmv: ?c.Sbmv,
        cblas_zhbmv: ?z.Sbmv,
        cblas_sspmv: ?s.Spmv,
        cblas_dspmv: ?d.Spmv,
        cblas_chpmv: ?c.Spmv,
        cblas_zhpmv: ?z.Spmv,
        cblas_sger: ?s.Ger,
        cblas_dger: ?d.Ger,
        cblas_cgeru: ?c.Ger,
        cblas_cgerc: ?c.Ger,
        cblas_zgeru: ?z.Ger,
        cblas_zgerc: ?z.Ger,
        cblas_ssyr: ?s.Syr,
        cblas_dsyr: ?d.Syr,
        cblas_cher: ?c.Syr,
        cblas_zher: ?z.Syr,
        cblas_sspr: ?s.Spr,
        cblas_dspr: ?d.Spr,
        cblas_chpr: ?c.Spr,
        cblas_zhpr: ?z.Spr,
        cblas_ssyr2: ?s.Syr2,
        cblas_dsyr2: ?d.Syr2,
        cblas_cher2: ?c.Syr2,
        cblas_zher2: ?z.Syr2,
        cblas_sspr2: ?s.Spr2,
        cblas_dspr2: ?d.Spr2,
        cblas_chpr2: ?c.Spr2,
        cblas_zhpr2: ?z.Spr2,

        // level 3
        cblas_sgemm: ?s.Gemm,
        cblas_dgemm: ?d.Gemm,
        cblas_cgemm: ?c.Gemm,
        cblas_zgemm: ?z.Gemm,
        cblas_ssymm: ?s.Symm,
        cblas_dsymm: ?d.Symm,
        cblas_csymm: ?c.Symm,
        cblas_zsymm: ?z.Symm,
        cblas_chemm: ?c.Symm,
        cblas_zhemm: ?z.Symm,
        cblas_ssyrk: ?s.Syrk,
        cblas_dsyrk: ?d.Syrk,
        cblas_csyrk: ?c.Syrk,
        cblas_zsyrk: ?z.Syrk,
        cblas_cherk: ?c.Herk,
        cblas_zherk: ?z.Herk,
        cblas_ssyr2k: ?s.Syr2k,
        cblas_dsyr2k: ?d.Syr2k,
        cblas_csyr2k: ?c.Syr2k,
        cblas_zsyr2k: ?z.Syr2k,
        cblas_cher2k: ?c.Her2k,
        cblas_zher2k: ?z.Her2k,
        cblas_strmm: ?s.Trmm,
        cblas_dtrmm: ?d.Trmm,
        cblas_ctrmm: ?c.Trmm,
        cblas_ztrmm: ?z.Trmm,
        cblas_strsm: ?s.Trmm,
        cblas_dtrsm: ?d.Trmm,
        cblas_ctrsm: ?c.Trmm,
        cblas_ztrsm: ?z.Trmm,

        // LAPACK
        spotrf_: ?s.Potrf,
        dpotrf_: ?d.Potrf,
        cpotrf_: ?c.Potrf,
        zpotrf_: ?z.Potrf,
    };
}

/// ILP64 names: the new interface only.
pub fn ilp64Names(comptime field: []const u8) []const [:0]const u8 {
    if (comptime fortran(field)) |base| return &.{ base ++ "$NEWLAPACK$ILP64", field ++ "$NEWLAPACK$ILP64" };
    return &.{field ++ "$NEWLAPACK$ILP64"};
}

/// 32-bit names: the new interface, then the classic one.
pub fn lp64Names(comptime field: []const u8) []const [:0]const u8 {
    if (comptime fortran(field)) |base| return &.{ base ++ "$NEWLAPACK", field ++ "$NEWLAPACK", field ++ "" };
    return &.{ field ++ "$NEWLAPACK", field ++ "" };
}

/// A LAPACK field's name without Fortran's trailing underscore. The macOS
/// 13.3 SDK's `libLAPACK.tbd` spells the new interface so (`dpotrf$NEWLAPACK`,
/// tests below); `load` tries the other spelling too, unseen on a Mac.
pub fn fortran(comptime field: []const u8) ?[]const u8 {
    return if (field[field.len - 1] == '_') field[0 .. field.len - 1] else null;
}

/// `BLASSetThreading` and `BLASGetThreading`, macOS 15+.
pub const ThreadApi = struct {
    BLASSetThreading: ?*const fn (Threading) callconv(.c) c_int,
    BLASGetThreading: ?*const fn () callconv(.c) Threading,
};

pub const Blas = union(enum) {
    ilp64: dll.Bound(Table(i64)),
    lp64: dll.Bound(Table(c_int)),
};

pub const AppleAccelerate = struct {
    blas: Blas,
    threads: ThreadApi,

    pub const ThreadError = error{ Unsupported, Rejected };

    /// The calling thread's model; null before macOS 15.
    pub fn threading(self: *const AppleAccelerate) ?Threading {
        const get = self.threads.BLASGetThreading orelse return null;
        return get();
    }

    /// Sets the calling thread's model (`BLASSetThreading`); other threads,
    /// and threads started later, keep the default (multi-threaded).
    pub fn setThreading(self: *const AppleAccelerate, model: Threading) ThreadError!void {
        const set = self.threads.BLASSetThreading orelse return error.Unsupported;
        if (set(model) != 0) return error.Rejected;
    }

    /// Makes the calling thread's calls single-threaded where macOS can
    /// (15+); a no-op before that, where only `VECLIB_MAXIMUM_THREADS` does.
    pub fn singleThreadHere(self: *const AppleAccelerate) void {
        self.setThreading(.single_threaded) catch {};
    }
};

pub const Options = struct {
    /// Run Accelerate on one thread: its threaded `gemv` rounds differently
    /// from its single-threaded one, so only one thread gives the same bits
    /// everywhere. Sets `VECLIB_MAXIMUM_THREADS=1` (`setMaxThreadsEnv`) and
    /// the loading thread's model; other threads need `singleThreadHere`.
    single_threaded: bool = true,
    /// Bind the 32-bit table even where ILP64 is there, as macOS before
    /// 13.3 would; for tests of the LP64 range checks.
    lp64: bool = false,
};

/// The table `load` binds from `lib` (`dll.bind`'s): ILP64 where `lib` has
/// any of it (macOS 13.3+) unless `lp64`, else the 32-bit one, its new names
/// first; null where it has neither.
pub fn bindBlas(lib: anytype, lp64: bool) ?Blas {
    const ilp = dll.bind(Table(i64), lib, ilp64Names);
    if (ilp.count() > 0 and !lp64) return .{ .ilp64 = ilp };
    const lp = dll.bind(Table(c_int), lib, lp64Names);
    if (lp.count() == 0) return null;
    return .{ .lp64 = lp };
}

/// Opens Accelerate and binds its CBLAS; null off macOS or without it. The
/// framework handle is kept open for the process: the table points into it.
pub fn load(options: Options) ?AppleAccelerate {
    if (comptime !available) return null;
    var fw = dll.openDarwinFramework("Accelerate") catch return null;
    const bound = bindBlas(&fw, options.lp64) orelse return null;
    const acc: AppleAccelerate = .{ .blas = bound, .threads = dll.bind(ThreadApi, &fw, dll.exactName).fns };
    if (options.single_threaded) {
        setMaxThreadsEnv(1) catch {};
        acc.singleThreadHere();
    }
    return acc;
}

pub const EnvError = error{Unsupported};

/// Sets `VECLIB_MAXIMUM_THREADS` to `n` for this process. Accelerate reads it
/// at the process's first BLAS call, so this counts only before that; later
/// it is ignored. `setenv` is not thread safe: call it before other threads
/// read the environment.
pub fn setMaxThreadsEnv(n: u32) EnvError!void {
    if (comptime !available) return error.Unsupported;
    var buf: [16:0]u8 = undefined;
    const v = std.mem.printSentinel(&buf, "{d}", .{n}, 0) catch unreachable;
    if (setenv("VECLIB_MAXIMUM_THREADS", v, 1) != 0) return error.Unsupported;
}

extern "c" fn setenv(name: [*:0]const u8, value: [*:0]const u8, overwrite: c_int) c_int;
