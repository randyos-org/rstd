//! Apple Accelerate's CBLAS, bound at run time: `load` opens the framework
//! through `dll.zig` and fills a `Table` by symbol name, so cross-compiles
//! need no SDK. Declarations are ours, from Apple's public documentation.
//!
//! Symbols, by macOS version:
//! - 13.3+: the new interface (LAPACK 3.9.1; 3.12.0 from macOS 26), as
//!   `cblas_dgemm$NEWLAPACK$ILP64` (64-bit integers, rstd's `isize`) and
//!   `cblas_dgemm$NEWLAPACK` (32-bit). `load` takes ILP64 when it has any.
//! - older: the classic `cblas_dgemm`, 32-bit integers.
//!
//! CBLAS, unlike rstd's interfaces: an `Order` comes first in level 2 and 3;
//! `cblas_i?amax` is 0-based; complex scalars go by pointer; complex `dot`
//! results come back through `cblas_?dot?_sub`'s last argument.
//!
//! Threads: Accelerate runs its own. `VECLIB_MAXIMUM_THREADS` caps them for
//! the process, read at its first BLAS call (`setMaxThreadsEnv`); macOS 15+
//! also has `BLASSetThreading`, per calling thread (`setThreading`). Either
//! one alone makes a call single-threaded; `load` sets both by default.

const std = @import("std");
const sysinfo = @import("builtin");
const dll = @import("../../dll.zig");
const types = @import("types.zig");
const blas = @import("__root__.zig");

const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const isComplex = types.isComplex;
const Float = types.Float;

/// Where `load` can find Accelerate.
pub const available = sysinfo.target.os.tag == .macos;
/// x86_64 macOS, whose Accelerate has quirks of its own (`quirks`).
const x86_64 = sysinfo.target.cpu.arch == .x86_64;

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
    };
}

/// The CBLAS routines rstd has, with `I` its integer (`i64` for ILP64,
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
    };
}

/// ILP64 names: the new interface only.
pub fn ilp64Names(comptime field: []const u8) []const [:0]const u8 {
    return &.{field ++ "$NEWLAPACK$ILP64"};
}

/// 32-bit names: the new interface, then the classic one.
pub fn lp64Names(comptime field: []const u8) []const [:0]const u8 {
    return &.{ field ++ "$NEWLAPACK", field ++ "" };
}

/// `BLASSetThreading` and `BLASGetThreading`, macOS 15+.
const ThreadApi = struct {
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

/// Opens Accelerate and binds its CBLAS; null off macOS or without it. The
/// framework handle is kept open for the process: the table points into it.
pub fn load(options: Options) ?AppleAccelerate {
    if (comptime !available) return null;
    var fw = dll.openDarwinFramework("Accelerate") catch return null;
    const ilp = dll.bind(Table(i64), &fw, ilp64Names);
    const bound: Blas = if (ilp.count() > 0 and !options.lp64) .{ .ilp64 = ilp } else blk: {
        const lp = dll.bind(Table(c_int), &fw, lp64Names);
        if (lp.count() == 0) return null;
        break :blk .{ .lp64 = lp };
    };
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

// ---- routing -----------------------------------------------------------------
//
// The public routines over the four classic types call the functions below
// after their own checks and quick returns, so the errors are rstd's on every
// platform and Accelerate never sees an argument the reference rejects. Each
// is false (or null) where rstd's own kernels run instead: off macOS, while
// off, for a routine, letter or quirk Accelerate does not match, and on the
// LP64 table for a size, stride or operand extent past `c_int`.

/// Whether the public routines over `T` can go to Accelerate.
pub fn covers(comptime T: type) bool {
    return available and (T == f32 or T == f64 or T == Complex32 or T == Complex64);
}

/// The framework, loaded once by `start` and never unloaded.
var loaded: ?AppleAccelerate = null;
var tried = false;
/// Whether the public routines call it.
var on = std.atomic.Value(bool).init(false);
/// Whether this thread's calls are single-threaded (`singleThreadHere`).
threadlocal var thread_set = false;
/// Calls made into Accelerate, counted in tests only.
pub var calls = std.atomic.Value(usize).init(0);

/// `blas.init`'s part: loads Accelerate the first time `enable` holds (it
/// `setenv`s, so only here), then routes to it if it loaded.
pub fn start(enable: bool) void {
    if (comptime !available) return;
    if (enable and !tried) {
        tried = true;
        loaded = load(.{});
    }
    on.store(enable and loaded != null, .release);
}

/// Routes to Accelerate or stops; true if it is in use after. Never loads
/// (only `start` does) nor unloads.
pub fn route(enable: bool) bool {
    const now = enable and loaded != null;
    on.store(now, .release);
    return now;
}

/// Whether the public routines call Accelerate now.
pub fn routing() bool {
    return on.load(.acquire);
}

/// The table for `T`'s routines, or null where rstd runs them. Sets the
/// calling thread single-threaded on its first call.
inline fn use(comptime T: type) ?*const Blas {
    if (comptime !covers(T)) return null;
    _ = blas.vtable();
    if (!on.load(.acquire)) return null;
    const acc = &loaded.?;
    if (!thread_set) {
        thread_set = true;
        acc.singleThreadHere();
    }
    return &acc.blas;
}

/// The CBLAS letter of `T`.
fn letter(comptime T: type) []const u8 {
    return switch (T) {
        f32 => "s",
        f64 => "d",
        Complex32 => "c",
        Complex64 => "z",
        else => unreachable,
    };
}

/// What `call` returns for `name`: its C result, an index as `isize`.
fn Ret(comptime name: []const u8) type {
    const R = @typeInfo(@typeInfo(@typeInfo(@FieldType(Table(i64), name)).optional.child).pointer.child).@"fn".return_type.?;
    return if (R == i64) isize else R;
}

/// Calls `name` with `args`, every `isize` narrowed to the table's integer;
/// null if the table lacks it or, on LP64, any argument or `extents` entry
/// (an operand's span in elements) is past `c_int`.
fn call(comptime name: []const u8, b: *const Blas, extents: []const isize, args: anytype) ?Ret(name) {
    switch (b.*) {
        inline else => |*bound, tag| {
            const f = @field(bound.fns, name) orelse return null;
            const I = if (tag == .ilp64) i64 else c_int;
            if (I != i64) for (extents) |e| if (e > std.math.maxInt(I)) return null;
            var conv: std.meta.ArgsTuple(@typeInfo(@TypeOf(f)).pointer.child) = undefined;
            inline for (0..args.len) |i| {
                conv[i] = if (@TypeOf(args[i]) == isize) std.math.cast(I, args[i]) orelse return null else args[i];
            }
            if (sysinfo.is_test) _ = calls.fetchAdd(1, .monotonic);
            const r = @call(.auto, f, conv);
            return if (@TypeOf(r) == I) @intCast(r) else r;
        },
    }
}

/// The span of `n` elements at stride `inc`.
fn span(n: isize, inc: isize) isize {
    return (n - 1) * @as(isize, @intCast(@abs(inc))) + 1;
}

/// The span of a `rows` by `cols` matrix at leading dimension `ld`.
fn mat(rows: isize, cols: isize, ld: isize) isize {
    return ld * (cols - 1) + rows;
}

/// A scalar as CBLAS takes it: by value if real, by pointer if complex.
inline fn scalar(comptime T: type, v: *const T) if (isComplex(T)) *const anyopaque else T {
    return if (comptime isComplex(T)) v else v.*;
}

/// CBLAS's `TRANS` for `t`; null for the letters it has none for (complex
/// `r`, and `o` to `d`). A real `r` and `c` are `n` and `t`.
fn transpose(comptime T: type, t: types.Trans) ?Transpose {
    if (t.conjX()) return null;
    if (comptime !isComplex(T)) return if (t.transposed()) .trans else .no_trans;
    return switch (t) {
        .n => .no_trans,
        .t => .trans,
        .c => .conj_trans,
        else => null,
    };
}

/// CBLAS's `UPLO`; null for `HEMVREV`'s `v` and `m`.
fn uplo(u: types.Uplo) ?Uplo {
    return switch (u) {
        .u => .upper,
        .l => .lower,
        .v, .m => null,
    };
}

fn diag(d: types.Diag) Diag {
    return if (d == .u) .unit else .non_unit;
}

fn side(s: types.Side) Side {
    return if (s == .l) .left else .right;
}

const col = Order.col_major;

// Quirks: where Accelerate's results differ from the reference's (macOS 27,
// M6; the special-value tests below found them), rstd runs the call.
// - complex `scal`, `axpy`, `axpby`: a scalar with a zero part. Accelerate
//   skips that part's products, so `(1.5, 0) * (inf, 1)` has a finite
//   imaginary part, where the reference's is NaN.
// - `axpby`: a zero `alpha`. Accelerate reads `x`, so NaN or an infinity
//   there reaches `y`.
// - `i?amax`: always (not routed). Accelerate does not take the first NaN.
// - `trmv`, `trsv`, `tbmv`, `tbsv`, `tpmv`, `tpsv`: NaN or an infinity in
//   `x`. Past small orders Accelerate multiplies it by the triangle's zeros
//   (`trmv`; complex `trsv` in `f32`), giving NaN.
// - `trmm`: NaN or an infinity in `B`, the same way.
// - `ger`, `geru`, `gerc`: NaN or an infinity in `x`. Accelerate does not
//   pass over a column whose `y` is zero, as the reference does.
// - `herk`, `her2k`, `her`, `hpr`, `her2`, `hpr2`: the reference sets the
//   diagonal's imaginary parts to zero; Accelerate scales them by `beta`
//   (level 3), and gives NaN where an infinity meets them (`her2`, `herk`,
//   `her2k`). They are zeroed after the call instead, but where `beta` is 1
//   and nothing is added, which leaves `C` as it was and is rstd's.
// - `herk`, `her2k`: NaN or an infinity in `C`'s triangle, with a nonzero
//   `beta`. Past small orders Accelerate scales `C` by `beta` as a complex
//   product, so a non-finite part reaches the other one.
// x86_64 (under Rosetta 2 on the M6; no Intel Mac was at hand):
// - `dsdot`, and `sdot`, `snrm2`, `scnrm2` built on it: always. Its
//   `cblas_dsdot` is off by `f32`'s epsilon, not `f64`'s.
// - `gemm`: NaN or an infinity in `A`. At an edge of `m` it gives NaN where
//   the reference's is infinite.

/// Zeroes the imaginary parts of the diagonal of the order `n` Hermitian
/// matrix at `a`, full (`ld`) or packed (`ld` 0) on the `lower` triangle.
fn realDiagonal(comptime T: type, n: isize, a: [*]T, ld: isize, lower: bool) void {
    if (comptime !isComplex(T)) return;
    const un: usize = @intCast(n);
    for (0..un) |j| {
        const at = if (ld != 0) j * @as(usize, @intCast(ld + 1)) else if (lower) j * (2 * un - j + 1) / 2 else j * (j + 1) / 2 + j;
        a[at].im = 0;
    }
}

/// `nonFinite` of the `lower` or upper triangle of the order `n` matrix at
/// `c`, leading dimension `ld`.
fn nonFiniteTri(comptime T: type, n: isize, c: [*]const T, ld: isize, lower: bool) bool {
    const un: usize = @intCast(n);
    for (0..un) |j| {
        const col_ = c + j * @as(usize, @intCast(ld));
        const bad = if (lower) nonFinite(T, n - @as(isize, @intCast(j)), col_ + j, 1) else nonFinite(T, @intCast(j + 1), col_, 1);
        if (bad) return true;
    }
    return false;
}

/// Whether complex `v` has a zero part.
fn zeroPart(comptime T: type, v: T) bool {
    if (comptime !isComplex(T)) return false;
    return v.re == 0 or v.im == 0;
}

/// Whether `n` elements of `x` at stride `inc`, from its lowest, hold NaN or
/// an infinity.
fn nonFinite(comptime T: type, n: isize, x: [*]const T, inc: isize) bool {
    const F = Float(T);
    const w: usize = if (comptime isComplex(T)) 2 else 1;
    const f: [*]const F = @ptrCast(x);
    const step = w * @abs(inc);
    // `v * 0` is NaN only for NaN and the infinities
    var acc: F = 0;
    for (0..@intCast(n)) |i| {
        for (0..w) |p| acc += f[i * step + p] * 0;
    }
    return acc != acc;
}

/// `nonFinite` of a `rows` by `cols` matrix at leading dimension `ld`.
fn nonFiniteMat(comptime T: type, rows: isize, cols: isize, b: [*]const T, ld: isize) bool {
    for (0..@intCast(cols)) |j| {
        if (nonFinite(T, rows, b + j * @as(usize, @intCast(ld)), 1)) return true;
    }
    return false;
}

// level 1

pub inline fn copy(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        return call("cblas_" ++ letter(T) ++ "copy", b, &.{ span(n, incx), span(n, incy) }, .{ n, x, incx, y, incy }) != null;
    }
}

pub inline fn swap(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        return call("cblas_" ++ letter(T) ++ "swap", b, &.{ span(n, incx), span(n, incy) }, .{ n, x, incx, y, incy }) != null;
    }
}

pub inline fn scal(comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        if (zeroPart(T, alpha)) return false;
        return call("cblas_" ++ letter(T) ++ "scal", b, &.{span(n, incx)}, .{ n, scalar(T, &alpha), x, incx }) != null;
    }
}

/// `x . y`, or `conj(x) . y`. `sdot` is `dsdot` rounded: rstd's contract sums
/// in `f64` (readme, One contract on every core).
pub inline fn dot(comptime T: type, comptime conj: bool, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) ?T {
    if (comptime !covers(T)) return null else {
        const b = use(T) orelse return null;
        const ext = &.{ span(n, incx), span(n, incy) };
        if (comptime T == f32) {
            if (comptime x86_64) return null;
            return @floatCast(call("cblas_dsdot", b, ext, .{ n, x, incx, y, incy }) orelse return null);
        }
        if (comptime T == f64) return call("cblas_ddot", b, ext, .{ n, x, incx, y, incy });
        var r: T = undefined;
        _ = call("cblas_" ++ letter(T) ++ "dot" ++ (if (conj) "c" else "u") ++ "_sub", b, ext, .{ n, x, incx, y, incy, &r }) orelse return null;
        return r;
    }
}

pub inline fn dsdot(n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize) ?f64 {
    const b = use(f32) orelse return null;
    if (comptime x86_64) return null;
    return call("cblas_dsdot", b, &.{ span(n, incx), span(n, incy) }, .{ n, x, incx, y, incy });
}

/// `s` and `c` as rstd's contract has them, the squares summed in `f64` by
/// `dsdot` and the root rounded once; `d` and `z` for `incx > 0`, the
/// reference's only stride.
pub inline fn nrm2(comptime T: type, n: isize, x: [*]T, incx: isize) ?Float(T) {
    if (comptime !covers(T)) return null else {
        const b = use(T) orelse return null;
        const f: [*]Float(T) = @ptrCast(x);
        if (comptime Float(T) == f32) {
            if (comptime x86_64) return null;
            const sq = if (comptime !isComplex(T))
                call("cblas_dsdot", b, &.{span(n, incx)}, .{ n, f, incx, f, incx }) orelse return null
            else if (incx == 1)
                call("cblas_dsdot", b, &.{2 * n}, .{ 2 * n, f, @as(isize, 1), f, @as(isize, 1) }) orelse return null
            else blk: {
                const re = call("cblas_dsdot", b, &.{span(n, 2 * incx)}, .{ n, f, 2 * incx, f, 2 * incx }) orelse return null;
                const im = call("cblas_dsdot", b, &.{span(n, 2 * incx)}, .{ n, f + 1, 2 * incx, f + 1, 2 * incx }) orelse return null;
                break :blk re + im;
            };
            return @floatCast(@sqrt(sq));
        }
        if (incx <= 0) return null;
        return call(if (comptime isComplex(T)) "cblas_dznrm2" else "cblas_dnrm2", b, &.{span(n, incx)}, .{ n, x, incx });
    }
}

/// For `incx > 0`; rstd's kernels give 0 otherwise.
pub inline fn asum(comptime T: type, n: isize, x: [*]T, incx: isize) ?Float(T) {
    if (comptime !covers(T)) return null else {
        const b = use(T) orelse return null;
        if (incx <= 0) return null;
        const name = comptime switch (T) {
            Complex32 => "cblas_scasum",
            Complex64 => "cblas_dzasum",
            else => "cblas_" ++ letter(T) ++ "asum",
        };
        return call(name, b, &.{span(n, incx)}, .{ n, x, incx });
    }
}

pub inline fn rot(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, c: Float(T), s: Float(T)) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        const name = comptime switch (T) {
            Complex32 => "cblas_csrot",
            Complex64 => "cblas_zdrot",
            else => "cblas_" ++ letter(T) ++ "rot",
        };
        return call(name, b, &.{ span(n, incx), span(n, incy) }, .{ n, x, incx, y, incy, c, s }) != null;
    }
}

pub inline fn rotm(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, param: *const [5]T) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        return call("cblas_" ++ letter(T) ++ "rotm", b, &.{ span(n, incx), span(n, incy) }, .{ n, x, incx, y, incy, param }) != null;
    }
}

pub inline fn axpy(comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        if (zeroPart(T, alpha)) return false;
        return call("cblas_" ++ letter(T) ++ "axpy", b, &.{ span(n, incx), span(n, incy) }, .{ n, scalar(T, &alpha), x, incx, y, incy }) != null;
    }
}

pub inline fn axpby(comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        if (zeroPart(T, alpha) or zeroPart(T, beta)) return false;
        if (comptime !isComplex(T)) if (alpha == 0) return false;
        return call("catlas_" ++ letter(T) ++ "axpby", b, &.{ span(n, incx), span(n, incy) }, .{ n, scalar(T, &alpha), x, incx, scalar(T, &beta), y, incy }) != null;
    }
}

// level 2

pub inline fn gemv(comptime T: type, trans: types.Trans, m: isize, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        const t = transpose(T, trans) orelse return false;
        const lenx = if (trans.transposed()) m else n;
        const leny = if (trans.transposed()) n else m;
        const ext = &.{ mat(m, n, lda), span(lenx, incx), span(leny, incy) };
        return call("cblas_" ++ letter(T) ++ "gemv", b, ext, .{ col, t, m, n, scalar(T, &alpha), a, lda, x, incx, scalar(T, &beta), y, incy }) != null;
    }
}

pub inline fn gbmv(comptime T: type, trans: types.Trans, m: isize, n: isize, kl: isize, ku: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        const t = transpose(T, trans) orelse return false;
        const lenx = if (trans.transposed()) m else n;
        const leny = if (trans.transposed()) n else m;
        const ext = &.{ mat(lda, n, lda), span(lenx, incx), span(leny, incy) };
        return call("cblas_" ++ letter(T) ++ "gbmv", b, ext, .{ col, t, m, n, kl, ku, scalar(T, &alpha), a, lda, x, incx, scalar(T, &beta), y, incy }) != null;
    }
}

/// `ger`, `geru`, and `gerc` for `conj`.
pub inline fn ger(comptime T: type, comptime conj: bool, m: isize, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        if (nonFinite(T, m, x, incx)) return false;
        const name = comptime "cblas_" ++ letter(T) ++ "ger" ++ if (!isComplex(T)) "" else if (conj) "c" else "u";
        return call(name, b, &.{ span(m, incx), span(n, incy), mat(m, n, lda) }, .{ col, m, n, scalar(T, &alpha), x, incx, y, incy, a, lda }) != null;
    }
}

/// `symv`, or `hemv` for `herm`.
pub inline fn symv(comptime T: type, comptime herm: bool, ul: types.Uplo, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T) or herm != isComplex(T)) return false else {
        const name = comptime "cblas_" ++ letter(T) ++ if (herm) "hemv" else "symv";
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        return call(name, b, &.{ mat(n, n, lda), span(n, incx), span(n, incy) }, .{ col, u, n, scalar(T, &alpha), a, lda, x, incx, scalar(T, &beta), y, incy }) != null;
    }
}

/// `sbmv` or `hbmv`, packed (`spmv`, `hpmv`) for `packed_`.
pub inline fn sbmv(comptime T: type, comptime herm: bool, comptime packed_: bool, ul: types.Uplo, n: isize, k: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T) or herm != isComplex(T)) return false else {
        const name = comptime "cblas_" ++ letter(T) ++ if (herm) (if (packed_) "hpmv" else "hbmv") else (if (packed_) "spmv" else "sbmv");
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const al = scalar(T, &alpha);
        const be = scalar(T, &beta);
        if (packed_) return call(name, b, &.{ @divExact(n * (n + 1), 2), span(n, incx), span(n, incy) }, .{ col, u, n, al, a, x, incx, be, y, incy }) != null;
        return call(name, b, &.{ mat(lda, n, lda), span(n, incx), span(n, incy) }, .{ col, u, n, k, al, a, lda, x, incx, be, y, incy }) != null;
    }
}

pub const TriKind = enum { trmv, trsv, tbmv, tbsv, tpmv, tpsv };

pub inline fn tri(comptime T: type, comptime kind: TriKind, ul: types.Uplo, trans: types.Trans, dg: types.Diag, n: isize, k: isize, a: [*]T, lda: isize, x: [*]T, incx: isize) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const t = transpose(T, trans) orelse return false;
        if (nonFinite(T, n, x, incx)) return false;
        const name = comptime "cblas_" ++ letter(T) ++ @tagName(kind);
        return switch (kind) {
            .trmv, .trsv => call(name, b, &.{ mat(n, n, lda), span(n, incx) }, .{ col, u, t, diag(dg), n, a, lda, x, incx }),
            .tbmv, .tbsv => call(name, b, &.{ mat(lda, n, lda), span(n, incx) }, .{ col, u, t, diag(dg), n, k, a, lda, x, incx }),
            .tpmv, .tpsv => call(name, b, &.{ @divExact(n * (n + 1), 2), span(n, incx) }, .{ col, u, t, diag(dg), n, a, x, incx }),
        } != null;
    }
}

/// `syr` or `her` (real `alpha`), packed (`spr`, `hpr`) for `packed_`.
pub inline fn syr(comptime T: type, comptime herm: bool, comptime packed_: bool, ul: types.Uplo, n: isize, alpha: if (herm) Float(T) else T, x: [*]T, incx: isize, a: [*]T, lda: isize) bool {
    if (comptime !covers(T) or herm != isComplex(T)) return false else {
        const name = comptime "cblas_" ++ letter(T) ++ if (herm) (if (packed_) "hpr" else "her") else (if (packed_) "spr" else "syr");
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const ok = if (packed_)
            call(name, b, &.{ span(n, incx), @divExact(n * (n + 1), 2) }, .{ col, u, n, alpha, x, incx, a }) != null
        else
            call(name, b, &.{ span(n, incx), mat(n, n, lda) }, .{ col, u, n, alpha, x, incx, a, lda }) != null;
        if (ok and herm) realDiagonal(T, n, a, if (packed_) 0 else lda, ul == .l);
        return ok;
    }
}

/// `syr2` or `her2`, packed (`spr2`, `hpr2`) for `packed_`.
pub inline fn syr2(comptime T: type, comptime herm: bool, comptime packed_: bool, ul: types.Uplo, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) bool {
    if (comptime !covers(T) or herm != isComplex(T)) return false else {
        const name = comptime "cblas_" ++ letter(T) ++ if (herm) (if (packed_) "hpr2" else "her2") else (if (packed_) "spr2" else "syr2");
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const al = scalar(T, &alpha);
        const ok = if (packed_)
            call(name, b, &.{ span(n, incx), span(n, incy), @divExact(n * (n + 1), 2) }, .{ col, u, n, al, x, incx, y, incy, a }) != null
        else
            call(name, b, &.{ span(n, incx), span(n, incy), mat(n, n, lda) }, .{ col, u, n, al, x, incx, y, incy, a, lda }) != null;
        if (ok and herm) realDiagonal(T, n, a, if (packed_) 0 else lda, ul == .l);
        return ok;
    }
}

// level 3

pub inline fn gemm(comptime T: type, transa: types.Trans, transb: types.Trans, m: isize, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b_: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        const ta = transpose(T, transa) orelse return false;
        const tb = transpose(T, transb) orelse return false;
        if (x86_64 and (if (transa.transposed()) nonFiniteMat(T, k, m, a, lda) else nonFiniteMat(T, m, k, a, lda))) return false;
        const ea = if (transa.transposed()) mat(k, m, lda) else mat(m, k, lda);
        const eb = if (transb.transposed()) mat(n, k, ldb) else mat(k, n, ldb);
        return call("cblas_" ++ letter(T) ++ "gemm", b, &.{ ea, eb, mat(m, n, ldc) }, .{ col, ta, tb, m, n, k, scalar(T, &alpha), a, lda, b_, ldb, scalar(T, &beta), c, ldc }) != null;
    }
}

/// `symm`, or `hemm` for `herm`.
pub inline fn symm(comptime T: type, comptime herm: bool, sd: types.Side, ul: types.Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b_: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) bool {
    if (comptime !covers(T) or (herm and !isComplex(T))) return false else {
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const ka = if (sd == .l) m else n;
        const name = comptime "cblas_" ++ letter(T) ++ if (herm) "hemm" else "symm";
        return call(name, b, &.{ mat(ka, ka, lda), mat(m, n, ldb), mat(m, n, ldc) }, .{ col, side(sd), u, m, n, scalar(T, &alpha), a, lda, b_, ldb, scalar(T, &beta), c, ldc }) != null;
    }
}

/// `syrk`, `herk` (real `alpha` and `beta`), `syr2k` or `her2k` (real
/// `beta`), as `herm` and `b_` say; `t` is a transposed `op`.
pub inline fn syrk(comptime T: type, comptime herm: bool, ul: types.Uplo, t: bool, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b_: ?[*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) bool {
    if (comptime !covers(T) or (herm and !isComplex(T))) return false else {
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const tr: Transpose = if (!t) .no_trans else if (herm) .conj_trans else .trans;
        if (comptime herm) {
            if (beta.re == 1 and (k == 0 or (alpha.re == 0 and alpha.im == 0))) return false;
            if (beta.re != 0 and nonFiniteTri(T, n, c, ldc, ul == .l)) return false;
        }
        const ea = if (t) mat(k, n, lda) else mat(n, k, lda);
        const ec = mat(n, n, ldc);
        const l = comptime letter(T);
        const ok = if (b_) |bb| blk: {
            const eb = if (t) mat(k, n, ldb) else mat(n, k, ldb);
            const ext = &.{ ea, eb, ec };
            if (comptime herm) break :blk call("cblas_" ++ l ++ "her2k", b, ext, .{ col, u, tr, n, k, scalar(T, &alpha), a, lda, bb, ldb, beta.re, c, ldc }) != null;
            break :blk call("cblas_" ++ l ++ "syr2k", b, ext, .{ col, u, tr, n, k, scalar(T, &alpha), a, lda, bb, ldb, scalar(T, &beta), c, ldc }) != null;
        } else if (comptime herm)
            call("cblas_" ++ l ++ "herk", b, &.{ ea, ec }, .{ col, u, tr, n, k, alpha.re, a, lda, beta.re, c, ldc }) != null
        else
            call("cblas_" ++ l ++ "syrk", b, &.{ ea, ec }, .{ col, u, tr, n, k, scalar(T, &alpha), a, lda, scalar(T, &beta), c, ldc }) != null;
        if (ok and herm) realDiagonal(T, n, c, ldc, ul == .l);
        return ok;
    }
}

/// `trmm`, or `trsm` for `solve`.
pub inline fn trmm(comptime T: type, comptime solve: bool, sd: types.Side, ul: types.Uplo, transa: types.Trans, dg: types.Diag, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b_: [*]T, ldb: isize) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const t = transpose(T, transa) orelse return false;
        const ka = if (sd == .l) m else n;
        if (!solve and nonFiniteMat(T, m, n, b_, ldb)) return false;
        const name = comptime "cblas_" ++ letter(T) ++ if (solve) "trsm" else "trmm";
        return call(name, b, &.{ mat(ka, ka, lda), mat(m, n, ldb) }, .{ col, side(sd), u, t, diag(dg), m, n, scalar(T, &alpha), a, lda, b_, ldb }) != null;
    }
}

// ---- tests -------------------------------------------------------------------

const zigstd = @import("../../zigstd.zig");
const check = @import("_test/blas1/check.zig");

/// The integer type of `tag`'s table.
fn Int(comptime tag: std.meta.Tag(Blas)) type {
    return if (tag == .ilp64) i64 else c_int;
}

fn fill(comptime F: type, rng: std.Random, out: []F) void {
    for (out) |*v| v.* = rng.float(F) * 2 - 1;
}

test "load binds the CBLAS table on macOS, and nothing elsewhere" {
    const acc = load(.{}) orelse return zigstd.testExpectTrue(!available);
    switch (acc.blas) {
        inline else => |b| {
            try zigstd.testExpectTrue(b.fns.cblas_dgemm != null and b.fns.cblas_ddot != null);
            try zigstd.testExpectTrue(b.fns.cblas_zdotu_sub != null and b.fns.cblas_isamax != null);
        },
    }
    if (acc.threads.BLASSetThreading != null) try zigstd.testExpectEqual(Threading.single_threaded, acc.threading().?);
}

test "the 32-bit table binds the new names first, and agrees with ILP64" {
    if (comptime !available) return error.SkipZigTest;
    const acc = load(.{}) orelse return error.SkipZigTest;
    if (acc.blas != .ilp64) return error.SkipZigTest;
    var lib = try dll.openDarwinFramework("Accelerate");
    const lp = dll.bind(Table(c_int), &lib, lp64Names);
    try zigstd.testExpectEqualStrings("cblas_ddot$NEWLAPACK", lp.name(lp64Names, .cblas_ddot).?);
    var x: [100]f64 = undefined;
    for (&x, 0..) |*v, i| v.* = @floatFromInt(i % 7);
    try zigstd.testExpectEqual(acc.blas.ilp64.fns.cblas_ddot.?(100, &x, 1, &x, 1), lp.fns.cblas_ddot.?(100, &x, 1, &x, 1));
}

test "ddot, dgemm, zdotu_sub and isamax agree with rstd" {
    const acc = load(.{}) orelse return error.SkipZigTest;
    // rstd's own kernels
    _ = blas.vtable();
    const was = routing();
    _ = route(false);
    defer _ = route(was);
    const gpa = zigstd.testAllocator;
    var prng: std.Random.DefaultPrng = .init(0x5eed);
    const rng = prng.random();
    switch (acc.blas) {
        inline else => |b, tag| {
            const I = Int(tag);
            const t = b.fns;

            // ddot, strided both ways
            var x: [600]f64 = undefined;
            var y: [600]f64 = undefined;
            fill(f64, rng, &x);
            fill(f64, rng, &y);
            const got = t.cblas_ddot.?(300, &x, 2, &y, -1);
            const want = try blas.dot(check.io, f64, 300, &x, 2, &y, -1);
            try zigstd.testExpectApproxEqRel(want, got, 1e-12);

            // zdotu_sub and zdotc_sub
            var zx: [257]Complex64 = undefined;
            var zy: [257]Complex64 = undefined;
            fill(f64, rng, @as(*[514]f64, @ptrCast(&zx)));
            fill(f64, rng, @as(*[514]f64, @ptrCast(&zy)));
            var zd: Complex64 = undefined;
            t.cblas_zdotu_sub.?(zx.len, &zx, 1, &zy, 1, &zd);
            const zw = try blas.dot(check.io, Complex64, zx.len, &zx, 1, &zy, 1);
            try zigstd.testExpectApproxEqAbs(zw.re, zd.re, 1e-12);
            try zigstd.testExpectApproxEqAbs(zw.im, zd.im, 1e-12);
            t.cblas_zdotc_sub.?(zx.len, &zx, 1, &zy, 1, &zd);
            const zc = try blas.dotc(check.io, Complex64, zx.len, &zx, 1, &zy, 1);
            try zigstd.testExpectApproxEqAbs(zc.re, zd.re, 1e-12);
            try zigstd.testExpectApproxEqAbs(zc.im, zd.im, 1e-12);

            // isamax: 0-based where rstd's is 1-based
            var sx: [1000]f32 = undefined;
            fill(f32, rng, &sx);
            sx[613] = -3;
            const ia: I = t.cblas_isamax.?(sx.len, &sx, 1);
            try zigstd.testExpectEqual(@as(I, 613), ia);
            try zigstd.testExpectEqual(try blas.index(check.io, f32, .amax, sx.len, &sx, 1), ia + 1);

            // dgemm, column- and row-major
            const m = 37;
            const n = 29;
            const k = 41;
            const a = try gpa.alloc(f64, m * k);
            defer gpa.free(a);
            const bm = try gpa.alloc(f64, k * n);
            defer gpa.free(bm);
            const c0 = try gpa.alloc(f64, m * n);
            defer gpa.free(c0);
            const c1 = try gpa.alloc(f64, m * n);
            defer gpa.free(c1);
            const c2 = try gpa.alloc(f64, m * n);
            defer gpa.free(c2);
            fill(f64, rng, a);
            fill(f64, rng, bm);
            fill(f64, rng, c0);
            @memcpy(c1, c0);
            // A^T, stored k by m
            t.cblas_dgemm.?(.col_major, .trans, .no_trans, m, n, k, 1.5, a.ptr, k, bm.ptr, k, -0.5, c0.ptr, m);
            try blas.gemm(gpa, check.io, f64, .t, .n, m, n, k, 1.5, a.ptr, k, bm.ptr, k, -0.5, c1.ptr, m);
            for (c0, c1) |g, w| try zigstd.testExpectApproxEqAbs(w, g, 1e-12 * k);
            // row-major: C^T = B^T A, the same numbers transposed
            t.cblas_dgemm.?(.row_major, .no_trans, .trans, n, m, k, 1, bm.ptr, k, a.ptr, k, 0, c2.ptr, m);
            try blas.gemm(gpa, check.io, f64, .t, .n, m, n, k, 1, a.ptr, k, bm.ptr, k, 0, c1.ptr, m);
            for (c1, c2) |w, g| try zigstd.testExpectApproxEqAbs(w, g, 1e-12 * k);
        },
    }
}

/// Each routine's result bits.
const Digest = struct { ddot: u64, dnrm2: u64, dasum: u64, dgemv_n: u64, dgemv_t: u64, dgemm: u64 };

const ReproData = struct {
    x: []f64,
    y: []f64,
    a: []f64,
    b: []f64,
    c: []f64,
    v: []f64,
    const n1 = 1_000_000;
    const n2 = 2000;
};

fn hash(s: []const f64) u64 {
    return std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(s));
}

fn digest(comptime I: type, t: *const Table(I), d: *const ReproData) Digest {
    const n1 = ReproData.n1;
    const n2 = ReproData.n2;
    var out: Digest = undefined;
    out.ddot = @bitCast(t.cblas_ddot.?(n1, d.x.ptr, 1, d.y.ptr, 1));
    out.dnrm2 = @bitCast(t.cblas_dnrm2.?(n1, d.x.ptr, 1));
    out.dasum = @bitCast(t.cblas_dasum.?(n1, d.x.ptr, 1));
    t.cblas_dgemv.?(.col_major, .no_trans, n2, n2, 1, d.a.ptr, n2, d.x.ptr, 1, 0, d.v.ptr, 1);
    out.dgemv_n = hash(d.v);
    t.cblas_dgemv.?(.col_major, .trans, n2, n2, 1, d.a.ptr, n2, d.x.ptr, 1, 0, d.v.ptr, 1);
    out.dgemv_t = hash(d.v);
    t.cblas_dgemm.?(.col_major, .no_trans, .no_trans, n2, n2, n2, 1, d.a.ptr, n2, d.b.ptr, n2, 0, d.c.ptr, n2);
    out.dgemm = hash(d.c);
    return out;
}

fn digestOnThread(comptime I: type, acc: *const AppleAccelerate, t: *const Table(I), d: *const ReproData, out: *Digest) void {
    acc.singleThreadHere();
    out.* = digest(I, t, d);
}

// Threaded, Accelerate's `gemv` with `A` untransposed rounds differently from
// one thread (macOS 27, M6); `ddot`, `dnrm2`, `dasum`, `A^T` `gemv` and `gemm`
// did not differ for any thread count tried. `load` pins one thread.
test "single-threaded results are bit-identical across runs and threads" {
    if (comptime !available or sysinfo.single_threaded) return error.SkipZigTest;
    const acc = load(.{}) orelse return error.SkipZigTest;
    const gpa = zigstd.testAllocator;
    var prng: std.Random.DefaultPrng = .init(0xacce1);
    const rng = prng.random();
    const n1 = ReproData.n1;
    const n2 = ReproData.n2;
    var d: ReproData = undefined;
    d.x = try gpa.alloc(f64, n1);
    defer gpa.free(d.x);
    d.y = try gpa.alloc(f64, n1);
    defer gpa.free(d.y);
    d.a = try gpa.alloc(f64, n2 * n2);
    defer gpa.free(d.a);
    d.b = try gpa.alloc(f64, n2 * n2);
    defer gpa.free(d.b);
    d.c = try gpa.alloc(f64, n2 * n2);
    defer gpa.free(d.c);
    d.v = try gpa.alloc(f64, n2);
    defer gpa.free(d.v);
    for ([_][]f64{ d.x, d.y, d.a, d.b }) |s| fill(f64, rng, s);

    switch (acc.blas) {
        inline else => |b, tag| {
            const I = Int(tag);
            const t = &b.fns;
            const first = digest(I, t, &d);
            try zigstd.testExpectEqual(first, digest(I, t, &d));
            var other: Digest = undefined;
            const th = try std.Thread.spawn(.{}, digestOnThread, .{ I, &acc, t, &d, &other });
            th.join();
            try zigstd.testExpectEqual(first, other);
            // the process-wide cap holds even where a thread asks for more
            if (acc.setThreading(.multi_threaded)) {
                defer acc.singleThreadHere();
                try zigstd.testExpectEqual(first, digest(I, t, &d));
            } else |_| {}
        },
    }
}

// ---- routing tests -------------------------------------------------------------

/// Fails unless Accelerate was called since `before` exactly when `want`.
fn expectCalled(want: bool, before: usize, what: []const u8) !void {
    if ((calls.load(.monotonic) > before) == want) return;
    check.say("{s}: Accelerate {s}called: ", .{ what, if (want) "not " else "" });
    return error.TestUnexpectedResult;
}

/// Counts calls into Accelerate routine by routine.
const Probe = struct {
    on: bool,
    last: usize,

    fn init(on_: bool) Probe {
        return .{ .on = on_, .last = calls.load(.monotonic) };
    }

    /// The routine just run called Accelerate if `routed` and `on`.
    fn saw(p: *Probe, routed: bool, what: []const u8) !void {
        try expectCalled(p.on and routed, p.last, what);
        p.last = calls.load(.monotonic);
    }
};

fn val(comptime T: type, re: Float(T), im: Float(T)) T {
    return if (comptime isComplex(T)) .init(re, im) else re;
}

/// Each public routine over `T` once, on small finite data: each one
/// Accelerate has calls it when `on`, and nothing else does.
fn everyRoutine(comptime T: type, on_: bool) !void {
    const al = zigstd.testAllocator;
    const io = check.io;
    const cx = comptime isComplex(T);
    const real = !cx;
    const al1 = val(T, 1.5, -0.5);
    const be1 = val(T, 0.5, 0.25);
    var a: [36]T = undefined;
    var bm: [36]T = undefined;
    var c: [36]T = undefined;
    var x: [6]T = undefined;
    var y: [6]T = undefined;
    check.fill(T, &a, 1);
    check.fill(T, &bm, 2);
    check.fill(T, &c, 3);
    check.fill(T, &x, 4);
    check.fill(T, &y, 5);
    // dominant diagonals, for the solves
    for (0..6) |i| a[i * 7] = val(T, 4, 0);
    for (0..6) |i| a[i] = val(T, 4, 0);
    const n: isize = 6;
    var p: Probe = .init(on_);

    try blas.copy(io, T, n, &x, 1, &y, 1);
    try p.saw(true, "copy");
    try blas.swap(io, T, n, &x, 1, &y, -1);
    try p.saw(true, "swap");
    try blas.scal(io, T, n, al1, &x, 1);
    try p.saw(true, "scal");
    _ = try blas.dot(io, T, n, &x, 1, &y, 1);
    try p.saw(!(x86_64 and T == f32), "dot");
    if (cx) {
        _ = try blas.dotc(io, T, n, &x, 1, &y, 1);
        try p.saw(true, "dotc");
        try blas.axpyc(io, T, n, al1, &x, 1, &y, 1);
        try p.saw(false, "axpyc");
        // `alpha` has a zero imaginary part (`quirks`)
        try blas.scalReal(io, T, n, 0.5, &x, 1);
        try p.saw(false, "scalReal");
    } else {
        var param = [5]T{ -1, 0.5, 0.25, -0.25, 0.75 };
        try blas.rotm(io, T, n, &x, 1, &y, 1, &param);
        try p.saw(true, "rotm");
    }
    if (T == f32) {
        _ = try blas.dsdot(io, n, &x, 1, &y, 1);
        try p.saw(!x86_64, "dsdot");
        _ = try blas.sdsdot(io, n, 0.5, &x, 1, &y, 1);
        try p.saw(!x86_64, "sdsdot");
    }
    _ = try blas.nrm2(io, T, n, &x, 1);
    try p.saw(!(x86_64 and Float(T) == f32), "nrm2");
    _ = try blas.nrm2(io, T, n, &x, -1);
    try p.saw(Float(T) == f32 and !x86_64, "nrm2, a negative stride");
    _ = try blas.asum(io, T, n, &x, 1);
    try p.saw(true, "asum");
    _ = try blas.sum(io, T, n, &x, 1);
    try p.saw(false, "sum");
    _ = try blas.index(io, T, .amax, n, &x, 1);
    try p.saw(false, "index amax");
    _ = try blas.index(io, T, .amin, n, &x, 1);
    try p.saw(false, "index amin");
    _ = try blas.extreme(io, T, .amax, n, &x, 1);
    try p.saw(false, "extreme amax");
    try blas.rot(io, T, n, &x, 1, &y, 1, 0.6, 0.8);
    try p.saw(true, "rot");
    try blas.axpy(io, T, n, al1, &x, 1, &y, 1);
    try p.saw(true, "axpy");
    try blas.axpby(io, T, n, al1, &x, 1, be1, &y, 1);
    try p.saw(true, "axpby");

    try blas.gemv(al, io, T, .n, n, n, al1, &a, n, &x, 1, be1, &y, 1);
    try p.saw(true, "gemv");
    try blas.gemv(al, io, T, if (cx) .c else .t, n, n, al1, &a, n, &x, 1, be1, &y, 1);
    try p.saw(true, "gemv transposed");
    try blas.gemv(al, io, T, .r, n, n, al1, &a, n, &x, 1, be1, &y, 1);
    try p.saw(real, "gemv r");
    try blas.gbmv(al, io, T, .t, n, n, 1, 1, al1, &a, 3, &x, 1, be1, &y, 1);
    try p.saw(true, "gbmv");
    try blas.ger(al, io, T, n, n, al1, &x, 1, &y, 1, &c, n);
    try p.saw(true, "ger");
    if (cx) {
        try blas.gerc(al, io, T, n, n, al1, &x, 1, &y, 1, &c, n);
        try p.saw(true, "gerc");
        try blas.hemv(al, io, T, .l, n, al1, &a, n, &x, 1, be1, &y, 1);
        try p.saw(true, "hemv");
        try blas.hemv(al, io, T, .v, n, al1, &a, n, &x, 1, be1, &y, 1);
        try p.saw(false, "hemv v");
        try blas.hbmv(al, io, T, .u, n, 1, al1, &a, 2, &x, 1, be1, &y, 1);
        try p.saw(true, "hbmv");
        try blas.hpmv(al, io, T, .l, n, al1, &c, &x, 1, be1, &y, 1);
        try p.saw(true, "hpmv");
        try blas.her(al, io, T, .u, n, 0.5, &x, 1, &c, n);
        try p.saw(true, "her");
        try blas.hpr(al, io, T, .l, n, 0.5, &x, 1, &c);
        try p.saw(true, "hpr");
        try blas.her2(al, io, T, .u, n, al1, &x, 1, &y, 1, &c, n);
        try p.saw(true, "her2");
        try blas.hpr2(al, io, T, .l, n, al1, &x, 1, &y, 1, &c);
        try p.saw(true, "hpr2");
    }
    try blas.symv(al, io, T, .u, n, al1, &a, n, &x, 1, be1, &y, 1);
    try p.saw(real, "symv");
    try blas.sbmv(al, io, T, .l, n, 1, al1, &a, 2, &x, 1, be1, &y, 1);
    try p.saw(real, "sbmv");
    try blas.spmv(al, io, T, .u, n, al1, &c, &x, 1, be1, &y, 1);
    try p.saw(real, "spmv");
    try blas.syr(al, io, T, .u, n, al1, &x, 1, &c, n);
    try p.saw(real, "syr");
    try blas.spr(al, io, T, .l, n, al1, &x, 1, &c);
    try p.saw(real, "spr");
    try blas.syr2(al, io, T, .u, n, al1, &x, 1, &y, 1, &c, n);
    try p.saw(real, "syr2");
    try blas.spr2(al, io, T, .l, n, al1, &x, 1, &y, 1, &c);
    try p.saw(real, "spr2");
    try blas.trmv(al, io, T, .u, .t, .n, n, &a, n, &x, 1);
    try p.saw(true, "trmv");
    try blas.trmv(al, io, T, .u, .r, .n, n, &a, n, &x, 1);
    try p.saw(real, "trmv r");
    try blas.trsv(al, io, T, .l, .n, .u, n, &a, n, &x, -1);
    try p.saw(true, "trsv");
    try blas.tbmv(al, io, T, .u, .n, .n, n, 1, &a, 2, &x, 1);
    try p.saw(true, "tbmv");
    try blas.tbsv(al, io, T, .l, .c, .n, n, 1, &a, 2, &x, 1);
    try p.saw(true, "tbsv");
    try blas.tpmv(al, io, T, .u, .n, .u, n, &a, &x, 1);
    try p.saw(true, "tpmv");
    try blas.tpsv(al, io, T, .l, .n, .u, n, &a, &x, 1);
    try p.saw(true, "tpsv");

    try blas.gemm(al, io, T, .n, if (cx) .c else .t, n, n, n, al1, &a, n, &bm, n, be1, &c, n);
    try p.saw(true, "gemm");
    try blas.gemm(al, io, T, .r, .n, n, n, n, al1, &a, n, &bm, n, be1, &c, n);
    try p.saw(real, "gemm r");
    try blas.symm(al, io, T, .l, .u, n, n, al1, &a, n, &bm, n, be1, &c, n);
    try p.saw(true, "symm");
    try blas.syrk(al, io, T, .u, .n, n, n, al1, &a, n, be1, &c, n);
    try p.saw(true, "syrk");
    try blas.syr2k(al, io, T, .l, .t, n, n, al1, &a, n, &bm, n, be1, &c, n);
    try p.saw(true, "syr2k");
    if (cx) {
        try blas.hemm(al, io, T, .r, .l, n, n, al1, &a, n, &bm, n, be1, &c, n);
        try p.saw(true, "hemm");
        try blas.herk(al, io, T, .l, .c, n, n, 1.5, &a, n, 0.5, &c, n);
        try p.saw(true, "herk");
        try blas.her2k(al, io, T, .u, .n, n, n, al1, &a, n, &bm, n, 0.5, &c, n);
        try p.saw(true, "her2k");
    }
    try blas.trmm(al, io, T, .r, .u, .t, .n, n, n, al1, &a, n, &c, n);
    try p.saw(true, "trmm");
    try blas.trsm(al, io, T, .l, .l, .n, .u, n, n, al1, &a, n, &c, n);
    try p.saw(true, "trsm");
}

test "the public routines call Accelerate on macOS where it has them, for s, d, c and z only" {
    _ = blas.vtable();
    const was = routing();
    defer _ = route(was);
    if (!route(true)) return zigstd.testExpectTrue(!available or loaded == null);
    inline for (.{ f32, f64, Complex32, Complex64 }) |T| try everyRoutine(T, true);
    // types it has nothing for
    var h: [4]f16 = @splat(1);
    const before = calls.load(.monotonic);
    _ = try blas.dot(check.io, f16, 4, &h, 1, &h, 1);
    inline for (check.x_real) |X| {
        var q: [4]X = @splat(1);
        _ = try blas.dot(check.io, X, 4, &q, 1, &q, 1);
    }
    try expectCalled(false, before, "f16 and XDouble dot");
}

test "the off switch runs rstd's kernels, and RSTD_BLAS_ACCELERATE=0 sets it at init" {
    _ = blas.vtable();
    const was = routing();
    defer _ = route(was);
    try zigstd.testExpectTrue(!blas.useAccelerate(false));
    inline for (.{ f32, f64, Complex32, Complex64 }) |T| try everyRoutine(T, false);
    if (!blas.useAccelerate(true)) return error.SkipZigTest;
    try everyRoutine(f64, true);

    blas.shutdown();
    defer {
        blas.shutdown();
        _ = blas.vtable();
    }
    var env: std.process.Environ.Map = .init(zigstd.testAllocator);
    defer env.deinit();
    try env.put("RSTD_BLAS_ACCELERATE", "0");
    try zigstd.testExpectTrue(blas.init(&env));
    try zigstd.testExpectTrue(!routing());
    try everyRoutine(Complex32, false);
    // loaded before, and never unloaded
    try zigstd.testExpectTrue(blas.useAccelerate(true));
    try everyRoutine(Complex32, true);
}

fn dotOnThread(out: *?Threading, before: usize, called: *bool) void {
    var x = [_]f64{ 1, 2, 3 };
    _ = blas.dot(check.io, f64, 3, &x, 1, &x, 1) catch {};
    called.* = calls.load(.monotonic) > before;
    out.* = loaded.?.threading();
}

test "each thread that calls in runs Accelerate single-threaded" {
    if (sysinfo.single_threaded) return error.SkipZigTest;
    _ = blas.vtable();
    if (!routing()) return error.SkipZigTest;
    var model: ?Threading = .multi_threaded;
    var called = false;
    const th = try std.Thread.spawn(.{}, dotOnThread, .{ &model, calls.load(.monotonic), &called });
    th.join();
    try zigstd.testExpectTrue(called);
    // null before macOS 15, where only `VECLIB_MAXIMUM_THREADS` does it
    if (model) |m| try zigstd.testExpectEqual(Threading.single_threaded, m);
}

test "on the 32-bit table, what does not fit c_int runs rstd's kernels" {
    // sizes past `c_int` need a 64-bit `isize`, as on every macOS
    if (comptime @bitSizeOf(isize) <= 32) return error.SkipZigTest;
    _ = blas.vtable();
    if (!routing()) return error.SkipZigTest;
    const lp = load(.{ .lp64 = true }) orelse return error.SkipZigTest;
    try zigstd.testExpectTrue(lp.blas == .lp64);
    const saved = loaded;
    loaded = lp;
    defer loaded = saved;
    const io = check.io;
    const al = zigstd.testAllocator;
    const big: isize = std.math.maxInt(c_int) + 1;
    var x = [_]f64{ 1, -2, 3 };
    var y = [_]f64{ 0.5, 4, -1 };

    var before = calls.load(.monotonic);
    try zigstd.testExpectEqual(@as(f64, -10.5), try blas.dot(io, f64, 3, &x, 1, &y, 1));
    try expectCalled(true, before, "dot");
    try zigstd.testExpectEqual(@as(isize, 3), try blas.index(io, f64, .amax, 3, &x, 1));
    // one element at a stride past `c_int`
    before = calls.load(.monotonic);
    try zigstd.testExpectEqual(@as(f64, 0.5), try blas.dot(io, f64, 1, &x, big, &y, 1));
    try zigstd.testExpectEqual(@as(f64, 1), try blas.asum(io, f64, 1, &x, big));
    try zigstd.testExpectEqual(@as(isize, 1), try blas.index(io, f64, .amax, 1, &x, big));
    // one column at a leading dimension past it
    var c = [_]f64{ 1, 1, 1 };
    try blas.gemv(al, io, f64, .n, 3, 1, 2, &x, big, &y, 1, 1, &c, 1);
    try zigstd.testExpectEqualSlices(f64, &.{ 2, -1, 4 }, &c);
    try blas.gemm(al, io, f64, .n, .n, 3, 1, 1, 1, &x, 3, &y, big, 0, &c, big);
    try zigstd.testExpectEqualSlices(f64, &.{ 0.5, -1, 1.5 }, &c);
    try expectCalled(false, before, "past c_int");
    // an operand whose span is past `c_int`, though each argument fits
    const half: isize = std.math.maxInt(c_int) / 2 + 1;
    try zigstd.testExpectTrue(call("cblas_ddot", &loaded.?.blas, &.{ span(3, half), 3 }, .{ @as(isize, 3), &x, half, &y, @as(isize, 1) }) == null);
    try zigstd.testExpectEqual(@as(f64, -10.5), call("cblas_ddot", &loaded.?.blas, &.{ 3, 3 }, .{ @as(isize, 3), &x, @as(isize, 1), &y, @as(isize, 1) }).?);
}

// ---- Accelerate against rstd's kernels on special values ----------------------

/// The routed level-1 routines.
const R1 = enum { copy, swap, scal, dot, dotc, dsdot, nrm2, asum, iamax, rot, rotm, axpy, axpby };

/// `alpha`, `beta` (`scal`, `axpy`, `axpby`), `c`, `s` (`rot`) and `rotm`'s
/// flag, by variant.
const scalars1 = [_]struct { a: [2]f64, b: [2]f64, cs: [2]f64, flag: f64 }{
    .{ .a = .{ 1.5, -0.5 }, .b = .{ 0.5, 0.25 }, .cs = .{ 0.6, 0.8 }, .flag = -1 },
    .{ .a = .{ 0, 0 }, .b = .{ 0.5, 0 }, .cs = .{ 0, 1 }, .flag = 0 },
    .{ .a = .{ 1.5, 0 }, .b = .{ 0, 0 }, .cs = .{ 1, 0 }, .flag = 1 },
    .{ .a = .{ 0, 0 }, .b = .{ 0, 0 }, .cs = .{ -0.8, 0.6 }, .flag = -2 },
    .{ .a = .{ 1, 0 }, .b = .{ 1, 0 }, .cs = .{ 0.6, -0.8 }, .flag = -1 },
    .{ .a = .{ 0, 1 }, .b = .{ 0, -0.5 }, .cs = .{ 0.8, 0.6 }, .flag = 0 },
};

fn variants1(r: R1) usize {
    return switch (r) {
        .scal, .rot, .rotm, .axpy, .axpby => scalars1.len,
        else => 1,
    };
}

fn has1(comptime T: type, r: R1) bool {
    return switch (r) {
        .dotc => isComplex(T),
        .dsdot => T == f32,
        .rotm => !isComplex(T),
        else => true,
    };
}

/// A level-1 run's value or index.
const Out1 = struct { v: [2]f128 = .{ 0, 0 }, i: isize = 0 };

fn out1(v: anytype) Out1 {
    if (comptime isComplex(@TypeOf(v))) return .{ .v = .{ v.re, v.im } };
    return .{ .v = .{ v, 0 } };
}

fn run1(comptime T: type, comptime r: R1, k: usize, n: usize, x: *check.Vector(T), y: *check.Vector(T)) !Out1 {
    const io = check.io;
    const ni: isize = @intCast(n);
    const s = scalars1[k];
    const F = Float(T);
    const alpha = val(T, @floatCast(s.a[0]), @floatCast(s.a[1]));
    const beta = val(T, @floatCast(s.b[0]), @floatCast(s.b[1]));
    const xb = x.base();
    const yb = y.base();
    switch (r) {
        .copy => try blas.copy(io, T, ni, xb, x.inc, yb, y.inc),
        .swap => try blas.swap(io, T, ni, xb, x.inc, yb, y.inc),
        .scal => try blas.scal(io, T, ni, alpha, xb, x.inc),
        .dot => return out1(try blas.dot(io, T, ni, xb, x.inc, yb, y.inc)),
        .dotc => if (comptime isComplex(T)) return out1(try blas.dotc(io, T, ni, xb, x.inc, yb, y.inc)),
        .dsdot => if (comptime T == f32) return out1(try blas.dsdot(io, ni, xb, x.inc, yb, y.inc)),
        .nrm2 => return out1(try blas.nrm2(io, T, ni, xb, x.inc)),
        .asum => return out1(try blas.asum(io, T, ni, xb, x.inc)),
        .iamax => return .{ .i = try blas.index(io, T, .amax, ni, xb, x.inc) },
        .rot => try blas.rot(io, T, ni, xb, x.inc, yb, y.inc, @floatCast(s.cs[0]), @floatCast(s.cs[1])),
        .rotm => if (comptime !isComplex(T)) {
            const param = [5]F{ @floatCast(s.flag), 0.75, -0.5, 0.25, 1.25 };
            try blas.rotm(io, T, ni, xb, x.inc, yb, y.inc, &param);
        },
        .axpy => try blas.axpy(io, T, ni, alpha, xb, x.inc, yb, y.inc),
        .axpby => try blas.axpby(io, T, ni, alpha, xb, x.inc, beta, yb, y.inc),
    }
    return .{};
}

/// `a` and `b` alike: NaN for NaN, the same infinity, else within `ulps` of
/// `scale`; with `lax`, any non-finite for any (complex infinities meet in
/// orders of their own).
fn alike(comptime F: type, a: f128, b: f128, ulps: f128, scale: f128, lax: bool) bool {
    if (!std.math.isFinite(a) or !std.math.isFinite(b)) {
        if (lax and !std.math.isFinite(a) and !std.math.isFinite(b)) return true;
        if (std.math.isNan(a) or std.math.isNan(b)) return std.math.isNan(a) and std.math.isNan(b);
        return a == b;
    }
    return @abs(a - b) <= ulps * (std.math.floatEps(F) * scale + std.math.floatTrueMin(F));
}

/// Element `i` of `v` as the routine read it (a zero stride reads one).
fn input(comptime T: type, v: *const check.Vector(T), i: usize) T {
    return if (v.inc == 0) v.was[v.n - 1] else v.was[i];
}

/// The case's results against rstd's: values within their accuracy, an
/// index exactly, vectors element by element, gaps untouched.
fn compare1(comptime T: type, comptime r: R1, n: usize, want: Out1, got: Out1, x1: *const check.Vector(T), y1: *const check.Vector(T), x2: *const check.Vector(T), y2: *const check.Vector(T)) !void {
    const F = Float(T);
    const cx = comptime isComplex(T);
    const nf: f128 = @floatFromInt(n);
    switch (r) {
        .iamax => if (want.i != got.i) {
            check.say("index {d}, rstd {d}: ", .{ got.i, want.i });
            return error.TestUnexpectedResult;
        },
        .dot, .dotc, .dsdot, .nrm2, .asum => {
            var scale: f128 = 0;
            for (0..n) |i| scale += switch (r) {
                .nrm2, .asum => check.mag(input(T, x1, i)),
                else => check.mag(input(T, x1, i)) * check.mag(input(T, y1, i)),
            };
            if (r == .nrm2) scale = check.mag(want.v[0]);
            // `sdot`, `snrm2` and `scnrm2` sum in `f64` and round once, as rstd does
            const single = (T == f32 and r == .dot) or (F == f32 and r == .nrm2);
            const ulps: f128 = if (single) 2 else nf + 4;
            const scale2 = if (single) check.mag(want.v[0]) + scale * 0x1p-40 else scale;
            const RF = if (r == .dsdot) f64 else F;
            if (!alike(RF, want.v[0], got.v[0], ulps, scale2, cx) or !alike(RF, want.v[1], got.v[1], ulps, scale2, cx)) {
                check.say("{e},{e}, rstd {e},{e}: ", .{ got.v[0], got.v[1], want.v[0], want.v[1] });
                return error.TestUnexpectedResult;
            }
        },
        else => {},
    }
    for ([_][2]*const check.Vector(T){ .{ x1, x2 }, .{ y1, y2 } }, 0..) |pair, which| {
        if (pair[0].inc == 0) continue;
        for (0..n) |i| {
            const scale = 4 * (check.mag(input(T, x1, i)) + check.mag(input(T, y1, i)));
            const w = pair[0].at(i);
            const g = pair[1].at(i);
            const ok = if (comptime cx)
                alike(F, w.re, g.re, 8, scale, true) and alike(F, w.im, g.im, 8, scale, true)
            else
                alike(F, w, g, 8, scale, false);
            if (!ok) {
                check.say("{s}[{d}] {any}, rstd {any}: ", .{ if (which == 0) "x" else "y", i, g, w });
                return error.TestUnexpectedResult;
            }
        }
        try pair[1].gaps();
    }
}

/// Room for four vectors of 1025 at stride 3, and their copies as they were.
var mem1: [4 * (2 * 1025 * 3 + 64) * @sizeOf(Complex64)]u8 align(16) = undefined;

/// One case through rstd's kernels and through Accelerate, compared.
fn case1(comptime T: type, comptime r: R1, k: usize, n: usize, inc: [2]isize, p: ?check.Pattern, on_y: bool) !void {
    var fba: std.heap.FixedBufferAllocator = .init(&mem1);
    const al = fba.allocator();
    const px = if (on_y) null else p;
    const py = if (on_y) p else null;
    var x1: check.Vector(T) = try .init(al, n, inc[0], 7 * n + 1, px);
    var y1: check.Vector(T) = try .init(al, n, inc[1], 7 * n + 2, py);
    var x2: check.Vector(T) = try .init(al, n, inc[0], 7 * n + 1, px);
    var y2: check.Vector(T) = try .init(al, n, inc[1], 7 * n + 2, py);
    _ = route(false);
    const want = try run1(T, r, k, n, &x1, &y1);
    _ = route(true);
    const got = try run1(T, r, k, n, &x2, &y2);
    try compare1(T, r, n, want, got, &x1, &y1, &x2, &y2);
}

/// Distinct failures seen, so one run reports each once.
const Seen = struct {
    keys: [256][96]u8 = undefined,
    lens: [256]usize = undefined,
    len: usize = 0,
    failed: bool = false,

    /// Records `key`; true the first time.
    fn first(s: *Seen, key: []const u8) bool {
        s.failed = true;
        for (s.keys[0..s.len], s.lens[0..s.len]) |*k, l| if (std.mem.eql(u8, k[0..l], key)) return false;
        if (s.len == s.keys.len) return false;
        const l = @min(key.len, 96);
        @memcpy(s.keys[s.len][0..l], key[0..l]);
        s.lens[s.len] = l;
        s.len += 1;
        return true;
    }
};

fn patternName(p: ?check.Pattern) []const u8 {
    const q = p orelse return "plain";
    return switch (q) {
        .one => |o| @tagName(o.kind),
        else => @tagName(q),
    };
}

/// Every size, stride and pattern of `r` over `T`; a failure is reported once
/// per variant, pattern kind and stride sign.
fn sweep1(comptime T: type, comptime r: R1, seen: *Seen) void {
    const two = switch (r) {
        .nrm2, .asum, .iamax, .scal => false,
        else => true,
    };
    const strides = [_][2]isize{ .{ 1, 1 }, .{ 2, 3 }, .{ -1, -1 }, .{ -3, 2 }, .{ 2, -3 }, .{ 0, 1 }, .{ 1, 0 } };
    const reads_x = r == .dot or r == .dotc or r == .dsdot or r == .copy or r == .axpy or r == .axpby;
    const reads_y = r == .dot or r == .dotc or r == .dsdot;
    for (0..variants1(r)) |k| for (strides) |inc| {
        // zero strides only where the vector is read, not written
        if (inc[0] == 0 and !reads_x) continue;
        if (inc[1] == 0 and !reads_y) continue;
        if (!two and inc[1] != 1) continue;
        for ([_]usize{ 1, 2, 5, 16, 33, 67, 300, 1025 }) |n| {
            var buf: [13]check.Pattern = undefined;
            const ps = check.patterns(n, &buf);
            for (0..ps.len + 1) |pi| for ([_]bool{ false, true }) |on_y| {
                const p: ?check.Pattern = if (pi == ps.len) null else ps[pi];
                if (on_y and (!two or p == null)) continue;
                case1(T, r, k, n, inc, p, on_y) catch {
                    var kb: [96]u8 = undefined;
                    const key = std.fmt.bufPrint(&kb, "{s} {s} v{d} {s} in {s} inc {d},{d}", .{ @tagName(r), @typeName(T), k, patternName(p), if (on_y) "y" else "x", std.math.sign(inc[0]), std.math.sign(inc[1]) }) catch "?";
                    if (seen.first(key)) check.say("n={d} inc={d},{d}: {s}\n", .{ n, inc[0], inc[1], key });
                };
            };
        }
    };
}

test "Accelerate's level-1 results match rstd's kernels on special values, strides and sizes" {
    _ = blas.vtable();
    const was = routing();
    defer _ = route(was);
    if (!route(true)) return error.SkipZigTest;
    var seen: Seen = .{};
    inline for (.{ f32, f64, Complex32, Complex64 }) |T| {
        inline for (comptime std.enums.values(R1)) |r| if (comptime has1(T, r)) sweep1(T, r, &seen);
    }
    if (seen.failed) return error.TestUnexpectedResult;
}

test "quirk cases run rstd's kernels, and the Hermitian diagonal comes out real" {
    _ = blas.vtable();
    if (!routing()) return error.SkipZigTest;
    const al = zigstd.testAllocator;
    const io = check.io;
    const inf = std.math.inf(f64);
    var z: [4]Complex64 = .{ .init(inf, 1), .init(1, 2), .init(-1, 0.5), .init(2, -1) };
    var w: [4]Complex64 = .{ .init(1, 1), .init(1, -1), .init(0.5, 0.5), .init(-2, 1) };
    var x = [_]f64{ 1, inf, -1, 2 };
    var y = [_]f64{ 1, 0, 2, 0.5 };
    var a: [16]f64 = @splat(1);
    for (0..4) |i| a[i * 5] = 2;
    var b: [16]f64 = @splat(0.5);
    b[5] = inf;
    var p: Probe = .init(true);

    try blas.scal(io, Complex64, 4, .init(1.5, 0), &z, 1);
    try p.saw(false, "complex scal, alpha with a zero part");
    try zigstd.testExpectTrue(std.math.isNan(z[0].im));
    try blas.axpy(io, Complex64, 4, .init(0, 1), &z, 1, &w, 1);
    try p.saw(false, "complex axpy, alpha with a zero part");
    try blas.axpby(io, f64, 4, 0, &x, 1, 2, &y, 1);
    try p.saw(false, "axpby, zero alpha");
    try zigstd.testExpectEqualSlices(f64, &.{ 2, 0, 4, 1 }, &y);
    try blas.axpby(io, Complex64, 4, .init(1, 1), &z, 1, .init(2, 0), &w, 1);
    try p.saw(false, "complex axpby, beta with a zero part");
    try blas.trmv(al, io, f64, .u, .n, .n, 4, &a, 4, &x, 1);
    try p.saw(false, "trmv, x not finite");
    try blas.ger(al, io, f64, 4, 4, 1, &x, 1, &y, 1, &a, 4);
    try p.saw(false, "ger, x not finite");
    try blas.trmm(al, io, f64, .l, .l, .n, .n, 4, 4, 1, &a, 4, &b, 4);
    try p.saw(false, "trmm, B not finite");
    // the same cases, finite, call it
    x = .{ 1, 3, -1, 2 };
    @memset(&a, 1);
    try blas.trmv(al, io, f64, .u, .n, .n, 4, &a, 4, &x, 1);
    try p.saw(true, "trmv");
    @memset(&b, 0.5);
    try blas.trmm(al, io, f64, .l, .l, .n, .n, 4, 4, 1, &a, 4, &b, 4);
    try p.saw(true, "trmm");

    // `herk` and `her2k` make the diagonal real; the reference never reads it
    var h: [12]Complex64 = undefined;
    var c: [9]Complex64 = undefined;
    for ([_]bool{ false, true }) |two| for ([_]bool{ false, true }) |special| {
        check.fill(Complex64, &h, 3);
        if (special) h[4] = .init(inf, 1);
        check.fill(Complex64, &c, 4);
        for (0..3) |i| c[i * 4].im = 3;
        if (two) {
            try blas.her2k(al, io, Complex64, .u, .n, 3, 4, .init(0.5, 0.25), &h, 3, &h, 3, 0.5, &c, 3);
        } else {
            try blas.herk(al, io, Complex64, .l, .n, 3, 4, 1.5, &h, 3, 0.5, &c, 3);
        }
        try p.saw(true, "herk and her2k");
        for (0..3) |i| try zigstd.testExpectEqual(@as(f64, 0), c[i * 4].im);
    };
    // and `her2`, packed too
    var hx: [3]Complex64 = .{ .init(inf, 1), .init(1, 2), .init(-1, 0.5) };
    check.fill(Complex64, &c, 5);
    try blas.her2(al, io, Complex64, .l, 3, .init(0.5, 0.25), &hx, 1, &hx, 1, &c, 3);
    try p.saw(true, "her2");
    for (0..3) |i| try zigstd.testExpectEqual(@as(f64, 0), c[i * 4].im);
    check.fill(Complex64, &c, 5);
    try blas.hpr2(al, io, Complex64, .u, 3, .init(0.5, 0.25), &hx, 1, &hx, 1, &c);
    try p.saw(true, "hpr2");
    for ([_]usize{ 0, 2, 5 }) |i| try zigstd.testExpectEqual(@as(f64, 0), c[i].im);
    // `beta` 1 and nothing added leaves `C` as it was, through rstd
    check.fill(Complex64, &c, 4);
    for (0..3) |i| c[i * 4].im = 3;
    try blas.herk(al, io, Complex64, .l, .n, 3, 0, 1.5, &h, 3, 1, &c, 3);
    try p.saw(false, "herk, beta 1 and k 0");
    for (0..3) |i| try zigstd.testExpectEqual(@as(f64, 3), c[i * 4].im);
    // a non-finite `C` with a nonzero `beta`, through rstd; with a zero one, Accelerate
    check.fill(Complex64, &h, 3);
    check.fill(Complex64, &c, 4);
    c[1] = .init(1, inf);
    try blas.herk(al, io, Complex64, .l, .n, 3, 4, 1.5, &h, 3, 0.5, &c, 3);
    try p.saw(false, "herk, C not finite");
    try zigstd.testExpectTrue(!std.math.isNan(c[1].re));
    c[1] = .init(1, inf);
    try blas.herk(al, io, Complex64, .l, .n, 3, 4, 1.5, &h, 3, 0, &c, 3);
    try p.saw(true, "herk, beta 0");
}
