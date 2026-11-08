//! Apple Accelerate's CBLAS, bound at run time through `dll.zig` with no SDK,
//! and the routing to it (design.md, Apple Accelerate). Sources:
//! <https://developer.apple.com/documentation/macos-release-notes/macos-13_3-release-notes>,
//! <https://developer.apple.com/documentation/accelerate/blas-library>,
//! <https://developer.apple.com/documentation/accelerate/blassetthreading(_:)>.

const std = @import("std");
const sysinfo = @import("builtin");
const dll = @import("../../dll.zig");
const types = @import("types.zig");
const blas = @import("__root__.zig");
const threading = @import("threading.zig");
const portable = @import("portable.zig");
const util = @import("util.zig");

const Io = std.Io;

const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const isComplex = types.isComplex;
const Float = types.Float;
const isZero = util.isZero;

/// Where `load` can find Accelerate.
pub const available = sysinfo.target.os.tag == .macos;
/// x86_64 macOS, whose Accelerate has quirks of its own (design.md).
pub const x86_64 = sysinfo.target.cpu.arch == .x86_64;

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

// ---- routing -----------------------------------------------------------------
//
// Called after the public routines' own checks, so the errors are rstd's;
// false (or null) where rstd's kernels run (design.md, Apple Accelerate).

/// Whether the public routines over `T` can go to Accelerate.
pub fn covers(comptime T: type) bool {
    return available and (T == f32 or T == f64 or T == Complex32 or T == Complex64);
}

/// The framework, loaded once by `start` and never unloaded.
pub var loaded: ?AppleAccelerate = null;
var tried = false;
/// Whether the public routines call it.
var on = std.atomic.Value(bool).init(false);
/// Whether this thread's calls are single-threaded (`singleThreadHere`).
threadlocal var thread_set = false;
/// Calls made into Accelerate, counted in tests only.
pub var calls = std.atomic.Value(usize).init(0);
/// Calls on a thread not single-threaded (macOS 15+), counted in tests only.
pub var unpinned = std.atomic.Value(usize).init(0);

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

/// The table for `T`'s routines, or null where rstd runs them.
inline fn use(comptime T: type) ?*const Blas {
    if (comptime !covers(T)) return null;
    _ = blas.vtable();
    if (!on.load(.acquire)) return null;
    return &loaded.?.blas;
}

/// Makes the calling thread's calls single-threaded, on its first.
fn pin() void {
    if (thread_set) return;
    thread_set = true;
    if (loaded) |*acc| acc.singleThreadHere();
}

/// What `call` returns for `name`: its C result, an index as `isize`.
fn Ret(comptime name: []const u8) type {
    const R = @typeInfo(@typeInfo(@typeInfo(@FieldType(Table(i64), name)).optional.child).pointer.child).@"fn".return_type.?;
    return if (R == i64) isize else R;
}

/// Calls `name` with `args`, every `isize` narrowed to the table's integer;
/// null if the table lacks it or, on LP64, any argument or `extents` entry
/// (an operand's span in elements) is past `c_int`. Pins the calling thread
/// (`pin`).
pub fn call(comptime name: []const u8, b: *const Blas, extents: []const isize, args: anytype) ?Ret(name) {
    switch (b.*) {
        inline else => |*bound, tag| {
            const f = @field(bound.fns, name) orelse return null;
            const I = if (tag == .ilp64) i64 else c_int;
            if (I != i64) for (extents) |e| if (e > std.math.maxInt(I)) return null;
            var conv: std.meta.ArgsTuple(@typeInfo(@TypeOf(f)).pointer.child) = undefined;
            inline for (0..args.len) |i| {
                conv[i] = if (@TypeOf(args[i]) == isize) std.math.cast(I, args[i]) orelse return null else args[i];
            }
            pin();
            if (sysinfo.is_test) {
                _ = calls.fetchAdd(1, .monotonic);
                if (loaded) |*acc| if (acc.threading()) |t| if (t != .single_threaded) {
                    _ = unpinned.fetchAdd(1, .monotonic);
                };
            }
            const r = @call(.auto, f, conv);
            return if (@TypeOf(r) == I) @intCast(r) else r;
        },
    }
}

/// The span of `n` elements at stride `inc`.
pub fn span(n: isize, inc: isize) isize {
    return (n - 1) * @as(isize, @intCast(@abs(inc))) + 1;
}

/// The span of a `rows` by `cols` matrix at leading dimension `ld`.
fn mat(rows: isize, cols: isize, ld: isize) isize {
    return ld * (cols - 1) + rows;
}

// ---- threads -----------------------------------------------------------------
//
// Accelerate runs single-threaded; a call splits by size alone into parts on
// rstd's workers, one Accelerate call each (threading.md, Vendor backends).

/// Level 1's parts start on multiples of this many elements, as rstd's own
/// do (`blas1/interface.zig`).
pub const vector_quantum = 1024;
/// Level 2's and 3's on multiples of this many rows or columns, as rstd's
/// level 2 (`blas2/interface.zig`): parts stay off each other's cache lines.
pub const matrix_quantum = 64;

/// `part(ctx, chunk)` for each of `family`'s chunks of `len` units, `work`
/// being its rule's measure: inline for one, else on `io`'s workers.
fn split(io: Io, comptime family: threading.Family, comptime T: type, work: u64, len: usize, quantum: usize, ctx: anytype, comptime part: fn (@TypeOf(ctx), threading.Chunk) void) Io.Cancelable!void {
    const p = threading.plan(family, T, work, len, .even, quantum);
    if (p.chunks.count == 1) return part(ctx, p.chunks.get(0));
    return threading.forEach(io, &p, ctx, part);
}

/// `split` for a reduction over `len` elements: the parts' sums, added in
/// chunk order.
fn fold(io: Io, comptime family: threading.Family, comptime T: type, comptime W: type, len: usize, ctx: anytype, comptime part: fn (@TypeOf(ctx), threading.Chunk) W) Io.Cancelable!W {
    const p = threading.plan(family, T, len, len, .even, vector_quantum);
    if (p.chunks.count == 1) return part(ctx, p.chunks.get(0));
    const S = struct {
        fn add(a: W, b: W) W {
            return if (comptime isComplex(W)) a.add(b) else a + b;
        }
    };
    return threading.reduce(io, &p, W, ctx, part, S.add);
}

/// What a reduction over `T` sums its parts in: `f64` for `f32`, so the
/// result rounds once, as rstd's own fold (`blas1/interface.zig`).
fn Sum(comptime T: type) type {
    const W = if (Float(T) == f32) f64 else Float(T);
    return if (isComplex(T)) std.math.Complex(W) else W;
}

fn rounded(comptime T: type, w: Sum(T)) T {
    if (comptime isComplex(T)) return .init(@floatCast(w.re), @floatCast(w.im));
    return @floatCast(w);
}

/// Whether the table has `name` and, on LP64, every entry of `ints` (each
/// integer argument's magnitude and each operand's span) fits `c_int`; then
/// so does any call over a part of the operands.
fn fits(comptime name: []const u8, b: *const Blas, ints: []const isize) bool {
    switch (b.*) {
        inline else => |*bound, tag| {
            if (@field(bound.fns, name) == null) return false;
            if (tag == .lp64) for (ints) |e| if (e > std.math.maxInt(c_int)) return false;
            return true;
        },
    }
}

/// `|v|`, as `fits` takes a stride.
fn mag(v: isize) isize {
    return @intCast(@abs(v));
}

/// The whole of `len` units, as one chunk.
fn whole(len: isize) threading.Chunk {
    return .{ .index = 0, .start = 0, .end = @intCast(len) };
}

/// Elements `c.start..c.end` of the `n` at `p`, stride `inc`, as CBLAS takes
/// a vector: from its lowest address, which is the last one for `inc < 0`.
fn sub(p: anytype, n: isize, inc: isize, c: threading.Chunk) @TypeOf(p) {
    const at = if (inc >= 0) c.start else @as(usize, @intCast(n)) - c.end;
    return p + at * @abs(inc);
}

/// Column `j` of a matrix at `p`, leading dimension `ld`.
fn column(p: anytype, j: usize, ld: isize) @TypeOf(p) {
    return p + j * @as(usize, @intCast(ld));
}

fn units(c: threading.Chunk) isize {
    return @intCast(c.end - c.start);
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

// Quirks: where Accelerate's results differ from the reference's, rstd runs
// the call (design.md, Apple Accelerate, lists them).

/// The level-1 routines `sends` covers; `dot` over `f32` is `sdot`.
pub const L1 = enum { copy, swap, scal, dot, dotc, dsdot, nrm2, asum, rot, rotm, axpy, axpby };

/// Which calls go to Accelerate: all, those whose every stride is 1, none.
const Strides = enum { all, unit, none };

/// `Strides` by letter.
const Letters = struct { s: Strides = .all, d: Strides = .all, c: Strides = .all, z: Strides = .all };

/// Where rstd's kernels are the faster, by speed alone: the arm64 flyoff of
/// 2026-10-29 (design.md, Apple Accelerate, has the rule and its rows);
/// strided there is 2. x86_64 has no data and keeps every call.
const sends: std.enums.EnumFieldStruct(L1, Letters, .{}) = if (x86_64) .{} else .{
    .nrm2 = .{ .s = .none, .d = .none, .c = .none, .z = .none },
    .dot = .{ .s = .none, .c = .unit, .z = .unit },
    .dotc = .{ .c = .unit, .z = .unit },
    .copy = .{ .s = .unit, .d = .unit },
    .swap = .{ .s = .unit, .d = .unit, .c = .unit, .z = .unit },
    .axpy = .{ .c = .unit, .z = .unit },
    .rot = .{ .c = .unit },
};

/// Whether `sends` lets `r` over `T` go to Accelerate; `unit` if every
/// stride is 1.
pub inline fn sent(comptime r: L1, comptime T: type, unit: bool) bool {
    return switch (comptime @field(@field(sends, @tagName(r)), portable.letter(T).?)) {
        .all => true,
        .unit => unit,
        .none => false,
    };
}

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

/// Whether some `j` has `x(j)`, and `y(j)` if given, zero: a column the
/// reference's rank updates pass over.
fn zeroColumn(comptime T: type, n: isize, x: [*]const T, incx: isize, y: ?[*]const T, incy: isize) bool {
    const un: usize = @intCast(n);
    for (0..un) |j| {
        if (!isZero(T, x[elem(j, un, incx)])) continue;
        const yv = y orelse return true;
        if (isZero(T, yv[elem(j, un, incy)])) return true;
    }
    return false;
}

/// Element `j`'s offset from the lowest of `n` at stride `inc`.
fn elem(j: usize, n: usize, inc: isize) usize {
    return (if (inc > 0) j else n - 1 - j) * @abs(inc);
}

/// Whether an entry on the diagonal of the order `n` matrix at `a` is too
/// small for its reciprocal: below the smallest normal (for complex, its
/// square root, as `1 / d` divides by `|d|^2`).
fn tinyDiagonal(comptime T: type, n: isize, a: [*]const T, lda: isize) bool {
    const F = Float(T);
    const lim = if (comptime isComplex(T)) @sqrt(std.math.floatMin(F)) else std.math.floatMin(F);
    const un: usize = @intCast(n);
    const ld: usize = @intCast(lda);
    for (0..un) |i| {
        const v = a[i * (ld + 1)];
        const m = if (comptime isComplex(T)) @max(@abs(v.re), @abs(v.im)) else @abs(v);
        if (m < lim) return true;
    }
    return false;
}

/// `nonFinite` of a `rows` by `cols` matrix at leading dimension `ld`.
fn nonFiniteMat(comptime T: type, rows: isize, cols: isize, b: [*]const T, ld: isize) bool {
    for (0..@intCast(cols)) |j| {
        if (nonFinite(T, rows, b + j * @as(usize, @intCast(ld)), 1)) return true;
    }
    return false;
}

// level 1

/// The element-wise level-1 routines `split` takes.
const Elem = enum { copy, swap, scal, rot, axpy };

/// One element-wise call: `x` and, but for `scal`, `y`, with its scalars.
fn Elementwise(comptime T: type, comptime r: Elem) type {
    return struct {
        const Self = @This();
        const name = "cblas_" ++ switch (r) {
            .rot => switch (T) {
                Complex32 => "csrot",
                Complex64 => "zdrot",
                else => portable.letter(T).? ++ "rot",
            },
            else => portable.letter(T).? ++ @tagName(r),
        };

        b: *const Blas,
        n: isize,
        x: [*]T,
        incx: isize,
        y: [*]T = undefined,
        incy: isize = 1,
        alpha: T = undefined,
        c: Float(T) = 0,
        s: Float(T) = 0,

        fn part(e: *const Self, ch: threading.Chunk) void {
            const len = units(ch);
            const x = sub(e.x, e.n, e.incx, ch);
            const y = if (r == .scal) e.y else sub(e.y, e.n, e.incy, ch);
            _ = switch (r) {
                .copy, .swap => call(name, e.b, &.{}, .{ len, x, e.incx, y, e.incy }),
                .scal => call(name, e.b, &.{}, .{ len, scalar(T, &e.alpha), x, e.incx }),
                .rot => call(name, e.b, &.{}, .{ len, x, e.incx, y, e.incy, e.c, e.s }),
                .axpy => call(name, e.b, &.{}, .{ len, scalar(T, &e.alpha), x, e.incx, y, e.incy }),
            } orelse unreachable;
        }

        /// The call over its parts; false where rstd runs it.
        fn run(e: *const Self, io: Io) Io.Cancelable!bool {
            if (!fits(name, e.b, &.{ e.n, mag(e.incx), mag(e.incy), span(e.n, e.incx), span(e.n, e.incy) })) return false;
            if (e.incy == 0 or (e.incx == 0 and (r == .swap or r == .scal or r == .rot))) return false;
            // a zero stride's elements depend on each other, as for rstd's own
            if (e.incx == 0) {
                e.part(whole(e.n));
                return true;
            }
            const len: usize = @intCast(e.n);
            try split(io, @field(threading.Family, @tagName(r)), T, len, len, vector_quantum, e, part);
            return true;
        }
    };
}

pub inline fn copy(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.copy, T, incx == 1 and incy == 1)) return false;
        const b = use(T) orelse return false;
        const e: Elementwise(T, .copy) = .{ .b = b, .n = n, .x = x, .incx = incx, .y = y, .incy = incy };
        return e.run(io);
    }
}

pub inline fn swap(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.swap, T, incx == 1 and incy == 1)) return false;
        const b = use(T) orelse return false;
        const e: Elementwise(T, .swap) = .{ .b = b, .n = n, .x = x, .incx = incx, .y = y, .incy = incy };
        return e.run(io);
    }
}

pub inline fn scal(io: Io, comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.scal, T, incx == 1)) return false;
        const b = use(T) orelse return false;
        if (zeroPart(T, alpha)) return false;
        const e: Elementwise(T, .scal) = .{ .b = b, .n = n, .x = x, .incx = incx, .alpha = alpha };
        return e.run(io);
    }
}

/// The reductions `fold` takes: `dot`, `dotc`, `dsdot` (`f32` only) and
/// `asum` (`x` only, `incx > 0`).
const Red = enum { dot, dotc, dsdot, asum };

/// One reduction, its parts summed in `Sum(T)`.
fn Reduction(comptime T: type, comptime r: Red) type {
    return struct {
        const Self = @This();
        /// `sdot` is `dsdot`: rstd's contract sums in `f64` (readme, One
        /// contract on every core).
        const name = if (r == .asum) switch (T) {
            Complex32 => "cblas_scasum",
            Complex64 => "cblas_dzasum",
            else => "cblas_" ++ portable.letter(T).? ++ "asum",
        } else switch (T) {
            f32 => "cblas_dsdot",
            f64 => "cblas_ddot",
            else => "cblas_" ++ portable.letter(T).? ++ (if (r == .dotc) "dotc" else "dotu") ++ "_sub",
        };
        const W = if (r == .asum) Sum(Float(T)) else Sum(T);

        b: *const Blas,
        n: isize,
        x: [*]T,
        incx: isize,
        y: [*]T = undefined,
        incy: isize = 1,

        fn part(e: *const Self, ch: threading.Chunk) W {
            const len = units(ch);
            const x = sub(e.x, e.n, e.incx, ch);
            if (r == .asum) return call(name, e.b, &.{}, .{ len, x, e.incx }).?;
            const y = sub(e.y, e.n, e.incy, ch);
            if (comptime !isComplex(T)) return call(name, e.b, &.{}, .{ len, x, e.incx, y, e.incy }).?;
            var v: T = undefined;
            _ = call(name, e.b, &.{}, .{ len, x, e.incx, y, e.incy, &v }).?;
            return .init(v.re, v.im);
        }

        /// The sum over its parts; null where rstd runs it.
        fn run(e: *const Self, io: Io) Io.Cancelable!?W {
            if (!fits(name, e.b, &.{ e.n, mag(e.incx), mag(e.incy), span(e.n, e.incx), span(e.n, e.incy) })) return null;
            // a zero stride runs serially, as rstd's own
            if (e.incx == 0 or e.incy == 0) return e.part(whole(e.n));
            return try fold(io, if (r == .asum) .asum else .dot, T, W, @intCast(e.n), e, part);
        }
    };
}

/// `x . y`, or `conj(x) . y`.
pub inline fn dot(io: Io, comptime T: type, comptime conj: bool, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Io.Cancelable!?T {
    if (comptime !covers(T)) return null else {
        if (!sent(if (conj) .dotc else .dot, T, incx == 1 and incy == 1)) return null;
        const b = use(T) orelse return null;
        if (comptime T == f32 and x86_64) return null;
        const e: Reduction(T, if (conj) .dotc else .dot) = .{ .b = b, .n = n, .x = x, .incx = incx, .y = y, .incy = incy };
        return rounded(T, try e.run(io) orelse return null);
    }
}

pub inline fn dsdot(io: Io, n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize) Io.Cancelable!?f64 {
    if (!sent(.dsdot, f32, incx == 1 and incy == 1)) return null;
    const b = use(f32) orelse return null;
    if (comptime x86_64) return null;
    const e: Reduction(f32, .dsdot) = .{ .b = b, .n = n, .x = x, .incx = incx, .y = y, .incy = incy };
    return e.run(io);
}

/// `s` and `c` as rstd's contract has them, the squares summed in `f64` by
/// `dsdot` and the root rounded once; `d` and `z` for `incx > 0`, the
/// reference's only stride.
pub inline fn nrm2(comptime T: type, n: isize, x: [*]T, incx: isize) ?Float(T) {
    if (comptime !covers(T)) return null else {
        if (!sent(.nrm2, T, incx == 1)) return null;
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
pub inline fn asum(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize) Io.Cancelable!?Float(T) {
    if (comptime !covers(T)) return null else {
        if (!sent(.asum, T, incx == 1)) return null;
        const b = use(T) orelse return null;
        if (incx <= 0) return null;
        const e: Reduction(T, .asum) = .{ .b = b, .n = n, .x = x, .incx = incx };
        return rounded(Float(T), try e.run(io) orelse return null);
    }
}

pub inline fn rot(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, c: Float(T), s: Float(T)) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.rot, T, incx == 1 and incy == 1)) return false;
        const b = use(T) orelse return false;
        const e: Elementwise(T, .rot) = .{ .b = b, .n = n, .x = x, .incx = incx, .y = y, .incy = incy, .c = c, .s = s };
        return e.run(io);
    }
}

pub inline fn rotm(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, param: *const [5]T) bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.rotm, T, incx == 1 and incy == 1)) return false;
        const b = use(T) orelse return false;
        if (incx == 0 or incy == 0) return false;
        return call("cblas_" ++ portable.letter(T).? ++ "rotm", b, &.{ span(n, incx), span(n, incy) }, .{ n, x, incx, y, incy, param }) != null;
    }
}

pub inline fn axpy(io: Io, comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.axpy, T, incx == 1 and incy == 1)) return false;
        const b = use(T) orelse return false;
        if (zeroPart(T, alpha)) return false;
        const e: Elementwise(T, .axpy) = .{ .b = b, .n = n, .x = x, .incx = incx, .y = y, .incy = incy, .alpha = alpha };
        return e.run(io);
    }
}

pub inline fn axpby(comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.axpby, T, incx == 1 and incy == 1)) return false;
        const b = use(T) orelse return false;
        if (zeroPart(T, alpha) or zeroPart(T, beta)) return false;
        if (comptime !isComplex(T)) if (alpha == 0) return false;
        if (incy == 0) return false;
        return call("catlas_" ++ portable.letter(T).? ++ "axpby", b, &.{ span(n, incx), span(n, incy) }, .{ n, scalar(T, &alpha), x, incx, scalar(T, &beta), y, incy }) != null;
    }
}

// level 2

/// `gemv`'s parts: rows of `y` without a transpose, columns with one, as
/// `gemv_thread.c` splits them.
fn Gemv(comptime T: type) type {
    return struct {
        const Self = @This();
        const name = "cblas_" ++ portable.letter(T).? ++ "gemv";

        b: *const Blas,
        t: Transpose,
        m: isize,
        n: isize,
        alpha: T,
        a: [*]T,
        lda: isize,
        x: [*]T,
        incx: isize,
        beta: T,
        y: [*]T,
        incy: isize,

        fn part(g: *const Self, ch: threading.Chunk) void {
            const tr = g.t != .no_trans;
            const m = if (tr) g.m else units(ch);
            const n = if (tr) units(ch) else g.n;
            const a = if (tr) column(g.a, ch.start, g.lda) else g.a + ch.start;
            const y = sub(g.y, if (tr) g.n else g.m, g.incy, ch);
            _ = call(name, g.b, &.{}, .{ col, g.t, m, n, scalar(T, &g.alpha), a, g.lda, g.x, g.incx, scalar(T, &g.beta), y, g.incy }).?;
        }
    };
}

pub inline fn gemv(io: Io, comptime T: type, trans: types.Trans, m: isize, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        const t = transpose(T, trans) orelse return false;
        const lenx = if (trans.transposed()) m else n;
        const leny = if (trans.transposed()) n else m;
        if (zeroPart(T, beta) and !isZero(T, beta) and nonFinite(T, leny, y, incy)) return false;
        const G = Gemv(T);
        if (!fits(G.name, b, &.{ m, n, lda, mag(incx), mag(incy), mat(m, n, lda), span(lenx, incx), span(leny, incy) })) return false;
        const g: G = .{ .b = b, .t = t, .m = m, .n = n, .alpha = alpha, .a = a, .lda = lda, .x = x, .incx = incx, .beta = beta, .y = y, .incy = incy };
        try split(io, .gemv, T, threading.product(&.{ @intCast(m), @intCast(n) }), @intCast(leny), matrix_quantum, &g, G.part);
        return true;
    }
}

pub inline fn gbmv(comptime T: type, trans: types.Trans, m: isize, n: isize, kl: isize, ku: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        const t = transpose(T, trans) orelse return false;
        const lenx = if (trans.transposed()) m else n;
        const leny = if (trans.transposed()) n else m;
        const ext = &.{ mat(lda, n, lda), span(lenx, incx), span(leny, incy) };
        return call("cblas_" ++ portable.letter(T).? ++ "gbmv", b, ext, .{ col, t, m, n, kl, ku, scalar(T, &alpha), a, lda, x, incx, scalar(T, &beta), y, incy }) != null;
    }
}

/// `ger`'s parts: columns of `A`, as `ger_thread.c` splits them.
fn Ger(comptime T: type, comptime conj: bool) type {
    return struct {
        const Self = @This();
        const name = "cblas_" ++ portable.letter(T).? ++ "ger" ++ if (!isComplex(T)) "" else if (conj) "c" else "u";

        b: *const Blas,
        m: isize,
        n: isize,
        alpha: T,
        x: [*]T,
        incx: isize,
        y: [*]T,
        incy: isize,
        a: [*]T,
        lda: isize,

        fn part(g: *const Self, ch: threading.Chunk) void {
            _ = call(name, g.b, &.{}, .{ col, g.m, units(ch), scalar(T, &g.alpha), g.x, g.incx, sub(g.y, g.n, g.incy, ch), g.incy, column(g.a, ch.start, g.lda), g.lda }).?;
        }
    };
}

/// `ger`, `geru`, and `gerc` for `conj`.
pub inline fn ger(io: Io, comptime T: type, comptime conj: bool, m: isize, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        if (nonFinite(T, m, x, incx)) return false;
        const G = Ger(T, conj);
        if (!fits(G.name, b, &.{ m, n, lda, mag(incx), mag(incy), span(m, incx), span(n, incy), mat(m, n, lda) })) return false;
        const g: G = .{ .b = b, .m = m, .n = n, .alpha = alpha, .x = x, .incx = incx, .y = y, .incy = incy, .a = a, .lda = lda };
        try split(io, .ger, T, threading.product(&.{ @intCast(m), @intCast(n) }), @intCast(n), matrix_quantum, &g, G.part);
        return true;
    }
}

/// `symv`, or `hemv` for `herm`.
pub inline fn symv(comptime T: type, comptime herm: bool, ul: types.Uplo, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T) or herm != isComplex(T)) return false else {
        const name = comptime "cblas_" ++ portable.letter(T).? ++ if (herm) "hemv" else "symv";
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        return call(name, b, &.{ mat(n, n, lda), span(n, incx), span(n, incy) }, .{ col, u, n, scalar(T, &alpha), a, lda, x, incx, scalar(T, &beta), y, incy }) != null;
    }
}

/// `sbmv` or `hbmv`, packed (`spmv`, `hpmv`) for `packed_`.
pub inline fn sbmv(comptime T: type, comptime herm: bool, comptime packed_: bool, ul: types.Uplo, n: isize, k: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T) or herm != isComplex(T)) return false else {
        const name = comptime "cblas_" ++ portable.letter(T).? ++ if (herm) (if (packed_) "hpmv" else "hbmv") else (if (packed_) "spmv" else "sbmv");
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
        const name = comptime "cblas_" ++ portable.letter(T).? ++ @tagName(kind);
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
        const name = comptime "cblas_" ++ portable.letter(T).? ++ if (herm) (if (packed_) "hpr" else "her") else (if (packed_) "spr" else "syr");
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        if (zeroColumn(T, n, x, incx, null, 0) and nonFinite(T, n, x, incx)) return false;
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
        const name = comptime "cblas_" ++ portable.letter(T).? ++ if (herm) (if (packed_) "hpr2" else "her2") else (if (packed_) "spr2" else "syr2");
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        if (zeroColumn(T, n, x, incx, y, incy) and (nonFinite(T, n, x, incx) or nonFinite(T, n, y, incy))) return false;
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

/// `gemm`'s parts: columns of `C` if it has as many as rows, else rows;
/// each takes the matching part of `op(B)` or `op(A)` and all of the other.
fn Gemm(comptime T: type) type {
    return struct {
        const Self = @This();
        const name = "cblas_" ++ portable.letter(T).? ++ "gemm";

        b: *const Blas,
        ta: Transpose,
        tb: Transpose,
        by_columns: bool,
        m: isize,
        n: isize,
        k: isize,
        alpha: T,
        a: [*]const T,
        lda: isize,
        b_: [*]const T,
        ldb: isize,
        beta: T,
        c: [*]T,
        ldc: isize,

        fn part(g: *const Self, ch: threading.Chunk) void {
            var m = g.m;
            var n = g.n;
            var a = g.a;
            var b_ = g.b_;
            var c = g.c;
            if (g.by_columns) {
                n = units(ch);
                b_ = if (g.tb == .no_trans) column(b_, ch.start, g.ldb) else b_ + ch.start;
                c = column(c, ch.start, g.ldc);
            } else {
                m = units(ch);
                a = if (g.ta == .no_trans) a + ch.start else column(a, ch.start, g.lda);
                c += ch.start;
            }
            _ = call(name, g.b, &.{}, .{ col, g.ta, g.tb, m, n, g.k, scalar(T, &g.alpha), a, g.lda, b_, g.ldb, scalar(T, &g.beta), c, g.ldc }).?;
        }
    };
}

pub inline fn gemm(io: Io, comptime T: type, transa: types.Trans, transb: types.Trans, m: isize, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b_: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        const ta = transpose(T, transa) orelse return false;
        const tb = transpose(T, transb) orelse return false;
        if (x86_64 and (if (transa.transposed()) nonFiniteMat(T, k, m, a, lda) else nonFiniteMat(T, m, k, a, lda))) return false;
        const ea = if (transa.transposed()) mat(k, m, lda) else mat(m, k, lda);
        const eb = if (transb.transposed()) mat(n, k, ldb) else mat(k, n, ldb);
        const G = Gemm(T);
        if (!fits(G.name, b, &.{ m, n, k, lda, ldb, ldc, ea, eb, mat(m, n, ldc) })) return false;
        const g: G = .{ .b = b, .ta = ta, .tb = tb, .by_columns = n >= m, .m = m, .n = n, .k = k, .alpha = alpha, .a = a, .lda = lda, .b_ = b_, .ldb = ldb, .beta = beta, .c = c, .ldc = ldc };
        const work = threading.product(&.{ @intCast(m), @intCast(n), @intCast(k) });
        try split(io, .gemm, T, work, @intCast(@max(m, n)), matrix_quantum, &g, G.part);
        return true;
    }
}

/// `symm`'s parts: columns of `B` and `C` for `S` on the left, rows on the
/// right, each with all of `S`.
fn Symm(comptime T: type, comptime herm: bool) type {
    return struct {
        const Self = @This();
        const name = "cblas_" ++ portable.letter(T).? ++ if (herm) "hemm" else "symm";

        b: *const Blas,
        sd: Side,
        u: Uplo,
        m: isize,
        n: isize,
        alpha: T,
        a: [*]const T,
        lda: isize,
        b_: [*]const T,
        ldb: isize,
        beta: T,
        c: [*]T,
        ldc: isize,

        fn part(g: *const Self, ch: threading.Chunk) void {
            const left = g.sd == .left;
            const m = if (left) g.m else units(ch);
            const n = if (left) units(ch) else g.n;
            const b_ = if (left) column(g.b_, ch.start, g.ldb) else g.b_ + ch.start;
            const c = if (left) column(g.c, ch.start, g.ldc) else g.c + ch.start;
            _ = call(name, g.b, &.{}, .{ col, g.sd, g.u, m, n, scalar(T, &g.alpha), g.a, g.lda, b_, g.ldb, scalar(T, &g.beta), c, g.ldc }).?;
        }
    };
}

/// `symm`, or `hemm` for `herm`.
pub inline fn symm(io: Io, comptime T: type, comptime herm: bool, sd: types.Side, ul: types.Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b_: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Io.Cancelable!bool {
    if (comptime !covers(T) or (herm and !isComplex(T))) return false else {
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const ka = if (sd == .l) m else n;
        const S = Symm(T, herm);
        if (!fits(S.name, b, &.{ m, n, lda, ldb, ldc, mat(ka, ka, lda), mat(m, n, ldb), mat(m, n, ldc) })) return false;
        const g: S = .{ .b = b, .sd = side(sd), .u = u, .m = m, .n = n, .alpha = alpha, .a = a, .lda = lda, .b_ = b_, .ldb = ldb, .beta = beta, .c = c, .ldc = ldc };
        // the work `blas3/interface.zig` thresholds
        const work = threading.product(&.{ 2, @intCast(m), @intCast(m), @intCast(n) });
        try split(io, .symm, T, work, @intCast(if (sd == .l) n else m), matrix_quantum, &g, S.part);
        return true;
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
        const l = comptime portable.letter(T).?;
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

/// `trmm`'s parts: columns of `B` for `A` on the left, rows on the right,
/// as `gemm_thread_n` and `gemm_thread_m` split them.
fn Trmm(comptime T: type, comptime solve: bool) type {
    return struct {
        const Self = @This();
        const name = "cblas_" ++ portable.letter(T).? ++ if (solve) "trsm" else "trmm";

        b: *const Blas,
        sd: Side,
        u: Uplo,
        t: Transpose,
        dg: Diag,
        m: isize,
        n: isize,
        alpha: T,
        a: [*]const T,
        lda: isize,
        b_: [*]T,
        ldb: isize,

        fn part(g: *const Self, ch: threading.Chunk) void {
            const left = g.sd == .left;
            const m = if (left) g.m else units(ch);
            const n = if (left) units(ch) else g.n;
            const b_ = if (left) column(g.b_, ch.start, g.ldb) else g.b_ + ch.start;
            _ = call(name, g.b, &.{}, .{ col, g.sd, g.u, g.t, g.dg, m, n, scalar(T, &g.alpha), g.a, g.lda, b_, g.ldb }).?;
        }
    };
}

/// `trmm`, or `trsm` for `solve`.
pub inline fn trmm(io: Io, comptime T: type, comptime solve: bool, sd: types.Side, ul: types.Uplo, transa: types.Trans, dg: types.Diag, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b_: [*]T, ldb: isize) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const t = transpose(T, transa) orelse return false;
        const ka = if (sd == .l) m else n;
        if (!solve and nonFiniteMat(T, m, n, b_, ldb)) return false;
        if (solve and dg == .n and tinyDiagonal(T, ka, a, lda)) return false;
        const R = Trmm(T, solve);
        if (!fits(R.name, b, &.{ m, n, lda, ldb, mat(ka, ka, lda), mat(m, n, ldb) })) return false;
        const g: R = .{ .b = b, .sd = side(sd), .u = u, .t = t, .dg = diag(dg), .m = m, .n = n, .alpha = alpha, .a = a, .lda = lda, .b_ = b_, .ldb = ldb };
        const work = threading.product(&.{ @intCast(m), @intCast(n) });
        try split(io, .trsm, T, work, @intCast(if (sd == .l) n else m), matrix_quantum, &g, R.part);
        return true;
    }
}

// LAPACK

/// `?potrf_`'s `INFO`, or null where rstd runs it. Accelerate's LAPACK runs
/// on its BLAS, single-threaded as `load` sets it.
pub inline fn potrf(comptime T: type, ul: types.Uplo, n: isize, a: [*]T, lda: isize) ?isize {
    if (comptime !covers(T) or (x86_64 and isComplex(T))) return null else {
        const b = use(T) orelse return null;
        const u: u8 = switch (ul) {
            .u => 'U',
            .l => 'L',
            .v, .m => return null,
        };
        switch (b.*) {
            inline else => |*bound, tag| {
                const f = @field(bound.fns, portable.letter(T).? ++ "potrf_") orelse return null;
                const I = if (tag == .ilp64) i64 else c_int;
                if (I != i64 and mat(n, n, lda) > std.math.maxInt(I)) return null;
                const nn = std.math.cast(I, n) orelse return null;
                const ld = std.math.cast(I, lda) orelse return null;
                var info: I = 0;
                if (sysinfo.is_test) _ = calls.fetchAdd(1, .monotonic);
                f(&u, &nn, @ptrCast(a), &ld, &info);
                return info;
            },
        }
    }
}
