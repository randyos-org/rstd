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

/// Where `load` can find Accelerate.
pub const available = sysinfo.os.tag == .macos;

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

pub const Accelerate = struct {
    blas: Blas,
    threads: ThreadApi,

    pub const ThreadError = error{ Unsupported, Rejected };

    /// The calling thread's model; null before macOS 15.
    pub fn threading(self: *const Accelerate) ?Threading {
        const get = self.threads.BLASGetThreading orelse return null;
        return get();
    }

    /// Sets the calling thread's model (`BLASSetThreading`); other threads,
    /// and threads started later, keep the default (multi-threaded).
    pub fn setThreading(self: *const Accelerate, model: Threading) ThreadError!void {
        const set = self.threads.BLASSetThreading orelse return error.Unsupported;
        if (set(model) != 0) return error.Rejected;
    }

    /// Makes the calling thread's calls single-threaded where macOS can
    /// (15+); a no-op before that, where only `VECLIB_MAXIMUM_THREADS` does.
    pub fn singleThreadHere(self: *const Accelerate) void {
        self.setThreading(.single_threaded) catch {};
    }
};

pub const Options = struct {
    /// Run Accelerate on one thread: its threaded `gemv` rounds differently
    /// from its single-threaded one, so only one thread gives the same bits
    /// everywhere. Sets `VECLIB_MAXIMUM_THREADS=1` (`setMaxThreadsEnv`) and
    /// the loading thread's model; other threads need `singleThreadHere`.
    single_threaded: bool = true,
};

/// Opens Accelerate and binds its CBLAS; null off macOS or without it. The
/// framework handle is kept open for the process: the table points into it.
pub fn load(options: Options) ?Accelerate {
    if (comptime !available) return null;
    var lib = dll.openFramework("Accelerate") catch return null;
    const ilp = dll.bind(Table(i64), &lib, ilp64Names);
    const blas: Blas = if (ilp.count() > 0) .{ .ilp64 = ilp } else blk: {
        const lp = dll.bind(Table(c_int), &lib, lp64Names);
        if (lp.count() == 0) return null;
        break :blk .{ .lp64 = lp };
    };
    const acc: Accelerate = .{ .blas = blas, .threads = dll.bind(ThreadApi, &lib, dll.exactName).fns };
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

// ---- tests -------------------------------------------------------------------

const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1/interface.zig");
const blas3 = @import("blas3/interface.zig");
const Complex64 = types.Complex64;

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
    var lib = try dll.openFramework("Accelerate");
    const lp = dll.bind(Table(c_int), &lib, lp64Names);
    try zigstd.testExpectEqualStrings("cblas_ddot$NEWLAPACK", lp.name(lp64Names, .cblas_ddot).?);
    var x: [100]f64 = undefined;
    for (&x, 0..) |*v, i| v.* = @floatFromInt(i % 7);
    try zigstd.testExpectEqual(acc.blas.ilp64.fns.cblas_ddot.?(100, &x, 1, &x, 1), lp.fns.cblas_ddot.?(100, &x, 1, &x, 1));
}

test "ddot, dgemm, zdotu_sub and isamax agree with rstd" {
    const acc = load(.{}) orelse return error.SkipZigTest;
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
            const want = try blas1.dot(f64, 300, &x, 2, &y, -1);
            try zigstd.testExpectApproxEqRel(want, got, 1e-12);

            // zdotu_sub and zdotc_sub
            var zx: [257]Complex64 = undefined;
            var zy: [257]Complex64 = undefined;
            fill(f64, rng, @as(*[514]f64, @ptrCast(&zx)));
            fill(f64, rng, @as(*[514]f64, @ptrCast(&zy)));
            var zd: Complex64 = undefined;
            t.cblas_zdotu_sub.?(zx.len, &zx, 1, &zy, 1, &zd);
            const zw = try blas1.dot(Complex64, zx.len, &zx, 1, &zy, 1);
            try zigstd.testExpectApproxEqAbs(zw.re, zd.re, 1e-12);
            try zigstd.testExpectApproxEqAbs(zw.im, zd.im, 1e-12);
            t.cblas_zdotc_sub.?(zx.len, &zx, 1, &zy, 1, &zd);
            const zc = try blas1.dotc(Complex64, zx.len, &zx, 1, &zy, 1);
            try zigstd.testExpectApproxEqAbs(zc.re, zd.re, 1e-12);
            try zigstd.testExpectApproxEqAbs(zc.im, zd.im, 1e-12);

            // isamax: 0-based where rstd's is 1-based
            var sx: [1000]f32 = undefined;
            fill(f32, rng, &sx);
            sx[613] = -3;
            const ia: I = t.cblas_isamax.?(sx.len, &sx, 1);
            try zigstd.testExpectEqual(@as(I, 613), ia);
            try zigstd.testExpectEqual(try blas1.index(f32, .amax, sx.len, &sx, 1), ia + 1);

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
            try blas3.gemm(gpa, f64, .t, .n, m, n, k, 1.5, a.ptr, k, bm.ptr, k, -0.5, c1.ptr, m);
            for (c0, c1) |g, w| try zigstd.testExpectApproxEqAbs(w, g, 1e-12 * k);
            // row-major: C^T = B^T A, the same numbers transposed
            t.cblas_dgemm.?(.row_major, .no_trans, .trans, n, m, k, 1, bm.ptr, k, a.ptr, k, 0, c2.ptr, m);
            try blas3.gemm(gpa, f64, .t, .n, m, n, k, 1, a.ptr, k, bm.ptr, k, 0, c1.ptr, m);
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

fn digestOnThread(comptime I: type, acc: *const Accelerate, t: *const Table(I), d: *const ReproData, out: *Digest) void {
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
