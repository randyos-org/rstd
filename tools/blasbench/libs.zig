//! The library races: `lib1`, `lib2`, `lib3`, and the loading of each library.

const std = @import("std");
const sysinfo = @import("builtin");
const Io = std.Io;
const rstd = @import("rstd");
const bench = @import("bench.zig");
const merge = @import("merge.zig");
const bb = @import("blas_bench.zig");
const kernels = @import("kernels.zig");
const ext = @import("ext.zig");

const blas = rstd.linalg.blas;
const lapack = rstd.linalg.lapack;
const Dll = rstd.dll.Dll;
const os_tag = sysinfo.target.os.tag;

const arm = bb.arm;
const clock = bb.clock;
const fillRandom = bb.fillRandom;
const Fn = bb.Fn;
const l1_timing = bb.l1_timing;
const l3_timing = bb.l3_timing;
const ob_asm = bb.ob_asm;
const Prec = bb.Prec;
const precs = bb.precs;

const asmFor = kernels.asmFor;
const asmGemm = kernels.asmGemm;
const Blocking = kernels.Blocking;
const Bufs = kernels.Bufs;
const candidate = kernels.candidate;
const checkGemm = kernels.checkGemm;
const checkGemv = kernels.checkGemv;
const checkL1 = kernels.checkL1;
const gemm_sizes = kernels.gemm_sizes;
const GemmCall = kernels.GemmCall;
const gemv_alpha = kernels.gemv_alpha;
const GemvCall = kernels.GemvCall;
const gemvCandidate = kernels.gemvCandidate;
const l1_cases = kernels.l1_cases;
const l1Bytes = kernels.l1Bytes;
const L1 = kernels.L1;
const L1Call = kernels.L1Call;
const l1Elem = kernels.l1Elem;
const l1Result = kernels.l1Result;
const l1Scalar = kernels.l1Scalar;
const l1Slot = kernels.l1Slot;
const l2_bytes = kernels.l2_bytes;
const reduces = kernels.reduces;
const rounded = kernels.rounded;

const ext_cols = ext.ext_cols;
const ext_l2 = ext.ext_l2;
const ext_l3 = ext.ext_l3;
const extTable = ext.extTable;

// ---- libraries -------------------------------------------------------------

pub const LibOptions = struct {
    /// Each library's path; `rstd`'s is unused.
    paths: std.EnumArray(Src, ?[]const u8) = .init(.{
        .rstd = null,
        .openblas = null,
        .blis = null,
        .accelerate = if (os_tag.isDarwin()) "/System/Library/Frameworks/Accelerate.framework/Accelerate" else null,
        .reference = null,
        .eigen = null,
        .aoclflame = null,
    }),
    /// Each library's thread count and rstd's ceiling; null for their own.
    threads: ?u31 = 1,
    rstd_io: RstdIo = .threaded,
    ob_core: ?[]const u8 = null,
    /// `--rstd-accelerate=0`: rstd's routines run its own kernels, not Accelerate.
    rstd_accelerate: bool = true,
};

/// Where a library contestant comes from: rstd's public routines, or a
/// library's `cblas_*` or Fortran-ABI routines. `reference` is
/// Reference-LAPACK's `liblapack` (its BLAS through its own dependency),
/// `eigen` the bench's shim (`nix/eigen_shim.cpp`), `aoclflame`
/// AOCL-libFLAME's LAPACK.
pub const Src = enum { rstd, openblas, blis, accelerate, reference, eigen, aoclflame };

/// Whether `src` is a library every table races (OpenBLAS, BLIS,
/// Accelerate): one that cannot load stops the bench, and only these load
/// in `alone`. The reference binaries are skipped when they cannot load.
fn raced(src: Src) bool {
    return switch (src) {
        .openblas, .blis, .accelerate => true,
        .rstd, .reference, .eigen, .aoclflame => false,
    };
}

/// Each library's handle; null when not given, and for `rstd`.
var lib_handles: std.EnumArray(Src, ?Dll) = .initFill(null);

pub fn loaded(src: Src) bool {
    if (src == .rstd) return !bb.alone;
    return lib_handles.get(src) != null;
}

/// The first error an rstd routine returned inside a kernel-ABI wrapper.
pub var lib_err: ?anyerror = null;

pub fn fail(e: anyerror) void {
    if (lib_err == null) lib_err = e;
}

/// rstd's routines take an allocator; a C program would hand them malloc's.
pub const rstd_alloc = std.heap.c_allocator;
/// The `Io` rstd's routines thread on.
pub const RstdIo = enum {
    /// The bench's `Io.Threaded`, limited to `--threads` less the caller.
    threaded,
    /// Serial: `Io.failing` runs every task inline.
    failing,
};

/// And an `Io`: `main` sets it.
pub var rstd_io: Io = .failing;

/// Which `Io` rstd ran on, for each library table's notes.
pub var io_note: []const u8 = "";
var io_note_buf: [96]u8 = undefined;

/// `--max-n`: rows of a larger `n` are skipped.
pub var n_max: usize = std.math.maxInt(usize);

/// The libraries' thread counts by environment, read when each loads, and
/// rstd's ceiling (`OPENBLAS_NUM_THREADS`), read when it selects its core.
pub fn setThreadEnv(threads: ?u31) !void {
    const t = threads orelse return;
    var buf: [16]u8 = undefined;
    const v = try std.mem.printSentinel(&buf, "{d}", .{t}, 0);
    for ([_][*:0]const u8{ "OPENBLAS_NUM_THREADS", "BLIS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS" }) |name| setEnv(name, v);
}

extern "c" fn setenv(name: [*:0]const u8, value: [*:0]const u8, overwrite: c_int) c_int;
pub extern "c" fn unsetenv(name: [*:0]const u8) c_int;
extern "kernel32" fn SetEnvironmentVariableA(name: [*:0]const u8, value: [*:0]const u8) callconv(.winapi) c_int;

fn setEnv(name: [*:0]const u8, value: [*:0]const u8) void {
    if (os_tag == .windows) _ = SetEnvironmentVariableA(name, value) else _ = setenv(name, value, 1);
}

/// The loader's reason for a failed open, where it gives one.
fn dlError() []const u8 {
    if (os_tag == .windows) return "";
    return if (std.c.dlerror()) |m| std.mem.span(m) else "";
}

/// `name` from `src`'s library: Accelerate's current LP64 entry
/// (`ACCELERATE_NEW_LAPACK`) where it has one, else the classic one.
fn lookup(comptime src: Src, comptime T: type, comptime name: [:0]const u8) ?T {
    const h = &(lib_handles.getPtr(src).* orelse return null);
    if (src == .accelerate) {
        if (h.lookup(T, name ++ "$NEWLAPACK")) |f| return f;
    }
    return h.lookup(T, name);
}

/// Fortran-ABI routine `base` (`dpotrf`) from `src`'s library: Accelerate's
/// new LP64 interface first (`dpotrf$NEWLAPACK`, or with the underscore:
/// `apple_accel.lp64Names`), else `base_`.
pub fn lookupF(comptime src: Src, comptime T: type, comptime base: []const u8) ?T {
    const h = &(lib_handles.getPtr(src).* orelse return null);
    if (src == .accelerate) {
        if (h.lookup(T, std.fmt.comptimePrint("{s}$NEWLAPACK", .{base}))) |f| return f;
        if (h.lookup(T, std.fmt.comptimePrint("{s}_$NEWLAPACK", .{base}))) |f| return f;
    }
    return h.lookup(T, std.fmt.comptimePrint("{s}_", .{base}));
}

/// Opens the libraries given, at `--threads` each (`setThreadEnv`), and
/// prints how each is configured, and rstd. OpenBLAS is required unless
/// `alone`; the reference binaries are skipped when they cannot load.
pub fn loadLibs(opts: LibOptions) !void {
    var buf: [64]u8 = undefined;
    if (opts.ob_core) |c| setEnv("OPENBLAS_CORETYPE", try std.mem.printSentinel(&buf, "{s}", .{c}, 0));
    const paths = opts.paths;
    if (paths.get(.openblas) == null and !bb.alone) return error.NeedOpenblasPath;
    var it = lib_handles.iterator();
    while (it.next()) |h| {
        const p = paths.get(h.key) orelse continue;
        // in `alone`, only what Accelerate is raced against
        if (bb.alone and !raced(h.key)) continue;
        h.value.* = Dll.open(p) catch |e| {
            std.debug.print("cannot load {s}: {t} {s}\n", .{ p, e, dlError() });
            if (raced(h.key)) return e;
            continue;
        };
    }

    const Str = *const fn () callconv(.c) [*:0]const u8;
    const Int = *const fn () callconv(.c) c_int;
    if (loaded(.openblas)) {
        const config = std.mem.span(lookup(.openblas, Str, "openblas_get_config").?());
        // `int` indices here; an ILP64 build (`USE64BITINT`) would need `long`
        if (std.mem.find(u8, config, "USE64BITINT") != null) return error.OpenblasIsIlp64;
        // Windows' CRT may have copied the environment already: by call too
        if (opts.threads) |t| lookup(.openblas, *const fn (c_int) callconv(.c) void, "openblas_set_num_threads").?(t);
        std.debug.print("openblas: {s}\nopenblas core: {s}, config: {s}, threads: {d}, lapack: {s}\n", .{
            paths.get(.openblas).?,
            lookup(.openblas, Str, "openblas_get_corename").?(),
            config,
            lookup(.openblas, Int, "openblas_get_num_threads").?(),
            lapackVersion(.openblas, &buf),
        });
    }
    if (loaded(.reference)) std.debug.print("reference: {s}, lapack {s}, one thread\n", .{ paths.get(.reference).?, lapackVersion(.reference, &buf) });
    if (loaded(.eigen)) {
        const v = if (lookup(.eigen, Str, "rstd_bench_eigen_version")) |f| std.mem.span(f()) else "?";
        std.debug.print("eigen: {s}, version {s}, one thread\n", .{ paths.get(.eigen).?, v });
    }
    if (loaded(.aoclflame)) std.debug.print("aocl-libflame: {s}, lapack {s}, threads by OMP_NUM_THREADS\n", .{ paths.get(.aoclflame).?, lapackVersion(.aoclflame, &buf) });

    if (loaded(.blis)) {
        // `dim_t` is 64 bits, `arch_t` an enum
        const Dim = *const fn () callconv(.c) i64;
        if (opts.threads) |t| lookup(.blis, *const fn (i64) callconv(.c) void, "bli_thread_set_num_threads").?(t);
        const arch = lookup(.blis, *const fn () callconv(.c) c_int, "bli_arch_query_id").?();
        std.debug.print("blis: {s}\nblis version: {s}, config: {s}, threads: {d}, cblas: {s}\n", .{
            paths.get(.blis).?,
            lookup(.blis, Str, "bli_info_get_version_str").?(),
            lookup(.blis, *const fn (c_int) callconv(.c) [*:0]const u8, "bli_arch_string").?(arch),
            lookup(.blis, Dim, "bli_thread_get_num_threads").?(),
            if (lookup(.blis, *const anyopaque, "cblas_dgemm") != null) "yes" else "no",
        });
    }

    if (loaded(.accelerate)) {
        // per thread: this one makes every call; more than one thread is
        // `VECLIB_MAXIMUM_THREADS`'s
        if (opts.threads == 1) {
            if (lookup(.accelerate, *const fn (c_uint) callconv(.c) c_int, "BLASSetThreading")) |f| {
                if (f(1) != 0) std.debug.print("warning: BLASSetThreading failed\n", .{});
            }
        }
        std.debug.print("accelerate: {s}, BLASGetThreading: {s}\n", .{ paths.get(.accelerate).?, accelThreading() });
    }
    if (bb.alone) {
        io_note = "rstd is not raced here.";
        return;
    }
    const routed = blas.useAccelerate(opts.rstd_accelerate);
    const io_name = switch (opts.rstd_io) {
        .threaded => "Io.Threaded",
        .failing => "Io.failing",
    };
    const ceiling = blas.config.threadCeiling();
    std.debug.print("rstd core: {s}, accelerate: {}, io: {s}, ceiling: {?d}, cpus: {d}\n", .{ blas.coreName(), routed, io_name, ceiling, blas.config.cpu_count });
    io_note = if (ceiling) |c|
        try std.fmt.bufPrint(&io_note_buf, "rstd runs on `{s}`, thread ceiling {d}.", .{ io_name, c })
    else
        try std.fmt.bufPrint(&io_note_buf, "rstd runs on `{s}`, threads up to the CPUs.", .{io_name});
    // rstd's own load of Accelerate makes this thread single-threaded
    if (loaded(.accelerate)) std.debug.print("accelerate after rstd's init: BLASGetThreading: {s}\n", .{accelThreading()});
}

/// `src`'s LAPACK version (`ilaver`), or `-` without one.
fn lapackVersion(comptime src: Src, buf: []u8) []const u8 {
    const f = lookupF(src, *const fn (*c_int, *c_int, *c_int) callconv(.c) void, "ilaver") orelse return "-";
    var v: [3]c_int = undefined;
    f(&v[0], &v[1], &v[2]);
    return std.fmt.bufPrint(buf, "{d}.{d}.{d}", .{ v[0], v[1], v[2] }) catch "?";
}

/// The calling thread's `BLASGetThreading`, as words.
fn accelThreading() []const u8 {
    const f = lookup(.accelerate, *const fn () callconv(.c) c_uint, "BLASGetThreading") orelse return "absent";
    return switch (f()) {
        0 => "multi-threaded",
        1 => "single-threaded",
        else => "?",
    };
}

/// `p`'s element as rstd's routines take it.
pub fn Num(comptime p: Prec) type {
    return if (p.complex) std.math.Complex(p.F) else p.F;
}

pub fn cast(comptime p: Prec, x: [*]p.F) [*]Num(p) {
    return @ptrCast(x);
}

/// cblas's `CBLAS_ORDER` and `CBLAS_TRANSPOSE`.
const col_major: c_int = 102;
const no_trans: c_int = 111;
const trans_t: c_int = 112;

/// Level-1 routine `r` over `p` from `src`, in the kernel ABI of its vtable
/// slot, so `checkL1` and `L1Call` take it as they take a kernel: `f` the
/// library's symbol, `k` the wrapper. `from` makes each library's struct
/// a type of its own; without it they would share one `f`.
fn L1Via(comptime src: Src, comptime r: L1, comptime p: Prec) type {
    const F = p.F;
    const T = Num(p);
    const Cx = blas.types.CComplex(F);
    const lib = src != .rstd;
    const s = p.name;
    // the real letter, for `scnrm2`, `csrot` and the like
    const re = if (F == f32) "s" else "d";
    const I = isize;
    const cx = p.complex;
    return switch (r) {
        .axpy => if (!cx) struct {
            const sym = "cblas_" ++ s ++ "axpy";
            const from = src;
            const Sig = fn (c_int, F, [*]F, c_int, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, _: I, _: I, a: F, x: [*]F, ix: I, y: [*]F, iy: I, _: [*]F, _: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), a, x, @intCast(ix), y, @intCast(iy)) else blas.axpy(rstd_io, T, n, a, x, ix, y, iy) catch |e| fail(e);
                return 0;
            }
        } else struct {
            const sym = "cblas_" ++ s ++ "axpy";
            const from = src;
            const Sig = fn (c_int, *const [2]F, [*]F, c_int, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, _: I, _: I, ar: F, ai: F, x: [*]F, ix: I, y: [*]F, iy: I, _: [*]F, _: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), &.{ ar, ai }, x, @intCast(ix), y, @intCast(iy)) else blas.axpy(rstd_io, T, n, .init(ar, ai), cast(p, x), ix, cast(p, y), iy) catch |e| fail(e);
                return 0;
            }
        },
        .dot => if (!cx) struct {
            const sym = "cblas_" ++ s ++ "dot";
            const from = src;
            const Sig = fn (c_int, [*]F, c_int, [*]F, c_int) callconv(.c) F;
            var f: *const Sig = undefined;
            fn k(n: I, x: [*]F, ix: I, y: [*]F, iy: I) callconv(.c) F {
                if (lib) return f(@intCast(n), x, @intCast(ix), y, @intCast(iy));
                return blas.dot(rstd_io, T, n, x, ix, y, iy) catch |e| {
                    fail(e);
                    return 0;
                };
            }
        } else struct {
            const sym = "cblas_" ++ s ++ "dotu_sub";
            const from = src;
            const Sig = fn (c_int, [*]F, c_int, [*]F, c_int, *[2]F) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, x: [*]F, ix: I, y: [*]F, iy: I) callconv(.c) Cx {
                if (lib) {
                    var o: [2]F = undefined;
                    f(@intCast(n), x, @intCast(ix), y, @intCast(iy), &o);
                    return .{ .re = o[0], .im = o[1] };
                }
                const v = blas.dot(rstd_io, T, n, cast(p, x), ix, cast(p, y), iy) catch |e| {
                    fail(e);
                    return .{ .re = 0, .im = 0 };
                };
                return .{ .re = v.re, .im = v.im };
            }
        },
        .scal => if (!cx) struct {
            const sym = "cblas_" ++ s ++ "scal";
            const from = src;
            const Sig = fn (c_int, F, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, _: I, _: I, a: F, x: [*]F, ix: I, _: [*]F, _: I, _: [*]F, _: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), a, x, @intCast(ix)) else blas.scal(rstd_io, T, n, a, x, ix) catch |e| fail(e);
                return 0;
            }
        } else struct {
            const sym = "cblas_" ++ s ++ "scal";
            const from = src;
            const Sig = fn (c_int, *const [2]F, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, _: I, _: I, ar: F, ai: F, x: [*]F, ix: I, _: [*]F, _: I, _: [*]F, _: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), &.{ ar, ai }, x, @intCast(ix)) else blas.scal(rstd_io, T, n, .init(ar, ai), cast(p, x), ix) catch |e| fail(e);
                return 0;
            }
        },
        .nrm2, .asum => struct {
            const sym = "cblas_" ++ (if (cx) re else "") ++ s ++ @tagName(r);
            const from = src;
            const Sig = fn (c_int, [*]F, c_int) callconv(.c) F;
            var f: *const Sig = undefined;
            fn k(n: I, x: [*]F, ix: I) callconv(.c) F {
                if (lib) return f(@intCast(n), x, @intCast(ix));
                return @field(blas, @tagName(r))(rstd_io, T, n, cast(p, x), ix) catch |e| {
                    fail(e);
                    return 0;
                };
            }
        },
        .iamax => struct {
            const sym = "cblas_i" ++ s ++ "amax";
            const from = src;
            // OpenBLAS returns `size_t`, Accelerate `int`: both from 0
            const Sig = fn (c_int, [*]F, c_int) callconv(.c) c_int;
            var f: *const Sig = undefined;
            fn k(n: I, x: [*]F, ix: I) callconv(.c) I {
                if (lib) return @as(I, f(@intCast(n), x, @intCast(ix))) + 1;
                return blas.index(rstd_io, T, .amax, n, cast(p, x), ix) catch |e| {
                    fail(e);
                    return 0;
                };
            }
        },
        .copy => struct {
            const sym = "cblas_" ++ s ++ "copy";
            const from = src;
            const Sig = fn (c_int, [*]F, c_int, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, x: [*]F, ix: I, y: [*]F, iy: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), x, @intCast(ix), y, @intCast(iy)) else blas.copy(rstd_io, T, n, cast(p, x), ix, cast(p, y), iy) catch |e| fail(e);
                return 0;
            }
        },
        .swap => if (!cx) struct {
            const sym = "cblas_" ++ s ++ "swap";
            const from = src;
            const Sig = fn (c_int, [*]F, c_int, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, _: I, _: I, _: F, x: [*]F, ix: I, y: [*]F, iy: I, _: [*]F, _: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), x, @intCast(ix), y, @intCast(iy)) else blas.swap(rstd_io, T, n, x, ix, y, iy) catch |e| fail(e);
                return 0;
            }
        } else struct {
            const sym = "cblas_" ++ s ++ "swap";
            const from = src;
            const Sig = fn (c_int, [*]F, c_int, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, _: I, _: I, _: F, _: F, x: [*]F, ix: I, y: [*]F, iy: I, _: [*]F, _: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), x, @intCast(ix), y, @intCast(iy)) else blas.swap(rstd_io, T, n, cast(p, x), ix, cast(p, y), iy) catch |e| fail(e);
                return 0;
            }
        },
        .rot => struct {
            const sym = "cblas_" ++ s ++ (if (cx) re else "") ++ "rot";
            const from = src;
            const Sig = fn (c_int, [*]F, c_int, [*]F, c_int, F, F) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, x: [*]F, ix: I, y: [*]F, iy: I, c: F, sn: F) callconv(.c) i32 {
                if (lib) f(@intCast(n), x, @intCast(ix), y, @intCast(iy), c, sn) else blas.rot(rstd_io, T, n, cast(p, x), ix, cast(p, y), iy, c, sn) catch |e| fail(e);
                return 0;
            }
        },
    };
}

/// `src`'s routine `r` as a kernel; null where the library lacks it.
fn libL1(comptime src: Src, comptime r: L1, comptime p: Prec) ?Fn(l1Slot(r, p)) {
    const W = L1Via(src, r, p);
    if (src != .rstd) W.f = lookup(src, *const W.Sig, W.sym) orelse return null;
    return &W.k;
}

/// `gemv` over `p` from `src` in the kernel ABI of `?gemv_n`/`?gemv_t`:
/// `y += alpha*op(A)*x` (`beta` 1).
fn GemvVia(comptime src: Src, comptime p: Prec, comptime t: []const u8) type {
    const F = p.F;
    const T = Num(p);
    const lib = src != .rstd;
    const I = isize;
    const tr = if (t[0] == 't') trans_t else no_trans;
    const rt: blas.Trans = if (t[0] == 't') .t else .n;
    return if (!p.complex) struct {
        const sym = "cblas_" ++ p.name ++ "gemv";
        const from = src;
        const Sig = fn (c_int, c_int, c_int, c_int, F, [*]F, c_int, [*]F, c_int, F, [*]F, c_int) callconv(.c) void;
        var f: *const Sig = undefined;
        fn k(m: I, n: I, _: I, a: F, ma: [*]F, lda: I, x: [*]F, ix: I, y: [*]F, iy: I, _: [*]F) callconv(.c) i32 {
            if (lib) f(col_major, tr, @intCast(m), @intCast(n), a, ma, @intCast(lda), x, @intCast(ix), 1, y, @intCast(iy)) else blas.gemv(rstd_alloc, rstd_io, T, rt, m, n, a, ma, lda, x, ix, 1, y, iy) catch |e| fail(e);
            return 0;
        }
    } else struct {
        const sym = "cblas_" ++ p.name ++ "gemv";
        const from = src;
        const Sig = fn (c_int, c_int, c_int, c_int, *const [2]F, [*]F, c_int, [*]F, c_int, *const [2]F, [*]F, c_int) callconv(.c) void;
        var f: *const Sig = undefined;
        fn k(m: I, n: I, _: I, ar: F, ai: F, ma: [*]F, lda: I, x: [*]F, ix: I, y: [*]F, iy: I, _: [*]F) callconv(.c) i32 {
            if (lib) f(col_major, tr, @intCast(m), @intCast(n), &.{ ar, ai }, ma, @intCast(lda), x, @intCast(ix), &.{ 1, 0 }, y, @intCast(iy)) else blas.gemv(rstd_alloc, rstd_io, T, rt, m, n, .init(ar, ai), cast(p, ma), lda, cast(p, x), ix, .init(1, 0), cast(p, y), iy) catch |e| fail(e);
            return 0;
        }
    };
}

fn libGemv(comptime src: Src, comptime p: Prec, comptime t: []const u8) ?Fn(p.name ++ "gemv_" ++ t) {
    const W = GemvVia(src, p, t);
    if (src != .rstd) W.f = lookup(src, *const W.Sig, W.sym) orelse return null;
    return &W.k;
}

/// Square `NN` gemm, `C += A*B`, column-major with leading dimension `n`.
fn GemmFn(comptime p: Prec) type {
    return *const fn (n: usize, a: [*]p.F, b: [*]p.F, c: [*]p.F) void;
}

fn GemmVia(comptime src: Src, comptime p: Prec) type {
    const F = p.F;
    const lib = src != .rstd;
    return if (!p.complex) struct {
        const sym = "cblas_" ++ p.name ++ "gemm";
        const from = src;
        const Sig = fn (c_int, c_int, c_int, c_int, c_int, c_int, F, [*]F, c_int, [*]F, c_int, F, [*]F, c_int) callconv(.c) void;
        var f: *const Sig = undefined;
        fn k(n: usize, a: [*]F, b: [*]F, c: [*]F) void {
            const ni: c_int = @intCast(n);
            if (lib) f(col_major, no_trans, no_trans, ni, ni, ni, 1, a, ni, b, ni, 1, c, ni) else blas.gemm(rstd_alloc, rstd_io, F, .n, .n, ni, ni, ni, 1, a, ni, b, ni, 1, c, ni) catch |e| fail(e);
        }
    } else struct {
        const sym = "cblas_" ++ p.name ++ "gemm";
        const from = src;
        const Sig = fn (c_int, c_int, c_int, c_int, c_int, c_int, *const [2]F, [*]F, c_int, [*]F, c_int, *const [2]F, [*]F, c_int) callconv(.c) void;
        var f: *const Sig = undefined;
        fn k(n: usize, a: [*]F, b: [*]F, c: [*]F) void {
            const ni: c_int = @intCast(n);
            const one: Num(p) = .init(1, 0);
            if (lib) f(col_major, no_trans, no_trans, ni, ni, ni, &.{ 1, 0 }, a, ni, b, ni, &.{ 1, 0 }, c, ni) else blas.gemm(rstd_alloc, rstd_io, Num(p), .n, .n, ni, ni, ni, one, cast(p, a), ni, cast(p, b), ni, one, cast(p, c), ni) catch |e| fail(e);
        }
    };
}

/// `src`'s gemm; null where the library has none.
fn libGemm(comptime src: Src, comptime p: Prec) ?GemmFn(p) {
    const W = GemmVia(src, p);
    if (src != .rstd) W.f = lookup(src, *const W.Sig, W.sym) orelse return null;
    return &W.k;
}

/// A library table's column: a contestant, and its library if it is one.
pub const Col = struct {
    label: []const u8,
    short: []const u8,
    lib: ?Src = null,
    /// Whether this build has it, for a column that is no library.
    on: bool = true,
    /// Its rows of `n` up to this only (a slow baseline).
    max_n: usize = std.math.maxInt(usize),

    /// rstd's own columns (rstd, `asm`, `vec x4`, `ob .S`) are out in `alone`.
    pub fn shown(c: Col) bool {
        return if (c.lib) |l| loaded(l) else c.on and !bb.alone;
    }
};

/// `lib1` and `lib2`'s contestants. Ratios are over the first, OpenBLAS.
pub const lib_cols = [_]Col{
    .{ .label = "openblas", .short = "ob", .lib = .openblas },
    .{ .label = "rstd", .short = "rstd" },
    .{ .label = "blis", .short = "blis", .lib = .blis },
    .{ .label = "accelerate", .short = "accel", .lib = .accelerate },
    .{ .label = "asm", .short = "asm", .on = arm },
    .{ .label = "vec x4", .short = "x4", .on = arm },
};
/// `lib3`'s: `ob .S` is OpenBLAS's vortex `.S` through the bench's loops.
const lib3_cols = lib_cols[0..4].* ++ [_]Col{.{ .label = "ob .S", .short = "ob .S", .on = ob_asm }};
pub const rstd_col = 1;

/// Column `o` of a race in the order checked: rstd's first, whose output
/// the libraries' are compared with.
pub fn rstdFirst(o: usize) usize {
    return switch (o) {
        0 => rstd_col,
        rstd_col => 0,
        else => o,
    };
}

/// rstd's largest difference from a library's result, in eps of the type
/// relative to the size of the terms, and what accumulation order explains.
pub const Diff = struct { eps: f64, bound: f64 };

/// What accumulation order explains for a sum of `n` terms: a few ulps
/// times `n`.
pub fn sumBound(n: usize) f64 {
    return 2 * (@as(f64, @floatFromInt(n)) + 8);
}

/// A table's header: `lead` and `rule`, then each shown contestant's rate,
/// its ratio over OpenBLAS, rstd's difference from each library, rstd's
/// error, and a hash of rstd's output.
pub fn printHeader(comptime lead: []const u8, comptime rule: []const u8, cols: []const Col) void {
    var n: usize = 0;
    std.debug.print("| " ++ lead ++ " |", .{});
    for (cols) |c| if (c.shown()) {
        std.debug.print(" {s} |", .{c.label});
        n += 1;
    };
    for (cols[1..]) |c| if (c.shown()) {
        std.debug.print(" {s}/ob |", .{c.short});
        n += 1;
    };
    if (bb.alone) {
        std.debug.print("\n| " ++ rule ++ " |", .{});
        for (0..n) |_| std.debug.print(" -: |", .{});
        std.debug.print("\n", .{});
        return;
    }
    for (cols) |c| if (c.lib != null and c.shown()) {
        std.debug.print(" rstd-{s} |", .{c.short});
        n += 1;
    };
    std.debug.print(" rstd err | rstd bits |\n| " ++ rule ++ " |", .{});
    for (0..n + 1) |_| std.debug.print(" -: |", .{});
    std.debug.print(" --- |\n", .{});
}

/// The column labels a merge (`merge.zig`) reads as rates, and each one's
/// short name in a ratio's header (`x/ob`).
pub const merge_labels = blk: {
    var l: []const merge.Label = &.{};
    for (lib_cols ++ ext_cols[4..] ++ lib3_cols[4..]) |c| l = l ++ [_]merge.Label{.{ .label = c.label, .short = c.short }};
    break :blk l[0..l.len].*;
};

/// Marks the start of a section's output for `merge.zig`, which reruns
/// noisy sections by name.
pub fn sectionMark(comptime name: []const u8) void {
    std.debug.print("\n<!-- section " ++ name ++ " -->\n", .{});
}

/// A row after its lead: `printHeader`'s columns; `-` where one is missing,
/// `!` on a difference past its bound.
pub fn printRow(cols: []const Col, rates: []const ?f64, diffs: []const ?Diff, mine: Rstd) void {
    for (cols, rates) |c, r| {
        if (!c.shown()) continue;
        if (r) |v| std.debug.print("{d:.1} | ", .{v}) else std.debug.print("- | ", .{});
    }
    for (cols[1..], rates[1..]) |c, r| {
        if (!c.shown()) continue;
        if (r != null and rates[0] != null) std.debug.print("{d:.2} | ", .{r.? / rates[0].?}) else std.debug.print("- | ", .{});
    }
    if (bb.alone) return std.debug.print("\n", .{});
    for (cols, diffs) |c, d| {
        if (c.lib == null or !c.shown()) continue;
        if (d) |v| std.debug.print("{d:.1}{s} | ", .{ v.eps, if (v.eps > v.bound) " !" else "" }) else std.debug.print("- | ", .{});
    }
    if (mine.err) |e| std.debug.print("{d:.1} | ", .{e}) else std.debug.print("- | ", .{});
    if (mine.bits) |h| std.debug.print("{x:0>8} |\n", .{@as(u32, @truncate(h))}) else std.debug.print("- |\n", .{});
}

/// rstd's error from the f64 reference, in eps of the type relative to the
/// size of the terms, and a hash of its output.
pub const Rstd = struct { err: ?f64 = null, bits: ?u64 = null };

/// A library table row's contestants, `Ctx`s entered in checking order
/// (`rstdFirst`) for `n` columns, and what `printRow` shows beside their rates.
pub fn Entrants(comptime Ctx: type, comptime n: usize) type {
    return struct {
        ctxs: [n]Ctx = undefined,
        /// Each column's place among `ctxs`; null for one not entered.
        at: [n]?usize = @splat(null),
        diffs: [n]?Diff = @splat(null),
        mine: Rstd = .{},
        count: usize = 0,

        /// Enters column `col`'s contestant.
        pub fn add(self: *@This(), col: usize, ctx: Ctx) *Ctx {
            self.ctxs[self.count] = ctx;
            self.at[col] = self.count;
            self.count += 1;
            return &self.ctxs[self.count - 1];
        }

        /// Races them; each column's rate, `work` (bytes or flops) per call
        /// in billions per second, null where it was not entered.
        pub fn race(self: *@This(), io: Io, timing: bench.Timing, work: f64) ![n]?f64 {
            var calls: [n]bench.Call = undefined;
            for (calls[0..self.count], self.ctxs[0..self.count]) |*c, *x| c.* = .of(x);
            var secs: [n]f64 = undefined;
            try bench.race(io, timing, calls[0..self.count], &secs);
            if (lib_err) |err| return err;
            var rates: [n]?f64 = undefined;
            for (self.at, &rates) |j, *rt| rt.* = if (j) |k| work / secs[k] / 1e9 else null;
            return rates;
        }

        /// The row after its lead.
        pub fn print(self: *const @This(), cols: []const Col, rates: []const ?f64) void {
            printRow(cols, rates, &self.diffs, self.mine);
        }
    };
}

/// What rstd's own table holds in slot `name`: converted `asm`, or
/// `portable`.
fn rstdRuns(comptime name: []const u8) []const u8 {
    if (bb.alone) return "-";
    return if (@field(blas.vtable().*, name).? == @field(bb.tables.portable, name).?) "portable" else "asm";
}

const diff_note = "`rstd-<lib>` is rstd's largest difference from that library's result, in eps of the type relative to the size of the terms; `!` marks one past what accumulation order explains: 0 for `copy`, `swap` and `i?amax`, 8 for the other element-wise routines, 2 (n + 8) for a sum of n terms. `rstd err` is rstd's own error from an f64 reference in the same units (gemm: at 256 entries), which tells whose result a difference comes from. `rstd bits` hashes rstd's output, the same on machines that run the same kernels";

pub fn libLevel1(io: Io, gpa: std.mem.Allocator) !void {
    sectionMark("lib1");
    try clock(io, "library level 1 start");
    std.debug.print(
        \\
        \\GB/s, best of {d} batches of at least {d} ms; ratios are over OpenBLAS. `rstd` is the public routine{s}. {s} {s}
        \\
        \\
    , .{ l1_timing.batches, @divTrunc(l1_timing.min_ns, std.time.ns_per_ms), if (arm) ", `asm` the converted kernel alone (`from`)" else "", io_note, diff_note });
    if (arm)
        printHeader("routine | p | n | inc | rstd runs | from", "------- | - | -: | --: | --------- | ----", &lib_cols)
    else
        printHeader("routine | p | n | inc | rstd runs", "------- | - | -: | --: | ---------", &lib_cols);
    inline for (precs) |p| {
        const len = l1_cases[3].bytes / @sizeOf(p.F);
        const b = try Bufs(p.F).init(gpa, len, 1);
        defer b.deinit(gpa);
        // rstd's output
        const kx = try gpa.alloc(p.F, b.x.len);
        defer gpa.free(kx);
        const ky = try gpa.alloc(p.F, b.y.len);
        defer gpa.free(ky);
        inline for (comptime std.meta.tags(L1)) |r| try raceLibL1(io, r, p, b, kx, ky);
    }
}

/// rstd's largest difference from a library's result: `b` holds the
/// library's output, `kx` and `ky` rstd's; `rs` and `ls` are reductions'.
fn diffL1(comptime r: L1, comptime p: Prec, n: usize, inc: usize, b: Bufs(p.F), kx: []const p.F, ky: []const p.F, rs: [2]f64, ls: [2]f64) Diff {
    const eps: f64 = std.math.floatEps(p.F);
    switch (r) {
        .iamax => return .{ .eps = if (rs[0] == ls[0]) 0 else std.math.inf(f64), .bound = 0 },
        .dot, .nrm2, .asum => {
            const s = l1Scalar(r, p, b, n, inc).scale;
            return .{ .eps = @max(@abs(rs[0] - ls[0]), @abs(rs[1] - ls[1])) / (eps * @max(s, 1e-300)), .bound = sumBound(n) };
        },
        .axpy, .scal, .copy, .swap, .rot => {
            var worst: f64 = 0;
            for (0..n) |k| {
                for (0..p.cs()) |q| {
                    const w = l1Elem(r, p, b, k, inc, q);
                    const at = k * inc * p.cs() + q;
                    const dx = @abs(@as(f64, kx[at]) - b.x[at]) / (eps * @max(w[2], 1e-300));
                    const dy = @abs(@as(f64, ky[at]) - b.y[at]) / (eps * @max(w[3], 1e-300));
                    worst = @max(worst, dx, dy);
                }
            }
            return .{ .eps = worst, .bound = if (r == .copy or r == .swap) 0 else 8 };
        },
    }
}

fn raceLibL1(io: Io, comptime r: L1, comptime p: Prec, b: Bufs(p.F), kx: []p.F, ky: []p.F) !void {
    const name = comptime l1Slot(r, p);
    const C = L1Call(r, p);
    const asm_k = if (arm and !bb.alone) asmFor(name) else null;
    const fns = [lib_cols.len]?Fn(name){
        libL1(.openblas, r, p),
        libL1(.rstd, r, p),
        libL1(.blis, r, p),
        libL1(.accelerate, r, p),
        if (asm_k) |a| a.f else null,
        if (arm) candidate(r, p, 4) else null,
    };
    for (l1_cases) |case| {
        const n = case.bytes / p.es() / case.inc;
        if (n > n_max) continue;
        const span = n * case.inc * p.cs();
        var es: Entrants(C, lib_cols.len) = .{};
        var ok = true;
        var rs: [2]f64 = .{ 0, 0 };
        for (0..lib_cols.len) |o| {
            const i = rstdFirst(o);
            if (!lib_cols[i].shown()) continue;
            const f = fns[i] orelse continue;
            const e = checkL1(r, p, f, n, case.inc, b);
            if (lib_err) |err| return err;
            if (e > 1) {
                std.debug.print("WRONG: {s} {s} {s} n={d} inc={d}: {d:.2} of tolerance\n", .{ @tagName(r), p.name, lib_cols[i].label, n, case.inc, e });
                ok = false;
            }
            const res: [2]f64 = if (comptime reduces(r)) l1Result(r, p, f, n, case.inc, b) else .{ 0, 0 };
            if (i == rstd_col) {
                rs = res;
                var h: std.hash.Wyhash = .init(0);
                if (comptime reduces(r)) h.update(std.mem.asBytes(&res)) else {
                    @memcpy(kx[0..span], b.x[0..span]);
                    @memcpy(ky[0..span], b.y[0..span]);
                    h.update(std.mem.sliceAsBytes(b.x[0..span]));
                    h.update(std.mem.sliceAsBytes(b.y[0..span]));
                }
                es.mine = .{ .err = e * @as(f64, if (comptime reduces(r)) @floatFromInt(n + 8) else 8), .bits = h.final() };
            } else if (lib_cols[i].lib != null) es.diffs[i] = diffL1(r, p, n, case.inc, b, kx, ky, rs, res);
            _ = es.add(i, .{ .f = f, .n = @intCast(n), .inc = @intCast(case.inc), .x = b.x.ptr, .y = b.y.ptr });
        }
        if (!ok) continue;
        b.reset();
        const rates = try es.race(io, l1_timing, l1Bytes(r, n, p.es()));
        std.debug.print("| {s} | {s} | {d} | {d} | {s} | ", .{ @tagName(r), p.name, n, case.inc, rstdRuns(name) });
        if (arm) std.debug.print("{s} | ", .{if (asm_k) |a| a.from else "none"});
        es.print(&lib_cols, &rates);
    }
}

pub fn libLevel2(io: Io, gpa: std.mem.Allocator) !void {
    sectionMark("lib2");
    try clock(io, "library level 2 start");
    std.debug.print(
        \\
        \\GFLOP/s, square, unit strides, `beta` 1, best of {d} batches of at least {d} ms; ratios are over OpenBLAS. {s} {s}
        \\
        \\
    , .{ l1_timing.batches, @divTrunc(l1_timing.min_ns, std.time.ns_per_ms), io_note, diff_note });
    if (arm)
        printHeader("routine | p | n | rstd runs | from", "------- | - | -: | --------- | ----", &lib_cols)
    else
        printHeader("routine | p | n | rstd runs", "------- | - | -: | ---------", &lib_cols);
    inline for (precs) |p| {
        inline for (.{ "n", "t" }) |t| try raceLibGemv(io, gpa, p, t);
    }
    try extTable(io, gpa, .l2, ext_l2);
}

/// The size of the terms of each float of `y0 + alpha*op(A)*x`, as
/// `checkGemv` has it.
fn gemvScale(comptime p: Prec, trans: bool, n: usize, a: []const p.F, x: []const p.F, y0: []const p.F, out: []f64) void {
    const cs = p.cs();
    const g = @abs(rounded(p, gemv_alpha));
    const al = if (p.complex) 2 * g else g;
    for (0..n) |i| {
        var mag: f64 = 0;
        for (0..n) |l| {
            const ai = if (trans) (l + i * n) * cs else (i + l * n) * cs;
            mag += magAt(p, a, ai) * magAt(p, x, l * cs);
        }
        for (0..cs) |q| out[i * cs + q] = @abs(y0[i * cs + q]) + al * mag;
    }
}

/// `|re| + |im|` of the element at float `i`.
fn magAt(comptime p: Prec, v: []const p.F, i: usize) f64 {
    return if (p.complex) @as(f64, @abs(v[i])) + @abs(v[i + 1]) else @abs(v[i]);
}

/// The largest `|want - got| / (eps * scale)` over `scale`'s floats.
fn diffScaled(comptime F: type, want: []const F, got: []const F, scale: []const f64) f64 {
    const eps: f64 = std.math.floatEps(F);
    var worst: f64 = 0;
    for (want[0..scale.len], got[0..scale.len], scale) |w, g, s| worst = @max(worst, @abs(@as(f64, w) - g) / (eps * @max(s, 1e-300)));
    return worst;
}

fn raceLibGemv(io: Io, gpa: std.mem.Allocator, comptime p: Prec, comptime t: []const u8) !void {
    const F = p.F;
    const name = p.name ++ "gemv_" ++ t;
    const C = GemvCall(p, name);
    const cs = p.cs();
    const max_n: usize = std.math.sqrt(l2_bytes[2] / p.es());
    const a = try gpa.alignedAlloc(F, .@"64", max_n * max_n * cs);
    defer gpa.free(a);
    const x = try gpa.alignedAlloc(F, .@"64", max_n * cs);
    defer gpa.free(x);
    const y = try gpa.alignedAlloc(F, .@"64", max_n * cs + 64);
    defer gpa.free(y);
    const y0 = try gpa.alloc(F, y.len);
    defer gpa.free(y0);
    const ky = try gpa.alloc(F, y.len);
    defer gpa.free(ky);
    const scale = try gpa.alloc(f64, max_n * cs);
    defer gpa.free(scale);
    const buf = try gpa.alignedAlloc(F, .@"64", 4096);
    defer gpa.free(buf);
    var rng = std.Random.DefaultPrng.init(2);
    const rand = rng.random();
    fillRandom(F, a, rand);
    fillRandom(F, x, rand);
    fillRandom(F, y0, rand);

    const asm_k = if (arm and !bb.alone) asmFor(name) else null;
    const fns = [lib_cols.len]?Fn(name){
        libGemv(.openblas, p, t),
        libGemv(.rstd, p, t),
        libGemv(.blis, p, t),
        libGemv(.accelerate, p, t),
        if (asm_k) |k| k.f else null,
        if (arm) gemvCandidate(p, t, 4) else null,
    };
    for (l2_bytes) |bytes| {
        const n: usize = std.math.sqrt(bytes / p.es());
        if (n > n_max) continue;
        const out = n * cs;
        gemvScale(p, t[0] == 't', n, a, x, y0, scale[0..out]);
        var es: Entrants(C, lib_cols.len) = .{};
        var ok = true;
        for (0..lib_cols.len) |o| {
            const i = rstdFirst(o);
            if (!lib_cols[i].shown()) continue;
            const f = fns[i] orelse continue;
            const ctx = es.add(i, .{ .f = f, .n = @intCast(n), .a = a.ptr, .x = x.ptr, .y = y.ptr, .buf = buf.ptr });
            @memcpy(y, y0);
            ctx.call();
            if (lib_err) |err| return err;
            const e = checkGemv(p, t[0] == 't', n, a, x, y, y0);
            if (e > 1) {
                std.debug.print("WRONG: gemv_{s} {s} {s} n={d}: {d:.2} of tolerance\n", .{ t, p.name, lib_cols[i].label, n, e });
                ok = false;
            }
            if (i == rstd_col) {
                @memcpy(ky, y);
                es.mine = .{ .err = e * @as(f64, @floatFromInt(n + 8)), .bits = std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(y[0..out])) };
            } else if (lib_cols[i].lib != null) es.diffs[i] = .{ .eps = diffScaled(F, ky, y, scale[0..out]), .bound = sumBound(n) };
        }
        if (!ok) continue;
        @memcpy(y, y0);
        const nf: f64 = @floatFromInt(n);
        const rates = try es.race(io, l1_timing, nf * nf * @as(f64, if (p.complex) 8 else 2));
        std.debug.print("| gemv_{s} | {s} | {d} | {s} | ", .{ t, p.name, n, rstdRuns(name) });
        if (arm) std.debug.print("{s} | ", .{if (asm_k) |k| k.from else "none"});
        es.print(&lib_cols, &rates);
    }
}

/// `GemmCall`, which is real only and needs the `.S`.
fn Loops(comptime p: Prec) type {
    return if (p.complex or !ob_asm) void else GemmCall(p);
}

/// One gemm contestant: a library's, or OpenBLAS's `.S` through the bench's
/// copy of `level3.c`'s loops (`GemmCall`).
fn LibGemmCall(comptime p: Prec) type {
    return struct {
        f: ?GemmFn(p),
        loops: ?Loops(p),
        n: usize,
        a: [*]p.F,
        b: [*]p.F,
        c: [*]p.F,
        fn call(s: *const @This()) void {
            if (s.f) |f| return f(s.n, s.a, s.b, s.c);
            if (comptime Loops(p) != void) s.loops.?.call();
        }
        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// The size of the terms of each float of `c0 + A*B`, square `n`:
/// `|c0| + |A| |B|`, a complex's size `|re| + |im|`, by rstd's dgemm.
/// `abs_a` and `abs_b` are scratch.
fn gemmScale(comptime p: Prec, n: usize, a: []const p.F, b: []const p.F, c0: []const p.F, abs_a: []f64, abs_b: []f64, out: []f64) !void {
    const cs = p.cs();
    const sq = n * n;
    for (0..sq) |i| {
        abs_a[i] = magAt(p, a, i * cs);
        abs_b[i] = magAt(p, b, i * cs);
        out[i] = magAt(p, c0, i * cs);
    }
    const ni: isize = @intCast(n);
    try blas.gemm(rstd_alloc, rstd_io, f64, .n, .n, ni, ni, ni, 1, abs_a.ptr, ni, abs_b.ptr, ni, 1, out.ptr, ni);
    // a complex's two floats share its scale, spread from the back
    if (p.complex) {
        var i = sq;
        while (i > 0) {
            i -= 1;
            out[2 * i + 1] = out[i];
            out[2 * i] = out[i];
        }
    }
}

pub fn libLevel3(io: Io, gpa: std.mem.Allocator) !void {
    sectionMark("lib3");
    try clock(io, "library level 3 start");
    std.debug.print(
        \\
        \\gemm NN, square, `C += A*B`: GFLOP/s, best of {d} batches of at least {d} ms; ratios are over OpenBLAS{s}. {s} {s}
        \\
        \\
    , .{ l3_timing.batches, @divTrunc(l3_timing.min_ns, std.time.ns_per_ms), if (ob_asm) ". `ob .S` is OpenBLAS's vortex `.S` through the bench's loops" else "", io_note, diff_note });
    printHeader("p | n", "- | -:", &lib3_cols);
    inline for (precs) |p| try raceLibGemm(io, gpa, p);
    try extTable(io, gpa, .l3, ext_l3);
}

fn raceLibGemm(io: Io, gpa: std.mem.Allocator, comptime p: Prec) !void {
    const F = p.F;
    const cs = p.cs();
    const max_n = gemm_sizes[gemm_sizes.len - 1];
    const len = max_n * max_n * cs;
    const a = try gpa.alignedAlloc(F, .@"64", len);
    defer gpa.free(a);
    const b = try gpa.alignedAlloc(F, .@"64", len);
    defer gpa.free(b);
    const c = try gpa.alignedAlloc(F, .@"64", len);
    defer gpa.free(c);
    const c0 = try gpa.alloc(F, len);
    defer gpa.free(c0);
    // rstd's output, and the size of each float's terms
    const kc = try gpa.alloc(F, len);
    defer gpa.free(kc);
    const scale = try gpa.alloc(f64, len);
    defer gpa.free(scale);
    const abs_a = try gpa.alloc(f64, max_n * max_n);
    defer gpa.free(abs_a);
    const abs_b = try gpa.alloc(f64, max_n * max_n);
    defer gpa.free(abs_b);
    var rng = std.Random.DefaultPrng.init(4);
    const rand = rng.random();
    for (a, b, c0) |*x, *y, *z| {
        x.* = @floatCast(rand.float(f64) * 2 - 1);
        y.* = @floatCast(rand.float(f64) * 2 - 1);
        z.* = @floatCast(rand.float(f64) * 2 - 1);
    }

    // the `.S` through the bench's loops, real only
    var sa: []align(64) F = &.{};
    var sb: []align(64) F = &.{};
    var blk: Blocking = undefined;
    if (comptime Loops(p) != void) {
        blk = .of(p, asmGemm(p));
        sa = try gpa.alignedAlloc(F, .@"64", blk.p * (blk.q + blk.um));
        sb = try gpa.alignedAlloc(F, .@"64", (blk.q + blk.um) * blk.r);
    }
    defer if (comptime Loops(p) != void) {
        gpa.free(sa);
        gpa.free(sb);
    };

    const fns = [lib3_cols.len]?GemmFn(p){ libGemm(.openblas, p), libGemm(.rstd, p), libGemm(.blis, p), libGemm(.accelerate, p), null };
    const LC = LibGemmCall(p);
    for (gemm_sizes) |n| {
        if (n > n_max) continue;
        const m = n * n * cs;
        // rstd's dgemm; nothing of rstd's runs in `alone`, where nothing is compared
        if (!bb.alone) try gemmScale(p, n, a[0..m], b[0..m], c0[0..m], abs_a, abs_b, scale[0..m]);
        var es: Entrants(LC, lib3_cols.len) = .{};
        for (0..lib3_cols.len) |o| {
            const i = rstdFirst(o);
            if (!lib3_cols[i].shown()) continue;
            const loops: ?Loops(p) = if (comptime Loops(p) != void) (if (i == 4) .{ .im = asmGemm(p), .blk = blk, .n = n, .a = a.ptr, .b = b.ptr, .c = c.ptr, .sa = sa.ptr, .sb = sb.ptr } else null) else null;
            if (fns[i] == null and loops == null) continue;
            const ctx = es.add(i, .{ .f = fns[i], .loops = loops, .n = n, .a = a.ptr, .b = b.ptr, .c = c.ptr });
            @memcpy(c[0..m], c0[0..m]);
            ctx.call();
            if (lib_err) |err| return err;
            const e = checkGemm(p, n, n, n, a[0..m], b[0..m], c[0..m], c0[0..m], 256);
            if (e > 1) {
                std.debug.print("WRONG: {s}gemm {s} n={d}: {d:.2} of tolerance\n", .{ p.name, lib3_cols[i].label, n, e });
                return error.WrongGemm;
            }
            if (i == rstd_col) {
                @memcpy(kc[0..m], c[0..m]);
                es.mine = .{ .err = e * @as(f64, @floatFromInt(n + 8)), .bits = std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(c[0..m])) };
            } else if (lib3_cols[i].lib != null) es.diffs[i] = .{ .eps = diffScaled(F, kc, c, scale[0..m]), .bound = sumBound(n) };
        }
        const nf: f64 = @floatFromInt(n);
        const rates = try es.race(io, l3_timing, nf * nf * nf * @as(f64, if (p.complex) 8 else 2));
        std.debug.print("| {s} | {d} | ", .{ p.name, n });
        es.print(&lib3_cols, &rates);
    }
}
