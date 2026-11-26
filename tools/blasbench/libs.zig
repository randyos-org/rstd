//! The library races: `lib1`, `lib2`, `lib3`, and the loading of each library.

const std = @import("std");
const sysinfo = @import("builtin");
const rstd = @import("rstd");
const bb = @import("blas_bench.zig");
const tb = @import("libs/table.zig");
const l1r = @import("libs/level1.zig");
const l2r = @import("libs/level2.zig");
const l3r = @import("libs/level3.zig");

const Io = std.Io;
const blas = rstd.math.kernels.blas;
const Dll = rstd.dll.Dll;
const os_tag = sysinfo.target.os.tag;
const Prec = bb.Prec;
pub const Col = tb.Col;
pub const Entrants = tb.Entrants;
pub const lib_cols = tb.lib_cols;
pub const merge_labels = tb.merge_labels;
pub const printHeader = tb.printHeader;
pub const printLead = tb.printLead;
pub const rstdFirst = tb.rstdFirst;
pub const rstd_col = tb.rstd_col;
pub const sectionMark = tb.sectionMark;
pub const sumBound = tb.sumBound;
pub const libLevel1 = l1r.libLevel1;
pub const libLevel2 = l2r.libLevel2;
pub const libLevel3 = l3r.libLevel3;

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
pub fn raced(src: Src) bool {
    return switch (src) {
        .openblas, .blis, .accelerate => true,
        .rstd, .reference, .eigen, .aoclflame => false,
    };
}

/// Each library's handle; null when not given, and for `rstd`.
pub var lib_handles: std.EnumArray(Src, ?Dll) = .initFill(null);

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
    /// The bench's `Io.Threaded`, with a program's default limits.
    threaded,
    /// Serial: `Io.failing` runs every task inline.
    failing,
};

/// And an `Io`: `main` sets it.
pub var rstd_io: Io = .failing;

/// Which `Io` rstd ran on, for each library table's notes.
pub var io_note: []const u8 = "";
pub var io_note_buf: [96]u8 = undefined;

/// `--openblas-prefix`: the prefix on each OpenBLAS symbol (`scipy_` in
/// SciPy's wheels).
pub var ob_prefix: []const u8 = "";

/// Whether the bench set `OPENBLAS_THREAD_TIMEOUT` to idle OpenBLAS's workers.
pub var ob_idle = false;

/// `--max-n`: rows of a larger `n` are skipped.
pub var n_max: usize = std.math.maxInt(usize);

/// `--only`: the routines the table sections race, comma separated; null for all.
pub var only: ?[]const u8 = null;

/// Whether `--only` leaves `name` in.
pub fn wanted(name: []const u8) bool {
    const list = only orelse return true;
    var it = std.mem.splitScalar(u8, list, ',');
    while (it.next()) |x| if (std.mem.eql(u8, x, name)) return true;
    return false;
}

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

pub fn setEnv(name: [*:0]const u8, value: [*:0]const u8) void {
    if (os_tag == .windows) _ = SetEnvironmentVariableA(name, value) else _ = setenv(name, value, 1);
}

/// The loader's reason for a failed open, where it gives one.
pub fn dlError() []const u8 {
    if (os_tag == .windows) return "";
    return if (std.c.dlerror()) |m| std.mem.span(m) else "";
}

/// `name` from `src`'s library: Accelerate's current LP64 entry
/// (`ACCELERATE_NEW_LAPACK`) where it has one, else the classic one.
pub fn lookup(comptime src: Src, comptime T: type, comptime name: [:0]const u8) ?T {
    const h = &(lib_handles.getPtr(src).* orelse return null);
    if (src == .accelerate) {
        if (h.lookup(T, name ++ "$NEWLAPACK")) |f| return f;
    }
    if (src == .openblas) return lookupPrefixed(h, T, name);
    return h.lookup(T, name);
}

/// `name` after `ob_prefix`, from OpenBLAS's library `h`.
pub fn lookupPrefixed(h: *Dll, comptime T: type, name: []const u8) ?T {
    var buf: [128]u8 = undefined;
    const full = std.mem.printSentinel(&buf, "{s}{s}", .{ ob_prefix, name }, 0) catch return null;
    return h.lookup(T, full);
}

/// `lookup`, or `error.MissingSymbol` after naming it: never a null call (#649).
fn need(comptime src: Src, comptime T: type, comptime name: [:0]const u8) error{MissingSymbol}!T {
    if (lookup(src, T, name)) |f| return f;
    const prefix = if (src == .openblas) ob_prefix else "";
    std.debug.print("{t} has no `{s}{s}`", .{ src, prefix, name });
    if (src == .openblas) std.debug.print("; is its symbol prefix (--openblas-prefix) right?", .{});
    std.debug.print("\n", .{});
    return error.MissingSymbol;
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
    if (src == .openblas) return lookupPrefixed(h, T, base ++ "_");
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
        // SciPy's wheels prefix every symbol
        if (ob_prefix.len == 0 and lookup(.openblas, Str, "openblas_get_config") == null) {
            if (lib_handles.getPtr(.openblas).*.?.lookup(Str, "scipy_openblas_get_config") != null) {
                ob_prefix = "scipy_";
                std.debug.print("openblas: symbols prefixed `scipy_`, as with --openblas-prefix=scipy_\n", .{});
            }
        }
        const config = std.mem.span((try need(.openblas, Str, "openblas_get_config"))());
        // `int` indices here; an ILP64 build (`USE64BITINT`) would need `long`
        if (std.mem.find(u8, config, "USE64BITINT") != null) return error.OpenblasIsIlp64;
        // Windows' CRT may have copied the environment already: by call too
        if (opts.threads) |t| (try need(.openblas, *const fn (c_int) callconv(.c) void, "openblas_set_num_threads"))(t);
        std.debug.print("openblas: {s}\nopenblas core: {s}, config: {s}, threads: {d}, lapack: {s}\n", .{
            paths.get(.openblas).?,
            (try need(.openblas, Str, "openblas_get_corename"))(),
            config,
            (try need(.openblas, Int, "openblas_get_num_threads"))(),
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
        if (opts.threads) |t| (try need(.blis, *const fn (i64) callconv(.c) void, "bli_thread_set_num_threads"))(t);
        const arch = (try need(.blis, *const fn () callconv(.c) c_int, "bli_arch_query_id"))();
        std.debug.print("blis: {s}\nblis version: {s}, config: {s}, threads: {d}, cblas: {s}\n", .{
            paths.get(.blis).?,
            (try need(.blis, Str, "bli_info_get_version_str"))(),
            (try need(.blis, *const fn (c_int) callconv(.c) [*:0]const u8, "bli_arch_string"))(arch),
            (try need(.blis, Dim, "bli_thread_get_num_threads"))(),
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
    // the race idles OpenBLAS's workers by a variable rstd reads too: rstd's
    // pool keeps upstream's default (#710)
    if (ob_idle) blas.config.thread_timeout = (blas.Config{}).thread_timeout;
    io_note = if (ceiling) |c|
        try std.mem.print(&io_note_buf, "rstd runs on `{s}`, thread ceiling {d}.", .{ io_name, c })
    else
        try std.mem.print(&io_note_buf, "rstd runs on `{s}`, threads up to the CPUs.", .{io_name});
    // rstd's own load of Accelerate makes this thread single-threaded
    if (loaded(.accelerate)) std.debug.print("accelerate after rstd's init: BLASGetThreading: {s}\n", .{accelThreading()});
}

/// `src`'s LAPACK version (`ilaver`), or `-` without one.
pub fn lapackVersion(comptime src: Src, buf: []u8) []const u8 {
    const f = lookupF(src, *const fn (*c_int, *c_int, *c_int) callconv(.c) void, "ilaver") orelse return "-";
    var v: [3]c_int = undefined;
    f(&v[0], &v[1], &v[2]);
    return std.mem.print(buf, "{d}.{d}.{d}", .{ v[0], v[1], v[2] }) catch "?";
}

/// The calling thread's `BLASGetThreading`, as words.
pub fn accelThreading() []const u8 {
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

/// rstd's view of the `n` elements `inc` apart that the bench passes from `x`.
pub fn view(comptime T: type, x: [*]T, n: isize, inc: isize) blas.Strided(T) {
    const count: usize = @intCast(@max(n, 0));
    const span = if (count == 0) 0 else (count - 1) * @abs(inc) + 1;
    return blas.Strided(T).init(x[0..span], count, inc) catch unreachable;
}

/// The bench's `rows` by `cols` matrix at `a`, leading dimension `ld`, as a slice.
pub fn matrix(comptime T: type, a: [*]T, rows: isize, cols: isize, ld: isize) []T {
    if (rows <= 0 or cols <= 0) return a[0..0];
    return a[0..@intCast((cols - 1) * ld + rows)];
}

/// cblas's `CBLAS_ORDER` and `CBLAS_TRANSPOSE`.
pub const col_major: c_int = 102;
pub const no_trans: c_int = 111;
pub const trans_t: c_int = 112;
