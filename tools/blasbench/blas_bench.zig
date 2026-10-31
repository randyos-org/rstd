//! rstd's BLAS kernels raced on Apple silicon, each in its vtable slot's
//! kernel ABI (`bench.zig`): the `vortex` object's kernel (converted
//! assembly; `armv8`'s where vortex has none of its own), the portable one,
//! and `@Vector` candidates (`blas.vector`). gemm races the Zig kernel
//! of each object against OpenBLAS's own vortex `.S`, assembled from a
//! checkout by the build. Prints Markdown table rows (`results-m6.md`).
//!
//!     zig build ssh-bench-aarch64-macos [-Dopenblas=<checkout>] [-- <section>...]
//!     zig build bench-native -Dopenblas-lib=<library> [-Dblis-lib=<library>] [-- <section>...]
//!     zig build qemu-bench-linux-<arch> [-Dvm-bench-arg=<section>...]
//!
//! Sections: `probe`, `l1`, `l2`, `l3`; all when none is named. arm64
//! only; `l3` only with the `.S` (`-Dopenblas`). `kern3`, only when named
//! and on every host: the selected core's gemm kernel alone against rstd's
//! whole `gemm` on one worker, every type (T21).
//!
//! `lib1`, `lib2`, `lib3` race rstd's public routines (`blas.dot`, ...)
//! against real OpenBLAS, BLIS and Accelerate, loaded by absolute path
//! (`results-m6-libraries.md`, `results-flyoff.md`), and give rstd's
//! largest difference from each library's result; on arm64 and x86_64,
//! Windows, Linux and macOS. `lapack` races rstd's LAPACK the same way.
//! `lib2`, `lib3` and `lapack` also hold the table-driven races
//! (`ext_l2`, `ext_l3`, `lapack_ops`), which take the reference binaries
//! too: reference LAPACK, Eigen and AOCL-libFLAME (`readme.md`). With the
//! others when none is named and `--openblas` is given. Options:
//! `--openblas=<library>` (`-Dopenblas-lib`), `--blis=<library>`,
//! `--accelerate=<path>` (macOS' framework by default),
//! `--reference=<liblapack>`, `--eigen=<shim>`, `--aocl-libflame=<library>`
//! (each skipped with a warning when it cannot load), `--threads=<n>` or
//! `--threads=default` (each library and rstd's ceiling; 1 by default),
//! `--rstd-io=failing` (rstd on `Io.failing`, serial, rather than its own
//! `Io.Threaded`), `--max-n=<n>` (library rows of `n` up to that only),
//! `--ob-core=<core>` (`OPENBLAS_CORETYPE`, which rstd reads too),
//! `--rstd-accelerate=0` (rstd's own kernels on macOS, not Accelerate),
//! `--accelerate-alone=1` (only the libraries, rstd never initialized: a
//! threaded run on macOS starts itself so for a fair Accelerate column).
//!
//!     blas_bench --merge <out.md> <log>...
//!
//! merges runs of one configuration (`merge.zig`, `flyoff.sh`).

const std = @import("std");
const sysinfo = @import("builtin");
const Io = std.Io;
const rstd = @import("rstd");
const bench_options = @import("bench_options");
const bench = @import("bench.zig");
const merge = @import("merge.zig");

const blas = rstd.linalg.blas;
const lapack = rstd.linalg.lapack;
const Vtable = blas.Vtable;
const Dll = rstd.dll.Dll;
const os_tag = sysinfo.target.os.tag;

/// The arm64 cores' own races: `probe`, `l1`, `l2`, and the `asm` and
/// `vec x4` library columns.
const arm = sysinfo.target.cpu.arch == .aarch64;
/// OpenBLAS's vortex gemm `.S` is linked (`-Dopenblas`): `l3`, `ob .S`.
const ob_asm = bench_options.ob_asm;

comptime {
    if (!arm and sysinfo.target.cpu.arch != .x86_64) @compileError("blas_bench runs on arm64 and x86_64");
    if (ob_asm and !arm) @compileError("the vortex `.S` is arm64");
}

/// The kernel type of vtable slot `name`.
fn Fn(comptime name: []const u8) type {
    return @typeInfo(@FieldType(Vtable, name)).optional.child;
}

const Tables = struct {
    vortex: Vtable,
    vortexm4: Vtable,
    armv8: Vtable,
    portable: Vtable,
};
var tables: Tables = undefined;

const Prec = struct {
    F: type,
    complex: bool,
    name: []const u8,

    fn cs(comptime p: Prec) usize {
        return if (p.complex) 2 else 1;
    }
    /// Bytes of one element.
    fn es(comptime p: Prec) usize {
        return p.cs() * @sizeOf(p.F);
    }
};

const precs = [_]Prec{
    .{ .F = f32, .complex = false, .name = "s" },
    .{ .F = f64, .complex = false, .name = "d" },
    .{ .F = f32, .complex = true, .name = "c" },
    .{ .F = f64, .complex = true, .name = "z" },
};

const l1_timing: bench.Timing = .{ .min_ns = 10 * std.time.ns_per_ms, .batches = 15 };
const l3_timing: bench.Timing = .{ .min_ns = 20 * std.time.ns_per_ms, .batches = 7 };

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const gpa = init.gpa;
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    // before anything touches rstd
    if (args.len > 1 and std.mem.eql(u8, args[1], "--merge")) return merge.run(io, gpa, args[2..], &merge_labels);

    if (comptime os_tag.isDarwin()) {
        // the performance cores; a `taskpolicy` clamp still overrides it
        if (pthread_set_qos_class_self_np(qos_class_user_interactive, 0) != 0)
            std.debug.print("warning: QoS not set\n", .{});
    }

    // `--name=value` options, then section names
    var opts: LibOptions = .{};
    var sections: std.ArrayList([]const u8) = .empty;
    for (args[1..]) |a| {
        if (std.mem.startsWith(u8, a, "--")) {
            const eq = std.mem.findScalar(u8, a, '=') orelse return error.BadOption;
            const v = a[eq + 1 ..];
            const k = std.meta.stringToEnum(Opt, a[2..eq]) orelse return error.BadOption;
            switch (k) {
                .openblas => opts.paths[0] = v,
                .blis => opts.paths[1] = v,
                .accelerate => opts.paths[2] = v,
                .reference => opts.paths[3] = v,
                .eigen => opts.paths[4] = v,
                .@"aocl-libflame" => opts.paths[5] = v,
                .threads => {
                    opts.threads = if (std.mem.eql(u8, v, "default")) null else std.fmt.parseInt(u31, v, 10) catch return error.BadOption;
                    if (opts.threads == 0) return error.BadOption;
                },
                .@"rstd-io" => opts.rstd_io = std.meta.stringToEnum(RstdIo, v) orelse return error.BadOption,
                .@"max-n" => n_max = std.fmt.parseInt(usize, v, 10) catch return error.BadOption,
                .@"ob-core" => opts.ob_core = v,
                .@"rstd-accelerate" => opts.rstd_accelerate = !std.mem.eql(u8, v, "0"),
                .@"accelerate-alone" => alone = !std.mem.eql(u8, v, "0"),
            }
        } else try sections.append(init.arena.allocator(), a);
    }

    // the parent's rstd set `VECLIB_MAXIMUM_THREADS=1`, which this process
    // inherited: Accelerate takes its own count unless `--threads` sets one
    if (comptime os_tag.isDarwin()) {
        if (alone and opts.threads == null) _ = unsetenv("VECLIB_MAXIMUM_THREADS");
    }
    // rstd's tables select nothing, but the alone run leaves rstd untouched
    if (!alone) {
        tables.portable = blas.cpu.generic.vtableFromCoreType(.generic);
        if (arm) {
            tables.vortex = blas.cpu.vtableFromCoreType(.vortex);
            tables.vortexm4 = blas.cpu.vtableFromCoreType(.vortexm4);
            tables.armv8 = blas.cpu.vtableFromCoreType(.armv8);
            if (tables.vortex.daxpy_k == tables.portable.daxpy_k) return error.VortexNotLinked;
        }
        std.debug.print("zig {s}, {t}-{t}, cpu model {s}, detected core {t}\n", .{ sysinfo.zig_version_string, sysinfo.target.cpu.arch, os_tag, sysinfo.target.cpu.model.name, blas.cpu.detect() });
    } else std.debug.print("zig {s}, {t}-{t}, cpu model {s}; the libraries alone, rstd never initialized\n", .{ sysinfo.zig_version_string, sysinfo.target.cpu.arch, os_tag, sysinfo.target.cpu.model.name });
    // before anything selects rstd's core, which reads its ceiling then
    try setThreadEnv(opts.threads);
    var threaded: Io.Threaded = .init(std.heap.c_allocator, if (opts.threads) |t| .{ .async_limit = .limited(t - 1), .concurrent_limit = .limited(t - 1) } else .{});
    defer threaded.deinit();
    if (opts.rstd_io == .threaded) rstd_io = threaded.io();
    const all = sections.items.len == 0;
    const want = struct {
        fn f(list: []const []const u8, name: []const u8) bool {
            for (list) |a| if (std.mem.eql(u8, a, name)) return true;
            return false;
        }
    }.f;
    if (!arm) {
        for ([_][]const u8{ "probe", "l1", "l2", "l3" }) |s| if (want(sections.items, s)) return error.ArmOnlySection;
    }
    if (!ob_asm and want(sections.items, "l3")) return error.NeedsOpenblasAsm;
    if (arm and !alone) {
        if (all or want(sections.items, "probe")) try probe(io);
        if (all or want(sections.items, "l1")) try level1(io, gpa);
        if (all or want(sections.items, "l2")) try level2(io, gpa);
    }
    if (ob_asm and !alone) {
        if (all or want(sections.items, "l3")) try level3(io, gpa);
    }
    // only when named: it times rstd alone, which the flyoff's runs don't want
    if (!alone and want(sections.items, "kern3")) try ownLevel3(io, gpa);
    const lib_sections = want(sections.items, "lib1") or want(sections.items, "lib2") or want(sections.items, "lib3") or want(sections.items, "lapack");
    if (lib_sections or (all and (opts.paths[0] != null or alone))) {
        try loadLibs(opts);
        if (all or want(sections.items, "lib1")) try libLevel1(io, gpa);
        if (all or want(sections.items, "lib2")) try libLevel2(io, gpa);
        if (all or want(sections.items, "lib3")) try libLevel3(io, gpa);
        if (all or want(sections.items, "lapack")) try libLapack(io, gpa);
        // rstd's init made Accelerate single-threaded here: its own threads
        // run in a process that never initializes rstd (T07)
        if (!alone and opts.threads != 1 and loaded(.accelerate)) try runAlone(io, args);
    }
    try clock(io, "end");
}

/// `--name=value` options.
const Opt = enum { openblas, blis, accelerate, reference, eigen, @"aocl-libflame", threads, @"rstd-io", @"max-n", @"ob-core", @"rstd-accelerate", @"accelerate-alone" };

/// `--accelerate-alone`: only the libraries race; nothing of rstd's runs, so
/// nothing sets `VECLIB_MAXIMUM_THREADS` or `BLASSetThreading`.
var alone = false;

/// This bench again with `args` and `--accelerate-alone=1`: the libraries'
/// tables, Accelerate on its own threads.
fn runAlone(io: Io, args: []const [:0]const u8) !void {
    var argv: [64][]const u8 = undefined;
    if (args.len + 1 > argv.len) return error.TooManyArguments;
    for (args, argv[0..args.len]) |a, *v| v.* = a;
    argv[args.len] = "--accelerate-alone=1";
    std.debug.print("\n{s}: the libraries again in a process without rstd, Accelerate on its own threads\n", .{args[0]});
    var child = try std.process.spawn(io, .{ .argv = argv[0 .. args.len + 1], .stdin = .ignore, .stdout = .inherit, .stderr = .inherit });
    switch (try child.wait(io)) {
        .exited => |code| if (code != 0) return error.AloneRunFailed,
        else => return error.AloneRunFailed,
    }
}

extern "c" fn pthread_set_qos_class_self_np(qos: c_uint, relpri: c_int) c_int;
const qos_class_user_interactive: c_uint = 0x21;

// ---- clock and peak --------------------------------------------------------

/// Eight dependent one-cycle adds per iteration: cycles per second.
fn clockLoop(iters: usize) void {
    var n = iters;
    var acc: u64 = 0;
    if (arm) asm volatile (
        \\1:
        \\ add %[a], %[a], #1
        \\ add %[a], %[a], #1
        \\ add %[a], %[a], #1
        \\ add %[a], %[a], #1
        \\ add %[a], %[a], #1
        \\ add %[a], %[a], #1
        \\ add %[a], %[a], #1
        \\ add %[a], %[a], #1
        \\ subs %[n], %[n], #1
        \\ b.ne 1b
        : [n] "+r" (n),
          [a] "+r" (acc),
        :
        : .{ .nzcv = true }) else asm volatile (
        \\1:
        \\ add $1, %[a]
        \\ add $1, %[a]
        \\ add $1, %[a]
        \\ add $1, %[a]
        \\ add $1, %[a]
        \\ add $1, %[a]
        \\ add $1, %[a]
        \\ add $1, %[a]
        \\ dec %[n]
        \\ jnz 1b
        : [n] "+r" (n),
          [a] "+r" (acc),
        :
        : .{ .cc = true });
}

/// Sixteen independent `fmla` chains per iteration, `.4s` or `.2d`.
fn peakLoop(comptime double: bool, iters: usize) void {
    var n = iters;
    asm volatile (peak_text[@intFromBool(double)]
        : [n] "+r" (n),
        :
        : .{ .nzcv = true, .v0 = true, .v1 = true, .v2 = true, .v3 = true, .v4 = true, .v5 = true, .v6 = true, .v7 = true, .v8 = true, .v9 = true, .v10 = true, .v11 = true, .v12 = true, .v13 = true, .v14 = true, .v15 = true, .v16 = true, .v17 = true });
}

/// `peakLoop`'s body, single then double: v0 to v15 accumulate v16 * v17.
const peak_text = blk: {
    var texts: [2][]const u8 = undefined;
    for (.{ "4s", "2d" }, 0..) |t, k| {
        var s: []const u8 = "";
        for (0..18) |i| s = s ++ std.fmt.comptimePrint(" movi v{d}.16b, #0\n", .{i});
        s = s ++ "1:\n";
        for (0..16) |i| s = s ++ std.fmt.comptimePrint(" fmla v{d}.{s}, v16.{s}, v17.{s}\n", .{ i, t, t, t });
        texts[k] = s ++ " subs %[n], %[n], #1\n b.ne 1b\n";
    }
    break :blk texts;
};

const LoopCtx = struct {
    iters: usize,
    f: *const fn (usize) void,
    fn run(o: *anyopaque) void {
        const c: *LoopCtx = @ptrCast(@alignCast(o));
        c.f(c.iters);
    }
};

/// GHz of the core this thread runs on, from the add chain.
fn ghz(io: Io) !f64 {
    var ctx: LoopCtx = .{ .iters = 1_000_000, .f = clockLoop };
    var out: [1]f64 = undefined;
    try bench.race(io, .{ .min_ns = 5 * std.time.ns_per_ms, .batches = 9 }, &.{.{ .ctx = &ctx, .run = LoopCtx.run }}, &out);
    return 8.0 * @as(f64, @floatFromInt(ctx.iters)) / out[0] / 1e9;
}

fn clock(io: Io, what: []const u8) !void {
    std.debug.print("clock ({s}): {d:.2} GHz\n", .{ what, try ghz(io) });
}

/// Peak `fmla` GFLOP/s, single and double, and the clock.
fn probe(io: Io) !void {
    var s: LoopCtx = .{ .iters = 1_000_000, .f = struct {
        fn f(n: usize) void {
            peakLoop(false, n);
        }
    }.f };
    var d: LoopCtx = .{ .iters = 1_000_000, .f = struct {
        fn f(n: usize) void {
            peakLoop(true, n);
        }
    }.f };
    var out: [2]f64 = undefined;
    try bench.race(io, l1_timing, &.{ .{ .ctx = &s, .run = LoopCtx.run }, .{ .ctx = &d, .run = LoopCtx.run } }, &out);
    const flops = 16.0 * 2.0 * @as(f64, @floatFromInt(s.iters));
    std.debug.print("peak fmla: {d:.1} GFLOP/s single, {d:.1} double\n", .{ flops * 4 / out[0] / 1e9, flops * 2 / out[1] / 1e9 });
    try clock(io, "probe");
}

// ---- level 1 ---------------------------------------------------------------

const L1 = enum { axpy, dot, scal, nrm2, asum, iamax, copy, swap, rot };

fn l1Slot(comptime r: L1, comptime p: Prec) []const u8 {
    return switch (r) {
        .dot => p.name ++ (if (p.complex) "dotu_k" else "dot_k"),
        .iamax => "i" ++ p.name ++ "amax_k",
        .rot => if (p.complex) p.name ++ (if (p.F == f32) "s" else "d") ++ "rot_k" else p.name ++ "rot_k",
        else => p.name ++ @tagName(r) ++ "_k",
    };
}

/// Bytes `r` moves for `n` elements at unit stride.
fn l1Bytes(r: L1, n: usize, es: usize) f64 {
    const vectors: f64 = switch (r) {
        .nrm2, .asum, .iamax => 1,
        .dot, .scal, .copy => 2,
        .axpy => 3,
        .swap, .rot => 4,
    };
    return vectors * @as(f64, @floatFromInt(n * es));
}

const axpy_alpha = 1e-6;
const scal_alpha = -1.0;
// |alpha| = 1, so repeated calls neither grow nor shrink x
const cscal_alpha = .{ 0.6, 0.8 };
const rot_cs = .{ 0.6, 0.8 };

/// `v` as the kernels over `p.F` get it.
fn rounded(comptime p: Prec, comptime v: comptime_float) f64 {
    return @as(p.F, v);
}

fn candidate(comptime r: L1, comptime p: Prec, comptime U: comptime_int) Fn(l1Slot(r, p)) {
    const K = if (p.complex) blas.vector.Complex(p.F, U) else blas.vector.Real(p.F, U);
    return &@field(K, if (r == .dot and p.complex) "dotu" else @tagName(r));
}

fn invoke(comptime r: L1, comptime p: Prec, f: Fn(l1Slot(r, p)), n: isize, inc: isize, x: [*]p.F, y: [*]p.F) void {
    const c = p.complex;
    switch (r) {
        .axpy => std.mem.doNotOptimizeAway(if (c) f(n, 0, 0, axpy_alpha, axpy_alpha, x, inc, y, inc, x, 0) else f(n, 0, 0, axpy_alpha, x, inc, y, inc, x, 0)),
        .dot => std.mem.doNotOptimizeAway(f(n, x, inc, y, inc)),
        .scal => std.mem.doNotOptimizeAway(if (c) f(n, 0, 0, cscal_alpha[0], cscal_alpha[1], x, inc, x, 0, x, 1) else f(n, 0, 0, scal_alpha, x, inc, x, 0, x, 1)),
        .nrm2, .asum, .iamax => std.mem.doNotOptimizeAway(f(n, x, inc)),
        .copy => std.mem.doNotOptimizeAway(f(n, x, inc, y, inc)),
        .swap => std.mem.doNotOptimizeAway(if (c) f(n, 0, 0, 0, 0, x, inc, y, inc, x, 0) else f(n, 0, 0, 0, x, inc, y, inc, x, 0)),
        .rot => std.mem.doNotOptimizeAway(f(n, x, inc, y, inc, rot_cs[0], rot_cs[1])),
    }
}

fn L1Call(comptime r: L1, comptime p: Prec) type {
    return struct {
        f: Fn(l1Slot(r, p)),
        n: isize,
        inc: isize,
        x: [*]p.F,
        y: [*]p.F,
        fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            invoke(r, p, s.f, s.n, s.inc, s.x, s.y);
        }
    };
}

/// x, y and the values they start from, with room past the longest vector.
fn Bufs(comptime F: type) type {
    return struct {
        x: []align(64) F,
        y: []align(64) F,
        x0: []align(64) F,
        y0: []align(64) F,

        const guard = 64;

        fn init(gpa: std.mem.Allocator, len: usize, seed: u64) !@This() {
            var b: @This() = undefined;
            inline for (.{ "x", "y", "x0", "y0" }) |name| @field(b, name) = try gpa.alignedAlloc(F, .@"64", len + guard);
            var rng = std.Random.DefaultPrng.init(seed);
            const rand = rng.random();
            for (b.x0, b.y0) |*a, *c| {
                a.* = @floatCast(rand.float(f64) * 2 - 1);
                c.* = @floatCast(rand.float(f64) * 2 - 1);
            }
            b.reset();
            return b;
        }

        fn deinit(b: @This(), gpa: std.mem.Allocator) void {
            inline for (.{ "x", "y", "x0", "y0" }) |name| gpa.free(@field(b, name));
        }

        fn reset(b: @This()) void {
            @memcpy(b.x, b.x0);
            @memcpy(b.y, b.y0);
        }
    };
}

/// Element `k` of a vector at stride `inc`, widened, as (re, im).
fn elem(comptime p: Prec, v: []const p.F, k: usize, inc: usize) [2]f64 {
    const at = k * inc * p.cs();
    return .{ v[at], if (p.complex) v[at + 1] else 0 };
}

/// Whether float `i` of a vector of `n` at stride `inc` belongs to it.
fn visited(comptime p: Prec, i: usize, n: usize, inc: usize) bool {
    const step = inc * p.cs();
    return i / step < n and i % step < p.cs();
}

/// What reduction `r` should give, (re, im), and the size of its terms.
fn l1Scalar(comptime r: L1, comptime p: Prec, b: Bufs(p.F), n: usize, inc: usize) struct { want: [2]f64, scale: f64 } {
    var re: f64 = 0;
    var im: f64 = 0;
    var mag: f64 = 0;
    for (0..n) |k| {
        const xv = elem(p, b.x0, k, inc);
        const yv = elem(p, b.y0, k, inc);
        switch (r) {
            .dot => {
                re += xv[0] * yv[0] - xv[1] * yv[1];
                im += xv[0] * yv[1] + xv[1] * yv[0];
                mag += (@abs(xv[0]) + @abs(xv[1])) * (@abs(yv[0]) + @abs(yv[1]));
            },
            .nrm2 => re += xv[0] * xv[0] + xv[1] * xv[1],
            .asum => re += @abs(xv[0]) + @abs(xv[1]),
            else => comptime unreachable,
        }
    }
    return switch (r) {
        .dot => .{ .want = .{ re, im }, .scale = mag },
        .nrm2 => .{ .want = .{ @sqrt(re), 0 }, .scale = @sqrt(re) },
        .asum => .{ .want = .{ re, 0 }, .scale = re },
        else => comptime unreachable,
    };
}

/// Part `q` of element `k` after `r` on fresh data: what x and y should
/// hold, then the size of the terms of each.
fn l1Elem(comptime r: L1, comptime p: Prec, b: Bufs(p.F), k: usize, inc: usize, q: usize) [4]f64 {
    const x0 = elem(p, b.x0, k, inc);
    const y0 = elem(p, b.y0, k, inc);
    var wx = x0[q];
    var wy = y0[q];
    var sx: f64 = @abs(x0[q]);
    var sy: f64 = @abs(y0[q]);
    switch (r) {
        .axpy => {
            const al = rounded(p, axpy_alpha);
            const a: [2]f64 = if (p.complex) .{ al, al } else .{ al, 0 };
            const t: [2]f64 = .{ a[0] * x0[0] - a[1] * x0[1], a[0] * x0[1] + a[1] * x0[0] };
            wy = y0[q] + t[q];
            sy = @abs(y0[q]) + @abs(a[0]) * (@abs(x0[0]) + @abs(x0[1])) + @abs(a[1]) * (@abs(x0[0]) + @abs(x0[1]));
        },
        .scal => {
            const a: [2]f64 = if (p.complex) .{ rounded(p, cscal_alpha[0]), rounded(p, cscal_alpha[1]) } else .{ scal_alpha, 0 };
            const t: [2]f64 = .{ a[0] * x0[0] - a[1] * x0[1], a[0] * x0[1] + a[1] * x0[0] };
            wx = t[q];
            sx = @abs(x0[0]) + @abs(x0[1]);
        },
        .copy => wy = x0[q],
        .swap => {
            wx = y0[q];
            wy = x0[q];
        },
        .rot => {
            const rc = rounded(p, rot_cs[0]);
            const rs = rounded(p, rot_cs[1]);
            wx = rc * x0[q] + rs * y0[q];
            wy = rc * y0[q] - rs * x0[q];
            sx = @abs(x0[q]) + @abs(y0[q]);
            sy = sx;
        },
        else => comptime unreachable,
    }
    return .{ wx, wy, sx, sy };
}

/// Whether `r` returns a value rather than writing its vectors.
fn reduces(r: L1) bool {
    return switch (r) {
        .dot, .nrm2, .asum, .iamax => true,
        .axpy, .scal, .copy, .swap, .rot => false,
    };
}

/// A reduction's result, (re, im), widened; `i?amax`'s index.
fn l1Result(comptime r: L1, comptime p: Prec, f: Fn(l1Slot(r, p)), n: usize, inc: usize, b: Bufs(p.F)) [2]f64 {
    const in: isize = @intCast(n);
    const ii: isize = @intCast(inc);
    return switch (r) {
        .dot => if (p.complex) blk: {
            const got = f(in, b.x.ptr, ii, b.y.ptr, ii);
            break :blk .{ got.re, got.im };
        } else .{ f(in, b.x.ptr, ii, b.y.ptr, ii), 0 },
        .nrm2, .asum => .{ f(in, b.x.ptr, ii), 0 },
        .iamax => .{ @floatFromInt(f(in, b.x.ptr, ii)), 0 },
        else => comptime unreachable,
    };
}

/// The largest error found as a multiple of its tolerance (above 1 is
/// wrong), after one call on fresh data. A routine that writes leaves its
/// output in `b`.
fn checkL1(comptime r: L1, comptime p: Prec, f: Fn(l1Slot(r, p)), n: usize, inc: usize, b: Bufs(p.F)) f64 {
    const F = p.F;
    const eps: f64 = std.math.floatEps(F);
    const nf: f64 = @floatFromInt(n);
    b.reset();
    var worst: f64 = 0;
    const span = b.x.len;

    // per element, against a tolerance on the size of its terms
    const Elem = struct {
        fn err(got: f64, want: f64, scale: f64, e: f64) f64 {
            if (got == want) return 0;
            return @abs(got - want) / (8 * e * @max(scale, 1e-300));
        }
    };
    if (comptime reduces(r)) {
        const got = l1Result(r, p, f, n, inc, b);
        if (r == .iamax) {
            worst = if (got[0] == @as(f64, @floatFromInt(firstMax(p, b.x0, n, inc)))) 0 else std.math.inf(f64);
        } else {
            const w = l1Scalar(r, p, b, n, inc);
            const tol = (nf + 8) * eps;
            worst = @max(@abs(got[0] - w.want[0]), @abs(got[1] - w.want[1])) / (tol * @max(w.scale, 1e-300));
        }
        // nothing written
        if (!std.mem.eql(F, b.x, b.x0) or !std.mem.eql(F, b.y, b.y0)) worst = std.math.inf(f64);
    } else {
        invoke(r, p, f, @intCast(n), @intCast(inc), b.x.ptr, b.y.ptr);
        for (0..n) |k| {
            const x1 = elem(p, b.x, k, inc);
            const y1 = elem(p, b.y, k, inc);
            for (0..p.cs()) |q| {
                const w = l1Elem(r, p, b, k, inc, q);
                worst = @max(worst, Elem.err(x1[q], w[0], w[2], eps), Elem.err(y1[q], w[1], w[3], eps));
            }
        }
        // every float not in the vectors as it was
        for (0..span) |i| {
            if (visited(p, i, n, inc)) continue;
            if (b.x[i] != b.x0[i] or b.y[i] != b.y0[i]) worst = std.math.inf(f64);
        }
    }
    if (std.math.isNan(worst)) return std.math.inf(f64);
    return worst;
}

/// The reference's `i?amax`: the first largest `|re| + |im|`, from one.
fn firstMax(comptime p: Prec, x: []const p.F, n: usize, inc: usize) isize {
    var best: p.F = -1;
    var at: usize = 0;
    for (0..n) |k| {
        const i = k * inc * p.cs();
        const v: p.F = if (p.complex) @abs(x[i]) + @abs(x[i + 1]) else @abs(x[i]);
        if (v > best) {
            best = v;
            at = k;
        }
    }
    return @intCast(at + 1);
}

const Case = struct { bytes: usize, inc: usize };

/// Bytes per vector: in L1, L1, L2, memory; then one strided case in L1.
const l1_cases = [_]Case{
    .{ .bytes = 4 << 10, .inc = 1 },
    .{ .bytes = 32 << 10, .inc = 1 },
    .{ .bytes = 1 << 20, .inc = 1 },
    .{ .bytes = 64 << 20, .inc = 1 },
    .{ .bytes = 32 << 10, .inc = 2 },
};

/// Which converted kernel a slot holds: vortex's own, else armv8's, else none.
fn asmFor(comptime name: []const u8) ?struct { f: Fn(name), from: []const u8 } {
    const port = @field(tables.portable, name).?;
    if (@field(tables.vortex, name).? != port) return .{ .f = @field(tables.vortex, name).?, .from = "vortex" };
    if (@field(tables.armv8, name).? != port) return .{ .f = @field(tables.armv8, name).?, .from = "armv8" };
    return null;
}

fn level1(io: Io, gpa: std.mem.Allocator) !void {
    try clock(io, "level 1 start");
    std.debug.print(
        \\
        \\GB/s, best of {d} batches of at least {d} ms; ratios are each over the converted kernel
        \\
        \\| routine | p | n | inc | asm | asm GB/s | portable | vec | vec x4 | portable/asm | vec/asm | vec x4/asm |
        \\| ------- | - | -: | --: | --- | -------: | -------: | --: | -----: | -----------: | ------: | ---------: |
        \\
    , .{ l1_timing.batches, @divTrunc(l1_timing.min_ns, std.time.ns_per_ms) });
    inline for (precs) |p| {
        const len = l1_cases[3].bytes / @sizeOf(p.F);
        const b = try Bufs(p.F).init(gpa, len, 1);
        defer b.deinit(gpa);
        inline for (comptime std.meta.tags(L1)) |r| try raceL1(io, r, p, b);
    }
}

fn raceL1(io: Io, comptime r: L1, comptime p: Prec, b: Bufs(p.F)) !void {
    const name = comptime l1Slot(r, p);
    const C = L1Call(r, p);
    const asm_k = asmFor(name);
    for (l1_cases) |case| {
        const n = case.bytes / p.es() / case.inc;
        var fns: [4]Fn(name) = undefined;
        var labels: [4][]const u8 = undefined;
        var count: usize = 0;
        if (asm_k) |a| {
            fns[count] = a.f;
            labels[count] = "asm";
            count += 1;
        }
        for ([_]Fn(name){ @field(tables.portable, name).?, candidate(r, p, 1), candidate(r, p, 4) }, [_][]const u8{ "portable", "vec", "vec x4" }) |f, l| {
            fns[count] = f;
            labels[count] = l;
            count += 1;
        }
        var ctxs: [4]C = undefined;
        var calls: [4]bench.Call = undefined;
        var ok = true;
        for (0..count) |i| {
            const e = checkL1(r, p, fns[i], n, case.inc, b);
            if (e > 1) {
                std.debug.print("WRONG: {s} {s} {s} n={d} inc={d}: {d:.2} of tolerance\n", .{ @tagName(r), p.name, labels[i], n, case.inc, e });
                ok = false;
            }
            ctxs[i] = .{ .f = fns[i], .n = @intCast(n), .inc = @intCast(case.inc), .x = b.x.ptr, .y = b.y.ptr };
            calls[i] = .{ .ctx = &ctxs[i], .run = C.run };
        }
        if (!ok) continue;
        b.reset();
        var secs: [4]f64 = undefined;
        try bench.race(io, l1_timing, calls[0..count], &secs);
        const bytes = l1Bytes(r, n, p.es());
        var gbs: [4]f64 = undefined;
        for (0..count) |i| gbs[i] = bytes / secs[i] / 1e9;
        const off: usize = if (asm_k != null) 1 else 0;
        std.debug.print("| {s} | {s} | {d} | {d} | {s} | ", .{ @tagName(r), p.name, n, case.inc, if (asm_k) |a| a.from else "none" });
        if (asm_k != null) std.debug.print("{d:.1} | ", .{gbs[0]}) else std.debug.print("- | ", .{});
        for (off..count) |i| std.debug.print("{d:.1} | ", .{gbs[i]});
        for (off..count) |i| {
            if (asm_k != null) std.debug.print("{d:.2} | ", .{gbs[i] / gbs[0]}) else std.debug.print("- | ", .{});
        }
        std.debug.print("\n", .{});
    }
}

// ---- level 2 ---------------------------------------------------------------

const gemv_alpha = 1e-3;

fn GemvCall(comptime p: Prec, comptime name: []const u8) type {
    return struct {
        f: Fn(name),
        n: isize,
        a: [*]p.F,
        x: [*]p.F,
        y: [*]p.F,
        buf: [*]p.F,
        fn call(s: *const @This()) void {
            const F = p.F;
            const al: F = gemv_alpha;
            std.mem.doNotOptimizeAway(if (p.complex)
                s.f(s.n, s.n, 0, al, al, s.a, s.n, s.x, 1, s.y, 1, s.buf)
            else
                s.f(s.n, s.n, 0, al, s.a, s.n, s.x, 1, s.y, 1, s.buf));
        }
        fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// Bytes of `A`: in L1, L2, memory.
const l2_bytes = [_]usize{ 32 << 10, 4 << 20, 64 << 20 };

fn level2(io: Io, gpa: std.mem.Allocator) !void {
    try clock(io, "level 2 start");
    std.debug.print(
        \\
        \\GFLOP/s, square, unit strides, best of {d} batches of at least {d} ms
        \\
        \\| routine | p | n | asm | asm GFLOP/s | portable | vec | vec x4 | portable/asm | vec/asm | vec x4/asm |
        \\| ------- | - | -: | --- | ----------: | -------: | --: | -----: | -----------: | ------: | ---------: |
        \\
    , .{ l1_timing.batches, @divTrunc(l1_timing.min_ns, std.time.ns_per_ms) });
    inline for (precs) |p| {
        inline for (.{ "n", "t" }) |t| try raceGemv(io, gpa, p, t);
    }
}

fn gemvCandidate(comptime p: Prec, comptime t: []const u8, comptime U: comptime_int) Fn(p.name ++ "gemv_" ++ t) {
    const K = if (p.complex) blas.vector.Complex(p.F, U) else blas.vector.Real(p.F, U);
    return &@field(K, "gemv_" ++ t);
}

fn raceGemv(io: Io, gpa: std.mem.Allocator, comptime p: Prec, comptime t: []const u8) !void {
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
    const buf = try gpa.alignedAlloc(F, .@"64", 4096);
    defer gpa.free(buf);
    var rng = std.Random.DefaultPrng.init(2);
    const rand = rng.random();
    for (a) |*v| v.* = @floatCast(rand.float(f64) * 2 - 1);
    for (x) |*v| v.* = @floatCast(rand.float(f64) * 2 - 1);
    for (y0) |*v| v.* = @floatCast(rand.float(f64) * 2 - 1);

    const asm_k = asmFor(name);
    for (l2_bytes) |bytes| {
        const n: usize = std.math.sqrt(bytes / p.es());
        var fns: [4]Fn(name) = undefined;
        var labels: [4][]const u8 = undefined;
        var count: usize = 0;
        if (asm_k) |k| {
            fns[count] = k.f;
            labels[count] = "asm";
            count += 1;
        }
        for ([_]Fn(name){ @field(tables.portable, name).?, gemvCandidate(p, t, 1), gemvCandidate(p, t, 4) }, [_][]const u8{ "portable", "vec", "vec x4" }) |f, l| {
            fns[count] = f;
            labels[count] = l;
            count += 1;
        }
        var ctxs: [4]C = undefined;
        var calls: [4]bench.Call = undefined;
        var ok = true;
        for (0..count) |i| {
            ctxs[i] = .{ .f = fns[i], .n = @intCast(n), .a = a.ptr, .x = x.ptr, .y = y.ptr, .buf = buf.ptr };
            calls[i] = .{ .ctx = &ctxs[i], .run = C.run };
            @memcpy(y, y0);
            ctxs[i].call();
            const e = checkGemv(p, t[0] == 't', n, a, x, y, y0);
            if (e > 1) {
                std.debug.print("WRONG: gemv_{s} {s} {s} n={d}: {d:.2} of tolerance\n", .{ t, p.name, labels[i], n, e });
                ok = false;
            }
        }
        if (!ok) continue;
        @memcpy(y, y0);
        var secs: [4]f64 = undefined;
        try bench.race(io, l1_timing, calls[0..count], &secs);
        const nf: f64 = @floatFromInt(n);
        const flops = nf * nf * @as(f64, if (p.complex) 8 else 2);
        var rate: [4]f64 = undefined;
        for (0..count) |i| rate[i] = flops / secs[i] / 1e9;
        const off: usize = if (asm_k != null) 1 else 0;
        std.debug.print("| gemv_{s} | {s} | {d} | {s} | ", .{ t, p.name, n, if (asm_k) |k| k.from else "none" });
        if (asm_k != null) std.debug.print("{d:.1} | ", .{rate[0]}) else std.debug.print("- | ", .{});
        for (off..count) |i| std.debug.print("{d:.1} | ", .{rate[i]});
        for (off..count) |i| {
            if (asm_k != null) std.debug.print("{d:.2} | ", .{rate[i] / rate[0]}) else std.debug.print("- | ", .{});
        }
        std.debug.print("\n", .{});
    }
}

/// `y0 + alpha*op(A)*x` against `y`, as in `checkL1`; `y` past `n` unchanged.
fn checkGemv(comptime p: Prec, trans: bool, n: usize, a: []const p.F, x: []const p.F, y: []const p.F, y0: []const p.F) f64 {
    const eps: f64 = std.math.floatEps(p.F);
    const cs = p.cs();
    const g = rounded(p, gemv_alpha);
    const al = [2]f64{ g, if (p.complex) g else 0 };
    var worst: f64 = 0;
    for (0..n) |i| {
        var re: f64 = 0;
        var im: f64 = 0;
        var mag: f64 = 0;
        for (0..n) |l| {
            // A(i, l) or A(l, i)
            const at = if (trans) (l + i * n) * cs else (i + l * n) * cs;
            const av = [2]f64{ a[at], if (p.complex) a[at + 1] else 0 };
            const xv = [2]f64{ x[l * cs], if (p.complex) x[l * cs + 1] else 0 };
            re += av[0] * xv[0] - av[1] * xv[1];
            im += av[0] * xv[1] + av[1] * xv[0];
            mag += (@abs(av[0]) + @abs(av[1])) * (@abs(xv[0]) + @abs(xv[1]));
        }
        const y_re: f64 = y0[i * cs];
        const y_im: f64 = if (p.complex) y0[i * cs + 1] else 0;
        const want = [2]f64{ y_re + al[0] * re - al[1] * im, y_im + al[0] * im + al[1] * re };
        for (0..cs) |q| {
            const got: f64 = y[i * cs + q];
            const scale = @abs(y0[i * cs + q]) + (@abs(al[0]) + @abs(al[1])) * mag;
            worst = @max(worst, @abs(got - want[q]) / ((@as(f64, @floatFromInt(n)) + 8) * eps * scale));
        }
    }
    for (n * cs..y.len) |i| {
        if (y[i] != y0[i]) worst = std.math.inf(f64);
    }
    if (std.math.isNan(worst)) return std.math.inf(f64);
    return worst;
}

// ---- level 3 ---------------------------------------------------------------

extern fn ob_sgemm_kernel(isize, isize, isize, f32, [*]f32, [*]f32, [*]f32, isize) callconv(.c) i32;
extern fn ob_sgemm_itcopy(isize, isize, [*]f32, isize, [*]f32) callconv(.c) i32;
extern fn ob_sgemm_oncopy(isize, isize, [*]f32, isize, [*]f32) callconv(.c) i32;
extern fn ob_dgemm_kernel(isize, isize, isize, f64, [*]f64, [*]f64, [*]f64, isize) callconv(.c) i32;
extern fn ob_dgemm_itcopy(isize, isize, [*]f64, isize, [*]f64) callconv(.c) i32;
extern fn ob_dgemm_oncopy(isize, isize, [*]f64, isize, [*]f64) callconv(.c) i32;

/// A gemm kernel with the copies that pack for it, for `NN` (`level3.c`
/// packs `A` with `itcopy` and `B` with `oncopy` there).
fn GemmImpl(comptime p: Prec) type {
    return struct {
        name: []const u8,
        kernel: Fn(p.name ++ "gemm_kernel"),
        itcopy: Fn(p.name ++ "gemm_itcopy"),
        oncopy: Fn(p.name ++ "gemm_oncopy"),
        um: usize,
        un: usize,

        fn of(name: []const u8, vt: *const Vtable) @This() {
            const pre = p.name ++ "gemm_";
            return .{
                .name = name,
                .kernel = @field(vt, pre ++ "kernel").?,
                .itcopy = @field(vt, pre ++ "itcopy").?,
                .oncopy = @field(vt, pre ++ "oncopy").?,
                .um = @intCast(@field(vt, pre ++ "unroll_m")),
                .un = @intCast(@field(vt, pre ++ "unroll_n")),
            };
        }
    };
}

/// OpenBLAS's vortex kernels: `sgemm_kernel_16x4.S` with `sgemm_tcopy_16.S`
/// and `sgemm_ncopy_4.S`; `dgemm_kernel_8x4.S` with `dgemm_tcopy_8.S` and
/// `dgemm_ncopy_4.S` (`KERNEL.NEOVERSEN1`, which `KERNEL.VORTEX` includes).
fn asmGemm(comptime p: Prec) GemmImpl(p) {
    return if (p.F == f32)
        .{ .name = "asm 16x4", .kernel = &ob_sgemm_kernel, .itcopy = &ob_sgemm_itcopy, .oncopy = &ob_sgemm_oncopy, .um = 16, .un = 4 }
    else
        .{ .name = "asm 8x4", .kernel = &ob_dgemm_kernel, .itcopy = &ob_dgemm_itcopy, .oncopy = &ob_dgemm_oncopy, .um = 8, .un = 4 };
}

/// `level3.c`'s blocking, from the vortex table, with `P` a whole number
/// of the kernel's row panels.
const Blocking = struct {
    p: usize,
    q: usize,
    r: usize,
    um: usize,
    un: usize,

    fn of(comptime p: Prec, im: GemmImpl(p)) Blocking {
        const pre = p.name ++ "gemm_";
        const bp: usize = @intCast(@field(tables.vortex, pre ++ "p"));
        return .{
            .p = @max(im.um, bp / im.um * im.um),
            .q = @intCast(@field(tables.vortex, pre ++ "q")),
            .r = @intCast(@field(tables.vortex, pre ++ "r")),
            .um = im.um,
            .un = im.un,
        };
    }

    fn split(rest: usize, max: usize, unroll: usize) usize {
        if (rest >= 2 * max) return max;
        if (rest > max) return @min(rest, (rest / 2 + unroll - 1) / unroll * unroll);
        return rest;
    }

    fn chunk(b: Blocking, rest: usize) usize {
        if (rest >= 3 * b.un) return 3 * b.un;
        if (rest > b.un) return b.un;
        return rest;
    }
};

/// `C += A*B`, square `n`, through `level3.c`'s loops as `blas3/driver.zig`
/// has them, with `im`'s kernel and copies.
fn GemmCall(comptime p: Prec) type {
    const F = p.F;
    return struct {
        im: GemmImpl(p),
        blk: Blocking,
        n: usize,
        a: [*]F,
        b: [*]F,
        c: [*]F,
        sa: [*]F,
        sb: [*]F,

        fn call(s: *const @This()) void {
            const n = s.n;
            const blk = s.blk;
            const L = isize;
            var js: usize = 0;
            while (js < n) : (js += blk.r) {
                const min_j = @min(blk.r, n - js);
                var ls: usize = 0;
                var min_l: usize = 0;
                while (ls < n) : (ls += min_l) {
                    min_l = Blocking.split(n - ls, blk.q, blk.um);
                    var min_i = Blocking.split(n, blk.p, blk.um);
                    const once = min_i == n;
                    _ = s.im.itcopy(@intCast(min_l), @intCast(min_i), s.a + ls * n, @intCast(n), s.sa);
                    var jjs = js;
                    while (jjs < js + min_j) {
                        const min_jj = blk.chunk(js + min_j - jjs);
                        const pb = if (once) s.sb else s.sb + min_l * (jjs - js);
                        _ = s.im.oncopy(@intCast(min_l), @intCast(min_jj), s.b + ls + jjs * n, @intCast(n), pb);
                        _ = s.im.kernel(@intCast(min_i), @intCast(min_jj), @intCast(min_l), 1, s.sa, pb, s.c + jjs * n, @as(L, @intCast(n)));
                        jjs += min_jj;
                    }
                    var is = min_i;
                    while (is < n) : (is += min_i) {
                        min_i = Blocking.split(n - is, blk.p, blk.um);
                        _ = s.im.itcopy(@intCast(min_l), @intCast(min_i), s.a + is + ls * n, @intCast(n), s.sa);
                        jjs = js;
                        while (jjs < js + min_j) {
                            const min_jj = blk.chunk(js + min_j - jjs);
                            _ = s.im.kernel(@intCast(min_i), @intCast(min_jj), @intCast(min_l), 1, s.sa, s.sb + min_l * (jjs - js), s.c + is + jjs * n, @as(L, @intCast(n)));
                            jjs += min_jj;
                        }
                    }
                }
            }
        }

        fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// The kernel alone on panels its own copies packed: `m` by `n`, depth `k`.
fn KernelCall(comptime p: Prec) type {
    const F = p.F;
    return struct {
        im: GemmImpl(p),
        m: usize,
        n: usize,
        k: usize,
        sa: [*]F,
        sb: [*]F,
        c: [*]F,

        fn call(s: *const @This()) void {
            _ = s.im.kernel(@intCast(s.m), @intCast(s.n), @intCast(s.k), 1, s.sa, s.sb, s.c, @intCast(s.m));
        }

        fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// `c0 + A*B` against `c` at `samples` entries (all when 0), as `checkL1`.
fn checkGemm(comptime F: type, m: usize, n: usize, k: usize, a: []const F, b: []const F, c: []const F, c0: []const F, samples: usize) f64 {
    const eps: f64 = std.math.floatEps(F);
    var rng = std.Random.DefaultPrng.init(3);
    const rand = rng.random();
    var worst: f64 = 0;
    const count = if (samples == 0) m * n else samples;
    for (0..count) |s| {
        const i = if (samples == 0) s % m else rand.uintLessThan(usize, m);
        const j = if (samples == 0) s / m else rand.uintLessThan(usize, n);
        var want: f64 = c0[i + j * m];
        var mag: f64 = @abs(want);
        for (0..k) |l| {
            want += @as(f64, a[i + l * m]) * b[l + j * k];
            mag += @abs(@as(f64, a[i + l * m]) * b[l + j * k]);
        }
        worst = @max(worst, @abs(c[i + j * m] - want) / ((@as(f64, @floatFromInt(k)) + 8) * eps * mag));
    }
    if (std.math.isNan(worst)) return std.math.inf(f64);
    return worst;
}

const gemm_sizes = [_]usize{ 64, 128, 256, 512, 1024, 2048 };

fn level3(io: Io, gpa: std.mem.Allocator) !void {
    try clock(io, "level 3 start");
    inline for (precs[0..2]) |p| try raceGemm(io, gpa, p);
    std.debug.print("\ncgemm, zgemm: no Zig kernel (only real gemm is ported), not raced\n", .{});
}

fn raceGemm(io: Io, gpa: std.mem.Allocator, comptime p: Prec) !void {
    const F = p.F;
    const impls = [_]GemmImpl(p){
        asmGemm(p),
        GemmImpl(p).of("zig vortex", &tables.vortex),
        GemmImpl(p).of("zig vortexm4", &tables.vortexm4),
        GemmImpl(p).of("zig portable", &tables.portable),
    };
    const max_n = gemm_sizes[gemm_sizes.len - 1];
    const a = try gpa.alignedAlloc(F, .@"64", max_n * max_n);
    defer gpa.free(a);
    const b = try gpa.alignedAlloc(F, .@"64", max_n * max_n);
    defer gpa.free(b);
    const c = try gpa.alignedAlloc(F, .@"64", max_n * max_n);
    defer gpa.free(c);
    const c0 = try gpa.alloc(F, max_n * max_n);
    defer gpa.free(c0);
    var rng = std.Random.DefaultPrng.init(4);
    const rand = rng.random();
    for (a, b, c0) |*x, *y, *z| {
        x.* = @floatCast(rand.float(f64) * 2 - 1);
        y.* = @floatCast(rand.float(f64) * 2 - 1);
        z.* = @floatCast(rand.float(f64) * 2 - 1);
    }

    var blks: [impls.len]Blocking = undefined;
    var sa_len: usize = 0;
    var sb_len: usize = 0;
    for (impls, &blks) |im, *blk| {
        blk.* = Blocking.of(p, im);
        sa_len = @max(sa_len, blk.p * (blk.q + blk.um));
        sb_len = @max(sb_len, (blk.q + blk.um) * blk.r);
    }
    var sas: [impls.len][]align(64) F = undefined;
    var sbs: [impls.len][]align(64) F = undefined;
    for (0..impls.len) |i| {
        sas[i] = try gpa.alignedAlloc(F, .@"64", sa_len);
        sbs[i] = try gpa.alignedAlloc(F, .@"64", sb_len);
    }
    defer {
        for (0..impls.len) |i| {
            gpa.free(sas[i]);
            gpa.free(sbs[i]);
        }
    }

    // the kernel alone, on panels in L1 and L2
    const km = 64;
    const kn = 24;
    const kk = 256;
    std.debug.print(
        \\
        \\{s}gemm kernel alone, {d} x {d}, depth {d}, its own copies' panels: GFLOP/s, best of {d} batches of at least {d} ms
        \\
        \\| kernel | tile | GFLOP/s | over asm |
        \\| ------ | ---- | ------: | -------: |
        \\
    , .{ p.name, km, kn, kk, l1_timing.batches, @divTrunc(l1_timing.min_ns, std.time.ns_per_ms) });
    {
        const KC = KernelCall(p);
        var ctxs: [impls.len]KC = undefined;
        var calls: [impls.len]bench.Call = undefined;
        for (impls, 0..) |im, i| {
            _ = im.itcopy(kk, km, a.ptr, km, sas[i].ptr);
            _ = im.oncopy(kk, kn, b.ptr, kk, sbs[i].ptr);
            ctxs[i] = .{ .im = im, .m = km, .n = kn, .k = kk, .sa = sas[i].ptr, .sb = sbs[i].ptr, .c = c.ptr };
            calls[i] = .{ .ctx = &ctxs[i], .run = KC.run };
            @memcpy(c[0 .. km * kn], c0[0 .. km * kn]);
            ctxs[i].call();
            const e = checkGemm(F, km, kn, kk, a[0 .. km * kk], b[0 .. kk * kn], c[0 .. km * kn], c0[0 .. km * kn], 0);
            if (e > 1) {
                std.debug.print("WRONG: {s}gemm kernel {s}: {d:.2} of tolerance\n", .{ p.name, im.name, e });
                return error.WrongKernel;
            }
        }
        var secs: [impls.len]f64 = undefined;
        try bench.race(io, l1_timing, &calls, &secs);
        const flops: f64 = 2.0 * km * kn * kk;
        for (impls, 0..) |im, i| std.debug.print("| {s} | {d} x {d} | {d:.1} | {d:.2} |\n", .{ im.name, im.um, im.un, flops / secs[i] / 1e9, secs[0] / secs[i] });
    }

    std.debug.print(
        \\
        \\{s}gemm NN, square, `C += A*B` through `level3.c`'s loops (vortex P {d}, Q {d}, R {d}): GFLOP/s, best of {d} batches of at least {d} ms
        \\
        \\| n |
    , .{ p.name, blks[0].p, blks[0].q, blks[0].r, l3_timing.batches, @divTrunc(l3_timing.min_ns, std.time.ns_per_ms) });
    for (impls) |im| std.debug.print(" {s} |", .{im.name});
    for (impls[1..]) |im| std.debug.print(" {s}/asm |", .{im.name});
    std.debug.print("\n| -: |", .{});
    for (0..2 * impls.len - 1) |_| std.debug.print(" -: |", .{});
    std.debug.print("\n", .{});

    const GC = GemmCall(p);
    for (gemm_sizes) |n| {
        var ctxs: [impls.len]GC = undefined;
        var calls: [impls.len]bench.Call = undefined;
        for (impls, 0..) |im, i| {
            ctxs[i] = .{ .im = im, .blk = blks[i], .n = n, .a = a.ptr, .b = b.ptr, .c = c.ptr, .sa = sas[i].ptr, .sb = sbs[i].ptr };
            calls[i] = .{ .ctx = &ctxs[i], .run = GC.run };
            @memcpy(c[0 .. n * n], c0[0 .. n * n]);
            ctxs[i].call();
            const e = checkGemm(F, n, n, n, a[0 .. n * n], b[0 .. n * n], c[0 .. n * n], c0[0 .. n * n], 256);
            if (e > 1) {
                std.debug.print("WRONG: {s}gemm {s} n={d}: {d:.2} of tolerance\n", .{ p.name, im.name, n, e });
                return error.WrongGemm;
            }
        }
        var secs: [impls.len]f64 = undefined;
        try bench.race(io, l3_timing, &calls, &secs);
        const nf: f64 = @floatFromInt(n);
        std.debug.print("| {d} |", .{n});
        for (0..impls.len) |i| std.debug.print(" {d:.1} |", .{2 * nf * nf * nf / secs[i] / 1e9});
        for (1..impls.len) |i| std.debug.print(" {d:.2} |", .{secs[0] / secs[i]});
        std.debug.print("\n", .{});
    }
}

// ---- rstd's own gemm --------------------------------------------------------

/// `kern3`: the selected core's gemm kernel alone, on a `P` by 64 block of
/// depth `Q` its own copies packed, against rstd's whole `gemm` on one
/// worker: how much of the kernel's rate the driver keeps (T21).
fn ownLevel3(io: Io, gpa: std.mem.Allocator) !void {
    try clock(io, "own level 3 start");
    inline for (precs) |p| try ownGemm(io, gpa, p);
}

/// The selected core's kernel alone on panels its own copies packed: `m` by
/// `n`, depth `k`.
fn OwnKernel(comptime p: Prec) type {
    const F = p.F;
    return struct {
        kernel: Fn(p.name ++ if (p.complex) "gemm_kernel_n" else "gemm_kernel"),
        m: usize,
        n: usize,
        k: usize,
        sa: [*]F,
        sb: [*]F,
        c: [*]F,

        fn call(s: *const @This()) void {
            if (p.complex)
                _ = s.kernel(@intCast(s.m), @intCast(s.n), @intCast(s.k), 1, 0, s.sa, s.sb, s.c, @intCast(s.m))
            else
                _ = s.kernel(@intCast(s.m), @intCast(s.n), @intCast(s.k), 1, s.sa, s.sb, s.c, @intCast(s.m));
        }

        fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// rstd's public `gemm`, `C += A*B`, square `n`, on `Io.failing` (one worker).
fn OwnGemm(comptime p: Prec) type {
    const F = p.F;
    return struct {
        n: usize,
        a: [*]F,
        b: [*]F,
        c: [*]F,

        fn call(s: *const @This()) void {
            const ni: isize = @intCast(s.n);
            const one: Num(p) = if (p.complex) .init(1, 0) else 1;
            blas.gemm(rstd_alloc, .failing, Num(p), .n, .n, ni, ni, ni, one, cast(p, s.a), ni, cast(p, s.b), ni, one, cast(p, s.c), ni) catch |e| fail(e);
        }

        fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// `checkGemm` for real or interleaved complex `p`.
fn checkOwn(comptime p: Prec, m: usize, n: usize, k: usize, a: []const p.F, b: []const p.F, c: []const p.F, c0: []const p.F, samples: usize) f64 {
    if (!p.complex) return checkGemm(p.F, m, n, k, a, b, c, c0, samples);
    const eps: f64 = std.math.floatEps(p.F);
    var rng = std.Random.DefaultPrng.init(3);
    const rand = rng.random();
    var worst: f64 = 0;
    for (0..samples) |_| {
        const i = rand.uintLessThan(usize, m);
        const j = rand.uintLessThan(usize, n);
        const o = 2 * (i + j * m);
        var re: f64 = c0[o];
        var im: f64 = c0[o + 1];
        var mag: f64 = @abs(re) + @abs(im);
        for (0..k) |l| {
            const ar: f64 = a[2 * (i + l * m)];
            const ai: f64 = a[2 * (i + l * m) + 1];
            const br: f64 = b[2 * (l + j * k)];
            const bi: f64 = b[2 * (l + j * k) + 1];
            re += ar * br - ai * bi;
            im += ar * bi + ai * br;
            mag += (@abs(ar) + @abs(ai)) * (@abs(br) + @abs(bi));
        }
        const err = @max(@abs(c[o] - re), @abs(c[o + 1] - im));
        worst = @max(worst, err / ((@as(f64, @floatFromInt(k)) + 8) * eps * mag));
    }
    if (std.math.isNan(worst)) return std.math.inf(f64);
    return worst;
}

const own_sizes = [_]usize{ 128, 256, 512, 1024 };

fn ownGemm(io: Io, gpa: std.mem.Allocator, comptime p: Prec) !void {
    const F = p.F;
    const cs = p.cs();
    const vt = blas.vtable();
    const pre = p.name ++ "gemm_";
    const um: usize = @intCast(@field(vt, pre ++ "unroll_m"));
    const un: usize = @intCast(@field(vt, pre ++ "unroll_n"));
    const bp: usize = @intCast(@field(vt, pre ++ "p"));
    // `P` a whole number of row panels, as `Blocking.of`
    const km = @max(um, bp / um * um);
    const kk: usize = @intCast(@field(vt, pre ++ "q"));
    const kn: usize = 64;
    const max_n = own_sizes[own_sizes.len - 1];
    const len = cs * @max(max_n * max_n, km * kk, kk * kn);
    const a = try gpa.alignedAlloc(F, .@"64", len);
    defer gpa.free(a);
    const b = try gpa.alignedAlloc(F, .@"64", len);
    defer gpa.free(b);
    const c = try gpa.alignedAlloc(F, .@"64", len);
    defer gpa.free(c);
    const c0 = try gpa.alloc(F, len);
    defer gpa.free(c0);
    const sa = try gpa.alignedAlloc(F, .@"64", cs * (km + um) * (kk + um));
    defer gpa.free(sa);
    const sb = try gpa.alignedAlloc(F, .@"64", cs * (kk + un) * (kn + un));
    defer gpa.free(sb);
    var rng = std.Random.DefaultPrng.init(4);
    const rand = rng.random();
    for (a, b, c0) |*x, *y, *z| {
        x.* = @floatCast(rand.float(f64) * 2 - 1);
        y.* = @floatCast(rand.float(f64) * 2 - 1);
        z.* = @floatCast(rand.float(f64) * 2 - 1);
    }
    const flop: f64 = if (p.complex) 8 else 2;

    std.debug.print(
        \\
        \\{s}gemm on rstd's {s} ({d} x {d} kernel, P {d}, Q {d}), one worker: GFLOP/s, best of {d} batches of at least {d} ms (kernel {d} ms)
        \\
        \\| call | GFLOP/s | of the kernel alone |
        \\| ---- | ------: | ------------------: |
        \\
    , .{ p.name, blas.coreName(), um, un, km, kk, l3_timing.batches, @divTrunc(l3_timing.min_ns, std.time.ns_per_ms), @divTrunc(l1_timing.min_ns, std.time.ns_per_ms) });

    _ = @field(vt, pre ++ "itcopy").?(@intCast(kk), @intCast(km), a.ptr, @intCast(km), sa.ptr);
    _ = @field(vt, pre ++ "oncopy").?(@intCast(kk), @intCast(kn), b.ptr, @intCast(kk), sb.ptr);
    var kc: OwnKernel(p) = .{ .kernel = @field(vt, if (p.complex) pre ++ "kernel_n" else pre ++ "kernel").?, .m = km, .n = kn, .k = kk, .sa = sa.ptr, .sb = sb.ptr, .c = c.ptr };
    @memcpy(c[0 .. cs * km * kn], c0[0 .. cs * km * kn]);
    kc.call();
    const ke = checkOwn(p, km, kn, kk, a[0 .. cs * km * kk], b[0 .. cs * kk * kn], c[0 .. cs * km * kn], c0[0 .. cs * km * kn], 256);
    if (ke > 1) {
        std.debug.print("WRONG: {s}gemm kernel: {d:.2} of tolerance\n", .{ p.name, ke });
        return error.WrongKernel;
    }
    var ksecs: [1]f64 = undefined;
    try bench.race(io, l1_timing, &.{.{ .ctx = &kc, .run = OwnKernel(p).run }}, &ksecs);
    const krate = flop * @as(f64, @floatFromInt(km * kn * kk)) / ksecs[0] / 1e9;
    std.debug.print("| kernel alone, {d} x {d}, depth {d} | {d:.1} | 1.00 |\n", .{ km, kn, kk, krate });

    for (own_sizes) |n| {
        var gc: OwnGemm(p) = .{ .n = n, .a = a.ptr, .b = b.ptr, .c = c.ptr };
        @memcpy(c[0 .. cs * n * n], c0[0 .. cs * n * n]);
        gc.call();
        const e = checkOwn(p, n, n, n, a[0 .. cs * n * n], b[0 .. cs * n * n], c[0 .. cs * n * n], c0[0 .. cs * n * n], 256);
        if (e > 1) {
            std.debug.print("WRONG: {s}gemm n={d}: {d:.2} of tolerance\n", .{ p.name, n, e });
            return error.WrongGemm;
        }
        var secs: [1]f64 = undefined;
        try bench.race(io, l3_timing, &.{.{ .ctx = &gc, .run = OwnGemm(p).run }}, &secs);
        const nf: f64 = @floatFromInt(n);
        const rate = flop * nf * nf * nf / secs[0] / 1e9;
        std.debug.print("| gemm, n = {d} | {d:.1} | {d:.2} |\n", .{ n, rate, rate / krate });
    }
}

// ---- libraries -------------------------------------------------------------

const LibOptions = struct {
    /// Each library's path by `Src` less one: OpenBLAS, BLIS, Accelerate,
    /// reference LAPACK, the Eigen shim, AOCL-libFLAME.
    paths: [lib_handles.len]?[]const u8 = .{ null, null, if (os_tag.isDarwin()) "/System/Library/Frameworks/Accelerate.framework/Accelerate" else null, null, null, null },
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
const Src = enum { rstd, openblas, blis, accelerate, reference, eigen, aoclflame };

/// Each library's handle, by `Src` less one; null when not given.
var lib_handles: [std.meta.tags(Src).len - 1]?Dll = @splat(null);

fn loaded(src: Src) bool {
    if (src == .rstd) return !alone;
    return lib_handles[@backingInt(src) - 1] != null;
}

/// The first error an rstd routine returned inside a kernel-ABI wrapper.
var lib_err: ?anyerror = null;

fn fail(e: anyerror) void {
    if (lib_err == null) lib_err = e;
}

/// rstd's routines take an allocator; a C program would hand them malloc's.
const rstd_alloc = std.heap.c_allocator;
/// The `Io` rstd's routines thread on.
const RstdIo = enum {
    /// The bench's `Io.Threaded`, limited to `--threads` less the caller.
    threaded,
    /// Serial: `Io.failing` runs every task inline.
    failing,
};

/// And an `Io`: `main` sets it.
var rstd_io: Io = .failing;

/// Which `Io` rstd ran on, for each library table's notes.
var io_note: []const u8 = "";
var io_note_buf: [96]u8 = undefined;

/// `--max-n`: rows of a larger `n` are skipped.
var n_max: usize = std.math.maxInt(usize);

/// The libraries' thread counts by environment, read when each loads, and
/// rstd's ceiling (`OPENBLAS_NUM_THREADS`), read when it selects its core.
fn setThreadEnv(threads: ?u31) !void {
    const t = threads orelse return;
    var buf: [16]u8 = undefined;
    const v = try std.mem.printSentinel(&buf, "{d}", .{t}, 0);
    for ([_][*:0]const u8{ "OPENBLAS_NUM_THREADS", "BLIS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS" }) |name| setEnv(name, v);
}

extern "c" fn setenv(name: [*:0]const u8, value: [*:0]const u8, overwrite: c_int) c_int;
extern "c" fn unsetenv(name: [*:0]const u8) c_int;
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
    const h = &(lib_handles[@backingInt(src) - 1] orelse return null);
    if (src == .accelerate) {
        if (h.lookup(T, name ++ "$NEWLAPACK")) |f| return f;
    }
    return h.lookup(T, name);
}

/// Fortran-ABI routine `base` (`dpotrf`) from `src`'s library: Accelerate's
/// new LP64 interface first (`dpotrf$NEWLAPACK`, or with the underscore:
/// `apple_accel.lp64Names`), else `base_`.
fn lookupF(comptime src: Src, comptime T: type, comptime base: []const u8) ?T {
    const h = &(lib_handles[@backingInt(src) - 1] orelse return null);
    if (src == .accelerate) {
        if (h.lookup(T, std.fmt.comptimePrint("{s}$NEWLAPACK", .{base}))) |f| return f;
        if (h.lookup(T, std.fmt.comptimePrint("{s}_$NEWLAPACK", .{base}))) |f| return f;
    }
    return h.lookup(T, std.fmt.comptimePrint("{s}_", .{base}));
}

/// Opens the libraries given, at `--threads` each (`setThreadEnv`), and
/// prints how each is configured, and rstd. OpenBLAS is required unless
/// `alone`; the reference binaries are skipped when they cannot load.
fn loadLibs(opts: LibOptions) !void {
    var buf: [64]u8 = undefined;
    if (opts.ob_core) |c| setEnv("OPENBLAS_CORETYPE", try std.mem.printSentinel(&buf, "{s}", .{c}, 0));
    const paths = opts.paths;
    if (paths[0] == null and !alone) return error.NeedOpenblasPath;
    for (paths, &lib_handles, 0..) |path, *h, i| {
        const p = path orelse continue;
        // in `alone`, only what Accelerate is raced against
        if (alone and i > 2) continue;
        h.* = Dll.open(p) catch |e| {
            std.debug.print("cannot load {s}: {t} {s}\n", .{ p, e, dlError() });
            if (i < 3) return e;
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
            paths[0].?,
            lookup(.openblas, Str, "openblas_get_corename").?(),
            config,
            lookup(.openblas, Int, "openblas_get_num_threads").?(),
            lapackVersion(.openblas, &buf),
        });
    }
    if (loaded(.reference)) std.debug.print("reference: {s}, lapack {s}, one thread\n", .{ paths[3].?, lapackVersion(.reference, &buf) });
    if (loaded(.eigen)) {
        const v = if (lookup(.eigen, Str, "rstd_bench_eigen_version")) |f| std.mem.span(f()) else "?";
        std.debug.print("eigen: {s}, version {s}, one thread\n", .{ paths[4].?, v });
    }
    if (loaded(.aoclflame)) std.debug.print("aocl-libflame: {s}, lapack {s}, threads by OMP_NUM_THREADS\n", .{ paths[5].?, lapackVersion(.aoclflame, &buf) });

    if (loaded(.blis)) {
        // `dim_t` is 64 bits, `arch_t` an enum
        const Dim = *const fn () callconv(.c) i64;
        if (opts.threads) |t| lookup(.blis, *const fn (i64) callconv(.c) void, "bli_thread_set_num_threads").?(t);
        const arch = lookup(.blis, *const fn () callconv(.c) c_int, "bli_arch_query_id").?();
        std.debug.print("blis: {s}\nblis version: {s}, config: {s}, threads: {d}, cblas: {s}\n", .{
            paths[1].?,
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
        std.debug.print("accelerate: {s}, BLASGetThreading: {s}\n", .{ paths[2].?, accelThreading() });
    }
    if (alone) {
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
fn Num(comptime p: Prec) type {
    return if (p.complex) std.math.Complex(p.F) else p.F;
}

fn cast(comptime p: Prec, x: [*]p.F) [*]Num(p) {
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
const Col = struct {
    label: []const u8,
    short: []const u8,
    lib: ?Src = null,
    /// Whether this build has it, for a column that is no library.
    on: bool = true,
    /// Its rows of `n` up to this only (a slow baseline).
    max_n: usize = std.math.maxInt(usize),

    /// rstd's own columns (rstd, `asm`, `vec x4`, `ob .S`) are out in `alone`.
    fn shown(c: Col) bool {
        return if (c.lib) |l| loaded(l) else c.on and !alone;
    }
};

/// `lib1` and `lib2`'s contestants. Ratios are over the first, OpenBLAS.
const lib_cols = [_]Col{
    .{ .label = "openblas", .short = "ob", .lib = .openblas },
    .{ .label = "rstd", .short = "rstd" },
    .{ .label = "blis", .short = "blis", .lib = .blis },
    .{ .label = "accelerate", .short = "accel", .lib = .accelerate },
    .{ .label = "asm", .short = "asm", .on = arm },
    .{ .label = "vec x4", .short = "x4", .on = arm },
};
/// `lib3`'s: `ob .S` is OpenBLAS's vortex `.S` through the bench's loops.
const lib3_cols = lib_cols[0..4].* ++ [_]Col{.{ .label = "ob .S", .short = "ob .S", .on = ob_asm }};
const rstd_col = 1;

/// Column `o` of a race in the order checked: rstd's first, whose output
/// the libraries' are compared with.
fn rstdFirst(o: usize) usize {
    return switch (o) {
        0 => rstd_col,
        rstd_col => 0,
        else => o,
    };
}

/// rstd's largest difference from a library's result, in eps of the type
/// relative to the size of the terms, and what accumulation order explains.
const Diff = struct { eps: f64, bound: f64 };

/// What accumulation order explains for a sum of `n` terms: a few ulps
/// times `n`.
fn sumBound(n: usize) f64 {
    return 2 * (@as(f64, @floatFromInt(n)) + 8);
}

/// A table's header: `lead` and `rule`, then each shown contestant's rate,
/// its ratio over OpenBLAS, rstd's difference from each library, rstd's
/// error, and a hash of rstd's output.
fn printHeader(comptime lead: []const u8, comptime rule: []const u8, cols: []const Col) void {
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
    if (alone) {
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
const merge_labels = blk: {
    var l: []const merge.Label = &.{};
    for (lib_cols ++ ext_cols[4..] ++ lib3_cols[4..]) |c| l = l ++ [_]merge.Label{.{ .label = c.label, .short = c.short }};
    break :blk l[0..l.len].*;
};

/// Marks the start of a section's output for `merge.zig`, which reruns
/// noisy sections by name.
fn sectionMark(comptime name: []const u8) void {
    std.debug.print("\n<!-- section " ++ name ++ " -->\n", .{});
}

/// A row after its lead: `printHeader`'s columns; `-` where one is missing,
/// `!` on a difference past its bound.
fn printRow(cols: []const Col, rates: []const ?f64, diffs: []const ?Diff, mine: Rstd) void {
    for (cols, rates) |c, r| {
        if (!c.shown()) continue;
        if (r) |v| std.debug.print("{d:.1} | ", .{v}) else std.debug.print("- | ", .{});
    }
    for (cols[1..], rates[1..]) |c, r| {
        if (!c.shown()) continue;
        if (r != null and rates[0] != null) std.debug.print("{d:.2} | ", .{r.? / rates[0].?}) else std.debug.print("- | ", .{});
    }
    if (alone) return std.debug.print("\n", .{});
    for (cols, diffs) |c, d| {
        if (c.lib == null or !c.shown()) continue;
        if (d) |v| std.debug.print("{d:.1}{s} | ", .{ v.eps, if (v.eps > v.bound) " !" else "" }) else std.debug.print("- | ", .{});
    }
    if (mine.err) |e| std.debug.print("{d:.1} | ", .{e}) else std.debug.print("- | ", .{});
    if (mine.bits) |h| std.debug.print("{x:0>8} |\n", .{@as(u32, @truncate(h))}) else std.debug.print("- |\n", .{});
}

/// rstd's error from the f64 reference, in eps of the type relative to the
/// size of the terms, and a hash of its output.
const Rstd = struct { err: ?f64 = null, bits: ?u64 = null };

/// What rstd's own table holds in slot `name`: converted `asm`, or
/// `portable`.
fn rstdRuns(comptime name: []const u8) []const u8 {
    if (alone) return "-";
    return if (@field(blas.vtable().*, name).? == @field(tables.portable, name).?) "portable" else "asm";
}

const diff_note = "`rstd-<lib>` is rstd's largest difference from that library's result, in eps of the type relative to the size of the terms; `!` marks one past what accumulation order explains: 0 for `copy`, `swap` and `i?amax`, 8 for the other element-wise routines, 2 (n + 8) for a sum of n terms. `rstd err` is rstd's own error from an f64 reference in the same units (gemm: at 256 entries), which tells whose result a difference comes from. `rstd bits` hashes rstd's output, the same on machines that run the same kernels";

fn libLevel1(io: Io, gpa: std.mem.Allocator) !void {
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
    const asm_k = if (arm and !alone) asmFor(name) else null;
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
        var ctxs: [lib_cols.len]C = undefined;
        var calls: [lib_cols.len]bench.Call = undefined;
        var at: [lib_cols.len]?usize = @splat(null);
        var diffs: [lib_cols.len]?Diff = @splat(null);
        var count: usize = 0;
        var ok = true;
        var rs: [2]f64 = .{ 0, 0 };
        var mine: Rstd = .{};
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
                mine = .{ .err = e * @as(f64, if (comptime reduces(r)) @floatFromInt(n + 8) else 8), .bits = h.final() };
            } else if (lib_cols[i].lib != null) diffs[i] = diffL1(r, p, n, case.inc, b, kx, ky, rs, res);
            ctxs[count] = .{ .f = f, .n = @intCast(n), .inc = @intCast(case.inc), .x = b.x.ptr, .y = b.y.ptr };
            at[i] = count;
            count += 1;
        }
        if (!ok) continue;
        for (0..count) |i| calls[i] = .{ .ctx = &ctxs[i], .run = C.run };
        b.reset();
        var secs: [lib_cols.len]f64 = undefined;
        try bench.race(io, l1_timing, calls[0..count], &secs);
        if (lib_err) |err| return err;
        const bytes = l1Bytes(r, n, p.es());
        var rates: [lib_cols.len]?f64 = undefined;
        for (at, &rates) |j, *rt| rt.* = if (j) |k| bytes / secs[k] / 1e9 else null;
        std.debug.print("| {s} | {s} | {d} | {d} | {s} | ", .{ @tagName(r), p.name, n, case.inc, rstdRuns(name) });
        if (arm) std.debug.print("{s} | ", .{if (asm_k) |a| a.from else "none"});
        printRow(&lib_cols, &rates, &diffs, mine);
    }
}

fn libLevel2(io: Io, gpa: std.mem.Allocator) !void {
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
    for (a) |*v| v.* = @floatCast(rand.float(f64) * 2 - 1);
    for (x) |*v| v.* = @floatCast(rand.float(f64) * 2 - 1);
    for (y0) |*v| v.* = @floatCast(rand.float(f64) * 2 - 1);

    const asm_k = if (arm and !alone) asmFor(name) else null;
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
        var ctxs: [lib_cols.len]C = undefined;
        var calls: [lib_cols.len]bench.Call = undefined;
        var at: [lib_cols.len]?usize = @splat(null);
        var diffs: [lib_cols.len]?Diff = @splat(null);
        var count: usize = 0;
        var ok = true;
        var mine: Rstd = .{};
        for (0..lib_cols.len) |o| {
            const i = rstdFirst(o);
            if (!lib_cols[i].shown()) continue;
            const f = fns[i] orelse continue;
            ctxs[count] = .{ .f = f, .n = @intCast(n), .a = a.ptr, .x = x.ptr, .y = y.ptr, .buf = buf.ptr };
            @memcpy(y, y0);
            ctxs[count].call();
            if (lib_err) |err| return err;
            const e = checkGemv(p, t[0] == 't', n, a, x, y, y0);
            if (e > 1) {
                std.debug.print("WRONG: gemv_{s} {s} {s} n={d}: {d:.2} of tolerance\n", .{ t, p.name, lib_cols[i].label, n, e });
                ok = false;
            }
            if (i == rstd_col) {
                @memcpy(ky, y);
                mine = .{ .err = e * @as(f64, @floatFromInt(n + 8)), .bits = std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(y[0..out])) };
            } else if (lib_cols[i].lib != null) diffs[i] = .{ .eps = diffScaled(F, ky, y, scale[0..out]), .bound = sumBound(n) };
            at[i] = count;
            count += 1;
        }
        if (!ok) continue;
        for (0..count) |i| calls[i] = .{ .ctx = &ctxs[i], .run = C.run };
        @memcpy(y, y0);
        var secs: [lib_cols.len]f64 = undefined;
        try bench.race(io, l1_timing, calls[0..count], &secs);
        if (lib_err) |err| return err;
        const nf: f64 = @floatFromInt(n);
        const flops = nf * nf * @as(f64, if (p.complex) 8 else 2);
        var rates: [lib_cols.len]?f64 = undefined;
        for (at, &rates) |j, *rt| rt.* = if (j) |k| flops / secs[k] / 1e9 else null;
        std.debug.print("| gemv_{s} | {s} | {d} | {s} | ", .{ t, p.name, n, rstdRuns(name) });
        if (arm) std.debug.print("{s} | ", .{if (asm_k) |k| k.from else "none"});
        printRow(&lib_cols, &rates, &diffs, mine);
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
        fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// `checkGemm` for square `n`, complex too.
fn checkGemmP(comptime p: Prec, n: usize, a: []const p.F, b: []const p.F, c: []const p.F, c0: []const p.F, samples: usize) f64 {
    if (!p.complex) return checkGemm(p.F, n, n, n, a, b, c, c0, samples);
    const eps: f64 = std.math.floatEps(p.F);
    var rng = std.Random.DefaultPrng.init(3);
    const rand = rng.random();
    var worst: f64 = 0;
    for (0..samples) |_| {
        const i = rand.uintLessThan(usize, n);
        const j = rand.uintLessThan(usize, n);
        const at = 2 * (i + j * n);
        var re: f64 = c0[at];
        var im: f64 = c0[at + 1];
        var mag: f64 = @abs(re) + @abs(im);
        for (0..n) |l| {
            const ai = 2 * (i + l * n);
            const bi = 2 * (l + j * n);
            const ar: f64 = a[ai];
            const aim: f64 = a[ai + 1];
            const br: f64 = b[bi];
            const bim: f64 = b[bi + 1];
            re += ar * br - aim * bim;
            im += ar * bim + aim * br;
            mag += (@abs(ar) + @abs(aim)) * (@abs(br) + @abs(bim));
        }
        const tol = (@as(f64, @floatFromInt(n)) + 8) * eps * mag;
        worst = @max(worst, @abs(c[at] - re) / tol, @abs(c[at + 1] - im) / tol);
    }
    if (std.math.isNan(worst)) return std.math.inf(f64);
    return worst;
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

fn libLevel3(io: Io, gpa: std.mem.Allocator) !void {
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
        if (!alone) try gemmScale(p, n, a[0..m], b[0..m], c0[0..m], abs_a, abs_b, scale[0..m]);
        var ctxs: [lib3_cols.len]LC = undefined;
        var calls: [lib3_cols.len]bench.Call = undefined;
        var at: [lib3_cols.len]?usize = @splat(null);
        var diffs: [lib3_cols.len]?Diff = @splat(null);
        var count: usize = 0;
        var mine: Rstd = .{};
        for (0..lib3_cols.len) |o| {
            const i = rstdFirst(o);
            if (!lib3_cols[i].shown()) continue;
            const loops: ?Loops(p) = if (comptime Loops(p) != void) (if (i == 4) .{ .im = asmGemm(p), .blk = blk, .n = n, .a = a.ptr, .b = b.ptr, .c = c.ptr, .sa = sa.ptr, .sb = sb.ptr } else null) else null;
            if (fns[i] == null and loops == null) continue;
            ctxs[count] = .{ .f = fns[i], .loops = loops, .n = n, .a = a.ptr, .b = b.ptr, .c = c.ptr };
            @memcpy(c[0..m], c0[0..m]);
            ctxs[count].call();
            if (lib_err) |err| return err;
            const e = checkGemmP(p, n, a[0..m], b[0..m], c[0..m], c0[0..m], 256);
            if (e > 1) {
                std.debug.print("WRONG: {s}gemm {s} n={d}: {d:.2} of tolerance\n", .{ p.name, lib3_cols[i].label, n, e });
                return error.WrongGemm;
            }
            if (i == rstd_col) {
                @memcpy(kc[0..m], c[0..m]);
                mine = .{ .err = e * @as(f64, @floatFromInt(n + 8)), .bits = std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(c[0..m])) };
            } else if (lib3_cols[i].lib != null) diffs[i] = .{ .eps = diffScaled(F, kc, c, scale[0..m]), .bound = sumBound(n) };
            at[i] = count;
            count += 1;
        }
        for (0..count) |i| calls[i] = .{ .ctx = &ctxs[i], .run = LC.run };
        var secs: [lib3_cols.len]f64 = undefined;
        try bench.race(io, l3_timing, calls[0..count], &secs);
        if (lib_err) |err| return err;
        const nf: f64 = @floatFromInt(n);
        const flops = nf * nf * nf * @as(f64, if (p.complex) 8 else 2);
        var rates: [lib3_cols.len]?f64 = undefined;
        for (at, &rates) |j, *rt| rt.* = if (j) |k| flops / secs[k] / 1e9 else null;
        std.debug.print("| {s} | {d} | ", .{ p.name, n });
        printRow(&lib3_cols, &rates, &diffs, mine);
    }
}

// ---- table-driven races ------------------------------------------------------
//
// The routines beyond `gemv` and `gemm` (`ext_l2`, `ext_l3`) and LAPACK
// (`lapack_ops`), one entry each, a struct of declarations:
//
// - `name`; `shape` (`Shape`); `fill` (`A`'s `Fill`, default `general`);
//   `prep` (`Prep`, default `none`); `restore` (what a call overwrites and
//   must find fresh again; default none, it accumulates); `out` (the `Buf`
//   checked, compared and hashed); `cols` (`out`'s columns where it is `b`
//   or `c` of an `mm` shape, default n); `level` (default its table's);
//   `band_rows` (a `band` shape's, default `2 band_k + 1`);
// - `flops(n)`, real; a complex call counts 4 times as many;
// - `sym(p)` (default `p.name ++ name`), `Sig(p)` and `lib(p, f, w)`: the
//   Fortran-ABI call (`lookupF`); a library without it shows `-`;
// - `rstd(p, w)`, absent until rstd has the routine; `has_rstd` gates one
//   written ahead of it;
// - `residual(p, w, i, j)`: `|want - got|` at output entry (i, j) and the
//   size of its terms, null off the stored part; without it only the
//   agreement with rstd is shown.
//
// Adding a routine is adding an entry to its table.

/// Where an entry's orders come from: `l2_bytes`, `gemm_sizes`,
/// `lapack_sizes`.
const Level = enum { l2, l3, lapack };

/// LAPACK's orders: `gemm_sizes` less the largest, where reference LAPACK
/// and `prepare` take many seconds a call.
const lapack_sizes = gemm_sizes[0 .. gemm_sizes.len - 1].*;

const Shape = enum {
    /// `A` n by n, `x` in `b` and `y` in `c`, n long; n from `l2_bytes`
    mv,
    /// `A` in band storage, `band_rows` by n; `x` and `y` as `mv`; n from
    /// `l2_bytes` over the rows
    band,
    /// `A`, `B` and `C` n by n
    mm,
};

/// `gbmv`'s and `gbtrf`'s `kl` and `ku`.
const band_k = 32;

const Fill = enum {
    /// every float in [-1, 1]
    general,
    /// `general` with 2n added to the diagonal's real parts: strictly
    /// diagonally dominant by rows and columns, so LU with partial pivoting
    /// exchanges no rows and every triangle is well conditioned
    dominant,
    /// `dominant` made Hermitian (symmetric if real), its diagonal real:
    /// positive definite
    spd,
};

const Prep = enum {
    none,
    /// `A` replaced by its Cholesky factor `U`, `A = Uᴴ U` (`prepare`)
    potrf,
    /// `A` replaced by `L` and `U` of `A = L U`, pivots 1 to n (`dominant`
    /// takes no row exchange; `prepare`)
    getrf,
};

/// `orig` is `A` before `prep`; `a0`, `b0` and `c0` what each call starts
/// from; `s` singular values, eigenvalues or `tau`.
const Buf = enum { a, a0, orig, b, b0, c, c0, s, u, vt };

const Restore = struct { a: bool = false, b: bool = false, c: bool = false };

const Cols = enum { one, n };

/// The table-driven races' contestants: the libraries' Fortran-ABI routines
/// and rstd's public ones. `max_n` holds in level 3 and LAPACK only.
const ext_cols = lib_cols[0..4].* ++ [_]Col{
    .{ .label = "reference", .short = "ref", .lib = .reference, .max_n = 1024 },
    .{ .label = "eigen", .short = "eigen", .lib = .eigen },
    .{ .label = "aocl-libflame", .short = "aocl", .lib = .aoclflame },
};

const ext_l2 = .{ Ger, Symv, Trsv, Gbmv };
const ext_l3 = .{ Gemm3, Symm, Syrk, Trmm, Trsm };
/// rstd's LAPACK, then the hooks: T34's SVD (faer), T111's QR and
/// eigensolver, T110's ReLAPACK set; each races the libraries now and
/// takes rstd's column once its routine lands.
const lapack_ops = .{ Potrf, Getrf, Getrs, Gesv, Trtri, Potri, Gesvd, Gesdd, Geqrf, Syevd, Sytrf, Gbtrf };

fn libLapack(io: Io, gpa: std.mem.Allocator) !void {
    sectionMark("lapack");
    try clock(io, "lapack start");
    try extTable(io, gpa, .lapack, lapack_ops);
}

fn extTable(io: Io, gpa: std.mem.Allocator, comptime lv: Level, comptime ops: anytype) !void {
    const t = if (lv == .l2) l1_timing else l3_timing;
    std.debug.print("\n{s} GFLOP/s, best of {d} batches of at least {d} ms; ratios are over OpenBLAS. {s} {s}\n\n", .{ ext_leads[@backingInt(lv)], t.batches, @divTrunc(t.min_ns, std.time.ns_per_ms), io_note, ext_note });
    printHeader("routine | p | n", "------- | - | -:", &ext_cols);
    inline for (ops) |E| {
        inline for (precs) |p| try raceExt(io, gpa, E, lv, p);
    }
}

const ext_leads = [_][]const u8{
    "Level 2 beyond `gemv`, each library's Fortran-ABI routine, unit strides, n as `gemv`'s (`A` of 32 KiB, 4 MiB, 64 MiB): `ger` `A += 1e-6 x yᵀ` (`geru` complex), `symv` `y += 1e-3 S x` (upper), `trsv` upper, no transpose, non-unit on a diagonally dominant `A`, `x` restored before each call (counted), `gbmv` `y += 1e-3 A x` with kl = ku = 32 (`A`'s band of the same bytes); rates by 2n², 2n², n², 2n (kl + ku + 1), four times as many complex.",
    "Level 3 beyond `gemm`, each library's Fortran-ABI routine, square n, `alpha` and `beta` 1: `gemm` NN, `symm` left upper, `syrk` upper NN, `trmm` and `trsm` left upper NN non-unit (`trsm` on a diagonally dominant `A`), `B` restored before each `trmm` and `trsm` call (an n² copy, counted); rates by 2n³, 2n³, n²(n + 1), n³, n³, four times as many complex; `reference` to n = 1024.",
    "LAPACK, square n, upper where it matters: rates by LAPACK's counts, `potrf` n³/3, `getrf` 2n³/3, `getrs` 2n³ (n right-hand sides), `gesv` 2n³/3 + 2n² (one), `trtri` n³/3, `potri` 2n³/3, and nominal ones: `gesvd` and `gesdd` 22n³ (thin `U` and `Vᵀ`, Golub and Van Loan's Golub-Reinsch count), `geqrf` 4n³/3, `syevd` (`heevd` complex) 9n³ with vectors, `sytrf` n³/3, `gbtrf` 4nk² (kl = ku = k = 32, at `gemv`'s band sizes); four times as many complex. `A` is diagonally dominant for `getrf`, `getrs`, `gesv` and `trtri` (no row exchanges, so pivots agree), Hermitian positive definite for `potrf`, `potri`, `syevd` and `sytrf`; `getrs` and `potri` start from a factorization made here in f64. Each call restores what it overwrites first (an n² copy, counted). `gesvd`, `gesdd`, `geqrf`, `syevd`, `sytrf` and `gbtrf` race the libraries until rstd has them.",
};

const ext_note = "`rstd err` is rstd's largest residual at the checked entries (all, or 256 sampled), in eps of the type relative to the size of its terms; `rstd-<lib>` is rstd's largest difference from that library's output in eps relative to the output's largest entry, `!` past 2 (terms + 8). `rstd bits` hashes rstd's output. Eigen runs one thread; `-` where a library lacks the routine";

/// A complex value, or a real one with 0 imaginary part, in f64.
const Pair = [2]f64;

fn cadd(a: Pair, b: Pair) Pair {
    return .{ a[0] + b[0], a[1] + b[1] };
}
fn csub(a: Pair, b: Pair) Pair {
    return .{ a[0] - b[0], a[1] - b[1] };
}
fn cmul(a: Pair, b: Pair) Pair {
    return .{ a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0] };
}
fn cdiv(a: Pair, b: Pair) Pair {
    const d = b[0] * b[0] + b[1] * b[1];
    return .{ (a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d };
}
fn cconj(a: Pair) Pair {
    return .{ a[0], -a[1] };
}
/// `|re| + |im|`, the size `checkGemv` gives a complex.
fn cmag(a: Pair) f64 {
    return @abs(a[0]) + @abs(a[1]);
}
fn cdist(a: Pair, b: Pair) f64 {
    return @max(@abs(a[0] - b[0]), @abs(a[1] - b[1]));
}
/// `v` rounded to `p`'s type, as the routines get it.
fn cnum(comptime p: Prec, v: f64) Pair {
    return .{ @as(p.F, @floatCast(v)), 0 };
}

/// Fortran's `CHARACTER*1` arguments, by reference; each call also passes
/// their hidden lengths (1) last, which a library built without them
/// ignores.
const fc = struct {
    const U: u8 = 'U';
    const L: u8 = 'L';
    const N: u8 = 'N';
    const S: u8 = 'S';
    const V: u8 = 'V';
};
const Ch = *const u8;
const In = *const c_int;
const one_c: c_int = 1;
const band_kc: c_int = band_k;

/// A Fortran scalar of `p` by reference: real, or (re, im).
fn Sc(comptime p: Prec) type {
    return if (p.complex) [2]p.F else p.F;
}
fn sc(comptime p: Prec, v: f64) Sc(p) {
    return if (p.complex) .{ @floatCast(v), 0 } else @floatCast(v);
}
/// The same scalar as rstd takes it.
fn num(comptime p: Prec, v: f64) Num(p) {
    return if (p.complex) .init(@floatCast(v), 0) else @floatCast(v);
}

fn decl(comptime E: type, comptime name: []const u8, comptime default: anytype) @TypeOf(default) {
    return if (@hasDecl(E, name)) @field(E, name) else default;
}
fn hasRstd(comptime E: type) bool {
    return decl(E, "has_rstd", @hasDecl(E, "rstd"));
}
fn bandRows(comptime E: type) usize {
    return decl(E, "band_rows", @as(usize, 2 * band_k + 1));
}
fn symOf(comptime E: type, comptime p: Prec) []const u8 {
    return if (@hasDecl(E, "sym")) E.sym(p) else p.name ++ E.name;
}

/// An entry's orders at level `lv`.
fn orders(comptime E: type, comptime lv: Level, comptime p: Prec) []const usize {
    return comptime switch (lv) {
        .l2 => blk: {
            var ns: [l2_bytes.len]usize = undefined;
            for (l2_bytes, &ns) |bytes, *n| n.* = switch (E.shape) {
                .band => bytes / (p.es() * bandRows(E)),
                .mv, .mm => std.math.sqrt(bytes / p.es()),
            };
            const out = ns;
            break :blk &out;
        },
        .l3 => &gemm_sizes,
        .lapack => &lapack_sizes,
    };
}

/// Elements of `A`, and of `B` and `C`, at order n.
fn lenA(comptime E: type, n: usize) usize {
    return switch (E.shape) {
        .mv, .mm => n * n,
        .band => bandRows(E) * n,
    };
}
fn lenV(comptime E: type, n: usize) usize {
    return switch (E.shape) {
        .mv, .band => n,
        .mm => n * n,
    };
}

/// Rows and columns of `E`'s output at order n.
fn outDims(comptime E: type, n: usize) [2]usize {
    return switch (E.out) {
        .a => .{ if (E.shape == .band) bandRows(E) else n, n },
        .b, .c => if (E.shape == .mm) .{ n, if (decl(E, "cols", Cols.n) == .one) 1 else n } else .{ n, 1 },
        .s => .{ n, 1 },
        else => @compileError("an output is a, b, c or s"),
    };
}

/// Floats of `E`'s output at order n; `s` is real.
fn outLen(comptime E: type, comptime p: Prec, n: usize) usize {
    const d = outDims(E, n);
    return d[0] * d[1] * (if (E.out == .s) 1 else p.cs());
}

/// Terms of one output entry, for the bounds.
fn terms(comptime E: type, n: usize) usize {
    return if (E.shape == .band) bandRows(E) else n;
}

/// A residual past this many eps is wrong, not rounding.
fn wrongLimit(comptime E: type, n: usize) f64 {
    return 64 * @as(f64, @floatFromInt(terms(E, n) + 8));
}

/// An entry's operands at its largest order, shared by its contestants.
fn Work(comptime p: Prec) type {
    return struct {
        const W = @This();
        const F = p.F;
        const cs = p.cs();

        n: usize = 0,
        a: []align(64) F,
        a0: []align(64) F,
        orig: []align(64) F,
        b: []align(64) F,
        b0: []align(64) F,
        c: []align(64) F,
        c0: []align(64) F,
        s: []align(64) F,
        u: []align(64) F,
        vt: []align(64) F,
        work: []align(64) F,
        rwork: []align(64) F,
        ipiv: []c_int,
        ipiv0: []c_int,
        rpiv: []isize,
        rpiv0: []isize,
        iwork: []c_int,
        info: c_int = 0,

        fn init(gpa: std.mem.Allocator, comptime E: type, comptime lv: Level, n: usize) !W {
            const la = lenA(E, n) * cs;
            const lb = lenV(E, n) * cs;
            const big = lv == .lapack;
            var w: W = undefined;
            w.n = 0;
            w.info = 0;
            inline for (.{ .{ "a", la }, .{ "a0", la }, .{ "orig", la }, .{ "b", lb }, .{ "b0", lb }, .{ "c", lb }, .{ "c0", lb } }) |f|
                @field(w, f[0]) = try gpa.alignedAlloc(F, .@"64", f[1]);
            w.s = try gpa.alignedAlloc(F, .@"64", n * cs);
            const sq = if (big) n * n * cs else 0;
            w.u = try gpa.alignedAlloc(F, .@"64", sq);
            w.vt = try gpa.alignedAlloc(F, .@"64", sq);
            w.work = try gpa.alignedAlloc(F, .@"64", if (big) lwork(n) * cs else 0);
            w.rwork = try gpa.alignedAlloc(F, .@"64", if (big) lrwork(n) else 0);
            w.iwork = try gpa.alloc(c_int, if (big) liwork(n) else 0);
            w.ipiv = try gpa.alloc(c_int, n);
            w.ipiv0 = try gpa.alloc(c_int, n);
            w.rpiv = try gpa.alloc(isize, n);
            w.rpiv0 = try gpa.alloc(isize, n);
            return w;
        }

        fn deinit(w: *W, gpa: std.mem.Allocator) void {
            inline for (.{ "a", "a0", "orig", "b", "b0", "c", "c0", "s", "u", "vt", "work", "rwork", "iwork", "ipiv", "ipiv0", "rpiv", "rpiv0" }) |f| gpa.free(@field(w, f));
        }

        /// Workspace for every routine's blocked path: past each one's
        /// documented need (`gesdd`'s 4n² + 7n, `syevd`'s 2n² + 6n + 1,
        /// `n nb` with `nb` up to 64).
        fn lwork(n: usize) usize {
            return 5 * n * n + 128 * n;
        }
        fn lrwork(n: usize) usize {
            return 5 * n * n + 64 * n;
        }
        fn liwork(n: usize) usize {
            return 8 * n + 64;
        }

        fn ni(w: *const W) c_int {
            return @intCast(w.n);
        }
        fn in(w: *const W) isize {
            return @intCast(w.n);
        }
        fn lw(w: *const W) c_int {
            return @intCast(lwork(w.n));
        }
        fn lrw(w: *const W) c_int {
            return @intCast(lrwork(w.n));
        }
        fn liw(w: *const W) c_int {
            return @intCast(liwork(w.n));
        }

        fn buf(w: *const W, comptime b: Buf) []const F {
            return @field(w, @tagName(b));
        }

        /// Element k of `b`, as (re, im); `s` is real.
        fn at(w: *const W, comptime b: Buf, k: usize) Pair {
            const v = w.buf(b);
            if (b == .s or !p.complex) return .{ v[k], 0 };
            return .{ v[2 * k], v[2 * k + 1] };
        }

        fn set(v: []F, k: usize, x: Pair) void {
            if (p.complex) {
                v[2 * k] = @floatCast(x[0]);
                v[2 * k + 1] = @floatCast(x[1]);
            } else v[k] = @floatCast(x[0]);
        }

        /// `E`'s output at the current order.
        fn out(w: *const W, comptime E: type) []const F {
            return w.buf(E.out)[0..outLen(E, p, w.n)];
        }

        /// Every operand fresh, before a check.
        fn reset(w: *W, comptime E: type) void {
            const la = lenA(E, w.n) * cs;
            const lb = lenV(E, w.n) * cs;
            @memcpy(w.a[0..la], w.a0[0..la]);
            @memcpy(w.b[0..lb], w.b0[0..lb]);
            @memcpy(w.c[0..lb], w.c0[0..lb]);
            @memcpy(w.ipiv[0..w.n], w.ipiv0[0..w.n]);
            @memcpy(w.rpiv[0..w.n], w.rpiv0[0..w.n]);
            w.info = 0;
        }

        /// What a call of `E` overwrites, fresh; inside its timing.
        fn restore(w: *W, comptime E: type) void {
            const r = decl(E, "restore", Restore{});
            if (r.a) {
                const la = lenA(E, w.n) * cs;
                @memcpy(w.a[0..la], w.a0[0..la]);
            }
            const lb = lenV(E, w.n) * cs;
            if (r.b) @memcpy(w.b[0..lb], w.b0[0..lb]);
            if (r.c) @memcpy(w.c[0..lb], w.c0[0..lb]);
        }

        /// The operands at order n: random, as `E`'s fill, and prepared.
        fn setup(w: *W, gpa: std.mem.Allocator, comptime E: type, n: usize) !void {
            w.n = n;
            var rng = std.Random.DefaultPrng.init(5);
            const r = rng.random();
            const la = lenA(E, n) * cs;
            const lb = lenV(E, n) * cs;
            for (w.orig[0..la]) |*v| v.* = @floatCast(r.float(f64) * 2 - 1);
            for (w.b0[0..lb]) |*v| v.* = @floatCast(r.float(f64) * 2 - 1);
            for (w.c0[0..lb]) |*v| v.* = @floatCast(r.float(f64) * 2 - 1);
            const fill = comptime decl(E, "fill", Fill.general);
            if (fill != .general) {
                if (E.shape == .band) @compileError("a band is filled general");
                if (fill == .spd) {
                    for (0..n) |j| for (0..j) |i| set(w.orig, j + i * n, cconj(w.at(.orig, i + j * n)));
                    for (0..n) |i| set(w.orig, i + i * n, .{ w.at(.orig, i + i * n)[0], 0 });
                }
                const two_n: F = @floatFromInt(2 * n);
                for (0..n) |i| w.orig[(i + i * n) * cs] += two_n;
            }
            @memcpy(w.a0[0..la], w.orig[0..la]);
            for (w.ipiv0[0..n], w.rpiv0[0..n], 1..) |*d, *e, k| {
                d.* = @intCast(k);
                e.* = @intCast(k);
            }
            const prep = comptime decl(E, "prep", Prep.none);
            if (prep != .none) try prepare(gpa, p, w, prep);
        }
    };
}

/// `prep`'s factorization of `orig` into `a0`, unblocked in f64: `potrf`
/// writes `U` over the upper triangle, `getrf` `L` and `U` over all of it.
fn prepare(gpa: std.mem.Allocator, comptime p: Prec, w: *Work(p), comptime what: Prep) !void {
    const n = w.n;
    const m = try gpa.alloc(Pair, n * n);
    defer gpa.free(m);
    for (m, 0..) |*v, k| v.* = w.at(.orig, k);
    switch (what) {
        .none => return,
        .potrf => for (0..n) |j| {
            for (0..j + 1) |i| {
                var s = m[i + j * n];
                for (0..i) |k| s = csub(s, cmul(cconj(m[k + i * n]), m[k + j * n]));
                m[i + j * n] = if (i < j) cdiv(s, m[i + i * n]) else .{ @sqrt(s[0]), 0 };
            }
        },
        .getrf => for (0..n) |k| {
            const piv = m[k + k * n];
            for (k + 1..n) |i| m[i + k * n] = cdiv(m[i + k * n], piv);
            for (k + 1..n) |j| {
                const u = m[k + j * n];
                for (k + 1..n) |i| m[i + j * n] = csub(m[i + j * n], cmul(m[i + k * n], u));
            }
        },
    }
    for (0..n) |j| for (0..n) |i| {
        if (what == .potrf and i > j) continue;
        Work(p).set(w.a0, i + j * n, m[i + j * n]);
    };
}

/// A contestant's largest residual over `E`'s output in eps of the type
/// relative to the size of its terms: every entry of a small output, else
/// 256 sampled; null without `residual`.
fn checkExt(comptime E: type, comptime p: Prec, w: *const Work(p)) ?f64 {
    if (!@hasDecl(E, "residual")) return null;
    const eps: f64 = std.math.floatEps(p.F);
    const d = outDims(E, w.n);
    var worst: f64 = 0;
    if (d[0] * d[1] <= 4096) {
        for (0..d[1]) |j| for (0..d[0]) |i| {
            const r = E.residual(p, w, i, j) orelse continue;
            worst = @max(worst, r[0] / (eps * @max(r[1], 1e-300)));
        };
    } else {
        var rng = std.Random.DefaultPrng.init(7);
        const rand = rng.random();
        var hits: usize = 0;
        var tries: usize = 0;
        while (hits < 256 and tries < 4096) : (tries += 1) {
            const r = E.residual(p, w, rand.uintLessThan(usize, d[0]), rand.uintLessThan(usize, d[1])) orelse continue;
            hits += 1;
            worst = @max(worst, r[0] / (eps * @max(r[1], 1e-300)));
        }
    }
    return if (std.math.isNan(worst)) std.math.inf(f64) else worst;
}

/// rstd's largest difference from `got` in eps relative to its own largest
/// entry.
fn diffNorm(comptime F: type, mine: []const F, got: []const F) f64 {
    var big: f64 = 0;
    var worst: f64 = 0;
    for (mine, got) |a, b| {
        big = @max(big, @abs(@as(f64, a)));
        worst = @max(worst, @abs(@as(f64, a) - b));
    }
    const r = worst / (std.math.floatEps(F) * @max(big, 1e-300));
    return if (std.math.isNan(r)) std.math.inf(f64) else r;
}

fn ExtCall(comptime E: type, comptime p: Prec) type {
    return struct {
        /// null: rstd's
        f: ?*const E.Sig(p),
        w: *Work(p),
        fn call(s: *const @This()) void {
            s.w.restore(E);
            if (s.f) |f| return E.lib(p, f, s.w);
            if (comptime hasRstd(E)) E.rstd(p, s.w) catch |e| fail(e);
        }
        fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

fn raceExt(io: Io, gpa: std.mem.Allocator, comptime E: type, comptime table_level: Level, comptime p: Prec) !void {
    const lv: Level = comptime decl(E, "level", table_level);
    const ns = comptime orders(E, lv, p);
    var w = try Work(p).init(gpa, E, lv, ns[ns.len - 1]);
    defer w.deinit(gpa);
    const kout = try gpa.alloc(p.F, outLen(E, p, ns[ns.len - 1]));
    defer gpa.free(kout);
    const Sig = *const E.Sig(p);
    const sym = comptime symOf(E, p);
    const fns = [ext_cols.len]?Sig{
        lookupF(.openblas, Sig, sym),
        null,
        lookupF(.blis, Sig, sym),
        lookupF(.accelerate, Sig, sym),
        lookupF(.reference, Sig, sym),
        lookupF(.eigen, Sig, sym),
        lookupF(.aoclflame, Sig, sym),
    };
    const EC = ExtCall(E, p);
    const timing = if (lv == .l2) l1_timing else l3_timing;
    for (ns) |n| {
        if (n > n_max) continue;
        try w.setup(gpa, E, n);
        const len = outLen(E, p, n);
        var ctxs: [ext_cols.len]EC = undefined;
        var calls: [ext_cols.len]bench.Call = undefined;
        var at: [ext_cols.len]?usize = @splat(null);
        var diffs: [ext_cols.len]?Diff = @splat(null);
        var count: usize = 0;
        var ok = true;
        var mine: Rstd = .{};
        var have_mine = false;
        for (0..ext_cols.len) |o| {
            const i = rstdFirst(o);
            const col = ext_cols[i];
            if (!col.shown() or (lv != .l2 and n > col.max_n)) continue;
            if (i == rstd_col) {
                if (!comptime hasRstd(E)) continue;
            } else if (fns[i] == null) continue;
            ctxs[count] = .{ .f = fns[i], .w = &w };
            w.reset(E);
            ctxs[count].call();
            if (lib_err) |err| return err;
            const e = checkExt(E, p, &w);
            if (w.info != 0 or (e != null and e.? > wrongLimit(E, n))) {
                std.debug.print("WRONG: {s} {s} {s} n={d}: info {d}, residual {?d:.1} eps\n", .{ E.name, p.name, col.label, n, w.info, e });
                ok = false;
            }
            const got = w.out(E);
            if (i == rstd_col) {
                @memcpy(kout[0..len], got);
                mine = .{ .err = e, .bits = std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(got)) };
                have_mine = true;
            } else if (have_mine and col.lib != null) diffs[i] = .{ .eps = diffNorm(p.F, kout[0..len], got), .bound = sumBound(terms(E, n)) };
            at[i] = count;
            count += 1;
        }
        if (!ok or count == 0) continue;
        for (0..count) |k| calls[k] = .{ .ctx = &ctxs[k], .run = EC.run };
        w.reset(E);
        var secs: [ext_cols.len]f64 = undefined;
        try bench.race(io, timing, calls[0..count], &secs);
        if (lib_err) |err| return err;
        const flops = E.flops(@floatFromInt(n)) * @as(f64, if (p.complex) 4 else 1);
        var rates: [ext_cols.len]?f64 = undefined;
        for (at, &rates) |j, *rt| rt.* = if (j) |k| flops / secs[k] / 1e9 else null;
        std.debug.print("| {s} | {s} | {d} | ", .{ E.name, p.name, n });
        printRow(&ext_cols, &rates, &diffs, mine);
    }
}

/// `A[i, k]` of a symmetric (`herm` false) or Hermitian matrix whose upper
/// triangle `b` stores, order n.
fn symAt(comptime p: Prec, w: *const Work(p), comptime b: Buf, comptime herm: bool, i: usize, k: usize) Pair {
    const n = w.n;
    if (i <= k) return w.at(b, i + k * n);
    const v = w.at(b, k + i * n);
    return if (herm) cconj(v) else v;
}

// ---- level 2 entries

const Ger = struct {
    pub const name = "ger";
    pub const shape: Shape = .mv;
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return 2 * n * n;
    }
    pub fn sym(comptime p: Prec) []const u8 {
        return p.name ++ (if (p.complex) "geru" else "ger");
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (In, In, *const Sc(p), [*]p.F, In, [*]p.F, In, [*]p.F, In) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const al = sc(p, axpy_alpha);
        f(&n, &n, &al, w.b.ptr, &one_c, w.c.ptr, &one_c, w.a.ptr, &n);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        try blas.ger(rstd_alloc, rstd_io, Num(p), w.in(), w.in(), num(p, axpy_alpha), cast(p, w.b.ptr), 1, cast(p, w.c.ptr), 1, cast(p, w.a.ptr), w.in());
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const k = i + j * w.n;
        const al = cnum(p, axpy_alpha);
        const x = w.at(.b0, i);
        const y = w.at(.c0, j);
        return .{ cdist(cadd(w.at(.a0, k), cmul(al, cmul(x, y))), w.at(.a, k)), cmag(w.at(.a0, k)) + cmag(al) * cmag(x) * cmag(y) };
    }
};

const Symv = struct {
    pub const name = "symv";
    pub const shape: Shape = .mv;
    pub const out: Buf = .c;
    pub fn flops(n: f64) f64 {
        return 2 * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, In, *const Sc(p), [*]p.F, In, [*]p.F, In, *const Sc(p), [*]p.F, In, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const al = sc(p, gemv_alpha);
        const one = sc(p, 1);
        f(&fc.U, &n, &al, w.a.ptr, &n, w.b.ptr, &one_c, &one, w.c.ptr, &one_c, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        try blas.symv(rstd_alloc, rstd_io, Num(p), .u, w.in(), num(p, gemv_alpha), cast(p, w.a.ptr), w.in(), cast(p, w.b.ptr), 1, num(p, 1), cast(p, w.c.ptr), 1);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, _: usize) ?[2]f64 {
        var s: Pair = .{ 0, 0 };
        var mag: f64 = 0;
        for (0..w.n) |k| {
            const a = symAt(p, w, .a0, false, i, k);
            const x = w.at(.b0, k);
            s = cadd(s, cmul(a, x));
            mag += cmag(a) * cmag(x);
        }
        const al = cnum(p, gemv_alpha);
        const y0 = w.at(.c0, i);
        return .{ cdist(cadd(y0, cmul(al, s)), w.at(.c, i)), cmag(y0) + cmag(al) * mag };
    }
};

const Trsv = struct {
    pub const name = "trsv";
    pub const shape: Shape = .mv;
    pub const fill: Fill = .dominant;
    pub const restore: Restore = .{ .b = true };
    pub const out: Buf = .b;
    pub fn flops(n: f64) f64 {
        return n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, Ch, Ch, In, [*]p.F, In, [*]p.F, In, usize, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        f(&fc.U, &fc.N, &fc.N, &n, w.a.ptr, &n, w.b.ptr, &one_c, 1, 1, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        try blas.trsv(rstd_alloc, rstd_io, Num(p), .u, .n, .n, w.in(), cast(p, w.a.ptr), w.in(), cast(p, w.b.ptr), 1);
    }
    /// `A x` against the right-hand side.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, _: usize) ?[2]f64 {
        var s: Pair = .{ 0, 0 };
        var mag: f64 = 0;
        for (i..w.n) |k| {
            const a = w.at(.a0, i + k * w.n);
            const x = w.at(.b, k);
            s = cadd(s, cmul(a, x));
            mag += cmag(a) * cmag(x);
        }
        const r = w.at(.b0, i);
        return .{ cdist(s, r), mag + cmag(r) };
    }
};

const Gbmv = struct {
    pub const name = "gbmv";
    pub const shape: Shape = .band;
    pub const out: Buf = .c;
    pub fn flops(n: f64) f64 {
        return 2 * n * @as(f64, @floatFromInt(bandRows(@This())));
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, In, In, In, In, *const Sc(p), [*]p.F, In, [*]p.F, In, *const Sc(p), [*]p.F, In, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lda: c_int = @intCast(bandRows(@This()));
        const al = sc(p, gemv_alpha);
        const one = sc(p, 1);
        f(&fc.N, &n, &n, &band_kc, &band_kc, &al, w.a.ptr, &lda, w.b.ptr, &one_c, &one, w.c.ptr, &one_c, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        try blas.gbmv(rstd_alloc, rstd_io, Num(p), .n, w.in(), w.in(), band_k, band_k, num(p, gemv_alpha), cast(p, w.a.ptr), @intCast(bandRows(@This())), cast(p, w.b.ptr), 1, num(p, 1), cast(p, w.c.ptr), 1);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, _: usize) ?[2]f64 {
        const lda = bandRows(@This());
        var s: Pair = .{ 0, 0 };
        var mag: f64 = 0;
        for (i -| band_k..@min(w.n, i + band_k + 1)) |k| {
            const a = w.at(.a0, band_k + i - k + k * lda);
            const x = w.at(.b0, k);
            s = cadd(s, cmul(a, x));
            mag += cmag(a) * cmag(x);
        }
        const al = cnum(p, gemv_alpha);
        const y0 = w.at(.c0, i);
        return .{ cdist(cadd(y0, cmul(al, s)), w.at(.c, i)), cmag(y0) + cmag(al) * mag };
    }
};

// ---- level 3 entries

/// `c0 + A B` at (i, j) against `c`, `A` read by `a(w, i, k)`.
fn productAt(comptime p: Prec, w: *const Work(p), i: usize, j: usize, comptime a: fn (*const Work(p), usize, usize) Pair) [2]f64 {
    var s = w.at(.c0, i + j * w.n);
    var mag = cmag(s);
    for (0..w.n) |k| {
        const x = a(w, i, k);
        const y = w.at(.b0, k + j * w.n);
        s = cadd(s, cmul(x, y));
        mag += cmag(x) * cmag(y);
    }
    return .{ cdist(s, w.at(.c, i + j * w.n)), mag };
}

/// `trmm`'s and `trsm`'s Fortran signature: side, uplo, transa, diag, m,
/// n, alpha, a, lda, b, ldb.
fn SigTri(comptime p: Prec) type {
    return fn (Ch, Ch, Ch, Ch, In, In, *const Sc(p), [*]p.F, In, [*]p.F, In, usize, usize, usize, usize) callconv(.c) void;
}

/// `trmm`'s and `trsm`'s call: left, upper, no transpose, non-unit, `alpha` 1.
fn libTri(comptime p: Prec, f: *const SigTri(p), w: *Work(p)) void {
    const n = w.ni();
    const one = sc(p, 1);
    f(&fc.L, &fc.U, &fc.N, &fc.N, &n, &n, &one, w.a.ptr, &n, w.b.ptr, &n, 1, 1, 1, 1);
}

/// Row i of upper triangular `A` (`a0`) times column j of `x`, and the
/// size of its terms.
fn triRow(comptime p: Prec, w: *const Work(p), comptime x: Buf, i: usize, j: usize) [2]Pair {
    const n = w.n;
    var s: Pair = .{ 0, 0 };
    var mag: f64 = 0;
    for (i..n) |k| {
        const a = w.at(.a0, i + k * n);
        const v = w.at(x, k + j * n);
        s = cadd(s, cmul(a, v));
        mag += cmag(a) * cmag(v);
    }
    return .{ s, .{ mag, 0 } };
}

const Gemm3 = struct {
    pub const name = "gemm";
    pub const shape: Shape = .mm;
    pub const out: Buf = .c;
    pub fn flops(n: f64) f64 {
        return 2 * n * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, Ch, In, In, In, *const Sc(p), [*]p.F, In, [*]p.F, In, *const Sc(p), [*]p.F, In, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const one = sc(p, 1);
        f(&fc.N, &fc.N, &n, &n, &n, &one, w.a.ptr, &n, w.b.ptr, &n, &one, w.c.ptr, &n, 1, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try blas.gemm(rstd_alloc, rstd_io, Num(p), .n, .n, n, n, n, num(p, 1), cast(p, w.a.ptr), n, cast(p, w.b.ptr), n, num(p, 1), cast(p, w.c.ptr), n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return productAt(p, w, i, j, struct {
            fn a(v: *const Work(p), r: usize, k: usize) Pair {
                return v.at(.a0, r + k * v.n);
            }
        }.a);
    }
};

const Symm = struct {
    pub const name = "symm";
    pub const shape: Shape = .mm;
    pub const out: Buf = .c;
    pub fn flops(n: f64) f64 {
        return 2 * n * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, Ch, In, In, *const Sc(p), [*]p.F, In, [*]p.F, In, *const Sc(p), [*]p.F, In, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const one = sc(p, 1);
        f(&fc.L, &fc.U, &n, &n, &one, w.a.ptr, &n, w.b.ptr, &n, &one, w.c.ptr, &n, 1, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try blas.symm(rstd_alloc, rstd_io, Num(p), .l, .u, n, n, num(p, 1), cast(p, w.a.ptr), n, cast(p, w.b.ptr), n, num(p, 1), cast(p, w.c.ptr), n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return productAt(p, w, i, j, struct {
            fn a(v: *const Work(p), r: usize, k: usize) Pair {
                return symAt(p, v, .a0, false, r, k);
            }
        }.a);
    }
};

const Syrk = struct {
    pub const name = "syrk";
    pub const shape: Shape = .mm;
    pub const out: Buf = .c;
    pub fn flops(n: f64) f64 {
        return n * n * (n + 1);
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, Ch, In, In, *const Sc(p), [*]p.F, In, *const Sc(p), [*]p.F, In, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const one = sc(p, 1);
        f(&fc.U, &fc.N, &n, &n, &one, w.a.ptr, &n, &one, w.c.ptr, &n, 1, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try blas.syrk(rstd_alloc, rstd_io, Num(p), .u, .n, n, n, num(p, 1), cast(p, w.a.ptr), n, num(p, 1), cast(p, w.c.ptr), n);
    }
    /// `c0 + A Aᵀ`, upper triangle.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        if (i > j) return null;
        const n = w.n;
        var s = w.at(.c0, i + j * n);
        var mag = cmag(s);
        for (0..n) |k| {
            const x = w.at(.a0, i + k * n);
            const y = w.at(.a0, j + k * n);
            s = cadd(s, cmul(x, y));
            mag += cmag(x) * cmag(y);
        }
        return .{ cdist(s, w.at(.c, i + j * n)), mag };
    }
};

const Trmm = struct {
    pub const name = "trmm";
    pub const shape: Shape = .mm;
    pub const restore: Restore = .{ .b = true };
    pub const out: Buf = .b;
    pub fn flops(n: f64) f64 {
        return n * n * n;
    }
    pub const Sig = SigTri;
    pub const lib = libTri;
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try blas.trmm(rstd_alloc, rstd_io, Num(p), .l, .u, .n, .n, n, n, num(p, 1), cast(p, w.a.ptr), n, cast(p, w.b.ptr), n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const r = triRow(p, w, .b0, i, j);
        return .{ cdist(r[0], w.at(.b, i + j * w.n)), r[1][0] };
    }
};

const Trsm = struct {
    pub const name = "trsm";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .dominant;
    pub const restore: Restore = .{ .b = true };
    pub const out: Buf = .b;
    pub fn flops(n: f64) f64 {
        return n * n * n;
    }
    pub const Sig = SigTri;
    pub const lib = libTri;
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try blas.trsm(rstd_alloc, rstd_io, Num(p), .l, .u, .n, .n, n, n, num(p, 1), cast(p, w.a.ptr), n, cast(p, w.b.ptr), n);
    }
    /// `A X` against the right-hand sides.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const r = triRow(p, w, .b, i, j);
        const b = w.at(.b0, i + j * w.n);
        return .{ cdist(r[0], b), r[1][0] + cmag(b) };
    }
};

// ---- LAPACK entries

/// `orig X` at (i, j), over k in `lo..hi`, against the identity; `X` read
/// by `x(w, k, j)`: the residual of an inverse.
fn inverseAt(comptime p: Prec, w: *const Work(p), i: usize, j: usize, lo: usize, hi: usize, comptime x: fn (*const Work(p), usize, usize) Pair) [2]f64 {
    var s: Pair = .{ 0, 0 };
    var mag: f64 = 0;
    for (lo..hi) |k| {
        const a = w.at(.orig, i + k * w.n);
        const v = x(w, k, j);
        s = cadd(s, cmul(a, v));
        mag += cmag(a) * cmag(v);
    }
    const id: f64 = if (i == j) 1 else 0;
    return .{ cdist(s, .{ id, 0 }), mag + id };
}

/// `orig X` at (i, j) against the right-hand sides `b0`; the solution is
/// in `b`.
fn solveAt(comptime p: Prec, w: *const Work(p), i: usize, j: usize) [2]f64 {
    const n = w.n;
    var s: Pair = .{ 0, 0 };
    var mag: f64 = 0;
    for (0..n) |k| {
        const a = w.at(.orig, i + k * n);
        const x = w.at(.b, k + j * n);
        s = cadd(s, cmul(a, x));
        mag += cmag(a) * cmag(x);
    }
    const r = w.at(.b0, i + j * n);
    return .{ cdist(s, r), mag + cmag(r) };
}

/// `uplo, n, a, lda, info`: `potrf`, `potri`.
fn SigUplo(comptime p: Prec) type {
    return fn (Ch, In, [*]p.F, In, *c_int, usize) callconv(.c) void;
}

fn libUplo(comptime p: Prec, f: *const SigUplo(p), w: *Work(p)) void {
    const n = w.ni();
    f(&fc.U, &n, w.a.ptr, &n, &w.info, 1);
}

const Potrf = struct {
    pub const name = "potrf";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .spd;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return n * n * n / 3;
    }
    pub const Sig = SigUplo;
    pub const lib = libUplo;
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        try lapack.potrf(rstd_alloc, rstd_io, Num(p), .u, w.in(), cast(p, w.a.ptr), w.in(), null);
    }
    /// `Uᴴ U` against `A`, upper triangle.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        if (i > j) return null;
        const n = w.n;
        var s: Pair = .{ 0, 0 };
        var mag: f64 = 0;
        for (0..i + 1) |k| {
            const x = cconj(w.at(.a, k + i * n));
            const y = w.at(.a, k + j * n);
            s = cadd(s, cmul(x, y));
            mag += cmag(x) * cmag(y);
        }
        return .{ cdist(s, w.at(.orig, i + j * n)), mag };
    }
};

const Getrf = struct {
    pub const name = "getrf";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .dominant;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return 2 * n * n * n / 3;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (In, In, [*]p.F, In, [*]c_int, *c_int) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        f(&n, &n, w.a.ptr, &n, w.ipiv.ptr, &w.info);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        try lapack.getrf(rstd_alloc, rstd_io, Num(p), w.in(), w.in(), cast(p, w.a.ptr), w.in(), w.rpiv.ptr, null);
        for (w.ipiv[0..w.n], w.rpiv[0..w.n]) |*d, s| d.* = @intCast(s);
    }
    /// `L U` against `P A`.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const n = w.n;
        // the row of `A` the exchanges, made in order, bring to row i
        var row = i;
        var k = n;
        while (k > 0) {
            k -= 1;
            const q: usize = @intCast(w.ipiv[k] - 1);
            if (row == k) row = q else if (row == q) row = k;
        }
        var s: Pair = .{ 0, 0 };
        var mag: f64 = 0;
        for (0..@min(i, j) + 1) |l| {
            const x: Pair = if (l == i) .{ 1, 0 } else w.at(.a, i + l * n);
            const y = w.at(.a, l + j * n);
            s = cadd(s, cmul(x, y));
            mag += cmag(x) * cmag(y);
        }
        return .{ cdist(s, w.at(.orig, row + j * n)), mag };
    }
};

const Getrs = struct {
    pub const name = "getrs";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .dominant;
    pub const prep: Prep = .getrf;
    pub const restore: Restore = .{ .b = true };
    pub const out: Buf = .b;
    pub fn flops(n: f64) f64 {
        return 2 * n * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, In, In, [*]p.F, In, [*]c_int, [*]p.F, In, *c_int, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        f(&fc.N, &n, &n, w.a.ptr, &n, w.ipiv.ptr, w.b.ptr, &n, &w.info, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try lapack.getrs(rstd_alloc, rstd_io, Num(p), .n, n, n, cast(p, w.a.ptr), n, w.rpiv.ptr, cast(p, w.b.ptr), n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return solveAt(p, w, i, j);
    }
};

const Gesv = struct {
    pub const name = "gesv";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .dominant;
    pub const restore: Restore = .{ .a = true, .b = true };
    pub const out: Buf = .b;
    pub const cols: Cols = .one;
    pub fn flops(n: f64) f64 {
        return 2 * n * n * n / 3 + 2 * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (In, In, [*]p.F, In, [*]c_int, [*]p.F, In, *c_int) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        f(&n, &one_c, w.a.ptr, &n, w.ipiv.ptr, w.b.ptr, &n, &w.info);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try lapack.gesv(rstd_alloc, rstd_io, Num(p), n, 1, cast(p, w.a.ptr), n, w.rpiv.ptr, cast(p, w.b.ptr), n, null);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return solveAt(p, w, i, j);
    }
};

const Trtri = struct {
    pub const name = "trtri";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .dominant;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return n * n * n / 3;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, Ch, In, [*]p.F, In, *c_int, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        f(&fc.U, &fc.N, &n, w.a.ptr, &n, &w.info, 1, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        try lapack.trtri(rstd_alloc, rstd_io, Num(p), .u, .n, w.in(), cast(p, w.a.ptr), w.in(), null);
    }
    /// `A X` against `I`, upper triangle.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        if (i > j) return null;
        return inverseAt(p, w, i, j, i, j + 1, struct {
            fn x(v: *const Work(p), k: usize, c: usize) Pair {
                return v.at(.a, k + c * v.n);
            }
        }.x);
    }
};

const Potri = struct {
    pub const name = "potri";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .spd;
    pub const prep: Prep = .potrf;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return 2 * n * n * n / 3;
    }
    pub const Sig = SigUplo;
    pub const lib = libUplo;
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        try lapack.potri(rstd_alloc, rstd_io, Num(p), .u, w.in(), cast(p, w.a.ptr), w.in(), null);
    }
    /// `A X` against `I`, `X` Hermitian from its upper triangle.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return inverseAt(p, w, i, j, 0, w.n, struct {
            fn x(v: *const Work(p), k: usize, c: usize) Pair {
                return symAt(p, v, .a, true, k, c);
            }
        }.x);
    }
};

/// `orig` against `U diag(s) Vᵀ`.
fn svdAt(comptime p: Prec, w: *const Work(p), i: usize, j: usize) [2]f64 {
    const n = w.n;
    var s: Pair = .{ 0, 0 };
    var mag: f64 = 0;
    for (0..n) |k| {
        const u = w.at(.u, i + k * n);
        const v = w.at(.vt, k + j * n);
        const sk = w.at(.s, k)[0];
        const x = cmul(u, v);
        s = cadd(s, .{ x[0] * sk, x[1] * sk });
        mag += cmag(u) * @abs(sk) * cmag(v);
    }
    return .{ cdist(s, w.at(.orig, i + j * n)), mag };
}

/// T34 (faer): thin `U` and `Vᵀ` (`S`).
const Gesvd = struct {
    pub const name = "gesvd";
    pub const shape: Shape = .mm;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .s;
    pub fn flops(n: f64) f64 {
        return 22 * n * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        const F = p.F;
        return if (p.complex)
            fn (Ch, Ch, In, In, [*]F, In, [*]F, [*]F, In, [*]F, In, [*]F, In, [*]F, *c_int, usize, usize) callconv(.c) void
        else
            fn (Ch, Ch, In, In, [*]F, In, [*]F, [*]F, In, [*]F, In, [*]F, In, *c_int, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lw = w.lw();
        if (p.complex)
            f(&fc.S, &fc.S, &n, &n, w.a.ptr, &n, w.s.ptr, w.u.ptr, &n, w.vt.ptr, &n, w.work.ptr, &lw, w.rwork.ptr, &w.info, 1, 1)
        else
            f(&fc.S, &fc.S, &n, &n, w.a.ptr, &n, w.s.ptr, w.u.ptr, &n, w.vt.ptr, &n, w.work.ptr, &lw, &w.info, 1, 1);
    }
    pub const has_rstd = @hasDecl(lapack, "gesvd");
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try lapack.gesvd(rstd_alloc, rstd_io, Num(p), .s, .s, n, n, cast(p, w.a.ptr), n, w.s.ptr, cast(p, w.u.ptr), n, cast(p, w.vt.ptr), n, null);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return svdAt(p, w, i, j);
    }
};

/// T34 (faer): `jobz` `S`.
const Gesdd = struct {
    pub const name = "gesdd";
    pub const shape: Shape = .mm;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .s;
    pub fn flops(n: f64) f64 {
        return 22 * n * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        const F = p.F;
        return if (p.complex)
            fn (Ch, In, In, [*]F, In, [*]F, [*]F, In, [*]F, In, [*]F, In, [*]F, [*]c_int, *c_int, usize) callconv(.c) void
        else
            fn (Ch, In, In, [*]F, In, [*]F, [*]F, In, [*]F, In, [*]F, In, [*]c_int, *c_int, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lw = w.lw();
        if (p.complex)
            f(&fc.S, &n, &n, w.a.ptr, &n, w.s.ptr, w.u.ptr, &n, w.vt.ptr, &n, w.work.ptr, &lw, w.rwork.ptr, w.iwork.ptr, &w.info, 1)
        else
            f(&fc.S, &n, &n, w.a.ptr, &n, w.s.ptr, w.u.ptr, &n, w.vt.ptr, &n, w.work.ptr, &lw, w.iwork.ptr, &w.info, 1);
    }
    pub const has_rstd = @hasDecl(lapack, "gesdd");
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try lapack.gesdd(rstd_alloc, rstd_io, Num(p), .s, n, n, cast(p, w.a.ptr), n, w.s.ptr, cast(p, w.u.ptr), n, cast(p, w.vt.ptr), n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return svdAt(p, w, i, j);
    }
};

/// T111 (faer): add `rstd` when its `geqrf` lands; `tau` in `s`.
const Geqrf = struct {
    pub const name = "geqrf";
    pub const shape: Shape = .mm;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return 4 * n * n * n / 3;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (In, In, [*]p.F, In, [*]p.F, [*]p.F, In, *c_int) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lw = w.lw();
        f(&n, &n, w.a.ptr, &n, w.s.ptr, w.work.ptr, &lw, &w.info);
    }
};

/// T111 (faer): `syevd`, `heevd` complex, with vectors; add `rstd` when its
/// eigensolver lands.
const Syevd = struct {
    pub const name = "syevd";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .spd;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .s;
    pub fn flops(n: f64) f64 {
        return 9 * n * n * n;
    }
    pub fn sym(comptime p: Prec) []const u8 {
        return p.name ++ (if (p.complex) "heevd" else "syevd");
    }
    pub fn Sig(comptime p: Prec) type {
        const F = p.F;
        return if (p.complex)
            fn (Ch, Ch, In, [*]F, In, [*]F, [*]F, In, [*]F, In, [*]c_int, In, *c_int, usize, usize) callconv(.c) void
        else
            fn (Ch, Ch, In, [*]F, In, [*]F, [*]F, In, [*]c_int, In, *c_int, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lw = w.lw();
        const liw = w.liw();
        if (p.complex) {
            const lrw = w.lrw();
            f(&fc.V, &fc.U, &n, w.a.ptr, &n, w.s.ptr, w.work.ptr, &lw, w.rwork.ptr, &lrw, w.iwork.ptr, &liw, &w.info, 1, 1);
        } else f(&fc.V, &fc.U, &n, w.a.ptr, &n, w.s.ptr, w.work.ptr, &lw, w.iwork.ptr, &liw, &w.info, 1, 1);
    }
    /// `A Z` against `Z diag(w)`.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const n = w.n;
        var s: Pair = .{ 0, 0 };
        var mag: f64 = 0;
        for (0..n) |k| {
            const a = w.at(.orig, i + k * n);
            const z = w.at(.a, k + j * n);
            s = cadd(s, cmul(a, z));
            mag += cmag(a) * cmag(z);
        }
        const ev = w.at(.s, j)[0];
        const z = w.at(.a, i + j * n);
        return .{ cdist(s, .{ z[0] * ev, z[1] * ev }), mag + @abs(ev) * cmag(z) };
    }
};

/// T110 (ReLAPACK): Bunch-Kaufman, symmetric (not Hermitian) complex too;
/// no residual, which would need the pivots' blocks.
const Sytrf = struct {
    pub const name = "sytrf";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .spd;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return n * n * n / 3;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, In, [*]p.F, In, [*]c_int, [*]p.F, In, *c_int, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lw = w.lw();
        f(&fc.U, &n, w.a.ptr, &n, w.ipiv.ptr, w.work.ptr, &lw, &w.info, 1);
    }
    pub const has_rstd = @hasDecl(lapack, "sytrf");
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        try lapack.sytrf(rstd_alloc, rstd_io, Num(p), .u, w.in(), cast(p, w.a.ptr), w.in(), w.rpiv.ptr, null);
    }
};

/// T110 (ReLAPACK): kl = ku = `band_k` in LAPACK's band storage, kl rows
/// for the fill, at `gemv`'s band sizes; no residual.
const Gbtrf = struct {
    pub const name = "gbtrf";
    pub const level: Level = .l2;
    pub const shape: Shape = .band;
    pub const band_rows: usize = 3 * band_k + 1;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return 4 * n * band_k * band_k;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (In, In, In, In, [*]p.F, In, [*]c_int, *c_int) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const ld: c_int = band_rows;
        f(&n, &n, &band_kc, &band_kc, w.a.ptr, &ld, w.ipiv.ptr, &w.info);
    }
    pub const has_rstd = @hasDecl(lapack, "gbtrf");
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        try lapack.gbtrf(rstd_alloc, rstd_io, Num(p), w.in(), w.in(), band_k, band_k, cast(p, w.a.ptr), band_rows, w.rpiv.ptr, null);
    }
};
