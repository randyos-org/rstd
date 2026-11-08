//! rstd's BLAS and LAPACK raced against OpenBLAS, BLIS, Accelerate and the reference binaries, and its arm64 kernels against
//! OpenBLAS's own `.S`, printing Markdown table rows. Sections, options and steps: `readme.md`.

const std = @import("std");
const sysinfo = @import("builtin");
const Io = std.Io;
const rstd = @import("rstd");
const bench_options = @import("bench_options");
const bench = @import("bench.zig");
const merge = @import("merge.zig");
const kernels = @import("kernels.zig");
const libs = @import("libs.zig");
const ext = @import("ext.zig");

const blas = rstd.linalg.blas;
const Vtable = blas.Vtable;
const os_tag = sysinfo.target.os.tag;

const level1 = kernels.level1;
const level2 = kernels.level2;
const level3 = kernels.level3;
const ownLevel3 = kernels.ownLevel3;

const libLevel1 = libs.libLevel1;
const libLevel2 = libs.libLevel2;
const libLevel3 = libs.libLevel3;
const LibOptions = libs.LibOptions;
const loaded = libs.loaded;
const loadLibs = libs.loadLibs;
const merge_labels = libs.merge_labels;
const RstdIo = libs.RstdIo;
const setThreadEnv = libs.setThreadEnv;
const unsetenv = libs.unsetenv;

const libLapack = ext.libLapack;

/// The arm64 cores' own races: `probe`, `l1`, `l2`, and the `asm` and
/// `vec x4` library columns.
pub const arm = sysinfo.target.cpu.arch == .aarch64;
/// OpenBLAS's vortex gemm `.S` is linked (`-Dopenblas`): `l3`, `ob .S`.
pub const ob_asm = bench_options.ob_asm;

comptime {
    if (!arm and sysinfo.target.cpu.arch != .x86_64) @compileError("blas_bench runs on arm64 and x86_64");
    if (ob_asm and !arm) @compileError("the vortex `.S` is arm64");
}

/// The kernel type of vtable slot `name`.
pub fn Fn(comptime name: []const u8) type {
    return @typeInfo(@FieldType(Vtable, name)).optional.child;
}

const Tables = struct {
    vortex: Vtable,
    vortexm4: Vtable,
    armv8: Vtable,
    portable: Vtable,
};
pub var tables: Tables = undefined;

pub const Prec = struct {
    F: type,
    complex: bool,
    name: []const u8,

    pub fn cs(comptime p: Prec) usize {
        return if (p.complex) 2 else 1;
    }
    /// Bytes of one element.
    pub fn es(comptime p: Prec) usize {
        return p.cs() * @sizeOf(p.F);
    }
};

pub const precs = [_]Prec{
    .{ .F = f32, .complex = false, .name = "s" },
    .{ .F = f64, .complex = false, .name = "d" },
    .{ .F = f32, .complex = true, .name = "c" },
    .{ .F = f64, .complex = true, .name = "z" },
};

/// `v` uniform in [-1, 1), from `rand`.
pub fn fillRandom(comptime F: type, v: []F, rand: std.Random) void {
    for (v) |*e| e.* = @floatCast(rand.float(f64) * 2 - 1);
}

pub const l1_timing: bench.Timing = .{ .min_ns = 10 * std.time.ns_per_ms, .batches = 15 };
pub const l3_timing: bench.Timing = .{ .min_ns = 20 * std.time.ns_per_ms, .batches = 7 };

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
                .@"max-n" => libs.n_max = std.fmt.parseInt(usize, v, 10) catch return error.BadOption,
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
    if (opts.rstd_io == .threaded) libs.rstd_io = threaded.io();
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
pub var alone = false;

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

pub fn clock(io: Io, what: []const u8) !void {
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
