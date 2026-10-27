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
//! only; `l3` only with the `.S` (`-Dopenblas`).
//!
//! `lib1`, `lib2`, `lib3` race rstd's public routines (`blas.dot`, ...)
//! against real OpenBLAS, BLIS and Accelerate, loaded by absolute path
//! (`results-m6-libraries.md`, `results-flyoff.md`), and give rstd's
//! largest difference from each library's result; on arm64 and x86_64,
//! Windows, Linux and macOS. With the others when none is named and
//! `--openblas` is given. Options: `--openblas=<library>`
//! (`-Dopenblas-lib`), `--blis=<library>`, `--accelerate=<path>` (macOS'
//! framework by default), `--threads=default` (else each library runs one
//! thread), `--ob-core=<core>` (`OPENBLAS_CORETYPE`, which rstd reads too).

const std = @import("std");
const sysinfo = @import("builtin");
const Io = std.Io;
const rstd = @import("rstd");
const bench_options = @import("bench_options");
const bench = @import("bench.zig");

const blas = rstd.linalg.blas;
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

    if (comptime os_tag.isDarwin()) {
        // the performance cores; a `taskpolicy` clamp still overrides it
        if (pthread_set_qos_class_self_np(qos_class_user_interactive, 0) != 0)
            std.debug.print("warning: QoS not set\n", .{});
    }

    tables.portable = blas.cpu.generic.vtableFromCoreType(.generic);
    if (arm) {
        tables.vortex = blas.cpu.vtableFromCoreType(.vortex);
        tables.vortexm4 = blas.cpu.vtableFromCoreType(.vortexm4);
        tables.armv8 = blas.cpu.vtableFromCoreType(.armv8);
        if (tables.vortex.daxpy_k == tables.portable.daxpy_k) return error.VortexNotLinked;
    }

    std.debug.print("zig {s}, {t}-{t}, cpu model {s}, detected core {t}\n", .{ sysinfo.zig_version_string, sysinfo.target.cpu.arch, os_tag, sysinfo.target.cpu.model.name, blas.cpu.detect() });

    // `--name=value` options, then section names
    var opts: LibOptions = .{};
    var sections: std.ArrayList([]const u8) = .empty;
    for (args[1..]) |a| {
        if (std.mem.startsWith(u8, a, "--")) {
            const eq = std.mem.findScalar(u8, a, '=') orelse return error.BadOption;
            const v = a[eq + 1 ..];
            const k = a[2..eq];
            if (std.mem.eql(u8, k, "openblas")) opts.openblas = v else if (std.mem.eql(u8, k, "blis")) opts.blis = v else if (std.mem.eql(u8, k, "accelerate")) opts.accelerate = v else if (std.mem.eql(u8, k, "threads")) opts.threads = v else if (std.mem.eql(u8, k, "ob-core")) opts.ob_core = v else return error.BadOption;
        } else try sections.append(init.arena.allocator(), a);
    }
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
    if (arm) {
        if (all or want(sections.items, "probe")) try probe(io);
        if (all or want(sections.items, "l1")) try level1(io, gpa);
        if (all or want(sections.items, "l2")) try level2(io, gpa);
    }
    if (ob_asm) {
        if (all or want(sections.items, "l3")) try level3(io, gpa);
    }
    const lib_sections = want(sections.items, "lib1") or want(sections.items, "lib2") or want(sections.items, "lib3");
    if (lib_sections or (all and opts.openblas != null)) {
        try loadLibs(opts);
        if (all or want(sections.items, "lib1")) try libLevel1(io, gpa);
        if (all or want(sections.items, "lib2")) try libLevel2(io, gpa);
        if (all or want(sections.items, "lib3")) try libLevel3(io, gpa);
    }
    try clock(io, "end");
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

// ---- libraries -------------------------------------------------------------

const LibOptions = struct {
    openblas: ?[]const u8 = null,
    blis: ?[]const u8 = null,
    accelerate: ?[]const u8 = if (os_tag.isDarwin()) "/System/Library/Frameworks/Accelerate.framework/Accelerate" else null,
    threads: []const u8 = "1",
    ob_core: ?[]const u8 = null,
};

/// Where a library contestant comes from: rstd's public routines, or a
/// library's `cblas_*`.
const Src = enum { rstd, openblas, blis, accelerate };

/// Each library's handle, by `Src` less one; null when not given.
var lib_handles: [3]?Dll = @splat(null);

fn loaded(src: Src) bool {
    return src == .rstd or lib_handles[@backingInt(src) - 1] != null;
}

/// The first error an rstd routine returned inside a kernel-ABI wrapper.
var lib_err: ?RstdError = null;

/// What the raced routines return: every level-1 error is a level-2 one.
const RstdError = blas.Level2Error || blas.Level3Error;

fn fail(e: RstdError) void {
    if (lib_err == null) lib_err = e;
}

/// rstd's routines take an allocator; a C program would hand them malloc's.
const rstd_alloc = std.heap.c_allocator;
/// And an `Io`: a serial one, as the libraries run one thread by default
/// (rstd's routines run on one thread with any `Io` so far).
const rstd_io = std.Io.failing;

extern "c" fn setenv(name: [*:0]const u8, value: [*:0]const u8, overwrite: c_int) c_int;
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

/// Opens the libraries given, one thread each unless `--threads=default`,
/// and prints how each is configured.
fn loadLibs(opts: LibOptions) !void {
    const single = !std.mem.eql(u8, opts.threads, "default");
    // read when the libraries load; Windows' CRT may have copied the
    // environment already, so the thread counts are set by call too
    if (single) {
        for ([_][*:0]const u8{ "OPENBLAS_NUM_THREADS", "BLIS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS" }) |v| setEnv(v, "1");
    }
    var buf: [64]u8 = undefined;
    if (opts.ob_core) |c| setEnv("OPENBLAS_CORETYPE", try std.mem.printSentinel(&buf, "{s}", .{c}, 0));
    const paths = [lib_handles.len]?[]const u8{ opts.openblas orelse return error.NeedOpenblasPath, opts.blis, opts.accelerate };
    for (paths, &lib_handles) |path, *h| {
        const p = path orelse continue;
        h.* = Dll.open(p) catch |e| {
            std.debug.print("cannot load {s}: {t} {s}\n", .{ p, e, dlError() });
            return e;
        };
    }

    const Str = *const fn () callconv(.c) [*:0]const u8;
    const Int = *const fn () callconv(.c) c_int;
    const config = std.mem.span(lookup(.openblas, Str, "openblas_get_config").?());
    // `int` indices here; an ILP64 build (`USE64BITINT`) would need `long`
    if (std.mem.find(u8, config, "USE64BITINT") != null) return error.OpenblasIsIlp64;
    if (single) lookup(.openblas, *const fn (c_int) callconv(.c) void, "openblas_set_num_threads").?(1);
    std.debug.print("openblas: {s}\nopenblas core: {s}, config: {s}, threads: {d}\n", .{
        paths[0].?,
        lookup(.openblas, Str, "openblas_get_corename").?(),
        config,
        lookup(.openblas, Int, "openblas_get_num_threads").?(),
    });

    if (loaded(.blis)) {
        // `dim_t` is 64 bits, `arch_t` an enum
        const Dim = *const fn () callconv(.c) i64;
        if (single) lookup(.blis, *const fn (i64) callconv(.c) void, "bli_thread_set_num_threads").?(1);
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
        // per thread: this one makes every call
        const Set = *const fn (c_uint) callconv(.c) c_int;
        const Get = *const fn () callconv(.c) c_uint;
        if (single) {
            if (lookup(.accelerate, Set, "BLASSetThreading")) |f| {
                if (f(1) != 0) std.debug.print("warning: BLASSetThreading failed\n", .{});
            }
        }
        const threading = if (lookup(.accelerate, Get, "BLASGetThreading")) |f| switch (f()) {
            0 => "multi-threaded",
            1 => "single-threaded",
            else => "?",
        } else "absent";
        std.debug.print("accelerate: {s}, BLASGetThreading: {s}\n", .{ paths[2].?, threading });
    }
    std.debug.print("rstd core: {s}, threads: 1\n", .{blas.coreName()});
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

    fn shown(c: Col) bool {
        return if (c.lib) |l| loaded(l) else c.on;
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
    for (cols) |c| if (c.lib != null and c.shown()) {
        std.debug.print(" rstd-{s} |", .{c.short});
        n += 1;
    };
    std.debug.print(" rstd err | rstd bits |\n| " ++ rule ++ " |", .{});
    for (0..n + 1) |_| std.debug.print(" -: |", .{});
    std.debug.print(" --- |\n", .{});
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
    for (cols, diffs) |c, d| {
        if (c.lib == null or !c.shown()) continue;
        if (d) |v| std.debug.print("{d:.1}{s} | ", .{ v.eps, if (v.eps > v.bound) " !" else "" }) else std.debug.print("- | ", .{});
    }
    std.debug.print("{d:.1} | {x:0>8} |\n", .{ mine.err, @as(u32, @truncate(mine.bits)) });
}

/// rstd's error from the f64 reference, in eps of the type relative to the
/// size of the terms, and a hash of its output.
const Rstd = struct { err: f64 = 0, bits: u64 = 0 };

/// What rstd's own table holds in slot `name`: converted `asm`, or
/// `portable`.
fn rstdRuns(comptime name: []const u8) []const u8 {
    return if (@field(blas.vtable().*, name).? == @field(tables.portable, name).?) "portable" else "asm";
}

const diff_note = "`rstd-<lib>` is rstd's largest difference from that library's result, in eps of the type relative to the size of the terms; `!` marks one past what accumulation order explains: 0 for `copy`, `swap` and `i?amax`, 8 for the other element-wise routines, 2 (n + 8) for a sum of n terms. `rstd err` is rstd's own error from an f64 reference in the same units (gemm: at 256 entries), which tells whose result a difference comes from. `rstd bits` hashes rstd's output, the same on machines that run the same kernels";

fn libLevel1(io: Io, gpa: std.mem.Allocator) !void {
    try clock(io, "library level 1 start");
    std.debug.print(
        \\
        \\GB/s, best of {d} batches of at least {d} ms; ratios are over OpenBLAS. `rstd` is the public routine{s}. {s}
        \\
        \\
    , .{ l1_timing.batches, @divTrunc(l1_timing.min_ns, std.time.ns_per_ms), if (arm) ", `asm` the converted kernel alone (`from`)" else "", diff_note });
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
    const asm_k = if (arm) asmFor(name) else null;
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
    try clock(io, "library level 2 start");
    std.debug.print(
        \\
        \\GFLOP/s, square, unit strides, `beta` 1, best of {d} batches of at least {d} ms; ratios are over OpenBLAS. {s}
        \\
        \\
    , .{ l1_timing.batches, @divTrunc(l1_timing.min_ns, std.time.ns_per_ms), diff_note });
    if (arm)
        printHeader("routine | p | n | rstd runs | from", "------- | - | -: | --------- | ----", &lib_cols)
    else
        printHeader("routine | p | n | rstd runs", "------- | - | -: | ---------", &lib_cols);
    inline for (precs) |p| {
        inline for (.{ "n", "t" }) |t| try raceLibGemv(io, gpa, p, t);
    }
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

    const asm_k = if (arm) asmFor(name) else null;
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
    try clock(io, "library level 3 start");
    std.debug.print(
        \\
        \\gemm NN, square, `C += A*B`: GFLOP/s, best of {d} batches of at least {d} ms; ratios are over OpenBLAS{s}. {s}
        \\
        \\
    , .{ l3_timing.batches, @divTrunc(l3_timing.min_ns, std.time.ns_per_ms), if (ob_asm) ". `ob .S` is OpenBLAS's vortex `.S` through the bench's loops" else "", diff_note });
    printHeader("p | n", "- | -:", &lib3_cols);
    inline for (precs) |p| try raceLibGemm(io, gpa, p);
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
        const m = n * n * cs;
        try gemmScale(p, n, a[0..m], b[0..m], c0[0..m], abs_a, abs_b, scale[0..m]);
        var ctxs: [lib3_cols.len]LC = undefined;
        var calls: [lib3_cols.len]bench.Call = undefined;
        var at: [lib3_cols.len]?usize = @splat(null);
        var diffs: [lib3_cols.len]?Diff = @splat(null);
        var count: usize = 0;
        var mine: Rstd = .{};
        for (0..lib3_cols.len) |o| {
            const i = rstdFirst(o);
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
