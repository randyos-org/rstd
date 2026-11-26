//! The kernel race `l3`: `gemm` kernels and their copies against OpenBLAS's.

const std = @import("std");
const rstd = @import("rstd");
const bench = @import("../bench.zig");
const bb = @import("../blas_bench.zig");

const Io = std.Io;
const blas = rstd.math.kernels.blas;
const Vtable = blas.Vtable;
const clock = bb.clock;
const Fn = bb.Fn;
const l1_timing = bb.l1_timing;
const l3_timing = bb.l3_timing;
const Prec = bb.Prec;
const precs = bb.precs;

// ---- level 3 ---------------------------------------------------------------

extern fn ob_sgemm_kernel(isize, isize, isize, f32, [*]f32, [*]f32, [*]f32, isize) callconv(.c) i32;
extern fn ob_sgemm_itcopy(isize, isize, [*]f32, isize, [*]f32) callconv(.c) i32;
extern fn ob_sgemm_oncopy(isize, isize, [*]f32, isize, [*]f32) callconv(.c) i32;
extern fn ob_dgemm_kernel(isize, isize, isize, f64, [*]f64, [*]f64, [*]f64, isize) callconv(.c) i32;
extern fn ob_dgemm_itcopy(isize, isize, [*]f64, isize, [*]f64) callconv(.c) i32;
extern fn ob_dgemm_oncopy(isize, isize, [*]f64, isize, [*]f64) callconv(.c) i32;

/// A gemm kernel with the copies that pack for it, for `NN` (`level3.c`
/// packs `A` with `itcopy` and `B` with `oncopy` there).
pub fn GemmImpl(comptime p: Prec) type {
    return struct {
        name: []const u8,
        kernel: Fn(p.name ++ "gemm_kernel"),
        itcopy: Fn(p.name ++ "gemm_itcopy"),
        oncopy: Fn(p.name ++ "gemm_oncopy"),
        um: usize,
        un: usize,

        pub fn of(name: []const u8, vt: *const Vtable) @This() {
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
pub fn asmGemm(comptime p: Prec) GemmImpl(p) {
    return if (p.F == f32)
        .{ .name = "asm 16x4", .kernel = &ob_sgemm_kernel, .itcopy = &ob_sgemm_itcopy, .oncopy = &ob_sgemm_oncopy, .um = 16, .un = 4 }
    else
        .{ .name = "asm 8x4", .kernel = &ob_dgemm_kernel, .itcopy = &ob_dgemm_itcopy, .oncopy = &ob_dgemm_oncopy, .um = 8, .un = 4 };
}

/// `level3.c`'s blocking, from the vortex table, with `P` a whole number
/// of the kernel's row panels.
pub const Blocking = struct {
    p: usize,
    q: usize,
    r: usize,
    um: usize,
    un: usize,

    pub fn of(comptime p: Prec, im: GemmImpl(p)) Blocking {
        const pre = p.name ++ "gemm_";
        const bp: usize = @intCast(@field(bb.tables.vortex, pre ++ "p"));
        return .{
            .p = @max(im.um, bp / im.um * im.um),
            .q = @intCast(@field(bb.tables.vortex, pre ++ "q")),
            .r = @intCast(@field(bb.tables.vortex, pre ++ "r")),
            .um = im.um,
            .un = im.un,
        };
    }

    pub fn split(rest: usize, max: usize, unroll: usize) usize {
        if (rest >= 2 * max) return max;
        if (rest > max) return @min(rest, (rest / 2 + unroll - 1) / unroll * unroll);
        return rest;
    }

    pub fn chunk(b: Blocking, rest: usize) usize {
        if (rest >= 3 * b.un) return 3 * b.un;
        if (rest > b.un) return b.un;
        return rest;
    }
};

/// `C += A*B`, square `n`, through `level3.c`'s loops as `blas3/driver.zig`
/// has them, with `im`'s kernel and copies.
pub fn GemmCall(comptime p: Prec) type {
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

        pub fn call(s: *const @This()) void {
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

        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// The kernel alone on panels its own copies packed: `m` by `n`, depth `k`.
pub fn KernelCall(comptime p: Prec) type {
    const F = p.F;
    return struct {
        im: GemmImpl(p),
        m: usize,
        n: usize,
        k: usize,
        sa: [*]F,
        sb: [*]F,
        c: [*]F,

        pub fn call(s: *const @This()) void {
            _ = s.im.kernel(@intCast(s.m), @intCast(s.n), @intCast(s.k), 1, s.sa, s.sb, s.c, @intCast(s.m));
        }

        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// `c0 + A*B` against `c` at `samples` entries (all when 0), as `checkL1`.
pub fn checkGemmReal(comptime F: type, m: usize, n: usize, k: usize, a: []const F, b: []const F, c: []const F, c0: []const F, samples: usize) f64 {
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

pub const gemm_sizes = [_]usize{ 64, 128, 256, 512, 1024, 2048 };

pub fn level3(io: Io, gpa: std.mem.Allocator) !void {
    try clock(io, "level 3 start");
    inline for (precs[0..2]) |p| try raceGemm(io, gpa, p);
    std.debug.print("\ncgemm, zgemm: no Zig kernel (only real gemm is ported), not raced\n", .{});
}

pub fn raceGemm(io: Io, gpa: std.mem.Allocator, comptime p: Prec) !void {
    const F = p.F;
    const impls = [_]GemmImpl(p){
        asmGemm(p),
        GemmImpl(p).of("zig vortex", &bb.tables.vortex),
        GemmImpl(p).of("zig vortexm4", &bb.tables.vortexm4),
        GemmImpl(p).of("zig portable", &bb.tables.portable),
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
            calls[i] = .of(&ctxs[i]);
            @memcpy(c[0 .. km * kn], c0[0 .. km * kn]);
            ctxs[i].call();
            const e = checkGemmReal(F, km, kn, kk, a[0 .. km * kk], b[0 .. kk * kn], c[0 .. km * kn], c0[0 .. km * kn], 0);
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
            calls[i] = .of(&ctxs[i]);
            @memcpy(c[0 .. n * n], c0[0 .. n * n]);
            ctxs[i].call();
            const e = checkGemmReal(F, n, n, n, a[0 .. n * n], b[0 .. n * n], c[0 .. n * n], c0[0 .. n * n], 256);
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
