//! `kern3`: the selected core's `gemm` kernel alone against rstd's whole `gemm`.

const std = @import("std");
const rstd = @import("rstd");
const bench = @import("../bench.zig");
const bb = @import("../blas_bench.zig");
const libs = @import("../libs.zig");
const k3 = @import("level3.zig");

const Io = std.Io;
const blas = rstd.math.kernels.blas;
const clock = bb.clock;
const Fn = bb.Fn;
const l1_timing = bb.l1_timing;
const l3_timing = bb.l3_timing;
const Prec = bb.Prec;
const precs = bb.precs;
const cast = libs.cast;
const fail = libs.fail;
const Num = libs.Num;
const rstd_alloc = libs.rstd_alloc;
const checkGemmReal = k3.checkGemmReal;

/// `kern3`: the selected core's gemm kernel alone, on a `P` by 64 block of
/// depth `Q` its own copies packed, against rstd's whole `gemm` on one
/// worker: how much of the kernel's rate the driver keeps (T21).
pub fn ownLevel3(io: Io, gpa: std.mem.Allocator) !void {
    try clock(io, "own level 3 start");
    inline for (precs) |p| try ownGemm(io, gpa, p);
}

/// The selected core's kernel alone on panels its own copies packed: `m` by
/// `n`, depth `k`.
pub fn OwnKernel(comptime p: Prec) type {
    const F = p.F;
    return struct {
        kernel: Fn(p.name ++ if (p.complex) "gemm_kernel_n" else "gemm_kernel"),
        m: usize,
        n: usize,
        k: usize,
        sa: [*]F,
        sb: [*]F,
        c: [*]F,

        pub fn call(s: *const @This()) void {
            if (p.complex)
                _ = s.kernel(@intCast(s.m), @intCast(s.n), @intCast(s.k), 1, 0, s.sa, s.sb, s.c, @intCast(s.m))
            else
                _ = s.kernel(@intCast(s.m), @intCast(s.n), @intCast(s.k), 1, s.sa, s.sb, s.c, @intCast(s.m));
        }

        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// rstd's public `gemm`, `C += A*B`, square `n`, on `Io.failing` (one worker).
pub fn OwnGemm(comptime p: Prec) type {
    const F = p.F;
    return struct {
        n: usize,
        a: [*]F,
        b: [*]F,
        c: [*]F,

        pub fn call(s: *const @This()) void {
            const one: Num(p) = if (p.complex) .init(1, 0) else 1;
            const sq = s.n * s.n;
            blas.gemm(rstd_alloc, .failing, Num(p), .n, .n, s.n, s.n, s.n, one, cast(p, s.a)[0..sq], s.n, cast(p, s.b)[0..sq], s.n, one, cast(p, s.c)[0..sq], s.n) catch |e| fail(e);
        }

        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// A gemm result checked against an f64 reference, for real or interleaved complex `p`.
pub fn checkGemm(comptime p: Prec, m: usize, n: usize, k: usize, a: []const p.F, b: []const p.F, c: []const p.F, c0: []const p.F, samples: usize) f64 {
    if (!p.complex) return checkGemmReal(p.F, m, n, k, a, b, c, c0, samples);
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

pub const own_sizes = [_]usize{ 128, 256, 512, 1024 };

pub fn ownGemm(io: Io, gpa: std.mem.Allocator, comptime p: Prec) !void {
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
    const ke = checkGemm(p, km, kn, kk, a[0 .. cs * km * kk], b[0 .. cs * kk * kn], c[0 .. cs * km * kn], c0[0 .. cs * km * kn], 256);
    if (ke > 1) {
        std.debug.print("WRONG: {s}gemm kernel: {d:.2} of tolerance\n", .{ p.name, ke });
        return error.WrongKernel;
    }
    var ksecs: [1]f64 = undefined;
    try bench.race(io, l1_timing, &.{.of(&kc)}, &ksecs);
    const krate = flop * @as(f64, @floatFromInt(km * kn * kk)) / ksecs[0] / 1e9;
    std.debug.print("| kernel alone, {d} x {d}, depth {d} | {d:.1} | 1.00 |\n", .{ km, kn, kk, krate });

    for (own_sizes) |n| {
        var gc: OwnGemm(p) = .{ .n = n, .a = a.ptr, .b = b.ptr, .c = c.ptr };
        @memcpy(c[0 .. cs * n * n], c0[0 .. cs * n * n]);
        gc.call();
        const e = checkGemm(p, n, n, n, a[0 .. cs * n * n], b[0 .. cs * n * n], c[0 .. cs * n * n], c0[0 .. cs * n * n], 256);
        if (e > 1) {
            std.debug.print("WRONG: {s}gemm n={d}: {d:.2} of tolerance\n", .{ p.name, n, e });
            return error.WrongGemm;
        }
        var secs: [1]f64 = undefined;
        try bench.race(io, l3_timing, &.{.of(&gc)}, &secs);
        const nf: f64 = @floatFromInt(n);
        const rate = flop * nf * nf * nf / secs[0] / 1e9;
        std.debug.print("| gemm, n = {d} | {d:.1} | {d:.2} |\n", .{ n, rate, rate / krate });
    }
}
