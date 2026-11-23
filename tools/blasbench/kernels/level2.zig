//! The kernel race `l2`: `gemv` kernels against OpenBLAS's.

const std = @import("std");
const rstd = @import("rstd");
const bench = @import("../bench.zig");
const bb = @import("../blas_bench.zig");
const kn = @import("../kernels.zig");

const Io = std.Io;
const blas = rstd.math.kernels.blas;
const clock = bb.clock;
const fillRandom = bb.fillRandom;
const Fn = bb.Fn;
const l1_timing = bb.l1_timing;
const Prec = bb.Prec;
const precs = bb.precs;
const Cf = kn.Cf;
const asmFor = kn.asmFor;
const elem = kn.elem;
const part = kn.part;
const rounded = kn.rounded;
const size = kn.size;

// ---- level 2 ---------------------------------------------------------------

pub const gemv_alpha = 1e-3;

pub fn GemvCall(comptime p: Prec, comptime name: []const u8) type {
    return struct {
        f: Fn(name),
        n: isize,
        a: [*]p.F,
        x: [*]p.F,
        y: [*]p.F,
        buf: [*]p.F,
        pub fn call(s: *const @This()) void {
            const F = p.F;
            const al: F = gemv_alpha;
            std.mem.doNotOptimizeAway(if (p.complex)
                s.f(s.n, s.n, 0, al, al, s.a, s.n, s.x, 1, s.y, 1, s.buf)
            else
                s.f(s.n, s.n, 0, al, s.a, s.n, s.x, 1, s.y, 1, s.buf));
        }
        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// Bytes of `A`: in L1, L2, memory.
pub const l2_bytes = [_]usize{ 32 << 10, 4 << 20, 64 << 20 };

pub fn level2(io: Io, gpa: std.mem.Allocator) !void {
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

pub fn gemvCandidate(comptime p: Prec, comptime t: []const u8, comptime U: comptime_int) Fn(p.name ++ "gemv_" ++ t) {
    const K = if (p.complex) blas.vector.Complex(p.F, U) else blas.vector.Real(p.F, U);
    return &@field(K, "gemv_" ++ t);
}

pub fn raceGemv(io: Io, gpa: std.mem.Allocator, comptime p: Prec, comptime t: []const u8) !void {
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
    fillRandom(F, a, rand);
    fillRandom(F, x, rand);
    fillRandom(F, y0, rand);

    const asm_k = asmFor(name, true);
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
        for ([_]Fn(name){ @field(bb.tables.portable, name).?, gemvCandidate(p, t, 1), gemvCandidate(p, t, 4) }, [_][]const u8{ "portable", "vec", "vec x4" }) |f, l| {
            fns[count] = f;
            labels[count] = l;
            count += 1;
        }
        var ctxs: [4]C = undefined;
        var calls: [4]bench.Call = undefined;
        var ok = true;
        for (0..count) |i| {
            ctxs[i] = .{ .f = fns[i], .n = @intCast(n), .a = a.ptr, .x = x.ptr, .y = y.ptr, .buf = buf.ptr };
            calls[i] = .of(&ctxs[i]);
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
pub fn checkGemv(comptime p: Prec, trans: bool, n: usize, a: []const p.F, x: []const p.F, y: []const p.F, y0: []const p.F) f64 {
    const eps: f64 = std.math.floatEps(p.F);
    const cs = p.cs();
    const g = rounded(p, gemv_alpha);
    const al: Cf = .init(g, if (p.complex) g else 0);
    var worst: f64 = 0;
    for (0..n) |i| {
        var sum: Cf = .init(0, 0);
        var mag: f64 = 0;
        for (0..n) |l| {
            // A(i, l) or A(l, i)
            const at = if (trans) (l + i * n) * cs else (i + l * n) * cs;
            const av = elem(p, a, at / cs, 1);
            const xv = elem(p, x, l, 1);
            sum = sum.add(av.mul(xv));
            mag += size(av) * size(xv);
        }
        const want = elem(p, y0, i, 1).add(al.mul(sum));
        for (0..cs) |q| {
            const got: f64 = y[i * cs + q];
            const scale = @abs(y0[i * cs + q]) + size(al) * mag;
            worst = @max(worst, @abs(got - part(want, q)) / ((@as(f64, @floatFromInt(n)) + 8) * eps * scale));
        }
    }
    for (n * cs..y.len) |i| {
        if (y[i] != y0[i]) worst = std.math.inf(f64);
    }
    if (std.math.isNan(worst)) return std.math.inf(f64);
    return worst;
}
