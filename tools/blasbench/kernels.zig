//! The kernel races: `l1` (arm64, and x86_64 against OpenBLAS's Haswell and Zen
//! kernels), `l2`, `l3` and `kern3` (arm64).

const std = @import("std");
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const obx86 = @import("obx86.zig");
const tb = @import("libs/table.zig");
const k2 = @import("kernels/level2.zig");
const k3 = @import("kernels/level3.zig");
const kg = @import("kernels/own_gemm.zig");

const Io = std.Io;
const blas = rstd.math.kernels.blas;
const clock = bb.clock;
const Fn = bb.Fn;
const l1_timing = bb.l1_timing;
const Prec = bb.Prec;
const precs = bb.precs;
pub const GemvCall = k2.GemvCall;
pub const checkGemv = k2.checkGemv;
pub const gemvCandidate = k2.gemvCandidate;
pub const gemv_alpha = k2.gemv_alpha;
pub const l2_bytes = k2.l2_bytes;
pub const level2 = k2.level2;
pub const Blocking = k3.Blocking;
pub const GemmCall = k3.GemmCall;
pub const asmGemm = k3.asmGemm;
pub const gemm_sizes = k3.gemm_sizes;
pub const level3 = k3.level3;
pub const checkGemm = kg.checkGemm;
pub const ownLevel3 = kg.ownLevel3;

// ---- level 1 ---------------------------------------------------------------

pub const L1 = enum { axpy, dot, scal, nrm2, asum, sum, iamax, copy, swap, rot };

/// `--routines`: the level-1 routines `l1` and `lib1` race, by name; all when null.
pub var l1_only: ?[]const u8 = null;

/// Whether `--routines` names `r`, or is absent.
pub fn l1Wanted(r: L1) bool {
    var it = std.mem.splitScalar(u8, l1_only orelse return true, ',');
    while (it.next()) |name| {
        if (std.mem.eql(u8, name, @tagName(r))) return true;
    }
    return false;
}

pub fn l1Slot(comptime r: L1, comptime p: Prec) []const u8 {
    return switch (r) {
        .dot => p.name ++ (if (p.complex) "dotu_k" else "dot_k"),
        .iamax => "i" ++ p.name ++ "amax_k",
        .rot => if (p.complex) p.name ++ (if (p.F == f32) "s" else "d") ++ "rot_k" else p.name ++ "rot_k",
        else => p.name ++ @tagName(r) ++ "_k",
    };
}

/// Bytes `r` moves for `n` elements at unit stride.
pub fn l1Bytes(r: L1, n: usize, es: usize) f64 {
    const vectors: f64 = switch (r) {
        .nrm2, .asum, .sum, .iamax => 1,
        .dot, .scal, .copy => 2,
        .axpy => 3,
        .swap, .rot => 4,
    };
    return vectors * @as(f64, @floatFromInt(n * es));
}

pub const axpy_alpha = 1e-6;
pub const scal_alpha = -1.0;
// |alpha| = 1, so repeated calls neither grow nor shrink x
pub const cscal_alpha = .{ 0.6, 0.8 };
pub const rot_cs = .{ 0.6, 0.8 };

/// `v` as the kernels over `p.F` get it.
pub fn rounded(comptime p: Prec, comptime v: comptime_float) f64 {
    return @as(p.F, v);
}

pub fn candidate(comptime r: L1, comptime p: Prec, comptime U: comptime_int) Fn(l1Slot(r, p)) {
    const K = if (p.complex) blas.vector.Complex(p.F, U) else blas.vector.Real(p.F, U);
    return &@field(K, if (r == .dot and p.complex) "dotu" else @tagName(r));
}

pub fn invoke(comptime r: L1, comptime p: Prec, f: Fn(l1Slot(r, p)), n: isize, inc: isize, x: [*]p.F, y: [*]p.F) void {
    const c = p.complex;
    switch (r) {
        .axpy => std.mem.doNotOptimizeAway(if (c) f(n, 0, 0, axpy_alpha, axpy_alpha, x, inc, y, inc, x, 0) else f(n, 0, 0, axpy_alpha, x, inc, y, inc, x, 0)),
        .dot => std.mem.doNotOptimizeAway(f(n, x, inc, y, inc)),
        .scal => std.mem.doNotOptimizeAway(if (c) f(n, 0, 0, cscal_alpha[0], cscal_alpha[1], x, inc, x, 0, x, 1) else f(n, 0, 0, scal_alpha, x, inc, x, 0, x, 1)),
        .nrm2, .asum, .sum, .iamax => std.mem.doNotOptimizeAway(f(n, x, inc)),
        .copy => std.mem.doNotOptimizeAway(f(n, x, inc, y, inc)),
        .swap => std.mem.doNotOptimizeAway(if (c) f(n, 0, 0, 0, 0, x, inc, y, inc, x, 0) else f(n, 0, 0, 0, x, inc, y, inc, x, 0)),
        .rot => std.mem.doNotOptimizeAway(f(n, x, inc, y, inc, rot_cs[0], rot_cs[1])),
    }
}

pub fn L1Call(comptime r: L1, comptime p: Prec) type {
    return struct {
        f: Fn(l1Slot(r, p)),
        n: isize,
        inc: isize,
        x: [*]p.F,
        y: [*]p.F,
        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            invoke(r, p, s.f, s.n, s.inc, s.x, s.y);
        }
    };
}

/// x, y and the values they start from, with room past the longest vector.
pub fn Bufs(comptime F: type) type {
    return struct {
        x: []align(64) F,
        y: []align(64) F,
        x0: []align(64) F,
        y0: []align(64) F,

        pub const guard = 64;

        pub fn init(gpa: std.mem.Allocator, len: usize, seed: u64) !@This() {
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

        pub fn deinit(b: @This(), gpa: std.mem.Allocator) void {
            inline for (.{ "x", "y", "x0", "y0" }) |name| gpa.free(@field(b, name));
        }

        pub fn reset(b: @This()) void {
            @memcpy(b.x, b.x0);
            @memcpy(b.y, b.y0);
        }
    };
}

/// A widened element, real ones with a zero imaginary part.
pub const Cf = std.math.Complex(f64);

/// Element `k` of a vector at stride `inc`, widened.
pub fn elem(comptime p: Prec, v: []const p.F, k: usize, inc: usize) Cf {
    const at = k * inc * p.cs();
    return .init(v[at], if (p.complex) v[at + 1] else 0);
}

/// Part `q` of `c`: 0 its real part, 1 its imaginary.
pub fn part(c: Cf, q: usize) f64 {
    return if (q == 0) c.re else c.im;
}

/// `|re| + |im|`, the size of `c`'s terms.
pub fn size(c: Cf) f64 {
    return @abs(c.re) + @abs(c.im);
}

/// Whether float `i` of a vector of `n` at stride `inc` belongs to it.
pub fn visited(comptime p: Prec, i: usize, n: usize, inc: usize) bool {
    const step = inc * p.cs();
    return i / step < n and i % step < p.cs();
}

/// What reduction `r` should give, (re, im), and the size of its terms.
pub fn l1Scalar(comptime r: L1, comptime p: Prec, b: Bufs(p.F), n: usize, inc: usize) struct { want: [2]f64, scale: f64 } {
    var sum: Cf = .init(0, 0);
    var mag: f64 = 0;
    for (0..n) |k| {
        const xv = elem(p, b.x0, k, inc);
        const yv = elem(p, b.y0, k, inc);
        switch (r) {
            .dot => {
                sum = sum.add(xv.mul(yv));
                mag += size(xv) * size(yv);
            },
            .nrm2 => sum.re += xv.re * xv.re + xv.im * xv.im,
            .asum => sum.re += size(xv),
            .sum => {
                sum.re += xv.re + xv.im;
                mag += size(xv);
            },
            else => comptime unreachable,
        }
    }
    return switch (r) {
        .dot => .{ .want = .{ sum.re, sum.im }, .scale = mag },
        .nrm2 => .{ .want = .{ @sqrt(sum.re), 0 }, .scale = @sqrt(sum.re) },
        .asum => .{ .want = .{ sum.re, 0 }, .scale = sum.re },
        .sum => .{ .want = .{ sum.re, 0 }, .scale = mag },
        else => comptime unreachable,
    };
}

/// Part `q` of element `k` after `r` on fresh data: what x and y should
/// hold, then the size of the terms of each.
pub fn l1Elem(comptime r: L1, comptime p: Prec, b: Bufs(p.F), k: usize, inc: usize, q: usize) [4]f64 {
    const x0 = elem(p, b.x0, k, inc);
    const y0 = elem(p, b.y0, k, inc);
    var wx = part(x0, q);
    var wy = part(y0, q);
    var sx: f64 = @abs(part(x0, q));
    var sy: f64 = @abs(part(y0, q));
    switch (r) {
        .axpy => {
            const al = rounded(p, axpy_alpha);
            const a: Cf = .init(al, if (p.complex) al else 0);
            wy = part(y0, q) + part(a.mul(x0), q);
            sy = @abs(part(y0, q)) + @abs(a.re) * size(x0) + @abs(a.im) * size(x0);
        },
        .scal => {
            const a: Cf = if (p.complex) .init(rounded(p, cscal_alpha[0]), rounded(p, cscal_alpha[1])) else .init(scal_alpha, 0);
            wx = part(a.mul(x0), q);
            sx = size(x0);
        },
        .copy => wy = part(x0, q),
        .swap => {
            wx = part(y0, q);
            wy = part(x0, q);
        },
        .rot => {
            const rc = rounded(p, rot_cs[0]);
            const rs = rounded(p, rot_cs[1]);
            wx = rc * part(x0, q) + rs * part(y0, q);
            wy = rc * part(y0, q) - rs * part(x0, q);
            sx = @abs(part(x0, q)) + @abs(part(y0, q));
            sy = sx;
        },
        else => comptime unreachable,
    }
    return .{ wx, wy, sx, sy };
}

/// Whether `r` returns a value rather than writing its vectors.
pub fn reduces(r: L1) bool {
    return switch (r) {
        .dot, .nrm2, .asum, .sum, .iamax => true,
        .axpy, .scal, .copy, .swap, .rot => false,
    };
}

/// A reduction's result, (re, im), widened; `i?amax`'s index.
pub fn l1Result(comptime r: L1, comptime p: Prec, f: Fn(l1Slot(r, p)), n: usize, inc: usize, b: Bufs(p.F)) [2]f64 {
    const in: isize = @intCast(n);
    const ii: isize = @intCast(inc);
    return switch (r) {
        .dot => if (p.complex) blk: {
            const got = f(in, b.x.ptr, ii, b.y.ptr, ii);
            break :blk .{ got.re, got.im };
        } else .{ f(in, b.x.ptr, ii, b.y.ptr, ii), 0 },
        .nrm2, .asum, .sum => .{ f(in, b.x.ptr, ii), 0 },
        .iamax => .{ @floatFromInt(f(in, b.x.ptr, ii)), 0 },
        else => comptime unreachable,
    };
}

/// The largest error found as a multiple of its tolerance (above 1 is
/// wrong), after one call on fresh data. A routine that writes leaves its
/// output in `b`.
pub fn checkL1(comptime r: L1, comptime p: Prec, f: Fn(l1Slot(r, p)), n: usize, inc: usize, b: Bufs(p.F)) f64 {
    const F = p.F;
    const eps: f64 = std.math.floatEps(F);
    const nf: f64 = @floatFromInt(n);
    b.reset();
    var worst: f64 = 0;
    const span = b.x.len;

    // per element, against a tolerance on the size of its terms
    const Elem = struct {
        pub fn err(got: f64, want: f64, scale: f64, e: f64) f64 {
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
                worst = @max(worst, Elem.err(part(x1, q), w[0], w[2], eps), Elem.err(part(y1, q), w[1], w[3], eps));
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
pub fn firstMax(comptime p: Prec, x: []const p.F, n: usize, inc: usize) isize {
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

pub const Case = struct { bytes: usize, inc: usize };

/// Bytes per vector: in L1, L1, L2, memory; then one strided case in L1.
pub const l1_cases = [_]Case{
    .{ .bytes = 4 << 10, .inc = 1 },
    .{ .bytes = 32 << 10, .inc = 1 },
    .{ .bytes = 1 << 20, .inc = 1 },
    .{ .bytes = 64 << 20, .inc = 1 },
    .{ .bytes = 32 << 10, .inc = 2 },
};

/// The reference kernel for a slot: on arm64 the converted one it holds,
/// vortex's (not with `vortex` false: `sum`, whose vortex slots hold the
/// `@Vector` kernels), else armv8's; on x86_64 OpenBLAS's own for the
/// selected core (`-Dopenblas`); else none.
pub fn asmFor(comptime name: []const u8, comptime vortex: bool) ?struct { f: Fn(name), from: []const u8 } {
    if (!bb.arm) {
        if (!bb.ob_x86) return null;
        const c = obx86.coreFor(blas.cpu.detect()) orelse return null;
        return .{ .f = obx86.kernel(name, c), .from = if (c == .zen) "ob zen" else "ob haswell" };
    }
    const port = @field(bb.tables.portable, name).?;
    if (vortex and @field(bb.tables.vortex, name).? != port) return .{ .f = @field(bb.tables.vortex, name).?, .from = "vortex" };
    if (@field(bb.tables.armv8, name).? != port) return .{ .f = @field(bb.tables.armv8, name).?, .from = "armv8" };
    return null;
}

pub fn level1(io: Io, gpa: std.mem.Allocator) !void {
    tb.sectionMark("l1");
    try clock(io, "level 1 start");
    // x86_64 adds the selected core's own kernel, `core`
    const head = if (bb.arm)
        \\| routine | p | n | inc | asm | asm GB/s | portable | vec | vec x4 | portable/asm | vec/asm | vec x4/asm |
        \\| ------- | - | -: | --: | --- | -------: | -------: | --: | -----: | -----------: | ------: | ---------: |
    else
        \\| routine | p | n | inc | asm | asm GB/s | core | portable | vec | vec x4 | core/asm | portable/asm | vec/asm | vec x4/asm |
        \\| ------- | - | -: | --: | --- | -------: | ---: | -------: | --: | -----: | -------: | -----------: | ------: | ---------: |
    ;
    std.debug.print(
        \\
        \\GB/s, best of {d} batches of at least {d} ms; ratios are each over the reference kernel (`asm`: the converted one on arm64, OpenBLAS's own on x86_64)
        \\
        \\{s}
        \\
    , .{ l1_timing.batches, @divTrunc(l1_timing.min_ns, std.time.ns_per_ms), head });
    inline for (precs) |p| {
        const len = l1_cases[3].bytes / @sizeOf(p.F);
        const b = try Bufs(p.F).init(gpa, len, 1);
        defer b.deinit(gpa);
        inline for (comptime std.meta.tags(L1)) |r| {
            if (l1Wanted(r)) try raceL1(io, r, p, b);
        }
    }
}

pub fn raceL1(io: Io, comptime r: L1, comptime p: Prec, b: Bufs(p.F)) !void {
    const name = comptime l1Slot(r, p);
    const C = L1Call(r, p);
    const asm_k = asmFor(name, r != .sum);
    for (l1_cases) |case| {
        const n = case.bytes / p.es() / case.inc;
        var fns: [5]Fn(name) = undefined;
        var labels: [5][]const u8 = undefined;
        var count: usize = 0;
        if (asm_k) |a| {
            fns[count] = a.f;
            labels[count] = "asm";
            count += 1;
        }
        if (!bb.arm) {
            fns[count] = @field(bb.tables.core, name).?;
            labels[count] = "core";
            count += 1;
        }
        for ([_]Fn(name){ @field(bb.tables.portable, name).?, candidate(r, p, 1), candidate(r, p, 4) }, [_][]const u8{ "portable", "vec", "vec x4" }) |f, l| {
            fns[count] = f;
            labels[count] = l;
            count += 1;
        }
        var ctxs: [5]C = undefined;
        var calls: [5]bench.Call = undefined;
        var ok = true;
        for (0..count) |i| {
            const e = checkL1(r, p, fns[i], n, case.inc, b);
            if (e > 1) {
                std.debug.print("WRONG: {s} {s} {s} n={d} inc={d}: {d:.2} of tolerance\n", .{ @tagName(r), p.name, labels[i], n, case.inc, e });
                ok = false;
            }
            ctxs[i] = .{ .f = fns[i], .n = @intCast(n), .inc = @intCast(case.inc), .x = b.x.ptr, .y = b.y.ptr };
            calls[i] = .of(&ctxs[i]);
        }
        if (!ok) continue;
        b.reset();
        var secs: [5]f64 = undefined;
        try bench.race(io, l1_timing, calls[0..count], &secs);
        const bytes = l1Bytes(r, n, p.es());
        var gbs: [5]f64 = undefined;
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
