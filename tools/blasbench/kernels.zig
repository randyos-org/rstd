//! The arm64 kernel races: `l1`, `l2`, `l3` and `kern3`.

const std = @import("std");
const Io = std.Io;
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const libs = @import("libs.zig");

const blas = rstd.linalg.blas;
const Vtable = blas.Vtable;

const clock = bb.clock;
const fillRandom = bb.fillRandom;
const Fn = bb.Fn;
const l1_timing = bb.l1_timing;
const l3_timing = bb.l3_timing;
const Prec = bb.Prec;
const precs = bb.precs;

const cast = libs.cast;
const fail = libs.fail;
const Num = libs.Num;
const rstd_alloc = libs.rstd_alloc;

// ---- level 1 ---------------------------------------------------------------

pub const L1 = enum { axpy, dot, scal, nrm2, asum, iamax, copy, swap, rot };

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
        .nrm2, .asum, .iamax => 1,
        .dot, .scal, .copy => 2,
        .axpy => 3,
        .swap, .rot => 4,
    };
    return vectors * @as(f64, @floatFromInt(n * es));
}

pub const axpy_alpha = 1e-6;
const scal_alpha = -1.0;
// |alpha| = 1, so repeated calls neither grow nor shrink x
const cscal_alpha = .{ 0.6, 0.8 };
const rot_cs = .{ 0.6, 0.8 };

/// `v` as the kernels over `p.F` get it.
pub fn rounded(comptime p: Prec, comptime v: comptime_float) f64 {
    return @as(p.F, v);
}

pub fn candidate(comptime r: L1, comptime p: Prec, comptime U: comptime_int) Fn(l1Slot(r, p)) {
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
const Cf = std.math.Complex(f64);

/// Element `k` of a vector at stride `inc`, widened.
fn elem(comptime p: Prec, v: []const p.F, k: usize, inc: usize) Cf {
    const at = k * inc * p.cs();
    return .init(v[at], if (p.complex) v[at + 1] else 0);
}

/// Part `q` of `c`: 0 its real part, 1 its imaginary.
fn part(c: Cf, q: usize) f64 {
    return if (q == 0) c.re else c.im;
}

/// `|re| + |im|`, the size of `c`'s terms.
fn size(c: Cf) f64 {
    return @abs(c.re) + @abs(c.im);
}

/// Whether float `i` of a vector of `n` at stride `inc` belongs to it.
fn visited(comptime p: Prec, i: usize, n: usize, inc: usize) bool {
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
            else => comptime unreachable,
        }
    }
    return switch (r) {
        .dot => .{ .want = .{ sum.re, sum.im }, .scale = mag },
        .nrm2 => .{ .want = .{ @sqrt(sum.re), 0 }, .scale = @sqrt(sum.re) },
        .asum => .{ .want = .{ sum.re, 0 }, .scale = sum.re },
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
        .dot, .nrm2, .asum, .iamax => true,
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
        .nrm2, .asum => .{ f(in, b.x.ptr, ii), 0 },
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
pub const l1_cases = [_]Case{
    .{ .bytes = 4 << 10, .inc = 1 },
    .{ .bytes = 32 << 10, .inc = 1 },
    .{ .bytes = 1 << 20, .inc = 1 },
    .{ .bytes = 64 << 20, .inc = 1 },
    .{ .bytes = 32 << 10, .inc = 2 },
};

/// Which converted kernel a slot holds: vortex's own, else armv8's, else none.
pub fn asmFor(comptime name: []const u8) ?struct { f: Fn(name), from: []const u8 } {
    const port = @field(bb.tables.portable, name).?;
    if (@field(bb.tables.vortex, name).? != port) return .{ .f = @field(bb.tables.vortex, name).?, .from = "vortex" };
    if (@field(bb.tables.armv8, name).? != port) return .{ .f = @field(bb.tables.armv8, name).?, .from = "armv8" };
    return null;
}

pub fn level1(io: Io, gpa: std.mem.Allocator) !void {
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
        for ([_]Fn(name){ @field(bb.tables.portable, name).?, candidate(r, p, 1), candidate(r, p, 4) }, [_][]const u8{ "portable", "vec", "vec x4" }) |f, l| {
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
            calls[i] = .of(&ctxs[i]);
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
    fillRandom(F, a, rand);
    fillRandom(F, x, rand);
    fillRandom(F, y0, rand);

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

        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// `c0 + A*B` against `c` at `samples` entries (all when 0), as `checkL1`.
fn checkGemmReal(comptime F: type, m: usize, n: usize, k: usize, a: []const F, b: []const F, c: []const F, c0: []const F, samples: usize) f64 {
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

fn raceGemm(io: Io, gpa: std.mem.Allocator, comptime p: Prec) !void {
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

// ---- rstd's own gemm --------------------------------------------------------

/// `kern3`: the selected core's gemm kernel alone, on a `P` by 64 block of
/// depth `Q` its own copies packed, against rstd's whole `gemm` on one
/// worker: how much of the kernel's rate the driver keeps (T21).
pub fn ownLevel3(io: Io, gpa: std.mem.Allocator) !void {
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

        pub fn run(o: *anyopaque) void {
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
