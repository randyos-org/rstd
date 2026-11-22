//! `Summation.fast`'s portable part: the real `dot` of `f32` added in `f32`
//! lanes, and the lane sums `kernels.stats` builds on. `asum`, `sum` and
//! complex `dot` take the core's own kernels instead (issue #584). Portable
//! `@Vector` code for the build's target, the same on every core.
//! Chunked as the plain sums are, the chunks' sums folded in chunk order in
//! `f64`, so the bits depend on neither the worker count nor the order (T198).

const std = @import("std");
const types = @import("../types.zig");
const util = @import("../util.zig");
const threading = @import("../threading.zig");
const shared = @import("../vector/util.zig");
const core = @import("interface.zig");

const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const mulAdd = shared.mulAdd;

/// Whether `.fast` adds `T` in its own lanes: `f32` and its `Complex`, where
/// `f32` lanes keep denormals (`util.vectorIeee`) and are not the x87's.
/// Elsewhere `.fast` is `.plain`.
pub fn serves(comptime T: type) bool {
    return Float(T) == f32 and util.vectorIeee(f32) and !util.onX87(f32);
}

const F = f32;
/// Lanes of a native vector.
const L = shared.lanes(F);
const N = @Vector(L, F);
/// Floats a step in each of two sums: four native vectors.
pub const W = 4 * L;
const V = @Vector(W, F);
/// Four strided floats.
const G = @Vector(4, F);

inline fn load(comptime T: type, p: [*]const F) T {
    return p[0..@typeInfo(T).vector.len].*;
}

/// Float `c` of element `i`, elements `w` floats and `inc` elements apart.
inline fn flt(comptime w: comptime_int, x: [*]const F, inc: isize, i: usize, comptime c: comptime_int) F {
    return util.at(F, @constCast(x), @as(i64, @intCast(i)) * inc * w + c).*;
}

/// The four floats of elements `i` on (four real ones, two complex).
inline fn gather(comptime w: comptime_int, x: [*]const F, inc: isize, i: usize) G {
    var g: [4]F = undefined;
    inline for (0..4) |k| g[k] = flt(w, x, inc, i + k / w, k % w);
    return g;
}

/// `v`'s upper half added to its lower until `n` lanes are left.
fn halved(comptime n: comptime_int, v: anytype) @Vector(n, F) {
    const m = @typeInfo(@TypeOf(v)).vector.len;
    if (m == n) return v;
    const h = m / 2;
    const lo = @shuffle(F, v, undefined, std.simd.iota(i32, h));
    const hi = @shuffle(F, v, undefined, std.simd.iota(i32, h) + @as(@Vector(h, i32), @splat(h)));
    return halved(n, lo + hi);
}

/// `map`'s `k` sums of `len` floats at `x`, `x` at an even float, so a lane's
/// parity is its float's (a complex element's `re` or `im`). `map.add(Vt,
/// acc, v)` adds a vector's terms to `acc`, `map.addOne(s, lane, f)` one
/// float's to lane `lane` of `s`. Float `i` of each `2W` into lane `i % W` of
/// accumulator `i / W % 2`, a last `W` into the first, then native vectors
/// into a third; the first two added and halved to `L` lanes, the third
/// added, halved to `keep` lanes, and the floats left added in order to lane
/// `i % keep`.
pub fn laneSums(comptime keep: comptime_int, map: anytype, len: usize, x: [*]const F) [@TypeOf(map).k][keep]F {
    const k = @TypeOf(map).k;
    var acc: [2][k]V = @splat(@splat(@splat(0)));
    var i: usize = 0;
    while (i + 2 * W <= len) : (i += 2 * W) {
        inline for (0..2) |q| map.add(V, &acc[q], load(V, x + i + q * W));
    }
    if (i + W <= len) {
        map.add(V, &acc[0], load(V, x + i));
        i += W;
    }
    var t: [k]N = @splat(@splat(0));
    while (i + L <= len) : (i += L) map.add(N, &t, load(N, x + i));
    var s: [k][keep]F = undefined;
    inline for (0..k) |j| s[j] = halved(keep, halved(L, acc[0][j] + acc[1][j]) + t[j]);
    while (i < len) : (i += 1) map.addOne(&s, i % keep, x[i]);
    return s;
}

/// `laneSums` of `n` strided elements of `w` floats: float `i` of each four
/// into lane `i % 4`; the lanes reduced (one sum) or halved to two, then the
/// rest in order, float `c` of an element to lane `c % keep`.
pub fn stridedSums(comptime keep: comptime_int, comptime w: comptime_int, map: anytype, n: usize, x: [*]const F, inc: isize) [@TypeOf(map).k][keep]F {
    const k = @TypeOf(map).k;
    var acc: [k]G = @splat(@splat(0));
    var i: usize = 0;
    while (i + 4 / w <= n) : (i += 4 / w) map.add(G, &acc, gather(w, x, inc, i));
    var s: [k][keep]F = undefined;
    inline for (0..k) |j| s[j] = if (keep == 1) .{@reduce(.Add, acc[j])} else halved(keep, acc[j]);
    while (i < n) : (i += 1) {
        inline for (0..w) |c| map.addOne(&s, c % keep, flt(w, x, inc, i, c));
    }
    return s;
}

/// The real dot of `len` floats, in `floats`' order, each product fused
/// where an instruction does it.
fn dotFloats(len: usize, x: [*]const F, y: [*]const F) F {
    var acc: [2]V = @splat(@splat(0));
    var i: usize = 0;
    while (i + 2 * W <= len) : (i += 2 * W) {
        inline for (0..2) |q| acc[q] = mulAdd(V, load(V, x + i + q * W), load(V, y + i + q * W), acc[q]);
    }
    if (i + W <= len) {
        acc[0] = mulAdd(V, load(V, x + i), load(V, y + i), acc[0]);
        i += W;
    }
    var t: N = @splat(0);
    while (i + L <= len) : (i += L) t = mulAdd(N, load(N, x + i), load(N, y + i), t);
    var s = halved(1, halved(L, acc[0] + acc[1]) + t)[0];
    while (i < len) : (i += 1) s = mulAdd(F, x[i], y[i], s);
    return s;
}

fn dotStrided(n: usize, x: [*]const F, incx: isize, y: [*]const F, incy: isize) F {
    var acc: G = @splat(0);
    var i: usize = 0;
    while (i + 4 <= n) : (i += 4) acc = mulAdd(G, gather(1, x, incx, i), gather(1, y, incy, i), acc);
    var s = @reduce(.Add, acc);
    while (i < n) : (i += 1) s = mulAdd(F, flt(1, x, incx, i, 0), flt(1, y, incy, i, 0), s);
    return s;
}

/// The vectors from element 0 at their strides.
const Ctx = struct {
    x: [*]const F,
    incx: isize,
    y: [*]const F,
    incy: isize,

    fn from(p: [*]const F, i: usize, inc: isize) [*]const F {
        return @ptrCast(util.at(F, @constCast(p), @as(i64, @intCast(i)) * inc));
    }

    /// A chunk's dot, widened for the fold.
    fn part(c: Ctx, ch: threading.Chunk) f64 {
        const n = ch.end - ch.start;
        const x = from(c.x, ch.start, c.incx);
        const y = from(c.y, ch.start, c.incy);
        return if (c.incx == 1 and c.incy == 1) dotFloats(n, x, y) else dotStrided(n, x, c.incx, y, c.incy);
    }
};

// rounded to `f64`, so the bits do not hang on an x87 register
fn merge(a: f64, b: f64) f64 {
    return util.rounded(f64, a + b);
}

/// The real `dot` of `f32` in `f32` lanes. Complex `.fast` and `asum`
/// and `sum` take the core's own kernels (`interface.zig`, issue #584).
pub fn dot(io: Io, n: usize, x: [*]const F, incx: isize, y: [*]const F, incy: isize) Io.Cancelable!F {
    comptime std.debug.assert(serves(F));
    if (n == 0) return 0;
    const ctx: Ctx = .{ .x = x, .incx = incx, .y = y, .incy = incy };
    if (threading.serial(.dot, F, n)) return @floatCast(ctx.part(.{ .index = 0, .start = 0, .end = n }));
    const p = threading.plan(.dot, F, n, n, .even, core.quantum);
    return @floatCast(try threading.reduce(io, &p, f64, ctx, Ctx.part, merge));
}

// ---- tests ----

/// The ops the tests take through the API's `.fast`.
const Op = enum { sum, asum, dot, dotc };

const testing = std.testing;
const api = @import("api.zig");
const strided = @import("../strided.zig");
const check = @import("../_test/check.zig");
const test_io = @import("../_test/io.zig");

const ConstStrided = strided.ConstStrided;
const Complex32 = types.Complex32;

/// A plain loop's result in `f64`, and the sum of its terms' magnitudes.
const Want = struct { re: f64 = 0, im: f64 = 0, mag: f64 = 0 };

fn want(comptime T: type, comptime op: Op, x: ConstStrided(T), y: ConstStrided(T)) Want {
    var r: Want = .{};
    for (0..x.n) |i| {
        const a = x.get(i);
        if (comptime !isComplex(T)) {
            const v: f64 = a;
            const t = switch (op) {
                .sum => v,
                .asum => @abs(v),
                .dot, .dotc => v * @as(f64, y.get(i)),
            };
            r.re += t;
            r.mag += @abs(t);
            continue;
        }
        const ar: f64 = a.re;
        const ai: f64 = if (op == .dotc) -a.im else a.im;
        switch (op) {
            .sum, .asum => {
                r.re += if (op == .sum) ar + ai else @abs(ar) + @abs(ai);
                r.mag += @abs(ar) + @abs(ai);
            },
            .dot, .dotc => {
                const b = y.get(i);
                const br: f64 = b.re;
                const bi: f64 = b.im;
                r.re += ar * br - ai * bi;
                r.im += ar * bi + ai * br;
                r.mag += @abs(ar * br) + @abs(ai * bi) + @abs(ar * bi) + @abs(ai * br);
            },
        }
    }
    return r;
}

/// `op` summed `.fast` through the API, as `re` and `im`.
fn got(comptime T: type, comptime op: Op, io: Io, x: ConstStrided(T), y: ConstStrided(T)) ![2]f64 {
    switch (op) {
        .sum => return .{ try api.sumWith(io, T, x, .fast), 0 },
        .asum => return .{ try api.asumWith(io, T, x, .fast), 0 },
        .dot, .dotc => {
            const d = if (op == .dotc) try api.dotcWith(io, T, x, y, .fast) else try api.dotWith(io, T, x, y, .fast);
            return if (comptime isComplex(T)) .{ d.re, d.im } else .{ d, 0 };
        },
    }
}

/// `.fast`'s error bound over `n` elements, in `f32` epsilons of the terms'
/// magnitudes: the longest chain of additions a term takes (`len / 4`
/// strided, fewer at unit strides; then the lanes and the rest) and a
/// product's rounding.
fn bound(comptime T: type, n: usize, mag: f64) f64 {
    const len = n * @as(usize, if (isComplex(T)) 2 else 1);
    return @as(f64, @floatFromInt(len / 4 + 2 * W + 10)) * std.math.floatEps(f32) * mag;
}

const ops = .{ Op.sum, Op.asum, Op.dot, Op.dotc };

test "fast: within the f32 lanes' error bound of a plain loop" {
    const al = testing.allocator;
    const n_max = 20000;
    const sizes = [_]usize{ 1, 2, 3, 7, W - 1, W, 2 * W - 1, 2 * W, 2 * W + 1, 100, 1000, 4097, n_max };
    const strides = [_][2]isize{ .{ 1, 1 }, .{ 2, 1 }, .{ 1, -3 }, .{ -2, -1 } };
    inline for (.{ f32, Complex32 }) |T| {
        const bx = try al.alloc(T, 3 * n_max);
        defer al.free(bx);
        const by = try al.alloc(T, 3 * n_max);
        defer al.free(by);
        check.fill(T, bx, 7);
        check.fill(T, by, 8);
        for (sizes) |n| for (strides) |s| {
            const x = try ConstStrided(T).init(bx[0 .. (n - 1) * @abs(s[0]) + 1], n, s[0]);
            const y = try ConstStrided(T).init(by[0 .. (n - 1) * @abs(s[1]) + 1], n, s[1]);
            inline for (ops) |op| {
                if (comptime op == .dotc and !isComplex(T)) continue;
                const w = want(T, op, x, y);
                const g = try got(T, op, check.io, x, y);
                const b = bound(T, n, w.mag);
                if (!(@abs(g[0] - w.re) <= b and @abs(g[1] - w.im) <= b)) {
                    std.debug.print("fast {t} {s} n={d} strides={any}: {any}, want {d} {d}, bound {e}\n", .{ op, @typeName(T), n, s, g, w.re, w.im, b });
                    return error.TestUnexpectedResult;
                }
            }
        };
    }
}

test "fast: f32 lanes drop what f64 keeps" {
    if (comptime !serves(f32)) return error.SkipZigTest;
    const n = 2 * W * 64;
    const tiny = std.math.ldexp(@as(f32, 1), -24);
    var x: [n]f32 = @splat(tiny);
    x[0] = 1;
    // lane 0 of the first sum drops the 63 others it meets after its 1
    const lanes: f32 = 1 + @as(f32, n - 64) * tiny;
    try testing.expectEqual(lanes, try api.sumWith(check.io, f32, .of(&x), .fast));
    try testing.expect(lanes < @as(f32, @floatCast(1 + @as(f64, n - 1) * tiny)));
}

test "fast: f64 sums as plain" {
    var x: [3000]f64 = undefined;
    var y: [3000]f64 = undefined;
    check.fill(f64, &x, 1);
    check.fill(f64, &y, 2);
    const io = check.io;
    try testing.expectEqual(try api.sumWith(io, f64, .of(&x), .plain), try api.sumWith(io, f64, .of(&x), .fast));
    try testing.expectEqual(try api.asumWith(io, f64, .of(&x), .plain), try api.asumWith(io, f64, .of(&x), .fast));
    try testing.expectEqual(try api.dotWith(io, f64, .of(&x), .of(&y), .plain), try api.dotWith(io, f64, .of(&x), .of(&y), .fast));
    var cx: [1500]types.Complex64 = undefined;
    var cy: [1500]types.Complex64 = undefined;
    check.fill(types.Complex64, &cx, 3);
    check.fill(types.Complex64, &cy, 4);
    const p = try api.dotcWith(io, types.Complex64, .of(&cx), .of(&cy), .plain);
    const f = try api.dotcWith(io, types.Complex64, .of(&cx), .of(&cy), .fast);
    try testing.expectEqual(p.re, f.re);
    try testing.expectEqual(p.im, f.im);
}

test "fast: one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    const n = (1 << 16) + 37;
    inline for (.{ f32, Complex32 }) |T| {
        const x = try al.alloc(T, 3 * n);
        defer al.free(x);
        const y = try al.alloc(T, n);
        defer al.free(y);
        check.fillMixed(T, x, 5);
        check.fillMixed(T, y, 6);
        const views = [_][2]ConstStrided(T){ .{ .of(x[0..n]), .of(y) }, .{ try .init(x, n, 3), try .init(y, n, -1) } };
        for (views) |v| inline for (ops) |op| {
            if (comptime op == .dotc and !isComplex(T)) continue;
            var outs: [test_io.worker_runs.len][2]f64 = undefined;
            for (&outs, 0..) |*o, c| {
                const io = runs.io(c);
                // the same chunks at every count: only the workers differ
                threading.test_serial_max = 4096;
                o.* = try got(T, op, io, v[0], v[1]);
            }
            for (outs[1..]) |o| for (0..2) |k| {
                try testing.expectEqual(@as(u64, @bitCast(outs[0][k])), @as(u64, @bitCast(o[k])));
            };
        };
    }
}
