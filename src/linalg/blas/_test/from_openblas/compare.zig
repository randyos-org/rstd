//! OpenBLAS's comparison programs, `test/compare_sgemv_{sbgemv,bgemv,shgemv}.c`
//! and `test/compare_sgemm_{sbgemm,bgemm,shgemm}.c` with `test_helpers.h`'s
//! `is_close`, v0.3.34-317: each half routine against `sgemv` or `sgemm` on
//! the unrounded values and against a plain loop over the rounded ones,
//! with the sources' tolerances. Random operands come from a seeded
//! generator, not `rand()`; conversions go through `sbtobf16` and `bf16to`,
//! as the C's `sbstobf16_` and `sbf16tos_`.
//!
//! The sizes are the sources' (every order to 256 for `gemv`; to 100 and
//! 256 for `gemm`, and `shgemm`'s scalar sweep on 3, 24, 55, 71, 128, 256);
//! `small` (the rigs, and every core) stops at 24. An order of 0 passes
//! leading dimensions of 1: rstd checks them against `max(1, rows)`, as the
//! reference BLAS does, where OpenBLAS lets 0 through.

const std = @import("std");
const ob = @import("__root__.zig");
const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");

const BF16 = types.BF16;
const HF16 = types.HF16;
const al = std.testing.allocator;

/// `test_helpers.h`'s `is_close`: `fabs(a - b) <= atol + rtol*fabs(b)`.
fn isClose(a: f32, b: f32, rtol: f64, atol: f64) bool {
    return @abs(@as(f64, a - b)) <= atol + rtol * @abs(@as(f64, b));
}

/// `sbstobf16_` of one value.
fn toBf(v: f32) !BF16 {
    var in = v;
    var out: BF16 = undefined;
    try blas.sbtobf16(ob.io, f32, 1, @ptrCast(&in), 1, @ptrCast(&out), 1);
    return out;
}

/// `float16to32`: `sbf16tos_` of one value.
fn fromBf(v: BF16) !f32 {
    var in = v;
    var out: f32 = undefined;
    try blas.bf16to(ob.io, f32, 1, @ptrCast(&in), 1, @ptrCast(&out), 1);
    return out;
}

/// `truncate_float32_to_bfloat16`: `v` through bfloat16 and back.
fn truncate(v: f32) !f32 {
    return fromBf(try toBf(v));
}

/// `v` in the half type `E`: `sbstobf16_` or a `_Float16` cast.
fn toHalf(comptime E: type, v: f32) !E {
    return if (E == BF16) toBf(v) else @floatCast(v);
}

/// A half value as `f32`: `float16to32` or a `float` cast.
fn fromHalf(comptime E: type, v: E) !f32 {
    return if (E == BF16) fromBf(v) else @floatCast(v);
}

/// `n` values at stride `inc` into `E`, as `toHalf` does each.
fn toHalves(comptime E: type, n: usize, src: []const f32, inc: usize, dst: []E) !void {
    if (E == BF16) return blas.sbtobf16(ob.io, f32, @intCast(n), @constCast(src.ptr), @intCast(inc), dst.ptr, @intCast(inc));
    for (0..n) |i| dst[i * inc] = @floatCast(src[i * inc]);
}

/// `n` values at stride `inc` out of `E`, as `fromHalf` does each.
fn fromHalves(comptime E: type, n: usize, src: []const E, inc: usize, dst: []f32) !void {
    if (E == BF16) return blas.bf16to(ob.io, f32, @intCast(n), @constCast(src.ptr), @intCast(inc), dst.ptr, @intCast(inc));
    for (0..n) |i| dst[i * inc] = @floatCast(src[i * inc]);
}

/// `((FLOAT) rand () / (FLOAT) RAND_MAX) + 0.5`.
fn rnd(r: std.Random) f32 {
    return r.float(f32) + 0.5;
}

fn alloc(comptime T: type, n: usize) ![]T {
    return al.alloc(T, @max(n, 1));
}

/// Counts failures, as the C's `ret`.
const Tally = struct {
    name: []const u8,
    failures: usize = 0,

    fn check(t: *Tally, ok: bool, comptime f: []const u8, args: anytype) void {
        if (ok) return;
        if (t.failures < 5) ob.say(f, args);
        t.failures += 1;
    }

    fn done(t: Tally) !void {
        if (t.failures == 0) return;
        ob.say("{s} FAILURES: {d}\n", .{ t.name, t.failures });
        return error.TestUnexpectedResult;
    }
};

// ---- gemv ----------------------------------------------------------------------

/// `compare_sgemv_<kind>.c`.
fn gemvAgainst(comptime kind: blas.HalfGemv, small: bool) !void {
    const In = if (kind == .shgemv) HF16 else BF16;
    const Out = if (kind == .bgemv) BF16 else f32;
    const loop: usize = if (small) 24 else 256;
    var prng = std.Random.DefaultPrng.init(std.hash.Wyhash.hash(0x5be3, @tagName(kind)));
    const r = prng.random();
    var t: Tally = .{ .name = @tagName(kind) };
    var beta: f32 = 0;
    while (beta < 3) : (beta += 1) {
        var alpha: f32 = 0;
        while (alpha < 3) : (alpha += 1) {
            for (0..2) |l| for (1..loop + 1) |x| {
                // `l` 1 tests strides other than one
                const inc = l + 1;
                const a = try alloc(f32, x * x);
                defer al.free(a);
                const b = try alloc(f32, x << @intCast(l));
                defer al.free(b);
                const c = try alloc(f32, x << @intCast(l));
                defer al.free(c);
                const aa = try alloc(In, x * x);
                defer al.free(aa);
                const bb = try alloc(In, x << @intCast(l));
                defer al.free(bb);
                const cc = try alloc(Out, x << @intCast(l));
                defer al.free(cc);
                const dd = try alloc(f32, x);
                defer al.free(dd);
                // `float16to32` of every element, once
                const wa = try alloc(f32, x * x);
                defer al.free(wa);
                const wb = try alloc(f32, x << @intCast(l));
                defer al.free(wb);
                for (0..x) |j| {
                    for (0..x) |i| a[j * x + i] = rnd(r);
                    b[j << @intCast(l)] = rnd(r);
                }
                try toHalves(In, x * x, a, 1, aa);
                try toHalves(In, x, b, inc, bb);
                try fromHalves(In, x * x, aa, 1, wa);
                try fromHalves(In, x, bb, inc, wb);
                for ([_]blas.Trans{ .n, .t }) |trans| {
                    @memset(c, 0);
                    @memset(cc, try toHalf(Out, 0));
                    @memset(dd, 0);
                    const xi: isize = @intCast(x);
                    const ii: isize = @intCast(inc);
                    try blas.gemv(al, ob.io, f32, trans, xi, xi, alpha, a.ptr, xi, b.ptr, ii, beta, c.ptr, ii);
                    try blas.gemvHalf(ob.io, kind, trans, xi, xi, try toHalf(Out, alpha), aa.ptr, xi, bb.ptr, ii, try toHalf(Out, beta), cc.ptr, ii);
                    for (dd) |*v| v.* *= beta;
                    for (0..x) |j| for (0..x) |i| {
                        if (trans == .n) {
                            dd[i] += alpha * wa[j * x + i] * wb[j << @intCast(l)];
                        } else {
                            dd[j] += alpha * wa[j * x + i] * wb[i << @intCast(l)];
                        }
                    };
                    for (0..x) |j| {
                        const got = try fromHalf(Out, cc[j << @intCast(l)]);
                        // `bgemv` against the others rounded to bfloat16
                        const vs_c = if (Out == BF16) try truncate(c[j << @intCast(l)]) else c[j << @intCast(l)];
                        const vs_d = if (Out == BF16) try truncate(dd[j]) else dd[j];
                        t.check(isClose(got, vs_c, 0.01, 0.001), "{s} {s} x={d} inc={d} alpha={d} beta={d} y[{d}]: {d}, sgemv {d}\n", .{ @tagName(kind), @tagName(trans), x, inc, alpha, beta, j, got, vs_c });
                        t.check(isClose(got, vs_d, 0.001, 0.0001), "{s} {s} x={d} inc={d} alpha={d} beta={d} y[{d}]: {d}, loop {d}\n", .{ @tagName(kind), @tagName(trans), x, inc, alpha, beta, j, got, vs_d });
                    }
                }
            };
        }
    }
    try t.done();
}

// ---- gemm ----------------------------------------------------------------------

/// One problem of `compare_sgemm_<kind>.c`: `m = n = k = x`, `A` and `B`
/// random, both in `f32` and rounded to `In`.
fn Operands(comptime In: type) type {
    return struct {
        const Self = @This();
        a: []f32,
        b: []f32,
        aa: []In,
        bb: []In,
        /// `float16to32` of `aa` and `bb`, once.
        wa: []f32,
        wb: []f32,

        fn init(x: usize, r: std.Random) !Self {
            const a = try alloc(f32, x * x);
            errdefer al.free(a);
            const b = try alloc(f32, x * x);
            errdefer al.free(b);
            const aa = try alloc(In, x * x);
            errdefer al.free(aa);
            const bb = try alloc(In, x * x);
            errdefer al.free(bb);
            const wa = try alloc(f32, x * x);
            errdefer al.free(wa);
            const wb = try alloc(f32, x * x);
            errdefer al.free(wb);
            for (a[0 .. x * x]) |*v| v.* = rnd(r);
            for (b[0 .. x * x]) |*v| v.* = rnd(r);
            try toHalves(In, x * x, a, 1, aa);
            try toHalves(In, x * x, b, 1, bb);
            try fromHalves(In, x * x, aa, 1, wa);
            try fromHalves(In, x * x, bb, 1, wb);
            return .{ .a = a, .b = b, .aa = aa, .bb = bb, .wa = wa, .wb = wb };
        }

        fn deinit(s: Self) void {
            inline for (.{ s.a, s.b, s.wa, s.wb }) |p| al.free(p);
            al.free(s.aa);
            al.free(s.bb);
        }

        /// The C's naive `DD[i * m + j] += AA * BB` for column `i`, row `j`.
        fn product(s: Self, x: usize, ta: bool, tb: bool, i: usize, j: usize, l: usize) f32 {
            const av = if (ta) s.wa[x * j + l] else s.wa[l * x + j];
            const bv = if (tb) s.wb[i + l * x] else s.wb[l + x * i];
            return av * bv;
        }
    };
}

const trans_pairs = [_][2]blas.Trans{ .{ .n, .n }, .{ .t, .n }, .{ .n, .t }, .{ .t, .t } };

/// The sources' orders: every one to `dense`, then `largest`.
fn orders(small: bool) [2]usize {
    return if (small) .{ 24, 24 } else .{ 100, 256 };
}

/// `compare_sgemm_sbgemm.c`, and the first part of `compare_sgemm_shgemm.c`,
/// whose loop check is as loose as its `sgemm` one.
fn gemmAgainst(comptime kind: blas.HalfGemm, small: bool) !void {
    const In = kind.In();
    const Out = kind.Out();
    const dense, const largest = orders(small);
    var prng = std.Random.DefaultPrng.init(std.hash.Wyhash.hash(0x5be7, @tagName(kind)));
    const r = prng.random();
    var t: Tally = .{ .name = @tagName(kind) };
    // `bgemm` takes `beta` 1 onto a `C` of 100
    const alpha: f32 = 1;
    const beta: f32 = if (kind == .bgemm) 1 else 0;
    const c0: f32 = if (kind == .bgemm) 100 else 0;
    const loop_rtol: f64, const loop_atol: f64 = if (kind == .sbgemm) .{ 0.001, 0.0001 } else .{ 0.01, 0.001 };
    for (0..largest + 1) |x| {
        if (x > dense and x != largest) continue;
        const ops = try Operands(In).init(x, r);
        defer ops.deinit();
        const c = try alloc(f32, x * x);
        defer al.free(c);
        const cc = try alloc(Out, x * x);
        defer al.free(cc);
        const dd = try alloc(f32, x * x);
        defer al.free(dd);
        const xi: isize = @intCast(x);
        const ld: isize = @max(1, xi);
        for (trans_pairs) |tr| {
            @memset(c, c0);
            @memset(cc, try toHalf(Out, c0));
            @memset(dd, c0);
            try blas.gemm(al, ob.io, f32, tr[0], tr[1], xi, xi, xi, alpha, ops.a.ptr, ld, ops.b.ptr, ld, beta, c.ptr, ld);
            try blas.gemmHalf(al, ob.io, kind, tr[0], tr[1], xi, xi, xi, try toHalf(Out, alpha), ops.aa.ptr, ld, ops.bb.ptr, ld, try toHalf(Out, beta), cc.ptr, ld);
            for (0..x) |i| for (0..x) |j| {
                for (0..x) |l| dd[i * x + j] += ops.product(x, tr[0] == .t, tr[1] == .t, i, j, l);
                const got = try fromHalf(Out, cc[i * x + j]);
                const vs_c = if (Out == BF16) try truncate(c[i * x + j]) else c[i * x + j];
                const vs_d = if (Out == BF16) try truncate(dd[i * x + j]) else dd[i * x + j];
                t.check(isClose(got, vs_c, 0.01, 0.001), "{s} {s}{s} x={d} C[{d},{d}]: {d}, sgemm {d}\n", .{ @tagName(kind), @tagName(tr[0]), @tagName(tr[1]), x, j, i, got, vs_c });
                t.check(isClose(got, vs_d, loop_rtol, loop_atol), "{s} {s}{s} x={d} C[{d},{d}]: {d}, loop {d}\n", .{ @tagName(kind), @tagName(tr[0]), @tagName(tr[1]), x, j, i, got, vs_d });
            };
        }
    }
    try t.done();
}

/// The second part of `compare_sgemm_shgemm.c`: `alpha` and `beta` 0, 0.5
/// and 1 on a zeroed `C`, against `sgemm` alone.
fn shgemmScalars(small: bool) !void {
    const xs: []const usize = if (small) &.{ 3, 24 } else &.{ 3, 24, 55, 71, 128, 256 };
    var prng = std.Random.DefaultPrng.init(0x5be9);
    const r = prng.random();
    var t: Tally = .{ .name = "shgemm" };
    for (xs) |x| {
        var alpha: f32 = 0;
        while (alpha <= 1) : (alpha += 0.5) {
            var beta: f32 = 0;
            while (beta <= 1) : (beta += 0.5) {
                const ops = try Operands(HF16).init(x, r);
                defer ops.deinit();
                const c = try alloc(f32, x * x);
                defer al.free(c);
                const cc = try alloc(f32, x * x);
                defer al.free(cc);
                const xi: isize = @intCast(x);
                for (trans_pairs) |tr| {
                    @memset(c, 0);
                    @memset(cc, 0);
                    try blas.gemm(al, ob.io, f32, tr[0], tr[1], xi, xi, xi, alpha, ops.a.ptr, xi, ops.b.ptr, xi, beta, c.ptr, xi);
                    try blas.gemmHalf(al, ob.io, .shgemm, tr[0], tr[1], xi, xi, xi, alpha, ops.aa.ptr, xi, ops.bb.ptr, xi, beta, cc.ptr, xi);
                    for (c, cc, 0..) |w, g, at| t.check(isClose(g, w, 0.01, 0.001), "shgemm {s}{s} x={d} alpha={d} beta={d} C[{d}]: {d}, sgemm {d}\n", .{ @tagName(tr[0]), @tagName(tr[1]), x, alpha, beta, at, g, w });
                }
            }
        }
    }
    try t.done();
}

/// Every program; `small` takes the rigs' sizes.
pub fn all(small: bool) !void {
    inline for (std.enums.values(blas.HalfGemv)) |kind| try gemvAgainst(kind, small);
    inline for (std.enums.values(blas.HalfGemm)) |kind| try gemmAgainst(kind, small);
    try shgemmScalars(small);
}

test "compare_sgemv_sbgemv, bgemv, shgemv: the half gemvs against sgemv and a plain loop" {
    inline for (std.enums.values(blas.HalfGemv)) |kind| try gemvAgainst(kind, ob.reduced());
}

test "compare_sgemm_sbgemm, bgemm, shgemm: the half gemms against sgemm and a plain loop" {
    inline for (std.enums.values(blas.HalfGemm)) |kind| try gemmAgainst(kind, ob.reduced());
}

test "compare_sgemm_shgemm: shgemm's alpha and beta sweep against sgemm" {
    try shgemmScalars(ob.reduced());
}
