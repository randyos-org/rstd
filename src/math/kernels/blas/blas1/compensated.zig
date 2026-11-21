//! Compensated summation (Neumaier's, "Rundungsfehleranalyse einiger
//! Verfahren zur Summation endlicher Summen", ZAMM 54, 1974) for `sum`,
//! `asum` and `dot`: portable code, chunked as the plain sums are, each
//! chunk's sum and compensation folded in chunk order, so the bits depend on
//! neither the worker count nor the order (T138).

const std = @import("std");
const types = @import("../types.zig");
const util = @import("../util.zig");
const threading = @import("../threading.zig");
const core = @import("interface.zig");

const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;

/// A running sum and its compensation, in `W` (`f64` for `f32`).
fn Acc(comptime W: type) type {
    return struct {
        s: W = 0,
        c: W = 0,

        const Self = @This();

        // the error-free sum needs `x` and `t` rounded to `W`, which the
        // x87 holds wider until stored (issue #476)
        inline fn add(a: *Self, x_wide: W) void {
            const x = util.rounded(W, x_wide);
            const t = util.rounded(W, a.s + x);
            if (@abs(a.s) >= @abs(x)) a.c += (a.s - t) + x else a.c += (x - t) + a.s;
            a.s = t;
        }

        fn merge(a: Self, b: Self) Self {
            var r = a;
            r.add(b.s);
            r.c += b.c;
            return r;
        }

        fn value(a: Self) W {
            return a.s + a.c;
        }
    };
}

fn Wide(comptime T: type) type {
    return if (Float(T) == f32) f64 else Float(T);
}

/// The real and imaginary accumulators; `im` unused for real `T`.
fn Pair(comptime T: type) type {
    return struct {
        re: Acc(Wide(T)) = .{},
        im: Acc(Wide(T)) = .{},

        fn merge(a: @This(), b: @This()) @This() {
            return .{ .re = a.re.merge(b.re), .im = a.im.merge(b.im) };
        }
    };
}

const Op = enum { sum, asum, dot, dotc };

/// The vectors from element 0 at their strides, and the op.
fn Ctx(comptime T: type, comptime op: Op) type {
    return struct {
        x: [*]const T,
        incx: isize,
        y: [*]const T,
        incy: isize,

        fn at(p: [*]const T, i: usize, inc: isize) T {
            return util.at(T, @constCast(p), @as(i64, @intCast(i)) * inc).*;
        }

        fn part(c: @This(), ch: threading.Chunk) Pair(T) {
            const W = Wide(T);
            var r: Pair(T) = .{};
            for (ch.start..ch.end) |i| {
                const a = at(c.x, i, c.incx);
                switch (op) {
                    .sum => if (comptime isComplex(T)) {
                        r.re.add(a.re);
                        r.re.add(a.im);
                    } else r.re.add(a),
                    .asum => if (comptime isComplex(T)) {
                        r.re.add(@abs(a.re));
                        r.re.add(@abs(a.im));
                    } else r.re.add(@abs(a)),
                    .dot, .dotc => {
                        const b = at(c.y, i, c.incy);
                        if (comptime isComplex(T)) {
                            const ar: W = a.re;
                            const ai: W = if (op == .dotc) -a.im else a.im;
                            const br: W = b.re;
                            const bi: W = b.im;
                            r.re.add(ar * br);
                            r.re.add(-(ai * bi));
                            r.im.add(ar * bi);
                            r.im.add(ai * br);
                        } else r.re.add(@as(W, a) * @as(W, b));
                    },
                }
            }
            return r;
        }
    };
}

fn reduce(io: Io, comptime T: type, comptime op: Op, comptime family: threading.Family, n: usize, ctx: Ctx(T, op)) Io.Cancelable!Pair(T) {
    const C = Ctx(T, op);
    // f32 and f64 thread, as the plain sums; wider types run in one chunk
    const F = Float(T);
    if ((F != f32 and F != f64) or types.isExtended(F) or threading.serial(family, T, n)) return C.part(ctx, .{ .index = 0, .start = 0, .end = n });
    const p = threading.plan(family, T, n, n, .even, core.quantum);
    return threading.reduce(io, &p, Pair(T), ctx, C.part, Pair(T).merge);
}

/// `sum`, compensated: the sum of `x`, of `re + im` for complex.
pub fn sum(io: Io, comptime T: type, n: usize, x: [*]const T, incx: isize) Io.Cancelable!Float(T) {
    if (n == 0) return 0;
    const r = try reduce(io, T, .sum, .asum, n, .{ .x = x, .incx = incx, .y = x, .incy = incx });
    return @floatCast(r.re.value());
}

/// `asum`, compensated: the sum of `|x|`, of `|re| + |im|` for complex.
pub fn asum(io: Io, comptime T: type, n: usize, x: [*]const T, incx: isize) Io.Cancelable!Float(T) {
    if (n == 0) return 0;
    const r = try reduce(io, T, .asum, .asum, n, .{ .x = x, .incx = incx, .y = x, .incy = incx });
    return @floatCast(r.re.value());
}

/// `dot` (`conj` for `dotc`), compensated, each product rounded once.
pub fn dot(io: Io, comptime T: type, comptime conj: bool, n: usize, x: [*]const T, incx: isize, y: [*]const T, incy: isize) Io.Cancelable!T {
    if (n == 0) return if (comptime isComplex(T)) .init(0, 0) else 0;
    const r = try reduce(io, T, if (conj) .dotc else .dot, .dot, n, .{ .x = x, .incx = incx, .y = y, .incy = incy });
    if (comptime isComplex(T)) return .init(@floatCast(r.re.value()), @floatCast(r.im.value()));
    return @floatCast(r.re.value());
}

// ---- tests ----

const testing = std.testing;
const api = @import("api.zig");
const strided = @import("../strided.zig");
const check = @import("../_test/check.zig");
const test_io = @import("../_test/io.zig");

test "compensated: the sums plain summation loses" {
    const io = check.io;
    const x = [_]f64{ 1e16, 1, -1e16, 1 };
    const ones = [_]f64{ 1, 1, 1, 1 };
    const v: strided.ConstStrided(f64) = .of(&x);
    try testing.expectEqual(@as(f64, 2), try api.sumWith(io, f64, v, .compensated));
    try testing.expectEqual(@as(f64, 2), try api.dotWith(io, f64, v, .of(&ones), .compensated));
    try testing.expectEqual(try api.sum(io, f64, v), try api.sumWith(io, f64, v, .plain));
    const neg = [_]f64{ -1e16, 1, 1e16, 1 };
    try testing.expectEqual(@as(f64, 2e16 + 2), try api.asumWith(io, f64, .of(&neg), .compensated));
    // a negative stride reads the view backward
    const back: strided.ConstStrided(f64) = try .init(&x, 4, -1);
    try testing.expectEqual(@as(f64, 2), try api.dotWith(io, f64, back, try .init(&ones, 4, -1), .compensated));
    // complex: (1 + 2i)(3 - i) + ..., and conj(x) . y
    const C = types.Complex64;
    const cx = [_]C{ .init(1, 2), .init(1e16, 0), .init(-1e16, 0) };
    const cy = [_]C{ .init(3, -1), .init(1, 0), .init(1, 0) };
    const d = try api.dotWith(io, C, .of(&cx), .of(&cy), .compensated);
    try testing.expectEqual(@as(f64, 5), d.re);
    try testing.expectEqual(@as(f64, 5), d.im);
    const dc = try api.dotcWith(io, C, .of(&cx), .of(&cy), .compensated);
    try testing.expectEqual(@as(f64, 1), dc.re);
    try testing.expectEqual(@as(f64, -7), dc.im);
    // f32 sums in f64, rounded once
    const xf = [_]f32{ 1e8, 1, -1e8, 1 };
    try testing.expectEqual(@as(f32, 2), try api.sumWith(io, f32, .of(&xf), .compensated));
}

test "compensated: one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    const al = testing.allocator;
    const saved = threading.test_serial_max;
    threading.test_serial_max = 1;
    defer threading.test_serial_max = saved;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    const n = 1 << 16;
    const x = try al.alloc(f64, n);
    defer al.free(x);
    check.fill(f64, x, 3);
    var got: [test_io.worker_runs.len]f64 = undefined;
    for (&got, 0..) |*g, c| g.* = try api.sumWith(runs.io(c), f64, .of(x), .compensated);
    for (got[1..]) |g| try testing.expectEqual(@as(u64, @bitCast(got[0])), @as(u64, @bitCast(g)));
}
