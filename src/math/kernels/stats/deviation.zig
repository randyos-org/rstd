//! `deviationSums`: `Σ(x − m)` and `Σ|x − m|²` in one pass, for the
//! corrected two-pass variance (Chan, Golub and LeVeque, "Algorithms for
//! computing the sample variance", The American Statistician 37, 1983):
//! `(Σd² − |Σd|²/n) / (n − 1)`. Chunked as the BLAS's sums, each chunk's
//! sums folded in chunk order, so the bits depend on neither the worker
//! count nor the order (#605).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const threading = @import("../blas/threading.zig");
const strided = @import("../blas/strided.zig");
const shared = @import("../blas/vector/util.zig");
const api = @import("../blas/blas1/api.zig");
const core = @import("../blas/blas1/interface.zig");
const fast = @import("../blas/blas1/fast.zig");
const compensated = @import("../blas/blas1/compensated.zig");

const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const ConstStrided = strided.ConstStrided;
const Summation = api.Summation;
const Wide = compensated.Wide;

/// `Σ(x − m)` and `Σ|x − m|²`.
pub fn DeviationSums(comptime T: type) type {
    return struct { sum: T, sum_sq: Float(T) };
}

/// `Σ(x − m)` and `Σ|x − m|²` over `x`, summed as `s` says: `.plain` in
/// `f64` for `f32` (`Float(T)` for wider types), `.compensated` with
/// Neumaier's compensation there, `.fast` in `f32` lanes for `f32` and
/// `Complex(f32)` (`fast.zig`'s), else as `.plain`. Taken in memory order,
/// whatever the stride's sign.
pub fn deviationSums(io: Io, comptime T: type, s: Summation, m: T, x: ConstStrided(T)) Io.Cancelable!DeviationSums(T) {
    const n = x.n;
    const inc: isize = if (n <= 1) 1 else @intCast(@abs(x.stride));
    const p: [*]const T = x.data.ptr;
    const r = switch (s) {
        .plain => try reduce(io, T, n, PlainCtx(T, false){ .x = p, .inc = inc, .m = m }),
        .compensated => try reduce(io, T, n, PlainCtx(T, true){ .x = p, .inc = inc, .m = m }),
        .fast => if (comptime fast.serves(T))
            try reduce(io, T, n, FastCtx(T){ .x = p, .inc = inc, .m = m })
        else
            try reduce(io, T, n, PlainCtx(T, false){ .x = p, .inc = inc, .m = m }),
    };
    const sum: T = if (comptime isComplex(T)) .init(@floatCast(r.re), @floatCast(r.im)) else @floatCast(r.re);
    return .{ .sum = sum, .sum_sq = @floatCast(r.sq) };
}

/// A chunk's sums: `Σd` (`re`, `im`) and `Σ|d|²` (`sq`), plain or
/// compensated, in `Wide(T)`.
fn Part(comptime T: type, comptime comp: bool) type {
    const W = Wide(T);
    return struct {
        re: Acc = .{},
        im: Acc = .{},
        sq: Acc = .{},

        const Self = @This();
        const Acc = if (comp) compensated.Acc(W) else Plain;

        const Plain = struct {
            s: W = 0,

            inline fn add(a: *Plain, v: W) void {
                a.s = util.rounded(W, a.s + v);
            }

            fn merge(a: Plain, b: Plain) Plain {
                return .{ .s = util.rounded(W, a.s + b.s) };
            }

            fn value(a: Plain) W {
                return a.s;
            }
        };

        fn merge(a: Self, b: Self) Self {
            return .{ .re = a.re.merge(b.re), .im = a.im.merge(b.im), .sq = a.sq.merge(b.sq) };
        }
    };
}

/// The folded result, in `Wide(T)`.
fn Result(comptime T: type) type {
    const W = Wide(T);
    return struct { re: W, im: W, sq: W };
}

/// `x` from element 0 at stride `inc`, and the shift `m`; `.plain` and
/// `.compensated` element by element.
fn PlainCtx(comptime T: type, comptime comp: bool) type {
    return struct {
        x: [*]const T,
        inc: isize,
        m: T,

        pub const P = Part(T, comp);

        pub fn part(c: @This(), ch: threading.Chunk) P {
            const W = Wide(T);
            var r: P = .{};
            for (ch.start..ch.end) |i| {
                const a = util.at(T, @constCast(c.x), @as(i64, @intCast(i)) * c.inc).*;
                if (comptime isComplex(T)) {
                    const dr = util.rounded(W, @as(W, a.re) - @as(W, c.m.re));
                    const di = util.rounded(W, @as(W, a.im) - @as(W, c.m.im));
                    r.re.add(dr);
                    r.im.add(di);
                    r.sq.add(dr * dr);
                    r.sq.add(di * di);
                } else {
                    const d = util.rounded(W, @as(W, a) - @as(W, c.m));
                    r.re.add(d);
                    r.sq.add(d * d);
                }
            }
            return r;
        }
    };
}

/// The deviations' terms for `fast.laneSums`: `d = x − m`, then `d` and
/// `d²` (fused where an instruction does it), `m`'s parts alternating over
/// a complex element's floats.
fn Deviations(comptime complex: bool) type {
    return struct {
        re: f32,
        im: f32,

        pub const k = 2;

        inline fn shift(d: @This(), comptime Vt: type) Vt {
            if (!complex) return @splat(d.re);
            const len = @typeInfo(Vt).vector.len;
            const even = std.simd.iota(u32, len) % @as(@Vector(len, u32), @splat(2)) == @as(@Vector(len, u32), @splat(0));
            return @select(f32, even, @as(Vt, @splat(d.re)), @as(Vt, @splat(d.im)));
        }

        pub inline fn add(d: @This(), comptime Vt: type, acc: *[k]Vt, v: Vt) void {
            const e = v - d.shift(Vt);
            acc[0] += e;
            acc[1] = shared.mulAdd(Vt, e, e, acc[1]);
        }

        pub inline fn addOne(d: @This(), s: anytype, lane: usize, v: f32) void {
            const e = v - (if (lane == 0) d.re else d.im);
            s[0][lane] += e;
            s[1][lane] = shared.mulAdd(f32, e, e, s[1][lane]);
        }
    };
}

/// `x` and `m` for `.fast`, in `f32` lanes; each chunk's sums widened to
/// `f64` for the fold.
fn FastCtx(comptime T: type) type {
    return struct {
        x: [*]const T,
        inc: isize,
        m: T,

        const complex = isComplex(T);
        const w = if (complex) 2 else 1;
        pub const P = Part(T, false);

        pub fn part(c: @This(), ch: threading.Chunk) P {
            const n = ch.end - ch.start;
            const x: [*]const f32 = @ptrCast(util.at(T, @constCast(c.x), @as(i64, @intCast(ch.start)) * c.inc));
            const map: Deviations(complex) = if (complex) .{ .re = c.m.re, .im = c.m.im } else .{ .re = c.m, .im = 0 };
            const s = if (c.inc == 1) fast.laneSums(w, map, w * n, x) else fast.stridedSums(w, w, map, n, x, c.inc);
            var r: P = .{};
            r.re.s = s[0][0];
            if (complex) r.im.s = s[0][1];
            r.sq.s = if (complex) @as(f64, s[1][0]) + @as(f64, s[1][1]) else s[1][0];
            return r;
        }
    };
}

fn reduce(io: Io, comptime T: type, n: usize, ctx: anytype) Io.Cancelable!Result(T) {
    const C = @TypeOf(ctx);
    const P = C.P;
    // f32 and f64 thread, as the BLAS's sums; wider types run in one chunk
    const F = Float(T);
    const one = (F != f32 and F != f64) or types.isExtended(F) or threading.serial(.dot, T, n);
    const r: P = if (n == 0) .{} else if (one) C.part(ctx, .{ .index = 0, .start = 0, .end = n }) else blk: {
        const p = threading.plan(.dot, T, n, n, .even, core.quantum);
        break :blk try threading.reduce(io, &p, P, ctx, C.part, P.merge);
    };
    return .{ .re = r.re.value(), .im = r.im.value(), .sq = r.sq.value() };
}
