//! Special values in random operands: the patterns a case puts in them, and
//! the comparison of a result against the wide reference where they appear.

const std = @import("std");
const types = @import("../types.zig");
const util = @import("../util.zig");
const check = @import("check.zig");
const reference = @import("reference.zig");

const isComplex = types.isComplex;
const FloatOf = types.Float;
const Wide = reference.Wide;

pub const Special = enum { nan, inf, ninf };

/// What a case puts in its random operands.
pub const Pattern = union(enum) {
    /// One special value at element `at`, in its real or imaginary part by
    /// `at`'s parity.
    one: struct { kind: Special, at: usize },
    /// By index mod 8: -0, +0, a denormal, a huge value, +Inf, -Inf, NaN, and
    /// one left random.
    mixed,
    /// Every element denormal.
    denormal,
    /// Every element near the square root of the smallest normal: squares
    /// underflow.
    tiny,
    /// Every element huge, though no product with a value of at most 2, or
    /// sum of 256 of them, overflows; squares do.
    huge,

    pub fn describe(p: Pattern) void {
        switch (p) {
            .one => |o| check.say("{s} at {d}: ", .{ @tagName(o.kind), o.at }),
            else => check.say("{s}: ", .{@tagName(p)}),
        }
    }
};

/// The patterns for a vector of `n`: each special value at the first, middle
/// and last element, then the whole-vector ones.
pub fn patterns(n: usize, buf: *[13]Pattern) []const Pattern {
    var len: usize = 0;
    for (std.enums.values(Special)) |kind| for ([_]usize{ 0, n / 2, n - 1 }, 0..) |at, k| {
        if ((k > 0 and at == 0) or (k == 2 and at == n / 2)) continue;
        buf[len] = .{ .one = .{ .kind = kind, .at = at } };
        len += 1;
    };
    for ([_]Pattern{ .mixed, .denormal, .tiny, .huge }) |p| {
        buf[len] = p;
        len += 1;
    }
    return buf[0..len];
}

/// A case's patterns for an operand of `len`: `patterns` natively, four on
/// the rigs, where emulation is slow.
pub fn patternsFor(len: usize, buf: *[13]Pattern) []const Pattern {
    if (!check.rig()) return patterns(len, buf);
    buf[0..4].* = .{ .{ .one = .{ .kind = .nan, .at = len / 2 } }, .{ .one = .{ .kind = .inf, .at = len - 1 } }, .mixed, .denormal };
    return buf[0..4];
}

pub fn special(comptime F: type, kind: Special) F {
    return switch (kind) {
        .nan => std.math.nan(F),
        .inf => std.math.inf(F),
        .ninf => -std.math.inf(F),
    };
}

/// The factors that take a value of at most 2 into each range.
fn denormalScale(comptime F: type) F {
    return comptime std.math.ldexp(@as(F, 1), std.math.floatExponentMin(F) - 4);
}
fn tinyScale(comptime F: type) F {
    return comptime std.math.ldexp(@as(F, 1), @divFloor(std.math.floatExponentMin(F), 2) - 4);
}
fn hugeScale(comptime F: type) F {
    return comptime std.math.ldexp(@as(F, 1), std.math.floatExponentMax(F) - 12);
}

/// Scalar where AltiVec flushes (`util.scalar`), so a denormal stays one.
fn scaled(comptime T: type, v: T, f: FloatOf(T)) T {
    const F = FloatOf(T);
    return if (comptime isComplex(T)) .init(util.scalar(F, v.re * f), util.scalar(F, v.im * f)) else util.scalar(T, v * f);
}

/// `f` in the real part of `v` for an even `i`, else the imaginary part.
fn put(comptime T: type, v: T, i: usize, f: FloatOf(T)) T {
    if (comptime !isComplex(T)) return f;
    return if (i % 2 == 0) .init(f, v.im) else .init(v.re, f);
}

/// Random element `v` at index `i` of an operand, shaped by `p`.
pub fn shape(comptime T: type, v: T, i: usize, p: Pattern) T {
    const F = FloatOf(T);
    return switch (p) {
        .one => |o| if (i == o.at) put(T, v, i, special(F, o.kind)) else v,
        .mixed => switch (i % 8) {
            0, 1 => scaled(T, v, if (i % 8 == 0) -0.0 else 0),
            2 => scaled(T, v, denormalScale(F)),
            3 => scaled(T, v, hugeScale(F)),
            4 => put(T, v, i / 8, std.math.inf(F)),
            5 => put(T, v, i / 8, -std.math.inf(F)),
            6 => put(T, v, i / 8, std.math.nan(F)),
            else => v,
        },
        .denormal => scaled(T, v, denormalScale(F)),
        .tiny => scaled(T, v, tinyScale(F)),
        .huge => scaled(T, v, hugeScale(F)),
    };
}

/// `got` against the reference `want` where special values may appear: NaN
/// for NaN, the infinity of `want`'s sign where it is one or past `T`'s
/// range, else within `ulps` units of `scale` plus `ulps` of `T`'s smallest
/// denormal, so a flushed denormal shows.
pub fn expectWide(comptime T: type, want: Wide(T), got: T, ulps: f128, scale: f128) !void {
    if (comptime isComplex(T)) {
        try component(FloatOf(T), want.re, got.re, ulps, scale);
        try component(FloatOf(T), want.im, got.im, ulps, scale);
    } else try component(T, want, got, ulps, scale);
}

/// `expectWide` of one float with `F`'s range and precision; `got` may be
/// wider (`dsdot`).
pub fn component(comptime F: type, want: f128, got: anytype, ulps: f128, scale: f128) !void {
    if (componentOk(F, want, got, ulps, scale)) return;
    const g: f128 = got;
    check.say("got {e}, want {e}: ", .{ g, want });
    return error.TestExpectedApproxEqAbs;
}

/// `expectWide` where complex infinities may meet in another order than the
/// reference's: a complex part that should be infinite or NaN may come out
/// either (Inf times a finite complex value is NaN in one order of the
/// multiplications, infinite in another).
pub fn expectWideAnyOrder(comptime T: type, want: Wide(T), got: T, ulps: f128, scale: f128) !void {
    if (comptime isComplex(T)) {
        const lax = !std.math.isFinite(want.re) and !std.math.isFinite(got.re);
        const lay = !std.math.isFinite(want.im) and !std.math.isFinite(got.im);
        if (!lax) try component(FloatOf(T), want.re, got.re, ulps, scale);
        if (!lay) try component(FloatOf(T), want.im, got.im, ulps, scale);
        return;
    }
    return expectWide(T, want, got, ulps, scale);
}

/// Whether `expectWide` would pass, saying nothing.
pub fn wideOk(comptime T: type, want: Wide(T), got: T, ulps: f128, scale: f128) bool {
    if (comptime isComplex(T)) return componentOk(FloatOf(T), want.re, got.re, ulps, scale) and componentOk(FloatOf(T), want.im, got.im, ulps, scale);
    return componentOk(T, want, got, ulps, scale);
}

fn componentOk(comptime F: type, want: f128, got: anytype, ulps: f128, scale: f128) bool {
    const g: f128 = got;
    if (std.math.isNan(want)) return std.math.isNan(g);
    if (@abs(want) > std.math.floatMax(F)) return std.math.isInf(g) and (g > 0) == (want > 0);
    return std.math.isFinite(g) and @abs(g - want) <= ulps * (std.math.floatEps(F) * scale + std.math.floatTrueMin(F));
}
