//! One case of each level-1 routine against the `f128` reference.

const std = @import("std");
const types = @import("../../../types.zig");
const util = @import("../../../util.zig");
const check = @import("../../check.zig");
const bitwise = @import("../../bitwise.zig");
const pattern = @import("../../pattern.zig");
const reference = @import("../../reference.zig");
const extremum = @import("../../../blas1/extremum.zig");
const ed = @import("../edge.zig");

const BlasLong = types.BlasLong;
const Float = check.FloatOf;
const isComplex = check.isComplex;
const Vector = check.Vector;
const Case = ed.Case;
const isZero = ed.isZero;
const lin = ed.lin;
const stride = ed.stride;
const vec = ed.vec;
const zero = ed.zero;

pub fn copyRun(comptime T: type, kern: anytype, c: Case) !void {
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    var y: Vector(T) = try vec(T, c, 1);
    defer y.deinit(c.al);
    _ = kern(@intCast(c.n), x.first(), stride(x), y.first(), stride(y));
    try x.unchanged();
    for (0..c.n) |i| if (!bitwise.same(T, x.was[i], y.at(i), true)) {
        check.say("element {d} not copied exactly: ", .{i});
        return error.TestUnexpectedResult;
    };
    try y.gaps();
}

pub fn swapRun(comptime T: type, kern: anytype, c: Case) !void {
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    var y: Vector(T) = try vec(T, c, 1);
    defer y.deinit(c.al);
    const n: BlasLong = @intCast(c.n);
    if (comptime isComplex(T)) {
        _ = kern(n, 0, 0, 0, 0, x.first(), stride(x), y.first(), stride(y), y.first(), 0);
    } else {
        _ = kern(n, 0, 0, 0, x.first(), stride(x), y.first(), stride(y), y.first(), 0);
    }
    for (0..c.n) |i| if (!bitwise.same(T, x.was[i], y.at(i), true) or !bitwise.same(T, y.was[i], x.at(i), true)) {
        check.say("element {d} not swapped exactly: ", .{i});
        return error.TestUnexpectedResult;
    };
    try x.gaps();
    try y.gaps();
}

pub fn Scal(comptime T: type, comptime K: type) type {
    return struct { kern: K, alpha: T };
}

/// `x := alpha*x` with the interface's flag; a zero keeps its sign.
pub fn scalRun(comptime T: type, ctx: anytype, c: Case) !void {
    errdefer check.say("alpha {any}: ", .{ctx.alpha});
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    const a = ctx.alpha;
    const n: BlasLong = @intCast(c.n);
    _ = @call(.auto, ctx.kern, .{ n, 0, 0 } ++ util.scalarArgs(T, a) ++ .{ x.first(), stride(x), x.first(), 0, x.first(), 1 });
    const wa = reference.wide(T, a);
    for (0..c.n) |i| {
        const want = reference.mulW(wa, x.w(i));
        try pattern.expectWide(T, want, x.at(i), reference.elem_ulps, reference.mag(wa) * reference.mag(x.w(i)));
        if (comptime !isComplex(T)) if (want == 0 and std.math.signbit(want) != std.math.signbit(x.at(i))) {
            check.say("element {d}: zero of the wrong sign: ", .{i});
            return error.TestUnexpectedResult;
        };
    }
    try x.gaps();
}

pub fn Axpy(comptime T: type, comptime K: type) type {
    return struct { kern: K, alpha: T, conj: bool };
}

/// `y += alpha*x`, `x` conjugated for `axpyc`.
pub fn axpyRun(comptime T: type, ctx: anytype, c: Case) !void {
    errdefer check.say("alpha {any}: ", .{ctx.alpha});
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    var y: Vector(T) = try vec(T, c, 1);
    defer y.deinit(c.al);
    const a = ctx.alpha;
    const n: BlasLong = @intCast(c.n);
    _ = @call(.auto, ctx.kern, .{ n, 0, 0 } ++ util.scalarArgs(T, a) ++ .{ x.first(), stride(x), y.first(), stride(y), y.first(), 0 });
    const wa = reference.wide(T, a);
    for (0..c.n) |i| {
        const xv = if (ctx.conj) reference.conjW(x.w(i)) else x.w(i);
        const want = reference.addW(y.w(i), reference.mulW(wa, xv));
        try pattern.expectWide(T, want, y.at(i), reference.elem_ulps, reference.mag(y.w(i)) + reference.mag(wa) * reference.mag(xv));
    }
    try x.unchanged();
    try y.gaps();
}

pub fn Axpby(comptime T: type, comptime K: type) type {
    return struct { kern: K, alpha: T, beta: T };
}

/// `y := alpha*x + beta*y`; a zero scalar's operand holds NaN, unread.
pub fn axpbyRun(comptime T: type, ctx: anytype, c: Case) !void {
    errdefer check.say("alpha {any}, beta {any}: ", .{ ctx.alpha, ctx.beta });
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    var y: Vector(T) = try vec(T, c, 1);
    defer y.deinit(c.al);
    const az = isZero(T, ctx.alpha);
    const bz = isZero(T, ctx.beta);
    if (az) x.poison();
    if (bz) y.poison();
    const a = ctx.alpha;
    const b = ctx.beta;
    const n: BlasLong = @intCast(c.n);
    _ = @call(.auto, ctx.kern, .{n} ++ util.scalarArgs(T, a) ++ .{ x.first(), stride(x) } ++ util.scalarArgs(T, b) ++ .{ y.first(), stride(y) });
    const wa = reference.wide(T, a);
    const wb = reference.wide(T, b);
    for (0..c.n) |i| {
        const ax = if (az) zero(T) else reference.mulW(wa, x.w(i));
        const by = if (bz) zero(T) else reference.mulW(wb, y.w(i));
        const scale = reference.mag(wa) * reference.mag(x.w(i)) + reference.mag(wb) * reference.mag(y.w(i));
        pattern.expectWide(T, reference.addW(ax, by), y.at(i), reference.elem_ulps, scale) catch |e| {
            check.say("y[{d}]: ", .{i});
            return e;
        };
    }
    try x.unchanged();
    try y.gaps();
}

pub fn Rot(comptime K: type) type {
    return struct { kern: K, c: f64, s: f64 };
}

/// `x := c*x + s*y`, `y := c*y - s*x`.
pub fn rotRun(comptime T: type, ctx: anytype, c: Case) !void {
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    var y: Vector(T) = try vec(T, c, 1);
    defer y.deinit(c.al);
    const F = Float(T);
    _ = ctx.kern(@intCast(c.n), x.first(), stride(x), y.first(), stride(y), @as(F, @floatCast(ctx.c)), @as(F, @floatCast(ctx.s)));
    for (0..c.n) |i| {
        const scale = @abs(ctx.c) * reference.mag(x.w(i)) + @abs(ctx.s) * reference.mag(y.w(i));
        try pattern.expectWide(T, lin(T, ctx.c, x.w(i), ctx.s, y.w(i)), x.at(i), reference.elem_ulps, scale);
        try pattern.expectWide(T, lin(T, ctx.c, y.w(i), -ctx.s, x.w(i)), y.at(i), reference.elem_ulps, scale);
    }
    try x.gaps();
    try y.gaps();
}

pub fn Rotm(comptime K: type) type {
    return struct { kern: K, flag: f64 };
}

/// `[x y] := H [x y]` for the flag's `H`; flag -2 changes nothing.
pub fn rotmRun(comptime T: type, ctx: anytype, c: Case) !void {
    errdefer check.say("flag {d}: ", .{ctx.flag});
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    var y: Vector(T) = try vec(T, c, 1);
    defer y.deinit(c.al);
    var param = [5]T{ @floatCast(ctx.flag), 0.75, -0.5, 0.25, 1.25 };
    _ = ctx.kern(@intCast(c.n), @ptrCast(x.base()), stride(x), @ptrCast(y.base()), stride(y), &param);
    if (ctx.flag == -2) {
        try x.unchanged();
        return y.unchanged();
    }
    // h11, h12, h21, h22
    const h: [4]f128 = if (ctx.flag < 0) .{ 0.75, 0.25, -0.5, 1.25 } else if (ctx.flag == 0) .{ 1, 0.25, -0.5, 1 } else .{ 0.75, 1, -1, 1.25 };
    for (0..c.n) |i| {
        const scale = reference.mag(x.w(i)) + reference.mag(y.w(i));
        try pattern.expectWide(T, lin(T, h[0], x.w(i), h[1], y.w(i)), x.at(i), reference.elem_ulps, scale);
        try pattern.expectWide(T, lin(T, h[2], x.w(i), h[3], y.w(i)), y.at(i), reference.elem_ulps, scale);
    }
    try x.gaps();
    try y.gaps();
}

pub fn Dot(comptime K: type) type {
    return struct { kern: K, conj: bool };
}

/// `x . y`, `x` conjugated for `dotc`; a `dsdot` result is `f64`. `y` never
/// takes `both`'s pattern: denormal times denormal is no result.
pub fn dotRun(comptime T: type, ctx: anytype, c: Case) !void {
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    var y: Vector(T) = try vec(T, c.alone(), 1);
    defer y.deinit(c.al);
    const got = ctx.kern(@intCast(c.n), x.first(), stride(x), y.first(), stride(y));
    var want = zero(T);
    var scale: f128 = 0;
    for (0..c.n) |i| {
        const xv = if (ctx.conj) reference.conjW(x.w(i)) else x.w(i);
        want = reference.addW(want, reference.mulW(xv, y.w(i)));
        scale += reference.mag(xv) * reference.mag(y.w(i));
    }
    const ulps = reference.dotUlps(c.n);
    if (comptime isComplex(T)) {
        try pattern.expectWide(T, want, .init(got.re, got.im), ulps, scale);
    } else {
        try pattern.component(T, want, got, ulps, scale);
    }
    try x.unchanged();
    try y.unchanged();
}

pub const Reduce = check.Reduce;

pub fn Reduction(comptime K: type) type {
    return struct { kern: K, op: Reduce };
}

/// `asum` (of `|re| + |im|`), `sum` (of `re + im`) or `nrm2`.
pub fn reduceRun(comptime T: type, ctx: anytype, c: Case) !void {
    if (comptime types.isExtended(Float(T))) if (ctx.op == .nrm2 and (c.p == .tiny or c.p == .huge)) return;
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    const got = ctx.kern(@intCast(c.n), x.first(), stride(x));
    var acc: f128 = 0;
    var mag: f128 = 0;
    for (0..c.n) |i| {
        const v = x.w(i);
        mag += reference.mag(v);
        acc += switch (ctx.op) {
            .asum => reference.abs1(v),
            .sum => if (comptime isComplex(T)) v.re + v.im else v,
            .nrm2 => if (comptime isComplex(T)) v.re * v.re + v.im * v.im else v * v,
        };
    }
    const want = if (ctx.op == .nrm2) reference.sqrtW(acc) else acc;
    const ulps = reference.sumUlps(c.n);
    try pattern.component(Float(T), want, got, ulps, if (ctx.op == .nrm2) want else mag);
    try x.unchanged();
}

pub fn Extreme(comptime K: type) type {
    return struct { kern: K, find: extremum.Find };
}

/// The 1-based index of the first NaN, else of the extreme element by a scan
/// that replaces only on a strict compare.
pub fn scan(comptime T: type, x: anytype, find: extremum.Find) usize {
    var best: usize = 0;
    var key: f128 = 0;
    for (0..x.n) |i| {
        const k = keyOf(T, x, find, i);
        if (std.math.isNan(k)) return i + 1;
        const better = if (find.max) k > key else k < key;
        if (best == 0 or better) {
            best = i + 1;
            key = k;
        }
    }
    return best;
}

pub fn keyOf(comptime T: type, x: anytype, find: extremum.Find, i: usize) f128 {
    const v = x.w(i);
    return if (comptime isComplex(T)) reference.abs1(v) else if (find.abs) @abs(v) else v;
}

/// The extreme's value, or its 1-based index, by the kernel's return type.
pub fn extremeRun(comptime T: type, ctx: anytype, c: Case) !void {
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    const got = ctx.kern(@intCast(c.n), x.first(), stride(x));
    try x.unchanged();
    const ref = scan(T, &x, ctx.find);
    const key = keyOf(T, &x, ctx.find, ref - 1);
    if (@TypeOf(got) == BlasLong) {
        if (got == @as(BlasLong, @intCast(ref))) return;
    } else {
        if ((std.math.isNan(key) and std.math.isNan(got)) or got == key) return;
        if (@abs(@as(f128, got) - key) <= reference.keyTol(Float(T)) * @abs(key)) return;
    }
    check.say("{any}, want {d}: ", .{ got, if (@TypeOf(got) == BlasLong) @as(f128, @floatFromInt(ref)) else keyOf(T, &x, ctx.find, ref - 1) });
    return error.TestUnexpectedResult;
}
