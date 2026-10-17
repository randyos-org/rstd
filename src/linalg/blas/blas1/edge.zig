//! Edge cases through every level-1 kernel, on every core this machine runs,
//! against an `f128` reference: the special values of `check.Pattern`, NaN
//! in every gap between strided elements and past the last one (which must
//! stay as it was), and strides negative and large.
//!
//! NaN and infinities reach every output they feed, and denormals do not
//! flush. Not propagation, on purpose, as the reference has it:
//! - `axpby` with a zero `beta` never reads `y`, with a zero `alpha` never
//!   reads `x` (`kernel/arm/axpby.c`); both zero store zeros;
//! - `rotm` with flag -2 leaves `x` and `y` alone;
//! - the extremes replace only on a strict `>` or `<`: NaN is never picked
//!   but as the first element, which nothing then replaces;
//! - `nrm2` is Inf for an infinity and NaN for a NaN, as LAPACK 3.10's; past
//!   `f64` it is the plain sum of squares (`nrm2/`), so tiny and huge
//!   operands are not its cases there;
//! - a zero `alpha` in `scal` follows its flag (`scal/`'s `checkFlag`); in
//!   `axpy` it is the interface's quick return, never a kernel call.
//!
//! Where a core's kernels do otherwise upstream (`issues.md`), `Quirk` says
//! what they do, and that is checked instead.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../types.zig");
const Vtable = @import("../Vtable.zig");
const check = @import("check.zig");
const extremum = @import("extremum.zig");

const Allocator = std.mem.Allocator;
const BlasLong = types.BlasLong;
const XDouble = types.XDouble;
const C = std.math.Complex;
const Float = check.FloatOf;
const isComplex = check.isComplex;
const Wide = check.Wide;
const Vector = check.Vector;

/// Sizes past every unroll, for strides up to 3.
const sizes = [_]usize{ 1, 2, 5, 16, 33, 67 };
const max_n = 67;
const strides = [_][2]isize{ .{ 1, 1 }, .{ 2, 3 } };
/// For routines whose kernels take negative strides: the interface passes the
/// far end and the stride as given (`rotm` the near end; its kernel finds the
/// far one).
const negative = [_][2]isize{ .{ -1, -1 }, .{ -3, 2 }, .{ 2, -3 } };
/// Strides of a page and more, on fewer sizes and one pattern.
const large = [_][2]isize{ .{ 4097, 1 }, .{ -65537, 65537 } };
const large_sizes = [_]usize{ 2, 5 };

/// What some cores' kernels do otherwise upstream, checked instead of the
/// reference (`issues.md`).
const Quirk = enum {
    /// The extremes replace their pick unless the new key is `<=` it (`>=`
    /// for a minimum), so a NaN is picked where it is met and the next
    /// element replaces it.
    nan_replaced,
    /// The extremes compare lane by lane, so what a NaN does depends on
    /// where it falls: any element's key or NaN, an index from 1 to `n + 1`
    /// (the interface clamps it).
    nan_lanes,
    /// Of equal extremes, any may be picked.
    ties,
    /// `nrm2` of an infinity, and no NaN, may be NaN.
    nrm2_inf_nan,
    /// `nrm2` of a NaN may be 0.
    nrm2_nan_zero,
    /// `nrm2` of denormals may be anything but NaN: 0, infinite, or off.
    nrm2_denormal,
    /// `nrm2` of denormals may be NaN too (`frsqrte` of a denormal sum).
    nrm2_denormal_nan,
    /// `axpby` with a zero `alpha` (`beta`) computes the elements past its
    /// blocks in full, from `x` (`y`) elements its blocks left behind, so a
    /// NaN or infinity anywhere in that operand may reach them as NaN; the
    /// operand holds its values here, not NaN.
    axpby_reads,
    /// `copy` may lose the sign of a zero.
    copy_zero_sign,
};

const Quirks = std.EnumSet(Quirk);

/// The quirks of `routine` (`isamax`, `zcopy`) on the core `check.eachCore`
/// is running; the portable kernels have none.
fn quirks(comptime routine: []const u8) Quirks {
    const core = check.core orelse return .empty;
    // the routine without its index letter and precision
    const at = comptime if (routine[0] == 'i' and routine.len > 4) 1 else 0;
    const prec = comptime routine[at];
    const name = comptime routine[at + 1 ..];
    const extreme = comptime std.mem.endsWith(u8, name, "max") or std.mem.endsWith(u8, name, "min");
    const is = struct {
        fn eq(comptime a: []const u8, comptime b: []const u8) bool {
            return comptime std.mem.eql(u8, a, b);
        }
    }.eq;
    var q: Quirks = .empty;
    switch (sysinfo.cpu.arch) {
        .x86_64 => {
            if (extreme) q.insert(.nan_lanes);
            if (is(name, "copy") and (core == .opteron or core == .opteron_sse3)) q.insert(.copy_zero_sign);
        },
        .aarch64, .aarch64_be => if (is(name, "amax")) q.insert(.nan_lanes),
        .x86 => if (is(name, "amin") or is(name, "min")) q.insert(.nan_replaced),
        // the real ones only
        .arm, .armeb, .thumb, .thumbeb => if ((is(name, "amin") or is(name, "min")) and (prec == 's' or prec == 'd')) q.insert(.nan_replaced),
        .sparc, .sparc64 => {
            if (extreme) q.insert(.ties);
            if (is(name, "nrm2")) q.setUnion(.initMany(&.{ .nrm2_inf_nan, .nrm2_denormal }));
        },
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => {
            if (extreme) q.insert(.nan_lanes);
            if (is(name, "nrm2")) q.setUnion(.initMany(&.{ .nrm2_inf_nan, .nrm2_denormal }));
            if (is(name, "nrm2") and (core == .ppcg4 or core == .ppc440)) q.insert(.nrm2_denormal_nan);
        },
        .loongarch64 => {
            if (extreme) q.insert(.nan_lanes);
            if (is(name, "nrm2")) q.setUnion(.initMany(&.{ .nrm2_inf_nan, .nrm2_nan_zero, .nrm2_denormal }));
            if (is(name, "axpby")) q.insert(.axpby_reads);
        },
        .mips64, .mips64el => if (core == .loongson3r3 or core == .loongson3r4) {
            if (extreme) q.insert(.ties);
            if (is(name, "nrm2")) q.setUnion(.initMany(&.{ .nrm2_inf_nan, .nrm2_denormal }));
        },
        else => {},
    }
    return q;
}

/// A case: size, strides, pattern, the operand that has it, and where the
/// operands' storage comes from.
const Case = struct { n: usize, inc: [2]isize, p: check.Pattern, on_y: bool, al: Allocator };

/// Operand `which` of `c`: 0 for `x`, 1 for `y`.
fn vec(comptime T: type, c: Case, which: u1) !Vector(T) {
    return .init(c.al, c.n, c.inc[which], 31 * c.n + which, if ((which == 1) == c.on_y) c.p else null);
}

fn stride(v: anytype) BlasLong {
    return @intCast(v.inc);
}

/// Whether `v`, when `zeroed`, holds an infinity or NaN.
fn nonFinite(comptime T: type, v: anytype, zeroed: bool) bool {
    if (!zeroed) return false;
    for (v.was) |e| if (!std.math.isFinite(check.abs1(check.wide(T, e)))) return true;
    return false;
}

fn isNanW(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) std.math.isNan(v.re) or std.math.isNan(v.im) else std.math.isNan(v);
}

fn zero(comptime T: type) Wide(T) {
    return check.wide(T, if (comptime isComplex(T)) T.init(0, 0) else 0);
}

fn isZero(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) v.re == 0 and v.im == 0 else v == 0;
}

/// `a*u + b*v` for real `a` and `b`.
fn lin(comptime T: type, a: f128, u: Wide(T), b: f128, v: Wide(T)) Wide(T) {
    if (comptime isComplex(T)) return .init(a * u.re + b * v.re, a * u.im + b * v.im);
    return a * u + b * v;
}

/// Storage for the small cases, so the rigs do not spend their time in the
/// testing allocator: four operands of the most room any takes.
var small_mem: [4 * ((max_n - 1) * 3 + 1 + check.tail + max_n) * @sizeOf(C(f128))]u8 align(16) = undefined;

/// Runs `run(T, ctx, case)` over every size, stride and pattern, the pattern
/// in `x` and, for `two`, then in `y`; with `neg`, also over negative strides,
/// and over large ones always. Names the case that fails.
fn sweep(comptime T: type, routine: []const u8, two: bool, neg: bool, ctx: anytype, comptime run: anytype) !void {
    var fba: std.heap.FixedBufferAllocator = .init(&small_mem);
    for ([_]bool{ false, true }) |big| for (if (big) &large else &(strides ++ negative)) |inc| {
        if (!neg and (inc[0] < 0 or inc[1] < 0)) continue;
        for (if (big) &large_sizes else &sizes) |n| {
            var buf: [13]check.Pattern = undefined;
            const ps: []const check.Pattern = if (big) &.{.mixed} else check.patterns(n, &buf);
            for (ps) |p| for ([_]bool{ false, true }) |on_y| {
                if (on_y and !two) continue;
                fba.reset();
                const c: Case = .{ .n = n, .inc = inc, .p = p, .on_y = on_y, .al = if (big) std.testing.allocator else fba.allocator() };
                run(T, ctx, c) catch |e| {
                    p.describe();
                    check.say("in {s}, {s} {s} n={d} inc={d},{d}: ", .{ if (on_y) "y" else "x", routine, @typeName(T), n, inc[0], inc[1] });
                    return e;
                };
            };
        }
    };
}

fn Copy(comptime K: type) type {
    return struct { kern: K, q: Quirks };
}

fn copyRun(comptime T: type, ctx: anytype, c: Case) !void {
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    var y: Vector(T) = try vec(T, c, 1);
    defer y.deinit(c.al);
    _ = ctx.kern(@intCast(c.n), x.first(), stride(x), y.first(), stride(y));
    try x.unchanged();
    for (0..c.n) |i| if (!check.same(T, x.was[i], y.at(i), !ctx.q.contains(.copy_zero_sign))) {
        check.say("element {d} not copied exactly: ", .{i});
        return error.TestUnexpectedResult;
    };
    try y.gaps();
}

fn swapRun(comptime T: type, kern: anytype, c: Case) !void {
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
    for (0..c.n) |i| if (!check.same(T, x.was[i], y.at(i), true) or !check.same(T, y.was[i], x.at(i), true)) {
        check.say("element {d} not swapped exactly: ", .{i});
        return error.TestUnexpectedResult;
    };
    try x.gaps();
    try y.gaps();
}

fn Scal(comptime T: type, comptime K: type) type {
    return struct { kern: K, alpha: T };
}

/// `x := alpha*x` with the interface's flag; a zero keeps its sign.
fn scalRun(comptime T: type, ctx: anytype, c: Case) !void {
    errdefer check.say("alpha {any}: ", .{ctx.alpha});
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    const a = ctx.alpha;
    const n: BlasLong = @intCast(c.n);
    if (comptime isComplex(T)) {
        _ = ctx.kern(n, 0, 0, a.re, a.im, x.first(), stride(x), x.first(), 0, x.first(), 1);
    } else {
        _ = ctx.kern(n, 0, 0, a, x.first(), stride(x), x.first(), 0, x.first(), 1);
    }
    const wa = check.wide(T, a);
    for (0..c.n) |i| {
        const want = check.mulW(wa, x.w(i));
        try check.expectWide(T, want, x.at(i), 4, check.mag(wa) * check.mag(x.w(i)));
        if (comptime !isComplex(T)) if (want == 0 and std.math.signbit(want) != std.math.signbit(x.at(i))) {
            check.say("element {d}: zero of the wrong sign: ", .{i});
            return error.TestUnexpectedResult;
        };
    }
    try x.gaps();
}

fn Axpy(comptime T: type, comptime K: type) type {
    return struct { kern: K, alpha: T, conj: bool };
}

/// `y += alpha*x`, `x` conjugated for `axpyc`.
fn axpyRun(comptime T: type, ctx: anytype, c: Case) !void {
    errdefer check.say("alpha {any}: ", .{ctx.alpha});
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    var y: Vector(T) = try vec(T, c, 1);
    defer y.deinit(c.al);
    const a = ctx.alpha;
    const n: BlasLong = @intCast(c.n);
    if (comptime isComplex(T)) {
        _ = ctx.kern(n, 0, 0, a.re, a.im, x.first(), stride(x), y.first(), stride(y), y.first(), 0);
    } else {
        _ = ctx.kern(n, 0, 0, a, x.first(), stride(x), y.first(), stride(y), y.first(), 0);
    }
    const wa = check.wide(T, a);
    for (0..c.n) |i| {
        const xv = if (ctx.conj) check.conjW(x.w(i)) else x.w(i);
        const want = check.addW(y.w(i), check.mulW(wa, xv));
        try check.expectWide(T, want, y.at(i), 4, check.mag(y.w(i)) + check.mag(wa) * check.mag(xv));
    }
    try x.unchanged();
    try y.gaps();
}

fn Axpby(comptime T: type, comptime K: type) type {
    return struct { kern: K, alpha: T, beta: T, q: Quirks };
}

/// `y := alpha*x + beta*y`; a zero scalar's operand holds NaN, unread.
fn axpbyRun(comptime T: type, ctx: anytype, c: Case) !void {
    errdefer check.say("alpha {any}, beta {any}: ", .{ ctx.alpha, ctx.beta });
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    var y: Vector(T) = try vec(T, c, 1);
    defer y.deinit(c.al);
    const az = isZero(T, ctx.alpha);
    const bz = isZero(T, ctx.beta);
    // with `axpby_reads` the unread operand holds its values, and an element
    // may be either sum
    const reads = ctx.q.contains(.axpby_reads);
    if (!reads) {
        if (az) x.poison();
        if (bz) y.poison();
    }
    const a = ctx.alpha;
    const b = ctx.beta;
    const n: BlasLong = @intCast(c.n);
    if (comptime isComplex(T)) {
        _ = ctx.kern(n, a.re, a.im, x.first(), stride(x), b.re, b.im, y.first(), stride(y));
    } else {
        _ = ctx.kern(n, a, x.first(), stride(x), b, y.first(), stride(y));
    }
    const wa = check.wide(T, a);
    const wb = check.wide(T, b);
    for (0..c.n) |i| {
        const ax = if (az) zero(T) else check.mulW(wa, x.w(i));
        const by = if (bz) zero(T) else check.mulW(wb, y.w(i));
        const scale = check.mag(wa) * check.mag(x.w(i)) + check.mag(wb) * check.mag(y.w(i));
        const full = check.addW(check.mulW(wa, x.w(i)), check.mulW(wb, y.w(i)));
        if (reads and (check.wideOk(T, full, y.at(i), 4, scale) or (nonFinite(T, &x, az) or nonFinite(T, &y, bz)) and isNanW(T, y.at(i)))) continue;
        check.expectWide(T, check.addW(ax, by), y.at(i), 4, scale) catch |e| {
            check.say("y[{d}]: ", .{i});
            return e;
        };
    }
    try x.unchanged();
    try y.gaps();
}

fn Rot(comptime K: type) type {
    return struct { kern: K, c: f64, s: f64 };
}

/// `x := c*x + s*y`, `y := c*y - s*x`.
fn rotRun(comptime T: type, ctx: anytype, c: Case) !void {
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    var y: Vector(T) = try vec(T, c, 1);
    defer y.deinit(c.al);
    const F = Float(T);
    _ = ctx.kern(@intCast(c.n), x.first(), stride(x), y.first(), stride(y), @as(F, @floatCast(ctx.c)), @as(F, @floatCast(ctx.s)));
    for (0..c.n) |i| {
        const scale = @abs(ctx.c) * check.mag(x.w(i)) + @abs(ctx.s) * check.mag(y.w(i));
        try check.expectWide(T, lin(T, ctx.c, x.w(i), ctx.s, y.w(i)), x.at(i), 4, scale);
        try check.expectWide(T, lin(T, ctx.c, y.w(i), -ctx.s, x.w(i)), y.at(i), 4, scale);
    }
    try x.gaps();
    try y.gaps();
}

fn Rotm(comptime K: type) type {
    return struct { kern: K, flag: f64 };
}

/// `[x y] := H [x y]` for the flag's `H`; flag -2 changes nothing.
fn rotmRun(comptime T: type, ctx: anytype, c: Case) !void {
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
        const scale = check.mag(x.w(i)) + check.mag(y.w(i));
        try check.expectWide(T, lin(T, h[0], x.w(i), h[1], y.w(i)), x.at(i), 4, scale);
        try check.expectWide(T, lin(T, h[2], x.w(i), h[3], y.w(i)), y.at(i), 4, scale);
    }
    try x.gaps();
    try y.gaps();
}

fn Dot(comptime K: type) type {
    return struct { kern: K, conj: bool };
}

/// `x . y`, `x` conjugated for `dotc`; a `dsdot` result is `f64`.
fn dotRun(comptime T: type, ctx: anytype, c: Case) !void {
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    var y: Vector(T) = try vec(T, c, 1);
    defer y.deinit(c.al);
    const got = ctx.kern(@intCast(c.n), x.first(), stride(x), y.first(), stride(y));
    var want = zero(T);
    var scale: f128 = 0;
    for (0..c.n) |i| {
        const xv = if (ctx.conj) check.conjW(x.w(i)) else x.w(i);
        want = check.addW(want, check.mulW(xv, y.w(i)));
        scale += check.mag(xv) * check.mag(y.w(i));
    }
    const ulps: f128 = @floatFromInt(4 * (c.n + 2));
    if (comptime isComplex(T)) {
        try check.expectWide(T, want, .init(got.re, got.im), ulps, scale);
    } else {
        try check.component(T, want, got, ulps, scale);
    }
    try x.unchanged();
    try y.unchanged();
}

const Reduce = enum { asum, sum, nrm2 };

fn Reduction(comptime K: type) type {
    return struct { kern: K, op: Reduce, q: Quirks };
}

/// `asum` (of `|re| + |im|`), `sum` (of `re + im`) or `nrm2`.
fn reduceRun(comptime T: type, ctx: anytype, c: Case) !void {
    if (comptime types.isExtended(Float(T))) if (ctx.op == .nrm2 and (c.p == .tiny or c.p == .huge)) return;
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    const got = ctx.kern(@intCast(c.n), x.first(), stride(x));
    var acc: f128 = 0;
    var mag: f128 = 0;
    for (0..c.n) |i| {
        const v = x.w(i);
        mag += check.mag(v);
        acc += switch (ctx.op) {
            .asum => check.abs1(v),
            .sum => if (comptime isComplex(T)) v.re + v.im else v,
            .nrm2 => if (comptime isComplex(T)) v.re * v.re + v.im * v.im else v * v,
        };
    }
    const want = if (ctx.op == .nrm2) check.sqrtW(acc) else acc;
    if (std.math.isInf(want) and std.math.isNan(got) and ctx.q.contains(.nrm2_inf_nan)) return;
    if (std.math.isNan(want) and got == 0 and ctx.q.contains(.nrm2_nan_zero)) return;
    if (want > 0 and want < std.math.floatMin(Float(T)) and ctx.q.contains(.nrm2_denormal)) {
        if (!std.math.isNan(got) or ctx.q.contains(.nrm2_denormal_nan)) return;
    }
    const ulps: f128 = @floatFromInt(2 * (c.n + 2));
    try check.component(Float(T), want, got, ulps, if (ctx.op == .nrm2) want else mag);
    try x.unchanged();
}

fn Extreme(comptime K: type) type {
    return struct { kern: K, find: extremum.Find, q: Quirks };
}

/// The 1-based index of the extreme element by the reference's scan, which
/// replaces only on a strict compare, or with `replaced` unless the compare
/// the other way holds (`nan_replaced`).
fn scan(comptime T: type, x: anytype, find: extremum.Find, replaced: bool) usize {
    var best: usize = 0;
    var key: f128 = 0;
    for (0..x.n) |i| {
        const k = keyOf(T, x, find, i);
        const better = if (replaced) !(if (find.max) k <= key else k >= key) else if (find.max) k > key else k < key;
        if (best == 0 or better) {
            best = i + 1;
            key = k;
        }
    }
    return best;
}

fn keyOf(comptime T: type, x: anytype, find: extremum.Find, i: usize) f128 {
    const v = x.w(i);
    return if (comptime isComplex(T)) check.abs1(v) else if (find.abs) @abs(v) else v;
}

/// The extreme's value, or its 1-based index, by the kernel's return type.
fn extremeRun(comptime T: type, ctx: anytype, c: Case) !void {
    var x: Vector(T) = try vec(T, c, 0);
    defer x.deinit(c.al);
    const got = ctx.kern(@intCast(c.n), x.first(), stride(x));
    try x.unchanged();
    var nan = false;
    for (0..c.n) |i| nan = nan or std.math.isNan(keyOf(T, &x, ctx.find, i));
    if (nan and ctx.q.contains(.nan_lanes)) return loose(T, &x, ctx.find, got);
    const ref = scan(T, &x, ctx.find, ctx.q.contains(.nan_replaced));
    const key = keyOf(T, &x, ctx.find, ref - 1);
    if (@TypeOf(got) == BlasLong) {
        const n: BlasLong = @intCast(c.n);
        if (got == @as(BlasLong, @intCast(ref))) return;
        // of equal keys, any
        if (ctx.q.contains(.ties) and got >= 1 and got <= n and keyOf(T, &x, ctx.find, @intCast(got - 1)) == key) return;
    } else {
        if (std.math.isNan(key) and std.math.isNan(got) or got == key) return;
        if (@abs(@as(f128, got) - key) <= 2 * std.math.floatEps(Float(T)) * @abs(key)) return;
    }
    check.say("{any}, want {d}: ", .{ got, if (@TypeOf(got) == BlasLong) @as(f128, @floatFromInt(ref)) else keyOf(T, &x, ctx.find, ref - 1) });
    return error.TestUnexpectedResult;
}

/// A `nan_lanes` result: NaN or any element's key, or an index from 1 to
/// `n + 1`.
fn loose(comptime T: type, x: anytype, find: extremum.Find, got: anytype) !void {
    if (@TypeOf(got) == BlasLong) {
        if (got >= 1 and got <= @as(BlasLong, @intCast(x.n)) + 1) return;
    } else {
        // or 0, where a maximum starts from
        if (std.math.isNan(got) or got == 0) return;
        for (0..x.n) |i| {
            const k = keyOf(T, x, find, i);
            if (got == k) return;
            if (@abs(@as(f128, got) - k) <= 2 * std.math.floatEps(Float(T)) * @abs(k)) return;
        }
    }
    check.say("{any} is no element's: ", .{got});
    return error.TestUnexpectedResult;
}

/// The first failure of each routine, so one run names every one.
var failed: ?anyerror = null;

fn note(r: anyerror!void) void {
    r catch |e| {
        check.say("\n", .{});
        failed = e;
    };
}

/// Every routine of `T`, whose slots begin with `l`.
fn slots(comptime l: []const u8, comptime T: type, vt: *const Vtable) void {
    const complex = comptime isComplex(T);
    const F = Float(T);
    const nan = std.math.nan(F);
    const one: T = if (complex) .init(1.5, -0.75) else 1.5;
    const inf: T = if (complex) .init(std.math.inf(F), 0.5) else std.math.inf(F);
    // NaN, and for complex a NaN beside a zero, which a test for zero by
    // `comiss` alone takes for zero
    const nan_a: T = if (complex) .init(nan, 0.5) else nan;
    const nan_b: T = if (complex) .init(0, nan) else nan;
    const zero_t: T = if (complex) .init(0, 0) else 0;
    const beta: T = if (complex) .init(-0.5, 0.25) else -0.5;
    const copy = @field(vt, l ++ "copy_k").?;
    note(sweep(T, l ++ "copy", false, true, Copy(@TypeOf(copy)){ .kern = copy, .q = quirks(l ++ "copy") }, copyRun));
    note(sweep(T, l ++ "swap", true, true, @field(vt, l ++ "swap_k").?, swapRun));
    const scal = @field(vt, l ++ "scal_k").?;
    for ([_]T{ one, inf, nan_a, nan_b }) |a| note(sweep(T, l ++ "scal", false, false, Scal(T, @TypeOf(scal)){ .kern = scal, .alpha = a }, scalRun));
    const axpy = @field(vt, l ++ "axpy_k").?;
    for ([_]T{ one, nan_a, nan_b }) |a| {
        note(sweep(T, l ++ "axpy", true, true, Axpy(T, @TypeOf(axpy)){ .kern = axpy, .alpha = a, .conj = false }, axpyRun));
        if (complex) note(sweep(T, l ++ "axpyc", true, true, Axpy(T, @TypeOf(axpy)){ .kern = @field(vt, l ++ "axpyc_k").?, .alpha = a, .conj = true }, axpyRun));
    }
    if (@hasField(Vtable, l ++ "axpby_k")) if (@field(vt, l ++ "axpby_k")) |k| {
        for ([_][2]T{ .{ one, beta }, .{ zero_t, beta }, .{ one, zero_t }, .{ zero_t, zero_t } }) |ab| {
            note(sweep(T, l ++ "axpby", true, true, Axpby(T, @TypeOf(k)){ .kern = k, .alpha = ab[0], .beta = ab[1], .q = quirks(l ++ "axpby") }, axpbyRun));
        }
    };
    const rot = @field(vt, if (complex) l ++ (if (l[0] == 'c') "s" else if (l[0] == 'z') "d" else "q") ++ "rot_k" else l ++ "rot_k").?;
    note(sweep(T, l ++ "rot", true, true, Rot(@TypeOf(rot)){ .kern = rot, .c = 0.6, .s = 0.8 }, rotRun));
    if (!complex) for ([_]f64{ -2, -1, 0, 1 }) |flag| {
        const k = @field(vt, l ++ "rotm_k").?;
        note(sweep(T, l ++ "rotm", true, true, Rotm(@TypeOf(k)){ .kern = k, .flag = flag }, rotmRun));
    };
    inline for (if (complex) .{ .{ "dotu", false }, .{ "dotc", true } } else .{.{ "dot", false }}) |d| {
        const k = @field(vt, l ++ d[0] ++ "_k").?;
        note(sweep(T, l ++ d[0], true, true, Dot(@TypeOf(k)){ .kern = k, .conj = d[1] }, dotRun));
    }
    inline for (comptime std.enums.values(Reduce)) |op| {
        const k = @field(vt, l ++ @tagName(op) ++ "_k").?;
        note(sweep(T, l ++ @tagName(op), false, op == .nrm2, Reduction(@TypeOf(k)){ .kern = k, .op = op, .q = quirks(l ++ @tagName(op)) }, reduceRun));
    }
    inline for (.{ "amax", "amin", "max", "min" }) |name| {
        if (comptime complex and name[0] != 'a') continue;
        const find: extremum.Find = .{ .abs = name[0] == 'a', .max = comptime std.mem.endsWith(u8, name, "max") };
        inline for (.{ "", "i" }) |i| {
            const k = @field(vt, i ++ l ++ name ++ "_k").?;
            note(sweep(T, i ++ l ++ name, false, false, Extreme(@TypeOf(k)){ .kern = k, .find = find, .q = quirks(i ++ l ++ name) }, extremeRun));
        }
    }
}

const vtable_types = .{ .{ "s", f32 }, .{ "d", f64 }, .{ "c", C(f32) }, .{ "z", C(f64) } };

fn checkAll(vt: *const Vtable) !void {
    failed = null;
    inline for (vtable_types) |p| slots(p[0], p[1], vt);
    note(sweep(f32, "dsdot", true, true, Dot(@TypeOf(vt.dsdot_k.?)){ .kern = vt.dsdot_k.?, .conj = false }, dotRun));
    if (failed) |e| return e;
}

test "portable level-1 kernels: special values, poisoned gaps, and negative and large strides" {
    var vt: Vtable = .{};
    @import("__root__.zig").installGeneric(&vt);
    try checkAll(&vt);
    failed = null;
    if (XDouble != f64) {
        slots("q", XDouble, &vt);
        slots("x", C(XDouble), &vt);
    }
    if (failed) |e| return e;
}

test "every linked core's level-1 kernels: special values, poisoned gaps, and negative and large strides" {
    try check.eachCore(checkAll);
}
