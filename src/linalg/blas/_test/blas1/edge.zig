//! Edge cases through every level-1 kernel on every core this machine runs,
//! against an `f128` reference: `pattern.Pattern`'s special values, NaN in the
//! gaps, negative and large strides; what may not propagate, `design.md`'s Verification.

const std = @import("std");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");
const util = @import("../../util.zig");
const check = @import("../check.zig");
const bitwise = @import("../bitwise.zig");
const pattern = @import("../pattern.zig");
const reference = @import("../reference.zig");
const extremum = @import("../../blas1/extremum.zig");

const Allocator = std.mem.Allocator;
const BlasLong = types.BlasLong;
const C = std.math.Complex;
const Float = check.FloatOf;
const isComplex = check.isComplex;
const Wide = reference.Wide;
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

/// A case: size, strides, pattern, the operand that has it, and where the
/// operands' storage comes from.
const Case = struct { n: usize, inc: [2]isize, p: pattern.Pattern, on_y: bool, al: Allocator };

/// Operand `which` of `c`: 0 for `x`, 1 for `y`.
fn vec(comptime T: type, c: Case, which: u1) !Vector(T) {
    return .init(c.al, c.n, c.inc[which], 31 * c.n + which, if ((which == 1) == c.on_y) c.p else null);
}

fn stride(v: anytype) BlasLong {
    return @intCast(v.inc);
}

fn zero(comptime T: type) Wide(T) {
    return reference.wide(T, if (comptime isComplex(T)) T.init(0, 0) else 0);
}

const isZero = util.isZero;

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
            var buf: [13]pattern.Pattern = undefined;
            const ps: []const pattern.Pattern = if (big) &.{.mixed} else pattern.patterns(n, &buf);
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

fn copyRun(comptime T: type, kern: anytype, c: Case) !void {
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
    for (0..c.n) |i| if (!bitwise.same(T, x.was[i], y.at(i), true) or !bitwise.same(T, y.was[i], x.at(i), true)) {
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

fn Axpby(comptime T: type, comptime K: type) type {
    return struct { kern: K, alpha: T, beta: T };
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
        const scale = @abs(ctx.c) * reference.mag(x.w(i)) + @abs(ctx.s) * reference.mag(y.w(i));
        try pattern.expectWide(T, lin(T, ctx.c, x.w(i), ctx.s, y.w(i)), x.at(i), reference.elem_ulps, scale);
        try pattern.expectWide(T, lin(T, ctx.c, y.w(i), -ctx.s, x.w(i)), y.at(i), reference.elem_ulps, scale);
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
        const scale = reference.mag(x.w(i)) + reference.mag(y.w(i));
        try pattern.expectWide(T, lin(T, h[0], x.w(i), h[1], y.w(i)), x.at(i), reference.elem_ulps, scale);
        try pattern.expectWide(T, lin(T, h[2], x.w(i), h[3], y.w(i)), y.at(i), reference.elem_ulps, scale);
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

const Reduce = check.Reduce;

fn Reduction(comptime K: type) type {
    return struct { kern: K, op: Reduce };
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

fn Extreme(comptime K: type) type {
    return struct { kern: K, find: extremum.Find };
}

/// The 1-based index of the first NaN, else of the extreme element by a scan
/// that replaces only on a strict compare.
fn scan(comptime T: type, x: anytype, find: extremum.Find) usize {
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

fn keyOf(comptime T: type, x: anytype, find: extremum.Find, i: usize) f128 {
    const v = x.w(i);
    return if (comptime isComplex(T)) reference.abs1(v) else if (find.abs) @abs(v) else v;
}

/// The extreme's value, or its 1-based index, by the kernel's return type.
fn extremeRun(comptime T: type, ctx: anytype, c: Case) !void {
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
    note(sweep(T, l ++ "copy", false, true, portable.slotPtr(vt, .from(l), .copy_k).*.?, copyRun));
    note(sweep(T, l ++ "swap", true, true, portable.slotPtr(vt, .from(l), .swap_k).*.?, swapRun));
    const scal = portable.slotPtr(vt, .from(l), .scal_k).*.?;
    for ([_]T{ one, inf, nan_a, nan_b }) |a| note(sweep(T, l ++ "scal", false, false, Scal(T, @TypeOf(scal)){ .kern = scal, .alpha = a }, scalRun));
    const axpy = portable.slotPtr(vt, .from(l), .axpy_k).*.?;
    // a zero alpha too, which the triangular drivers pass (`-x(j)`): `0 * x`
    for ([_]T{ one, zero_t, nan_a, nan_b }) |a| {
        note(sweep(T, l ++ "axpy", true, true, Axpy(T, @TypeOf(axpy)){ .kern = axpy, .alpha = a, .conj = false }, axpyRun));
        if (complex) note(sweep(T, l ++ "axpyc", true, true, Axpy(T, @TypeOf(axpy)){ .kern = portable.slotPtr(vt, .from(l), .axpyc_k).*.?, .alpha = a, .conj = true }, axpyRun));
    }
    if (@hasField(Vtable, portable.slotName(.from(l), .axpby_k))) if (portable.slotPtr(vt, .from(l), .axpby_k).*) |k| {
        for ([_][2]T{ .{ one, beta }, .{ zero_t, beta }, .{ one, zero_t }, .{ zero_t, zero_t } }) |ab| {
            note(sweep(T, l ++ "axpby", true, true, Axpby(T, @TypeOf(k)){ .kern = k, .alpha = ab[0], .beta = ab[1] }, axpbyRun));
        }
    };
    const rot = portable.slotPtr(vt, .from(l), .rot_k).*.?;
    note(sweep(T, l ++ "rot", true, true, Rot(@TypeOf(rot)){ .kern = rot, .c = 0.6, .s = 0.8 }, rotRun));
    if (!complex) for ([_]f64{ -2, -1, 0, 1 }) |flag| {
        const k = portable.slotPtr(vt, .from(l), .rotm_k).*.?;
        note(sweep(T, l ++ "rotm", true, true, Rotm(@TypeOf(k)){ .kern = k, .flag = flag }, rotmRun));
    };
    const prec: portable.Prec = comptime .from(l);
    inline for (if (complex) .{ .{ portable.Op.dotu_k, false }, .{ portable.Op.dotc_k, true } } else .{.{ portable.Op.dot_k, false }}) |d| {
        const k = portable.slotPtr(vt, prec, d[0]).*.?;
        note(sweep(T, portable.slotName(prec, d[0]), true, true, Dot(@TypeOf(k)){ .kern = k, .conj = d[1] }, dotRun));
    }
    inline for (comptime std.enums.values(Reduce)) |op| {
        const k = portable.slotPtr(vt, prec, op.op()).*.?;
        note(sweep(T, portable.slotName(prec, op.op()), false, op == .nrm2, Reduction(@TypeOf(k)){ .kern = k, .op = op }, reduceRun));
    }
    inline for (check.extremes) |e| {
        if (comptime complex and !e[2].abs) continue;
        inline for (.{ e[0], e[1] }) |op| {
            const k = portable.slotPtr(vt, prec, op).*.?;
            note(sweep(T, portable.slotName(prec, op), false, false, Extreme(@TypeOf(k)){ .kern = k, .find = e[2] }, extremeRun));
        }
    }
}

const vtable_types = check.sd_letters;

fn checkAll(vt: *const Vtable) !void {
    failed = null;
    inline for (vtable_types) |p| slots(p[0], p[1], vt);
    note(sweep(f32, "dsdot", true, true, Dot(@TypeOf(vt.dsdot_k.?)){ .kern = vt.dsdot_k.?, .conj = false }, dotRun));
    if (failed) |e| return e;
}

const XDouble = types.XDouble;

test "portable level-1 kernels: special values, poisoned gaps, and negative and large strides" {
    const vt = check.portableVtable();
    try checkAll(vt);
    failed = null;
    if (XDouble != f64) {
        slots("q", XDouble, vt);
        slots("x", C(XDouble), vt);
    }
    if (failed) |e| return e;
}

test "every linked core's level-1 kernels: special values, poisoned gaps, and negative and large strides" {
    try check.eachCore(checkAll);
}

// ---- the same bits at any address ----------------------------------------------

/// Room for both `shifted` operands at the largest shift, 64-byte aligned so
/// each shift is the same misalignment on every run.
var shift_mem: [2 * (offset_n + 8)]C(f64) align(64) = undefined;
const offset_n = 300;

/// `kern` (a reduction of one or two vectors, by `two`) on the same data
/// shifted by 0 to 7 floats in each operand, any alignment a float may have:
/// the same bits each time, as `threading.md`'s reproducibility needs.
fn shifted(comptime T: type, routine: []const u8, comptime two: bool, kern: anytype) !void {
    const F = Float(T);
    const mem: [*]F = @ptrCast(&shift_mem);
    const w = @sizeOf(T) / @sizeOf(F);
    var vals: [2][offset_n]T = undefined;
    check.fill(T, &vals[0], 7);
    check.fill(T, &vals[1], 8);
    for ([_]usize{ 1, 3, 7, 16, 33, 67, 100, offset_n }) |n| {
        var first: ?[2]f64 = null;
        for (0..8) |ox| for (0..if (two) 8 else 1) |oy| {
            const x = mem + ox;
            const y = mem + w * (offset_n + 8) + oy;
            @memcpy(x[0 .. w * n], @as([*]const F, @ptrCast(&vals[0]))[0 .. w * n]);
            @memcpy(y[0 .. w * n], @as([*]const F, @ptrCast(&vals[1]))[0 .. w * n]);
            const got = if (two) kern(@intCast(n), x, 1, y, 1) else kern(@intCast(n), x, 1);
            const bits: [2]f64 = if (@typeInfo(@TypeOf(got)) == .@"struct") .{ got.re, got.im } else .{ got, 0 };
            if (first) |f| {
                if (!bitwise.same(f64, f[0], bits[0], true) or !bitwise.same(f64, f[1], bits[1], true)) {
                    check.say("{s} {s} n={d}: shifted by {d},{d} gave {any}, unshifted {any}: ", .{ routine, @typeName(T), n, ox, oy, bits, f });
                    return error.TestUnexpectedResult;
                }
            } else first = bits;
        };
    }
}

/// The element-wise kernels `shiftedWise` runs.
const Wise = enum { scal, axpy, axpby, rot };

/// `kern` (an element-wise kernel, by `op`, with `alpha` `a`) on `shifted`'s
/// data and shifts: the same bits in `x` and `y` each time.
fn shiftedWise(comptime T: type, routine: []const u8, comptime op: Wise, kern: anytype, a: T) !void {
    const F = Float(T);
    const complex = comptime isComplex(T);
    const mem: [*]F = @ptrCast(&shift_mem);
    const w = @sizeOf(T) / @sizeOf(F);
    var vals: [2][offset_n]T = undefined;
    check.fill(T, &vals[0], 7);
    check.fill(T, &vals[1], 8);
    const b: T = if (complex) .init(-0.5, 0.25) else -0.5;
    var first: [2][w * offset_n]F = undefined;
    for ([_]usize{ 1, 3, 7, 16, 33, 67, 100, offset_n }) |n| {
        const m = w * n;
        for (0..8) |ox| for (0..8) |oy| {
            const x = mem + ox;
            const y = mem + w * (offset_n + 8) + oy;
            @memcpy(x[0..m], @as([*]const F, @ptrCast(&vals[0]))[0..m]);
            @memcpy(y[0..m], @as([*]const F, @ptrCast(&vals[1]))[0..m]);
            const ni: BlasLong = @intCast(n);
            _ = switch (op) {
                .scal => @call(.auto, kern, .{ ni, 0, 0 } ++ util.scalarArgs(T, a) ++ .{ x, 1, x, 0, x, 1 }),
                .axpy => @call(.auto, kern, .{ ni, 0, 0 } ++ util.scalarArgs(T, a) ++ .{ x, 1, y, 1, y, 0 }),
                .axpby => @call(.auto, kern, .{ni} ++ util.scalarArgs(T, a) ++ .{ x, 1 } ++ util.scalarArgs(T, b) ++ .{ y, 1 }),
                .rot => kern(ni, x, 1, y, 1, @as(F, 0.6), @as(F, 0.8)),
            };
            if (ox == 0 and oy == 0) {
                @memcpy(first[0][0..m], x[0..m]);
                @memcpy(first[1][0..m], y[0..m]);
            } else if (!std.mem.eql(u8, std.mem.sliceAsBytes(first[0][0..m]), std.mem.sliceAsBytes(x[0..m])) or
                !std.mem.eql(u8, std.mem.sliceAsBytes(first[1][0..m]), std.mem.sliceAsBytes(y[0..m])))
            {
                check.say("{s} {s} n={d} alpha {any}: shifted by {d},{d} gave other bits than unshifted: ", .{ routine, @typeName(T), n, a, ox, oy });
                return error.TestUnexpectedResult;
            }
        };
    }
}

fn checkShifted(vt: *const Vtable) !void {
    failed = null;
    inline for (vtable_types) |p| {
        const l = p[0];
        const T = p[1];
        const prec: portable.Prec = comptime .from(l);
        inline for (comptime std.enums.values(Reduce)) |r| note(shifted(T, portable.slotName(prec, r.op()), false, portable.slotPtr(vt, prec, r.op()).*.?));
        inline for (if (comptime isComplex(T)) .{ portable.Op.dotu_k, portable.Op.dotc_k } else .{portable.Op.dot_k}) |r| note(shifted(T, portable.slotName(prec, r), true, portable.slotPtr(vt, prec, r).*.?));
        const a: T = if (comptime isComplex(T)) .init(1.5, -0.75) else 1.5;
        // complex kernels branch on a zero part (`scalReal` passes one)
        const scals: []const T = if (comptime isComplex(T)) &.{ a, .init(1.5, 0), .init(0, -0.75) } else &.{a};
        for (scals) |s| note(shiftedWise(T, l ++ "scal", .scal, portable.slotPtr(vt, .from(l), .scal_k).*.?, s));
        inline for (if (comptime isComplex(T)) .{ portable.Op.axpy_k, portable.Op.axpyc_k } else .{portable.Op.axpy_k}) |r| note(shiftedWise(T, portable.slotName(prec, r), .axpy, portable.slotPtr(vt, prec, r).*.?, a));
        if (@hasField(Vtable, portable.slotName(.from(l), .axpby_k))) if (portable.slotPtr(vt, .from(l), .axpby_k).*) |k| note(shiftedWise(T, l ++ "axpby", .axpby, k, a));
        note(shiftedWise(T, l ++ "rot", .rot, portable.slotPtr(vt, .from(l), .rot_k).*.?, a));
    }
    note(shifted(f32, "dsdot", true, vt.dsdot_k.?));
    if (failed) |e| return e;
}

test "portable level-1 reductions and element-wise kernels: the same bits whatever the address" {
    const vt = check.portableVtable();
    try checkShifted(vt);
}

test "every linked core's level-1 reductions and element-wise kernels: the same bits whatever the address" {
    try check.eachCore(checkShifted);
}

// ---- one contract on every core ----------------------------------------------

/// Long enough for every `@Vector` kernel's blocks (`4 * W` keys, 256 for
/// `f32` on AVX-512) to hold several.
const contract_n = 600;

/// Room for a `contract` operand.
var contract_mem: [contract_n]C(f64) = undefined;

/// Where `contract` puts special values in `x`: NaN at `nan`, and again at
/// `nan2`, infinities at `inf`, by element; null for none. `part` puts them
/// in the imaginary part of a complex element; `big` makes every element
/// with no special value `floatMax / 2` in each part, so a complex key
/// overflows.
const Specials = struct { nan: ?usize = null, nan2: ?usize = null, inf: ?usize = null, part: bool = false, big: bool = false };

/// Index and value routines (`i?amax`, `?amax`, ...) and `nrm2` against the
/// portable kernels on random data of every length to `contract_n`, with
/// NaN and infinities at the first, a middle and the last element, inside and
/// past the vector blocks: the same index, NaN or infinity where they give
/// one, and values within their accuracy class.
fn contract(vt: *const Vtable) !void {
    const ref = check.portableVtable();
    failed = null;
    inline for (vtable_types) |p| note(contractOf(p[0], p[1], vt, ref));
    if (failed) |e| return e;
}

fn contractOf(comptime l: []const u8, comptime T: type, vt: *const Vtable, ref: *const Vtable) !void {
    const F = Float(T);
    const complex = comptime isComplex(T);
    const w = if (complex) 2 else 1;
    const x: [*]F = @ptrCast(&contract_mem);
    for ([_]usize{ 1, 2, 9, 63, 257, 300, contract_n }) |n| {
        const mid = n / 2;
        const cases = [_]Specials{
            .{},
            .{ .nan = 0 },
            .{ .nan = mid },
            .{ .nan = n - 1 },
            .{ .nan = mid, .part = true },
            .{ .nan = mid, .nan2 = n - 1, .inf = 0 },
            .{ .inf = mid },
            .{ .inf = n - 1, .nan = n - 1, .part = true },
            .{ .big = true },
            .{ .big = true, .inf = n - 1, .part = true },
        };
        for (cases) |sp| {
            check.fill(F, x[0 .. w * n], n + 3);
            for (0..n) |i| {
                const at = w * i + @intFromBool(sp.part and complex);
                if (sp.big) for (0..w) |c| {
                    x[w * i + c] = std.math.copysign(std.math.floatMax(F) / 2, x[w * i + c]);
                };
                if (sp.inf == i) x[at] = -std.math.inf(F);
                if (sp.nan == i or sp.nan2 == i) x[at] = std.math.nan(F);
            }
            contractCase(l, T, vt, ref, n, x) catch |e| {
                check.say("{s} n={d} {any}: ", .{ @typeName(T), n, sp });
                return e;
            };
        }
    }
}

fn contractCase(comptime l: []const u8, comptime T: type, vt: *const Vtable, ref: *const Vtable, n: usize, x: [*]Float(T)) !void {
    const F = Float(T);
    const complex = comptime isComplex(T);
    const ni: BlasLong = @intCast(n);
    const prec: portable.Prec = comptime .from(l);
    inline for (check.extremes) |e| {
        if (comptime complex and !e[2].abs) continue;
        const i_got = portable.slotPtr(vt, prec, e[1]).*.?(ni, x, 1);
        const i_want = portable.slotPtr(ref, prec, e[1]).*.?(ni, x, 1);
        if (i_got != i_want) {
            check.say("{s}: {d}, the portable kernel {d}: ", .{ portable.slotName(prec, e[1]), i_got, i_want });
            return error.TestUnexpectedResult;
        }
        const got = portable.slotPtr(vt, prec, e[0]).*.?(ni, x, 1);
        const want = portable.slotPtr(ref, prec, e[0]).*.?(ni, x, 1);
        if (!(got == want or (std.math.isNan(got) and std.math.isNan(want)))) {
            check.say("{s}: {e}, the portable kernel {e}: ", .{ portable.slotName(prec, e[0]), got, want });
            return error.TestUnexpectedResult;
        }
    }
    const got = portable.slotPtr(vt, .from(l), .nrm2_k).*.?(ni, x, 1);
    const want = portable.slotPtr(ref, .from(l), .nrm2_k).*.?(ni, x, 1);
    const w: f128 = if (complex) 2 else 1;
    const ulps: f128 = if (F == f32) 2 else @as(f128, @floatFromInt(n)) * w + 4;
    if (std.math.isNan(want) != std.math.isNan(got) or std.math.isInf(want) != std.math.isInf(got) or
        (std.math.isFinite(want) and @abs(@as(f128, got) - want) > ulps * std.math.floatEps(F) * @as(f128, want)))
    {
        check.say("{s}nrm2: {e}, the portable kernel {e}: ", .{ l, got, want });
        return error.TestUnexpectedResult;
    }
}

test "every linked core keeps the portable kernels' contract: the extremes' NaN and index, nrm2's class" {
    try check.eachCore(contract);
}
