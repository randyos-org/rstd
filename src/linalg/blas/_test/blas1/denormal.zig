//! Denormal results through every level-1 kernel, on every core this machine
//! runs: the operands denormal and the scalars near 1, so a flushed input or
//! result shows against the `f128` reference. ARMv7 NEON and pre-VSX AltiVec
//! flush `f32` lanes (`util.vectorIeee`).

const std = @import("std");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("check.zig");
const extremum = @import("../../blas1/extremum.zig");
const util = @import("../../util.zig");

const Allocator = std.mem.Allocator;
const BlasLong = types.BlasLong;
const C = std.math.Complex;
const Float = check.FloatOf;
const isComplex = check.isComplex;
const Wide = check.Wide;
const Vector = check.Vector;

/// Past every kernel's blocks, the extremes' `4 * W` keys too.
const sizes = [_]usize{ 1, 7, 33, 67, 300 };
const max_n = 300;
const strides = [_][2]isize{ .{ 1, 1 }, .{ 2, 3 } };

/// Two operands at the largest size and stride.
var mem: [2 * ((max_n - 1) * 3 + 1 + check.tail + max_n) * @sizeOf(C(f64))]u8 align(16) = undefined;

/// `x` denormal, `y` denormal too or, for `y_denormal` false, near 1.
fn operands(comptime T: type, al: Allocator, n: usize, inc: [2]isize, y_denormal: bool) ![2]Vector(T) {
    return .{
        try .init(al, n, inc[0], 31 * n, .denormal),
        try .init(al, n, inc[1], 31 * n + 1, if (y_denormal) .denormal else null),
    };
}

fn stride(v: anytype) BlasLong {
    return @intCast(v.inc);
}

/// Runs `run(T, ctx, al, n, inc)` over every size and stride; names the case
/// that fails.
fn sweep(comptime T: type, routine: []const u8, ctx: anytype, comptime run: anytype) !void {
    var fba: std.heap.FixedBufferAllocator = .init(&mem);
    for (strides) |inc| for (sizes) |n| {
        fba.reset();
        run(T, ctx, fba.allocator(), n, inc) catch |e| {
            check.say("{s} {s} n={d} inc={d},{d}: ", .{ routine, @typeName(T), n, inc[0], inc[1] });
            return e;
        };
    };
}

fn moveRun(comptime T: type, ctx: anytype, al: Allocator, n: usize, inc: [2]isize) !void {
    var v = try operands(T, al, n, inc, true);
    const nl: BlasLong = @intCast(n);
    if (ctx.swap) {
        if (comptime isComplex(T)) {
            _ = ctx.kern(nl, 0, 0, 0, 0, v[0].first(), stride(v[0]), v[1].first(), stride(v[1]), v[1].first(), 0);
        } else {
            _ = ctx.kern(nl, 0, 0, 0, v[0].first(), stride(v[0]), v[1].first(), stride(v[1]), v[1].first(), 0);
        }
    } else _ = ctx.kern(nl, v[0].first(), stride(v[0]), v[1].first(), stride(v[1]));
    for (0..n) |i| {
        const x_ok = if (ctx.swap) check.same(T, v[1].was[i], v[0].at(i), true) else true;
        if (!x_ok or !check.same(T, v[0].was[i], v[1].at(i), true)) {
            check.say("element {d} not moved exactly: ", .{i});
            return error.TestUnexpectedResult;
        }
    }
}

fn scalRun(comptime T: type, ctx: anytype, al: Allocator, n: usize, inc: [2]isize) !void {
    var v = try operands(T, al, n, inc, true);
    const x = &v[0];
    const a = ctx.alpha;
    if (comptime isComplex(T)) {
        _ = ctx.kern(@intCast(n), 0, 0, a.re, a.im, x.first(), stride(x.*), x.first(), 0, x.first(), 1);
    } else {
        _ = ctx.kern(@intCast(n), 0, 0, a, x.first(), stride(x.*), x.first(), 0, x.first(), 1);
    }
    const wa = check.wide(T, a);
    for (0..n) |i| try check.expectWide(T, check.mulW(wa, x.w(i)), x.at(i), 4, check.mag(wa) * check.mag(x.w(i)));
}

/// `axpy` (`axpyc` for `conj`) and, with a `beta`, `axpby`.
fn axpyRun(comptime T: type, ctx: anytype, al: Allocator, n: usize, inc: [2]isize) !void {
    var v = try operands(T, al, n, inc, true);
    const x = &v[0];
    const y = &v[1];
    const a = ctx.alpha;
    const b = ctx.beta;
    const nl: BlasLong = @intCast(n);
    if (ctx.axpby) {
        if (comptime isComplex(T)) {
            _ = ctx.kern(nl, a.re, a.im, x.first(), stride(x.*), b.re, b.im, y.first(), stride(y.*));
        } else {
            _ = ctx.kern(nl, a, x.first(), stride(x.*), b, y.first(), stride(y.*));
        }
    } else if (comptime isComplex(T)) {
        _ = ctx.kern(nl, 0, 0, a.re, a.im, x.first(), stride(x.*), y.first(), stride(y.*), y.first(), 0);
    } else {
        _ = ctx.kern(nl, 0, 0, a, x.first(), stride(x.*), y.first(), stride(y.*), y.first(), 0);
    }
    const wa = check.wide(T, a);
    const wb = check.wide(T, b);
    for (0..n) |i| {
        const xv = if (ctx.conj) check.conjW(x.w(i)) else x.w(i);
        const want = check.addW(check.mulW(wa, xv), check.mulW(wb, y.w(i)));
        const scale = check.mag(wa) * check.mag(xv) + check.mag(wb) * check.mag(y.w(i));
        try check.expectWide(T, want, y.at(i), 4, scale);
    }
}

fn lin(comptime T: type, a: f128, u: Wide(T), b: f128, w: Wide(T)) Wide(T) {
    return check.addW(check.scaleW(a, u), check.scaleW(b, w));
}

/// `rot` by `c` 0.6 and `s` 0.8, or `rotm` for a `flag`.
fn rotRun(comptime T: type, ctx: anytype, al: Allocator, n: usize, inc: [2]isize) !void {
    var v = try operands(T, al, n, inc, true);
    const x = &v[0];
    const y = &v[1];
    const F = Float(T);
    // h11, h12, h21, h22
    var h: [4]f128 = .{ 0.6, 0.8, -0.8, 0.6 };
    if (comptime @hasField(@TypeOf(ctx), "flag")) {
        const flag = ctx.flag;
        var param = [5]F{ @floatCast(flag), 0.75, -0.5, 0.25, 1.25 };
        h = if (flag < 0) .{ 0.75, 0.25, -0.5, 1.25 } else if (flag == 0) .{ 1, 0.25, -0.5, 1 } else .{ 0.75, 1, -1, 1.25 };
        _ = ctx.kern(@intCast(n), @ptrCast(x.base()), stride(x.*), @ptrCast(y.base()), stride(y.*), &param);
    } else {
        _ = ctx.kern(@intCast(n), x.first(), stride(x.*), y.first(), stride(y.*), @as(F, 0.6), @as(F, 0.8));
    }
    for (0..n) |i| {
        const scale = check.mag(x.w(i)) + check.mag(y.w(i));
        try check.expectWide(T, lin(T, h[0], x.w(i), h[1], y.w(i)), x.at(i), 4, scale);
        try check.expectWide(T, lin(T, h[2], x.w(i), h[3], y.w(i)), y.at(i), 4, scale);
    }
}

/// `x . y`, `x` conjugated for `conj`, `y` near 1.
fn dotRun(comptime T: type, ctx: anytype, al: Allocator, n: usize, inc: [2]isize) !void {
    var v = try operands(T, al, n, inc, false);
    const got = ctx.kern(@intCast(n), v[0].first(), stride(v[0]), v[1].first(), stride(v[1]));
    var want = check.wide(T, if (comptime isComplex(T)) T.init(0, 0) else 0);
    var scale: f128 = 0;
    for (0..n) |i| {
        const xv = if (ctx.conj) check.conjW(v[0].w(i)) else v[0].w(i);
        want = check.addW(want, check.mulW(xv, v[1].w(i)));
        scale += check.mag(xv) * check.mag(v[1].w(i));
    }
    const ulps: f128 = @floatFromInt(4 * (n + 2));
    if (comptime isComplex(T)) return check.expectWide(T, want, .init(got.re, got.im), ulps, scale);
    try check.component(T, want, got, ulps, scale);
}

const Reduce = enum { asum, sum, nrm2 };

fn reduceRun(comptime T: type, ctx: anytype, al: Allocator, n: usize, inc: [2]isize) !void {
    var v = try operands(T, al, n, inc, true);
    const x = &v[0];
    const got = ctx.kern(@intCast(n), x.first(), stride(x.*));
    var acc: f128 = 0;
    var mag: f128 = 0;
    for (0..n) |i| {
        const w = x.w(i);
        mag += check.mag(w);
        acc += switch (ctx.op) {
            .asum => check.abs1(w),
            .sum => if (comptime isComplex(T)) w.re + w.im else w,
            .nrm2 => if (comptime isComplex(T)) w.re * w.re + w.im * w.im else w * w,
        };
    }
    const want = if (ctx.op == .nrm2) check.sqrtW(acc) else acc;
    try check.component(Float(T), want, got, @floatFromInt(2 * (n + 2)), if (ctx.op == .nrm2) want else mag);
}

/// The extreme's value, or its 1-based index, by the kernel's return type.
fn extremeRun(comptime T: type, ctx: anytype, al: Allocator, n: usize, inc: [2]isize) !void {
    var v = try operands(T, al, n, inc, true);
    const x = &v[0];
    const got = ctx.kern(@intCast(n), x.first(), stride(x.*));
    var best: usize = 0;
    var key: f128 = 0;
    for (0..n) |i| {
        const w = x.w(i);
        const k = if (comptime isComplex(T)) check.abs1(w) else if (ctx.find.abs) @abs(w) else w;
        if (i == 0 or (if (ctx.find.max) k > key else k < key)) {
            best = i;
            key = k;
        }
    }
    const ok = if (@TypeOf(got) == BlasLong) got == best + 1 else got == key or
        @abs(@as(f128, got) - key) <= 2 * std.math.floatEps(Float(T)) * @abs(key);
    if (ok) return;
    check.say("{any}, want index {d}, key {e}: ", .{ got, best + 1, key });
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
    const one: T = if (complex) .init(1.5, -0.75) else 1.5;
    const beta: T = if (complex) .init(-0.5, 0.25) else -0.5;
    const zero: T = if (complex) .init(0, 0) else 0;
    const unit: T = if (complex) .init(1, 0) else 1;
    inline for (.{ .{ "copy", false }, .{ "swap", true } }) |m| {
        note(sweep(T, l ++ m[0], .{ .kern = @field(vt, l ++ m[0] ++ "_k").?, .swap = m[1] }, moveRun));
    }
    note(sweep(T, l ++ "scal", .{ .kern = @field(vt, l ++ "scal_k").?, .alpha = one }, scalRun));
    inline for (if (complex) .{ "axpy", "axpyc" } else .{"axpy"}) |r| {
        // a zero alpha too: `y + 0*x` is `y`
        for ([_]T{ one, zero }) |a| {
            const ctx = .{ .kern = @field(vt, l ++ r ++ "_k").?, .alpha = a, .beta = unit, .conj = r.len == 5, .axpby = false };
            note(sweep(T, l ++ r, ctx, axpyRun));
        }
    }
    if (@hasField(Vtable, l ++ "axpby_k")) if (@field(vt, l ++ "axpby_k")) |k| {
        for ([_][2]T{ .{ one, beta }, .{ zero, beta }, .{ one, zero } }) |ab| {
            note(sweep(T, l ++ "axpby", .{ .kern = k, .alpha = ab[0], .beta = ab[1], .conj = false, .axpby = true }, axpyRun));
        }
    };
    const rot = l ++ (if (complex) (if (l[0] == 'c') "s" else "d") else "") ++ "rot_k";
    note(sweep(T, rot, .{ .kern = @field(vt, rot).? }, rotRun));
    if (!complex) for ([_]f64{ -1, 0, 1 }) |flag| {
        note(sweep(T, l ++ "rotm", .{ .kern = @field(vt, l ++ "rotm_k").?, .flag = flag }, rotRun));
    };
    inline for (if (complex) .{ .{ "dotu", false }, .{ "dotc", true } } else .{.{ "dot", false }}) |d| {
        note(sweep(T, l ++ d[0], .{ .kern = @field(vt, l ++ d[0] ++ "_k").?, .conj = d[1] }, dotRun));
    }
    inline for (comptime std.enums.values(Reduce)) |op| {
        note(sweep(T, l ++ @tagName(op), .{ .kern = @field(vt, l ++ @tagName(op) ++ "_k").?, .op = op }, reduceRun));
    }
    inline for (.{ "amax", "amin", "max", "min" }) |name| {
        if (comptime complex and name[0] != 'a') continue;
        const find: extremum.Find = .{ .abs = name[0] == 'a', .max = comptime std.mem.endsWith(u8, name, "max") };
        inline for (.{ "", "i" }) |i| {
            note(sweep(T, i ++ l ++ name, .{ .kern = @field(vt, i ++ l ++ name ++ "_k").?, .find = find }, extremeRun));
        }
    }
}

fn checkAll(vt: *const Vtable) !void {
    failed = null;
    inline for (.{ .{ "s", f32 }, .{ "d", f64 }, .{ "c", C(f32) }, .{ "z", C(f64) } }) |p| slots(p[0], p[1], vt);
    note(sweep(f32, "dsdot", .{ .kern = vt.dsdot_k.?, .conj = false }, dotRun));
    if (failed) |e| return e;
}

test "portable level-1 kernels: denormal results" {
    var vt: Vtable = .{};
    @import("../../blas1/__root__.zig").installGeneric(&vt);
    try checkAll(&vt);
}

test "every linked core's level-1 kernels: denormal results" {
    try check.eachCore(checkAll);
}

test "where util.vectorIeee holds, f32 lanes keep denormals" {
    if (comptime !util.vectorIeee(f32)) return error.SkipZigTest;
    const tiny = std.math.floatTrueMin(f32);
    var a: @Vector(4, f32) = .{ tiny, 3 * tiny, -tiny, 0x1p-130 };
    var b: @Vector(4, f32) = @splat(2);
    std.mem.doNotOptimizeAway(&a);
    std.mem.doNotOptimizeAway(&b);
    const want: @Vector(4, f32) = .{ 3 * tiny, 9 * tiny, -3 * tiny, 0x3p-130 };
    if (!@reduce(.And, a * b + a == want)) return error.TestUnexpectedResult;
}
