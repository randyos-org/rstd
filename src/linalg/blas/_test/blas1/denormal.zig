//! Denormal results through every level-1 kernel on every core this machine
//! runs, against the `f128` reference; ARMv7 NEON and pre-VSX AltiVec flush
//! `f32` lanes (`util.vectorIeee`).

const std = @import("std");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");
const check = @import("../check.zig");
const bitwise = @import("../bitwise.zig");
const pattern = @import("../pattern.zig");
const reference = @import("../reference.zig");
const extremum = @import("../../blas1/extremum.zig");
const util = @import("../../util.zig");

const Allocator = std.mem.Allocator;
const BlasLong = types.BlasLong;
const C = std.math.Complex;
const Float = check.FloatOf;
const isComplex = check.isComplex;
const Wide = reference.Wide;
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
        const x_ok = if (ctx.swap) bitwise.same(T, v[1].was[i], v[0].at(i), true) else true;
        if (!x_ok or !bitwise.same(T, v[0].was[i], v[1].at(i), true)) {
            check.say("element {d} not moved exactly: ", .{i});
            return error.TestUnexpectedResult;
        }
    }
}

fn scalRun(comptime T: type, ctx: anytype, al: Allocator, n: usize, inc: [2]isize) !void {
    var v = try operands(T, al, n, inc, true);
    const x = &v[0];
    const a = ctx.alpha;
    _ = @call(.auto, ctx.kern, .{ @as(BlasLong, @intCast(n)), 0, 0 } ++ util.scalarArgs(T, a) ++ .{ x.first(), stride(x.*), x.first(), 0, x.first(), 1 });
    const wa = reference.wide(T, a);
    for (0..n) |i| try pattern.expectWide(T, reference.mulW(wa, x.w(i)), x.at(i), reference.elem_ulps, reference.mag(wa) * reference.mag(x.w(i)));
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
        _ = @call(.auto, ctx.kern, .{nl} ++ util.scalarArgs(T, a) ++ .{ x.first(), stride(x.*) } ++ util.scalarArgs(T, b) ++ .{ y.first(), stride(y.*) });
    } else {
        _ = @call(.auto, ctx.kern, .{ nl, 0, 0 } ++ util.scalarArgs(T, a) ++ .{ x.first(), stride(x.*), y.first(), stride(y.*), y.first(), 0 });
    }
    const wa = reference.wide(T, a);
    const wb = reference.wide(T, b);
    for (0..n) |i| {
        const xv = if (ctx.conj) reference.conjW(x.w(i)) else x.w(i);
        const want = reference.addW(reference.mulW(wa, xv), reference.mulW(wb, y.w(i)));
        const scale = reference.mag(wa) * reference.mag(xv) + reference.mag(wb) * reference.mag(y.w(i));
        try pattern.expectWide(T, want, y.at(i), reference.elem_ulps, scale);
    }
}

fn lin(comptime T: type, a: f128, u: Wide(T), b: f128, w: Wide(T)) Wide(T) {
    return reference.addW(reference.scaleW(a, u), reference.scaleW(b, w));
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
        const scale = reference.mag(x.w(i)) + reference.mag(y.w(i));
        try pattern.expectWide(T, lin(T, h[0], x.w(i), h[1], y.w(i)), x.at(i), reference.elem_ulps, scale);
        try pattern.expectWide(T, lin(T, h[2], x.w(i), h[3], y.w(i)), y.at(i), reference.elem_ulps, scale);
    }
}

/// `x . y`, `x` conjugated for `conj`, `y` near 1.
fn dotRun(comptime T: type, ctx: anytype, al: Allocator, n: usize, inc: [2]isize) !void {
    var v = try operands(T, al, n, inc, false);
    const got = ctx.kern(@intCast(n), v[0].first(), stride(v[0]), v[1].first(), stride(v[1]));
    var want = reference.wide(T, if (comptime isComplex(T)) T.init(0, 0) else 0);
    var scale: f128 = 0;
    for (0..n) |i| {
        const xv = if (ctx.conj) reference.conjW(v[0].w(i)) else v[0].w(i);
        want = reference.addW(want, reference.mulW(xv, v[1].w(i)));
        scale += reference.mag(xv) * reference.mag(v[1].w(i));
    }
    const ulps = reference.dotUlps(n);
    if (comptime isComplex(T)) return pattern.expectWide(T, want, .init(got.re, got.im), ulps, scale);
    try pattern.component(T, want, got, ulps, scale);
}

const Reduce = check.Reduce;

fn reduceRun(comptime T: type, ctx: anytype, al: Allocator, n: usize, inc: [2]isize) !void {
    var v = try operands(T, al, n, inc, true);
    const x = &v[0];
    const got = ctx.kern(@intCast(n), x.first(), stride(x.*));
    var acc: f128 = 0;
    var mag: f128 = 0;
    for (0..n) |i| {
        const w = x.w(i);
        mag += reference.mag(w);
        acc += switch (ctx.op) {
            .asum => reference.abs1(w),
            .sum => if (comptime isComplex(T)) w.re + w.im else w,
            .nrm2 => if (comptime isComplex(T)) w.re * w.re + w.im * w.im else w * w,
        };
    }
    const want = if (ctx.op == .nrm2) reference.sqrtW(acc) else acc;
    try pattern.component(Float(T), want, got, reference.sumUlps(n), if (ctx.op == .nrm2) want else mag);
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
        const k = if (comptime isComplex(T)) reference.abs1(w) else if (ctx.find.abs) @abs(w) else w;
        if (i == 0 or (if (ctx.find.max) k > key else k < key)) {
            best = i;
            key = k;
        }
    }
    const ok = if (@TypeOf(got) == BlasLong) got == best + 1 else got == key or
        @abs(@as(f128, got) - key) <= reference.keyTol(Float(T)) * @abs(key);
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
    const prec: portable.Prec = comptime .from(l);
    inline for (.{ .{ portable.Op.copy_k, false }, .{ portable.Op.swap_k, true } }) |m| {
        note(sweep(T, portable.slotName(prec, m[0]), .{ .kern = portable.slotPtr(vt, prec, m[0]).*.?, .swap = m[1] }, moveRun));
    }
    note(sweep(T, l ++ "scal", .{ .kern = portable.slotPtr(vt, .from(l), .scal_k).*.?, .alpha = one }, scalRun));
    inline for (if (complex) .{ portable.Op.axpy_k, portable.Op.axpyc_k } else .{portable.Op.axpy_k}) |r| {
        // a zero alpha too: `y + 0*x` is `y`
        for ([_]T{ one, zero }) |a| {
            const ctx = .{ .kern = portable.slotPtr(vt, prec, r).*.?, .alpha = a, .beta = unit, .conj = r == .axpyc_k, .axpby = false };
            note(sweep(T, portable.slotName(prec, r), ctx, axpyRun));
        }
    }
    if (@hasField(Vtable, portable.slotName(.from(l), .axpby_k))) if (portable.slotPtr(vt, .from(l), .axpby_k).*) |k| {
        for ([_][2]T{ .{ one, beta }, .{ zero, beta }, .{ one, zero } }) |ab| {
            note(sweep(T, l ++ "axpby", .{ .kern = k, .alpha = ab[0], .beta = ab[1], .conj = false, .axpby = true }, axpyRun));
        }
    };
    note(sweep(T, portable.slotName(.from(l), .rot_k), .{ .kern = portable.slotPtr(vt, .from(l), .rot_k).*.? }, rotRun));
    if (!complex) for ([_]f64{ -1, 0, 1 }) |flag| {
        note(sweep(T, l ++ "rotm", .{ .kern = portable.slotPtr(vt, .from(l), .rotm_k).*.?, .flag = flag }, rotRun));
    };
    inline for (if (complex) .{ .{ portable.Op.dotu_k, false }, .{ portable.Op.dotc_k, true } } else .{.{ portable.Op.dot_k, false }}) |d| {
        note(sweep(T, portable.slotName(prec, d[0]), .{ .kern = portable.slotPtr(vt, prec, d[0]).*.?, .conj = d[1] }, dotRun));
    }
    inline for (comptime std.enums.values(Reduce)) |op| {
        note(sweep(T, portable.slotName(prec, op.op()), .{ .kern = portable.slotPtr(vt, prec, op.op()).*.?, .op = op }, reduceRun));
    }
    inline for (check.extremes) |e| {
        if (comptime complex and !e[2].abs) continue;
        inline for (.{ e[0], e[1] }) |op| {
            note(sweep(T, portable.slotName(prec, op), .{ .kern = portable.slotPtr(vt, prec, op).*.?, .find = e[2] }, extremeRun));
        }
    }
}

fn checkAll(vt: *const Vtable) !void {
    failed = null;
    inline for (check.sd_letters) |p| slots(p[0], p[1], vt);
    note(sweep(f32, "dsdot", .{ .kern = vt.dsdot_k.?, .conj = false }, dotRun));
    if (failed) |e| return e;
}

test "portable level-1 kernels: denormal results" {
    const vt = check.portableVtable();
    try checkAll(vt);
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

test "where util.vectorIeee fails, f32 lanes flush denormals" {
    if (comptime util.vectorIeee(f32)) return error.SkipZigTest;
    const tiny = std.math.floatTrueMin(f32);
    // denormal inputs, and a denormal product of normal ones
    var a: @Vector(4, f32) = .{ tiny, 3 * tiny, -tiny, 0x1p-126 };
    var b: @Vector(4, f32) = .{ 2, 2, 2, 0.5 };
    std.mem.doNotOptimizeAway(&a);
    std.mem.doNotOptimizeAway(&b);
    if (!@reduce(.And, a * b == @as(@Vector(4, f32), @splat(0)))) return error.TestUnexpectedResult;
}
