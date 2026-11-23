//! Edge cases through every level-1 kernel on every core this machine runs,
//! against an `f128` reference: `pattern.Pattern`'s special values, NaN in the
//! gaps, negative and large strides; what may not propagate, `design.md`'s Verification.
//! Denormal results too, unflushed, each in tests of its own; ARMv7 NEON and
//! pre-VSX AltiVec flush `f32` lanes (`util.vectorIeee`).

const std = @import("std");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");
const util = @import("../../util.zig");
const check = @import("../check.zig");
const pattern = @import("../pattern.zig");
const reference = @import("../reference.zig");
const er = @import("edge/runs.zig");
const es = @import("edge/shifted.zig");
const ec = @import("edge/contract.zig");

const Allocator = std.mem.Allocator;
const BlasLong = types.BlasLong;
const C = std.math.Complex;
const Float = check.FloatOf;
const isComplex = check.isComplex;
const Wide = reference.Wide;
const Vector = check.Vector;
const Axpby = er.Axpby;
const Axpy = er.Axpy;
const Dot = er.Dot;
const Extreme = er.Extreme;
const Reduce = er.Reduce;
const Reduction = er.Reduction;
const Rot = er.Rot;
const Rotm = er.Rotm;
const Scal = er.Scal;
const axpbyRun = er.axpbyRun;
const axpyRun = er.axpyRun;
const copyRun = er.copyRun;
const dotRun = er.dotRun;
const extremeRun = er.extremeRun;
const reduceRun = er.reduceRun;
const rotRun = er.rotRun;
const rotmRun = er.rotmRun;
const scalRun = er.scalRun;
const swapRun = er.swapRun;
const checkShifted = es.checkShifted;
const contract = ec.contract;

/// Sizes past every unroll, for strides up to 3.
pub const sizes = [_]usize{ 1, 2, 5, 16, 33, 67 };
pub const max_n = 67;
pub const strides = [_][2]isize{ .{ 1, 1 }, .{ 2, 3 } };
/// For routines whose kernels take negative strides: the interface passes the
/// far end and the stride as given (`rotm` the near end; its kernel finds the
/// far one).
pub const negative = [_][2]isize{ .{ -1, -1 }, .{ -3, 2 }, .{ 2, -3 } };
/// Strides of a page and more, on fewer sizes and one pattern.
pub const large = [_][2]isize{ .{ 4097, 1 }, .{ -65537, 65537 } };
pub const large_sizes = [_]usize{ 2, 5 };
/// The denormal sweep's: past every kernel's blocks, the extremes' `4 * W`
/// keys too.
pub const denormal_sizes = [_]usize{ 1, 7, 33, 67, 300 };

/// What a sweep puts in the operands: `pattern.Pattern`'s special values one
/// operand at a time, or denormals in every operand (`dot`'s `y` near 1).
pub const Mode = enum { specials, denormals };

/// A case: size, strides, pattern, the operand that has it (`both`: every
/// operand), and where the operands' storage comes from.
pub const Case = struct {
    n: usize,
    inc: [2]isize,
    p: pattern.Pattern,
    on_y: bool,
    both: bool = false,
    al: Allocator,

    /// `c` with the pattern in its one operand only.
    pub fn alone(c: Case) Case {
        var d = c;
        d.both = false;
        return d;
    }
};

/// Operand `which` of `c`: 0 for `x`, 1 for `y`.
pub fn vec(comptime T: type, c: Case, which: u1) !Vector(T) {
    return .init(c.al, c.n, c.inc[which], 31 * c.n + which, if (c.both or (which == 1) == c.on_y) c.p else null);
}

pub fn stride(v: anytype) BlasLong {
    return @intCast(v.inc);
}

pub fn zero(comptime T: type) Wide(T) {
    return reference.wide(T, if (comptime isComplex(T)) T.init(0, 0) else 0);
}

pub const isZero = util.isZero;

/// `a*u + b*v` for real `a` and `b`.
pub fn lin(comptime T: type, a: f128, u: Wide(T), b: f128, v: Wide(T)) Wide(T) {
    if (comptime isComplex(T)) return .init(a * u.re + b * v.re, a * u.im + b * v.im);
    return a * u + b * v;
}

/// Storage for the small cases, so the rigs do not spend their time in the
/// testing allocator: four operands of the most room any takes.
pub var small_mem: [4 * ((max_n - 1) * 3 + 1 + check.tail + max_n) * @sizeOf(C(f128))]u8 align(16) = undefined;

/// Runs `run(T, ctx, case)` over every size, stride and pattern, the pattern
/// in `x` and, for `two`, then in `y`; with `neg`, also over negative strides,
/// and over large ones always. For `.denormals`, over `denormal_sizes` and
/// the positive strides, denormals in both operands for `two`. Names the
/// case that fails.
pub fn sweep(comptime T: type, mode: Mode, routine: []const u8, two: bool, neg: bool, ctx: anytype, comptime run: anytype) !void {
    var fba: std.heap.FixedBufferAllocator = .init(&small_mem);
    if (mode == .denormals) {
        for (strides) |inc| for (denormal_sizes) |n| {
            fba.reset();
            const c: Case = .{ .n = n, .inc = inc, .p = .denormal, .on_y = false, .both = two, .al = if (n > max_n) std.testing.allocator else fba.allocator() };
            run(T, ctx, c) catch |e| {
                check.say("denormal {s} {s} n={d} inc={d},{d}: ", .{ routine, @typeName(T), n, inc[0], inc[1] });
                return e;
            };
        };
        return;
    }
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

/// The first failure of each routine, so one run names every one (`check.keep`).
pub var failed: ?anyerror = null;

/// Every routine of `T`, whose slots begin with `l`.
pub fn slots(comptime l: []const u8, comptime T: type, mode: Mode, vt: *const Vtable) void {
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
    check.keep(&failed, sweep(T, mode, l ++ "copy", false, true, portable.slotPtr(vt, .from(l), .copy_k).*.?, copyRun));
    check.keep(&failed, sweep(T, mode, l ++ "swap", true, true, portable.slotPtr(vt, .from(l), .swap_k).*.?, swapRun));
    const scal = portable.slotPtr(vt, .from(l), .scal_k).*.?;
    for ([_]T{ one, inf, nan_a, nan_b }) |a| check.keep(&failed, sweep(T, mode, l ++ "scal", false, false, Scal(T, @TypeOf(scal)){ .kern = scal, .alpha = a }, scalRun));
    const axpy = portable.slotPtr(vt, .from(l), .axpy_k).*.?;
    // a zero alpha too, which the triangular drivers pass (`-x(j)`): `0 * x`
    for ([_]T{ one, zero_t, nan_a, nan_b }) |a| {
        check.keep(&failed, sweep(T, mode, l ++ "axpy", true, true, Axpy(T, @TypeOf(axpy)){ .kern = axpy, .alpha = a, .conj = false }, axpyRun));
        if (complex) check.keep(&failed, sweep(T, mode, l ++ "axpyc", true, true, Axpy(T, @TypeOf(axpy)){ .kern = portable.slotPtr(vt, .from(l), .axpyc_k).*.?, .alpha = a, .conj = true }, axpyRun));
    }
    if (@hasField(Vtable, portable.slotName(.from(l), .axpby_k))) if (portable.slotPtr(vt, .from(l), .axpby_k).*) |k| {
        for ([_][2]T{ .{ one, beta }, .{ zero_t, beta }, .{ one, zero_t }, .{ zero_t, zero_t } }) |ab| {
            check.keep(&failed, sweep(T, mode, l ++ "axpby", true, true, Axpby(T, @TypeOf(k)){ .kern = k, .alpha = ab[0], .beta = ab[1] }, axpbyRun));
        }
    };
    const rot = portable.slotPtr(vt, .from(l), .rot_k).*.?;
    check.keep(&failed, sweep(T, mode, l ++ "rot", true, true, Rot(@TypeOf(rot)){ .kern = rot, .c = 0.6, .s = 0.8 }, rotRun));
    if (!complex) for ([_]f64{ -2, -1, 0, 1 }) |flag| {
        const k = portable.slotPtr(vt, .from(l), .rotm_k).*.?;
        check.keep(&failed, sweep(T, mode, l ++ "rotm", true, true, Rotm(@TypeOf(k)){ .kern = k, .flag = flag }, rotmRun));
    };
    const prec: portable.Prec = comptime .from(l);
    inline for (if (complex) .{ .{ portable.Op.dotu_k, false }, .{ portable.Op.dotc_k, true } } else .{.{ portable.Op.dot_k, false }}) |d| {
        const k = portable.slotPtr(vt, prec, d[0]).*.?;
        check.keep(&failed, sweep(T, mode, portable.slotName(prec, d[0]), true, true, Dot(@TypeOf(k)){ .kern = k, .conj = d[1] }, dotRun));
    }
    inline for (comptime std.enums.values(Reduce)) |op| {
        const k = portable.slotPtr(vt, prec, op.op()).*.?;
        check.keep(&failed, sweep(T, mode, portable.slotName(prec, op.op()), false, op == .nrm2, Reduction(@TypeOf(k)){ .kern = k, .op = op }, reduceRun));
    }
    inline for (check.extremes) |e| {
        if (comptime complex and !e[2].abs) continue;
        inline for (.{ e[0], e[1] }) |op| {
            const k = portable.slotPtr(vt, prec, op).*.?;
            check.keep(&failed, sweep(T, mode, portable.slotName(prec, op), false, false, Extreme(@TypeOf(k)){ .kern = k, .find = e[2] }, extremeRun));
        }
    }
}

pub const vtable_types = check.sd_letters;

pub fn checkMode(vt: *const Vtable, mode: Mode) !void {
    failed = null;
    inline for (vtable_types) |p| slots(p[0], p[1], mode, vt);
    check.keep(&failed, sweep(f32, mode, "dsdot", true, true, Dot(@TypeOf(vt.dsdot_k.?)){ .kern = vt.dsdot_k.?, .conj = false }, dotRun));
    if (failed) |e| return e;
}

pub fn checkAll(vt: *const Vtable) !void {
    return checkMode(vt, .specials);
}

pub fn checkDenormals(vt: *const Vtable) !void {
    return checkMode(vt, .denormals);
}

pub const XDouble = types.XDouble;

test "portable level-1 kernels: special values, poisoned gaps, and negative and large strides" {
    const vt = check.portableVtable();
    try checkAll(vt);
    failed = null;
    if (XDouble != f64) {
        slots("q", XDouble, .specials, vt);
        slots("x", C(XDouble), .specials, vt);
    }
    if (failed) |e| return e;
}

test "every linked core's level-1 kernels: special values, poisoned gaps, and negative and large strides" {
    try check.eachCore(checkAll);
}

// ---- denormals --------------------------------------------------------------------

test "portable level-1 kernels: denormal results" {
    const vt = check.portableVtable();
    try checkDenormals(vt);
}

test "every linked core's level-1 kernels: denormal results" {
    try check.eachCore(checkDenormals);
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

test "portable level-1 reductions and element-wise kernels: the same bits whatever the address" {
    const vt = check.portableVtable();
    try checkShifted(vt);
}

test "every linked core's level-1 reductions and element-wise kernels: the same bits whatever the address" {
    try check.eachCore(checkShifted);
}

test "every linked core keeps the portable kernels' contract: the extremes' NaN and index, nrm2's class" {
    try check.eachCore(contract);
}
