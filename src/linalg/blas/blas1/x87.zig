//! The sums of OpenBLAS's x87 `q`/`x` kernels (`kernel/x86_64/dot.S`,
//! `zdot.S`, `asum.S`, `zasum.S`, `nrm2.S`, `znrm2.S`), which the portable
//! reductions take for floats wider than `f64`. At `f80` they give the
//! assembly's bits: the same accumulators, filled and folded in its order,
//! with nothing fused. No overflow guard in `nrm2`, as in the assembly.
//! `kernel/x86/`'s are in the same order. The gemv sums are in
//! `blas2/gemv/x87.zig`.

const std = @import("std");
const types = @import("../types.zig");
const util = @import("../util.zig");
const check = @import("check.zig");
const zigstd = @import("../../../zigstd.zig");

const BlasLong = types.BlasLong;

/// `dot.S`: four accumulators by `i % 4`, the tail into the first.
pub fn dot(comptime F: type, n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) F {
    var a: [4]F = @splat(0);
    var i: BlasLong = 0;
    const n4 = n & -4;
    while (i < n4) : (i += 4) {
        inline for (0..4) |c| a[c] += util.at(F, x, (i + c) * inc_x).* * util.at(F, y, (i + c) * inc_y).*;
    }
    while (i < n) : (i += 1) a[0] += util.at(F, x, i * inc_x).* * util.at(F, y, i * inc_y).*;
    return (a[3] + a[1]) + (a[2] + a[0]);
}

/// `zdot.S`: each of an element's four products in its own sum, combined at
/// the end; `conj` conjugates `x`. Real and imaginary parts.
pub fn zdot(comptime F: type, comptime conj: bool, n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) [2]F {
    var a: [4]F = @splat(0);
    var i: BlasLong = 0;
    while (i < n) : (i += 1) {
        const xr = util.at(F, x, i * inc_x * 2).*;
        const xi = util.at(F, x, i * inc_x * 2 + 1).*;
        const yr = util.at(F, y, i * inc_y * 2).*;
        const yi = util.at(F, y, i * inc_y * 2 + 1).*;
        a[0] += xr * yr;
        a[1] += xr * yi;
        a[2] += xi * yr;
        a[3] += xi * yi;
    }
    return if (conj) .{ a[0] + a[3], a[1] - a[2] } else .{ a[0] - a[3], a[1] + a[2] };
}

/// What `sum` adds up per float.
pub const Term = enum { abs, square };

/// `asum.S`, `nrm2.S` and their complex forms, before `nrm2`'s square root:
/// eight floats a trip into accumulators by float `% 4`, then each remaining
/// element's float `c` (of `w`) into accumulator `c`.
pub fn sum(comptime F: type, comptime w: comptime_int, comptime term: Term, n: BlasLong, x: [*]F, inc_x: BlasLong) F {
    const per = 8 / w;
    var a: [4]F = @splat(0);
    var i: BlasLong = 0;
    const main = n & -per;
    while (i < main) : (i += per) {
        inline for (0..8) |f| a[f % 4] += of(F, term, util.at(F, x, (i + f / w) * inc_x * w + f % w).*);
    }
    while (i < n) : (i += 1) {
        inline for (0..w) |c| a[c] += of(F, term, util.at(F, x, i * inc_x * w + c).*);
    }
    return a[3] + ((a[2] + a[0]) + a[1]);
}

inline fn of(comptime F: type, comptime term: Term, v: F) F {
    return switch (term) {
        .abs => @abs(v),
        .square => v * v,
    };
}

// ---- tests -----------------------------------------------------------------

/// The same sums written out plainly: accumulator by position, then folded.
const Plain = struct {
    fn dot(comptime F: type, n: usize, x: []const F, ix: usize, y: []const F, iy: usize) F {
        var a: [4]F = @splat(0);
        for (0..n) |i| a[if (i < n / 4 * 4) i % 4 else 0] += x[i * ix] * y[i * iy];
        return (a[3] + a[1]) + (a[2] + a[0]);
    }

    fn zdot(comptime F: type, conj: bool, n: usize, x: []const F, ix: usize, y: []const F, iy: usize) [2]F {
        var a: [4]F = @splat(0);
        for (0..n) |i| for (0..2) |p| for (0..2) |q| {
            a[2 * p + q] += x[2 * i * ix + p] * y[2 * i * iy + q];
        };
        return if (conj) .{ a[0] + a[3], a[1] - a[2] } else .{ a[0] - a[3], a[1] + a[2] };
    }

    fn sum(comptime F: type, w: usize, term: Term, n: usize, x: []const F, ix: usize) F {
        var a: [4]F = @splat(0);
        const main = n / (8 / w) * (8 / w);
        for (0..n) |i| for (0..w) |c| {
            const v = x[i * ix * w + c];
            const t = if (term == .abs) @abs(v) else v * v;
            a[if (i < main) (i * w + c) % 4 else c] += t;
        };
        return a[3] + ((a[2] + a[0]) + a[1]);
    }
};

fn expectBits(comptime F: type, want: F, got: F, what: []const u8, n: usize, inc: usize) !void {
    const U = @Int(.unsigned, @bitSizeOf(F));
    if (@as(U, @bitCast(want)) == @as(U, @bitCast(got))) return;
    check.say("{s} {s} n={d} inc={d}: got {e}, want {e}\n", .{ what, @typeName(F), n, inc, got, want });
    return error.TestExpectedEqual;
}

test "the x87 sums fold in the assembly's order, bit for bit" {
    inline for (.{ f80, f128 }) |F| {
        var x: [check.capacity * 2]F = undefined;
        var y: [check.capacity * 2]F = undefined;
        for (check.sizes) |n| for ([_]usize{ 1, 3 }) |inc| {
            check.fill(F, &x, n + 11);
            check.fill(F, &y, n + 12);
            const len: BlasLong = @intCast(n);
            const li: BlasLong = @intCast(inc);
            try expectBits(F, Plain.dot(F, n, &x, inc, &y, 1), dot(F, len, &x, li, &y, 1), "dot", n, inc);
            inline for (.{ false, true }) |conj| {
                const want = Plain.zdot(F, conj, n, &x, inc, &y, 2);
                const got = zdot(F, conj, len, &x, li, &y, 2);
                try expectBits(F, want[0], got[0], "zdot re", n, inc);
                try expectBits(F, want[1], got[1], "zdot im", n, inc);
            }
            inline for (.{ 1, 2 }) |w| inline for (.{ Term.abs, Term.square }) |term| {
                try expectBits(F, Plain.sum(F, w, term, n, &x, inc), sum(F, w, term, len, &x, li), @tagName(term), n, inc);
            };
        };
    }
}

test "the x87 sums at f80 carry the full significand" {
    // 1 + 2^-63 squared is 1 + 2^-62 in f80 and 1 in f64
    var x = [_]f80{ 1 + std.math.ldexp(@as(f80, 1), -63), 0, 0, 0, 0 };
    try zigstd.testExpectEqual(1 + std.math.ldexp(@as(f80, 1), -62), dot(f80, 5, &x, 1, &x, 1));
    try zigstd.testExpectEqual(1 + std.math.ldexp(@as(f80, 1), -62), sum(f80, 1, .square, 5, &x, 1));
}

test "the x87 sums give the values traced by hand from the .S files" {
    // 2^64 + 1 rounds to 2^64 in f80, so each value below depends on which
    // accumulator a term lands in and on the fold order. Traced once from
    // kernel/x86_64/*.S: dot `faddp %st, %st(k)` per i % 4, tail into st(0),
    // folded (a3 + a1) + (a2 + a0); zdot a0..a3 = xr*yr, xr*yi, xi*yr, xi*yi,
    // re = a0 - a3 (`fsubp` is fsubr in AT&T), im = a1 + a2, conj a0 + a3 and
    // a1 - a2; asum/nrm2 eight floats a trip by float % 4, tail floats into
    // a0 (zasum/znrm2: re to a0, im to a1), folded a3 + ((a2 + a0) + a1).
    const F = f80;
    const t = std.math.ldexp(@as(F, 1), 64);
    const h = std.math.ldexp(@as(F, 1), 32);
    var ones: [8]F = @splat(1);
    // dot: a0 = 1 + 1 + 1 (tail), a1 = t, a2 = 1, a3 = -t: (0) + (1 + 3)
    var dx = [_]F{ 1, t, 1, -t, 1, 1 };
    try zigstd.testExpectEqual(@as(F, 4), dot(F, 6, &dx, 1, &ones, 1));
    // zdot of (t, t), (1, -1) with (1, 1), (1, 1): a0 = t, a1 = t, a2 = t - 1, a3 = t - 1
    var zx = [_]F{ t, t, 1, -1 };
    const u = zdot(F, false, 2, &zx, 1, &ones, 1);
    try zigstd.testExpectEqual(@as(F, 1), u[0]);
    try zigstd.testExpectEqual(2 * t, u[1]); // 2t - 1 rounds to even
    const c = zdot(F, true, 2, &zx, 1, &ones, 1);
    try zigstd.testExpectEqual(2 * t, c[0]);
    try zigstd.testExpectEqual(@as(F, 1), c[1]);
    // asum: a0 = 1, a1 = 1, a2 = t: (t + 1) + 1 = t; a tail 1 into a0 gives
    // (t + 2) + 1 = t + 4 (ties to even), into a1 it would give t + 2
    var ax = [_]F{ 1, 1, t, 0, 0, 0, 0, 0, 1 };
    try zigstd.testExpectEqual(t, sum(F, 1, .abs, 8, &ax, 1));
    try zigstd.testExpectEqual(t + 4, sum(F, 1, .abs, 9, &ax, 1));
    // zasum: four complex a trip; the tail's imaginary part goes to a1
    var zax = [_]F{ 1, 1, t, 0, 0, 0, 0, 0, 0, 1 };
    try zigstd.testExpectEqual(t, sum(F, 2, .abs, 4, &zax, 1));
    try zigstd.testExpectEqual(t + 2, sum(F, 2, .abs, 5, &zax, 1));
    // nrm2 and znrm2 before the root: the same with squares
    var nx = [_]F{ 1, 1, h, 0, 0, 0, 0, 0, 1 };
    try zigstd.testExpectEqual(t, sum(F, 1, .square, 8, &nx, 1));
    try zigstd.testExpectEqual(t + 4, sum(F, 1, .square, 9, &nx, 1));
    var znx = [_]F{ 1, 1, h, 0, 0, 0, 0, 0, 0, 1 };
    try zigstd.testExpectEqual(t + 2, sum(F, 2, .square, 5, &znx, 1));
}
