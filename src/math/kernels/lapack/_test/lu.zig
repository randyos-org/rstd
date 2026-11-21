//! Tests of `../lu.zig`'s `getrf`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const portable = @import("../../blas/portable.zig");
const Vtable = @import("../../blas/Vtable.zig");
const ops = @import("../ops.zig");
const lu = @import("../lu.zig");

const Allocator = std.mem.Allocator;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const getrf = lu.getrf;

const testing = std.testing;
const check = @import("../../blas/_test/check.zig");
const guard = @import("../../blas/_test/guard.zig");
const reference = @import("../../blas/_test/reference.zig");
const lcheck = @import("check.zig");
const fixtures = @import("fixtures.zig");

const W = lcheck.W;
const test_types = lcheck.test_types;
const isNanT = lcheck.isNanT;
const fromF64 = lcheck.fromF64;

pub const nan = std.math.nan(f64);
const inf = std.math.inf(f64);

/// `vt` with every precision's `q` at `2 * unroll_n + 1`, so the blocks are
/// capped at an odd width just past the base case's; static, as a vtable is
/// large for a rig's stack.
var shrunk: Vtable = undefined;

pub fn shrink(vt: *const Vtable) *const Vtable {
    shrunk = vt.*;
    inline for ([_]portable.Gemm{ .s, .d, .q, .c, .z, .x }) |g| {
        const un = @max(1, portable.paramPtr(&shrunk, g, .unroll_n).*);
        portable.paramPtr(&shrunk, g, .q).* = 2 * un + 1;
    }
    return &shrunk;
}

/// `m x n` shapes: square, tall, wide, a row and a column, across both
/// blockings' block edges.
const shapes = [_][2]usize{
    .{ 1, 1 },   .{ 2, 2 },   .{ 3, 3 },   .{ 7, 7 },   .{ 1, 9 },   .{ 9, 1 },   .{ 17, 17 },   .{ 33, 33 },
    .{ 40, 25 }, .{ 25, 40 }, .{ 50, 50 }, .{ 70, 70 }, .{ 90, 61 }, .{ 61, 90 }, .{ 130, 130 },
};

fn factorCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, shrink, test_types, factorOn);
}

fn factorOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
    _ = v;
    const al = testing.allocator;
    for (shapes) |s| {
        if (check.rig() and s[0] * s[1] > 70 * 70) continue;
        var p = try fixtures.LuProblem(T).init(al, s[0], s[1], s[0] * 31 + s[1]);
        defer p.deinit(al);
        var poison: guard.Poison = .{ .child = al };
        const zero = p.factorIt(poison.allocator()) catch |e| {
            check.say("getrf {s} {d}x{d}\n", .{ @typeName(T), s[0], s[1] });
            return e;
        };
        try testing.expectEqual(@as(?usize, null), zero);
        p.expectFactor(al) catch |e| {
            check.say("{s}, blocking {s}\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
            return e;
        };
    }
}

test "getrf: P L U = A, square, tall and wide, on every linked core" {
    try check.eachCore(factorCases);
}

/// `getrf` on `p` must report the zero pivot `want`, its factorization
/// complete.
fn expectSingular(comptime T: type, al: Allocator, p: *fixtures.LuProblem(T), want: usize, what: []const u8) !void {
    var poison: guard.Poison = .{ .child = al };
    const zero = try p.factorIt(poison.allocator());
    if (zero != want) {
        check.say("{s} {s} {d}x{d}: zero pivot {?d}, want {d}\n", .{ what, @typeName(T), p.m, p.n, zero, want });
        return error.TestUnexpectedResult;
    }
    p.expectFactor(al) catch |e| {
        check.say("{s} {s} {d}x{d}, singular at {d}\n", .{ what, @typeName(T), p.m, p.n, want });
        return e;
    };
}

/// Columns made zero: the first and last, and either side of block edges on
/// both blockings.
const zero_at = [_]usize{ 0, 1, 2, 5, 12, 13, 17, 25, 32, 33, 69 };

fn singularCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, shrink, test_types, singularOn);
}

fn singularOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
    _ = v;
    const al = testing.allocator;
    const n: usize = if (vi == 0) 40 else 70;
    for ([_]usize{ n, n + 9 }) |m| {
        var p = try fixtures.LuProblem(T).init(al, m, n, 5 + m);
        defer p.deinit(al);
        const aw = try al.dupe(W(T), p.aw);
        defer al.free(aw);
        for (zero_at) |k| {
            if (k >= n) continue;
            // a zero column `k`, then a second one after it
            @memcpy(p.aw, aw);
            for (0..m) |i| p.aw[i + k * m] = reference.wide(T, scalar(T, 0));
            if (k + 3 < n) for (0..m) |i| {
                p.aw[i + (k + 3) * m] = reference.wide(T, scalar(T, 0));
            };
            p.store();
            try expectSingular(T, al, &p, k, "zero column");
        }
        // rank one, `x yᵀ` with `x` of powers of two and `y(0) = 1`: every
        // step is exact, so `U(1, 1)` (0-based) is zero
        for (0..m) |i| for (0..n) |j| {
            const xi = std.math.ldexp(@as(f64, if (i % 3 == 1) -1 else 1), -@as(i32, @intCast((i * 7) % 5)));
            const yj: f64 = if (j == 0) 1 else @floatFromInt(@as(i32, @intCast(j % 7)) - 3);
            p.aw[i + j * m] = reference.wide(T, fromF64(T, xi * yj));
        };
        p.store();
        try expectSingular(T, al, &p, 1, "rank one");
    }
}

test "getrf: the first zero pivot, the factorization completed, on every linked core" {
    try check.eachCore(singularCases);
}

/// A NaN anywhere is a pivot's or reaches one, and is never zero: no error,
/// even where it reaches a zero column, as the reference's `INFO`.
fn nanCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, shrink, test_types, nanOn);
}

fn nanOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
    _ = v;
    _ = vi;
    const al = testing.allocator;
    const n = 40;
    var p = try fixtures.LuProblem(T).init(al, n, n, 11);
    defer p.deinit(al);
    for ([_][2]usize{ .{ 0, 0 }, .{ 7, 3 }, .{ 39, 20 }, .{ 0, 39 } }) |at| {
        p.store();
        p.a[at[0] + at[1] * p.lda] = fromF64(T, nan);
        // a zero column past the NaN's
        const z = @min(at[1] + 2, n - 1);
        if (z != at[1]) for (0..n) |i| {
            p.a[i + z * p.lda] = scalar(T, 0);
        };
        var poison: guard.Poison = .{ .child = al };
        if (try p.factorIt(poison.allocator())) |k| {
            check.say("getrf {s}: NaN at ({d},{d}): zero pivot {d}\n", .{ @typeName(T), at[0], at[1], k });
            return error.TestUnexpectedResult;
        }
        // the last pivot is NaN
        try testing.expect(isNanT(T, p.a[(n - 1) * (p.lda + 1)]));
    }
}

test "getrf: NaN is never a zero pivot, on every linked core" {
    try check.eachCore(nanCases);
}

/// An untouched pivot, which no row index can be.
pub const unset = std.math.maxInt(usize);

/// The `m x 1` column `col` factored: its pivot (0-based) and `L`.
fn column(comptime T: type, comptime m: usize, col: [m]T) !struct { usize, [m]T } {
    var a = col;
    var ipiv: [1]usize = .{unset};
    _ = try getrf(testing.failing_allocator, check.io, T, m, 1, &a, m, &ipiv);
    return .{ ipiv[0], a };
}

test "getrf: the first largest |re| + |im| pivots, the first NaN before it" {
    inline for (test_types) |T| {
        const c = struct {
            fn f(re: f64, im: f64) T {
                return if (comptime isComplex(T)) .init(@floatCast(re), @floatCast(im)) else @floatCast(re);
            }
        }.f;
        // ties keep the first; no swap for the diagonal's
        try testing.expectEqual(@as(usize, 0), (try column(T, 3, .{ c(1, 0), c(-1, 0), c(1, 0) }))[0]);
        try testing.expectEqual(@as(usize, 1), (try column(T, 4, .{ c(-1, 0), c(2, 0), c(-2, 0), c(2, 0) }))[0]);
        if (comptime isComplex(T)) {
            // 7, 5, 7, 7: the first; and 3 - 4i over 5, the same modulus
            try testing.expectEqual(@as(usize, 0), (try column(T, 4, .{ c(3, 4), c(5, 0), c(-4, 3), c(0, -7) }))[0]);
            try testing.expectEqual(@as(usize, 2), (try column(T, 3, .{ c(1, 0), c(5, 0), c(3, -4) }))[0]);
        }
        // the first NaN wins, and fills L
        const r = try column(T, 4, .{ c(1, 0), c(nan, 0), c(5, 0), c(nan, 0) });
        try testing.expectEqual(@as(usize, 1), r[0]);
        for (r[1][1..]) |x| try testing.expect(if (comptime isComplex(T)) std.math.isNan(x.re) else std.math.isNan(x));
        // an infinite pivot: finite elements over it are zero, infinite ones NaN
        const s = try column(T, 3, .{ c(1, 0), c(inf, 0), c(-inf, 0) });
        try testing.expectEqual(@as(usize, 1), s[0]);
        try testing.expectEqual(c(0, 0), s[1][1]);
        try testing.expect(if (comptime isComplex(T)) std.math.isNan(s[1][2].re) else std.math.isNan(s[1][2]));
    }
}

test "getrf: a denormal pivot is swapped in and divided by" {
    inline for (test_types) |T| {
        const F = Float(T);
        const tm = std.math.floatTrueMin(F);
        // [[4tm, 1], [8tm, 3]]: rows swapped, L(2, 1) = 1/2, U(2, 2) = -1/2
        var a: [4]T = .{ scalar(T, 4 * tm), scalar(T, 8 * tm), scalar(T, 1), scalar(T, 3) };
        var ipiv: [2]usize = @splat(unset);
        try testing.expectEqual(@as(?usize, null), try getrf(testing.failing_allocator, check.io, T, 2, 2, &a, 2, &ipiv));
        try testing.expectEqual([2]usize{ 1, 1 }, ipiv);
        try testing.expectEqual([4]T{ scalar(T, 8 * tm), scalar(T, 0.5), scalar(T, 3), scalar(T, -0.5) }, a);
        if (comptime isComplex(T)) {
            // 4tm over 8tm i: -i/2
            const r = try column(T, 2, .{ .init(4 * tm, 0), .init(0, 8 * tm) });
            try testing.expectEqual(@as(usize, 1), r[0]);
            try testing.expectEqual(T.init(0, -0.5), r[1][1]);
        }
    }
}

test "getrf: an exactly singular U, its first zero pivot returned" {
    inline for (.{ f64, types.Complex64 }) |T| {
        const al = testing.failing_allocator;
        var ipiv: [3]usize = @splat(unset);
        // a zero matrix's first pivot is zero, a zero row's `U(0, 0)` too
        var z: [9]T = @splat(scalar(T, 0));
        try testing.expectEqual(@as(?usize, 0), try getrf(al, check.io, T, 3, 3, &z, 3, &ipiv));
        try testing.expectEqual([3]usize{ 0, 1, 2 }, ipiv);
        var row: [3]T = .{ scalar(T, 0), scalar(T, 2), scalar(T, 3) };
        try testing.expectEqual(@as(?usize, 0), try getrf(al, check.io, T, 1, 3, &row, 1, &ipiv));
        // [[1, 2], [2, 4]]: rows swapped, then U(1, 1) = 2 - 1/2 4 = 0
        var a: [4]T = .{ scalar(T, 1), scalar(T, 2), scalar(T, 2), scalar(T, 4) };
        try testing.expectEqual(@as(?usize, 1), try getrf(al, check.io, T, 2, 2, &a, 2, &ipiv));
        try testing.expectEqual([4]T{ scalar(T, 2), scalar(T, 0.5), scalar(T, 4), scalar(T, 0) }, a);
        a = .{ scalar(T, 1), scalar(T, 2), scalar(T, 2), scalar(T, 5) };
        try testing.expectEqual(@as(?usize, null), try getrf(al, check.io, T, 2, 2, &a, 2, &ipiv));
        try testing.expectEqual([4]T{ scalar(T, 2), scalar(T, 0.5), scalar(T, 5), scalar(T, -0.5) }, a);
    }
}

test "getrf: a slice too short or leading dimension too small, nothing touched when empty" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    var ipiv: [2]usize = @splat(unset);
    try testing.expectError(error.LeadingDimensionTooSmall, getrf(fa, check.io, f64, 2, 2, &a, 1, &ipiv));
    try testing.expectError(error.LeadingDimensionTooSmall, getrf(fa, check.io, f64, 0, 2, &a, 0, &ipiv));
    try testing.expectError(error.SliceTooShort, getrf(fa, check.io, f64, 2, 3, &a, 2, &ipiv));
    try testing.expectError(error.SliceTooShort, getrf(fa, check.io, f64, 2, 2, &a, std.math.maxInt(usize) / 2, &ipiv));
    try testing.expectError(error.SliceTooShort, getrf(fa, check.io, f64, 2, 2, &a, 2, ipiv[0..1]));
    try testing.expectEqual(@as(?usize, null), try getrf(fa, check.io, f64, 0, 2, &a, 1, &ipiv));
    try testing.expectEqual(@as(?usize, null), try getrf(fa, check.io, f64, 2, 0, &a, 2, &ipiv));
    for (a) |v| try testing.expect(std.math.isNan(v));
    try testing.expectEqual([2]usize{ unset, unset }, ipiv);
}
