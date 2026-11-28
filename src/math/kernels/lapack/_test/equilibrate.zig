//! Tests: `geequ` (Reference-LAPACK 3.12.0 `TESTING/LIN/dchkeq.f`'s GEEQU
//! cases) and `laqge`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const lcheck = @import("check.zig");
const equilibrate = @import("../equilibrate.zig");

const testing = std.testing;
const Float = types.Float;
const isComplex = types.isComplex;

const nsz = 5;

/// `|got - want| / want`.
fn rel(comptime F: type, got: F, want: F) F {
    return @abs((got - want) / want);
}

/// `10^e` by products, and `1 / 10^-e` for negative `e`, as `dchkeq.f`'s
/// `POW` and `RPOW`.
fn pow10(comptime F: type, e: i32) F {
    var p: F = 1;
    for (0..@abs(e)) |_| p *= 10;
    return if (e < 0) 1 / p else p;
}

test "lapack geequ: dchkeq's powers of ten, every m x n up to 5 x 5" {
    inline for (lcheck.test_types) |T| {
        const F = Float(T);
        var worst: F = 0;
        for (0..nsz + 1) |n| for (0..nsz + 1) |m| {
            var a: [nsz * nsz]T = undefined;
            // A(i, j) = (-10)^(i + j), 1-based, inside the m x n corner
            for (0..nsz) |j| for (0..nsz) |i| {
                const e: i32 = @intCast(i + j + 2);
                const v: f64 = if (i < m and j < n) std.math.pow(f64, -10, @floatFromInt(e)) else 0;
                a[i + j * nsz] = lcheck.fromF64(T, v);
            };
            var r: [nsz]F = undefined;
            var c: [nsz]F = undefined;
            const s = try equilibrate.geequ(T, m, n, &a, nsz, &r, &c);
            if (m == 0 or n == 0) {
                try testing.expectEqual(@as(F, 1), s.rowcnd);
                try testing.expectEqual(@as(F, 1), s.colcnd);
                try testing.expectEqual(@as(F, 0), s.amax);
                continue;
            }
            const mi: i32 = @intCast(m);
            const ni: i32 = @intCast(n);
            worst = @max(worst, rel(F, s.rowcnd, pow10(F, 1 - mi)), rel(F, s.colcnd, pow10(F, 1 - ni)), rel(F, s.amax, pow10(F, mi + ni)));
            for (0..m) |i| worst = @max(worst, rel(F, r[i], pow10(F, -(@as(i32, @intCast(i)) + 1 + ni))));
            for (0..n) |j| worst = @max(worst, rel(F, c[j], pow10(F, ni - 1 - @as(i32, @intCast(j)))));
        };
        try testing.expect(worst / std.math.floatEps(F) < lcheck.thresh);
    }
}

test "lapack geequ: a zero row or column; bad arguments" {
    inline for (lcheck.test_types) |T| {
        const F = Float(T);
        var a: [nsz * nsz]T = @splat(lcheck.fromF64(T, 1));
        var r: [nsz]F = undefined;
        var c: [nsz]F = undefined;
        for (0..nsz) |j| a[3 + j * nsz] = lcheck.fromF64(T, 0);
        try testing.expectError(error.ZeroRow, equilibrate.geequ(T, nsz, nsz, &a, nsz, &r, &c));
        for (0..nsz) |j| a[3 + j * nsz] = lcheck.fromF64(T, 1);
        for (0..nsz) |i| a[i + 3 * nsz] = lcheck.fromF64(T, 0);
        try testing.expectError(error.ZeroColumn, equilibrate.geequ(T, nsz, nsz, &a, nsz, &r, &c));
        try testing.expectError(error.LeadingDimensionTooSmall, equilibrate.geequ(T, nsz, nsz, &a, nsz - 1, &r, &c));
        try testing.expectError(error.SliceTooShort, equilibrate.geequ(T, nsz, nsz, &a, nsz, r[0..4], &c));
        try testing.expectError(error.SliceTooShort, equilibrate.geequ(T, nsz, nsz, a[0..20], nsz, &r, &c));
    }
}

/// `laqge` on a 2 x 2 of ones with `r = {2, 4}`, `c = {8, 16}`.
fn laqgeCase(comptime T: type, s: equilibrate.Scaling(Float(T)), want: equilibrate.Equilibration) !void {
    const F = Float(T);
    var a: [4]T = @splat(lcheck.fromF64(T, if (comptime isComplex(T)) 1 else -1));
    if (comptime isComplex(T)) {
        for (&a) |*x| x.im = -3;
    }
    const a0 = a;
    const r = [_]F{ 2, 4 };
    const c = [_]F{ 8, 16 };
    try testing.expectEqual(want, try equilibrate.laqge(T, 2, 2, &a, 2, &r, &c, s));
    for (0..2) |j| for (0..2) |i| {
        const k: F = (if (want == .columns or want == .both) c[j] else 1) * (if (want == .rows or want == .both) r[i] else 1);
        const x = a0[i + 2 * j];
        const w: T = if (comptime isComplex(T)) .init(k * x.re, k * x.im) else k * x;
        try testing.expectEqual(w, a[i + 2 * j]);
    };
}

test "lapack laqge: rows, columns, both or none by the 0.1 ratios and amax's range" {
    inline for (lcheck.test_types) |T| {
        const F = Float(T);
        const tiny = std.math.floatMin(F) / std.math.floatEps(F);
        try laqgeCase(T, .{ .rowcnd = 0.5, .colcnd = 0.5, .amax = 1 }, .none);
        try laqgeCase(T, .{ .rowcnd = 0.1, .colcnd = 0.1, .amax = 1 }, .none);
        try laqgeCase(T, .{ .rowcnd = 0.5, .colcnd = 0.01, .amax = 1 }, .columns);
        try laqgeCase(T, .{ .rowcnd = 0.01, .colcnd = 0.5, .amax = 1 }, .rows);
        try laqgeCase(T, .{ .rowcnd = 0.01, .colcnd = 0.01, .amax = 1 }, .both);
        // amax outside [small, 1 / small] scales rows even when rowcnd is fine
        try laqgeCase(T, .{ .rowcnd = 0.5, .colcnd = 0.5, .amax = tiny / 2 }, .rows);
        try laqgeCase(T, .{ .rowcnd = 0.5, .colcnd = 0.5, .amax = 2 / tiny }, .rows);
        var a: [4]T = undefined;
        const r = [_]F{ 1, 1 };
        try testing.expectEqual(equilibrate.Equilibration.none, try equilibrate.laqge(T, 0, 2, &a, 1, &r, &r, .{ .rowcnd = 0, .colcnd = 0, .amax = 0 }));
        try testing.expectError(error.SliceTooShort, equilibrate.laqge(T, 2, 2, &a, 2, r[0..1], &r, .{ .rowcnd = 0, .colcnd = 0, .amax = 0 }));
    }
}
