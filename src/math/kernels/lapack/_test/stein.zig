//! Tests: `stein` on `stebz`'s eigenvalues, as Reference-LAPACK 3.12.0's
//! `TESTING/EIG/dchkst.f` checks `DSTEIN`: `dstt22.f`'s residual and
//! orthogonality ratios under `sep.in`'s threshold (60), random, split,
//! graded and Wilkinson matrices, real and complex `Z`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const ops = @import("../ops.zig");
const bisect = @import("../stebz.zig");
const inverse = @import("../stein.zig");
const stebz_tests = @import("stebz.zig");
const random = @import("../laruv.zig");

const testing = std.testing;
const Float = types.Float;
const Kind = stebz_tests.Kind;

/// `sep.in`'s `THRESH`.
const thresh = 60;

fn fail(comptime T: type, kind: Kind, n: usize, what: []const u8, value: anytype) error{TestUnexpectedResult} {
    check.say("stein {s} {t} n {d}: {s} {any}\n", .{ @typeName(T), kind, n, what, value });
    return error.TestUnexpectedResult;
}

fn case(comptime T: type, kind: Kind, n: usize) !void {
    const F = Float(T);
    const al = testing.allocator;
    const d = try al.alloc(F, n);
    defer al.free(d);
    const e = try al.alloc(F, n -| 1);
    defer al.free(e);
    stebz_tests.tridiagonal(F, kind, d, e);
    const w = try al.alloc(F, n);
    defer al.free(w);
    const iblock = try al.alloc(usize, n);
    defer al.free(iblock);
    const isplit = try al.alloc(usize, n);
    defer al.free(isplit);
    const found = try bisect.stebz(al, F, .all, .by_block, d, e, 2 * std.math.floatMin(F), w, iblock, isplit);
    const m = found.count;
    const z = try al.alloc(T, n * m);
    defer al.free(z);
    @memset(z, lcheck.poisoned(T));
    const ifail = try al.alloc(usize, m);
    defer al.free(ifail);
    const failed = try inverse.stein(al, check.io, T, d, e, w[0..m], iblock[0..m], isplit, z, n, ifail);
    if (failed != 0) return fail(T, kind, n, "unconverged", ifail[0..failed]);

    const ulp = std.math.floatEps(F);
    const nf: F = @floatFromInt(n);
    var tnorm: F = 0;
    for (0..n) |i| tnorm = @max(tnorm, @abs(d[i]) + (if (i > 0) @abs(e[i - 1]) else 0) + (if (i + 1 < n) @abs(e[i]) else 0));
    // `‖T z - w z‖₁ / (n ulp ‖T‖₁)`, and `‖I - Zᵀ Z‖₁ / (n ulp)`
    var r1: F = 0;
    var r2: F = 0;
    for (0..m) |j| {
        const zj = z[j * n ..][0..n];
        var s: F = 0;
        for (0..n) |i| {
            if (comptime types.isComplex(T)) if (zj[i].im != 0) return fail(T, kind, n, "imaginary part at", .{ i, j });
            var t = (d[i] - w[j]) * ops.re(T, zj[i]);
            if (i > 0) t += e[i - 1] * ops.re(T, zj[i - 1]);
            if (i + 1 < n) t += e[i] * ops.re(T, zj[i + 1]);
            s += @abs(t);
        }
        r1 = @max(r1, s / (nf * ulp * @max(tnorm, std.math.floatMin(F))));
        var o: F = 0;
        for (0..m) |k| {
            var dot: F = 0;
            for (0..n) |i| dot += ops.re(T, zj[i]) * ops.re(T, z[k * n + i]);
            o += @abs((if (j == k) @as(F, 1) else 0) - dot);
        }
        r2 = @max(r2, o / (nf * ulp));
    }
    if (!(r1 < thresh)) return fail(T, kind, n, "residual", r1);
    if (!(r2 < thresh)) return fail(T, kind, n, "orthogonality", r2);
}

test "stein: dstt22's ratios on stebz's eigenvalues, every type" {
    inline for (lcheck.test_types) |T| {
        for (std.meta.tags(Kind)) |kind| for ([_]usize{ 1, 2, 5, 23, 60 }) |n| try case(T, kind, n);
    }
}

test "stein: empty sizes, bad arguments" {
    const al = testing.allocator;
    const d = [_]f64{ 3, 1, 2 };
    const e = [_]f64{ 1, 1 };
    var z: [9]f64 = undefined;
    var ifail: [3]usize = undefined;
    const isplit = [_]usize{3};
    try testing.expectEqual(@as(usize, 0), try inverse.stein(al, check.io, f64, &d, &e, &.{}, &.{}, &isplit, z[0..0], 3, &ifail));
    const ee = testing.expectError;
    try ee(error.TooManyEigenvalues, inverse.stein(al, check.io, f64, d[0..1], e[0..0], &.{ 1, 2 }, &.{ 0, 0 }, &isplit, &z, 1, &ifail));
    try ee(error.EigenvaluesOutOfOrder, inverse.stein(al, check.io, f64, &d, &e, &.{ 2, 1 }, &.{ 0, 0 }, &isplit, &z, 3, &ifail));
    try ee(error.BlocksOutOfOrder, inverse.stein(al, check.io, f64, &d, &e, &.{ 1, 2 }, &.{ 1, 0 }, &.{ 1, 3 }, &z, 3, &ifail));
    try ee(error.SplitOutOfOrder, inverse.stein(al, check.io, f64, &d, &e, &.{1}, &.{0}, &.{4}, &z, 3, &ifail));
    try ee(error.LeadingDimensionTooSmall, inverse.stein(al, check.io, f64, &d, &e, &.{1}, &.{0}, &isplit, &z, 2, &ifail));
}

/// A row of `laruv.mm` (or a seed) as the 48-bit number its digits spell.
fn digits(row: anytype) u64 {
    var v: u64 = 0;
    for (row) |dgt| v = v * 4096 + dgt;
    return v;
}

test "laruv: MM's rows are its first row's powers mod 2^48, and the draws follow them" {
    const mask: u64 = (1 << 48) - 1;
    const a = digits(random.mm[0]);
    var p = a;
    for (random.mm, 0..) |row, i| {
        if (digits(row) != p) return fail(f64, .random, i, "MM row", row);
        p = (p *% a) & mask;
    }
    var seed: [4]u32 = @splat(1);
    const s0 = digits(seed);
    var x: [3]f64 = undefined;
    random.laruv(f64, &seed, &x);
    for (x, 0..) |v, i| try testing.expectEqual(@as(f64, @floatFromInt((s0 *% digits(random.mm[i])) & mask)) / (1 << 48), v);
    try testing.expectEqual((s0 *% digits(random.mm[2])) & mask, digits(seed));
}
