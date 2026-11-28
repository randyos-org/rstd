//! Tests: `syevx` by Reference-LAPACK 3.12.0's `TESTING/EIG/ddrvst.f`
//! `DSYEVX` cases: `dsyt22.f`'s `‖Zᴴ A Z - W‖₁ / (n ‖A‖₁ ulp)` and
//! `‖I - Zᴴ Z‖₁ / (n ulp)`, the eigenvalues against `syevd`'s within
//! `dstech.f`'s tolerance, all under `sep.in`'s threshold (60): every range,
//! `jobz`, `uplo` and `abstol`, random, tiny and huge matrices (the latter
//! two scaled by `ISCALE`).

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const ops = @import("../ops.zig");
const norm = @import("../norm.zig");
const evd = @import("../evd.zig");
const interface = @import("../interface.zig");
const selected = @import("../syevx.zig");

const testing = std.testing;
const Float = types.Float;
const Mat = ops.Mat;
const Uplo = types.Uplo;
const SpectrumRange = interface.SpectrumRange;

/// `sep.in`'s `THRESH`.
const thresh = 60;

/// `A`'s kinds: random Hermitian; the same scaled below `ISCALE`'s `rmin`
/// and above its `rmax`; two random diagonal blocks, so `T` splits and
/// `stebz` returns the eigenvalues by block, out of order.
const Kind = enum { random, tiny, huge, blocks };

fn fail(comptime T: type, kind: Kind, uplo: Uplo, n: usize, what: []const u8, value: anytype) error{TestUnexpectedResult} {
    check.say("syevx {s} {t} uplo {t} n {d}: {s} {any}\n", .{ @typeName(T), kind, uplo, n, what, value });
    return error.TestUnexpectedResult;
}

fn matrix(comptime T: type, kind: Kind, a: Mat(T)) void {
    const F = Float(T);
    lcheck.hermitian(T, a, a.rows * 41 + 3);
    const safmin = std.math.floatMin(F);
    const smlnum = safmin / std.math.floatEps(F);
    const s: F = switch (kind) {
        .random => 1,
        .tiny => @sqrt(smlnum) / 64,
        .huge => @min(@sqrt(1 / smlnum), 1 / @sqrt(@sqrt(safmin))) * 64,
        .blocks => 1,
    };
    const h = a.rows / 2;
    if (kind == .blocks) for (0..a.cols) |j| for (0..a.rows) |i| {
        if ((i < h) != (j < h)) a.at(i, j).* = ops.zero(T);
    };
    for (0..a.cols) |j| for (0..a.rows) |i| {
        a.at(i, j).* = ops.mulReal(T, a.get(i, j), s);
    };
}

fn case(comptime T: type, kind: Kind, uplo: Uplo, n: usize) !void {
    const F = Float(T);
    const al = testing.allocator;
    const asav = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, asav);
    matrix(T, kind, asav);
    const a = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, a);
    // every eigenvalue, by `syevd`
    const want = try al.alloc(F, n);
    defer al.free(want);
    a.copyFrom(asav);
    try evd.syevd(al, check.io, T, .none, uplo, n, a.slice(), a.ld, want);
    const ulp = std.math.floatEps(F);
    var mx: F = 0;
    for (want) |v| mx = @max(mx, @abs(v));
    const tol = thresh * ulp * mx + std.math.floatMin(F) / ulp;
    const anorm = @max(try norm.lange(al, check.io, T, .one, n, n, asav.slice(), asav.ld), std.math.floatMin(F));
    const nf: F = @floatFromInt(@max(n, 1));

    const w = try al.alloc(F, n);
    defer al.free(w);
    const z = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, z);
    const ifail = try al.alloc(usize, n);
    defer al.free(ifail);
    const first = n / 4;
    const last = n -| (1 + n / 3);
    const runs = [_]struct { range: SpectrumRange(F), from: usize, count: usize }{
        .{ .range = .all, .from = 0, .count = n },
        .{ .range = .{ .indices = .{ .first = first, .last = last } }, .from = first, .count = last - first + 1 },
        .{ .range = .{ .values = .{ .lo = if (first > 0) 0.5 * (want[first - 1] + want[first]) else want[0] - 1, .hi = if (last + 1 < n) 0.5 * (want[last] + want[last + 1]) else want[n - 1] + 1 } }, .from = first, .count = last - first + 1 },
    };
    for (runs, 0..) |run, ri| {
        if (n < 3 and ri > 0) continue;
        for ([_]F{ 0, 2 * std.math.floatMin(F) }) |abstol| for ([_]interface.VectorJob{ .none, .compute }) |jobz| {
            a.copyFrom(asav);
            const got = try selected.syevx(al, check.io, T, jobz, run.range, uplo, n, a.slice(), a.ld, abstol, w, z.slice(), z.ld, ifail);
            if (got.count != run.count or got.unconverged != 0 or got.values_unconverged or got.incomplete) return fail(T, kind, uplo, n, "found", got);
            for (w[0..got.count], want[run.from..][0..got.count]) |x, y| if (!(@abs(x - y) <= tol)) return fail(T, kind, uplo, n, "eigenvalue", .{ x, y });
            if (jobz == .none) continue;
            const m = got.count;
            const zm = z.sub(0, 0, n, m);
            // `Zᴴ A Z - W`, `I - Zᴴ Z`
            const az = try ops.alloc(T, al, n, m);
            defer ops.free(T, al, az);
            lcheck.mul(T, false, asav, false, zm, az);
            const r = try ops.alloc(T, al, m, m);
            defer ops.free(T, al, r);
            lcheck.mul(T, true, zm, false, az, r);
            for (0..m) |i| r.at(i, i).* = ops.sub(T, r.get(i, i), ops.scalar(T, w[i]));
            const r1 = try norm.lange(al, check.io, T, .one, m, m, r.slice(), r.ld) / (nf * anorm * ulp);
            if (!(r1 < thresh)) return fail(T, kind, uplo, n, "Zᴴ A Z - W", r1);
            lcheck.mul(T, true, zm, false, zm, r);
            for (0..m) |i| r.at(i, i).* = ops.sub(T, r.get(i, i), ops.one(T));
            const r2 = try norm.lange(al, check.io, T, .one, m, m, r.slice(), r.ld) / (nf * ulp);
            if (!(r2 < thresh)) return fail(T, kind, uplo, n, "I - Zᴴ Z", r2);
        };
    }
}

test "syevx: ddrvst's ratios, every range, jobz, uplo and abstol, scaled matrices, every type" {
    inline for (lcheck.test_types) |T| for (std.meta.tags(Kind)) |kind| for ([_]Uplo{ .u, .l }) |uplo| for ([_]usize{ 1, 2, 5, 23, 40 }) |n| {
        try case(T, kind, uplo, n);
    };
}

test "syevx: one element, a range missing it, empty, bad arguments" {
    const al = testing.allocator;
    var a = [_]f64{ 3, 1, 1, 2 };
    var w: [2]f64 = undefined;
    var z: [4]f64 = undefined;
    var ifail: [2]usize = undefined;
    const one = try selected.syevx(al, check.io, f64, .compute, .all, .u, 1, a[0..1], 1, 0, &w, &z, 1, &ifail);
    try testing.expectEqual(@as(usize, 1), one.count);
    try testing.expectEqual(@as(f64, 3), w[0]);
    try testing.expectEqual(@as(f64, 1), z[0]);
    const none = try selected.syevx(al, check.io, f64, .none, .{ .values = .{ .lo = 3, .hi = 4 } }, .u, 1, a[0..1], 1, 0, &w, &z, 1, &ifail);
    try testing.expectEqual(@as(usize, 0), none.count);
    const empty = try selected.syevx(al, check.io, f64, .compute, .all, .l, 0, a[0..0], 1, 0, w[0..0], z[0..0], 1, ifail[0..0]);
    try testing.expectEqual(@as(usize, 0), empty.count);
    const ee = testing.expectError;
    try ee(error.UnsupportedUplo, selected.syevx(al, check.io, f64, .none, .all, .v, 2, &a, 2, 0, &w, &z, 2, &ifail));
    try ee(error.EmptyInterval, selected.syevx(al, check.io, f64, .none, .{ .values = .{ .lo = 1, .hi = 1 } }, .u, 2, &a, 2, 0, &w, &z, 2, &ifail));
    try ee(error.IndexOutOfRange, selected.syevx(al, check.io, f64, .none, .{ .indices = .{ .first = 1, .last = 2 } }, .u, 2, &a, 2, 0, &w, &z, 2, &ifail));
    try ee(error.SliceTooShort, selected.syevx(al, check.io, f64, .compute, .all, .u, 2, &a, 2, 0, &w, z[0..3], 2, &ifail));
    try ee(error.LeadingDimensionTooSmall, selected.syevx(al, check.io, f64, .compute, .all, .u, 2, &a, 2, 0, &w, &z, 1, &ifail));
}
