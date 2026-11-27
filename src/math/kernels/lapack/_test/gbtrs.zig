//! Tests: `gbtrs` and `gbsv`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const ops = @import("../ops.zig");
const gbtrf_factor = @import("../gbtrf.zig");
const gbtrs_solve = @import("../gbtrs.zig");
const drivers = @import("../drivers.zig");

const testing = std.testing;
const Float = types.Float;
const Trans = types.Trans;

/// A random band `A` (`n x n`, `kl`, `ku`), full and in band storage
/// (`ldab = 2 kl + ku + 1`).
fn Band(comptime T: type) type {
    return struct {
        n: usize,
        kl: usize,
        ku: usize,
        full: []T,
        ab: []T,

        fn init(n: usize, kl: usize, ku: usize, seed: u64) !@This() {
            const al = testing.allocator;
            const ld = 2 * kl + ku + 1;
            const full = try al.alloc(T, n * n);
            @memset(full, ops.zero(T));
            const ab = try al.alloc(T, ld * n);
            @memset(ab, ops.zero(T));
            const vals = try al.alloc(T, n * n);
            defer al.free(vals);
            check.fill(T, vals, seed);
            for (0..n) |j| for (0..n) |i| {
                if (i + ku < j or i > j + kl) continue;
                full[i + j * n] = vals[i + j * n];
                ab[kl + ku + i - j + j * ld] = vals[i + j * n];
            };
            return .{ .n = n, .kl = kl, .ku = ku, .full = full, .ab = ab };
        }

        fn deinit(b: @This()) void {
            testing.allocator.free(b.full);
            testing.allocator.free(b.ab);
        }
    };
}

/// `‖op(A) X - B‖` over the entries.
fn residual(comptime T: type, trans: Trans, n: usize, nrhs: usize, a: []const T, x: []const T, b: []const T) Float(T) {
    var worst: Float(T) = 0;
    for (0..nrhs) |c| for (0..n) |i| {
        var s = ops.zero(T);
        for (0..n) |k| {
            var e = if (trans.transposed()) a[k + i * n] else a[i + k * n];
            if (trans.conjA()) e = ops.conj(T, e);
            s = ops.add(T, s, ops.mul(T, e, x[k + c * n]));
        }
        worst = @max(worst, ops.abs(T, ops.sub(T, s, b[i + c * n])));
    };
    return worst;
}

test "gbtrs: op(A) X = B from gbtrf's factor, every trans and band shape" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| for ([_][3]usize{ .{ 1, 0, 0 }, .{ 5, 1, 2 }, .{ 12, 3, 1 }, .{ 9, 2, 0 }, .{ 9, 0, 3 } }) |s| {
        const n = s[0];
        const kl = s[1];
        const ku = s[2];
        const band = try Band(T).init(n, kl, ku, n * 11 + kl * 3 + ku);
        defer band.deinit();
        const ld = 2 * kl + ku + 1;
        const ipiv = try al.alloc(usize, n);
        defer al.free(ipiv);
        try testing.expectEqual(null, try gbtrf_factor.gbtrf(T, n, n, kl, ku, band.ab, ld, ipiv));
        const nrhs = 3;
        const b = try al.alloc(T, n * nrhs);
        defer al.free(b);
        check.fill(T, b, n + 5);
        for ([_]Trans{ .n, .t, .r, .c }) |trans| {
            const x = try al.dupe(T, b);
            defer al.free(x);
            try gbtrs_solve.gbtrs(al, check.io, T, trans, n, kl, ku, nrhs, band.ab, ld, ipiv, x, n);
            var big: Float(T) = 0;
            for (band.full) |v| big = @max(big, ops.abs(T, v));
            var xbig: Float(T) = 0;
            for (x) |v| xbig = @max(xbig, ops.abs(T, v));
            const tol = 200 * @as(Float(T), @floatFromInt(n * (kl + ku + 1))) * std.math.floatEps(Float(T)) * @max(big, 1) * @max(xbig, 1);
            const res = residual(T, trans, n, nrhs, band.full, x, b);
            if (!(res <= tol)) {
                check.say("gbtrs {s} {t} n={d} kl={d} ku={d}: residual {e}, tol {e}\n", .{ @typeName(T), trans, n, kl, ku, res, tol });
                return error.TestUnexpectedResult;
            }
        }
    };
}

test "gbsv: A X = B; a singular band; bad arguments" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| {
        const n = 10;
        const kl = 2;
        const ku = 1;
        const band = try Band(T).init(n, kl, ku, 77);
        defer band.deinit();
        var ipiv: [n]usize = undefined;
        var b: [n * 2]T = undefined;
        check.fill(T, &b, 78);
        var x = b;
        try drivers.gbsv(al, check.io, T, n, kl, ku, 2, band.ab, 2 * kl + ku + 1, &ipiv, &x, n);
        try testing.expect(residual(T, .n, n, 2, band.full, &x, &b) <= 1e-12 * 100);
    }
    var z: [8]f64 = @splat(0);
    var piv: [2]usize = undefined;
    var b = [_]f64{ 1, 2 };
    try testing.expectError(error.Singular, drivers.gbsv(al, check.io, f64, 2, 1, 1, 1, &z, 4, &piv, &b, 2));
    try testing.expectEqualSlices(f64, &.{ 1, 2 }, &b);
    try testing.expectError(error.LeadingDimensionTooSmall, drivers.gbsv(al, check.io, f64, 2, 1, 1, 1, &z, 3, &piv, &b, 2));
    var badpiv = [_]usize{ 5, 0 };
    try testing.expectError(error.PivotOutOfRange, gbtrs_solve.gbtrs(al, check.io, f64, .n, 2, 1, 1, 1, &z, 4, &badpiv, &b, 2));
}
