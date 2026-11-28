//! Tests: `ggevReal`.

const std = @import("std");
const check = @import("../../blas/_test/check.zig");
const ggev = @import("../ggev_real.zig");
const complex_test = @import("ggev_complex.zig");

const testing = std.testing;

fn ggevCase(comptime T: type, n: usize, a0: []const T, b0: []const T) !void {
    const al = testing.allocator;
    const Z = std.math.Complex(T);
    const a = try al.dupe(T, a0);
    defer al.free(a);
    const b = try al.dupe(T, b0);
    defer al.free(b);
    const w = try al.alloc(T, 3 * n);
    defer al.free(w);
    const vl = try al.alloc(T, n * n);
    defer al.free(vl);
    const vr = try al.alloc(T, n * n);
    defer al.free(vr);
    try ggev.ggevReal(al, check.io, T, .compute, .compute, n, a, n, b, n, w[0..n], w[n .. 2 * n], w[2 * n ..], vl, n, vr, n);
    // complex copies of the pencil, for the residuals
    const za = try al.alloc(Z, 2 * n * n);
    defer al.free(za);
    for (a0, b0, 0..) |x, y, k| {
        za[k] = Z.init(x, 0);
        za[n * n + k] = Z.init(y, 0);
    }
    const zv = try al.alloc(Z, n);
    defer al.free(zv);
    const tol = 300 * @as(T, @floatFromInt(n)) * std.math.floatEps(T);
    var j: usize = 0;
    while (j < n) {
        const pair = w[n + j] != 0;
        // a pair's eigenvalues conjugate (`alpha` and `beta` each scale apart)
        if (pair) {
            try testing.expect(j + 1 < n and w[n + j] > 0 and w[n + j + 1] < 0);
            const im1 = w[n + j] / w[2 * n + j];
            const im2 = w[n + j + 1] / w[2 * n + j + 1];
            const re1 = w[j] / w[2 * n + j];
            const re2 = w[j + 1] / w[2 * n + j + 1];
            try testing.expect(@abs(re1 - re2) + @abs(im1 + im2) <= 10 * tol * (1 + @abs(re1) + @abs(im1)));
        }
        for (0..@as(usize, if (pair) 2 else 1)) |k| {
            const alpha = Z.init(w[j + k], w[n + j + k]);
            const beta = Z.init(w[2 * n + j + k], 0);
            for ([_][]const T{ vr, vl }, [_]bool{ false, true }) |v, left| {
                // a pair's `v(:, j) ± i v(:, j + 1)`
                var nrm: T = 0;
                for (zv, 0..) |*e, i| {
                    const re = v[i + j * n];
                    const im: T = if (pair) v[i + (j + 1) * n] else 0;
                    e.* = Z.init(re, if (k == 0) im else -im);
                    nrm = std.math.hypot(nrm, std.math.hypot(re, im));
                }
                try testing.expectApproxEqAbs(@as(T, 1), nrm, tol);
                const res = complex_test.residual(Z, n, za[0 .. n * n], za[n * n ..], alpha, beta, zv, left);
                if (!(res <= tol)) check.say("ggevReal {s} n={d} j={d} k={d} pair={} left={}: residual {e}, tol {e}\n", .{ @typeName(T), n, j, k, pair, left, res, tol });
                try testing.expect(res <= tol);
            }
        }
        j += if (pair) 2 else 1;
    }
}

test "ggevReal: real and complex-pair eigenvectors of random pencils" {
    const al = testing.allocator;
    inline for (.{ f32, f64 }) |T| for ([_]usize{ 1, 2, 3, 6, 12 }) |n| {
        const a = try al.alloc(T, n * n);
        defer al.free(a);
        const b = try al.alloc(T, n * n);
        defer al.free(b);
        check.fill(T, a, n * 5 + 3);
        check.fill(T, b, n * 5 + 4);
        try ggevCase(T, n, a, b);
    };
    // a rotation pencil: a pair only; B singular: an infinite eigenvalue
    const c = @cos(@as(f64, 0.5));
    const sn = @sin(@as(f64, 0.5));
    try ggevCase(f64, 2, &.{ c, sn, -sn, c }, &.{ 1, 0, 0, 1 });
    var a: [16]f64 = undefined;
    var b: [16]f64 = undefined;
    check.fill(f64, &a, 51);
    check.fill(f64, &b, 52);
    for (0..4) |i| b[i] = 0;
    try ggevCase(f64, 4, &a, &b);
}

test "ggevReal: a defective pencil's vectors stay finite" {
    const al = testing.allocator;
    const n = 4;
    var a: [n * n]f64 = undefined;
    var b: [n * n]f64 = undefined;
    for (0..n) |j| for (0..n) |i| {
        a[i + j * n] = if (i == j) 2 else if (i < j) 0.5 else 0;
        b[i + j * n] = if (i == j) 1 else 0;
    };
    var w: [3 * n]f64 = undefined;
    var vl: [n * n]f64 = undefined;
    var vr: [n * n]f64 = undefined;
    try ggev.ggevReal(al, check.io, f64, .compute, .compute, n, &a, n, &b, n, w[0..n], w[n .. 2 * n], w[2 * n ..], &vl, n, &vr, n);
    for (vr, vl) |x, y| try testing.expect(std.math.isFinite(x) and std.math.isFinite(y));
}
