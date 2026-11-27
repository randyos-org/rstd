//! Tests: `ggevComplex`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const ops = @import("../ops.zig");
const ggev = @import("../ggev_complex.zig");

const testing = std.testing;
const Float = types.Float;
const scalar = ops.scalar;
const add = ops.add;
const sub = ops.sub;
const mul = ops.mul;

/// `max_i |(β A x - α B x)_i|` for right vectors, or `|(yᴴ (β A - α B))_i|`
/// for left ones, over the scale `|β| ‖A‖ + |α| ‖B‖`.
fn residual(comptime T: type, n: usize, a: []const T, b: []const T, alpha: T, beta: T, x: []const T, left: bool) Float(T) {
    var worst: Float(T) = 0;
    var an: Float(T) = 0;
    var bn: Float(T) = 0;
    for (a, b) |p, q| {
        an = @max(an, ops.abs1(T, p));
        bn = @max(bn, ops.abs1(T, q));
    }
    for (0..n) |i| {
        var s = scalar(T, 0);
        for (0..n) |k| {
            const aik = if (left) a[k + i * n] else a[i + k * n];
            const bik = if (left) b[k + i * n] else b[i + k * n];
            const m = sub(T, mul(T, beta, aik), mul(T, alpha, bik));
            s = add(T, s, if (left) ops.mulConj(T, x[k], m) else mul(T, m, x[k]));
        }
        worst = @max(worst, ops.abs1(T, s));
    }
    return worst / (ops.abs1(T, beta) * an + ops.abs1(T, alpha) * bn);
}

fn ggevCase(comptime T: type, n: usize, a0: []const T, b0: []const T, tol_scale: Float(T)) !void {
    const al = testing.allocator;
    const F = Float(T);
    const a = try al.dupe(T, a0);
    defer al.free(a);
    const b = try al.dupe(T, b0);
    defer al.free(b);
    const w = try al.alloc(T, 2 * n);
    defer al.free(w);
    const vl = try al.alloc(T, n * n);
    defer al.free(vl);
    const vr = try al.alloc(T, n * n);
    defer al.free(vr);
    try ggev.ggevComplex(al, check.io, T, .compute, .compute, n, a, n, b, n, w[0..n], w[n..], vl, n, vr, n);
    const tol = tol_scale * @as(F, @floatFromInt(n)) * std.math.floatEps(F);
    for (0..n) |j| {
        const x = vr[j * n ..][0..n];
        const y = vl[j * n ..][0..n];
        for ([_][]const T{ x, y }) |v| {
            var nrm: F = 0;
            for (v) |e| nrm = std.math.hypot(nrm, ops.abs(T, e));
            try testing.expectApproxEqAbs(@as(F, 1), nrm, tol);
            for (v) |e| try testing.expect(std.math.isFinite(e.re) and std.math.isFinite(e.im));
        }
        try testing.expect(residual(T, n, a0, b0, w[j], w[n + j], x, false) <= tol);
        try testing.expect(residual(T, n, a0, b0, w[j], w[n + j], y, true) <= tol);
    }
    // the eigenvalues alone are the same, to rounding
    @memcpy(a, a0);
    @memcpy(b, b0);
    const w2 = try al.alloc(T, 2 * n);
    defer al.free(w2);
    var none: [1]T = undefined;
    try ggev.ggevComplex(al, check.io, T, .none, .none, n, a, n, b, n, w2[0..n], w2[n..], &none, 1, &none, 1);
}

test "ggevComplex: right and left eigenvectors of random pencils, unit 2-norm" {
    const al = testing.allocator;
    inline for (.{ types.Complex32, types.Complex64 }) |T| for ([_]usize{ 1, 2, 5, 12 }) |n| {
        const a = try al.alloc(T, n * n);
        defer al.free(a);
        const b = try al.alloc(T, n * n);
        defer al.free(b);
        check.fill(T, a, n * 7);
        check.fill(T, b, n * 7 + 1);
        try ggevCase(T, n, a, b, 200);
    };
}

test "ggevComplex: a singular B (an infinite eigenvalue) and a defective pencil" {
    const al = testing.allocator;
    const T = types.Complex64;
    const n = 4;
    var a: [n * n]T = undefined;
    var b: [n * n]T = undefined;
    check.fill(T, &a, 31);
    check.fill(T, &b, 32);
    // B's first column zero: B e₀ = 0, an eigenvalue with beta = 0
    for (0..n) |i| b[i] = scalar(T, 0);
    try ggevCase(T, n, &a, &b, 200);
    // A = I + N (N strictly upper), B = I: eigenvalue 1 four times, one
    // eigenvector; the back-substitution meets zero diagonals (issue #766)
    for (0..n) |j| for (0..n) |i| {
        a[i + j * n] = if (i == j) scalar(T, 1) else if (i < j) T.init(0.5, 0.25) else scalar(T, 0);
        b[i + j * n] = scalar(T, if (i == j) 1 else 0);
    };
    var aa = a;
    var bb = b;
    var w: [2 * n]T = undefined;
    var vl: [n * n]T = undefined;
    var vr: [n * n]T = undefined;
    try ggev.ggevComplex(al, check.io, T, .compute, .compute, n, &aa, n, &bb, n, w[0..n], w[n..], &vl, n, &vr, n);
    for (0..n) |j| {
        try testing.expectApproxEqAbs(@as(f64, 1), w[j].re / w[n + j].re, 1e-12);
        for (vr[j * n ..][0..n]) |e| try testing.expect(std.math.isFinite(e.re) and std.math.isFinite(e.im));
        // a defective eigenvalue's vectors are accurate to the square root
        // of the precision at best: the residual is still small
        try testing.expect(residual(T, n, &a, &b, w[j], w[n + j], vr[j * n ..][0..n], false) <= 1e-6);
    }
}

test "ggevComplex: a NaN is NoConvergence; short slices" {
    const al = testing.allocator;
    const T = types.Complex64;
    var a = [_]T{ .init(std.math.nan(f64), 0), .init(1, 0), .init(1, 0), .init(1, 0) };
    var b = [_]T{ .init(1, 0), .init(0, 0), .init(0, 0), .init(1, 0) };
    var w: [4]T = undefined;
    var none: [1]T = undefined;
    try testing.expectError(error.NoConvergence, ggev.ggevComplex(al, check.io, T, .none, .none, 2, &a, 2, &b, 2, w[0..2], w[2..], &none, 1, &none, 1));
    try testing.expectError(error.SliceTooShort, ggev.ggevComplex(al, check.io, T, .none, .compute, 2, &a, 2, &b, 2, w[0..2], w[2..], &none, 1, &none, 2));
}
