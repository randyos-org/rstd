//! Tests: `trcon` and `pocon` against the exact `1 / (‖A‖ ‖A⁻¹‖)`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const ops = @import("../ops.zig");
const norm = @import("../norm.zig");
const cholesky = @import("../cholesky.zig");
const estimates = @import("../condition_estimates.zig");

const testing = std.testing;
const Float = types.Float;
const Uplo = types.Uplo;

/// `A⁻¹` of the full `n x n` `a` by Gauss-Jordan without pivoting (the
/// test matrices are diagonally dominant), packed.
fn inverse(comptime T: type, al: std.mem.Allocator, n: usize, a: []const T) ![]T {
    const w = try al.dupe(T, a);
    defer al.free(w);
    const x = try al.alloc(T, n * n);
    for (0..n) |j| for (0..n) |i| {
        x[i + j * n] = ops.scalar(T, if (i == j) 1 else 0);
    };
    for (0..n) |k| {
        const p = w[k + k * n];
        for (0..n) |j| {
            w[k + j * n] = ops.div(T, w[k + j * n], p);
            x[k + j * n] = ops.div(T, x[k + j * n], p);
        }
        for (0..n) |i| if (i != k) {
            const f = w[i + k * n];
            for (0..n) |j| {
                w[i + j * n] = ops.sub(T, w[i + j * n], ops.mul(T, f, w[k + j * n]));
                x[i + j * n] = ops.sub(T, x[i + j * n], ops.mul(T, f, x[k + j * n]));
            }
        };
    }
    return x;
}

/// The estimate is at least the exact value (`lacn2` under-estimates
/// `‖A⁻¹‖`) and within a factor of 10 of it.
fn near(comptime F: type, est: F, exact: F) !void {
    if (!(est >= exact * (1 - 1e-3) and est <= 10 * exact)) {
        check.say("condition estimate {e}, exact {e}\n", .{ est, exact });
        return error.TestUnexpectedResult;
    }
}

test "trcon: within lacn2's reach of the exact value, every type, uplo, diag and norm" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| for ([_]usize{ 1, 3, 12 }) |n| for ([_]Uplo{ .u, .l }) |uplo| for ([_]types.Diag{ .n, .u }) |diag| {
        const F = Float(T);
        const a = try al.alloc(T, n * n);
        defer al.free(a);
        check.fill(T, a, n * 13 + 5);
        // the full triangle `trcon` sees, its diagonal shifted away from 0
        const full = try al.alloc(T, n * n);
        defer al.free(full);
        for (0..n) |j| for (0..n) |i| {
            const in = if (uplo == .u) i < j else i > j;
            full[i + j * n] = if (i == j) (if (diag == .u) ops.scalar(T, 1) else ops.add(T, a[i + j * n], ops.scalar(T, 2))) else if (in) ops.mulReal(T, a[i + j * n], 0.3) else ops.zero(T);
            if (in or (i == j and diag == .n)) a[i + j * n] = full[i + j * n];
        };
        const inv = try inverse(T, al, n, full);
        defer al.free(inv);
        for ([_]norm.NormType{ .one, .infinity }) |kind| {
            const exact = 1 / (try norm.lange(al, check.io, T, kind, n, n, full, n) * try norm.lange(al, check.io, T, kind, n, n, inv, n));
            try near(F, try estimates.trcon(al, check.io, T, kind, uplo, diag, n, a, n), exact);
        }
    };
}

test "pocon: within lacn2's reach of the exact value; edge cases" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| for ([_]usize{ 1, 4, 15 }) |n| for ([_]Uplo{ .u, .l }) |uplo| {
        const F = Float(T);
        // `A = M Mᴴ + n I`, its one-norm, its factor
        const m = try al.alloc(T, n * n);
        defer al.free(m);
        check.fill(T, m, n * 17 + 2);
        const a = try al.alloc(T, n * n);
        defer al.free(a);
        for (0..n) |j| for (0..n) |i| {
            var t = ops.zero(T);
            for (0..n) |k| t = ops.add(T, t, ops.mul(T, m[i + k * n], ops.conj(T, m[j + k * n])));
            a[i + j * n] = if (i == j) ops.add(T, t, ops.scalar(T, @floatFromInt(n))) else t;
        };
        const anorm = try norm.lange(al, check.io, T, .one, n, n, a, n);
        const inv = try inverse(T, al, n, a);
        defer al.free(inv);
        const exact = 1 / (anorm * try norm.lange(al, check.io, T, .one, n, n, inv, n));
        const f = try al.dupe(T, a);
        defer al.free(f);
        try cholesky.potrf(al, check.io, T, uplo, n, f, n);
        try near(F, try estimates.pocon(al, check.io, T, uplo, n, f, n, anorm), exact);
    };
    var a = [_]f64{ 2, 0, 0, 3 };
    try testing.expectEqual(1, try estimates.pocon(testing.allocator, check.io, f64, .u, 0, a[0..0], 1, 1));
    try testing.expectEqual(0, try estimates.pocon(testing.allocator, check.io, f64, .u, 2, &a, 2, 0));
    try testing.expectError(error.NormOutOfRange, estimates.pocon(testing.allocator, check.io, f64, .u, 2, &a, 2, -1));
    try testing.expectError(error.UnsupportedNorm, estimates.trcon(testing.allocator, check.io, f64, .frobenius, .u, .n, 2, &a, 2));
    // a zero triangle: rcond 0
    var z = [_]f64{ 0, 0, 0, 0 };
    try testing.expectEqual(0, try estimates.trcon(testing.allocator, check.io, f64, .one, .u, .n, 2, &z, 2));
}
