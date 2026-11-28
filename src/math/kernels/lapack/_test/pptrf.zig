//! Tests: `pptrf` and `pptrs`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const ops = @import("../ops.zig");
const cholesky = @import("../cholesky.zig");
const packed_cholesky = @import("../pptrf.zig");

const testing = std.testing;
const Float = types.Float;
const Uplo = types.Uplo;

/// `uplo`'s triangle of the `n x n` `a`, packed by columns.
fn pack(comptime T: type, uplo: Uplo, n: usize, a: []const T, ap: []T) void {
    var k: usize = 0;
    for (0..n) |j| {
        const rows = if (uplo == .u) [2]usize{ 0, j + 1 } else [2]usize{ j, n };
        for (rows[0]..rows[1]) |i| {
            ap[k] = a[i + j * n];
            k += 1;
        }
    }
}

test "pptrf: potrf's factor in packed storage; pptrs solves A X = B" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| for ([_]Uplo{ .u, .l }) |uplo| for ([_]usize{ 1, 2, 7, 30 }) |n| {
        // `A = M Mᴴ + n I`
        const m = try al.alloc(T, n * n);
        defer al.free(m);
        check.fill(T, m, n * 3 + @as(u64, if (uplo == .u) 1 else 0));
        const a = try al.alloc(T, n * n);
        defer al.free(a);
        for (0..n) |j| for (0..n) |i| {
            var t = ops.zero(T);
            for (0..n) |k| t = ops.add(T, t, ops.mul(T, m[i + k * n], ops.conj(T, m[j + k * n])));
            a[i + j * n] = if (i == j) ops.add(T, t, ops.scalar(T, @floatFromInt(n))) else t;
        };
        const ap = try al.alloc(T, n * (n + 1) / 2);
        defer al.free(ap);
        pack(T, uplo, n, a, ap);
        try packed_cholesky.pptrf(al, check.io, T, uplo, n, ap);
        // against `potrf`'s, packed
        const f = try al.dupe(T, a);
        defer al.free(f);
        try cholesky.potrf(al, check.io, T, uplo, n, f, n);
        const fp = try al.alloc(T, ap.len);
        defer al.free(fp);
        pack(T, uplo, n, f, fp);
        var big: Float(T) = 0;
        for (a) |v| big = @max(big, ops.abs(T, v));
        const eps = std.math.floatEps(Float(T));
        for (ap, fp) |x, y| try testing.expect(ops.abs(T, ops.sub(T, x, y)) <= 50 * @as(Float(T), @floatFromInt(n)) * eps * @sqrt(big));
        const nrhs = 3;
        const b = try al.alloc(T, n * nrhs);
        defer al.free(b);
        check.fill(T, b, n + 9);
        const x = try al.dupe(T, b);
        defer al.free(x);
        try packed_cholesky.pptrs(al, check.io, T, uplo, n, nrhs, ap, x, n);
        var worst: Float(T) = 0;
        for (0..nrhs) |c| for (0..n) |i| {
            var s = ops.zero(T);
            for (0..n) |k| s = ops.add(T, s, ops.mul(T, a[i + k * n], x[k + c * n]));
            worst = @max(worst, ops.abs(T, ops.sub(T, s, b[i + c * n])));
        };
        try testing.expect(worst <= 100 * @as(Float(T), @floatFromInt(n * n)) * eps * big);
    };
}

test "pptrf: a minor not positive definite, NaN, bad arguments" {
    const al = testing.allocator;
    // `diag(1, -1)`: the second pivot fails and is stored
    for ([_]Uplo{ .u, .l }) |uplo| {
        var ap = [_]f64{ 1, 0, -1 };
        try testing.expectError(error.NotPositiveDefinite, packed_cholesky.pptrf(al, check.io, f64, uplo, 2, &ap));
        try testing.expectEqualSlices(f64, &.{ 1, 0, -1 }, &ap);
        var nan = [_]f64{ std.math.nan(f64), 0, 1 };
        try testing.expectError(error.NotPositiveDefinite, packed_cholesky.pptrf(al, check.io, f64, uplo, 2, &nan));
    }
    var c = [_]types.Complex64{ .init(4, 3), .init(1, 1), .init(2, 0) };
    try packed_cholesky.pptrf(al, check.io, types.Complex64, .u, 2, &c);
    // the diagonal comes out real; its imaginary part is never read
    try testing.expectEqual(@as(f64, 2), c[0].re);
    try testing.expectEqual(@as(f64, 0), c[0].im);
    var short = [_]f64{ 1, 0 };
    try testing.expectError(error.SliceTooShort, packed_cholesky.pptrf(al, check.io, f64, .u, 2, &short));
    try testing.expectError(error.UnsupportedUplo, packed_cholesky.pptrf(al, check.io, f64, .v, 1, &short));
    var b = [_]f64{1};
    try testing.expectError(error.SliceTooShort, packed_cholesky.pptrs(al, check.io, f64, .l, 2, 1, &short, &b, 2));
}
