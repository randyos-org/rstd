//! Tests: `sytri` and `hetri`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const ops = @import("../ops.zig");
const sytrf_factor = @import("../sytrf.zig");
const sytri_invert = @import("../sytri.zig");

const testing = std.testing;
const Float = types.Float;
const Uplo = types.Uplo;
const Pivot = sytrf_factor.Pivot;

test "sytri, hetri: A⁻¹ A = I from sytrf's and hetrf's factors, 2 x 2 blocks too" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| inline for (.{ false, true }) |herm| {
        if (herm and !comptime types.isComplex(T)) continue;
        for ([_]Uplo{ .u, .l }) |uplo| for ([_]usize{ 1, 2, 6, 25 }) |n| for ([_]bool{ false, true }) |zero_diag| {
            // symmetric (or Hermitian), indefinite; a zero diagonal forces
            // 2 x 2 blocks
            const a = try al.alloc(T, n * n);
            defer al.free(a);
            check.fill(T, a, n * 5 + @as(u64, @intFromBool(zero_diag)));
            for (0..n) |j| for (0..j + 1) |i| {
                a[j + i * n] = if (herm) ops.conj(T, a[i + j * n]) else a[i + j * n];
            };
            for (0..n) |i| a[i + i * n] = if (zero_diag and n > 1) ops.zero(T) else if (herm) ops.scalar(T, ops.re(T, a[i + i * n])) else a[i + i * n];
            const f = try al.dupe(T, a);
            defer al.free(f);
            const ipiv = try al.alloc(Pivot, n);
            defer al.free(ipiv);
            const factor = if (herm) sytrf_factor.hetrf else sytrf_factor.sytrf;
            try testing.expectEqual(@as(?usize, null), try factor(al, check.io, T, uplo, n, f, n, ipiv));
            if (herm) try sytri_invert.hetri(al, check.io, T, uplo, n, f, n, ipiv) else try sytri_invert.sytri(al, check.io, T, uplo, n, f, n, ipiv);
            // the inverse's other triangle from its `uplo` one
            for (0..n) |j| for (0..n) |i| {
                const stored = if (uplo == .u) i <= j else i >= j;
                if (!stored) f[i + j * n] = if (herm) ops.conj(T, f[j + i * n]) else f[j + i * n];
            };
            if (herm) for (0..n) |i| try testing.expectEqual(@as(Float(T), 0), f[i + i * n].im);
            var worst: Float(T) = 0;
            var big: Float(T) = 0;
            var ibig: Float(T) = 0;
            for (a, f) |x, y| {
                big = @max(big, ops.abs(T, x));
                ibig = @max(ibig, ops.abs(T, y));
            }
            for (0..n) |j| for (0..n) |i| {
                var s = ops.zero(T);
                for (0..n) |k| s = ops.add(T, s, ops.mul(T, f[i + k * n], a[k + j * n]));
                worst = @max(worst, ops.abs(T, ops.sub(T, s, if (i == j) ops.one(T) else ops.zero(T))));
            };
            const tol = 200 * @as(Float(T), @floatFromInt(n * n)) * std.math.floatEps(Float(T)) * big * ibig;
            if (!(worst <= tol)) {
                check.say("sytri {s} herm={} {t} n={d} zero_diag={}: |A⁻¹ A - I| {e}, tol {e}\n", .{ @typeName(T), herm, uplo, n, zero_diag, worst, tol });
                return error.TestUnexpectedResult;
            }
        };
    };
}

test "sytri: a zero 1 x 1 block, bad pivots and shapes" {
    const al = testing.allocator;
    var a = [_]f64{ 0, 0, 0, 1 };
    const one_one = [_]Pivot{ .of(0, .one), .of(1, .one) };
    try testing.expectError(error.Singular, sytri_invert.sytri(al, check.io, f64, .l, 2, &a, 2, &one_one));
    try testing.expectEqualSlices(f64, &.{ 0, 0, 0, 1 }, &a);
    const bad = [_]Pivot{ .of(0, .one), .of(2, .one) };
    try testing.expectError(error.PivotOutOfRange, sytri_invert.sytri(al, check.io, f64, .u, 2, &a, 2, &bad));
    try testing.expectError(error.SliceTooShort, sytri_invert.sytri(al, check.io, f64, .u, 2, &a, 2, one_one[0..1]));
}
