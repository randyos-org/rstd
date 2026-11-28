//! Tests: `lansy`, `lanhe`, `lantr` against `lange` on the full matrix.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const ops = @import("../ops.zig");
const norm = @import("../norm.zig");
const structured = @import("../norm_structured.zig");

const testing = std.testing;
const Float = types.Float;
const Uplo = types.Uplo;
const NormType = norm.NormType;

const norms = [_]NormType{ .max_abs, .one, .infinity, .frobenius };

fn close(comptime F: type, got: F, want: F) !void {
    const tol = 16 * std.math.floatEps(F) * @max(want, 1);
    if (!(@abs(got - want) <= tol)) {
        check.say("structured norm: {e}, lange {e}\n", .{ got, want });
        return error.TestUnexpectedResult;
    }
}

test "lansy, lanhe: the full matrix's lange, both triangles, every norm" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| inline for (.{ false, true }) |herm| {
        if (herm and !comptime types.isComplex(T)) continue;
        for ([_]usize{ 1, 4, 9 }) |n| for ([_]Uplo{ .u, .l }) |uplo| {
            const a = try al.alloc(T, n * n);
            defer al.free(a);
            check.fill(T, a, n * 3 + 1);
            // the full matrix the triangle stands for
            const full = try al.alloc(T, n * n);
            defer al.free(full);
            for (0..n) |j| for (0..n) |i| {
                const in = if (uplo == .u) i <= j else i >= j;
                const v = if (in) a[i + j * n] else a[j + i * n];
                full[i + j * n] = if (!in and herm) ops.conj(T, v) else v;
                if (herm and i == j) full[i + j * n] = ops.scalar(T, ops.re(T, v));
            };
            for (norms) |kind| {
                const got = if (herm) try structured.lanhe(al, check.io, T, kind, uplo, n, a, n) else try structured.lansy(al, check.io, T, kind, uplo, n, a, n);
                try close(Float(T), got, try norm.lange(al, check.io, T, kind, n, n, full, n));
            }
        };
    };
}

test "lantr: the trapezoid's lange, unit and not, tall and wide" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| for ([_][2]usize{ .{ 1, 1 }, .{ 5, 3 }, .{ 3, 5 }, .{ 6, 6 } }) |s| {
        const m = s[0];
        const n = s[1];
        for ([_]Uplo{ .u, .l }) |uplo| for ([_]types.Diag{ .n, .u }) |diag| {
            const a = try al.alloc(T, m * n);
            defer al.free(a);
            check.fill(T, a, m * 7 + n);
            const full = try al.alloc(T, m * n);
            defer al.free(full);
            for (0..n) |j| for (0..m) |i| {
                const in = if (uplo == .u) i <= j else i >= j;
                full[i + j * m] = if (i == j and diag == .u) ops.scalar(T, 1) else if (in) a[i + j * m] else ops.zero(T);
            };
            for (norms) |kind| try close(Float(T), try structured.lantr(al, check.io, T, kind, uplo, diag, m, n, a, m), try norm.lange(al, check.io, T, kind, m, n, full, m));
        };
    };
}

test "lansy, lantr: NaN, empty, bad shapes" {
    const al = testing.allocator;
    var a = [_]f64{ 1, std.math.nan(f64), 2, 3 };
    try testing.expect(std.math.isNan(try structured.lansy(al, check.io, f64, .one, .l, 2, &a, 2)));
    try testing.expect(std.math.isNan(try structured.lantr(al, check.io, f64, .max_abs, .l, .n, 2, 2, &a, 2)));
    try testing.expectEqual(0, try structured.lansy(al, check.io, f64, .frobenius, .u, 0, a[0..0], 1));
    try testing.expectEqual(0, try structured.lantr(al, check.io, f64, .one, .u, .u, 0, 3, a[0..0], 1));
    try testing.expectError(error.LeadingDimensionTooSmall, structured.lansy(al, check.io, f64, .one, .u, 2, &a, 1));
    try testing.expectError(error.SliceTooShort, structured.lantr(al, check.io, f64, .one, .u, .n, 2, 3, &a, 2));
}
