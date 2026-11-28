//! `?PBTRS`: solves with `pbtrf`'s band Cholesky factor, Reference-LAPACK
//! 3.12.0 (`lapack-netlib/SRC/dpbtrs.f`, `zpbtrs.f`, OpenBLAS v0.3.34-317;
//! modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;
const Uplo = types.Uplo;

/// What `pbtrs` returns: a `uplo` that is neither triangle, a slice too
/// short or a leading dimension too small (`ldab` under `kd + 1`), and the
/// allocator's and cancellation.
pub const PbtrsError = interface.UploError || interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `?PBTRS` (`dpbtrs.f`, `zpbtrs.f`): solves `A X = B` for the `n x nrhs`
/// `X`, in place of `B`, given `pbtrf`'s `Uᴴ U` or `L Lᴴ` of the positive
/// definite band `A` (`kd` off-diagonals) in `ab`: two `tbsv` solves per
/// column of `B`.
pub fn pbtrs(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, kd: usize, nrhs: usize, ab: []const T, ldab: usize, b: []T, ldb: usize) PbtrsError!void {
    comptime ops.requireFloat(T, "pbtrs");
    if (uplo.isRev()) return error.UnsupportedUplo;
    if (ldab < kd + 1) return error.LeadingDimensionTooSmall;
    try util.holds(ab.len, kd + 1, n, ldab);
    try util.holds(b.len, n, nrhs, ldb);
    if (n == 0 or nrhs == 0) return;
    const upper = uplo == .u;
    const adjoint: Trans = if (comptime types.isComplex(T)) .c else .t;
    for (0..nrhs) |j| {
        const x = b.ptr + j * ldb;
        // `Uᴴ Y = B`, then `U X = Y`; or `L Y = B`, then `Lᴴ X = Y`
        try blas_calls.tbsv(allocator, io, T, upper, if (upper) adjoint else .n, false, n, kd, ab.ptr, ldab, x, 1);
        try blas_calls.tbsv(allocator, io, T, upper, if (upper) .n else adjoint, false, n, kd, ab.ptr, ldab, x, 1);
    }
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const relapack_pbtrf = @import("relapack_pbtrf.zig");

test "pbtrs: A X = B from pbtrf's band factor, both triangles, every type" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| for ([_][2]usize{ .{ 1, 0 }, .{ 6, 1 }, .{ 15, 3 }, .{ 9, 0 } }) |s| for ([_]Uplo{ .u, .l }) |uplo| {
        const F = types.Float(T);
        const n = s[0];
        const kd = s[1];
        const ld = kd + 1;
        // a diagonally dominant Hermitian band, full and packed
        const full = try al.alloc(T, n * n);
        defer al.free(full);
        @memset(full, ops.zero(T));
        const vals = try al.alloc(T, n * n);
        defer al.free(vals);
        check.fill(T, vals, n * 7 + kd);
        for (0..n) |j| for (j + 1..@min(n, j + kd + 1)) |i| {
            full[i + j * n] = vals[i + j * n];
            full[j + i * n] = ops.conj(T, vals[i + j * n]);
        };
        for (0..n) |i| {
            var r: F = 1;
            for (0..n) |k| r += ops.abs(T, full[i + k * n]);
            full[i + i * n] = ops.scalar(T, r);
        }
        const ab = try al.alloc(T, ld * n);
        defer al.free(ab);
        @memset(ab, ops.zero(T));
        for (0..n) |j| for (0..n) |i| {
            if (uplo == .u and i <= j and j - i <= kd) ab[kd + i - j + j * ld] = full[i + j * n];
            if (uplo == .l and i >= j and i - j <= kd) ab[i - j + j * ld] = full[i + j * n];
        };
        try relapack_pbtrf.pbtrf(al, check.io, T, uplo, n, kd, ab, ld);
        const nrhs = 2;
        const b = try al.alloc(T, n * nrhs);
        defer al.free(b);
        check.fill(T, b, n + 3);
        const x = try al.dupe(T, b);
        defer al.free(x);
        try pbtrs(al, check.io, T, uplo, n, kd, nrhs, ab, ld, x, n);
        var worst: F = 0;
        for (0..nrhs) |c| for (0..n) |i| {
            var t = ops.zero(T);
            for (0..n) |k| t = ops.add(T, t, ops.mul(T, full[i + k * n], x[k + c * n]));
            worst = @max(worst, ops.abs(T, ops.sub(T, t, b[i + c * n])));
        };
        try testing.expect(worst <= 100 * @as(F, @floatFromInt(n * (kd + 1))) * std.math.floatEps(F) * @as(F, @floatFromInt(2 * kd + 2)));
    };
    var ab = [_]f64{ 1, 2 };
    var b = [_]f64{ 1, 1 };
    try testing.expectError(error.LeadingDimensionTooSmall, pbtrs(al, check.io, f64, .u, 2, 1, 1, &ab, 1, &b, 2));
    try testing.expectError(error.SliceTooShort, pbtrs(al, check.io, f64, .u, 2, 1, 1, &ab, 2, &b, 2));
}
