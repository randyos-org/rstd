//! LAPACK's simple drivers, a factorization and then its solve:
//! Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dposv.f`, `dsysv.f`,
//! `zhesv.f`, OpenBLAS v0.3.34-317; modified BSD).

const std = @import("std");
const options = @import("lapack_options");
const util = @import("../blas/util.zig");
const types = @import("../blas/types.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const cholesky = @import("cholesky.zig");
const relapack = @import("relapack.zig");
const sytrf_factor = @import("sytrf.zig");
const sytrs_solve = @import("sytrs.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Pivot = sytrf_factor.Pivot;

/// The `potrf` the package exports (ReLAPACK's with `-Drelapack`).
const potrf = if (options.relapack) relapack.potrf else cholesky.potrf;

/// What `posv` returns: `potrf`'s (`NotPositiveDefinite` when a leading
/// minor is not, `B` untouched then), which hold `potrs`'s.
pub const PosvError = cholesky.PotrfError;

/// `?POSV` (`dposv.f`, `zposv.f`): solves `A X = B` for the symmetric
/// (Hermitian) positive definite `n x n` `A` in `a`'s `uplo` triangle and
/// the `n x nrhs` `B`: `potrf` leaves the factor in `a`, then `potrs`
/// leaves `X` in `b`. Every shape is checked before either runs.
pub fn posv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize) PosvError!void {
    comptime ops.requireFloat(T, "posv");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    try util.holds(b.len, n, nrhs, ldb);
    try potrf(allocator, io, T, uplo, n, a, lda);
    try cholesky.potrs(allocator, io, T, uplo, n, nrhs, a, lda, b, ldb);
}

/// What `sysv` and `hesv` return: `sytrf`'s and `sytrs`'s shapes, and
/// `Singular` for an exactly zero block of `D` (the reference's `INFO >
/// 0`), the factor written and `B` untouched.
pub const SysvError = sytrf_factor.SytrfError || sytrs_solve.SytrsError || interface.SingularError;

/// `?SYSV` (`dsysv.f`, `zsysv.f`): solves `A X = B` for the symmetric
/// (complex: complex symmetric) `n x n` `A` in `a`'s `uplo` triangle:
/// `sytrf`'s `U D Uᵀ` or `L D Lᵀ` with its `ipiv` (`n`), then `sytrs`. The
/// reference takes `?SYTRS2` when its workspace allows, the same solution
/// in other rounding.
pub fn sysv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, nrhs: usize, a: []T, lda: usize, ipiv: []Pivot, b: []T, ldb: usize) SysvError!void {
    comptime ops.requireFloat(T, "sysv");
    return indefinite(allocator, io, T, false, uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

/// `?HESV` (`zhesv.f`, `chesv.f`): `sysv` for the Hermitian `A`, by `hetrf`
/// and `hetrs`. Complex only.
pub fn hesv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, nrhs: usize, a: []T, lda: usize, ipiv: []Pivot, b: []T, ldb: usize) SysvError!void {
    comptime ops.requireComplex(T, "hesv");
    return indefinite(allocator, io, T, true, uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

fn indefinite(allocator: Allocator, io: Io, comptime T: type, comptime hermitian: bool, uplo: Uplo, n: usize, nrhs: usize, a: []T, lda: usize, ipiv: []Pivot, b: []T, ldb: usize) SysvError!void {
    try interface.holdsTriangle(uplo, n, a.len, lda);
    try util.holds(b.len, n, nrhs, ldb);
    if (ipiv.len < n) return error.SliceTooShort;
    const factor = if (hermitian) sytrf_factor.hetrf else sytrf_factor.sytrf;
    if (try factor(allocator, io, T, uplo, n, a, lda, ipiv) != null) return error.Singular;
    const solve = if (hermitian) sytrs_solve.hetrs else sytrs_solve.sytrs;
    try solve(T, uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");

/// `‖A X - B‖` over the entries, `A` Hermitian from its full storage.
fn residual(comptime T: type, n: usize, nrhs: usize, a: []const T, x: []const T, b: []const T) types.Float(T) {
    var worst: types.Float(T) = 0;
    for (0..nrhs) |j| for (0..n) |i| {
        var s = ops.zero(T);
        for (0..n) |k| s = ops.add(T, s, ops.mul(T, a[i + k * n], x[k + j * n]));
        worst = @max(worst, ops.abs(T, ops.sub(T, s, b[i + j * n])));
    };
    return worst;
}

test "posv: A X = B for positive definite A, both triangles; not positive definite" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| for ([_]Uplo{ .u, .l }) |uplo| for ([_][2]usize{ .{ 1, 1 }, .{ 5, 3 }, .{ 40, 2 } }) |s| {
        const n = s[0];
        const nrhs = s[1];
        // `A = M Mᴴ + n I`
        const m = try al.alloc(T, n * n);
        defer al.free(m);
        check.fill(T, m, n * 3 + nrhs);
        const a = try al.alloc(T, n * n);
        defer al.free(a);
        for (0..n) |j| for (0..n) |i| {
            var t = ops.zero(T);
            for (0..n) |k| t = ops.add(T, t, ops.mul(T, m[i + k * n], ops.conj(T, m[j + k * n])));
            a[i + j * n] = if (i == j) ops.add(T, t, ops.scalar(T, @floatFromInt(n))) else t;
        };
        const b = try al.alloc(T, n * nrhs);
        defer al.free(b);
        check.fill(T, b, n + 7);
        const f = try al.dupe(T, a);
        defer al.free(f);
        const x = try al.dupe(T, b);
        defer al.free(x);
        try posv(al, check.io, T, uplo, n, nrhs, f, n, x, n);
        var big: types.Float(T) = 0;
        for (a) |v| big = @max(big, ops.abs(T, v));
        const tol = 100 * @as(types.Float(T), @floatFromInt(n * n)) * std.math.floatEps(types.Float(T)) * big;
        try testing.expect(residual(T, n, nrhs, a, x, b) <= tol);
    };
    // `diag(1, -1)`: the second minor fails, `B` is left as it was
    var a = [_]f64{ 1, 0, 0, -1 };
    var b = [_]f64{ 3, 4 };
    try testing.expectError(error.NotPositiveDefinite, posv(al, check.io, f64, .u, 2, 1, &a, 2, &b, 2));
    try testing.expectEqualSlices(f64, &.{ 3, 4 }, &b);
    try testing.expectError(error.SliceTooShort, posv(al, check.io, f64, .l, 2, 2, &a, 2, &b, 2));
}

test "sysv, hesv: A X = B for indefinite A, both triangles; a singular A" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| inline for (.{ false, true }) |herm| {
        if (herm and !comptime types.isComplex(T)) continue;
        for ([_]Uplo{ .u, .l }) |uplo| for ([_][2]usize{ .{ 1, 1 }, .{ 6, 2 }, .{ 37, 3 } }) |s| {
            const n = s[0];
            const nrhs = s[1];
            // symmetric (or Hermitian) and indefinite
            const a = try al.alloc(T, n * n);
            defer al.free(a);
            check.fill(T, a, n * 5 + nrhs);
            for (0..n) |j| for (0..j + 1) |i| {
                a[j + i * n] = if (herm) ops.conj(T, a[i + j * n]) else a[i + j * n];
            };
            if (herm) for (0..n) |i| {
                a[i + i * n] = ops.scalar(T, ops.re(T, a[i + i * n]));
            };
            const b = try al.alloc(T, n * nrhs);
            defer al.free(b);
            check.fill(T, b, n + 11);
            const f = try al.dupe(T, a);
            defer al.free(f);
            const x = try al.dupe(T, b);
            defer al.free(x);
            const ipiv = try al.alloc(Pivot, n);
            defer al.free(ipiv);
            if (herm) try hesv(al, check.io, T, uplo, n, nrhs, f, n, ipiv, x, n) else try sysv(al, check.io, T, uplo, n, nrhs, f, n, ipiv, x, n);
            var big: types.Float(T) = 0;
            for (a) |v| big = @max(big, ops.abs(T, v));
            var xbig: types.Float(T) = 0;
            for (x) |v| xbig = @max(xbig, ops.abs(T, v));
            const tol = 1000 * @as(types.Float(T), @floatFromInt(n * n)) * std.math.floatEps(types.Float(T)) * big * @max(xbig, 1);
            const res = residual(T, n, nrhs, a, x, b);
            if (!(res <= tol)) {
                check.say("sysv {s} herm={} n={d}: residual {e}, tol {e}\n", .{ @typeName(T), herm, n, res, tol });
                return error.TestUnexpectedResult;
            }
        };
    };
    var z = [_]f64{ 0, 0, 0, 0 };
    var b = [_]f64{ 1, 2 };
    var ipiv: [2]Pivot = undefined;
    try testing.expectError(error.Singular, sysv(al, check.io, f64, .l, 2, 1, &z, 2, &ipiv, &b, 2));
    try testing.expectEqualSlices(f64, &.{ 1, 2 }, &b);
    try testing.expectError(error.SliceTooShort, sysv(al, check.io, f64, .l, 2, 1, &z, 2, ipiv[0..1], &b, 2));
}
