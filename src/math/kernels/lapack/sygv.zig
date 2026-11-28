//! `sygv`, `hegv`, `sygvd`, `hegvd`: Reference-LAPACK's `xSYGV`, `xHEGV`,
//! `xSYGVD`, `xHEGVD` (`lapack-netlib/SRC/dsygv.f` and the rest, OpenBLAS
//! v0.3.34-317): `potrf` of `B`, `sygst` (`hegst`) to standard form, `syev`
//! or `syevd`, the eigenvectors brought back by `trsm` or `trmm`.

const std = @import("std");
const types = @import("../blas/types.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const relapack = @import("relapack.zig");
const evd = @import("evd.zig");
const lapack = @import("__root__.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Trans = types.Trans;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;

pub const GeneralizedProblem = relapack.GeneralizedProblem;
pub const VectorJob = interface.VectorJob;

/// What the four return: a reversed `uplo`, a slice too short (`a`, `b`, or
/// `w` under `n`) or a leading dimension too small; `NotPositiveDefinite`
/// when `B` is not (the reference's `INFO > N`), `B` then part factored;
/// `NoConvergence` from the eigensolver, also for a NaN or infinity.
pub const SygvError = interface.TriangleError || error{ NotPositiveDefinite, NoConvergence } || Allocator.Error || Io.Cancelable;

/// `?SYGV`: the eigenvalues of the symmetric-definite `problem` (`A x = λ B
/// x`, `A B x = λ x` or `B A x = λ x`; `A` symmetric, `B` symmetric positive
/// definite, both from their `uplo` triangles) into `w` (nondecreasing); with
/// `jobz` `.compute` the eigenvectors over `A`, normalized as the
/// reference's (`Xᵀ B X = I`, or `Xᵀ B⁻¹ X = I` for `B A x = λ x`). `B`
/// takes its Cholesky factor. Real only; QR iteration.
pub fn sygv(allocator: Allocator, io: Io, comptime T: type, problem: GeneralizedProblem, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, b: []T, ldb: usize, w: []Float(T)) SygvError!void {
    comptime ops.requireReal(T, "sygv");
    return solve(allocator, io, T, false, problem, jobz, uplo, n, a, lda, b, ldb, w);
}

/// `?SYGVD`: `sygv` by divide and conquer (`syevd`).
pub fn sygvd(allocator: Allocator, io: Io, comptime T: type, problem: GeneralizedProblem, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, b: []T, ldb: usize, w: []Float(T)) SygvError!void {
    comptime ops.requireReal(T, "sygvd");
    return solve(allocator, io, T, true, problem, jobz, uplo, n, a, lda, b, ldb, w);
}

/// `?HEGV`: `sygv` for Hermitian `A` and Hermitian positive definite `B`
/// (`Xᴴ B X = I`). Complex only.
pub fn hegv(allocator: Allocator, io: Io, comptime T: type, problem: GeneralizedProblem, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, b: []T, ldb: usize, w: []Float(T)) SygvError!void {
    comptime ops.requireComplex(T, "hegv");
    return solve(allocator, io, T, false, problem, jobz, uplo, n, a, lda, b, ldb, w);
}

/// `?HEGVD`: `hegv` by divide and conquer (`heevd`).
pub fn hegvd(allocator: Allocator, io: Io, comptime T: type, problem: GeneralizedProblem, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, b: []T, ldb: usize, w: []Float(T)) SygvError!void {
    comptime ops.requireComplex(T, "hegvd");
    return solve(allocator, io, T, true, problem, jobz, uplo, n, a, lda, b, ldb, w);
}

fn solve(allocator: Allocator, io: Io, comptime T: type, comptime divide: bool, problem: GeneralizedProblem, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, b: []T, ldb: usize, w: []Float(T)) SygvError!void {
    try interface.holdsTriangle(uplo, n, a.len, lda);
    try interface.holdsTriangle(uplo, n, b.len, ldb);
    if (w.len < n) return error.SliceTooShort;
    if (n == 0) return;
    lapack.potrf(allocator, io, T, uplo, n, b, ldb) catch |e| return ops.narrow(SygvError, e);
    const gst = if (comptime isComplex(T)) relapack.hegst else relapack.sygst;
    try gst(allocator, io, T, problem, uplo, n, a, lda, b, ldb);
    const ev = if (divide) evd.syevd else evd.syev;
    try ev(allocator, io, T, jobz, uplo, n, a, lda, w);
    if (jobz == .none) return;
    const ct: Trans = if (comptime isComplex(T)) .c else .t;
    const upper = uplo == .u;
    switch (problem) {
        // `x = L⁻ᴴ y` or `U⁻¹ y`
        .ax_equals_lambda_bx, .abx_equals_lambda_x => try relapack.trsm(allocator, io, T, true, .l, uplo, if (upper) .n else ct, .n, n, n, scalar(T, 1), b.ptr, ldb, a.ptr, lda),
        // `x = L y` or `Uᴴ y`
        .bax_equals_lambda_x => try relapack.trsm(allocator, io, T, false, .l, uplo, if (upper) ct else .n, .n, n, n, scalar(T, 1), b.ptr, ldb, a.ptr, lda),
    }
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");

/// `y = M x` for the `n x n` `m`.
fn matVec(comptime T: type, n: usize, m: []const T, x: []const T, y: []T) void {
    for (0..n) |i| {
        var s = scalar(T, 0);
        for (0..n) |k| s = ops.add(T, s, ops.mul(T, m[i + k * n], x[k]));
        y[i] = s;
    }
}

fn maxAbs1(comptime T: type, v: []const T) Float(T) {
    var r: Float(T) = 0;
    for (v) |x| r = @max(r, ops.abs1(T, x));
    return r;
}

fn sygvCases(comptime T: type, comptime divide: bool) !void {
    const al = testing.allocator;
    const herm = comptime isComplex(T);
    const f = if (herm) (if (divide) hegvd else hegv) else (if (divide) sygvd else sygv);
    const F = Float(T);
    for ([_]usize{ 1, 4, 17 }) |n| for ([_]Uplo{ .l, .u }) |uplo| for ([_]GeneralizedProblem{ .ax_equals_lambda_bx, .abx_equals_lambda_x, .bax_equals_lambda_x }) |problem| {
        const a0 = try al.alloc(T, n * n);
        defer al.free(a0);
        const m = try al.alloc(T, n * n);
        defer al.free(m);
        check.fill(T, a0, n * 3 + 1);
        check.fill(T, m, n * 5 + 2);
        const b0 = try al.alloc(T, n * n);
        defer al.free(b0);
        // A Hermitian; B = M Mᴴ + n I
        for (0..n) |j| for (0..n) |i| {
            if (i > j) a0[i + j * n] = ops.conj(T, a0[j + i * n]);
            if (i == j) a0[i + i * n] = scalar(T, ops.re(T, a0[i + i * n]));
            var s = scalar(T, if (i == j) @floatFromInt(n) else 0);
            for (0..n) |k| s = ops.add(T, s, ops.mul(T, m[i + k * n], ops.conj(T, m[j + k * n])));
            b0[i + j * n] = s;
        };
        const a = try al.dupe(T, a0);
        defer al.free(a);
        const b = try al.dupe(T, b0);
        defer al.free(b);
        const w = try al.alloc(F, n);
        defer al.free(w);
        try f(al, check.io, T, problem, .compute, uplo, n, a, n, b, n, w);
        for (1..n) |i| try testing.expect(w[i - 1] <= w[i]);
        const t1 = try al.alloc(T, n);
        defer al.free(t1);
        const t2 = try al.alloc(T, n);
        defer al.free(t2);
        const eps = std.math.floatEps(F);
        const scale = maxAbs1(T, a0) + maxAbs1(T, b0) * (@abs(w[0]) + @abs(w[n - 1]));
        for (0..n) |j| {
            const x = a[j * n ..][0..n];
            const lam = scalar(T, w[j]);
            // left `A x` (`A B x`, `B A x`); right `λ B x` (`λ x`)
            switch (problem) {
                .ax_equals_lambda_bx => {
                    matVec(T, n, a0, x, t1);
                    matVec(T, n, b0, x, t2);
                    for (t2) |*v| v.* = ops.mul(T, lam, v.*);
                },
                .abx_equals_lambda_x => {
                    matVec(T, n, b0, x, t2);
                    matVec(T, n, a0, t2, t1);
                    for (t2, x) |*v, xi| v.* = ops.mul(T, lam, xi);
                },
                .bax_equals_lambda_x => {
                    matVec(T, n, a0, x, t2);
                    matVec(T, n, b0, t2, t1);
                    for (t2, x) |*v, xi| v.* = ops.mul(T, lam, xi);
                },
            }
            var r: F = 0;
            for (t1, t2) |p, q| r = @max(r, ops.abs1(T, ops.sub(T, p, q)));
            const bound = 200 * @as(F, @floatFromInt(n)) * eps * scale * @max(maxAbs1(T, x), 1) * (1 + maxAbs1(T, b0));
            if (!(r <= bound)) {
                check.say("{s} n={d} {t} {t}: eigenpair {d} residual {e}, bound {e}\n", .{ @typeName(T), n, uplo, problem, j, r, bound });
                return error.TestUnexpectedResult;
            }
        }
    };
}

test "sygv, sygvd, hegv, hegvd: every problem type's eigenpairs, both triangles" {
    inline for (.{ f32, f64, types.Complex64 }) |T| inline for (.{ false, true }) |divide| try sygvCases(T, divide);
}

test "sygv: B not positive definite; short slices" {
    const al = testing.allocator;
    var a = [_]f64{ 1, 0, 0, 1 };
    var b = [_]f64{ -1, 0, 0, 1 };
    var w: [2]f64 = undefined;
    try testing.expectError(error.NotPositiveDefinite, sygv(al, check.io, f64, .ax_equals_lambda_bx, .none, .l, 2, &a, 2, &b, 2, &w));
    try testing.expectError(error.SliceTooShort, sygv(al, check.io, f64, .ax_equals_lambda_bx, .none, .l, 2, &a, 2, &b, 2, w[0..1]));
}
