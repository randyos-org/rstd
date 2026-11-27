//! `?GGES` without eigenvalue ordering: the generalized Schur form from
//! `ggev_complex.zig`'s `reduce` and the QZ, as faer v0.24.4's
//! `compute_gevd_generic` builds it before its eigenvectors
//! (`faer/src/linalg/gevd/mod.rs`, MIT).

const std = @import("std");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const ggev_complex = @import("ggev_complex.zig");
const qz_complex = @import("qz_complex.zig");
const qz_real = @import("qz_real.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;

pub const VectorJob = ggev_complex.VectorJob;
pub const GgevError = ggev_complex.GgevError;

/// `?GGES`, complex: `(A, B) = (Q S Zᴴ, Q T Zᴴ)` with `S` and `T` upper
/// triangular (in `a` and `b`; `T`'s diagonal real and nonnegative), the
/// eigenvalues `alpha[j] / beta[j]` their diagonals, and the Schur vectors
/// `Q` in `vsl` and `Z` in `vsr` as `jobvsl` and `jobvsr` ask.
pub fn ggesComplex(allocator: Allocator, io: Io, comptime T: type, jobvsl: VectorJob, jobvsr: VectorJob, n: usize, a: []T, lda: usize, b: []T, ldb: usize, alpha: []T, beta: []T, vsl: []T, ldvsl: usize, vsr: []T, ldvsr: usize) GgevError!void {
    comptime ops.requireComplex(T, "ggesComplex");
    if (alpha.len < n or beta.len < n) return error.SliceTooShort;
    const left = jobvsl == .compute;
    const right = jobvsr == .compute;
    try check(T, n, a, lda, b, ldb, left, vsl, ldvsl, right, vsr, ldvsr);
    if (n == 0) return;
    try ggev_complex.reduce(allocator, io, T, left, right, n, a, lda, b, ldb, vsl, ldvsl, vsr, ldvsr);
    qz_complex.hgeqzComplex(T, .schur_form, if (left) .accumulate else .none, if (right) .accumulate else .none, n, 0, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr) catch |e| return ops.narrow(GgevError, e);
}

/// `?GGES`, real: as `ggesComplex` with the real generalized Schur form, `S`
/// quasi-triangular (a 2 x 2 block per complex pair, `T`'s block diagonal)
/// and the eigenvalues `(alphar[j] + i alphai[j]) / beta[j]`.
pub fn ggesReal(allocator: Allocator, io: Io, comptime T: type, jobvsl: VectorJob, jobvsr: VectorJob, n: usize, a: []T, lda: usize, b: []T, ldb: usize, alphar: []T, alphai: []T, beta: []T, vsl: []T, ldvsl: usize, vsr: []T, ldvsr: usize) GgevError!void {
    comptime ops.requireReal(T, "ggesReal");
    if (alphar.len < n or alphai.len < n or beta.len < n) return error.SliceTooShort;
    const left = jobvsl == .compute;
    const right = jobvsr == .compute;
    try check(T, n, a, lda, b, ldb, left, vsl, ldvsl, right, vsr, ldvsr);
    if (n == 0) return;
    try ggev_complex.reduce(allocator, io, T, left, right, n, a, lda, b, ldb, vsl, ldvsl, vsr, ldvsr);
    qz_real.hgeqzReal(io, T, .schur_form, if (left) .accumulate else .none, if (right) .accumulate else .none, n, 0, n, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr) catch |e| return ops.narrow(GgevError, e);
}

fn check(comptime T: type, n: usize, a: []T, lda: usize, b: []T, ldb: usize, left: bool, vsl: []T, ldvsl: usize, right: bool, vsr: []T, ldvsr: usize) GgevError!void {
    try util.holds(a.len, n, n, lda);
    try util.holds(b.len, n, n, ldb);
    if (left) try util.holds(vsl.len, n, n, ldvsl);
    if (right) try util.holds(vsr.len, n, n, ldvsr);
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const types = @import("../blas/types.zig");
const blas_check = @import("../blas/_test/check.zig");
const pencil = @import("_test/pencil.zig");

fn ggesCase(comptime T: type, n: usize, seed: u64) !void {
    const al = testing.allocator;
    const complex = comptime types.isComplex(T);
    const F = types.Float(T);
    const a0 = try al.alloc(T, n * n);
    defer al.free(a0);
    const b0 = try al.alloc(T, n * n);
    defer al.free(b0);
    blas_check.fill(T, a0, seed);
    blas_check.fill(T, b0, seed + 1);
    const a = try al.dupe(T, a0);
    defer al.free(a);
    const b = try al.dupe(T, b0);
    defer al.free(b);
    const q = try al.alloc(T, n * n);
    defer al.free(q);
    const z = try al.alloc(T, n * n);
    defer al.free(z);
    const w = try al.alloc(T, 3 * n);
    defer al.free(w);
    const ld = @max(n, 1);
    if (complex) {
        try ggesComplex(al, blas_check.io, T, .compute, .compute, n, a, ld, b, ld, w[0..n], w[n .. 2 * n], q, ld, z, ld);
    } else try ggesReal(al, blas_check.io, T, .compute, .compute, n, a, ld, b, ld, w[0..n], w[n .. 2 * n], w[2 * n ..], q, ld, z, ld);
    const tol = 60 * @as(F, @floatFromInt(n * n)) * std.math.floatEps(F);
    for (0..n) |j| for (j + 1..n) |i| {
        try testing.expect(ops.isZero(T, b[i + j * n]));
        // `S` triangular, or below a real pair's 2 x 2 block
        const pair_row = !complex and i == j + 1 and w[n + j] != 0;
        if (!pair_row) try testing.expect(ops.isZero(T, a[i + j * n]));
    };
    try testing.expect(pencil.unitary(T, n, q) <= tol);
    try testing.expect(pencil.unitary(T, n, z) <= tol);
    try testing.expect(pencil.residual(T, n, q, a, z, a0) <= tol * 4);
    try testing.expect(pencil.residual(T, n, q, b, z, b0) <= tol * 4);
}

test "ggesComplex, ggesReal: (A, B) = (Q S Zᴴ, Q T Zᴴ), the forms as promised" {
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]usize{ 0, 1, 2, 5, 13 }) |n| try ggesCase(T, n, n * 9 + 2);
}
