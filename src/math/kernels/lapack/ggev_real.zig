//! `?GGEV` for real pencils: faer v0.24.4's `compute_gevd_generic` for real
//! types (`faer/src/linalg/gevd/mod.rs`, MIT): the steps of
//! `ggev_complex.zig`'s `reduce`, the real QZ, then `ggev_real_vectors.zig`.

const std = @import("std");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const ggev_complex = @import("ggev_complex.zig");
const qz_real = @import("qz_real.zig");
const vectors = @import("ggev_real_vectors.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;

pub const VectorJob = ggev_complex.VectorJob;
pub const GgevError = ggev_complex.GgevError;

/// `?GGEV`, real: the generalized eigenvalues `(alphar[j] + i alphai[j]) /
/// beta[j]` of `(A, B)` (`n x n`; both overwritten), a complex pair
/// together with the positive imaginary part first, and as `jobvl` and
/// `jobvr` ask the left (`yᴴ A = λ yᴴ B`) and right (`A x = λ B x`)
/// eigenvectors: a real eigenvalue's in its column, a pair's `v(:, j) ± i
/// v(:, j + 1)` in two, each of unit 2-norm. `beta[j] = 0` is infinite.
pub fn ggevReal(allocator: Allocator, io: Io, comptime T: type, jobvl: VectorJob, jobvr: VectorJob, n: usize, a: []T, lda: usize, b: []T, ldb: usize, alphar: []T, alphai: []T, beta: []T, vl: []T, ldvl: usize, vr: []T, ldvr: usize) GgevError!void {
    comptime ops.requireReal(T, "ggevReal");
    try util.holds(a.len, n, n, lda);
    try util.holds(b.len, n, n, ldb);
    if (alphar.len < n or alphai.len < n or beta.len < n) return error.SliceTooShort;
    const left = jobvl == .compute;
    const right = jobvr == .compute;
    if (left) try util.holds(vl.len, n, n, ldvl);
    if (right) try util.holds(vr.len, n, n, ldvr);
    if (n == 0) return;
    try ggev_complex.reduce(allocator, io, T, left, right, n, a, lda, b, ldb, vl, ldvl, vr, ldvr);
    const job: qz_real.SchurJob = if (left or right) .schur_form else .eigenvalues;
    qz_real.hgeqzReal(io, T, job, if (left) .accumulate else .none, if (right) .accumulate else .none, n, 0, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr) catch |e| return ops.narrow(GgevError, e);
    if (left or right) try vectors.eigenvectors(allocator, T, n, a, lda, b, ldb, if (left) vl else null, ldvl, if (right) vr else null, ldvr);
}
