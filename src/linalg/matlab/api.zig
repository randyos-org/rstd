//! Friendly, MATLAB-command-named bindings on top of the low-level
//! Fortran-derived kernels in this directory - e.g. `eig` here is the
//! dispatch logic the original REPL's `EIG` command used to pick between
//! the Hermitian (`htridi`+`imtql2`+`htribk`) and general
//! (`corth`+`comqr3`) eigensolvers, not a new algorithm. Pure aliases
//! (`lu`, `svd`, `qr`, `hilb`) just repoint a REPL name at the underlying
//! kernel; see that kernel's doc comment for the actual contract.
//!
//! Functions here that need scratch space sized at runtime (`eig`,
//! `cond`, `rank`, `pinv`) take a `std.mem.Allocator` and return
//! caller-owned slices; everything else works directly on caller-supplied
//! buffers like the low-level kernels do, so this module stays usable
//! without an allocator where the original REPL command didn't need one
//! either.

const std = @import("std");
const blas1 = @import("blas1.zig");
const lu_mod = @import("lu.zig");
const cholesky = @import("cholesky.zig");
const generators = @import("generators.zig");
const tridiagonal = @import("tridiagonal.zig");
const eigen_symmetric = @import("eigen_symmetric.zig");
const hessenberg = @import("hessenberg.zig");
const eigen_general = @import("eigen_general.zig");
const svd_mod = @import("svd.zig");
const qr_mod = @import("qr.zig");

const Cx = blas1.Cx;
const ix = blas1.ix;

/// LU: factor `A` in place (packed LINPACK/LAPACK-style: unit multipliers
/// below the diagonal, `U` on and above it) - REPL name for `wgefa`.
pub const lu = lu_mod.wgefa;

/// HILB: REPL name for `hilber` (the closed-form inverse Hilbert matrix).
pub const hilb = generators.hilber;

/// SVD: REPL name for `wsvdc`.
pub const svd = svd_mod.wsvdc;

/// QR: REPL name for `wqrdc`. Pair with `qrsl` (== `wqrsl`) to extract `Q`,
/// solve least squares, etc. - see `qr.zig`'s tests for a worked example.
pub const qr = qr_mod.wqrdc;
pub const qrsl = qr_mod.wqrsl;

/// RREF: REPL name for `generators.rref`.
pub const rref = generators.rref;

/// DET: determinant of the `n`x`n` matrix `ar`/`ai` (column-major, leading
/// dimension `lda`), via LU factorization. `A` is destroyed; `ipvt`
/// (length `n`) is scratch.
pub fn det(ar: []f64, ai: []f64, lda: usize, n: usize, ipvt: []usize) lu_mod.Determinant {
    _ = lu_mod.wgefa(ar, ai, lda, n, ipvt);
    var result: lu_mod.Determinant = undefined;
    var no_workspace: [0]f64 = .{};
    lu_mod.wgedi(ar, ai, lda, n, ipvt, &no_workspace, &no_workspace, 10, &result);
    return result;
}

/// INV: invert the `n`x`n` matrix `ar`/`ai` in place, via LU factorization
/// with condition estimation. `ipvt` (length `n`) and `zr`/`zi` (length
/// `n`) are scratch. Returns the reciprocal condition number (`wgeco`'s
/// `rcond`) so callers can judge accuracy - a value so small that
/// `1.0 + rcond == 1.0` means `A` is numerically singular and the
/// "inverse" is not to be trusted.
pub fn inv(ar: []f64, ai: []f64, lda: usize, n: usize, ipvt: []usize, zr: []f64, zi: []f64) f64 {
    const rcond = lu_mod.wgeco(ar, ai, lda, n, ipvt, zr, zi);
    // wgeco's own use of zr/zi (an approximate null vector) is done by the
    // time we get here, so it doubles as wgedi's scratch workspace.
    lu_mod.wgedi(ar, ai, lda, n, ipvt, zr, zi, 1, null);
    return rcond;
}

/// CHOL: Cholesky-factor the Hermitian positive-definite matrix `ar`/`ai`
/// in place into upper-triangular `R` (`ctrans(R)*R == A`), zeroing the
/// strict lower triangle so the result is a clean triangular matrix
/// (unlike raw `wpofa`, which leaves it untouched). Returns 0 on success,
/// or the 1-based order of the leading minor that failed to be positive
/// definite (the matrix is left only partially factored in that case).
pub fn chol(ar: []f64, ai: []f64, lda: usize, n: usize) usize {
    const info = cholesky.wpofa(ar, ai, lda, n);
    if (info != 0) return info;
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        blas1.wset(n - j, Cx.init(0.0, 0.0), ar[ix(lda, j + 1, j)..], ai[ix(lda, j + 1, j)..], 1);
    }
    return 0;
}

/// Eigenvalues, and optionally eigenvectors, of an `n`x`n` complex matrix.
/// All slices are owned by the caller and must be freed with the same
/// allocator used to create them (via `deinit`).
pub const EigResult = struct {
    values_r: []f64,
    values_i: []f64,
    /// Column-major `n`x`n`, present iff `eig` was asked to compute
    /// vectors. For a non-Hermitian `A` these are unnormalized, matching
    /// `comqr3`'s (`job == 2`) EISPACK COMQR2 behavior.
    vectors_r: ?[]f64,
    vectors_i: ?[]f64,

    pub fn deinit(self: EigResult, allocator: std.mem.Allocator) void {
        allocator.free(self.values_r);
        allocator.free(self.values_i);
        if (self.vectors_r) |v| allocator.free(v);
        if (self.vectors_i) |v| allocator.free(v);
    }
};

pub const EigError = error{EigenvalueDidNotConverge} || std.mem.Allocator.Error;

fn isHermitian(ar: []const f64, ai: []const f64, n: usize) bool {
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        var i: usize = 1;
        while (i <= j) : (i += 1) {
            if (ar[ix(n, i, j)] != ar[ix(n, j, i)] or ai[ix(n, i, j)] != -ai[ix(n, j, i)]) return false;
        }
    }
    return true;
}

/// EIG: eigenvalues (and, if `want_vectors`, eigenvectors) of the `n`x`n`
/// complex matrix `ar`/`ai` (column-major, leading dimension `n` - not
/// mutated). Automatically dispatches to the Hermitian eigensolver
/// (`htridi`+`imtql2`+`htribk`, real eigenvalues) or the general one
/// (`corth`+`comqr3`), exactly as the original REPL's `EIG` command did.
pub fn eig(allocator: std.mem.Allocator, ar: []const f64, ai: []const f64, n: usize, want_vectors: bool) EigError!EigResult {
    const work_r = try allocator.alloc(f64, n * n);
    defer allocator.free(work_r);
    const work_i = try allocator.alloc(f64, n * n);
    defer allocator.free(work_i);
    @memcpy(work_r, ar[0 .. n * n]);
    @memcpy(work_i, ai[0 .. n * n]);

    const values_r = try allocator.alloc(f64, n);
    errdefer allocator.free(values_r);
    const values_i = try allocator.alloc(f64, n);
    errdefer allocator.free(values_i);

    if (isHermitian(ar, ai, n)) {
        @memset(values_i, 0.0);
        const e = try allocator.alloc(f64, n);
        defer allocator.free(e);
        const e2 = try allocator.alloc(f64, n);
        defer allocator.free(e2);
        const taur = try allocator.alloc(f64, n);
        defer allocator.free(taur);
        const taui = try allocator.alloc(f64, n);
        defer allocator.free(taui);

        tridiagonal.htridi(n, n, work_r, work_i, values_r, e, e2, taur, taui);

        if (!want_vectors) {
            var no_vectors: [0]f64 = .{};
            const ierr = eigen_symmetric.imtql2(n, values_r, e, &no_vectors, n, 0);
            if (ierr != 0) return error.EigenvalueDidNotConverge;
            return .{ .values_r = values_r, .values_i = values_i, .vectors_r = null, .vectors_i = null };
        }

        const zr = try allocator.alloc(f64, n * n);
        errdefer allocator.free(zr);
        const zi = try allocator.alloc(f64, n * n);
        errdefer allocator.free(zi);
        @memset(zr, 0.0);
        @memset(zi, 0.0);
        var k: usize = 1;
        while (k <= n) : (k += 1) zr[ix(n, k, k)] = 1.0;

        const ierr = eigen_symmetric.imtql2(n, values_r, e, zr, n, 1);
        if (ierr != 0) return error.EigenvalueDidNotConverge;
        tridiagonal.htribk(n, n, work_r, work_i, taur, taui, n, zr, zi);
        return .{ .values_r = values_r, .values_i = values_i, .vectors_r = zr, .vectors_i = zi };
    }

    const ortr = try allocator.alloc(f64, n);
    defer allocator.free(ortr);
    const orti = try allocator.alloc(f64, n);
    defer allocator.free(orti);
    @memset(ortr, 0.0);
    @memset(orti, 0.0);

    hessenberg.corth(n, n, 1, n, work_r, work_i, ortr, orti);

    const zr = try allocator.alloc(f64, n * n);
    errdefer allocator.free(zr);
    const zi = try allocator.alloc(f64, n * n);
    errdefer allocator.free(zi);
    @memset(zr, 0.0);
    @memset(zi, 0.0);

    const job: usize = if (want_vectors) 2 else 0;
    const ierr = eigen_general.comqr3(n, n, 1, n, ortr, orti, work_r, work_i, values_r, values_i, zr, zi, job);
    if (ierr != 0) {
        allocator.free(zr);
        allocator.free(zi);
        return error.EigenvalueDidNotConverge;
    }
    if (!want_vectors) {
        allocator.free(zr);
        allocator.free(zi);
        return .{ .values_r = values_r, .values_i = values_i, .vectors_r = null, .vectors_i = null };
    }
    return .{ .values_r = values_r, .values_i = values_i, .vectors_r = zr, .vectors_i = zi };
}

fn singularValues(allocator: std.mem.Allocator, ar: []const f64, ai: []const f64, n: usize, p: usize) ![]f64 {
    const work_r = try allocator.alloc(f64, n * p);
    defer allocator.free(work_r);
    const work_i = try allocator.alloc(f64, n * p);
    defer allocator.free(work_i);
    @memcpy(work_r, ar[0 .. n * p]);
    @memcpy(work_i, ai[0 .. n * p]);

    const mm = @min(n + 1, p);
    const sr = try allocator.alloc(f64, mm);
    errdefer allocator.free(sr);
    const si = try allocator.alloc(f64, mm);
    defer allocator.free(si);
    const er = try allocator.alloc(f64, p);
    defer allocator.free(er);
    const ei = try allocator.alloc(f64, p);
    defer allocator.free(ei);
    const workr = try allocator.alloc(f64, n);
    defer allocator.free(workr);
    const worki = try allocator.alloc(f64, n);
    defer allocator.free(worki);
    var dummy: [1]f64 = .{0.0};

    _ = svd_mod.wsvdc(work_r, work_i, n, n, p, sr, si, er, ei, &dummy, &dummy, 1, &dummy, &dummy, 1, workr, worki, 0);
    return sr;
}

/// COND: ratio of the largest to smallest singular value of the `n`x`p`
/// matrix `ar`/`ai` (2-norm condition number). `A` is not mutated.
pub fn cond(allocator: std.mem.Allocator, ar: []const f64, ai: []const f64, n: usize, p: usize) !f64 {
    const s = try singularValues(allocator, ar, ai, n, p);
    defer allocator.free(s);
    const smallest = s[@min(n, p) - 1];
    if (smallest == 0.0) return std.math.inf(f64);
    return s[0] / smallest;
}

/// RANK: number of singular values of the `n`x`p` matrix `ar`/`ai` that
/// exceed `tol` (pass `null` to use the same default the original used:
/// `max(n, p) * eps * largest_singular_value`). `A` is not mutated.
pub fn rank(allocator: std.mem.Allocator, ar: []const f64, ai: []const f64, n: usize, p: usize, tol: ?f64) !usize {
    const s = try singularValues(allocator, ar, ai, n, p);
    defer allocator.free(s);
    const mn = @min(n, p);
    const threshold = tol orelse @as(f64, @floatFromInt(@max(n, p))) * std.math.floatEps(f64) * s[0];
    var k: usize = 0;
    for (s[0..mn]) |sv| {
        if (sv <= threshold) break;
        k += 1;
    }
    return k;
}

/// PINV: Moore-Penrose pseudoinverse of the `n`x`p` matrix `ar`/`ai`, via
/// SVD (`pinv(A) = V * pinv(S) * ctrans(U)`, treating singular values at or
/// below `tol` as zero - `null` uses the same default as `rank`). Writes
/// the `p`x`n` result into caller-supplied `result_r`/`result_i`
/// (column-major, leading dimension `p`). `A` is not mutated.
pub fn pinv(
    allocator: std.mem.Allocator,
    ar: []const f64,
    ai: []const f64,
    n: usize,
    p: usize,
    tol: ?f64,
    result_r: []f64,
    result_i: []f64,
) !void {
    const work_r = try allocator.alloc(f64, n * p);
    defer allocator.free(work_r);
    const work_i = try allocator.alloc(f64, n * p);
    defer allocator.free(work_i);
    @memcpy(work_r, ar[0 .. n * p]);
    @memcpy(work_i, ai[0 .. n * p]);

    const mn = @min(n, p);
    const mm = @min(n + 1, p);
    const sr = try allocator.alloc(f64, mm);
    defer allocator.free(sr);
    const si = try allocator.alloc(f64, mm);
    defer allocator.free(si);
    const er = try allocator.alloc(f64, p);
    defer allocator.free(er);
    const ei = try allocator.alloc(f64, p);
    defer allocator.free(ei);
    const ur = try allocator.alloc(f64, n * n);
    defer allocator.free(ur);
    const ui = try allocator.alloc(f64, n * n);
    defer allocator.free(ui);
    const vr = try allocator.alloc(f64, p * p);
    defer allocator.free(vr);
    const vi = try allocator.alloc(f64, p * p);
    defer allocator.free(vi);
    const workr = try allocator.alloc(f64, n);
    defer allocator.free(workr);
    const worki = try allocator.alloc(f64, n);
    defer allocator.free(worki);

    _ = svd_mod.wsvdc(work_r, work_i, n, n, p, sr, si, er, ei, ur, ui, n, vr, vi, p, workr, worki, 11);

    const threshold = tol orelse @as(f64, @floatFromInt(@max(n, p))) * std.math.floatEps(f64) * sr[0];
    var k: usize = 0;
    for (sr[0..mn]) |sv| {
        if (sv <= threshold) break;
        k += 1;
    }

    var i: usize = 1;
    while (i <= p) : (i += 1) {
        var j: usize = 1;
        while (j <= n) : (j += 1) {
            var re: f64 = 0.0;
            var im: f64 = 0.0;
            var l: usize = 1;
            while (l <= k) : (l += 1) {
                const v = Cx.init(vr[ix(p, i, l)], vi[ix(p, i, l)]);
                const u_conj = Cx.init(ur[ix(n, j, l)], -ui[ix(n, j, l)]);
                const term = blas1.wmul(v, u_conj);
                re += term.re / sr[l - 1];
                im += term.im / sr[l - 1];
            }
            result_r[ix(p, i, j)] = re;
            result_i[ix(p, i, j)] = im;
        }
    }
}

const testing = std.testing;

test "det, inv and chol agree on a small real matrix" {
    var ar = [_]f64{ 4.0, 2.0, 2.0, 3.0 }; // [[4,2],[2,3]], det=8, SPD
    var ai: [4]f64 = @splat(0.0);
    var ipvt: [2]usize = @splat(0);

    var ar_for_det = ar;
    var ai_for_det = ai;
    const d = det(&ar_for_det, &ai_for_det, 2, 2, &ipvt);
    try testing.expectApproxEqAbs(8.0, d.mantissa.re * std.math.pow(f64, 10.0, d.exponent10), 1e-9);

    var ar_for_inv = ar;
    var ai_for_inv = ai;
    var zr: [2]f64 = @splat(0.0);
    var zi: [2]f64 = @splat(0.0);
    const rcond = inv(&ar_for_inv, &ai_for_inv, 2, 2, &ipvt, &zr, &zi);
    try testing.expect(rcond > 0.0);
    // inv([[4,2],[2,3]]) = [[3,-2],[-2,4]] / 8.
    try testing.expectApproxEqAbs(3.0 / 8.0, ar_for_inv[ix(2, 1, 1)], 1e-9);
    try testing.expectApproxEqAbs(-2.0 / 8.0, ar_for_inv[ix(2, 1, 2)], 1e-9);
    try testing.expectApproxEqAbs(4.0 / 8.0, ar_for_inv[ix(2, 2, 2)], 1e-9);

    const info = chol(&ar, &ai, 2, 2);
    try testing.expectEqual(@as(usize, 0), info);
    try testing.expectApproxEqAbs(0.0, ar[ix(2, 2, 1)], 1e-12); // lower triangle zeroed
}

test "eig dispatches to the Hermitian path for a symmetric matrix" {
    // A = [[2,1],[1,2]]: eigenvalues 1 and 3.
    const ar = [_]f64{ 2.0, 1.0, 1.0, 2.0 };
    const ai: [4]f64 = @splat(0.0);
    const result = try eig(testing.allocator, &ar, &ai, 2, true);
    defer result.deinit(testing.allocator);

    const lo = @min(result.values_r[0], result.values_r[1]);
    const hi = @max(result.values_r[0], result.values_r[1]);
    try testing.expectApproxEqAbs(1.0, lo, 1e-9);
    try testing.expectApproxEqAbs(3.0, hi, 1e-9);
    try testing.expect(result.vectors_r != null);
}

test "eig dispatches to the general path for a non-symmetric matrix" {
    const ar = [_]f64{ 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 };
    const ai: [9]f64 = @splat(0.0);
    const result = try eig(testing.allocator, &ar, &ai, 3, false);
    defer result.deinit(testing.allocator);

    try testing.expect(result.vectors_r == null);
    var trace: f64 = 0.0;
    for (result.values_r) |v| trace += v;
    try testing.expectApproxEqAbs(16.0, trace, 1e-8);
}

test "cond, rank and pinv agree for a known rank-deficient matrix" {
    // A = [[1,2],[2,4]] (rank 1, columns parallel).
    const ar = [_]f64{ 1.0, 2.0, 2.0, 4.0 };
    const ai: [4]f64 = @splat(0.0);

    const c = try cond(testing.allocator, &ar, &ai, 2, 2);
    try testing.expect(c > 1e10 or std.math.isInf(c));

    const r = try rank(testing.allocator, &ar, &ai, 2, 2, null);
    try testing.expectEqual(@as(usize, 1), r);

    var pr: [4]f64 = @splat(0.0);
    var pi: [4]f64 = @splat(0.0);
    try pinv(testing.allocator, &ar, &ai, 2, 2, null, &pr, &pi);
    // A * pinv(A) * A should reconstruct A (pseudoinverse identity).
    var apa_full: [4]f64 = @splat(0.0);
    for (0..2) |row| {
        for (0..2) |col| {
            var s: f64 = 0.0;
            for (0..2) |m| {
                var ap: f64 = 0.0;
                for (0..2) |l| ap += ar[ix(2, row + 1, l + 1)] * pr[ix(2, l + 1, m + 1)];
                s += ap * ar[ix(2, m + 1, col + 1)];
            }
            apa_full[ix(2, row + 1, col + 1)] = s;
        }
    }
    for (0..4) |idx| try testing.expectApproxEqAbs(ar[idx], apa_full[idx], 1e-6);
}
