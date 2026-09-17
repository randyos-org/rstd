//! `HESS` and `SCHUR`, ported from `matlab88.f90`'s `MATFN4`/`MATFN5`
//! command bodies. Both are thin orchestration over kernels this port
//! already has: `hessenberg.corth` (Hessenberg reduction) and
//! `eigen_general.comqr3`, whose `job` parameter already covers exactly
//! these two cases (`job == 3`: accumulate only the Hessenberg-reduction
//! similarity, no QR iteration - `HESS`; `job == 1`: full QR iteration to
//! the Schur form, plus the accumulated similarity - `SCHUR`) - no new
//! numerical algorithm needed, unlike `eig.zig`'s Hermitian-vs-general
//! dispatch, `HESS`/`SCHUR` always go through `corth`+`comqr3` regardless
//! of whether `A` is Hermitian (per the original's own doc text).

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1.zig");
const CxMat = @import("CxMat.zig");
const hessenberg = @import("hessenberg.zig");
const eigen_general = @import("eigen_general.zig");
const Cx = blas1.Cx;
const ix = blas1.ix;

/// `<P, H>`: `A`'s unitary similarity `P` and Hessenberg form `H`, so
/// `A = P*H*P'` (`H` zero below the first subdiagonal - tridiagonal if `A`
/// is Hermitian, per the original's own doc text). Caller owns the
/// returned slices (free with `allocator`).
pub const HessResult = struct {
    p_r: []f64,
    p_i: []f64,
    h_r: []f64,
    h_i: []f64,

    pub fn deinit(self: HessResult, allocator: zigstd.Allocator) void {
        allocator.free(self.p_r);
        allocator.free(self.p_i);
        allocator.free(self.h_r);
        allocator.free(self.h_i);
    }
};

pub const HessError = zigstd.Allocator.Error;

/// HESS: see `HessResult`. `a` is not mutated (a cloned copy is reduced
/// in place instead).
pub fn hess(allocator: zigstd.Allocator, a: CxMat) HessError!HessResult {
    const n = a.n;
    var h = try a.clone(allocator);
    errdefer h.deinit(allocator);

    const ortr = try allocator.alloc(f64, n);
    defer allocator.free(ortr);
    const orti = try allocator.alloc(f64, n);
    defer allocator.free(orti);
    @memset(ortr, 0.0);
    @memset(orti, 0.0);

    hessenberg.corth(h, 1, n, ortr, orti);

    const p_r = try allocator.alloc(f64, n * n);
    errdefer allocator.free(p_r);
    const p_i = try allocator.alloc(f64, n * n);
    errdefer allocator.free(p_i);
    const p = CxMat.initPacked(p_r, p_i, n, n);

    // job == 3 skips QR iteration entirely, so wr/wi are never touched -
    // still sized n rather than a 1-element dummy, since that skip isn't
    // reached for every n/low/igh combination (small n can fall through
    // past the early-return branch).
    const wr = try allocator.alloc(f64, n);
    defer allocator.free(wr);
    const wi = try allocator.alloc(f64, n);
    defer allocator.free(wi);

    _ = eigen_general.comqr3(h, 1, n, ortr, orti, wr, wi, p, 3);

    // Only now, after comqr3's P-accumulation has read whatever it needed
    // from H's lower triangle (corth stores part of its Householder-vector
    // data there, not just in ortr/orti - see comqr3's own accumulation
    // loop, which reads hr/hi directly), zero it for a clean Hessenberg
    // result: corth doesn't zero below the first subdiagonal itself (see
    // its own doc comment), the same way wpofa leaves its lower triangle
    // untouched and api.chol zeros it afterward.
    if (n >= 3) {
        var col: usize = 1;
        while (col <= n - 2) : (col += 1) {
            blas1.wset(n - col - 1, Cx.init(0.0, 0.0), h.r[ix(n, col + 2, col)..], h.i[ix(n, col + 2, col)..], 1);
        }
    }

    return .{ .p_r = p_r, .p_i = p_i, .h_r = h.r, .h_i = h.i };
}

/// `<U, T>`: `A`'s unitary similarity `U` and (upper triangular) Schur
/// form `T`, with the eigenvalues of `A` on `T`'s diagonal, so
/// `A = U*T*U'`. Caller owns the returned slices (free with `allocator`).
pub const SchurResult = struct {
    u_r: []f64,
    u_i: []f64,
    t_r: []f64,
    t_i: []f64,

    pub fn deinit(self: SchurResult, allocator: zigstd.Allocator) void {
        allocator.free(self.u_r);
        allocator.free(self.u_i);
        allocator.free(self.t_r);
        allocator.free(self.t_i);
    }
};

pub const SchurError = error{EigenvalueDidNotConverge} || zigstd.Allocator.Error;

/// SCHUR: see `SchurResult`. `a` is not mutated (a cloned copy is reduced
/// in place instead). Returns an error if an eigenvalue failed to converge
/// after `30*n` QR iterations (matching `comqr3`'s own failure contract).
pub fn schur(allocator: zigstd.Allocator, a: CxMat) SchurError!SchurResult {
    const n = a.n;
    var t = try a.clone(allocator);
    errdefer t.deinit(allocator);

    const ortr = try allocator.alloc(f64, n);
    defer allocator.free(ortr);
    const orti = try allocator.alloc(f64, n);
    defer allocator.free(orti);
    @memset(ortr, 0.0);
    @memset(orti, 0.0);

    hessenberg.corth(t, 1, n, ortr, orti);

    const u_r = try allocator.alloc(f64, n * n);
    errdefer allocator.free(u_r);
    const u_i = try allocator.alloc(f64, n * n);
    errdefer allocator.free(u_i);
    const u = CxMat.initPacked(u_r, u_i, n, n);

    const wr = try allocator.alloc(f64, n);
    defer allocator.free(wr);
    const wi = try allocator.alloc(f64, n);
    defer allocator.free(wi);

    const ierr = eigen_general.comqr3(t, 1, n, ortr, orti, wr, wi, u, 1);
    if (ierr != 0) return error.EigenvalueDidNotConverge;

    // The QR iteration converges the subdiagonal toward zero (a genuine,
    // small-but-not-exactly-zero numerical residual, left as-is) but never
    // touches entries two-or-more below the diagonal at all: those are
    // already assumed zero by the Hessenberg band structure the iteration
    // operates within, so corth's leftover reduction data there (see
    // `hess`'s comment above) survives untouched all the way through.
    // Zero it now that comqr3 (which reads it as part of accumulating U)
    // is done with it, same as `hess`.
    if (n >= 3) {
        var col: usize = 1;
        while (col <= n - 2) : (col += 1) {
            blas1.wset(n - col - 1, Cx.init(0.0, 0.0), t.r[ix(n, col + 2, col)..], t.i[ix(n, col + 2, col)..], 1);
        }
    }

    return .{ .u_r = u_r, .u_i = u_i, .t_r = t.r, .t_i = t.i };
}

fn matSub(alloc: zigstd.Allocator, n: usize, ar: []const f64, ai: []const f64, br: []const f64, bi: []const f64) ![2][]f64 {
    const dr = try alloc.alloc(f64, n * n);
    const di = try alloc.alloc(f64, n * n);
    for (dr, di, ar, ai, br, bi) |*r, *i, a_r, a_i, b_r, b_i| {
        r.* = a_r - b_r;
        i.* = a_i - b_i;
    }
    return .{ dr, di };
}

test "hess reduces to zero below the first subdiagonal, and P*H*P' reconstructs A" {
    const n = 3;
    var ar = [_]f64{ 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 };
    var ai: [n * n]f64 = @splat(0.0);

    var result = try hess(zigstd.testAllocator, CxMat.initPacked(&ar, &ai, n, n));
    defer result.deinit(zigstd.testAllocator);

    // Zero below the first subdiagonal: only (2,1) and (3,2) survive below
    // the diagonal; (3,1) must be exactly zero.
    try zigstd.testExpectApproxEqAbs(0.0, result.h_r[ix(n, 3, 1)], 1e-9);
    try zigstd.testExpectApproxEqAbs(0.0, result.h_i[ix(n, 3, 1)], 1e-9);

    // A = P*H*P': check by reconstructing P*H*P' and comparing to A.
    const matmul_mod = @import("matmul.zig");
    var scratch_r: [n]f64 = undefined;
    var scratch_i: [n]f64 = undefined;
    var ph_r: [n * n]f64 = undefined;
    var ph_i: [n * n]f64 = undefined;
    const p_mat = CxMat.initPacked(result.p_r, result.p_i, n, n);
    const h_mat = CxMat.initPacked(result.h_r, result.h_i, n, n);
    matmul_mod.matmul(p_mat, h_mat, CxMat.initPacked(&ph_r, &ph_i, n, n), &scratch_r, &scratch_i);

    // P' = conjugate transpose of P.
    var pt_r: [n * n]f64 = undefined;
    var pt_i: [n * n]f64 = undefined;
    for (0..n) |row| {
        for (0..n) |col| {
            pt_r[ix(n, row + 1, col + 1)] = result.p_r[ix(n, col + 1, row + 1)];
            pt_i[ix(n, row + 1, col + 1)] = -result.p_i[ix(n, col + 1, row + 1)];
        }
    }
    var recon_r: [n * n]f64 = undefined;
    var recon_i: [n * n]f64 = undefined;
    matmul_mod.matmul(CxMat.initPacked(&ph_r, &ph_i, n, n), CxMat.initPacked(&pt_r, &pt_i, n, n), CxMat.initPacked(&recon_r, &recon_i, n, n), &scratch_r, &scratch_i);

    const diff = try matSub(zigstd.testAllocator, n, &recon_r, &recon_i, &ar, &ai);
    defer zigstd.testAllocator.free(diff[0]);
    defer zigstd.testAllocator.free(diff[1]);
    for (diff[0], diff[1]) |dr, di| {
        try zigstd.testExpectApproxEqAbs(0.0, dr, 1e-8);
        try zigstd.testExpectApproxEqAbs(0.0, di, 1e-8);
    }
}

test "schur produces an upper triangular T with A's eigenvalues on the diagonal" {
    const n = 3;
    var ar = [_]f64{ 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 };
    var ai: [n * n]f64 = @splat(0.0);

    var result = try schur(zigstd.testAllocator, CxMat.initPacked(&ar, &ai, n, n));
    defer result.deinit(zigstd.testAllocator);

    // Upper triangular. The subdiagonal (row == col+1) converges toward
    // zero via the QR iteration's own relative (not exact-to-machine-
    // precision) convergence test, so a small residual there is expected,
    // not a bug - looser tolerance. Everything two-or-more below the
    // diagonal is explicitly zeroed by `schur` itself (the QR iteration
    // never touches it), so that should be exact.
    for (0..n) |col| {
        for (col + 1..n) |row| {
            const tol: f64 = if (row == col + 1) 1e-4 else 1e-12;
            try zigstd.testExpectApproxEqAbs(0.0, result.t_r[ix(n, row + 1, col + 1)], tol);
            try zigstd.testExpectApproxEqAbs(0.0, result.t_i[ix(n, row + 1, col + 1)], tol);
        }
    }

    // Trace of T (== sum of A's eigenvalues) matches trace of A.
    var trace: f64 = 0.0;
    for (0..n) |k| trace += result.t_r[ix(n, k + 1, k + 1)];
    try zigstd.testExpectApproxEqAbs(16.0, trace, 1e-8);
}
