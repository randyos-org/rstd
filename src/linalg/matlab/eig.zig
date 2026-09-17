//! `EIG`'s dispatch logic, kept as its own file (rather than folded into
//! `CxMat.zig`, whose `eig` method delegates here) to keep files under
//! this project's line-count guideline: pick the Hermitian eigensolver
//! (`tridiagonal.htridi` + `eigen_symmetric.imtql2` + `tridiagonal.htribk`)
//! or the general one (`hessenberg.corth` + `eigen_general.comqr3`),
//! exactly as the original REPL's `EIG` command did.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1.zig");
const CxMat = @import("CxMat.zig");
const tridiagonal = @import("tridiagonal.zig");
const eigen_symmetric = @import("eigen_symmetric.zig");
const hessenberg = @import("hessenberg.zig");
const eigen_general = @import("eigen_general.zig");

const ix = blas1.ix;

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

    pub fn deinit(self: EigResult, allocator: zigstd.Allocator) void {
        allocator.free(self.values_r);
        allocator.free(self.values_i);
        if (self.vectors_r) |v| allocator.free(v);
        if (self.vectors_i) |v| allocator.free(v);
    }
};

pub const EigError = error{EigenvalueDidNotConverge} || zigstd.Allocator.Error;

fn isHermitian(a: CxMat) bool {
    var j: usize = 1;
    while (j <= a.n) : (j += 1) {
        var i: usize = 1;
        while (i <= j) : (i += 1) {
            if (a.r[a.idx(i, j)] != a.r[a.idx(j, i)] or a.i[a.idx(i, j)] != -a.i[a.idx(j, i)]) return false;
        }
    }
    return true;
}

/// EIG: eigenvalues (and, if `want_vectors`, eigenvectors) of the complex
/// matrix `a` (not mutated - a cloned copy is reduced in place instead).
/// Automatically dispatches to the Hermitian eigensolver
/// (`htridi`+`imtql2`+`htribk`, real eigenvalues) or the general one
/// (`corth`+`comqr3`), exactly as the original REPL's `EIG` command did.
pub fn eig(allocator: zigstd.Allocator, a: CxMat, want_vectors: bool) EigError!EigResult {
    const n = a.n;
    var work = try a.clone(allocator);
    defer work.deinit(allocator);

    const values_r = try allocator.alloc(f64, n);
    errdefer allocator.free(values_r);
    const values_i = try allocator.alloc(f64, n);
    errdefer allocator.free(values_i);

    if (isHermitian(a)) {
        @memset(values_i, 0.0);
        const e = try allocator.alloc(f64, n);
        defer allocator.free(e);
        const e2 = try allocator.alloc(f64, n);
        defer allocator.free(e2);
        const taur = try allocator.alloc(f64, n);
        defer allocator.free(taur);
        const taui = try allocator.alloc(f64, n);
        defer allocator.free(taui);

        tridiagonal.htridi(work, values_r, e, e2, taur, taui);

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
        tridiagonal.htribk(work, taur, taui, CxMat.initPacked(zr, zi, n, n));
        return .{ .values_r = values_r, .values_i = values_i, .vectors_r = zr, .vectors_i = zi };
    }

    const ortr = try allocator.alloc(f64, n);
    defer allocator.free(ortr);
    const orti = try allocator.alloc(f64, n);
    defer allocator.free(orti);
    @memset(ortr, 0.0);
    @memset(orti, 0.0);

    hessenberg.corth(work, 1, n, ortr, orti);

    const zr = try allocator.alloc(f64, n * n);
    errdefer allocator.free(zr);
    const zi = try allocator.alloc(f64, n * n);
    errdefer allocator.free(zi);
    @memset(zr, 0.0);
    @memset(zi, 0.0);

    const job: usize = if (want_vectors) 2 else 0;
    const ierr = eigen_general.comqr3(work, 1, n, ortr, orti, values_r, values_i, CxMat.initPacked(zr, zi, n, n), job);
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


test "eig dispatches to the Hermitian path for a symmetric matrix" {
    // A = [[2,1],[1,2]]: eigenvalues 1 and 3.
    var ar = [_]f64{ 2.0, 1.0, 1.0, 2.0 };
    var ai: [4]f64 = @splat(0.0);
    const result = try eig(zigstd.testAllocator, CxMat.initPacked(&ar, &ai, 2, 2), true);
    defer result.deinit(zigstd.testAllocator);

    const lo = @min(result.values_r[0], result.values_r[1]);
    const hi = @max(result.values_r[0], result.values_r[1]);
    try zigstd.testExpectApproxEqAbs(1.0, lo, 1e-9);
    try zigstd.testExpectApproxEqAbs(3.0, hi, 1e-9);
    try zigstd.testExpectTrue(result.vectors_r != null);
}

test "eig dispatches to the general path for a non-symmetric matrix" {
    var ar = [_]f64{ 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 };
    var ai: [9]f64 = @splat(0.0);
    const result = try eig(zigstd.testAllocator, CxMat.initPacked(&ar, &ai, 3, 3), false);
    defer result.deinit(zigstd.testAllocator);

    try zigstd.testExpectTrue(result.vectors_r == null);
    var trace: f64 = 0.0;
    for (result.values_r) |v| trace += v;
    try zigstd.testExpectApproxEqAbs(16.0, trace, 1e-8);
}
