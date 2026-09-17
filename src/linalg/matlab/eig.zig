//! `EIG`'s dispatch logic, split out of `api.zig` to keep files under this
//! project's line-count guideline: pick the Hermitian eigensolver
//! (`tridiagonal.htridi` + `eigen_symmetric.imtql2` + `tridiagonal.htribk`)
//! or the general one (`hessenberg.corth` + `eigen_general.comqr3`),
//! exactly as the original REPL's `EIG` command did.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1.zig");
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
pub fn eig(allocator: zigstd.Allocator, ar: []const f64, ai: []const f64, n: usize, want_vectors: bool) EigError!EigResult {
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


test "eig dispatches to the Hermitian path for a symmetric matrix" {
    // A = [[2,1],[1,2]]: eigenvalues 1 and 3.
    const ar = [_]f64{ 2.0, 1.0, 1.0, 2.0 };
    const ai: [4]f64 = @splat(0.0);
    const result = try eig(zigstd.testAllocator, &ar, &ai, 2, true);
    defer result.deinit(zigstd.testAllocator);

    const lo = @min(result.values_r[0], result.values_r[1]);
    const hi = @max(result.values_r[0], result.values_r[1]);
    try zigstd.testExpectApproxEqAbs(1.0, lo, 1e-9);
    try zigstd.testExpectApproxEqAbs(3.0, hi, 1e-9);
    try zigstd.testExpectTrue(result.vectors_r != null);
}

test "eig dispatches to the general path for a non-symmetric matrix" {
    const ar = [_]f64{ 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 };
    const ai: [9]f64 = @splat(0.0);
    const result = try eig(zigstd.testAllocator, &ar, &ai, 3, false);
    defer result.deinit(zigstd.testAllocator);

    try zigstd.testExpectTrue(result.vectors_r == null);
    var trace: f64 = 0.0;
    for (result.values_r) |v| trace += v;
    try zigstd.testExpectApproxEqAbs(16.0, trace, 1e-8);
}
