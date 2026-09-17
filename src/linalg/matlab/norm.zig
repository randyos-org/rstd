//! `NORM`, ported from `matlab88.f90`'s `MATFN5` command body. Split into
//! separate vector/matrix functions rather than one shape-dispatching
//! `norm()`, matching this port's existing convention (see `matrixfun.zig`'s
//! doc comment on why `sin`/`cos`/... got the same treatment).
//!
//! Vector: `NORM(V,P) = (sum(|V(I)|**P))**(1/P)`, `NORM(V) = NORM(V,2)`,
//! `NORM(V,'INF') = MAX(ABS(V(I)))`.
//! Matrix: `NORM(X)`/`NORM(X,2)` is the largest singular value,
//! `NORM(X,1)` the 1-norm, `NORM(X,'INF')` the infinity-norm,
//! `NORM(X,'FRO')` the Frobenius norm.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1.zig");
const svd = @import("svd.zig");
const ix = blas1.ix;

/// `NORM(V,2)` / `NORM(V)`: Euclidean norm. Alias of `blas1.wnrm2`, which
/// is already the overflow-safe incremental-`pythag` algorithm the
/// original used for this.
pub const vecNorm2 = blas1.wnrm2;

/// `NORM(V,P)`: general p-norm, `(sum(|v_i|**p))**(1/p)`. For `p == 2`,
/// prefer `vecNorm2` (`wnrm2`), which is overflow-safe; this general form
/// isn't (matches the original, which didn't special-case overflow safety
/// for arbitrary `p` either).
pub fn vecNormP(n: usize, xr: []const f64, xi: []const f64, incx: isize, p: f64) f64 {
    var acc: f64 = 0.0;
    var i_x: isize = if (incx < 0) (-@as(isize, @intCast(n)) + 1) * incx else 0;
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const u: usize = @intCast(i_x);
        acc += std.math.pow(f64, blas1.pythag(xr[u], xi[u]), p);
        i_x += incx;
    }
    return std.math.pow(f64, acc, 1.0 / p);
}

/// `NORM(V,'INF')`: max magnitude over `v`'s entries.
pub fn vecNormInf(n: usize, xr: []const f64, xi: []const f64, incx: isize) f64 {
    var m: f64 = 0.0;
    var i_x: isize = if (incx < 0) (-@as(isize, @intCast(n)) + 1) * incx else 0;
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const u: usize = @intCast(i_x);
        m = @max(m, blas1.pythag(xr[u], xi[u]));
        i_x += incx;
    }
    return m;
}

/// `NORM(X)` / `NORM(X,2)`: largest singular value of `ar`/`ai` (`n`x`p`,
/// column-major, leading dimension `n` - destroyed, matching `wsvdc`).
/// `sr`/`si` (length `min(n+1,p)`), `er`/`ei` (length `p`), and
/// `workr`/`worki` (length `n`) are `wsvdc`'s own scratch - see
/// `svdapi.zig`'s `singularValues` (the same job-0, no-vectors call shape
/// this mirrors) for the exact sizes. Only `sr[0]` (the largest singular
/// value) is read back here.
pub fn matNorm2(
    ar: []f64,
    ai: []f64,
    n: usize,
    p: usize,
    sr: []f64,
    si: []f64,
    er: []f64,
    ei: []f64,
    workr: []f64,
    worki: []f64,
) usize {
    var no_vecs: [1]f64 = .{0.0};
    return svd.wsvdc(ar, ai, n, n, p, sr, si, er, ei, &no_vecs, &no_vecs, 1, &no_vecs, &no_vecs, 1, workr, worki, 0);
}

/// `NORM(X,1)`: max column sum of `|re|+|im|` (matches `lu.wgeco`'s own
/// `anorm` computation, the same convention this port already uses
/// internally for a matrix 1-norm).
pub fn matNorm1(ar: []const f64, ai: []const f64, lda: usize, m: usize, n: usize) f64 {
    var result: f64 = 0.0;
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        result = @max(result, blas1.wasum(m, ar[ix(lda, 1, j)..], ai[ix(lda, 1, j)..], 1));
    }
    return result;
}

/// `NORM(X,'INF')`: max row sum of `|re|+|im|`. A row is strided (not
/// contiguous) in this column-major layout, so this doesn't take
/// `wasum`'s `@Vector` fast path - matching every other row-wise access in
/// this port.
pub fn matNormInf(ar: []const f64, ai: []const f64, lda: usize, m: usize, n: usize) f64 {
    var result: f64 = 0.0;
    var i: usize = 1;
    while (i <= m) : (i += 1) {
        result = @max(result, blas1.wasum(n, ar[ix(lda, i, 1)..], ai[ix(lda, i, 1)..], @intCast(lda)));
    }
    return result;
}

/// `NORM(X,'FRO')`: Frobenius norm, `sqrt(sum(|x_ij|**2))` - the whole
/// `m`x`n` matrix read as one flat length-`m*n` vector and passed straight
/// to `wnrm2` (a column-major matrix's storage already *is* one flat
/// contiguous run).
pub fn matNormFro(ar: []const f64, ai: []const f64, m: usize, n: usize) f64 {
    return blas1.wnrm2(m * n, ar[0 .. m * n], ai[0 .. m * n], 1);
}

test "vecNorm2 matches wnrm2 and the 3-4-5 triangle" {
    const xr = [_]f64{ 3.0, 0.0 };
    const xi = [_]f64{ 0.0, 4.0 };
    try zigstd.testExpectApproxEqAbs(5.0, vecNorm2(2, &xr, &xi, 1), 1e-12);
}

test "vecNormP with p=2 matches vecNorm2, and p=1 is the sum of moduli" {
    // moduli: pythag(3,4)=5, pythag(1,0)=1, pythag(-2,0)=2.
    const xr = [_]f64{ 3.0, 1.0, -2.0 };
    const xi = [_]f64{ 4.0, 0.0, 0.0 };
    try zigstd.testExpectApproxEqAbs(vecNorm2(3, &xr, &xi, 1), vecNormP(3, &xr, &xi, 1, 2.0), 1e-9);
    try zigstd.testExpectApproxEqAbs(8.0, vecNormP(3, &xr, &xi, 1, 1.0), 1e-9);
}

test "vecNormInf finds the largest-magnitude entry" {
    // moduli: 1, 2, pythag(0.5,9)=sqrt(81.25)=9.0138782...
    const xr = [_]f64{ 1.0, -2.0, 0.5 };
    const xi = [_]f64{ 0.0, 0.0, 9.0 };
    try zigstd.testExpectApproxEqAbs(@sqrt(81.25), vecNormInf(3, &xr, &xi, 1), 1e-9);
}

test "matNorm1/matNormInf/matNormFro match hand computation for a real 2x2" {
    // A = [[1, -3], [4, 2]] (column-major: col1=[1,4], col2=[-3,2]).
    const ar = [_]f64{ 1.0, 4.0, -3.0, 2.0 };
    const ai: [4]f64 = @splat(0.0);

    // 1-norm: max column sum of abs = max(1+4, 3+2) = 5.
    try zigstd.testExpectApproxEqAbs(5.0, matNorm1(&ar, &ai, 2, 2, 2), 1e-12);
    // inf-norm: max row sum of abs = max(1+3, 4+2) = 6.
    try zigstd.testExpectApproxEqAbs(6.0, matNormInf(&ar, &ai, 2, 2, 2), 1e-12);
    // Frobenius: sqrt(1+16+9+4) = sqrt(30).
    try zigstd.testExpectApproxEqAbs(@sqrt(30.0), matNormFro(&ar, &ai, 2, 2), 1e-9);
}

test "matNorm2 (largest singular value) matches a known case" {
    // A diagonal matrix's singular values are the absolute diagonal
    // entries, so the largest is trivially known.
    var ar = [_]f64{ 3.0, 0.0, 0.0, -5.0 };
    var ai: [4]f64 = @splat(0.0);
    var sr: [2]f64 = undefined;
    var si: [2]f64 = undefined;
    var er: [2]f64 = undefined;
    var ei: [2]f64 = undefined;
    var workr: [2]f64 = undefined;
    var worki: [2]f64 = undefined;

    const info = matNorm2(&ar, &ai, 2, 2, &sr, &si, &er, &ei, &workr, &worki);
    try zigstd.testExpectEqual(@as(usize, 0), info);
    try zigstd.testExpectApproxEqAbs(5.0, sr[0], 1e-9);
}
