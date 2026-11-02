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
const scalar = @import("scalar.zig");
const CxMat = @import("CxMat.zig");
const svdapi = @import("svdapi.zig");
const blas = @import("../blas/__root__.zig");
const level1 = scalar.level1;
const Io = std.Io;

/// `NORM(V,2)` / `NORM(V)`: Euclidean norm, `hypot` of `dnrm2` over the
/// real and imaginary parts (M2: BLAS 1, which scales where WNRM2's fast
/// path lost digits, M1's E8). `incx` as BLAS's: `xr`, `xi` start at the
/// storage's first element, a negative stride walking it backward.
pub fn vecNorm2(io: Io, n: usize, xr: []const f64, xi: []const f64, incx: isize) blas.Level1Error!f64 {
    const ni: isize = @intCast(n);
    const re = try blas.nrm2(io, f64, ni, @constCast(xr.ptr), incx);
    return std.math.hypot(re, try blas.nrm2(io, f64, ni, @constCast(xi.ptr), incx));
}

/// `NORM(V,P)`: general p-norm, `(sum(|v_i|**p))**(1/p)`. For `p == 2`,
/// prefer `vecNorm2`, which is overflow-safe; this general form
/// isn't (matches the original, which didn't special-case overflow safety
/// for arbitrary `p` either).
pub fn vecNormP(n: usize, xr: []const f64, xi: []const f64, incx: isize, p: f64) f64 {
    var acc: f64 = 0.0;
    var i_x: isize = if (incx < 0) (-@as(isize, @intCast(n)) + 1) * incx else 0;
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const u: usize = @intCast(i_x);
        acc += std.math.pow(f64, std.math.hypot(xr[u], xi[u]), p);
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
        m = @max(m, std.math.hypot(xr[u], xi[u]));
        i_x += incx;
    }
    return m;
}

/// `NORM(X)` / `NORM(X,2)`: largest singular value of `self` (not
/// mutated), by `gesvd`; 0 for an empty `self`.
pub fn matNorm2(self: CxMat, allocator: std.mem.Allocator, io: std.Io) svdapi.Error!f64 {
    if (self.m == 0 or self.n == 0) return 0;
    const s = try svdapi.singularValues(self, allocator, io);
    defer allocator.free(s);
    return s[0];
}

/// `dasum` of `n` entries of `r` and of `i`, `inc` apart: the sum of
/// `|re| + |im|`.
fn asum(io: Io, n: usize, r: []f64, i: []f64, inc: usize) Io.Cancelable!f64 {
    const ni: isize = @intCast(n);
    const si: isize = @intCast(inc);
    return try level1(blas.asum(io, f64, ni, r.ptr, si)) + try level1(blas.asum(io, f64, ni, i.ptr, si));
}

/// `NORM(X,1)`: max column sum of `|re|+|im|` (WGECO's `anorm`), by
/// `dasum`.
pub fn matNorm1(self: CxMat, io: Io) Io.Cancelable!f64 {
    var result: f64 = 0.0;
    for (1..self.n + 1) |j| {
        const c = self.col(j);
        result = @max(result, try asum(io, self.m, c.r, c.i, 1));
    }
    return result;
}

/// `NORM(X,'INF')`: max row sum of `|re|+|im|`, by `dasum` along each row.
pub fn matNormInf(self: CxMat, io: Io) Io.Cancelable!f64 {
    var result: f64 = 0.0;
    for (1..self.m + 1) |i| {
        const row = self.rowView(i);
        result = @max(result, try asum(io, self.n, row.r, row.i, self.lda));
    }
    return result;
}

/// `NORM(X,'FRO')`: Frobenius norm, `sqrt(sum(|x_ij|**2))`, `hypot` over
/// the columns' `vecNorm2`.
pub fn matNormFro(self: CxMat, io: Io) Io.Cancelable!f64 {
    var result: f64 = 0.0;
    for (1..self.n + 1) |j| {
        const c = self.col(j);
        result = std.math.hypot(result, try level1(vecNorm2(io, self.m, c.r, c.i, 1)));
    }
    return result;
}

const zigstd = @import("../../zigstd.zig");

const test_io = std.testing.io;

test "vecNorm2 matches the 3-4-5 triangle" {
    const xr = [_]f64{ 3.0, 0.0 };
    const xi = [_]f64{ 0.0, 4.0 };
    try zigstd.testExpectApproxEqAbs(5.0, try vecNorm2(test_io, 2, &xr, &xi, 1), 1e-15);
}

test "vecNorm2 scales: huge, tiny, partly underflowing and zero vectors (M1's E8)" {
    // squares of 1e300 overflow and of 1e-200 underflow; at 1e-160 the
    // squares are subnormal, where WNRM2's fast path lost digits
    for ([_]f64{ 1e300, 1e-200, 1e-160 }) |x| {
        const xr = [_]f64{ x, x };
        const xi = [_]f64{ 0, x };
        const want = @sqrt(3.0) * x;
        try zigstd.testExpectApproxEqAbs(want, try vecNorm2(test_io, 2, &xr, &xi, 1), want * 4 * std.math.floatEps(f64));
    }
    const zr: [4]f64 = @splat(0.0);
    try zigstd.testExpectEqual(@as(f64, 0.0), try vecNorm2(test_io, 4, &zr, &zr, 1));
}

test "vecNormP with p=2 matches vecNorm2, and p=1 is the sum of moduli" {
    // moduli: hypot(3,4)=5, hypot(1,0)=1, hypot(-2,0)=2.
    const xr = [_]f64{ 3.0, 1.0, -2.0 };
    const xi = [_]f64{ 4.0, 0.0, 0.0 };
    try zigstd.testExpectApproxEqAbs(try vecNorm2(test_io, 3, &xr, &xi, 1), vecNormP(3, &xr, &xi, 1, 2.0), 1e-9);
    try zigstd.testExpectApproxEqAbs(8.0, vecNormP(3, &xr, &xi, 1, 1.0), 1e-9);
}

test "vecNormInf finds the largest-magnitude entry" {
    // moduli: 1, 2, hypot(0.5,9)=sqrt(81.25)=9.0138782...
    const xr = [_]f64{ 1.0, -2.0, 0.5 };
    const xi = [_]f64{ 0.0, 0.0, 9.0 };
    try zigstd.testExpectApproxEqAbs(@sqrt(81.25), vecNormInf(3, &xr, &xi, 1), 1e-9);
}

test "matNorm1/matNormInf/matNormFro match hand computation for a real 2x2" {
    // A = [[1, -3], [4, 2]] (column-major: col1=[1,4], col2=[-3,2]).
    var ar = [_]f64{ 1.0, 4.0, -3.0, 2.0 };
    var ai: [4]f64 = @splat(0.0);
    const a = CxMat.initPacked(&ar, &ai, 2, 2);

    // 1-norm: max column sum of abs = max(1+4, 3+2) = 5.
    try zigstd.testExpectApproxEqAbs(5.0, try matNorm1(a, test_io), 1e-15);
    // inf-norm: max row sum of abs = max(1+3, 4+2) = 6.
    try zigstd.testExpectApproxEqAbs(6.0, try matNormInf(a, test_io), 1e-15);
    // Frobenius: sqrt(1+16+9+4) = sqrt(30).
    try zigstd.testExpectApproxEqAbs(@sqrt(30.0), try matNormFro(a, test_io), 1e-15);
}

test "matNorm2 (largest singular value) matches a known case" {
    // A diagonal matrix's singular values are the absolute diagonal
    // entries, so the largest is trivially known.
    var ar = [_]f64{ 3.0, 0.0, 0.0, -5.0 };
    var ai: [4]f64 = @splat(0.0);
    const s = try matNorm2(CxMat.initPacked(&ar, &ai, 2, 2), std.testing.allocator, std.testing.io);
    try zigstd.testExpectApproxEqAbs(5.0, s, 1e-15);
}
