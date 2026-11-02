//! The matlab port's LU-based operations, `det` and `inv` here and the
//! divisions in `divide.zig`, over rstd's LAPACK (`getrf`, `getrs`,
//! `getri`, `gecon`; M2) in place of matlab88.f90's WGEFA, WGESL, WGECO and
//! WGEDI. They convert at the boundary (`CxMat.zig`): a real matrix goes to
//! the real routines in place on its `r`, a complex one through one
//! interleaved copy. WGEDI's determinant scaling stays (M1's E1,
//! `linalg/design.md`): a mantissa times `10^exponent` over `getrf`'s
//! diagonal, where a plain product of it overflows.

const std = @import("std");
const scalar = @import("scalar.zig");
const blas = @import("../blas/__root__.zig");
const CxMat = @import("CxMat.zig");
const lapack = @import("../lapack/__root__.zig");
const Cx = scalar.Cx;
const Complex = CxMat.Complex;
const Allocator = std.mem.Allocator;
const Io = std.Io;

/// What the LU operations return: `Singular` for an exact zero pivot where
/// the result needs `U`'s inverse.
pub const Error = error{Singular} || Allocator.Error || Io.Cancelable;

pub const Determinant = struct { mantissa: Cx, exponent10: f64 };

/// `getrf`'s `P L U` of `a` (`n x n`): in place on `a.r` for a real `a`,
/// else (or for `complex`) in `z`, an interleaved copy. A zero pivot is not
/// an error here: the factors are complete (`singular`), as `getrf` leaves
/// them.
pub const Lu = struct {
    a: CxMat,
    ipiv: []isize,
    z: ?[]Complex,
    singular: bool,

    pub fn init(a: CxMat, complex: bool, allocator: Allocator, io: Io) Error!Lu {
        const n = a.n;
        const ipiv = try allocator.alloc(isize, n);
        errdefer allocator.free(ipiv);
        const z: ?[]Complex = if (complex or !a.isReal()) try a.toPacked(Complex, allocator) else null;
        errdefer if (z) |x| allocator.free(x);
        var singular = false;
        const ni: isize = @intCast(n);
        const r = if (z) |x|
            lapack.getrf(allocator, io, Complex, ni, ni, x.ptr, ni, ipiv.ptr, null)
        else
            lapack.getrf(allocator, io, f64, ni, ni, a.r.ptr, @intCast(a.lda), ipiv.ptr, null);
        r catch |e| switch (e) {
            error.Singular => singular = true,
            error.OutOfMemory, error.Canceled => |x| return x,
            // the arguments come from `a`
            else => unreachable,
        };
        return .{ .a = a, .ipiv = ipiv, .z = z, .singular = singular };
    }

    pub fn deinit(f: Lu, allocator: Allocator) void {
        allocator.free(f.ipiv);
        if (f.z) |x| allocator.free(x);
    }

    /// `op(A) X = B` in place of `b` (`n x p`), by `getrs`; a complex `b`
    /// needs `f` factored `complex`.
    pub fn solve(f: Lu, trans: lapack.Trans, b: CxMat, allocator: Allocator, io: Io) Error!void {
        const ni: isize = @intCast(f.a.n);
        const p: isize = @intCast(b.n);
        if (f.z) |x| {
            const zb = try b.toPacked(Complex, allocator);
            defer allocator.free(zb);
            try solved(lapack.getrs(allocator, io, Complex, trans, ni, p, x.ptr, ni, f.ipiv.ptr, zb.ptr, ni));
            b.fromPacked(Complex, zb);
        } else try solved(lapack.getrs(allocator, io, f64, trans, ni, p, f.a.r.ptr, @intCast(f.a.lda), f.ipiv.ptr, b.r.ptr, @intCast(b.lda)));
    }

    fn solved(r: lapack.GetrsError!void) Error!void {
        r catch |e| switch (e) {
            error.OutOfMemory, error.Canceled => |x| return x,
            // the arguments come from the factors
            else => unreachable,
        };
    }

    /// `U(k, k)`, 1-based.
    fn diag(f: Lu, k: usize) Cx {
        if (f.z) |x| {
            const e = x[(k - 1) * (f.a.n + 1)];
            return Cx.init(e.re, e.im);
        }
        return Cx.init(f.a.r[f.a.idx(k, k)], 0);
    }
};

/// DET: the determinant of the square `self` (destroyed), as WGEDI's job 10
/// scales it over `getrf`'s pivots and diagonal.
pub fn det(self: CxMat, allocator: Allocator, io: Io) Error!Determinant {
    const f = try Lu.init(self, false, allocator, io);
    defer f.deinit(allocator);
    var d = Cx.init(1.0, 0.0);
    var exp10: f64 = 0.0;
    for (1..self.n + 1) |i| {
        if (f.ipiv[i - 1] != i) d = Cx.init(-d.re, -d.im);
        d = f.diag(i).mul(d);
        if (blas.cabs1(Cx, d) == 0.0) break;
        while (blas.cabs1(Cx, d) < 1.0) {
            d = Cx.init(d.re * 10.0, d.im * 10.0);
            exp10 -= 1.0;
        }
        while (blas.cabs1(Cx, d) >= 10.0) {
            d = Cx.init(d.re / 10.0, d.im / 10.0);
            exp10 += 1.0;
        }
    }
    return .{ .mantissa = d, .exponent10 = exp10 };
}

/// The 1-norm of the square `self`, as `gecon` wants it before `getrf`.
pub fn norm1(self: CxMat, allocator: Allocator, io: Io) Allocator.Error!f64 {
    const ni: isize = @intCast(self.n);
    if (self.isReal()) return lapack.lange(allocator, io, f64, .o, ni, ni, self.r.ptr, @intCast(self.lda)) catch |e| switch (e) {
        error.OutOfMemory => |x| x,
        else => unreachable,
    };
    const z = try self.toPacked(Complex, allocator);
    defer allocator.free(z);
    return lapack.lange(allocator, io, Complex, .o, ni, ni, z.ptr, ni) catch |e| switch (e) {
        error.OutOfMemory => |x| x,
        else => unreachable,
    };
}

/// `gecon`'s 1-norm `rcond` of the factored `f`, `anorm` its matrix's 1-norm;
/// 0 where the estimate fails (an infinite or NaN inverse norm), as
/// numerically singular.
pub fn rcond(f: Lu, allocator: Allocator, io: Io, anorm: f64) Error!f64 {
    var rc: f64 = 0;
    const ni: isize = @intCast(f.a.n);
    const r = if (f.z) |x|
        lapack.gecon(allocator, io, Complex, .o, ni, x.ptr, ni, anorm, &rc)
    else
        lapack.gecon(allocator, io, f64, .o, ni, f.a.r.ptr, @intCast(f.a.lda), anorm, &rc);
    r catch |e| switch (e) {
        error.EstimateFailed => return 0,
        error.OutOfMemory, error.Canceled => |x| return x,
        else => unreachable,
    };
    return rc;
}

/// INV: `self` (square) inverted in place, by `getrf`, `gecon` and `getri`.
/// Returns the reciprocal 1-norm condition estimate (WGECO's `rcond`): so
/// small that `1 + rcond == 1` means the inverse is not to be trusted. An
/// exact zero pivot is `error.Singular`, `self` then destroyed.
pub fn inv(self: CxMat, allocator: Allocator, io: Io) Error!f64 {
    const anorm = try norm1(self, allocator, io);
    const f = try Lu.init(self, false, allocator, io);
    defer f.deinit(allocator);
    if (f.singular) return error.Singular;
    const rc = try rcond(f, allocator, io, anorm);
    const ni: isize = @intCast(self.n);
    const r = if (f.z) |x|
        lapack.getri(allocator, io, Complex, ni, x.ptr, ni, f.ipiv.ptr, null)
    else
        lapack.getri(allocator, io, f64, ni, self.r.ptr, @intCast(self.lda), f.ipiv.ptr, null);
    r catch |e| switch (e) {
        error.Singular, error.OutOfMemory, error.Canceled => |x| return x,
        else => unreachable,
    };
    if (f.z) |x| self.fromPacked(Complex, x);
    return rc;
}

const zigstd = @import("../../zigstd.zig");
const testing = std.testing;

test "det and inv of a diagonal matrix, real and complex" {
    // A = diag(2, 3, -4), then with an imaginary part
    for ([_]f64{ 0, 1 }) |im| {
        var ar = [_]f64{ 2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, -4.0 };
        var ai = [_]f64{ im, 0, 0, 0, 0, 0, 0, 0, 0 };
        var br = ar;
        var bi = ai;
        const d = try det(CxMat.initPacked(&ar, &ai, 3, 3), testing.allocator, testing.io);
        // (2 + i im) 3 (-4)
        try zigstd.testExpectApproxEqAbs(-24.0, d.mantissa.re * std.math.pow(f64, 10.0, d.exponent10), 1e-9);
        try zigstd.testExpectApproxEqAbs(-12.0 * im, d.mantissa.im * std.math.pow(f64, 10.0, d.exponent10), 1e-9);
        const rc = try inv(CxMat.initPacked(&br, &bi, 3, 3), testing.allocator, testing.io);
        try zigstd.testExpectTrue(rc > 0.0);
        // 1 / (2 + i im), 1/3, -1/4
        const q = 4 + im * im;
        try zigstd.testExpectApproxEqAbs(2 / q, br[0], 1e-12);
        try zigstd.testExpectApproxEqAbs(-im / q, bi[0], 1e-12);
        try zigstd.testExpectApproxEqAbs(1.0 / 3.0, br[scalar.ix(3, 2, 2)], 1e-12);
        try zigstd.testExpectApproxEqAbs(-0.25, br[scalar.ix(3, 3, 3)], 1e-12);
    }
}

test "det keeps its range where the pivots' product overflows (M1's E1); a singular inv fails" {
    // det = -3e300 exactly representable, but 1e200 * 1e200 overflows on the way
    var ar = [_]f64{ 1e200, 0, 0, 0, 1e200, 0, 0, 0, -3e-100 };
    var ai: [9]f64 = @splat(0.0);
    var product: f64 = 1;
    for (1..4) |k| product *= ar[scalar.ix(3, k, k)];
    try zigstd.testExpectTrue(std.math.isInf(product));
    const d = try det(CxMat.initPacked(&ar, &ai, 3, 3), testing.allocator, testing.io);
    try zigstd.testExpectApproxEqAbs(-3.0, d.mantissa.re, 1e-12);
    try zigstd.testExpectEqual(@as(f64, 0), d.mantissa.im);
    try zigstd.testExpectEqual(@as(f64, 300), d.exponent10);
    // [[1, 2], [2, 4]]: a zero pivot, det 0, inv Singular
    var sr = [_]f64{ 1, 2, 2, 4 };
    var si: [4]f64 = @splat(0.0);
    var tr = sr;
    var ti = si;
    const z = try det(CxMat.initPacked(&sr, &si, 2, 2), testing.allocator, testing.io);
    try zigstd.testExpectEqual(@as(f64, 0), z.mantissa.re);
    try testing.expectError(error.Singular, inv(CxMat.initPacked(&tr, &ti, 2, 2), testing.allocator, testing.io));
}
