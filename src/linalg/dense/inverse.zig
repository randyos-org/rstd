//! Determinant, inverse and solve.

const sysinfo = @import("builtin");
const std = @import("std");
const Allocator = std.mem.Allocator;
const Io = std.Io;
const blas = @import("../../math/kernels/blas/__root__.zig");
const ext = @import("../../math/kernels/blas/ext/interface.zig");
const l3 = @import("../../math/kernels/blas/blas3/core.zig");
const l2 = @import("../../math/kernels/blas/blas2/core.zig");
const l1 = @import("../../math/kernels/blas/blas1/interface.zig");
const lapack = @import("../../math/kernels/lapack/__root__.zig");
const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const ops = @import("../../math/kernels/lapack/ops.zig");

const View = matrix.View;
const ConstView = matrix.ConstView;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const eMul = matrix.eMul;
const eMulAdd = matrix.eMulAdd;

const eigen = @import("eigen.zig");
const isHermitian = eigen.isHermitian;
const util = @import("util.zig");
const ElementOf = util.ElementOf;
const Error = util.Error;
const Path = util.Path;
const constOf = util.constOf;
const copyInto = util.copyInto;
const copyOf = util.copyOf;
const fixed = util.fixed;
const narrow = util.narrow;
const pathOf = util.pathOf;
const viewOf = util.viewOf;

pub const DetError = error{NotSquare} || Allocator.Error || Io.Cancelable;

/// The determinant (MATLAB's `det`); 1 for 0 by 0, 0 when singular.
/// `linalg.matrix.det` up to the inline limit; else the product of `getrf`'s
/// pivots, renormalized by powers of 2 as it goes so that it overflows or
/// underflows only if the result does (design.md, M1's E1). The scaling is
/// exact, so where the plain product stays in range it has its bits.
pub fn det(allocator: Allocator, io: Io, a: anytype) DetError!ElementOf(@TypeOf(a)) {
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    return detWith(pathOf(&.{n}), allocator, io, ca);
}

pub fn detWith(path: Path, allocator: Allocator, io: Io, ca: anytype) (Allocator.Error || Io.Cancelable)!@TypeOf(ca).Element {
    const T = @TypeOf(ca).Element;
    const n = ca.rows;
    if (n == 0) return ops.one(T);
    if (path == .fixed) return fixed(T, "det", n, .{ca});
    const f = try copyOf(T, allocator, ca);
    defer f.deinit(allocator);
    const ipiv = try allocator.alloc(usize, n);
    defer allocator.free(ipiv);
    if (lapack.getrf(allocator, io, T, n, n, f.data, n, ipiv) catch |e| return narrow(Error, e)) |_| return ops.zero(T);
    var d = ops.one(T);
    var exp: i64 = 0;
    for (ipiv, 0..) |p, i| {
        d = eMul(d, f.at(i, i));
        if (p != i) d = ops.neg(T, d);
        const x = std.math.frexp(if (comptime isComplex(T)) @max(@abs(d.re), @abs(d.im)) else d);
        exp += x.exponent;
        d = scaleBy(d, -x.exponent);
    }
    return scaleBy(d, @intCast(std.math.clamp(exp, -1 << 16, 1 << 16)));
}

/// `x 2^e`.
pub fn scaleBy(x: anytype, e: i32) @TypeOf(x) {
    if (comptime !isComplex(@TypeOf(x))) return std.math.ldexp(x, e);
    return .{ .re = std.math.ldexp(x.re, e), .im = std.math.ldexp(x.im, e) };
}

pub const InvError = error{ NotSquare, Singular } || Allocator.Error || Io.Cancelable;

/// `a⁻¹` (MATLAB's `inv`), owned by the caller: `linalg.matrix.inv` up to the
/// inline limit, else `getrf` and `getri`; `error.Singular` with no result.
pub fn inv(allocator: Allocator, io: Io, a: anytype) InvError!Mat(ElementOf(@TypeOf(a))) {
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const r = try Mat(ElementOf(@TypeOf(a))).init(allocator, ca.rows, ca.cols);
    errdefer r.deinit(allocator);
    try invWith(pathOf(&.{ca.rows}), allocator, io, r.view(), ca);
    return r;
}

/// `out = a⁻¹`, `out` untouched on an error.
pub fn invInto(allocator: Allocator, io: Io, out: anytype, a: anytype) (InvError || error{DimensionMismatch})!void {
    const o = viewOf(out);
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    if (o.rows != ca.rows or o.cols != ca.cols) return error.DimensionMismatch;
    return invWith(pathOf(&.{ca.rows}), allocator, io, o, ca);
}

pub fn invWith(path: Path, allocator: Allocator, io: Io, o: anytype, a: ConstView(@TypeOf(o).Element)) (error{Singular} || Allocator.Error || Io.Cancelable)!void {
    const T = @TypeOf(o).Element;
    const n = a.rows;
    if (n == 0) return;
    if (path == .fixed) return fixed(error{Singular}!void, "inv", n, .{ o, a });
    const f = try copyOf(T, allocator, a);
    defer f.deinit(allocator);
    const ipiv = try allocator.alloc(usize, n);
    defer allocator.free(ipiv);
    if (lapack.getrf(allocator, io, T, n, n, f.data, n, ipiv) catch |e| return narrow(Error, e)) |_| return error.Singular;
    lapack.getri(allocator, io, T, n, f.data, n, ipiv) catch |e| return narrow(error{Singular} || Error, e);
    copyInto(T, o.data, o.ld, f.constView());
}

pub const SolveError = error{ DimensionMismatch, Singular } || Allocator.Error || Io.Cancelable;

/// `a \ b` (MATLAB's `mldivide`), owned by the caller: square `a` by
/// `linalg.matrix.solve`'s LU within the inline limit, else by `structureOf`,
/// `error.Singular` with no result; rectangular `a` by `gelsy`'s least-norm
/// least squares, rank cut at `max(m, n) eps` (MATLAB's QR tolerance).
pub fn solve(allocator: Allocator, io: Io, a: anytype, b: anytype) SolveError!Mat(ElementOf(@TypeOf(a))) {
    const ca = constOf(a);
    const cb = constOf(b);
    if (cb.rows != ca.rows) return error.DimensionMismatch;
    const r = try Mat(ElementOf(@TypeOf(a))).init(allocator, ca.cols, cb.cols);
    errdefer r.deinit(allocator);
    try solveWith(pathOf(&.{ ca.rows, cb.cols }), allocator, io, r.view(), ca, cb);
    return r;
}

/// `out = a \ b`, `out` untouched on an error.
pub fn solveInto(allocator: Allocator, io: Io, out: anytype, a: anytype, b: anytype) SolveError!void {
    const o = viewOf(out);
    const ca = constOf(a);
    const cb = constOf(b);
    if (cb.rows != ca.rows or o.rows != ca.cols or o.cols != cb.cols) return error.DimensionMismatch;
    return solveWith(pathOf(&.{ ca.rows, cb.cols }), allocator, io, o, ca, cb);
}

pub const Structure = enum { upper, lower, hermitian, general };

/// `solve`'s route past the inline limit, from `a`'s values alone (MATLAB's
/// order): triangular `trtrs`; Hermitian with a positive real diagonal
/// `potrf` and `potrs`, LU if that fails; else `gesv`.
pub fn structureOf(a: anytype) Structure {
    const T = @TypeOf(a).Element;
    const n = a.rows;
    var upper = true;
    var lower = true;
    for (0..n) |j| for (0..n) |i| {
        if (ops.isZero(T, a.at(i, j))) continue;
        if (i > j) upper = false;
        if (i < j) lower = false;
    };
    if (upper) return .upper;
    if (lower) return .lower;
    for (0..n) |j| {
        const d = a.at(j, j);
        const positive = if (comptime isComplex(T)) d.re > 0 else d > 0;
        if (!positive) return .general;
    }
    return if (isHermitian(a)) .hermitian else .general;
}

pub fn solveWith(path: Path, allocator: Allocator, io: Io, o: anytype, a: ConstView(@TypeOf(o).Element), b: ConstView(@TypeOf(o).Element)) SolveError!void {
    const T = @TypeOf(o).Element;
    const n = a.cols;
    if (a.rows != n) return leastSquares(allocator, io, o, a, b);
    if (n == 0) return;
    if (path == .fixed) return fixed(error{Singular}!void, "solve", n, .{ o, a, b });
    const x = try copyOf(T, allocator, b);
    defer x.deinit(allocator);
    const ld = @max(n, 1);
    switch (structureOf(a)) {
        .upper, .lower => |s| {
            const uplo: lapack.Uplo = if (s == .upper) .u else .l;
            lapack.trtrs(allocator, io, T, uplo, .n, .n, n, b.cols, a.data, a.ld, x.data, ld) catch |e| return narrow(error{Singular} || Error, e);
        },
        .hermitian, .general => |s| {
            const f = try copyOf(T, allocator, a);
            defer f.deinit(allocator);
            if (s == .hermitian) cholesky: {
                lapack.potrf(allocator, io, T, .u, n, f.data, ld) catch |e| switch (e) {
                    error.NotPositiveDefinite => {
                        copyInto(T, f.data, n, a);
                        break :cholesky;
                    },
                    else => return narrow(Error, e),
                };
                lapack.potrs(allocator, io, T, .u, n, b.cols, f.data, ld, x.data, ld) catch |e| return narrow(Error, e);
                return copyInto(T, o.data, o.ld, x.constView());
            }
            const ipiv = try allocator.alloc(usize, n);
            defer allocator.free(ipiv);
            lapack.gesv(allocator, io, T, n, b.cols, f.data, ld, ipiv, x.data, ld) catch |e| return narrow(error{Singular} || Error, e);
        },
    }
    copyInto(T, o.data, o.ld, x.constView());
}

/// `gelsy`'s least-norm least-squares solution, for rectangular `a`.
pub fn leastSquares(allocator: Allocator, io: Io, o: anytype, a: ConstView(@TypeOf(o).Element), b: ConstView(@TypeOf(o).Element)) (Allocator.Error || Io.Cancelable)!void {
    const T = @TypeOf(o).Element;
    const F = Real(T);
    const m = a.rows;
    const n = a.cols;
    if (m == 0 or n == 0) {
        for (0..o.cols) |j| @memset(o.col(j), ops.zero(T));
        return;
    }
    const f = try copyOf(T, allocator, a);
    defer f.deinit(allocator);
    const ldx = @max(m, n);
    const x = try allocator.alloc(T, ldx * b.cols);
    defer allocator.free(x);
    copyInto(T, x, ldx, b);
    const jpvt = try allocator.alloc(usize, n);
    defer allocator.free(jpvt);
    @memset(jpvt, 0);
    const rcond = @as(F, @floatFromInt(ldx)) * std.math.floatEps(F);
    _ = lapack.gelsy(allocator, io, T, m, n, b.cols, f.data, m, x, ldx, jpvt, rcond) catch |e| return narrow(Error, e);
    for (0..o.cols) |j| @memcpy(o.col(j), x[j * ldx ..][0..n]);
}
