//! `rstd.linalg`'s run-time-shape functions (design.md, "Operations over
//! views"), re-exported flat at `rstd.linalg`. Each operand is a `View` or
//! `ConstView`, a `Mat`, a pointer to a `StaticMat` or a slice (a column),
//! in any mix; a vector is a view of one column or one row. Shapes at or
//! below the inline limit take the fixed-size forms' operation order, so a
//! `Mat` and a `StaticMat` holding the same numbers give the same bits; past
//! it, the BLAS and LAPACK. Outputs must not overlap inputs.

const sysinfo = @import("builtin");
const std = @import("std");
const Allocator = std.mem.Allocator;
const Io = std.Io;
const blas = @import("blas/__root__.zig");
const lapack = @import("lapack/__root__.zig");
const vector = @import("vector.zig");
const matrix = @import("matrix.zig");
const ops = @import("lapack/ops.zig");

const View = matrix.View;
const ConstView = matrix.ConstView;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const eMul = matrix.eMul;
const eMulAdd = matrix.eMulAdd;

// ---- vectors ----------------------------------------------------------------

pub const VectorError = error{DimensionMismatch} || Io.Cancelable;

/// `xᴴ y` (`Σ conj(x[i]) y[i]`, MATLAB's `dot`); `linalg.vector.dot`'s order
/// up to the inline limit, else the BLAS's `dot` or `dotc`.
pub fn dot(io: Io, x: anytype, y: anytype) VectorError!ElementOf(@TypeOf(x)) {
    const sx = try strided(constOf(x));
    const sy = try strided(constOf(y));
    if (sx.n != sy.n) return error.DimensionMismatch;
    return dotWith(pathOf(&.{sx.n}), io, sx, sy);
}

fn dotWith(path: Path, io: Io, sx: anytype, sy: @TypeOf(sx)) Io.Cancelable!@TypeOf(sx.ptr[0]) {
    const T = @TypeOf(sx.ptr[0]);
    if (sx.n == 0) return ops.zero(T);
    if (path == .fixed) return fixed(T, "dot", sx.n, .{ sx, sy });
    const f = if (comptime isComplex(T)) blas.dotc else blas.dot;
    return f(io, T, toIsize(sx.n), sx.ptr, toIsize(sx.inc), sy.ptr, toIsize(sy.inc)) catch |e| return narrow(Io.Cancelable, e);
}

/// The 2-norm (MATLAB's `norm` of a vector; `linalg.vector`'s `norm`
/// normalizes): `linalg.vector.length` up to the inline limit, else `nrm2`.
pub fn norm(io: Io, x: anytype) VectorError!Real(ElementOf(@TypeOf(x))) {
    const s = try strided(constOf(x));
    return normWith(pathOf(&.{s.n}), io, s);
}

fn normWith(path: Path, io: Io, s: anytype) Io.Cancelable!Real(@TypeOf(s.ptr[0])) {
    const T = @TypeOf(s.ptr[0]);
    if (s.n == 0) return 0;
    if (path == .fixed) return fixed(Real(T), "length", s.n, .{s});
    return blas.nrm2(io, T, toIsize(s.n), s.ptr, toIsize(s.inc)) catch |e| return narrow(Io.Cancelable, e);
}

/// `Σ x[i]`: in order up to the inline limit, else the BLAS's `sum`, the
/// real and imaginary parts apart for complex.
pub fn sum(io: Io, x: anytype) VectorError!ElementOf(@TypeOf(x)) {
    const s = try strided(constOf(x));
    return sumWith(pathOf(&.{s.n}), io, s);
}

fn sumWith(path: Path, io: Io, s: anytype) Io.Cancelable!@TypeOf(s.ptr[0]) {
    const T = @TypeOf(s.ptr[0]);
    if (s.n == 0) return ops.zero(T);
    if (path == .fixed) return fixed(T, "sum", s.n, .{s});
    if (comptime !isComplex(T)) return blas.sum(io, T, toIsize(s.n), s.ptr, toIsize(s.inc)) catch |e| return narrow(Io.Cancelable, e);
    const F = Real(T);
    const parts: [*]F = @ptrCast(s.ptr);
    const inc = toIsize(2 * s.inc);
    return .{
        .re = blas.sum(io, F, toIsize(s.n), parts, inc) catch |e| return narrow(Io.Cancelable, e),
        .im = blas.sum(io, F, toIsize(s.n), parts + 1, inc) catch |e| return narrow(Io.Cancelable, e),
    };
}

/// The largest element, or the first NaN; ties keep the first; 0 for an
/// empty vector, as the BLAS's `max`. Real only.
pub fn max(io: Io, x: anytype) VectorError!ElementOf(@TypeOf(x)) {
    return extreme(.max, io, x);
}

/// The smallest element, or the first NaN, as `max`.
pub fn min(io: Io, x: anytype) VectorError!ElementOf(@TypeOf(x)) {
    return extreme(.min, io, x);
}

fn extreme(comptime which: blas.Extreme, io: Io, x: anytype) VectorError!ElementOf(@TypeOf(x)) {
    const T = ElementOf(@TypeOf(x));
    if (comptime isComplex(T)) @compileError("max and min order real elements only");
    const s = try strided(constOf(x));
    return extremeWith(which, pathOf(&.{s.n}), io, s);
}

fn extremeWith(comptime which: blas.Extreme, path: Path, io: Io, s: anytype) Io.Cancelable!@TypeOf(s.ptr[0]) {
    const T = @TypeOf(s.ptr[0]);
    if (s.n == 0) return 0;
    if (path == .fixed) return fixed(T, "extreme", s.n, .{ which, s });
    return blas.extreme(io, T, which, toIsize(s.n), s.ptr, toIsize(s.inc)) catch |e| return narrow(Io.Cancelable, e);
}

/// `y += a x`, by the BLAS's `axpy`.
pub fn axpy(io: Io, a: anytype, x: anytype, y: anytype) VectorError!void {
    const T = ElementOf(@TypeOf(y));
    const sx = try strided(constOf(x));
    const sy = try strided(viewOf(y).constView());
    if (sx.n != sy.n) return error.DimensionMismatch;
    const alpha: T = a;
    blas.axpy(io, T, toIsize(sx.n), alpha, sx.ptr, toIsize(sx.inc), sy.ptr, toIsize(sy.inc)) catch |e| return narrow(Io.Cancelable, e);
}

/// `x *= a`, element by element as `linalg.matrix.scale` (design.md,
/// Dispatch: element-wise forms are inline at every size).
pub fn scale(a: anytype, x: anytype) void {
    const v = viewOf(x);
    const alpha: @TypeOf(v).Element = a;
    for (0..v.cols) |j| for (v.col(j)) |*e| {
        e.* = eMul(e.*, alpha);
    };
}

// ---- products and transposes ------------------------------------------------

pub const MulError = error{DimensionMismatch} || Allocator.Error || Io.Cancelable;

/// `a b` (MATLAB's `*`), owned by the caller.
pub fn mul(allocator: Allocator, io: Io, a: anytype, b: anytype) MulError!Mat(ElementOf(@TypeOf(a))) {
    const ca = constOf(a);
    const cb = constOf(b);
    if (ca.cols != cb.rows) return error.DimensionMismatch;
    const r = try Mat(ElementOf(@TypeOf(a))).init(allocator, ca.rows, cb.cols);
    errdefer r.deinit(allocator);
    try mulInto(allocator, io, r.view(), ca, cb);
    return r;
}

/// `out = a b`: up to the inline limit in every dimension, `linalg.matrix.mul`'s
/// order (each column of `out` the columns of `a` scaled by `b`'s and
/// summed in order); else `gemv` for one column, `gemm` otherwise.
pub fn mulInto(allocator: Allocator, io: Io, out: anytype, a: anytype, b: anytype) MulError!void {
    const o = viewOf(out);
    const ca = constOf(a);
    const cb = constOf(b);
    if (ca.cols != cb.rows or o.rows != ca.rows or o.cols != cb.cols) return error.DimensionMismatch;
    return mulWith(pathOf(&.{ ca.rows, cb.cols, ca.cols }), allocator, io, o, ca, cb);
}

fn mulWith(path: Path, allocator: Allocator, io: Io, o: anytype, ca: ConstView(@TypeOf(o).Element), cb: ConstView(@TypeOf(o).Element)) (Allocator.Error || Io.Cancelable)!void {
    const T = @TypeOf(o).Element;
    const m = ca.rows;
    const n = cb.cols;
    const k = ca.cols;
    if (path == .fixed) {
        for (0..n) |j| for (0..m) |i| {
            var acc = ops.zero(T);
            if (k > 0) {
                acc = eMul(ca.at(i, 0), cb.at(0, j));
                for (1..k) |p| acc = eMulAdd(ca.at(i, p), cb.at(p, j), acc);
            }
            o.data[i + j * o.ld] = acc;
        };
        return;
    }
    const pa: [*]T = @constCast(ca.data.ptr);
    const pb: [*]T = @constCast(cb.data.ptr);
    const one = ops.one(T);
    if (n == 1) {
        blas.gemv(allocator, io, T, .n, toIsize(m), toIsize(k), one, pa, toIsize(ca.ld), pb, 1, ops.zero(T), o.data.ptr, 1) catch |e| return narrow(Error, e);
        return;
    }
    blas.gemm(allocator, io, T, .n, .n, toIsize(m), toIsize(n), toIsize(k), one, pa, toIsize(ca.ld), pb, toIsize(cb.ld), ops.zero(T), o.data.ptr, toIsize(o.ld)) catch |e| return narrow(Error, e);
}

pub const TransposeError = error{DimensionMismatch} || Io.Cancelable;

/// `aᵀ` (MATLAB's `a.'`), owned by the caller.
pub fn transpose(allocator: Allocator, io: Io, a: anytype) (Allocator.Error || Io.Cancelable)!Mat(ElementOf(@TypeOf(a))) {
    return transposed(.t, allocator, io, a);
}

/// `aᴴ` (MATLAB's `a'`); `aᵀ` for real `T`.
pub fn ctranspose(allocator: Allocator, io: Io, a: anytype) (Allocator.Error || Io.Cancelable)!Mat(ElementOf(@TypeOf(a))) {
    return transposed(.c, allocator, io, a);
}

fn transposed(comptime trans: blas.Trans, allocator: Allocator, io: Io, a: anytype) (Allocator.Error || Io.Cancelable)!Mat(ElementOf(@TypeOf(a))) {
    const ca = constOf(a);
    const r = try Mat(ElementOf(@TypeOf(a))).init(allocator, ca.cols, ca.rows);
    errdefer r.deinit(allocator);
    transposeWith(trans, io, r.view(), ca) catch |e| return narrow(Io.Cancelable, e);
    return r;
}

/// `out = aᵀ`, by the BLAS's `omatcopy`.
pub fn transposeInto(io: Io, out: anytype, a: anytype) TransposeError!void {
    return transposeWith(.t, io, viewOf(out), constOf(a));
}

/// `out = aᴴ`, by the BLAS's `omatcopy`.
pub fn ctransposeInto(io: Io, out: anytype, a: anytype) TransposeError!void {
    return transposeWith(.c, io, viewOf(out), constOf(a));
}

fn transposeWith(comptime trans: blas.Trans, io: Io, o: anytype, a: anytype) TransposeError!void {
    const T = @TypeOf(a).Element;
    if (o.rows != a.cols or o.cols != a.rows) return error.DimensionMismatch;
    if (a.rows == 0 or a.cols == 0) return;
    const one = ops.one(T);
    const t: blas.Trans = if (comptime isComplex(T)) trans else .t;
    blas.omatcopy(io, T, t, toIsize(a.rows), toIsize(a.cols), one, a.data.ptr, toIsize(a.ld), o.data.ptr, toIsize(o.ld)) catch |e| return narrow(Io.Cancelable, e);
}

// ---- determinant, inverse and solve -----------------------------------------

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

fn detWith(path: Path, allocator: Allocator, io: Io, ca: anytype) (Allocator.Error || Io.Cancelable)!@TypeOf(ca).Element {
    const T = @TypeOf(ca).Element;
    const n = ca.rows;
    if (n == 0) return ops.one(T);
    if (path == .fixed) return fixed(T, "det", n, .{ca});
    const f = try copyOf(T, allocator, ca);
    defer f.deinit(allocator);
    const ipiv = try allocator.alloc(isize, n);
    defer allocator.free(ipiv);
    lapack.getrf(allocator, io, T, toIsize(n), toIsize(n), f.data.ptr, toIsize(n), ipiv.ptr, null) catch |e| switch (e) {
        error.Singular => return ops.zero(T),
        else => return narrow(Error, e),
    };
    var d = ops.one(T);
    var exp: i64 = 0;
    for (ipiv, 0..) |p, i| {
        d = eMul(d, f.at(i, i));
        if (p != i + 1) d = ops.neg(T, d);
        const x = std.math.frexp(if (comptime isComplex(T)) @max(@abs(d.re), @abs(d.im)) else d);
        exp += x.exponent;
        d = scaleBy(d, -x.exponent);
    }
    return scaleBy(d, @intCast(std.math.clamp(exp, -1 << 16, 1 << 16)));
}

/// `x 2^e`.
fn scaleBy(x: anytype, e: i32) @TypeOf(x) {
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

fn invWith(path: Path, allocator: Allocator, io: Io, o: anytype, a: ConstView(@TypeOf(o).Element)) (error{Singular} || Allocator.Error || Io.Cancelable)!void {
    const T = @TypeOf(o).Element;
    const n = a.rows;
    if (n == 0) return;
    if (path == .fixed) return fixed(error{Singular}!void, "inv", n, .{ o, a });
    const f = try copyOf(T, allocator, a);
    defer f.deinit(allocator);
    const ipiv = try allocator.alloc(isize, n);
    defer allocator.free(ipiv);
    lapack.getrf(allocator, io, T, toIsize(n), toIsize(n), f.data.ptr, toIsize(n), ipiv.ptr, null) catch |e| return narrow(error{Singular} || Error, e);
    lapack.getri(allocator, io, T, toIsize(n), f.data.ptr, toIsize(n), ipiv.ptr, null) catch |e| return narrow(error{Singular} || Error, e);
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

const Structure = enum { upper, lower, hermitian, general };

/// `solve`'s route past the inline limit, from `a`'s values alone (MATLAB's
/// order): triangular `trtrs`; Hermitian with a positive real diagonal
/// `potrf` and `potrs`, LU if that fails; else `gesv`.
fn structureOf(a: anytype) Structure {
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

fn solveWith(path: Path, allocator: Allocator, io: Io, o: anytype, a: ConstView(@TypeOf(o).Element), b: ConstView(@TypeOf(o).Element)) SolveError!void {
    const T = @TypeOf(o).Element;
    const n = a.cols;
    if (a.rows != n) return leastSquares(allocator, io, o, a, b);
    if (n == 0) return;
    if (path == .fixed) return fixed(error{Singular}!void, "solve", n, .{ o, a, b });
    const x = try copyOf(T, allocator, b);
    defer x.deinit(allocator);
    const nn = toIsize(n);
    const k = toIsize(b.cols);
    const ld = toIsize(@max(n, 1));
    switch (structureOf(a)) {
        .upper, .lower => |s| {
            const uplo: lapack.Uplo = if (s == .upper) .u else .l;
            lapack.trtrs(allocator, io, T, uplo, .n, .n, nn, k, a.data.ptr, toIsize(a.ld), x.data.ptr, ld, null) catch |e| return narrow(error{Singular} || Error, e);
        },
        .hermitian, .general => |s| {
            const f = try copyOf(T, allocator, a);
            defer f.deinit(allocator);
            if (s == .hermitian) cholesky: {
                lapack.potrf(allocator, io, T, .u, nn, f.data.ptr, ld, null) catch |e| switch (e) {
                    error.NotPositiveDefinite => {
                        copyInto(T, f.data, n, a);
                        break :cholesky;
                    },
                    else => return narrow(Error, e),
                };
                lapack.potrs(allocator, io, T, .u, nn, k, f.data.ptr, ld, x.data.ptr, ld) catch |e| return narrow(Error, e);
                return copyInto(T, o.data, o.ld, x.constView());
            }
            const ipiv = try allocator.alloc(isize, n);
            defer allocator.free(ipiv);
            lapack.gesv(allocator, io, T, nn, k, f.data.ptr, ld, ipiv.ptr, x.data.ptr, ld, null) catch |e| return narrow(error{Singular} || Error, e);
        },
    }
    copyInto(T, o.data, o.ld, x.constView());
}

/// `gelsy`'s least-norm least-squares solution, for rectangular `a`.
fn leastSquares(allocator: Allocator, io: Io, o: anytype, a: ConstView(@TypeOf(o).Element), b: ConstView(@TypeOf(o).Element)) (Allocator.Error || Io.Cancelable)!void {
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
    const jpvt = try allocator.alloc(isize, n);
    defer allocator.free(jpvt);
    @memset(jpvt, 0);
    const rcond = @as(F, @floatFromInt(ldx)) * std.math.floatEps(F);
    var got_rank: isize = undefined;
    lapack.gelsy(allocator, io, T, toIsize(m), toIsize(n), toIsize(b.cols), f.data.ptr, toIsize(m), x.ptr, toIsize(ldx), jpvt.ptr, rcond, &got_rank) catch |e| return narrow(Error, e);
    for (0..o.cols) |j| @memcpy(o.col(j), x[j * ldx ..][0..n]);
}

// ---- factorizations ---------------------------------------------------------

/// `lu`'s result, `P A = L U` for the `m` by `n` `A`, `getrf`'s factors.
pub fn Lu(comptime T: type) type {
    return struct {
        /// `L` below the diagonal (its unit diagonal implied), `U` on and
        /// above it.
        factors: Mat(T),
        /// Row `i` was interchanged with row `ipiv[i]`, for `i` from 0 up.
        ipiv: []usize,
        /// The first `i` with `U(i, i)` exactly 0; a singular matrix still
        /// factors completely (design.md, decision 5).
        zero_pivot: ?usize,

        const Self = @This();

        pub fn deinit(f: Self, allocator: Allocator) void {
            f.factors.deinit(allocator);
            allocator.free(f.ipiv);
        }

        /// `L`, `m` by `min(m, n)`, unit lower trapezoidal.
        pub fn l(f: Self, allocator: Allocator) Allocator.Error!Mat(T) {
            const k = @min(f.factors.rows, f.factors.cols);
            const r = try Mat(T).zeros(allocator, f.factors.rows, k);
            for (0..k) |j| {
                r.data[j + j * r.rows] = ops.one(T);
                for (j + 1..r.rows) |i| r.data[i + j * r.rows] = f.factors.at(i, j);
            }
            return r;
        }

        /// `U`, `min(m, n)` by `n`, upper trapezoidal.
        pub fn u(f: Self, allocator: Allocator) Allocator.Error!Mat(T) {
            const k = @min(f.factors.rows, f.factors.cols);
            const r = try Mat(T).zeros(allocator, k, f.factors.cols);
            for (0..r.cols) |j| for (0..@min(j + 1, k)) |i| {
                r.data[i + j * k] = f.factors.at(i, j);
            };
            return r;
        }
    };
}

/// `P A = L U` with partial pivoting (MATLAB's `lu`), by `getrf`.
pub fn lu(allocator: Allocator, io: Io, a: anytype) (Allocator.Error || Io.Cancelable)!Lu(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const k = @min(ca.rows, ca.cols);
    const f = try copyOf(T, allocator, ca);
    errdefer f.deinit(allocator);
    const ipiv = try allocator.alloc(usize, k);
    errdefer allocator.free(ipiv);
    const piv = try allocator.alloc(isize, k);
    defer allocator.free(piv);
    var zero: usize = 0;
    var zero_pivot: ?usize = null;
    lapack.getrf(allocator, io, T, toIsize(ca.rows), toIsize(ca.cols), f.data.ptr, toIsize(@max(ca.rows, 1)), piv.ptr, &zero) catch |e| switch (e) {
        error.Singular => zero_pivot = zero - 1,
        else => return narrow(Error, e),
    };
    for (ipiv, piv) |*d, p| d.* = @intCast(p - 1);
    return .{ .factors = f, .ipiv = ipiv, .zero_pivot = zero_pivot };
}

pub const CholError = error{ NotSquare, NotPositiveDefinite } || Allocator.Error || Io.Cancelable;

/// The upper `R` with `A = Rᴴ R` (MATLAB's `chol`), by `potrf` on `a`'s
/// upper triangle (the lower is not read); zeros below the diagonal.
pub fn chol(allocator: Allocator, io: Io, a: anytype) CholError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    const r = try copyOf(T, allocator, ca);
    errdefer r.deinit(allocator);
    lapack.potrf(allocator, io, T, .u, toIsize(n), r.data.ptr, toIsize(@max(n, 1)), null) catch |e| return narrow(error{NotPositiveDefinite} || Error, e);
    for (0..n) |j| @memset(r.data[j * n + j + 1 .. (j + 1) * n], ops.zero(T));
    return r;
}

/// `qr`'s shapes for the `m` by `n` `A`, `k = min(m, n)`: `.economy` `Q` `m`
/// by `k` and `R` `k` by `n`; `.full` `Q` `m` by `m` and `R` `m` by `n`.
pub const QrMode = enum { economy, full };

/// `qr`'s result.
pub fn Qr(comptime T: type) type {
    return struct {
        /// Orthonormal columns.
        q: Mat(T),
        /// Upper trapezoidal.
        r: Mat(T),

        pub fn deinit(f: @This(), allocator: Allocator) void {
            f.q.deinit(allocator);
            f.r.deinit(allocator);
        }
    };
}

/// `A = Q R` (MATLAB's `qr`, unpivoted), by `geqrf` and `orgqr` (`ungqr`);
/// `R`'s diagonal keeps LAPACK's signs.
pub fn qr(allocator: Allocator, io: Io, a: anytype, mode: QrMode) (Allocator.Error || Io.Cancelable)!Qr(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const m = ca.rows;
    const n = ca.cols;
    const k = @min(m, n);
    const qc = if (mode == .full) m else k;
    const f = try copyOf(T, allocator, ca);
    defer f.deinit(allocator);
    const tau = try allocator.alloc(T, k);
    defer allocator.free(tau);
    const ld = toIsize(@max(m, 1));
    lapack.geqrf(allocator, io, T, toIsize(m), toIsize(n), f.data.ptr, ld, tau.ptr) catch |e| return narrow(Error, e);
    const r = try Mat(T).zeros(allocator, qc, n);
    errdefer r.deinit(allocator);
    for (0..n) |j| for (0..@min(j + 1, k)) |i| {
        r.data[i + j * qc] = f.at(i, j);
    };
    const q = try Mat(T).zeros(allocator, m, qc);
    errdefer q.deinit(allocator);
    for (0..k) |j| @memcpy(q.data[j * m ..][0..m], f.data[j * m ..][0..m]);
    if (m > 0) lapack.orgqr(allocator, io, T, toIsize(m), toIsize(qc), toIsize(k), q.data.ptr, ld, tau.ptr) catch |e| return narrow(Error, e);
    return .{ .q = q, .r = r };
}

// ---- norms, SVD and eigenvalues ---------------------------------------------

/// `matrixNorm`'s norms: the largest column sum of `|a(i, j)|`, the largest
/// row sum, Frobenius, and the largest singular value.
pub const MatrixNorm = enum { one, inf, fro, two };

pub const MatrixNormError = error{NoConvergence} || Allocator.Error || Io.Cancelable;

/// A matrix norm (MATLAB's `norm(A, p)`); 0 when empty. `lange`, or `.two`
/// the first singular value from `gesvd`.
pub fn matrixNorm(allocator: Allocator, io: Io, a: anytype, which: MatrixNorm) MatrixNormError!Real(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const m = ca.rows;
    const n = ca.cols;
    if (m == 0 or n == 0) return 0;
    const kind: lapack.Norm = switch (which) {
        .one => .o,
        .inf => .i,
        .fro => .f,
        .two => {
            const s = try singularValues(allocator, io, ca);
            defer allocator.free(s);
            return s[0];
        },
    };
    return lapack.lange(allocator, io, T, kind, toIsize(m), toIsize(n), ca.data.ptr, toIsize(ca.ld)) catch |e| return narrow(Allocator.Error, e);
}

/// `svd`'s outputs for the `m` by `n` `A`, `k = min(m, n)`: `.values` `s`
/// alone; `.economy` `U` `m` by `k` and `V` `n` by `k`; `.full` `U` `m` by
/// `m` and `V` `n` by `n`.
pub const SvdMode = enum { values, economy, full };

/// `svd`'s result, `A = U Σ Vᴴ`.
pub fn Svd(comptime T: type) type {
    return struct {
        /// The `k` singular values, nonincreasing.
        s: []Real(T),
        /// Orthonormal columns; null for `.values`.
        u: ?Mat(T),
        v: ?Mat(T),

        pub fn deinit(f: @This(), allocator: Allocator) void {
            allocator.free(f.s);
            if (f.u) |u| u.deinit(allocator);
            if (f.v) |v| v.deinit(allocator);
        }
    };
}

/// `svd`'s, `pinv`'s, `rank`'s and `cond`'s: `NoConvergence` also for a NaN
/// in `A`.
pub const SvdError = error{NoConvergence} || Allocator.Error || Io.Cancelable;

/// `InvalidA` (a NaN) is `NoConvergence`.
fn svdError(e: anyerror) SvdError {
    return if (e == error.InvalidA) error.NoConvergence else narrow(SvdError, e);
}

/// `A = U Σ Vᴴ` (MATLAB's `svd`): the values alone by `gesvd` (QR
/// iteration, which keeps a graded matrix's small values), with vectors by
/// `gesdd` (divide and conquer), so the two may differ in the last bits.
/// Empty `A` with `.full`: `U` and `V` are identities.
pub fn svd(allocator: Allocator, io: Io, a: anytype, mode: SvdMode) SvdError!Svd(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    if (mode == .values) return .{ .s = try singularValues(allocator, io, ca), .u = null, .v = null };
    const m = ca.rows;
    const n = ca.cols;
    const k = @min(m, n);
    const uc = if (mode == .full) m else k;
    const vc = if (mode == .full) n else k;
    const s = try allocator.alloc(Real(T), k);
    errdefer allocator.free(s);
    const u = try Mat(T).zeros(allocator, m, uc);
    errdefer u.deinit(allocator);
    const v = try Mat(T).zeros(allocator, n, vc);
    errdefer v.deinit(allocator);
    if (k == 0) {
        for (0..uc) |i| u.data[i * (m + 1)] = ops.one(T);
        for (0..vc) |i| v.data[i * (n + 1)] = ops.one(T);
        return .{ .s = s, .u = u, .v = v };
    }
    const f = try copyOf(T, allocator, ca);
    defer f.deinit(allocator);
    const vt = try allocator.alloc(T, vc * n);
    defer allocator.free(vt);
    lapack.gesdd(allocator, io, T, if (mode == .full) .a else .s, toIsize(m), toIsize(n), f.data.ptr, toIsize(m), s.ptr, u.data.ptr, toIsize(m), vt.ptr, toIsize(vc)) catch |e| return svdError(e);
    for (0..n) |j| for (0..vc) |i| {
        v.data[j + i * n] = ops.conj(T, vt[i + j * vc]);
    };
    return .{ .s = s, .u = u, .v = v };
}

/// `a`'s `min(m, n)` singular values by `gesvd`, owned by the caller.
fn singularValues(allocator: Allocator, io: Io, a: anytype) SvdError![]Real(@TypeOf(a).Element) {
    const T = @TypeOf(a).Element;
    const m = a.rows;
    const n = a.cols;
    const s = try allocator.alloc(Real(T), @min(m, n));
    errdefer allocator.free(s);
    if (s.len == 0) return s;
    const f = try copyOf(T, allocator, a);
    defer f.deinit(allocator);
    var none: [1]T = undefined;
    lapack.gesvd(allocator, io, T, .n, .n, toIsize(m), toIsize(n), f.data.ptr, toIsize(m), s.ptr, &none, 1, &none, 1, null) catch |e| return svdError(e);
    return s;
}

/// The spacing of floats at `x` (MATLAB's `eps(x)`).
fn epsOf(x: anytype) @TypeOf(x) {
    const F = @TypeOf(x);
    if (x == 0) return std.math.floatTrueMin(F);
    const e = std.math.frexp(x).exponent;
    return @max(std.math.ldexp(std.math.floatEps(F), e - 1), std.math.floatTrueMin(F));
}

/// How many of the nonincreasing `s` exceed `tol`, by default
/// `max(m, n) eps(σ₁)` (MATLAB's `rank` and `pinv`).
fn rankOf(comptime F: type, s: []const F, m: usize, n: usize, tol: ?F) usize {
    if (s.len == 0) return 0;
    const t = tol orelse (@as(F, @floatFromInt(@max(m, n))) * epsOf(s[0]));
    var r: usize = 0;
    while (r < s.len and s[r] > t) r += 1;
    return r;
}

/// The Moore-Penrose pseudoinverse (MATLAB's `pinv`), `n` by `m`, owned by
/// the caller: `V Σ⁺ Uᴴ` over the singular values above `tol`, by default
/// `max(m, n) eps(σ₁)`.
pub fn pinv(allocator: Allocator, io: Io, a: anytype, tol: ?Real(ElementOf(@TypeOf(a)))) SvdError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const m = ca.rows;
    const n = ca.cols;
    const x = try Mat(T).zeros(allocator, n, m);
    errdefer x.deinit(allocator);
    const f = try svd(allocator, io, ca, .economy);
    defer f.deinit(allocator);
    const r = rankOf(Real(T), f.s, m, n, tol);
    if (r == 0) return x;
    const u = f.u.?;
    const v = f.v.?;
    for (f.s[0..r], 0..) |sj, j| for (v.data[j * n ..][0..n]) |*e| {
        e.* = if (comptime isComplex(T)) T.init(e.re / sj, e.im / sj) else e.* / sj;
    };
    const t: blas.Trans = if (comptime isComplex(T)) .c else .t;
    blas.gemm(allocator, io, T, .n, t, toIsize(n), toIsize(m), toIsize(r), ops.one(T), v.data.ptr, toIsize(n), u.data.ptr, toIsize(m), ops.zero(T), x.data.ptr, toIsize(n)) catch |e| return narrow(Error, e);
    return x;
}

/// The number of singular values above `tol`, by default `max(m, n)
/// eps(σ₁)` (MATLAB's `rank`); 0 when empty.
pub fn rank(allocator: Allocator, io: Io, a: anytype, tol: ?Real(ElementOf(@TypeOf(a)))) SvdError!usize {
    const ca = constOf(a);
    const s = try singularValues(allocator, io, ca);
    defer allocator.free(s);
    return rankOf(Real(ElementOf(@TypeOf(a))), s, ca.rows, ca.cols, tol);
}

/// The 2-norm condition number `σ₁ / σₖ` (MATLAB's `cond`); infinite when
/// `σₖ` is 0, 0 when empty.
pub fn cond(allocator: Allocator, io: Io, a: anytype) SvdError!Real(ElementOf(@TypeOf(a))) {
    const s = try singularValues(allocator, io, constOf(a));
    defer allocator.free(s);
    if (s.len == 0) return 0;
    const last = s[s.len - 1];
    return if (last == 0) std.math.inf(@TypeOf(last)) else s[0] / last;
}

/// `eig`'s result: Hermitian `A` (decided by its values alone, exactly)
/// has real eigenvalues, nondecreasing, and orthonormal eigenvectors;
/// otherwise they are complex, a real `A`'s conjugate pairs together with
/// the positive imaginary part first, and each eigenvector (`A v = λ v`) of
/// unit 2-norm with its largest entry real (design.md, decision 6).
pub fn Eig(comptime T: type) type {
    const C = std.math.Complex(Real(T));
    return union(enum) {
        hermitian: struct { values: []Real(T), vectors: ?Mat(T) },
        general: struct { values: []C, vectors: ?Mat(C) },

        pub fn deinit(e: @This(), allocator: Allocator) void {
            switch (e) {
                inline else => |r| {
                    allocator.free(r.values);
                    if (r.vectors) |v| v.deinit(allocator);
                },
            }
        }
    };
}

pub const EigError = error{ NotSquare, NoConvergence } || Allocator.Error || Io.Cancelable;

/// Eigenvalues, and with `want_vectors` eigenvectors (MATLAB's `eig`):
/// `syevd` (`heevd`) for Hermitian `A`, else `geev`. A NaN in `A` gives
/// `error.NoConvergence`.
pub fn eig(allocator: Allocator, io: Io, a: anytype, want_vectors: bool) EigError!Eig(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const C = std.math.Complex(Real(T));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    const ld = toIsize(@max(n, 1));
    const jobz: lapack.Jobz = if (want_vectors) .v else .n;
    const f = try copyOf(T, allocator, ca);
    if (isHermitian(ca)) {
        errdefer f.deinit(allocator);
        const w = try allocator.alloc(Real(T), n);
        errdefer allocator.free(w);
        lapack.syevd(allocator, io, T, jobz, .u, toIsize(n), f.data.ptr, ld, w.ptr) catch |e| return svdError(e);
        if (!want_vectors) f.deinit(allocator);
        return .{ .hermitian = .{ .values = w, .vectors = if (want_vectors) f else null } };
    }
    defer f.deinit(allocator);
    const w = try allocator.alloc(C, n);
    errdefer allocator.free(w);
    const vr: ?Mat(T) = if (want_vectors) try Mat(T).init(allocator, n, n) else null;
    var none: [1]T = undefined;
    const pr = if (vr) |x| x.data.ptr else &none;
    if (comptime isComplex(T)) {
        errdefer if (vr) |x| x.deinit(allocator);
        lapack.geev(allocator, io, T, .n, jobz, toIsize(n), f.data.ptr, ld, w.ptr, &none, 1, pr, ld, null) catch |e| return svdError(e);
        return .{ .general = .{ .values = w, .vectors = vr } };
    }
    defer if (vr) |x| x.deinit(allocator);
    const wr = try allocator.alloc(T, 2 * n);
    defer allocator.free(wr);
    lapack.geev(allocator, io, T, .n, jobz, toIsize(n), f.data.ptr, ld, .{ .re = wr.ptr, .im = wr[n..].ptr }, &none, 1, pr, ld, null) catch |e| return svdError(e);
    for (w, wr[0..n], wr[n..]) |*d, re, im| d.* = .init(re, im);
    const x = vr orelse return .{ .general = .{ .values = w, .vectors = null } };
    // `v(:, j) ± i v(:, j + 1)` for a pair, `wi(j) > 0` first
    const v = try Mat(C).init(allocator, n, n);
    var j: usize = 0;
    while (j < n) : (j += 1) {
        const c = x.data[j * n ..][0..n];
        if (wr[n + j] == 0) {
            for (v.data[j * n ..][0..n], c) |*d, re| d.* = .init(re, 0);
            continue;
        }
        const ci = x.data[(j + 1) * n ..][0..n];
        for (v.data[j * n ..][0..n], v.data[(j + 1) * n ..][0..n], c, ci) |*p, *q, re, im| {
            p.* = .init(re, im);
            q.* = .init(re, -im);
        }
        j += 1;
    }
    return .{ .general = .{ .values = w, .vectors = v } };
}

/// `a == aᴴ` exactly, the diagonal real.
fn isHermitian(a: anytype) bool {
    const T = @TypeOf(a).Element;
    for (0..a.rows) |j| {
        if (comptime isComplex(T)) if (a.at(j, j).im != 0) return false;
        for (j + 1..a.rows) |i| {
            const x = a.at(i, j);
            const y = ops.conj(T, a.at(j, i));
            const same = if (comptime isComplex(T)) x.re == y.re and x.im == y.im else x == y;
            if (!same) return false;
        }
    }
    return true;
}

// ---- plumbing ---------------------------------------------------------------

/// The inline limit `L` per dimension, `rstd_linalg_options.inline_limit`
/// (design.md, Dispatch).
const inline_limit = matrix.inline_limit;

/// The largest inline size built: `L`, and 16 in tests, which run each path
/// at every size 1 to 16 (design.md, Dispatch).
const fixed_max = if (sysinfo.is_test) @max(inline_limit, 16) else inline_limit;

/// The inline kernels or the BLAS and LAPACK; by size alone outside tests.
const Path = enum { fixed, blas };

/// `.fixed` when every dimension is at or below the inline limit.
fn pathOf(dims: []const usize) Path {
    for (dims) |d| if (d > inline_limit) return .blas;
    return .fixed;
}

/// The element type of an operand.
fn ElementOf(comptime X: type) type {
    switch (@typeInfo(X)) {
        .pointer => |p| switch (p.size) {
            .slice => return p.child,
            .one => return if (@typeInfo(p.child) == .array) @typeInfo(p.child).array.child else p.child.Element,
            else => {},
        },
        .@"struct" => if (@hasField(X, "data") and @typeInfo(@FieldType(X, "data")) == .pointer) return X.Element,
        else => {},
    }
    @compileError("expects a View, ConstView, Mat, pointer to a StaticMat or slice, found " ++ @typeName(X));
}

/// `x` read-only.
fn constOf(x: anytype) ConstView(ElementOf(@TypeOf(x))) {
    const T = ElementOf(@TypeOf(x));
    return switch (@typeInfo(@TypeOf(x))) {
        .pointer => |p| if (p.size == .slice)
            column(T, true, x)
        else if (@typeInfo(p.child) == .array) column(T, true, x[0..]) else x.constView(),
        else => x.constView(),
    };
}

/// `x` writable.
fn viewOf(x: anytype) View(ElementOf(@TypeOf(x))) {
    const X = @TypeOf(x);
    const T = ElementOf(X);
    if (X == View(T)) return x;
    if (X == ConstView(T)) @compileError("expects a writable operand, found a ConstView");
    return switch (@typeInfo(X)) {
        .pointer => |p| if (p.size == .slice)
            column(T, false, x)
        else if (@typeInfo(p.child) == .array) column(T, false, x[0..]) else x.view(),
        else => x.view(),
    };
}

fn column(comptime T: type, comptime read_only: bool, s: if (read_only) []const T else []T) if (read_only) ConstView(T) else View(T) {
    return .{ .data = s, .rows = s.len, .cols = 1, .ld = @max(s.len, 1) };
}

/// A vector's elements: `n` of them `inc` apart from `ptr`.
fn Strided(comptime T: type) type {
    return struct {
        ptr: [*]T,
        n: usize,
        inc: usize,

        fn at(s: @This(), i: usize) T {
            return s.ptr[i * s.inc];
        }

        /// The first `c` elements as an array, for the fixed-size forms.
        fn gather(s: @This(), comptime c: usize) [c]T {
            var r: [c]T = undefined;
            for (&r, 0..) |*x, i| x.* = s.at(i);
            return r;
        }
    };
}

/// `linalg.vector`'s forms at length `c` and `linalg.matrix`'s at `c` by `c`,
/// one non-inline function per size so that each gets its own comptime
/// branch quota.
fn Fixed(comptime c: usize) type {
    return struct {
        fn load(a: anytype) matrix.StaticMat(c, c, @TypeOf(a).Element) {
            var m: matrix.StaticMat(c, c, @TypeOf(a).Element) = undefined;
            for (&m.data, 0..) |*mc, j| @memcpy(mc, a.col(j));
            return m;
        }

        fn det(a: anytype) @TypeOf(a).Element {
            return matrix.any_size.det(load(a));
        }

        fn inv(o: anytype, a: ConstView(@TypeOf(o).Element)) error{Singular}!void {
            const r = try matrix.any_size.inv(load(a));
            for (&r.data, 0..) |*rc, j| @memcpy(o.col(j), rc);
        }

        /// `matrix.solve`'s LU, factored once for all of `b`'s columns.
        fn solve(o: anytype, a: ConstView(@TypeOf(o).Element), b: ConstView(@TypeOf(o).Element)) error{Singular}!void {
            const T = @TypeOf(o).Element;
            var f = load(a).data;
            var perm: [c]usize = undefined;
            var odd: bool = undefined;
            try matrix.luInPlace(c, T, &f, &perm, &odd);
            for (0..b.cols) |j| {
                const x = matrix.luSolve(c, T, &f, &perm, b.col(j)[0..c].*);
                @memcpy(o.col(j), &x);
            }
        }

        fn dot(x: anytype, y: @TypeOf(x)) @TypeOf(x.ptr[0]) {
            return vector.dot(x.gather(c), y.gather(c));
        }

        fn length(x: anytype) Real(@TypeOf(x.ptr[0])) {
            return vector.length(x.gather(c));
        }

        fn sum(x: anytype) @TypeOf(x.ptr[0]) {
            return vector.sum(x.gather(c));
        }

        fn extreme(comptime which: blas.Extreme, x: anytype) @TypeOf(x.ptr[0]) {
            return if (which == .max) vector.max(x.gather(c)) else vector.min(x.gather(c));
        }
    };
}

/// `Fixed(n).name(args)`, `n` from 1 to `fixed_max`.
inline fn fixed(comptime R: type, comptime name: []const u8, n: usize, args: anytype) R {
    switch (n) {
        inline 1...fixed_max => |c| return @call(.auto, @field(Fixed(c), name), args),
        else => unreachable,
    }
}

/// `v` as a vector: one column, or one row (stride `ld`).
fn strided(v: anytype) error{DimensionMismatch}!Strided(@TypeOf(v).Element) {
    const p: [*]@TypeOf(v).Element = @constCast(v.data.ptr);
    if (v.cols == 1 or v.rows == 0 or v.cols == 0) return .{ .ptr = p, .n = if (v.cols == 0) 0 else v.rows, .inc = 1 };
    if (v.rows == 1) return .{ .ptr = p, .n = v.cols, .inc = v.ld };
    return error.DimensionMismatch;
}

fn toIsize(n: usize) isize {
    return @intCast(n);
}

/// What the BLAS and LAPACK calls here pass on (`ops.narrow`): the arguments
/// come from views, so their argument errors are bugs.
const Error = Allocator.Error || Io.Cancelable;
const narrow = ops.narrow;

/// `a` into `dst`, column stride `ld`.
fn copyInto(comptime T: type, dst: []T, ld: usize, a: ConstView(T)) void {
    for (0..a.cols) |j| @memcpy(dst[j * ld ..][0..a.rows], a.col(j));
}

/// A packed copy of `a`, for LAPACK to overwrite.
fn copyOf(comptime T: type, allocator: Allocator, a: ConstView(T)) Allocator.Error!Mat(T) {
    const r = try Mat(T).init(allocator, a.rows, a.cols);
    copyInto(T, r.data, @max(a.rows, 1), a);
    return r;
}

// ---- tests ------------------------------------------------------------------

const zigstd = @import("../zigstd.zig");
const bitwise = @import("blas/_test/bitwise.zig");
const testing = std.testing;
const sameBits = bitwise.sameBits;
const C64 = std.math.Complex(f64);
const test_io = Io.failing;

fn randomElem(comptime T: type, rnd: std.Random) T {
    if (comptime isComplex(T)) return .init(rnd.float(f64) * 2 - 1, rnd.float(f64) * 2 - 1);
    return std.math.ldexp(rnd.float(T) * 2 - 1, rnd.intRangeAtMost(i32, -4, 4));
}

fn fill(comptime T: type, rnd: std.Random, s: []T) void {
    for (s) |*x| x.* = randomElem(T, rnd);
}

/// `mulInto` on a `Mat` against `linalg.matrix.mul` on `StaticMat`s of the
/// same numbers, `R` by `K` times `K` by `C`.
fn expectStaticBits(comptime T: type, comptime r: usize, comptime k: usize, comptime c: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    for (0..20) |_| {
        var a: matrix.StaticMat(r, k, T) = undefined;
        var b: matrix.StaticMat(k, c, T) = undefined;
        fill(T, rnd, a.view().data);
        fill(T, rnd, b.view().data);
        const want = matrix.mul(a, b);
        const got = try mul(al, test_io, &a, &b);
        defer got.deinit(al);
        var w = want;
        try testing.expect(sameBits(T, w.view().data, got.data));
        // a column: `M v`
        const v: [k]T = b.data[0];
        var col: [r]T = undefined;
        try mulInto(al, test_io, &col, &a, &v);
        try testing.expect(sameBits(T, &col, &@as([r]T, w.data[0])));
    }
}

test "up to the inline limit, a Mat multiplies to a StaticMat's bits" {
    var prng: std.Random.DefaultPrng = .init(0x7e57);
    const shapes = .{ .{ 1, 1, 1 }, .{ 2, 3, 4 }, .{ 4, 4, 4 }, .{ 3, 8, 5 }, .{ 8, 8, 8 } };
    inline for (.{ f32, f64, C64 }) |T| inline for (shapes) |s| try expectStaticBits(T, s[0], s[1], s[2], prng.random());
}

/// `Σ_p a(i, p) b(p, j)` in `f64` (or its complex), against `got` within
/// `k` roundings of `T`.
fn expectNearProduct(comptime T: type, a: ConstView(T), b: ConstView(T), got: ConstView(T)) !void {
    const eps: f64 = std.math.floatEps(Real(T));
    for (0..got.cols) |j| for (0..got.rows) |i| {
        var re: f64 = 0;
        var im: f64 = 0;
        var mag: f64 = 0;
        for (0..a.cols) |p| {
            const x = a.at(i, p);
            const y = b.at(p, j);
            if (comptime isComplex(T)) {
                re += x.re * y.re - x.im * y.im;
                im += x.re * y.im + x.im * y.re;
                mag += (@abs(x.re) + @abs(x.im)) * (@abs(y.re) + @abs(y.im));
            } else {
                re += @as(f64, x) * y;
                mag += @abs(@as(f64, x) * y);
            }
        }
        const g = got.at(i, j);
        const tol = 2 * @as(f64, @floatFromInt(a.cols + 1)) * eps * mag;
        if (comptime isComplex(T)) {
            try testing.expect(@abs(g.re - re) <= tol and @abs(g.im - im) <= tol);
        } else try testing.expect(@abs(@as(f64, g) - re) <= tol);
    };
}

fn expectBlasProduct(comptime T: type, m: usize, n: usize, k: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    const a = try Mat(T).init(al, m, k);
    defer a.deinit(al);
    const b = try Mat(T).init(al, k, n);
    defer b.deinit(al);
    fill(T, rnd, a.data);
    fill(T, rnd, b.data);
    // undefined storage may hold NaNs: `beta = 0` must not read them
    const out = try Mat(T).init(al, m, n);
    defer out.deinit(al);
    const nan = std.math.nan(Real(T));
    @memset(out.data, if (comptime isComplex(T)) T.init(nan, nan) else nan);
    try mulInto(al, test_io, out, a, b);
    try expectNearProduct(T, a.constView(), b.constView(), out.constView());
}

test "past the inline limit, mul is the BLAS's gemm or gemv" {
    var prng: std.Random.DefaultPrng = .init(0xb1a5);
    inline for (.{ f32, f64, C64 }) |T| {
        try expectBlasProduct(T, 9, 7, 5, prng.random());
        try expectBlasProduct(T, 3, 4, 20, prng.random());
        try expectBlasProduct(T, 40, 33, 17, prng.random());
        try expectBlasProduct(T, 30, 1, 12, prng.random());
    }
}

fn expectVectorOps(comptime T: type, n: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    // `x` is a row of a 3-row matrix, stride 3; `y` a column
    const xs = try al.alloc(T, 3 * n);
    defer al.free(xs);
    fill(T, rnd, xs);
    const xv = try ConstView(T).init(xs, 1, n, 3);
    const y = try al.alloc(T, n);
    defer al.free(y);
    fill(T, rnd, y);
    var d_re: f64 = 0;
    var d_im: f64 = 0;
    var s_re: f64 = 0;
    var s_im: f64 = 0;
    var sq: f64 = 0;
    var mag: f64 = 0;
    for (0..n) |i| {
        const x = xv.at(0, i);
        if (comptime isComplex(T)) {
            d_re += x.re * y[i].re + x.im * y[i].im;
            d_im += x.re * y[i].im - x.im * y[i].re;
            s_re += x.re;
            s_im += x.im;
            sq += x.re * x.re + x.im * x.im;
            mag += (@abs(x.re) + @abs(x.im)) * (@abs(y[i].re) + @abs(y[i].im));
        } else {
            d_re += @as(f64, x) * y[i];
            s_re += x;
            sq += @as(f64, x) * x;
            mag += @abs(@as(f64, x) * y[i]);
        }
    }
    const eps: f64 = std.math.floatEps(Real(T));
    const tol = 2 * @as(f64, @floatFromInt(n + 1)) * eps;
    const d = try dot(test_io, xv, y);
    const s = try sum(test_io, xv);
    const nrm: f64 = try norm(test_io, xv);
    try testing.expect(@abs(nrm - @sqrt(sq)) <= tol * @sqrt(sq));
    if (comptime isComplex(T)) {
        try testing.expect(@abs(d.re - d_re) <= tol * mag and @abs(d.im - d_im) <= tol * mag);
        try testing.expect(@abs(s.re - s_re) <= tol * @as(f64, @floatFromInt(n)) and @abs(s.im - s_im) <= tol * @as(f64, @floatFromInt(n)));
    } else {
        try testing.expect(@abs(@as(f64, d) - d_re) <= tol * mag);
        try testing.expect(@abs(@as(f64, s) - s_re) <= tol * @as(f64, @floatFromInt(n)));
        var hi = xv.at(0, 0);
        var lo = hi;
        for (1..n) |i| {
            hi = @max(hi, xv.at(0, i));
            lo = @min(lo, xv.at(0, i));
        }
        try testing.expectEqual(hi, try max(test_io, xv));
        try testing.expectEqual(lo, try min(test_io, xv));
    }
}

test "dot, norm, sum, max and min, inline and on the BLAS, on rows and columns" {
    var prng: std.Random.DefaultPrng = .init(0xd07);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ 1, 3, 8, 9, 100 }) |n| try expectVectorOps(T, n, prng.random());
}

test "up to the inline limit, vector reductions keep linalg.vector's bits" {
    var prng: std.Random.DefaultPrng = .init(0x5a5e);
    var a: [6]f64 = undefined;
    var b: [6]f64 = undefined;
    fill(f64, prng.random(), &a);
    fill(f64, prng.random(), &b);
    try testing.expectEqual(vector.dot(a, b), try dot(test_io, &a, &b));
    try testing.expectEqual(vector.length(a), try norm(test_io, &a));
    try testing.expectEqual(vector.sum(a), try sum(test_io, &a));
    var z: [5]C64 = undefined;
    fill(C64, prng.random(), &z);
    try testing.expectEqual(vector.dot(z, z), try dot(test_io, &z, &z));
}

test "transposes, axpy and scale" {
    const al = testing.allocator;
    var a = matrix.fromRows(matrix.StaticMat(2, 3, f64), .{ .{ 1, 2, 3 }, .{ 4, 5, 6 } });
    const t = try transpose(al, test_io, &a);
    defer t.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 1, 2, 3, 4, 5, 6 }, t.data);
    var too_small: [4]f64 = undefined;
    try testing.expectError(error.DimensionMismatch, transposeInto(test_io, &too_small, &a));

    var z: matrix.StaticMat(1, 2, C64) = .{ .data = .{ .{.init(1, 2)}, .{.init(3, -4)} } };
    const zh = try ctranspose(al, test_io, &z);
    defer zh.deinit(al);
    try testing.expectEqual(C64.init(1, -2), zh.data[0]);
    try testing.expectEqual(C64.init(3, 4), zh.data[1]);

    var y = [_]f64{ 1, 1, 1 };
    try axpy(test_io, 2, &[_]f64{ 1, 2, 3 }, &y);
    try testing.expectEqualSlices(f64, &.{ 3, 5, 7 }, &y);
    try testing.expectError(error.DimensionMismatch, axpy(test_io, 2, &[_]f64{ 1, 2 }, &y));
    scale(-1, &a);
    try testing.expectEqual(@as(f64, -6), a.data[2][1]);
    try testing.expectEqual(matrix.scale(matrix.fromCols(matrix.M2f32, .{ .{ 0.1, 0.2 }, .{ 0.3, 0.4 } }), 0.7).data[1][0], blk: {
        var s = matrix.fromCols(matrix.M2f32, .{ .{ 0.1, 0.2 }, .{ 0.3, 0.4 } });
        scale(0.7, &s);
        break :blk s.data[1][0];
    });
}

test "shape errors" {
    const al = testing.allocator;
    var a: matrix.StaticMat(2, 3, f64) = .{ .data = @splat(@splat(1)) };
    try testing.expectError(error.DimensionMismatch, dot(test_io, &a, &a));
    try testing.expectError(error.DimensionMismatch, dot(test_io, &[_]f64{1}, &[_]f64{ 1, 2 }));
    try testing.expectError(error.DimensionMismatch, mul(al, test_io, &a, &a));
    var out: matrix.M3f64 = undefined;
    try testing.expectError(error.DimensionMismatch, mulInto(al, test_io, &out, &a, &a));
    try testing.expectEqual(@as(f64, 0), try sum(test_io, &[_]f64{}));
    try testing.expectError(error.NotSquare, det(al, test_io, &a));
    try testing.expectError(error.NotSquare, inv(al, test_io, &a));
    try testing.expectError(error.NotSquare, chol(al, test_io, &a));
    try testing.expectError(error.DimensionMismatch, solve(al, test_io, &out, &a));
    try testing.expectError(error.DimensionMismatch, invInto(al, test_io, &a, &out));
    var x: [2]f64 = undefined;
    try testing.expectError(error.DimensionMismatch, solveInto(al, test_io, &x, &out, &[_]f64{ 1, 2, 3 }));
}

/// `det`, `inv` and `solve` on a `Mat` against `linalg.matrix`'s on the
/// `StaticMat` of the same numbers.
fn expectFixedBits(comptime T: type, comptime c: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    for (0..10) |_| {
        var a: matrix.StaticMat(c, c, T) = undefined;
        var b: matrix.StaticMat(c, 3, T) = undefined;
        fill(T, rnd, a.view().data);
        fill(T, rnd, b.view().data);
        const m = try copyOf(T, al, a.constView());
        defer m.deinit(al);
        try testing.expect(sameBits(T, &[_]T{matrix.det(a)}, &[_]T{try det(al, test_io, m)}));
        var want_inv = try matrix.inv(a);
        const got_inv = try inv(al, test_io, m);
        defer got_inv.deinit(al);
        try testing.expect(sameBits(T, want_inv.view().data, got_inv.data));
        var want_x = try matrix.solve(a, b);
        const got_x = try solve(al, test_io, m, &b);
        defer got_x.deinit(al);
        try testing.expect(sameBits(T, want_x.view().data, got_x.data));
    }
}

test "up to the inline limit, det, inv and solve keep a StaticMat's bits" {
    var prng: std.Random.DefaultPrng = .init(0x50f7);
    inline for (.{ f32, f64, C64 }) |T| inline for (.{ 1, 2, 3, 4, 6, 8 }) |c| try expectFixedBits(T, c, prng.random());
}

fn absOf(x: anytype) f64 {
    return if (comptime isComplex(@TypeOf(x))) @abs(x.re) + @abs(x.im) else @abs(@as(f64, x));
}

fn maxAbs(comptime T: type, a: ConstView(T)) f64 {
    var r: f64 = 0;
    for (0..a.cols) |j| for (a.col(j)) |x| {
        r = @max(r, absOf(x));
    };
    return r;
}

/// `max |a(i, j) - b(i, j)|`, `|re| + |im|` for complex.
fn maxDiff(comptime T: type, a: ConstView(T), b: ConstView(T)) f64 {
    var r: f64 = 0;
    for (0..a.cols) |j| for (0..a.rows) |i| {
        const x = a.at(i, j);
        const y = b.at(i, j);
        const d = if (comptime isComplex(T)) @abs(x.re - y.re) + @abs(x.im - y.im) else @abs(@as(f64, x) - y);
        r = @max(r, d);
    };
    return r;
}

/// 100 `n` roundings of `T`.
fn tolOf(comptime T: type, n: usize) f64 {
    return 100 * @as(f64, @floatFromInt(@max(n, 1))) * std.math.floatEps(Real(T));
}

fn randomMat(comptime T: type, rnd: std.Random, m: usize, n: usize) !Mat(T) {
    const r = try Mat(T).init(testing.allocator, m, n);
    fill(T, rnd, r.data);
    return r;
}

/// `a(i, i) += s`.
fn shiftDiag(comptime T: type, a: Mat(T), s: f64) void {
    for (0..@min(a.rows, a.cols)) |i| {
        const d = &a.data[i + i * a.rows];
        d.* = if (comptime isComplex(T)) T.init(d.re + s, d.im) else d.* + @as(T, @floatCast(s));
    }
}

/// `Mᴴ M + n I`, exactly Hermitian.
fn hermitianPd(comptime T: type, rnd: std.Random, n: usize) !Mat(T) {
    const al = testing.allocator;
    const m = try randomMat(T, rnd, n, n);
    defer m.deinit(al);
    const mh = try ctranspose(al, test_io, m);
    defer mh.deinit(al);
    const a = try mul(al, test_io, mh, m);
    for (0..n) |j| {
        for (j + 1..n) |i| a.data[i + j * n] = ops.conj(T, a.data[j + i * n]);
        if (comptime isComplex(T)) a.data[j + j * n].im = 0;
    }
    shiftDiag(T, a, @floatFromInt(n));
    return a;
}

/// `‖a x - b‖ <= tol (n ‖a‖ ‖x‖ + ‖b‖)`, max norms.
fn expectSolves(comptime T: type, a: ConstView(T), x: ConstView(T), b: ConstView(T)) !void {
    const al = testing.allocator;
    const ax = try mul(al, test_io, a, x);
    defer ax.deinit(al);
    const n: f64 = @floatFromInt(a.cols);
    const bound = tolOf(T, a.cols) * (n * maxAbs(T, a) * maxAbs(T, x) + maxAbs(T, b));
    try testing.expect(maxDiff(T, ax.constView(), b) <= bound);
}

fn expectSolveRoutes(comptime T: type, n: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    const b = try randomMat(T, rnd, n, 3);
    defer b.deinit(al);
    const general = try randomMat(T, rnd, n, n);
    defer general.deinit(al);
    shiftDiag(T, general, 2);
    const upper = try randomMat(T, rnd, n, n);
    defer upper.deinit(al);
    const lower = try randomMat(T, rnd, n, n);
    defer lower.deinit(al);
    for (0..n) |j| for (0..n) |i| {
        if (i > j) upper.data[i + j * n] = ops.zero(T);
        if (i < j) lower.data[i + j * n] = ops.zero(T);
    };
    shiftDiag(T, upper, 4);
    shiftDiag(T, lower, 4);
    const pd = try hermitianPd(T, rnd, n);
    defer pd.deinit(al);
    // Hermitian, positive diagonal, indefinite: `potrf` fails, LU solves
    const indefinite = try Mat(T).zeros(al, n, n);
    defer indefinite.deinit(al);
    shiftDiag(T, indefinite, 1);
    indefinite.data[1] = if (comptime isComplex(T)) T.init(2, 0) else 2;
    indefinite.data[n] = indefinite.data[1];
    try testing.expectError(error.NotPositiveDefinite, chol(al, test_io, indefinite));
    for ([_]Mat(T){ general, upper, lower, pd, indefinite }) |a| {
        const x = try solve(al, test_io, a, b);
        defer x.deinit(al);
        try expectSolves(T, a.constView(), x.constView(), b.constView());
        // one right-hand side: `trsv` and `gesv`'s single column
        const x1 = try Mat(T).init(al, n, 1);
        defer x1.deinit(al);
        try solveInto(al, test_io, x1, a, b.constView().col(0));
        try expectSolves(T, a.constView(), x1.constView(), (try b.constView().sub(0, 0, n, 1)));
    }
}

test "past the inline limit, solve routes by structure, each to a solution" {
    var prng: std.Random.DefaultPrng = .init(0x501e);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ 9, 30 }) |n| try expectSolveRoutes(T, n, prng.random());
    // a small square with many right-hand sides is past the limit too
    const al = testing.allocator;
    const a = try randomMat(f64, prng.random(), 3, 3);
    defer a.deinit(al);
    shiftDiag(f64, a, 2);
    const b = try randomMat(f64, prng.random(), 3, 20);
    defer b.deinit(al);
    const x = try solve(al, test_io, a, b);
    defer x.deinit(al);
    try expectSolves(f64, a.constView(), x.constView(), b.constView());
}

test "solve and inv: error.Singular leaves out untouched" {
    const al = testing.allocator;
    inline for (.{ 3, 10 }) |n| {
        // ones: Hermitian, not positive definite, LU's second pivot exactly 0
        const ones = try Mat(f64).init(al, n, n);
        defer ones.deinit(al);
        @memset(ones.data, 1);
        // zero: triangular, `trtrs`
        const zero = try Mat(f64).zeros(al, n, n);
        defer zero.deinit(al);
        for ([_]Mat(f64){ ones, zero }) |a| {
            var out: [n * n]f64 = @splat(7);
            const v = try View(f64).init(&out, n, n, n);
            try testing.expectError(error.Singular, solveInto(al, test_io, v, a, a));
            try testing.expectError(error.Singular, invInto(al, test_io, v, a));
            try testing.expect(std.mem.allEqual(f64, &out, 7));
            try testing.expectEqual(@as(f64, 0), try det(al, test_io, a));
        }
    }
}

fn expectLeastSquares(comptime T: type, m: usize, n: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    const a = try randomMat(T, rnd, m, n);
    defer a.deinit(al);
    const x0 = try randomMat(T, rnd, n, 2);
    defer x0.deinit(al);
    const b = try mul(al, test_io, a, x0);
    defer b.deinit(al);
    const x = try solve(al, test_io, a, b);
    defer x.deinit(al);
    try expectSolves(T, a.constView(), x.constView(), b.constView());
    if (m > n) {
        // full column rank: the only solution
        try testing.expect(maxDiff(T, x.constView(), x0.constView()) <= tolOf(T, m) * 10 * maxAbs(T, x0.constView()));
        return;
    }
    // least norm: `x = aᴴ (a aᴴ)⁻¹ b`
    const ah = try ctranspose(al, test_io, a);
    defer ah.deinit(al);
    const g = try mul(al, test_io, a, ah);
    defer g.deinit(al);
    const y = try solve(al, test_io, g, b);
    defer y.deinit(al);
    const want = try mul(al, test_io, ah, y);
    defer want.deinit(al);
    try testing.expect(maxDiff(T, x.constView(), want.constView()) <= tolOf(T, n) * 10 * maxAbs(T, want.constView()));
}

test "rectangular solve is gelsy's least-norm least squares" {
    var prng: std.Random.DefaultPrng = .init(0x9e15);
    inline for (.{ f32, f64, C64 }) |T| {
        try expectLeastSquares(T, 15, 9, prng.random());
        try expectLeastSquares(T, 4, 12, prng.random());
        try expectLeastSquares(T, 3, 2, prng.random());
    }
    const al = testing.allocator;
    const empty = try Mat(f64).init(al, 0, 4);
    defer empty.deinit(al);
    const x = try solve(al, test_io, empty, try empty.constView().sub(0, 0, 0, 2));
    defer x.deinit(al);
    try testing.expect(x.rows == 4 and x.cols == 2 and std.mem.allEqual(f64, x.data, 0));
}

test "det past the inline limit: exact cases, scaled past the range, against inv" {
    const al = testing.allocator;
    // `diag(1e200, 1e200, -3e-100, 1, ...)` (design.md, M1's E1)
    const e1 = try Mat(f64).zeros(al, 10, 10);
    defer e1.deinit(al);
    shiftDiag(f64, e1, 1);
    e1.data[0] = 1e200;
    e1.data[11] = 1e200;
    e1.data[22] = -3e-100;
    const d = try det(al, test_io, e1);
    try testing.expect(@abs(d + 3e300) <= 1e-14 * 3e300);
    // `2 I` and a row swap, exact
    const two = try Mat(f64).zeros(al, 12, 12);
    defer two.deinit(al);
    shiftDiag(f64, two, 2);
    try testing.expectEqual(@as(f64, 4096), try det(al, test_io, two));
    const swap = try Mat(f64).zeros(al, 9, 9);
    defer swap.deinit(al);
    shiftDiag(f64, swap, 1);
    swap.data[0] = 0;
    swap.data[10] = 0;
    swap.data[1] = 1;
    swap.data[9] = 1;
    try testing.expectEqual(@as(f64, -1), try det(al, test_io, swap));
    var zi: [81]C64 = @splat(.init(0, 0));
    for (0..9) |i| zi[i * 10] = .init(0, 1);
    try testing.expectEqual(C64.init(0, 1), try det(al, test_io, try ConstView(C64).init(&zi, 9, 9, 9)));
    try testing.expectEqual(@as(f64, 1), try det(al, test_io, try ConstView(f64).init(&.{}, 0, 0, 1)));

    var prng: std.Random.DefaultPrng = .init(0xde7);
    inline for (.{ f64, C64 }) |T| {
        const a = try randomMat(T, prng.random(), 12, 12);
        defer a.deinit(al);
        shiftDiag(T, a, 1);
        const ai = try inv(al, test_io, a);
        defer ai.deinit(al);
        const p = eMul(try det(al, test_io, a), try det(al, test_io, ai));
        try testing.expect(absOf(if (comptime isComplex(T)) T.init(p.re - 1, p.im) else p - 1) <= 1e-10);
    }
}

test "inv past the inline limit: a a⁻¹ is I" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x1a7);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ 9, 25 }) |n| {
        const a = try randomMat(T, prng.random(), n, n);
        defer a.deinit(al);
        shiftDiag(T, a, 2);
        const ai = try inv(al, test_io, a);
        defer ai.deinit(al);
        const id = try Mat(T).zeros(al, n, n);
        defer id.deinit(al);
        shiftDiag(T, id, 1);
        try expectSolves(T, a.constView(), ai.constView(), id.constView());
    };
}

/// `P A = L U`, rows swapped by `ipiv` in order.
fn expectLu(comptime T: type, a: Mat(T)) !Lu(T) {
    const al = testing.allocator;
    const f = try lu(al, test_io, a);
    errdefer f.deinit(al);
    const l = try f.l(al);
    defer l.deinit(al);
    const u = try f.u(al);
    defer u.deinit(al);
    const p = try mul(al, test_io, l, u);
    defer p.deinit(al);
    const pa = try copyOf(T, al, a.constView());
    defer pa.deinit(al);
    for (f.ipiv, 0..) |r, i| {
        if (r == i) continue;
        for (0..pa.cols) |j| std.mem.swap(T, &pa.data[i + j * pa.rows], &pa.data[r + j * pa.rows]);
    }
    try testing.expect(maxDiff(T, pa.constView(), p.constView()) <= tolOf(T, a.cols) * maxAbs(T, a.constView()));
    return f;
}

test "lu: P A = L U, singular too" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x10);
    inline for (.{ f32, f64, C64 }) |T| for ([_][2]usize{ .{ 10, 7 }, .{ 7, 10 }, .{ 12, 12 }, .{ 3, 2 } }) |s| {
        const a = try randomMat(T, prng.random(), s[0], s[1]);
        defer a.deinit(al);
        const f = try expectLu(T, a);
        defer f.deinit(al);
        try testing.expectEqual(@as(?usize, null), f.zero_pivot);
    };
    // a zero first column: `U(0, 0)` is 0, and the rest still factors
    const a = try randomMat(f64, prng.random(), 9, 9);
    defer a.deinit(al);
    @memset(a.data[0..9], 0);
    const f = try expectLu(f64, a);
    defer f.deinit(al);
    try testing.expectEqual(@as(?usize, 0), f.zero_pivot);
}

test "chol: Rᴴ R is A, R upper" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0xc401);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ 3, 20 }) |n| {
        const a = try hermitianPd(T, prng.random(), n);
        defer a.deinit(al);
        const r = try chol(al, test_io, a);
        defer r.deinit(al);
        for (0..n) |j| for (j + 1..n) |i| try testing.expect(ops.isZero(T, r.at(i, j)));
        const rh = try ctranspose(al, test_io, r);
        defer rh.deinit(al);
        const p = try mul(al, test_io, rh, r);
        defer p.deinit(al);
        try testing.expect(maxDiff(T, p.constView(), a.constView()) <= tolOf(T, n) * maxAbs(T, a.constView()));
    };
}

fn expectQr(comptime T: type, a: Mat(T), mode: QrMode) !void {
    const al = testing.allocator;
    const f = try qr(al, test_io, a, mode);
    defer f.deinit(al);
    const m = a.rows;
    const qc = if (mode == .full) m else @min(m, a.cols);
    try testing.expect(f.q.rows == m and f.q.cols == qc and f.r.rows == qc and f.r.cols == a.cols);
    for (0..a.cols) |j| for (@min(j + 1, qc)..qc) |i| try testing.expect(ops.isZero(T, f.r.at(i, j)));
    const p = try mul(al, test_io, f.q, f.r);
    defer p.deinit(al);
    try testing.expect(maxDiff(T, p.constView(), a.constView()) <= tolOf(T, m) * @max(maxAbs(T, a.constView()), 1));
    const qh = try ctranspose(al, test_io, f.q);
    defer qh.deinit(al);
    const qhq = try mul(al, test_io, qh, f.q);
    defer qhq.deinit(al);
    const id = try Mat(T).zeros(al, qc, qc);
    defer id.deinit(al);
    shiftDiag(T, id, 1);
    try testing.expect(maxDiff(T, qhq.constView(), id.constView()) <= tolOf(T, m));
}

test "qr: Q R is A, Q orthonormal, both modes" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x9a);
    inline for (.{ f32, f64, C64 }) |T| for ([_][2]usize{ .{ 12, 7 }, .{ 7, 12 }, .{ 5, 5 }, .{ 5, 0 } }) |s| {
        const a = try randomMat(T, prng.random(), s[0], s[1]);
        defer a.deinit(al);
        try expectQr(T, a, .economy);
        try expectQr(T, a, .full);
    };
}

test "matrixNorm: one, inf, fro and two" {
    const al = testing.allocator;
    const a = matrix.fromRows(matrix.M2f64, .{ .{ 1, -2 }, .{ 3, 4 } });
    try testing.expectEqual(@as(f64, 6), try matrixNorm(al, test_io, &a, .one));
    try testing.expectEqual(@as(f64, 7), try matrixNorm(al, test_io, &a, .inf));
    try testing.expectApproxEqRel(@sqrt(30.0), try matrixNorm(al, test_io, &a, .fro), 1e-15);
    try testing.expectApproxEqRel(@sqrt(15 + @sqrt(125.0)), try matrixNorm(al, test_io, &a, .two), 1e-14);
    const z: matrix.StaticMat(1, 2, C64) = .{ .data = .{ .{.init(3, 4)}, .{.init(0, 0)} } };
    try testing.expectApproxEqRel(@as(f64, 5), try matrixNorm(al, test_io, &z, .two), 1e-15);
    try testing.expectEqual(@as(f64, 0), try matrixNorm(al, test_io, try ConstView(f64).init(&.{}, 0, 3, 1), .two));
}

/// `x(:, j) *= d[j]`, real `d`.
fn scaleCols(comptime T: type, x: Mat(T), d: []const Real(T)) void {
    for (d, 0..) |dj, j| for (x.data[j * x.rows ..][0..x.rows]) |*e| {
        e.* = if (comptime isComplex(T)) T.init(e.re * dj, e.im * dj) else e.* * dj;
    };
}

/// `qᴴ q` is `I`.
fn expectOrthonormal(comptime T: type, q: Mat(T)) !void {
    const al = testing.allocator;
    const qh = try ctranspose(al, test_io, q);
    defer qh.deinit(al);
    const qhq = try mul(al, test_io, qh, q);
    defer qhq.deinit(al);
    const id = try Mat(T).zeros(al, q.cols, q.cols);
    defer id.deinit(al);
    shiftDiag(T, id, 1);
    try testing.expect(maxDiff(T, qhq.constView(), id.constView()) <= tolOf(T, q.rows));
}

fn expectSvd(comptime T: type, a: Mat(T), mode: SvdMode) !void {
    const al = testing.allocator;
    const m = a.rows;
    const n = a.cols;
    const k = @min(m, n);
    const f = try svd(al, test_io, a, mode);
    defer f.deinit(al);
    try testing.expectEqual(k, f.s.len);
    for (f.s, 0..) |x, i| try testing.expect(x >= 0 and (i == 0 or x <= f.s[i - 1]));
    const values = try svd(al, test_io, a, .values);
    defer values.deinit(al);
    for (f.s, values.s) |x, y| try testing.expect(@abs(x - y) <= tolOf(T, @max(m, n)) * @max(values.s[0], 1));
    if (mode == .values) return testing.expect(f.u == null and f.v == null);
    const u = f.u.?;
    const v = f.v.?;
    const uc = if (mode == .full) m else k;
    const vc = if (mode == .full) n else k;
    try testing.expect(u.rows == m and u.cols == uc and v.rows == n and v.cols == vc);
    try expectOrthonormal(T, u);
    try expectOrthonormal(T, v);
    // `U(:, 1:k) Σ V(:, 1:k)ᴴ`
    const us = try copyOf(T, al, try u.constView().sub(0, 0, m, k));
    defer us.deinit(al);
    scaleCols(T, us, f.s);
    const vh = try ctranspose(al, test_io, try v.constView().sub(0, 0, n, k));
    defer vh.deinit(al);
    const p = try mul(al, test_io, us, vh);
    defer p.deinit(al);
    try testing.expect(maxDiff(T, p.constView(), a.constView()) <= tolOf(T, @max(m, n)) * @max(maxAbs(T, a.constView()), 1));
}

test "svd: U Σ Vᴴ is A, U and V orthonormal, every mode" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x5fd);
    inline for (.{ f32, f64, C64 }) |T| for ([_][2]usize{ .{ 12, 7 }, .{ 7, 12 }, .{ 5, 5 }, .{ 5, 0 }, .{ 0, 3 }, .{ 1, 1 } }) |s| {
        const a = try randomMat(T, prng.random(), s[0], s[1]);
        defer a.deinit(al);
        inline for (.{ .values, .economy, .full }) |mode| try expectSvd(T, a, mode);
    };
    // past `gesdd`'s divide and conquer threshold (128)
    const big = try randomMat(f64, prng.random(), 140, 130);
    defer big.deinit(al);
    try expectSvd(f64, big, .economy);
    // known values: `diag(3, -4)` padded with a zero row
    const d = matrix.fromRows(matrix.StaticMat(3, 2, f64), .{ .{ 3, 0 }, .{ 0, -4 }, .{ 0, 0 } });
    const f = try svd(al, test_io, &d, .full);
    defer f.deinit(al);
    try testing.expectApproxEqRel(@as(f64, 4), f.s[0], 1e-15);
    try testing.expectApproxEqRel(@as(f64, 3), f.s[1], 1e-15);
}

/// `a V = V Λ`, `V` orthonormal, `Λ` nondecreasing.
fn expectHermitianEig(comptime T: type, a: Mat(T)) !void {
    const al = testing.allocator;
    const n = a.rows;
    const e = try eig(al, test_io, a, true);
    defer e.deinit(al);
    const h = e.hermitian;
    const v = h.vectors.?;
    for (h.values[1..], h.values[0 .. n - 1]) |x, y| try testing.expect(x >= y);
    try expectOrthonormal(T, v);
    const av = try mul(al, test_io, a, v);
    defer av.deinit(al);
    const vl = try copyOf(T, al, v.constView());
    defer vl.deinit(al);
    scaleCols(T, vl, h.values);
    const bound = tolOf(T, n) * @as(f64, @floatFromInt(n)) * maxAbs(T, a.constView());
    try testing.expect(maxDiff(T, av.constView(), vl.constView()) <= bound);
    const values = try eig(al, test_io, a, false);
    defer values.deinit(al);
    try testing.expect(values.hermitian.vectors == null);
    for (h.values, values.hermitian.values) |x, y| try testing.expect(@abs(x - y) <= bound);
}

/// `a v = λ v` per column, each `v` of unit norm.
fn expectGeneralEig(comptime T: type, a: Mat(T)) !void {
    const C = std.math.Complex(Real(T));
    const al = testing.allocator;
    const n = a.rows;
    const e = try eig(al, test_io, a, true);
    defer e.deinit(al);
    const g = e.general;
    const v = g.vectors.?;
    const ac = try Mat(C).init(al, n, n);
    defer ac.deinit(al);
    for (ac.data, a.data) |*d, x| d.* = if (comptime isComplex(T)) x else C.init(x, 0);
    const av = try mul(al, test_io, ac, v);
    defer av.deinit(al);
    const vl = try copyOf(C, al, v.constView());
    defer vl.deinit(al);
    for (g.values, 0..) |l, j| for (vl.data[j * n ..][0..n]) |*x| {
        x.* = x.mul(l);
    };
    const bound = tolOf(T, n) * @as(f64, @floatFromInt(n)) * maxAbs(T, a.constView());
    try testing.expect(maxDiff(C, av.constView(), vl.constView()) <= bound);
    for (0..n) |j| {
        var s: f64 = 0;
        for (v.constView().col(j)) |x| s += x.re * x.re + x.im * x.im;
        try testing.expect(@abs(s - 1) <= tolOf(T, n));
    }
    const values = try eig(al, test_io, a, false);
    defer values.deinit(al);
    try testing.expect(values.general.vectors == null);
    for (g.values, values.general.values) |x, y| try testing.expect(x.re == y.re and x.im == y.im);
}

test "eig: Hermitian input real, the rest complex" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0xe16);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ 1, 5, 20 }) |n| {
        const h = try hermitianPd(T, prng.random(), n);
        defer h.deinit(al);
        try expectHermitianEig(T, h);
        const a = try randomMat(T, prng.random(), n, n);
        defer a.deinit(al);
        // exactly Hermitian at random is impossible past 1 by 1
        if (n > 1) try expectGeneralEig(T, a);
    };
    // a rotation: `±i`, the positive imaginary part first, vectors conjugate
    const r = matrix.fromRows(matrix.M2f64, .{ .{ 0, -1 }, .{ 1, 0 } });
    const e = try eig(al, test_io, &r, true);
    defer e.deinit(al);
    const w = e.general.values;
    try testing.expect(@abs(w[0].re) + @abs(w[0].im - 1) <= 1e-15 and @abs(w[1].re) + @abs(w[1].im + 1) <= 1e-15);
    const v = e.general.vectors.?;
    for (v.constView().col(0), v.constView().col(1)) |x, y| try testing.expect(x.re == y.re and x.im == -y.im);
    // symmetric: real, ascending
    const s = matrix.fromRows(matrix.M2f64, .{ .{ 2, 1 }, .{ 1, 2 } });
    const es = try eig(al, test_io, &s, false);
    defer es.deinit(al);
    try testing.expectApproxEqRel(@as(f64, 1), es.hermitian.values[0], 1e-15);
    try testing.expectApproxEqRel(@as(f64, 3), es.hermitian.values[1], 1e-15);
    // empty is Hermitian; not square fails
    const z = try eig(al, test_io, try ConstView(f64).init(&.{}, 0, 0, 1), true);
    defer z.deinit(al);
    try testing.expect(z.hermitian.values.len == 0 and z.hermitian.vectors.?.data.len == 0);
    try testing.expectError(error.NotSquare, eig(al, test_io, try ConstView(f64).init(&.{ 1, 2 }, 2, 1, 2), false));
}

test "a NaN gives error.NoConvergence" {
    const al = testing.allocator;
    const nan = std.math.nan(f64);
    const g = matrix.fromRows(matrix.M2f64, .{ .{ 1, nan }, .{ 3, 4 } });
    const s = matrix.fromRows(matrix.M2f64, .{ .{ nan, 1 }, .{ 1, 4 } });
    inline for (.{ &g, &s }) |a| {
        try testing.expectError(error.NoConvergence, eig(al, test_io, a, true));
        try testing.expectError(error.NoConvergence, svd(al, test_io, a, .values));
        try testing.expectError(error.NoConvergence, svd(al, test_io, a, .economy));
        try testing.expectError(error.NoConvergence, pinv(al, test_io, a, null));
        try testing.expectError(error.NoConvergence, rank(al, test_io, a, null));
        try testing.expectError(error.NoConvergence, cond(al, test_io, a));
    }
}

/// The four Moore-Penrose conditions: `a x a = a`, `x a x = x`, `a x` and
/// `x a` Hermitian.
fn expectPinv(comptime T: type, a: Mat(T)) !void {
    const al = testing.allocator;
    const x = try pinv(al, test_io, a, null);
    defer x.deinit(al);
    try testing.expect(x.rows == a.cols and x.cols == a.rows);
    const ax = try mul(al, test_io, a, x);
    defer ax.deinit(al);
    const xa = try mul(al, test_io, x, a);
    defer xa.deinit(al);
    const axa = try mul(al, test_io, ax, a);
    defer axa.deinit(al);
    const xax = try mul(al, test_io, xa, x);
    defer xax.deinit(al);
    const axh = try ctranspose(al, test_io, ax);
    defer axh.deinit(al);
    const xah = try ctranspose(al, test_io, xa);
    defer xah.deinit(al);
    const c = try cond(al, test_io, a);
    const t = tolOf(T, @max(a.rows, a.cols)) * @max(1, if (std.math.isInf(c)) 1 else c);
    const na = @max(maxAbs(T, a.constView()), 1);
    const nx = @max(maxAbs(T, x.constView()), 1);
    try testing.expect(maxDiff(T, axa.constView(), a.constView()) <= t * na);
    try testing.expect(maxDiff(T, xax.constView(), x.constView()) <= t * nx);
    try testing.expect(maxDiff(T, axh.constView(), ax.constView()) <= t);
    try testing.expect(maxDiff(T, xah.constView(), xa.constView()) <= t);
}

test "pinv, rank and cond" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x919);
    inline for (.{ f32, f64, C64 }) |T| {
        for ([_][2]usize{ .{ 12, 7 }, .{ 7, 12 }, .{ 9, 9 } }) |s| {
            const a = try randomMat(T, prng.random(), s[0], s[1]);
            defer a.deinit(al);
            try expectPinv(T, a);
            try testing.expectEqual(@min(s[0], s[1]), try rank(al, test_io, a, null));
        }
        // rank 3: `x yᴴ`, 10 by 8
        const x = try randomMat(T, prng.random(), 10, 3);
        defer x.deinit(al);
        const y = try randomMat(T, prng.random(), 3, 8);
        defer y.deinit(al);
        const a = try mul(al, test_io, x, y);
        defer a.deinit(al);
        try testing.expectEqual(@as(usize, 3), try rank(al, test_io, a, null));
        try expectPinv(T, a);
        try testing.expect(try cond(al, test_io, a) > 1e-3 / std.math.floatEps(Real(T)));
    }
    // `rank`'s tolerance, and `cond` the ratio of the extremes
    const d = matrix.fromRows(matrix.M3f64, .{ .{ 4, 0, 0 }, .{ 0, 1e-3, 0 }, .{ 0, 0, -1e-9 } });
    try testing.expectEqual(@as(usize, 3), try rank(al, test_io, &d, null));
    try testing.expectEqual(@as(usize, 2), try rank(al, test_io, &d, 1e-6));
    try testing.expectApproxEqRel(@as(f64, 4e9), try cond(al, test_io, &d), 1e-14);
    // `pinv(diag(4, 1e-3, -1e-9), 1e-6)` drops the last value
    const p = try pinv(al, test_io, &d, 1e-6);
    defer p.deinit(al);
    try testing.expectApproxEqRel(@as(f64, 0.25), p.at(0, 0), 1e-15);
    try testing.expectApproxEqRel(@as(f64, 1e3), p.at(1, 1), 1e-14);
    try testing.expectEqual(@as(f64, 0), p.at(2, 2));
    // zero and empty
    const z: matrix.M3f64 = .{ .data = @splat(@splat(0)) };
    try testing.expectEqual(@as(usize, 0), try rank(al, test_io, &z, null));
    try testing.expect(std.math.isInf(try cond(al, test_io, &z)));
    const pz = try pinv(al, test_io, &z, null);
    defer pz.deinit(al);
    try testing.expect(maxAbs(f64, pz.constView()) == 0);
    const e = try ConstView(f64).init(&.{}, 0, 3, 1);
    try testing.expectEqual(@as(usize, 0), try rank(al, test_io, e, null));
    try testing.expectEqual(@as(f64, 0), try cond(al, test_io, e));
    const pe = try pinv(al, test_io, e, null);
    defer pe.deinit(al);
    try testing.expect(pe.rows == 3 and pe.cols == 0);
}

/// A wider `T` for references: `f64` for `f32`, `f128` for `f64`.
fn Wide(comptime T: type) type {
    const W = if (Real(T) == f32) f64 else f128;
    return if (isComplex(T)) std.math.Complex(W) else W;
}

fn widen(x: anytype) Wide(@TypeOf(x)) {
    return if (comptime isComplex(@TypeOf(x))) .init(x.re, x.im) else x;
}

/// `|re| + |im|` of a wide value, as `f64`.
fn wideAbs(x: anytype) f64 {
    return @floatCast(if (comptime isComplex(@TypeOf(x))) @abs(x.re) + @abs(x.im) else @abs(x));
}

fn wideDiff(x: anytype, ref: Wide(@TypeOf(x))) f64 {
    const w = widen(x);
    return wideAbs(if (comptime isComplex(@TypeOf(x))) w.sub(ref) else w - ref);
}

/// `acc + x y`, unfused.
fn wideMulAdd(acc: anytype, x: @TypeOf(acc), y: @TypeOf(acc)) @TypeOf(acc) {
    return if (comptime isComplex(@TypeOf(acc))) acc.add(x.mul(y)) else acc + x * y;
}

/// `det(a)` by LU with partial pivoting in the wide type.
fn wideDet(comptime T: type, a: ConstView(T)) !Wide(T) {
    const W = Wide(T);
    const n = a.rows;
    const f = try testing.allocator.alloc(W, n * n);
    defer testing.allocator.free(f);
    for (0..n) |j| for (0..n) |i| {
        f[i + j * n] = widen(a.at(i, j));
    };
    var d = ops.one(W);
    for (0..n) |k| {
        var p = k;
        for (k + 1..n) |i| {
            if (wideAbs(f[i + k * n]) > wideAbs(f[p + k * n])) p = i;
        }
        if (p != k) {
            d = ops.neg(W, d);
            for (0..n) |j| std.mem.swap(W, &f[k + j * n], &f[p + j * n]);
        }
        const piv = f[k + k * n];
        d = if (comptime isComplex(W)) d.mul(piv) else d * piv;
        for (k + 1..n) |i| {
            const l = ops.neg(W, if (comptime isComplex(W)) f[i + k * n].div(piv) else f[i + k * n] / piv);
            for (k + 1..n) |j| f[i + j * n] = wideMulAdd(f[i + j * n], l, f[k + j * n]);
        }
    }
    return d;
}

fn expectPaths(comptime T: type, n: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    const W = Wide(T);
    const a = try randomMat(T, rnd, n, n);
    defer a.deinit(al);
    shiftDiag(T, a, 2);
    const b = try randomMat(T, rnd, n, 3);
    defer b.deinit(al);
    const x = b.constView().col(0);
    const y = b.constView().col(1);
    const sx = try strided(constOf(x));
    const sy = try strided(constOf(y));
    // the references, and the sums of magnitudes that bound their errors
    var dot_ref = ops.zero(W);
    var sum_ref = ops.zero(W);
    var sq: Real(W) = 0;
    var dot_mag: f64 = 0;
    var sum_mag: f64 = 0;
    for (x, y) |xi, yi| {
        dot_ref = wideMulAdd(dot_ref, ops.conj(W, widen(xi)), widen(yi));
        sum_ref = wideMulAdd(sum_ref, widen(xi), ops.one(W));
        const w = widen(xi);
        sq += if (comptime isComplex(W)) w.re * w.re + w.im * w.im else w * w;
        dot_mag += wideAbs(widen(xi)) * wideAbs(widen(yi));
        sum_mag += wideAbs(widen(xi));
    }
    const prod_ref = try al.alloc(W, n * 3);
    defer al.free(prod_ref);
    for (0..3) |j| for (0..n) |i| {
        var acc = ops.zero(W);
        for (0..n) |p| acc = wideMulAdd(acc, widen(a.at(i, p)), widen(b.at(p, j)));
        prod_ref[i + j * n] = acc;
    };
    const det_ref = try wideDet(T, a.constView());
    const tol = tolOf(T, n);
    const nf: f64 = @floatFromInt(n);
    for ([_]Path{ .fixed, .blas }) |path| {
        try testing.expect(wideDiff(try dotWith(path, test_io, sx, sy), dot_ref) <= tol * dot_mag);
        try testing.expect(wideDiff(try sumWith(path, test_io, sx), sum_ref) <= tol * sum_mag);
        try testing.expect(wideDiff(try normWith(path, test_io, sx), @sqrt(sq)) <= tol * @as(f64, @floatCast(@sqrt(sq))));
        if (comptime !isComplex(T)) {
            try testing.expectEqual(std.mem.max(T, x), try extremeWith(.max, path, test_io, sx));
            try testing.expectEqual(std.mem.min(T, x), try extremeWith(.min, path, test_io, sx));
        }
        // `gemm` and, one column, `gemv`
        const o = try Mat(T).init(al, n, 3);
        defer o.deinit(al);
        const mul_bound = tol * nf * maxAbs(T, a.constView()) * maxAbs(T, b.constView());
        for ([_]usize{ 3, 1 }) |k| {
            try mulWith(path, al, test_io, try o.view().sub(0, 0, n, k), a.constView(), try b.constView().sub(0, 0, n, k));
            for (0..k) |j| for (0..n) |i| try testing.expect(wideDiff(o.at(i, j), prod_ref[i + j * n]) <= mul_bound);
        }
        try testing.expect(wideDiff(try detWith(path, al, test_io, a.constView()), det_ref) <= tol * nf * wideAbs(det_ref));
        const ai = try Mat(T).init(al, n, n);
        defer ai.deinit(al);
        try invWith(path, al, test_io, ai.view(), a.constView());
        const p = try mul(al, test_io, a, ai);
        defer p.deinit(al);
        const id = try Mat(T).zeros(al, n, n);
        defer id.deinit(al);
        shiftDiag(T, id, 1);
        try testing.expect(maxDiff(T, p.constView(), id.constView()) <= tol * nf * maxAbs(T, a.constView()) * maxAbs(T, ai.constView()));
        try solveWith(path, al, test_io, o.view(), a.constView(), b.constView());
        try expectSolves(T, a.constView(), o.constView(), b.constView());
    }
}

test "each path at every size 1 to 16, against a wide reference" {
    var prng: std.Random.DefaultPrng = .init(0xd15);
    inline for (.{ f32, f64, C64 }) |T| for (1..17) |n| try expectPaths(T, n, prng.random());
}

test {
    zigstd.testRefAllDecls(@This());
}
