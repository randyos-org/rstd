//! `rstd.linalg`'s run-time-shape functions (design.md, "Operations over
//! views"), re-exported flat at `rstd.linalg`. Each operand is a `View` or
//! `ConstView`, a `Mat`, a pointer to a `StaticMat` or a slice (a column),
//! in any mix; a vector is a view of one column or one row. Shapes at or
//! below the inline limit take the fixed-size forms' operation order, so a
//! `Mat` and a `StaticMat` holding the same numbers give the same bits; past
//! it, the BLAS and LAPACK. Outputs must not overlap inputs.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Io = std.Io;
const blas = @import("blas/__root__.zig");
const lapack = @import("lapack/__root__.zig");
const vector = @import("vector.zig");
const matrix = @import("matrix.zig");

const View = matrix.View;
const ConstView = matrix.ConstView;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const eZero = matrix.eZero;
const eOne = matrix.eOne;
const eMul = matrix.eMul;
const eMulAdd = matrix.eMulAdd;
const eNeg = matrix.eNeg;
const eConj = matrix.eConj;
const eIsZero = matrix.eIsZero;

/// The inline limit `L` per dimension (design.md, Dispatch); M11 makes it
/// `rstd_linalg_options.inline_limit`.
const inline_limit = 8;

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
            return matrix.det(load(a));
        }

        fn inv(o: anytype, a: ConstView(@TypeOf(o).Element)) error{Singular}!void {
            const r = try matrix.inv(load(a));
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

/// `v` as a vector: one column, or one row (stride `ld`).
fn strided(v: anytype) error{DimensionMismatch}!Strided(@TypeOf(v).Element) {
    const p: [*]@TypeOf(v).Element = @constCast(v.data.ptr);
    if (v.cols == 1 or v.rows == 0 or v.cols == 0) return .{ .ptr = p, .n = if (v.cols == 0) 0 else v.rows, .inc = 1 };
    if (v.rows == 1) return .{ .ptr = p, .n = v.cols, .inc = v.ld };
    return error.DimensionMismatch;
}

fn len(n: usize) isize {
    return @intCast(n);
}

/// A BLAS call's error, of which only `Canceled` can arise here: the
/// arguments come from views, so the BLAS's argument errors are bugs.
fn canceled(e: anyerror) Io.Cancelable {
    return switch (e) {
        error.Canceled => error.Canceled,
        else => unreachable,
    };
}

/// As `canceled`, for the calls that may also allocate.
fn allocOrCanceled(e: anyerror) (Allocator.Error || Io.Cancelable) {
    return switch (e) {
        error.Canceled => error.Canceled,
        error.OutOfMemory => error.OutOfMemory,
        else => unreachable,
    };
}

pub const VectorError = error{DimensionMismatch} || Io.Cancelable;

/// `xᴴ y` (`Σ conj(x[i]) y[i]`, MATLAB's `dot`); `linalg.vector.dot`'s order
/// up to the inline limit, else the BLAS's `dot` or `dotc`.
pub fn dot(io: Io, x: anytype, y: anytype) VectorError!ElementOf(@TypeOf(x)) {
    const T = ElementOf(@TypeOf(x));
    const sx = try strided(constOf(x));
    const sy = try strided(constOf(y));
    if (sx.n != sy.n) return error.DimensionMismatch;
    switch (sx.n) {
        0 => return eZero(T),
        inline 1...inline_limit => |c| return Fixed(c).dot(sx, sy),
        else => {},
    }
    const f = if (comptime isComplex(T)) blas.dotc else blas.dot;
    return f(io, T, len(sx.n), sx.ptr, len(sx.inc), sy.ptr, len(sy.inc)) catch |e| return canceled(e);
}

/// The 2-norm (MATLAB's `norm` of a vector; `linalg.vector`'s `norm`
/// normalizes): `linalg.vector.length` up to the inline limit, else `nrm2`.
pub fn norm(io: Io, x: anytype) VectorError!Real(ElementOf(@TypeOf(x))) {
    const T = ElementOf(@TypeOf(x));
    const s = try strided(constOf(x));
    switch (s.n) {
        0 => return 0,
        inline 1...inline_limit => |c| return Fixed(c).length(s),
        else => {},
    }
    return blas.nrm2(io, T, len(s.n), s.ptr, len(s.inc)) catch |e| return canceled(e);
}

/// `Σ x[i]`: in order up to the inline limit, else the BLAS's `sum`, the
/// real and imaginary parts apart for complex.
pub fn sum(io: Io, x: anytype) VectorError!ElementOf(@TypeOf(x)) {
    const T = ElementOf(@TypeOf(x));
    const s = try strided(constOf(x));
    switch (s.n) {
        0 => return eZero(T),
        inline 1...inline_limit => |c| return Fixed(c).sum(s),
        else => {},
    }
    if (comptime !isComplex(T)) return blas.sum(io, T, len(s.n), s.ptr, len(s.inc)) catch |e| return canceled(e);
    const F = Real(T);
    const parts: [*]F = @ptrCast(s.ptr);
    const inc = len(2 * s.inc);
    return .{
        .re = blas.sum(io, F, len(s.n), parts, inc) catch |e| return canceled(e),
        .im = blas.sum(io, F, len(s.n), parts + 1, inc) catch |e| return canceled(e),
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
    switch (s.n) {
        0 => return 0,
        inline 1...inline_limit => |c| return Fixed(c).extreme(which, s),
        else => {},
    }
    return blas.extreme(io, T, which, len(s.n), s.ptr, len(s.inc)) catch |e| return canceled(e);
}

/// `y += a x`, by the BLAS's `axpy`.
pub fn axpy(io: Io, a: anytype, x: anytype, y: anytype) VectorError!void {
    const T = ElementOf(@TypeOf(y));
    const sx = try strided(constOf(x));
    const sy = try strided(viewOf(y).constView());
    if (sx.n != sy.n) return error.DimensionMismatch;
    const alpha: T = a;
    blas.axpy(io, T, len(sx.n), alpha, sx.ptr, len(sx.inc), sy.ptr, len(sy.inc)) catch |e| return canceled(e);
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
    const T = ElementOf(@TypeOf(a));
    const o = viewOf(out);
    const ca = constOf(a);
    const cb = constOf(b);
    if (ca.cols != cb.rows or o.rows != ca.rows or o.cols != cb.cols) return error.DimensionMismatch;
    const m = ca.rows;
    const n = cb.cols;
    const k = ca.cols;
    if (m <= inline_limit and n <= inline_limit and k <= inline_limit) {
        for (0..n) |j| for (0..m) |i| {
            var acc = eZero(T);
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
    const one: T = if (comptime isComplex(T)) T.init(1, 0) else 1;
    if (n == 1) {
        blas.gemv(allocator, io, T, .n, len(m), len(k), one, pa, len(ca.ld), pb, 1, eZero(T), o.data.ptr, 1) catch |e| return allocOrCanceled(e);
        return;
    }
    blas.gemm(allocator, io, T, .n, .n, len(m), len(n), len(k), one, pa, len(ca.ld), pb, len(cb.ld), eZero(T), o.data.ptr, len(o.ld)) catch |e| return allocOrCanceled(e);
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
    transposeWith(trans, io, r.view(), ca) catch |e| return switch (e) {
        error.DimensionMismatch => unreachable,
        error.Canceled => error.Canceled,
    };
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
    const one: T = if (comptime isComplex(T)) T.init(1, 0) else 1;
    const t: blas.Trans = if (comptime isComplex(T)) trans else .t;
    blas.omatcopy(io, T, t, len(a.rows), len(a.cols), one, a.data.ptr, len(a.ld), o.data.ptr, len(o.ld)) catch |e| return canceled(e);
}

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

/// As `allocOrCanceled`, passing `Singular` through.
fn singularOr(e: anyerror) (error{Singular} || Allocator.Error || Io.Cancelable) {
    return switch (e) {
        error.Singular => error.Singular,
        else => allocOrCanceled(e),
    };
}

pub const DetError = error{NotSquare} || Allocator.Error || Io.Cancelable;

/// The determinant (MATLAB's `det`); 1 for 0 by 0, 0 when singular.
/// `linalg.matrix.det` up to the inline limit; else the product of `getrf`'s
/// pivots, renormalized by powers of 2 as it goes so that it overflows or
/// underflows only if the result does (design.md, M1's E1). The scaling is
/// exact, so where the plain product stays in range it has its bits.
pub fn det(allocator: Allocator, io: Io, a: anytype) DetError!ElementOf(@TypeOf(a)) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    switch (n) {
        0 => return eOne(T),
        inline 1...inline_limit => |c| return Fixed(c).det(ca),
        else => {},
    }
    const f = try copyOf(T, allocator, ca);
    defer f.deinit(allocator);
    const ipiv = try allocator.alloc(isize, n);
    defer allocator.free(ipiv);
    lapack.getrf(allocator, io, T, len(n), len(n), f.data.ptr, len(n), ipiv.ptr, null) catch |e| switch (e) {
        error.Singular => return eZero(T),
        else => return allocOrCanceled(e),
    };
    var d = eOne(T);
    var exp: i64 = 0;
    for (ipiv, 0..) |p, i| {
        d = eMul(d, f.at(i, i));
        if (p != i + 1) d = eNeg(d);
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
    try invWith(allocator, io, r.view(), ca);
    return r;
}

/// `out = a⁻¹`, `out` untouched on an error.
pub fn invInto(allocator: Allocator, io: Io, out: anytype, a: anytype) (InvError || error{DimensionMismatch})!void {
    const o = viewOf(out);
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    if (o.rows != ca.rows or o.cols != ca.cols) return error.DimensionMismatch;
    return invWith(allocator, io, o, ca);
}

fn invWith(allocator: Allocator, io: Io, o: anytype, a: ConstView(@TypeOf(o).Element)) (error{Singular} || Allocator.Error || Io.Cancelable)!void {
    const T = @TypeOf(o).Element;
    const n = a.rows;
    switch (n) {
        0 => return,
        inline 1...inline_limit => |c| return Fixed(c).inv(o, a),
        else => {},
    }
    const f = try copyOf(T, allocator, a);
    defer f.deinit(allocator);
    const ipiv = try allocator.alloc(isize, n);
    defer allocator.free(ipiv);
    lapack.getrf(allocator, io, T, len(n), len(n), f.data.ptr, len(n), ipiv.ptr, null) catch |e| return singularOr(e);
    lapack.getri(allocator, io, T, len(n), f.data.ptr, len(n), ipiv.ptr, null) catch |e| return singularOr(e);
    copyInto(T, o.data, o.ld, f.constView());
}

pub const SolveError = error{ DimensionMismatch, Singular } || Allocator.Error || Io.Cancelable;

/// `a \ b` (MATLAB's `mldivide`), owned by the caller, `a` `m` by `n` and `b`
/// `m` by `k`. Square `a` up to the inline limit with `k` within it too:
/// `linalg.matrix.solve`'s LU. Past it, by `a`'s values alone (MATLAB's
/// order): triangular by `trtrs`; Hermitian with a positive real diagonal
/// by `potrf` and `potrs`, LU if that fails; else `gesv`. Square `a` gives
/// `error.Singular` with no result. Rectangular `a`: the least-squares
/// solution of least norm by `gelsy`, its rank cut at `max(m, n) eps`, as
/// MATLAB's QR tolerance.
pub fn solve(allocator: Allocator, io: Io, a: anytype, b: anytype) SolveError!Mat(ElementOf(@TypeOf(a))) {
    const ca = constOf(a);
    const cb = constOf(b);
    if (cb.rows != ca.rows) return error.DimensionMismatch;
    const r = try Mat(ElementOf(@TypeOf(a))).init(allocator, ca.cols, cb.cols);
    errdefer r.deinit(allocator);
    try solveWith(allocator, io, r.view(), ca, cb);
    return r;
}

/// `out = a \ b`, `out` untouched on an error.
pub fn solveInto(allocator: Allocator, io: Io, out: anytype, a: anytype, b: anytype) SolveError!void {
    const o = viewOf(out);
    const ca = constOf(a);
    const cb = constOf(b);
    if (cb.rows != ca.rows or o.rows != ca.cols or o.cols != cb.cols) return error.DimensionMismatch;
    return solveWith(allocator, io, o, ca, cb);
}

const Structure = enum { upper, lower, hermitian, general };

/// `solve`'s route past the inline limit, from `a`'s values alone.
fn structureOf(a: anytype) Structure {
    const T = @TypeOf(a).Element;
    const n = a.rows;
    var upper = true;
    var lower = true;
    for (0..n) |j| for (0..n) |i| {
        if (eIsZero(a.at(i, j))) continue;
        if (i > j) upper = false;
        if (i < j) lower = false;
    };
    if (upper) return .upper;
    if (lower) return .lower;
    for (0..n) |j| {
        const d = a.at(j, j);
        const positive = if (comptime isComplex(T)) d.im == 0 and d.re > 0 else d > 0;
        if (!positive) return .general;
        for (j + 1..n) |i| {
            const x = a.at(i, j);
            const y = eConj(a.at(j, i));
            const same = if (comptime isComplex(T)) x.re == y.re and x.im == y.im else x == y;
            if (!same) return .general;
        }
    }
    return .hermitian;
}

fn solveWith(allocator: Allocator, io: Io, o: anytype, a: ConstView(@TypeOf(o).Element), b: ConstView(@TypeOf(o).Element)) SolveError!void {
    const T = @TypeOf(o).Element;
    const n = a.cols;
    if (a.rows != n) return leastSquares(allocator, io, o, a, b);
    if (n <= inline_limit and b.cols <= inline_limit) switch (n) {
        0 => return,
        inline 1...inline_limit => |c| return Fixed(c).solve(o, a, b),
        else => unreachable,
    };
    const x = try copyOf(T, allocator, b);
    defer x.deinit(allocator);
    const nn = len(n);
    const k = len(b.cols);
    const ld = len(@max(n, 1));
    switch (structureOf(a)) {
        .upper, .lower => |s| {
            const uplo: lapack.Uplo = if (s == .upper) .u else .l;
            lapack.trtrs(allocator, io, T, uplo, .n, .n, nn, k, a.data.ptr, len(a.ld), x.data.ptr, ld, null) catch |e| return singularOr(e);
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
                    else => return allocOrCanceled(e),
                };
                lapack.potrs(allocator, io, T, .u, nn, k, f.data.ptr, ld, x.data.ptr, ld) catch |e| return allocOrCanceled(e);
                return copyInto(T, o.data, o.ld, x.constView());
            }
            const ipiv = try allocator.alloc(isize, n);
            defer allocator.free(ipiv);
            lapack.gesv(allocator, io, T, nn, k, f.data.ptr, ld, ipiv.ptr, x.data.ptr, ld, null) catch |e| return singularOr(e);
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
        for (0..o.cols) |j| @memset(o.col(j), eZero(T));
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
    var rank: isize = undefined;
    lapack.gelsy(allocator, io, T, len(m), len(n), len(b.cols), f.data.ptr, len(m), x.ptr, len(ldx), jpvt.ptr, rcond, &rank) catch |e| return allocOrCanceled(e);
    for (0..o.cols) |j| @memcpy(o.col(j), x[j * ldx ..][0..n]);
}

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
                r.data[j + j * r.rows] = eOne(T);
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
    lapack.getrf(allocator, io, T, len(ca.rows), len(ca.cols), f.data.ptr, len(@max(ca.rows, 1)), piv.ptr, &zero) catch |e| switch (e) {
        error.Singular => zero_pivot = zero - 1,
        else => return allocOrCanceled(e),
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
    lapack.potrf(allocator, io, T, .u, len(n), r.data.ptr, len(@max(n, 1)), null) catch |e| return switch (e) {
        error.NotPositiveDefinite => error.NotPositiveDefinite,
        else => allocOrCanceled(e),
    };
    for (0..n) |j| @memset(r.data[j * n + j + 1 .. (j + 1) * n], eZero(T));
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
    const ld = len(@max(m, 1));
    lapack.geqrf(allocator, io, T, len(m), len(n), f.data.ptr, ld, tau.ptr) catch |e| return allocOrCanceled(e);
    const r = try Mat(T).zeros(allocator, qc, n);
    errdefer r.deinit(allocator);
    for (0..n) |j| for (0..@min(j + 1, k)) |i| {
        r.data[i + j * qc] = f.at(i, j);
    };
    const q = try Mat(T).zeros(allocator, m, qc);
    errdefer q.deinit(allocator);
    for (0..k) |j| @memcpy(q.data[j * m ..][0..m], f.data[j * m ..][0..m]);
    if (m > 0) lapack.orgqr(allocator, io, T, len(m), len(qc), len(k), q.data.ptr, ld, tau.ptr) catch |e| return allocOrCanceled(e);
    return .{ .q = q, .r = r };
}

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
            const f = try copyOf(T, allocator, ca);
            defer f.deinit(allocator);
            const s = try allocator.alloc(Real(T), @min(m, n));
            defer allocator.free(s);
            var none: [1]T = undefined;
            lapack.gesvd(allocator, io, T, .n, .n, len(m), len(n), f.data.ptr, len(m), s.ptr, &none, 1, &none, 1, null) catch |e| return switch (e) {
                error.NoConvergence => error.NoConvergence,
                else => allocOrCanceled(e),
            };
            return s[0];
        },
    };
    return lapack.lange(allocator, io, T, kind, len(m), len(n), ca.data.ptr, len(ca.ld)) catch |e| return switch (e) {
        error.OutOfMemory => error.OutOfMemory,
        else => unreachable,
    };
}

const zigstd = @import("../zigstd.zig");
const testing = std.testing;
const C64 = std.math.Complex(f64);
const test_io = Io.failing;

test "views: init's checks, sub, col; StaticMat and Mat views share storage" {
    var buf = [_]f64{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    try testing.expectError(error.DimensionMismatch, View(f64).init(&buf, 3, 2, 2));
    try testing.expectError(error.SliceTooShort, View(f64).init(&buf, 3, 4, 3));
    const v = try View(f64).init(&buf, 3, 3, 3);
    try testing.expectEqual(@as(f64, 8), v.at(1, 2));
    const s = try v.sub(1, 1, 2, 2);
    try testing.expectEqual(@as(f64, 5), s.at(0, 0));
    try testing.expectEqual(@as(f64, 9), s.at(1, 1));
    try testing.expectEqualSlices(f64, &.{ 5, 6 }, s.col(0));
    try testing.expectError(error.DimensionMismatch, v.sub(2, 0, 2, 1));

    var sm = matrix.fromCols(matrix.M2f64, .{ .{ 1, 2 }, .{ 3, 4 } });
    sm.view().data[3] = 40;
    try testing.expectEqual(@as(f64, 40), sm.data[1][1]);
    try testing.expectEqual(@as(f64, 3), sm.constView().at(0, 1));
    const m = try Mat(f64).zeros(testing.allocator, 2, 3);
    defer m.deinit(testing.allocator);
    m.view().data[5] = 7;
    try testing.expectEqual(@as(f64, 7), m.at(1, 2));
}

fn randomElem(comptime T: type, rnd: std.Random) T {
    if (comptime isComplex(T)) return .init(rnd.float(f64) * 2 - 1, rnd.float(f64) * 2 - 1);
    return std.math.ldexp(rnd.float(T) * 2 - 1, rnd.intRangeAtMost(i32, -4, 4));
}

fn fill(comptime T: type, rnd: std.Random, s: []T) void {
    for (s) |*x| x.* = randomElem(T, rnd);
}

fn sameBits(comptime T: type, a: []const T, b: []const T) bool {
    return std.mem.eql(u8, std.mem.sliceAsBytes(a), std.mem.sliceAsBytes(b));
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
        for (j + 1..n) |i| a.data[i + j * n] = eConj(a.data[j + i * n]);
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
        if (i > j) upper.data[i + j * n] = eZero(T);
        if (i < j) lower.data[i + j * n] = eZero(T);
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
        for (0..n) |j| for (j + 1..n) |i| try testing.expect(eIsZero(r.at(i, j)));
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
    for (0..a.cols) |j| for (@min(j + 1, qc)..qc) |i| try testing.expect(eIsZero(f.r.at(i, j)));
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

test {
    zigstd.testRefAllDecls(@This());
}
