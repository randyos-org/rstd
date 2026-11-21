//! Dot products, norms, sums, extremes, `axpy` and `scale` over vectors.

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

const util = @import("util.zig");
const ElementOf = util.ElementOf;
const Path = util.Path;
const constOf = util.constOf;
const fixed = util.fixed;
const narrow = util.narrow;
const pathOf = util.pathOf;
const strided = util.strided;
const toIsize = util.toIsize;
const viewOf = util.viewOf;

pub const VectorError = error{DimensionMismatch} || Io.Cancelable;

/// `xᴴ y` (`Σ conj(x[i]) y[i]`, MATLAB's `dot`); `linalg.vector.dot`'s order
/// up to the inline limit, else the BLAS's `dot` or `dotc`.
pub fn dot(io: Io, x: anytype, y: anytype) VectorError!ElementOf(@TypeOf(x)) {
    const sx = try strided(constOf(x));
    const sy = try strided(constOf(y));
    if (sx.n != sy.n) return error.DimensionMismatch;
    return dotWith(pathOf(&.{sx.n}), io, sx, sy);
}

pub fn dotWith(path: Path, io: Io, sx: anytype, sy: @TypeOf(sx)) Io.Cancelable!@TypeOf(sx.ptr[0]) {
    const T = @TypeOf(sx.ptr[0]);
    if (sx.n == 0) return ops.zero(T);
    if (path == .fixed) return fixed(T, "dot", sx.n, .{ sx, sy });
    const f = if (comptime isComplex(T)) l1.dotc else l1.dot;
    return f(io, T, toIsize(sx.n), sx.ptr, toIsize(sx.inc), sy.ptr, toIsize(sy.inc)) catch |e| return narrow(Io.Cancelable, e);
}

/// The 2-norm (MATLAB's `norm` of a vector; `linalg.vector`'s `norm`
/// normalizes): `linalg.vector.length` up to the inline limit, else `nrm2`.
pub fn norm(io: Io, x: anytype) VectorError!Real(ElementOf(@TypeOf(x))) {
    const s = try strided(constOf(x));
    return normWith(pathOf(&.{s.n}), io, s);
}

pub fn normWith(path: Path, io: Io, s: anytype) Io.Cancelable!Real(@TypeOf(s.ptr[0])) {
    const T = @TypeOf(s.ptr[0]);
    if (s.n == 0) return 0;
    if (path == .fixed) return fixed(Real(T), "length", s.n, .{s});
    return l1.nrm2(io, T, toIsize(s.n), s.ptr, toIsize(s.inc)) catch |e| return narrow(Io.Cancelable, e);
}

/// `Σ x[i]`: in order up to the inline limit, else the BLAS's `sum`, the
/// real and imaginary parts apart for complex.
pub fn sum(io: Io, x: anytype) VectorError!ElementOf(@TypeOf(x)) {
    const s = try strided(constOf(x));
    return sumWith(pathOf(&.{s.n}), io, s);
}

pub fn sumWith(path: Path, io: Io, s: anytype) Io.Cancelable!@TypeOf(s.ptr[0]) {
    const T = @TypeOf(s.ptr[0]);
    if (s.n == 0) return ops.zero(T);
    if (path == .fixed) return fixed(T, "sum", s.n, .{s});
    if (comptime !isComplex(T)) return l1.sum(io, T, toIsize(s.n), s.ptr, toIsize(s.inc)) catch |e| return narrow(Io.Cancelable, e);
    const F = Real(T);
    const parts: [*]F = @ptrCast(s.ptr);
    const inc = toIsize(2 * s.inc);
    return .{
        .re = l1.sum(io, F, toIsize(s.n), parts, inc) catch |e| return narrow(Io.Cancelable, e),
        .im = l1.sum(io, F, toIsize(s.n), parts + 1, inc) catch |e| return narrow(Io.Cancelable, e),
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

pub fn extreme(comptime which: blas.Extreme, io: Io, x: anytype) VectorError!ElementOf(@TypeOf(x)) {
    const T = ElementOf(@TypeOf(x));
    if (comptime isComplex(T)) @compileError("max and min order real elements only");
    const s = try strided(constOf(x));
    return extremeWith(which, pathOf(&.{s.n}), io, s);
}

pub fn extremeWith(comptime which: blas.Extreme, path: Path, io: Io, s: anytype) Io.Cancelable!@TypeOf(s.ptr[0]) {
    const T = @TypeOf(s.ptr[0]);
    if (s.n == 0) return 0;
    if (path == .fixed) return fixed(T, "extreme", s.n, .{ which, s });
    return l1.extreme(io, T, which, toIsize(s.n), s.ptr, toIsize(s.inc)) catch |e| return narrow(Io.Cancelable, e);
}

/// `y += a x`, by the BLAS's `axpy`.
pub fn axpy(io: Io, a: anytype, x: anytype, y: anytype) VectorError!void {
    const T = ElementOf(@TypeOf(y));
    const sx = try strided(constOf(x));
    const sy = try strided(viewOf(y).constView());
    if (sx.n != sy.n) return error.DimensionMismatch;
    const alpha: T = a;
    l1.axpy(io, T, toIsize(sx.n), alpha, sx.ptr, toIsize(sx.inc), sy.ptr, toIsize(sy.inc)) catch |e| return narrow(Io.Cancelable, e);
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
