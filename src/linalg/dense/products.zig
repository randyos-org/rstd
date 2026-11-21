//! Products and transposes.

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
const Error = util.Error;
const Path = util.Path;
const constOf = util.constOf;
const narrow = util.narrow;
const pathOf = util.pathOf;
const toIsize = util.toIsize;
const viewOf = util.viewOf;

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

pub fn mulWith(path: Path, allocator: Allocator, io: Io, o: anytype, ca: ConstView(@TypeOf(o).Element), cb: ConstView(@TypeOf(o).Element)) (Allocator.Error || Io.Cancelable)!void {
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
        l2.gemv(allocator, io, T, .n, toIsize(m), toIsize(k), one, pa, toIsize(ca.ld), pb, 1, ops.zero(T), o.data.ptr, 1) catch |e| return narrow(Error, e);
        return;
    }
    l3.gemm(allocator, io, T, .n, .n, toIsize(m), toIsize(n), toIsize(k), one, pa, toIsize(ca.ld), pb, toIsize(cb.ld), ops.zero(T), o.data.ptr, toIsize(o.ld)) catch |e| return narrow(Error, e);
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

pub fn transposed(comptime trans: blas.Trans, allocator: Allocator, io: Io, a: anytype) (Allocator.Error || Io.Cancelable)!Mat(ElementOf(@TypeOf(a))) {
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

pub fn transposeWith(comptime trans: blas.Trans, io: Io, o: anytype, a: anytype) TransposeError!void {
    const T = @TypeOf(a).Element;
    if (o.rows != a.cols or o.cols != a.rows) return error.DimensionMismatch;
    if (a.rows == 0 or a.cols == 0) return;
    const one = ops.one(T);
    const t: blas.Trans = if (comptime isComplex(T)) trans else .t;
    ext.omatcopy(io, T, t, toIsize(a.rows), toIsize(a.cols), one, a.data.ptr, toIsize(a.ld), o.data.ptr, toIsize(o.ld)) catch |e| return narrow(Io.Cancelable, e);
}
