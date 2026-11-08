//! The BLAS calls the LAPACK routines make, on `Mat` views and plain
//! pointers. Level 3 goes through the drivers on the interface's worker rule,
//! so Accelerate never runs it; level 1 and 2 through the public routines.
//! Every threaded call keeps the BLAS's rule: the same bits for any worker
//! count.

const std = @import("std");
const blas = @import("../blas/__root__.zig");
const types = @import("../blas/types.zig");
const driver = @import("../blas/blas3/driver.zig");
const threaded = @import("../blas/blas3/threaded.zig");
const threading = @import("../blas/threading.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const Mat = ops.Mat;
const Error = ops.Error;
const scalar = ops.scalar;
const zero = ops.zero;
const add = ops.add;
const sub = ops.sub;
const mul = ops.mul;
const conj = ops.conj;
const mulReal = ops.mulReal;
const abs1 = ops.abs1;
const isZero = ops.isZero;
const narrow = ops.narrow;
const table = interface.table;

/// How a level-3 operand enters: as is, transposed, conjugated, or
/// conjugate transposed (real: `r` is `n`, `c` is `t`).
pub const Op = enum {
    n,
    t,
    r,
    c,

    fn transposed(op: Op) bool {
        return op == .t or op == .c;
    }

    /// Only for a complex `T`.
    fn conjugated(op: Op, comptime T: type) bool {
        return isComplex(T) and (op == .r or op == .c);
    }

    /// `n` or `r` by `conj`, `t` or `c` by `conj` and `transpose`.
    pub fn of(transpose: bool, conj_: bool) Op {
        return if (transpose) (if (conj_) .c else .t) else (if (conj_) .r else .n);
    }
};

/// `C := alpha op(A) op(B) + beta C`, `C` being `m x n`.
pub fn gemm(al: Allocator, io: Io, comptime T: type, ta: Op, tb: Op, alpha: T, a: Mat(T), b: Mat(T), beta: T, c: Mat(T)) Error!void {
    const k = if (ta.transposed()) a.rows else a.cols;
    if (c.rows == 0 or c.cols == 0) return;
    const w = threaded.workers(.gemm, T, threading.product(&.{ c.rows, c.cols, k }));
    try driver.gemm(al, io, T, table(T), ta.transposed(), tb.transposed(), .of(ta.conjugated(T), tb.conjugated(T)), c.rows, c.cols, k, alpha, a.ptr, a.ld, b.ptr, b.ld, beta, c.ptr, c.ld, w);
}

/// `B := alpha op(A) B` (`solve`: `alpha op(A)⁻¹ B`) from the left, `A`
/// triangular.
pub fn tri(al: Allocator, io: Io, comptime T: type, comptime solve: bool, upper: bool, op: Op, unit: bool, alpha: T, a: Mat(T), b: Mat(T)) Error!void {
    if (b.rows == 0 or b.cols == 0) return;
    const w = threaded.workers(.trsm, T, threading.product(&.{ b.rows, b.cols }));
    const f = if (solve) driver.trsm else driver.trmm;
    try f(al, io, T, table(T), true, upper, op.transposed(), op.conjugated(T), unit, b.rows, b.cols, alpha, a.ptr, a.ld, b.ptr, b.ld, w);
}

/// `B := alpha B op(A)` (`solve`: `alpha B op(A)⁻¹`) from the right, `A`
/// triangular.
pub fn triRight(al: Allocator, io: Io, comptime T: type, comptime solve: bool, upper: bool, op: Op, unit: bool, alpha: T, a: Mat(T), b: Mat(T)) Error!void {
    if (b.rows == 0 or b.cols == 0) return;
    const w = threaded.workers(.trsm, T, threading.product(&.{ b.rows, b.cols }));
    const f = if (solve) driver.trsm else driver.trmm;
    try f(al, io, T, table(T), false, upper, op.transposed(), op.conjugated(T), unit, b.rows, b.cols, alpha, a.ptr, a.ld, b.ptr, b.ld, w);
}

/// The upper triangle of `C := alpha Aᴴ A + beta C` (real: `Aᵀ A`), `A`
/// being `k x n`; a complex diagonal comes out real.
pub fn gramUpper(al: Allocator, io: Io, comptime T: type, alpha: Float(T), a: Mat(T), beta: Float(T), c: Mat(T)) Error!void {
    const n = a.cols;
    if (n == 0) return;
    const w = threaded.workers(.syrk, T, threading.product(&.{ n + 1, n, a.rows }));
    try driver.syrk(al, io, T, table(T), comptime isComplex(T), true, true, n, a.rows, scalar(T, alpha), a.ptr, a.ld, null, 0, scalar(T, beta), c.ptr, c.ld, w);
}

/// `y := alpha op(A) x + beta y` (`blas.gemv`; `trans` as its letters), `A`
/// being `m x n`. Unlike `gemv`, `y` is scaled by `beta` even for an empty
/// product.
pub fn gemv(al: Allocator, io: Io, comptime T: type, trans: Trans, alpha: T, a: Mat(T), x: [*]T, incx: usize, beta: T, y: [*]T, incy: usize) Error!void {
    const leny = if (trans.transposed()) a.cols else a.rows;
    if (a.rows == 0 or a.cols == 0) {
        for (0..leny) |i| y[i * incy] = if (isZero(T, beta)) zero(T) else mul(T, beta, y[i * incy]);
        return;
    }
    const t: Trans = if (comptime isComplex(T)) trans else switch (trans) {
        .n, .r, .o, .s => .n,
        .t, .c, .u, .d => .t,
    };
    blas.gemv(al, io, T, t, @intCast(a.rows), @intCast(a.cols), alpha, a.ptr, @intCast(a.ld), x, @intCast(incx), beta, y, @intCast(incy)) catch |e| return narrow(Error, e);
}

/// `A += alpha x yᵀ` (`blas.ger`), unconjugated.
pub fn ger(al: Allocator, io: Io, comptime T: type, alpha: T, x: [*]T, incx: usize, y: [*]T, incy: usize, a: Mat(T)) Error!void {
    if (a.rows == 0 or a.cols == 0) return;
    blas.ger(al, io, T, @intCast(a.rows), @intCast(a.cols), alpha, x, @intCast(incx), y, @intCast(incy), a.ptr, @intCast(a.ld)) catch |e| return narrow(Error, e);
}

/// `‖x‖₂` (`blas.nrm2`).
pub fn nrm2(io: Io, comptime T: type, n: usize, x: [*]T, inc: usize) Float(T) {
    if (n == 0) return 0;
    return blas.nrm2(io, T, @intCast(n), x, @intCast(inc)) catch unreachable;
}

/// `?SCAL` (`?DSCAL` for complex) by a real `s`.
pub fn scal(io: Io, comptime T: type, n: usize, s: Float(T), x: [*]T, inc: usize) Error!void {
    if (n == 0) return;
    const r = if (comptime isComplex(T)) blas.scalReal(io, T, @intCast(n), s, x, @intCast(inc)) else blas.scal(io, T, @intCast(n), s, x, @intCast(inc));
    r catch |e| return narrow(Error, e);
}

/// `?SCAL` by a complex `s`.
pub fn scalComplex(io: Io, comptime T: type, n: usize, s: T, x: [*]T, inc: usize) Error!void {
    if (n == 0) return;
    blas.scal(io, T, @intCast(n), s, x, @intCast(inc)) catch |e| return narrow(Error, e);
}

/// `?AXPY`.
pub fn axpy(io: Io, comptime T: type, n: usize, alpha: T, x: [*]const T, incx: usize, y: [*]T, incy: usize) Error!void {
    if (n == 0) return;
    blas.axpy(io, T, @intCast(n), alpha, @constCast(x), @intCast(incx), y, @intCast(incy)) catch |e| return narrow(Error, e);
}

/// `?DOT`, `?DOTU`, or with `conj_` `?DOTC` (`conj(x)ᵀ y`).
pub fn dot(io: Io, comptime T: type, conj_: bool, n: usize, x: [*]const T, incx: usize, y: [*]const T, incy: usize) Error!T {
    if (n == 0) return zero(T);
    const r = if (comptime isComplex(T)) (if (conj_) blas.dotc(io, T, @intCast(n), @constCast(x), @intCast(incx), @constCast(y), @intCast(incy)) else blas.dot(io, T, @intCast(n), @constCast(x), @intCast(incx), @constCast(y), @intCast(incy))) else blas.dot(io, T, @intCast(n), @constCast(x), @intCast(incx), @constCast(y), @intCast(incy));
    return r catch |e| return narrow(Error, e);
}

/// The first index of the largest `abs1` (`I?AMAX`); NaN is passed over.
pub fn iamax(comptime T: type, n: usize, x: [*]const T, inc: usize) usize {
    var best: usize = 0;
    var m: Float(T) = -1;
    for (0..n) |i| {
        const v = abs1(T, x[i * inc]);
        if (v > m) {
            m = v;
            best = i;
        }
    }
    return best;
}

/// `?SWAP` through rstd's BLAS.
pub fn swap(io: Io, comptime T: type, n: usize, x: [*]T, incx: usize, y: [*]T, incy: usize) Error!void {
    if (n == 0) return;
    blas.swap(io, T, @intCast(n), x, @intCast(incx), y, @intCast(incy)) catch |e| return narrow(Error, e);
}

/// `DROT` through rstd's BLAS.
pub fn rot(io: Io, comptime F: type, n: usize, x: [*]F, incx: usize, y: [*]F, incy: usize, c: F, s: F) Error!void {
    if (n == 0) return;
    blas.rot(io, F, @intCast(n), x, @intCast(incx), y, @intCast(incy), c, s) catch |e| return narrow(Error, e);
}

/// `ZROT` (`zrot.f`): `x := c x + s y`, `y := c y - conj(s) x`, `c` real.
pub fn zrot(comptime T: type, n: usize, x: [*]T, incx: usize, y: [*]T, incy: usize, c: Float(T), s: T) void {
    for (0..n) |i| {
        const xi = x[i * incx];
        const yi = y[i * incy];
        x[i * incx] = add(T, mulReal(T, xi, c), mul(T, s, yi));
        y[i * incy] = sub(T, mulReal(T, yi, c), mul(T, conj(T, s), xi));
    }
}

/// `?TRSV` through rstd's BLAS.
pub fn trsv(al: Allocator, io: Io, comptime T: type, upper: bool, op: Trans, unit: bool, n: usize, a: [*]const T, lda: usize, x: [*]T) Error!void {
    blas.trsv(al, io, T, if (upper) .u else .l, op, if (unit) .u else .n, @intCast(n), @constCast(a), @intCast(lda), x, 1) catch |e| return narrow(Error, e);
}
