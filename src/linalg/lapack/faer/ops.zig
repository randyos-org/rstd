//! What the faer ports share: a column-major matrix view, element arithmetic
//! over a real or `Complex` `T`, and the BLAS calls that stand in for faer's
//! own `matmul`, `triangular::matmul` and `triangular_solve`. Level 3 goes
//! through the drivers on the interface's worker rule (as `../lu.zig` calls
//! them, so Accelerate never runs them); level 1 and 2 through the public
//! routines. Every threaded call keeps the BLAS's rule: the same bits for any
//! worker count.

const std = @import("std");
const blas = @import("../../blas/__root__.zig");
const types = @import("../../blas/types.zig");
const util = @import("../../blas/util.zig");
const driver = @import("../../blas/blas3/driver.zig");
const threaded = @import("../../blas/blas3/threaded.zig");
const threading = @import("../../blas/threading.zig");
const chol = @import("../cholesky.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
pub const Float = types.Float;
pub const isComplex = types.isComplex;
pub const Trans = types.Trans;

/// What the internal routines return.
pub const Error = Allocator.Error || Io.Cancelable;

/// Tests only: caps the iterative solvers (`bidiag_svd.qrAlgorithm`, both
/// `lahqr`s) at this many iterations, so `NoConvergence` can be reached.
pub var test_iter_limit: ?usize = null;

/// `limit`, or `test_iter_limit` when set.
pub fn iterLimit(limit: usize) usize {
    return test_iter_limit orelse limit;
}

/// A column-major view: element `(i, j)` at `ptr[i + j * ld]`.
pub fn Mat(comptime T: type) type {
    return struct {
        ptr: [*]T,
        rows: usize,
        cols: usize,
        ld: usize,

        const Self = @This();

        pub fn at(m: Self, i: usize, j: usize) *T {
            return &m.ptr[i + j * m.ld];
        }

        pub fn get(m: Self, i: usize, j: usize) T {
            return m.ptr[i + j * m.ld];
        }

        /// The `r x c` block at `(i, j)`.
        pub fn sub(m: Self, i: usize, j: usize, r: usize, c: usize) Self {
            return .{ .ptr = m.ptr + i + j * m.ld, .rows = r, .cols = c, .ld = m.ld };
        }

        pub fn col(m: Self, j: usize) [*]T {
            return m.ptr + j * m.ld;
        }

        pub fn fill(m: Self, v: T) void {
            for (0..m.cols) |j| @memset(m.col(j)[0..m.rows], v);
        }

        /// Zero, with ones on the diagonal.
        pub fn identity(m: Self) void {
            m.fill(zero(T));
            for (0..@min(m.rows, m.cols)) |i| m.at(i, i).* = one(T);
        }

        pub fn copyFrom(m: Self, src: Self) void {
            for (0..m.cols) |j| @memcpy(m.col(j)[0..m.rows], src.col(j)[0..m.rows]);
        }
    };
}

/// An owned `rows x cols` matrix, `ld = max(rows, 1)`.
pub fn alloc(comptime T: type, al: Allocator, rows: usize, cols: usize) Allocator.Error!Mat(T) {
    const ld = @max(rows, 1);
    const buf = try al.alloc(T, ld * cols);
    return .{ .ptr = buf.ptr, .rows = rows, .cols = cols, .ld = ld };
}

pub fn free(comptime T: type, al: Allocator, m: Mat(T)) void {
    al.free(m.ptr[0 .. m.ld * m.cols]);
}

// ---- element arithmetic -----------------------------------------------------------

pub fn zero(comptime T: type) T {
    return chol.scalar(T, 0);
}

pub fn one(comptime T: type) T {
    return chol.scalar(T, 1);
}

pub const scalar = chol.scalar;

pub fn re(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) x.re else x;
}

pub fn add(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) .init(x.re + y.re, x.im + y.im) else x + y;
}

pub fn sub(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) .init(x.re - y.re, x.im - y.im) else x - y;
}

pub fn neg(comptime T: type, x: T) T {
    return if (comptime isComplex(T)) .init(-x.re, -x.im) else -x;
}

pub fn mul(comptime T: type, x: T, y: T) T {
    if (comptime !isComplex(T)) return x * y;
    return .init(x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re);
}

pub fn conj(comptime T: type, x: T) T {
    return if (comptime isComplex(T)) .init(x.re, -x.im) else x;
}

pub fn mulReal(comptime T: type, x: T, r: Float(T)) T {
    return if (comptime isComplex(T)) .init(x.re * r, x.im * r) else x * r;
}

/// `|x|`, by `hypot` for complex.
pub fn abs(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) std.math.hypot(x.re, x.im) else @abs(x);
}

/// The principal square root. Complex: std's for `f32` and `f64`, which has
/// no other type (issues.md), and the same musl `csqrt` for the rest (T124).
pub fn sqrt(comptime T: type, x: T) T {
    if (comptime !isComplex(T)) return @sqrt(x);
    const F = Float(T);
    if (F == f32 or F == f64) return std.math.complex.sqrt(x);
    // past this, x + hypot(x, y) may overflow
    const threshold = comptime std.math.floatMax(F) / (1 + std.math.sqrt2);
    var xr = x.re;
    var xi = x.im;
    if (xr == 0 and xi == 0) return .init(0, xi);
    if (std.math.isInf(xi)) return .init(std.math.inf(F), xi);
    if (std.math.isNan(xr)) return .init(xr, (xi - xi) / (xi - xi));
    if (std.math.isInf(xr)) {
        return if (std.math.signbit(xr)) .init(@abs(xi - xi), std.math.copysign(xr, xi)) else .init(xr, std.math.copysign(xi - xi, xi));
    }
    const scale = @abs(xr) >= threshold or @abs(xi) >= threshold;
    if (scale) {
        xr *= 0.25;
        xi *= 0.25;
    }
    var r: T = undefined;
    if (xr >= 0) {
        const t = @sqrt((xr + std.math.hypot(xr, xi)) * 0.5);
        r = .init(t, xi / (2 * t));
    } else {
        const t = @sqrt((-xr + std.math.hypot(xr, xi)) * 0.5);
        r = .init(@abs(xi) / (2 * t), std.math.copysign(t, xi));
    }
    return if (scale) .init(2 * r.re, 2 * r.im) else r;
}

/// `|re| + |im|` (faer's `abs1`, the reference's `CABS1`); `|x|` for real.
pub fn abs1(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) @abs(x.re) + @abs(x.im) else @abs(x);
}

pub fn abs2(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) x.re * x.re + x.im * x.im else x * x;
}

/// `1 / x`, complex by Smith's algorithm (`util.cdiv`).
pub fn recip(comptime T: type, x: T) T {
    if (comptime !isComplex(T)) return 1 / x;
    const q = util.cdiv(Float(T), 1, 0, x.re, x.im);
    return .init(q[0], q[1]);
}

pub fn isZero(comptime T: type, x: T) bool {
    return if (comptime isComplex(T)) x.re == 0 and x.im == 0 else x == 0;
}

pub fn isFinite(comptime T: type, x: T) bool {
    return if (comptime isComplex(T)) std.math.isFinite(x.re) and std.math.isFinite(x.im) else std.math.isFinite(x);
}

/// `DLASCL` type `G`: `A := A * (to / from)`, in steps that neither
/// overflow nor underflow (`lapack-netlib/SRC/dlascl.f`). `from` is finite
/// and nonzero.
pub fn lascl(comptime T: type, from: Float(T), to: Float(T), a: Mat(T)) void {
    const F = Float(T);
    const smlnum = std.math.floatMin(F);
    const bignum = 1 / smlnum;
    var cfromc = from;
    var ctoc = to;
    while (true) {
        const cfrom1 = cfromc * smlnum;
        var f: F = undefined;
        var done = true;
        if (cfrom1 == cfromc) {
            f = ctoc / cfromc;
        } else {
            const cto1 = ctoc / bignum;
            if (cto1 == ctoc) {
                f = ctoc;
                cfromc = 1;
            } else if (@abs(cfrom1) > @abs(ctoc) and ctoc != 0) {
                f = smlnum;
                done = false;
                cfromc = cfrom1;
            } else if (@abs(cto1) > @abs(cfromc)) {
                f = bignum;
                done = false;
                ctoc = cto1;
            } else {
                f = ctoc / cfromc;
                if (f == 1) return;
            }
        }
        for (0..a.cols) |j| for (a.col(j)[0..a.rows]) |*x| {
            x.* = mulReal(T, x.*, f);
        };
        if (done) return;
    }
}

/// `xLANGE('M')`: the largest `|a(i, j)|`, NaN if any is NaN.
pub fn normMax(comptime T: type, a: Mat(T)) Float(T) {
    var x: Float(T) = 0;
    for (0..a.cols) |j| for (0..a.rows) |i| {
        const v = abs(T, a.get(i, j));
        if (std.math.isNan(v)) return v;
        x = @max(x, v);
    };
    return x;
}

/// The scale `dgesvd.f`, `dgeev.f` and `dgels.f` bring a matrix of
/// max-norm `anrm` to: `SMLNUM` below it, `BIGNUM` above, else none.
pub fn scaleTo(comptime F: type, anrm: F, smlnum: F) ?F {
    const bignum = 1 / smlnum;
    return if (anrm > 0 and anrm < smlnum) smlnum else if (anrm > bignum) bignum else null;
}

// ---- BLAS ------------------------------------------------------------------------

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
    try driver.gemm(al, io, T, chol.table(T), ta.transposed(), tb.transposed(), .of(ta.conjugated(T), tb.conjugated(T)), c.rows, c.cols, k, alpha, a.ptr, a.ld, b.ptr, b.ld, beta, c.ptr, c.ld, w);
}

/// `B := alpha op(A) B` (`solve`: `alpha op(A)⁻¹ B`) from the left, `A`
/// triangular.
pub fn tri(al: Allocator, io: Io, comptime T: type, comptime solve: bool, upper: bool, op: Op, unit: bool, alpha: T, a: Mat(T), b: Mat(T)) Error!void {
    if (b.rows == 0 or b.cols == 0) return;
    const w = threaded.workers(.trsm, T, threading.product(&.{ b.rows, b.cols }));
    const f = if (solve) driver.trsm else driver.trmm;
    try f(al, io, T, chol.table(T), true, upper, op.transposed(), op.conjugated(T), unit, b.rows, b.cols, alpha, a.ptr, a.ld, b.ptr, b.ld, w);
}

/// `B := alpha B op(A)` (`solve`: `alpha B op(A)⁻¹`) from the right, `A`
/// triangular.
pub fn triRight(al: Allocator, io: Io, comptime T: type, comptime solve: bool, upper: bool, op: Op, unit: bool, alpha: T, a: Mat(T), b: Mat(T)) Error!void {
    if (b.rows == 0 or b.cols == 0) return;
    const w = threaded.workers(.trsm, T, threading.product(&.{ b.rows, b.cols }));
    const f = if (solve) driver.trsm else driver.trmm;
    try f(al, io, T, chol.table(T), false, upper, op.transposed(), op.conjugated(T), unit, b.rows, b.cols, alpha, a.ptr, a.ld, b.ptr, b.ld, w);
}

/// The upper triangle of `C := alpha Aᴴ A + beta C` (real: `Aᵀ A`), `A`
/// being `k x n`; a complex diagonal comes out real.
pub fn gramUpper(al: Allocator, io: Io, comptime T: type, alpha: Float(T), a: Mat(T), beta: Float(T), c: Mat(T)) Error!void {
    const n = a.cols;
    if (n == 0) return;
    const w = threaded.workers(.syrk, T, threading.product(&.{ n + 1, n, a.rows }));
    try driver.syrk(al, io, T, chol.table(T), comptime isComplex(T), true, true, n, a.rows, scalar(T, alpha), a.ptr, a.ld, null, 0, scalar(T, beta), c.ptr, c.ld, w);
}

/// A public level-1 or level-2 routine's error, given arguments it accepts.
fn narrow(e: anyerror) Error {
    return switch (e) {
        error.OutOfMemory => error.OutOfMemory,
        error.Canceled => error.Canceled,
        else => unreachable,
    };
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
    blas.gemv(al, io, T, t, @intCast(a.rows), @intCast(a.cols), alpha, a.ptr, @intCast(a.ld), x, @intCast(incx), beta, y, @intCast(incy)) catch |e| return narrow(e);
}

/// `A += alpha x yᵀ` (`blas.ger`), unconjugated.
pub fn ger(al: Allocator, io: Io, comptime T: type, alpha: T, x: [*]T, incx: usize, y: [*]T, incy: usize, a: Mat(T)) Error!void {
    if (a.rows == 0 or a.cols == 0) return;
    blas.ger(al, io, T, @intCast(a.rows), @intCast(a.cols), alpha, x, @intCast(incx), y, @intCast(incy), a.ptr, @intCast(a.ld)) catch |e| return narrow(e);
}

/// `‖x‖₂` (`blas.nrm2`).
pub fn nrm2(io: Io, comptime T: type, n: usize, x: [*]T, inc: usize) Float(T) {
    if (n == 0) return 0;
    return blas.nrm2(io, T, @intCast(n), x, @intCast(inc)) catch unreachable;
}

test {
    _ = Mat(f64);
}
