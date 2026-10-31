//! Matrix norms and condition estimates, from Reference-LAPACK 3.12.0 as
//! OpenBLAS v0.3.34-317 ships it (`lapack-netlib/SRC/`, modified BSD):
//! `xLANGE` with `xLASSQ`, `xLACN2` (Higham's estimator), `xLATRS`,
//! `xRSCL` and `xGECON`. Neither faer nor libflame has its own: faer's
//! reductions have no operator norms or estimator, libflame's are the
//! reference through f2c.
//!
//! Real and complex share one loop; where `zlatrs.f` departs from
//! `dlatrs.f` (half-scaled bounds, `|re| + |im|` sizes), the complex branch
//! says so. Complex division is Smith's (`util.cdiv`) where the reference
//! takes `ZLADIV`. Vector products go through rstd's BLAS (`dot`, `axpy`,
//! `trsv`, `scal`); the bits do not depend on the worker count.

const std = @import("std");
const blas = @import("../blas/__root__.zig");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const chol = @import("cholesky.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;

/// What the internal routines return.
pub const Error = Allocator.Error || Io.Cancelable;

/// `NORM`: the largest `|a(i, j)|` (`m`), the one norm (`o`, LAPACK's `1`
/// or `O`), the infinity norm (`i`) or the Frobenius norm (`f`, `F` or
/// `E`). `|.|` is the modulus for complex.
pub const Norm = enum { m, o, i, f };

// ---- element helpers ---------------------------------------------------------------

/// `|re| + |im|` (`CABS1`), `|x|` for real.
pub fn abs1(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) @abs(x.re) + @abs(x.im) else @abs(x);
}

/// `x / y`, complex by Smith's algorithm (`util.cdiv`), for `ZLADIV`.
pub fn div(comptime T: type, x: T, y: T) T {
    if (comptime !isComplex(T)) return x / y;
    const q = util.cdiv(Float(T), x.re, x.im, y.re, y.im);
    return .init(q[0], q[1]);
}

/// A public level-1 or level-2 routine's error, given arguments it accepts.
pub fn narrow(e: anyerror) Error {
    return switch (e) {
        error.OutOfMemory => error.OutOfMemory,
        error.Canceled => error.Canceled,
        else => unreachable,
    };
}

/// `?SCAL` (`?DSCAL` for complex) by a real `s`.
pub fn scal(io: Io, comptime T: type, n: usize, s: Float(T), x: [*]T, inc: usize) Error!void {
    if (n == 0) return;
    const r = if (comptime isComplex(T)) blas.scalReal(io, T, @intCast(n), s, x, @intCast(inc)) else blas.scal(io, T, @intCast(n), s, x, @intCast(inc));
    r catch |e| return narrow(e);
}

/// `?SCAL` by a complex `s`.
pub fn scalComplex(io: Io, comptime T: type, n: usize, s: T, x: [*]T, inc: usize) Error!void {
    if (n == 0) return;
    blas.scal(io, T, @intCast(n), s, x, @intCast(inc)) catch |e| return narrow(e);
}

/// `?AXPY`.
pub fn axpy(io: Io, comptime T: type, n: usize, alpha: T, x: [*]const T, incx: usize, y: [*]T, incy: usize) Error!void {
    if (n == 0) return;
    blas.axpy(io, T, @intCast(n), alpha, @constCast(x), @intCast(incx), y, @intCast(incy)) catch |e| return narrow(e);
}

/// `?DOT`, `?DOTU`, or with `conj` `?DOTC` (`conj(x)ᵀ y`).
pub fn dot(io: Io, comptime T: type, conj: bool, n: usize, x: [*]const T, incx: usize, y: [*]const T, incy: usize) Error!T {
    if (n == 0) return ops.zero(T);
    const r = if (comptime isComplex(T)) (if (conj) blas.dotc(io, T, @intCast(n), @constCast(x), @intCast(incx), @constCast(y), @intCast(incy)) else blas.dot(io, T, @intCast(n), @constCast(x), @intCast(incx), @constCast(y), @intCast(incy))) else blas.dot(io, T, @intCast(n), @constCast(x), @intCast(incx), @constCast(y), @intCast(incy));
    return r catch |e| return narrow(e);
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

// ---- xLASET -----------------------------------------------------------------------

/// What `laset` returns.
pub const LasetError = error{ InvalidM, InvalidN, InvalidLda };

/// `?LASET` (`dlaset.f`, `zlaset.f`): the `m x n` `A`'s off-diagonal
/// entries set to `alpha` and its diagonal to `beta`; with `uplo`, only that
/// strict triangle and the diagonal, the rest untouched (null: all of it).
pub fn laset(comptime T: type, uplo: ?types.Uplo, m: isize, n: isize, alpha: T, beta: T, a: [*]T, lda: isize) LasetError!void {
    comptime chol.floatType(T, "laset");
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    const rows: usize = @intCast(m);
    const la: usize = @intCast(lda);
    for (0..@intCast(n)) |j| {
        const from = if (uplo == .l) @min(j, rows) else 0;
        const to = if (uplo == .u) @min(j + 1, rows) else rows;
        for (from..to) |i| a[i + j * la] = if (i == j) beta else alpha;
    }
}

// ---- xLANGE, xLASSQ ----------------------------------------------------------------

/// What `lange` returns.
pub const LangeError = error{ InvalidM, InvalidN, InvalidLda } || Allocator.Error;

/// `?LANGE` (`dlange.f`, `zlange.f`): `norm` of the `m x n` `A`; 0 when
/// empty. A NaN entry gives NaN. The reference checks no argument; here a
/// negative size or short `lda` is an error. `i` takes `m` reals of
/// workspace.
pub fn lange(allocator: Allocator, io: Io, comptime T: type, norm: Norm, m: isize, n: isize, a: [*]const T, lda: isize) LangeError!Float(T) {
    comptime chol.floatType(T, "lange");
    _ = io;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (m == 0 or n == 0) return 0;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    return matNorm(allocator, T, norm, @intCast(m), @intCast(n), a, @intCast(lda));
}

/// `lange` after its checks.
pub fn matNorm(al: Allocator, comptime T: type, norm: Norm, m: usize, n: usize, a: [*]const T, lda: usize) Allocator.Error!Float(T) {
    const F = Float(T);
    if (m == 0 or n == 0) return 0;
    var value: F = 0;
    switch (norm) {
        .m => for (0..n) |j| for (0..m) |i| {
            const t = ops.abs(T, a[i + j * lda]);
            if (value < t or std.math.isNan(t)) value = t;
        },
        .o => for (0..n) |j| {
            var sum: F = 0;
            for (0..m) |i| sum += ops.abs(T, a[i + j * lda]);
            if (value < sum or std.math.isNan(sum)) value = sum;
        },
        .i => {
            const work = try al.alloc(F, m);
            defer al.free(work);
            @memset(work, 0);
            for (0..n) |j| for (0..m) |i| {
                work[i] += ops.abs(T, a[i + j * lda]);
            };
            for (work) |t| if (value < t or std.math.isNan(t)) {
                value = t;
            };
        },
        .f => {
            var scale: F = 0;
            var sum: F = 1;
            for (0..n) |j| lassq(T, m, a + j * lda, 1, &scale, &sum);
            value = scale * @sqrt(sum);
        },
    }
    return value;
}

/// Blue's thresholds and scalings (`la_constants.f90`'s `tsml`, `tbig`,
/// `ssml`, `sbig`) from Fortran's `digits`, `minexponent`, `maxexponent`.
fn Blue(comptime F: type) type {
    return struct {
        const digits = std.math.floatFractionalBits(F) + 1;
        const minexp = std.math.floatExponentMin(F) + 1;
        const maxexp = std.math.floatExponentMax(F) + 1;
        const tsml = std.math.ldexp(@as(F, 1), -@divFloor(1 - minexp, 2));
        const tbig = std.math.ldexp(@as(F, 1), @divFloor(maxexp - digits + 1, 2));
        const ssml = std.math.ldexp(@as(F, 1), -@divFloor(minexp - digits, 2));
        const sbig = std.math.ldexp(@as(F, 1), -(-@divFloor(-(maxexp + digits - 1), 2)));
    };
}

/// `?LASSQ` (`dlassq.f90`, `zlassq.f90`, Blue's algorithm): updates
/// `scale` and `sumsq` so `scale² sumsq` gains `Σ |x(i)|²` (complex: both
/// parts' squares). A NaN in either leaves both as they are.
pub fn lassq(comptime T: type, n: usize, x: [*]const T, inc: usize, scale: *Float(T), sumsq: *Float(T)) void {
    const F = Float(T);
    const c = Blue(F);
    if (std.math.isNan(scale.*) or std.math.isNan(sumsq.*)) return;
    if (sumsq.* == 0) scale.* = 1;
    if (scale.* == 0) {
        scale.* = 1;
        sumsq.* = 0;
    }
    if (n == 0) return;
    var notbig = true;
    var asml: F = 0;
    var amed: F = 0;
    var abig: F = 0;
    for (0..n) |i| {
        const v = x[i * inc];
        const parts: []const F = if (comptime isComplex(T)) &.{ v.re, v.im } else &.{v};
        for (parts) |p| {
            const ax = @abs(p);
            if (ax > c.tbig) {
                abig += (ax * c.sbig) * (ax * c.sbig);
                notbig = false;
            } else if (ax < c.tsml) {
                if (notbig) asml += (ax * c.ssml) * (ax * c.ssml);
            } else {
                amed += ax * ax;
            }
        }
    }
    // the sum so far joins the accumulator of its size
    if (sumsq.* > 0) {
        const ax = scale.* * @sqrt(sumsq.*);
        if (ax > c.tbig) {
            if (scale.* > 1) {
                scale.* *= c.sbig;
                abig += scale.* * (scale.* * sumsq.*);
            } else {
                abig += scale.* * (scale.* * (c.sbig * (c.sbig * sumsq.*)));
            }
        } else if (ax < c.tsml) {
            if (notbig) {
                if (scale.* < 1) {
                    scale.* *= c.ssml;
                    asml += scale.* * (scale.* * sumsq.*);
                } else {
                    asml += scale.* * (scale.* * (c.ssml * (c.ssml * sumsq.*)));
                }
            }
        } else {
            amed += scale.* * (scale.* * sumsq.*);
        }
    }
    if (abig > 0) {
        if (amed > 0 or std.math.isNan(amed)) abig += (amed * c.sbig) * c.sbig;
        scale.* = 1 / c.sbig;
        sumsq.* = abig;
    } else if (asml > 0) {
        if (amed > 0 or std.math.isNan(amed)) {
            const m = @sqrt(amed);
            const s = @sqrt(asml) / c.ssml;
            const ymin = if (s > m) m else s;
            const ymax = if (s > m) s else m;
            scale.* = 1;
            sumsq.* = ymax * ymax * (1 + (ymin / ymax) * (ymin / ymax));
        } else {
            scale.* = 1 / c.ssml;
            sumsq.* = asml;
        }
    } else {
        scale.* = 1;
        sumsq.* = amed;
    }
}

// ---- xLACN2 ------------------------------------------------------------------------

/// `ITMAX` (`dlacn2.f`): the most power-method steps.
const itmax = 5;

/// `?LACN2` (`dlacn2.f`, `zlacn2.f`; N. J. Higham, "FORTRAN codes for
/// estimating the one-norm of a real or complex matrix", ACM TOMS 14(4),
/// 1988): a lower bound on `‖B‖₁`, with `v` (`n`) holding `B w` for the
/// `w` it found, `‖v‖₁ = est`. The reference's reverse communication is a
/// callback: `apply(ctx, adjoint, x)` sets `x := B x` (its `KASE = 1`) or,
/// with `adjoint`, `x := Bᴴ x` (`KASE = 2`); any error it returns ends the
/// estimate. `x` (`n`) is workspace.
pub fn lacn2(al: Allocator, comptime T: type, comptime E: type, n: usize, v: [*]T, x: [*]T, ctx: anytype, comptime apply: fn (@TypeOf(ctx), bool, [*]T) E!void) (E || Allocator.Error)!Float(T) {
    const F = Float(T);
    const nf: F = @floatFromInt(n);
    const isgn = try al.alloc(bool, if (comptime isComplex(T)) 0 else n);
    defer al.free(isgn);
    for (0..n) |i| x[i] = ops.scalar(T, 1 / nf);
    try apply(ctx, false, x);
    if (n == 1) {
        v[0] = x[0];
        return ops.abs(T, v[0]);
    }
    var est = sum1(T, n, x);
    signs(T, n, x, isgn);
    try apply(ctx, true, x);
    var j = imax1(T, n, x);
    var iter: usize = 2;
    while (true) {
        @memset(x[0..n], ops.zero(T));
        x[j] = ops.one(T);
        try apply(ctx, false, x);
        @memcpy(v[0..n], x[0..n]);
        const estold = est;
        est = sum1(T, n, v);
        if (comptime !isComplex(T)) {
            // repeated signs: converged
            for (0..n) |i| {
                if ((x[i] >= 0) != isgn[i]) break;
            } else break;
        }
        if (est <= estold) break;
        signs(T, n, x, isgn);
        try apply(ctx, true, x);
        const jlast = j;
        j = imax1(T, n, x);
        const last = if (comptime isComplex(T)) ops.abs(T, x[jlast]) else x[jlast];
        if (last != ops.abs(T, x[j]) and iter < itmax) {
            iter += 1;
            continue;
        }
        break;
    }
    // the alternating-sign test vector
    var altsgn: F = 1;
    for (0..n) |i| {
        x[i] = ops.scalar(T, altsgn * (1 + @as(F, @floatFromInt(i)) / (nf - 1)));
        altsgn = -altsgn;
    }
    try apply(ctx, false, x);
    const temp = 2 * (sum1(T, n, x) / (3 * nf));
    if (temp > est) {
        @memcpy(v[0..n], x[0..n]);
        est = temp;
    }
    return est;
}

/// `DASUM`, or `DZSUM1` (moduli) for complex.
fn sum1(comptime T: type, n: usize, x: [*]const T) Float(T) {
    var s: Float(T) = 0;
    for (0..n) |i| s += ops.abs(T, x[i]);
    return s;
}

/// `IDAMAX`, or `IZMAX1` (by modulus) for complex.
fn imax1(comptime T: type, n: usize, x: [*]const T) usize {
    var best: usize = 0;
    var m: Float(T) = -1;
    for (0..n) |i| {
        const a = ops.abs(T, x[i]);
        if (a > m) {
            m = a;
            best = i;
        }
    }
    return best;
}

/// `x := sign(x)`: real `±1` (NaN: -1), kept in `isgn`; complex `x / |x|`,
/// or 1 at or below the safe minimum.
fn signs(comptime T: type, n: usize, x: [*]T, isgn: []bool) void {
    const F = Float(T);
    for (0..n) |i| {
        if (comptime isComplex(T)) {
            const a = ops.abs(T, x[i]);
            x[i] = if (a > std.math.floatMin(F)) .init(x[i].re / a, x[i].im / a) else ops.one(T);
        } else {
            isgn[i] = x[i] >= 0;
            x[i] = if (isgn[i]) 1 else -1;
        }
    }
}

// ---- xRSCL, xLATRS -----------------------------------------------------------------

/// `?RSCL`, `?DRSCL` (`drscl.f`, `zdrscl.f`): `x := x / sa` by safe
/// multiplications, `sa` real, never forming an overflowing `1 / sa`.
pub fn rscl(io: Io, comptime T: type, n: usize, sa: Float(T), x: [*]T, inc: usize) Error!void {
    const F = Float(T);
    if (n == 0) return;
    const smlnum = std.math.floatMin(F);
    const bignum = 1 / smlnum;
    var cden = sa;
    var cnum: F = 1;
    while (true) {
        const cden1 = cden * smlnum;
        const cnum1 = cnum / bignum;
        var mul: F = undefined;
        var done = false;
        if (@abs(cden1) > @abs(cnum) and cnum != 0) {
            mul = smlnum;
            cden = cden1;
        } else if (@abs(cnum1) > @abs(cden)) {
            mul = bignum;
            cnum = cnum1;
        } else {
            mul = cnum / cden;
            done = true;
        }
        try scal(io, T, n, mul, x, inc);
        if (done) return;
    }
}

/// `TRANS` for `latrs`: `n`, `t`, or `c` (complex only).
pub const Op = enum { n, t, c };

/// `?LATRS` (`dlatrs.f`, `zlatrs.f`): solves `op(A) x = scale b` for the
/// `n x n` triangular `A` with `scale <= 1` chosen so `x` cannot overflow;
/// returns `scale`. `cnorm` (`n`) holds the off-diagonal column norms
/// (`abs1` sums), computed unless `normin`. A zero on the diagonal gives a
/// null vector with `scale = 0`. `trsv` takes over where the bounds allow.
pub fn latrs(al: Allocator, io: Io, comptime T: type, upper: bool, op: Op, unit: bool, normin: bool, n: usize, a: [*]const T, lda: usize, x: [*]T, cnorm: [*]Float(T)) Error!Float(T) {
    const F = Float(T);
    const cx = comptime isComplex(T);
    // zlatrs halves its bounds, so a sum of two parts cannot overflow
    const half: F = if (cx) 0.5 else 1;
    var scale: F = 1;
    if (n == 0) return scale;
    const notran = op == .n;
    const smlnum = std.math.floatMin(F) / std.math.floatEps(F);
    const bignum = 1 / smlnum;
    const overflow = std.math.floatMax(F);
    const at = struct {
        fn f(p: [*]const T, ld: usize, i: usize, j: usize) T {
            return p[i + j * ld];
        }
    }.f;
    if (!normin) {
        for (0..n) |j| {
            const lo, const len = if (upper) .{ @as(usize, 0), j } else .{ j + 1, n - j - 1 };
            var s: F = 0;
            for (0..len) |i| s += abs1(T, at(a, lda, lo + i, j));
            cnorm[j] = s;
        }
    }
    // scale the column norms by `tscal` if the largest exceeds `bignum`
    var tmax = cnorm[iamax(F, n, cnorm, 1)];
    var tscal: F = 1;
    if (tmax > bignum * half) {
        if (tmax <= overflow) {
            tscal = half / (smlnum * tmax);
            for (0..n) |j| cnorm[j] *= tscal;
        } else {
            // a column's sum overflowed: rescale from its largest entry
            tmax = 0;
            for (0..n) |j| {
                const lo, const len = if (upper) .{ @as(usize, 0), j } else .{ j + 1, n - j - 1 };
                for (0..len) |i| {
                    const v = at(a, lda, lo + i, j);
                    tmax = if (cx) @max(tmax, @abs(v.re), @abs(v.im)) else @max(tmax, @abs(v));
                }
            }
            if (tmax <= overflow) {
                tscal = 1 / (smlnum * tmax);
                for (0..n) |j| {
                    if (cnorm[j] <= overflow) {
                        cnorm[j] *= tscal;
                        continue;
                    }
                    if (cx) tscal *= 2;
                    cnorm[j] = 0;
                    const lo, const len = if (upper) .{ @as(usize, 0), j } else .{ j + 1, n - j - 1 };
                    for (0..len) |i| cnorm[j] += tscal * abs2(T, at(a, lda, lo + i, j));
                    if (cx) tscal *= 0.5;
                }
            } else {
                try trsv(al, io, T, upper, op, unit, n, a, lda, x);
                return scale;
            }
        }
    }
    // bound the growth of the solution
    var xmax: F = if (cx) blk: {
        var m: F = 0;
        for (0..n) |j| m = @max(m, abs2(T, x[j]));
        break :blk m;
    } else @abs(x[iamax(T, n, x, 1)]);
    var xbnd = xmax;
    // the order of `j`: forward for lower `A` (upper `op(A)`)
    const forward = upper != notran;
    const step = struct {
        fn idx(fwd: bool, nn: usize, k: usize) usize {
            return if (fwd) k else nn - 1 - k;
        }
    }.idx;
    var grow: F = 0;
    if (tscal == 1) grow: {
        if (notran) {
            if (!unit) {
                grow = half / @max(xbnd, smlnum);
                xbnd = grow;
                for (0..n) |k| {
                    const j = step(forward, n, k);
                    if (grow <= smlnum) break :grow;
                    const tjj = abs1(T, at(a, lda, j, j));
                    if (!cx or tjj >= smlnum) xbnd = @min(xbnd, @min(1, tjj) * grow) else xbnd = 0;
                    grow = if (tjj + cnorm[j] >= smlnum) grow * (tjj / (tjj + cnorm[j])) else 0;
                }
                grow = xbnd;
            } else {
                grow = @min(1, half / @max(xbnd, smlnum));
                for (0..n) |k| {
                    const j = step(forward, n, k);
                    if (grow <= smlnum) break :grow;
                    grow *= 1 / (1 + cnorm[j]);
                }
            }
        } else {
            if (!unit) {
                grow = half / @max(xbnd, smlnum);
                xbnd = grow;
                for (0..n) |k| {
                    const j = step(forward, n, k);
                    if (grow <= smlnum) break :grow;
                    const xj = 1 + cnorm[j];
                    grow = @min(grow, xbnd / xj);
                    const tjj = abs1(T, at(a, lda, j, j));
                    if (!cx or tjj >= smlnum) {
                        if (xj > tjj) xbnd *= tjj / xj;
                    } else xbnd = 0;
                }
                grow = @min(grow, xbnd);
            } else {
                grow = @min(1, half / @max(xbnd, smlnum));
                for (0..n) |k| {
                    const j = step(forward, n, k);
                    if (grow <= smlnum) break :grow;
                    grow /= 1 + cnorm[j];
                }
            }
        }
    }
    if (grow * tscal > smlnum) {
        // no overflow: the BLAS's solve
        try trsv(al, io, T, upper, op, unit, n, a, lda, x);
    } else {
        if (xmax > bignum * half) {
            scale = (bignum * half) / xmax;
            try scal(io, T, n, scale, x, 1);
            xmax = bignum;
        } else if (cx) {
            xmax *= 2;
        }
        const conj = cx and op == .c;
        for (0..n) |k| {
            const j = step(forward, n, k);
            if (notran) {
                // x(j) / A(j, j), then the update of the rest by column j
                var xj = abs1(T, x[j]);
                skip: {
                    const tjjs = if (!unit) ops.mulReal(T, at(a, lda, j, j), tscal) else blk: {
                        if (tscal == 1) break :skip;
                        break :blk ops.scalar(T, tscal);
                    };
                    const tjj = abs1(T, tjjs);
                    if (tjj > smlnum) {
                        if (tjj < 1 and xj > tjj * bignum) {
                            const rec = 1 / xj;
                            try scal(io, T, n, rec, x, 1);
                            scale *= rec;
                            xmax *= rec;
                        }
                        x[j] = div(T, x[j], tjjs);
                        xj = abs1(T, x[j]);
                    } else if (tjj > 0) {
                        if (xj > tjj * bignum) {
                            var rec = (tjj * bignum) / xj;
                            if (cnorm[j] > 1) rec /= cnorm[j];
                            try scal(io, T, n, rec, x, 1);
                            scale *= rec;
                            xmax *= rec;
                        }
                        x[j] = div(T, x[j], tjjs);
                        xj = abs1(T, x[j]);
                    } else {
                        @memset(x[0..n], ops.zero(T));
                        x[j] = ops.one(T);
                        xj = 1;
                        scale = 0;
                        xmax = 0;
                    }
                }
                if (xj > 1) {
                    var rec = 1 / xj;
                    if (cnorm[j] > (bignum - xmax) * rec) {
                        rec *= 0.5;
                        try scal(io, T, n, rec, x, 1);
                        scale *= rec;
                    }
                } else if (xj * cnorm[j] > bignum - xmax) {
                    try scal(io, T, n, 0.5, x, 1);
                    scale *= 0.5;
                }
                const alpha = ops.neg(T, ops.mulReal(T, x[j], tscal));
                if (upper) {
                    if (j > 0) {
                        try axpy(io, T, j, alpha, a + j * lda, 1, x, 1);
                        xmax = abs1(T, x[iamax(T, j, x, 1)]);
                    }
                } else if (j + 1 < n) {
                    try axpy(io, T, n - j - 1, alpha, a + j + 1 + j * lda, 1, x + j + 1, 1);
                    xmax = abs1(T, x[j + 1 + iamax(T, n - j - 1, x + j + 1, 1)]);
                }
            } else {
                // x(j) - op(A)(j, :) x, then the division
                var xj = abs1(T, x[j]);
                var uscal = ops.scalar(T, tscal);
                var rec = 1 / @max(xmax, 1);
                const diag = struct {
                    fn f(p: [*]const T, ld: usize, jj: usize, u: bool, c: bool, ts: F) T {
                        if (u) return ops.scalar(T, ts);
                        const d = p[jj + jj * ld];
                        return ops.mulReal(T, if (c) ops.conj(T, d) else d, ts);
                    }
                }.f;
                if (cnorm[j] > (bignum - xj) * rec) {
                    rec *= 0.5;
                    const tjjs = diag(a, lda, j, unit, conj, tscal);
                    const tjj = abs1(T, tjjs);
                    if (tjj > 1) {
                        rec = @min(1, rec * tjj);
                        uscal = div(T, uscal, tjjs);
                    }
                    if (rec < 1) {
                        try scal(io, T, n, rec, x, 1);
                        scale *= rec;
                        xmax *= rec;
                    }
                }
                const lo, const len = if (upper) .{ @as(usize, 0), j } else .{ j + 1, n - j - 1 };
                var sumj = ops.zero(T);
                if (eq(T, uscal, ops.one(T))) {
                    sumj = try dot(io, T, conj, len, a + lo + j * lda, 1, x + lo, 1);
                } else {
                    for (0..len) |i| {
                        const aij = at(a, lda, lo + i, j);
                        sumj = ops.add(T, sumj, ops.mul(T, ops.mul(T, if (conj) ops.conj(T, aij) else aij, uscal), x[lo + i]));
                    }
                }
                if (eq(T, uscal, ops.scalar(T, tscal))) {
                    x[j] = ops.sub(T, x[j], sumj);
                    xj = abs1(T, x[j]);
                    skip: {
                        if (unit and tscal == 1) break :skip;
                        const tjjs = diag(a, lda, j, unit, conj, tscal);
                        const tjj = abs1(T, tjjs);
                        if (tjj > smlnum) {
                            if (tjj < 1 and xj > tjj * bignum) {
                                const r = 1 / xj;
                                try scal(io, T, n, r, x, 1);
                                scale *= r;
                                xmax *= r;
                            }
                            x[j] = div(T, x[j], tjjs);
                        } else if (tjj > 0) {
                            if (xj > tjj * bignum) {
                                const r = (tjj * bignum) / xj;
                                try scal(io, T, n, r, x, 1);
                                scale *= r;
                                xmax *= r;
                            }
                            x[j] = div(T, x[j], tjjs);
                        } else {
                            @memset(x[0..n], ops.zero(T));
                            x[j] = ops.one(T);
                            scale = 0;
                            xmax = 0;
                        }
                    }
                } else {
                    x[j] = ops.sub(T, div(T, x[j], diag(a, lda, j, unit, conj, tscal)), sumj);
                }
                xmax = @max(xmax, abs1(T, x[j]));
            }
        }
        scale /= tscal;
    }
    if (tscal != 1) for (0..n) |j| {
        cnorm[j] *= 1 / tscal;
    };
    return scale;
}

/// `CABS2` (`zlatrs.f`): `|re / 2| + |im / 2|`; `|x|` for real.
fn abs2(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) @abs(x.re / 2) + @abs(x.im / 2) else @abs(x);
}

fn eq(comptime T: type, x: T, y: T) bool {
    return if (comptime isComplex(T)) x.re == y.re and x.im == y.im else x == y;
}

/// `?TRSV` through rstd's BLAS.
fn trsv(al: Allocator, io: Io, comptime T: type, upper: bool, op: Op, unit: bool, n: usize, a: [*]const T, lda: usize, x: [*]T) Error!void {
    const t: Trans = switch (op) {
        .n => .n,
        .t => .t,
        .c => .c,
    };
    blas.trsv(al, io, T, if (upper) .u else .l, t, if (unit) .u else .n, @intCast(n), @constCast(a), @intCast(lda), x, 1) catch |e| return narrow(e);
}

// ---- xGECON ------------------------------------------------------------------------

/// What `gecon` returns: `InvalidNorm` for other than `o` or `i`;
/// `InvalidAnorm` for a negative, infinite or NaN `anorm` (the reference's
/// `INFO = -5`, the last two after its checks); `EstimateFailed` for the
/// reference's `INFO = 1`, an estimate that is NaN, infinite, or of a zero
/// inverse norm.
pub const GeconError = error{ InvalidNorm, InvalidN, InvalidLda, InvalidAnorm, EstimateFailed } || Allocator.Error || Io.Cancelable;

/// `?GECON` (`dgecon.f`, `zgecon.f`): an estimate of the reciprocal
/// condition number `1 / (‖A‖ ‖A⁻¹‖)` in the one (`o`) or infinity (`i`)
/// norm, from `getrf`'s `L U` in `A` (`n x n`) and `anorm`, `A`'s norm
/// before factoring. `‖A⁻¹‖` comes from `lacn2` over solves with `L` and
/// `U` (`latrs`). Written to `rcond`: 1 for `n = 0`, 0 for `anorm = 0` or
/// a solve that would overflow, the estimate otherwise; NaN `anorm` writes
/// NaN and returns `InvalidAnorm`, and `EstimateFailed` leaves the failed
/// value. Takes `5 n` elements of workspace.
pub fn gecon(allocator: Allocator, io: Io, comptime T: type, norm: Norm, n: isize, a: [*]const T, lda: isize, anorm: Float(T), rcond: *Float(T)) GeconError!void {
    comptime chol.floatType(T, "gecon");
    const F = Float(T);
    const onenrm = switch (norm) {
        .o => true,
        .i => false,
        .m, .f => return error.InvalidNorm,
    };
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (anorm < 0) return error.InvalidAnorm;
    rcond.* = 0;
    if (n == 0) {
        rcond.* = 1;
        return;
    }
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (anorm == 0) return;
    if (std.math.isNan(anorm)) {
        rcond.* = anorm;
        return error.InvalidAnorm;
    }
    if (anorm > std.math.floatMax(F)) return error.InvalidAnorm;
    const un: usize = @intCast(n);
    const work = try allocator.alloc(T, 2 * un);
    defer allocator.free(work);
    const cnorm = try allocator.alloc(F, 2 * un);
    defer allocator.free(cnorm);
    var ctx: Gecon(T) = .{ .al = allocator, .io = io, .onenrm = onenrm, .n = un, .a = a, .lda = @intCast(lda), .cnorm = cnorm.ptr };
    const ainvnm = lacn2(allocator, T, Gecon(T).Stop, un, work.ptr + un, work.ptr, &ctx, Gecon(T).apply) catch |e| switch (e) {
        error.Overflow => return,
        error.OutOfMemory, error.Canceled => |x| return x,
    };
    if (ainvnm == 0) return error.EstimateFailed;
    rcond.* = (1 / ainvnm) / anorm;
    if (std.math.isNan(rcond.*) or rcond.* > std.math.floatMax(F)) return error.EstimateFailed;
}

/// `gecon`'s `B = A⁻¹` for `lacn2`: two `latrs` solves, `L` then `U` (or
/// their transposes in reverse), and the reference's early exit (`GO TO
/// 20`) as `error.Overflow`.
fn Gecon(comptime T: type) type {
    return struct {
        al: Allocator,
        io: Io,
        onenrm: bool,
        n: usize,
        a: [*]const T,
        lda: usize,
        cnorm: [*]Float(T),
        normin: bool = false,

        const Self = @This();
        const Stop = Error || error{Overflow};

        fn apply(s: *Self, adjoint: bool, x: [*]T) Stop!void {
            const F = Float(T);
            const op: Op = if (comptime isComplex(T)) .c else .t;
            const lower_norms = s.cnorm;
            const upper_norms = s.cnorm + s.n;
            var sl: F = undefined;
            var su: F = undefined;
            if (adjoint != s.onenrm) {
                sl = try latrs(s.al, s.io, T, false, .n, true, s.normin, s.n, s.a, s.lda, x, lower_norms);
                su = try latrs(s.al, s.io, T, true, .n, false, s.normin, s.n, s.a, s.lda, x, upper_norms);
            } else {
                su = try latrs(s.al, s.io, T, true, op, false, s.normin, s.n, s.a, s.lda, x, upper_norms);
                sl = try latrs(s.al, s.io, T, false, op, true, s.normin, s.n, s.a, s.lda, x, lower_norms);
            }
            const scale = sl * su;
            s.normin = true;
            if (scale != 1) {
                const ix = iamax(T, s.n, x, 1);
                if (scale < abs1(T, x[ix]) * std.math.floatMin(F) or scale == 0) return error.Overflow;
                try rscl(s.io, T, s.n, scale, x, 1);
            }
        }
    };
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/blas1/check.zig");
const lu = @import("lu.zig");

const test_types = chol.test_types;

/// `THRESH` of the reference's linear-equation tests (`TESTING/dtest.in`).
pub const thresh = 30;

/// `T` from parts.
pub fn cplx(comptime T: type, r: f64, i: f64) T {
    return if (comptime isComplex(T)) .init(@floatCast(r), @floatCast(i)) else @floatCast(r);
}

test "lapack norm: lange's four norms on a known 2 x 3, NaN, empty, bad arguments" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        // lda 3, columns (1, -2), (3 + 4i, 0), (-6, 0); complex (3, 4) has modulus 5
        const a = [_]T{ cplx(T, 1, 0), cplx(T, -2, 0), cplx(T, 9, 9), cplx(T, 3, if (isComplex(T)) 4 else 0), cplx(T, 0, 0), cplx(T, 9, 9), cplx(T, -6, 0), cplx(T, 0, 0) };
        const c2: F = if (comptime isComplex(T)) 5 else 3;
        try testing.expectEqual(@as(F, 6), try lange(al, check.io, T, .m, 2, 3, &a, 3));
        try testing.expectEqual(@as(F, 6), try lange(al, check.io, T, .o, 2, 3, &a, 3));
        try testing.expectEqual(1 + c2 + 6, try lange(al, check.io, T, .i, 2, 3, &a, 3));
        const fro = try lange(al, check.io, T, .f, 2, 3, &a, 3);
        try testing.expectApproxEqRel(@sqrt(41 + c2 * c2), fro, 4 * std.math.floatEps(F));
        var b = a;
        b[1] = chol.poisoned(T);
        for ([_]Norm{ .m, .o, .i, .f }) |nm| try testing.expect(std.math.isNan(try lange(al, check.io, T, nm, 2, 3, &b, 3)));
        try testing.expectEqual(@as(F, 0), try lange(testing.failing_allocator, check.io, T, .i, 0, 3, &b, 1));
        try testing.expectEqual(@as(F, 0), try lange(testing.failing_allocator, check.io, T, .f, 2, 0, &b, 2));
    }
    var x = [_]f64{0};
    try testing.expectError(error.InvalidM, lange(al, check.io, f64, .m, -1, -1, &x, 0));
    try testing.expectError(error.InvalidN, lange(al, check.io, f64, .m, 1, -1, &x, 0));
    try testing.expectError(error.InvalidLda, lange(al, check.io, f64, .m, 2, 1, &x, 1));
}

test "lapack norm: lange's Frobenius norm neither overflows nor underflows (Blue's lassq)" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        // sqrt(20) big stays finite
        const big = std.math.floatMax(F) / 16;
        const tiny = std.math.floatTrueMin(F) * 16;
        for ([_]F{ big, tiny, 1 }) |s| {
            var a: [8]T = undefined;
            for (&a, 0..) |*v, i| v.* = ops.scalar(T, s * @as(F, @floatFromInt(i % 2 + 1)));
            // 4 (1 + 4) s squared
            const want = @sqrt(@as(F, 20)) * s;
            const got = try lange(al, check.io, T, .f, 4, 2, &a, 4);
            try testing.expect(std.math.isFinite(got) and got > 0);
            // denormal entries carry few bits
            const tol: F = if (s == tiny) 0.15 else 8 * std.math.floatEps(F);
            try testing.expectApproxEqRel(want, got, tol);
        }
        // a big and a small entry together
        var c = [_]T{ ops.scalar(T, big), ops.scalar(T, tiny), ops.scalar(T, 1) };
        try testing.expectEqual(big, try lange(al, check.io, T, .f, 3, 1, &c, 3));
    }
}

/// `DGET06` (`TESTING/LIN/dget06.f`): how far an estimate `rcond` is from
/// the true `rcondc`, as a ratio less `1 - eps`.
pub fn ratio(comptime F: type, rcond: F, rcondc: F) F {
    const eps = std.math.floatEps(F);
    if (rcond > 0) return if (rcondc > 0) @max(rcond, rcondc) / @min(rcond, rcondc) - (1 - eps) else rcond / eps;
    return if (rcondc > 0) rcondc / eps else 0;
}

/// `gecon`'s estimate against `1 / (‖A‖ ‖A⁻¹‖)` from `getrs`'s inverse,
/// as the reference's `dchkge.f` checks it, in both norms.
fn geconCase(comptime T: type, al: Allocator, n: usize, grade: Float(T), seed: u64) !void {
    const F = Float(T);
    const a = try al.alloc(T, n * n);
    defer al.free(a);
    check.fill(T, a, seed);
    // row i scaled by grade^i
    for (0..n) |i| {
        const s: F = @floatCast(std.math.pow(f64, @floatCast(grade), @floatFromInt(i)));
        for (0..n) |j| a[i + j * n] = ops.mulReal(T, a[i + j * n], s);
    }
    const ipiv = try al.alloc(isize, n);
    defer al.free(ipiv);
    const inv = try al.alloc(T, n * n);
    defer al.free(inv);
    const f = try al.alloc(T, n * n);
    defer al.free(f);
    const ni: isize = @intCast(n);
    @memcpy(f, a);
    try lu.getrf(al, check.io, T, ni, ni, f.ptr, ni, ipiv.ptr, null);
    @memset(inv, ops.zero(T));
    for (0..n) |i| inv[i * (n + 1)] = ops.one(T);
    try lu.getrs(al, check.io, T, .n, ni, ni, f.ptr, ni, ipiv.ptr, inv.ptr, ni);
    for ([_]Norm{ .o, .i }) |nm| {
        const anorm = try lange(al, check.io, T, nm, ni, ni, a.ptr, ni);
        const ainvnm = try lange(al, check.io, T, nm, ni, ni, inv.ptr, ni);
        const rcondc = (1 / anorm) / ainvnm;
        var rcond: F = undefined;
        try gecon(al, check.io, T, nm, ni, f.ptr, ni, anorm, &rcond);
        const r = ratio(F, rcond, rcondc);
        if (!(r < thresh)) {
            check.say("gecon {s} n {d} grade {d} norm {s}: rcond {e} against {e}, ratio {d}\n", .{ @typeName(T), n, @as(f64, @floatCast(grade)), @tagName(nm), @as(f64, @floatCast(rcond)), @as(f64, @floatCast(rcondc)), @as(f64, @floatCast(r)) });
            return error.TestUnexpectedResult;
        }
    }
}

test "lapack gecon: estimates within the reference's threshold of 1 / (|A| |inv(A)|), both norms" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 3, 5, 10, 50 }) |n| {
            // the last grade reaches a condition near 1 / sqrt(eps)
            const steep: F = @floatCast(std.math.pow(f64, @floatCast(@sqrt(std.math.floatEps(F))), 1 / @as(f64, @floatFromInt(@max(n, 2) - 1))));
            for ([_]F{ 1, 0.5, steep }, 0..) |g, k| try geconCase(T, al, n, g, n * 7 + k);
        }
    }
}

test "lapack gecon: n = 0, anorm = 0, 1 x 1, a zero pivot, NaN and infinite anorm, an overflowing estimate, bad arguments" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        var rcond: F = -1;
        var a = [_]T{ cplx(T, 2, 0), cplx(T, 1, 0), cplx(T, 0, 0), cplx(T, 3, 0) };
        try gecon(al, check.io, T, .o, 0, &a, 1, 1, &rcond);
        try testing.expectEqual(@as(F, 1), rcond);
        try gecon(al, check.io, T, .o, 2, &a, 2, 0, &rcond);
        try testing.expectEqual(@as(F, 0), rcond);
        // 1 x 1: exact
        try gecon(al, check.io, T, .i, 1, &a, 1, 2, &rcond);
        try testing.expectEqual(@as(F, 1), rcond);
        // U singular: the overflowing solve stops at 0
        var s = [_]T{ cplx(T, 2, 0), cplx(T, 0.5, 0), cplx(T, 1, 0), cplx(T, 0, 0) };
        for ([_]Norm{ .o, .i }) |nm| {
            rcond = -1;
            try gecon(al, check.io, T, nm, 2, &s, 2, 3, &rcond);
            try testing.expectEqual(@as(F, 0), rcond);
        }
        try testing.expectError(error.InvalidAnorm, gecon(al, check.io, T, .o, 2, &a, 2, std.math.nan(F), &rcond));
        try testing.expect(std.math.isNan(rcond));
        try testing.expectError(error.InvalidAnorm, gecon(al, check.io, T, .o, 2, &a, 2, std.math.inf(F), &rcond));
        try testing.expectEqual(@as(F, 0), rcond);
        // an anorm so small the estimate overflows: the reference's INFO = 1
        try testing.expectError(error.EstimateFailed, gecon(al, check.io, T, .o, 1, &a, 1, std.math.floatTrueMin(F), &rcond));
        try testing.expect(std.math.isInf(rcond));
    }
    var x = [_]f64{1};
    var r: f64 = 0;
    try testing.expectError(error.InvalidNorm, gecon(al, check.io, f64, .f, -1, &x, 0, -1, &r));
    try testing.expectError(error.InvalidN, gecon(al, check.io, f64, .o, -1, &x, 0, -1, &r));
    try testing.expectError(error.InvalidLda, gecon(al, check.io, f64, .o, 2, &x, 1, -1, &r));
    try testing.expectError(error.InvalidAnorm, gecon(al, check.io, f64, .o, 1, &x, 1, -1, &r));
}

/// `latrs` on an `n x n` triangle with diagonal `d`: `scale` in [0, 1], `x`
/// finite, and `‖op(A) x - scale b‖ / (‖A‖ ‖x‖ n eps)` below the threshold.
fn latrsCase(comptime T: type, al: Allocator, upper: bool, op: Op, unit: bool, n: usize, d: Float(T), seed: u64) !void {
    const F = Float(T);
    const a = try al.alloc(T, n * n);
    defer al.free(a);
    check.fill(T, a, seed);
    for (0..n) |j| for (0..n) |i| {
        if (if (upper) i > j else i < j) a[i + j * n] = chol.poisoned(T);
    };
    for (0..n) |i| a[i * (n + 1)] = if (unit) chol.poisoned(T) else ops.scalar(T, d);
    const b = try al.alloc(T, n);
    defer al.free(b);
    check.fill(T, b, seed + 1);
    const x = try al.alloc(T, n);
    defer al.free(x);
    @memcpy(x, b);
    const cnorm = try al.alloc(F, n);
    defer al.free(cnorm);
    const scale = try latrs(al, check.io, T, upper, op, unit, false, n, a.ptr, n, x.ptr, cnorm.ptr);
    try testing.expect(scale >= 0 and scale <= 1);
    var xmax: F = 0;
    for (x) |v| xmax = @max(xmax, abs1(T, v));
    try testing.expect(std.math.isFinite(xmax));
    var amax: F = 0;
    var rmax: F = 0;
    for (0..n) |i| {
        var s = ops.mulReal(T, b[i], -scale);
        for (0..n) |j| {
            const r, const c = if (op == .n) .{ i, j } else .{ j, i };
            if (if (upper) r > c else r < c) continue;
            var e = if (r == c) (if (unit) ops.one(T) else ops.scalar(T, d)) else a[r + c * n];
            if (op == .c) e = ops.conj(T, e);
            amax = @max(amax, abs1(T, e));
            s = ops.add(T, s, ops.mul(T, e, x[j]));
        }
        rmax = @max(rmax, abs1(T, s));
    }
    const eps = std.math.floatEps(F);
    const nf: F = @floatFromInt(n);
    const denom = amax * xmax * nf * eps;
    if (!(rmax <= thresh * denom)) {
        check.say("latrs {s} n {d} upper {} op {s} unit {} d {e}: residual {e}, bound {e}\n", .{ @typeName(T), n, upper, @tagName(op), unit, @as(f64, @floatCast(d)), @as(f64, @floatCast(rmax)), @as(f64, @floatCast(denom)) });
        return error.TestUnexpectedResult;
    }
}

test "lapack norm: latrs solves with scaling, ordinary and near-overflow triangles, every form" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const op_list: []const Op = if (comptime isComplex(T)) &.{ .n, .t, .c } else &.{ .n, .t };
        for ([_]bool{ true, false }) |upper| for (op_list) |op| for ([_]bool{ false, true }) |unit| {
            for ([_]usize{ 1, 2, 7, 40 }, 0..) |n, k| {
                // a diagonal of 3 (trsv), and one so small the solution would overflow
                for ([_]F{ 3, std.math.floatMin(F) * 64 }) |d| try latrsCase(T, al, upper, op, unit, n, d, n * 5 + k);
            }
        };
    }
}

test "lapack norm: laset fills the diagonal and the chosen strict triangle, wide and tall" {
    inline for (.{ f64, types.Complex64 }) |T| {
        const x = chol.scalar(T, -1);
        const al_ = chol.scalar(T, 2);
        const be = chol.scalar(T, 5);
        for ([_][2]usize{ .{ 3, 3 }, .{ 2, 4 }, .{ 4, 2 } }) |mn| for ([_]?types.Uplo{ null, .u, .l }) |uplo| {
            const m = mn[0];
            const n = mn[1];
            var a: [16]T = @splat(x);
            try laset(T, uplo, @intCast(m), @intCast(n), al_, be, &a, @intCast(m));
            for (0..n) |j| for (0..m) |i| {
                const inside = uplo == null or (uplo == .u and i <= j) or (uplo == .l and i >= j);
                const want = if (!inside) x else if (i == j) be else al_;
                try testing.expectEqual(want, a[i + j * m]);
            };
        };
    }
    var a: [4]f64 = @splat(0);
    try testing.expectError(error.InvalidM, laset(f64, null, -1, 1, 0, 0, &a, 1));
    try testing.expectError(error.InvalidN, laset(f64, null, 1, -1, 0, 0, &a, 1));
    try testing.expectError(error.InvalidLda, laset(f64, null, 2, 2, 0, 0, &a, 1));
}
