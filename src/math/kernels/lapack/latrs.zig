//! `xLATRS` and `xRSCL`: a triangular solve scaled against overflow, Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dlatrs.f`, `zlatrs.f`, `drscl.f`, OpenBLAS v0.3.34-317; modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const strided = @import("../blas/strided.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const Error = ops.Error;
const Strided = strided.Strided;

/// `?RSCL`, `?DRSCL` (`drscl.f`, `zdrscl.f`): `x := x / sa` by safe
/// multiplications, `sa` real, never forming an overflowing `1 / sa`.
pub fn rscl(io: Io, comptime T: type, sa: Float(T), x: Strided(T)) Error!void {
    const F = Float(T);
    const n = x.n;
    if (n == 0) return;
    // elementwise: memory order serves either stride sign
    const inc: usize = @abs(x.stride);
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
        try blas_calls.scal(io, T, n, mul, x.data.ptr, inc);
        if (done) return;
    }
}

/// `?LATRS` (`dlatrs.f`, `zlatrs.f`): solves `op(A) x = scale b` (`op` `n`,
/// `t`, or `c` for complex) for the `n x n` triangular `A` with `scale <= 1`
/// chosen so `x` cannot overflow; returns `scale`. `cnorm` (`n`) holds the
/// off-diagonal column norms (`abs1` sums), computed unless `normin`. A zero
/// on the diagonal gives a null vector with `scale = 0`. `trsv` takes over
/// where the bounds allow. `x` and `cnorm` take `n` each.
pub fn latrs(al: Allocator, io: Io, comptime T: type, upper: bool, op: Trans, unit: bool, normin: bool, n: usize, a: []const T, lda: usize, x: []T, cnorm: []Float(T)) LatrsError!Float(T) {
    switch (op) {
        .n, .t, .c => {},
        .r, .o, .u, .s, .d => return error.UnsupportedTranspose,
    }
    try util.holds(a.len, n, n, lda);
    if (x.len < n or cnorm.len < n) return error.SliceTooShort;
    return scaledSolve(al, io, T, upper, op, unit, normin, n, a.ptr, lda, x.ptr, cnorm.ptr);
}

/// What `latrs` returns: an `op` other than `n`, `t` or `c`, a slice too
/// short or a leading dimension too small, and the allocator's and
/// cancellation.
pub const LatrsError = interface.TransError || interface.ShapeError || Error;

/// `latrs` after its checks.
pub fn scaledSolve(al: Allocator, io: Io, comptime T: type, upper: bool, op: Trans, unit: bool, normin: bool, n: usize, a: [*]const T, lda: usize, x: [*]T, cnorm: [*]Float(T)) Error!Float(T) {
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
    const am: ops.Mat(T) = .{ .ptr = @constCast(a), .rows = n, .cols = n, .ld = lda };
    if (!normin) {
        for (0..n) |j| {
            const lo, const len = if (upper) .{ @as(usize, 0), j } else .{ j + 1, n - j - 1 };
            var s: F = 0;
            for (0..len) |i| s += ops.abs1(T, am.get(lo + i, j));
            cnorm[j] = s;
        }
    }
    // scale the column norms by `tscal` if the largest exceeds `bignum`
    var tmax = cnorm[blas_calls.iamax(F, n, cnorm, 1)];
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
                    const v = am.get(lo + i, j);
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
                    for (0..len) |i| cnorm[j] += tscal * cabs2(T, am.get(lo + i, j));
                    if (cx) tscal *= 0.5;
                }
            } else {
                try blas_calls.trsv(al, io, T, upper, op, unit, n, a, lda, x);
                return scale;
            }
        }
    }
    // bound the growth of the solution
    var xmax: F = if (cx) blk: {
        var m: F = 0;
        for (0..n) |j| m = @max(m, cabs2(T, x[j]));
        break :blk m;
    } else @abs(x[blas_calls.iamax(T, n, x, 1)]);
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
                    const tjj = ops.abs1(T, am.get(j, j));
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
                    const tjj = ops.abs1(T, am.get(j, j));
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
        try blas_calls.trsv(al, io, T, upper, op, unit, n, a, lda, x);
    } else {
        if (xmax > bignum * half) {
            scale = (bignum * half) / xmax;
            try blas_calls.scal(io, T, n, scale, x, 1);
            xmax = bignum;
        } else if (cx) {
            xmax *= 2;
        }
        const conj = cx and op == .c;
        for (0..n) |k| {
            const j = step(forward, n, k);
            if (notran) {
                // x(j) / A(j, j), then the update of the rest by column j
                var xj = ops.abs1(T, x[j]);
                skip: {
                    const tjjs = if (!unit) ops.mulReal(T, am.get(j, j), tscal) else blk: {
                        if (tscal == 1) break :skip;
                        break :blk ops.scalar(T, tscal);
                    };
                    try divideScaled(io, T, n, x, j, xj, tjjs, cnorm[j], bignum, smlnum, &scale, &xmax);
                    xj = ops.abs1(T, x[j]);
                }
                if (xj > 1) {
                    var rec = 1 / xj;
                    if (cnorm[j] > (bignum - xmax) * rec) {
                        rec *= 0.5;
                        try blas_calls.scal(io, T, n, rec, x, 1);
                        scale *= rec;
                    }
                } else if (xj * cnorm[j] > bignum - xmax) {
                    try blas_calls.scal(io, T, n, 0.5, x, 1);
                    scale *= 0.5;
                }
                const alpha = ops.neg(T, ops.mulReal(T, x[j], tscal));
                if (upper) {
                    if (j > 0) {
                        try blas_calls.axpy(io, T, j, alpha, a + j * lda, 1, x, 1);
                        xmax = ops.abs1(T, x[blas_calls.iamax(T, j, x, 1)]);
                    }
                } else if (j + 1 < n) {
                    try blas_calls.axpy(io, T, n - j - 1, alpha, a + j + 1 + j * lda, 1, x + j + 1, 1);
                    xmax = ops.abs1(T, x[j + 1 + blas_calls.iamax(T, n - j - 1, x + j + 1, 1)]);
                }
            } else {
                // x(j) - op(A)(j, :) x, then the division
                var xj = ops.abs1(T, x[j]);
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
                    const tjj = ops.abs1(T, tjjs);
                    if (tjj > 1) {
                        rec = @min(1, rec * tjj);
                        uscal = ops.div(T, uscal, tjjs);
                    }
                    if (rec < 1) {
                        try blas_calls.scal(io, T, n, rec, x, 1);
                        scale *= rec;
                        xmax *= rec;
                    }
                }
                const lo, const len = if (upper) .{ @as(usize, 0), j } else .{ j + 1, n - j - 1 };
                var sumj = ops.zero(T);
                if (eq(T, uscal, ops.one(T))) {
                    sumj = try blas_calls.dot(io, T, conj, len, a + lo + j * lda, 1, x + lo, 1);
                } else {
                    for (0..len) |i| {
                        const aij = am.get(lo + i, j);
                        sumj = ops.add(T, sumj, ops.mul(T, ops.mul(T, if (conj) ops.conj(T, aij) else aij, uscal), x[lo + i]));
                    }
                }
                if (eq(T, uscal, ops.scalar(T, tscal))) {
                    x[j] = ops.sub(T, x[j], sumj);
                    xj = ops.abs1(T, x[j]);
                    skip: {
                        if (unit and tscal == 1) break :skip;
                        const tjjs = diag(a, lda, j, unit, conj, tscal);
                        try divideScaled(io, T, n, x, j, xj, tjjs, null, bignum, smlnum, &scale, &xmax);
                    }
                } else {
                    x[j] = ops.sub(T, ops.div(T, x[j], diag(a, lda, j, unit, conj, tscal)), sumj);
                }
                xmax = @max(xmax, ops.abs1(T, x[j]));
            }
        }
        scale /= tscal;
    }
    if (tscal != 1) for (0..n) |j| {
        cnorm[j] *= 1 / tscal;
    };
    return scale;
}

/// `x(j) / tjjs`, `xj` the `abs1` of `x(j)`, with `x` (`n`) scaled first so
/// the quotient cannot overflow, as `xlatrs.f`'s two copies: by `1 / xj` for
/// `|tjjs|` below 1, by `tjj bignum / xj` (and `1 / cnorm` past 1, when
/// given) for a tiny one; a zero `tjjs` gives `x = e(j)`, `scale = 0`.
fn divideScaled(io: Io, comptime T: type, n: usize, x: [*]T, j: usize, xj: Float(T), tjjs: T, cnorm: ?Float(T), bignum: Float(T), smlnum: Float(T), scale: *Float(T), xmax: *Float(T)) Error!void {
    const tjj = ops.abs1(T, tjjs);
    if (tjj > smlnum) {
        if (tjj < 1 and xj > tjj * bignum) {
            const rec = 1 / xj;
            try blas_calls.scal(io, T, n, rec, x, 1);
            scale.* *= rec;
            xmax.* *= rec;
        }
        x[j] = ops.div(T, x[j], tjjs);
    } else if (tjj > 0) {
        if (xj > tjj * bignum) {
            var rec = (tjj * bignum) / xj;
            if (cnorm) |c| {
                if (c > 1) rec /= c;
            }
            try blas_calls.scal(io, T, n, rec, x, 1);
            scale.* *= rec;
            xmax.* *= rec;
        }
        x[j] = ops.div(T, x[j], tjjs);
    } else {
        @memset(x[0..n], ops.zero(T));
        x[j] = ops.one(T);
        scale.* = 0;
        xmax.* = 0;
    }
}

/// `CABS2` (`zlatrs.f`): `|re / 2| + |im / 2|`; `|x|` for real.
fn cabs2(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) @abs(x.re / 2) + @abs(x.im / 2) else @abs(x);
}

fn eq(comptime T: type, x: T, y: T) bool {
    return if (comptime isComplex(T)) x.re == y.re and x.im == y.im else x == y;
}
