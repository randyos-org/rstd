//! Complex/real BLAS level-1 style kernels ported from `matlab88.f90` (the
//! 1988 Argonne National Laboratory MATLAB source; see `matlab88_doc.txt` /
//! `matlab88_help.txt` in this directory). Names match the original Fortran
//! (case-insensitive there), lowercased to Zig convention.
//!
//! Complex vectors/matrices keep the original split real/imaginary array
//! layout (`xr`/`xi` rather than an interleaved `Complex(f64)`) on purpose:
//! it matches `AR`/`AI` in the Fortran, and it is what lets the stride-1
//! fast paths below vectorize with `@Vector` instead of needing shuffles to
//! pull real/imaginary lanes apart. A single complex scalar (a multiplier,
//! a temporary, a return value) is still `Cx = std.math.Complex(f64)`.
//!
//! An explicit `incx`/`incy` stride (may be negative, as in reference BLAS)
//! mirrors the Fortran calling convention: the slice already starts at the
//! first accessed element and the stride selects successive ones - passing
//! `AR(K,J)` as an actual argument in Fortran hands the callee a new base
//! address into the same array. Non-unit strides fall back to a plain
//! scalar loop; only `incx == incy == 1` (the common case in this port)
//! takes the SIMD path.
//!
//! The original's `FLOP(...)` wrapper (a flop-counter for the REPL's status
//! line) is dropped throughout this port; it had no effect on results.

const std = @import("std");

const rsimd = @import("../../simd/__root__.zig");
const mulAdd = rsimd.mulAdd;

pub const Cx = std.math.Complex(f64);

const vec_len = rsimd.zigvec.f64_simd_lanes;
const Vf64 = rsimd.zigvec.SimdVecF64;

inline fn loadVec(s: []const f64, i: usize) Vf64 {
    return s[i..][0..vec_len].*;
}

inline fn storeVec(s: []f64, i: usize, v: Vf64) void {
    s[i..][0..vec_len].* = v;
}

/// CABS1: `|re| + |im|`, the cheap "1-norm" pivoting/singularity metric
/// LINPACK uses in place of the true complex magnitude. `anytype` (rather
/// than a fixed `f64`) so the same body also works as `reduce.mapReduce`'s
/// per-chunk map in `wasum` below, called once against `f64` scalars (the
/// tail) and once against `Vf64` chunks (the SIMD body).
pub fn cabs1(re: anytype, im: @TypeOf(re)) @TypeOf(re) {
    return @abs(re) + @abs(im);
}

/// 0-based offset of the first element accessed by a Fortran-style
/// 1-based loop of `n` steps with increment `inc` (negative increments walk
/// backward from the end, as in the reference BLAS).
fn stride0(n: usize, inc: isize) isize {
    return if (inc < 0) (-@as(isize, @intCast(n)) + 1) * inc else 0;
}

/// 0-based flat offset of 1-based matrix element `(i, j)` in a column-major
/// array with leading dimension `lda` - the Zig equivalent of Fortran
/// handing a subroutine the address of `A(i, j)`.
pub fn ix(lda: usize, i: usize, j: usize) usize {
    return (i - 1) + (j - 1) * lda;
}

/// PYTHAG: `sqrt(a*a + b*b)` without spurious overflow/underflow.
pub fn pythag(a: f64, b: f64) f64 {
    return std.math.hypot(a, b);
}

/// WMUL: `c = a*b`.
pub fn wmul(a: Cx, b: Cx) Cx {
    return a.mul(b);
}

/// WDIV: `c = a/b`, scaled (Smith's algorithm) to avoid overflow.
pub fn wdiv(a: Cx, b: Cx) Cx {
    const s = @abs(b.re) + @abs(b.im);
    const ars = a.re / s;
    const ais = a.im / s;
    const brs = b.re / s;
    const bis = b.im / s;
    const d = mulAdd(f64, bis, bis, brs * brs);
    return Cx.init(
        mulAdd(f64, ais, bis, ars * brs) / d,
        mulAdd(f64, -ars, bis, ais * brs) / d,
    );
}

/// WSIGN: if `y != 0`, `z = x*y/|y|`; if `y == 0`, `z = x`.
pub fn wsign(x: Cx, y: Cx) Cx {
    const t = pythag(y.re, y.im);
    if (t == 0.0) return x;
    return wmul(Cx.init(y.re / t, y.im / t), x);
}

/// WLOG: `y = log(x)`. Delegates to `std.math.complex.log`, except this
/// port's branch-cut convention for the negative real axis is preserved
/// explicitly: for `x.im == 0.0` (which also matches `-0.0`) and
/// `x.re < 0.0`, the angle is forced to `+pi` rather than left to
/// `atan2`'s signed-zero handling (which would give `-pi` for a `-0.0`
/// imaginary part).
pub fn wlog(x: Cx) Cx {
    var r = std.math.complex.log(x);
    if (x.im == 0.0 and x.re < 0.0) r.im = @abs(r.im);
    return r;
}

/// WATAN: `y = atan(x) = (i/2)*log((i+x)/(i-x))`.
pub fn watan(x: Cx) Cx {
    if (x.im == 0.0) return Cx.init(std.math.atan2(x.re, 1.0), 0.0);
    const t = wdiv(Cx.init(x.re, 1.0 + x.im), Cx.init(-x.re, 1.0 - x.im));
    const l = wlog(t);
    return Cx.init(-l.im / 2.0, l.re / 2.0);
}

/// WNRM2: Euclidean norm of `x`. Tries a fast SIMD sum-of-squares first
/// (unit stride only) and only falls back to `wnrm2Safe`'s sequential
/// `pythag`-chain (the original's own algorithm) when that visibly
/// failed - matches `linalg.vector.length`'s gated design, and for the
/// same reason: chaining `pythag` unconditionally costs real per-call
/// overhead (NaN/Inf checks, a compensated-arithmetic correction term)
/// that a lane-parallel sum-of-squares avoids for the overwhelmingly
/// common case of ordinary-magnitude components.
pub fn wnrm2(n: usize, xr: []const f64, xi: []const f64, incx: isize) f64 {
    if (n == 0) return 0.0;
    if (incx == 1) {
        const sumSquares = struct {
            fn f(re: anytype, im: @TypeOf(re)) @TypeOf(re) {
                return mulAdd(@TypeOf(re), im, im, re * re);
            }
        }.f;
        const sq = rsimd.reduce.mapReduce(.Add, f64, sumSquares, .{ xr[0..n], xi[0..n] });
        // sq == 0.0 is ambiguous by itself: genuinely all-zero input, or
        // every squared term underflowed - wasum (sum of non-negative
        // |re|+|im| terms, no cancellation possible) tells them apart
        // cheaply, and already has its own SIMD fast path too.
        const looks_fine = std.math.isFinite(sq) and (sq != 0.0 or wasum(n, xr, xi, incx) == 0.0);
        if (looks_fine) return @sqrt(sq);
    }
    return wnrm2Safe(n, xr, xi, incx);
}

/// The overflow/underflow-safe fallback `wnrm2` uses when its fast path
/// fails, and directly for non-unit strides (`mapReduce` has no stride
/// support - see `wasum`'s/`wdot`'s doc comments for the same reasoning).
fn wnrm2Safe(n: usize, xr: []const f64, xi: []const f64, incx: isize) f64 {
    var s: f64 = 0.0;
    if (n == 0) return s;
    var i_x = stride0(n, incx);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const u: usize = @intCast(i_x);
        s = pythag(s, xr[u]);
        s = pythag(s, xi[u]);
        i_x += incx;
    }
    return s;
}

/// WASUM: 1-norm of `x` (sum of `|re|+|im|`).
pub fn wasum(n: usize, xr: []const f64, xi: []const f64, incx: isize) f64 {
    if (n == 0) return 0.0;
    if (incx == 1) return rsimd.reduce.mapReduce(.Add, f64, cabs1, .{ xr[0..n], xi[0..n] });
    var s: f64 = 0.0;
    var i_x = stride0(n, incx);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const u: usize = @intCast(i_x);
        s += cabs1(xr[u], xi[u]);
        i_x += incx;
    }
    return s;
}

/// IWAMAX: 1-based index of the element of largest `|re|+|im|` (0 if `n==0`).
/// Kept 1-based (unlike the rest of this port) because call sites lean on
/// Fortran's `iwamax(...) + k - 1` idiom to turn a local index into a global
/// one; renumbering here would just move the `+1` elsewhere.
pub fn iwamax(n: usize, xr: []const f64, xi: []const f64, incx: isize) usize {
    if (n == 0) return 0;
    var k: usize = 1;
    var s: f64 = 0.0;
    var i_x = stride0(n, incx);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const u: usize = @intCast(i_x);
        const p = cabs1(xr[u], xi[u]);
        if (p > s) {
            k = i + 1;
            s = p;
        }
        i_x += incx;
    }
    return k;
}

/// WRSCAL: `x *= s` for a real scalar `s`.
pub fn wrscal(n: usize, s: f64, xr: []f64, xi: []f64, incx: isize) void {
    if (n == 0) return;
    if (incx == 1) {
        const sv: Vf64 = @splat(s);
        var i: usize = 0;
        while (i + vec_len <= n) : (i += vec_len) {
            storeVec(xr, i, loadVec(xr, i) * sv);
            storeVec(xi, i, loadVec(xi, i) * sv);
        }
        while (i < n) : (i += 1) {
            xr[i] *= s;
            xi[i] *= s;
        }
        return;
    }
    var i_x = stride0(n, incx);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const u: usize = @intCast(i_x);
        xr[u] *= s;
        xi[u] *= s;
        i_x += incx;
    }
}

/// WSCAL: `x *= s` for a complex scalar `s`.
pub fn wscal(n: usize, s: Cx, xr: []f64, xi: []f64, incx: isize) void {
    if (n == 0) return;
    if (incx == 1) {
        const sr: Vf64 = @splat(s.re);
        const si: Vf64 = @splat(s.im);
        var i: usize = 0;
        while (i + vec_len <= n) : (i += vec_len) {
            const xrv = loadVec(xr, i);
            const xiv = loadVec(xi, i);
            storeVec(xr, i, mulAdd(Vf64, -si, xiv, sr * xrv));
            storeVec(xi, i, mulAdd(Vf64, si, xrv, sr * xiv));
        }
        while (i < n) : (i += 1) {
            const t = wmul(s, Cx.init(xr[i], xi[i]));
            xr[i] = t.re;
            xi[i] = t.im;
        }
        return;
    }
    var i_x = stride0(n, incx);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const u: usize = @intCast(i_x);
        const t = wmul(s, Cx.init(xr[u], xi[u]));
        xr[u] = t.re;
        xi[u] = t.im;
        i_x += incx;
    }
}

/// WAXPY: `y += s*x`.
pub fn waxpy(
    n: usize,
    s: Cx,
    xr: []const f64,
    xi: []const f64,
    incx: isize,
    yr: []f64,
    yi: []f64,
    incy: isize,
) void {
    if (n == 0 or (s.re == 0.0 and s.im == 0.0)) return;
    if (incx == 1 and incy == 1) {
        const sr: Vf64 = @splat(s.re);
        const si: Vf64 = @splat(s.im);
        var i: usize = 0;
        while (i + vec_len <= n) : (i += vec_len) {
            const xrv = loadVec(xr, i);
            const xiv = loadVec(xi, i);
            storeVec(yr, i, mulAdd(Vf64, -si, xiv, mulAdd(Vf64, sr, xrv, loadVec(yr, i))));
            storeVec(yi, i, mulAdd(Vf64, si, xrv, mulAdd(Vf64, sr, xiv, loadVec(yi, i))));
        }
        while (i < n) : (i += 1) {
            yr[i] = mulAdd(f64, -s.im, xi[i], mulAdd(f64, s.re, xr[i], yr[i]));
            yi[i] = mulAdd(f64, s.im, xr[i], mulAdd(f64, s.re, xi[i], yi[i]));
        }
        return;
    }
    var i_x = stride0(n, incx);
    var i_y = stride0(n, incy);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const xu: usize = @intCast(i_x);
        const yu: usize = @intCast(i_y);
        yr[yu] = mulAdd(f64, -s.im, xi[xu], mulAdd(f64, s.re, xr[xu], yr[yu]));
        yi[yu] = mulAdd(f64, s.im, xr[xu], mulAdd(f64, s.re, xi[xu], yi[yu]));
        i_x += incx;
        i_y += incy;
    }
}

const DotKind = enum { unconjugated_re, unconjugated_im, conjugated_re, conjugated_im };

fn wdot(
    comptime kind: DotKind,
    n: usize,
    xr: []const f64,
    xi: []const f64,
    incx: isize,
    yr: []const f64,
    yi: []const f64,
    incy: isize,
) f64 {
    if (n == 0) return 0.0;
    if (incx == 1 and incy == 1) {
        const mapFn = struct {
            fn f(xrv: anytype, xiv: @TypeOf(xrv), yrv: @TypeOf(xrv), yiv: @TypeOf(xrv)) @TypeOf(xrv) {
                const T = @TypeOf(xrv);
                return switch (kind) {
                    .unconjugated_re => mulAdd(T, -xiv, yiv, xrv * yrv),
                    .unconjugated_im => mulAdd(T, xiv, yrv, xrv * yiv),
                    .conjugated_re => mulAdd(T, xiv, yiv, xrv * yrv),
                    .conjugated_im => mulAdd(T, -xiv, yrv, xrv * yiv),
                };
            }
        }.f;
        return rsimd.reduce.mapReduce(.Add, f64, mapFn, .{ xr[0..n], xi[0..n], yr[0..n], yi[0..n] });
    }
    var s: f64 = 0.0;
    var i_x = stride0(n, incx);
    var i_y = stride0(n, incy);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const xu: usize = @intCast(i_x);
        const yu: usize = @intCast(i_y);
        s += switch (kind) {
            .unconjugated_re => mulAdd(f64, -xi[xu], yi[yu], xr[xu] * yr[yu]),
            .unconjugated_im => mulAdd(f64, xi[xu], yr[yu], xr[xu] * yi[yu]),
            .conjugated_re => mulAdd(f64, xi[xu], yi[yu], xr[xu] * yr[yu]),
            .conjugated_im => mulAdd(f64, -xi[xu], yr[yu], xr[xu] * yi[yu]),
        };
        i_x += incx;
        i_y += incy;
    }
    return s;
}

/// WDOTUR: `re(sum(x_k * y_k))` (unconjugated dot product).
pub fn wdotur(n: usize, xr: []const f64, xi: []const f64, incx: isize, yr: []const f64, yi: []const f64, incy: isize) f64 {
    return wdot(.unconjugated_re, n, xr, xi, incx, yr, yi, incy);
}

/// WDOTUI: `im(sum(x_k * y_k))` (unconjugated dot product).
pub fn wdotui(n: usize, xr: []const f64, xi: []const f64, incx: isize, yr: []const f64, yi: []const f64, incy: isize) f64 {
    return wdot(.unconjugated_im, n, xr, xi, incx, yr, yi, incy);
}

/// WDOTCR: `re(sum(conj(x_k) * y_k))` (conjugated dot product).
pub fn wdotcr(n: usize, xr: []const f64, xi: []const f64, incx: isize, yr: []const f64, yi: []const f64, incy: isize) f64 {
    return wdot(.conjugated_re, n, xr, xi, incx, yr, yi, incy);
}

/// WDOTCI: `im(sum(conj(x_k) * y_k))` (conjugated dot product).
pub fn wdotci(n: usize, xr: []const f64, xi: []const f64, incx: isize, yr: []const f64, yi: []const f64, incy: isize) f64 {
    return wdot(.conjugated_im, n, xr, xi, incx, yr, yi, incy);
}

/// WCOPY: `y = x`.
pub fn wcopy(n: usize, xr: []const f64, xi: []const f64, incx: isize, yr: []f64, yi: []f64, incy: isize) void {
    if (n == 0) return;
    if (incx == 1 and incy == 1) {
        @memcpy(yr[0..n], xr[0..n]);
        @memcpy(yi[0..n], xi[0..n]);
        return;
    }
    var i_x = stride0(n, incx);
    var i_y = stride0(n, incy);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        yr[@intCast(i_y)] = xr[@intCast(i_x)];
        yi[@intCast(i_y)] = xi[@intCast(i_x)];
        i_x += incx;
        i_y += incy;
    }
}

/// WSET: fill `y` with the constant `val`.
pub fn wset(n: usize, val: Cx, yr: []f64, yi: []f64, incy: isize) void {
    if (n == 0) return;
    if (incy == 1) {
        @memset(yr[0..n], val.re);
        @memset(yi[0..n], val.im);
        return;
    }
    var i_y = stride0(n, incy);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const u: usize = @intCast(i_y);
        yr[u] = val.re;
        yi[u] = val.im;
        i_y += incy;
    }
}

/// WSWAP: exchange `x` and `y`.
pub fn wswap(n: usize, xr: []f64, xi: []f64, incx: isize, yr: []f64, yi: []f64, incy: isize) void {
    if (n == 0) return;
    var i_x = stride0(n, incx);
    var i_y = stride0(n, incy);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const xu: usize = @intCast(i_x);
        const yu: usize = @intCast(i_y);
        std.mem.swap(f64, &xr[xu], &yr[yu]);
        std.mem.swap(f64, &xi[xu], &yi[yu]);
        i_x += incx;
        i_y += incy;
    }
}

/// RSET: fill `y` with the real constant `val`.
pub fn rset(n: usize, val: f64, y: []f64, incy: isize) void {
    if (n == 0) return;
    if (incy == 1) {
        @memset(y[0..n], val);
        return;
    }
    var i_y = stride0(n, incy);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        y[@intCast(i_y)] = val;
        i_y += incy;
    }
}

/// RSWAP: exchange real vectors `x` and `y`.
pub fn rswap(n: usize, x: []f64, incx: isize, y: []f64, incy: isize) void {
    if (n == 0) return;
    var i_x = stride0(n, incx);
    var i_y = stride0(n, incy);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        std.mem.swap(f64, &x[@intCast(i_x)], &y[@intCast(i_y)]);
        i_x += incx;
        i_y += incy;
    }
}

/// RROT: apply the plane rotation `(c, s)` to real vectors `x` and `y`.
pub fn rrot(n: usize, x: []f64, incx: isize, y: []f64, incy: isize, c: f64, s: f64) void {
    if (n == 0) return;
    if (incx == 1 and incy == 1) {
        const cv: Vf64 = @splat(c);
        const sv: Vf64 = @splat(s);
        var i: usize = 0;
        while (i + vec_len <= n) : (i += vec_len) {
            const xv = loadVec(x, i);
            const yv = loadVec(y, i);
            storeVec(x, i, mulAdd(Vf64, sv, yv, cv * xv));
            storeVec(y, i, mulAdd(Vf64, -sv, xv, cv * yv));
        }
        while (i < n) : (i += 1) {
            const temp = mulAdd(f64, s, y[i], c * x[i]);
            y[i] = mulAdd(f64, -s, x[i], c * y[i]);
            x[i] = temp;
        }
        return;
    }
    var i_x = stride0(n, incx);
    var i_y = stride0(n, incy);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const xu: usize = @intCast(i_x);
        const yu: usize = @intCast(i_y);
        const temp = mulAdd(f64, s, y[yu], c * x[xu]);
        y[yu] = mulAdd(f64, -s, x[xu], c * y[yu]);
        x[xu] = temp;
        i_x += incx;
        i_y += incy;
    }
}

pub const Givens = struct { r: f64, c: f64, s: f64 };

/// RROTG: construct a Givens plane rotation `(c, s)` such that
/// `c*da + s*db == r` and `-s*da + c*db == 0`.
///
/// The original also returns a compact `z` encoding from which `(c, s)`
/// could be reconstructed later; every call site in this port immediately
/// overwrites its own copy of `db` without ever reading `z` back, so it is
/// omitted here.
pub fn rrotg(da: f64, db: f64) Givens {
    var rho = db;
    if (@abs(da) > @abs(db)) rho = da;
    var c: f64 = 1.0;
    var s: f64 = 0.0;
    const r = std.math.copysign(pythag(da, db), rho);
    if (r != 0.0) {
        c = da / r;
        s = db / r;
    }
    return .{ .r = r, .c = c, .s = s };
}

const zigstd = @import("../../zigstd.zig");

test "wdiv is the inverse of wmul" {
    const a = Cx.init(3.0, -2.0);
    const b = Cx.init(-1.5, 4.0);
    const prod = wmul(a, b);
    const back = wdiv(prod, b);
    try zigstd.testExpectApproxEqAbs(a.re, back.re, 1e-12);
    try zigstd.testExpectApproxEqAbs(a.im, back.im, 1e-12);
}

test "waxpy/wdot round trip over a SIMD-width-crossing vector" {
    const n = vec_len * 2 + 3;
    var xr: [n]f64 = undefined;
    var xi: [n]f64 = undefined;
    var yr: [n]f64 = undefined;
    var yi: [n]f64 = undefined;
    for (0..n) |k| {
        xr[k] = @floatFromInt(k + 1);
        xi[k] = @floatFromInt(2 * k);
        yr[k] = 1.0;
        yi[k] = 0.0;
    }
    waxpy(n, Cx.init(2.0, 0.0), &xr, &xi, 1, &yr, &yi, 1);
    for (0..n) |k| {
        try zigstd.testExpectApproxEqAbs(1.0 + 2.0 * xr[k], yr[k], 1e-12);
        try zigstd.testExpectApproxEqAbs(2.0 * xi[k], yi[k], 1e-12);
    }

    const dot_re = wdotcr(n, &xr, &xi, 1, &xr, &xi, 1);
    const expected: f64 = blk: {
        var acc: f64 = 0.0;
        for (0..n) |k| acc += xr[k] * xr[k] + xi[k] * xi[k];
        break :blk acc;
    };
    try zigstd.testExpectApproxEqAbs(expected, dot_re, 1e-9);
}

test "iwamax finds the largest-magnitude element" {
    const xr = [_]f64{ 1.0, -2.0, 0.5 };
    const xi = [_]f64{ 0.0, 0.0, 9.0 };
    try zigstd.testExpectEqual(@as(usize, 3), iwamax(3, &xr, &xi, 1));
}

test "wnrm2 matches a hand-computed value for ordinary magnitudes" {
    const xr = [_]f64{ 3.0, 0.0 };
    const xi = [_]f64{ 0.0, 4.0 };
    try zigstd.testExpectApproxEqAbs(5.0, wnrm2(2, &xr, &xi, 1), 1e-12);
}

test "wnrm2's fast path falls back instead of overflowing to inf on huge components" {
    // Squaring 1e300 overflows f64 to +inf; the fallback's incremental
    // pythag chain never squares anything that large without hypot's own
    // overflow-safe rescale first.
    const xr = [_]f64{ 1e300, 1e300 };
    const xi: [2]f64 = @splat(0.0);
    const want = @sqrt(2.0) * 1e300;
    const got = wnrm2(2, &xr, &xi, 1);
    try zigstd.testExpectApproxEqAbs(want, got, want * 1e-12);
    try zigstd.testExpectTrue(std.math.isFinite(got));
}

test "wnrm2's fast path falls back instead of underflowing to 0 on tiny components" {
    // Squaring 1e-200 underflows f64 to exactly 0, losing the vector
    // entirely if trusted directly; wasum's all-zero check catches it.
    const xr = [_]f64{ 1e-200, 1e-200 };
    const xi: [2]f64 = @splat(0.0);
    const want = @sqrt(2.0) * 1e-200;
    const got = wnrm2(2, &xr, &xi, 1);
    try zigstd.testExpectApproxEqAbs(want, got, want * 1e-12);
    try zigstd.testExpectTrue(got > 0);
}

test "wnrm2 returns exactly 0 for the genuine zero vector, not a false-positive fallback" {
    const xr: [4]f64 = @splat(0.0);
    const xi: [4]f64 = @splat(0.0);
    try zigstd.testExpectEqual(@as(f64, 0.0), wnrm2(4, &xr, &xi, 1));
}

test "wnrm2's fast and safe paths agree on a SIMD-width-crossing vector" {
    const n = vec_len * 2 + 3;
    var xr: [n]f64 = undefined;
    var xi: [n]f64 = undefined;
    for (0..n) |k| {
        xr[k] = @floatFromInt(k + 1);
        xi[k] = @floatFromInt(2 * k);
    }
    try zigstd.testExpectApproxEqAbs(wnrm2(n, &xr, &xi, 1), wnrm2Safe(n, &xr, &xi, 1), 1e-9);
}
