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

pub const Cx = std.math.Complex(f64);

const vec_len: usize = std.simd.suggestVectorLength(f64) orelse 4;
const Vf64 = @Vector(vec_len, f64);

inline fn loadVec(s: []const f64, i: usize) Vf64 {
    return s[i..][0..vec_len].*;
}

inline fn storeVec(s: []f64, i: usize, v: Vf64) void {
    s[i..][0..vec_len].* = v;
}

/// CABS1: `|re| + |im|`, the cheap "1-norm" pivoting/singularity metric
/// LINPACK uses in place of the true complex magnitude.
pub fn cabs1(re: f64, im: f64) f64 {
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
    const d = brs * brs + bis * bis;
    return Cx.init((ars * brs + ais * bis) / d, (ais * brs - ars * bis) / d);
}

/// WSIGN: if `y != 0`, `z = x*y/|y|`; if `y == 0`, `z = x`.
pub fn wsign(x: Cx, y: Cx) Cx {
    const t = pythag(y.re, y.im);
    if (t == 0.0) return x;
    return wmul(Cx.init(y.re / t, y.im / t), x);
}

/// WSQRT: `y = sqrt(x)` with `re(y) >= 0` and `sign(im(y)) == sign(im(x))`.
pub fn wsqrt(x: Cx) Cx {
    const tr = x.re;
    const ti = x.im;
    var s = @sqrt(0.5 * (pythag(tr, ti) + @abs(tr)));
    var yr: f64 = 0.0;
    var yi: f64 = 0.0;
    if (tr >= 0.0) yr = s;
    if (ti < 0.0) s = -s;
    if (tr <= 0.0) yi = s;
    if (tr < 0.0) yr = 0.5 * (ti / yi);
    if (tr > 0.0) yi = 0.5 * (ti / yr);
    return Cx.init(yr, yi);
}

/// WLOG: `y = log(x)`.
pub fn wlog(x: Cx) Cx {
    const r = pythag(x.re, x.im);
    var t = std.math.atan2(x.im, x.re);
    if (x.im == 0.0 and x.re < 0.0) t = @abs(t);
    return Cx.init(@log(r), t);
}

/// WATAN: `y = atan(x) = (i/2)*log((i+x)/(i-x))`.
pub fn watan(x: Cx) Cx {
    if (x.im == 0.0) return Cx.init(std.math.atan2(x.re, 1.0), 0.0);
    const t = wdiv(Cx.init(x.re, 1.0 + x.im), Cx.init(-x.re, 1.0 - x.im));
    const l = wlog(t);
    return Cx.init(-l.im / 2.0, l.re / 2.0);
}

/// WNRM2: Euclidean norm of `x`. Left scalar-only: the original accumulates
/// via repeated `PYTHAG(s, ...)` specifically to stay overflow-safe, which
/// is an inherently sequential refinement rather than a lane-parallel sum.
pub fn wnrm2(n: usize, xr: []const f64, xi: []const f64, incx: isize) f64 {
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
    if (incx == 1) {
        var acc: Vf64 = @splat(0.0);
        var i: usize = 0;
        while (i + vec_len <= n) : (i += vec_len) {
            acc += @abs(loadVec(xr, i)) + @abs(loadVec(xi, i));
        }
        var s: f64 = @reduce(.Add, acc);
        while (i < n) : (i += 1) s += cabs1(xr[i], xi[i]);
        return s;
    }
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
            storeVec(xr, i, sr * xrv - si * xiv);
            storeVec(xi, i, sr * xiv + si * xrv);
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
            storeVec(yr, i, loadVec(yr, i) + sr * xrv - si * xiv);
            storeVec(yi, i, loadVec(yi, i) + sr * xiv + si * xrv);
        }
        while (i < n) : (i += 1) {
            yr[i] += s.re * xr[i] - s.im * xi[i];
            yi[i] += s.re * xi[i] + s.im * xr[i];
        }
        return;
    }
    var i_x = stride0(n, incx);
    var i_y = stride0(n, incy);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const xu: usize = @intCast(i_x);
        const yu: usize = @intCast(i_y);
        yr[yu] += s.re * xr[xu] - s.im * xi[xu];
        yi[yu] += s.re * xi[xu] + s.im * xr[xu];
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
        var acc: Vf64 = @splat(0.0);
        var i: usize = 0;
        while (i + vec_len <= n) : (i += vec_len) {
            const xrv = loadVec(xr, i);
            const xiv = loadVec(xi, i);
            const yrv = loadVec(yr, i);
            const yiv = loadVec(yi, i);
            acc += switch (kind) {
                .unconjugated_re => xrv * yrv - xiv * yiv,
                .unconjugated_im => xrv * yiv + xiv * yrv,
                .conjugated_re => xrv * yrv + xiv * yiv,
                .conjugated_im => xrv * yiv - xiv * yrv,
            };
        }
        var s: f64 = @reduce(.Add, acc);
        while (i < n) : (i += 1) {
            s += switch (kind) {
                .unconjugated_re => xr[i] * yr[i] - xi[i] * yi[i],
                .unconjugated_im => xr[i] * yi[i] + xi[i] * yr[i],
                .conjugated_re => xr[i] * yr[i] + xi[i] * yi[i],
                .conjugated_im => xr[i] * yi[i] - xi[i] * yr[i],
            };
        }
        return s;
    }
    var s: f64 = 0.0;
    var i_x = stride0(n, incx);
    var i_y = stride0(n, incy);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const xu: usize = @intCast(i_x);
        const yu: usize = @intCast(i_y);
        s += switch (kind) {
            .unconjugated_re => xr[xu] * yr[yu] - xi[xu] * yi[yu],
            .unconjugated_im => xr[xu] * yi[yu] + xi[xu] * yr[yu],
            .conjugated_re => xr[xu] * yr[yu] + xi[xu] * yi[yu],
            .conjugated_im => xr[xu] * yi[yu] - xi[xu] * yr[yu],
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
            storeVec(x, i, cv * xv + sv * yv);
            storeVec(y, i, cv * yv - sv * xv);
        }
        while (i < n) : (i += 1) {
            const temp = c * x[i] + s * y[i];
            y[i] = c * y[i] - s * x[i];
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
        const temp = c * x[xu] + s * y[yu];
        y[yu] = c * y[yu] - s * x[xu];
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

const testing = std.testing;

test "wdiv is the inverse of wmul" {
    const a = Cx.init(3.0, -2.0);
    const b = Cx.init(-1.5, 4.0);
    const prod = wmul(a, b);
    const back = wdiv(prod, b);
    try testing.expectApproxEqAbs(a.re, back.re, 1e-12);
    try testing.expectApproxEqAbs(a.im, back.im, 1e-12);
}

test "wsqrt squares back to the input" {
    const x = Cx.init(-3.0, 4.0);
    const r = wsqrt(x);
    const sq = wmul(r, r);
    try testing.expectApproxEqAbs(x.re, sq.re, 1e-12);
    try testing.expectApproxEqAbs(x.im, sq.im, 1e-12);
    try testing.expect(r.re >= 0.0);
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
        try testing.expectApproxEqAbs(1.0 + 2.0 * xr[k], yr[k], 1e-12);
        try testing.expectApproxEqAbs(2.0 * xi[k], yi[k], 1e-12);
    }

    const dot_re = wdotcr(n, &xr, &xi, 1, &xr, &xi, 1);
    const expected: f64 = blk: {
        var acc: f64 = 0.0;
        for (0..n) |k| acc += xr[k] * xr[k] + xi[k] * xi[k];
        break :blk acc;
    };
    try testing.expectApproxEqAbs(expected, dot_re, 1e-9);
}

test "iwamax finds the largest-magnitude element" {
    const xr = [_]f64{ 1.0, -2.0, 0.5 };
    const xi = [_]f64{ 0.0, 0.0, 9.0 };
    try testing.expectEqual(@as(usize, 3), iwamax(3, &xr, &xi, 1));
}
