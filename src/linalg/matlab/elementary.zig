//! Elementary complex scalar functions ported from `matlab88.f90`'s
//! `MATFN2` (the `SIN`/`COS`/`EXP`/`POWER` command bodies) and its `ROUND`
//! utility. `blas1.zig` already carries `wsqrt`/`wlog`/`watan` (needed
//! internally by the eigen/SVD/QR kernels); this file rounds out the rest
//! of the elementary functions the original exposed as REPL commands.
//!
//! These operate on a single `Cx` scalar. For applying one across a whole
//! vector or matrix, see `api.zig`'s elementwise (`sin`, `cos`, ...) and
//! matrix-function (`sinm`, `cosm`, ...) wrappers.

const std = @import("std");
const blas1 = @import("blas1.zig");
const Cx = blas1.Cx;

/// SIN: complex sine, `sin(x+iy) = sin(x)*cosh(y) + i*cos(x)*sinh(y)`.
pub fn wsin(x: Cx) Cx {
    return Cx.init(@sin(x.re) * std.math.cosh(x.im), @cos(x.re) * std.math.sinh(x.im));
}

/// COS: complex cosine, `cos(x+iy) = cos(x)*cosh(y) - i*sin(x)*sinh(y)`.
pub fn wcos(x: Cx) Cx {
    return Cx.init(@cos(x.re) * std.math.cosh(x.im), -(@sin(x.re) * std.math.sinh(x.im)));
}

/// EXP: complex exponential.
pub fn wexp(x: Cx) Cx {
    const r = @exp(x.re);
    return Cx.init(r * @cos(x.im), r * @sin(x.im));
}

/// POWER: `x^y = exp(y * log(x))` (principal branch).
pub fn wpow(x: Cx, y: Cx) Cx {
    return wexp(blas1.wmul(blas1.wlog(x), y));
}

/// ROUND: round `x` to the nearest integer via exact binary decomposition
/// of its magnitude rather than naive `x + 0.5` truncation, so rounding
/// stays exact even when `x` is too large for `+0.5` to change its nearest
/// representable value (returns `x` unchanged in that case).
pub fn round(x: f64) f64 {
    const h: f64 = 1.0e9;
    var z = @abs(x);
    var y = z + 1.0;
    if (y == z) return x;

    y = 0.0;
    var e: f64 = h;
    while (e < z) e *= 2.0;
    while (e > h) {
        if (e <= z) {
            y += e;
            z -= e;
        }
        e /= 2.0;
    }
    z = @floor(z + 0.5);
    y += z;
    if (x < 0.0) y = -y;
    return y;
}

const testing = std.testing;

test "wsin/wcos/wexp agree with the real axis and Euler's identity" {
    const x = Cx.init(1.3, 0.0);
    try testing.expectApproxEqAbs(@sin(1.3), wsin(x).re, 1e-12);
    try testing.expectApproxEqAbs(@cos(1.3), wcos(x).re, 1e-12);
    try testing.expectApproxEqAbs(@exp(1.3), wexp(x).re, 1e-12);

    // e^(i*pi) = -1.
    const ipi = Cx.init(0.0, std.math.pi);
    const r = wexp(ipi);
    try testing.expectApproxEqAbs(-1.0, r.re, 1e-9);
    try testing.expectApproxEqAbs(0.0, r.im, 1e-9);
}

test "wpow matches repeated multiplication for an integer exponent" {
    const x = Cx.init(1.5, -0.5);
    const cubed = blas1.wmul(blas1.wmul(x, x), x);
    const viaPow = wpow(x, Cx.init(3.0, 0.0));
    try testing.expectApproxEqAbs(cubed.re, viaPow.re, 1e-9);
    try testing.expectApproxEqAbs(cubed.im, viaPow.im, 1e-9);
}

test "round matches ordinary rounding for everyday magnitudes" {
    try testing.expectEqual(@as(f64, 3.0), round(2.6));
    try testing.expectEqual(@as(f64, -3.0), round(-2.6));
    try testing.expectEqual(@as(f64, 2.0), round(2.4));
    try testing.expectEqual(@as(f64, 0.0), round(0.0));
}
