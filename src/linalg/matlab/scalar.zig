//! `matlab88.f90`'s complex scalars (WDIV, WSIGN, WLOG, WATAN), `ix` and
//! `level1`; level 1 BLAS is rstd's since M2, and FLOP counting is gone.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const ops = @import("../lapack/ops.zig");

const mulAdd = @import("../vector.zig").mulAdd;

pub const Cx = std.math.Complex(f64);

/// 0-based flat offset of 1-based matrix element `(i, j)` in a column-major
/// array with leading dimension `lda`; the Zig equivalent of Fortran
/// handing a subroutine the address of `A(i, j)`.
pub fn ix(lda: usize, i: usize, j: usize) usize {
    return (i - 1) + (j - 1) * lda;
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
    const t = std.math.hypot(y.re, y.im);
    if (t == 0.0) return x;
    return Cx.init(y.re / t, y.im / t).mul(x);
}

/// WLOG: `std.math.complex.log`, with the angle `+pi` on the negative real
/// axis even where `x.im` is `-0.0`.
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

/// A level-1 BLAS result, its strides chosen by the caller: only
/// `error.Canceled` remains.
pub fn level1(r: anytype) std.Io.Cancelable!@typeInfo(@TypeOf(r)).error_union.payload {
    return r catch |e| ops.narrow(std.Io.Cancelable, e);
}

// ---- tests ------------------------------------------------------------------

test "wdiv is the inverse of mul" {
    const a = Cx.init(3.0, -2.0);
    const b = Cx.init(-1.5, 4.0);
    const prod = a.mul(b);
    const back = wdiv(prod, b);
    try zigstd.testExpectApproxEqAbs(a.re, back.re, 1e-12);
    try zigstd.testExpectApproxEqAbs(a.im, back.im, 1e-12);
}
