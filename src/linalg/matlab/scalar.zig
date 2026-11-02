//! The matlab port's complex scalar helpers from `matlab88.f90` (`WDIV`,
//! `WSIGN`, `WLOG`, `WATAN`), its indexing (`ix`) and its BLAS error adapter
//! (`level1`). Names match the original Fortran, lowercased. Its vector
//! kernels (`WNRM2`, `WASUM`, `IWAMAX`, `WSCAL`, `WAXPY`, `WSWAP`, `RSWAP`)
//! are rstd's BLAS level 1 since M2; `PYTHAG`, `WMUL` and `CABS1` are
//! `std.math.hypot`, `Complex.mul` and the BLAS's `cabs1`.
//!
//! The original's `FLOP(...)` wrapper (a flop-counter for the REPL's status
//! line) is dropped throughout this port; it had no effect on results.

const std = @import("std");

const mulAdd = @import("../vector.zig").mulAdd;

pub const Cx = std.math.Complex(f64);

/// 0-based flat offset of 1-based matrix element `(i, j)` in a column-major
/// array with leading dimension `lda` - the Zig equivalent of Fortran
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

/// A level-1 BLAS result, its strides chosen by the caller: only
/// `error.Canceled` remains.
pub fn level1(r: anytype) std.Io.Cancelable!@typeInfo(@TypeOf(r)).error_union.payload {
    return r catch |e| switch (e) {
        error.Canceled => error.Canceled,
        error.InvalidIncx, error.InvalidIncy => unreachable,
    };
}

const zigstd = @import("../../zigstd.zig");

test "wdiv is the inverse of mul" {
    const a = Cx.init(3.0, -2.0);
    const b = Cx.init(-1.5, 4.0);
    const prod = a.mul(b);
    const back = wdiv(prod, b);
    try zigstd.testExpectApproxEqAbs(a.re, back.re, 1e-12);
    try zigstd.testExpectApproxEqAbs(a.im, back.im, 1e-12);
}
