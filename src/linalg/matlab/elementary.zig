//! POWER, from `matlab88.f90`'s `MATFN2`; SIN, COS, EXP and SQRT are
//! `std.math.complex` calls in `matrixfun.zig`.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const scalar = @import("scalar.zig");
const Cx = scalar.Cx;

/// POWER: `x^y = exp(y log(x))`, principal branch, by `std.math.complex.pow`;
/// the REPL's `**` operator itself is not ported.
pub fn wpow(x: Cx, y: Cx) Cx {
    return std.math.complex.pow(x, y);
}

// ---- tests ------------------------------------------------------------------

test "wpow matches repeated multiplication for an integer exponent" {
    const x = Cx.init(1.5, -0.5);
    const cubed = x.mul(x).mul(x);
    const viaPow = wpow(x, Cx.init(3.0, 0.0));
    try zigstd.testExpectApproxEqAbs(cubed.re, viaPow.re, 1e-9);
    try zigstd.testExpectApproxEqAbs(cubed.im, viaPow.im, 1e-9);
}
