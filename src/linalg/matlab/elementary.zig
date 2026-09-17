//! Elementary complex scalar functions ported from `matlab88.f90`'s
//! `MATFN2` (the `POWER` command body - `SIN`/`COS`/`EXP`/`SQRT` are
//! `std.math.complex` calls directly, wired in `matrixfun.zig`).
//!
//! Operates on a single `Cx` scalar. For applying one across a whole
//! vector or matrix, see `__root__.zig`'s elementwise (`sin`, `cos`, ...)
//! and `CxMat`'s matrix-function (`sinm`, `cosm`, ...) methods.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1.zig");
const Cx = blas1.Cx;

/// POWER: `x^y = exp(y * log(x))` (principal branch). Same formula as
/// `std.math.complex.pow`, inheriting `exp`/`log`'s own overflow handling.
/// A real function, not `pub const wpow = std.math.complex.pow;`: not
/// currently needed as a first-class value anywhere, but keeping it typed
/// concretely matches `blas1.wlog`/`watan`'s convention and costs nothing.
/// Not yet wired to a REPL-level entry point (the original exposed this via
/// the `**` operator, which - like `*` matrix multiplication - is REPL
/// operator-dispatch logic this port hasn't ported yet, not a named
/// subroutine); callable directly meanwhile.
pub fn wpow(x: Cx, y: Cx) Cx {
    return std.math.complex.pow(x, y);
}

test "wpow matches repeated multiplication for an integer exponent" {
    const x = Cx.init(1.5, -0.5);
    const cubed = blas1.wmul(blas1.wmul(x, x), x);
    const viaPow = wpow(x, Cx.init(3.0, 0.0));
    try zigstd.testExpectApproxEqAbs(cubed.re, viaPow.re, 1e-9);
    try zigstd.testExpectApproxEqAbs(cubed.im, viaPow.im, 1e-9);
}
