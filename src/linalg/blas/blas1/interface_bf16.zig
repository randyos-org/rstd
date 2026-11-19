//! The level-1 interface's bfloat16 routines (`interface.zig`).

const std = @import("std");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const bf = @import("../ext/bf16/__root__.zig");
const interface = @import("interface.zig");

const Io = std.Io;
const BF16 = types.BF16;
const first = util.first;
const Error = interface.Error;
const spanX = interface.spanX;
const spanXY = interface.spanXY;

/// `SBDOT` (`interface/bf16dot.c`): the dot of two bfloat16 vectors, 0 for
/// `n <= 0`, returned in `f32`.
pub fn sbdot(io: Io, n: isize, x: [*]BF16, incx: isize, y: [*]BF16, incy: isize) Error!f32 {
    _ = io;
    if (n <= 0) return 0;
    try spanXY(BF16, BF16, n, incx, incy);
    return blas.kernel(.sbdot_k)(n, first(BF16, x, n, incx), incx, first(BF16, y, n, incy), incy);
}

/// `BSCAL` (`interface/scal.c` built `BGEMM`): `x := alpha*x` over bfloat16.
/// Nothing for `incx <= 0`, and nothing to do for `alpha == 1`; like `scal`
/// it passes the kernel flag 1.
pub fn bscal(io: Io, n: isize, alpha: BF16, x: [*]BF16, incx: isize) Error!void {
    _ = io;
    if (incx <= 0 or n <= 0) return;
    try spanX(BF16, n, incx);
    if (bf.scalarToF32(alpha) == 1) return;
    _ = blas.kernel(.bscal_k)(n, 0, 0, alpha, x, incx, x, 0, x, 1);
}

/// `SBSTOBF16` and `SBDTOBF16` (`interface/tobf16.c`): `out := bf16(in)`,
/// rounding to nearest even; `f32` or `f64` in. The errors name `in` as
/// `x` and `out` as `y`.
pub fn sbtobf16(io: Io, comptime F: type, n: isize, in: [*]F, inc_in: isize, out: [*]BF16, inc_out: isize) Error!void {
    _ = io;
    if (n <= 0) return;
    const kern = if (comptime F == f32)
        blas.kernel(.sbstobf16_k)
    else if (comptime F == f64)
        blas.kernel(.sbdtobf16_k)
    else
        @compileError("sbtobf16: " ++ @typeName(F) ++ " is not f32 or f64");
    try spanXY(F, BF16, n, inc_in, inc_out);
    kern(n, first(F, in, n, inc_in), inc_in, first(BF16, out, n, inc_out), inc_out);
}

/// `SBF16TOS` and `DBF16TOD` (`interface/bf16to.c`): `out := in` widened to
/// `f32` or `f64`; errors as `sbtobf16`'s.
pub fn bf16to(io: Io, comptime F: type, n: isize, in: [*]BF16, inc_in: isize, out: [*]F, inc_out: isize) Error!void {
    _ = io;
    if (n <= 0) return;
    const kern = if (comptime F == f32)
        blas.kernel(.sbf16tos_k)
    else if (comptime F == f64)
        blas.kernel(.dbf16tod_k)
    else
        @compileError("bf16to: " ++ @typeName(F) ++ " is not f32 or f64");
    try spanXY(BF16, F, n, inc_in, inc_out);
    kern(n, first(BF16, in, n, inc_in), inc_in, first(F, out, n, inc_out), inc_out);
}
