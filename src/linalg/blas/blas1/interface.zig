//! The level-1 interface: the checks and pointer set-up OpenBLAS's
//! `interface/*.c` do before calling a kernel. Each takes the element type `T`
//! explicitly, any float type or its `Complex`; `n` and strides are `isize`.
//! A type with vtable slots (`f32`, `f64`, `XDouble` and their `Complex`)
//! runs the core's kernels, any other the portable ones built for it
//! (`portable.kernel`).
//!
//! A negative stride names the vector's last element first: the interface
//! moves the pointer to it, and the kernel walks back from there. Kernels
//! assume `n >= 1`; every routine here returns before calling one otherwise.
//!
//! A vector that could not fit the address space returns an `ArgError`
//! naming its stride, `x`'s checked first; there is no xerbla.

const std = @import("std");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const check = @import("check.zig");
const portable = @import("../portable.zig");

const BlasLong = types.BlasLong;
const BF16 = types.BF16;
const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

const isComplex = types.isComplex;
const Float = types.Float;

/// `T`'s elements as the kernels take them: floats, two per complex.
fn raw(comptime T: type, p: [*]T) [*]Float(T) {
    return @ptrCast(p);
}

/// Bad arguments, by the reference's names.
pub const ArgError = error{ InvalidIncx, InvalidIncy };

/// Where a vector of `n` with stride `inc` starts: its last element when
/// `inc < 0`. The vector must have passed `spanX`/`spanXY`.
fn first(comptime T: type, p: [*]T, n: BlasLong, inc: BlasLong) [*]T {
    return if (inc < 0) @ptrCast(util.at(T, p, -(n - 1) * inc)) else p;
}

/// Refuses an `x` that could not exist in memory (`util.fits`), so no index
/// overflows.
fn spanX(comptime T: type, n: isize, incx: isize) ArgError!void {
    if (!util.fits(n, incx, @sizeOf(T))) return error.InvalidIncx;
}

/// `spanX` for `x` of `X`, then `y` of `Y`.
fn spanXY(comptime X: type, comptime Y: type, n: isize, incx: isize, incy: isize) ArgError!void {
    try spanX(X, n, incx);
    if (!util.fits(n, incy, @sizeOf(Y))) return error.InvalidIncy;
}

// ---- copy, swap, scal --------------------------------------------------------

/// `?COPY` (`interface/copy.c`): `y := x`.
pub fn copy(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) ArgError!void {
    if (n <= 0) return;
    try spanXY(T, T, n, incx, incy);
    _ = portable.fetch(T, .copy_k, blas.vtable)(n, raw(T, first(T, x, n, incx)), incx, raw(T, first(T, y, n, incy)), incy);
}

/// `?SWAP` (`interface/swap.c`, `zswap.c`): `x <-> y`.
pub fn swap(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) ArgError!void {
    if (n <= 0) return;
    try spanXY(T, T, n, incx, incy);
    const px = raw(T, first(T, x, n, incx));
    const py = raw(T, first(T, y, n, incy));
    if (comptime isComplex(T)) {
        _ = portable.fetch(T, .swap_k, blas.vtable)(n, 0, 0, 0, 0, px, incx, py, incy, py, 0);
    } else {
        _ = portable.fetch(T, .swap_k, blas.vtable)(n, 0, 0, 0, px, incx, py, incy, py, 0);
    }
}

/// `?SCAL` (`interface/scal.c`, `zscal.c`): `x := alpha*x`. Nothing for
/// `incx <= 0`, and nothing to do for `alpha == 1`. With the flag the
/// interface passes (1), a zero `alpha` leaves NaN and infinity as NaN, as
/// `alpha*x` would, and finite elements 0.
pub fn scal(comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize) ArgError!void {
    if (n <= 0 or incx <= 0) return;
    try spanX(T, n, incx);
    const px = raw(T, x);
    if (comptime isComplex(T)) {
        if (alpha.re == 1 and alpha.im == 0) return;
        _ = portable.fetch(T, .scal_k, blas.vtable)(n, 0, 0, alpha.re, alpha.im, px, incx, px, 0, px, 1);
    } else {
        if (alpha == 1) return;
        _ = portable.fetch(T, .scal_k, blas.vtable)(n, 0, 0, alpha, px, incx, px, 0, px, 1);
    }
}

/// `CSSCAL`/`ZDSCAL`: complex `x` scaled by a real `alpha`.
pub fn scalReal(comptime T: type, n: isize, alpha: Float(T), x: [*]T, incx: isize) ArgError!void {
    if (comptime !isComplex(T)) @compileError("scalReal: " ++ @typeName(T) ++ " is not complex");
    return scal(T, n, .init(alpha, 0), x, incx);
}

// ---- dot ---------------------------------------------------------------------

/// `?DOT`/`?DOTU` (`interface/dot.c`, `zdot.c`): `x . y`, 0 for `n <= 0`.
pub fn dot(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) ArgError!T {
    return dotImpl(T, false, n, x, incx, y, incy);
}

/// `?DOTC`: `conj(x) . y`. Complex only.
pub fn dotc(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) ArgError!T {
    if (comptime !isComplex(T)) @compileError("dotc: " ++ @typeName(T) ++ " is not complex");
    return dotImpl(T, true, n, x, incx, y, incy);
}

fn dotImpl(comptime T: type, comptime conj: bool, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) ArgError!T {
    if (n <= 0) return zero(T);
    try spanXY(T, T, n, incx, incy);
    const px = raw(T, first(T, x, n, incx));
    const py = raw(T, first(T, y, n, incy));
    if (comptime isComplex(T)) {
        return types.complexOf(portable.fetch(T, if (conj) .dotc_k else .dotu_k, blas.vtable)(n, px, incx, py, incy));
    }
    return portable.fetch(T, .dot_k, blas.vtable)(n, px, incx, py, incy);
}

/// `DSDOT` (`interface/dsdot.c`): the dot of two `f32` vectors, summed and
/// returned in `f64`.
pub fn dsdot(n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize) ArgError!f64 {
    if (n <= 0) return 0;
    try spanXY(f32, f32, n, incx, incy);
    return blas.kernel(.dsdot_k)(n, first(f32, x, n, incx), incx, first(f32, y, n, incy), incy);
}

/// `SDSDOT` (`interface/sdsdot.c`): `alpha + x . y`, summed in `f64` and
/// returned as `f32`; `alpha` alone for `n <= 0`.
pub fn sdsdot(n: isize, alpha: f32, x: [*]f32, incx: isize, y: [*]f32, incy: isize) ArgError!f32 {
    if (n <= 0) return alpha;
    return @floatCast(try dsdot(n, x, incx, y, incy) + alpha);
}

// ---- norms and sums ----------------------------------------------------------

/// `?NRM2`/`SCNRM2`/`DZNRM2` (`interface/nrm2.c`): the Euclidean norm. For
/// `incx == 0` it is `sqrt(n) * |x[0]|`, worked out here as the reference
/// does; a negative stride is the same vector backwards.
pub fn nrm2(comptime T: type, n: isize, x: [*]T, incx: isize) ArgError!Float(T) {
    const F = Float(T);
    if (n <= 0) return 0;
    // the count and its root in the accumulator type: `n` may be past `F`'s range
    const A = types.Acc(F);
    const fn_: A = @floatFromInt(n);
    if (comptime !isComplex(T)) {
        if (n == 1) return @abs(x[0]);
        if (incx == 0) return @floatCast(@sqrt(fn_) * @abs(@as(A, x[0])));
    } else if (incx == 0) {
        const fr: A = @abs(x[0].re);
        const fi: A = @abs(x[0].im);
        const lo = @min(fr, fi);
        const hi = @max(fr, fi);
        if (hi == 0) return 0;
        if (hi == lo) return @floatCast(@sqrt(fn_) * @sqrt(@as(A, 2)) * hi);
        return @floatCast(@sqrt(fn_) * hi * @sqrt(1 + (lo / hi) * (lo / hi)));
    }
    try spanX(T, n, incx);
    return portable.fetch(T, .nrm2_k, blas.vtable)(n, raw(T, first(T, x, n, incx)), incx);
}

/// `?ASUM`/`SCASUM`/`DZASUM` (`interface/asum.c`): the sum of `|x|`, of
/// `|re| + |im|` for complex. 0 for `n <= 0`; the kernel gives 0 for
/// `incx <= 0`.
pub fn asum(comptime T: type, n: isize, x: [*]T, incx: isize) ArgError!Float(T) {
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    return portable.fetch(T, .asum_k, blas.vtable)(n, raw(T, x), incx);
}

/// `?SUM` (`interface/sum.c`, an OpenBLAS extension): the sum of `x`, of
/// `re + im` for complex.
pub fn sum(comptime T: type, n: isize, x: [*]T, incx: isize) ArgError!Float(T) {
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    return portable.fetch(T, .sum_k, blas.vtable)(n, raw(T, x), incx);
}

// ---- extremes ----------------------------------------------------------------

/// Which extreme: by magnitude (`|re| + |im|` for complex) or by value.
pub const Extreme = enum { amax, amin, max, min };

/// `I?AMAX`, `I?AMIN`, `I?MAX`, `I?MIN` (`interface/imax.c`): the 1-based
/// index of the first extreme element; 0 for `n <= 0` (and, from the
/// kernel, for `incx <= 0`). `max` and `min` are real only.
pub fn index(comptime T: type, comptime which: Extreme, n: isize, x: [*]T, incx: isize) ArgError!isize {
    comptime realOnly(T, which);
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    const ret: BlasLong = portable.fetch(T, @field(portable.Op, "i" ++ @tagName(which) ++ "_k"), blas.vtable)(n, raw(T, x), incx);
    return @intCast(@min(ret, n));
}

/// `?AMAX`, `?AMIN`, `?MAX`, `?MIN` (`interface/max.c`, OpenBLAS extensions):
/// the extreme value itself, `|x|` for `amax`/`amin`. 0 for `n <= 0`; for
/// `incx == 0`, `x[0]` (by magnitude for `amax`/`amin`).
pub fn extreme(comptime T: type, comptime which: Extreme, n: isize, x: [*]T, incx: isize) ArgError!Float(T) {
    comptime realOnly(T, which);
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    if (incx == 0) {
        if (comptime isComplex(T)) return @abs(x[0].re) + @abs(x[0].im);
        return if (which == .amax or which == .amin) @abs(x[0]) else x[0];
    }
    return portable.fetch(T, @field(portable.Op, @tagName(which) ++ "_k"), blas.vtable)(n, raw(T, x), incx);
}

fn realOnly(comptime T: type, comptime which: Extreme) void {
    if (isComplex(T) and (which == .max or which == .min)) @compileError(@tagName(which) ++ ": " ++ @typeName(T) ++ " is not real");
}

// ---- rotations ---------------------------------------------------------------

/// `?ROT`, `CSROT`/`ZDROT` (`interface/rot.c`, `zrot.c`): the plane rotation
/// `x := c*x + s*y`, `y := c*y - s*x`, with real `c` and `s`.
pub fn rot(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, c: Float(T), s: Float(T)) ArgError!void {
    if (n <= 0) return;
    try spanXY(T, T, n, incx, incy);
    const px = raw(T, first(T, x, n, incx));
    const py = raw(T, first(T, y, n, incy));
    _ = portable.fetch(T, .rot_k, blas.vtable)(n, px, incx, py, incy, c, s);
}

/// `?ROTM` (`interface/rotm.c`): the modified rotation `param` describes
/// (`param[0]` the flag). The reference passes everything to the kernel,
/// which handles negative strides itself; nothing for `n <= 0`.
pub fn rotm(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, param: *const [5]T) ArgError!void {
    if (comptime isComplex(T)) @compileError("rotm: " ++ @typeName(T) ++ " is not real");
    if (n <= 0) return;
    try spanXY(T, T, n, incx, incy);
    _ = portable.fetch(T, .rotm_k, blas.vtable)(n, x, incx, y, incy, @constCast(param));
}

// ---- axpy, axpby -------------------------------------------------------------

/// `?AXPY` (`interface/axpy.c`, `zaxpy.c`): `y := alpha*x + y`. Nothing to
/// do for `alpha == 0`; both strides zero is `y[0] += n*alpha*x[0]`.
pub fn axpy(comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize) ArgError!void {
    return axpyImpl(T, false, n, alpha, x, incx, y, incy);
}

/// `?AXPYC`: `y := alpha*conj(x) + y`. Complex only.
pub fn axpyc(comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize) ArgError!void {
    if (comptime !isComplex(T)) @compileError("axpyc: " ++ @typeName(T) ++ " is not complex");
    return axpyImpl(T, true, n, alpha, x, incx, y, incy);
}

fn axpyImpl(comptime T: type, comptime conj: bool, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize) ArgError!void {
    if (n <= 0) return;
    if (comptime isComplex(T)) {
        if (alpha.re == 0 and alpha.im == 0) return;
        if (incx == 0 and incy == 0) {
            // the reference does this unconjugated for `axpyc` too
            // the count in the accumulator type: `n` may be past `F`'s range
            const A = types.Acc(Float(T));
            const cn: A = @floatFromInt(n);
            y[0].re += @floatCast(cn * (alpha.re * x[0].re - alpha.im * x[0].im));
            y[0].im += @floatCast(cn * (alpha.im * x[0].re + alpha.re * x[0].im));
            return;
        }
        try spanXY(T, T, n, incx, incy);
        const px = raw(T, first(T, x, n, incx));
        const py = raw(T, first(T, y, n, incy));
        _ = portable.fetch(T, if (conj) .axpyc_k else .axpy_k, blas.vtable)(n, 0, 0, alpha.re, alpha.im, px, incx, py, incy, py, 0);
    } else {
        if (alpha == 0) return;
        if (incx == 0 and incy == 0) {
            y[0] += @floatCast(@as(types.Acc(T), @floatFromInt(n)) * alpha * x[0]);
            return;
        }
        try spanXY(T, T, n, incx, incy);
        const px = first(T, x, n, incx);
        const py = first(T, y, n, incy);
        _ = portable.fetch(T, .axpy_k, blas.vtable)(n, 0, 0, alpha, px, incx, py, incy, py, 0);
    }
}

/// `?AXPBY` (`interface/axpby.c`, `zaxpby.c`, an OpenBLAS extension):
/// `y := alpha*x + beta*y`.
pub fn axpby(comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) ArgError!void {
    if (n <= 0) return;
    try spanXY(T, T, n, incx, incy);
    const px = raw(T, first(T, x, n, incx));
    const py = raw(T, first(T, y, n, incy));
    if (comptime isComplex(T)) {
        _ = portable.fetch(T, .axpby_k, blas.vtable)(n, alpha.re, alpha.im, px, incx, beta.re, beta.im, py, incy);
    } else {
        _ = portable.fetch(T, .axpby_k, blas.vtable)(n, alpha, px, incx, beta, py, incy);
    }
}

// ---- bfloat16 ----------------------------------------------------------------

/// `SBDOT` (`interface/bf16dot.c`): the dot of two bfloat16 vectors, 0 for
/// `n <= 0`, returned in `f32`.
pub fn sbdot(n: isize, x: [*]BF16, incx: isize, y: [*]BF16, incy: isize) ArgError!f32 {
    if (n <= 0) return 0;
    try spanXY(BF16, BF16, n, incx, incy);
    return blas.kernel(.sbdot_k)(n, first(BF16, x, n, incx), incx, first(BF16, y, n, incy), incy);
}

/// `BSCAL` (`interface/scal.c` built `BGEMM`): `x := alpha*x` over bfloat16.
/// Nothing for `incx <= 0`, and nothing to do for `alpha == 1`; like `scal`
/// it passes the kernel flag 1.
pub fn bscal(n: isize, alpha: BF16, x: [*]BF16, incx: isize) ArgError!void {
    if (incx <= 0 or n <= 0) return;
    try spanX(BF16, n, incx);
    if (@import("../ext/bf16/__root__.zig").scalarToF32(alpha) == 1) return;
    _ = blas.kernel(.bscal_k)(n, 0, 0, alpha, x, incx, x, 0, x, 1);
}

/// `SBSTOBF16` and `SBDTOBF16` (`interface/tobf16.c`): `out := bf16(in)`,
/// rounding to nearest even; `f32` or `f64` in. The errors name `in` as
/// `x` and `out` as `y`.
pub fn sbtobf16(comptime F: type, n: isize, in: [*]F, inc_in: isize, out: [*]BF16, inc_out: isize) ArgError!void {
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
pub fn bf16to(comptime F: type, n: isize, in: [*]BF16, inc_in: isize, out: [*]F, inc_out: isize) ArgError!void {
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

// ---- rotation set-up ---------------------------------------------------------

/// `?ROTG` (`interface/rotg.c`): the rotation `[c s; -s c]` that zeroes `b`
/// against `a`. `a` becomes `r`, `b` the reconstruction value `z`.
/// Computed in `XDouble`, as the reference does in `long double` on x86 (x87
/// where `XDouble` is), or in `T` where that is wider; a narrower input is
/// rounded once on the way out.
pub fn rotg(comptime T: type, a: *T, b: *T, c: *T, s: *T) void {
    if (comptime isComplex(T)) @compileError("rotg: use zrotg for " ++ @typeName(T));
    const X = if (@bitSizeOf(T) > @bitSizeOf(XDouble)) T else XDouble;
    const da: X = a.*;
    const db: X = b.*;
    const safmin: X = std.math.floatMin(T);
    const ada = @abs(da);
    const adb = @abs(db);
    if (adb == 0) {
        c.* = 1;
        s.* = 0;
        b.* = 0;
        return;
    }
    if (ada == 0) {
        c.* = 0;
        s.* = 1;
        a.* = b.*;
        b.* = 1;
        return;
    }
    const sc = @min(@max(safmin, @max(ada, adb)), 1 / safmin);
    const sigma = std.math.copysign(@as(X, 1), if (ada > adb) da else db);
    const aa = da / sc;
    const bb = db / sc;
    const r = sigma * sc * @sqrt(aa * aa + bb * bb);
    const cx = da / r;
    const sx = db / r;
    var z: X = 1;
    if (ada > adb) z = sx;
    if (ada <= adb and cx != 0) z = 1 / cx;
    c.* = @floatCast(cx);
    s.* = @floatCast(sx);
    a.* = @floatCast(r);
    b.* = @floatCast(z);
}

/// `CROTG`/`ZROTG`: the rotation with real `c` and complex `s` that zeroes
/// `b` against `a`, `[c s; -conj(s) c] [a; b] = [r; 0]`; `a` becomes `r`.
/// Follows LAPACK's `zlartg`, scaling against overflow and underflow.
/// `interface/zrotg.c` transcribes the same algorithm but its scaled branches
/// are wrong (square roots where `zlartg` has squared norms, `a` rescaled by
/// the wrong factor, a multiply where a divide belongs), so results there
/// are wrong in the last bits and sometimes worse; this is the sound
/// algorithm, not a copy of `zrotg.c`, and is tested against a wide-precision
/// reference over all its branches.
pub fn zrotg(comptime T: type, a: *T, b: T, c: *Float(T), s: *T) void {
    if (comptime !isComplex(T)) @compileError("zrotg: " ++ @typeName(T) ++ " is not complex");
    const F = Float(T);
    const safmin = std.math.floatMin(F);
    const safmax = 1 / safmin;
    const rtmin = @sqrt(safmin);
    const f = a.*;
    const g = b;
    if (g.re == 0 and g.im == 0) {
        c.* = 1;
        s.* = .init(0, 0);
        return;
    }
    if (f.re == 0 and f.im == 0) {
        c.* = 0;
        if (g.re == 0 or g.im == 0) {
            const d = @abs(if (g.re == 0) g.im else g.re);
            s.* = scale(T, g.conjugate(), 1 / d);
            a.* = .init(d, 0);
            return;
        }
        const g1 = @max(@abs(g.re), @abs(g.im));
        if (g1 > rtmin and g1 < @sqrt(safmax / 2)) {
            const d = @sqrt(abssq(T, g));
            s.* = scale(T, g.conjugate(), 1 / d);
            a.* = .init(d, 0);
        } else {
            const u = @min(safmax, @max(safmin, g1));
            const gs = scale(T, g, 1 / u);
            const d = @sqrt(abssq(T, gs));
            s.* = scale(T, gs.conjugate(), 1 / d);
            a.* = .init(d * u, 0);
        }
        return;
    }
    const f1 = @max(@abs(f.re), @abs(f.im));
    const g1 = @max(@abs(g.re), @abs(g.im));
    var rtmax = @sqrt(safmax / 4);
    // unscaled, or scaled by `u` (and `f` by `v`, as `w = v/u`)
    var u: F = 1;
    var w: F = 1;
    var fs = f;
    var gs = g;
    if (!(f1 > rtmin and f1 < rtmax and g1 > rtmin and g1 < rtmax)) {
        u = @min(safmax, @max(safmin, @max(f1, g1)));
        gs = scale(T, g, 1 / u);
        if (f1 / u < rtmin) {
            const v = @min(safmax, @max(safmin, f1));
            w = v / u;
            fs = scale(T, f, 1 / v);
        } else {
            fs = scale(T, f, 1 / u);
        }
    }
    const f2 = abssq(T, fs);
    const h2 = f2 * w * w + abssq(T, gs);
    var r: T = undefined;
    if (f2 >= h2 * safmin) {
        c.* = @sqrt(f2 / h2);
        r = scale(T, fs, 1 / c.*);
        rtmax *= 2;
        s.* = if (f2 > rtmin and h2 < rtmax)
            gs.conjugate().mul(scale(T, fs, 1 / @sqrt(f2 * h2)))
        else
            gs.conjugate().mul(scale(T, r, 1 / h2));
    } else {
        // intermediates might overflow
        const d = @sqrt(f2 * h2);
        c.* = f2 / d;
        r = if (c.* >= safmin) scale(T, fs, 1 / c.*) else scale(T, fs, h2 / d);
        s.* = gs.conjugate().mul(scale(T, fs, 1 / d));
    }
    c.* *= w;
    a.* = scale(T, r, u);
}

fn scale(comptime T: type, z: T, f: Float(T)) T {
    return .init(z.re * f, z.im * f);
}

fn abssq(comptime T: type, z: T) Float(T) {
    return z.re * z.re + z.im * z.im;
}

/// `?ROTMG` (`interface/rotmg.c`): the modified rotation that zeroes the
/// second component of `(sqrt(d1)*x1, sqrt(d2)*y1)`, as `param` for `rotm`.
/// `d1`, `d2` and `x1` are updated.
pub fn rotmg(comptime T: type, d1: *T, d2: *T, x1: *T, y1: T, param: *[5]T) void {
    if (comptime isComplex(T)) @compileError("rotmg: " ++ @typeName(T) ++ " is not real");
    // the reference's bounds, 4096 and its square; a narrow type (`f16`)
    // keeps `d1`, `d2` within its range by 16 and 256 instead
    const narrow = @bitSizeOf(T) < 32;
    const gam: T = if (narrow) 16 else 4096;
    const gamsq: T = if (narrow) 256 else 16777216;
    const rgamsq: T = if (narrow) 1.0 / 256.0 else 5.9604645e-8;
    var h11: T = 0;
    var h21: T = 0;
    var h12: T = 0;
    var h22: T = 0;
    var flag: T = -1;
    if (d2.* == 0 or y1 == 0) {
        param[0] = -2;
        return;
    }
    if (d1.* < 0) {
        d1.* = 0;
        d2.* = 0;
        x1.* = 0;
    } else if ((d1.* == 0 or x1.* == 0) and d2.* > 0) {
        flag = 1;
        h12 = 1;
        h21 = -1;
        x1.* = y1;
        std.mem.swap(T, d1, d2);
    } else {
        const p2 = d2.* * y1;
        const p1 = d1.* * x1.*;
        const q2 = p2 * y1;
        const q1 = p1 * x1.*;
        if (@abs(q1) > @abs(q2)) {
            h11 = 1;
            h22 = 1;
            h21 = -y1 / x1.*;
            h12 = p2 / p1;
            const u = 1 - h12 * h21;
            if (u <= 0) {
                // only through rounding, as in the reference
                h11 = 0;
                h12 = 0;
                h21 = 0;
                h22 = 0;
                d1.* = 0;
                d2.* = 0;
                x1.* = 0;
            } else {
                flag = 0;
                d1.* /= u;
                d2.* /= u;
                x1.* *= u;
            }
        } else if (q2 < 0) {
            d1.* = 0;
            d2.* = 0;
            x1.* = 0;
        } else {
            flag = 1;
            h21 = -1;
            h12 = 1;
            h11 = p1 / p2;
            h22 = x1.* / y1;
            const u = 1 + h11 * h22;
            const t = d2.* / u;
            d2.* = d1.* / u;
            d1.* = t;
            x1.* = y1 * u;
        }
        while (d1.* <= rgamsq and d1.* != 0) {
            flag = -1;
            d1.* *= gam * gam;
            x1.* /= gam;
            h11 /= gam;
            h12 /= gam;
        }
        while (@abs(d1.*) > gamsq) {
            flag = -1;
            d1.* /= gam * gam;
            x1.* *= gam;
            h11 *= gam;
            h12 *= gam;
        }
        while (@abs(d2.*) <= rgamsq and d2.* != 0) {
            flag = -1;
            d2.* *= gam * gam;
            h21 /= gam;
            h22 /= gam;
        }
        while (@abs(d2.*) > gamsq) {
            flag = -1;
            d2.* /= gam * gam;
            h21 *= gam;
            h22 *= gam;
        }
    }
    if (flag < 0) {
        param[1] = h11;
        param[2] = h21;
        param[3] = h12;
        param[4] = h22;
    } else if (flag == 0) {
        param[2] = h21;
        param[3] = h12;
    } else {
        param[1] = h11;
        param[4] = h22;
    }
    param[0] = flag;
}

// ---- tests -------------------------------------------------------------------

const reals = .{ f32, f64 } ++ check.x_real;
const all = .{ f32, f64, Complex32, Complex64 } ++ check.x_real ++ check.x_complex;

fn zero(comptime T: type) T {
    return if (comptime isComplex(T)) .init(0, 0) else 0;
}

fn same(comptime T: type, want: []const T, got: []const T) !void {
    try check.expectCloseSlices(T, want, got, check.tolerance(Float(T), 4));
}

fn sameInt(want: isize, got: isize) !void {
    if (want == got) return;
    check.say("got {d}, want {d}: ", .{ got, want });
    return error.TestExpectedEqual;
}

test "interface: n <= 0 changes nothing and returns zero" {
    inline for (all) |T| {
        var x: [4]T = undefined;
        var y: [4]T = undefined;
        check.fill(T, &x, 1);
        check.fill(T, &y, 2);
        const x0 = x;
        const y0 = y;
        for ([_]isize{ 0, -1 }) |n| {
            try copy(T, n, &x, 1, &y, 1);
            try swap(T, n, &x, 1, &y, 1);
            try scal(T, n, zero(T), &x, 1);
            // the reference has no extended axpby
            if (comptime !types.isExtended(Float(T))) try axpby(T, n, zero(T), &x, 1, zero(T), &y, 1);
            try rot(T, n, &x, 1, &y, 1, 0, 1);
            try same(T, &x0, &x);
            try same(T, &y0, &y);
            try check.expectClose(T, zero(T), try dot(T, n, &x, 1, &y, 1), 0);
            try check.expectClose(Float(T), 0, try nrm2(T, n, &x, 1), 0);
            try check.expectClose(Float(T), 0, try asum(T, n, &x, 1), 0);
            try check.expectClose(Float(T), 0, try sum(T, n, &x, 1), 0);
            try check.expectClose(Float(T), 0, try extreme(T, .amax, n, &x, 1), 0);
            try sameInt(0, try index(T, .amax, n, &x, 1));
        }
    }
    try check.expectClose(f32, 2.5, try sdsdot(0, 2.5, undefined, 1, undefined, 1), 0);
}

test "interface: a vector past the address space is its stride's error, x's first" {
    const testing = std.testing;
    // two elements of 2 bytes or more at this stride overflow
    const big = std.math.maxInt(isize) / 2;
    const bad = [_]isize{ big, -big, std.math.minInt(isize) };
    inline for (all) |T| {
        var x: [2]T = undefined;
        var y: [2]T = undefined;
        check.fill(T, &x, 1);
        check.fill(T, &y, 2);
        const two: T = if (comptime isComplex(T)) .init(2, 0) else 2;
        const extended = comptime types.isExtended(Float(T));
        for (bad) |b| {
            const cases = [_]struct { ix: isize, iy: isize, err: ArgError }{
                .{ .ix = b, .iy = 1, .err = error.InvalidIncx },
                .{ .ix = 1, .iy = b, .err = error.InvalidIncy },
                .{ .ix = b, .iy = b, .err = error.InvalidIncx },
            };
            for (cases) |c| {
                try testing.expectError(c.err, copy(T, 2, &x, c.ix, &y, c.iy));
                try testing.expectError(c.err, swap(T, 2, &x, c.ix, &y, c.iy));
                try testing.expectError(c.err, dot(T, 2, &x, c.ix, &y, c.iy));
                try testing.expectError(c.err, rot(T, 2, &x, c.ix, &y, c.iy, 0.6, 0.8));
                try testing.expectError(c.err, axpy(T, 2, two, &x, c.ix, &y, c.iy));
                if (comptime !extended) try testing.expectError(c.err, axpby(T, 2, two, &x, c.ix, two, &y, c.iy));
                if (comptime isComplex(T)) {
                    try testing.expectError(c.err, dotc(T, 2, &x, c.ix, &y, c.iy));
                    try testing.expectError(c.err, axpyc(T, 2, two, &x, c.ix, &y, c.iy));
                } else {
                    try testing.expectError(c.err, rotm(T, 2, &x, c.ix, &y, c.iy, &.{ -1, 2, 3, -1, 0.5 }));
                }
            }
            // `scal` and `bscal` skip a negative stride before checking it
            if (b > 0) try testing.expectError(error.InvalidIncx, scal(T, 2, two, &x, b));
            try testing.expectError(error.InvalidIncx, nrm2(T, 2, &x, b));
            try testing.expectError(error.InvalidIncx, asum(T, 2, &x, b));
            try testing.expectError(error.InvalidIncx, sum(T, 2, &x, b));
            try testing.expectError(error.InvalidIncx, index(T, .amax, 2, &x, b));
            try testing.expectError(error.InvalidIncx, extreme(T, .amax, 2, &x, b));
            // the early returns come first
            try copy(T, 0, &x, b, &y, b);
            try testing.expectEqual(zero(T), try dot(T, 0, &x, b, &y, b));
            if (comptime !isComplex(T)) try rotm(T, 0, &x, b, &y, b, &.{ -1, 2, 3, -1, 0.5 });
        }
        // nothing was written
        var x0: [2]T = undefined;
        var y0: [2]T = undefined;
        check.fill(T, &x0, 1);
        check.fill(T, &y0, 2);
        try same(T, &x0, &x);
        try same(T, &y0, &y);
    }
    var f = [_]f32{ 1, 2 };
    var d = [_]f64{ 1, 2 };
    var h = [_]BF16{ b16(1), b16(2) };
    for (bad) |b| {
        try testing.expectError(error.InvalidIncx, dsdot(2, &f, b, &f, 1));
        try testing.expectError(error.InvalidIncy, dsdot(2, &f, 1, &f, b));
        try testing.expectError(error.InvalidIncx, sdsdot(2, 0.5, &f, b, &f, b));
        try testing.expectError(error.InvalidIncx, sbdot(2, &h, b, &h, 1));
        try testing.expectError(error.InvalidIncy, sbdot(2, &h, 1, &h, b));
        try testing.expectError(error.InvalidIncx, sbtobf16(f32, 2, &f, b, &h, 1));
        try testing.expectError(error.InvalidIncy, sbtobf16(f64, 2, &d, 1, &h, b));
        try testing.expectError(error.InvalidIncx, bf16to(f32, 2, &h, b, &f, 1));
        try testing.expectError(error.InvalidIncy, bf16to(f64, 2, &h, 1, &d, b));
        if (b > 0) try testing.expectError(error.InvalidIncx, bscal(2, b16(2), &h, b));
    }
}

test "interface: a negative stride walks the vector from its end" {
    inline for (all) |T| {
        const n = 7;
        var x: [n * 2]T = undefined;
        var y: [n * 3]T = undefined;
        check.fill(T, &x, 3);
        check.fill(T, &y, 4);
        // `xr` is x as stride -2 sees it, stored forwards at stride 1
        var xr: [n]T = undefined;
        for (0..n) |i| xr[i] = x[(n - 1 - i) * 2];
        var y1 = y;
        var y2 = y;
        try copy(T, n, &x, -2, &y1, 3);
        try copy(T, n, &xr, 1, &y2, 3);
        try same(T, &y2, &y1);
        try check.expectClose(T, try dot(T, n, &xr, 1, &y, 3), try dot(T, n, &x, -2, &y, 3), check.tolerance(Float(T), 4 * n));
        try check.expectClose(Float(T), try nrm2(T, n, &xr, 1), try nrm2(T, n, &x, -2), check.tolerance(Float(T), 4 * n));
    }
}

test "interface: scal does nothing for incx <= 0 or alpha 1; a zero alpha keeps NaN" {
    inline for (all) |T| {
        var x: [4]T = undefined;
        check.fill(T, &x, 5);
        const x0 = x;
        const two: T = if (comptime isComplex(T)) .init(2, 0) else 2;
        const one: T = if (comptime isComplex(T)) .init(1, 0) else 1;
        try scal(T, 4, two, &x, 0);
        try scal(T, 4, two, &x, -1);
        try scal(T, 4, one, &x, 1);
        try same(T, &x0, &x);
        const nan = std.math.nan(Float(T));
        x[1] = if (comptime isComplex(T)) .init(nan, nan) else nan;
        try scal(T, 4, zero(T), &x, 1);
        const at1 = if (comptime isComplex(T)) x[1].re else x[1];
        if (!std.math.isNan(at1)) return error.TestExpectedNan;
        x[1] = zero(T);
        try same(T, &@as([4]T, @splat(zero(T))), &x);
    }
    var z = [_]Complex64{ .init(1, 2), .init(-3, 4) };
    try scalReal(Complex64, 2, 2, &z, 1);
    try same(Complex64, &.{ .init(2, 4), .init(-6, 8) }, &z);
}

test "interface: dotc conjugates x, and dsdot/sdsdot sum in f64" {
    var x = [_]Complex32{ .init(1, 2), .init(3, -1) };
    var y = [_]Complex32{ .init(2, 1), .init(-1, 1) };
    // conj(1+2i)(2+i) + conj(3-i)(-1+i) = (4-3i) + (-4+2i)
    try check.expectClose(Complex32, .init(0, -1), try dotc(Complex32, 2, &x, 1, &y, 1), 1e-6);
    // (1+2i)(2+i) + (3-i)(-1+i) = 5i + (-2+4i)
    try check.expectClose(Complex32, .init(-2, 9), try dot(Complex32, 2, &x, 1, &y, 1), 1e-6);
    // 1e8 + 1 - 1e8 is 1 in f64 but 0 in f32
    var a = [_]f32{ 1e4, 1, -1e4 };
    var b = [_]f32{ 1e4, 1, 1e4 };
    try check.expectClose(f64, 1, try dsdot(3, &a, 1, &b, 1), 0);
    try check.expectClose(f32, 1.5, try sdsdot(3, 0.5, &a, 1, &b, 1), 0);
}

test "interface: nrm2 of a zero stride is sqrt(n) * |x[0]|" {
    var x = [_]f64{-3};
    try check.expectClose(f64, 3 * @sqrt(4.0), try nrm2(f64, 4, &x, 0), 1e-12);
    var z = [_]Complex64{.init(3, -4)};
    try check.expectClose(f64, 5 * @sqrt(4.0), try nrm2(Complex64, 4, &z, 0), 1e-12);
    var w = [_]Complex64{.init(2, -2)};
    try check.expectClose(f64, @sqrt(9.0) * @sqrt(2.0) * 2, try nrm2(Complex64, 9, &w, 0), 1e-12);
}

test "interface: extremes, by magnitude and by value" {
    var x = [_]f64{ 3, -7, 7, 1, -2 };
    try sameInt(2, try index(f64, .amax, 5, &x, 1));
    try sameInt(4, try index(f64, .amin, 5, &x, 1));
    try sameInt(3, try index(f64, .max, 5, &x, 1));
    try sameInt(2, try index(f64, .min, 5, &x, 1));
    try sameInt(0, try index(f64, .amax, 5, &x, 0));
    try check.expectClose(f64, 7, try extreme(f64, .amax, 5, &x, 1), 0);
    try check.expectClose(f64, 1, try extreme(f64, .amin, 5, &x, 1), 0);
    try check.expectClose(f64, 7, try extreme(f64, .max, 5, &x, 1), 0);
    try check.expectClose(f64, -7, try extreme(f64, .min, 5, &x, 1), 0);
    try check.expectClose(f64, -7, try extreme(f64, .max, 5, x[1..], 0), 0);
    try check.expectClose(f64, 7, try extreme(f64, .amax, 5, x[1..], 0), 0);
    var z = [_]Complex32{ .init(1, 1), .init(-3, 2), .init(0, 4) };
    try sameInt(2, try index(Complex32, .amax, 3, &z, 1));
    try check.expectClose(f32, 5, try extreme(Complex32, .amax, 3, &z, 1), 0);
}

test "interface: rot, rotm and axpby match a plain loop" {
    inline for (reals) |T| {
        const n = 9;
        var x: [n * 2]T = undefined;
        var y: [n]T = undefined;
        check.fill(T, &x, 6);
        check.fill(T, &y, 7);
        var wx = x;
        var wy = y;
        for (0..n) |i| {
            const a = wx[(n - 1 - i) * 2];
            const b = wy[i];
            wx[(n - 1 - i) * 2] = 0.6 * a + 0.8 * b;
            wy[i] = 0.6 * b - 0.8 * a;
        }
        try rot(T, n, &x, -2, &y, 1, 0.6, 0.8);
        try same(T, &wx, &x);
        try same(T, &wy, &y);

        const param = [5]T{ -1, 2, 3, -1, 0.5 };
        var mx = x;
        var my = y;
        for (0..n) |i| {
            const a = mx[i * 2];
            const b = my[i];
            mx[i * 2] = 2 * a + -1 * b;
            my[i] = 3 * a + 0.5 * b;
        }
        try rotm(T, n, &x, 2, &y, 1, &param);
        try same(T, &mx, &x);
        try same(T, &my, &y);

        if (comptime types.isExtended(T)) continue;
        var ay = y;
        for (0..n) |i| ay[i] = 1.5 * x[i * 2] - 0.5 * y[i];
        try axpby(T, n, 1.5, &x, 2, -0.5, &y, 1);
        try same(T, &ay, &y);
    }
}

test "interface: axpy returns for n <= 0 or zero alpha, and collapses zero strides" {
    var x = [_]f64{ 1, 2 };
    var y = [_]f64{ 5, 6 };
    try axpy(f64, 0, 1, &x, 1, &y, 1);
    try axpy(f64, -3, 1, &x, 1, &y, 1);
    try axpy(f64, 2, 0, &x, 1, &y, 1);
    try same(f64, &.{ 5, 6 }, &y);
    var cx = [_]Complex64{.init(1, 1)};
    var cy = [_]Complex64{.init(5, 6)};
    try axpy(Complex64, 1, .init(0, 0), &cx, 1, &cy, 1);
    try same(Complex64, &.{.init(5, 6)}, &cy);

    var x1 = [_]f64{3};
    var y1 = [_]f64{1};
    try axpy(f64, 4, 2, &x1, 0, &y1, 0);
    try same(f64, &.{25}, &y1);
    var cz = [_]Complex64{.init(0, 0)};
    try axpy(Complex64, 2, .init(2, 3), &cx, 0, &cz, 0);
    // 2 * ((2+3i)(1+i)) = 2 * (-1+5i)
    try same(Complex64, &.{.init(-2, 10)}, &cz);
}

test "interface: axpy with negative strides, and axpyc" {
    var x = [_]f64{ 1, 2, 3 };
    var y = [_]f64{ 0, 0, 0 };
    try axpy(f64, 3, 1, &x, -1, &y, 1);
    try same(f64, &.{ 3, 2, 1 }, &y);
    // both negative: the same pairing as both positive
    var y2 = [_]f64{ 0, 0, 0 };
    try axpy(f64, 3, 1, &x, -1, &y2, -1);
    try same(f64, &x, &y2);
    var cx = [_]Complex32{ .init(1, 0), .init(2, 0), .init(3, 0) };
    var cy = [_]Complex32{ .init(0, 0), .init(0, 0), .init(0, 0) };
    try axpy(Complex32, 3, .init(1, 0), &cx, 1, &cy, -1);
    try same(Complex32, &.{ .init(3, 0), .init(2, 0), .init(1, 0) }, &cy);

    var zx = [_]Complex64{.init(1, 1)};
    var zy = [_]Complex64{.init(0, 0)};
    try axpyc(Complex64, 1, .init(2, 3), &zx, 1, &zy, 1);
    // (2+3i) * conj(1+i) = 5 + i
    try same(Complex64, &.{.init(5, 1)}, &zy);

    var qx = [_]XDouble{ 1, 2, 3 };
    var qy = [_]XDouble{ 1, 1, 1 };
    try axpy(XDouble, 3, 2, &qx, 1, &qy, 1);
    try same(XDouble, &.{ 3, 5, 7 }, &qy);
}

test "interface: rotg zeroes b, and rotmg's param zeroes the second component" {
    inline for (reals) |T| {
        var a: T = 3;
        var b: T = 4;
        var c: T = undefined;
        var s: T = undefined;
        rotg(T, &a, &b, &c, &s);
        try check.expectClose(T, 5, a, check.tolerance(T, 4));
        try check.expectClose(T, 0.6, c, check.tolerance(T, 4));
        try check.expectClose(T, 0.8, s, check.tolerance(T, 4));
        // |a| <= |b| and c != 0: z = 1/c
        try check.expectClose(T, 1 / @as(T, 0.6), b, check.tolerance(T, 8));

        var d1: T = 2;
        var d2: T = 3;
        var x1: T = 1.5;
        const y1: T = 0.5;
        const d1_0 = d1;
        const d2_0 = d2;
        const x1_0 = x1;
        var param: [5]T = undefined;
        rotmg(T, &d1, &d2, &x1, y1, &param);
        // apply H to (x1, y1): the second component vanishes
        var hx = [1]T{x1_0};
        var hy = [1]T{y1};
        try rotm(T, 1, &hx, 1, &hy, 1, &param);
        try check.expectClose(T, 0, hy[0], check.tolerance(T, 8));
        try check.expectClose(T, x1, hx[0], check.tolerance(T, 8));
        // and the weighted length is kept: d1*x1^2 + d2*y1^2
        try check.expectClose(T, d1_0 * x1_0 * x1_0 + d2_0 * y1 * y1, d1 * x1 * x1, check.tolerance(T, 16));
    }
}

test "interface: zrotg zeroes b against a over every branch" {
    inline for (.{ Complex32, Complex64 } ++ check.x_complex) |C| {
        const F = Float(C);
        const eps: f64 = if (F == f32) 1e-5 else 1e-12;
        // magnitudes that take the unscaled and both scaled paths, with a or b zero,
        // b real or imaginary only, and a tiny against a huge b (and the reverse)
        const big: f64 = if (F == f32) 1e30 else 1e300;
        const tiny: f64 = if (F == f32) 1e-30 else 1e-300;
        const cases = [_][4]f64{
            .{ 1, 2, 3, -1 },
            .{ 1 * tiny, 2 * tiny, 3 * tiny, -1 * tiny },
            .{ 1 * big, 2 * big, 3 * big, -1 * big },
            .{ 0, 0, 3, -4 },
            .{ 0, 0, 3 * big, -4 * big },
            .{ 0, 0, 3 * tiny, -4 * tiny },
            .{ 0, 0, 0, 5 },
            .{ 0, 0, -5, 0 },
            .{ 1, 2, 0, 0 },
            .{ tiny, tiny, big, big },
            .{ big, big, tiny, tiny },
            .{ 1, -2, 0, 3 },
            .{ -2, 1, 7, 0 },
        };
        for (cases) |k| {
            const a0: C = .init(@floatCast(k[0]), @floatCast(k[1]));
            const b: C = .init(@floatCast(k[2]), @floatCast(k[3]));
            var a = a0;
            var c: F = undefined;
            var s: C = undefined;
            zrotg(C, &a, b, &c, &s);
            // [c s; -conj(s) c] [a0; b] = [r; 0], |r| = |(a0, b)|, c real, c^2 + |s|^2 = 1
            const norm = @max(@max(@abs(k[0]), @abs(k[1])), @max(@abs(k[2]), @abs(k[3])));
            const top = scale(C, a0, c).add(s.mul(b));
            const bottom = C.init(-s.re, s.im).mul(a0).add(scale(C, b, c));
            try check.expectClose(C, a, top, @floatCast(eps * norm));
            try check.expectClose(C, .init(0, 0), bottom, @floatCast(eps * norm));
            const W = check.Wide(F);
            try check.expectClose(W, 1, @as(W, c) * c + abssq(C, s), eps);
            // |r| against the norm of (a0, b), formed in f64 on scaled parts
            const inv = 1 / norm;
            const want = @sqrt(k[0] * inv * k[0] * inv + k[1] * inv * k[1] * inv + k[2] * inv * k[2] * inv + k[3] * inv * k[3] * inv);
            const got = @sqrt(@as(W, a.re) * inv * a.re * inv + @as(W, a.im) * inv * a.im * inv);
            try check.expectClose(W, want, got, eps);
        }
    }
}

/// Level 1 over `T`, which has no vtable slots, against an `f128` reference:
/// `x` at stride 2 and `y` at stride -1.
fn anyFloat(comptime T: type) !void {
    const W = check.Wide(T);
    const n = 37;
    const ulps = 2 * (n + 2);
    var x: [2 * n]T = undefined;
    var y: [n]T = undefined;
    check.fill(T, &x, 11);
    check.fill(T, &y, 12);
    var wx: [n]W = undefined;
    var wy: [n]W = undefined;
    for (0..n) |i| {
        wx[i] = check.wide(T, x[2 * i]);
        wy[i] = check.wide(T, y[n - 1 - i]);
    }
    const alpha: T = if (comptime isComplex(T)) .init(1.25, -0.5) else 1.25;
    const beta: T = if (comptime isComplex(T)) .init(-0.75, 0.25) else -0.75;
    const wa = check.wide(T, alpha);
    const wb = check.wide(T, beta);

    // reductions
    var d = zero(W);
    var dc = zero(W);
    var dscale: f128 = 0;
    var ssq: f128 = 0;
    var s1: f128 = 0;
    var best: usize = 0;
    for (0..n) |i| {
        d = check.addW(d, check.mulW(wx[i], wy[i]));
        dc = check.addW(dc, check.mulW(check.conjW(wx[i]), wy[i]));
        dscale += check.abs1(wx[i]) * check.abs1(wy[i]);
        ssq += if (comptime isComplex(T)) wx[i].re * wx[i].re + wx[i].im * wx[i].im else wx[i] * wx[i];
        s1 += check.abs1(wx[i]);
        if (check.abs1(wx[i]) > check.abs1(wx[best])) best = i;
    }
    try check.expectNear(T, d, try dot(T, n, &x, 2, &y, -1), ulps, dscale);
    if (comptime isComplex(T)) try check.expectNear(T, dc, try dotc(T, n, &x, 2, &y, -1), ulps, dscale);
    const norm = @sqrt(ssq);
    try check.expectNear(Float(T), norm, try nrm2(T, n, &x, 2), ulps, norm);
    try check.expectNear(Float(T), s1, try asum(T, n, &x, 2), ulps, s1);
    try check.expectNear(Float(T), check.abs1(wx[best]), try extreme(T, .amax, n, &x, 2), 4, s1);
    // the first of near-equal keys may differ; its key must be the largest
    const at = try index(T, .amax, n, &x, 2);
    try check.expectNear(Float(T), check.abs1(wx[best]), @floatCast(check.abs1(wx[@intCast(at - 1)])), 4, check.abs1(wx[best]));

    // updates, element by element
    var y2 = y;
    try axpy(T, n, alpha, &x, 2, &y2, -1);
    for (0..n) |i| {
        const want = check.addW(check.mulW(wa, wx[i]), wy[i]);
        try check.expectNear(T, want, y2[n - 1 - i], 8, check.abs1(wa) * check.abs1(wx[i]) + check.abs1(wy[i]));
    }
    y2 = y;
    try axpby(T, n, alpha, &x, 2, beta, &y2, -1);
    for (0..n) |i| {
        const want = check.addW(check.mulW(wa, wx[i]), check.mulW(wb, wy[i]));
        try check.expectNear(T, want, y2[n - 1 - i], 8, check.abs1(wa) * check.abs1(wx[i]) + check.abs1(wb) * check.abs1(wy[i]));
    }
    var x2 = x;
    try scal(T, n, alpha, &x2, 2);
    for (0..n) |i| try check.expectNear(T, check.mulW(wa, wx[i]), x2[2 * i], 8, check.abs1(wa) * check.abs1(wx[i]));
    const c: Float(T) = 0.6;
    const sn: Float(T) = 0.8;
    x2 = x;
    y2 = y;
    try rot(T, n, &x2, 2, &y2, -1, c, sn);
    const wc = check.wide(Float(T), c);
    const ws = check.wide(Float(T), sn);
    for (0..n) |i| {
        const mag = check.abs1(wx[i]) + check.abs1(wy[i]);
        const rx = if (comptime isComplex(T)) W.init(wc * wx[i].re + ws * wy[i].re, wc * wx[i].im + ws * wy[i].im) else wc * wx[i] + ws * wy[i];
        const ry = if (comptime isComplex(T)) W.init(wc * wy[i].re - ws * wx[i].re, wc * wy[i].im - ws * wx[i].im) else wc * wy[i] - ws * wx[i];
        try check.expectNear(T, rx, x2[2 * i], 8, mag);
        try check.expectNear(T, ry, y2[n - 1 - i], 8, mag);
    }
    x2 = x;
    y2 = y;
    try swap(T, n, &x2, 2, &y2, -1);
    try copy(T, n, &x2, 2, &y2, -1);
    for (0..n) |i| {
        try check.expectNear(T, wy[i], x2[2 * i], 0, 0);
        try check.expectNear(T, wy[i], y2[n - 1 - i], 0, 0);
    }
    if (comptime isComplex(T)) return;

    var ps = [5]T{ -1, 2, 3, -1, 0.5 };
    x2 = x;
    y2 = y;
    try rotm(T, n, &x2, 2, &y2, -1, &ps);
    for (0..n) |i| {
        const mag = 2 * @abs(wx[i]) + @abs(wy[i]);
        try check.expectNear(T, 2 * wx[i] - wy[i], x2[2 * i], 8, mag);
        try check.expectNear(T, 3 * wx[i] + 0.5 * wy[i], y2[n - 1 - i], 8, mag);
    }
    try check.expectNear(T, sumOf(T, &wx), try sum(T, n, &x, 2), ulps, s1);
    var ga: T = 3;
    var gb: T = 4;
    var gc: T = undefined;
    var gs: T = undefined;
    rotg(T, &ga, &gb, &gc, &gs);
    try check.expectNear(T, 5, ga, 4, 5);
    try check.expectNear(T, 0.6, gc, 4, 1);
    try check.expectNear(T, 0.8, gs, 4, 1);
}

fn sumOf(comptime T: type, w: []const check.Wide(T)) f128 {
    var acc: f128 = 0;
    for (w) |v| acc += v;
    return acc;
}

test "interface: level 1 at any float type runs the portable kernels built for it" {
    inline for (check.any_floats) |F| {
        inline for (.{ F, std.math.Complex(F) }) |T| anyFloat(T) catch |e| {
            check.say("{s}\n", .{@typeName(T)});
            return e;
        };
    }
}

test "interface: f16 reductions accumulate wide, and counts past its range stay finite" {
    const H = f16;
    const C = std.math.Complex(H);
    const al = std.testing.allocator;
    // past 2048 an f16 sum of ones stalls; past 65504 it overflows
    const n = 5000;
    const ones = try al.alloc(H, n);
    defer al.free(ones);
    @memset(ones, 1);
    try check.expectClose(H, 5000, try asum(H, n, ones.ptr, 1), 0);
    try check.expectClose(H, 5000, try sum(H, n, ones.ptr, 1), 0);
    try check.expectClose(H, 5000, try dot(H, n, ones.ptr, 1, ones.ptr, 1), 0);
    const big = try al.alloc(H, 4000);
    defer al.free(big);
    for (big, 0..) |*v, i| v.* = if (i < 2000) 40 else -40;
    try check.expectClose(H, 0, try sum(H, 4000, big.ptr, 1), 0);
    const cones = try al.alloc(C, 3000);
    defer al.free(cones);
    @memset(cones, .init(1, 0));
    try check.expectClose(C, .init(3000, 0), try dot(C, 3000, cones.ptr, 1, cones.ptr, 1), 0);
    try check.expectClose(C, .init(3000, 0), try dotc(C, 3000, cones.ptr, 1, cones.ptr, 1), 0);
    // gemv's transposed sums, real and complex
    var y = [_]H{0};
    try @import("../blas2/interface.zig").gemv(al, H, .t, n, 1, 1, ones.ptr, n, ones.ptr, 1, 0, &y, 1);
    try check.expectClose(H, 5000, y[0], 0);
    var cy = [_]C{.init(0, 0)};
    try @import("../blas2/interface.zig").gemv(al, C, .c, 3000, 1, .init(1, 0), cones.ptr, 3000, cones.ptr, 1, .init(0, 0), &cy, 1);
    try check.expectClose(C, .init(3000, 0), cy[0], 0);
    // a zero stride: sqrt(n) |x[0]| with n past f16's range
    var x = [_]H{0.01};
    const want: f64 = @sqrt(70000.0) * @as(f64, @as(H, 0.01));
    try check.expectClose(f64, want, try nrm2(H, 70000, &x, 0), 4e-3);
    var cx = [_]C{.init(0.01, 0)};
    try check.expectClose(f64, want, try nrm2(C, 70000, &cx, 0), 4e-3);
    var yz = [_]H{0};
    try axpy(H, 70000, 1, &x, 0, &yz, 0);
    try check.expectClose(f64, 700, yz[0], 1);
}

test "interface: rotmg in f16 rescales within its range" {
    const H = f16;
    // d1 small and large enough to take each rescaling loop
    for ([_][4]H{ .{ 1.0 / 1024.0, 1, 1, 0.5 }, .{ 1000, 1, 1, 0.5 }, .{ 2, 3, 1.5, 0.5 }, .{ 1, 1.0 / 2048.0, 0.25, 4 } }) |c| {
        var d1 = c[0];
        var d2 = c[1];
        var x1 = c[2];
        const y1 = c[3];
        var param: [5]H = @splat(0);
        rotmg(H, &d1, &d2, &x1, y1, &param);
        for ([_]H{ d1, d2, x1 }) |v| try std.testing.expect(std.math.isFinite(v) and v != 0);
        var hx = [1]H{c[2]};
        var hy = [1]H{y1};
        try rotm(H, 1, &hx, 1, &hy, 1, &param);
        try check.expectClose(f64, 0, hy[0], 1e-2 * @abs(@as(f64, c[2])));
        try check.expectClose(f64, x1, hx[0], 4e-3 * @abs(@as(f64, x1)));
        // the weighted length d1*x1^2 + d2*y1^2 is kept
        const was = @as(f64, c[0]) * c[2] * c[2] + @as(f64, c[1]) * y1 * y1;
        try check.expectClose(f64, was, @as(f64, d1) * x1 * x1, 1e-2 * was);
    }
}

const bf = @import("../ext/bf16/__root__.zig");

fn b16(v: f32) BF16 {
    return bf.toBf16(v);
}

test "interface: sbdot, the conversions and bscal do nothing for n <= 0" {
    var x = [_]BF16{ b16(1), b16(2) };
    var f = [_]f32{ 7, 8 };
    var d = [_]f64{ 7, 8 };
    for ([_]isize{ 0, -3 }) |n| {
        try std.testing.expectEqual(@as(f32, 0), try sbdot(n, &x, 1, &x, 1));
        try sbtobf16(f32, n, &f, 1, &x, 1);
        try sbtobf16(f64, n, &d, 1, &x, 1);
        try bf16to(f32, n, &x, 1, &f, 1);
        try bf16to(f64, n, &x, 1, &d, 1);
        try bscal(n, b16(0), &x, 1);
    }
    try std.testing.expectEqual(@as(f32, 7), f[0]);
    try std.testing.expectEqual(@as(f64, 8), d[1]);
    try std.testing.expectEqual(@as(f32, 2), bf.fromBf16(x[1]));
}

test "interface: a negative stride reads and writes bfloat16 vectors from their ends" {
    const n = 5;
    var f: [n * 2]f32 = undefined;
    check.fill(f32, &f, 21);
    var h: [n * 3]BF16 = @splat(b16(0));
    // stride -3 puts element 0 last: the last element of the buffer region
    try sbtobf16(f32, n, &f, 2, &h, -3);
    for (0..n) |i| {
        try std.testing.expectEqual(@backingInt(b16(f[i * 2])), @backingInt(h[(n - 1 - i) * 3]));
    }
    var back: [n]f64 = undefined;
    try bf16to(f64, n, &h, -3, &back, 1);
    for (0..n) |i| try std.testing.expectEqual(@as(f64, bf.fromBf16(b16(f[i * 2]))), back[i]);
    // a reversed read of the same storage dots against itself forwards
    var fwd: [n]BF16 = undefined;
    for (0..n) |i| fwd[i] = h[(n - 1 - i) * 3];
    try std.testing.expectEqual(try sbdot(n, &fwd, 1, &fwd, 1), try sbdot(n, &h, -3, &h, -3));
    var want: f64 = 0;
    for (0..n) |i| want += @as(f64, bf.fromBf16(fwd[i])) * bf.fromBf16(fwd[i]);
    try check.expectClose(f32, @floatCast(want), try sbdot(n, &fwd, 1, &fwd, 1), check.tolerance(f32, 4 * n));
}

test "interface: bscal skips alpha 1 and a non-positive stride, and scales the rest" {
    var x = [_]BF16{ b16(1.5), b16(-2), b16(4) };
    try bscal(3, b16(1), &x, 1);
    try bscal(3, b16(2), &x, 0);
    try bscal(3, b16(2), &x, -1);
    try std.testing.expectEqual(@as(f32, -2), bf.fromBf16(x[1]));
    try bscal(3, b16(2), &x, 1);
    try std.testing.expectEqual(@as(f32, 3), bf.fromBf16(x[0]));
    try std.testing.expectEqual(@as(f32, -4), bf.fromBf16(x[1]));
    try std.testing.expectEqual(@as(f32, 8), bf.fromBf16(x[2]));
    // flag 1: a zero alpha zeroes the finite elements and keeps NaN
    x[2] = b16(std.math.nan(f32));
    try bscal(3, b16(0), &x, 1);
    try std.testing.expectEqual(@as(f32, 0), bf.fromBf16(x[0]));
    try std.testing.expect(std.math.isNan(bf.fromBf16(x[2])));
}

test "interface: rotg handles scaled and special inputs" {
    inline for (reals) |T| {
        const big: T = std.math.floatMax(T) / 4;
        const tiny: T = std.math.floatMin(T) * 4;
        for ([_][2]T{ .{ big, big }, .{ tiny, tiny }, .{ big, tiny }, .{ tiny, big }, .{ -3, 4 }, .{ 3, -4 } }) |ab| {
            var a = ab[0];
            var b = ab[1];
            var c: T = undefined;
            var s: T = undefined;
            rotg(T, &a, &b, &c, &s);
            try std.testing.expect(std.math.isFinite(a));
            const W = check.Wide(T);
            try check.expectClose(W, 1, @as(W, c) * c + @as(W, s) * s, 1e-6);
            // the rotation zeroes b: -s*a0 + c*b0 = 0, relative to the norm
            const norm = @max(@abs(ab[0]), @abs(ab[1]));
            const resid = -@as(W, s) * ab[0] + @as(W, c) * ab[1];
            try std.testing.expect(@abs(resid) <= 1e-6 * norm);
        }
        var a: T = 0;
        var b: T = 0;
        var c: T = undefined;
        var s: T = undefined;
        rotg(T, &a, &b, &c, &s);
        try std.testing.expectEqual(@as(T, 1), c);
        try std.testing.expectEqual(@as(T, 0), s);
        a = std.math.nan(T);
        b = 1;
        rotg(T, &a, &b, &c, &s);
        try std.testing.expect(std.math.isNan(a));
    }
}

test "portable kernels: special values and ranges" {
    const nrm2k = @import("nrm2/__root__.zig");
    const scalk = @import("scal/__root__.zig");
    const iamaxk = @import("iamax/__root__.zig");
    const dotk = @import("dot/__root__.zig");
    // nrm2 neither overflows nor underflows
    inline for (.{ .{ f32, nrm2k.snrm2_generic }, .{ f64, nrm2k.dnrm2_generic } }) |p| {
        const T = p[0];
        const big: T = std.math.floatMax(T) / 2;
        const tiny: T = std.math.floatMin(T);
        var x = [_]T{ big, big, big };
        try check.expectClose(f64, @as(f64, big) * @sqrt(3.0), p[1](3, &x, 1), @as(f64, big) * 1e-6);
        x = .{ tiny, tiny, tiny };
        try check.expectClose(f64, @as(f64, tiny) * @sqrt(3.0), p[1](3, &x, 1), @as(f64, tiny) * 1e-6);
        x = .{ 0, 0, 0 };
        try std.testing.expectEqual(@as(T, 0), p[1](3, &x, 1));
    }
    // scal with a zero alpha: the flag makes NaN and Inf stay NaN
    var v = [_]f64{ 1, std.math.inf(f64), std.math.nan(f64), -2 };
    _ = scalk.dscal_generic(4, 0, 0, 0, &v, 1, &v, 1, &v, 1);
    try std.testing.expectEqual(@as(f64, 0), v[0]);
    try std.testing.expect(std.math.isNan(v[1]) and std.math.isNan(v[2]));
    try std.testing.expectEqual(@as(f64, 0), v[3]);
    v = .{ 1, std.math.inf(f64), std.math.nan(f64), -2 };
    _ = scalk.dscal_generic(4, 0, 0, 0, &v, 1, &v, 1, &v, 0);
    for (v) |e| try std.testing.expectEqual(@as(f64, 0), e);
    // iamax: first of equal maxima, 1-based; a NaN first element is never replaced
    var ties = [_]f64{ 1, -3, 3, 2 };
    try std.testing.expectEqual(@as(BlasLong, 2), iamaxk.idamax_generic(4, &ties, 1));
    var nan_first = [_]f64{ std.math.nan(f64), 5, 7 };
    try std.testing.expectEqual(@as(BlasLong, 1), iamaxk.idamax_generic(3, &nan_first, 1));
    // dsdot accumulates in double: a sum f32 would lose
    var xs: [3]f32 = .{ 16777216, 1, 1 };
    try std.testing.expectEqual(@as(f64, 16777216.0 * 16777216.0 + 2), dotk.dsdot_generic(3, &xs, 1, &xs, 1));
}
