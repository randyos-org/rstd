//! The level-1 interface: the checks and pointer set-up OpenBLAS's
//! `interface/*.c` do before calling a kernel. Each takes the element type `T`
//! (`f32`, `f64`, `XDouble`, or their `Complex`) explicitly; `n` and strides
//! are `BlasInt`, as in the Fortran interface.
//!
//! A negative stride names the vector's last element first: the interface
//! moves the pointer to it, and the kernel walks back from there. Kernels
//! assume `n >= 1`; every routine here returns before calling one otherwise.

const std = @import("std");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");
const check = @import("check.zig");

const BlasInt = types.BlasInt;
const BlasLong = types.BlasLong;
const BF16 = types.BF16;
const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

fn isComplex(comptime T: type) bool {
    return T == Complex32 or T == Complex64 or T == ComplexXD;
}

/// The float type of `T`'s parts.
fn Float(comptime T: type) type {
    return if (isComplex(T)) @FieldType(T, "re") else T;
}

/// The reference's precision letter for `T` (`f64` before `XDouble`, which
/// may be the same type).
fn letter(comptime T: type) []const u8 {
    if (T == f32 or T == Complex32) return if (T == f32) "s" else "c";
    if (T == f64 or T == Complex64) return if (T == f64) "d" else "z";
    if (T == XDouble or T == ComplexXD) return if (T == XDouble) "q" else "x";
    @compileError("not a BLAS element type: " ++ @typeName(T));
}

/// The kernel for `T` in the slot named `prefix ++ letter(T) ++ base ++ "_k"`.
fn slot(comptime T: type, comptime prefix: []const u8, comptime base: []const u8) @typeInfo(@FieldType(Vtable, prefix ++ letter(T) ++ base ++ "_k")).optional.child {
    return blas.kernel(@field(blas.Slot, prefix ++ letter(T) ++ base ++ "_k"));
}

/// `T`'s elements as the kernels take them: floats, two per complex.
fn raw(comptime T: type, p: [*]T) [*]Float(T) {
    return @ptrCast(p);
}

/// Where a vector of `n` with stride `inc` starts: its last element when
/// `inc < 0`.
fn first(comptime T: type, p: [*]T, n: BlasLong, inc: BlasLong) [*]T {
    return if (inc < 0) @ptrCast(util.at(T, p, -(n - 1) * inc)) else p;
}

fn long(v: BlasInt) BlasLong {
    return @intCast(v);
}

// ---- copy, swap, scal --------------------------------------------------------

/// `?COPY` (`interface/copy.c`): `y := x`.
pub fn copy(comptime T: type, n: BlasInt, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt) void {
    if (n <= 0) return;
    const kn = long(n);
    _ = slot(T, "", "copy")(kn, raw(T, first(T, x, kn, incx)), incx, raw(T, first(T, y, kn, incy)), incy);
}

/// `?SWAP` (`interface/swap.c`, `zswap.c`): `x <-> y`.
pub fn swap(comptime T: type, n: BlasInt, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt) void {
    if (n <= 0) return;
    const kn = long(n);
    const px = raw(T, first(T, x, kn, incx));
    const py = raw(T, first(T, y, kn, incy));
    if (comptime isComplex(T)) {
        _ = slot(T, "", "swap")(kn, 0, 0, 0, 0, px, incx, py, incy, py, 0);
    } else {
        _ = slot(T, "", "swap")(kn, 0, 0, 0, px, incx, py, incy, py, 0);
    }
}

/// `?SCAL` (`interface/scal.c`, `zscal.c`): `x := alpha*x`. Nothing for
/// `incx <= 0`, and nothing to do for `alpha == 1`. With the flag the
/// interface passes (1), a zero `alpha` leaves NaN and infinity as NaN, as
/// `alpha*x` would, and finite elements 0.
pub fn scal(comptime T: type, n: BlasInt, alpha: T, x: [*]T, incx: BlasInt) void {
    if (n <= 0 or incx <= 0) return;
    const px = raw(T, x);
    if (comptime isComplex(T)) {
        if (alpha.re == 1 and alpha.im == 0) return;
        _ = slot(T, "", "scal")(long(n), 0, 0, alpha.re, alpha.im, px, incx, px, 0, px, 1);
    } else {
        if (alpha == 1) return;
        _ = slot(T, "", "scal")(long(n), 0, 0, alpha, px, incx, px, 0, px, 1);
    }
}

/// `CSSCAL`/`ZDSCAL`: complex `x` scaled by a real `alpha`.
pub fn scalReal(comptime T: type, n: BlasInt, alpha: Float(T), x: [*]T, incx: BlasInt) void {
    if (comptime !isComplex(T)) @compileError("scalReal: " ++ @typeName(T) ++ " is not complex");
    scal(T, n, .init(alpha, 0), x, incx);
}

// ---- dot ---------------------------------------------------------------------

/// `?DOT`/`?DOTU` (`interface/dot.c`, `zdot.c`): `x . y`, 0 for `n <= 0`.
pub fn dot(comptime T: type, n: BlasInt, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt) T {
    return dotImpl(T, false, n, x, incx, y, incy);
}

/// `?DOTC`: `conj(x) . y`. Complex only.
pub fn dotc(comptime T: type, n: BlasInt, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt) T {
    if (comptime !isComplex(T)) @compileError("dotc: " ++ @typeName(T) ++ " is not complex");
    return dotImpl(T, true, n, x, incx, y, incy);
}

fn dotImpl(comptime T: type, comptime conj: bool, n: BlasInt, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt) T {
    if (n <= 0) return zero(T);
    const kn = long(n);
    const px = raw(T, first(T, x, kn, incx));
    const py = raw(T, first(T, y, kn, incy));
    if (comptime isComplex(T)) {
        return slot(T, "", if (conj) "dotc" else "dotu")(kn, px, incx, py, incy).toComplex();
    }
    return slot(T, "", "dot")(kn, px, incx, py, incy);
}

/// `DSDOT` (`interface/dsdot.c`): the dot of two `f32` vectors, summed and
/// returned in `f64`.
pub fn dsdot(n: BlasInt, x: [*]f32, incx: BlasInt, y: [*]f32, incy: BlasInt) f64 {
    if (n <= 0) return 0;
    const kn = long(n);
    return blas.kernel(.dsdot_k)(kn, first(f32, x, kn, incx), incx, first(f32, y, kn, incy), incy);
}

/// `SDSDOT` (`interface/sdsdot.c`): `alpha + x . y`, summed in `f64` and
/// returned as `f32`; `alpha` alone for `n <= 0`.
pub fn sdsdot(n: BlasInt, alpha: f32, x: [*]f32, incx: BlasInt, y: [*]f32, incy: BlasInt) f32 {
    if (n <= 0) return alpha;
    return @floatCast(dsdot(n, x, incx, y, incy) + alpha);
}

// ---- norms and sums ----------------------------------------------------------

/// `?NRM2`/`SCNRM2`/`DZNRM2` (`interface/nrm2.c`): the Euclidean norm. For
/// `incx == 0` it is `sqrt(n) * |x[0]|`, worked out here as the reference
/// does; a negative stride is the same vector backwards.
pub fn nrm2(comptime T: type, n: BlasInt, x: [*]T, incx: BlasInt) Float(T) {
    const F = Float(T);
    if (n <= 0) return 0;
    const fn_: F = @floatFromInt(n);
    if (comptime !isComplex(T)) {
        if (n == 1) return @abs(x[0]);
        if (incx == 0) return @sqrt(fn_) * @abs(x[0]);
    } else if (incx == 0) {
        const fr = @abs(x[0].re);
        const fi = @abs(x[0].im);
        const lo = @min(fr, fi);
        const hi = @max(fr, fi);
        if (hi == 0) return hi;
        if (hi == lo) return @sqrt(fn_) * @sqrt(@as(F, 2)) * hi;
        return @sqrt(fn_) * hi * @sqrt(1 + (lo / hi) * (lo / hi));
    }
    const kn = long(n);
    return slot(T, "", "nrm2")(kn, raw(T, first(T, x, kn, incx)), incx);
}

/// `?ASUM`/`SCASUM`/`DZASUM` (`interface/asum.c`): the sum of `|x|`, of
/// `|re| + |im|` for complex. 0 for `n <= 0`; the kernel gives 0 for
/// `incx <= 0`.
pub fn asum(comptime T: type, n: BlasInt, x: [*]T, incx: BlasInt) Float(T) {
    if (n <= 0) return 0;
    return slot(T, "", "asum")(long(n), raw(T, x), incx);
}

/// `?SUM` (`interface/sum.c`, an OpenBLAS extension): the sum of `x`, of
/// `re + im` for complex.
pub fn sum(comptime T: type, n: BlasInt, x: [*]T, incx: BlasInt) Float(T) {
    if (n <= 0) return 0;
    return slot(T, "", "sum")(long(n), raw(T, x), incx);
}

// ---- extremes ----------------------------------------------------------------

/// Which extreme: by magnitude (`|re| + |im|` for complex) or by value.
pub const Extreme = enum { amax, amin, max, min };

/// `I?AMAX`, `I?AMIN`, `I?MAX`, `I?MIN` (`interface/imax.c`): the 1-based
/// index of the first extreme element; 0 for `n <= 0` (and, from the
/// kernel, for `incx <= 0`). `max` and `min` are real only.
pub fn index(comptime T: type, comptime which: Extreme, n: BlasInt, x: [*]T, incx: BlasInt) BlasInt {
    if (n <= 0) return 0;
    const ret: BlasLong = slot(T, "i", @tagName(which))(long(n), raw(T, x), incx);
    return @intCast(@min(ret, long(n)));
}

/// `?AMAX`, `?AMIN`, `?MAX`, `?MIN` (`interface/max.c`, OpenBLAS extensions):
/// the extreme value itself, `|x|` for `amax`/`amin`. 0 for `n <= 0`; for
/// `incx == 0`, `x[0]` (by magnitude for `amax`/`amin`).
pub fn extreme(comptime T: type, comptime which: Extreme, n: BlasInt, x: [*]T, incx: BlasInt) Float(T) {
    if (n <= 0) return 0;
    if (incx == 0) {
        if (comptime isComplex(T)) return @abs(x[0].re) + @abs(x[0].im);
        return if (which == .amax or which == .amin) @abs(x[0]) else x[0];
    }
    return slot(T, "", @tagName(which))(long(n), raw(T, x), incx);
}

// ---- rotations ---------------------------------------------------------------

/// `?ROT`, `CSROT`/`ZDROT` (`interface/rot.c`, `zrot.c`): the plane rotation
/// `x := c*x + s*y`, `y := c*y - s*x`, with real `c` and `s`.
pub fn rot(comptime T: type, n: BlasInt, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt, c: Float(T), s: Float(T)) void {
    if (n <= 0) return;
    const kn = long(n);
    const px = raw(T, first(T, x, kn, incx));
    const py = raw(T, first(T, y, kn, incy));
    const kern = if (comptime T == Complex32)
        blas.kernel(.csrot_k)
    else if (comptime T == Complex64)
        blas.kernel(.zdrot_k)
    else if (comptime T == ComplexXD)
        blas.kernel(.xqrot_k)
    else
        slot(T, "", "rot");
    _ = kern(kn, px, incx, py, incy, c, s);
}

/// `?ROTM` (`interface/rotm.c`): the modified rotation `param` describes
/// (`param[0]` the flag). The reference passes everything to the kernel,
/// which handles `n <= 0` and negative strides itself.
pub fn rotm(comptime T: type, n: BlasInt, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt, param: *const [5]T) void {
    if (comptime isComplex(T)) @compileError("rotm: " ++ @typeName(T) ++ " is not real");
    _ = slot(T, "", "rotm")(long(n), x, incx, y, incy, @constCast(param));
}

// ---- axpy, axpby -------------------------------------------------------------

/// `?AXPY` (`interface/axpy.c`, `zaxpy.c`): `y := alpha*x + y`. Nothing to
/// do for `alpha == 0`; both strides zero is `y[0] += n*alpha*x[0]`.
pub fn axpy(comptime T: type, n: BlasInt, alpha: T, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt) void {
    axpyImpl(T, false, n, alpha, x, incx, y, incy);
}

/// `?AXPYC`: `y := alpha*conj(x) + y`. Complex only.
pub fn axpyc(comptime T: type, n: BlasInt, alpha: T, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt) void {
    if (comptime !isComplex(T)) @compileError("axpyc: " ++ @typeName(T) ++ " is not complex");
    axpyImpl(T, true, n, alpha, x, incx, y, incy);
}

fn axpyImpl(comptime T: type, comptime conj: bool, n: BlasInt, alpha: T, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt) void {
    if (n <= 0) return;
    const kn = long(n);
    if (comptime isComplex(T)) {
        if (alpha.re == 0 and alpha.im == 0) return;
        if (incx == 0 and incy == 0) {
            // the reference does this unconjugated for `axpyc` too
            const cn: Float(T) = @floatFromInt(kn);
            y[0].re += cn * (alpha.re * x[0].re - alpha.im * x[0].im);
            y[0].im += cn * (alpha.im * x[0].re + alpha.re * x[0].im);
            return;
        }
        const px = raw(T, first(T, x, kn, incx));
        const py = raw(T, first(T, y, kn, incy));
        _ = slot(T, "", if (conj) "axpyc" else "axpy")(kn, 0, 0, alpha.re, alpha.im, px, incx, py, incy, py, 0);
    } else {
        if (alpha == 0) return;
        if (incx == 0 and incy == 0) {
            y[0] += @as(T, @floatFromInt(kn)) * alpha * x[0];
            return;
        }
        const px = first(T, x, kn, incx);
        const py = first(T, y, kn, incy);
        _ = slot(T, "", "axpy")(kn, 0, 0, alpha, px, incx, py, incy, py, 0);
    }
}

/// `?AXPBY` (`interface/axpby.c`, `zaxpby.c`, an OpenBLAS extension):
/// `y := alpha*x + beta*y`.
pub fn axpby(comptime T: type, n: BlasInt, alpha: T, x: [*]T, incx: BlasInt, beta: T, y: [*]T, incy: BlasInt) void {
    if (n <= 0) return;
    const kn = long(n);
    const px = raw(T, first(T, x, kn, incx));
    const py = raw(T, first(T, y, kn, incy));
    if (comptime isComplex(T)) {
        _ = slot(T, "", "axpby")(kn, alpha.re, alpha.im, px, incx, beta.re, beta.im, py, incy);
    } else {
        _ = slot(T, "", "axpby")(kn, alpha, px, incx, beta, py, incy);
    }
}

// ---- rotation set-up ---------------------------------------------------------

/// `?ROTG` (`interface/rotg.c`): the rotation `[c s; -s c]` that zeroes `b`
/// against `a`. `a` becomes `r`, `b` the reconstruction value `z`.
pub fn rotg(comptime T: type, a: *T, b: *T, c: *T, s: *T) void {
    if (comptime isComplex(T)) @compileError("rotg: use zrotg for " ++ @typeName(T));
    const safmin = std.math.floatMin(T);
    const ada = @abs(a.*);
    const adb = @abs(b.*);
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
    const sigma: T = std.math.copysign(@as(T, 1), if (ada > adb) a.* else b.*);
    const aa = a.* / sc;
    const bb = b.* / sc;
    const r = sigma * sc * @sqrt(aa * aa + bb * bb);
    c.* = a.* / r;
    s.* = b.* / r;
    var z: T = 1;
    if (ada > adb) z = s.*;
    if (ada <= adb and c.* != 0) z = 1 / c.*;
    a.* = r;
    b.* = z;
}

/// `CROTG`/`ZROTG`: the rotation with real `c` and complex `s` that zeroes
/// `b` against `a`, `[c s; -conj(s) c] [a; b] = [r; 0]`; `a` becomes `r`.
/// Follows LAPACK's `zlartg`, scaling against overflow and underflow.
/// `interface/zrotg.c` transcribes the same algorithm but slips in the
/// scaled branches (a sign in `s`, square roots for squared norms, `|a|^2`
/// for `a`), so this ports `zlartg` itself.
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
    const gam: T = 4096;
    const gamsq: T = 16777216;
    const rgamsq: T = 5.9604645e-8;
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

const reals = .{ f32, f64 };
const all = .{ f32, f64, Complex32, Complex64 };

fn zero(comptime T: type) T {
    return if (comptime isComplex(T)) .init(0, 0) else 0;
}

fn same(comptime T: type, want: []const T, got: []const T) !void {
    try check.expectCloseSlices(T, want, got, check.tolerance(Float(T), 4));
}

fn sameInt(want: BlasInt, got: BlasInt) !void {
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
        for ([_]BlasInt{ 0, -1 }) |n| {
            copy(T, n, &x, 1, &y, 1);
            swap(T, n, &x, 1, &y, 1);
            scal(T, n, zero(T), &x, 1);
            axpby(T, n, zero(T), &x, 1, zero(T), &y, 1);
            rot(T, n, &x, 1, &y, 1, 0, 1);
            try same(T, &x0, &x);
            try same(T, &y0, &y);
            try check.expectClose(T, zero(T), dot(T, n, &x, 1, &y, 1), 0);
            try check.expectClose(Float(T), 0, nrm2(T, n, &x, 1), 0);
            try check.expectClose(Float(T), 0, asum(T, n, &x, 1), 0);
            try check.expectClose(Float(T), 0, sum(T, n, &x, 1), 0);
            try check.expectClose(Float(T), 0, extreme(T, .amax, n, &x, 1), 0);
            try sameInt(0, index(T, .amax, n, &x, 1));
        }
    }
    try check.expectClose(f32, 2.5, sdsdot(0, 2.5, undefined, 1, undefined, 1), 0);
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
        copy(T, n, &x, -2, &y1, 3);
        copy(T, n, &xr, 1, &y2, 3);
        try same(T, &y2, &y1);
        try check.expectClose(T, dot(T, n, &xr, 1, &y, 3), dot(T, n, &x, -2, &y, 3), check.tolerance(Float(T), 4 * n));
        try check.expectClose(Float(T), nrm2(T, n, &xr, 1), nrm2(T, n, &x, -2), check.tolerance(Float(T), 4 * n));
    }
}

test "interface: scal does nothing for incx <= 0 or alpha 1; a zero alpha keeps NaN" {
    inline for (all) |T| {
        var x: [4]T = undefined;
        check.fill(T, &x, 5);
        const x0 = x;
        const two: T = if (comptime isComplex(T)) .init(2, 0) else 2;
        const one: T = if (comptime isComplex(T)) .init(1, 0) else 1;
        scal(T, 4, two, &x, 0);
        scal(T, 4, two, &x, -1);
        scal(T, 4, one, &x, 1);
        try same(T, &x0, &x);
        const nan = std.math.nan(Float(T));
        x[1] = if (comptime isComplex(T)) .init(nan, nan) else nan;
        scal(T, 4, zero(T), &x, 1);
        const at1 = if (comptime isComplex(T)) x[1].re else x[1];
        if (!std.math.isNan(at1)) return error.TestExpectedNan;
        x[1] = zero(T);
        try same(T, &@as([4]T, @splat(zero(T))), &x);
    }
    var z = [_]Complex64{ .init(1, 2), .init(-3, 4) };
    scalReal(Complex64, 2, 2, &z, 1);
    try same(Complex64, &.{ .init(2, 4), .init(-6, 8) }, &z);
}

test "interface: dotc conjugates x, and dsdot/sdsdot sum in f64" {
    var x = [_]Complex32{ .init(1, 2), .init(3, -1) };
    var y = [_]Complex32{ .init(2, 1), .init(-1, 1) };
    // conj(1+2i)(2+i) + conj(3-i)(-1+i) = (4-3i) + (-4+2i)
    try check.expectClose(Complex32, .init(0, -1), dotc(Complex32, 2, &x, 1, &y, 1), 1e-6);
    // (1+2i)(2+i) + (3-i)(-1+i) = 5i + (-2+4i)
    try check.expectClose(Complex32, .init(-2, 9), dot(Complex32, 2, &x, 1, &y, 1), 1e-6);
    // 1e8 + 1 - 1e8 is 1 in f64 but 0 in f32
    var a = [_]f32{ 1e4, 1, -1e4 };
    var b = [_]f32{ 1e4, 1, 1e4 };
    try check.expectClose(f64, 1, dsdot(3, &a, 1, &b, 1), 0);
    try check.expectClose(f32, 1.5, sdsdot(3, 0.5, &a, 1, &b, 1), 0);
}

test "interface: nrm2 of a zero stride is sqrt(n) * |x[0]|" {
    var x = [_]f64{-3};
    try check.expectClose(f64, 3 * @sqrt(4.0), nrm2(f64, 4, &x, 0), 1e-12);
    var z = [_]Complex64{.init(3, -4)};
    try check.expectClose(f64, 5 * @sqrt(4.0), nrm2(Complex64, 4, &z, 0), 1e-12);
    var w = [_]Complex64{.init(2, -2)};
    try check.expectClose(f64, @sqrt(9.0) * @sqrt(2.0) * 2, nrm2(Complex64, 9, &w, 0), 1e-12);
}

test "interface: extremes, by magnitude and by value" {
    var x = [_]f64{ 3, -7, 7, 1, -2 };
    try sameInt(2, index(f64, .amax, 5, &x, 1));
    try sameInt(4, index(f64, .amin, 5, &x, 1));
    try sameInt(3, index(f64, .max, 5, &x, 1));
    try sameInt(2, index(f64, .min, 5, &x, 1));
    try sameInt(0, index(f64, .amax, 5, &x, 0));
    try check.expectClose(f64, 7, extreme(f64, .amax, 5, &x, 1), 0);
    try check.expectClose(f64, 1, extreme(f64, .amin, 5, &x, 1), 0);
    try check.expectClose(f64, 7, extreme(f64, .max, 5, &x, 1), 0);
    try check.expectClose(f64, -7, extreme(f64, .min, 5, &x, 1), 0);
    try check.expectClose(f64, -7, extreme(f64, .max, 5, x[1..], 0), 0);
    try check.expectClose(f64, 7, extreme(f64, .amax, 5, x[1..], 0), 0);
    var z = [_]Complex32{ .init(1, 1), .init(-3, 2), .init(0, 4) };
    try sameInt(2, index(Complex32, .amax, 3, &z, 1));
    try check.expectClose(f32, 5, extreme(Complex32, .amax, 3, &z, 1), 0);
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
        rot(T, n, &x, -2, &y, 1, 0.6, 0.8);
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
        rotm(T, n, &x, 2, &y, 1, &param);
        try same(T, &mx, &x);
        try same(T, &my, &y);

        var ay = y;
        for (0..n) |i| ay[i] = 1.5 * x[i * 2] - 0.5 * y[i];
        axpby(T, n, 1.5, &x, 2, -0.5, &y, 1);
        try same(T, &ay, &y);
    }
}

test "interface: axpy returns for n <= 0 or zero alpha, and collapses zero strides" {
    var x = [_]f64{ 1, 2 };
    var y = [_]f64{ 5, 6 };
    axpy(f64, 0, 1, &x, 1, &y, 1);
    axpy(f64, -3, 1, &x, 1, &y, 1);
    axpy(f64, 2, 0, &x, 1, &y, 1);
    try same(f64, &.{ 5, 6 }, &y);
    var cx = [_]Complex64{.init(1, 1)};
    var cy = [_]Complex64{.init(5, 6)};
    axpy(Complex64, 1, .init(0, 0), &cx, 1, &cy, 1);
    try same(Complex64, &.{.init(5, 6)}, &cy);

    var x1 = [_]f64{3};
    var y1 = [_]f64{1};
    axpy(f64, 4, 2, &x1, 0, &y1, 0);
    try same(f64, &.{25}, &y1);
    var cz = [_]Complex64{.init(0, 0)};
    axpy(Complex64, 2, .init(2, 3), &cx, 0, &cz, 0);
    // 2 * ((2+3i)(1+i)) = 2 * (-1+5i)
    try same(Complex64, &.{.init(-2, 10)}, &cz);
}

test "interface: axpy with negative strides, and axpyc" {
    var x = [_]f64{ 1, 2, 3 };
    var y = [_]f64{ 0, 0, 0 };
    axpy(f64, 3, 1, &x, -1, &y, 1);
    try same(f64, &.{ 3, 2, 1 }, &y);
    // both negative: the same pairing as both positive
    var y2 = [_]f64{ 0, 0, 0 };
    axpy(f64, 3, 1, &x, -1, &y2, -1);
    try same(f64, &x, &y2);
    var cx = [_]Complex32{ .init(1, 0), .init(2, 0), .init(3, 0) };
    var cy = [_]Complex32{ .init(0, 0), .init(0, 0), .init(0, 0) };
    axpy(Complex32, 3, .init(1, 0), &cx, 1, &cy, -1);
    try same(Complex32, &.{ .init(3, 0), .init(2, 0), .init(1, 0) }, &cy);

    var zx = [_]Complex64{.init(1, 1)};
    var zy = [_]Complex64{.init(0, 0)};
    axpyc(Complex64, 1, .init(2, 3), &zx, 1, &zy, 1);
    // (2+3i) * conj(1+i) = 5 + i
    try same(Complex64, &.{.init(5, 1)}, &zy);

    var qx = [_]XDouble{ 1, 2, 3 };
    var qy = [_]XDouble{ 1, 1, 1 };
    axpy(XDouble, 3, 2, &qx, 1, &qy, 1);
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
        rotm(T, 1, &hx, 1, &hy, 1, &param);
        try check.expectClose(T, 0, hy[0], check.tolerance(T, 8));
        try check.expectClose(T, x1, hx[0], check.tolerance(T, 8));
        // and the weighted length is kept: d1*x1^2 + d2*y1^2
        try check.expectClose(T, d1_0 * x1_0 * x1_0 + d2_0 * y1 * y1, d1 * x1 * x1, check.tolerance(T, 16));
    }
}

test "interface: zrotg zeroes b against a, unscaled and scaled" {
    const C = Complex64;
    // ordinary, then tiny and huge enough to take the scaled branches
    for ([_]f64{ 1, 1e-300, 1e300 }) |mag| {
        const a0: C = .init(1 * mag, 2 * mag);
        const b: C = .init(3 * mag, -1 * mag);
        var a = a0;
        var c: f64 = undefined;
        var s: C = undefined;
        zrotg(C, &a, b, &c, &s);
        // [c s; -conj(s) c] [a0; b] = [r; 0], with |r| = |(a0, b)|
        const top = scale(C, a0, c).add(s.mul(b));
        const bottom = C.init(-s.re, s.im).mul(a0).add(scale(C, b, c));
        try check.expectClose(C, a, top, 1e-12 * mag);
        try check.expectClose(C, .init(0, 0), bottom, 1e-12 * mag);
        try check.expectClose(f64, @sqrt(15.0), @sqrt(abssq(C, scale(C, a, 1 / mag))), 1e-12);
        try check.expectClose(f64, 1, c * c + abssq(C, s), 1e-12);
    }
}
