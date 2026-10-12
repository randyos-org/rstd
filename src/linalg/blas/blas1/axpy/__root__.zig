//! `AXPY`: `y := alpha*x + y` for `f32`, `f64`, `XDouble` and their complex
//! counterparts, plus the conjugated complex form (`CAXPYC`/`ZAXPYC`/`XAXPYC`).
//!
//! Layers, top down (mirroring `interface/{,z}axpy.c` -> `gotoblas_t` -> `kernel/`):
//!
//!   * `axpy`/`taxpy` (and `axpyc`/`taxpyc`): the BLAS interface. Infers or
//!     takes the element type, applies the interface's early-outs and
//!     negative-stride pointer adjustment.
//!   * `saxpy`/`daxpy`/`qaxpy`/`caxpy`/`zaxpy`/`xaxpy` (+ `*c`): typed kernel
//!     wrappers over the vtable, no interface semantics.
//!   * `blas.kernel(.saxpy_k)` and the rest: the vtable slots, in the
//!     reference's kernel ABI (padding `dummy*` arguments, interleaved
//!     complex, `i32` result).
//!   * `*_generic`: the portable kernels, which the slots hold until `init`
//!     and wherever a core has none of its own. A core's own are in `<arch>/`,
//!     installed by `kernels/<arch>.zig`.
//!
//! Complex values are `std.math.Complex(F)` rather than the reference's
//! `da_r`/`da_i` over interleaved floats, so `util.at` needs no doubled stride.

const std = @import("std");
const zigstd = @import("../../../../zigstd.zig");
const util = @import("../../util.zig");
const types = @import("../../types.zig");
const cpu = @import("../../cpu/__root__.zig");
const blas = @import("../../__root__.zig");
const Vtable = @import("../../Vtable.zig");
const rsimd = @import("../../../../simd/__root__.zig");
const fma = @import("../../fma.zig");

const BlasInt = types.BlasInt;
const BlasLong = types.BlasLong;
const XDouble = types.XDouble;
const Complex = std.math.Complex;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

// ---- kernel ABI and drivers -------------------------------------------------

pub const kernel = @import("kernel.zig");
pub const RealKernel = kernel.RealKernel;
pub const ComplexKernel = kernel.ComplexKernel;
pub const realEntry = kernel.realEntry;
pub const complexEntry = kernel.complexEntry;
pub const Real = kernel.Real;
pub const Cplx = kernel.Cplx;

/// Portable block kernel: `y[0..n1] += alpha * x[0..n1]` on `@Vector`s, the
/// `V_SIMD` branch of `daxpy_kernel_8`/`saxpy_kernel_16`.
fn SimdBlock(comptime T: type) type {
    return struct {
        pub fn run(n1: BlasLong, alpha: T, x: [*]T, y: [*]T) void {
            const lanes = comptime rsimd.vectorLen(T);
            const V = @Vector(lanes, T);
            const valpha: V = @splat(alpha);
            var i: usize = 0;
            const end: usize = @intCast(n1);
            while (i < end) : (i += lanes) {
                const vx: V = x[i..][0..lanes].*;
                const vy: V = y[i..][0..lanes].*;
                y[i..][0..lanes].* = rsimd.mulAdd(V, valpha, vx, vy);
            }
        }
    };
}

/// Whether `T` gets the `@Vector` block kernel in the portable path.
fn hasSimdBlock(comptime T: type) bool {
    return T == f32 or T == f64;
}

fn GenericReal(comptime T: type) type {
    return Real(T, 128 / @sizeOf(T), if (hasSimdBlock(T)) SimdBlock(T) else null, fma.fuse);
}

fn GenericCplx(comptime F: type, comptime conj: bool) type {
    return Cplx(F, conj, 16, null, fma.fuse);
}

// ---- generic implementations (vtable targets) ------------------------------

pub const saxpy_generic = realEntry(f32, GenericReal(f32));
pub const daxpy_generic = realEntry(f64, GenericReal(f64));
pub const qaxpy_generic = realEntry(XDouble, GenericReal(XDouble));
pub const caxpy_generic = complexEntry(f32, GenericCplx(f32, false));
pub const caxpyc_generic = complexEntry(f32, GenericCplx(f32, true));
pub const zaxpy_generic = complexEntry(f64, GenericCplx(f64, false));
pub const zaxpyc_generic = complexEntry(f64, GenericCplx(f64, true));
pub const xaxpy_generic = complexEntry(XDouble, GenericCplx(XDouble, false));
pub const xaxpyc_generic = complexEntry(XDouble, GenericCplx(XDouble, true));

/// Points every axpy slot at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.saxpy_k = saxpy_generic;
    vt.daxpy_k = daxpy_generic;
    vt.qaxpy_k = qaxpy_generic;
    vt.caxpy_k = caxpy_generic;
    vt.caxpyc_k = caxpyc_generic;
    vt.zaxpy_k = zaxpy_generic;
    vt.zaxpyc_k = zaxpyc_generic;
    vt.xaxpy_k = xaxpy_generic;
    vt.xaxpyc_k = xaxpyc_generic;
}

// ---- typed kernel wrappers -------------------------------------------------
// Kernel semantics: strides may be negative, but then the pointers must
// already name the first element visited (`taxpy` does that adjustment).

pub fn saxpy(n: BlasLong, da: f32, x: [*]f32, incx: BlasLong, y: [*]f32, incy: BlasLong) void {
    _ = blas.kernel(.saxpy_k)(n, 0, 0, da, x, incx, y, incy, y, 0);
}

pub fn daxpy(n: BlasLong, da: f64, x: [*]f64, incx: BlasLong, y: [*]f64, incy: BlasLong) void {
    _ = blas.kernel(.daxpy_k)(n, 0, 0, da, x, incx, y, incy, y, 0);
}

pub fn qaxpy(n: BlasLong, da: XDouble, x: [*]XDouble, incx: BlasLong, y: [*]XDouble, incy: BlasLong) void {
    _ = blas.kernel(.qaxpy_k)(n, 0, 0, da, x, incx, y, incy, y, 0);
}

pub fn caxpy(n: BlasLong, da: Complex32, x: [*]Complex32, incx: BlasLong, y: [*]Complex32, incy: BlasLong) void {
    _ = blas.kernel(.caxpy_k)(n, 0, 0, da.re, da.im, @ptrCast(x), incx, @ptrCast(y), incy, @ptrCast(y), 0);
}

pub fn caxpyc(n: BlasLong, da: Complex32, x: [*]Complex32, incx: BlasLong, y: [*]Complex32, incy: BlasLong) void {
    _ = blas.kernel(.caxpyc_k)(n, 0, 0, da.re, da.im, @ptrCast(x), incx, @ptrCast(y), incy, @ptrCast(y), 0);
}

pub fn zaxpy(n: BlasLong, da: Complex64, x: [*]Complex64, incx: BlasLong, y: [*]Complex64, incy: BlasLong) void {
    _ = blas.kernel(.zaxpy_k)(n, 0, 0, da.re, da.im, @ptrCast(x), incx, @ptrCast(y), incy, @ptrCast(y), 0);
}

pub fn zaxpyc(n: BlasLong, da: Complex64, x: [*]Complex64, incx: BlasLong, y: [*]Complex64, incy: BlasLong) void {
    _ = blas.kernel(.zaxpyc_k)(n, 0, 0, da.re, da.im, @ptrCast(x), incx, @ptrCast(y), incy, @ptrCast(y), 0);
}

pub fn xaxpy(n: BlasLong, da: ComplexXD, x: [*]ComplexXD, incx: BlasLong, y: [*]ComplexXD, incy: BlasLong) void {
    _ = blas.kernel(.xaxpy_k)(n, 0, 0, da.re, da.im, @ptrCast(x), incx, @ptrCast(y), incy, @ptrCast(y), 0);
}

pub fn xaxpyc(n: BlasLong, da: ComplexXD, x: [*]ComplexXD, incx: BlasLong, y: [*]ComplexXD, incy: BlasLong) void {
    _ = blas.kernel(.xaxpyc_k)(n, 0, 0, da.re, da.im, @ptrCast(x), incx, @ptrCast(y), incy, @ptrCast(y), 0);
}

// ---- interface -------------------------------------------------------------

/// `SAXPY`/`DAXPY`/`CAXPY`/`ZAXPY` (and the `Q`/`X` extended forms), unified:
/// `DY := DA*DX + DY`. The precision is `da`'s type instead of a name prefix.
///
/// Follows the Fortran interface (`SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)`)
/// in argument names, order and meaning; see `taxpy`.
pub fn axpy(
    n: BlasInt,
    da: anytype,
    dx: [*]@TypeOf(da),
    incx: BlasInt,
    dy: [*]@TypeOf(da),
    incy: BlasInt,
) void {
    taxpy(@TypeOf(da), n, da, dx, incx, dy, incy);
}

/// `axpy` with the element type given explicitly.
///
/// Interface semantics (`interface/{,z}axpy.c`): returns for `n <= 0` or
/// `da == 0`; a negative stride names the vector's last element first, so the
/// pointer is moved to the start before calling the kernel; both strides zero
/// collapses to `y[0] += n*da*x[0]`.
///
/// This is the interface signature, not the kernel one (see `daxpy_k`).
pub fn taxpy(
    comptime T: type,
    n: BlasInt,
    da: T,
    dx: [*]T,
    incx: BlasInt,
    dy: [*]T,
    incy: BlasInt,
) void {
    taxpyImpl(T, false, n, da, dx, incx, dy, incy);
}

/// `axpy` with `x` conjugated: `CAXPYC`/`ZAXPYC`/`XAXPYC`. Complex only.
pub fn axpyc(
    n: BlasInt,
    da: anytype,
    dx: [*]@TypeOf(da),
    incx: BlasInt,
    dy: [*]@TypeOf(da),
    incy: BlasInt,
) void {
    taxpyc(@TypeOf(da), n, da, dx, incx, dy, incy);
}

pub fn taxpyc(
    comptime T: type,
    n: BlasInt,
    da: T,
    dx: [*]T,
    incx: BlasInt,
    dy: [*]T,
    incy: BlasInt,
) void {
    if (comptime !isComplex(T)) @compileError("taxpyc: " ++ @typeName(T) ++ " is not complex");
    taxpyImpl(T, true, n, da, dx, incx, dy, incy);
}

fn isComplex(comptime T: type) bool {
    return T == Complex32 or T == Complex64 or T == ComplexXD;
}

fn taxpyImpl(
    comptime T: type,
    comptime conj: bool,
    n: BlasInt,
    da: T,
    dx: [*]T,
    incx: BlasInt,
    dy: [*]T,
    incy: BlasInt,
) void {
    comptime types.assertBlasType(T, @src());

    if (n <= 0) return;

    const cplx = comptime isComplex(T);
    if (cplx) {
        if (da.re == 0 and da.im == 0) return;
    } else {
        if (da == 0) return;
    }

    const kn: BlasLong = n;
    const kincx: BlasLong = incx;
    const kincy: BlasLong = incy;

    if (kincx == 0 and kincy == 0) {
        // The reference applies this to the `CONJ` build too, unconjugated.
        if (cplx) {
            const F = @FieldType(T, "re");
            const cn: F = @floatFromInt(kn);
            dy[0].re += cn * (da.re * dx[0].re - da.im * dx[0].im);
            dy[0].im += cn * (da.im * dx[0].re + da.re * dx[0].im);
        } else {
            dy[0] += @as(T, @floatFromInt(kn)) * da * dx[0];
        }
        return;
    }

    const x: [*]T = if (kincx < 0) @ptrCast(util.at(T, dx, -(kn - 1) * kincx)) else dx;
    const y: [*]T = if (kincy < 0) @ptrCast(util.at(T, dy, -(kn - 1) * kincy)) else dy;

    if (comptime T == f32) {
        saxpy(kn, da, x, kincx, y, kincy);
    } else if (comptime T == f64) {
        daxpy(kn, da, x, kincx, y, kincy);
    } else if (comptime T == XDouble) {
        qaxpy(kn, da, x, kincx, y, kincy);
    } else if (comptime T == Complex32) {
        if (conj) caxpyc(kn, da, x, kincx, y, kincy) else caxpy(kn, da, x, kincx, y, kincy);
    } else if (comptime T == Complex64) {
        if (conj) zaxpyc(kn, da, x, kincx, y, kincy) else zaxpy(kn, da, x, kincx, y, kincy);
    } else if (comptime T == ComplexXD) {
        if (conj) xaxpyc(kn, da, x, kincx, y, kincy) else xaxpy(kn, da, x, kincx, y, kincy);
    } else {
        @compileError("axpy: no implementation for " ++ @typeName(T) ++ " (have: f32, f64, XDouble, Complex32, Complex64, ComplexXD)");
    }
}

// ---- tests: interface -------------------------------------------------------

test "axpy, taxpy and the typed wrappers agree" {
    var xs = [_]f64{ 1, 2, 3 };
    var y0 = [_]f64{ 1, 1, 1 };
    var y1 = y0;
    var y2 = y0;
    axpy(3, @as(f64, 2), &xs, 1, &y0, 1);
    taxpy(f64, 3, 2, &xs, 1, &y1, 1);
    daxpy(3, 2, &xs, 1, &y2, 1);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 3, 5, 7 }, &y0);
    try zigstd.testExpectEqualSlices(f64, &y0, &y1);
    try zigstd.testExpectEqualSlices(f64, &y0, &y2);
}

test "interface: n <= 0 and zero alpha leave y untouched" {
    var x = [_]f64{ 1, 2 };
    var y = [_]f64{ 5, 6 };
    axpy(0, @as(f64, 1), &x, 1, &y, 1);
    axpy(-3, @as(f64, 1), &x, 1, &y, 1);
    axpy(2, @as(f64, 0), &x, 1, &y, 1);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 5, 6 }, &y);

    var cx = [_]Complex64{Complex64.init(1, 1)};
    var cy = [_]Complex64{Complex64.init(5, 6)};
    axpy(1, Complex64.init(0, 0), &cx, 1, &cy, 1);
    try zigstd.testExpectEqual(Complex64.init(5, 6), cy[0]);
}

test "interface: zero strides collapse to n*alpha*x[0]" {
    var x = [_]f64{3};
    var y = [_]f64{1};
    axpy(4, @as(f64, 2), &x, 0, &y, 0);
    try zigstd.testExpectEqual(@as(f64, 25), y[0]);

    var cx = [_]Complex64{Complex64.init(1, 1)};
    var cy = [_]Complex64{Complex64.init(0, 0)};
    axpy(2, Complex64.init(2, 3), &cx, 0, &cy, 0);
    // 2 * ((2+3i)(1+i)) = 2 * (-1+5i)
    try zigstd.testExpectEqual(Complex64.init(-2, 10), cy[0]);
}

test "interface: negative stride names the last element first" {
    var x = [_]f64{ 1, 2, 3 };
    var y = [_]f64{ 0, 0, 0 };
    axpy(3, @as(f64, 1), &x, -1, &y, 1);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 3, 2, 1 }, &y);

    // Both negative: same element pairing as both positive.
    var y2 = [_]f64{ 0, 0, 0 };
    axpy(3, @as(f64, 1), &x, -1, &y2, -1);
    try zigstd.testExpectEqualSlices(f64, &x, &y2);

    var cx = [_]Complex32{ Complex32.init(1, 0), Complex32.init(2, 0), Complex32.init(3, 0) };
    var cy = [_]Complex32{ Complex32.init(0, 0), Complex32.init(0, 0), Complex32.init(0, 0) };
    axpy(3, Complex32.init(1, 0), &cx, 1, &cy, -1);
    try zigstd.testExpectEqualSlices(Complex32, &[_]Complex32{ Complex32.init(3, 0), Complex32.init(2, 0), Complex32.init(1, 0) }, &cy);
}

test "interface: axpyc conjugates x" {
    var x = [_]Complex64{Complex64.init(1, 1)};
    var y = [_]Complex64{Complex64.init(0, 0)};
    axpyc(1, Complex64.init(2, 3), &x, 1, &y, 1);
    // (2+3i) * conj(1+i) = (2+3i)(1-i) = 5 + i
    try zigstd.testExpectEqual(Complex64.init(5, 1), y[0]);
}

test "interface: extended-precision types" {
    var x = [_]XDouble{ 1, 2, 3 };
    var y = [_]XDouble{ 1, 1, 1 };
    axpy(3, @as(XDouble, 2), &x, 1, &y, 1);
    try zigstd.testExpectEqualSlices(XDouble, &[_]XDouble{ 3, 5, 7 }, &y);
}

// ---- tests: portable kernels vs. a scalar reference ------------------------

fn fill(comptime T: type, buf: []T, seed: u32) void {
    var prng = std.Random.DefaultPrng.init(seed);
    const r = prng.random();
    for (buf) |*v| {
        if (comptime isComplex(T)) {
            v.* = T.init(@floatCast(r.float(f64) * 4 - 2), @floatCast(r.float(f64) * 4 - 2));
        } else {
            v.* = @floatCast(r.float(f64) * 4 - 2);
        }
    }
}

fn FloatOf(comptime T: type) type {
    return if (isComplex(T)) @FieldType(T, "re") else T;
}

fn expectClose(comptime T: type, want: []const T, got: []const T) !void {
    const tol: FloatOf(T) = if (FloatOf(T) == f32) 1e-5 else 1e-12;
    try zigstd.testExpectEqual(want.len, got.len);
    for (want, got) |w, g| {
        if (comptime isComplex(T)) {
            try zigstd.testExpectApproxEqAbs(w.re, g.re, tol);
            try zigstd.testExpectApproxEqAbs(w.im, g.im, tol);
        } else {
            try zigstd.testExpectApproxEqAbs(w, g, tol);
        }
    }
}

/// Checks `kern` (vtable ABI) against a plain loop for contiguous and
/// strided data at sizes that straddle every block and unroll boundary.
pub fn checkReal(comptime T: type, kern: RealKernel(T)) !void {
    const sizes = [_]usize{ 0, 1, 3, 4, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 257 };
    const alpha: T = 1.5;
    for (sizes) |n| {
        var xb: [800]T = undefined;
        var yb: [800]T = undefined;
        var want: [800]T = undefined;
        fill(T, &xb, @intCast(n + 1));
        fill(T, &yb, @intCast(n + 7));
        inline for (.{ [2]usize{ 1, 1 }, [2]usize{ 2, 1 }, [2]usize{ 1, 3 }, [2]usize{ 2, 2 } }) |inc| {
            want = yb;
            var y = yb;
            for (0..n) |i| want[i * inc[1]] += alpha * xb[i * inc[0]];
            _ = kern(@intCast(n), 0, 0, alpha, &xb, inc[0], &y, inc[1], &y, 0);
            try expectClose(T, &want, &y);
        }
    }
}

pub fn checkComplex(comptime F: type, comptime conj: bool, kern: ComplexKernel(F)) !void {
    const C = Complex(F);
    const sizes = [_]usize{ 0, 1, 3, 4, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 129 };
    const alpha = C.init(1.5, -0.75);
    for (sizes) |n| {
        var xb: [400]C = undefined;
        var yb: [400]C = undefined;
        var want: [400]C = undefined;
        fill(C, &xb, @intCast(n + 1));
        fill(C, &yb, @intCast(n + 7));
        inline for (.{ [2]usize{ 1, 1 }, [2]usize{ 2, 1 }, [2]usize{ 1, 3 } }) |inc| {
            want = yb;
            var y = yb;
            for (0..n) |i| {
                const xv = xb[i * inc[0]];
                const yp = &want[i * inc[1]];
                if (conj) {
                    yp.re += alpha.re * xv.re + alpha.im * xv.im;
                    yp.im -= alpha.re * xv.im - alpha.im * xv.re;
                } else {
                    yp.re += alpha.re * xv.re - alpha.im * xv.im;
                    yp.im += alpha.re * xv.im + alpha.im * xv.re;
                }
            }
            _ = kern(@intCast(n), 0, 0, alpha.re, alpha.im, @ptrCast(&xb), inc[0], @ptrCast(&y), inc[1], @ptrCast(&y), 0);
            try expectClose(C, &want, &y);
        }
    }
}

test "generic kernels match a scalar loop (all types, contiguous and strided)" {
    try checkReal(f32, saxpy_generic);
    try checkReal(f64, daxpy_generic);
    try checkReal(XDouble, qaxpy_generic);
    try checkComplex(f32, false, caxpy_generic);
    try checkComplex(f32, true, caxpyc_generic);
    try checkComplex(f64, false, zaxpy_generic);
    try checkComplex(f64, true, zaxpyc_generic);
    try checkComplex(XDouble, false, xaxpy_generic);
    try checkComplex(XDouble, true, xaxpyc_generic);
}

test "a core's kernels fill every axpy slot, for the host's core type" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    cpu.kernel.install(&vt, cpu.detect());
    try zigstd.testExpectTrue(vt.saxpy_k != null and vt.daxpy_k != null and vt.qaxpy_k != null);
    try zigstd.testExpectTrue(vt.caxpy_k != null and vt.caxpyc_k != null);
    try zigstd.testExpectTrue(vt.zaxpy_k != null and vt.zaxpyc_k != null);
    try zigstd.testExpectTrue(vt.xaxpy_k != null and vt.xaxpyc_k != null);
    try checkReal(f32, vt.saxpy_k.?);
    try checkReal(f64, vt.daxpy_k.?);
    try checkComplex(f32, false, vt.caxpy_k.?);
    try checkComplex(f32, true, vt.caxpyc_k.?);
    try checkComplex(f64, false, vt.zaxpy_k.?);
    try checkComplex(f64, true, vt.zaxpyc_k.?);
}

test "every core's kernels match a scalar loop, where this machine runs them" {
    inline for (@typeInfo(cpu.CoreType).@"enum".field_names) |name| {
        const core = @field(cpu.CoreType, name);
        if (!@hasDecl(cpu.arch, "canRun") or cpu.arch.canRun(core)) {
            var vt: Vtable = .{};
            installGeneric(&vt);
            cpu.kernel.install(&vt, core);
            try checkReal(f32, vt.saxpy_k.?);
            try checkReal(f64, vt.daxpy_k.?);
            try checkComplex(f32, false, vt.caxpy_k.?);
            try checkComplex(f32, true, vt.caxpyc_k.?);
            try checkComplex(f64, false, vt.zaxpy_k.?);
            try checkComplex(f64, true, vt.zaxpyc_k.?);
        }
    }
}

test {
    zigstd.testRefAllDecls(@This());
}
