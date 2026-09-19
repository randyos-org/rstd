//! Shared pieces of the axpy kernels: the reference's kernel ABI and the
//! `CNAME` drivers of `kernel/x86_64/{s,d,c,z}axpy.c`. Nothing here depends on
//! the BLAS root, so the per-core kernel objects (`blaskernel_*.zig`) can use
//! it too.
//!
//! A driver is an optional block kernel for the contiguous part plus the
//! reference's scalar tail and strided loops. `fuse` says whether the
//! compiler flags of the core being built contract multiply-adds (see
//! `fma.zig`); it is applied per statement the way clang's contraction is:
//! only within a single expression, so the strided loop's separate
//! `m = da * x; y += m` statements stay un-fused.

const std = @import("std");
const util = @import("../../util.zig");
const types = @import("../../types.zig");

const BlasLong = types.BlasLong;
const Complex = std.math.Complex;

// ---- kernel ABI ------------------------------------------------------------

/// Vtable slot type for the real kernels (`saxpy_k`, `daxpy_k`, `qaxpy_k`).
pub fn RealKernel(comptime T: type) type {
    return *const fn (BlasLong, BlasLong, BlasLong, T, [*]T, BlasLong, [*]T, BlasLong, [*]T, BlasLong) callconv(.c) i32;
}

/// Vtable slot type for the complex kernels; `F` is the component type.
pub fn ComplexKernel(comptime F: type) type {
    return *const fn (BlasLong, BlasLong, BlasLong, F, F, [*]F, BlasLong, [*]F, BlasLong, [*]F, BlasLong) callconv(.c) i32;
}

/// Wraps `Impl.call(n, alpha, x, incx, y, incy)` in the real kernel ABI.
pub fn realEntry(comptime T: type, comptime Impl: type) RealKernel(T) {
    return &struct {
        fn f(n: BlasLong, _: BlasLong, _: BlasLong, da: T, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong, _: [*]T, _: BlasLong) callconv(.c) i32 {
            Impl.call(n, da, x, incx, y, incy);
            return 0;
        }
    }.f;
}

/// Wraps `Impl.call` in the complex kernel ABI (interleaved `F` pairs).
pub fn complexEntry(comptime F: type, comptime Impl: type) ComplexKernel(F) {
    return &struct {
        fn f(n: BlasLong, _: BlasLong, _: BlasLong, da_r: F, da_i: F, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F, _: BlasLong) callconv(.c) i32 {
            Impl.call(n, Complex(F).init(da_r, da_i), @ptrCast(x), incx, @ptrCast(y), incy);
            return 0;
        }
    }.f;
}

// ---- drivers ---------------------------------------------------------------

/// `a*b + c`, contracted or not.
pub inline fn fmadd(comptime fuse: bool, comptime T: type, a: T, b: T, c: T) T {
    return if (fuse) @mulAdd(T, a, b, c) else a * b + c;
}

/// Real kernel driver. `mask` is the reference's `n & -mask` block size
/// (16 for `f64`, 32 for `f32`); `Block.run(n1, alpha, x, y)` handles that
/// many contiguous elements, or `Block` is null for a scalar-only kernel.
pub fn Real(comptime T: type, comptime mask: comptime_int, comptime Block: ?type, comptime fuse: bool) type {
    return struct {
        pub fn call(n: BlasLong, alpha: T, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong) void {
            if (n <= 0) return;

            if (incx == 1 and incy == 1) {
                var i: BlasLong = 0;
                if (comptime Block) |B| {
                    i = n & -mask;
                    if (i != 0) B.run(i, alpha, x, y);
                }
                while (i < n) : (i += 1) {
                    const idx: usize = @intCast(i);
                    y[idx] = fmadd(fuse, T, alpha, x[idx], y[idx]);
                }
                return;
            }

            const n1: BlasLong = n & -4;
            var i: BlasLong = 0;
            var ix: BlasLong = 0;
            var iy: BlasLong = 0;
            while (i < n1) : (i += 4) {
                const m1 = alpha * util.at(T, x, ix).*;
                const m2 = alpha * util.at(T, x, ix + incx).*;
                const m3 = alpha * util.at(T, x, ix + 2 * incx).*;
                const m4 = alpha * util.at(T, x, ix + 3 * incx).*;

                util.at(T, y, iy).* += m1;
                util.at(T, y, iy + incy).* += m2;
                util.at(T, y, iy + 2 * incy).* += m3;
                util.at(T, y, iy + 3 * incy).* += m4;

                ix += incx * 4;
                iy += incy * 4;
            }
            while (i < n) : (i += 1) {
                const yp = util.at(T, y, iy);
                yp.* = fmadd(fuse, T, alpha, util.at(T, x, ix).*, yp.*);
                ix += incx;
                iy += incy;
            }
        }
    };
}

/// Complex kernel driver; `conj` selects the `CONJ` (`*axpyc`) formula.
/// `Block.run(n1, alpha, x, y)` (complex pointers) handles `n & -mask`
/// contiguous elements, or null for scalar only.
pub fn Cplx(comptime F: type, comptime conj: bool, comptime mask: comptime_int, comptime Block: ?type, comptime fuse: bool) type {
    const C = Complex(F);
    return struct {
        // One C statement per component. Where clang contracts, it fuses the
        // first product of the parenthesized sum/difference and leaves the
        // outer `y +=` alone (`fmuladd(ar, xr, -(ai*xi))`, then a plain add).
        inline fn madd(y: *C, alpha: C, x: C) void {
            if (comptime !conj) {
                y.re += fmadd(fuse, F, alpha.re, x.re, -(alpha.im * x.im));
                y.im += fmadd(fuse, F, alpha.re, x.im, alpha.im * x.re);
            } else {
                y.re += fmadd(fuse, F, alpha.re, x.re, alpha.im * x.im);
                y.im -= fmadd(fuse, F, alpha.re, x.im, -(alpha.im * x.re));
            }
        }

        pub fn call(n: BlasLong, alpha: C, x: [*]C, incx: BlasLong, y: [*]C, incy: BlasLong) void {
            if (n <= 0) return;

            var i: BlasLong = 0;
            if (incx == 1 and incy == 1) {
                if (comptime Block) |B| {
                    i = n & -mask;
                    if (i != 0) B.run(i, alpha, x, y);
                }
                while (i < n) : (i += 1) {
                    const idx: usize = @intCast(i);
                    madd(&y[idx], alpha, x[idx]);
                }
                return;
            }

            var ix: BlasLong = 0;
            var iy: BlasLong = 0;
            while (i < n) : (i += 1) {
                madd(util.at(C, y, iy), alpha, util.at(C, x, ix).*);
                ix += incx;
                iy += incy;
            }
        }
    };
}
