//! `GEMV` kernels: `y += alpha * op(A) * x`, `A` column-major `m` by `n`. The
//! vtable wrappers, the portable kernels (`*_generic`), and tests of every
//! core's kernels against plain arithmetic.
//!
//! Argument names are copied from the OpenBLAS kernel sources (see the
//! `CNAME(...)` definitions under `kernel/<arch>/*.c`). Several are named
//! `dummy*` in the reference: kernel-ABI padding slots this shape does not use.

const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");

const BlasLong = types.BlasLong;
const BlasInt = types.BlasInt;
const HF16 = types.HF16;
const BF16 = types.BF16;
const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

pub fn shgemv_n(m: BlasLong, n: BlasLong, alpha: f32, a: [*]HF16, lda: BlasLong, x: [*]HF16, inc_x: BlasLong, beta: f32, y: [*]f32, inc_y: BlasLong) i32 {
    return blas.kernel(.shgemv_n)(m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);
}
pub fn shgemv_t(m: BlasLong, n: BlasLong, alpha: f32, a: [*]HF16, lda: BlasLong, x: [*]HF16, inc_x: BlasLong, beta: f32, y: [*]f32, inc_y: BlasLong) i32 {
    return blas.kernel(.shgemv_t)(m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);
}
pub fn bgemv_n(m: BlasLong, n: BlasLong, alpha: BF16, a: [*]BF16, lda: BlasLong, x: [*]BF16, inc_x: BlasLong, beta: BF16, y: [*]BF16, inc_y: BlasLong) i32 {
    return blas.kernel(.bgemv_n)(m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);
}
pub fn bgemv_t(m: BlasLong, n: BlasLong, alpha: BF16, a: [*]BF16, lda: BlasLong, x: [*]BF16, inc_x: BlasLong, beta: BF16, y: [*]BF16, inc_y: BlasLong) i32 {
    return blas.kernel(.bgemv_t)(m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);
}
pub fn sbgemv_n(m: BlasLong, n: BlasLong, alpha: f32, a: [*]BF16, lda: BlasLong, x: [*]BF16, inc_x: BlasLong, beta: f32, y: [*]f32, inc_y: BlasLong) i32 {
    return blas.kernel(.sbgemv_n)(m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);
}
pub fn sbgemv_t(m: BlasLong, n: BlasLong, alpha: f32, a: [*]BF16, lda: BlasLong, x: [*]BF16, inc_x: BlasLong, beta: f32, y: [*]f32, inc_y: BlasLong) i32 {
    return blas.kernel(.sbgemv_t)(m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);
}
pub fn sgemv_n(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.sgemv_n)(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn sgemv_t(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.sgemv_t)(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn dgemv_n(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.dgemv_n)(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn dgemv_t(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.dgemv_t)(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn qgemv_n(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.qgemv_n)(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn qgemv_t(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.qgemv_t)(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_n(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgemv_n)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_t(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgemv_t)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_r(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgemv_r)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_c(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgemv_c)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_o(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgemv_o)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_u(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgemv_u)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_s(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgemv_s)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn cgemv_d(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, a: [*]f32, lda: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgemv_d)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_n(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgemv_n)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_t(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgemv_t)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_r(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgemv_r)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_c(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgemv_c)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_o(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgemv_o)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_u(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgemv_u)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_s(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgemv_s)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn zgemv_d(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, a: [*]f64, lda: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgemv_d)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_n(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgemv_n)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_t(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgemv_t)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_r(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgemv_r)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_c(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgemv_c)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_o(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgemv_o)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_u(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgemv_u)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_s(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgemv_s)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}
pub fn xgemv_d(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, a: [*]XDouble, lda: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgemv_d)(m, n, dummy1, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../../blas1/check.zig");
const half = @import("half.zig");

/// `kernel/arm/gemv_n.c` and `gemv_t.c`: `y += alpha*A*x` (`n`, `y` of `m`)
/// or `y += alpha*A^T*x` (`t`, `y` of `n`). The interface has already scaled
/// `y` by `beta`.
fn Gemv(comptime F: type) type {
    return struct {
        fn n(m: BlasLong, cols: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(.c) i32 {
            var j: BlasLong = 0;
            while (j < cols) : (j += 1) {
                const tv = alpha * util.at(F, x, j * inc_x).*;
                var i: BlasLong = 0;
                while (i < m) : (i += 1) util.at(F, y, i * inc_y).* += tv * util.at(F, a, i + j * lda).*;
            }
            return 0;
        }
        fn t(m: BlasLong, cols: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(.c) i32 {
            var j: BlasLong = 0;
            while (j < cols) : (j += 1) {
                var acc: F = 0;
                var i: BlasLong = 0;
                while (i < m) : (i += 1) acc += util.at(F, a, i + j * lda).* * util.at(F, x, i * inc_x).*;
                util.at(F, y, j * inc_y).* += alpha * acc;
            }
            return 0;
        }
    };
}

/// `kernel/arm/zgemv_n.c` and `zgemv_t.c`: `y += alpha*op(A)*x` with `A`
/// conjugated for `conj` (`CONJ`) and `x` for `xconj` (`XCONJ`): `r` is
/// `conj(A)*x`, `o` is `A*conj(x)`, `c` is `A^H*x`, and so on.
fn Zgemv(comptime F: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool) type {
    return struct {
        fn get(p: [*]F, i: BlasLong, c: bool) std.math.Complex(F) {
            const re = util.at(F, p, 2 * i).*;
            const im = util.at(F, p, 2 * i + 1).*;
            return .init(re, if (c) -im else im);
        }
        fn k(m: BlasLong, cols: BlasLong, _: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(.c) i32 {
            const alpha: std.math.Complex(F) = .init(alpha_r, alpha_i);
            var j: BlasLong = 0;
            while (j < cols) : (j += 1) {
                if (trans) {
                    var acc: std.math.Complex(F) = .init(0, 0);
                    var i: BlasLong = 0;
                    while (i < m) : (i += 1) acc = acc.add(get(a, i + j * lda, conj).mul(get(x, i * inc_x, xconj)));
                    const v = alpha.mul(acc);
                    util.at(F, y, 2 * j * inc_y).* += v.re;
                    util.at(F, y, 2 * j * inc_y + 1).* += v.im;
                } else {
                    const tv = alpha.mul(get(x, j * inc_x, xconj));
                    var i: BlasLong = 0;
                    while (i < m) : (i += 1) {
                        const v = tv.mul(get(a, i + j * lda, conj));
                        util.at(F, y, 2 * i * inc_y).* += v.re;
                        util.at(F, y, 2 * i * inc_y + 1).* += v.im;
                    }
                }
            }
            return 0;
        }
    };
}

/// The complex variants by suffix: transposed, `A` conjugated, `x` conjugated.
const variants = .{
    .{ "n", false, false, false }, .{ "t", true, false, false },
    .{ "r", false, true, false },  .{ "c", true, true, false },
    .{ "o", false, false, true },  .{ "u", true, false, true },
    .{ "s", false, true, true },   .{ "d", true, true, true },
};

/// Points every gemv slot at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    half.installGeneric(vt);
    vt.sgemv_n = &Gemv(f32).n;
    vt.sgemv_t = &Gemv(f32).t;
    vt.dgemv_n = &Gemv(f64).n;
    vt.dgemv_t = &Gemv(f64).t;
    vt.qgemv_n = &Gemv(XDouble).n;
    vt.qgemv_t = &Gemv(XDouble).t;
    inline for (.{ .{ "c", f32 }, .{ "z", f64 }, .{ "x", XDouble } }) |p| {
        inline for (variants) |v| {
            @field(vt, p[0] ++ "gemv_" ++ v[0]) = &Zgemv(p[1], v[1], v[2], v[3]).k;
        }
    }
}

// ---- tests -----------------------------------------------------------------

const C = std.math.Complex;
/// Shapes across the kernels' blocks; small enough for the rigs' stacks.
const sizes = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 }, .{ 3, 2 }, .{ 4, 4 }, .{ 7, 5 }, .{ 8, 8 }, .{ 9, 16 }, .{ 16, 9 }, .{ 17, 17 }, .{ 33, 4 }, .{ 4, 33 } };
const strides = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 } };
const cap = 33 * 3 + 2; // a vector's span at stride 3

/// `kern` against plain arithmetic; `T` is the element type and the flags
/// are the variant's.
fn checkGemv(comptime T: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    const alpha: T = if (complex) .init(1.25, -0.5) else 1.25;
    for (sizes) |mn| for (strides) |inc| {
        const m = mn[0];
        const n = mn[1];
        const lda = m + 1;
        var a: [400]T = undefined; // lda * n at most 18 * 17
        var x: [cap]T = undefined;
        var y: [cap]T = undefined;
        check.fill(T, &a, m * 7 + n);
        check.fill(T, &x, m + 3);
        check.fill(T, &y, n + 5);
        var want = y;
        const rows = if (trans) n else m; // y's length
        const inner = if (trans) m else n; // x's length
        for (0..rows) |r| {
            var acc: T = if (complex) .init(0, 0) else 0;
            for (0..inner) |q| {
                const e = if (trans) a[q + r * lda] else a[r + q * lda];
                const xv = x[q * inc[0]];
                if (complex) {
                    const ae = if (conj) e.conjugate() else e;
                    const xe = if (xconj) xv.conjugate() else xv;
                    acc = acc.add(ae.mul(xe));
                } else acc += e * xv;
            }
            if (complex) want[r * inc[1]] = want[r * inc[1]].add(alpha.mul(acc)) else want[r * inc[1]] += alpha * acc;
        }
        if (complex) {
            _ = kern(@intCast(m), @intCast(n), 0, alpha.re, alpha.im, @ptrCast(&a), @intCast(lda), @ptrCast(&x), @intCast(inc[0]), @ptrCast(&y), @intCast(inc[1]), check.work(F));
        } else {
            _ = kern(@intCast(m), @intCast(n), 0, alpha, &a, @intCast(lda), &x, @intCast(inc[0]), &y, @intCast(inc[1]), check.work(F));
        }
        check.expectCloseSlices(T, &want, &y, check.tolerance(F, 4 * inner)) catch |e| {
            check.say("m={d} n={d} trans={} conj={} xconj={}: ", .{ m, n, trans, conj, xconj });
            return check.case(e, "gemv", T, rows, inc);
        };
    };
}

fn checkAll(vt: *const Vtable) !void {
    try checkGemv(f32, false, false, false, vt.sgemv_n.?);
    try checkGemv(f32, true, false, false, vt.sgemv_t.?);
    try checkGemv(f64, false, false, false, vt.dgemv_n.?);
    try checkGemv(f64, true, false, false, vt.dgemv_t.?);
    inline for (.{ .{ "c", f32 }, .{ "z", f64 } }) |p| {
        inline for (variants) |v| {
            try checkGemv(C(p[1]), v[1], v[2], v[3], @field(vt, p[0] ++ "gemv_" ++ v[0]).?);
        }
    }
}

// the half-precision kernels and their tests
comptime {
    _ = half;
}

test "portable gemv kernels match plain arithmetic" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkGemv(XDouble, false, false, false, vt.qgemv_n.?);
    try checkGemv(XDouble, true, false, false, vt.qgemv_t.?);
    inline for (variants) |v| try checkGemv(C(XDouble), v[1], v[2], v[3], @field(vt, "xgemv_" ++ v[0]).?);
}

test "every linked core's gemv kernels match plain arithmetic, where this machine runs them" {
    try check.eachCore(checkAll);
}
