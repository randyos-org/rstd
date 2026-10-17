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
const x87 = @import("x87.zig");
const sym = @import("../symmetric.zig");
const guard = @import("../guard.zig");

/// `kernel/arm/gemv_n.c` and `gemv_t.c`: `y += alpha*A*x` (`n`, `y` of `m`)
/// or `y += alpha*A^T*x` (`t`, `y` of `n`). The interface has already scaled
/// `y` by `beta`. Wider than `f64`, the x87 `qgemv_n.S`/`qgemv_t.S` sums
/// instead (`x87.zig`).
fn Gemv(comptime F: type) type {
    return struct {
        fn n(m: BlasLong, cols: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
            if (comptime types.isExtended(F)) {
                x87.gemv(F, false, x87.colBlock(false, .unknown), m, cols, alpha, a, lda, x, inc_x, y, inc_y);
                return 0;
            }
            var j: BlasLong = 0;
            while (j < cols) : (j += 1) {
                const tv = alpha * util.at(F, x, j * inc_x).*;
                var i: BlasLong = 0;
                while (i < m) : (i += 1) util.at(F, y, i * inc_y).* += tv * util.at(F, a, i + j * lda).*;
            }
            return 0;
        }
        fn t(m: BlasLong, cols: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
            if (comptime types.isExtended(F)) {
                x87.gemv(F, true, x87.rowBlock(false), m, cols, alpha, a, lda, x, inc_x, y, inc_y);
                return 0;
            }
            var j: BlasLong = 0;
            while (j < cols) : (j += 1) {
                const A = types.Acc(F);
                var acc: A = 0;
                var i: BlasLong = 0;
                while (i < m) : (i += 1) acc += @as(A, util.at(F, a, i + j * lda).*) * util.at(F, x, i * inc_x).*;
                util.at(F, y, j * inc_y).* += @floatCast(alpha * acc);
            }
            return 0;
        }
    };
}

/// `kernel/arm/zgemv_n.c` and `zgemv_t.c`: `y += alpha*op(A)*x` with `A`
/// conjugated for `conj` (`CONJ`) and `x` for `xconj` (`XCONJ`): `r` is
/// `conj(A)*x`, `o` is `A*conj(x)`, `c` is `A^H*x`, and so on. Wider than
/// `f64`, the x87 `xgemv_n.S`/`xgemv_t.S` sums instead (`x87.zig`).
fn Zgemv(comptime F: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool) type {
    const W = std.math.Complex(types.Acc(F));
    return struct {
        fn wide(v: std.math.Complex(F)) W {
            return .init(v.re, v.im);
        }
        fn get(p: [*]F, i: BlasLong, c: bool) std.math.Complex(F) {
            const re = util.at(F, p, 2 * i).*;
            const im = util.at(F, p, 2 * i + 1).*;
            return .init(re, if (c) -im else im);
        }
        fn k(m: BlasLong, cols: BlasLong, _: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
            if (comptime types.isExtended(F)) {
                const p = comptime if (trans) x87.rowBlock(true) else x87.colBlock(true, .unknown);
                x87.zgemv(F, trans, conj, xconj, p, m, cols, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y);
                return 0;
            }
            const alpha: std.math.Complex(F) = .init(alpha_r, alpha_i);
            var j: BlasLong = 0;
            while (j < cols) : (j += 1) {
                if (trans) {
                    var acc: W = .init(0, 0);
                    var i: BlasLong = 0;
                    while (i < m) : (i += 1) acc = acc.add(wide(get(a, i + j * lda, conj)).mul(wide(get(x, i * inc_x, xconj))));
                    const v = wide(alpha).mul(acc);
                    util.at(F, y, 2 * j * inc_y).* += @floatCast(v.re);
                    util.at(F, y, 2 * j * inc_y + 1).* += @floatCast(v.im);
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

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    if (!types.isComplex(T)) return struct {
        pub const gemv_n = &Gemv(F).n;
        pub const gemv_t = &Gemv(F).t;
    };
    return struct {
        pub const gemv_n = &Zgemv(F, false, false, false).k;
        pub const gemv_t = &Zgemv(F, true, false, false).k;
        pub const gemv_r = &Zgemv(F, false, true, false).k;
        pub const gemv_c = &Zgemv(F, true, true, false).k;
        pub const gemv_o = &Zgemv(F, false, false, true).k;
        pub const gemv_u = &Zgemv(F, true, false, true).k;
        pub const gemv_s = &Zgemv(F, false, true, true).k;
        pub const gemv_d = &Zgemv(F, true, true, true).k;
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
/// Shapes across the kernels' blocks, then bigger ones through their
/// multi-block loops (storage is static: `guard.operand`).
const sizes = [_][2]usize{
    .{ 1, 1 },   .{ 2, 3 },    .{ 3, 2 },    .{ 4, 4 },   .{ 7, 5 },   .{ 8, 8 }, .{ 9, 16 }, .{ 16, 9 }, .{ 17, 17 }, .{ 33, 4 }, .{ 4, 33 },
    .{ 65, 65 }, .{ 100, 37 }, .{ 37, 100 }, .{ 129, 3 }, .{ 3, 129 },
};
const strides = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 } };

/// `kern` against plain arithmetic; `T` is the element type and the flags
/// are the variant's. Its buffer is the interface's, exactly (`sym.gemvLen`),
/// 32-byte aligned as a stack one is, with nothing written outside it.
fn checkGemv(comptime T: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool, kern: anytype) !void {
    try checkGemvAt(T, trans, conj, xconj, kern, &sizes, &strides, null);
}

/// Past the x87 kernels' 4096-row blocks along the sum (`t` forms) and past
/// every column block (`n` forms, `DTB_ENTRIES` and `x87.colBlock`); on the
/// rigs only the shape that crosses the form's own block, at unit stride.
const big_sizes = [_][2]usize{ .{ 4100, 5 }, .{ 5, 4100 } };

fn checkGemvBig(comptime T: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool, kern: anytype) !void {
    const rig = check.rig();
    const shapes: []const [2]usize = if (!rig) &big_sizes else if (trans) big_sizes[0..1] else big_sizes[1..2];
    const incs: []const [2]usize = if (rig) strides[0..1] else &strides;
    try checkGemvAt(T, trans, conj, xconj, kern, shapes, incs, std.testing.allocator);
}

/// `checkGemv` over `shapes` and `incs`, its operands from `al` where given,
/// else static.
fn checkGemvAt(comptime T: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool, kern: anytype, shapes: []const [2]usize, incs: []const [2]usize, al: ?std.mem.Allocator) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    const cs = if (complex) 2 else 1;
    const alpha: T = if (complex) .init(1.25, -0.5) else 1.25;
    for (shapes) |mn| for (incs) |inc| {
        const m = mn[0];
        const n = mn[1];
        const lda = m + 1;
        const rows = if (trans) n else m; // y's length
        const inner = if (trans) m else n; // x's length
        const ops = try Operands(T).get(al, lda * n, guard.span(inner, inc[0]), guard.span(rows, inc[1]));
        defer ops.free(al);
        const a = ops.a;
        const x = ops.x;
        const y = ops.y;
        const want = ops.want;
        check.fill(T, a, m * 7 + n);
        check.fill(T, x, m + 3);
        check.fill(T, y, n + 5);
        @memcpy(want, y);
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
        const bytes = sym.gemvLen(F, cs, m, n) * @sizeOf(F);
        // x86_64 `zgemv_t.S` needs a buffer free of NaN, which its callers clear
        const buf = guard.open(F, bytes, 32, if (complex and trans) .junk else .nan);
        if (complex) {
            _ = kern(@intCast(m), @intCast(n), 0, alpha.re, alpha.im, @ptrCast(a.ptr), @intCast(lda), @ptrCast(x.ptr), @intCast(inc[0]), @ptrCast(y.ptr), @intCast(inc[1]), buf);
        } else {
            _ = kern(@intCast(m), @intCast(n), 0, alpha, a.ptr, @intCast(lda), x.ptr, @intCast(inc[0]), y.ptr, @intCast(inc[1]), buf);
        }
        check.expectCloseSlices(T, want, y, check.tolerance(F, 4 * inner)) catch |e| {
            check.say("m={d} n={d} trans={} conj={} xconj={}: ", .{ m, n, trans, conj, xconj });
            return check.case(e, "gemv", T, rows, inc);
        };
        guard.intact(bytes, 32) catch |e| {
            check.say("m={d} n={d} trans={} conj={} xconj={}: ", .{ m, n, trans, conj, xconj });
            return check.case(e, "gemv buffer", T, rows, inc);
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

/// `checkGemv`'s operands: static, or from an allocator for big shapes.
fn Operands(comptime T: type) type {
    return struct {
        a: []T,
        x: []T,
        y: []T,
        want: []T,

        fn get(al: ?std.mem.Allocator, na: usize, nx: usize, ny: usize) !@This() {
            const g = al orelse return .{ .a = guard.operand(T, .a, na), .x = guard.operand(T, .x, nx), .y = guard.operand(T, .y, ny), .want = guard.operand(T, .want, ny) };
            const a = try g.alloc(T, na);
            errdefer g.free(a);
            const x = try g.alloc(T, nx);
            errdefer g.free(x);
            const y = try g.alloc(T, ny);
            errdefer g.free(y);
            return .{ .a = a, .x = x, .y = y, .want = try g.alloc(T, ny) };
        }

        fn free(self: @This(), al: ?std.mem.Allocator) void {
            const g = al orelse return;
            g.free(self.want);
            g.free(self.y);
            g.free(self.x);
            g.free(self.a);
        }
    };
}

/// `checkAll` past the blocks.
fn checkAllBig(vt: *const Vtable) !void {
    try checkGemvBig(f32, false, false, false, vt.sgemv_n.?);
    try checkGemvBig(f32, true, false, false, vt.sgemv_t.?);
    try checkGemvBig(f64, false, false, false, vt.dgemv_n.?);
    try checkGemvBig(f64, true, false, false, vt.dgemv_t.?);
    inline for (.{ .{ "c", f32 }, .{ "z", f64 } }) |p| {
        inline for (variants) |v| {
            try checkGemvBig(C(p[1]), v[1], v[2], v[3], @field(vt, p[0] ++ "gemv_" ++ v[0]).?);
        }
    }
}

// the half-precision and x87 kernels and their tests
comptime {
    _ = half;
    _ = x87;
}

test "portable gemv kernels match plain arithmetic" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkGemv(XDouble, false, false, false, vt.qgemv_n.?);
    try checkGemv(XDouble, true, false, false, vt.qgemv_t.?);
    inline for (variants) |v| try checkGemv(C(XDouble), v[1], v[2], v[3], @field(vt, "xgemv_" ++ v[0]).?);
}

test "portable gemv kernels at the extended type match plain arithmetic" {
    const E = check.Ext;
    try checkGemv(E, false, false, false, &Gemv(E).n);
    try checkGemv(E, true, false, false, &Gemv(E).t);
    inline for (variants) |v| try checkGemv(C(E), v[1], v[2], v[3], &Zgemv(E, v[1], v[2], v[3]).k);
}

test "every linked core's gemv kernels match plain arithmetic, where this machine runs them" {
    try check.eachCore(checkAll);
}

test "portable gemv kernels past the x87 row block and every column block" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAllBig(&vt);
    try checkGemvBig(XDouble, false, false, false, vt.qgemv_n.?);
    try checkGemvBig(XDouble, true, false, false, vt.qgemv_t.?);
    inline for (variants) |v| try checkGemvBig(C(XDouble), v[1], v[2], v[3], @field(vt, "xgemv_" ++ v[0]).?);
}

test "every linked core's gemv kernels past the x87 row block and every column block" {
    try check.eachCore(checkAllBig);
}
