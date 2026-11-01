//! `GEMV` kernels: `y += alpha * op(A) * x`, `A` column-major `m` by `n`. The
//! vtable wrappers, the portable kernels (`*_generic`), the extended types'
//! x87 sums, and tests of every core's kernels against plain arithmetic.
//!
//! The half-precision gemvs (`shgemv` `f16` in and `f32` out, `sbgemv` `bf16`
//! in and `f32` out, `bgemv` `bf16` throughout) compute
//! `y = alpha*op(A)*x + beta*y`; their kernels are `ext/bf16/kernel.zig`'s
//! `Gemv`, the same on every core. See the readme.
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
const sysinfo = @import("builtin");
const gemm_params = @import("../../blas3/gemm/params.zig");
const bf16_kernel = @import("../../ext/bf16/kernel.zig");
const vector = @import("../../vector.zig");

/// `kernel/arm/gemv_n.c` and `gemv_t.c`: `y += alpha*A*x` (`n`, `y` of `m`)
/// or `y += alpha*A^T*x` (`t`, `y` of `n`). The interface has already scaled
/// `y` by `beta`. Wider than `f64`, the x87 `qgemv_n.S`/`qgemv_t.S` sums
/// instead (`extendedGemv`).
fn Gemv(comptime F: type) type {
    return struct {
        fn n(m: BlasLong, cols: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
            if (comptime types.isExtended(F)) {
                extendedGemv(F, false, extendedColBlock(false, .unknown), m, cols, alpha, a, lda, x, inc_x, y, inc_y);
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
                extendedGemv(F, true, extendedRowBlock(false), m, cols, alpha, a, lda, x, inc_x, y, inc_y);
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
/// `f64`, the x87 `xgemv_n.S`/`xgemv_t.S` sums instead (`extendedZgemv`).
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
                const p = comptime if (trans) extendedRowBlock(true) else extendedColBlock(true, .unknown);
                extendedZgemv(F, trans, conj, xconj, p, m, cols, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y);
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

// ---- extended types ------------------------------------------------------

// The sums of OpenBLAS's x87 `q`/`x` gemv kernels (`kernel/x86_64/qgemv_n.S`,
// `qgemv_t.S`, `xgemv_n.S`, `xgemv_t.S`; on x86, `kernel/x86/`'s, which differ
// only in their blocks), which the portable gemv takes for floats wider than
// `f64`. At `f80` they give the assembly's bits: the same blocks and
// accumulators, filled and folded in its order, with nothing fused.

const X86CoreType = gemm_params.CoreOf(.x86);

const x86 = sysinfo.target.cpu.arch == .x86;

/// Rows a block of `qgemv_t.S` (`complex`: `xgemv_t.S`): 4096 on x86_64,
/// 1000 and 400 on x86.
fn extendedRowBlock(comptime complex: bool) BlasLong {
    return if (!x86) 4096 else if (complex) 400 else 1000;
}

/// Columns a block of `qgemv_n.S` (`complex`: `xgemv_n.S`): 32 on x86_64. On
/// x86, `core`'s `DTB_DEFAULT_ENTRIES`, halved for `xgemv_n` on `PENTIUM4` and
/// `ATHLON` (`qgemv_n`'s 32 on `ATHLON` and `OPTERON` is their DTB anyway).
fn extendedColBlock(comptime complex: bool, comptime core: X86CoreType) BlasLong {
    if (!x86) return 32;
    const dtb = gemm_params.dtbFor(.x86, core);
    const half = complex and (core == .northwood or core == .prescott or core == .athlon);
    return if (half) @divExact(dtb, 2) else dtb;
}

/// `qgemv_n.S` and `qgemv_t.S`: per block of `p` along the sum, each element
/// of `y` in its whole fours is one sequential sum; the one to three left take
/// four accumulators by index % 4, eight a trip, the tail into the first,
/// folded ((a0 + a1) + a2) + a3. Then `y += sum * alpha`, block by block.
fn extendedGemv(comptime F: type, comptime trans: bool, comptime p: BlasLong, m: BlasLong, n: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) void {
    const outs = if (trans) n else m;
    const inner = if (trans) m else n;
    const seq = outs & -4;
    // `A`'s step along the sum and along `y`
    const step = if (trans) 1 else lda;
    const across = if (trans) lda else 1;
    var s: BlasLong = 0;
    while (s < inner) : (s += p) {
        const len = @min(p, inner - s);
        var o: BlasLong = 0;
        while (o < outs) : (o += 1) {
            const ao = s * step + o * across;
            var c: [4]F = @splat(0);
            var k: BlasLong = 0;
            if (o < seq) {
                while (k < len) : (k += 1) c[0] += util.at(F, a, ao + k * step).* * util.at(F, x, (s + k) * inc_x).*;
            } else {
                const main = len & -8;
                while (k < main) : (k += 8) {
                    inline for (0..8) |r| c[r % 4] += util.at(F, a, ao + (k + r) * step).* * util.at(F, x, (s + k + r) * inc_x).*;
                }
                while (k < len) : (k += 1) c[0] += util.at(F, a, ao + k * step).* * util.at(F, x, (s + k) * inc_x).*;
                c[0] = ((c[0] + c[1]) + c[2]) + c[3];
            }
            util.at(F, y, o * inc_y).* += c[0] * alpha;
        }
    }
}

/// `xgemv_n.S` and `xgemv_t.S`: per block of `p` along the sum, each element
/// of `y` sums `ar*xr`, `ai*xr` (subtracted for `conj`), `ar*xi` and `ai*xi`
/// apart, combines them as the variant's signs say, and adds `alpha` times it.
fn extendedZgemv(comptime F: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool, comptime p: BlasLong, m: BlasLong, n: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) void {
    const outs = if (trans) n else m;
    const inner = if (trans) m else n;
    const step = if (trans) 1 else lda;
    const across = if (trans) lda else 1;
    var s: BlasLong = 0;
    while (s < inner) : (s += p) {
        const len = @min(p, inner - s);
        var o: BlasLong = 0;
        while (o < outs) : (o += 1) {
            var c: [4]F = @splat(0);
            var k: BlasLong = 0;
            while (k < len) : (k += 1) {
                const e = 2 * (s * step + o * across + k * step);
                const ar = util.at(F, a, e).*;
                const ai = util.at(F, a, e + 1).*;
                const xr = util.at(F, x, 2 * (s + k) * inc_x).*;
                const xi = util.at(F, x, 2 * (s + k) * inc_x + 1).*;
                c[0] += ar * xr;
                if (conj) c[1] -= ai * xr else c[1] += ai * xr;
                c[2] += ar * xi;
                c[3] += ai * xi;
            }
            const re = if (conj == xconj) c[0] - c[3] else c[0] + c[3];
            const im = if (xconj) c[1] - c[2] else c[1] + c[2];
            util.at(F, y, 2 * o * inc_y).* += re * alpha_r - im * alpha_i;
            util.at(F, y, 2 * o * inc_y + 1).* += im * alpha_r + re * alpha_i;
        }
    }
}

/// The `n`-form kernels over `XDouble` with `p` columns a block.
fn ExtendedBlocked(comptime p: BlasLong) type {
    const F = XDouble;
    return struct {
        fn q(m: BlasLong, n: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
            extendedGemv(F, false, p, m, n, alpha, a, lda, x, inc_x, y, inc_y);
            return 0;
        }
        fn Z(comptime conj: bool, comptime xconj: bool) type {
            return struct {
                fn k(m: BlasLong, n: BlasLong, _: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
                    extendedZgemv(F, false, conj, xconj, p, m, n, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y);
                    return 0;
                }
            };
        }
    };
}

/// On x86, points the extended `n`-form slots at kernels with `core`'s column
/// block (`extendedColBlock`); the portable ones have the unknown core's.
pub fn installExtendedCore(vt: *Vtable, core: X86CoreType) void {
    if (comptime !x86 or !types.isExtended(XDouble)) return;
    switch (core) {
        inline else => |c| {
            vt.qgemv_n = &ExtendedBlocked(extendedColBlock(false, c)).q;
            const Z = ExtendedBlocked(extendedColBlock(true, c)).Z;
            inline for (.{ .{ "n", false, false }, .{ "r", true, false }, .{ "o", false, true }, .{ "s", true, true } }) |v| {
                @field(vt, "xgemv_" ++ v[0]) = &Z(v[1], v[2]).k;
            }
        },
    }
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
    bf16_kernel.installGenericGemv(vt);
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

const check = @import("../../_test/blas1/check.zig");
const sym = @import("../symmetric.zig");
const guard = @import("../../_test/blas2/guard.zig");
const zigstd = @import("../../../../zigstd.zig");
const fma = @import("../../fma.zig");
const x86cpu = @import("../../cpu/x86.zig");

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
/// every column block (`n` forms, `DTB_ENTRIES` and `extendedColBlock`); on the
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
        // NaN: no kernel lets what the buffer held reach `y` (`buffers.md`)
        const buf = guard.open(F, bytes, 32, .nan);
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

/// Room for `shifted`'s `A`, `x` and `y` at the largest shift, 64-byte
/// aligned so each shift is the same misalignment on every run.
var shift_mem: [3][(100 + 3) * 13 + 4]C(f64) align(64) = undefined;
const shift_shapes = [_][2]usize{ .{ 7, 5 }, .{ 33, 9 }, .{ 100, 13 } };

/// `kern` on the same `A`, `x` and `y`, each shifted by 0 to 3 floats, any
/// alignment a float may have: the same bits in `y` each time, as
/// `threading.md`'s reproducibility needs.
fn shifted(comptime T: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    const w = if (complex) 2 else 1;
    const alpha: T = if (complex) .init(1.25, -0.5) else 1.25;
    var vals: [3][(100 + 3) * 13]T = undefined;
    var first: [100]T = undefined;
    for (shift_shapes) |mn| {
        const m = mn[0];
        const n = mn[1];
        const lda = m + 3;
        const rows = if (trans) n else m;
        const inner = if (trans) m else n;
        const lens = [3]usize{ lda * n, inner, rows };
        for (0..3) |k| check.fill(T, vals[k][0..lens[k]], m + 5 * k);
        for (0..64) |o| {
            var p: [3][*]F = undefined;
            for (0..3) |k| {
                p[k] = @as([*]F, @ptrCast(&shift_mem[k])) + (o >> @intCast(2 * k)) % 4;
                @memcpy(p[k][0 .. w * lens[k]], @as([*]const F, @ptrCast(&vals[k]))[0 .. w * lens[k]]);
            }
            const bytes = sym.gemvLen(F, w, m, n) * @sizeOf(F);
            const buf = guard.open(F, bytes, 32, .nan);
            if (complex) {
                _ = kern(@intCast(m), @intCast(n), 0, alpha.re, alpha.im, p[0], @intCast(lda), p[1], 1, p[2], 1, buf);
            } else {
                _ = kern(@intCast(m), @intCast(n), 0, alpha, p[0], @intCast(lda), p[1], 1, p[2], 1, buf);
            }
            const got = @as([*]align(1) T, @ptrCast(p[2]))[0..rows];
            if (o == 0) {
                @memcpy(first[0..rows], got);
                continue;
            }
            for (0..rows) |r| if (!check.same(T, first[r], got[r], true)) {
                check.say("{s} m={d} n={d} trans={} conj={} xconj={}, A, x, y shifted by {d},{d},{d}: y[{d}] {any}, unshifted {any}\n", .{ @typeName(T), m, n, trans, conj, xconj, o % 4, (o >> 2) % 4, o >> 4, r, got[r], first[r] });
                return error.TestUnexpectedResult;
            };
        }
    }
}

/// Every gemv kernel, the first failure of each named.
fn checkShifted(vt: *const Vtable) !void {
    var failed: ?anyerror = null;
    shifted(f32, false, false, false, vt.sgemv_n.?) catch |e| {
        failed = e;
    };
    shifted(f32, true, false, false, vt.sgemv_t.?) catch |e| {
        failed = e;
    };
    shifted(f64, false, false, false, vt.dgemv_n.?) catch |e| {
        failed = e;
    };
    shifted(f64, true, false, false, vt.dgemv_t.?) catch |e| {
        failed = e;
    };
    inline for (.{ .{ "c", f32 }, .{ "z", f64 } }) |p| {
        inline for (variants) |v| {
            shifted(C(p[1]), v[1], v[2], v[3], @field(vt, p[0] ++ "gemv_" ++ v[0]).?) catch |e| {
                failed = e;
            };
        }
    }
    if (failed) |e| return e;
}

test "portable gemv kernels: the same bits whatever the address" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkShifted(&vt);
}

test "every linked core's gemv kernels: the same bits whatever the address" {
    try check.eachCore(checkShifted);
}

// the half-precision gemvs, `ext/bf16/kernel.zig`'s on every core

const half_sizes = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 }, .{ 3, 2 }, .{ 4, 4 }, .{ 7, 5 }, .{ 8, 8 }, .{ 9, 16 }, .{ 16, 9 }, .{ 17, 17 }, .{ 33, 4 }, .{ 4, 33 } };
const half_strides = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 }, .{ 3, 2 } };
const half_alphas = [_]f32{ 1.25, -0.5 };
const half_betas = [_]f32{ 0, 1, 0.5 };
const half_cap = 33 * 3 + 2; // a vector's span at stride 3

fn fillAs(comptime T: type, buf: []T, seed: u64) void {
    var f: [400]f32 = undefined;
    check.fill(f32, f[0..buf.len], seed);
    for (buf, f[0..buf.len]) |*d, s| d.* = bf16_kernel.narrow(T, s);
}

fn checkHalf(comptime E: type, comptime Y: type, comptime trans: bool, kern: anytype) !void {
    for (half_sizes) |mn| for (half_strides) |inc| for (half_alphas) |alpha| for (half_betas) |beta| {
        const m = mn[0];
        const n = mn[1];
        const lda = m + 1;
        var a: [400]E = undefined; // lda * n at most 18 * 17
        var x: [half_cap]E = undefined;
        var y: [half_cap]Y = undefined;
        fillAs(E, &a, m * 7 + n);
        fillAs(E, &x, m + 3);
        fillAs(Y, &y, n + 5);
        const rows = if (trans) n else m;
        const inner = if (trans) m else n;
        // With beta 0, y must never be read: poison it.
        if (beta == 0) for (&y) |*v| {
            v.* = bf16_kernel.narrow(Y, std.math.nan(f32));
        };
        const orig = y;
        _ = kern(@intCast(m), @intCast(n), bf16_kernel.narrow(Y, alpha), &a, @intCast(lda), &x, @intCast(inc[0]), bf16_kernel.narrow(Y, beta), &y, @intCast(inc[1]));
        var touched: [half_cap]bool = @splat(false);
        for (0..rows) |r| {
            var acc: f64 = 0;
            for (0..inner) |q| {
                const e = if (trans) a[q + r * lda] else a[r + q * lda];
                acc += @as(f64, bf16_kernel.wide(E, e)) * bf16_kernel.wide(E, x[q * inc[0]]);
            }
            const yi = r * inc[1];
            touched[yi] = true;
            const prior: f64 = if (beta == 0) 0 else bf16_kernel.wide(Y, orig[yi]);
            const want = @as(f64, alpha) * acc + @as(f64, beta) * prior;
            const got: f64 = bf16_kernel.wide(Y, y[yi]);
            const tol: f64 = if (Y == BF16)
                0.005 * @abs(want) + 1e-4
            else
                check.tolerance(f32, 4 * inner + 4);
            if (!(@abs(want - got) <= tol)) {
                check.say("got {e}, want {e} within {e}: m={d} n={d} trans={} alpha={d} beta={d} row={d}: ", .{ got, want, tol, m, n, trans, alpha, beta, r });
                return check.case(error.TestExpectedApproxEqAbs, "half gemv", E, rows, inc);
            }
        }
        for (0..half_cap) |i| {
            if (touched[i]) continue;
            if (!std.mem.eql(u8, std.mem.asBytes(&y[i]), std.mem.asBytes(&orig[i]))) {
                check.say("untouched y[{d}] changed: m={d} n={d} trans={}: ", .{ i, m, n, trans });
                return check.case(error.TestUnexpectedResult, "half gemv", E, rows, inc);
            }
        }
    };
}

fn checkAllHalf(vt: *const Vtable) !void {
    inline for (bf16_kernel.families) |f| {
        try checkHalf(f[1], f[2], false, @field(vt, f[0] ++ "_n").?);
        try checkHalf(f[1], f[2], true, @field(vt, f[0] ++ "_t").?);
    }
}

test "half-precision gemv: the portable kernels match plain arithmetic" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAllHalf(&vt);
}

test "half-precision gemv: every linked core's kernels match plain arithmetic, where this machine runs them" {
    try check.eachCore(checkAllHalf);
}

/// `x` rounded to `f32`: each step of the models is, as on any IEEE unit, so
/// the x87 gives the same bits wherever values spill.
inline fn r32(x: f32) f32 {
    return util.rounded(f32, x);
}

/// `y` after the kernel, by a loop of a different shape (columns outermost for
/// `n`, every row's accumulator kept in an array), bit for bit.
fn halfModel(comptime E: type, comptime Y: type, comptime trans: bool, comptime fuse: bool, m: usize, n: usize, alpha: Y, a: []const E, lda: usize, x: []const E, incx: usize, beta: Y, y: []Y, incy: usize) void {
    const rows = if (trans) n else m;
    const inner = if (trans) m else n;
    var acc: [40]f32 = @splat(0);
    for (0..inner) |q| for (0..rows) |r| {
        const av = bf16_kernel.wide(E, if (trans) a[q + r * lda] else a[r + q * lda]);
        const xv = bf16_kernel.wide(E, x[q * incx]);
        acc[r] = if (fuse) @mulAdd(f32, av, xv, acc[r]) else r32(r32(av * xv) + acc[r]);
    };
    for (0..rows) |r| {
        const af = bf16_kernel.wide(Y, alpha);
        const bf = bf16_kernel.wide(Y, beta);
        const yi = r * incy;
        y[yi] = bf16_kernel.narrow(Y, if (bf == 0) r32(af * acc[r]) else if (fuse) @mulAdd(f32, af, acc[r], bf * bf16_kernel.wide(Y, y[yi])) else r32(r32(af * acc[r]) + r32(bf * bf16_kernel.wide(Y, y[yi]))));
    }
}

fn checkHalfBits(comptime E: type, comptime Y: type, comptime trans: bool, comptime fuse: bool) !void {
    return checkHalfBitsOf(E, Y, trans, fuse, bf16_kernel.Gemv(E, Y, trans, fuse).k);
}

fn checkHalfBitsOf(comptime E: type, comptime Y: type, comptime trans: bool, comptime fuse: bool, kern: anytype) !void {
    for (half_sizes) |mn| for (half_strides) |inc| for ([_]f32{ 1.25, -0.5, 0 }) |alpha| for (half_betas) |beta| {
        const m = mn[0];
        const n = mn[1];
        const lda = m + 2;
        var a: [400]E = undefined;
        var x: [half_cap]E = undefined;
        var y: [half_cap]Y = undefined;
        fillAs(E, &a, m * 5 + n);
        fillAs(E, &x, m + 9);
        fillAs(Y, &y, n + 1);
        var want = y;
        halfModel(E, Y, trans, fuse, m, n, bf16_kernel.narrow(Y, alpha), &a, lda, &x, inc[0], bf16_kernel.narrow(Y, beta), &want, inc[1]);
        _ = kern(@intCast(m), @intCast(n), bf16_kernel.narrow(Y, alpha), &a, @intCast(lda), &x, @intCast(inc[0]), bf16_kernel.narrow(Y, beta), &y, @intCast(inc[1]));
        if (!std.mem.eql(u8, std.mem.asBytes(&want), std.mem.asBytes(&y))) {
            check.say("m={d} n={d} trans={} fuse={} alpha={d} beta={d}: ", .{ m, n, trans, fuse, alpha, beta });
            return check.case(error.TestUnexpectedResult, "half gemv bits", E, m, inc);
        }
    };
}

test "half-precision gemv: the kernels are the model, bit for bit, fused or not" {
    inline for (bf16_kernel.families) |f| inline for (.{ false, true }) |trans| {
        try checkHalfBits(f[1], f[2], trans, false);
        if (comptime fma.fuseFor(f32)) try checkHalfBits(f[1], f[2], trans, true);
    };
}

test "half-precision gemv: the @Vector kernels are the model, bit for bit" {
    if (comptime !vector.ported.builds(.sbgemv_n)) return error.SkipZigTest;
    inline for (bf16_kernel.families) |f| {
        try checkHalfBitsOf(f[1], f[2], false, false, vector.ported.entry(@field(blas.Slot, f[0] ++ "_n")));
        try checkHalfBitsOf(f[1], f[2], true, false, vector.ported.entry(@field(blas.Slot, f[0] ++ "_t")));
    }
}

test "half-precision gemv: an empty operand leaves beta*y" {
    inline for (.{ false, true }) |trans| {
        const gen = bf16_kernel.Gemv(BF16, f32, trans, false).k;
        var a: [4]BF16 = @splat(bf16_kernel.narrow(BF16, 1));
        var x: [4]BF16 = @splat(bf16_kernel.narrow(BF16, 1));
        var y0: [4]f32 = .{ 2, 4, 6, 8 };
        _ = gen(if (trans) 0 else 4, if (trans) 4 else 0, 3, &a, 4, &x, 1, 0.5, &y0, 1);
        try std.testing.expectEqualSlices(f32, &.{ 1, 2, 3, 4 }, &y0);
    }
}

test "half-precision gemv: the row blocks of the non-transposed kernel sum each row as a row-by-row loop does" {
    const al = std.testing.allocator;
    inline for (bf16_kernel.families) |f| {
        const E = f[1];
        const Y = f[2];
        for ([_][2]usize{ .{ 255, 3 }, .{ 256, 2 }, .{ 257, 5 }, .{ 600, 7 } }) |mn| for ([_][2]usize{ .{ 1, 1 }, .{ 3, 2 } }) |inc| {
            const m = mn[0];
            const n = mn[1];
            const lda = m + 1;
            const a = try al.alloc(E, lda * n);
            defer al.free(a);
            const x = try al.alloc(E, n * inc[0]);
            defer al.free(x);
            const y = try al.alloc(Y, m * inc[1]);
            defer al.free(y);
            const want = try al.alloc(Y, m * inc[1]);
            defer al.free(want);
            var prng = std.Random.DefaultPrng.init(m + n);
            const r = prng.random();
            for (a) |*v| v.* = bf16_kernel.narrow(E, r.float(f32) * 4 - 2);
            for (x) |*v| v.* = bf16_kernel.narrow(E, r.float(f32) * 4 - 2);
            for (y, want) |*v, *w| {
                v.* = bf16_kernel.narrow(Y, r.float(f32) * 4 - 2);
                w.* = v.*;
            }
            const alpha: f32 = 1.25;
            const beta: f32 = 0.5;
            for (0..m) |i| {
                var acc: f32 = 0;
                for (0..n) |j| acc = r32(r32(bf16_kernel.wide(E, a[i + j * lda]) * bf16_kernel.wide(E, x[j * inc[0]])) + acc);
                const w = &want[i * inc[1]];
                w.* = bf16_kernel.narrow(Y, r32(r32(bf16_kernel.wide(Y, bf16_kernel.narrow(Y, alpha)) * acc) + r32(bf16_kernel.wide(Y, bf16_kernel.narrow(Y, beta)) * bf16_kernel.wide(Y, w.*))));
            }
            _ = bf16_kernel.Gemv(E, Y, false, false).k(@intCast(m), @intCast(n), bf16_kernel.narrow(Y, alpha), a.ptr, @intCast(lda), x.ptr, @intCast(inc[0]), bf16_kernel.narrow(Y, beta), y.ptr, @intCast(inc[1]));
            try std.testing.expect(std.mem.eql(u8, std.mem.sliceAsBytes(want), std.mem.sliceAsBytes(y)));
        };
    }
}

// the extended types' x87 sums

/// The same sums written out plainly: by element of `y`, then by block.
const PlainExtended = struct {
    fn gemv(comptime F: type, trans: bool, p: usize, m: usize, n: usize, alpha: F, a: []const F, lda: usize, x: []const F, ix: usize, y: []F, iy: usize) void {
        const outs = if (trans) n else m;
        const inner = if (trans) m else n;
        for (0..outs) |o| {
            var b: usize = 0;
            while (b < inner) : (b += p) {
                const len = @min(p, inner - b);
                const seq = o < outs / 4 * 4;
                var acc: [4]F = @splat(0);
                for (0..len) |k| {
                    const e = if (trans) a[b + k + o * lda] else a[o + (b + k) * lda];
                    acc[if (seq or k >= len / 8 * 8) 0 else k % 4] += e * x[(b + k) * ix];
                }
                const sum = if (seq) acc[0] else ((acc[0] + acc[1]) + acc[2]) + acc[3];
                y[o * iy] += sum * alpha;
            }
        }
    }

    fn zgemv(comptime F: type, trans: bool, conj: bool, xconj: bool, p: usize, m: usize, n: usize, alpha: [2]F, a: []const F, lda: usize, x: []const F, ix: usize, y: []F, iy: usize) void {
        const outs = if (trans) n else m;
        const inner = if (trans) m else n;
        for (0..outs) |o| {
            var b: usize = 0;
            while (b < inner) : (b += p) {
                var c: [4]F = @splat(0);
                for (b..@min(b + p, inner)) |i| {
                    const e = 2 * (if (trans) i + o * lda else o + i * lda);
                    // accumulator 2q + r takes part r of `A` times part q of `x`
                    for (0..2) |q| for (0..2) |r| {
                        const t = a[e + r] * x[2 * i * ix + q];
                        if (conj and r == 1 and q == 0) c[2 * q + r] -= t else c[2 * q + r] += t;
                    };
                }
                const re = if (conj == xconj) c[0] - c[3] else c[0] + c[3];
                const im = if (xconj) c[1] - c[2] else c[1] + c[2];
                y[2 * o * iy] += re * alpha[0] - im * alpha[1];
                y[2 * o * iy + 1] += im * alpha[0] + re * alpha[1];
            }
        }
    }
};

fn expectExtendedBits(comptime F: type, want: []const F, got: []const F) !void {
    const U = @Int(.unsigned, @bitSizeOf(F));
    for (want, got, 0..) |w, g, i| if (@as(U, @bitCast(w)) != @as(U, @bitCast(g))) {
        check.say("y[{d}]: got {e}, want {e}\n", .{ i, g, w });
        return error.TestExpectedEqual;
    };
}

const extended_variants = .{
    .{ false, false }, .{ true, false }, .{ false, true }, .{ true, true },
};

/// One shape at one block size, real and every complex variant, against `PlainExtended`.
fn checkExtendedShape(comptime F: type, comptime trans: bool, comptime p: BlasLong, m: usize, n: usize, inc: [2]usize) !void {
    const al = std.testing.allocator;
    const lda = m + 1;
    const outs = if (trans) n else m;
    const inner = if (trans) m else n;
    const a = try al.alloc(F, 2 * lda * n);
    defer al.free(a);
    const x = try al.alloc(F, 2 * ((inner - 1) * inc[0] + 1));
    defer al.free(x);
    const y = try al.alloc(F, 2 * ((outs - 1) * inc[1] + 1));
    defer al.free(y);
    const want = try al.alloc(F, y.len);
    defer al.free(want);
    check.fill(F, a, m * 7 + n);
    check.fill(F, x, m + 3);
    const alpha: F = 1.25;
    const ml: BlasLong = @intCast(m);
    const nl: BlasLong = @intCast(n);
    const ll: BlasLong = @intCast(lda);
    const ix: BlasLong = @intCast(inc[0]);
    const iy: BlasLong = @intCast(inc[1]);
    errdefer check.say("{s} trans={} p={d} m={d} n={d} inc={d},{d}: ", .{ @typeName(F), trans, p, m, n, inc[0], inc[1] });

    check.fill(F, y, n + 5);
    @memcpy(want, y);
    PlainExtended.gemv(F, trans, @intCast(p), m, n, alpha, a, lda, x, inc[0], want, inc[1]);
    extendedGemv(F, trans, p, ml, nl, alpha, a.ptr, ll, x.ptr, ix, y.ptr, iy);
    try expectExtendedBits(F, want, y);

    inline for (extended_variants) |v| {
        errdefer check.say("conj={} xconj={} ", .{ v[0], v[1] });
        check.fill(F, y, n + 6);
        @memcpy(want, y);
        PlainExtended.zgemv(F, trans, v[0], v[1], @intCast(p), m, n, .{ alpha, -0.5 }, a, lda, x, inc[0], want, inc[1]);
        extendedZgemv(F, trans, v[0], v[1], p, ml, nl, alpha, -0.5, a.ptr, ll, x.ptr, ix, y.ptr, iy);
        try expectExtendedBits(F, want, y);
    }
}

test "extended gemv: the sums follow the assembly's order, bit for bit" {
    const shapes = [_][2]usize{
        .{ 1, 1 },   .{ 2, 3 },   .{ 3, 2 },   .{ 4, 4 },   .{ 5, 7 }, .{ 7, 5 }, .{ 8, 9 }, .{ 9, 8 },
        .{ 17, 13 }, .{ 13, 17 }, .{ 37, 70 }, .{ 70, 37 },
    };
    inline for (.{ f80, f128 }) |F| {
        for (shapes) |s| for ([_][2]usize{ .{ 1, 1 }, .{ 2, 3 } }) |inc| {
            inline for (.{ comptime extendedRowBlock(false), comptime extendedRowBlock(true), 20 }) |p| try checkExtendedShape(F, true, p, s[0], s[1], inc);
            inline for (.{ comptime extendedColBlock(false, .unknown), comptime extendedColBlock(true, .unknown), 20 }) |p| try checkExtendedShape(F, false, p, s[0], s[1], inc);
        };
        // past one block of rows, the transposed kernels' blocking
        try checkExtendedShape(F, true, extendedRowBlock(false), 4100, 6, .{ 2, 3 });
        try checkExtendedShape(F, true, extendedRowBlock(true), 4100, 6, .{ 1, 1 });
        try checkExtendedShape(F, false, extendedColBlock(false, .unknown), 6, 4100, .{ 2, 3 });
    }
}

test "extended gemv: the sums at f80 carry the full significand" {
    // 1 + 2^-63 squared is 1 + 2^-62 in f80 and 1 in f64
    const e = 1 + std.math.ldexp(@as(f80, 1), -63);
    const want = 1 + std.math.ldexp(@as(f80, 1), -62);
    var a = [_]f80{ e, 0 };
    var y = [_]f80{ 0, 0 };
    extendedGemv(f80, true, 4096, 1, 1, 1, &a, 1, &a, 1, &y, 1);
    try zigstd.testExpectEqual(want, y[0]);
    y[0] = 0;
    extendedZgemv(f80, false, false, false, 32, 1, 1, 1, 0, &a, 1, &a, 1, &y, 1);
    try zigstd.testExpectEqual(want, y[0]);
}

test "extended gemv: the sums give the values traced by hand from the .S files" {
    // 2^64 + 1 rounds to 2^64 in f80 (ties to even), so each value below
    // depends on which accumulator a term lands in, on the fold order and on
    // where blocks end. Traced once from kernel/x86_64/{q,x}gemv_{n,t}.S;
    // the plain loop of the portable gemv gives the value in parentheses.
    const F = f80;
    const t = std.math.ldexp(@as(F, 1), 64);
    var ones: [16]F = @splat(1);

    // qgemv_t: columns in fours, each one sum (`faddp` into st(1)..st(4)); the
    // n % 4 left go into four sums by row % 4, eight a trip (.L65), the tail
    // into the first (.L71), folded ((a0 + a1) + a2) + a3 (.L70). Each column
    // is t, 1, 0, 0, 0, 1, 0, 0, 1 and y starts at -t: one sum gives t, so 0;
    // column 4's a0 = t + 1 = t, a1 = 2 give t + 2, so 2 (0). A tail into a1
    // would make a1 = 3 and give t + 4.
    const col = [9]F{ t, 1, 0, 0, 0, 1, 0, 0, 1 };
    var at: [45]F = undefined;
    for (0..5) |j| @memcpy(at[9 * j ..][0..9], &col);
    var y: [5]F = @splat(-t);
    extendedGemv(F, true, extendedRowBlock(false), 9, 5, 1, &at, 9, &ones, 1, &y, 1);
    try zigstd.testExpectEqualSlices(F, &.{ 0, 0, 0, 0, 2 }, &y);

    // qgemv_n: rows in fours, each one sum over the columns; the m % 4 left as
    // qgemv_t's leftover columns. Each row is `col`: 0 for rows 0 to 3, 2 for
    // row 4. The portable `y += (alpha*x)*a` gives 3 for every row.
    var an: [45]F = undefined;
    for (0..9) |j| an[5 * j ..][0..5].* = @splat(col[j]);
    y = @splat(-t);
    extendedGemv(F, false, extendedColBlock(false, .unknown), 5, 9, 1, &an, 5, &ones, 1, &y, 1);
    try zigstd.testExpectEqualSlices(F, &.{ 0, 0, 0, 0, 2 }, &y);

    // xgemv_n and xgemv_t: four sums c1..c4 = ar*xr, ai*xr, ar*xi, ai*xi
    // (`CONJ`: c2 -= ai*xr, `fsubrp`), then re = c1 - c4, im = c2 + c3; `CONJ`
    // re = c1 + c4; `XCONJ` re = c1 + c4, im = c2 - c3; both, re = c1 - c4,
    // im = c2 - c3 (`fsubp` is fsubr in AT&T). A = (t, t), (1, 1) and
    // x = (1, 1), (1, 0): c1 = t + 1 = t, c2 = t (`CONJ`: -t - 1 = -t),
    // c3 = c4 = t. With alpha = 1, y = alpha * (re, im).
    var az = [_]F{ t, t, 1, 1 };
    var xz = [_]F{ 1, 1, 1, 0 };
    const want = [4][2]F{
        .{ 0, 2 * t }, // n, t (1, 2t)
        .{ 2 * t, 0 }, // r, c (2t, -1)
        .{ 2 * t, 0 }, // o, u (2t, 1)
        .{ 0, -2 * t }, // s, d (1, -2t)
    };
    inline for (extended_variants, want) |v, w| {
        var yz: [2]F = @splat(0);
        extendedZgemv(F, true, v[0], v[1], extendedRowBlock(true), 2, 1, 1, 0, &az, 2, &xz, 1, &yz, 1);
        try zigstd.testExpectEqualSlices(F, &w, &yz);
        yz = @splat(0);
        extendedZgemv(F, false, v[0], v[1], extendedColBlock(true, .unknown), 1, 2, 1, 0, &az, 1, &xz, 1, &yz, 1);
        try zigstd.testExpectEqualSlices(F, &w, &yz);
    }
}

test "extended gemv: the sums add each block to y as it ends" {
    // y starts at -t; t at the first index of the sum and 1 just past the
    // first block. ExtendedBlocked: -t + t, then + 1, so 1; one sum: t + 1 = t, so 0.
    const al = std.testing.allocator;
    const F = f80;
    const t = std.math.ldexp(@as(F, 1), 64);
    inline for (.{ false, true }) |complex| {
        const w = if (complex) 2 else 1;
        // the transposed kernels: rows in blocks, four columns
        {
            const p = comptime extendedRowBlock(complex);
            const m: usize = @intCast(p + 1);
            const a = try al.alloc(F, w * m * 4);
            defer al.free(a);
            const x = try al.alloc(F, w * m);
            defer al.free(x);
            @memset(a, 0);
            for (0..m) |i| {
                x[w * i] = 1;
                if (complex) x[w * i + 1] = 0;
            }
            for (0..4) |j| {
                a[w * j * m] = t;
                a[w * (j * m + m - 1)] = 1;
            }
            var y: [8]F = @splat(0);
            for (0..4) |j| y[w * j] = -t;
            if (complex) {
                extendedZgemv(F, true, false, false, p, @intCast(m), 4, 1, 0, a.ptr, @intCast(m), x.ptr, 1, &y, 1);
                try zigstd.testExpectEqualSlices(F, &.{ 1, 0, 1, 0, 1, 0, 1, 0 }, &y);
            } else {
                extendedGemv(F, true, p, @intCast(m), 4, 1, a.ptr, @intCast(m), x.ptr, 1, &y, 1);
                try zigstd.testExpectEqualSlices(F, &.{ 1, 1, 1, 1 }, y[0..4]);
            }
        }
        // the plain kernels: columns in blocks, four rows
        {
            const p = comptime extendedColBlock(complex, .unknown);
            const n: usize = @intCast(p + 1);
            const a = try al.alloc(F, w * 4 * n);
            defer al.free(a);
            const x = try al.alloc(F, w * n);
            defer al.free(x);
            @memset(a, 0);
            for (0..n) |j| {
                x[w * j] = 1;
                if (complex) x[w * j + 1] = 0;
            }
            for (0..4) |i| {
                a[w * i] = t;
                a[w * (i + 4 * (n - 1))] = 1;
            }
            var y: [8]F = @splat(0);
            for (0..4) |i| y[w * i] = -t;
            if (complex) {
                extendedZgemv(F, false, false, false, p, 4, @intCast(n), 1, 0, a.ptr, 4, x.ptr, 1, &y, 1);
                try zigstd.testExpectEqualSlices(F, &.{ 1, 0, 1, 0, 1, 0, 1, 0 }, &y);
            } else {
                extendedGemv(F, false, p, 4, @intCast(n), 1, a.ptr, 4, x.ptr, 1, &y, 1);
                try zigstd.testExpectEqualSlices(F, &.{ 1, 1, 1, 1 }, y[0..4]);
            }
        }
    }
}

test "extended gemv: the blocks by target and core" {
    if (x86) {
        try zigstd.testExpectEqual(@as(BlasLong, 1000), extendedRowBlock(false));
        try zigstd.testExpectEqual(@as(BlasLong, 400), extendedRowBlock(true));
        try zigstd.testExpectEqual(@as(BlasLong, 256), extendedColBlock(false, .core2));
        try zigstd.testExpectEqual(@as(BlasLong, 256), extendedColBlock(true, .core2));
        try zigstd.testExpectEqual(@as(BlasLong, 64), extendedColBlock(false, .northwood));
        try zigstd.testExpectEqual(@as(BlasLong, 32), extendedColBlock(true, .northwood));
        try zigstd.testExpectEqual(@as(BlasLong, 32), extendedColBlock(false, .athlon));
        try zigstd.testExpectEqual(@as(BlasLong, 16), extendedColBlock(true, .athlon));
        try zigstd.testExpectEqual(@as(BlasLong, 32), extendedColBlock(true, .opteron));
    } else {
        try zigstd.testExpectEqual(@as(BlasLong, 4096), extendedRowBlock(false));
        try zigstd.testExpectEqual(@as(BlasLong, 4096), extendedRowBlock(true));
        try zigstd.testExpectEqual(@as(BlasLong, 32), extendedColBlock(false, .core2));
        try zigstd.testExpectEqual(@as(BlasLong, 32), extendedColBlock(true, .athlon));
    }
}

/// x86's `qgemv_n.S` and `xgemv_n.S` columns a block (`P`) by core:
/// `getarch.c`'s `DTB_DEFAULT_ENTRIES`, which `xgemv_n.S` halves for
/// `PENTIUM4` (northwood, prescott) and `ATHLON`; `qgemv_n.S` takes 32 on
/// `ATHLON` and `OPTERON`, their own anyway.
const x87_col_blocks = [_]struct { core: X86CoreType, q: usize, x: usize }{
    .{ .core = .katmai, .q = 64, .x = 64 },
    .{ .core = .northwood, .q = 64, .x = 32 },
    .{ .core = .prescott, .q = 64, .x = 32 },
    .{ .core = .athlon, .q = 32, .x = 16 },
    .{ .core = .opteron, .q = 32, .x = 32 },
    .{ .core = .barcelona, .q = 48, .x = 48 },
    .{ .core = .bobcat, .q = 40, .x = 40 },
    .{ .core = .core2, .q = 256, .x = 256 },
    .{ .core = .nano, .q = 128, .x = 128 },
};

/// A core's vtable, static: a vtable is large for a rig's stack.
var x87_core_vt: Vtable = undefined;

/// `y := -t + A*x` on four rows, `t` = 2^64 in column 0 and 1 in column
/// `at`, the rest 0, `x` ones: a block ending before `at` gives 1, one
/// holding both gives `t + 1 = t`, so 0.
fn x87BlockProbe(comptime complex: bool, at: usize, n: usize) !XDouble {
    const al = std.testing.allocator;
    const F = XDouble;
    const w = if (complex) 2 else 1;
    const t = std.math.ldexp(@as(F, 1), 64);
    const a = try al.alloc(F, w * 4 * n);
    defer al.free(a);
    const x = try al.alloc(F, w * n);
    defer al.free(x);
    @memset(a, 0);
    @memset(x, 0);
    for (0..n) |j| x[w * j] = 1;
    for (0..4) |i| {
        a[w * i] = t;
        a[w * (i + 4 * at)] = 1;
    }
    var y: [8]F = @splat(0);
    for (0..4) |i| y[w * i] = -t;
    const ln: BlasLong = @intCast(n);
    if (complex) {
        _ = x87_core_vt.xgemv_n.?(4, ln, 0, 1, 0, a.ptr, 4, x.ptr, 1, &y, 1, undefined);
    } else {
        _ = x87_core_vt.qgemv_n.?(4, ln, 0, 1, a.ptr, 4, x.ptr, 1, &y, 1, undefined);
    }
    for (1..4) |i| try zigstd.testExpectEqual(y[0], y[w * i]);
    return y[0];
}

test "extended gemv: each x86 core's qgemv_n and xgemv_n slots take its column block" {
    if (comptime !x86 or !types.isExtended(XDouble)) return error.SkipZigTest;
    for (x87_col_blocks) |c| {
        x87_core_vt = x86cpu.vtableFromCoreType(c.core);
        inline for (.{ false, true }) |complex| {
            const p = if (complex) c.x else c.q;
            errdefer check.say("{s} {s}gemv_n, block {d}: ", .{ @tagName(c.core), if (complex) "x" else "q", p });
            // the 1 last in the first block, then first in the second
            try zigstd.testExpectEqual(@as(XDouble, 0), try x87BlockProbe(complex, p - 1, 2 * p + 1));
            try zigstd.testExpectEqual(@as(XDouble, 1), try x87BlockProbe(complex, p, 2 * p + 1));
        }
    }
}

test "x86_64's dgemv_n vector block and SSE2 dgemv_t blocks give the reference's arithmetic, bit for bit" {
    if (sysinfo.target.cpu.arch != .x86_64) return error.SkipZigTest;
    const driver = @import("kernel.zig");
    const xk = @import("x86_64/kernel.zig");
    const own = @import("x86_64/uk/dgemv_t_4.zig");
    var a: [4 * 40]f64 = undefined;
    var x: [40]f64 = undefined;
    check.fill(f64, &a, 11);
    check.fill(f64, &x, 12);
    const ap = [4][*]f64{ a[0..].ptr, a[40..].ptr, a[80..].ptr, a[120..].ptr };
    const alpha = [1]f64{1.25};
    // Skylake-X's `dgemv_kernel_4x4` and `4x2` against a row-by-row model
    inline for (.{ false, true }) |wide| inline for (.{ false, true }) |fuse| {
        const B = xk.VectorBlockWith(.dgemv_n, wide, fuse);
        var n: usize = 4;
        while (n <= 40) : (n += 4) {
            inline for (.{ 4, 2 }) |cols| {
                var y: [40]f64 = undefined;
                check.fill(f64, &y, n);
                var want = y;
                for (0..n) |i| {
                    var sum = fma.madd(fuse, f64, ap[0][i], x[0], ap[1][i] * x[1]);
                    for (2..cols) |c| sum = fma.madd(fuse, f64, ap[c][i], x[c], sum);
                    want[i] = fma.madd(fuse, f64, sum, alpha[0], want[i]);
                }
                if (cols == 4) B.dgemv_kernel_4x4(@intCast(n), &ap, &x, &y, &alpha) else B.dgemv_kernel_4x2(@intCast(n), &ap, &x, &y, &alpha);
                try std.testing.expectEqualSlices(f64, &want, &y);
            }
        }
    };
    // the SSE2 forms of `dgemv_t_4.c`'s `kernel_4x2` and `4x1` against its assembly
    if (comptime !std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .sse3)) return;
    const S = driver.Fallback(f64, fma.fuse);
    var n: usize = 2;
    while (n <= 40) : (n += 2) {
        var want: [2]f64 = undefined;
        var got: [2]f64 = undefined;
        own.dgemv_kernel_4x2(@intCast(n), ap[0], ap[1], &x, &want);
        S.t4x2(@intCast(n), ap[0], ap[1], &x, &got);
        try std.testing.expectEqualSlices(f64, &want, &got);
        own.dgemv_kernel_4x1(@intCast(n), ap[2], &x, &want);
        S.t4x1(@intCast(n), ap[2], &x, &got);
        try std.testing.expectEqual(want[0], got[0]);
    }
}

test "x86_64's Skylake-X sgemv (Haswell's blocks and its own) matches plain arithmetic, where this build has AVX2 and FMA" {
    if (sysinfo.target.cpu.arch != .x86_64) return error.SkipZigTest;
    const f = sysinfo.target.cpu.features;
    if (comptime !std.Target.x86.featureSetHasAll(f, .{ .avx2, .fma })) return error.SkipZigTest;
    const S = @import("x86_64/kernel.zig").Skylake;
    // `n`: the few-column path, the column-block path with its 2- and
    // 1-column tails (`n` past 48), and 4 columns, which takes the second
    try checkGemv(f32, false, false, false, S.entry(.sgemv_n));
    try checkGemvAt(f32, false, false, false, S.entry(.sgemv_n), &.{ .{ 37, 51 }, .{ 70, 55 }, .{ 20, 4 }, .{ 129, 49 } }, &strides, std.testing.allocator);
    // `t`: `checkGemv`'s `lda = m + 1` takes the column-block path, and
    // `lda == m` every shape of up to 8 rows across its blocks
    try checkGemv(f32, true, false, false, S.entry(.sgemv_t));
    try checkSkylakeT(S.entry(.sgemv_t));
}

/// `kern` with `lda == m` against plain arithmetic, for `m` up to 8.
fn checkSkylakeT(kern: anytype) !void {
    var a: [8 * 70]f32 = undefined;
    var x: [16]f32 = undefined;
    var y: [210]f32 = undefined;
    for (1..9) |m| for ([_]usize{ 1, 2, 3, 5, 7, 9, 15, 16, 17, 31, 33, 35, 47, 64, 69 }) |n| for (strides) |inc| {
        check.fill(f32, &a, m * 3 + n);
        check.fill(f32, &x, m);
        check.fill(f32, &y, n + 1);
        var want = y;
        for (0..n) |j| {
            var acc: f64 = 0;
            for (0..m) |i| acc += @as(f64, a[j * m + i]) * x[i * inc[0]];
            want[j * inc[1]] = @floatCast(@as(f64, want[j * inc[1]]) + 1.25 * acc);
        }
        var buf: [160]f32 = @splat(std.math.nan(f32));
        _ = kern(@intCast(m), @intCast(n), 0, 1.25, &a, @intCast(m), &x, @intCast(inc[0]), &y, @intCast(inc[1]), &buf);
        check.expectCloseSlices(f32, &want, &y, check.tolerance(f32, 4 * m + 4)) catch |e| {
            check.say("sgemv_t m={d} n={d}: ", .{ m, n });
            return e;
        };
    };
}
