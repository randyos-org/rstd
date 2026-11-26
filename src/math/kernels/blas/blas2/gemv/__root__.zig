//! `GEMV` kernels, `y += alpha * op(A) * x` on column-major `A`: the portable
//! ones (`*_generic`; the extended types' are `extended.zig`'s x87 sums) and
//! tests of every core's against plain arithmetic. The half-precision gemvs'
//! kernels are `ext/bf16/kernel.zig`'s `Gemv`, the same on every core.

const types = @import("../../types.zig");
const portable = @import("../../portable.zig");
const std = @import("std");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const bf16_kernel = @import("../../ext/bf16/kernel.zig");
const extended = @import("extended.zig");

const BlasLong = types.BlasLong;
const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

/// `kernel/arm/gemv_n.c` and `gemv_t.c`: `y += alpha*A*x` (`n`, `y` of `m`)
/// or `y += alpha*A^T*x` (`t`, `y` of `n`). The interface has already scaled
/// `y` by `beta`. Wider than `f64`, the x87 `qgemv_n.S`/`qgemv_t.S` sums
/// instead (`extendedGemv`).
pub fn Gemv(comptime F: type) type {
    return struct {
        pub fn n(m: BlasLong, cols: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
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
        pub fn t(m: BlasLong, cols: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
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
pub fn Zgemv(comptime F: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool) type {
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
        pub fn k(m: BlasLong, cols: BlasLong, _: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
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

// ---- extended types: `extended.zig` ----------------------------------------

const X86CoreType = extended.X86CoreType;
const extendedRowBlock = extended.extendedRowBlock;
const extendedColBlock = extended.extendedColBlock;
const extendedGemv = extended.extendedGemv;
const extendedZgemv = extended.extendedZgemv;
pub const installExtendedCore = extended.installExtendedCore;

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

/// The complex variants by op: transposed, `A` conjugated, `x`
/// conjugated, as `types.Trans` says.
pub const variants = blk: {
    const all = std.enums.values(types.Trans);
    var v: [all.len]struct { portable.Op, bool, bool, bool } = undefined;
    for (all, &v) |t, *e| e.* = .{ portable.gemvOp(t), t.transposed(), t.conjA(), t.conjX() };
    break :blk v;
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
            portable.slotPtr(vt, .from(p[0]), v[0]).* = &Zgemv(p[1], v[1], v[2], v[3]).k;
        }
    }
}

test {
    _ = extended;
}
