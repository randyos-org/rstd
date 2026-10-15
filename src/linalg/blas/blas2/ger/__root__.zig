//! `GER` kernels: the rank-1 update `A += alpha*x*y^T` (conjugated for the
//! complex `c`, `v`, `d` forms). The vtable wrappers, the portable kernels,
//! and tests of every core's kernels against plain arithmetic.
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

pub fn sbger_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, a: [*]f32, lda: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.sbger_k)(m, n, dummy1, alpha, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn sger_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, a: [*]f32, lda: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.sger_k)(m, n, dummy1, alpha, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn dger_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, a: [*]f64, lda: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.dger_k)(m, n, dummy1, alpha, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn qger_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, a: [*]XDouble, lda: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.qger_k)(m, n, dummy1, alpha, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn cgeru_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, a: [*]f32, lda: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgeru_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn cgerc_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, a: [*]f32, lda: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgerc_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn cgerv_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, a: [*]f32, lda: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgerv_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn cgerd_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f32, alpha_i: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, a: [*]f32, lda: BlasLong, buffer: [*]f32) i32 {
    return blas.kernel(.cgerd_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn zgeru_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, a: [*]f64, lda: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgeru_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn zgerc_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, a: [*]f64, lda: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgerc_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn zgerv_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, a: [*]f64, lda: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgerv_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn zgerd_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: f64, alpha_i: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, a: [*]f64, lda: BlasLong, buffer: [*]f64) i32 {
    return blas.kernel(.zgerd_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn xgeru_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, a: [*]XDouble, lda: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgeru_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn xgerc_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, a: [*]XDouble, lda: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgerc_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn xgerv_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, a: [*]XDouble, lda: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgerv_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}
pub fn xgerd_k(m: BlasLong, n: BlasLong, dummy1: BlasLong, alpha_r: XDouble, alpha_i: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, a: [*]XDouble, lda: BlasLong, buffer: [*]XDouble) i32 {
    return blas.kernel(.xgerd_k)(m, n, dummy1, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../../blas1/check.zig");

/// `kernel/generic/ger.c`: the rank-1 update `A += alpha*x*y^T`, `A`
/// column-major `m` by `n`.
fn Ger(comptime F: type) type {
    return struct {
        fn k(m: BlasLong, n: BlasLong, _: BlasLong, alpha: F, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, a: [*]F, lda: BlasLong, _: [*]F) callconv(.c) i32 {
            var j: BlasLong = 0;
            while (j < n) : (j += 1) {
                const t = alpha * util.at(F, y, j * inc_y).*;
                var i: BlasLong = 0;
                while (i < m) : (i += 1) util.at(F, a, i + j * lda).* += t * util.at(F, x, i * inc_x).*;
            }
            return 0;
        }
    };
}

/// `kernel/generic/zger.c`: `A += alpha*x*y^T` with `y` conjugated for `conj`
/// (`gerc`, `gerd`) and `x` for `xconj` (`gerv`, `gerd`).
fn Zger(comptime F: type, comptime conj: bool, comptime xconj: bool) type {
    return struct {
        fn get(p: [*]F, i: BlasLong, c: bool) std.math.Complex(F) {
            const im = util.at(F, p, 2 * i + 1).*;
            return .init(util.at(F, p, 2 * i).*, if (c) -im else im);
        }
        fn k(m: BlasLong, n: BlasLong, _: BlasLong, alpha_r: F, alpha_i: F, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, a: [*]F, lda: BlasLong, _: [*]F) callconv(.c) i32 {
            const alpha: std.math.Complex(F) = .init(alpha_r, alpha_i);
            var j: BlasLong = 0;
            while (j < n) : (j += 1) {
                const t = alpha.mul(get(y, j * inc_y, conj));
                var i: BlasLong = 0;
                while (i < m) : (i += 1) {
                    const v = t.mul(get(x, i * inc_x, xconj));
                    util.at(F, a, 2 * (i + j * lda)).* += v.re;
                    util.at(F, a, 2 * (i + j * lda) + 1).* += v.im;
                }
            }
            return 0;
        }
    };
}

/// The complex variants by suffix: `y` conjugated, `x` conjugated.
const variants = .{ .{ "u", false, false }, .{ "c", true, false }, .{ "v", false, true }, .{ "d", true, true } };

/// Points every ger slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.sger_k = &Ger(f32).k;
    vt.dger_k = &Ger(f64).k;
    vt.qger_k = &Ger(XDouble).k;
    inline for (.{ .{ "c", f32 }, .{ "z", f64 }, .{ "x", XDouble } }) |p| {
        inline for (variants) |v| @field(vt, p[0] ++ "ger" ++ v[0] ++ "_k") = &Zger(p[1], v[1], v[2]).k;
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
fn checkGer(comptime T: type, comptime conj: bool, comptime xconj: bool, kern: anytype) !void {
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
        var want = a;
        for (0..n) |j| for (0..m) |i| {
            const xv = x[i * inc[0]];
            const yv = y[j * inc[1]];
            if (complex) {
                const xe = if (xconj) xv.conjugate() else xv;
                const ye = if (conj) yv.conjugate() else yv;
                want[i + j * lda] = want[i + j * lda].add(alpha.mul(xe).mul(ye));
            } else want[i + j * lda] += alpha * xv * yv;
        };
        if (complex) {
            _ = kern(@intCast(m), @intCast(n), 0, alpha.re, alpha.im, @ptrCast(&x), @intCast(inc[0]), @ptrCast(&y), @intCast(inc[1]), @ptrCast(&a), @intCast(lda), check.work(F));
        } else {
            _ = kern(@intCast(m), @intCast(n), 0, alpha, &x, @intCast(inc[0]), &y, @intCast(inc[1]), &a, @intCast(lda), check.work(F));
        }
        check.expectCloseSlices(T, &want, &a, check.tolerance(F, 4)) catch |e| {
            check.say("m={d} n={d}: ", .{ m, n });
            return check.case(e, "ger", T, m, inc);
        };
    };
}

fn checkAll(vt: *const Vtable) !void {
    try checkGer(f32, false, false, vt.sger_k.?);
    try checkGer(f64, false, false, vt.dger_k.?);
    inline for (.{ .{ "c", f32 }, .{ "z", f64 } }) |p| {
        inline for (variants) |v| try checkGer(C(p[1]), v[1], v[2], @field(vt, p[0] ++ "ger" ++ v[0] ++ "_k").?);
    }
}

test "portable ger kernels match plain arithmetic" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkGer(XDouble, false, false, vt.qger_k.?);
    inline for (variants) |v| try checkGer(C(XDouble), v[1], v[2], @field(vt, "xger" ++ v[0] ++ "_k").?);
}

test "every linked core's ger kernels match plain arithmetic, where this machine runs them" {
    try check.eachCore(checkAll);
}
