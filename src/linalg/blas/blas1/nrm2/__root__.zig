//! `NRM2`: the Euclidean norm of `x`. The vtable wrappers, the portable
//! kernels (`*_generic`), and tests of every core's kernels against a plain loop.
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

pub fn sbnrm2_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.sbnrm2_k)(n, x, inc_x);
}
pub fn snrm2_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.snrm2_k)(n, x, inc_x);
}
pub fn dnrm2_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas.kernel(.dnrm2_k)(n, x, inc_x);
}
pub fn qnrm2_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas.kernel(.qnrm2_k)(n, x, inc_x);
}
pub fn cnrm2_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.cnrm2_k)(n, x, inc_x);
}
pub fn znrm2_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas.kernel(.znrm2_k)(n, x, inc_x);
}
pub fn xnrm2_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas.kernel(.xnrm2_k)(n, x, inc_x);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../check.zig");
const x87 = @import("../x87.zig");

/// `kernel/arm/nrm2.c` and `znrm2.c`: the Euclidean norm, scaled by the
/// largest magnitude so far so the squares cannot overflow; an element is `w`
/// floats. Only a zero stride is refused: a negative one walks backwards.
/// Wider than `f64`, x86_64 `nrm2.S`'s and `znrm2.S`'s plain sum of squares
/// (`x87.zig`), unscaled: it overflows where they do.
fn Nrm2(comptime F: type, comptime w: comptime_int) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong) callconv(types.kernelCc(F)) F {
            if (n <= 0 or inc_x == 0) return 0;
            if (comptime types.isExtended(F)) return @sqrt(x87.sum(F, w, .square, n, x, inc_x));
            if (w == 1 and n == 1) return @abs(x[0]);
            var scale: types.Acc(F) = 0;
            var ssq: types.Acc(F) = 1;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                inline for (0..w) |c| {
                    const v = util.at(F, x, i * inc_x * w + c).*;
                    if (v != 0) {
                        const a: types.Acc(F) = @abs(v);
                        if (scale < a) {
                            ssq = 1 + ssq * (scale / a) * (scale / a);
                            scale = a;
                        } else {
                            ssq += (a / scale) * (a / scale);
                        }
                    }
                }
            }
            return @floatCast(scale * @sqrt(ssq));
        }
    };
}

pub const snrm2_generic = &Nrm2(f32, 1).k;
pub const dnrm2_generic = &Nrm2(f64, 1).k;
pub const qnrm2_generic = &Nrm2(XDouble, 1).k;
pub const cnrm2_generic = &Nrm2(f32, 2).k;
pub const znrm2_generic = &Nrm2(f64, 2).k;
pub const xnrm2_generic = &Nrm2(XDouble, 2).k;

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    const w = if (types.isComplex(T)) 2 else 1;
    return struct {
        pub const nrm2_k = &Nrm2(F, w).k;
    };
}

/// Points every nrm2 slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.snrm2_k = snrm2_generic;
    vt.dnrm2_k = dnrm2_generic;
    vt.qnrm2_k = qnrm2_generic;
    vt.cnrm2_k = cnrm2_generic;
    vt.znrm2_k = znrm2_generic;
    vt.xnrm2_k = xnrm2_generic;
}

// ---- tests -----------------------------------------------------------------

/// `kern` against a plain loop in `f64`; `T` is the element type.
fn checkNrm2(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    for (check.sizes) |n| {
        for ([_]usize{ 1, 2, 3 }) |inc| {
            var x: [check.capacity]T = undefined;
            check.fill(T, &x, n + 1);
            const W = check.Wide(F);
            var acc: W = 0;
            for (0..n) |i| {
                const v = x[i * inc];
                acc += if (comptime complex) @as(W, v.re) * v.re + @as(W, v.im) * v.im else @as(W, v) * v;
            }
            const want = @sqrt(acc);
            const got = kern(@intCast(n), @as([*]F, @ptrCast(&x)), @intCast(inc));
            check.expectCloseRef(F, want, got, check.tolerance(F, 2 * n) * @max(1, want)) catch |e| return check.case(e, "nrm2", T, n, .{ inc, 0 });
        }
    }
}

fn checkAll(vt: *const Vtable) !void {
    try checkNrm2(f32, vt.snrm2_k.?);
    try checkNrm2(f64, vt.dnrm2_k.?);
    try checkNrm2(std.math.Complex(f32), vt.cnrm2_k.?);
    try checkNrm2(std.math.Complex(f64), vt.znrm2_k.?);
}

test "portable nrm2 kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkNrm2(XDouble, vt.qnrm2_k.?);
    try checkNrm2(std.math.Complex(XDouble), vt.xnrm2_k.?);
}

test "portable nrm2 kernels at the extended type match a plain loop" {
    try checkNrm2(check.Ext, &Nrm2(check.Ext, 1).k);
    try checkNrm2(std.math.Complex(check.Ext), &Nrm2(check.Ext, 2).k);
}

test "every linked core's nrm2 kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}

test "nrm2 past 64 bits is the plain sum of squares, overflow included; at 64 bits or less, scaled" {
    // by width, so a 64-bit `c_longdouble` of its own takes the scaled path
    inline for (.{ f32, f64, c_longdouble, f80, f128 }) |F| {
        const big = if (comptime types.isExtended(F)) std.math.ldexp(@as(F, 1), 16000) else std.math.floatMax(F) / 2;
        var x = [_]F{ big, big };
        const got = Nrm2(F, 1).k(2, &x, 1);
        if (comptime types.isExtended(F)) {
            try std.testing.expect(std.math.isInf(got));
        } else {
            try std.testing.expectApproxEqRel(big * @sqrt(@as(F, 2)), got, 1e-6);
        }
    }
}
