//! `ASUM`: the sum of `|x|`. The vtable wrappers, the portable kernels
//! (`*_generic`), and tests of every core's kernels against a plain loop.
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

pub fn sbasum_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.sbasum_k)(n, x, inc_x);
}
pub fn sasum_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.sasum_k)(n, x, inc_x);
}
pub fn dasum_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas.kernel(.dasum_k)(n, x, inc_x);
}
pub fn qasum_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas.kernel(.qasum_k)(n, x, inc_x);
}
pub fn casum_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.casum_k)(n, x, inc_x);
}
pub fn zasum_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas.kernel(.zasum_k)(n, x, inc_x);
}
pub fn xasum_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas.kernel(.xasum_k)(n, x, inc_x);
}

// ---- portable kernels ------------------------------------------------------

const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

/// `kernel/arm/asum.c` and `zasum.c`: the sum of `|x|`, of `|re| + |im|` for
/// complex; an element is `w` floats. Wider than `f64`, x86_64 `asum.S`'s
/// and `zasum.S`'s x87 sum (`util.extendedSum`).
fn Asum(comptime F: type, comptime w: comptime_int) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong) callconv(types.kernelCc(F)) F {
            if (n <= 0 or inc_x <= 0) return 0;
            if (comptime types.isExtended(F)) return util.extendedSum(F, w, .abs, n, x, inc_x);
            var sum: types.Acc(F) = 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                inline for (0..w) |c| sum += @abs(util.at(F, x, i * inc_x * w + c).*);
            }
            return @floatCast(sum);
        }
    };
}

pub const sasum_generic = &Asum(f32, 1).k;
pub const dasum_generic = &Asum(f64, 1).k;
pub const qasum_generic = &Asum(XDouble, 1).k;
pub const casum_generic = &Asum(f32, 2).k;
pub const zasum_generic = &Asum(f64, 2).k;
pub const xasum_generic = &Asum(XDouble, 2).k;

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    const w = if (types.isComplex(T)) 2 else 1;
    return struct {
        pub const asum_k = &Asum(F, w).k;
    };
}

/// Points every asum slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.sasum_k = sasum_generic;
    vt.dasum_k = dasum_generic;
    vt.qasum_k = qasum_generic;
    vt.casum_k = casum_generic;
    vt.zasum_k = zasum_generic;
    vt.xasum_k = xasum_generic;
}

// ---- tests -----------------------------------------------------------------

const zigstd = @import("../../../../zigstd.zig");
const std = @import("std");
const check = @import("../../_test/blas1/check.zig");

/// `kern` against a plain loop in `f64`; `T` is the element type.
fn checkAsum(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    for (check.sizes) |n| {
        for ([_]usize{ 1, 2, 3 }) |inc| {
            var x: [check.capacity]T = undefined;
            check.fill(T, &x, n + 1);
            const W = check.Wide(F);
            var want: W = 0;
            for (0..n) |i| {
                const v = x[i * inc];
                want += if (comptime check.isComplex(T)) @abs(@as(W, v.re)) + @abs(@as(W, v.im)) else @abs(@as(W, v));
            }
            const got = kern(@intCast(n), @as([*]F, @ptrCast(&x)), @intCast(inc));
            check.expectCloseRef(F, want, got, check.tolerance(F, 2 * n)) catch |e| return check.case(e, "asum", T, n, .{ inc, 0 });
        }
    }
}

fn checkAll(vt: *const Vtable) !void {
    try checkAsum(f32, vt.sasum_k.?);
    try checkAsum(f64, vt.dasum_k.?);
    try checkAsum(std.math.Complex(f32), vt.casum_k.?);
    try checkAsum(std.math.Complex(f64), vt.zasum_k.?);
}

test "portable asum kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkAsum(XDouble, vt.qasum_k.?);
    try checkAsum(std.math.Complex(XDouble), vt.xasum_k.?);
}

test "portable asum kernels at the extended type match a plain loop" {
    try checkAsum(check.Ext, &Asum(check.Ext, 1).k);
    try checkAsum(std.math.Complex(check.Ext), &Asum(check.Ext, 2).k);
}

test "every linked core's asum kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}

test "the extended asum sums fold in asum.S's and zasum.S's order, bit for bit" {
    inline for (.{ f80, f128 }) |F| {
        var x: [check.capacity * 2]F = undefined;
        for (check.sizes) |n| for ([_]usize{ 1, 3 }) |inc| {
            check.fill(F, &x, n + 11);
            const len: BlasLong = @intCast(n);
            const li: BlasLong = @intCast(inc);
            inline for (.{ 1, 2 }) |w| {
                try check.expectBits(F, check.plainExtendedSum(F, w, .abs, n, &x, inc), util.extendedSum(F, w, .abs, len, &x, li), "asum", n, inc);
            }
        };
    }
}

test "the extended asum sums give the values traced by hand from asum.S and zasum.S" {
    // 2^64 + 1 rounds to 2^64 in f80, so each value below depends on which
    // accumulator a term lands in and on the fold order. Traced once from
    // kernel/x86_64/: eight floats a trip by float % 4, tail floats into a0
    // (zasum: re to a0, im to a1), folded a3 + ((a2 + a0) + a1).
    const F = f80;
    const t = std.math.ldexp(@as(F, 1), 64);
    // a0 = 1, a1 = 1, a2 = t: (t + 1) + 1 = t; a tail 1 into a0 gives
    // (t + 2) + 1 = t + 4 (ties to even), into a1 it would give t + 2
    var ax = [_]F{ 1, 1, t, 0, 0, 0, 0, 0, 1 };
    try zigstd.testExpectEqual(t, util.extendedSum(F, 1, .abs, 8, &ax, 1));
    try zigstd.testExpectEqual(t + 4, util.extendedSum(F, 1, .abs, 9, &ax, 1));
    // zasum: four complex a trip; the tail's imaginary part goes to a1
    var zax = [_]F{ 1, 1, t, 0, 0, 0, 0, 0, 0, 1 };
    try zigstd.testExpectEqual(t, util.extendedSum(F, 2, .abs, 4, &zax, 1));
    try zigstd.testExpectEqual(t + 2, util.extendedSum(F, 2, .abs, 5, &zax, 1));
}
