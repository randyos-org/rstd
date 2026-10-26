//! `SWAP`: `x <-> y`. The vtable wrappers, the portable kernels (`*_generic`),
//! and tests of every core's kernels against a plain loop.
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

pub fn sbswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas.kernel(.sbswap_k)(n, dummy0, dummy1, dummy3, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn sswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas.kernel(.sswap_k)(n, dummy0, dummy1, dummy3, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn dswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, dummy: [*]f64, dummy2: BlasLong) i32 {
    return blas.kernel(.dswap_k)(n, dummy0, dummy1, dummy3, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn qswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, dummy: [*]XDouble, dummy2: BlasLong) i32 {
    return blas.kernel(.qswap_k)(n, dummy0, dummy1, dummy3, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn cswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f32, dummy4: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas.kernel(.cswap_k)(n, dummy0, dummy1, dummy3, dummy4, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn zswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f64, dummy4: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, dummy: [*]f64, dummy2: BlasLong) i32 {
    return blas.kernel(.zswap_k)(n, dummy0, dummy1, dummy3, dummy4, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn xswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: XDouble, dummy4: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, dummy: [*]XDouble, dummy2: BlasLong) i32 {
    return blas.kernel(.xswap_k)(n, dummy0, dummy1, dummy3, dummy4, x, inc_x, y, inc_y, dummy, dummy2);
}

// ---- portable kernels ------------------------------------------------------

const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

/// `kernel/arm/swap.c` and `zswap.c`: exchanges `x` and `y`; an element is
/// `w` floats.
fn swapElements(comptime F: type, comptime w: comptime_int, n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) void {
    var i: BlasLong = 0;
    while (i < n) : (i += 1) {
        inline for (0..w) |c| {
            const px = util.at(F, x, i * inc_x * w + c);
            const py = util.at(F, y, i * inc_y * w + c);
            const t = px.*;
            px.* = py.*;
            py.* = t;
        }
    }
}

fn Swap(comptime F: type) type {
    return struct {
        fn real(n: BlasLong, _: BlasLong, _: BlasLong, _: F, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F, _: BlasLong) callconv(types.kernelCc(F)) i32 {
            swapElements(F, 1, n, x, inc_x, y, inc_y);
            return 0;
        }
        fn complex(n: BlasLong, _: BlasLong, _: BlasLong, _: F, _: F, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F, _: BlasLong) callconv(types.kernelCc(F)) i32 {
            swapElements(F, 2, n, x, inc_x, y, inc_y);
            return 0;
        }
    };
}

pub const sswap_generic = &Swap(f32).real;
pub const dswap_generic = &Swap(f64).real;
pub const qswap_generic = &Swap(XDouble).real;
pub const cswap_generic = &Swap(f32).complex;
pub const zswap_generic = &Swap(f64).complex;
pub const xswap_generic = &Swap(XDouble).complex;

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    return struct {
        pub const swap_k = if (types.isComplex(T)) &Swap(F).complex else &Swap(F).real;
    };
}

/// Points every swap slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.sswap_k = sswap_generic;
    vt.dswap_k = dswap_generic;
    vt.qswap_k = qswap_generic;
    vt.cswap_k = cswap_generic;
    vt.zswap_k = zswap_generic;
    vt.xswap_k = xswap_generic;
}

// ---- tests -----------------------------------------------------------------

const std = @import("std");
const check = @import("../../_test/blas1/check.zig");

/// `check.zero_strides` and both zero, which a kernel that loads a block
/// before storing it gets wrong.
const strides = check.zero_strides ++ [_][2]isize{.{ 0, 0 }};

/// `kern` against a plain loop; `T` is the element type, real or complex.
fn checkSwap(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    for (check.sizes) |n| {
        for (strides) |inc| {
            var x: [check.capacity]T = undefined;
            var y: [check.capacity]T = undefined;
            check.fill(T, &x, n + 1);
            check.fill(T, &y, n + 7);
            var want_x = x;
            var want_y = y;
            for (0..n) |i| std.mem.swap(T, &want_x[check.pos(n, inc[0], i)], &want_y[check.pos(n, inc[1], i)]);
            const px: [*]F = @ptrCast(&x[check.pos(n, inc[0], 0)]);
            const py: [*]F = @ptrCast(&y[check.pos(n, inc[1], 0)]);
            if (comptime complex) {
                _ = kern(@intCast(n), 0, 0, 0, 0, px, inc[0], py, inc[1], py, 0);
            } else {
                _ = kern(@intCast(n), 0, 0, 0, px, inc[0], py, inc[1], py, 0);
            }
            check.expectCloseSlices(T, &want_x, &x, 0) catch |e| return check.case(e, "swap x", T, n, inc);
            check.expectCloseSlices(T, &want_y, &y, 0) catch |e| return check.case(e, "swap y", T, n, inc);
        }
    }
}

fn checkAll(vt: *const Vtable) !void {
    var failed: ?anyerror = null;
    check.keep(&failed, checkSwap(f32, vt.sswap_k.?));
    check.keep(&failed, checkSwap(f64, vt.dswap_k.?));
    check.keep(&failed, checkSwap(std.math.Complex(f32), vt.cswap_k.?));
    check.keep(&failed, checkSwap(std.math.Complex(f64), vt.zswap_k.?));
    if (failed) |e| return e;
}

test "portable swap kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkSwap(XDouble, vt.qswap_k.?);
    try checkSwap(std.math.Complex(XDouble), vt.xswap_k.?);
}

test "portable swap kernels at the extended type match a plain loop" {
    try checkSwap(check.Ext, &Swap(check.Ext).real);
    try checkSwap(std.math.Complex(check.Ext), &Swap(check.Ext).complex);
}

test "every linked core's swap kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}
