//! `ROT`: the plane rotation `x := c*x + s*y`, `y := c*y - s*x`. The vtable
//! wrappers, the portable kernels (`*_generic`), and tests of every core's
//! kernels against a plain loop.
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

pub fn sbrot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, c: f32, s: f32) i32 {
    return blas.kernel(.sbrot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn srot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, c: f32, s: f32) i32 {
    return blas.kernel(.srot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn drot_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, c: f64, s: f64) i32 {
    return blas.kernel(.drot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn qrot_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, c: XDouble, s: XDouble) i32 {
    return blas.kernel(.qrot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn csrot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, c: f32, s: f32) i32 {
    return blas.kernel(.csrot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn zdrot_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, c: f64, s: f64) i32 {
    return blas.kernel(.zdrot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn xqrot_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, c: XDouble, s: XDouble) i32 {
    return blas.kernel(.xqrot_k)(n, x, inc_x, y, inc_y, c, s);
}

// ---- portable kernels ------------------------------------------------------

const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

/// `kernel/arm/rot.c` and `zrot.c`: the plane rotation `x := c*x + s*y`,
/// `y := c*y - s*x`; an element is `w` floats, each rotated alike (`csrot`
/// and `zdrot` take a real `c` and `s`).
fn Rot(comptime F: type, comptime w: comptime_int) type {
    return struct {
        fn k(n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, c: F, s: F) callconv(types.kernelCc(F)) i32 {
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                inline for (0..w) |j| {
                    const px = util.at(F, x, i * inc_x * w + j);
                    const py = util.at(F, y, i * inc_y * w + j);
                    const t = c * px.* + s * py.*;
                    py.* = c * py.* - s * px.*;
                    px.* = t;
                }
            }
            return 0;
        }
    };
}

pub const srot_generic = &Rot(f32, 1).k;
pub const drot_generic = &Rot(f64, 1).k;
pub const qrot_generic = &Rot(XDouble, 1).k;
pub const csrot_generic = &Rot(f32, 2).k;
pub const zdrot_generic = &Rot(f64, 2).k;
pub const xqrot_generic = &Rot(XDouble, 2).k;

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    const w = if (types.isComplex(T)) 2 else 1;
    return struct {
        pub const rot_k = &Rot(F, w).k;
    };
}

/// Points every rot slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.srot_k = srot_generic;
    vt.drot_k = drot_generic;
    vt.qrot_k = qrot_generic;
    vt.csrot_k = csrot_generic;
    vt.zdrot_k = zdrot_generic;
    vt.xqrot_k = xqrot_generic;
}

// ---- tests -----------------------------------------------------------------

const std = @import("std");
const sysinfo = @import("builtin");
const check = @import("../../_test/blas1/check.zig");
const fma = @import("../../fma.zig");
const driver = @import("kernel.zig");

/// `check.zero_strides` and both zero, which a kernel that loads a block
/// before storing it gets wrong.
const strides = check.zero_strides ++ [_][2]isize{.{ 0, 0 }};

/// `kern` against a plain loop; `T` is the element type, real or complex.
fn checkRot(comptime T: type, kern: anytype) !void {
    const F = check.FloatOf(T);
    const c: F = 0.8;
    const s: F = -0.6;
    for (check.sizes) |n| {
        for (strides) |inc| {
            var x: [check.capacity]T = undefined;
            var y: [check.capacity]T = undefined;
            check.fill(T, &x, n + 1);
            check.fill(T, &y, n + 7);
            var want_x = x;
            var want_y = y;
            for (0..n) |i| {
                const ix = check.pos(n, inc[0], i);
                const iy = check.pos(n, inc[1], i);
                const a = want_x[ix];
                const b = want_y[iy];
                if (comptime check.isComplex(T)) {
                    want_x[ix] = .init(c * a.re + s * b.re, c * a.im + s * b.im);
                    want_y[iy] = .init(c * b.re - s * a.re, c * b.im - s * a.im);
                } else {
                    want_x[ix] = c * a + s * b;
                    want_y[iy] = c * b - s * a;
                }
            }
            const px: [*]F = @ptrCast(&x[check.pos(n, inc[0], 0)]);
            _ = kern(@intCast(n), px, inc[0], @ptrCast(&y[check.pos(n, inc[1], 0)]), inc[1], c, s);
            const tol = check.tolerance(F, 2);
            check.expectStrided(T, &want_x, &x, tol, n, inc) catch |e| return check.case(e, "rot x", T, n, inc);
            check.expectStrided(T, &want_y, &y, tol, n, inc) catch |e| return check.case(e, "rot y", T, n, inc);
        }
    }
}

fn checkAll(vt: *const Vtable) !void {
    var failed: ?anyerror = null;
    check.keep(&failed, checkRot(f32, vt.srot_k.?));
    check.keep(&failed, checkRot(f64, vt.drot_k.?));
    check.keep(&failed, checkRot(std.math.Complex(f32), vt.csrot_k.?));
    check.keep(&failed, checkRot(std.math.Complex(f64), vt.zdrot_k.?));
    if (failed) |e| return e;
}

test "portable rot kernels match a plain loop" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkRot(XDouble, vt.qrot_k.?);
    try checkRot(std.math.Complex(XDouble), vt.xqrot_k.?);
}

test "portable rot kernels at the extended type match a plain loop" {
    try checkRot(check.Ext, &Rot(check.Ext, 1).k);
    try checkRot(std.math.Complex(check.Ext), &Rot(check.Ext, 2).k);
}

test "every linked core's rot kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}

test "x86_64's rot drivers rotate every element as the reference's blocks do, bit for bit" {
    // the cores these drivers serve fuse; so must the test's own C loops
    if (sysinfo.target.cpu.arch != .x86_64 or !fma.fuse) return error.SkipZigTest;
    const index = @import("x86_64/uk/__root__.zig");
    const Vector = @import("x86_64/kernel.zig").VectorBlock;
    inline for (.{ .srot_k, .drot_k }) |slot| {
        const F = driver.Elem(slot);
        // a vector block (haswell) and the C one (cooperlake)
        inline for (.{ "haswell", "cooperlake" }) |core| {
            const kern = driver.entryOf(slot, core, comptime index.find(slot, core).?, Vector);
            const c: F = 0.8;
            const s: F = -0.6;
            for (check.sizes) |n| for ([_]isize{ 1, 2 }) |inc| {
                var x: [check.capacity]F = undefined;
                var y: [check.capacity]F = undefined;
                check.fill(F, &x, n + 3);
                check.fill(F, &y, n + 5);
                var want_x = x;
                var want_y = y;
                for (0..n) |i| {
                    const k = i * @as(usize, @intCast(inc));
                    driver.rotate(F, true, &want_x[k], &want_y[k], c, s);
                }
                _ = kern(@intCast(n), &x, inc, &y, inc, c, s);
                try std.testing.expectEqualSlices(F, &want_x, &x);
                try std.testing.expectEqualSlices(F, &want_y, &y);
            };
        }
    }
}
