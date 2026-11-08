//! `AXPY`: `y := alpha*x + y`, and `AXPYC` with `x` conjugated. The portable
//! kernels (`*_generic`) and tests of every core's kernels against a plain
//! loop; the interface is `blas1.interface.axpy`. Complex values are
//! interleaved floats in the kernels, as in the reference.

const std = @import("std");
const sysinfo = @import("builtin");
const zigstd = @import("../../../../zigstd.zig");
const util = @import("../../util.zig");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const fma = @import("../../fma.zig");
const corespec = @import("../../cpu/corespec.zig");
pub const kernel = @import("driver.zig");

const BlasLong = types.BlasLong;
const XDouble = types.XDouble;
const Complex = std.math.Complex;

// ---- kernel ABI and drivers -------------------------------------------------

pub const RealKernel = kernel.RealKernel;
pub const ComplexKernel = kernel.ComplexKernel;
pub const realEntry = kernel.realEntry;
pub const complexEntry = kernel.complexEntry;
pub const Real = kernel.Real;
pub const Cplx = kernel.Cplx;

/// Portable block kernel: `y[0..n1] += alpha * x[0..n1]` on `@Vector`s, the
/// `V_SIMD` branch of `daxpy_kernel_8`/`saxpy_kernel_16`.
fn SimdBlock(comptime T: type) type {
    return struct {
        pub fn run(n1: BlasLong, alpha: T, x: [*]T, y: [*]T) void {
            const lanes = comptime std.simd.suggestVectorLength(T) orelse 1;
            const V = @Vector(lanes, T);
            const valpha: V = @splat(alpha);
            var i: usize = 0;
            const end: usize = @intCast(n1);
            while (i < end) : (i += lanes) {
                const vx: V = x[i..][0..lanes].*;
                const vy: V = y[i..][0..lanes].*;
                y[i..][0..lanes].* = fma.madd(fma.fuse, V, valpha, vx, vy);
            }
        }
    };
}

/// Whether `T` gets the `@Vector` block kernel in the portable path: not
/// where its lanes flush denormals.
fn hasSimdBlock(comptime T: type) bool {
    return (T == f32 or T == f64) and util.vectorIeee(T);
}

/// The blocks of `kernel/x86_64/daxpy.c` and `saxpy.c` (`n & -16`, `n & -32`).
fn GenericReal(comptime T: type) type {
    return Real(T, 128 / @sizeOf(T), if (hasSimdBlock(T)) SimdBlock(T) else null, fma.fuseFor(T));
}

fn GenericCplx(comptime F: type, comptime conj: bool) type {
    return Cplx(F, conj, 16, null, fma.fuseFor(F));
}

// ---- generic implementations (vtable targets) ------------------------------

pub const saxpy_generic = realEntry(f32, GenericReal(f32));
pub const daxpy_generic = realEntry(f64, GenericReal(f64));
pub const qaxpy_generic = realEntry(XDouble, GenericReal(XDouble));
pub const caxpy_generic = complexEntry(f32, GenericCplx(f32, false));
pub const caxpyc_generic = complexEntry(f32, GenericCplx(f32, true));
pub const zaxpy_generic = complexEntry(f64, GenericCplx(f64, false));
pub const zaxpyc_generic = complexEntry(f64, GenericCplx(f64, true));
pub const xaxpy_generic = complexEntry(XDouble, GenericCplx(XDouble, false));
pub const xaxpyc_generic = complexEntry(XDouble, GenericCplx(XDouble, true));

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    return struct {
        pub const axpy_k = if (types.isComplex(T)) complexEntry(F, GenericCplx(F, false)) else realEntry(F, GenericReal(F));
        pub const axpyc_k = complexEntry(F, GenericCplx(F, true));
    };
}

/// Points every axpy slot at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.saxpy_k = saxpy_generic;
    vt.daxpy_k = daxpy_generic;
    vt.qaxpy_k = qaxpy_generic;
    vt.caxpy_k = caxpy_generic;
    vt.caxpyc_k = caxpyc_generic;
    vt.zaxpy_k = zaxpy_generic;
    vt.zaxpyc_k = zaxpyc_generic;
    vt.xaxpy_k = xaxpy_generic;
    vt.xaxpyc_k = xaxpyc_generic;
}

const isComplex = check.isComplex;

// ---- tests: portable kernels vs. a scalar reference ------------------------

const check = @import("../../_test/check.zig");
const fill = check.fill;

const FloatOf = check.FloatOf;

/// `expectClose`, or with a zero stride `check.expectStrided`'s bound.
fn expectStrided(comptime T: type, want: []const T, got: []const T, n: usize, inc: [2]isize) !void {
    if (inc[0] != 0 and inc[1] != 0) return expectClose(T, want, got);
    check.expectStrided(T, want, got, 0, n, inc) catch |e| return check.case(e, "axpy", T, n, inc);
}

fn expectClose(comptime T: type, want: []const T, got: []const T) !void {
    // by width: past 64 bits nothing fuses, so the kernels give the loop's bits
    const tol: FloatOf(T) = switch (@typeInfo(FloatOf(T)).float.bits) {
        32 => 1e-5,
        64 => 1e-12,
        else => 0,
    };
    try zigstd.testExpectEqual(want.len, got.len);
    for (want, got) |w, g| {
        if (comptime isComplex(T)) {
            try zigstd.testExpectApproxEqAbs(w.re, g.re, tol);
            try zigstd.testExpectApproxEqAbs(w.im, g.im, tol);
        } else {
            try zigstd.testExpectApproxEqAbs(w, g, tol);
        }
    }
}

/// Checks `kern` (vtable ABI) against a plain loop for contiguous and
/// strided data at sizes that straddle every block and unroll boundary.
pub fn checkReal(comptime T: type, kern: RealKernel(T)) !void {
    const sizes = [_]usize{ 0, 1, 3, 4, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 257 };
    const alpha: T = 1.5;
    for (sizes) |n| {
        var xb: [800]T = undefined;
        var yb: [800]T = undefined;
        var want: [800]T = undefined;
        fill(T, &xb, @intCast(n + 1));
        fill(T, &yb, @intCast(n + 7));
        for (check.zero_strides) |inc| {
            want = yb;
            var y = yb;
            for (0..n) |i| want[check.pos(n, inc[1], i)] += alpha * xb[check.pos(n, inc[0], i)];
            const py: [*]T = @ptrCast(&y[check.pos(n, inc[1], 0)]);
            _ = kern(@intCast(n), 0, 0, alpha, @ptrCast(&xb[check.pos(n, inc[0], 0)]), inc[0], py, inc[1], py, 0);
            try expectStrided(T, &want, &y, n, inc);
        }
    }
}

pub fn checkComplex(comptime F: type, comptime conj: bool, kern: ComplexKernel(F)) !void {
    const C = Complex(F);
    const sizes = [_]usize{ 0, 1, 3, 4, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 129 };
    const alpha = C.init(1.5, -0.75);
    for (sizes) |n| {
        var xb: [400]C = undefined;
        var yb: [400]C = undefined;
        var want: [400]C = undefined;
        fill(C, &xb, @intCast(n + 1));
        fill(C, &yb, @intCast(n + 7));
        for (check.zero_strides) |inc| {
            want = yb;
            var y = yb;
            for (0..n) |i| {
                const xv = xb[check.pos(n, inc[0], i)];
                const yp = &want[check.pos(n, inc[1], i)];
                if (conj) {
                    yp.re += alpha.re * xv.re + alpha.im * xv.im;
                    yp.im -= alpha.re * xv.im - alpha.im * xv.re;
                } else {
                    yp.re += alpha.re * xv.re - alpha.im * xv.im;
                    yp.im += alpha.re * xv.im + alpha.im * xv.re;
                }
            }
            const py: [*]F = @ptrCast(&y[check.pos(n, inc[1], 0)]);
            _ = kern(@intCast(n), 0, 0, alpha.re, alpha.im, @ptrCast(&xb[check.pos(n, inc[0], 0)]), inc[0], py, inc[1], py, 0);
            try expectStrided(C, &want, &y, n, inc);
        }
    }
}

const cpu = @import("../../cpu/__root__.zig");

test "portable axpy kernels match a scalar loop (all types, contiguous and strided)" {
    try checkReal(f32, saxpy_generic);
    try checkReal(f64, daxpy_generic);
    try checkReal(XDouble, qaxpy_generic);
    try checkComplex(f32, false, caxpy_generic);
    try checkComplex(f32, true, caxpyc_generic);
    try checkComplex(f64, false, zaxpy_generic);
    try checkComplex(f64, true, zaxpyc_generic);
    try checkComplex(XDouble, false, xaxpy_generic);
    try checkComplex(XDouble, true, xaxpyc_generic);
}

test "portable axpy kernels at the extended type give a plain unfused loop's bits" {
    const E = check.Ext;
    try checkReal(E, realEntry(E, GenericReal(E)));
    try checkComplex(E, false, complexEntry(E, GenericCplx(E, false)));
    try checkComplex(E, true, complexEntry(E, GenericCplx(E, true)));
}

test "a core's axpy kernels fill every axpy slot, for the host's core type" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    cpu.kernel.install(&vt, cpu.detect());
    try zigstd.testExpectTrue(vt.saxpy_k != null and vt.daxpy_k != null and vt.qaxpy_k != null);
    try zigstd.testExpectTrue(vt.caxpy_k != null and vt.caxpyc_k != null);
    try zigstd.testExpectTrue(vt.zaxpy_k != null and vt.zaxpyc_k != null);
    try zigstd.testExpectTrue(vt.xaxpy_k != null and vt.xaxpyc_k != null);
    try checkReal(f32, vt.saxpy_k.?);
    try checkReal(f64, vt.daxpy_k.?);
    try checkComplex(f32, false, vt.caxpy_k.?);
    try checkComplex(f32, true, vt.caxpyc_k.?);
    try checkComplex(f64, false, vt.zaxpy_k.?);
    try checkComplex(f64, true, vt.zaxpyc_k.?);
}

fn checkAll(vt: *const Vtable) !void {
    var failed: ?anyerror = null;
    check.keep(&failed, checkReal(f32, vt.saxpy_k.?));
    check.keep(&failed, checkReal(f64, vt.daxpy_k.?));
    check.keep(&failed, checkComplex(f32, false, vt.caxpy_k.?));
    check.keep(&failed, checkComplex(f32, true, vt.caxpyc_k.?));
    check.keep(&failed, checkComplex(f64, false, vt.zaxpy_k.?));
    check.keep(&failed, checkComplex(f64, true, vt.zaxpyc_k.?));
    if (failed) |e| return e;
}

test "every linked core's axpy kernels match a scalar loop, where this machine runs them" {
    try check.eachCore(checkAll);
}

/// `(1+e)(1+e) - (1+2e)` is `e*e` fused, 0 rounded: scalar (`n` 1) and
/// vector (`n` 64) paths.
fn checkFusedReal(comptime T: type, kern: RealKernel(T)) !void {
    const e = std.math.floatEps(T);
    for ([_]usize{ 1, 64 }) |n| {
        var x: [64]T = @splat(1 + e);
        var y: [64]T = @splat(-(1 + 2 * e));
        _ = kern(@intCast(n), 0, 0, 1 + e, &x, 1, &y, 1, &y, 0);
        for (y[0..n]) |v| try zigstd.testExpectEqual(e * e, v);
    }
}

fn checkFused(vt: *const Vtable) !void {
    if (comptime sysinfo.target.cpu.arch.isMips()) {
        const c = check.core orelse return;
        if (!(corespec.mips(c).fuse_fma orelse false)) return;
        try checkFusedReal(f32, vt.saxpy_k.?);
        try checkFusedReal(f64, vt.daxpy_k.?);
    }
}

test "the Loongson 3A cores fuse axpy, as their flags do" {
    if (comptime !sysinfo.target.cpu.arch.isMips()) return error.SkipZigTest;
    // every mips rig is QEMU, which rounds `madd.fmt` (issues.md)
    if (check.rig()) return error.SkipZigTest;
    try check.eachCore(checkFused);
}

test {
    zigstd.testRefAllDecls(@This());
}
