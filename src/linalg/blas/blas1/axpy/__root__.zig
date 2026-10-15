//! `AXPY`: `y := alpha*x + y`, and `AXPYC` with `x` conjugated. The portable
//! kernels (`*_generic`) and tests of every core's kernels against a plain
//! loop; the interface is `blas1.interface.axpy`. Complex values are
//! interleaved floats in the kernels, as in the reference.

const std = @import("std");
const zigstd = @import("../../../../zigstd.zig");
const util = @import("../../util.zig");
const types = @import("../../types.zig");
const cpu = @import("../../cpu/__root__.zig");
const blas = @import("../../__root__.zig");
const Vtable = @import("../../Vtable.zig");
const rsimd = @import("../../../../simd/__root__.zig");
const fma = @import("../../fma.zig");
const registry = @import("../../asm_registry.zig");

const BlasLong = types.BlasLong;
const XDouble = types.XDouble;
const Complex = std.math.Complex;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

// ---- kernel ABI and drivers -------------------------------------------------

pub const kernel = @import("kernel.zig");
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
            const lanes = comptime rsimd.vectorLen(T);
            const V = @Vector(lanes, T);
            const valpha: V = @splat(alpha);
            var i: usize = 0;
            const end: usize = @intCast(n1);
            while (i < end) : (i += lanes) {
                const vx: V = x[i..][0..lanes].*;
                const vy: V = y[i..][0..lanes].*;
                y[i..][0..lanes].* = rsimd.mulAdd(V, valpha, vx, vy);
            }
        }
    };
}

/// Whether `T` gets the `@Vector` block kernel in the portable path.
fn hasSimdBlock(comptime T: type) bool {
    return T == f32 or T == f64;
}

fn GenericReal(comptime T: type) type {
    return Real(T, 128 / @sizeOf(T), if (hasSimdBlock(T)) SimdBlock(T) else null, fma.fuse);
}

fn GenericCplx(comptime F: type, comptime conj: bool) type {
    return Cplx(F, conj, 16, null, fma.fuse);
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

fn isComplex(comptime T: type) bool {
    return T == Complex32 or T == Complex64 or T == ComplexXD;
}

// ---- tests: portable kernels vs. a scalar reference ------------------------

fn fill(comptime T: type, buf: []T, seed: u32) void {
    var prng = std.Random.DefaultPrng.init(seed);
    const r = prng.random();
    for (buf) |*v| {
        if (comptime isComplex(T)) {
            v.* = T.init(@floatCast(r.float(f64) * 4 - 2), @floatCast(r.float(f64) * 4 - 2));
        } else {
            v.* = @floatCast(r.float(f64) * 4 - 2);
        }
    }
}

fn FloatOf(comptime T: type) type {
    return if (isComplex(T)) @FieldType(T, "re") else T;
}

fn expectClose(comptime T: type, want: []const T, got: []const T) !void {
    const tol: FloatOf(T) = if (FloatOf(T) == f32) 1e-5 else 1e-12;
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
        inline for (.{ [2]usize{ 1, 1 }, [2]usize{ 2, 1 }, [2]usize{ 1, 3 }, [2]usize{ 2, 2 } }) |inc| {
            want = yb;
            var y = yb;
            for (0..n) |i| want[i * inc[1]] += alpha * xb[i * inc[0]];
            _ = kern(@intCast(n), 0, 0, alpha, &xb, inc[0], &y, inc[1], &y, 0);
            try expectClose(T, &want, &y);
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
        inline for (.{ [2]usize{ 1, 1 }, [2]usize{ 2, 1 }, [2]usize{ 1, 3 } }) |inc| {
            want = yb;
            var y = yb;
            for (0..n) |i| {
                const xv = xb[i * inc[0]];
                const yp = &want[i * inc[1]];
                if (conj) {
                    yp.re += alpha.re * xv.re + alpha.im * xv.im;
                    yp.im -= alpha.re * xv.im - alpha.im * xv.re;
                } else {
                    yp.re += alpha.re * xv.re - alpha.im * xv.im;
                    yp.im += alpha.re * xv.im + alpha.im * xv.re;
                }
            }
            _ = kern(@intCast(n), 0, 0, alpha.re, alpha.im, @ptrCast(&xb), inc[0], @ptrCast(&y), inc[1], @ptrCast(&y), 0);
            try expectClose(C, &want, &y);
        }
    }
}

test "generic kernels match a scalar loop (all types, contiguous and strided)" {
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

test "a core's kernels fill every axpy slot, for the host's core type" {
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

test "every linked core's kernels match a scalar loop, where this machine runs them" {
    inline for (@typeInfo(cpu.CoreType).@"enum".field_names) |name| {
        const core = @field(cpu.CoreType, name);
        if (registry.isLinked(core) and (!@hasDecl(cpu.arch, "canRun") or cpu.arch.canRun(core))) {
            var vt: Vtable = .{};
            installGeneric(&vt);
            cpu.kernel.install(&vt, core);
            try checkReal(f32, vt.saxpy_k.?);
            try checkReal(f64, vt.daxpy_k.?);
            try checkComplex(f32, false, vt.caxpy_k.?);
            try checkComplex(f32, true, vt.caxpyc_k.?);
            try checkComplex(f64, false, vt.zaxpy_k.?);
            try checkComplex(f64, true, vt.zaxpyc_k.?);
        }
    }
}

test {
    zigstd.testRefAllDecls(@This());
}
