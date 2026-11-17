//! pffft's vector macro layer (`src/simd/pf_float.h`, `pf_sse1_float.h`,
//! `pf_neon_float.h`, `pf_scalar_float.h`'s 4xScalar mode and
//! `pf_fma_x86.h`; marton78/pffft efc1dd2d, UCAR and BSD-3, LICENSE) as
//! `@Vector(4, f32)` code. SSE on x86, NEON on arm and aarch64, the
//! 4xScalar mode elsewhere: every target runs pffft's 4-lane algorithm, so
//! the unordered layout is the same everywhere.

const std = @import("std");
const sysinfo = @import("builtin");
const vector = @import("../../linalg/vector.zig");
const matrix = @import("../../linalg/matrix.zig");

pub const V = @Vector(4, f32);
/// pffft's `v4sf *`: callers' slices hold only f32 alignment.
pub const Ptr = [*]align(4) V;
pub const ConstPtr = [*]align(4) const V;

const cpu = sysinfo.cpu;

/// Whether pffft's `VMADD` and `VMSUB` fuse here: SSE with FMA3
/// (`PFFFT_X86_HAVE_FMA`), NEON on aarch64, NEON with VFPv4 on arm
/// (`__ARM_FEATURE_FMA`); never in the 4xScalar mode.
pub const fused = switch (cpu.arch) {
    .x86, .x86_64 => std.Target.x86.featureSetHas(cpu.features, .sse) and std.Target.x86.featureSetHas(cpu.features, .fma),
    .aarch64, .aarch64_be => std.Target.aarch64.featureSetHas(cpu.features, .neon),
    .arm, .armeb, .thumb, .thumbeb => std.Target.arm.featureSetHas(cpu.features, .neon) and std.Target.arm.featureSetHas(cpu.features, .vfp4d16sp),
    else => false,
};

/// pffft's `M_PI`, `M_SQRT2` and `M_SQRT1_2`.
pub const m_pi = 3.14159265358979323846;
pub const m_sqrt2 = 1.41421356237309504880;
pub const m_sqrt1_2 = 0.707106781186547524400844362104849039;

/// A C double literal cast to float, `(float)x`: rounded twice, as pffft's.
pub fn lit(comptime x: comptime_float) f32 {
    return @floatCast(@as(f64, x));
}

/// `LD_PS1`.
pub inline fn ld(f: f32) V {
    return @splat(f);
}

/// `VMADD(a, b, c)`: `a b + c`.
pub inline fn madd(a: V, b: V, c: V) V {
    return if (fused) @mulAdd(V, a, b, c) else a * b + c;
}

/// `VMSUB(a, b, c)`: `c - a b` (`_mm_fnmadd_ps`, `vfmsq_f32`).
pub inline fn msub(a: V, b: V, c: V) V {
    return if (fused) @mulAdd(V, -a, b, c) else c - a * b;
}

/// `SVMUL`.
pub inline fn smul(f: f32, v: V) V {
    return ld(f) * v;
}

/// `SVMADD`.
pub inline fn smadd(f: f32, v: V, c: V) V {
    return madd(ld(f), v, c);
}

/// `SVMSUB`.
pub inline fn smsub(f: f32, v: V, c: V) V {
    return msub(ld(f), v, c);
}

/// `VCPLXMUL`: `(ar, ai) *= (br, bi)`.
pub inline fn cplxMul(ar: *V, ai: *V, br: V, bi: V) void {
    const tmp = ar.* * bi;
    ar.* = msub(ai.*, bi, ar.* * br);
    ai.* = madd(ai.*, br, tmp);
}

/// `VCPLXMULCONJ`: `(ar, ai) *= conj(br, bi)`.
pub inline fn cplxMulConj(ar: *V, ai: *V, br: V, bi: V) void {
    const tmp = ai.* * bi;
    ai.* = msub(ar.*, bi, ai.* * br);
    ar.* = madd(ar.*, br, tmp);
}

/// `INTERLEAVE2`.
pub inline fn interleave2(a: V, b: V) [2]V {
    return vector.interleave(a, b);
}

/// `UNINTERLEAVE2`.
pub inline fn uninterleave2(a: V, b: V) [2]V {
    return vector.uninterleave(a, b);
}

/// `VSWAPHL`.
pub inline fn swapHL(a: V, b: V) V {
    return vector.swapHalves(a, b);
}

/// `VTRANSPOSE4`, in place.
pub inline fn transpose4(x0: *V, x1: *V, x2: *V, x3: *V) void {
    const t = matrix.transposeVecs(4, f32, .{ x0.*, x1.*, x2.*, x3.* });
    x0.* = t[0];
    x1.* = t[1];
    x2.* = t[2];
    x3.* = t[3];
}

/// The floats behind `p`, for pffft's `v4sf_union` lane accesses.
pub inline fn floats(p: Ptr) [*]f32 {
    return @ptrCast(p);
}

pub inline fn constFloats(p: ConstPtr) [*]const f32 {
    return @ptrCast(p);
}
