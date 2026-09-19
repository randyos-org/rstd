//! Whether a multiply-add is fused, decided the way the reference's compile
//! flags decide it.
//!
//! OpenBLAS builds each core's kernels with that core's flags, and only some
//! of them enable FMA (`-mfma`, or a `-march=` that includes it); the fusion
//! itself is then the C compiler's contraction of `a*b + c`. Here a build
//! states its policy once and every multiply-add asks this file:
//!
//!   * the root source file may declare `pub const blas_fuse_fma: bool`; the
//!     per-core kernel objects do, since their CPU features can include FMA
//!     only so an assembler accepts a hand-written microkernel while the
//!     core's flags leave contraction off;
//!   * otherwise the answer follows the target: fused where the hardware has
//!     an FMA unit, so `@mulAdd` never degrades into a libm `fma()` call.
//!
//! `rsimd.mulAdd` delegates here.

const std = @import("std");
const sysinfo = @import("builtin");
const root = @import("root");

/// Whether the target has a hardware fused multiply-add.
pub const hardware = switch (sysinfo.cpu.arch) {
    // x86 FMA3 is an optional extension, not baseline; must check.
    .x86_64, .x86 => std.Target.x86.featureSetHas(sysinfo.cpu.features, .fma),
    // AArch64 scalar and NEON FMA is mandatory since ARMv8-A.
    .aarch64, .aarch64_be => true,
    // 32-bit ARM FMA needs VFPv4 or later, not guaranteed on the ARMv7 baseline.
    .arm, .armeb, .thumb, .thumbeb => std.Target.arm.featureSetHasAny(sysinfo.cpu.features, .{
        .vfp4, .vfp4d16, .vfp4d16sp, .vfp4sp,
    }),
    // The base PowerPC FPU has had fmadd/fmsub since POWER1.
    .powerpc, .powerpcle, .powerpc64, .powerpc64le => true,
    // FMA is part of the RISC-V F/D extensions themselves.
    .riscv32, .riscv32be, .riscv64, .riscv64be => std.Target.riscv.featureSetHasAny(sysinfo.cpu.features, .{
        .f, .d,
    }),
    else => false,
};

/// The build's policy: the root's `blas_fuse_fma` if it declares one.
pub const fuse: bool = if (@hasDecl(root, "blas_fuse_fma")) root.blas_fuse_fma else hardware;

/// `a*b + c`, fused when `fuse` says so (`@mulAdd` is strictly single-rounded,
/// so it must not be used on targets without an FMA unit).
pub inline fn mulAdd(comptime T: type, a: T, b: T, c: T) T {
    return if (fuse) @mulAdd(T, a, b, c) else a * b + c;
}

test "mulAdd agrees with a*b+c on exact values" {
    try std.testing.expectEqual(@as(f64, 14.0), mulAdd(f64, 2.0, 3.0, 8.0));
    try std.testing.expectEqual(@as(f32, 14.0), mulAdd(f32, 2.0, 3.0, 8.0));
}
