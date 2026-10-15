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
//!     an FMA unit, so `@mulAdd` never degrades into a libm `fma()` call;
//!   * never for floats wider than `f64` (`fuseFor`): `@mulAdd(f80)` is a
//!     soft `fmal` call, and the reference's x87 `q`/`x` kernels never fuse.
//!
//! Contraction is taken the way a C compiler does it within one expression
//! (clang's `-ffp-contract=on`, and gcc's `fast` for the single-expression
//! statements the reference kernels use): `a*b + c*d` is `fma(a, b, c*d)`,
//! `a*b - c*d` is `fma(a, b, -(c*d))`, a longer chain fuses each product into
//! the running sum from the left, and a product kept in its own statement, or
//! converted (`float` to `double`) before the add, stays rounded.
//!
//! `rsimd.mulAdd` delegates here.

const std = @import("std");
const sysinfo = @import("builtin");
const root = @import("root");
const types = @import("types.zig");

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
    // LoongArch has fmadd.s/fmadd.d in the base F/D ISA.
    .loongarch32, .loongarch64 => std.Target.loongarch.featureSetHasAny(sysinfo.cpu.features, .{ .f, .d }),
    // z/Architecture has the BFP multiply-and-add (MAEBR/MADBR) in its base.
    .s390x => true,
    // Before release 6, MADD.fmt rounds the product, so it equals a*b + c and
    // gcc's use of it never changes a result; R6 adds the fused MADDF.fmt.
    .mips, .mipsel, .mips64, .mips64el => std.Target.mips.featureSetHasAny(sysinfo.cpu.features, .{
        .mips32r6, .mips64r6,
    }),
    // No scalar fused multiply-add: LLVM has none for SPARC (gcc's needs
    // -mfmaf, which OpenBLAS never passes), and WebAssembly's only exists as
    // `relaxed_simd` vector madd (see `vector_madd`).
    else => false,
};

/// Whether a vector multiply-add (`v_muladd_*` of `intrin.h`, MSA `fmadd`)
/// fuses on this target: where scalar fusion does, for MIPS with MSA (whose
/// `fmadd` is fused whatever the scalar ISA), and for WebAssembly with
/// `relaxed_simd` (`wasm_f32x4_relaxed_madd`; the spec lets an engine fuse or
/// not, taken here as fused, as on engines with hardware FMA).
pub const vector_madd = hardware or switch (sysinfo.cpu.arch) {
    .mips, .mipsel, .mips64, .mips64el => std.Target.mips.featureSetHas(sysinfo.cpu.features, .msa),
    .wasm32, .wasm64 => std.Target.wasm.featureSetHas(sysinfo.cpu.features, .relaxed_simd),
    else => false,
};

/// The build's policy: the root's `blas_fuse_fma` if it declares one.
pub const fuse: bool = if (@hasDecl(root, "blas_fuse_fma")) root.blas_fuse_fma else hardware;

/// The same for vector multiply-adds; a root's `blas_fuse_fma_vector` wins,
/// then its `blas_fuse_fma`.
pub const fuse_vector: bool = if (@hasDecl(root, "blas_fuse_fma_vector"))
    root.blas_fuse_fma_vector
else if (@hasDecl(root, "blas_fuse_fma")) root.blas_fuse_fma else vector_madd;

/// Whether `T`, a float or a vector of them, can fuse at all: not past `f64`.
pub fn fusible(comptime T: type) bool {
    const E = switch (@typeInfo(T)) {
        .vector => |v| v.child,
        else => T,
    };
    return !types.isExtended(E);
}

/// `fuse` for multiply-adds in `T`.
pub fn fuseFor(comptime T: type) bool {
    return fuse and fusible(T);
}

/// `a*b + c`, fused when `fuseFor` says so (`@mulAdd` is strictly
/// single-rounded, so it must not be used on targets without an FMA unit).
pub inline fn mulAdd(comptime T: type, a: T, b: T, c: T) T {
    return if (comptime fuseFor(T)) @mulAdd(T, a, b, c) else a * b + c;
}

test "mulAdd agrees with a*b+c on exact values" {
    try std.testing.expectEqual(@as(f64, 14.0), mulAdd(f64, 2.0, 3.0, 8.0));
    try std.testing.expectEqual(@as(f32, 14.0), mulAdd(f32, 2.0, 3.0, 8.0));
}

test "extended floats never fuse" {
    inline for (.{ f80, f128, @Vector(2, f80) }) |T| try std.testing.expect(!fuseFor(T));
    // by width, not identity: a 64-bit `c_longdouble` of its own fuses as f64 does
    try std.testing.expectEqual(fuse and @typeInfo(c_longdouble).float.bits <= 64, fuseFor(c_longdouble));
    // (1 + 2^-63)(1 - 2^-63) - 1 is -2^-126 fused, 0 rounded first
    const e: f80 = std.math.ldexp(@as(f80, 1), -63);
    try std.testing.expectEqual(@as(f80, 0), mulAdd(f80, 1 + e, 1 - e, -1));
}
