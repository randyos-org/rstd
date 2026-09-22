//! BLAS kernel registry for MIPS: what `build/blaskern.zig` and `_kernel.zig`
//! ask of a registry. See `x86.zig` for what one provides.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const axpy = @import("../blas1/axpy/mips64/kernel.zig");

pub const Core = coretypes.MipsCoreType;

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch.isMIPS();
}

/// Whether `core`'s kernels run on `cpu`: a core's code is for one word size,
/// the Loongson 3A kernels are release 2 code, which release 6 does not run,
/// and the SIMD kernels need the release that has the SIMD architecture.
pub fn runsOn(cpu: std.Target.Cpu, core: Core) bool {
    const wide = cpu.arch.isMIPS64();
    return switch (core) {
        .mips24k, .mips1004k, .p5600 => !wide and (core != .p5600 or cpu.has(.mips, .mips32r5)),
        .i6400, .i6500, .p6600 => wide and cpu.has(.mips, .mips64r6),
        .loongson3r3 => wide and !cpu.has(.mips, .mips64r6),
        // release 2 code with the SIMD architecture, whose objects do not link with release 6's
        .loongson3r4 => wide and !cpu.has(.mips, .mips64r6),
        .mips64_generic, .unknown => wide,
    };
}

/// Whether any routine has a kernel for `core`.
pub fn hasKernels(comptime arch: std.Target.Cpu.Arch, comptime core: Core) bool {
    _ = arch;
    return axpy.has(core);
}

/// Exports every kernel `core` has.
pub fn exportCore(comptime core: Core) void {
    axpy.exportCore(core);
}
