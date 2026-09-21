//! BLAS kernel registry for MIPS: what `build/blaskern.zig` and `_kernel.zig`
//! ask of a registry. See `x86.zig` for what one provides.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const axpy = @import("../blas1/axpy/mips64/kernel.zig");

pub const Core = coretypes.MipsCoreType;

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch == .mips64 or arch == .mips64el;
}

/// Whether `core`'s kernels run on `cpu`: they are release 2 code, which
/// release 6 does not run.
pub fn runsOn(cpu: std.Target.Cpu, core: Core) bool {
    _ = core;
    return !cpu.has(.mips, .mips64r6);
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
