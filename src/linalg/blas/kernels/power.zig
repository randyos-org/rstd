//! BLAS kernel registry for POWER: what `build/blaskern.zig` and `_kernel.zig`
//! ask of a registry. See `x86.zig` for what one provides.

const std = @import("std");
const sysinfo = @import("builtin");
const coretypes = @import("../cpu/_coretypes.zig");
const axpy = @import("../blas1/axpy/power/kernel.zig");

pub const Core = coretypes.PowerCoreType;

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch.isPowerPC();
}

/// Whether any routine has a kernel for `core` on `arch`.
pub fn hasKernels(comptime arch: std.Target.Cpu.Arch, comptime core: Core) bool {
    return axpy.has(arch, core);
}

/// Exports every kernel `core` has, for the architecture being compiled.
pub fn exportCore(comptime core: Core) void {
    axpy.exportCore(sysinfo.cpu.arch, core);
}
