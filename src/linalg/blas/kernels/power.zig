//! BLAS kernel registry for power: no kernels are ported yet, so no per-core
//! object is built and every core keeps the portable kernels. See
//! `x86.zig` for what a registry provides and how to fill one in.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");

pub const Core = coretypes.PowerCoreType;

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    _ = arch;
    return false;
}

/// Whether any routine has a kernel for `core`.
pub fn hasKernels(comptime core: Core) bool {
    _ = core;
    return false;
}

/// Exports every kernel `core` has.
pub fn exportCore(comptime core: Core) void {
    _ = core;
}
