//! BLAS kernel registry for architectures without a registry of their own: no kernels are ported yet, so no per-core
//! object is built and every core keeps the portable kernels. See
//! `x86.zig` for what a registry provides and how to fill one in.

const std = @import("std");
const Vtable = @import("../Vtable.zig");
const coretypes = @import("../cpu/_coretypes.zig");

pub const Core = coretypes.GenericCoreType;

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    _ = arch;
    return false;
}

/// Whether any routine has a kernel for `core`.
pub fn hasKernels(comptime arch: std.Target.Cpu.Arch, comptime core: Core) bool {
    _ = arch;
    _ = core;
    return false;
}

/// Exports every kernel `core` has.
pub fn exportCore(comptime core: Core) void {
    _ = core;
}

/// No kernels; the portable ones stay.
pub fn install(vt: *Vtable, core: Core) void {
    _ = vt;
    _ = core;
}
