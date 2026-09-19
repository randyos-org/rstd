//! x86_64 BLAS kernel registry (`kernels/<arch>.zig`): which routines have kernels for which cores,
//! and how to export them. The per-core kernel object root (`_kernel.zig`)
//! exports through it, and the build (`build/blaskern.zig`) asks it which
//! cores need an object at all.
//!
//! To add a routine, import its `kernel_x86_64.zig` here and add it to
//! `hasKernels` and `exportCore`.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const axpy = @import("../blas1/axpy/kernel_x86_64.zig");

pub const Core = coretypes.X86CoreType;

/// Whether this registry has kernels for `arch` (the 32-bit x86 has none).
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch == .x86_64;
}

/// Whether any routine has a kernel for `core`.
pub fn hasKernels(comptime core: Core) bool {
    return axpy.has(core, "daxpy") or axpy.has(core, "saxpy") or
        axpy.has(core, "caxpy") or axpy.has(core, "zaxpy");
}

/// Exports every kernel `core` has.
pub fn exportCore(comptime core: Core) void {
    axpy.exportCore(core);
}
