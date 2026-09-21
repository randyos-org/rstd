//! x86 BLAS kernel registry (`kernels/<arch>.zig`), for both x86_64 and the
//! 32-bit x86: which routines have kernels for which cores, and how to export
//! them. The per-core kernel object root (`_kernel.zig`)
//! exports through it, and the build (`build/blaskern.zig`) asks it which
//! cores need an object at all.
//!
//! To add a routine, import its `x86_64/kernel.zig` and `x86/kernel.zig` here
//! and add it to `hasKernels` and `exportCore`.

const std = @import("std");
const sysinfo = @import("builtin");
const coretypes = @import("../cpu/_coretypes.zig");
const axpy = @import("../blas1/axpy/x86_64/kernel.zig");
const axpy32 = @import("../blas1/axpy/x86/kernel.zig");

pub const Core = coretypes.X86CoreType;

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch == .x86_64 or arch == .x86;
}

/// Whether any routine has a kernel for `core` on `arch`.
pub fn hasKernels(comptime arch: std.Target.Cpu.Arch, comptime core: Core) bool {
    return switch (arch) {
        .x86_64 => axpy.has(core, "daxpy") or axpy.has(core, "saxpy") or
            axpy.has(core, "caxpy") or axpy.has(core, "zaxpy"),
        .x86 => axpy32.has(core),
        else => false,
    };
}

/// Exports every kernel `core` has, for the architecture being compiled.
pub fn exportCore(comptime core: Core) void {
    switch (sysinfo.cpu.arch) {
        .x86_64 => axpy.exportCore(core),
        .x86 => axpy32.exportCore(core),
        else => {},
    }
}
