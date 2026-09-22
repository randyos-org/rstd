//! BLAS kernel registry for WebAssembly: what `build/blaskern.zig` and
//! `_kernel.zig` ask of a registry. See `x86.zig` for what one provides.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const axpy = @import("../blas1/axpy/wasm/kernel.zig");

pub const Core = coretypes.WasmCoreType;

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch.isWasm();
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
