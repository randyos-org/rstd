//! BLAS kernel registry for wasm: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu/wasm.zig` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `modules` to wire it in.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const registry = @import("../asm_registry.zig");
const Vtable = @import("../Vtable.zig");

pub const Core = coretypes.WasmCoreType;

/// Every routine with kernels here: each routine's converted assembly (its
/// `asm/__root__.zig`), and `axpy`'s own module, which adds the kernels
/// upstream writes in C.
const modules = .{
    @import("../blas1/axpy/wasm/kernel.zig"),
};

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch.isWasm();
}

/// Whether any routine has a kernel for `core`.
pub fn hasKernels(comptime arch: std.Target.Cpu.Arch, comptime core: Core) bool {
    _ = arch;
    inline for (modules) |m| {
        if (comptime m.has(core)) return true;
    }
    return false;
}

/// Exports every kernel `core` has.
pub fn exportCore(comptime core: Core) void {
    inline for (modules) |m| m.exportCore(core);
}

/// Points the vtable at the kernels `core` has; the portable ones stay where
/// it has none.
pub fn install(vt: *Vtable, core: Core) void {
    switch (core) {
        inline else => |c| {
            inline for (modules) |m| m.install(vt, c);
        },
    }
}
