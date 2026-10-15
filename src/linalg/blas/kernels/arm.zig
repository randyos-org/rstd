//! BLAS kernel registry for arm: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu/arm.zig` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `modules` to wire it in.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const registry = @import("../asm_registry.zig");
const Vtable = @import("../Vtable.zig");

pub const Core = coretypes.ArmCoreType;

/// Every routine with kernels here: each routine's converted assembly (its
/// `asm/__root__.zig`), and `axpy`'s own module, which adds the kernels
/// upstream writes in C.
const modules = .{
    @import("../blas1/axpy/arm/kernel.zig"),
    registry.Kernels(Core, @import("../blas1/amax/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/amin/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/asum/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/dot/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/iamax/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/iamin/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/imax/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/imin/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/max/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/min/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/nrm2/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/rot/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/sum/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/swap/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas2/gemv/arm/asm/__root__.zig")),
};

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch == .arm;
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

/// Points the vtable at the kernels `core` has (see `asm_registry.install`).
pub fn install(vt: *Vtable, core: Core) void {
    registry.install(vt, modules, core);
}
