//! BLAS kernel registry for power: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu/power.zig` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `modules` to wire it in.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const registry = @import("../asm_registry.zig");
const Vtable = @import("../Vtable.zig");

pub const Core = coretypes.PowerCoreType;

/// Every routine with kernels here: each routine's converted assembly (its
/// `asm/__root__.zig`), and `axpy`'s own module, which adds the kernels
/// upstream writes in C.
const modules = .{
    @import("../blas1/axpy/power/kernel.zig"),
    registry.Kernels(Core, @import("../blas1/amax/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/amin/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/asum/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/copy/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/dot/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/iamax/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/iamin/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/imax/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/imin/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/max/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/min/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/nrm2/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/rot/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/scal/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/sum/power/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/swap/power/asm/__root__.zig")),
};

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch.isPowerPC();
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
