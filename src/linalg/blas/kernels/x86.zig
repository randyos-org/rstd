//! BLAS kernel registry for x86: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu/x86.zig` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `converted` to wire it in.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const corespec = @import("../cpu/_corespec.zig");
const registry = @import("../asm_registry.zig");
const Vtable = @import("../Vtable.zig");
const vector = @import("../vector.zig");

pub const Core = coretypes.X86CoreType;

/// Every routine with kernels here: each routine's converted assembly (its
/// `asm/__root__.zig`), and `axpy`'s own module, which adds the kernels
/// upstream writes in C.
const converted = .{
    @import("../blas1/axpy/x86/kernel.zig"),
    registry.Kernels(Core, @import("../blas1/asum/x86/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/copy/x86/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/dot/x86/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/rot/x86/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/scal/x86/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/sum/x86/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/swap/x86/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas2/gemv/x86/asm/__root__.zig")),
    @import("../blas3/gemm/kernel.zig").Registry(Core),
};

/// Where the `@Vector` kernels (`vector.zig`) go: the `routed` slots.
const vector_rule: vector.Rule(Core) = .{ .routed = null };

/// The `@Vector` kernels go over those.
const modules = converted ++ .{vector.Registry(Core, vector_rule, converted)};

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch == .x86;
}

/// Whether any routine has a kernel for `core`; none for a core another's
/// object serves (`corespec.objectCore`).
pub fn hasKernels(comptime arch: std.Target.Cpu.Arch, comptime core: Core) bool {
    if (corespec.objectCore(arch, core) != core) return false;
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
