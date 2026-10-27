//! BLAS kernel registry for arm64: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu/arm64.zig` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `converted` to wire it in.

const std = @import("std");
const sysinfo = @import("builtin");
const coretypes = @import("../cpu/_coretypes.zig");
const registry = @import("../asm_registry.zig");
const Vtable = @import("../Vtable.zig");
const vector = @import("../vector.zig");

pub const Core = coretypes.Arm64CoreType;

/// Every routine with kernels here: each routine's converted assembly (its
/// `asm/__root__.zig`), and `axpy`'s own module, which adds the kernels
/// upstream writes in C.
const converted = .{
    @import("../blas1/axpy/arm64/kernel.zig"),
    registry.Kernels(Core, @import("../blas1/asum/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/copy/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/dot/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/rot/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/scal/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/sum/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/swap/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas2/gemv/arm64/asm/__root__.zig")),
    @import("../blas3/gemm/kernel.zig").Registry(Core),
};

/// The `@Vector` kernels go over those.
const modules = converted ++ .{vector.Registry(Core, usesVector)};

/// Where the `@Vector` kernels (`vector.zig`) go: every slot of theirs but
/// `rot` on `vortex` and `vortexm4`, where they match or beat the converted
/// kernels (`tools/blasbench/results-m6.md`), and on every core each slot no
/// converted kernel fills, mostly C upstream, for the portable scalar loop.
/// On Apple's systems, whose Accelerate is to serve the rest, only the
/// `routed` slots.
fn usesVector(comptime core: Core, comptime role: []const u8) bool {
    if (vector.routed(role)) return true;
    if (sysinfo.target.os.tag.isDarwin()) return false;
    if ((core == .vortex or core == .vortexm4) and !std.mem.endsWith(u8, role, "rot")) return true;
    inline for (converted) |m| {
        if (@hasDecl(m, "hasRole") and m.hasRole(core, role)) return false;
    }
    return true;
}

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch == .aarch64;
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
