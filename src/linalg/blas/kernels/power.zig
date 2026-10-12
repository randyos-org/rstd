//! BLAS kernel registry for POWER: what `build/blaskern.zig` and `_kernel.zig`
//! ask of a registry. See `x86.zig` for what one provides.

const std = @import("std");
const sysinfo = @import("builtin");
const coretypes = @import("../cpu/_coretypes.zig");
const axpy = @import("../blas1/axpy/power/kernel.zig");
const amax = @import("../blas1/amax/power/kernel.zig");
const amin = @import("../blas1/amin/power/kernel.zig");
const asum = @import("../blas1/asum/power/kernel.zig");
const axpby = @import("../blas1/axpby/power/kernel.zig");
const copy = @import("../blas1/copy/power/kernel.zig");
const dot = @import("../blas1/dot/power/kernel.zig");
const iamax = @import("../blas1/iamax/power/kernel.zig");
const iamin = @import("../blas1/iamin/power/kernel.zig");
const imax = @import("../blas1/imax/power/kernel.zig");
const imin = @import("../blas1/imin/power/kernel.zig");
const max = @import("../blas1/max/power/kernel.zig");
const min = @import("../blas1/min/power/kernel.zig");
const nrm2 = @import("../blas1/nrm2/power/kernel.zig");
const rot = @import("../blas1/rot/power/kernel.zig");
const rotm = @import("../blas1/rotm/power/kernel.zig");
const scal = @import("../blas1/scal/power/kernel.zig");
const sum = @import("../blas1/sum/power/kernel.zig");
const swap = @import("../blas1/swap/power/kernel.zig");

/// Every routine registered above but `axpy`, each with the uniform
/// `has(core)` / `exportCore(core)` pair. `axpy`'s take the target `arch`
/// as well (its 32-bit-only cores have no 64-bit kernels), so it is called
/// on its own. Add a routine here to wire it in.
const modules = .{
    amax, amin, asum, axpby, copy, dot, iamax, iamin, imax,
    imin, max,  min,  nrm2,  rot,  rotm, scal, sum,   swap,
};

pub const Core = coretypes.PowerCoreType;

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch.isPowerPC();
}

/// Whether `core`'s kernels run on `cpu`: the 32-bit cores have no 64-bit code,
/// and the vector ones need their instructions.
pub fn runsOn(cpu: std.Target.Cpu, core: Core) bool {
    return switch (core) {
        .power8, .power9 => cpu.has(.powerpc, .vsx),
        .power10 => cpu.has(.powerpc, .vsx) and cpu.has(.powerpc, .isa_v31_instructions),
        .ppc440, .ppc440fp2, .ppcg4 => !cpu.arch.isPowerPC64(),
        else => true,
    };
}

/// Whether any routine has a kernel for `core` on `arch`.
pub fn hasKernels(comptime arch: std.Target.Cpu.Arch, comptime core: Core) bool {
    if (comptime axpy.has(arch, core)) return true;
    inline for (modules) |m| {
        if (comptime m.has(core)) return true;
    }
    return false;
}

/// Exports every kernel `core` has, for the architecture being compiled.
pub fn exportCore(comptime core: Core) void {
    axpy.exportCore(sysinfo.cpu.arch, core);
    inline for (modules) |m| {
        if (comptime m.has(core)) m.exportCore(core);
    }
}
