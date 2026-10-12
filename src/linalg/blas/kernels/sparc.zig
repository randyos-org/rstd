//! BLAS kernel registry for SPARC: what `build/blaskern.zig` and `_kernel.zig`
//! ask of a registry. See `x86.zig` for what one provides.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const axpy = @import("../blas1/axpy/sparc/kernel.zig");
const amax = @import("../blas1/amax/sparc/kernel.zig");
const amin = @import("../blas1/amin/sparc/kernel.zig");
const asum = @import("../blas1/asum/sparc/kernel.zig");
const axpby = @import("../blas1/axpby/sparc/kernel.zig");
const copy = @import("../blas1/copy/sparc/kernel.zig");
const dot = @import("../blas1/dot/sparc/kernel.zig");
const iamax = @import("../blas1/iamax/sparc/kernel.zig");
const iamin = @import("../blas1/iamin/sparc/kernel.zig");
const imax = @import("../blas1/imax/sparc/kernel.zig");
const imin = @import("../blas1/imin/sparc/kernel.zig");
const max = @import("../blas1/max/sparc/kernel.zig");
const min = @import("../blas1/min/sparc/kernel.zig");
const nrm2 = @import("../blas1/nrm2/sparc/kernel.zig");
const rot = @import("../blas1/rot/sparc/kernel.zig");
const rotm = @import("../blas1/rotm/sparc/kernel.zig");
const scal = @import("../blas1/scal/sparc/kernel.zig");
const sum = @import("../blas1/sum/sparc/kernel.zig");
const swap = @import("../blas1/swap/sparc/kernel.zig");

/// Every routine registered above, each with the uniform `has(core)` /
/// `exportCore(core)` pair `hasKernels`/`exportCore` need. Add a routine
/// here (not to those functions) to wire it in.
const modules = .{
    axpy, amax, amin, asum, axpby, copy, dot, iamax, iamin,
    imax, imin, max,  min,  nrm2,  rot,  rotm, scal, sum,  swap,
};

pub const Core = coretypes.SparcCoreType;

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch.isSPARC();
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
    inline for (modules) |m| {
        if (comptime m.has(core)) m.exportCore(core);
    }
}
