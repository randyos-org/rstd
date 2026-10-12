//! x86_64 BLAS kernel registry (`kernels/<arch>.zig`): which routines have
//! kernels for which cores, and how to export them. The per-core kernel
//! object root (`_kernel.zig`) exports through it, and the build
//! (`build/blaskern.zig`) asks it which cores need an object at all.
//!
//! 32-bit x86 is the separate `kernels/x86.zig`; the two share nothing
//! below their `CoreType`, which `cpu/_x86.zig`'s shared detector defines.
//!
//! To add a routine, import its `x86_64/kernel.zig` here and add it to
//! `hasKernels` and `exportCore`.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const axpy = @import("../blas1/axpy/x86_64/kernel.zig");
const amax = @import("../blas1/amax/x86_64/kernel.zig");
const amin = @import("../blas1/amin/x86_64/kernel.zig");
const asum = @import("../blas1/asum/x86_64/kernel.zig");
const axpby = @import("../blas1/axpby/x86_64/kernel.zig");
const copy = @import("../blas1/copy/x86_64/kernel.zig");
const dot = @import("../blas1/dot/x86_64/kernel.zig");
const iamax = @import("../blas1/iamax/x86_64/kernel.zig");
const iamin = @import("../blas1/iamin/x86_64/kernel.zig");
const imax = @import("../blas1/imax/x86_64/kernel.zig");
const imin = @import("../blas1/imin/x86_64/kernel.zig");
const max = @import("../blas1/max/x86_64/kernel.zig");
const min = @import("../blas1/min/x86_64/kernel.zig");
const nrm2 = @import("../blas1/nrm2/x86_64/kernel.zig");
const rot = @import("../blas1/rot/x86_64/kernel.zig");
const rotm = @import("../blas1/rotm/x86_64/kernel.zig");
const scal = @import("../blas1/scal/x86_64/kernel.zig");
const sum = @import("../blas1/sum/x86_64/kernel.zig");
const swap = @import("../blas1/swap/x86_64/kernel.zig");

/// Every routine registered above, each with the uniform `has(core)` /
/// `exportCore(core)` pair `hasKernels`/`exportCore` need. Add a routine
/// here (not to those functions) to wire it in.
const modules = .{
    axpy, amax, amin, asum, axpby, copy, dot, iamax, iamin,
    imax, imin, max,  min,  nrm2,  rot,  rotm, scal, sum,  swap,
};

pub const Core = coretypes.X86CoreType;

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch == .x86_64;
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
