//! BLAS kernel registry for MIPS: what `build/blaskern.zig` and `_kernel.zig`
//! ask of a registry. See `x86.zig` for what one provides.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const axpy = @import("../blas1/axpy/mips64/kernel.zig");
const amax = @import("../blas1/amax/mips64/kernel.zig");
const amin = @import("../blas1/amin/mips64/kernel.zig");
const asum = @import("../blas1/asum/mips64/kernel.zig");
const axpby = @import("../blas1/axpby/mips64/kernel.zig");
const copy = @import("../blas1/copy/mips64/kernel.zig");
const dot = @import("../blas1/dot/mips64/kernel.zig");
const iamax = @import("../blas1/iamax/mips64/kernel.zig");
const iamin = @import("../blas1/iamin/mips64/kernel.zig");
const imax = @import("../blas1/imax/mips64/kernel.zig");
const imin = @import("../blas1/imin/mips64/kernel.zig");
const max = @import("../blas1/max/mips64/kernel.zig");
const min = @import("../blas1/min/mips64/kernel.zig");
const nrm2 = @import("../blas1/nrm2/mips64/kernel.zig");
const rot = @import("../blas1/rot/mips64/kernel.zig");
const rotm = @import("../blas1/rotm/mips64/kernel.zig");
const scal = @import("../blas1/scal/mips64/kernel.zig");
const sum = @import("../blas1/sum/mips64/kernel.zig");
const swap = @import("../blas1/swap/mips64/kernel.zig");

/// Every routine registered above, each with the uniform `has(core)` /
/// `exportCore(core)` pair `hasKernels`/`exportCore` need. Add a routine
/// here (not to those functions) to wire it in.
const modules = .{
    axpy, amax, amin, asum, axpby, copy, dot, iamax, iamin,
    imax, imin, max,  min,  nrm2,  rot,  rotm, scal, sum,  swap,
};

pub const Core = coretypes.MipsCoreType;

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch.isMIPS();
}

/// Whether `core`'s kernels run on `cpu`: a core's code is for one word size,
/// the Loongson 3A kernels are release 2 code, which release 6 does not run,
/// and the SIMD kernels need the release that has the SIMD architecture.
pub fn runsOn(cpu: std.Target.Cpu, core: Core) bool {
    const wide = cpu.arch.isMIPS64();
    return switch (core) {
        .mips24k, .mips1004k, .p5600 => !wide and (core != .p5600 or cpu.has(.mips, .mips32r5)),
        .i6400, .i6500, .p6600 => wide and cpu.has(.mips, .mips64r6),
        .loongson3r3 => wide and !cpu.has(.mips, .mips64r6),
        // release 2 code with the SIMD architecture, whose objects do not link with release 6's
        .loongson3r4 => wide and !cpu.has(.mips, .mips64r6),
        .mips64_generic, .unknown => wide,
    };
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
