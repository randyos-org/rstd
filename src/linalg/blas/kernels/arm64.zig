//! BLAS kernel registry for arm64: what `build/blaskern.zig` and `_kernel.zig`
//! ask of a registry. See `x86.zig` for what one provides.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const axpy = @import("../blas1/axpy/arm64/kernel.zig");
const amax = @import("../blas1/amax/arm64/kernel.zig");
const amin = @import("../blas1/amin/arm64/kernel.zig");
const asum = @import("../blas1/asum/arm64/kernel.zig");
const axpby = @import("../blas1/axpby/arm64/kernel.zig");
const copy = @import("../blas1/copy/arm64/kernel.zig");
const dot = @import("../blas1/dot/arm64/kernel.zig");
const iamax = @import("../blas1/iamax/arm64/kernel.zig");
const iamin = @import("../blas1/iamin/arm64/kernel.zig");
const imax = @import("../blas1/imax/arm64/kernel.zig");
const imin = @import("../blas1/imin/arm64/kernel.zig");
const max = @import("../blas1/max/arm64/kernel.zig");
const min = @import("../blas1/min/arm64/kernel.zig");
const nrm2 = @import("../blas1/nrm2/arm64/kernel.zig");
const rot = @import("../blas1/rot/arm64/kernel.zig");
const rotm = @import("../blas1/rotm/arm64/kernel.zig");
const scal = @import("../blas1/scal/arm64/kernel.zig");
const sum = @import("../blas1/sum/arm64/kernel.zig");
const swap = @import("../blas1/swap/arm64/kernel.zig");

/// Every routine registered above, each with the uniform `has(core)` /
/// `exportCore(core)` pair `hasKernels`/`exportCore` need. Add a routine
/// here (not to those functions) to wire it in.
const modules = .{
    axpy, amax, amin, asum, axpby, copy, dot, iamax, iamin,
    imax, imin, max,  min,  nrm2,  rot,  rotm, scal, sum,  swap,
};

pub const Core = coretypes.Arm64CoreType;

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch == .aarch64;
}

/// Whether `core`'s kernels run on `cpu`: the SVE cores' need the extension.
pub fn runsOn(cpu: std.Target.Cpu, core: Core) bool {
    return switch (core) {
        .a64fx, .neoversev1 => cpu.has(.aarch64, .sve),
        else => true,
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
