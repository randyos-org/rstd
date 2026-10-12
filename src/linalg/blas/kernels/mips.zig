//! BLAS kernel registry for mips: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu/mips.zig` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `modules` to wire it in.

const std = @import("std");
const sysinfo = @import("builtin");
const coretypes = @import("../cpu/_coretypes.zig");
const registry = @import("../asm_registry.zig");
const Vtable = @import("../Vtable.zig");

pub const Core = coretypes.MipsCoreType;

/// Every routine with kernels here: each routine's converted assembly (its
/// `asm/__root__.zig`), and `axpy`'s own module, which adds the kernels
/// upstream writes in C.
const modules = .{
    @import("../blas1/axpy/mips64/kernel.zig"),
    registry.Kernels(Core, @import("../blas1/amax/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/amin/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/asum/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/copy/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/dot/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/iamax/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/iamin/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/imax/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/imin/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/max/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/min/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/nrm2/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/scal/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/sum/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/swap/mips64/asm/__root__.zig")),
};

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
    inline for (modules) |m| m.exportCore(core);
}

/// Points the vtable at the kernels `core` has; the portable ones stay where
/// it has none, and for a core whose object is not linked for this target.
pub fn install(vt: *Vtable, core: Core) void {
    switch (core) {
        inline else => |c| {
            if (comptime !runsOn(sysinfo.cpu, c)) return;
            inline for (modules) |m| m.install(vt, c);
        },
    }
}
