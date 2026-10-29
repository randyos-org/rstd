//! BLAS kernel registry for x86_64: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu/x86_64.zig` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `converted` to wire it in.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const registry = @import("../asm_registry.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");
const vector = @import("../vector.zig");

const Slot = util.Slot;

pub const Core = coretypes.X86CoreType;

/// Every routine with kernels here: each routine's converted assembly (its
/// `asm/__root__.zig`), and for `axpy` and `dot` their own modules, which add
/// the kernels upstream writes in C.
const converted = .{
    @import("../blas1/axpy/x86_64/kernel.zig"),
    registry.Kernels(Core, @import("../blas1/asum/x86_64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/copy/x86_64/asm/__root__.zig")),
    @import("../blas1/dot/x86_64/kernel.zig"),
    @import("../blas1/rot/x86_64/kernel.zig"),
    @import("../blas1/scal/x86_64/kernel.zig"),
    registry.Kernels(Core, @import("../blas1/sum/x86_64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/swap/x86_64/asm/__root__.zig")),
    @import("../blas2/gemv/x86_64/kernel.zig"),
    registry.Kernels(Core, @import("../blas2/hemv/x86_64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas2/symv/x86_64/asm/__root__.zig")),
    @import("../blas3/gemm/kernel.zig").Registry(Core),
};

/// The `@Vector` kernels go over those.
const modules = converted ++ .{vector.Registry(Core, usesVector)};

/// Where the `@Vector` kernels (`vector.zig`) go: the `routed` slots, those
/// routed for their alignment (`vector.peels`) only on the cores whose
/// upstream assembly they replace, which is not converted (`skips`).
fn usesVector(comptime core: Core, comptime slot: Slot) bool {
    if (!vector.routed(slot)) return false;
    if (!vector.peels(slot)) return true;
    inline for (converted) |m| {
        if (@hasDecl(m, "skips") and m.skips(core, slot)) return true;
    }
    return false;
}

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
    inline for (modules) |m| m.exportCore(core);
}

/// Points the vtable at the kernels `core` has (see `asm_registry.install`).
pub fn install(vt: *Vtable, core: Core) void {
    registry.install(vt, modules, core);
}

test "the @Vector kernels take every slot whose peeling assembly was left out" {
    // the routines' own indexes, so a module that forgets `skips` is caught
    const peeling_asm = .{
        registry.Kernels(Core, @import("../blas1/asum/x86_64/asm/__root__.zig")),
        registry.Kernels(Core, @import("../blas1/sum/x86_64/asm/__root__.zig")),
        registry.Kernels(Core, @import("../blas1/dot/x86_64/asm/__root__.zig")),
        registry.Kernels(Core, @import("../blas2/gemv/x86_64/asm/__root__.zig")),
    };
    @setEvalBranchQuota(100_000);
    inline for (comptime std.meta.tags(Core)) |core| {
        inline for (vector.peeling) |slot| {
            const left_out = comptime for (peeling_asm) |a| {
                if (a.skips(core, slot)) break true;
            } else false;
            if (left_out and !usesVector(core, slot)) {
                std.debug.print("{t} {t}\n", .{ core, slot });
                return error.TestUnexpectedResult;
            }
        }
    }
}
