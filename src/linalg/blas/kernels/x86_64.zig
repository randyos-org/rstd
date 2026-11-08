//! BLAS kernel registry for x86_64: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu.Arch` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `converted` to wire it in.

const std = @import("std");
const coretypes = @import("../cpu/coretypes.zig");
const registry = @import("../registry.zig");
const vector = @import("../vector.zig");

pub const Core = coretypes.X86CoreType;

/// Every routine with kernels here: each routine's converted assembly (its
/// `asm/__root__.zig`), and a routine's own module (`kernel.zig`) where it
/// adds the kernels upstream writes in C or intrinsics.
const converted = .{
    @import("../blas1/axpy/x86_64/kernel.zig"),
    @import("../blas1/asum/x86_64/kernel.zig"),
    registry.Kernels(Core, @import("../blas1/copy/x86_64/asm/__root__.zig")),
    @import("../blas1/dot/x86_64/kernel.zig"),
    @import("../blas1/rot/x86_64/kernel.zig"),
    @import("../blas1/scal/x86_64/kernel.zig"),
    registry.Kernels(Core, @import("../blas1/sum/x86_64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/swap/x86_64/asm/__root__.zig")),
    @import("../blas2/gemv/x86_64/kernel.zig"),
    @import("../blas2/ger/x86_64/kernel.zig"),
    registry.Kernels(Core, @import("../blas2/hemv/x86_64/asm/__root__.zig")),
    @import("../blas2/symv/x86_64/kernel.zig"),
    @import("../blas3/kernels.zig").Registry(Core),
    @import("../blas3/gemm/small_x86_64.zig"),
};

/// Where the `@Vector` kernels (`vector.zig`) go: the `routed` slots, those
/// routed for their alignment (`vector.peels`) only on the cores whose
/// upstream assembly they replace, which is not converted (`skips`).
const vector_rule: registry.Rule(Core) = .{ .routed = null };

const Vector = registry.Vector(Core, vector_rule, converted);

/// The `@Vector` kernels go over those.
const modules = converted ++ .{Vector};

/// The registry (`registry.Arch`).
pub const Registry = registry.Arch(Core, &.{.x86_64}, modules, false);

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
            if (left_out and !Vector.uses(core, slot)) {
                std.debug.print("{t} {t}\n", .{ core, slot });
                return error.TestUnexpectedResult;
            }
        }
    }
}
