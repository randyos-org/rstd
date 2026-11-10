//! BLAS kernel registry for loongarch: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu.Arch` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `converted` to wire it in.

const coretypes = @import("../cpu/coretypes.zig");
const registry = @import("../registry.zig");

pub const Core = coretypes.LoongarchCoreType;

/// Every routine with kernels here: each routine's converted assembly (its
/// `asm/__root__.zig`).
const converted = .{
    registry.Kernels(Core, @import("../blas1/axpy/loongarch64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/asum/loongarch64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/axpby/loongarch64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/copy/loongarch64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/dot/loongarch64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/rot/loongarch64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/scal/loongarch64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/sum/loongarch64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/swap/loongarch64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas2/gemv/loongarch64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas2/symv/loongarch64/asm/__root__.zig")),
    @import("../blas3/kernels.zig").Registry(Core),
};

/// Where the `@Vector` kernels (`vector.zig`) go: the `routed` slots.
const vector_rule: registry.Rule(Core) = .{ .routed = null };

/// The `@Vector` kernels go over those.
const modules = converted ++ .{registry.Vector(Core, vector_rule, converted)};

/// The registry (`registry.Arch`).
pub const Registry = registry.Arch(Core, &.{.loongarch64}, modules, false);
