//! BLAS kernel registry for arm: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu.Arch` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `converted` to wire it in.

const coretypes = @import("../cpu/coretypes.zig");
const registry = @import("../registry.zig");

pub const Core = coretypes.ArmCoreType;

/// Every routine with kernels here: each routine's converted assembly (its
/// `asm/__root__.zig`), and `axpy`'s own module, which adds the kernels
/// upstream writes in C.
const converted = .{
    @import("../blas1/axpy/arm/kernel.zig"),
    registry.Kernels(Core, @import("../blas1/asum/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/dot/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/rot/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/sum/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/swap/arm/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas2/gemv/arm/asm/__root__.zig")),
    @import("../blas3/kernels.zig").Registry(Core),
};

/// Where the `@Vector` kernels (`vector.zig`) go: the `routed` slots, but
/// `s` and `c` on the NEON cores, whose lanes flush `f32` (`registry.flushing`).
const vector_rule: registry.Rule(Core) = .{ .routed = null };

/// The `@Vector` kernels go over those.
const modules = converted ++ .{registry.Vector(Core, vector_rule, converted)};

/// The registry (`registry.Arch`).
pub const Registry = registry.Arch(Core, &.{.arm}, modules, false);
