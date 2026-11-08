//! BLAS kernel registry for arm64: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu.Arch` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `converted` to wire it in.

const coretypes = @import("../cpu/coretypes.zig");
const registry = @import("../registry.zig");

pub const Core = coretypes.Arm64CoreType;

/// Every routine with kernels here: each routine's converted assembly (its
/// `asm/__root__.zig`), and a routine's own module (`kernel.zig`) where it
/// adds the kernels upstream writes in C or intrinsics.
const converted = .{
    @import("../blas1/axpy/arm64/kernel.zig"),
    registry.Kernels(Core, @import("../blas1/asum/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/copy/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/dot/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/rot/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/scal/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/sum/arm64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/swap/arm64/asm/__root__.zig")),
    @import("../blas2/gemv/arm64/kernel.zig"),
    @import("../blas2/symv/arm64/kernel.zig"),
    @import("../blas3/kernels.zig").Registry(Core),
};

/// Where the `@Vector` kernels (`vector.zig`) go: the `routed` slots, and
/// off Apple each other slot of theirs no converted kernel fills or that
/// they displace on `vortex` and `vortexm4` (`registry.displacing`), whose
/// objects then leave the converted kernels out.
const vector_rule: registry.Rule(Core) = .unfilled;

/// The `@Vector` kernels go over those.
const modules = converted ++ .{registry.Vector(Core, vector_rule, converted)};

/// The registry (`registry.Arch`).
pub const Registry = registry.Arch(Core, &.{.aarch64}, modules, false);
