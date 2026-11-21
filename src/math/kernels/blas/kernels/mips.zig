//! BLAS kernel registry for mips: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu.Arch` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `converted` to wire it in.

const coretypes = @import("../cpu/coretypes.zig");
const registry = @import("../registry.zig");

pub const Core = coretypes.MipsCoreType;

/// Every routine with kernels here: each routine's converted assembly (its
/// `asm/__root__.zig`), and `axpy`'s own module, which adds the kernels
/// upstream writes in C.
const converted = .{
    @import("../blas1/axpy/mips64/kernel.zig"),
    registry.Kernels(Core, @import("../blas1/asum/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/copy/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/dot/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/scal/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/sum/mips64/asm/__root__.zig")),
    registry.Kernels(Core, @import("../blas1/swap/mips64/asm/__root__.zig")),
    @import("../blas3/kernels.zig").Registry(Core),
};

/// Where the `@Vector` kernels (`vector.zig`) go: the `routed` slots.
const vector_rule: registry.Rule(Core) = .{ .routed = null };

/// The `@Vector` kernels go over those.
const modules = converted ++ .{registry.Vector(Core, vector_rule, converted)};

/// The registry (`registry.Arch`).
pub const Registry = registry.Arch(Core, &.{ .mips, .mipsel, .mips64, .mips64el }, modules, false);
