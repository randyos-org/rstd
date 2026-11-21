//! BLAS kernel registry for riscv: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu.Arch` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `modules` to wire it in.

const coretypes = @import("../cpu/coretypes.zig");
const registry = @import("../registry.zig");

pub const Core = coretypes.RiscvCoreType;

/// Every routine with kernels here, all upstream's vector intrinsics: `axpy`
/// in `asm`, the rest `@Vector` ports (`vector.ported`), and the gemm tile.
const converted = .{
    @import("../blas1/axpy/riscv64/kernel.zig"),
    @import("../blas1/axpby/riscv64/kernel.zig"),
    @import("../blas1/rotm/riscv64/kernel.zig"),
    @import("../blas2/gemv/riscv64/kernel.zig"),
    @import("../blas2/hemv/riscv64/kernel.zig"),
    @import("../blas2/symv/riscv64/kernel.zig"),
    @import("../blas3/kernels.zig").Registry(Core),
};

/// Where the `@Vector` kernels (`vector.zig`) go: on the vector cores, the
/// `routed` slots and every other slot of theirs no port fills. The cores
/// without the vector extension (`riscv64_generic`, `u74`, `c910v`, whose
/// is pre-ratification) have no object.
const vector_rule: registry.Rule(Core) = .{ .unfilled_on = &.{ .riscv64_zvl128b, .riscv64_zvl256b, .x280 } };

/// The `@Vector` kernels go over those.
const modules = converted ++ .{registry.Vector(Core, vector_rule, converted)};

/// The registry (`registry.Arch`).
pub const Registry = registry.Arch(Core, &.{.riscv64}, modules, false);
