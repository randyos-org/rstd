//! BLAS kernel registry for riscv: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu/riscv.zig` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `modules` to wire it in.

const std = @import("std");
const coretypes = @import("../cpu/_coretypes.zig");
const registry = @import("../asm_registry.zig");
const Vtable = @import("../Vtable.zig");
const vector = @import("../vector.zig");

pub const Core = coretypes.RiscvCoreType;

/// Every routine with kernels here, all upstream's vector intrinsics: `axpy`
/// in `asm`, the rest `@Vector` ports (`vector.ported`), and the gemm tile.
const converted = .{
    @import("../blas1/axpy/riscv64/kernel.zig"),
    @import("../blas1/axpby/riscv64/kernel.zig"),
    @import("../blas1/dot/riscv64/kernel.zig"),
    @import("../blas1/rotm/riscv64/kernel.zig"),
    @import("../blas2/gemv/riscv64/kernel.zig"),
    @import("../blas2/hemv/riscv64/kernel.zig"),
    @import("../blas2/symv/riscv64/kernel.zig"),
    @import("../blas3/gemm/kernel.zig").Registry(Core),
};

/// Where the `@Vector` kernels (`vector.zig`) go: on the vector cores, the
/// `routed` slots and every other slot of theirs no port fills. The cores
/// without the vector extension (`riscv64_generic`, `u74`, `c910v`, whose
/// is pre-ratification) have no object.
const vector_rule: vector.Rule(Core) = .{ .unfilled_on = &.{ .riscv64_zvl128b, .riscv64_zvl256b, .x280 } };

/// The `@Vector` kernels go over those.
const modules = converted ++ .{vector.Registry(Core, vector_rule, converted)};

/// Whether this registry has kernels for `arch`.
pub fn supports(arch: std.Target.Cpu.Arch) bool {
    return arch == .riscv64;
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
