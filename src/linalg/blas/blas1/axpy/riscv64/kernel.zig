//! RISC-V axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`), and how they get installed into
//! the vtable (`install`): cores without an entry keep the portable kernels
//! from `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`,
//! not here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! The reference's kernels here are C over the vector extension's
//! intrinsics, `axpy_rvv.c` and `zaxpy_rvv.c` for the cores whose
//! `KERNEL.<CORE>` names them; they are in `rvv.zig` as the assembly the
//! intrinsics stand for. The cores that use the plain C loops
//! (`riscv64_generic`, `u74`) keep the portable kernels.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const Vtable = @import("../../../Vtable.zig");
const util = @import("../../../util.zig");
const kernel = @import("../kernel.zig");
const rvv = @import("rvv.zig");

const Slot = util.Slot;

pub const Core = coretypes.RiscvCoreType;

/// Whether `core` uses the `_rvv` kernels.
fn usesRvv(comptime core: Core) bool {
    return switch (core) {
        .riscv64_zvl128b, .x280 => true,
        else => false,
    };
}

/// Whether `slot` has a kernel of its own on `core`.
pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
    _ = slot;
    return usesRvv(core);
}

/// Whether `core` has any axpy kernel of its own.
pub fn has(comptime core: Core) bool {
    return usesRvv(core);
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (comptime !usesRvv(core)) return;
    registry.exportKernel(&rvv.Real(f32).call, .saxpy_k, core);
    registry.exportKernel(&rvv.Real(f64).call, .daxpy_k, core);
    registry.exportKernel(&rvv.Cplx(f32, false).call, .caxpy_k, core);
    registry.exportKernel(&rvv.Cplx(f32, true).call, .caxpyc_k, core);
    registry.exportKernel(&rvv.Cplx(f64, false).call, .zaxpy_k, core);
    registry.exportKernel(&rvv.Cplx(f64, true).call, .zaxpyc_k, core);
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installSlots(vt, core, &kernel.slots, hasSlot);
}
