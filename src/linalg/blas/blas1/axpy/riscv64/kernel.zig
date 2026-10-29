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
const kernel = @import("../kernel.zig");
const rvv = @import("rvv.zig");

pub const Core = coretypes.RiscvCoreType;

/// The hand ports and where they serve.
pub const ports = [_]registry.Port(Core){
    .{ .module = kernel.Ported(rvv), .cores = &.{ .riscv64_zvl128b, .x280 }, .slots = &.{ .saxpy_k, .daxpy_k, .caxpy_k, .caxpyc_k, .zaxpy_k, .zaxpyc_k }, .written = .assembly },
};

const P = registry.Ports(Core, "riscv64", &ports);

pub const hasSlot = P.hasSlot;
pub const has = P.has;

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    inline for (kernel.slots) |slot| P.exportOne(slot, core);
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installSlots(vt, core, &kernel.slots, hasSlot);
}
