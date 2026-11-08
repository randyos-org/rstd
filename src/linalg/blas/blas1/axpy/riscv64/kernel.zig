//! RISC-V axpy kernels per core (`rvv.zig`: `axpy_rvv.c`, `zaxpy_rvv.c` as the assembly their intrinsics stand for); tests in `../__root__.zig`.

const coretypes = @import("../../../cpu/coretypes.zig");
const registry = @import("../../../registry.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../driver.zig");
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
