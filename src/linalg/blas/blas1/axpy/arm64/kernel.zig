//! AArch64 axpy kernels per core (`asm/`, `sve.zig` for A64FX and Neoverse V1, `thunderx.zig`'s block); tests in `../__root__.zig`.

const coretypes = @import("../../../cpu/coretypes.zig");
const registry = @import("../../../registry.zig");
const util = @import("../../../util.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../driver.zig");
const uk = @import("../../../uk_block.zig");
const sve = @import("sve.zig");
const thunderx = @import("thunderx.zig");

const Slot = util.Slot;
const R = registry.Kernels(coretypes.Arm64CoreType, @import("asm/__root__.zig"));

pub const Core = coretypes.Arm64CoreType;

/// `daxpy_thunderx.c`'s block (`thunderx.zig`) under the shared driver.
const Thunderx = struct {
    pub fn entry(comptime slot: Slot) kernel.Kernel(slot) {
        return uk.driven(kernel, slot, .{
            .mask = 32,
            .block = uk.Block.assembly,
            .kernel = "daxpy_kernel_8",
            .module = thunderx,
        }, uk.NoVector);
    }
};

/// The hand ports and where they serve.
pub const ports = [_]registry.Port(Core){
    .{ .module = kernel.Ported(sve), .cores = &.{ .a64fx, .neoversev1 }, .slots = &.{ .saxpy_k, .daxpy_k }, .written = .assembly },
    .{ .module = Thunderx, .cores = &.{.thunderx}, .slots = &.{.daxpy_k}, .written = .scalar },
};

const P = registry.Ports(Core, "arm64", &ports);

/// Whether `slot` has a kernel of its own on `core`.
pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
    return P.hasSlot(core, slot) or R.hasSlot(core, slot);
}

/// Whether `core` has any axpy kernel of its own.
pub fn has(comptime core: Core) bool {
    inline for (kernel.slots) |slot| {
        if (comptime hasSlot(core, slot)) return true;
    }
    return false;
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    inline for (kernel.slots) |slot| {
        if (comptime P.hasSlot(core, slot)) {
            P.exportOne(slot, core);
        } else if (comptime R.hasSlot(core, slot)) {
            R.exportOne(slot, core);
        }
    }
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installSlots(vt, core, &kernel.slots, hasSlot);
}
