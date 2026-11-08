//! 32-bit x86 axpy kernels per core (`asm/`; C ones stay portable until ported); tests in `../__root__.zig`.

const coretypes = @import("../../../cpu/coretypes.zig");
const registry = @import("../../../registry.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../driver.zig");

const R = registry.Kernels(coretypes.X86CoreType, @import("asm/__root__.zig"));

pub const Core = coretypes.X86CoreType;
pub const has = R.has;
pub const hasSlot = R.hasSlot;
pub const exportCore = R.exportCore;

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installSlots(vt, core, &kernel.slots, hasSlot);
}
