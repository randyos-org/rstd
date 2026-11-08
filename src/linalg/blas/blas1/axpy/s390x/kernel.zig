//! s390x axpy kernels per core (`zarch/{s,d,c,z}axpy.c`: `../driver.zig` around a `uk/` block); tests in `../__root__.zig`.

const coretypes = @import("../../../cpu/coretypes.zig");
const registry = @import("../../../registry.zig");
const Vtable = @import("../../../Vtable.zig");
const util = @import("../../../util.zig");
const kernel = @import("../driver.zig");
const uk = @import("../../../uk_block.zig");
const ukidx = @import("uk/__root__.zig");

const Slot = util.Slot;

pub const Core = coretypes.S390xCoreType;

/// The entry for `slot` on `core`, if it has one.
fn entry(comptime slot: Slot, comptime core: Core) ?ukidx.Entry {
    return ukidx.find(slot, @tagName(core));
}

/// Whether `slot` has a kernel of its own on `core`.
pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
    return entry(slot, core) != null;
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
        if (comptime entry(slot, core)) |e| {
            registry.exportKernel(uk.driven(kernel, slot, e, uk.NoVector), slot, core);
        }
    }
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installSlots(vt, core, &kernel.slots, hasSlot);
}
