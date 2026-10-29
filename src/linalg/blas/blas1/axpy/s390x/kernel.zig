//! s390x axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`), and how they get installed into
//! the vtable (`install`): cores without an entry keep the portable kernels
//! from `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`,
//! not here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! The reference's `zarch/{s,d,c,z}axpy.c` are the shared C driver around a
//! block that is inline assembly on the vector facility: the block is in
//! `uk/`, as written, and composed with the driver in `../kernel.zig`.
//! `uk/__root__.zig` says which core uses which. The cores whose
//! `KERNEL.<CORE>` names the plain C loops keep the portable kernels.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const Vtable = @import("../../../Vtable.zig");
const util = @import("../../../util.zig");
const kernel = @import("../kernel.zig");
const uk = @import("../../uk_block.zig");
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
