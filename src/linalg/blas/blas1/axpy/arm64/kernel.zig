//! AArch64 axpy kernels, one set per core, built into the per-core kernel objects
//! (`src/linalg/blas/_kernel.zig`), and how they get installed into the
//! vtable (`install`): cores without an entry keep the portable kernels from
//! `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`, not
//! here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! The kernels that are assembly upstream are in `asm/`, one Zig file per
//! source and precision with its `#if`s and macros resolved (`asm/__root__.zig`
//! says which core uses which). `A64FX` and `NEOVERSEV1` use `axpy_sve.c`
//! (`sve.zig`); `THUNDERX` uses `daxpy_thunderx.c`, its block (`thunderx.zig`)
//! under the shared driver. Other cores keep the portable kernels.

const sysinfo = @import("builtin");
const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const util = @import("../../../util.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../kernel.zig");
const uk = @import("../../uk_block.zig");
const sve = @import("sve.zig");
const thunderx = @import("thunderx.zig");

const Slot = util.Slot;
const R = registry.Kernels(coretypes.Arm64CoreType, @import("asm/__root__.zig"));

pub const Core = coretypes.Arm64CoreType;

const thunderx_entry = .{
    .mask = 32,
    .block = uk.Block.assembly,
    .kernel = "daxpy_kernel_8",
    .module = thunderx,
};

fn usesSve(comptime core: Core, comptime slot: Slot) bool {
    return (core == .a64fx or core == .neoversev1) and !kernel.form(slot).complex;
}

fn usesThunderx(comptime core: Core, comptime slot: Slot) bool {
    return core == .thunderx and slot == .daxpy_k;
}

/// Whether `slot` has a kernel of its own on `core`.
pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
    return R.hasSlot(core, slot) or usesSve(core, slot) or usesThunderx(core, slot);
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
        if (comptime R.hasSlot(core, slot)) {
            R.exportOne(slot, core);
        } else if (comptime usesSve(core, slot)) {
            registry.exportKernel(&sve.Real(kernel.form(slot).F).call, slot, core);
        } else if (comptime usesThunderx(core, slot)) {
            registry.exportKernel(uk.driven(kernel, slot, thunderx_entry, uk.NoVector), slot, core);
        }
    }
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installSlots(vt, core, &kernel.slots, hasSlot);
}
