//! POWER axpy kernels, one set per core, built into the per-core kernel objects
//! (`src/linalg/blas/_kernel.zig`), and how they get installed into the
//! vtable (`install`): cores without an entry keep the portable kernels from
//! `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`, not
//! here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! The kernels that are assembly upstream are in `asm/`, one Zig file per
//! source and precision with its `#if`s and macros resolved (`asm/__root__.zig`
//! says which core uses which). `POWER8` and up use the C kernels of
//! `c_kernels.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const util = @import("../../../util.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../kernel.zig");
const c_kernels = @import("c_kernels.zig");

const Slot = util.Slot;
const R = registry.Kernels(coretypes.PowerCoreType, @import("asm/__root__.zig"));

pub const Core = coretypes.PowerCoreType;

/// Which C kernels `core` uses (`KERNEL.POWER8`, `.POWER9`, `.POWER10`); the
/// assembly `caxpy` alternatives of the first two are for older compilers.
fn cSource(comptime core: Core) ?c_kernels.Source {
    return switch (core) {
        .power8, .power9 => .power8,
        .power10 => .power10,
        else => null,
    };
}

/// Exports every axpy kernel `core` has.
pub fn exportCore(comptime core: Core) void {
    inline for (kernel.slots) |slot| {
        if (comptime cSource(core)) |source| {
            registry.exportKernel(c_kernels.entry(slot, source), slot, core);
        } else if (comptime R.hasSlot(core, slot)) {
            R.exportOne(slot, core);
        }
    }
}

/// Whether `core` has any axpy kernel. (The 32-bit cores have no 64-bit
/// code; `kernels/power.zig`'s `runsOn` keeps them off 64-bit targets.)
pub fn has(comptime core: Core) bool {
    return cSource(core) != null or R.has(core);
}

/// Whether `slot` has a kernel on `core` for this target.
pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
    return cSource(core) != null or R.hasSlot(core, slot);
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installSlots(vt, core, &kernel.slots, hasSlot);
}
