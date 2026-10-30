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

/// The C kernels of `source` (`c_kernels.zig`) as a port.
fn CKernels(comptime source: c_kernels.Source) type {
    return struct {
        pub fn entry(comptime slot: Slot) kernel.Kernel(slot) {
            return c_kernels.entry(slot, source);
        }
    };
}

/// The hand ports and where they serve: the C kernels of `KERNEL.POWER8`,
/// `.POWER9` and `.POWER10`, whose blocks are converted inline assembly
/// (`uk/`) or AltiVec intrinsics on vectors. The assembly `caxpy`
/// alternatives of the first two are for older compilers.
pub const ports = [_]registry.Port(Core){
    .{ .module = CKernels(.power8), .cores = &.{ .power8, .power9 }, .slots = &.{ .saxpy_k, .caxpy_k, .caxpyc_k }, .written = .vector },
    .{ .module = CKernels(.power8), .cores = &.{ .power8, .power9 }, .slots = &.{ .daxpy_k, .zaxpy_k, .zaxpyc_k }, .written = .assembly },
    .{ .module = CKernels(.power10), .cores = &.{.power10}, .slots = &.{ .saxpy_k, .daxpy_k, .caxpy_k, .caxpyc_k, .zaxpy_k, .zaxpyc_k }, .written = .assembly },
};

const P = registry.Ports(Core, "power", &ports);

/// Exports every axpy kernel `core` has.
pub fn exportCore(comptime core: Core) void {
    inline for (kernel.slots) |slot| {
        if (comptime P.hasSlot(core, slot)) {
            P.exportOne(slot, core);
        } else if (comptime R.hasSlot(core, slot)) {
            R.exportOne(slot, core);
        }
    }
}

/// Whether `core` has any axpy kernel. (The 32-bit cores have no 64-bit
/// code; `kernels/power.zig`'s `runsOn` keeps them off 64-bit targets.)
pub fn has(comptime core: Core) bool {
    return P.has(core) or R.has(core);
}

/// Whether `slot` has a kernel on `core` for this target.
pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
    return P.hasSlot(core, slot) or R.hasSlot(core, slot);
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installSlots(vt, core, &kernel.slots, hasSlot);
}
