//! MIPS axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`), for both word sizes, and how they
//! get installed into the vtable (`install`): cores without an entry keep
//! the portable kernels from `axpy.installGeneric`. Per-core tests live in
//! `axpy/__root__.zig`, not here -- they need the full axpy interface
//! (`checkReal`, ...), which imports outside `src/linalg/blas/` and so
//! cannot appear in this file (see `kernels/__root__.zig`).
//!
//! What a core uses is its `KERNEL.<CORE>` choice, made in the reference by
//! `NO_MSA` where a core may or may not have the SIMD architecture: the SIMD
//! cores use `mips/{s,d,c,z}axpy_msa.c` (`msa.zig`), the Loongson 3A R3 the
//! whole-function `.S` kernels in `asm/` (`asm/__root__.zig`). The cores with
//! neither keep the portable kernels. On MIPS64 the architecture's `KERNEL`
//! then resets the complex ones to `arm/zaxpy.c` (it is read after the core's
//! file), so only `saxpy` and `daxpy` use the SIMD kernels there.

const sysinfo = @import("builtin");
const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const util = @import("../../../util.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../kernel.zig");
const msa = @import("../mips/msa.zig");

const Slot = util.Slot;
const R = registry.Kernels(coretypes.MipsCoreType, @import("asm/__root__.zig"));

pub const Core = coretypes.MipsCoreType;

/// Whether `core` uses the SIMD kernel for `slot`.
fn usesMsa(comptime core: Core, comptime slot: Slot) bool {
    if (kernel.form(slot).complex and sysinfo.target.cpu.arch.isMips64()) return false;
    return switch (core) {
        .loongson3r4, .p5600, .i6400, .i6500, .p6600 => true,
        else => false,
    };
}

/// Whether `core` uses the SIMD kernel for any slot.
fn usesAnyMsa(comptime core: Core) bool {
    inline for (kernel.slots) |s| if (usesMsa(core, s)) return true;
    return false;
}

/// Whether `slot` has a kernel of its own on `core`.
pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
    return usesMsa(core, slot) or R.hasSlot(core, slot);
}

/// Whether `core` has any axpy kernel of its own.
pub fn has(comptime core: Core) bool {
    return usesAnyMsa(core) or R.has(core);
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (comptime usesAnyMsa(core)) {
        inline for (.{
            .{ Slot.saxpy_k, &msa.Real(f32).call },
            .{ Slot.daxpy_k, &msa.Real(f64).call },
            .{ Slot.caxpy_k, &msa.Cplx(f32, false).call },
            .{ Slot.caxpyc_k, &msa.Cplx(f32, true).call },
            .{ Slot.zaxpy_k, &msa.Cplx(f64, false).call },
            .{ Slot.zaxpyc_k, &msa.Cplx(f64, true).call },
        }) |e| {
            if (comptime usesMsa(core, e[0])) registry.exportKernel(e[1], e[0], core);
        }
    }
    R.exportCore(core);
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installSlots(vt, core, &kernel.slots, hasSlot);
}
