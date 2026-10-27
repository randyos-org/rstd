//! WebAssembly axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`), and how they get installed into
//! the vtable (`install`): cores without an entry keep the portable kernels
//! from `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`,
//! not here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! Only `saxpy` and `daxpy` have a kernel of their own (`simd.zig`, from
//! `kernel/wasm/axpy.c`, for the core with SIMD128); the complex kernels
//! upstream are the plain C loops, which the portable kernels are.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const Vtable = @import("../../../Vtable.zig");
const util = @import("../../../util.zig");
const kernel = @import("../kernel.zig");
const simd = @import("simd.zig");

const Slot = util.Slot;

pub const Core = coretypes.WasmCoreType;

/// Whether `slot` has a kernel of its own on `core`.
pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
    return core == .wasm128_generic and !kernel.form(slot).complex;
}

/// Whether `core` has any axpy kernel of its own.
pub fn has(comptime core: Core) bool {
    return core == .wasm128_generic;
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (comptime core != .wasm128_generic) return;
    registry.exportKernel(&simd.Real(f32).call, .saxpy_k, core);
    registry.exportKernel(&simd.Real(f64).call, .daxpy_k, core);
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installSlots(vt, core, &kernel.slots, hasSlot);
}
