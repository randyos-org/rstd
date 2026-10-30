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
const kernel = @import("../kernel.zig");
const simd = @import("simd.zig");

pub const Core = coretypes.WasmCoreType;

/// The hand ports and where they serve.
pub const ports = [_]registry.Port(Core){
    .{ .module = kernel.Ported(simd), .cores = &.{.wasm128_generic}, .slots = &.{ .saxpy_k, .daxpy_k }, .written = .vector },
};

const P = registry.Ports(Core, "wasm", &ports);

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
