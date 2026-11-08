//! WebAssembly axpy kernels: `saxpy` and `daxpy` on SIMD128 (`simd.zig`, `kernel/wasm/axpy.c`); tests in `../__root__.zig`.

const coretypes = @import("../../../cpu/coretypes.zig");
const registry = @import("../../../registry.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../driver.zig");
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
