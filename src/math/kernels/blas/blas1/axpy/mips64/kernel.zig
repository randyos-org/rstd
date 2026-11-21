//! MIPS64 axpy kernels per core (`msa.zig` for `saxpy` and `daxpy` where the core has MSA, `asm/` for Loongson 3A R3); tests in `../__root__.zig`.

const sysinfo = @import("builtin");
const coretypes = @import("../../../cpu/coretypes.zig");
const registry = @import("../../../registry.zig");
const util = @import("../../../util.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../driver.zig");
const msa = @import("../mips/msa.zig");

const Slot = util.Slot;
const R = registry.Kernels(coretypes.MipsCoreType, @import("asm/__root__.zig"));

pub const Core = coretypes.MipsCoreType;

/// The hand ports and where they serve: the SIMD kernels, but for the
/// complex ones on MIPS64.
pub const ports = [_]registry.Port(Core){
    .{ .module = kernel.Ported(msa), .cores = &.{ .loongson3r4, .p5600, .i6400, .i6500, .p6600 }, .slots = &.{ .saxpy_k, .daxpy_k }, .written = .vector },
    .{ .module = kernel.Ported(msa), .cores = &.{ .loongson3r4, .p5600, .i6400, .i6500, .p6600 }, .slots = &.{ .caxpy_k, .caxpyc_k, .zaxpy_k, .zaxpyc_k }, .written = .vector, .arch = "mips" },
};

const P = registry.Ports(Core, if (sysinfo.target.cpu.arch.isMips64()) "mips64" else "mips", &ports);

/// Whether `slot` has a kernel of its own on `core`.
pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
    return P.hasSlot(core, slot) or R.hasSlot(core, slot);
}

/// Whether `core` has any axpy kernel of its own.
pub fn has(comptime core: Core) bool {
    return P.has(core) or R.has(core);
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    inline for (kernel.slots) |slot| P.exportOne(slot, core);
    R.exportCore(core);
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installSlots(vt, core, &kernel.slots, hasSlot);
}
