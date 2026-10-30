//! arm64 gemv kernels: the converted `.S` kernels (`asm/`), and the
//! half-precision ports on the BF16 cores: `sbgemv_n_neon.c`,
//! `sbgemv_t_bfdot.c` and `bgemv_n_sve_v3x4.c` as `vector.Half`, which keeps
//! the portable kernels' bits where upstream's BFDOT and BFMLAL round as
//! they do. The other gemv slots take `vector.zig`'s kernels through the
//! registry's `vector_rule`. Tests live in `gemv/__root__.zig`.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const util = @import("../../../util.zig");
const vector = @import("../../../vector.zig");
const Vtable = @import("../../../Vtable.zig");

const Slot = util.Slot;

pub const Core = coretypes.Arm64CoreType;

/// The ports and where they serve (`KERNEL.<CORE>`'s `SBGEMV?KERNEL`, `BGEMV?KERNEL`).
pub const ports = [_]registry.Port(Core){
    .{ .module = vector.ported, .cores = &.{ .neoversen2, .neoversev1, .neoversev2 }, .slots = &.{ .sbgemv_n, .sbgemv_t, .bgemv_n, .bgemv_t }, .written = .vector },
};

const P = registry.Ports(Core, "arm64", &ports);
const A = registry.Kernels(Core, @import("asm/__root__.zig"));

pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
    return A.hasSlot(core, slot) or P.hasSlot(core, slot);
}

pub fn has(comptime core: Core) bool {
    return A.has(core) or P.has(core);
}

pub fn exportCore(comptime core: Core) void {
    A.exportCore(core);
    P.exportCore(core);
}

pub fn install(vt: *Vtable, comptime core: Core) void {
    A.install(vt, core);
    P.install(vt, core);
}
