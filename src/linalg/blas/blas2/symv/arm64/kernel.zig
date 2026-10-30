//! arm64 `symv` kernels: the `@Vector` port (`../kernel.zig`) of
//! `kernel/arm64/symv_L_asimd_4x4.c` and `symv_U_asimd_4x4.c` (NEON
//! intrinsics). The SVE cores' `symv_[LU]_sve_v1x4.c` are not ported: no
//! tier-2 rig runs those cores, so they keep the portable kernels. Tests
//! live in `symv/__root__.zig`.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const util = @import("../../../util.zig");
const kernel = @import("../kernel.zig");

const Slot = util.Slot;

pub const Core = coretypes.Arm64CoreType;

const Asimd = struct {
    pub fn entry(comptime slot: Slot) util.SlotFn(slot) {
        return kernel.arm64(slot);
    }
};

/// The ports and where they serve (`KERNEL.<CORE>`'s `?SYMV_?_KERNEL`).
pub const ports = [_]registry.Port(Core){
    .{ .module = Asimd, .cores = &.{ .neoversen1, .ampere1, .vortex, .vortexm4 }, .slots = &.{ .ssymv_U, .ssymv_L, .dsymv_U, .dsymv_L }, .written = .vector },
};

const P = registry.Ports(Core, "arm64", &ports);

pub const hasSlot = P.hasSlot;
pub const has = P.has;
pub const exportCore = P.exportCore;
pub const install = P.install;
