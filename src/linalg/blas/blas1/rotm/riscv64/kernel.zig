//! RISC-V `rotm` kernels: the `@Vector` ports (`vector.ported`) of OpenBLAS
//! `kernel/riscv64/rotm_rvv.c`, for the cores whose `KERNEL.<CORE>` names them. Per-core tests live
//! in `rotm/__root__.zig`.

const coretypes = @import("../../../cpu/coretypes.zig");
const registry = @import("../../../registry.zig");
const vector = @import("../../../vector.zig");

pub const Core = coretypes.RiscvCoreType;

/// The ports and where they serve.
pub const ports = [_]registry.Port(Core){
    .{ .module = vector.ported, .cores = &.{ .riscv64_zvl128b, .riscv64_zvl256b, .x280 }, .slots = &.{ .srotm_k, .drotm_k }, .written = .vector },
};

const P = registry.Ports(Core, "riscv64", &ports);

pub const hasSlot = P.hasSlot;
pub const has = P.has;
pub const exportCore = P.exportCore;
pub const install = P.install;
