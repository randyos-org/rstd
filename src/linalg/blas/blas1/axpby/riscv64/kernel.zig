//! RISC-V `axpby` kernels: the `@Vector` ports (`vector.ported`) of OpenBLAS
//! `kernel/riscv64/axpby_rvv.c`, `zaxpby_rvv.c` and `axpby_vector_v2.c`, for the cores whose
//! `KERNEL.<CORE>` names them. zvl256b's complex `axpby` is plain C upstream (`../arm/zaxpby.c`)
//! and stays portable. Per-core tests live in `axpby/__root__.zig`.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const vector = @import("../../../vector.zig");

pub const Core = coretypes.RiscvCoreType;

/// The ports and where they serve.
pub const ports = [_]registry.Port(Core){
    .{ .module = vector.ported, .cores = &.{ .riscv64_zvl128b, .x280 }, .slots = &.{ .saxpby_k, .daxpby_k, .caxpby_k, .zaxpby_k }, .written = .vector },
    .{ .module = vector.ported, .cores = &.{.riscv64_zvl256b}, .slots = &.{ .saxpby_k, .daxpby_k }, .written = .vector },
};

const P = registry.Ports(Core, "riscv64", &ports);

pub const hasSlot = P.hasSlot;
pub const has = P.has;
pub const exportCore = P.exportCore;
pub const install = P.install;
