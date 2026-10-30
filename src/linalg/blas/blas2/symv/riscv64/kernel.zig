//! RISC-V `symv` kernels: the `@Vector` ports (`vector.ported`) of OpenBLAS
//! `kernel/riscv64/symv_U_rvv.c`, `symv_L_rvv.c`, `zsymv_U_rvv.c`, `zsymv_L_rvv.c`,
//! `symv_U_vector.c` and `symv_L_vector.c`, for the cores whose `KERNEL.<CORE>` names them.
//! zvl256b's complex `symv` is `../generic/zsymv_k.c` upstream, which the portable kernel composes
//! from the core's own gemv. Per-core tests live in `symv/__root__.zig`.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const vector = @import("../../../vector.zig");

pub const Core = coretypes.RiscvCoreType;

/// The ports and where they serve.
pub const ports = [_]registry.Port(Core){
    .{ .module = vector.ported, .cores = &.{ .riscv64_zvl128b, .riscv64_zvl256b, .x280 }, .slots = &.{ .ssymv_U, .ssymv_L, .dsymv_U, .dsymv_L }, .written = .vector },
    .{ .module = vector.ported, .cores = &.{ .riscv64_zvl128b, .x280 }, .slots = &.{ .csymv_U, .csymv_L, .zsymv_U, .zsymv_L }, .written = .vector },
};

const P = registry.Ports(Core, "riscv64", &ports);

pub const hasSlot = P.hasSlot;
pub const has = P.has;
pub const exportCore = P.exportCore;
pub const install = P.install;
