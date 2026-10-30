//! RISC-V `dsdot` kernels: the `@Vector` ports (`vector.ported`) of OpenBLAS
//! `kernel/riscv64/dot_rvv.c` with `DSDOT`, for the cores whose `KERNEL.<CORE>` names them. The
//! other `dot` slots take `vector.zig`'s kernels through the registry's `vector_rule`; zvl256b's
//! `dsdot` is plain C upstream (`../generic/dot.c`). Per-core tests live in `dot/__root__.zig`.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const vector = @import("../../../vector.zig");

pub const Core = coretypes.RiscvCoreType;

/// The ports and where they serve.
pub const ports = [_]registry.Port(Core){
    .{ .module = vector.ported, .cores = &.{ .riscv64_zvl128b, .x280 }, .slots = &.{.dsdot_k}, .written = .vector },
};

const P = registry.Ports(Core, "riscv64", &ports);

pub const hasSlot = P.hasSlot;
pub const has = P.has;
pub const exportCore = P.exportCore;
pub const install = P.install;
