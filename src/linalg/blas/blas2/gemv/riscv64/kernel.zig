//! RISC-V half-precision `gemv` kernels: the `@Vector` ports (`vector.ported`) of OpenBLAS
//! `kernel/riscv64/sbgemv_n_vector.c` and `sbgemv_t_vector.c` (`sbgemv`, and `shgemv` under
//! `HFLOAT16`), for the cores whose `KERNEL.<CORE>` names them. They keep the portable kernels'
//! bits (`vector.Half`). The other `gemv` slots take `vector.zig`'s kernels through the registry's
//! `vector_rule`. Per-core tests live in `gemv/__root__.zig`.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const vector = @import("../../../vector.zig");

pub const Core = coretypes.RiscvCoreType;

/// The ports and where they serve.
pub const ports = [_]registry.Port(Core){
    .{ .module = vector.ported, .cores = &.{.riscv64_zvl256b}, .slots = &.{ .sbgemv_n, .sbgemv_t, .shgemv_n, .shgemv_t }, .written = .vector },
};

const P = registry.Ports(Core, "riscv64", &ports);

pub const hasSlot = P.hasSlot;
pub const has = P.has;
pub const exportCore = P.exportCore;
pub const install = P.install;
