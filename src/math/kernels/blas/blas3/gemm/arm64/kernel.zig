//! arm64 gemm kernels of their own: `vortexm4`'s SME2 direct real `gemm`
//! (`KERNEL.VORTEXM4`'s clang build, `sme_sgemm_kernel.c` and
//! `sme_dgemm_kernel.c`; `../sme2/kernel.zig`).

const coretypes = @import("../../../cpu/coretypes.zig");
const registry = @import("../../../registry.zig");
const sme2 = @import("../sme2/kernel.zig");

pub const Core = coretypes.Arm64CoreType;

/// The ports and where they serve.
pub const ports = [_]registry.Port(Core){
    .{ .module = sme2, .cores = &.{.vortexm4}, .slots = &.{ .sme_sgemm_kernel, .sme_dgemm_kernel }, .written = .assembly },
};

const P = registry.Ports(Core, "arm64", &ports);
pub const hasSlot = P.hasSlot;
pub const has = P.has;
pub const exportCore = P.exportCore;
pub const install = P.install;
