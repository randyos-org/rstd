//! 32-bit ARM CPU detection.
//!
//! Note the reference has no `dynamic_arm.c` -- 32-bit ARM is the one
//! supported architecture with no runtime dispatch at all, so the kernel
//! set is fixed at build time from `TARGET` (ARMV5/ARMV6/ARMV7/NEON). This
//! exists so every architecture presents the same shape to the dispatcher,
//! and reports what `AT_HWCAP` says rather than inventing a model table.

const std = @import("std");
const coretypes = @import("_coretypes.zig");
const sysinfo = @import("builtin");
const auxv = @import("_auxv.zig");

pub const CoreType = coretypes.ArmCoreType;
pub const kernel = @import("../kernels/arm.zig");

pub const HWCAP_VFP: usize = 1 << 6;
pub const HWCAP_NEON: usize = 1 << 12;
pub const HWCAP_VFPv3: usize = 1 << 13;
pub const HWCAP_VFPv4: usize = 1 << 16;

pub fn detect() CoreType {
    const caps = auxv.hwcap();
    if (caps & HWCAP_VFPv4 != 0) return .vfpv4;
    if (caps & HWCAP_NEON != 0) return .neon;
    if (caps & HWCAP_VFPv3 != 0) return .armv7;
    if (caps & HWCAP_VFP != 0) return .armv6;
    return .armv5;
}
const Vtable = @import("../Vtable.zig");
const blas1 = @import("../blas1/__root__.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");
pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    blas1.installGeneric(&vt);
    kernel.install(&vt, core_type);
    gemm_params.installCore(&vt, .arm, core_type, .fromConfig(blas.config));
    return vt;
}

test "detect runs on 32-bit ARM targets" {
    if (comptime !sysinfo.cpu.arch.isArm()) return error.SkipZigTest;
    _ = detect();
}
