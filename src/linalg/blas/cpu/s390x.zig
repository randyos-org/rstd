//! IBM Z (s390x) CPU detection, mirroring `driver/others/dynamic_zarch.c`.
//!
//! Selection is purely capability-based, straight out of `AT_HWCAP`. As
//! the reference notes, vector registers may not be used on a z13 or newer
//! unless the OS kernel advertises support -- it has to save and restore
//! them across context switches -- so the hwcap bit, not the CPU model, is
//! the thing that decides.

const std = @import("std");
const coretypes = @import("_coretypes.zig");
const sysinfo = @import("builtin");
const auxv = @import("_auxv.zig");

pub const CoreType = coretypes.S390xCoreType;
pub const kernel = @import("../kernels/s390x.zig");

/// Vector Facility (SIMD for double).
pub const HWCAP_S390_VX: usize = 2048;
/// Vector-Enhancements Facility 1 (float SIMD).
pub const HWCAP_S390_VXE: usize = 8192;

pub fn detect() CoreType {
    // `LD_HWCAP_MASK` masks bits out of what the kernel reports, as in the reference
    const caps = auxv.hwcap() & (blas.config.ld_hwcap_mask orelse ~@as(u64, 0));
    if (caps & HWCAP_S390_VXE != 0) return .z14;
    if (caps & HWCAP_S390_VX != 0) return .z13;
    return .zarch_generic;
}

const Vtable = @import("../Vtable.zig");
const blas1 = @import("../blas1/__root__.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");
pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    blas1.installGeneric(&vt);
    kernel.install(&vt, core_type);
    gemm_params.installCore(&vt, .s390x, core_type, .fromConfig(blas.config));
    return vt;
}

test "detect runs on s390x targets" {
    if (comptime sysinfo.cpu.arch != .s390x) return error.SkipZigTest;
    _ = detect();
}
