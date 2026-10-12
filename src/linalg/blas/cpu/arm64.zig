//! AArch64 CPU detection, mirroring `driver/others/dynamic_arm64.c`.
//!
//! Unlike x86, AArch64 cannot identify itself unaided: `MIDR_EL1` is a
//! system register that faults in user mode unless the kernel emulates the
//! read, which Linux advertises through `AT_HWCAP & HWCAP_CPUID`. So the
//! reference has three separate paths -- Darwin (Apple silicon, answered
//! by `sysctl`), Linux/Android (`MIDR_EL1`, with `/proc/cpuinfo` and
//! `/sys/.../midr_el1` as fallbacks), and everything else (generic).

const std = @import("std");
const coretypes = @import("_coretypes.zig");
const zigstd = @import("../../../zigstd.zig");
const sysinfo = @import("builtin");
const auxv = @import("_auxv.zig");

pub const CoreType = coretypes.Arm64CoreType;
pub const kernel = @import("../kernels/arm64.zig");

pub const HWCAP_CPUID: usize = 1 << 11;
pub const HWCAP_SVE: usize = 1 << 22;
pub const HWCAP2_SME: usize = 1 << 23;

pub const Midr = struct {
    implementer: u32,
    variant: u32,
    architecture: u32,
    part: u32,
    revision: u32,
};

/// Reads `MIDR_EL1`. Only legal once `HWCAP_CPUID` says the kernel will
/// emulate it -- calling this unguarded on a kernel without that support
/// raises SIGILL.
pub fn readMidr() Midr {
    const raw = asm volatile ("mrs %[out], midr_el1"
        : [out] "=r" (-> u64),
    );
    return .{
        .implementer = @truncate((raw >> 24) & 0xff),
        .variant = @truncate((raw >> 20) & 0xf),
        .architecture = @truncate((raw >> 16) & 0xf),
        .part = @truncate((raw >> 4) & 0xfff),
        .revision = @truncate(raw & 0xf),
    };
}

pub fn supportSve() bool {
    return auxv.hwcap() & HWCAP_SVE != 0;
}

pub fn supportSme() bool {
    if (sysinfo.os.tag.isDarwin()) return darwinSysctlFlag("hw.optional.arm.FEAT_SME");
    return auxv.hwcap2() & HWCAP2_SME != 0;
}

extern "c" fn sysctlbyname(
    name: [*:0]const u8,
    oldp: ?*anyopaque,
    oldlenp: ?*usize,
    newp: ?*const anyopaque,
    newlen: usize,
) c_int;

fn darwinSysctlFlag(comptime name: [:0]const u8) bool {
    if (!sysinfo.link_libc) return false;
    var value: u64 = 0;
    var len: usize = @sizeOf(u64);
    if (sysctlbyname(name.ptr, &value, &len, null, 0) != 0) return false;
    return value != 0;
}

pub fn detect() CoreType {
    // Apple silicon: no MIDR access, and `sysctl` answers the only
    // question the reference asks here.
    if (sysinfo.os.tag.isDarwin()) {
        return if (supportSme()) .vortexm4 else .vortex;
    }

    // Anywhere that is not Linux/Android has no way to ask, so the
    // reference returns NULL and lands on the generic kernels.
    if (sysinfo.os.tag != .linux) return .armv8;

    // Without CPUID emulation the reference falls back to parsing
    // /proc/cpuinfo or /sys/devices/system/cpu/cpu0/regs/identification/
    // midr_el1. Reading MIDR_EL1 anyway would fault, so stay generic.
    if (auxv.hwcap() & HWCAP_CPUID == 0) {
        return if (supportSve()) .armv8sve else .armv8;
    }

    const midr = readMidr();
    const core: CoreType = switch (midr.implementer) {
        0x41 => switch (midr.part) { // ARM
            0xd03 => .cortexa53,
            0xd04 => .cortexa53, // Cortex-A35
            0xd05 => .cortexa55,
            0xd07 => .cortexa57,
            0xd08 => .cortexa72,
            0xd09 => .cortexa73,
            0xd0c => .neoversen1,
            0xd40 => .neoversev1,
            0xd49 => .neoversen2,
            0xd4f => .neoversev2,
            else => if (supportSve()) .armv8sve else .armv8,
        },
        0x42, 0x43 => switch (midr.part) { // Broadcom / Cavium
            0x0a1 => .thunderx,
            0x0af => .thunderx2t99,
            0x0b8 => .thunderx3t110,
            else => .armv8,
        },
        0x46 => .a64fx, // Fujitsu
        0x48 => .tsv110, // HiSilicon
        0x50 => .emag8180, // Ampere (APM)
        0x51 => .falkor, // Qualcomm
        0x61 => if (supportSme()) .vortexm4 else .vortex, // Apple
        0xc0 => .neoversen1, // Ampere Computing
        else => if (supportSve()) .armv8sve else .armv8,
    };

    // The reference gives SME-capable parts the SME kernels regardless of
    // which microarchitecture matched above.
    if (core == .armv8 and supportSme()) return .armv9sme;
    return core;
}

const Vtable = @import("../Vtable.zig");
const blas1 = @import("../blas1/__root__.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");
pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    blas1.installGeneric(&vt);
    kernel.install(&vt, core_type);
    gemm_params.installCore(&vt, .arm64, core_type, .fromConfig(blas.config));
    return vt;
}

test "detect runs and returns a concrete core type" {
    if (comptime sysinfo.cpu.arch != .aarch64 and sysinfo.cpu.arch != .aarch64_be) {
        return error.SkipZigTest;
    }
    try zigstd.testExpectTrue(detect() != .unknown);
}
