//! LoongArch CPU detection, mirroring
//! `driver/others/dynamic_loongarch64.c`.
//!
//! Two sources: `AT_HWCAP` for the vector extensions (LSX 128-bit, LASX
//! 256-bit), and the `cpucfg` instruction for the processor identity.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const sysinfo = @import("builtin");
const auxv = @import("../_auxv.zig");

pub const CoreType = coretypes.LoongarchCoreType;

pub const HWCAP_LSX: usize = 1 << 4;
pub const HWCAP_LASX: usize = 1 << 5;

/// `cpucfg` is unprivileged, like x86's `cpuid`.
pub fn cpucfg(word: u32) u32 {
    return asm volatile ("cpucfg %[out], %[in]"
        : [out] "=r" (-> u32),
        : [in] "r" (word),
    );
}

pub fn supportLsx() bool {
    return auxv.hwcap() & HWCAP_LSX != 0;
}

pub fn supportLasx() bool {
    return auxv.hwcap() & HWCAP_LASX != 0;
}

pub fn detect() CoreType {
    if (supportLasx()) return .la464;
    if (supportLsx()) return .la264;
    return .la64_generic;
}

test "detect runs on LoongArch targets" {
    if (comptime !sysinfo.cpu.arch.isLoongArch()) return error.SkipZigTest;
    _ = detect();
}
