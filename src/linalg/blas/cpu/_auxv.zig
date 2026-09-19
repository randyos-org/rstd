//! `getauxval` access, the mechanism every non-x86 detector in
//! `driver/others/dynamic_*.c` leans on.
//!
//! x86 is the easy case -- `cpuid` is an unprivileged instruction, so
//! detection needs no OS cooperation at all. Every other architecture has
//! to ask the kernel: AArch64 needs `AT_HWCAP & HWCAP_CPUID` before it may
//! even read `MIDR_EL1`, PowerPC reads `AT_PLATFORM`, s390x and LoongArch
//! read capability bits directly out of `AT_HWCAP`, and so on. That means
//! none of them can detect anything on a freestanding target until we
//! provide an equivalent ourselves.

const std = @import("std");
const sysinfo = @import("builtin");

pub const AT_HWCAP: usize = 16;
pub const AT_HWCAP2: usize = 26;
pub const AT_PLATFORM: usize = 15;

/// Returns 0 when there is no auxiliary vector to read -- either because
/// the OS is not Linux, or because we are freestanding and nobody handed
/// us one. Callers treat 0 as "no capabilities known", which lands them on
/// the generic kernels.
pub fn get(tag: usize) usize {
    if (sysinfo.os.tag != .linux) return 0;
    return std.os.linux.getauxval(tag);
}

pub fn has(tag: usize, bit: usize) bool {
    return get(tag) & bit != 0;
}

pub fn hwcap() usize {
    return get(AT_HWCAP);
}

pub fn hwcap2() usize {
    return get(AT_HWCAP2);
}

test "auxv access is safe to call on any target" {
    // On non-Linux this is a compile-time-constant 0; on Linux without an
    // auxiliary vector it is a runtime 0. Either way it must not fault.
    _ = hwcap();
    _ = hwcap2();
}
