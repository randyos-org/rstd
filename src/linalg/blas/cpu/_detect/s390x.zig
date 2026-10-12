//! IBM Z (s390x) CPU detection, mirroring `driver/others/dynamic_zarch.c`.
//!
//! Selection is purely capability-based, straight out of `AT_HWCAP`. As
//! the reference notes, vector registers may not be used on a z13 or newer
//! unless the OS kernel advertises support -- it has to save and restore
//! them across context switches -- so the hwcap bit, not the CPU model, is
//! the thing that decides.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const sysinfo = @import("builtin");
const auxv = @import("../_auxv.zig");

pub const CoreType = coretypes.S390xCoreType;

/// Vector Facility (SIMD for double).
pub const HWCAP_S390_VX: usize = 2048;
/// Vector-Enhancements Facility 1 (float SIMD).
pub const HWCAP_S390_VXE: usize = 8192;

pub fn detect() CoreType {
    return detectMasked(~@as(u64, 0));
}

/// `detect`, with the bits of `mask` (`LD_HWCAP_MASK`) only, as in the reference.
pub fn detectMasked(mask: u64) CoreType {
    const caps = auxv.hwcap() & mask;
    if (caps & HWCAP_S390_VXE != 0) return .z14;
    if (caps & HWCAP_S390_VX != 0) return .z13;
    return .zarch_generic;
}

test "detect runs on s390x targets" {
    if (comptime sysinfo.cpu.arch != .s390x) return error.SkipZigTest;
    _ = detect();
}
