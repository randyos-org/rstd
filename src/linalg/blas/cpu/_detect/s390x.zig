//! IBM Z (s390x) CPU detection, mirroring `driver/others/dynamic_zarch.c`.
//!
//! Selection is purely capability-based, straight out of `AT_HWCAP`. As
//! the reference notes, vector registers may not be used on a z13 or newer
//! unless the OS kernel advertises support -- it has to save and restore
//! them across context switches -- so the hwcap bit, not the CPU model, is
//! the thing that decides.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const select = @import("../_select.zig");
const sysinfo = @import("builtin");
const auxv = @import("../_auxv.zig");

pub const CoreType = coretypes.S390xCoreType;
pub const Decision = select.Decision(CoreType);

/// Vector Facility (SIMD for double).
pub const HWCAP_S390_VX: usize = 2048;
/// Vector-Enhancements Facility 1 (float SIMD).
pub const HWCAP_S390_VXE: usize = 8192;

pub const Caps = struct {
    /// `AT_HWCAP`, under `LD_HWCAP_MASK`
    hwcap: usize = 0,
};

pub fn probe() Caps {
    return probeMasked(~@as(u64, 0));
}

/// `probe`, with the bits of `mask` (`LD_HWCAP_MASK`) only, as in the reference.
pub fn probeMasked(mask: u64) Caps {
    return .{ .hwcap = auxv.hwcap() & @as(usize, @truncate(mask)) };
}

/// z14's vector enhancements need the Vector Facility too.
pub fn runnable(core: CoreType, c: Caps) bool {
    return switch (core) {
        .zarch_generic => true,
        .z13 => c.hwcap & HWCAP_S390_VX != 0,
        .z14 => c.hwcap & HWCAP_S390_VX != 0 and c.hwcap & HWCAP_S390_VXE != 0,
    };
}

pub fn decide(c: Caps) Decision {
    if (c.hwcap & HWCAP_S390_VXE != 0) return select.guard(@This(), .{ .core = .z14 }, c);
    if (c.hwcap & HWCAP_S390_VX != 0) return .{ .core = .z13 };
    return .{ .core = .zarch_generic };
}

pub fn detect() CoreType {
    return decide(probe()).core;
}

pub fn detectNoted() Decision {
    return decide(probe());
}

/// `detect`, under `LD_HWCAP_MASK`.
pub fn detectMasked(mask: u64) CoreType {
    return decide(probeMasked(mask)).core;
}

test "the vector bits pick the core" {
    try std.testing.expectEqual(CoreType.zarch_generic, decide(.{}).core);
    try std.testing.expectEqual(CoreType.z13, decide(.{ .hwcap = HWCAP_S390_VX }).core);
    try std.testing.expectEqual(CoreType.z14, decide(.{ .hwcap = HWCAP_S390_VX | HWCAP_S390_VXE }).core);
    // VXE without VX is no machine; z14 kernels do not run on it
    try std.testing.expectEqual(CoreType.zarch_generic, decide(.{ .hwcap = HWCAP_S390_VXE }).core);
}

test "LD_HWCAP_MASK hides capabilities" {
    const all = HWCAP_S390_VX | HWCAP_S390_VXE | 1;
    const z14 = decide(.{ .hwcap = all & 0xffff }).core;
    try std.testing.expectEqual(CoreType.z14, z14);
    try std.testing.expectEqual(CoreType.z13, decide(.{ .hwcap = all & ~HWCAP_S390_VXE }).core);
    try std.testing.expectEqual(CoreType.zarch_generic, decide(.{ .hwcap = all & ~(HWCAP_S390_VX | HWCAP_S390_VXE) }).core);
}

test "detect runs on s390x targets" {
    if (comptime sysinfo.cpu.arch != .s390x) return error.SkipZigTest;
    _ = detect();
    _ = detectMasked(0);
}
