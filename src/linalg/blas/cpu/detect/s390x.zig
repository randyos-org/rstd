//! IBM Z (s390x) CPU detection, mirroring `driver/others/dynamic_zarch.c`:
//! by `AT_HWCAP` alone, as vector registers are usable only where the OS
//! saves them across context switches.

const std = @import("std");
const coretypes = @import("../coretypes.zig");
const select = @import("../select.zig");
const sysinfo = @import("builtin");
const osinfo = @import("../osinfo.zig");

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
    const hw = if (comptime sysinfo.target.os.tag == .freestanding) facilities() else osinfo.hwcap();
    return .{ .hwcap = hw & @as(usize, @truncate(mask)) };
}

/// With no OS to ask, the facility list (`STFLE`, unprivileged) as hwcap bits:
/// facility 129 is the Vector Facility, 135 Vector-Enhancements 1.
fn facilities() usize {
    if (comptime sysinfo.target.cpu.arch != .s390x) return 0;
    var list: [3]u64 = @splat(0);
    _ = asm volatile ("stfle 0(%[p])"
        : [n] "={r0}" (-> u64),
        : [p] "a" (&list),
          [in] "{r0}" (@as(u64, list.len - 1)),
        : .{ .memory = true, .ps = true });
    const has = struct {
        fn bit(l: [3]u64, n: u8) bool {
            return (l[n / 64] >> @intCast(63 - n % 64)) & 1 != 0;
        }
    }.bit;
    return (if (has(list, 129)) HWCAP_S390_VX else 0) | (if (has(list, 135)) HWCAP_S390_VXE else 0);
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
    if (comptime sysinfo.target.cpu.arch != .s390x) return error.SkipZigTest;
    try std.testing.expect(select.canRun(@This(), select.detect(@This(), null)));
    _ = select.detect(@This(), 0);
}
