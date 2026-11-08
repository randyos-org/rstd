//! LoongArch CPU detection, mirroring `driver/others/dynamic_loongarch64.c`:
//! `AT_HWCAP` for LSX and LASX, `cpucfg` word 0 (the PRID) for the series.

const std = @import("std");
const coretypes = @import("../coretypes.zig");
const select = @import("../select.zig");
const sysinfo = @import("builtin");
const osinfo = @import("../osinfo.zig");

pub const CoreType = coretypes.LoongarchCoreType;
pub const Decision = select.Decision(CoreType);

pub const HWCAP_LSX: usize = 1 << 4;
pub const HWCAP_LASX: usize = 1 << 5;
const CPUCFG2_LSX: u32 = 1 << 6;
const CPUCFG2_LASX: u32 = 1 << 7;

const PRID_SERIES_MASK: u32 = 0xf000;
const PRID_SERIES_LA264: u32 = 0xa000;
const PRID_SERIES_LA364: u32 = 0xb000;
const PRID_SERIES_LA464: u32 = 0xc000;
const PRID_SERIES_LA664: u32 = 0xd000;

pub const Feature = enum {
    /// 128-bit vectors
    lsx,
    /// 256-bit vectors
    lasx,
};

pub const Caps = struct {
    prid: u32 = 0,
    features: std.EnumSet(Feature) = .empty,

    pub fn has(c: Caps, f: Feature) bool {
        return c.features.contains(f);
    }

    /// Every one of `fs`, the set built at compile time.
    pub fn hasAll(c: Caps, comptime fs: []const Feature) bool {
        return c.features.supersetOf(comptime std.EnumSet(Feature).initMany(fs));
    }
};

/// `cpucfg` is unprivileged, like x86's `cpuid`.
pub fn cpucfg(word: u32) u32 {
    return asm volatile ("cpucfg %[out], %[in]"
        : [out] "=r" (-> u32),
        : [in] "r" (word),
    );
}

pub fn supportLsx() bool {
    return osinfo.hwcap() & HWCAP_LSX != 0;
}

pub fn supportLasx() bool {
    return osinfo.hwcap() & HWCAP_LASX != 0;
}

pub fn probe() Caps {
    if (comptime !sysinfo.target.cpu.arch.isLoongarch()) return .{};
    switch (comptime sysinfo.target.os.tag) {
        // no OS to ask: `cpucfg` word 2 says what the CPU has
        .freestanding, .uefi => {
            const cfg2 = cpucfg(2);
            return .{ .prid = cpucfg(0), .features = .init(.{ .lsx = cfg2 & CPUCFG2_LSX != 0, .lasx = cfg2 & CPUCFG2_LASX != 0 }) };
        },
        else => return .{ .prid = cpucfg(0), .features = .init(.{ .lsx = supportLsx(), .lasx = supportLasx() }) },
    }
}

pub fn runnable(core: CoreType, c: Caps) bool {
    return switch (core) {
        .la264, .loongson2k1000 => c.has(.lsx),
        .la464, .loongson3r5 => c.hasAll(&.{ .lsx, .lasx }),
        .la64_generic, .loongsongeneric => true,
        .unknown => false,
    };
}

/// The series the PRID names picks the ceiling, the OS-enabled extensions the
/// core within it; an unknown series is generic.
pub fn decide(c: Caps) Decision {
    const d: Decision = switch (c.prid & PRID_SERIES_MASK) {
        PRID_SERIES_LA464, PRID_SERIES_LA664 => if (c.has(.lasx))
            .{ .core = .la464 }
        else if (c.has(.lsx))
            .{ .core = .la264 }
        else
            .{ .core = .la64_generic },
        PRID_SERIES_LA264, PRID_SERIES_LA364 => if (c.has(.lsx))
            .{ .core = .la264 }
        else
            .{ .core = .la64_generic },
        else => .{ .core = .la64_generic },
    };
    return select.guard(@This(), d, c);
}

test "series and OS extensions pick the core" {
    const t = struct {
        fn core(prid: u32, lsx: bool, lasx: bool) CoreType {
            return decide(.{ .prid = prid, .features = .init(.{ .lsx = lsx, .lasx = lasx }) }).core;
        }
    };
    try std.testing.expectEqual(CoreType.la464, t.core(0x14c010, true, true)); // 3A5000
    try std.testing.expectEqual(CoreType.la464, t.core(0x14d010, true, true)); // 3A6000
    try std.testing.expectEqual(CoreType.la264, t.core(0x14c010, true, false));
    try std.testing.expectEqual(CoreType.la64_generic, t.core(0x14c010, false, false));
    try std.testing.expectEqual(CoreType.la264, t.core(0x14a010, true, false)); // 2K1000LA
    try std.testing.expectEqual(CoreType.la264, t.core(0x14b010, true, true));
    try std.testing.expectEqual(CoreType.la64_generic, t.core(0x14a010, false, false));
    try std.testing.expectEqual(CoreType.la64_generic, t.core(0x140010, true, true));
}

test "detect runs on LoongArch targets" {
    if (comptime !sysinfo.target.cpu.arch.isLoongarch()) return error.SkipZigTest;
    try std.testing.expect(select.canRun(@This(), select.detect(@This(), null)));
}
