//! LoongArch CPU detection, mirroring `driver/others/dynamic_loongarch64.c`.
//!
//! Two sources: `AT_HWCAP` for the vector extensions (LSX 128-bit, LASX
//! 256-bit), and `cpucfg` word 0 (the PRID) for the processor series.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const select = @import("../_select.zig");
const sysinfo = @import("builtin");
const auxv = @import("../_auxv.zig");

pub const CoreType = coretypes.LoongarchCoreType;
pub const Decision = select.Decision(CoreType);

pub const HWCAP_LSX: usize = 1 << 4;
pub const HWCAP_LASX: usize = 1 << 5;

const PRID_SERIES_MASK: u32 = 0xf000;
const PRID_SERIES_LA264: u32 = 0xa000;
const PRID_SERIES_LA364: u32 = 0xb000;
const PRID_SERIES_LA464: u32 = 0xc000;
const PRID_SERIES_LA664: u32 = 0xd000;

pub const Caps = struct {
    prid: u32 = 0,
    lsx: bool = false,
    lasx: bool = false,
};

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

pub fn probe() Caps {
    if (comptime !sysinfo.cpu.arch.isLoongArch()) return .{};
    return .{ .prid = cpucfg(0), .lsx = supportLsx(), .lasx = supportLasx() };
}

pub fn runnable(core: CoreType, c: Caps) bool {
    return switch (core) {
        .la264, .loongson2k1000 => c.lsx,
        .la464, .loongson3r5 => c.lsx and c.lasx,
        .la64_generic, .loongsongeneric => true,
        .unknown => false,
    };
}

/// The series the PRID names picks the ceiling, the OS-enabled extensions the
/// core within it; an unknown series is generic.
pub fn decide(c: Caps) Decision {
    const d: Decision = switch (c.prid & PRID_SERIES_MASK) {
        PRID_SERIES_LA464, PRID_SERIES_LA664 => if (c.lasx)
            .{ .core = .la464 }
        else if (c.lsx)
            .{ .core = .la264 }
        else
            .{ .core = .la64_generic },
        PRID_SERIES_LA264, PRID_SERIES_LA364 => if (c.lsx)
            .{ .core = .la264 }
        else
            .{ .core = .la64_generic },
        else => .{ .core = .la64_generic },
    };
    return select.guard(@This(), d, c);
}

pub fn detect() CoreType {
    return decide(probe()).core;
}

pub fn detectNoted() Decision {
    return decide(probe());
}

test "series and OS extensions pick the core" {
    const t = struct {
        fn core(prid: u32, lsx: bool, lasx: bool) CoreType {
            return decide(.{ .prid = prid, .lsx = lsx, .lasx = lasx }).core;
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
    if (comptime !sysinfo.cpu.arch.isLoongArch()) return error.SkipZigTest;
    _ = detect();
}
