//! MIPS / MIPS64 CPU detection, mirroring `driver/others/dynamic_mips64.c`.
//!
//! Loongson's `cpucfg` instruction when `AT_HWCAP` advertises it, otherwise the
//! model name in `/proc/cpuinfo`. The reference tests the hwcap bit the wrong way
//! round (it runs `cpucfg` where the kernel says it is absent, and reads
//! `/proc/cpuinfo` where it is present), and its `cpucfg` word has both
//! registers zeroed; both are done as the manuals say here.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const select = @import("../_select.zig");
const sysinfo = @import("builtin");
const auxv = @import("../_auxv.zig");

pub const CoreType = coretypes.MipsCoreType;
pub const Decision = select.Decision(CoreType);

pub const HWCAP_LOONGSON_CPUCFG: usize = 1 << 14;

const LOONGSON_CFG1: u32 = 0x1;
const MMI_MASK: u32 = 0x00000010;
const MSA_MASK: u32 = 0x00000020;

pub const Feature = enum {
    /// Loongson multimedia instructions
    mmi,
    /// MIPS SIMD architecture
    msa,
};

pub const Caps = struct {
    features: std.EnumSet(Feature) = .empty,

    pub fn has(c: Caps, f: Feature) bool {
        return c.features.contains(f);
    }
};

/// Whether this machine can run `core`'s kernels (for tests): the Loongson
/// 3A R4's need the SIMD architecture, which no release 2 machine here has.
pub fn canRun(core: CoreType) bool {
    return core != .loongson3r4 or sysinfo.cpu.has(.mips, .msa);
}

/// `CPUCFG`'s encoding as the Linux kernel writes it (`read_cpucfg`): the
/// register holding the word number in bits 25..21, the result in 15..11.
pub fn cpucfgWord(reg: u5, res: u5) u32 {
    return 0xc8080118 | (@as(u32, reg) << 21) | (@as(u32, res) << 11);
}

/// Loongson's configuration-word read, with `$4` in and `$2` out
/// (`cpucfgWord(4, 2)`).
pub fn cpucfg(word: u32) u32 {
    return asm volatile (
        \\.set push
        \\.set noreorder
        \\.insn
        \\.word 0xc8881118
        \\.set pop
        : [out] "={$2}" (-> u32),
        : [in] "{$4}" (word),
    );
}

/// Caps from CPUCFG word 1.
pub fn capsFromCfg1(cfg: u32) Caps {
    return .{ .features = .init(.{ .mmi = cfg & MMI_MASK != 0, .msa = cfg & MSA_MASK != 0 }) };
}

/// Caps from the first "model name" line of `/proc/cpuinfo`: 3A/3B 3000 are the
/// release 3 parts (MMI), 3A/3B 4000 the release 4 ones (MMI and MSA).
pub fn capsFromCpuinfo(text: []const u8) Caps {
    const name = auxv.cpuinfoField(text, "model name") orelse return .{};
    const has = struct {
        fn any(n: []const u8, a: []const u8, b: []const u8) bool {
            return std.mem.indexOf(u8, n, a) != null or std.mem.indexOf(u8, n, b) != null;
        }
    }.any;
    if (has(name, "Loongson-3A3000", "Loongson-3B3000")) return .{ .features = .initOne(.mmi) };
    if (has(name, "Loongson-3A4000", "Loongson-3B4000")) return .{ .features = .initMany(&.{ .mmi, .msa }) };
    return .{};
}

pub fn probe() Caps {
    if (comptime !sysinfo.cpu.arch.isMIPS()) return .{};
    if (sysinfo.os.tag != .linux) return .{};
    if (auxv.hwcap() & HWCAP_LOONGSON_CPUCFG != 0) return capsFromCfg1(cpucfg(LOONGSON_CFG1));
    var buf: [4096]u8 = undefined;
    return capsFromCpuinfo(auxv.readFile("/proc/cpuinfo", &buf));
}

pub fn runnable(core: CoreType, c: Caps) bool {
    return switch (core) {
        .loongson3r4 => c.has(.msa),
        .loongson3r3 => c.has(.mmi) or c.has(.msa),
        .unknown => false,
        else => true,
    };
}

pub fn decide(c: Caps) Decision {
    if (c.has(.msa)) return .{ .core = .loongson3r4 };
    if (c.has(.mmi)) return .{ .core = .loongson3r3 };
    return .{ .core = .mips64_generic };
}

pub fn detect() CoreType {
    return decide(probe()).core;
}

pub fn detectNoted() Decision {
    return decide(probe());
}

test "the cpucfg instruction word" {
    const w = cpucfgWord(4, 2);
    try std.testing.expectEqual(@as(u32, 0xc8881118), w);
    try std.testing.expectEqual(@as(u32, 4), w >> 21 & 0x1f);
    try std.testing.expectEqual(@as(u32, 2), w >> 11 & 0x1f);
    try std.testing.expectEqual(@as(u32, 0xc8080118), w & ~@as(u32, (0x1f << 21) | (0x1f << 11)));
}

test "cpucfg word 1 picks the core" {
    try std.testing.expectEqual(CoreType.loongson3r4, decide(capsFromCfg1(MSA_MASK | MMI_MASK)).core);
    try std.testing.expectEqual(CoreType.loongson3r3, decide(capsFromCfg1(MMI_MASK)).core);
    try std.testing.expectEqual(CoreType.mips64_generic, decide(capsFromCfg1(0)).core);
}

test "cpuinfo model names pick the core" {
    try std.testing.expectEqual(CoreType.loongson3r3, decide(capsFromCpuinfo("system type\t: x\nmodel name\t: Loongson-3A3000 @ 1.4GHz\n")).core);
    try std.testing.expectEqual(CoreType.loongson3r3, decide(capsFromCpuinfo("model name\t: Loongson-3B3000\n")).core);
    try std.testing.expectEqual(CoreType.loongson3r4, decide(capsFromCpuinfo("model name\t: Loongson-3A4000 @ 1.8GHz\n")).core);
    try std.testing.expectEqual(CoreType.mips64_generic, decide(capsFromCpuinfo("model name\t: MIPS P5600\n")).core);
    try std.testing.expectEqual(CoreType.mips64_generic, decide(capsFromCpuinfo("")).core);
}

test "the Loongson cores need their SIMD" {
    try std.testing.expect(runnable(.loongson3r3, .{ .features = .initOne(.mmi) }));
    try std.testing.expect(!runnable(.loongson3r4, .{ .features = .initOne(.mmi) }));
    try std.testing.expect(runnable(.loongson3r4, .{ .features = .initMany(&.{ .mmi, .msa }) }));
    try std.testing.expect(runnable(.mips64_generic, .{}));
}

test "detect runs on MIPS targets" {
    if (comptime !sysinfo.cpu.arch.isMIPS()) return error.SkipZigTest;
    _ = detect();
}
