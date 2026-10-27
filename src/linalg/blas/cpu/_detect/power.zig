//! PowerPC CPU detection, mirroring `driver/others/dynamic_power.c`.
//!
//! The reference uses GCC's `__builtin_cpu_is`, which glibc backs with
//! `AT_PLATFORM` from the auxiliary vector, `__builtin_cpu_supports` for
//! POWER10 (ISA 3.1 and MMA), and on AIX or without the builtins the processor
//! version register. Here: `AT_PLATFORM` first; then the PVR where the OS
//! lets a program read it (Linux emulates `mfpvr`); then the ISA bits of
//! `AT_HWCAP` and `AT_HWCAP2`. A CPU none of them names gets the portable
//! kernels, where the reference guesses POWER8. The older cores (POWER3 to 5,
//! 970, Cell, G4, 440) are static targets upstream; they are detected here so
//! no CPU runs another generation's instructions.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const select = @import("../_select.zig");
const sysinfo = @import("builtin");
const auxv = @import("../_auxv.zig");

pub const CoreType = coretypes.PowerCoreType;
pub const Decision = select.Decision(CoreType);

/// The `AT_PLATFORM` names Linux gives (`arch/powerpc/kernel/cputable.c`) that
/// map to a core; tags match the strings.
pub const Platform = enum {
    power11,
    power10,
    power9,
    power8,
    @"power7+",
    power7,
    power6x,
    power6,
    @"power5+",
    power5,
    ppc970,
    @"ppc-cell-be",
    power4,
    power3,
    ppc7450,
    ppc7400,
    ppc440,
    ppc440gp,
};

pub const PPC_FEATURE_64: usize = 0x40000000;
pub const PPC_FEATURE_HAS_ALTIVEC: usize = 0x10000000;
pub const PPC_FEATURE_POWER4: usize = 0x00080000;
pub const PPC_FEATURE_POWER5: usize = 0x00040000;
pub const PPC_FEATURE_POWER5_PLUS: usize = 0x00020000;
pub const PPC_FEATURE_CELL: usize = 0x00010000;
pub const PPC_FEATURE_ARCH_2_05: usize = 0x00001000;
pub const PPC_FEATURE_ARCH_2_06: usize = 0x00000100;
pub const PPC_FEATURE2_ARCH_2_07: usize = 0x80000000;
pub const PPC_FEATURE2_ARCH_3_00: usize = 0x00800000;
pub const PPC_FEATURE2_ARCH_3_1: usize = 0x00040000;
pub const PPC_FEATURE2_MMA: usize = 0x00020000;

/// The core family a CPU belongs to, from any source.
pub const Generation = enum { unknown, power3, power4, ppc970, cell, power5, power6, power8, power9, power10, g4, p440 };

pub const Caps = struct {
    platform: ?Platform = null,
    /// from the PVR
    generation: Generation = .unknown,
    hwcap: usize = 0,
    hwcap2: usize = 0,
};

/// The generation for a PVR value: the reference's table (first match) and the
/// older processors' version numbers.
pub fn fromPvr(pvr: u32) Generation {
    const T = struct { mask: u32, value: u32, gen: Generation };
    const table = [_]T{
        .{ .mask = 0xffffffff, .value = 0x0f000001, .gen = .power5 }, // POWER6 in P5+ mode
        .{ .mask = 0xffff0000, .value = 0x003e0000, .gen = .power6 }, // POWER6X
        .{ .mask = 0xffff0000, .value = 0x003f0000, .gen = .power6 }, // POWER7
        .{ .mask = 0xffff0000, .value = 0x004a0000, .gen = .power6 }, // POWER7+
        .{ .mask = 0xffff0000, .value = 0x004b0000, .gen = .power8 }, // POWER8E
        .{ .mask = 0xffff0000, .value = 0x004c0000, .gen = .power8 }, // POWER8NVL
        .{ .mask = 0xffff0000, .value = 0x004d0000, .gen = .power8 }, // POWER8
        .{ .mask = 0xffffefff, .value = 0x004e0200, .gen = .power9 }, // POWER9 DD2.0
        .{ .mask = 0xffffefff, .value = 0x004e0201, .gen = .power9 }, // POWER9 DD2.1
        .{ .mask = 0xffff0000, .value = 0x004e0000, .gen = .power9 }, // POWER9 DD2.2 and later
        .{ .mask = 0xffff0000, .value = 0x00800000, .gen = .power10 },
        .{ .mask = 0xffff0000, .value = 0x00820000, .gen = .power10 }, // POWER11
        // not in the reference's table
        .{ .mask = 0xffff0000, .value = 0x003a0000, .gen = .power5 },
        .{ .mask = 0xffff0000, .value = 0x003b0000, .gen = .power5 }, // POWER5+
        .{ .mask = 0xffff0000, .value = 0x00390000, .gen = .ppc970 },
        .{ .mask = 0xffff0000, .value = 0x003c0000, .gen = .ppc970 }, // 970FX
        .{ .mask = 0xffff0000, .value = 0x00440000, .gen = .ppc970 }, // 970MP
        .{ .mask = 0xffff0000, .value = 0x00450000, .gen = .ppc970 }, // 970GX
        .{ .mask = 0xffff0000, .value = 0x00700000, .gen = .cell },
        .{ .mask = 0xffff0000, .value = 0x00350000, .gen = .power4 },
        .{ .mask = 0xffff0000, .value = 0x00380000, .gen = .power4 }, // POWER4+
        .{ .mask = 0xffff0000, .value = 0x00400000, .gen = .power3 },
        .{ .mask = 0xffff0000, .value = 0x00410000, .gen = .power3 }, // POWER3+
        .{ .mask = 0xffff0000, .value = 0x000c0000, .gen = .g4 }, // 7400
        .{ .mask = 0xffff0000, .value = 0x800c0000, .gen = .g4 }, // 7410
        .{ .mask = 0xfff80000, .value = 0x80000000, .gen = .g4 }, // 7450 to 7448 (0x8000 to 0x8004)
    };
    for (table) |t| if (pvr & t.mask == t.value) {
        // the 7450 family is 0x8000 to 0x8004 only
        if (t.mask == 0xfff80000 and pvr >> 16 > 0x8004) return .unknown;
        return t.gen;
    };
    return .unknown;
}

fn fromPlatform(p: Platform) Generation {
    return switch (p) {
        .power11, .power10 => .power10,
        .power9 => .power9,
        .power8 => .power8,
        .@"power7+", .power7, .power6x, .power6 => .power6,
        .@"power5+", .power5 => .power5,
        .ppc970 => .ppc970,
        .@"ppc-cell-be" => .cell,
        .power4 => .power4,
        .power3 => .power3,
        .ppc7450, .ppc7400 => .g4,
        .ppc440, .ppc440gp => .p440,
    };
}

/// The server generation `AT_HWCAP` and `AT_HWCAP2` vouch for.
fn fromHwcap(c: Caps) Generation {
    const h2 = c.hwcap2;
    const h = c.hwcap;
    if (h2 & PPC_FEATURE2_ARCH_3_1 != 0 and h2 & PPC_FEATURE2_MMA != 0) return .power10;
    if (h2 & PPC_FEATURE2_ARCH_3_00 != 0) return .power9;
    if (h2 & PPC_FEATURE2_ARCH_2_07 != 0) return .power8;
    if (h & (PPC_FEATURE_ARCH_2_06 | PPC_FEATURE_ARCH_2_05) != 0) return .power6;
    if (h & (PPC_FEATURE_POWER5 | PPC_FEATURE_POWER5_PLUS) != 0) return .power5;
    if (h & PPC_FEATURE_CELL != 0) return .cell;
    if (h & PPC_FEATURE_POWER4 != 0) return if (h & PPC_FEATURE_HAS_ALTIVEC != 0) .ppc970 else .power4;
    return .unknown;
}

/// `AT_PLATFORM` (a NUL-terminated string such as "power9"), if it is one we
/// know. A chip in compatibility mode reports the mode it is in.
pub fn platform() ?Platform {
    const name = auxv.platform() orelse return null;
    return std.meta.stringToEnum(Platform, name);
}

fn mfpvr() u32 {
    return asm volatile ("mfpvr %[out]"
        : [out] "=r" (-> u32),
    );
}

pub fn probe() Caps {
    if (comptime !sysinfo.target.cpu.arch.isPowerpc()) return .{};
    // no OS to ask: the PVR, readable in supervisor state
    if (comptime sysinfo.target.os.tag == .freestanding) return .{ .generation = fromPvr(mfpvr()) };
    var c: Caps = .{ .platform = platform(), .hwcap = auxv.hwcap(), .hwcap2 = auxv.hwcap2() };
    if (sysinfo.target.os.tag == .linux and c.platform == null) {
        c.generation = fromPvr(mfpvr());
    }
    return c;
}

/// The generation the names, the PVR and the ISA bits settle on. POWER10 needs
/// ISA 3.1 and MMA from the OS (`__builtin_cpu_supports`), else it is POWER9.
fn generationOf(c: Caps) Generation {
    const named = if (c.platform) |p| fromPlatform(p) else c.generation;
    const g = if (named == .unknown) fromHwcap(c) else named;
    if (g == .power10 and c.hwcap2 != 0 and fromHwcap(c) != .power10) return .power9;
    return g;
}

/// Server levels, each running the code of the ones below it: ISA 1 (POWER3),
/// 2.00 (POWER4, 970 and Cell, which add AltiVec), 2.02 (POWER5), 2.05 and
/// 2.06 (POWER6, 7), 2.07, 3.00, 3.1.
fn level(g: Generation) u8 {
    return switch (g) {
        .unknown, .g4, .p440 => 0,
        .power3 => 1,
        .power4, .ppc970, .cell => 2,
        .power5 => 3,
        .power6 => 4,
        .power8 => 5,
        .power9 => 6,
        .power10 => 7,
    };
}

fn hasAltivec(c: Caps, g: Generation) bool {
    return c.hwcap & PPC_FEATURE_HAS_ALTIVEC != 0 or switch (g) {
        .ppc970, .cell, .g4, .power6, .power8, .power9, .power10 => true,
        else => false,
    };
}

/// Whether the kernels of `core` run: a server core needs the CPU's level to be
/// at least its own (and AltiVec for the 970's), G4 needs AltiVec, the 440 is
/// only the 440.
pub fn runnable(core: CoreType, c: Caps) bool {
    const g = generationOf(c);
    const have = @max(level(g), level(fromHwcap(c)));
    return switch (core) {
        .unknown => false,
        .power3 => have >= 1,
        .power4 => have >= 2,
        .ppc970, .cell => have >= 2 and hasAltivec(c, g),
        .power5 => have >= 3,
        .power6 => have >= 4,
        .power8 => have >= 5,
        .power9 => have >= 6,
        .power10 => have >= 7,
        .ppcg4 => g == .g4 or (have == 0 and c.hwcap & PPC_FEATURE_HAS_ALTIVEC != 0),
        .ppc440, .ppc440fp2 => g == .p440,
    };
}

pub fn decide(c: Caps) Decision {
    const d: Decision = .{ .core = switch (generationOf(c)) {
        .power10 => .power10,
        .power9 => .power9,
        .power8 => .power8,
        .power6 => .power6,
        .power5 => .power5,
        .power4 => .power4,
        .power3 => .power3,
        .ppc970 => .ppc970,
        .cell => .cell,
        .g4 => .ppcg4,
        .p440 => .ppc440,
        .unknown => .unknown,
    } };
    if (d.core == .unknown) return d;
    return select.guard(@This(), d, c);
}

/// Whether the kernels of `core` run on this CPU.
pub fn supportsCore(core: CoreType) bool {
    return runnable(core, probe());
}

pub fn detect() CoreType {
    return decide(probe()).core;
}

pub fn detectNoted() Decision {
    return decide(probe());
}

test "AT_PLATFORM names the core" {
    const both = PPC_FEATURE2_ARCH_3_1 | PPC_FEATURE2_MMA;
    const t = struct {
        fn core(p: Platform) CoreType {
            return decide(.{ .platform = p }).core;
        }
    };
    try std.testing.expectEqual(CoreType.power6, t.core(.power7));
    try std.testing.expectEqual(CoreType.power6, t.core(.@"power7+"));
    try std.testing.expectEqual(CoreType.power6, t.core(.power6x));
    try std.testing.expectEqual(CoreType.power5, t.core(.@"power5+"));
    try std.testing.expectEqual(CoreType.power4, t.core(.power4));
    try std.testing.expectEqual(CoreType.power3, t.core(.power3));
    try std.testing.expectEqual(CoreType.ppc970, t.core(.ppc970));
    try std.testing.expectEqual(CoreType.cell, t.core(.@"ppc-cell-be"));
    try std.testing.expectEqual(CoreType.ppcg4, t.core(.ppc7450));
    try std.testing.expectEqual(CoreType.ppc440, t.core(.ppc440gp));
    try std.testing.expectEqual(CoreType.power8, decide(.{ .platform = .power8, .hwcap2 = PPC_FEATURE2_ARCH_2_07 }).core);
    try std.testing.expectEqual(CoreType.power9, decide(.{ .platform = .power9, .hwcap2 = PPC_FEATURE2_ARCH_3_00 }).core);
    try std.testing.expectEqual(CoreType.power10, decide(.{ .platform = .power10, .hwcap2 = both | PPC_FEATURE2_ARCH_3_00 }).core);
    try std.testing.expectEqual(CoreType.power10, decide(.{ .platform = .power11, .hwcap2 = both | PPC_FEATURE2_ARCH_3_00 }).core);
}

test "POWER10 needs ISA 3.1 and MMA, else POWER9" {
    const d = decide(.{ .platform = .power10, .hwcap2 = PPC_FEATURE2_ARCH_3_00 | PPC_FEATURE2_ARCH_3_1 });
    try std.testing.expectEqual(CoreType.power9, d.core);
}

test "a CPU the platform string does not name is told by its ISA bits" {
    try std.testing.expectEqual(CoreType.power8, decide(.{ .hwcap2 = PPC_FEATURE2_ARCH_2_07 }).core);
    try std.testing.expectEqual(CoreType.power9, decide(.{ .hwcap2 = PPC_FEATURE2_ARCH_2_07 | PPC_FEATURE2_ARCH_3_00 }).core);
    try std.testing.expectEqual(CoreType.power10, decide(.{ .hwcap2 = PPC_FEATURE2_ARCH_3_1 | PPC_FEATURE2_MMA | PPC_FEATURE2_ARCH_3_00 }).core);
    try std.testing.expectEqual(CoreType.power6, decide(.{ .hwcap = PPC_FEATURE_ARCH_2_06 }).core);
    try std.testing.expectEqual(CoreType.power5, decide(.{ .hwcap = PPC_FEATURE_POWER5_PLUS }).core);
    try std.testing.expectEqual(CoreType.ppc970, decide(.{ .hwcap = PPC_FEATURE_POWER4 | PPC_FEATURE_HAS_ALTIVEC | PPC_FEATURE_64 }).core);
    try std.testing.expectEqual(CoreType.power4, decide(.{ .hwcap = PPC_FEATURE_POWER4 | PPC_FEATURE_64 }).core);
}

test "a CPU nothing names gets the portable kernels" {
    try std.testing.expectEqual(CoreType.unknown, decide(.{}).core);
    // a 603 or 750: no ISA bits the cores need
    try std.testing.expectEqual(CoreType.unknown, decide(.{ .hwcap = 0x0c000000 }).core);
}

test "the PVR table" {
    try std.testing.expectEqual(Generation.power6, fromPvr(0x003f0201)); // POWER7
    try std.testing.expectEqual(Generation.power8, fromPvr(0x004d0200));
    try std.testing.expectEqual(Generation.power8, fromPvr(0x004b0201));
    try std.testing.expectEqual(Generation.power9, fromPvr(0x004e1202));
    try std.testing.expectEqual(Generation.power9, fromPvr(0x004e0200));
    try std.testing.expectEqual(Generation.power10, fromPvr(0x00800200));
    try std.testing.expectEqual(Generation.power10, fromPvr(0x00820100));
    try std.testing.expectEqual(Generation.ppc970, fromPvr(0x00390202));
    try std.testing.expectEqual(Generation.ppc970, fromPvr(0x00440100)); // 970MP
    try std.testing.expectEqual(Generation.power5, fromPvr(0x003a0201));
    try std.testing.expectEqual(Generation.power5, fromPvr(0x003b0300));
    try std.testing.expectEqual(Generation.g4, fromPvr(0x800c1104)); // 7410
    try std.testing.expectEqual(Generation.g4, fromPvr(0x80040202)); // 7448
    try std.testing.expectEqual(Generation.unknown, fromPvr(0x80050000));
    try std.testing.expectEqual(Generation.unknown, fromPvr(0x00080200)); // 750
    try std.testing.expectEqual(CoreType.power6, decide(.{ .generation = .power6 }).core);
    try std.testing.expectEqual(CoreType.ppc970, decide(.{ .generation = .ppc970 }).core);
    try std.testing.expectEqual(CoreType.ppcg4, decide(.{ .generation = .g4 }).core);
    try std.testing.expectEqual(CoreType.power9, decide(.{ .generation = .power10, .hwcap2 = PPC_FEATURE2_ARCH_3_00 }).core);
    try std.testing.expectEqual(CoreType.power10, decide(.{ .generation = .power10 }).core);
}

test "runnable checks the ISA level" {
    const g5: Caps = .{ .platform = .ppc970, .hwcap = PPC_FEATURE_POWER4 | PPC_FEATURE_HAS_ALTIVEC | PPC_FEATURE_64 };
    try std.testing.expect(runnable(.ppc970, g5));
    try std.testing.expect(runnable(.power4, g5));
    try std.testing.expect(!runnable(.power5, g5));
    try std.testing.expect(!runnable(.power6, g5));
    const p5: Caps = .{ .generation = .power5 };
    try std.testing.expect(runnable(.power5, p5));
    try std.testing.expect(!runnable(.power6, p5));
    // POWER5 has no AltiVec
    try std.testing.expect(!runnable(.ppc970, p5));
    const g4: Caps = .{ .platform = .ppc7450, .hwcap = PPC_FEATURE_HAS_ALTIVEC };
    try std.testing.expect(runnable(.ppcg4, g4));
    try std.testing.expect(!runnable(.power4, g4));
    try std.testing.expect(!runnable(.ppc440, g4));
    try std.testing.expect(runnable(.ppc440, .{ .platform = .ppc440 }));
    try std.testing.expect(!runnable(.power6, .{}));
    try std.testing.expect(!runnable(.power3, .{}));
    // every decision is runnable
    for (std.enums.values(Platform)) |p| {
        const c: Caps = .{ .platform = p };
        try std.testing.expect(runnable(decide(c).core, c));
    }
}

test "detect runs on PowerPC targets" {
    if (comptime !sysinfo.target.cpu.arch.isPowerpc()) return error.SkipZigTest;
    const core = detect();
    if (core != .unknown) try std.testing.expect(supportsCore(core));
}
