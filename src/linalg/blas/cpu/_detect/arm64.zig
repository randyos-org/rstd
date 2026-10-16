//! AArch64 CPU detection, mirroring `driver/others/dynamic_arm64.c`.
//!
//! Unlike x86, AArch64 cannot identify itself unaided: `MIDR_EL1` is a
//! system register that faults in user mode unless the kernel emulates the
//! read, which Linux and FreeBSD advertise through `AT_HWCAP & HWCAP_CPUID`.
//! So there are separate paths: Darwin (Apple silicon, answered by `sysctl`),
//! Linux/Android (`MIDR_EL1`, with `/sys/.../midr_el1` for every CPU when the
//! kernel lacks the emulation), Windows (the register's copy in the registry),
//! and everything else (generic). `modelCore` is the part table.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const select = @import("../_select.zig");
const sysinfo = @import("builtin");
const auxv = @import("../_auxv.zig");

pub const CoreType = coretypes.Arm64CoreType;
pub const Decision = select.Decision(CoreType);

pub const HWCAP_CPUID: usize = 1 << 11;
pub const HWCAP_SVE: usize = 1 << 22;
pub const HWCAP2_SVE2: usize = 1 << 1;
pub const HWCAP2_BF16: usize = 1 << 14;
pub const HWCAP2_SME: usize = 1 << 23;

pub const Feature = enum {
    /// Apple silicon, which has no `MIDR_EL1` and answers by `sysctl`
    apple,
    /// `implementer` and `part` are valid
    midr,
    sve,
    sve2,
    sme,
    bf16,
};

/// What the decisions read.
pub const Caps = struct {
    features: std.EnumSet(Feature) = .empty,
    implementer: u32 = 0,
    part: u32 = 0,

    pub fn has(c: Caps, f: Feature) bool {
        return c.features.contains(f);
    }

    /// Every one of `fs`, the set built at compile time.
    pub fn hasAll(c: Caps, comptime fs: []const Feature) bool {
        return c.features.supersetOf(comptime std.EnumSet(Feature).initMany(fs));
    }
};

pub const Midr = struct {
    implementer: u32,
    variant: u32,
    architecture: u32,
    part: u32,
    revision: u32,

    pub fn fromRaw(raw: u64) Midr {
        return .{
            .implementer = @truncate((raw >> 24) & 0xff),
            .variant = @truncate((raw >> 20) & 0xf),
            .architecture = @truncate((raw >> 16) & 0xf),
            .part = @truncate((raw >> 4) & 0xfff),
            .revision = @truncate(raw & 0xf),
        };
    }
};

/// Reads `MIDR_EL1`. Only legal once `HWCAP_CPUID` says the kernel will
/// emulate it -- calling this unguarded on a kernel without that support
/// raises SIGILL.
pub fn readMidr() Midr {
    const raw = asm volatile ("mrs %[out], midr_el1"
        : [out] "=r" (-> u64),
    );
    return .fromRaw(raw);
}

// ---- probes ----------------------------------------------------------------

/// An ID register, by its `S3_0_C0_Cn_m` name so the assembler needs no
/// feature (SVE, SME) to accept it.
fn idRegister(comptime name: []const u8) u64 {
    return asm volatile ("mrs %[out], " ++ name
        : [out] "=r" (-> u64),
    );
}

/// What the ID registers say, with no OS to ask: `ID_AA64PFR0_EL1.SVE`,
/// `ID_AA64ZFR0_EL1.SVEver`, `ID_AA64ISAR1_EL1.BF16`, `ID_AA64PFR1_EL1.SME`.
fn idRegisters() Caps {
    const sve = (idRegister("s3_0_c0_c4_0") >> 32) & 0xf != 0;
    const m = readMidr();
    return .{
        .features = .init(.{
            .midr = true,
            .sve = sve,
            // ZFR0 reads as zero without SVE
            .sve2 = sve and idRegister("s3_0_c0_c4_4") & 0xf != 0,
            .bf16 = (idRegister("s3_0_c0_c6_1") >> 44) & 0xf != 0,
            .sme = (idRegister("s3_0_c0_c4_1") >> 24) & 0xf != 0,
        }),
        .implementer = m.implementer,
        .part = m.part,
    };
}

/// A hex number as the sysfs files print it ("0x00000000411fd070"), the low 32
/// bits, as `strtoul` into an int gives.
pub fn parseHex(text: []const u8) ?u32 {
    var t = std.mem.trim(u8, text, " \t\r\n\x00");
    if (std.mem.startsWith(u8, t, "0x") or std.mem.startsWith(u8, t, "0X")) t = t[2..];
    var end: usize = 0;
    while (end < t.len and std.ascii.isHex(t[end])) end += 1;
    if (end == 0) return null;
    return @truncate(std.fmt.parseInt(u64, t[0..end], 16) catch return null);
}

/// How many CPUs `/sys/devices/system/cpu/possible` ("0-7", or "0") names.
pub fn cpuCountFromPossible(text: []const u8) usize {
    const t = std.mem.trim(u8, text, " \t\r\n\x00");
    if (!std.mem.startsWith(u8, t, "0-")) return 1;
    var end: usize = 2;
    while (end < t.len and std.ascii.isDigit(t[end])) end += 1;
    return (std.fmt.parseInt(usize, t[2..end], 10) catch 0) + 1;
}

/// The per-CPU sysfs pass for kernels without the `MIDR_EL1` emulation, one CPU
/// at a time: the implementer of the last CPU, and the part of the last CPU
/// whose capacity is 1000 or more (the big cores; the scale tops out at 1024),
/// else of the last CPU. `dynamic_arm64.c` reads `cpu_capability` as hex, a
/// file Linux does not have; `cpuid_arm64.c` reads `cpu_capacity` as decimal,
/// which is what Linux has. Its counts of low, mid and high CPUs are unused.
/// The implementer and part a sysfs scan settles on.
pub const Part = struct { implementer: u32, part: u32 };

pub const SysfsScan = struct {
    implementer: u32 = 0,
    last_part: u32 = 0,
    fast_part: u32 = 0,
    count: usize = 0,

    /// One CPU's `midr_el1` and `cpu_capacity` file text (empty if absent).
    /// False when the MIDR cannot be read, where the reference gives up.
    pub fn add(s: *SysfsScan, midr_text: []const u8, capacity_text: []const u8) bool {
        const midr = parseHex(midr_text) orelse return false;
        s.implementer = (midr >> 24) & 0xff;
        s.last_part = (midr >> 4) & 0xfff;
        const cap = std.mem.trim(u8, capacity_text, " \t\r\n\x00");
        if (std.fmt.parseInt(u32, cap, 10)) |c| {
            if (c >= 1000) s.fast_part = s.last_part;
        } else |_| {}
        s.count += 1;
        return true;
    }

    pub fn result(s: SysfsScan) ?Part {
        if (s.count == 0) return null;
        return .{ .implementer = s.implementer, .part = if (s.fast_part != 0) s.fast_part else s.last_part };
    }
};

/// All CPUs' MIDR from sysfs; null if a CPU's file is missing.
fn midrFromSysfs() ?Part {
    var buf: [64]u8 = undefined;
    var cap: [16]u8 = undefined;
    const possible = auxv.readFile("/sys/devices/system/cpu/possible", &buf);
    // without `possible`, every CPU until the first without a MIDR file
    const n = if (possible.len != 0) cpuCountFromPossible(possible) else 4096;
    var scan: SysfsScan = .{};
    for (0..n) |i| {
        var path: [96]u8 = undefined;
        const midr_path = std.mem.printSentinel(&path, "/sys/devices/system/cpu/cpu{d}/regs/identification/midr_el1", .{i}, 0) catch return null;
        const midr = auxv.readFile(midr_path, &buf);
        if (midr.len == 0 and possible.len == 0) break;
        const cap_path = std.mem.printSentinel(&path, "/sys/devices/system/cpu/cpu{d}/cpu_capacity", .{i}, 0) catch return null;
        if (!scan.add(midr, auxv.readFile(cap_path, &cap))) return null;
    }
    return scan.result();
}

/// `MIDR_EL1`'s copy in the registry, Windows' only way to the CPU's identity.
fn midrFromRegistry() ?u64 {
    if (comptime sysinfo.os.tag != .windows) return null;
    const advapi = struct {
        extern "advapi32" fn RegGetValueW(hkey: usize, subkey: [*:0]const u16, value: [*:0]const u16, flags: u32, kind: ?*u32, data: ?*anyopaque, size: ?*u32) callconv(.winapi) i32;
    };
    const HKEY_LOCAL_MACHINE: usize = 0x80000002;
    const RRF_RT_REG_QWORD: u32 = 0x40;
    var raw: u64 = 0;
    var size: u32 = @sizeOf(u64);
    const rc = advapi.RegGetValueW(
        HKEY_LOCAL_MACHINE,
        std.unicode.utf8ToUtf16LeStringLiteral("HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0"),
        std.unicode.utf8ToUtf16LeStringLiteral("CP 4000"),
        RRF_RT_REG_QWORD,
        null,
        &raw,
        &size,
    );
    return if (rc == 0) raw else null;
}

pub fn probe() Caps {
    var c: Caps = .{};
    switch (comptime sysinfo.os.tag) {
        .driverkit, .ios, .maccatalyst, .macos, .tvos, .visionos, .watchos => {
            c.features.insert(.apple);
            c.features.setPresent(.sme, (auxv.sysctlInt("hw.optional.arm.FEAT_SME") orelse 0) != 0);
            return c;
        },
        .windows => {
            const PF = std.os.windows.PF;
            c.features = .init(.{
                .sve = auxv.processorFeature(@backingInt(PF.ARM_SVE_INSTRUCTIONS_AVAILABLE)),
                .sve2 = auxv.processorFeature(@backingInt(PF.ARM_SVE2_INSTRUCTIONS_AVAILABLE)),
                .sme = auxv.processorFeature(@backingInt(PF.ARM_SME_INSTRUCTIONS_AVAILABLE)),
                .bf16 = auxv.processorFeature(@backingInt(PF.ARM_V86_BF16_INSTRUCTIONS_AVAILABLE)),
            });
            if (midrFromRegistry()) |raw| {
                const m = Midr.fromRaw(raw);
                c.features.insert(.midr);
                c.implementer = m.implementer;
                c.part = m.part;
            }
            return c;
        },
        // no OS to ask: the ID registers, readable at EL1 and up
        .freestanding, .uefi => return idRegisters(),
        else => {},
    }
    const hw = auxv.hwcap();
    const hw2 = auxv.hwcap2();
    c.features = .init(.{
        .sve = hw & HWCAP_SVE != 0,
        .sve2 = hw2 & HWCAP2_SVE2 != 0,
        .sme = hw2 & HWCAP2_SME != 0,
        .bf16 = hw2 & HWCAP2_BF16 != 0,
    });
    if (hw & HWCAP_CPUID != 0) {
        const m = readMidr();
        c.features.insert(.midr);
        c.implementer = m.implementer;
        c.part = m.part;
    } else if (sysinfo.os.tag == .linux) {
        if (midrFromSysfs()) |m| {
            c.features.insert(.midr);
            c.implementer = m.implementer;
            c.part = m.part;
        }
    }
    return c;
}

// ---- core type selection ---------------------------------------------------

const neoversen1_fallback = "Your OS does not support SVE instructions. OpenBLAS is using Neoverse N1 kernels as a fallback, which may give poorer performance.";

/// The SVE cores want SVE in the OS, else Neoverse N1's kernels (`dynamic_arm64.c`).
fn sveOr(c: Caps, core: CoreType) Decision {
    if (!c.has(.sve)) return .{ .core = .neoversen1, .note = neoversen1_fallback };
    return .{ .core = core };
}

/// `get_coretype`'s `switch (implementer)`: the part table. Null where the
/// reference falls out of the switch to the SME and SVE tests.
pub fn modelCore(c: Caps) ?Decision {
    if (!c.has(.midr)) return null;
    switch (c.implementer) {
        0x41 => switch (c.part) { // Arm
            0xd03 => return .{ .core = .cortexa53 },
            0xd07 => return .{ .core = .cortexa57 },
            0xd08 => return .{ .core = .cortexa72 },
            0xd09, 0xd0a, 0xd0b => return .{ .core = .cortexa73 }, // A73, A75, A76
            0xd0c => return .{ .core = .neoversen1 },
            0xd49 => return sveOr(c, .neoversen2),
            0xd40 => return sveOr(c, .neoversev1),
            0xd4f, 0xd83, 0xd84, 0xd85, 0xd87 => return sveOr(c, .neoversev2), // V2, V3AE, V3, X925, A725
            0xd05 => return .{ .core = .cortexa55 },
            else => {},
        },
        0x42 => if (c.part == 0x516) return .{ .core = .thunderx2t99 }, // Broadcom Vulcan
        0x43 => switch (c.part) { // Cavium
            0x0a1 => return .{ .core = .thunderx },
            0x0af => return .{ .core = .thunderx2t99 },
            0x0b8 => return .{ .core = .thunderx3t110 },
            else => {},
        },
        0x46 => if (c.part == 0x001) return .{ .core = .a64fx }, // Fujitsu
        0x48 => if (c.part == 0xd01) return .{ .core = .tsv110 }, // HiSilicon
        0x50 => if (c.part == 0x000) return .{ .core = .emag8180 }, // Applied Micro
        0xc0 => switch (c.part) { // Ampere
            0xac3, 0xac4 => return .{ .core = .neoversen1 },
            else => {},
        },
        0x51 => if (c.part == 0xc00) return .{ .core = .falkor }, // Qualcomm
        0x61 => return .{ .core = if (c.has(.sme)) .vortexm4 else .vortex }, // Apple
        else => {},
    }
    return null;
}

/// The core for `c`: Darwin's two, else the part table, else the SME and SVE
/// tests, else generic ARMv8, which the reference reports when it cannot tell.
pub fn decide(c: Caps) Decision {
    return select.guard(@This(), pick(c), c);
}

fn pick(c: Caps) Decision {
    if (c.has(.apple)) return .{ .core = if (c.has(.sme)) .vortexm4 else .vortex };
    if (modelCore(c)) |d| return d;
    if (c.has(.sme)) return .{ .core = if (c.has(.sve)) .armv9sme else .vortexm4 };
    if (c.has(.sve)) return .{ .core = .armv8sve };
    return .{ .core = .armv8 };
}

/// Whether the kernels of `core` run on this CPU.
pub fn supportsCore(core: CoreType) bool {
    return runnable(core, probe());
}

/// Whether the kernels of `core` run on a CPU with `c`'s extensions. The ISA
/// level a core is built for (`v8_2a`...) cannot be told apart by the probes; a
/// part table that names a core only names one the CPU's generation has.
pub fn runnable(core: CoreType, c: Caps) bool {
    return switch (core) {
        .a64fx, .armv8sve, .cortexa510, .cortexa710, .cortexx2 => c.has(.sve),
        .neoversev1 => c.hasAll(&.{ .sve, .bf16 }),
        .neoversen2, .neoversev2 => c.hasAll(&.{ .sve, .sve2, .bf16 }),
        .armv9sme => c.hasAll(&.{ .sve, .sve2, .sme }),
        .vortexm4 => c.has(.sme) or c.has(.apple),
        .unknown => false,
        else => true,
    };
}

/// The core to fall back to when the CPU or OS lacks what `core` needs: its own
/// class first (Neoverse to Neoverse N1, as the reference does without SVE),
/// then generic.
pub fn parent(core: CoreType) ?CoreType {
    return switch (core) {
        .neoversev2 => .neoversen2,
        .neoversen2 => .neoversev1,
        .neoversev1, .a64fx, .armv9sme => .armv8sve,
        .armv8sve, .vortex => .neoversen1,
        .vortexm4 => .vortex,
        .cortexa55 => .cortexa53,
        .armv8, .unknown => null,
        else => .armv8,
    };
}

pub fn detect() CoreType {
    return decide(probe()).core;
}

pub fn detectNoted() Decision {
    return decide(probe());
}

// ---- tests -----------------------------------------------------------------

fn arm(part: u32) Caps {
    return .{ .features = .initOne(.midr), .implementer = 0x41, .part = part };
}

test "Arm parts" {
    const sve: Caps = .{ .features = .initMany(&.{ .midr, .sve, .sve2, .bf16 }), .implementer = 0x41 };
    var c = arm(0xd03);
    try std.testing.expectEqual(CoreType.cortexa53, decide(c).core);
    c.part = 0xd05;
    try std.testing.expectEqual(CoreType.cortexa55, decide(c).core);
    c.part = 0xd07;
    try std.testing.expectEqual(CoreType.cortexa57, decide(c).core);
    c.part = 0xd08;
    try std.testing.expectEqual(CoreType.cortexa72, decide(c).core);
    for ([_]u32{ 0xd09, 0xd0a, 0xd0b }) |p| {
        c.part = p;
        try std.testing.expectEqual(CoreType.cortexa73, decide(c).core);
    }
    c.part = 0xd0c;
    try std.testing.expectEqual(CoreType.neoversen1, decide(c).core);
    var s = sve;
    s.part = 0xd40;
    try std.testing.expectEqual(CoreType.neoversev1, decide(s).core);
    s.part = 0xd49;
    try std.testing.expectEqual(CoreType.neoversen2, decide(s).core);
    for ([_]u32{ 0xd4f, 0xd83, 0xd84, 0xd85, 0xd87 }) |p| {
        s.part = p;
        try std.testing.expectEqual(CoreType.neoversev2, decide(s).core);
    }
}

test "an SVE part on an OS without SVE gets Neoverse N1 and a warning" {
    for ([_]u32{ 0xd40, 0xd49, 0xd4f }) |p| {
        const d = decide(arm(p));
        try std.testing.expectEqual(CoreType.neoversen1, d.core);
        try std.testing.expect(d.note != null);
    }
}

test "other implementers" {
    const t = struct {
        fn core(imp: u32, part: u32) CoreType {
            return decide(.{ .features = .initMany(&.{ .midr, .sve, .sve2, .bf16 }), .implementer = imp, .part = part }).core;
        }
    };
    try std.testing.expectEqual(CoreType.thunderx2t99, t.core(0x42, 0x516)); // Vulcan
    try std.testing.expectEqual(CoreType.thunderx, t.core(0x43, 0x0a1));
    try std.testing.expectEqual(CoreType.thunderx2t99, t.core(0x43, 0x0af));
    try std.testing.expectEqual(CoreType.thunderx3t110, t.core(0x43, 0x0b8));
    try std.testing.expectEqual(CoreType.a64fx, t.core(0x46, 0x001));
    try std.testing.expectEqual(CoreType.tsv110, t.core(0x48, 0xd01));
    try std.testing.expectEqual(CoreType.emag8180, t.core(0x50, 0x000));
    try std.testing.expectEqual(CoreType.neoversen1, t.core(0xc0, 0xac3));
    try std.testing.expectEqual(CoreType.neoversen1, t.core(0xc0, 0xac4));
    try std.testing.expectEqual(CoreType.falkor, t.core(0x51, 0xc00));
    try std.testing.expectEqual(CoreType.vortex, t.core(0x61, 0x022));
    // a part the table lacks: generic
    try std.testing.expectEqual(CoreType.armv8sve, t.core(0x41, 0xd0d));
    try std.testing.expectEqual(CoreType.armv8sve, t.core(0xc0, 0x001));
    try std.testing.expectEqual(CoreType.armv8sve, t.core(0x99, 0x001));
}

test "an unknown part takes SME, then SVE, as the reference does" {
    var c: Caps = .{ .features = .initMany(&.{ .midr, .sve, .sve2 }), .implementer = 0x99, .part = 1 };
    try std.testing.expectEqual(CoreType.armv8sve, decide(c).core);
    c.features.insert(.sme);
    try std.testing.expectEqual(CoreType.armv9sme, decide(c).core);
    c.features.remove(.sve);
    c.features.remove(.sve2);
    try std.testing.expectEqual(CoreType.vortexm4, decide(c).core);
    // without the identity register the same tests apply
    try std.testing.expectEqual(CoreType.armv8sve, decide(.{ .features = .initOne(.sve) }).core);
    try std.testing.expectEqual(CoreType.armv8, decide(.{}).core);
}

test "Apple silicon" {
    try std.testing.expectEqual(CoreType.vortex, decide(.{ .features = .initOne(.apple) }).core);
    try std.testing.expectEqual(CoreType.vortexm4, decide(.{ .features = .initMany(&.{ .apple, .sme }) }).core);
    // Apple's own implementer code, as Linux (Asahi) reports it
    try std.testing.expectEqual(CoreType.vortex, decide(.{ .features = .initOne(.midr), .implementer = 0x61 }).core);
}

test "a Neoverse part short of what its core needs stays in its class" {
    // V2 with SVE but without SVE2 or BF16 in the hwcaps
    const v2: Caps = .{ .features = .initMany(&.{ .midr, .sve }), .implementer = 0x41, .part = 0xd4f };
    try std.testing.expectEqual(CoreType.armv8sve, decide(v2).core);
    // V1 without BF16: generic SVE, not A64FX's or a static target's kernels
    const v1: Caps = .{ .features = .initMany(&.{ .midr, .sve }), .implementer = 0x41, .part = 0xd40 };
    try std.testing.expectEqual(CoreType.armv8sve, decide(v1).core);
    // with SVE2 and BF16 it is V2
    const v2_full: Caps = .{ .features = .initMany(&.{ .midr, .sve, .sve2, .bf16 }), .implementer = 0x41, .part = 0xd4f };
    try std.testing.expectEqual(CoreType.neoversev2, decide(v2_full).core);
    for (std.enums.values(CoreType)) |c| try std.testing.expect(c == .unknown or c == .armv8 or parent(c) != null);
}

test "no CPU is given a core it cannot run" {
    // every combination of the extensions the probes see, on every implementer
    // and part the table names
    const parts = [_][2]u32{
        .{ 0x41, 0xd03 }, .{ 0x41, 0xd05 }, .{ 0x41, 0xd07 }, .{ 0x41, 0xd0c }, .{ 0x41, 0xd40 },
        .{ 0x41, 0xd49 }, .{ 0x41, 0xd4f }, .{ 0x42, 0x516 }, .{ 0x46, 0x001 }, .{ 0x61, 0 },
        .{ 0x99, 0 },
    };
    var mask: u32 = 0;
    while (mask < 16) : (mask += 1) {
        for (parts) |p| {
            const c: Caps = .{
                .features = .init(.{ .midr = true, .sve = mask & 1 != 0, .sve2 = mask & 2 != 0, .sme = mask & 4 != 0, .bf16 = mask & 8 != 0 }),
                .implementer = p[0],
                .part = p[1],
            };
            try std.testing.expect(runnable(decide(c).core, c));
        }
    }
}

test "MIDR_EL1 fields" {
    const m = Midr.fromRaw(0x414fd0c1); // Neoverse N1 r3p1
    try std.testing.expectEqual(@as(u32, 0x41), m.implementer);
    try std.testing.expectEqual(@as(u32, 0xd0c), m.part);
    try std.testing.expectEqual(@as(u32, 4), m.variant);
    try std.testing.expectEqual(@as(u32, 1), m.revision);
}

test "sysfs parsing" {
    try std.testing.expectEqual(@as(?u32, 0x411fd070), parseHex("0x00000000411fd070\n"));
    try std.testing.expectEqual(@as(?u32, 0x1000), parseHex("1000"));
    try std.testing.expectEqual(@as(?u32, null), parseHex("\n"));
    try std.testing.expectEqual(@as(usize, 8), cpuCountFromPossible("0-7\n"));
    try std.testing.expectEqual(@as(usize, 1), cpuCountFromPossible("0\n"));
    try std.testing.expectEqual(@as(usize, 1), cpuCountFromPossible(""));
}

test "without the MIDR emulation the last CPU's identity decides, unless a big one is named" {
    // little Cortex-A55s then big Cortex-A76s, as sysfs prints them
    const a55 = "0x00000000411fd050\n";
    const a76 = "0x00000000414fd0b0\n";
    var plain: SysfsScan = .{};
    for ([_][]const u8{ a55, a55, a76, a76 }) |m| try std.testing.expect(plain.add(m, ""));
    try std.testing.expectEqual(@as(u32, 0x41), plain.result().?.implementer);
    try std.testing.expectEqual(@as(u32, 0xd0b), plain.result().?.part);
    // capacities: the CPU rated 1000 or more wins even when it is not last
    var rated: SysfsScan = .{};
    try std.testing.expect(rated.add(a76, "1024\n"));
    try std.testing.expect(rated.add(a55, "446\n"));
    try std.testing.expect(rated.add(a55, "446\n"));
    try std.testing.expectEqual(@as(u32, 0xd0b), rated.result().?.part);
    // "446" is decimal: as hex (the reference's reading) it would pass 1000
    var little: SysfsScan = .{};
    try std.testing.expect(little.add(a76, "446\n"));
    try std.testing.expect(little.add(a55, "446\n"));
    try std.testing.expectEqual(@as(u32, 0xd05), little.result().?.part);
    // an unreadable MIDR stops the pass
    var bad: SysfsScan = .{};
    try std.testing.expect(!bad.add("", ""));
    try std.testing.expect(bad.result() == null);
}

test "detect runs and returns a concrete core type" {
    if (comptime sysinfo.cpu.arch != .aarch64 and sysinfo.cpu.arch != .aarch64_be) {
        return error.SkipZigTest;
    }
    try std.testing.expect(detect() != .unknown);
}
