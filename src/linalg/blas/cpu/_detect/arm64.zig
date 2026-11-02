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

pub const HWCAP_FCMA: usize = 1 << 14;
pub const HWCAP_CPUID: usize = 1 << 11;
pub const HWCAP_SVE: usize = 1 << 22;
pub const HWCAP2_SVE2: usize = 1 << 1;
pub const HWCAP2_BF16: usize = 1 << 14;
pub const HWCAP2_SME: usize = 1 << 23;
pub const HWCAP2_SME_F64F64: usize = 1 << 25;

pub const Feature = enum {
    /// Apple silicon, which has no `MIDR_EL1` and answers by `sysctl`
    apple,
    /// `implementer` and `part` are valid
    midr,
    sve,
    sve2,
    sme,
    /// SME's `f64` outer products (FEAT_SME_F64F64), optional in SME; the
    /// `armv9sme` and `vortexm4` objects are built with it
    sme_f64f64,
    bf16,
    /// ARMv8.3's complex `FCMLA` and `FCADD`, which LLVM may give the
    /// `v8_3a` cores' Zig kernels
    fcma,
};

pub const Features = std.EnumSet(Feature);

/// What the decisions read.
pub const Caps = struct {
    features: Features = .empty,
    implementer: u32 = 0,
    part: u32 = 0,

    pub fn has(c: Caps, f: Feature) bool {
        return c.features.contains(f);
    }

    /// Every one of `fs`, the set built at compile time.
    pub fn hasAll(c: Caps, comptime fs: []const Feature) bool {
        return c.features.supersetOf(comptime Features.initMany(fs));
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

/// The ID registers a bare-metal probe reads (EL1 and up). `ZFR0` and `SMFR0`
/// read as zero without SVE and SME.
pub const IdRegs = struct {
    midr: u64 = 0,
    /// `ID_AA64PFR0_EL1`: SVE [35:32]
    pfr0: u64 = 0,
    /// `ID_AA64PFR1_EL1`: SME [27:24]
    pfr1: u64 = 0,
    /// `ID_AA64ZFR0_EL1`: SVEver [3:0]
    zfr0: u64 = 0,
    /// `ID_AA64ISAR1_EL1`: BF16 [47:44], FCMA [19:16]
    isar1: u64 = 0,
    /// `ID_AA64SMFR0_EL1`: F64F64 [48]
    smfr0: u64 = 0,
};

/// What the ID registers say.
pub fn idCaps(r: IdRegs) Caps {
    const sve = (r.pfr0 >> 32) & 0xf != 0;
    const sme = (r.pfr1 >> 24) & 0xf != 0;
    const m = Midr.fromRaw(r.midr);
    return .{
        .features = .init(.{
            .midr = true,
            .sve = sve,
            .sve2 = sve and r.zfr0 & 0xf != 0,
            .sme = sme,
            .sme_f64f64 = sme and (r.smfr0 >> 48) & 1 != 0,
            .bf16 = (r.isar1 >> 44) & 0xf != 0,
            .fcma = (r.isar1 >> 16) & 0xf != 0,
        }),
        .implementer = m.implementer,
        .part = m.part,
    };
}

/// The ID registers, with no OS to ask.
fn idRegisters() Caps {
    return idCaps(.{
        .midr = idRegister("midr_el1"),
        .pfr0 = idRegister("s3_0_c0_c4_0"),
        .pfr1 = idRegister("s3_0_c0_c4_1"),
        .zfr0 = idRegister("s3_0_c0_c4_4"),
        .isar1 = idRegister("s3_0_c0_c6_1"),
        .smfr0 = idRegister("s3_0_c0_c4_5"),
    });
}

/// What the hwcaps say (Linux, FreeBSD).
pub fn hwcapFeatures(hw: usize, hw2: usize) Features {
    const sme = hw2 & HWCAP2_SME != 0;
    return .init(.{
        .sve = hw & HWCAP_SVE != 0,
        .sve2 = hw2 & HWCAP2_SVE2 != 0,
        .sme = sme,
        .sme_f64f64 = sme and hw2 & HWCAP2_SME_F64F64 != 0,
        .bf16 = hw2 & HWCAP2_BF16 != 0,
        .fcma = hw & HWCAP_FCMA != 0,
    });
}

/// What Darwin's `hw.optional.arm.FEAT_SME` and `FEAT_SME_F64F64` say (null:
/// no such sysctl, as on macOS before SME).
pub fn darwinFeatures(sme: ?u64, sme_f64f64: ?u64) Features {
    const has_sme = (sme orelse 0) != 0;
    return .init(.{ .apple = true, .sme = has_sme, .sme_f64f64 = has_sme and (sme_f64f64 orelse 0) != 0 });
}

const PF = std.os.windows.PF;

/// What `IsProcessorFeaturePresent` says, `present` answering for it. Windows
/// before the SME flags answers false, so no SME core there.
pub fn windowsFeatures(present: *const fn (PF) bool) Features {
    const sme = present(.ARM_SME_INSTRUCTIONS_AVAILABLE);
    return .init(.{
        .sve = present(.ARM_SVE_INSTRUCTIONS_AVAILABLE),
        .sve2 = present(.ARM_SVE2_INSTRUCTIONS_AVAILABLE),
        .sme = sme,
        .sme_f64f64 = sme and present(.ARM_SME_F64F64_INSTRUCTIONS_AVAILABLE),
        .bf16 = present(.ARM_V86_BF16_INSTRUCTIONS_AVAILABLE),
        // no flag of its own: JSCVT, which ARMv8.3 makes mandatory with it
        // (LLVM's `v8.3a` implies both `complxnum` and `jsconv`)
        .fcma = present(.ARM_V83_JSCVT_INSTRUCTIONS_AVAILABLE),
    });
}

fn processorFeature(f: PF) bool {
    return auxv.processorFeature(@backingInt(f));
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

/// The implementer and part a sysfs scan settles on.
pub const Part = struct { implementer: u32, part: u32 };

/// The per-CPU sysfs pass for kernels without the `MIDR_EL1` emulation, one CPU
/// at a time: the implementer of the last CPU, and the part of the last CPU
/// whose capacity is 1000 or more (the big cores; Linux v6.17 prints it in
/// decimal, the biggest core at 1024: `drivers/base/topology.c:226`,
/// `drivers/base/arch_topology.c:254`), else of the last CPU.
/// `dynamic_arm64.c` reads `cpu_capability` as hex, a file Linux does not
/// have; `cpuid_arm64.c` reads `cpu_capacity` as decimal. Its counts of low,
/// mid and high CPUs are unused.
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

/// All CPUs' MIDR from sysfs, `os.read(path, buf)` reading each file (empty
/// if absent); null if a CPU's file is missing.
pub fn midrFromSysfs(os: anytype) ?Part {
    var buf: [64]u8 = undefined;
    var cap: [16]u8 = undefined;
    const possible = os.read("/sys/devices/system/cpu/possible", &buf);
    // without `possible`, every CPU until the first without a MIDR file
    const n = if (possible.len != 0) cpuCountFromPossible(possible) else 4096;
    var scan: SysfsScan = .{};
    for (0..n) |i| {
        var path: [96]u8 = undefined;
        const midr_path = std.mem.printSentinel(&path, "/sys/devices/system/cpu/cpu{d}/regs/identification/midr_el1", .{i}, 0) catch return null;
        const midr = os.read(midr_path, &buf);
        if (midr.len == 0 and possible.len == 0) break;
        const cap_path = std.mem.printSentinel(&path, "/sys/devices/system/cpu/cpu{d}/cpu_capacity", .{i}, 0) catch return null;
        if (!scan.add(midr, os.read(cap_path, &cap))) return null;
    }
    return scan.result();
}

/// What the hwcap OSes (Linux, FreeBSD) say: `hw` and `hw2`, then the
/// identity from `os.midr()` (`readMidr`, legal only under `HWCAP_CPUID`) or,
/// where `os.sysfs` holds, `midrFromSysfs(os)`.
pub fn hwcapCaps(hw: usize, hw2: usize, os: anytype) Caps {
    var c: Caps = .{ .features = hwcapFeatures(hw, hw2) };
    const part: ?Part = if (hw & HWCAP_CPUID != 0) blk: {
        const m: Midr = os.midr();
        break :blk .{ .implementer = m.implementer, .part = m.part };
    } else if (@TypeOf(os).sysfs) midrFromSysfs(os) else null;
    if (part) |p| {
        c.features.insert(.midr);
        c.implementer = p.implementer;
        c.part = p.part;
    }
    return c;
}

/// This process's OS, for `hwcapCaps`.
const LiveOs = struct {
    const sysfs = sysinfo.target.os.tag == .linux;

    fn midr(_: LiveOs) Midr {
        return readMidr();
    }

    fn read(_: LiveOs, path: [*:0]const u8, buf: []u8) []const u8 {
        return auxv.readFile(path, buf);
    }
};

/// What Windows says: `present` for `IsProcessorFeaturePresent`, `midr` the
/// registry's `CP 4000` (null if unread).
pub fn windowsCaps(present: *const fn (PF) bool, midr: ?u64) Caps {
    var c: Caps = .{ .features = windowsFeatures(present) };
    if (midr) |raw| {
        const m = Midr.fromRaw(raw);
        c.features.insert(.midr);
        c.implementer = m.implementer;
        c.part = m.part;
    }
    return c;
}

/// The key holding the first CPU's cached ID registers, as `CP <encoding>`
/// values; `CP 4000` is `MIDR_EL1` (`ARM64_SYSREG(3, 0, 0, 0, 0)` in
/// `winnt.h`; Zig's `lib/std/zig/system/windows.zig` and LLVM's
/// `lib/TargetParser/Host.cpp` read it so).
const cpu0_key = "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0";

/// A `REG_QWORD` under `HKEY_LOCAL_MACHINE`; null if absent or another type.
fn registryQword(comptime subkey: []const u8, comptime value: []const u8) ?u64 {
    if (comptime sysinfo.target.os.tag != .windows) return null;
    const advapi = struct {
        extern "advapi32" fn RegGetValueW(hkey: usize, subkey: [*:0]const u16, value: [*:0]const u16, flags: u32, kind: ?*u32, data: ?*anyopaque, size: ?*u32) callconv(.winapi) i32;
    };
    const HKEY_LOCAL_MACHINE: usize = 0x80000002;
    const RRF_RT_REG_QWORD: u32 = 0x40;
    var raw: u64 = 0;
    var size: u32 = @sizeOf(u64);
    const rc = advapi.RegGetValueW(
        HKEY_LOCAL_MACHINE,
        std.unicode.utf8ToUtf16LeStringLiteral(subkey),
        std.unicode.utf8ToUtf16LeStringLiteral(value),
        RRF_RT_REG_QWORD,
        null,
        &raw,
        &size,
    );
    return if (rc == 0) raw else null;
}

pub fn probe() Caps {
    switch (comptime sysinfo.target.os.tag) {
        .driverkit, .ios, .maccatalyst, .macos, .tvos, .visionos, .watchos => {
            return .{ .features = darwinFeatures(auxv.sysctlInt("hw.optional.arm.FEAT_SME"), auxv.sysctlInt("hw.optional.arm.FEAT_SME_F64F64")) };
        },
        .windows => return windowsCaps(&processorFeature, registryQword(cpu0_key, "CP 4000")),
        // no OS to ask: the ID registers, readable at EL1 and up
        .freestanding, .uefi => return idRegisters(),
        else => return hwcapCaps(auxv.hwcap(), auxv.hwcap2(), LiveOs{}),
    }
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
    // the reference gives SME without SVE `vortexm4` (`dynamic_arm64.c:543-547`);
    // its object is built for the M4, so other CPUs take the SVE test (T18)
    if (c.has(.sme) and c.has(.sve)) return .{ .core = .armv9sme };
    if (c.has(.sve)) return .{ .core = .armv8sve };
    return .{ .core = .armv8 };
}

/// Whether the kernels of `core` run on this CPU.
pub fn supportsCore(core: CoreType) bool {
    return runnable(core, probe());
}

/// Whether the kernels of `core` run on a CPU with `c`'s extensions. The ISA
/// level a core is built for (`v8_2a`...) cannot be told apart by the probes; a
/// part table that names a core only names one the CPU's generation has. But
/// the kernels of the cores built for ARMv8.3 and up hold `FCMLA` (the
/// conjugating complex `trsm` kernels), so those cores want `fcma`; the M1
/// and later have it. The SME cores' objects are built with SME's optional F64F64,
/// as the reference's are, so they want it too. Apple's cores are built for
/// Apple's CPUs (`apple_m1`, `apple_m4`), so they want one.
pub fn runnable(core: CoreType, c: Caps) bool {
    return switch (core) {
        .a64fx, .armv8sve => c.has(.sve),
        .cortexa510, .cortexa710, .cortexx2 => c.hasAll(&.{ .sve, .fcma }),
        .neoversev1 => c.hasAll(&.{ .sve, .bf16, .fcma }),
        .neoversen2, .neoversev2 => c.hasAll(&.{ .sve, .sve2, .bf16, .fcma }),
        .armv9sme => c.hasAll(&.{ .sve, .sve2, .sme, .sme_f64f64, .fcma }),
        .vortexm4 => isApple(c) and c.hasAll(&.{ .sme, .sme_f64f64 }) and (c.has(.fcma) or c.has(.apple)),
        .vortex => isApple(c) and (c.has(.fcma) or c.has(.apple)),
        .thunderx3t110, .ampere1 => c.has(.fcma) or c.has(.apple),
        .unknown => false,
        else => true,
    };
}

/// Apple silicon: Darwin's answer, or Apple's implementer code (Linux, bare metal).
fn isApple(c: Caps) bool {
    return c.has(.apple) or (c.has(.midr) and c.implementer == 0x61);
}

/// `runnable` cannot see the ISA level a core is built for, so a forced core's
/// stand-in must need nothing the detected core does not (`select.choose`).
pub const level_unseen = true;

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
    const sve: Caps = .{ .features = .initMany(&.{ .midr, .sve, .sve2, .bf16, .fcma }), .implementer = 0x41 };
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
            return decide(.{ .features = .initMany(&.{ .midr, .sve, .sve2, .bf16, .fcma }), .implementer = imp, .part = part }).core;
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
    var c: Caps = .{ .features = .initMany(&.{ .midr, .sve, .sve2, .fcma }), .implementer = 0x99, .part = 1 };
    try std.testing.expectEqual(CoreType.armv8sve, decide(c).core);
    c.features.insert(.sme);
    // SME without F64F64: the SME cores' objects would fault
    try std.testing.expectEqual(CoreType.armv8sve, decide(c).core);
    c.features.insert(.sme_f64f64);
    try std.testing.expectEqual(CoreType.armv9sme, decide(c).core);
    c.features.remove(.sve);
    c.features.remove(.sve2);
    // SME without SVE: the reference's `vortexm4` is Apple's alone (T18)
    try std.testing.expectEqual(CoreType.armv8, decide(c).core);
    try std.testing.expect(!runnable(.vortexm4, c));
    try std.testing.expect(!runnable(.vortex, c));
    // without the identity register the same tests apply
    try std.testing.expectEqual(CoreType.armv8sve, decide(.{ .features = .initOne(.sve) }).core);
    try std.testing.expectEqual(CoreType.armv8, decide(.{}).core);
}

test "Apple silicon" {
    try std.testing.expectEqual(CoreType.vortex, decide(.{ .features = .initOne(.apple) }).core);
    try std.testing.expectEqual(CoreType.vortexm4, decide(.{ .features = .initMany(&.{ .apple, .sme, .sme_f64f64 }) }).core);
    try std.testing.expectEqual(CoreType.vortex, decide(.{ .features = .initMany(&.{ .apple, .sme }) }).core);
    try std.testing.expect(!runnable(.vortexm4, .{ .features = .initOne(.apple) }));
    // Apple's own implementer code, as Linux (Asahi) reports it
    try std.testing.expectEqual(CoreType.vortex, decide(.{ .features = .initMany(&.{ .midr, .fcma }), .implementer = 0x61 }).core);
    try std.testing.expectEqual(CoreType.vortexm4, decide(.{ .features = .initMany(&.{ .midr, .fcma, .sme, .sme_f64f64 }), .implementer = 0x61 }).core);
}

test "an ARMv8.3 core without the complex instructions takes generic kernels" {
    // a CPU qemu or a hypervisor names ThunderX3 while hiding FCMA
    const c: Caps = .{ .features = .initOne(.midr), .implementer = 0x43, .part = 0x0b8 };
    try std.testing.expectEqual(CoreType.armv8, decide(c).core);
    try std.testing.expect(!runnable(.vortex, .{}));
    try std.testing.expect(runnable(.vortex, .{ .features = .initOne(.apple) }));
}

test "a Neoverse part short of what its core needs stays in its class" {
    // V2 with SVE but without SVE2 or BF16 in the hwcaps
    const v2: Caps = .{ .features = .initMany(&.{ .midr, .sve }), .implementer = 0x41, .part = 0xd4f };
    try std.testing.expectEqual(CoreType.armv8sve, decide(v2).core);
    // V1 without BF16: generic SVE, not A64FX's or a static target's kernels
    const v1: Caps = .{ .features = .initMany(&.{ .midr, .sve }), .implementer = 0x41, .part = 0xd40 };
    try std.testing.expectEqual(CoreType.armv8sve, decide(v1).core);
    // with SVE2 and BF16 it is V2
    const v2_full: Caps = .{ .features = .initMany(&.{ .midr, .sve, .sve2, .bf16, .fcma }), .implementer = 0x41, .part = 0xd4f };
    try std.testing.expectEqual(CoreType.neoversev2, decide(v2_full).core);
    for (std.enums.values(CoreType)) |c| try std.testing.expect(c == .unknown or c == .armv8 or parent(c) != null);
}

/// Every set of the extensions the probes see, on every implementer and part
/// the table names, then without the identity register.
fn everyCaps(comptime body: fn (Caps) anyerror!void) !void {
    const parts = [_][2]u32{
        .{ 0x41, 0xd03 }, .{ 0x41, 0xd05 }, .{ 0x41, 0xd07 }, .{ 0x41, 0xd0c }, .{ 0x41, 0xd40 },
        .{ 0x41, 0xd49 }, .{ 0x41, 0xd4f }, .{ 0x42, 0x516 }, .{ 0x43, 0x0b8 }, .{ 0x46, 0x001 },
        .{ 0x61, 0 },     .{ 0xc0, 0xac3 }, .{ 0x99, 0 },
    };
    const n = @typeInfo(Feature).@"enum".field_names.len;
    var mask: usize = 0;
    while (mask < 1 << n) : (mask += 1) {
        var fs: Features = .empty;
        for (std.enums.values(Feature), 0..) |f, i| fs.setPresent(f, (mask >> @intCast(i)) & 1 != 0);
        if (!fs.contains(.midr)) {
            try body(.{ .features = fs });
            continue;
        }
        for (parts) |p| try body(.{ .features = fs, .implementer = p[0], .part = p[1] });
    }
}

test "no CPU is given a core it cannot run" {
    try everyCaps(struct {
        fn f(c: Caps) !void {
            try std.testing.expect(runnable(decide(c).core, c));
        }
    }.f);
}

test "no SME core is chosen without F64F64" {
    try everyCaps(struct {
        fn f(c: Caps) !void {
            const core = decide(c).core;
            if (!c.has(.sme_f64f64)) try std.testing.expect(core != .armv9sme and core != .vortexm4);
        }
    }.f);
    // nor by name
    const Linked = struct {
        pub const arch: std.Target.Cpu.Arch = .aarch64;
        pub const fallback: ?CoreType = null;
        pub fn ok(_: @This(), _: CoreType) bool {
            return true;
        }
    };
    const NoLog = struct {
        pub fn warn(_: @This(), _: i32, comptime _: []const u8, _: anytype) void {}
    };
    const c: Caps = .{ .features = .initMany(&.{ .sve, .sve2, .sme, .bf16, .fcma }) };
    for ([_][]const u8{ "armv9sme", "vortexm4" }) |name| {
        const r = select.choose(@This(), c, name, Linked{}, NoLog{});
        try std.testing.expect(r.core != .armv9sme and r.core != .vortexm4);
    }
}

test "a forced core the CPU cannot run falls back within what detection found" {
    const Linked = struct {
        pub const arch: std.Target.Cpu.Arch = .aarch64;
        pub const fallback: ?CoreType = null;
        pub fn ok(_: @This(), _: CoreType) bool {
            return true;
        }
    };
    const NoLog = struct {
        pub fn warn(_: @This(), _: i32, comptime _: []const u8, _: anytype) void {}
    };
    const D = @This();
    // a Cortex-A72 (Raspberry Pi 4) told it is a V1: no ARMv8.2 Neoverse N1
    const a72: Caps = .{ .features = .initOne(.midr), .implementer = 0x41, .part = 0xd08 };
    const r = select.choose(D, a72, "neoversev1", Linked{}, NoLog{});
    try std.testing.expect(select.subsetOf(CoreType, r.core, .cortexa72));
    try everyCaps(struct {
        fn f(c: Caps) !void {
            const det = decide(c).core;
            for (std.enums.values(CoreType)) |want| {
                if (want == .unknown) continue;
                const got = select.choose(D, c, @tagName(want), Linked{}, NoLog{}).core;
                try std.testing.expect(runnable(got, c));
                if (got != want) try std.testing.expect(select.subsetOf(CoreType, got, det));
            }
        }
    }.f);
}

test "the hwcaps decode" {
    const f = hwcapFeatures(HWCAP_SVE | HWCAP_FCMA | HWCAP_CPUID, HWCAP2_SVE2 | HWCAP2_SME | HWCAP2_SME_F64F64 | HWCAP2_BF16);
    try std.testing.expect(f.eql(.initMany(&.{ .sve, .sve2, .sme, .sme_f64f64, .bf16, .fcma })));
    try std.testing.expect(hwcapFeatures(0, HWCAP2_SME).eql(.initOne(.sme)));
    // F64F64 counts only with SME
    try std.testing.expect(hwcapFeatures(0, HWCAP2_SME_F64F64).eql(.empty));
    try std.testing.expect(hwcapFeatures(0, 0).eql(.empty));
}

test "the ID registers decode" {
    // Neoverse V2's MIDR; SVE, SME2 and SMFR0 with FA64, I16I64 and F64F64
    const regs: IdRegs = .{
        .midr = 0x410fd4f0,
        .pfr0 = 1 << 32,
        .pfr1 = 2 << 24,
        .zfr0 = 1,
        .isar1 = (1 << 44) | (1 << 16),
        .smfr0 = (1 << 63) | (0xf << 52) | (1 << 48),
    };
    const c = idCaps(regs);
    try std.testing.expect(c.features.eql(.initMany(&.{ .midr, .sve, .sve2, .sme, .sme_f64f64, .bf16, .fcma })));
    try std.testing.expectEqual(@as(u32, 0x41), c.implementer);
    try std.testing.expectEqual(@as(u32, 0xd4f), c.part);
    // SMFR0 without its F64F64 bit, and SMFR0's bits without SME
    var r = regs;
    r.smfr0 &= ~@as(u64, 1 << 48);
    try std.testing.expect(!idCaps(r).has(.sme_f64f64));
    try std.testing.expect(idCaps(r).has(.sme));
    r = regs;
    r.pfr1 = 0;
    try std.testing.expect(!idCaps(r).has(.sme) and !idCaps(r).has(.sme_f64f64));
    // ZFR0 counts only with SVE
    try std.testing.expect(idCaps(.{ .zfr0 = 1 }).features.eql(.initOne(.midr)));
}

test "Darwin's sysctls decode" {
    try std.testing.expect(darwinFeatures(1, 1).eql(.initMany(&.{ .apple, .sme, .sme_f64f64 })));
    try std.testing.expect(darwinFeatures(1, 0).eql(.initMany(&.{ .apple, .sme })));
    try std.testing.expect(darwinFeatures(1, null).eql(.initMany(&.{ .apple, .sme })));
    try std.testing.expect(darwinFeatures(null, 1).eql(.initOne(.apple)));
    try std.testing.expect(darwinFeatures(null, null).eql(.initOne(.apple)));
}

/// An `IsProcessorFeaturePresent` that answers true for `on` alone.
fn fakePf(comptime on: []const PF) *const fn (PF) bool {
    return &struct {
        fn f(p: PF) bool {
            return std.mem.findScalar(PF, on, p) != null;
        }
    }.f;
}

test "Windows' processor features decode" {
    const all = windowsFeatures(fakePf(&.{
        .ARM_SVE_INSTRUCTIONS_AVAILABLE,        .ARM_SVE2_INSTRUCTIONS_AVAILABLE,     .ARM_SME_INSTRUCTIONS_AVAILABLE,
        .ARM_SME_F64F64_INSTRUCTIONS_AVAILABLE, .ARM_V86_BF16_INSTRUCTIONS_AVAILABLE, .ARM_V83_JSCVT_INSTRUCTIONS_AVAILABLE,
    }));
    try std.testing.expect(all.eql(.initMany(&.{ .sve, .sve2, .sme, .sme_f64f64, .bf16, .fcma })));
    try std.testing.expect(windowsFeatures(fakePf(&.{.ARM_SME_INSTRUCTIONS_AVAILABLE})).eql(.initOne(.sme)));
    // F64F64 counts only with SME; a Windows without the flags knows neither
    try std.testing.expect(windowsFeatures(fakePf(&.{.ARM_SME_F64F64_INSTRUCTIONS_AVAILABLE})).eql(.empty));
    try std.testing.expect(windowsFeatures(fakePf(&.{})).eql(.empty));
}

test "MIDR_EL1 fields" {
    const m = Midr.fromRaw(0x414fd0c1); // Neoverse N1 r4p1
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

// ---- recorded machines -------------------------------------------------------
//
// The OS paths no machine here runs (T47), fed what real machines report.
// MIDR values and ID registers are QEMU v10.1.0's models of Arm's TRMs
// (`target/arm/cpu64.c`, `target/arm/tcg/cpu64.c`); implementer and part
// codes are Linux v7.3's `arch/arm64/include/asm/cputype.h`.

/// `MIDR_CPU_MODEL(imp, part)` (Linux v7.3 `cputype.h`): r0p0, for parts with
/// no recorded revision; detection reads only the implementer and part.
fn cpuModel(imp: u64, part: u64) u64 {
    return (imp << 24) | (0xf << 16) | (part << 4);
}

/// Each feature's hwcap word and bit, and its `/proc/cpuinfo` name (Linux
/// v7.3 `arch/arm64/kernel/cpuinfo.c`, `hwcap_str`).
const hwcap_names = .{
    .{ Feature.sve, "sve", 0, HWCAP_SVE },
    .{ Feature.sve2, "sve2", 1, HWCAP2_SVE2 },
    .{ Feature.sme, "sme", 1, HWCAP2_SME },
    .{ Feature.sme_f64f64, "smef64f64", 1, HWCAP2_SME_F64F64 },
    .{ Feature.bf16, "bf16", 1, HWCAP2_BF16 },
    .{ Feature.fcma, "fcma", 0, HWCAP_FCMA },
};

/// The hwcaps a kernel reports for `fs`, with `HWCAP_CPUID` if `cpuid`.
fn hwcapsOf(fs: Features, cpuid: bool) [2]usize {
    var hw: [2]usize = .{ if (cpuid) HWCAP_CPUID else 0, 0 };
    inline for (hwcap_names) |n| {
        if (fs.contains(n[0])) hw[n[2]] |= n[3];
    }
    return hw;
}

/// The features a `/proc/cpuinfo` block's `Features` line names.
fn cpuinfoFeatures(text: []const u8) !Features {
    var fs: Features = .empty;
    var words = std.mem.tokenizeScalar(u8, auxv.cpuinfoField(text, "Features") orelse return error.TestUnexpectedResult, ' ');
    while (words.next()) |w| {
        inline for (hwcap_names) |n| {
            if (std.mem.eql(u8, w, n[1])) fs.insert(n[0]);
        }
    }
    return fs;
}

/// `MIDR_EL1` from a `/proc/cpuinfo` block; arm64 prints `CPU architecture:
/// 8` for the register's 0xf (`cpuinfo.c:277-282`).
fn cpuinfoMidr(text: []const u8) !u64 {
    const get = struct {
        fn f(s: []const u8, key: []const u8) !u64 {
            return std.fmt.parseInt(u64, auxv.cpuinfoField(s, key) orelse return error.TestUnexpectedResult, 0);
        }
    }.f;
    const imp = try get(text, "CPU implementer");
    const variant = try get(text, "CPU variant");
    const part = try get(text, "CPU part");
    const rev = try get(text, "CPU revision");
    return (imp << 24) | (variant << 20) | (0xf << 16) | (part << 4) | rev;
}

/// A Linux whose sysfs holds `midrs` (one per CPU) and, if any, `capacity`;
/// `midr` answers for CPU `cpu`, as the kernel's emulated `mrs` does for the
/// CPU it runs on. Files print as the kernel prints them: `possible` as a
/// list ("0-3"), `midr_el1` as `0x%016llx` (`cpuinfo.c:334`), `cpu_capacity`
/// in decimal (`drivers/base/arch_topology.c`).
const FakeLinux = struct {
    const sysfs = true;
    const dir = "/sys/devices/system/cpu/";

    midrs: []const u64,
    capacity: []const u32 = &.{},
    cpu: usize = 0,

    fn midr(f: FakeLinux) Midr {
        return .fromRaw(f.midrs[f.cpu]);
    }

    fn read(f: FakeLinux, path: []const u8, buf: []u8) []const u8 {
        const out = struct {
            fn print(b: []u8, comptime fmt: []const u8, args: anytype) []const u8 {
                return std.fmt.bufPrint(b, fmt, args) catch b[0..0];
            }
        }.print;
        if (std.mem.eql(u8, path, dir ++ "possible")) {
            if (f.midrs.len == 1) return out(buf, "0\n", .{});
            return out(buf, "0-{d}\n", .{f.midrs.len - 1});
        }
        if (!std.mem.startsWith(u8, path, dir ++ "cpu")) return buf[0..0];
        const rest = path[dir.len + 3 ..];
        const slash = std.mem.findScalar(u8, rest, '/') orelse return buf[0..0];
        const i = std.fmt.parseInt(usize, rest[0..slash], 10) catch return buf[0..0];
        if (i >= f.midrs.len) return buf[0..0];
        const file = rest[slash + 1 ..];
        if (std.mem.eql(u8, file, "regs/identification/midr_el1")) return out(buf, "0x{x:0>16}\n", .{f.midrs[i]});
        if (std.mem.eql(u8, file, "cpu_capacity") and f.capacity.len != 0) return out(buf, "{d}\n", .{f.capacity[i]});
        return buf[0..0];
    }
};

/// `want` on both identity paths: the emulated `mrs` (on each CPU) and,
/// without `HWCAP_CPUID`, the sysfs pass.
fn expectLinux(want: CoreType, midrs: []const u64, capacity: []const u32, fs: Features) !void {
    for ([_]bool{ true, false }) |cpuid| {
        const hw = hwcapsOf(fs, cpuid);
        for (0..midrs.len) |cpu| {
            const c = hwcapCaps(hw[0], hw[1], FakeLinux{ .midrs = midrs, .capacity = capacity, .cpu = cpu });
            try std.testing.expect(c.has(.midr));
            try std.testing.expectEqual(want, decide(c).core);
            // the sysfs pass reads every CPU itself
            if (!cpuid) break;
        }
    }
}

test "recorded Raspberry Pis pick OpenBLAS's core with and without the MIDR emulation" {
    // Pi 3 (Cortex-A53 r0p4): Arch Linux ARM forum, 2018-02-21,
    // <https://archlinuxarm.org/forum/viewtopic.php?p=57811>
    const pi3 = "processor\t: 0\nFeatures\t: fp asimd evtstrm crc32 cpuid\nCPU implementer\t: 0x41\n" ++
        "CPU architecture: 8\nCPU variant\t: 0x0\nCPU part\t: 0xd03\nCPU revision\t: 4\n";
    // Pi 4 (Cortex-A72 r0p3), Linux 6.12.75+rpt-rpi-v8:
    // <https://www.t2linux.com/hardware/board/Raspberry/Pi4/>
    const pi4 = "processor\t: 0\nFeatures\t: fp asimd evtstrm crc32 cpuid\nCPU implementer\t: 0x41\n" ++
        "CPU architecture: 8\nCPU variant\t: 0x0\nCPU part\t: 0xd08\nCPU revision\t: 3\n";
    // Pi 5 (Cortex-A76 r4p1): CNX Software's review, 2023-11-05,
    // <https://www.cnx-software.com/2023/11/05/raspberry-pi-5-review-raspberry-pi-os-bookworm-benchmarks-power-consumption/>
    const pi5 = "processor\t: 0\nFeatures\t: fp asimd evtstrm aes pmull sha1 sha2 crc32 atomics fphp asimdhp " ++
        "cpuid asimdrdm lrcpc dcpop asimddp\nCPU implementer\t: 0x41\nCPU architecture: 8\nCPU variant\t: 0x4\n" ++
        "CPU part\t: 0xd0b\nCPU revision\t: 1\n";
    try std.testing.expectEqual(@as(u64, 0x410fd034), try cpuinfoMidr(pi3));
    try std.testing.expectEqual(@as(u64, 0x410fd083), try cpuinfoMidr(pi4));
    try std.testing.expectEqual(@as(u64, 0x414fd0b1), try cpuinfoMidr(pi5));
    // four identical CPUs each; identical CPUs all rate 1024
    const none: []const u32 = &.{};
    const equal: [4]u32 = @splat(1024);
    for ([_][]const u32{ none, &equal }) |cap| {
        try expectLinux(.cortexa53, &@as([4]u64, @splat(try cpuinfoMidr(pi3))), cap, try cpuinfoFeatures(pi3));
        try expectLinux(.cortexa72, &@as([4]u64, @splat(try cpuinfoMidr(pi4))), cap, try cpuinfoFeatures(pi4));
        // the reference gives A76 the A73 kernels (`dynamic_arm64.c:444`)
        try expectLinux(.cortexa73, &@as([4]u64, @splat(try cpuinfoMidr(pi5))), cap, try cpuinfoFeatures(pi5));
    }
}

test "recorded Neoverse N1 and N2 pick OpenBLAS's core" {
    // QEMU's `neoverse-n1` (r4p1) and `neoverse-n2` (r0p3); the kernel's
    // hwcaps are these ID register fields
    const n1 = idCaps(.{ .midr = 0x414fd0c1, .pfr0 = 0x1100000010111112, .pfr1 = 0x20, .isar1 = 0x100001 });
    const n2 = idCaps(.{ .midr = 0x410fd493, .pfr0 = 0x1201111120111112, .pfr1 = 0x221, .zfr0 = 0x0000110100110021, .isar1 = 0x0011111101211052 });
    try std.testing.expect(n1.features.eql(.initOne(.midr)));
    try std.testing.expect(n2.features.eql(.initMany(&.{ .midr, .sve, .sve2, .bf16, .fcma })));
    try expectLinux(.neoversen1, &.{0x414fd0c1}, &.{}, n1.features);
    try expectLinux(.neoversen2, &.{0x410fd493}, &.{}, n2.features);
    // a kernel without SVE: N1's kernels and the reference's warning
    var no_sve = n2.features;
    no_sve.remove(.sve);
    no_sve.remove(.sve2);
    try expectLinux(.neoversen1, &.{0x410fd493}, &.{}, no_sve);
    const hw = hwcapsOf(no_sve, true);
    try std.testing.expect(decide(hwcapCaps(hw[0], hw[1], FakeLinux{ .midrs = &.{0x410fd493} })).note != null);
}

test "a recorded Apple M1 under Asahi picks vortex" {
    // `apple,icestorm` CPUs 0-3 then `apple,firestorm` 4-7 (Linux v7.3
    // `arch/arm64/boot/dts/apple/t8103.dtsi`); FCMA is ARMv8.3's, and LLVM's
    // `apple-m1` is ARMv8.4-A
    const m1 = [_]u64{
        cpuModel(0x61, 0x022), cpuModel(0x61, 0x022), cpuModel(0x61, 0x022), cpuModel(0x61, 0x022),
        cpuModel(0x61, 0x023), cpuModel(0x61, 0x023), cpuModel(0x61, 0x023), cpuModel(0x61, 0x023),
    };
    try expectLinux(.vortex, &m1, &.{}, .initOne(.fcma));
    // the sysfs pass settles on the last CPU, a Firestorm
    try std.testing.expectEqual(@as(u32, 0x023), midrFromSysfs(FakeLinux{ .midrs = &m1 }).?.part);
}

test "a sysfs without the files gives no identity" {
    const hw = hwcapsOf(.empty, false);
    const Empty = struct {
        const sysfs = true;
        fn midr(_: @This()) Midr {
            unreachable;
        }
        fn read(_: @This(), _: []const u8, buf: []u8) []const u8 {
            return buf[0..0];
        }
    };
    const c = hwcapCaps(hw[0], hw[1], Empty{});
    try std.testing.expect(!c.has(.midr));
    try std.testing.expectEqual(CoreType.armv8, decide(c).core);
    // FreeBSD and the rest have no sysfs pass
    const NoSysfs = struct {
        const sysfs = false;
        fn midr(_: @This()) Midr {
            unreachable;
        }
        fn read(_: @This(), _: []const u8, _: []u8) []const u8 {
            unreachable;
        }
    };
    try std.testing.expect(!hwcapCaps(hw[0], hw[1], NoSysfs{}).has(.midr));
}

test "recorded Windows arm64 machines pick OpenBLAS's core" {
    // OpenBLAS's run-time detection gives Windows ARMV8 (`dynamic_arm64.c:361`);
    // with the registry's MIDR, the part table and SVE test apply as on Linux.
    // Snapdragon X Elite: Oryon (0x51, 0x001), ARMv8.7-A, so JSCVT and BF16
    const elite = windowsCaps(fakePf(&.{ .ARM_V83_JSCVT_INSTRUCTIONS_AVAILABLE, .ARM_V86_BF16_INSTRUCTIONS_AVAILABLE }), cpuModel(0x51, 0x001));
    try std.testing.expectEqual(@as(u32, 0x001), elite.part);
    try std.testing.expectEqual(CoreType.armv8, decide(elite).core);
    // Microsoft SQ1 and SQ2: Kryo 495, ARMv8.2-A Gold (0x804) and Silver (0x805)
    for ([_]u64{ 0x804, 0x805 }) |p| {
        try std.testing.expectEqual(CoreType.armv8, decide(windowsCaps(fakePf(&.{}), cpuModel(0x51, p))).core);
    }
    // Azure Cobalt 100 (GitHub's `windows-11-arm`): Microsoft's (0x6d) N2
    // r0p0; the reference has no 0x6d, so its SVE test decides
    const cobalt = cpuModel(0x6d, 0xd49);
    const sve = windowsCaps(fakePf(&.{
        .ARM_SVE_INSTRUCTIONS_AVAILABLE,      .ARM_SVE2_INSTRUCTIONS_AVAILABLE,
        .ARM_V86_BF16_INSTRUCTIONS_AVAILABLE, .ARM_V83_JSCVT_INSTRUCTIONS_AVAILABLE,
    }), cobalt);
    try std.testing.expectEqual(@as(u32, 0x6d), sve.implementer);
    try std.testing.expectEqual(CoreType.armv8sve, decide(sve).core);
    // a Windows without the SVE flags
    const no_sve = windowsCaps(fakePf(&.{ .ARM_V86_BF16_INSTRUCTIONS_AVAILABLE, .ARM_V83_JSCVT_INSTRUCTIONS_AVAILABLE }), cobalt);
    try std.testing.expectEqual(CoreType.armv8, decide(no_sve).core);
    // the registry unread
    const unread = windowsCaps(fakePf(&.{}), null);
    try std.testing.expect(!unread.has(.midr));
    try std.testing.expectEqual(CoreType.armv8, decide(unread).core);
}

test "the registry read works on Windows" {
    if (comptime sysinfo.target.os.tag != .windows) return error.SkipZigTest;
    // a `REG_QWORD` every Windows 10 and 11 has
    try std.testing.expect(registryQword("SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion", "InstallTime") != null);
    // `CP 4000` on arm64 alone; x86's key has none
    const midr = registryQword(cpu0_key, "CP 4000");
    try std.testing.expect((midr != null) == (sysinfo.target.cpu.arch == .aarch64));
}

test "detect runs and returns a concrete core type" {
    if (comptime sysinfo.target.cpu.arch != .aarch64 and sysinfo.target.cpu.arch != .aarch64_be) {
        return error.SkipZigTest;
    }
    try std.testing.expect(detect() != .unknown);
}
