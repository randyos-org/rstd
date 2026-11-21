//! AArch64 CPU detection, mirroring `driver/others/dynamic_arm64.c`; each
//! OS's probe is in design.md (Core detection). `modelCore` is the part table.

const std = @import("std");
const coretypes = @import("../coretypes.zig");
const select = @import("../select.zig");
const sysinfo = @import("builtin");
const osinfo = @import("../osinfo.zig");

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
/// emulate it; calling this unguarded on a kernel without that support
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
    return osinfo.processorFeature(@backingInt(f));
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
        const midr = osinfo.parseHex(midr_text) orelse return false;
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
    const n = if (possible.len != 0) osinfo.cpuCountFromPossible(possible) else 4096;
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
        return osinfo.readFile(path, buf);
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
pub const cpu0_key = "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0";

pub fn probe() Caps {
    switch (comptime sysinfo.target.os.tag) {
        .driverkit, .ios, .maccatalyst, .macos, .tvos, .visionos, .watchos => {
            return .{ .features = darwinFeatures(osinfo.sysctlInt("hw.optional.arm.FEAT_SME"), osinfo.sysctlInt("hw.optional.arm.FEAT_SME_F64F64")) };
        },
        .windows => return windowsCaps(&processorFeature, osinfo.registryQword(cpu0_key, "CP 4000")),
        // no OS to ask: the ID registers, readable at EL1 and up
        .freestanding, .uefi => return idRegisters(),
        else => return hwcapCaps(osinfo.hwcap(), osinfo.hwcap2(), LiveOs{}),
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
    // the reference gives SME without SVE `vortexm4` (`dynamic_arm64.c:544-548`);
    // its object is built for the M4, so other CPUs take the SVE test (T18)
    if (c.has(.sme) and c.has(.sve)) return .{ .core = .armv9sme };
    if (c.has(.sve)) return .{ .core = .armv8sve };
    return .{ .core = .armv8 };
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
