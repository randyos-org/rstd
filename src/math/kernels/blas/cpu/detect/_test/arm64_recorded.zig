//! `detect/arm64.zig`'s tests on recorded machines.

const std = @import("std");
const select = @import("../../select.zig");
const sysinfo = @import("builtin");
const osinfo = @import("../../osinfo.zig");
const arm64 = @import("../arm64.zig");
const decisions = @import("arm64.zig");

const CoreType = arm64.CoreType;
const HWCAP_FCMA = arm64.HWCAP_FCMA;
const HWCAP_CPUID = arm64.HWCAP_CPUID;
const HWCAP_SVE = arm64.HWCAP_SVE;
const HWCAP2_SVE2 = arm64.HWCAP2_SVE2;
const HWCAP2_BF16 = arm64.HWCAP2_BF16;
const HWCAP2_SME = arm64.HWCAP2_SME;
const HWCAP2_SME_F64F64 = arm64.HWCAP2_SME_F64F64;
const Feature = arm64.Feature;
const Features = arm64.Features;
const Midr = arm64.Midr;
const idCaps = arm64.idCaps;
const midrFromSysfs = arm64.midrFromSysfs;
const hwcapCaps = arm64.hwcapCaps;
const windowsCaps = arm64.windowsCaps;
const cpu0_key = arm64.cpu0_key;
const decide = arm64.decide;
const fakePf = decisions.fakePf;

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
    var words = std.mem.tokenizeScalar(u8, osinfo.cpuinfoField(text, "Features") orelse return error.TestUnexpectedResult, ' ');
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
            return std.fmt.parseInt(u64, osinfo.cpuinfoField(s, key) orelse return error.TestUnexpectedResult, 0);
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
    pub const sysfs = true;
    const dir = "/sys/devices/system/cpu/";

    midrs: []const u64,
    capacity: []const u32 = &.{},
    cpu: usize = 0,

    pub fn midr(f: FakeLinux) Midr {
        return .fromRaw(f.midrs[f.cpu]);
    }

    pub fn read(f: FakeLinux, path: []const u8, buf: []u8) []const u8 {
        const out = struct {
            fn print(b: []u8, comptime fmt: []const u8, args: anytype) []const u8 {
                return std.mem.print(b, fmt, args) catch b[0..0];
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
        // the reference gives A76 the A73 kernels (`dynamic_arm64.c:445`)
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
        pub const sysfs = true;
        pub fn midr(_: @This()) Midr {
            unreachable;
        }
        pub fn read(_: @This(), _: []const u8, buf: []u8) []const u8 {
            return buf[0..0];
        }
    };
    const c = hwcapCaps(hw[0], hw[1], Empty{});
    try std.testing.expect(!c.has(.midr));
    try std.testing.expectEqual(CoreType.armv8, decide(c).core);
    // FreeBSD and the rest have no sysfs pass
    const NoSysfs = struct {
        pub const sysfs = false;
        pub fn midr(_: @This()) Midr {
            unreachable;
        }
        pub fn read(_: @This(), _: []const u8, _: []u8) []const u8 {
            unreachable;
        }
    };
    try std.testing.expect(!hwcapCaps(hw[0], hw[1], NoSysfs{}).has(.midr));
}

test "recorded Windows arm64 machines pick OpenBLAS's core" {
    // OpenBLAS's run-time detection gives Windows ARMV8 (`dynamic_arm64.c:362`);
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
    try std.testing.expect(osinfo.registryQword("SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion", "InstallTime") != null);
    // `CP 4000` on arm64 alone; x86's key has none
    const midr = osinfo.registryQword(cpu0_key, "CP 4000");
    try std.testing.expect((midr != null) == (sysinfo.target.cpu.arch == .aarch64));
}

test "detect runs and returns a concrete core type" {
    if (comptime sysinfo.target.cpu.arch != .aarch64 and sysinfo.target.cpu.arch != .aarch64_be) {
        return error.SkipZigTest;
    }
    try std.testing.expect(select.detect(arm64, null) != .unknown);
}
