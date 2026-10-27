//! 32-bit ARM CPU detection, mirroring `cpuid_arm.c`.
//!
//! The reference has no `dynamic_arm.c` -- 32-bit ARM is the one supported
//! architecture with no runtime dispatch at all, so the kernel set is fixed at
//! build time from `TARGET`, which `getarch` finds from `/proc/cpuinfo` (the CPU
//! part, the model name and the feature list). This does the same at run time,
//! from the part there and `AT_HWCAP`'s VFP and NEON bits, which are the feature
//! list's source.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const select = @import("../_select.zig");
const sysinfo = @import("builtin");
const auxv = @import("../_auxv.zig");

pub const CoreType = coretypes.ArmCoreType;
pub const Decision = select.Decision(CoreType);

pub const HWCAP_VFP: usize = 1 << 6;
pub const HWCAP_NEON: usize = 1 << 12;
pub const HWCAP_VFPv3: usize = 1 << 13;
pub const HWCAP_VFPv3D16: usize = 1 << 14;
pub const HWCAP_VFPv4: usize = 1 << 16;
/// 32 double registers; the VFPv3 kernels use d16 to d31.
pub const HWCAP_VFPD32: usize = 1 << 19;

pub const Caps = struct {
    hwcap: usize = 0,
    /// "CPU part" of `/proc/cpuinfo` (0xc09 Cortex-A9, 0xc0f Cortex-A15); 0 if unknown
    part: u32 = 0,
};

/// The "CPU part" of `/proc/cpuinfo` text (the first, as `cpuid_arm.c` reads it).
pub fn partOf(text: []const u8) u32 {
    const v = auxv.cpuinfoField(text, "CPU part") orelse return 0;
    if (!std.mem.startsWith(u8, v, "0x")) return 0;
    return std.fmt.parseInt(u32, v[2..], 16) catch 0;
}

pub fn probe() Caps {
    if (comptime !sysinfo.target.cpu.arch.isArm()) return .{};
    switch (comptime sysinfo.target.os.tag) {
        // no OS to ask: the ID registers, readable at PL1 and up
        .freestanding, .uefi => return idRegisters(),
        else => {},
    }
    var buf: [4096]u8 = undefined;
    return .{ .hwcap = auxv.hwcap(), .part = partOf(auxv.readFile("/proc/cpuinfo", &buf)) };
}

fn cp15(comptime reg: []const u8) u32 {
    return asm volatile ("mrc p15, 0, %[out], " ++ reg
        : [out] "=r" (-> u32),
    );
}

/// The VFP system registers `vmrs` reads, by number.
const VfpReg = enum(u4) { mvfr1 = 6, mvfr0 = 7, fpexc = 8 };

/// `vmrs r0, <reg>` by its encoding: it assembles where the target has no
/// VFP, as under the soft-float ABI, whose code may still use a VFP the CPU
/// has.
fn vmrs(comptime reg: VfpReg) u32 {
    const inst = if (comptime sysinfo.target.cpu.arch.isThumb()) ".inst.w" else ".inst";
    return asm volatile (std.fmt.comptimePrint("{s} 0x{x}", .{ inst, 0xeef00a10 | (@as(u32, @backingInt(reg)) << 16) })
        : [out] "={r0}" (-> u32),
    );
}

/// `MIDR`'s part, and the hwcaps Linux's `vfp_init` derives from `MVFR0` and
/// `MVFR1`; none where `CPACR` or `FPEXC.EN` leave the FPU off.
fn idRegisters() Caps {
    const part = (cp15("c0, c0, 0") >> 4) & 0xfff;
    if ((cp15("c1, c0, 2") >> 20) & 0x3 == 0 or vmrs(.fpexc) & (1 << 30) == 0) return .{ .part = part };
    const mvfr0 = vmrs(.mvfr0);
    const mvfr1 = vmrs(.mvfr1);
    var hw: usize = 0;
    // single or double precision; VFPv3 is level 2 in either
    if (mvfr0 & 0xff0 != 0) hw |= HWCAP_VFP;
    if (@max((mvfr0 >> 4) & 0xf, (mvfr0 >> 8) & 0xf) >= 2) {
        hw |= HWCAP_VFPv3 | @as(usize, if (mvfr0 & 0xf == 1) HWCAP_VFPv3D16 else HWCAP_VFPD32);
    }
    if (mvfr1 & 0x000fff00 == 0x00011100) hw |= HWCAP_NEON;
    if (mvfr1 & 0xf0000000 == 0x10000000) hw |= HWCAP_VFPv4;
    return .{ .hwcap = hw, .part = part };
}

test "partOf reads the first CPU part line" {
    try std.testing.expectEqual(@as(u32, 0xc09), partOf("processor\t: 0\nCPU architecture: 7\nCPU part\t: 0xc09\nCPU part\t: 0xc0f\n"));
    try std.testing.expectEqual(@as(u32, 0), partOf("Hardware\t: BCM2835\n"));
}

/// Whether the kernels of `core` run with these capabilities: the ARMv7 ones use
/// VFPv3 with 32 registers (NEON implies it), the ARMv6 ones VFP.
pub fn runnable(core: CoreType, c: Caps) bool {
    const vfp3 = c.hwcap & (HWCAP_VFPv3 | HWCAP_VFPv4 | HWCAP_NEON) != 0 and c.hwcap & (HWCAP_VFPD32 | HWCAP_NEON) != 0;
    return switch (core) {
        .armv5 => true,
        .armv6 => c.hwcap & (HWCAP_VFP | HWCAP_VFPv3 | HWCAP_VFPv4 | HWCAP_NEON) != 0,
        .armv7 => vfp3,
        .neon, .cortexa9, .cortexa15 => c.hwcap & HWCAP_NEON != 0,
        .vfpv4 => vfp3 and c.hwcap & HWCAP_VFPv4 != 0,
    };
}

/// `cpuid_arm.c`'s `detect`, from `c`: the Cortex-A9 and A15 parts (ARMv7's
/// kernels with NEON ones for sgemm) where NEON is there, other VFPv3 CPUs
/// ARMv7, VFP-only ones ARMv6, anything else ARMv5, the generic kernels.
pub fn decide(c: Caps) Decision {
    if (c.part == 0xc09 and runnable(.cortexa9, c)) return .{ .core = .cortexa9 };
    if (c.part == 0xc0f and runnable(.cortexa15, c)) return .{ .core = .cortexa15 };
    if (runnable(.armv7, c)) return .{ .core = .armv7 };
    if (runnable(.armv6, c)) return .{ .core = .armv6 };
    return .{ .core = .armv5 };
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

test "decide follows the VFP and NEON capability bits" {
    // Cortex-A9 with NEON, VFPv3-D32
    const a9: Caps = .{ .hwcap = HWCAP_VFP | HWCAP_NEON | HWCAP_VFPv3 | HWCAP_VFPD32, .part = 0xc09 };
    try std.testing.expectEqual(CoreType.cortexa9, decide(a9).core);
    var a15 = a9;
    a15.part = 0xc0f;
    try std.testing.expectEqual(CoreType.cortexa15, decide(a15).core);
    // a Cortex-A9 without NEON stays ARMv7
    try std.testing.expectEqual(CoreType.armv7, decide(.{ .hwcap = HWCAP_VFP | HWCAP_VFPv3 | HWCAP_VFPD32, .part = 0xc09 }).core);
    // Cortex-A7/A15 with VFPv4
    try std.testing.expectEqual(CoreType.armv7, decide(.{ .hwcap = HWCAP_VFP | HWCAP_NEON | HWCAP_VFPv3 | HWCAP_VFPv4 | HWCAP_VFPD32 }).core);
    // Tegra 2: VFPv3 with 16 registers only; the VFPv3 kernels would fault
    try std.testing.expectEqual(CoreType.armv6, decide(.{ .hwcap = HWCAP_VFP | HWCAP_VFPv3 | HWCAP_VFPv3D16 }).core);
    // ARM11
    try std.testing.expectEqual(CoreType.armv6, decide(.{ .hwcap = HWCAP_VFP }).core);
    // no FPU, or no auxiliary vector
    try std.testing.expectEqual(CoreType.armv5, decide(.{}).core);
}

test "no CPU is given a core it cannot run" {
    const bits = [_]usize{ HWCAP_VFP, HWCAP_NEON, HWCAP_VFPv3, HWCAP_VFPv3D16, HWCAP_VFPv4, HWCAP_VFPD32 };
    var mask: usize = 0;
    while (mask < 1 << bits.len) : (mask += 1) {
        var hw: usize = 0;
        for (bits, 0..) |b, i| {
            if ((mask >> @intCast(i)) & 1 != 0) hw |= b;
        }
        try std.testing.expect(runnable(decide(.{ .hwcap = hw }).core, .{ .hwcap = hw }));
    }
}

test "detect runs on 32-bit ARM targets" {
    if (comptime !sysinfo.target.cpu.arch.isArm()) return error.SkipZigTest;
    try std.testing.expect(supportsCore(detect()));
}
