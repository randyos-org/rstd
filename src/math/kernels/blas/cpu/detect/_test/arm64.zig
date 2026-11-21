//! `detect/arm64.zig`'s tests: the decisions and the probes' decoding.

const std = @import("std");
const select = @import("../../select.zig");
const arm64 = @import("../arm64.zig");

const PF = std.os.windows.PF;
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
const Caps = arm64.Caps;
const Midr = arm64.Midr;
const IdRegs = arm64.IdRegs;
const idCaps = arm64.idCaps;
const hwcapFeatures = arm64.hwcapFeatures;
const darwinFeatures = arm64.darwinFeatures;
const windowsFeatures = arm64.windowsFeatures;
const SysfsScan = arm64.SysfsScan;
const decide = arm64.decide;
const runnable = arm64.runnable;
const parent = arm64.parent;

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
    const c: Caps = .{ .features = .initMany(&.{ .sve, .sve2, .sme, .bf16, .fcma }) };
    for ([_][]const u8{ "armv9sme", "vortexm4" }) |name| {
        const r = select.choose(arm64, c, name, select.testing.AllLinked(.aarch64, CoreType){}, select.testing.NoLog{});
        try std.testing.expect(r.core != .armv9sme and r.core != .vortexm4);
    }
}

test "a forced core the CPU cannot run falls back within what detection found" {
    const D = arm64;
    // a Cortex-A72 (Raspberry Pi 4) told it is a V1: no ARMv8.2 Neoverse N1
    const a72: Caps = .{ .features = .initOne(.midr), .implementer = 0x41, .part = 0xd08 };
    const r = select.choose(D, a72, "neoversev1", select.testing.AllLinked(.aarch64, CoreType){}, select.testing.NoLog{});
    try std.testing.expect(select.subsetOf(CoreType, r.core, .cortexa72));
    try everyCaps(struct {
        fn f(c: Caps) !void {
            const det = decide(c).core;
            for (std.enums.values(CoreType)) |want| {
                if (want == .unknown) continue;
                const got = select.choose(D, c, @tagName(want), select.testing.AllLinked(.aarch64, CoreType){}, select.testing.NoLog{}).core;
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
pub fn fakePf(comptime on: []const PF) *const fn (PF) bool {
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
