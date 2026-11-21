//! `detect/x86.zig`'s tests.

const std = @import("std");
const corespec = @import("../../corespec.zig");
const select = @import("../../select.zig");
const sysinfo = @import("builtin");
const x86 = @import("../x86.zig");

const CoreType = x86.CoreType;
const Vendor = x86.Vendor;
const Caps = x86.Caps;
const vendorOf = x86.vendorOf;
const probe = x86.probe;
const runnable = x86.runnable;
const nehalem_fallback = x86.nehalem_fallback;
const sandybridge_fallback = x86.sandybridge_fallback;
const haswell_fallback = x86.haswell_fallback;
const decideWidth = x86.decideWidth;

/// A CPU as the tests describe it: CPUID.1:EAX and the extensions the OS lets it use.
fn cpu(v: Vendor, eax1: u32, level: enum { sse3, sse41, avx, avx2, avx512, bf16, amx }) Caps {
    const n = @backingInt(level);
    return .{
        .vendor = v,
        .eax1 = eax1,
        .max_ext_leaf = 0x80000008,
        .features = .init(.{
            .sse3 = true,
            .ssse3 = true,
            .sse4_1 = n >= 1,
            .avx = n >= 2,
            // the Bulldozer family has FMA3 and FMA4 below AVX2
            .fma = n >= 3 or (v == .amd and n == 2),
            .fma4 = v == .amd and n == 2,
            .avx2 = n >= 3,
            .avx512 = n >= 4,
            .avx512_bf16 = n >= 5,
            .amx_bf16 = n >= 6,
        }),
    };
}

fn expectCore(want: CoreType, c: Caps) !void {
    try std.testing.expectEqual(want, decideWidth(c, true).core);
}

test "Intel models map to their cores" {
    // CPUID.1:EAX values of real parts
    try expectCore(.haswell, cpu(.intel, 0x306c3, .avx2)); // Haswell
    try expectCore(.haswell, cpu(.intel, 0x406f1, .avx2)); // Broadwell-E
    try expectCore(.haswell, cpu(.intel, 0x506e3, .avx2)); // Skylake client
    try expectCore(.haswell, cpu(.intel, 0x906ea, .avx2)); // Coffee Lake
    try expectCore(.skylakex, cpu(.intel, 0x50654, .avx512)); // Skylake-SP
    try expectCore(.cooperlake, cpu(.intel, 0x5065b, .bf16)); // model 5 with bf16
    try expectCore(.skylakex, cpu(.intel, 0x606a6, .avx512)); // Ice Lake SP
    try expectCore(.skylakex, cpu(.intel, 0x706e5, .avx512)); // Ice Lake client
    try expectCore(.skylakex, cpu(.intel, 0x806c1, .avx512)); // Tiger Lake
    try expectCore(.sapphirerapids, cpu(.intel, 0x806f8, .amx)); // Sapphire Rapids
    try expectCore(.cooperlake, cpu(.intel, 0x806f8, .bf16)); // ... without AMX
    try expectCore(.sapphirerapids, cpu(.intel, 0xc06f2, .amx)); // Emerald Rapids
    try expectCore(.cooperlake, cpu(.intel, 0xc06f2, .bf16)); // ... without AMX
    try expectCore(.skylakex, cpu(.intel, 0xc06f2, .avx512));
    try expectCore(.cooperlake, cpu(.intel, 0x90672, .bf16)); // Alder Lake with bf16
    try expectCore(.haswell, cpu(.intel, 0x90672, .avx2)); // Alder Lake
    try expectCore(.sandybridge, cpu(.intel, 0x206a7, .avx)); // Sandy Bridge
    try expectCore(.sandybridge, cpu(.intel, 0x306a9, .avx)); // Ivy Bridge
    try expectCore(.nehalem, cpu(.intel, 0x106a5, .sse41)); // Nehalem
    try expectCore(.nehalem, cpu(.intel, 0x206c2, .sse41)); // Westmere
    try expectCore(.penryn, cpu(.intel, 0x10676, .sse41)); // Penryn
    try expectCore(.dunnington, cpu(.intel, 0x106d1, .sse41)); // Dunnington
    try expectCore(.core2, cpu(.intel, 0x6f6, .sse3)); // Conroe (exmodel 0, model 15)
    try expectCore(.atom, cpu(.intel, 0x106c2, .sse3)); // Diamondville
    try expectCore(.prescott, cpu(.intel, 0xf49, .sse3)); // Prescott
}

test "32-bit-only models stay 32-bit, and 64-bit builds move them to Prescott" {
    const katmai = cpu(.intel, 0x673, .sse3); // Pentium III
    try std.testing.expectEqual(CoreType.katmai, decideWidth(katmai, false).core);
    try std.testing.expectEqual(CoreType.prescott, decideWidth(katmai, true).core);
    const northwood = cpu(.intel, 0xf24, .sse3);
    try std.testing.expectEqual(CoreType.northwood, decideWidth(northwood, false).core);
    const banias = cpu(.intel, 0x6d8, .sse3); // Pentium M
    try std.testing.expectEqual(CoreType.banias, decideWidth(banias, false).core);
    const coppermine = cpu(.intel, 0x686, .sse3);
    try std.testing.expectEqual(CoreType.coppermine, decideWidth(coppermine, false).core);
}

test "an OS that hides AVX state gets the old kernels, with the reference's warning" {
    const sandy = decideWidth(cpu(.intel, 0x206a7, .sse41), true);
    try std.testing.expectEqual(CoreType.nehalem, sandy.core);
    try std.testing.expectEqualStrings(nehalem_fallback, sandy.note.?);
    const hsw_avx = decideWidth(cpu(.intel, 0x306c3, .avx), true);
    try std.testing.expectEqual(CoreType.sandybridge, hsw_avx.core);
    try std.testing.expectEqualStrings(sandybridge_fallback, hsw_avx.note.?);
    const hsw_none = decideWidth(cpu(.intel, 0x306c3, .sse41), true);
    try std.testing.expectEqual(CoreType.nehalem, hsw_none.core);
    try std.testing.expectEqualStrings(nehalem_fallback, hsw_none.note.?);
    // a real Sandy or Ivy Bridge is no fallback: no warning
    for ([_]u32{ 0x206a7, 0x206d7, 0x306a9, 0x306e4 }) |eax| {
        const d = decideWidth(cpu(.intel, eax, .avx), true);
        try std.testing.expectEqual(CoreType.sandybridge, d.core);
        try std.testing.expectEqual(@as(?[]const u8, null), d.note);
    }
    // Skylake-X without the ZMM state keeps AVX2
    const skx = decideWidth(cpu(.intel, 0x50654, .avx2), true);
    try std.testing.expectEqual(CoreType.haswell, skx.core);
    try std.testing.expectEqualStrings(haswell_fallback, skx.note.?);
    // a part without a note
    try std.testing.expect(decideWidth(cpu(.intel, 0x306c3, .avx2), true).note == null);
}

test "AMD families" {
    try expectCore(.opteron, blk: {
        var c = cpu(.amd, 0x20f12, .sse3); // K8
        c.features.remove(.sse3);
        c.ecx1 = 0;
        c.features.remove(.ssse3);
        break :blk c;
    });
    var k8e = cpu(.amd, 0x40f12, .sse3);
    k8e.ecx1 = 1;
    try expectCore(.opteron_sse3, k8e);
    try expectCore(.barcelona, cpu(.amd, 0x100f22, .sse3)); // K10
    try expectCore(.barcelona, cpu(.amd, 0x300f10, .sse3)); // Llano
    try expectCore(.bobcat, cpu(.amd, 0x500f20, .sse3)); // Bobcat
    try expectCore(.bobcat, cpu(.amd, 0x700f01, .sse3)); // Jaguar
    try expectCore(.bulldozer, cpu(.amd, 0x600f12, .avx)); // Zambezi
    try expectCore(.piledriver, cpu(.amd, 0x600f20, .avx)); // Vishera
    try expectCore(.piledriver, cpu(.amd, 0x610f01, .avx)); // Trinity (exmodel 1, model 0)
    try expectCore(.steamroller, cpu(.amd, 0x630f01, .avx)); // Kaveri
    try expectCore(.excavator, cpu(.amd, 0x660f51, .avx)); // family 15h model 65h: A12-9800, PRO A12-8870 (InstLatx64 dumps)
    try expectCore(.excavator, cpu(.amd, 0x660f00, .avx)); // exmodel 6, model 0
    try expectCore(.zen, cpu(.amd, 0x800f82, .avx2)); // Zen+
    try expectCore(.zen, cpu(.amd, 0x870f10, .avx2)); // Zen 2
    try expectCore(.zen, cpu(.amd, 0xa20f10, .avx2)); // Zen 3 (no AVX-512)
    try expectCore(.skylakex, cpu(.amd, 0xa10f11, .avx512)); // Zen 4 reports 19h
    try expectCore(.cooperlake, cpu(.amd, 0xa20f12, .bf16));
    try expectCore(.zen, cpu(.hygon, 0x900f01, .avx2)); // Dhyana
}

test "AMD without AVX uses Barcelona" {
    const d = decideWidth(cpu(.amd, 0x600f20, .sse3), true);
    try std.testing.expectEqual(CoreType.barcelona, d.core);
    try std.testing.expect(d.note != null);
    try expectCore(.barcelona, cpu(.amd, 0x870f10, .sse3));
    // AVX without AVX2 or FMA cannot run Zen's kernels; Sandy Bridge's do
    try expectCore(.sandybridge, cpu(.amd, 0x870f10, .avx));
}

test "Athlon needs 3DNow! and its extension" {
    var c = cpu(.amd, 0x680, .sse3); // family 6
    c.features.remove(.sse3);
    c.features.remove(.ssse3);
    c.edx_ext1 = (1 << 30) | (1 << 31);
    try std.testing.expectEqual(CoreType.athlon, decideWidth(c, false).core);
    // 64-bit builds have none
    try std.testing.expectEqual(CoreType.opteron, decideWidth(c, true).core);
    c.edx_ext1 = 1 << 31;
    try std.testing.expectEqual(CoreType.katmai, decideWidth(c, false).core);
}

test "VIA and Zhaoxin" {
    try expectCore(.nano, cpu(.centaur, 0x6f2, .sse3)); // Nano: family 6 model 15 stepping 2
    try expectCore(.nehalem, cpu(.centaur, 0x6fe, .sse41)); // stepping 14
    try expectCore(.dunnington, cpu(.centaur, 0x507f0, .sse41)); // family 7, exmodel 5
    try expectCore(.zen, cpu(.centaur, 0x507f0, .avx2));
    try expectCore(.zen, cpu(.zhaoxin, 0x107f0, .avx2));
    try expectCore(.nehalem, cpu(.zhaoxin, 0x107f0, .sse41));
    try expectCore(.dunnington, cpu(.zhaoxin, 0x507f0, .sse41));
}

test "an unknown model takes the richest core the CPU and OS support" {
    // Raptor Lake: exmodel 0xb is in no table
    try expectCore(.haswell, cpu(.intel, 0xb06a2, .avx2));
    try expectCore(.skylakex, cpu(.intel, 0xb06a2, .avx512));
    try expectCore(.cooperlake, cpu(.intel, 0xb06a2, .bf16));
    try expectCore(.sapphirerapids, cpu(.intel, 0xb06a2, .amx));
    // Granite Rapids: exmodel 0xa, model 13
    try expectCore(.sapphirerapids, cpu(.intel, 0xa06d1, .amx));
    try expectCore(.sandybridge, cpu(.unknown, 0, .avx));
    try expectCore(.prescott, cpu(.unknown, 0, .sse3));
    // no SSE3 at all: a 64-bit build still gets a core it can run
    var old = cpu(.unknown, 0, .sse3);
    old.features.remove(.sse3);
    old.features.remove(.ssse3);
    try expectCore(.opteron, old);
    try std.testing.expectEqual(CoreType.katmai, decideWidth(old, false).core);
}

test "FMA3 missing under AVX2 does not give the Haswell kernels" {
    var c = cpu(.intel, 0x306c3, .avx2);
    c.features.remove(.fma);
    try expectCore(.sandybridge, c);
}

test "vendorOf reads the leaf 0 strings" {
    const w = struct {
        fn v(s: *const [12]u8, max: u32) Vendor {
            return vendorOf(std.mem.readInt(u32, s[0..4], .little), std.mem.readInt(u32, s[4..8], .little), std.mem.readInt(u32, s[8..12], .little), max);
        }
    };
    try std.testing.expectEqual(Vendor.intel, w.v("GenuineIntel", 0x16));
    try std.testing.expectEqual(Vendor.amd, w.v("AuthenticAMD", 0x10));
    try std.testing.expectEqual(Vendor.hygon, w.v("HygonGenuine", 0xd));
    try std.testing.expectEqual(Vendor.centaur, w.v("CentaurHauls", 0x6));
    try std.testing.expectEqual(Vendor.zhaoxin, w.v("  Shanghai  ", 0xd));
    // unknown strings: Intel when the highest leaf is 0 or has bit 8 or 10
    try std.testing.expectEqual(Vendor.intel, w.v("TransmetaCPU", 0));
    try std.testing.expectEqual(Vendor.intel, w.v("TransmetaCPU", 0x500));
    try std.testing.expectEqual(Vendor.unknown, w.v("TransmetaCPU", 0x10));
}

test "the display model adds the extended field" {
    const skx = cpu(.intel, 0x50654, .avx512);
    try std.testing.expectEqual(@as(u32, 0x55), skx.displayModel());
    const zen2 = cpu(.amd, 0x870f10, .avx2);
    try std.testing.expectEqual(@as(u32, 0x71), zen2.displayModel());
}

test "every core the decisions can return is runnable on a CPU that has what it needs" {
    // the richest 32-bit CPU with 3DNow! runs every core but the FMA4 ones
    var full = cpu(.intel, 0x806f8, .amx);
    full.is64 = false;
    full.edx_ext1 = (1 << 30) | (1 << 31);
    for (std.enums.values(CoreType)) |core| {
        if (core == .unknown) continue;
        const s = corespec.x86(core);
        var fma4 = false;
        for (s.features) |f| fma4 = fma4 or std.mem.eql(u8, f, "fma4");
        try std.testing.expectEqual(!fma4, runnable(core, full));
    }
}

test "64-bit builds run no 32-bit-only core, and Athlon needs 3DNow!" {
    var c = cpu(.intel, 0x806f8, .amx);
    c.edx_ext1 = (1 << 30) | (1 << 31);
    c.is64 = true;
    for ([_]CoreType{ .katmai, .coppermine, .northwood, .banias, .athlon, .p5, .p6, .viac3, .yonah }) |k|
        try std.testing.expect(!runnable(k, c));
    c.is64 = false;
    try std.testing.expect(runnable(.athlon, c));
    c.edx_ext1 = 0;
    try std.testing.expect(!runnable(.athlon, c));
    try std.testing.expect(runnable(.katmai, c));
}

test "feature probes execute without faulting" {
    if (comptime !sysinfo.target.cpu.arch.isX86()) return error.SkipZigTest;
    const c = probe();
    // Nothing to assert about which CPU this is; the point is that every
    // cpuid/xgetbv path runs to completion and reports a coherent result.
    try std.testing.expect(c.family() > 0);
    if (c.has(.avx512)) try std.testing.expect(c.has(.avx2));
    if (c.has(.avx2)) try std.testing.expect(c.has(.avx));
}

test "detect returns a concrete core type the host can run" {
    if (comptime !sysinfo.target.cpu.arch.isX86()) return error.SkipZigTest;
    try std.testing.expect(select.detect(x86, null) != .unknown);
    try std.testing.expect(select.canRun(x86, select.detect(x86, null)));
}
