//! x86 / x86_64 CPU detection, mirroring `driver/others/dynamic.c`, shared
//! by both widths (`cpuid` is the same in either): `x86.zig` and `x86_64.zig`
//! add each width's kernel registry and vtable.
//!
//! This is the only architecture that can detect itself with no help from
//! the OS: `cpuid` and `xgetbv` are unprivileged, so everything here works
//! identically on Windows, Linux, macOS and freestanding. The probes are thin
//! (`probe` fills a `Caps`); `modelCore` and `decide` are pure functions of it,
//! which is what the tests feed.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const corespec = @import("../_corespec.zig");
const select = @import("../_select.zig");
const sysinfo = @import("builtin");

pub const CoreType = coretypes.X86CoreType;
pub const Decision = select.Decision(CoreType);

pub const Vendor = enum { unknown, intel, amd, centaur, hygon, zhaoxin };

/// The coarse SIMD capability tiers that actually change which kernel we
/// want, as opposed to the fine-grained microarchitecture identity above.
/// The reference distinguishes both: two chips can share a tier (Zen and
/// Haswell are both AVX2) and still want different cache blocking.
/// `_unresolved` is only `cachedTier`'s pre-probe state; never returned.
pub const Tier = enum(u8) { _unresolved, baseline, avx2, avx512 };

/// What the decisions read, probed once. The `support_*` flags include the OS
/// having enabled the state (XCR0), as `dynamic.c` does.
pub const Caps = struct {
    vendor: Vendor = .unknown,
    /// CPUID.1:EAX (stepping, model, family, extended model and family)
    eax1: u32 = 0,
    /// CPUID.1:ECX, for SSE3 (the Opteron split)
    ecx1: u32 = 0,
    /// CPUID.80000001h:EDX when the leaf exists; bits 30 and 31 are 3DNow!
    edx_ext1: u32 = 0,
    /// CPUID.80000000h:EAX, the highest extended leaf
    max_ext_leaf: u32 = 0,
    sse3: bool = false,
    ssse3: bool = false,
    sse4_1: bool = false,
    avx: bool = false,
    avx2: bool = false,
    fma3: bool = false,
    fma4: bool = false,
    avx512: bool = false,
    avx512_bf16: bool = false,
    amx_bf16: bool = false,
    /// a 64-bit build, which has no 32-bit-only cores
    is64: bool = @sizeOf(usize) == 8,

    pub fn family(c: Caps) u32 {
        return (c.eax1 >> 8) & 0xf;
    }
    pub fn exfamily(c: Caps) u32 {
        return (c.eax1 >> 20) & 0xff;
    }
    pub fn model(c: Caps) u32 {
        return (c.eax1 >> 4) & 0xf;
    }
    pub fn exmodel(c: Caps) u32 {
        return (c.eax1 >> 16) & 0xf;
    }
    pub fn stepping(c: Caps) u32 {
        return c.eax1 & 0xf;
    }
    /// The display family and model: the extended fields added in as the
    /// manuals say (family 0xf adds the family, 6 and 0xf the model).
    pub fn displayFamily(c: Caps) u32 {
        return if (c.family() == 0xf) c.family() + c.exfamily() else c.family();
    }
    pub fn displayModel(c: Caps) u32 {
        return if (c.family() == 6 or c.family() == 0xf) c.model() + (c.exmodel() << 4) else c.model();
    }
    /// 3DNow! and its extension, which the Athlon test wants both of.
    pub fn has3dnow(c: Caps) bool {
        return c.max_ext_leaf & 0xffff >= 1 and c.edx_ext1 & (1 << 30) != 0 and c.edx_ext1 & (1 << 31) != 0;
    }
};

/// The old summary, kept for callers that want the display family and model.
pub const Features = struct {
    vendor: Vendor = .unknown,
    family: u32 = 0,
    model: u32 = 0,
    stepping: u32 = 0,
    avx: bool = false,
    avx2: bool = false,
    avx512: bool = false,
    avx512_bf16: bool = false,
    amx_bf16: bool = false,

    pub fn tier(f: Features) Tier {
        if (f.avx512) return .avx512;
        if (f.avx2) return .avx2;
        return .baseline;
    }
};

// ---- raw instruction access ------------------------------------------------

pub const CpuidResult = struct { eax: u32, ebx: u32, ecx: u32, edx: u32 };

pub fn cpuid(leaf: u32, subleaf: u32) CpuidResult {
    var eax: u32 = undefined;
    var ebx: u32 = undefined;
    var ecx: u32 = undefined;
    var edx: u32 = undefined;
    asm volatile ("cpuid"
        : [_] "={eax}" (eax),
          [_] "={ebx}" (ebx),
          [_] "={ecx}" (ecx),
          [_] "={edx}" (edx),
        : [_] "{eax}" (leaf),
          [_] "{ecx}" (subleaf),
    );
    return .{ .eax = eax, .ebx = ebx, .ecx = ecx, .edx = edx };
}

/// Reads XCR0. Only meaningful once CPUID has reported OSXSAVE -- executing
/// `xgetbv` without that support faults.
pub fn xgetbv() u32 {
    return asm volatile (
        \\xor %%ecx, %%ecx
        \\xgetbv
        : [_] "={eax}" (-> u32),
        :
        : .{ .edx = true, .ecx = true });
}

/// `cpuid(leaf, sub)`, or zeros when the CPU has no such leaf: below the
/// reference's reading, which takes whatever the highest leaf returns.
fn cpuidChecked(max_leaf: u32, leaf: u32, sub: u32) CpuidResult {
    if (max_leaf < leaf) return .{ .eax = 0, .ebx = 0, .ecx = 0, .edx = 0 };
    return cpuid(leaf, sub);
}

// ---- feature detection -----------------------------------------------------

// The CPUID leaf 0 vendor strings `vendor` distinguishes; tags match the
/// strings.
pub const VendorId = enum { GenuineIntel, AuthenticAMD, CentaurHauls, HygonGenuine, @"  Shanghai  " };

/// The vendor for CPUID leaf 0's `ebx`, `edx`, `ecx` and `eax` (`get_vendor`).
/// A string it does not know is Intel's when the highest leaf is 0 or has
/// bits 8 and 10 set, as the reference reads it.
pub fn vendorOf(ebx: u32, edx: u32, ecx: u32, max_leaf: u32) Vendor {
    var name: [12]u8 = undefined;
    std.mem.writeInt(u32, name[0..4], ebx, .little);
    std.mem.writeInt(u32, name[4..8], edx, .little);
    std.mem.writeInt(u32, name[8..12], ecx, .little);
    const id = std.meta.stringToEnum(VendorId, &name) orelse
        return if (max_leaf == 0 or max_leaf & 0x500 != 0) .intel else .unknown;
    return switch (id) {
        .GenuineIntel => .intel,
        .AuthenticAMD => .amd,
        .CentaurHauls => .centaur,
        .HygonGenuine => .hygon,
        .@"  Shanghai  " => .zhaoxin,
    };
}

pub fn vendor() Vendor {
    const l = cpuid(0, 0);
    return vendorOf(l.ebx, l.edx, l.ecx, l.eax);
}

/// Probes the running CPU: every `cpuid` and `xgetbv` the decisions need.
pub fn probe() Caps {
    var c: Caps = .{};
    const leaf0 = cpuid(0, 0);
    const max = leaf0.eax;
    c.vendor = vendorOf(leaf0.ebx, leaf0.edx, leaf0.ecx, max);

    const leaf1 = cpuidChecked(max, 1, 0);
    c.eax1 = leaf1.eax;
    c.ecx1 = leaf1.ecx;
    c.sse3 = leaf1.ecx & (1 << 0) != 0;
    c.ssse3 = leaf1.ecx & (1 << 9) != 0;
    c.sse4_1 = leaf1.ecx & (1 << 19) != 0;

    // 3DNow! lives in the extended leaves, which only AMD and clones have
    const ext_max = cpuid(0x80000000, 0).eax;
    // a CPU without extended leaves can answer 0x80000000 with a basic leaf
    if (ext_max & 0x80000000 != 0) {
        c.max_ext_leaf = ext_max;
        if (ext_max >= 0x80000001) c.edx_ext1 = cpuid(0x80000001, 0).edx;
    }

    // OS state: XSAVE, OSXSAVE and AVX reported, and XCR0 enabling XMM and YMM
    const state_bits = (1 << 26) | (1 << 27) | (1 << 28);
    const xcr0: u32 = if (leaf1.ecx & state_bits == state_bits) xgetbv() else 0;
    c.avx = xcr0 & 0b110 == 0b110;
    if (!c.avx) return c;
    c.fma3 = leaf1.ecx & (1 << 12) != 0;
    if (c.max_ext_leaf >= 0x80000001) c.fma4 = cpuid(0x80000001, 0).ecx & (1 << 16) != 0;

    const leaf7 = cpuidChecked(max, 7, 0);
    c.avx2 = leaf7.ebx & (1 << 5) != 0;
    // AVX512VL and the OS saving opmask/ZMM state. macOS saves the ZMM state
    // lazily, on the first AVX-512 instruction, so XCR0 does not show it yet.
    const zmm_saved = sysinfo.os.tag.isDarwin() or xcr0 & 0xe0 == 0xe0;
    c.avx512 = c.avx2 and leaf7.ebx & (1 << 31) != 0 and zmm_saved;
    if (!c.avx512) return c;
    c.avx512_bf16 = max >= 7 and cpuid(7, 1).eax & (1 << 5) != 0;

    // CPUID.7.0:EDX AMX-TILE (24) and AMX-BF16 (22); CPUID.D.0:EAX[17:18] state
    const amx_tile = 0x01000000;
    const amx_bf16 = 0x00400000;
    const amx_state = 0x00060000;
    c.amx_bf16 = leaf7.edx & amx_tile != 0 and leaf7.edx & amx_bf16 != 0 and
        cpuidChecked(max, 0xd, 0).eax & amx_state == amx_state;
    return c;
}

/// Whether the host can run what the reference builds for `core`, from `c`:
/// every extension in the core's spec (`_corespec.zig`). The named models there
/// are the AVX-512 cores.
pub fn runnable(core: CoreType, c: Caps) bool {
    switch (core) {
        .unknown => return false,
        // 64-bit builds have no 32-bit-only cores (`gotoblas_dynamic_init`)
        .katmai, .coppermine, .northwood, .banias, .p5, .p6, .viac3, .yonah => if (c.is64) return false,
        // its kernels use 3DNow! (`prefetchw`)
        .athlon => if (c.is64 or !c.has3dnow()) return false,
        else => {},
    }
    const s = corespec.x86(core);
    if (s.model != null) {
        if (!(c.avx512 and c.fma3)) return false;
        if (core == .cooperlake and !c.avx512_bf16) return false;
        if (core == .sapphirerapids and !(c.avx512_bf16 and c.amx_bf16)) return false;
    }
    // the extensions `_corespec.zig` names, by their LLVM feature names
    const Feature = enum { sse3, ssse3, sse4_1, avx, avx2, fma, fma4 };
    for (s.features) |f| {
        const ok = switch (std.meta.stringToEnum(Feature, f) orelse return false) {
            .sse3 => c.sse3,
            .ssse3 => c.ssse3,
            .sse4_1 => c.sse4_1,
            .avx => c.avx,
            .avx2 => c.avx2,
            .fma => c.fma3,
            .fma4 => c.fma4,
        };
        if (!ok) return false;
    }
    return true;
}

/// Whether this machine can run `core`'s kernels (for tests).
pub fn supportsCore(core: CoreType) bool {
    return runnable(core, probe());
}

pub fn features() Features {
    const c = probe();
    return .{
        .vendor = c.vendor,
        .family = c.displayFamily(),
        .model = c.displayModel(),
        .stepping = c.stepping(),
        .avx = c.avx,
        .avx2 = c.avx2,
        .avx512 = c.avx512,
        .avx512_bf16 = c.avx512_bf16,
        .amx_bf16 = c.amx_bf16,
    };
}

// ---- core type selection ---------------------------------------------------

const nehalem_fallback = "Your OS does not support AVX instructions. OpenBLAS is using Nehalem kernels as a fallback, which may give poorer performance.";
const sandybridge_fallback = "Your OS does not support AVX2 instructions. OpenBLAS is using Sandybridge kernels as a fallback, which may give poorer performance.";
const haswell_fallback = "Your OS does not support AVX512VL instructions. OpenBLAS is using Haswell kernels as a fallback, which may give poorer performance.";
const barcelona_fallback = "Your OS does not support AVX instructions. OpenBLAS is using Barcelona kernels as a fallback, which may give poorer performance.";

/// The rungs of the AVX ladder, richest first: where a model's descent starts.
const Rung = enum { sapphire, cooper, skx, haswell, sandy };

/// The descent `dynamic.c` repeats for each model: the richest core of the
/// model's generation that the OS lets the CPU use, then the next one down, with
/// the warning each step has. `avx2_note` is the models whose Haswell step warns.
fn ladder(c: Caps, start: Rung, avx2_note: bool) Decision {
    const s = @intFromEnum(start);
    if (s <= @intFromEnum(Rung.sapphire) and c.amx_bf16) return .{ .core = .sapphirerapids };
    if (s <= @intFromEnum(Rung.cooper) and c.avx512_bf16) return .{ .core = .cooperlake };
    if (s <= @intFromEnum(Rung.skx) and c.avx512) return .{ .core = .skylakex };
    if (s <= @intFromEnum(Rung.haswell) and c.avx2) return .{ .core = .haswell, .note = if (avx2_note) haswell_fallback else null };
    // Sandy Bridge on a Sandy or Ivy Bridge row is no fallback
    if (c.avx) return .{ .core = .sandybridge, .note = if (s < @intFromEnum(Rung.sandy)) sandybridge_fallback else null };
    return .{ .core = .nehalem, .note = nehalem_fallback };
}

/// What an Intel family 6 model is: a core, or a ladder to descend.
const Row = struct {
    exmodel: u8,
    models: []const u8,
    to: union(enum) { core: CoreType, ladder: struct { start: Rung, avx2_note: bool = false } },
};

/// `get_coretype`'s switch for Intel family 6, by extended model.
const intel6 = [_]Row{
    .{ .exmodel = 0, .models = &.{ 0, 1, 2, 3, 4, 5, 6, 7 }, .to = .{ .core = .katmai } },
    .{ .exmodel = 0, .models = &.{ 8, 10, 11 }, .to = .{ .core = .coppermine } },
    .{ .exmodel = 0, .models = &.{ 9, 13, 14 }, .to = .{ .core = .banias } },
    .{ .exmodel = 0, .models = &.{15}, .to = .{ .core = .core2 } },
    .{ .exmodel = 1, .models = &.{6}, .to = .{ .core = .core2 } },
    .{ .exmodel = 1, .models = &.{7}, .to = .{ .core = .penryn } },
    .{ .exmodel = 1, .models = &.{13}, .to = .{ .core = .dunnington } },
    .{ .exmodel = 1, .models = &.{ 10, 11, 14, 15 }, .to = .{ .core = .nehalem } },
    .{ .exmodel = 1, .models = &.{12}, .to = .{ .core = .atom } },
    // Clarkdale/Arrandale, Westmere
    .{ .exmodel = 2, .models = &.{ 5, 12, 14, 15 }, .to = .{ .core = .nehalem } },
    // Sandy Bridge, Sandy Bridge-E
    .{ .exmodel = 2, .models = &.{ 10, 13 }, .to = .{ .ladder = .{ .start = .sandy } } },
    // Ivy Bridge
    .{ .exmodel = 3, .models = &.{ 10, 14 }, .to = .{ .ladder = .{ .start = .sandy } } },
    // Haswell, Broadwell
    .{ .exmodel = 3, .models = &.{ 12, 13, 15 }, .to = .{ .ladder = .{ .start = .haswell } } },
    // Bay Trail
    .{ .exmodel = 3, .models = &.{7}, .to = .{ .core = .atom } },
    // Haswell, Broadwell, Skylake
    .{ .exmodel = 4, .models = &.{ 5, 6, 7, 14, 15 }, .to = .{ .ladder = .{ .start = .haswell } } },
    // Braswell, Avoton
    .{ .exmodel = 4, .models = &.{ 12, 13 }, .to = .{ .core = .nehalem } },
    // Broadwell-DE, Skylake
    .{ .exmodel = 5, .models = &.{ 6, 14 }, .to = .{ .ladder = .{ .start = .haswell } } },
    // Skylake-X, Cascade Lake, Cooper Lake
    .{ .exmodel = 5, .models = &.{5}, .to = .{ .ladder = .{ .start = .cooper, .avx2_note = true } } },
    // Knights Landing
    .{ .exmodel = 5, .models = &.{7}, .to = .{ .ladder = .{ .start = .haswell, .avx2_note = true } } },
    // Apollo Lake, Denverton
    .{ .exmodel = 5, .models = &.{ 12, 15 }, .to = .{ .core = .nehalem } },
    // Cannon Lake
    .{ .exmodel = 6, .models = &.{6}, .to = .{ .ladder = .{ .start = .haswell } } },
    // Ice Lake SP
    .{ .exmodel = 6, .models = &.{ 10, 12 }, .to = .{ .ladder = .{ .start = .cooper } } },
    // Goldmont Plus
    .{ .exmodel = 7, .models = &.{10}, .to = .{ .core = .nehalem } },
    // Ice Lake
    .{ .exmodel = 7, .models = &.{ 13, 14 }, .to = .{ .ladder = .{ .start = .skx, .avx2_note = true } } },
    // Tiger Lake
    .{ .exmodel = 8, .models = &.{ 12, 13 }, .to = .{ .ladder = .{ .start = .skx, .avx2_note = true } } },
    // Kaby Lake, Coffee Lake
    .{ .exmodel = 8, .models = &.{14}, .to = .{ .ladder = .{ .start = .haswell } } },
    // Sapphire Rapids
    .{ .exmodel = 8, .models = &.{15}, .to = .{ .ladder = .{ .start = .sapphire } } },
    // Alder Lake
    .{ .exmodel = 9, .models = &.{ 7, 10 }, .to = .{ .ladder = .{ .start = .cooper } } },
    .{ .exmodel = 9, .models = &.{14}, .to = .{ .ladder = .{ .start = .haswell } } },
    .{ .exmodel = 10, .models = &.{ 5, 6 }, .to = .{ .ladder = .{ .start = .haswell } } },
    .{ .exmodel = 10, .models = &.{7}, .to = .{ .ladder = .{ .start = .skx } } },
};

fn zen(c: Caps, start: Rung) Decision {
    if (start == .cooper and c.avx512_bf16) return .{ .core = .cooperlake };
    if (start != .haswell and c.avx512) return .{ .core = .skylakex };
    if (c.avx) return .{ .core = .zen };
    return .{ .core = .barcelona, .note = barcelona_fallback };
}

/// The AMD (and Hygon) tests, in `get_coretype`'s order. Null where it
/// returns NULL or falls through to the later vendors.
fn amd(c: Caps) ?Decision {
    const family = c.family();
    const exfamily = c.exfamily();
    const model = c.model();
    const exmodel = c.exmodel();
    if (family <= 0xe) {
        // verify 3dnow and 3dnowext before claiming an Athlon
        return if (c.has3dnow()) .{ .core = .athlon } else null;
    }
    if (family != 0xf) return null;
    const need_avx = Decision{ .core = .barcelona, .note = barcelona_fallback };
    switch (exfamily) {
        // K8, family 11h
        0, 2 => {
            if (exmodel == 6 and model == 11) { // QEMU virtual cpu
                if (c.avx512_bf16) return .{ .core = .cooperlake };
                if (c.avx512) return .{ .core = .skylakex };
                return if (c.avx2) .{ .core = .zen } else .{ .core = .barcelona };
            }
            return .{ .core = if (c.ecx1 & 1 != 0) .opteron_sse3 else .opteron };
        },
        // Bobcat (14h), Jaguar (16h)
        5, 7 => return .{ .core = .bobcat },
        // Bulldozer family (15h)
        6 => {
            const core: CoreType = switch (model) {
                1 => .bulldozer,
                2, 3 => .piledriver,
                5 => .excavator,
                0, 8 => switch (exmodel) {
                    1 => .piledriver, // Trinity
                    3 => .steamroller,
                    6 => .excavator,
                    else => return null,
                },
                else => return null,
            };
            return if (c.avx) .{ .core = core } else need_avx;
        },
        // Zen (17h), Hygon Dhyana
        8, 9 => return zen(c, .haswell),
        // Zen 4, 5 (19h, 1Ah)
        10, 11 => return zen(c, .cooper),
        // K10 (10h) and Llano (12h): the reference returns NULL and lands on
        // Prescott; Barcelona is their core
        1, 3 => return .{ .core = .barcelona },
        else => return null,
    }
}

fn centaur(c: Caps) ?Decision {
    switch (c.family()) {
        0x6 => {
            if (c.model() == 0xf and c.stepping() < 0xe) return .{ .core = .nano };
            return .{ .core = .nehalem };
        },
        0x7 => switch (c.exmodel()) {
            5, 6 => return .{ .core = if (c.avx2) .zen else .dunnington },
            else => return .{ .core = .nehalem },
        },
        else => return if (c.family() >= 0x8) .{ .core = .nehalem } else null,
    }
}

fn zhaoxin(c: Caps) Decision {
    if (c.avx2) return .{ .core = .zen };
    if (c.family() == 0x7 and c.exmodel() == 5) return .{ .core = .dunnington };
    return .{ .core = .nehalem };
}

/// `get_coretype`: the model tables. Null where the reference returns NULL, so
/// the capability ladder decides.
pub fn modelCore(c: Caps) ?Decision {
    switch (c.vendor) {
        .intel => switch (c.family()) {
            0x6 => for (intel6) |r| {
                if (r.exmodel != c.exmodel()) continue;
                if (std.mem.indexOfScalar(u8, r.models, @intCast(c.model())) == null) continue;
                return switch (r.to) {
                    .core => |k| .{ .core = k },
                    .ladder => |l| ladder(c, l.start, l.avx2_note),
                };
            },
            0xf => return .{ .core = if (c.model() <= 0x2) .northwood else .prescott },
            else => {},
        },
        .amd, .hygon => if (amd(c)) |d| return d,
        .centaur => if (centaur(c)) |d| return d,
        .zhaoxin => return zhaoxin(c),
        .unknown => {},
    }
    return null;
}

/// The capability ladder `gotoblas_dynamic_init` falls back on for a CPU the
/// model tables do not know: the richest core the CPU and OS support.
fn capsCore(c: Caps, is64: bool) CoreType {
    if (c.avx512_bf16) return .cooperlake;
    if (c.avx512) return .skylakex;
    if (c.avx2) return .haswell;
    if (c.avx) return .sandybridge;
    // 32-bit builds end at Katmai, the reference's pick; 64-bit at Prescott. A
    // CPU with no SSE3 gets the most that it can run.
    if (c.sse3) return .prescott;
    return if (is64) .opteron else .katmai;
}

/// The core for `c`: the model tables, else the capability ladder, with the
/// sanity step that 64-bit builds have no 32-bit-only cores. Never a core `c`
/// cannot run: where the pick needs more than `c` has (a 64-bit pick of Prescott
/// on a CPU without SSE3, FMA missing under an AVX2 pick, 3DNow! cores), the
/// closest core by ISA stands in.
pub fn decideWidth(c: Caps, is64: bool) Decision {
    var caps = c;
    caps.is64 = is64;
    var d = modelCore(caps) orelse Decision{ .core = capsCore(caps, is64) };
    if (is64) switch (d.core) {
        .katmai, .coppermine, .northwood, .banias, .athlon => d = .{ .core = .prescott },
        else => {},
    };
    return select.guard(@This(), d, caps);
}

/// `decideWidth` for the word size `c` says.
pub fn decide(c: Caps) Decision {
    return decideWidth(c, c.is64);
}

/// Picks a core type for the running CPU (see `decideWidth`).
pub fn detect() CoreType {
    return detectNoted().core;
}

/// `detect`, with the warning the reference gives when a fallback was taken.
pub fn detectNoted() Decision {
    return decide(probe());
}

/// The coarse tier, for callers that only care about vector width and FMA
/// rather than the specific microarchitecture.
pub fn detectTier() Tier {
    const c = probe();
    if (c.avx512) return .avx512;
    if (c.avx2) return .avx2;
    return .baseline;
}

/// Module-level so every kernel and element type shares one probe. Racing
/// writers store the same value, so `.monotonic` suffices.
var cached_tier = std.atomic.Value(Tier).init(._unresolved);

/// `detectTier`, probed once and cached; use this on hot paths (cpuid is slow).
pub fn cachedTier() Tier {
    const c = cached_tier.load(.monotonic);
    if (c != ._unresolved) return c;
    const t = detectTier();
    cached_tier.store(t, .monotonic);
    return t;
}

// ---- tests -----------------------------------------------------------------

/// A CPU as the tests describe it: CPUID.1:EAX and the extensions the OS lets it use.
fn cpu(v: Vendor, eax1: u32, level: enum { sse3, sse41, avx, avx2, avx512, bf16, amx }) Caps {
    var c: Caps = .{ .vendor = v, .eax1 = eax1 };
    const n = @intFromEnum(level);
    c.sse3 = true;
    c.ssse3 = true;
    c.sse4_1 = n >= 1;
    c.avx = n >= 2;
    // the Bulldozer family has FMA3 and FMA4 below AVX2
    c.fma3 = n >= 3 or (v == .amd and n == 2);
    c.fma4 = v == .amd and n == 2;
    c.avx2 = n >= 3;
    c.avx512 = n >= 4;
    c.avx512_bf16 = n >= 5;
    c.amx_bf16 = n >= 6;
    c.max_ext_leaf = 0x80000008;
    return c;
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
        c.sse3 = false;
        c.ecx1 = 0;
        c.ssse3 = false;
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
    try expectCore(.excavator, cpu(.amd, 0x660f51, .avx)); // Carrizo (model 5)
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
    c.sse3 = false;
    c.ssse3 = false;
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
    try expectCore(.sandybridge, cpu(.unknown, 0, .avx));
    try expectCore(.prescott, cpu(.unknown, 0, .sse3));
    // no SSE3 at all: a 64-bit build still gets a core it can run
    var old = cpu(.unknown, 0, .sse3);
    old.sse3 = false;
    old.ssse3 = false;
    try expectCore(.opteron, old);
    try std.testing.expectEqual(CoreType.katmai, decideWidth(old, false).core);
}

test "FMA3 missing under AVX2 does not give the Haswell kernels" {
    var c = cpu(.intel, 0x306c3, .avx2);
    c.fma3 = false;
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
    // unknown strings: Intel when the highest leaf is 0 or has bits 8 and 10
    try std.testing.expectEqual(Vendor.intel, w.v("TransmetaCPU", 0));
    try std.testing.expectEqual(Vendor.intel, w.v("TransmetaCPU", 0x500));
    try std.testing.expectEqual(Vendor.unknown, w.v("TransmetaCPU", 0x10));
}

test "display family and model add the extended fields" {
    const skx = cpu(.intel, 0x50654, .avx512);
    try std.testing.expectEqual(@as(u32, 6), skx.displayFamily());
    try std.testing.expectEqual(@as(u32, 0x55), skx.displayModel());
    const zen2 = cpu(.amd, 0x870f10, .avx2);
    try std.testing.expectEqual(@as(u32, 0x17), zen2.displayFamily());
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
    if (comptime !sysinfo.cpu.arch.isX86()) return error.SkipZigTest;
    const f = features();
    // Nothing to assert about which CPU this is -- the point is that every
    // cpuid/xgetbv path runs to completion and reports a coherent result.
    try std.testing.expect(f.family > 0);
    if (f.avx512) try std.testing.expect(f.avx2);
    if (f.avx2) try std.testing.expect(f.avx);
}

test "cachedTier matches detectTier and is stable across calls" {
    if (comptime !sysinfo.cpu.arch.isX86()) return error.SkipZigTest;
    try std.testing.expectEqual(detectTier(), cachedTier());
    try std.testing.expectEqual(cachedTier(), cachedTier());
}

test "detect returns a concrete core type the host can run" {
    if (comptime !sysinfo.cpu.arch.isX86()) return error.SkipZigTest;
    try std.testing.expect(detect() != .unknown);
    try std.testing.expect(supportsCore(detect()));
    switch (detectTier()) {
        ._unresolved => unreachable,
        .baseline, .avx2, .avx512 => {},
    }
}
