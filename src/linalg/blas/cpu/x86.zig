//! x86 / x86_64 CPU detection, mirroring `driver/others/dynamic.c`.
//!
//! This is the only architecture that can detect itself with no help from
//! the OS: `cpuid` and `xgetbv` are unprivileged, so everything here works
//! identically on Windows, Linux, macOS and freestanding.

const std = @import("std");
const coretypes = @import("_coretypes.zig");
const corespec = @import("_corespec.zig");
const zigstd = @import("../../../zigstd.zig");
const sysinfo = @import("builtin");

pub const CoreType = coretypes.X86CoreType;
pub const kernel = @import("../kernels/x86.zig");

pub const Vendor = enum { unknown, intel, amd, centaur, hygon, zhaoxin };

/// The coarse SIMD capability tiers that actually change which kernel we
/// want, as opposed to the fine-grained microarchitecture identity above.
/// The reference distinguishes both: two chips can share a tier (Zen and
/// Haswell are both AVX2) and still want different cache blocking.
/// `_unresolved` is only `cachedTier`'s pre-probe state; never returned.
pub const Tier = enum(u8) { _unresolved, baseline, avx2, avx512 };

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

// ---- feature detection -----------------------------------------------------

/// Mirrors `support_avx()`: the CPU must report XSAVE/OSXSAVE/AVX *and* the
/// OS must have enabled saving the AVX register state. CPUID alone is not
/// enough -- a hypervisor or an older OS can leave the state disabled, and
/// using AVX anyway faults.
pub fn supportAvx() bool {
    const leaf1 = cpuid(1, 0);
    const bits = (1 << 26) | (1 << 27) | (1 << 28); // XSAVE | OSXSAVE | AVX
    if (leaf1.ecx & bits != bits) return false;
    return xgetbv() & 0b110 == 0b110;
}

/// Mirrors `support_avx2()`: AVX plus CPUID.7.0:EBX[5].
pub fn supportAvx2() bool {
    if (!supportAvx()) return false;
    return cpuid(7, 0).ebx & (1 << 5) != 0;
}

/// FMA3: CPUID.1:ECX[12], on top of `supportAvx` (which checks the OS state).
pub fn supportFma3() bool {
    if (!supportAvx()) return false;
    return cpuid(1, 0).ecx & (1 << 12) != 0;
}

/// FMA4 (AMD Bulldozer only): CPUID.80000001h:ECX[16], on top of `supportAvx`.
pub fn supportFma4() bool {
    if (!supportAvx()) return false;
    if (cpuid(0x80000000, 0).eax < 0x80000001) return false;
    return cpuid(0x80000001, 0).ecx & (1 << 16) != 0;
}

/// Mirrors `support_avx512()`: AVX512VL reported, plus the OS saving the
/// opmask/ZMM_Hi256/Hi16_ZMM state (XCR0 bits 5-7).
pub fn supportAvx512() bool {
    if (!supportAvx()) return false;
    const leaf7 = cpuid(7, 0);
    if (leaf7.ebx & (1 << 5) == 0) return false; // no AVX2 => no AVX-512 either
    if (leaf7.ebx & (1 << 31) == 0) return false; // AVX512VL
    return xgetbv() & 0xe0 == 0xe0;
}

/// SSE3, SSSE3 and SSE4.1: CPUID.1:ECX bits 0, 9 and 19.
fn supportSse(comptime bit: u5) bool {
    return cpuid(1, 0).ecx & (1 << bit) != 0;
}

/// Whether the host can run what the reference builds for `core`: every
/// extension in the core's spec (`_corespec.zig`). The named models there are
/// the AVX-512 cores.
pub fn supportsCore(core: CoreType) bool {
    const s = corespec.x86(core);
    if (s.model != null and !(supportAvx512() and supportFma3())) return false;
    for (s.features) |f| {
        const ok = if (std.mem.eql(u8, f, "sse3")) supportSse(0) else if (std.mem.eql(u8, f, "ssse3")) supportSse(9) else if (std.mem.eql(u8, f, "sse4_1")) supportSse(19) else if (std.mem.eql(u8, f, "avx")) supportAvx() else if (std.mem.eql(u8, f, "avx2")) supportAvx2() else if (std.mem.eql(u8, f, "fma")) supportFma3() else if (std.mem.eql(u8, f, "fma4")) supportFma4() else false;
        if (!ok) return false;
    }
    return true;
}

/// Mirrors `support_avx512_bf16()`: CPUID.7.1:EAX[5].
pub fn supportAvx512Bf16() bool {
    if (!supportAvx512()) return false;
    return cpuid(7, 1).eax & (1 << 5) != 0;
}

/// Mirrors `support_amx_bf16()`: AMX-TILE and AMX-BF16 reported in
/// CPUID.7.0:EDX, and the AMX state enabled in CPUID.D.0:EAX[17:18].
pub fn supportAmxBf16() bool {
    if (!supportAvx512()) return false;
    const leaf7 = cpuid(7, 0);
    const amx_tile = 0x01000000;
    const amx_bf16 = 0x00400000;
    if (leaf7.edx & amx_tile == 0 or leaf7.edx & amx_bf16 == 0) return false;
    const amx_enabled = 0x00060000;
    return cpuid(0xd, 0).eax & amx_enabled == amx_enabled;
}

pub fn vendor() Vendor {
    const leaf0 = cpuid(0, 0);
    var name: [12]u8 = undefined;
    std.mem.writeInt(u32, name[0..4], leaf0.ebx, .little);
    std.mem.writeInt(u32, name[4..8], leaf0.edx, .little);
    std.mem.writeInt(u32, name[8..12], leaf0.ecx, .little);

    if (std.mem.eql(u8, &name, "GenuineIntel")) return .intel;
    if (std.mem.eql(u8, &name, "AuthenticAMD")) return .amd;
    if (std.mem.eql(u8, &name, "CentaurHauls")) return .centaur;
    if (std.mem.eql(u8, &name, "HygonGenuine")) return .hygon;
    if (std.mem.eql(u8, &name, "  Shanghai  ")) return .zhaoxin;
    return .unknown;
}

pub fn features() Features {
    var f: Features = .{ .vendor = vendor() };

    const leaf1 = cpuid(1, 0);
    // Family/model use the extended-field encoding: family 0xf adds the
    // extended family, and families 6 and 0xf add the extended model.
    var family = (leaf1.eax >> 8) & 0xf;
    var model = (leaf1.eax >> 4) & 0xf;
    if (family == 0xf) family += (leaf1.eax >> 20) & 0xff;
    if (family == 6 or family == 0xf) model += ((leaf1.eax >> 16) & 0xf) << 4;

    f.family = family;
    f.model = model;
    f.stepping = leaf1.eax & 0xf;

    f.avx = supportAvx();
    f.avx2 = supportAvx2();
    f.avx512 = supportAvx512();
    f.avx512_bf16 = supportAvx512Bf16();
    f.amx_bf16 = supportAmxBf16();
    return f;
}

// ---- core type selection ---------------------------------------------------

/// Picks a core type for the running CPU.
///
/// The reference does this in two stages: a large vendor + family/model
/// table (`get_coretype`), then a capability-based fallback in
/// `gotoblas_dynamic_init` for anything the table does not recognize. The
/// fallback is what actually runs on any CPU newer than the table, so it
/// is the part that matters most; the model table below covers the modern
/// entries and is meant to grow rather than ever be exhaustive.
pub fn detect() CoreType {
    const f = features();

    switch (f.vendor) {
        .intel => {
            if (f.family == 6) {
                switch (f.model) {
                    // Sapphire Rapids / Emerald Rapids
                    0x8f, 0xcf => if (f.avx512) return .sapphirerapids,
                    // Cooper Lake
                    0x55 => if (f.avx512_bf16) return .cooperlake else if (f.avx512) return .skylakex,
                    // Skylake-X / Cascade Lake / Ice Lake family
                    0x6a, 0x6c, 0x7d, 0x7e, 0x9d => if (f.avx512) return .skylakex,
                    // Haswell / Broadwell / client Skylake and successors
                    0x3c, 0x3f, 0x45, 0x46, 0x3d, 0x47, 0x4f, 0x56 => return .haswell,
                    0x4e, 0x5e, 0x8e, 0x9e, 0xa5, 0xa6, 0xa7 => return .haswell,
                    // Sandy Bridge / Ivy Bridge
                    0x2a, 0x2d, 0x3a, 0x3e => return .sandybridge,
                    // Nehalem / Westmere
                    0x1a, 0x1e, 0x1f, 0x2e, 0x25, 0x2c, 0x2f => return .nehalem,
                    // Penryn / Core2
                    0x17, 0x1d => return .penryn,
                    0x0f, 0x16 => return .core2,
                    // Atom
                    0x1c, 0x26, 0x27, 0x35, 0x36 => return .atom,
                    else => {},
                }
            }
        },
        .amd => {
            // Zen and later all report family >= 0x17.
            if (f.family >= 0x17) return .zen;
            switch (f.family) {
                0x15 => return .bulldozer, // Bulldozer/Piledriver/Steamroller/Excavator
                0x14, 0x16 => return .bobcat,
                0x10, 0x11, 0x12 => return .barcelona,
                0x0f => return .opteron_sse3,
                0x06 => return .athlon,
                else => {},
            }
        },
        .centaur, .zhaoxin => return .nano,
        .hygon => return .zen,
        .unknown => {},
    }

    // Capability-based fallback, matching `gotoblas_dynamic_init`.
    if (f.avx512_bf16) return .cooperlake;
    if (f.avx512) return .skylakex;
    if (f.avx2) return .haswell;
    if (f.avx) return .sandybridge;
    return .prescott;
}

/// The coarse tier, for callers that only care about vector width and FMA
/// rather than the specific microarchitecture.
pub fn detectTier() Tier {
    if (supportAvx512()) return .avx512;
    if (supportAvx2()) return .avx2;
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

const Vtable = @import("../Vtable.zig");
const axpy = @import("../blas1/axpy/__root__.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");
pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    axpy.install(&vt, core_type);
    gemm_params.installCore(&vt, (if (@sizeOf(usize) == 8) .x86_64 else .x86), core_type, .fromConfig(blas.config));
    return vt;
}

test "feature probes execute without faulting" {
    if (comptime !sysinfo.cpu.arch.isX86()) return error.SkipZigTest;
    const f = features();
    // Nothing to assert about which CPU this is -- the point is that every
    // cpuid/xgetbv path runs to completion and reports a coherent result.
    try zigstd.testExpectTrue(f.family > 0);
    if (f.avx512) try zigstd.testExpectTrue(f.avx2);
    if (f.avx2) try zigstd.testExpectTrue(f.avx);
}

test "cachedTier matches detectTier and is stable across calls" {
    if (comptime !sysinfo.cpu.arch.isX86()) return error.SkipZigTest;
    try zigstd.testExpectEqual(detectTier(), cachedTier());
    try zigstd.testExpectEqual(cachedTier(), cachedTier());
}

test "detect returns a concrete core type and a matching tier" {
    if (comptime !sysinfo.cpu.arch.isX86()) return error.SkipZigTest;
    try zigstd.testExpectTrue(detect() != .unknown);
    switch (detectTier()) {
        ._unresolved => unreachable,
        .baseline, .avx2, .avx512 => {},
    }
}
