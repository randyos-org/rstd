//! The core type enums of every architecture, in one file free of rstd
//! imports so that code which only needs to name a core (the per-core kernel
//! objects, which are compiled on their own) does not pull in the detection
//! code. Each `<arch>.zig` re-exports its own as `CoreType`.

const std = @import("std");

// ---- x86 ----
/// Core types the reference can select between, in the order of its
/// `corename[]` table (`driver/others/dynamic.c`). Names are matched
/// case-insensitively against `OPENBLAS_CORETYPE`.
pub const X86CoreType = enum {
    unknown,
    katmai,
    coppermine,
    northwood,
    prescott,
    banias,
    atom,
    core2,
    penryn,
    dunnington,
    nehalem,
    athlon,
    opteron,
    opteron_sse3,
    barcelona,
    nano,
    sandybridge,
    bobcat,
    bulldozer,
    piledriver,
    haswell,
    steamroller,
    excavator,
    zen,
    skylakex,
    cooperlake,
    sapphirerapids,
    // Not detected at run time (`dynamic.c` has no such core); the targets of
    // static builds, named as in the reference's `TARGET`.
    p5,
    p6,
    viac3,
    yonah,
};

// ---- arm64 ----
/// Mirrors the reference's `corename[]` ordering.
pub const Arm64CoreType = enum {
    armv8,
    cortexa53,
    cortexa57,
    cortexa72,
    cortexa73,
    falkor,
    thunderx,
    thunderx2t99,
    tsv110,
    emag8180,
    neoversen1,
    neoversev1,
    neoversev2,
    neoversen2,
    thunderx3t110,
    cortexa55,
    armv8sve,
    a64fx,
    armv9sme,
    vortex,
    vortexm4,
    unknown,
    // Not detected at run time (`dynamic_arm64.c` has no such core); the
    // targets of static builds, named as in the reference's `TARGET`.
    cortexa76,
    cortexa510,
    cortexa710,
    cortexx1,
    cortexx2,
    ampere1,
    ft2000,
};

// ---- arm ----
pub const ArmCoreType = enum {
    armv5,
    armv6,
    armv7,
    neon,
    vfpv4,
    // static targets of the reference, which `cpuid_arm.c` finds by part
    cortexa9,
    cortexa15,
};

// ---- power ----
/// Mirrors the reference's `corename[]` ordering.
pub const PowerCoreType = enum {
    unknown,
    power6,
    power8,
    power9,
    power10,
    // Not detected at run time (`dynamic_power.c` has no such core); the
    // targets of static builds, named as in the reference's `TARGET`.
    power3,
    power4,
    power5,
    ppc970,
    cell,
    ppc440,
    ppc440fp2,
    ppcg4,
};

// ---- riscv ----
/// Vector kernels are selected by guaranteed vector register length
/// (`VLEN`), not by microarchitecture.
pub const RiscvCoreType = enum {
    riscv64_generic,
    riscv64_zvl128b,
    riscv64_zvl256b,
    // Not detected at run time (`dynamic_riscv64.c` has no such core); the
    // targets of static builds, named as in the reference's `TARGET`.
    c910v,
    u74,
    x280,
};

// ---- s390x ----
pub const S390xCoreType = enum {
    zarch_generic,
    z13,
    z14,
};

// ---- mips ----
/// Mirrors the reference's `corename[]` ordering.
pub const MipsCoreType = enum {
    mips64_generic,
    loongson3r3,
    loongson3r4,
    unknown,
    // Not detected at run time (`dynamic_mips64.c` has no such core); the
    // targets of static builds, named as in the reference's `TARGET`. The
    // first three are 32-bit, the last three MIPS64 release 6.
    mips24k,
    mips1004k,
    p5600,
    i6400,
    i6500,
    p6600,
};

// ---- loongarch ----
/// Mirrors the reference's `corename[]` ordering.
pub const LoongarchCoreType = enum {
    la64_generic,
    la264,
    la464,
    loongsongeneric,
    loongson2k1000,
    loongson3r5,
    unknown,
};

// ---- sparc ----
/// The reference's two SPARC targets, named as in its `TARGET`: `SPARC` is
/// V9 and `SPARCV7` the generic one; neither is detected at run time.
pub const SparcCoreType = enum {
    sparcv7,
    sparc,
};

// ---- wasm ----
/// The reference has one target for WebAssembly, `WASM128_GENERIC`; the other
/// core type is a module compiled without SIMD128.
pub const WasmCoreType = enum {
    wasm_generic,
    wasm128_generic,
};

// ---- generic ----
pub const GenericCoreType = enum { generic };

/// The cores of `Core` no detector returns: the reference's static-only
/// `TARGET`s, which only `OPENBLAS_CORETYPE` names (`select.staticOnly`).
pub fn staticOnly(comptime Core: type) []const Core {
    return switch (Core) {
        X86CoreType => &.{ .p5, .p6, .viac3, .yonah },
        Arm64CoreType => &.{ .cortexa76, .cortexa510, .cortexa710, .cortexx1, .cortexx2, .ampere1, .ft2000 },
        RiscvCoreType => &.{ .c910v, .u74, .x280 },
        MipsCoreType => &.{ .mips24k, .mips1004k, .p5600, .i6400, .i6500, .p6600 },
        else => &.{},
    };
}

/// Every family's core type.
pub const all = .{
    X86CoreType,   Arm64CoreType, ArmCoreType,       PowerCoreType,
    RiscvCoreType, MipsCoreType,  LoongarchCoreType, S390xCoreType,
    SparcCoreType, WasmCoreType,  GenericCoreType,
};

/// The architecture families: one detector and kernel registry each, but
/// x86 and x86_64 share a detector and core type.
pub const Family = enum { x86_64, x86, arm64, arm, power, riscv, s390x, mips, loongarch, sparc, wasm, generic };

/// The family of `arch`.
pub fn familyOf(arch: std.Target.Cpu.Arch) Family {
    return switch (arch) {
        .x86_64 => .x86_64,
        .x86 => .x86,
        .aarch64, .aarch64_be => .arm64,
        .arm, .armeb, .thumb, .thumbeb => .arm,
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => .power,
        .riscv32, .riscv32be, .riscv64, .riscv64be => .riscv,
        .s390x => .s390x,
        .mips, .mipsel, .mips64, .mips64el => .mips,
        .loongarch32, .loongarch64 => .loongarch,
        .sparc, .sparc64 => .sparc,
        .wasm32, .wasm64 => .wasm,
        else => .generic,
    };
}

/// The core type of `family`.
pub fn CoreType(comptime family: Family) type {
    return switch (family) {
        .x86_64, .x86 => X86CoreType,
        .arm64 => Arm64CoreType,
        .arm => ArmCoreType,
        .power => PowerCoreType,
        .riscv => RiscvCoreType,
        .s390x => S390xCoreType,
        .mips => MipsCoreType,
        .loongarch => LoongarchCoreType,
        .sparc => SparcCoreType,
        .wasm => WasmCoreType,
        .generic => GenericCoreType,
    };
}
