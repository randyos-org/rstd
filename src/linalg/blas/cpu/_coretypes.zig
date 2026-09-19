//! The core type enums of every architecture, in one import-free file so that
//! code which only needs to name a core (the per-core kernel objects, which are
//! compiled on their own) does not pull in the detection code. Each `<arch>.zig`
//! re-exports its own as `CoreType`.

// ---- x86.zig ----
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
};

// ---- arm64.zig ----
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
};

// ---- arm.zig ----
pub const ArmCoreType = enum {
    armv5,
    armv6,
    armv7,
    neon,
    vfpv4,
};

// ---- power.zig ----
/// Mirrors the reference's `corename[]` ordering.
pub const PowerCoreType = enum {
    unknown,
    power6,
    power8,
    power9,
    power10,
};

// ---- riscv.zig ----
/// Vector kernels are selected by guaranteed vector register length
/// (`VLEN`), not by microarchitecture.
pub const RiscvCoreType = enum {
    riscv64_generic,
    riscv64_zvl128b,
    riscv64_zvl256b,
};

// ---- s390x.zig ----
pub const S390xCoreType = enum {
    zarch_generic,
    z13,
    z14,
};

// ---- mips.zig ----
/// Mirrors the reference's `corename[]` ordering.
pub const MipsCoreType = enum {
    mips64_generic,
    loongson3r3,
    loongson3r4,
    unknown,
};

// ---- loongarch.zig ----
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

// ---- generic.zig ----
pub const GenericCoreType = enum { generic };
