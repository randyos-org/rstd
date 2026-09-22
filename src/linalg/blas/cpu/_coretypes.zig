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

// ---- riscv.zig ----
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
// ---- sparc.zig ----
/// The reference's two SPARC targets, named as in its `TARGET`: `SPARC` is
/// V9 and `SPARCV7` the generic one; neither is detected at run time.
pub const SparcCoreType = enum {
    sparcv7,
    sparc,
};

// ---- wasm.zig ----
/// The reference has one target for WebAssembly, `WASM128_GENERIC`; the other
/// core type is a module compiled without SIMD128.
pub const WasmCoreType = enum {
    wasm_generic,
    wasm128_generic,
};

pub const GenericCoreType = enum { generic };
