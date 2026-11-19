//! What the CPUs that run each core have (`corespec.zig`), from the hardware.

const coretypes = @import("coretypes.zig");
const corespec = @import("corespec.zig");

const Spec = corespec.Spec;
const no_altivec = corespec.no_altivec;
const mips_simd = corespec.mips_simd;

/// RV64GC as features: Zig 0.17 has no such model.
const rv64gc = [_][]const u8{ "m", "a", "f", "d", "c", "zicsr", "zifencei" };

/// What every CPU of `core`'s kind has: the LLVM model of the CPU it is named
/// for, or of the weakest part it covers (else the generic one), plus what
/// `runnable` requires. The detectors' last resorts can still pick a core for
/// a weaker CPU (`katmai` for any 32-bit x86 without SSE3, as the reference
/// does; design.md, Core detection and selection); that is theirs, not this.
pub fn machine(core: anytype) Spec {
    return switch (@TypeOf(core)) {
        coretypes.X86CoreType => switch (core) {
            .unknown => .{},
            .p5 => .{ .model = "pentium" },
            .p6 => .{ .model = "pentiumpro" },
            .katmai, .coppermine => .{ .model = "pentium3" },
            .northwood => .{ .model = "pentium4" },
            .banias => .{ .model = "pentium_m" },
            .yonah => .{ .model = "yonah" },
            .viac3 => .{ .model = "c3" },
            .athlon => .{ .model = "athlon" },
            .opteron => .{ .model = "k8" },
            .opteron_sse3 => .{ .model = "k8_sse3" },
            // also unknown CPUs with SSE3, and AMD's without AVX (QEMU's)
            .prescott, .barcelona => .{ .model = "prescott" },
            // VIA's Nano has no model of its own
            .core2, .nano => .{ .model = "core2" },
            .atom => .{ .model = "bonnell" },
            // Nehalem's: also VIA's and Zhaoxin's with SSE4.1 alone
            .penryn, .dunnington, .nehalem => .{ .model = "penryn" },
            .sandybridge => .{ .model = "sandybridge" },
            .bobcat => .{ .model = "btver1" },
            .bulldozer => .{ .model = "bdver1" },
            .piledriver => .{ .model = "bdver2" },
            .steamroller => .{ .model = "bdver3" },
            .excavator => .{ .model = "bdver4" },
            .haswell => .{ .model = "haswell" },
            // also Zhaoxin's and VIA's AVX2 parts: none of AMD's own extensions
            .zen => .{ .model = "znver1", .without = &.{ "sse4a", "clzero", "mwaitx" } },
            .skylakex => .{ .model = "skylake_avx512" },
            .cooperlake => .{ .model = "cooperlake" },
            .sapphirerapids => .{ .model = "sapphirerapids" },
        },
        coretypes.Arm64CoreType => switch (core) {
            .armv8, .ft2000, .unknown => .{ .features = &.{"neon"} },
            .cortexa53 => .{ .model = "cortex_a53" },
            .cortexa55 => .{ .model = "cortex_a55" },
            .cortexa57 => .{ .model = "cortex_a57" },
            .cortexa72 => .{ .model = "cortex_a72" },
            // also the A75 and A76 parts
            .cortexa73 => .{ .model = "cortex_a73" },
            .cortexa76 => .{ .model = "cortex_a76" },
            .cortexa510 => .{ .model = "cortex_a510" },
            .cortexa710 => .{ .model = "cortex_a710" },
            .cortexx1 => .{ .model = "cortex_x1" },
            .cortexx2 => .{ .model = "cortex_x2" },
            .falkor => .{ .model = "falkor" },
            .thunderx => .{ .model = "thunderx" },
            .thunderx2t99 => .{ .model = "thunderx2t99" },
            .thunderx3t110 => .{ .model = "thunderx3t110" },
            .tsv110 => .{ .model = "tsv110" },
            .emag8180 => .{ .model = "emag" },
            // also Ampere's first parts and the SVE cores without SVE
            .neoversen1 => .{ .model = "neoverse_n1" },
            .neoversev1 => .{ .model = "neoverse_v1" },
            .neoversen2 => .{ .model = "neoverse_n2" },
            .neoversev2 => .{ .model = "neoverse_v2" },
            .a64fx => .{ .model = "a64fx" },
            .ampere1 => .{ .model = "ampere1" },
            // SVE; SVE2, SME and SME's optional F64F64 (`runnable`)
            .armv8sve => .{ .features = &.{ "neon", "sve" } },
            .armv9sme => .{ .features = &.{ "v9a", "sme", "sme_f64f64" } },
            // Apple's: macOS and Linux run on the M1 and later
            .vortex => .{ .model = "apple_m1" },
            .vortexm4 => .{ .model = "apple_m4" },
        },
        // `runnable`: VFP; VFPv3 with 32 registers; NEON (with VFPv4 on the A15)
        coretypes.ArmCoreType => switch (core) {
            .armv5 => .{ .features = &.{"v5te"} },
            .armv6 => .{ .features = &.{ "v6", "vfp2" } },
            .armv7 => .{ .features = &.{ "v7a", "vfp3" } },
            .neon => .{ .features = &.{ "v7a", "neon" } },
            .vfpv4 => .{ .features = &.{ "v7a", "vfp4" } },
            .cortexa9 => .{ .model = "cortex_a9", .features = &.{"neon"} },
            .cortexa15 => .{ .model = "cortex_a15", .features = &.{ "neon", "vfp4" } },
        },
        // LLVM's POWER3 to POWER5 have AltiVec, which the CPUs do not
        coretypes.PowerCoreType => switch (core) {
            .unknown => .{},
            .power3 => .{ .model = "pwr3", .without = &no_altivec },
            .power4 => .{ .model = "pwr4", .without = &no_altivec },
            .power5 => .{ .model = "pwr5", .without = &no_altivec },
            .power6 => .{ .model = "pwr6" },
            .power8 => .{ .model = "pwr8" },
            .power9 => .{ .model = "pwr9" },
            .power10 => .{ .model = "pwr10" },
            .ppc970, .cell => .{ .model = "970" },
            .ppc440, .ppc440fp2 => .{ .model = "440" },
            .ppcg4 => .{ .model = "g4" },
        },
        // `runnable`: by guaranteed vector length only
        coretypes.RiscvCoreType => switch (core) {
            .riscv64_generic, .c910v, .u74 => .{ .features = &rv64gc },
            .riscv64_zvl128b => .{ .features = &(rv64gc ++ [_][]const u8{"v"}) },
            .riscv64_zvl256b => .{ .features = &(rv64gc ++ [_][]const u8{ "v", "zvl256b" }) },
            .x280 => .{ .features = &(rv64gc ++ [_][]const u8{ "v", "zvl512b" }) },
        },
        coretypes.MipsCoreType => switch (core) {
            .mips64_generic, .unknown => .{ .model = "mips64" },
            .loongson3r3 => .{ .model = "mips64r2" },
            .loongson3r4 => .{ .model = "mips64r5", .features = &.{"msa"} },
            .mips24k, .mips1004k => .{ .model = "mips32r2" },
            .p5600 => .{ .model = "p5600", .features = &mips_simd },
            .i6400 => .{ .model = "i6400", .features = &mips_simd },
            .i6500 => .{ .model = "i6500", .features = &mips_simd },
            .p6600 => .{ .model = "mips64r6", .features = &.{"msa"} },
        },
        coretypes.LoongarchCoreType => switch (core) {
            .la64_generic, .loongsongeneric, .unknown => .{ .model = "la64v1_0" },
            .la264, .loongson2k1000 => .{ .model = "la64v1_0", .features = &.{"lsx"} },
            .la464, .loongson3r5 => .{ .model = "la464" },
        },
        coretypes.S390xCoreType => switch (core) {
            .zarch_generic => .{ .model = "arch8" },
            .z13 => .{ .model = "z13" },
            .z14 => .{ .model = "z14" },
        },
        coretypes.SparcCoreType => switch (core) {
            .sparcv7 => .{ .model = "v7" },
            .sparc => .{ .model = "v9" },
        },
        coretypes.WasmCoreType => switch (core) {
            .wasm_generic => .{ .model = "mvp" },
            .wasm128_generic => .{ .model = "mvp", .features = &.{"simd128"} },
        },
        else => .{},
    };
}
