//! What each core type is compiled with, taken from how OpenBLAS builds its
//! per-core kernels (`DYNAMIC_ARCH`): the ISA extensions the core's flags
//! enable (`Makefile.<arch>`, the `HAVE_*` tables in `getarch.c`) and the
//! explicit overrides in `kernel/Makefile`.
//!
//! This is the single description of that. The build (`build/blaskern.zig`)
//! reads it to give each kernel object its CPU features, grouping cores whose
//! specs are equal into one object, and the runtime reads the same data, for
//! instance to check that the host can run the kernels a core type selects.
//! It is plain data, so both can use it.
//!
//! Names are those of `std.Target.<arch>`'s models and features; a test checks
//! every one resolves. Only the ISA is described: the reference also passes
//! `-mtune=`, which only affects scheduling and has no counterpart here.
//! The non-x86 entries follow the primary (current-compiler) branch of each
//! Makefile and have not been exercised yet; no kernels are built for those
//! architectures so far.

const std = @import("std");
const coretypes = @import("_coretypes.zig");

pub const Spec = struct {
    /// A named model where the reference passes `-march=`/`-mcpu=` for a core;
    /// otherwise the baseline plus `features`.
    model: ?[]const u8 = null,
    /// Features added on top of the model or baseline.
    features: []const []const u8 = &.{},
    /// Whether the reference's flags for the core contract multiply-adds
    /// (see `fma.zig`); null follows the hardware.
    fuse_fma: ?bool = null,

    /// Whether two specs build the same code, so one object can serve both.
    pub fn eql(a: Spec, b: Spec) bool {
        if ((a.model == null) != (b.model == null)) return false;
        if (a.model) |m| if (!std.mem.eql(u8, m, b.model.?)) return false;
        if (a.fuse_fma != b.fuse_fma) return false;
        if (a.features.len != b.features.len) return false;
        for (a.features) |fa| {
            for (b.features) |fb| {
                if (std.mem.eql(u8, fa, fb)) break;
            } else return false;
        }
        return true;
    }
};

// ---- x86 -------------------------------------------------------------------
// Flags: `Makefile.x86_64` (`HAVE_*` -> -msse3 -mssse3 -msse4.1 -mavx -mavx2),
// `kernel/Makefile` (`-march=` for the AVX-512 cores; `-mavx2 -mfma` for HASWELL
// and ZEN, the clang convention followed here). FMA3/FMA4 on the AMD cores is
// listed only so an assembler accepts their hand-written microkernels; the
// reference's flags do not enable contraction for them (`fuse_fma` false).

const x86_sse41 = [_][]const u8{ "sse3", "ssse3", "sse4_1" };

pub fn x86(core: coretypes.X86CoreType) Spec {
    return switch (core) {
        .nehalem => .{ .features = &x86_sse41, .fuse_fma = false },
        .sandybridge => .{ .features = &(x86_sse41 ++ [_][]const u8{"avx"}), .fuse_fma = false },
        .bulldozer => .{ .features = &.{ "sse3", "avx", "fma4" }, .fuse_fma = false },
        .piledriver, .steamroller, .excavator => .{ .features = &.{ "sse3", "sse4_1", "avx", "fma", "fma4" }, .fuse_fma = false },
        .haswell, .zen => .{ .features = &(x86_sse41 ++ [_][]const u8{ "avx", "avx2", "fma" }), .fuse_fma = true },
        .skylakex => .{ .model = "skylake_avx512", .fuse_fma = true },
        .cooperlake => .{ .model = "cooperlake", .fuse_fma = true },
        .sapphirerapids => .{ .model = "sapphirerapids", .fuse_fma = true },
        // Older cores: what `getarch.c` lists as `HAVE_*`.
        .prescott, .opteron_sse3, .barcelona => .{ .features = &.{"sse3"}, .fuse_fma = false },
        .core2, .atom, .nano, .bobcat => .{ .features = &.{ "sse3", "ssse3" }, .fuse_fma = false },
        .penryn, .dunnington => .{ .features = &x86_sse41, .fuse_fma = false },
        .yonah => .{ .features = &.{"sse3"}, .fuse_fma = false },
        .opteron, .katmai, .coppermine, .northwood, .banias, .athlon, .p5, .p6, .viac3, .unknown => .{ .fuse_fma = false },
    };
}

// ---- aarch64 (`Makefile.arm64`, `kernel/Makefile`) --------------------------

pub fn arm64(core: coretypes.Arm64CoreType) Spec {
    return switch (core) {
        .armv8, .cortexa53, .cortexa57, .cortexa72, .cortexa73, .falkor, .thunderx, .emag8180, .ft2000, .unknown => .{},
        .cortexa76, .cortexx1 => .{ .features = &.{"v8_2a"} },
        .cortexa510, .cortexa710, .cortexx2 => .{ .features = &.{ "v8_4a", "sve" } },
        .ampere1 => .{ .features = &.{ "v8_6a", "aes", "sha2", "crc", "fullfp16", "sha3", "rand" } },
        .cortexa55, .tsv110, .neoversen1 => .{ .features = &.{"v8_2a"} },
        .thunderx2t99 => .{ .features = &.{"v8_1a"} },
        .thunderx3t110, .vortex => .{ .features = &.{"v8_3a"} },
        .vortexm4 => .{ .features = &.{ "v8_4a", "sme", "sme_f64f64" } },
        .armv8sve => .{ .features = &.{"sve"} },
        .a64fx => .{ .features = &.{ "v8_2a", "sve" } },
        .neoversev1 => .{ .features = &.{ "v8_4a", "sve", "bf16" } },
        .neoversen2 => .{ .features = &.{ "v8_5a", "sve", "sve2", "bf16" } },
        .neoversev2 => .{ .model = "neoverse_v2" },
        .armv9sme => .{ .features = &.{ "v9a", "sve2", "sme", "sme_f64f64" } },
    };
}

// ---- 32-bit arm (`Makefile.arm`) --------------------------------------------

pub fn arm(core: coretypes.ArmCoreType) Spec {
    return switch (core) {
        .armv5 => .{},
        .armv6 => .{ .features = &.{ "v6", "vfp2" } },
        .armv7 => .{ .features = &.{ "v7a", "vfp3" } },
        .neon => .{ .features = &.{ "v7a", "neon" } },
        .vfpv4 => .{ .features = &.{ "v7a", "vfp4" } },
        // `-mfpu=neon` for their `HAVE_NEON`
        .cortexa9, .cortexa15 => .{ .features = &.{ "v7a", "neon" } },
    };
}

// ---- power (`Makefile.power`) -----------------------------------------------

pub fn power(core: coretypes.PowerCoreType) Spec {
    return switch (core) {
        .unknown => .{},
        .power6 => .{ .model = "pwr6" },
        .power8 => .{ .model = "pwr8", .features = &.{"vsx"} },
        .power9 => .{ .model = "pwr9", .features = &.{"vsx"} },
        .power10 => .{ .model = "pwr10", .features = &.{"vsx"} },
        // `-mcpu` of the static targets in `Makefile.power`
        .power3 => .{ .model = "pwr3" },
        .power4 => .{ .model = "pwr4" },
        .power5 => .{ .model = "pwr5" },
        .ppc970, .cell => .{ .model = "970" },
        .ppc440, .ppc440fp2 => .{ .model = "440" },
        .ppcg4 => .{ .model = "g4" },
    };
}

// ---- riscv64 (`Makefile.riscv64`; `RISCV64_OPT` is `rv64imafdcv`) ------------

pub fn riscv(core: coretypes.RiscvCoreType) Spec {
    return switch (core) {
        .riscv64_generic => .{ .features = &.{ "m", "a", "f", "d", "c" } },
        .riscv64_zvl128b => .{ .features = &.{ "m", "a", "f", "d", "c", "v", "zvl128b" } },
        .riscv64_zvl256b => .{ .features = &.{ "m", "a", "f", "d", "c", "v", "zvl256b" } },
        // `Makefile.riscv64`: the T-Head core has the pre-ratification vector
        // extension, which this toolchain cannot assemble (it has no
        // `xtheadvector`), SiFive's U74 none, and SiFive's X280 the ratified one
        .c910v => .{ .features = &.{ "m", "a", "f", "d", "c" } },
        .u74 => .{ .features = &.{ "m", "a", "f", "d", "c" } },
        .x280 => .{ .features = &.{ "m", "a", "f", "d", "c", "v", "zvl512b" } },
    };
}

// ---- mips (`Makefile.mips64`: `MSA_FLAGS = -mmsa -mfp64`) --------------------

pub fn mips(core: coretypes.MipsCoreType) Spec {
    return switch (core) {
        // MSA needs the revision-5 feature bit explicitly; GCC's `-mmsa` implies
        // it, LLVM's integrated assembler does not.
        .loongson3r4 => .{ .features = &.{ "msa", "fp64", "mips64r5" } },
        // `Makefile.mips`, `Makefile.mips64`: release 5 and release 6 cores with
        // the SIMD architecture; the 24K and 1004K have none
        .mips24k, .mips1004k => .{ .model = "mips32r2" },
        .p5600 => .{ .model = "p5600", .features = &.{ "msa", "fp64" } },
        .i6400, .i6500, .p6600 => .{ .model = "mips64r6", .features = &.{ "msa", "fp64" } },
        // `-march=mips64r2 -mtune=loongson3a`; the 4-operand `madd.s` needs it
        .loongson3r3 => .{ .features = &.{"mips64r2"} },
        .mips64_generic, .unknown => .{},
    };
}

// ---- loongarch --------------------------------------------------------------

pub fn loongarch(core: coretypes.LoongarchCoreType) Spec {
    return switch (core) {
        .la264, .loongson2k1000 => .{ .features = &.{"lsx"} },
        .la464, .loongson3r5 => .{ .features = &.{ "lsx", "lasx" } },
        .la64_generic, .loongsongeneric, .unknown => .{},
    };
}

// ---- s390x (`Makefile.zarch`: `-march=z13|z14 -mzvector`) ---------------------

pub fn s390x(core: coretypes.S390xCoreType) Spec {
    return switch (core) {
        .z13 => .{ .model = "z13" },
        .z14 => .{ .model = "z14" },
        .zarch_generic => .{},
    };
}

// ---- sparc (`Makefile.sparc`: `-mcpu=v9`) -----------------------------------

pub fn sparc(core: coretypes.SparcCoreType) Spec {
    return switch (core) {
        .sparc => .{ .model = "v9" },
        .sparcv7 => .{},
    };
}

// ---- wasm (`Makefile.wasm`: `-msimd128`) ---------------------------------------

pub fn wasm(core: coretypes.WasmCoreType) Spec {
    return switch (core) {
        .wasm128_generic => .{ .features = &.{"simd128"} },
        .wasm_generic => .{},
    };
}

/// The spec of any architecture's core type.
pub fn of(core: anytype) Spec {
    return switch (@TypeOf(core)) {
        coretypes.X86CoreType => x86(core),
        coretypes.Arm64CoreType => arm64(core),
        coretypes.ArmCoreType => arm(core),
        coretypes.PowerCoreType => power(core),
        coretypes.RiscvCoreType => riscv(core),
        coretypes.MipsCoreType => mips(core),
        coretypes.LoongarchCoreType => loongarch(core),
        coretypes.S390xCoreType => s390x(core),
        coretypes.SparcCoreType => sparc(core),
        coretypes.WasmCoreType => wasm(core),
        else => .{},
    };
}

// ---- what the build makes objects for ---------------------------------------
// `build/blaskern.zig` reads these instead of the registries, which would make
// its compile depend on every kernel source.

/// The core type of `arch`'s registry (`kernels/__root__.zig`'s `ForArch`).
pub fn CoreType(comptime arch: std.Target.Cpu.Arch) type {
    return switch (arch) {
        .x86_64, .x86 => coretypes.X86CoreType,
        .aarch64, .aarch64_be => coretypes.Arm64CoreType,
        .arm, .armeb, .thumb, .thumbeb => coretypes.ArmCoreType,
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => coretypes.PowerCoreType,
        .riscv32, .riscv32be, .riscv64, .riscv64be => coretypes.RiscvCoreType,
        .s390x => coretypes.S390xCoreType,
        .mips, .mipsel, .mips64, .mips64el => coretypes.MipsCoreType,
        .loongarch32, .loongarch64 => coretypes.LoongarchCoreType,
        .sparc, .sparc64 => coretypes.SparcCoreType,
        .wasm32, .wasm64 => coretypes.WasmCoreType,
        else => coretypes.GenericCoreType,
    };
}

/// The detector of `arch` (`detect()`): light, unlike `cpu/<arch>.zig`, which
/// adds the kernels.
pub fn Detector(comptime arch: std.Target.Cpu.Arch) type {
    return switch (arch) {
        .x86_64, .x86 => @import("_detect/x86.zig"),
        .aarch64, .aarch64_be => @import("_detect/arm64.zig"),
        .arm, .armeb, .thumb, .thumbeb => @import("_detect/arm.zig"),
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => @import("_detect/power.zig"),
        .riscv32, .riscv32be, .riscv64, .riscv64be => @import("_detect/riscv.zig"),
        .s390x => @import("_detect/s390x.zig"),
        .mips, .mipsel, .mips64, .mips64el => @import("_detect/mips.zig"),
        .loongarch32, .loongarch64 => @import("_detect/loongarch.zig"),
        .sparc, .sparc64 => @import("_detect/sparc.zig"),
        .wasm32, .wasm64 => @import("_detect/wasm.zig"),
        else => @import("_detect/generic.zig"),
    };
}

/// The cores of `arch` with kernels: each registry's `hasKernels`, as compiled
/// for that architecture. Derived from the registries; re-derive rather than
/// edit (the test in `kernels/__root__.zig` fails with the right list).
pub fn kernelCores(comptime arch: std.Target.Cpu.Arch) []const CoreType(arch) {
    return switch (arch) {
        .x86_64 => &.{ .prescott, .atom, .core2, .penryn, .dunnington, .nehalem, .opteron, .opteron_sse3, .barcelona, .nano, .sandybridge, .bobcat, .bulldozer, .piledriver, .haswell, .steamroller, .excavator, .zen, .skylakex, .cooperlake, .sapphirerapids },
        .x86 => &.{ .katmai, .coppermine, .northwood, .prescott, .banias, .atom, .core2, .penryn, .dunnington, .nehalem, .athlon, .opteron, .opteron_sse3, .barcelona, .nano, .sandybridge, .bobcat, .bulldozer, .piledriver, .haswell, .zen, .p5, .p6, .viac3, .yonah },
        .aarch64 => &.{ .armv8, .cortexa53, .cortexa57, .cortexa72, .cortexa73, .falkor, .thunderx, .thunderx2t99, .tsv110, .emag8180, .neoversen1, .neoversev1, .neoversev2, .neoversen2, .thunderx3t110, .cortexa55, .armv8sve, .a64fx, .armv9sme, .vortex, .vortexm4, .cortexa76, .cortexa510, .cortexa710, .cortexx1, .cortexx2, .ampere1, .ft2000 },
        .arm => &.{ .armv6, .armv7, .cortexa9, .cortexa15 },
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => &.{ .power6, .power8, .power9, .power10, .power3, .power4, .power5, .ppc970, .cell, .ppc440, .ppc440fp2, .ppcg4 },
        .mips, .mipsel, .mips64, .mips64el => &.{ .loongson3r3, .loongson3r4, .p5600, .i6400, .i6500, .p6600 },
        .loongarch64 => &.{ .la64_generic, .la264, .la464 },
        .riscv64 => &.{ .riscv64_zvl128b, .x280 },
        .s390x => &.{ .z13, .z14 },
        .sparc, .sparc64 => &.{ .sparcv7, .sparc },
        .wasm32, .wasm64 => &.{.wasm128_generic},
        else => &.{},
    };
}

/// Whether `core`'s kernels run on `cpu`.
pub fn runsOn(cpu: std.Target.Cpu, core: anytype) bool {
    return switch (@TypeOf(core)) {
        // the SVE cores need the extension
        coretypes.Arm64CoreType => switch (core) {
            .a64fx, .neoversev1, .cortexa510, .cortexa710, .cortexx2 => cpu.has(.aarch64, .sve),
            else => true,
        },
        // a core's code is for one word size; the Loongson 3A kernels are
        // release 2 code, which release 6 does not run; the SIMD kernels need
        // the release that has the SIMD architecture
        coretypes.MipsCoreType => runsOnMips(cpu, core),
        // the 32-bit cores have no 64-bit code; the vector ones need their
        // instructions
        coretypes.PowerCoreType => switch (core) {
            .power8, .power9 => cpu.has(.powerpc, .vsx),
            .power10 => cpu.has(.powerpc, .vsx) and cpu.has(.powerpc, .isa_v31_instructions),
            .ppc440, .ppc440fp2, .ppcg4 => !cpu.arch.isPowerPC64(),
            else => true,
        },
        else => true,
    };
}

fn runsOnMips(cpu: std.Target.Cpu, core: coretypes.MipsCoreType) bool {
    const wide = cpu.arch.isMIPS64();
    return switch (core) {
        .mips24k, .mips1004k, .p5600 => !wide and (core != .p5600 or cpu.has(.mips, .mips32r5)),
        .i6400, .i6500, .p6600 => wide and cpu.has(.mips, .mips64r6),
        .loongson3r3, .loongson3r4 => wide and !cpu.has(.mips, .mips64r6),
        .mips64_generic, .unknown => wide,
    };
}

// ---- tests -----------------------------------------------------------------

/// Every model and feature named for `Core` exists in `Target` (a
/// `std.Target.<arch>` namespace).
fn check(comptime Target: type, comptime Core: type) !void {
    inline for (@typeInfo(Core).@"enum".field_names) |core_name| {
        const s = comptime of(@field(Core, core_name));
        if (s.model) |m| try std.testing.expect(@hasDecl(Target.cpu, m));
        for (s.features) |name| try std.testing.expect(std.meta.stringToEnum(Target.Feature, name) != null);
    }
}

test "every model and feature named in a spec exists for its architecture" {
    try check(std.Target.x86, coretypes.X86CoreType);
    try check(std.Target.aarch64, coretypes.Arm64CoreType);
    try check(std.Target.arm, coretypes.ArmCoreType);
    try check(std.Target.powerpc, coretypes.PowerCoreType);
    try check(std.Target.riscv, coretypes.RiscvCoreType);
    try check(std.Target.mips, coretypes.MipsCoreType);
    try check(std.Target.loongarch, coretypes.LoongarchCoreType);
    try check(std.Target.s390x, coretypes.S390xCoreType);
    try check(std.Target.sparc, coretypes.SparcCoreType);
    try check(std.Target.wasm, coretypes.WasmCoreType);
}

test "cores that share flags compare equal, others do not" {
    try std.testing.expect(x86(.haswell).eql(x86(.zen)));
    try std.testing.expect(x86(.piledriver).eql(x86(.excavator)));
    try std.testing.expect(!x86(.haswell).eql(x86(.sandybridge)));
    try std.testing.expect(!x86(.skylakex).eql(x86(.cooperlake)));
}
