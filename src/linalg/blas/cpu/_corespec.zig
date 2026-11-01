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
    /// Features taken out of the model: AltiVec on the cores without VSX,
    /// whose lanes flush denormals and which LLVM's auto-vectorizer would
    /// otherwise use for scalar loops (`util.vectorIeee`).
    without: []const []const u8 = &.{},
    /// Whether the reference's flags for the core contract multiply-adds
    /// (see `fma.zig`); null follows the hardware.
    fuse_fma: ?bool = null,

    /// Whether two specs build the same code, so one object can serve both.
    pub fn eql(a: Spec, b: Spec) bool {
        if ((a.model == null) != (b.model == null)) return false;
        if (a.model) |m| if (!std.mem.eql(u8, m, b.model.?)) return false;
        if (a.fuse_fma != b.fuse_fma) return false;
        return sameNames(a.features, b.features) and sameNames(a.without, b.without);
    }

    fn sameNames(a: []const []const u8, b: []const []const u8) bool {
        if (a.len != b.len) return false;
        for (a) |fa| {
            for (b) |fb| {
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
        .nehalem, .penryn, .dunnington => .{ .features = &x86_sse41, .fuse_fma = false },
        .sandybridge => .{ .features = &(x86_sse41 ++ [_][]const u8{"avx"}), .fuse_fma = false },
        .bulldozer => .{ .features = &.{ "sse3", "avx", "fma4" }, .fuse_fma = false },
        .piledriver, .steamroller, .excavator => .{ .features = &.{ "sse3", "sse4_1", "avx", "fma", "fma4" }, .fuse_fma = false },
        .haswell, .zen => .{ .features = &(x86_sse41 ++ [_][]const u8{ "avx", "avx2", "fma" }), .fuse_fma = true },
        .skylakex => .{ .model = "skylake_avx512", .fuse_fma = true },
        .cooperlake => .{ .model = "cooperlake", .fuse_fma = true },
        .sapphirerapids => .{ .model = "sapphirerapids", .fuse_fma = true },
        // Older cores: what `getarch.c` lists as `HAVE_*`.
        .prescott, .opteron_sse3, .barcelona, .yonah => .{ .features = &.{"sse3"}, .fuse_fma = false },
        .core2, .atom, .nano, .bobcat => .{ .features = &.{ "sse3", "ssse3" }, .fuse_fma = false },
        .opteron, .katmai, .coppermine, .northwood, .banias, .athlon, .p5, .p6, .viac3, .unknown => .{ .fuse_fma = false },
    };
}

// ---- aarch64 (`Makefile.arm64`, `kernel/Makefile`) --------------------------
// `kernel/Makefile`'s `-march=armv8.4-a+sve+bf16` for NEOVERSEN2, NEOVERSEV1,
// ARMV8SVE and A64FX is not followed; each takes its `Makefile.arm64` flags.

pub fn arm64(core: coretypes.Arm64CoreType) Spec {
    return switch (core) {
        .armv8, .cortexa53, .cortexa57, .cortexa72, .cortexa73, .falkor, .thunderx, .emag8180, .ft2000, .unknown => .{},
        .cortexa76, .cortexx1, .cortexa55, .tsv110, .neoversen1 => .{ .features = &.{"v8_2a"} },
        .cortexa510, .cortexa710, .cortexx2 => .{ .features = &.{ "v8_4a", "sve" } },
        .ampere1 => .{ .features = &.{ "v8_6a", "aes", "sha2", "crc", "fullfp16", "sha3", "rand" } },
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
        // Cortex-A9 and A15: `-mfpu=neon` for their `HAVE_NEON`
        .neon, .cortexa9, .cortexa15 => .{ .features = &.{ "v7a", "neon" } },
        .vfpv4 => .{ .features = &.{ "v7a", "vfp4" } },
    };
}

// ---- power (`Makefile.power`) -----------------------------------------------

const no_altivec = [_][]const u8{"altivec"};

pub fn power(core: coretypes.PowerCoreType) Spec {
    return switch (core) {
        .unknown => .{},
        .power6 => .{ .model = "pwr6", .without = &no_altivec },
        .power8 => .{ .model = "pwr8", .features = &.{"vsx"} },
        .power9 => .{ .model = "pwr9", .features = &.{"vsx"} },
        .power10 => .{ .model = "pwr10", .features = &.{"vsx"} },
        // `Makefile.power` gives the static targets no `-mcpu`; these are their
        // CPUs' LLVM models. LLVM's POWER3 to POWER5 have AltiVec too
        .power3 => .{ .model = "pwr3", .without = &no_altivec },
        .power4 => .{ .model = "pwr4", .without = &no_altivec },
        .power5 => .{ .model = "pwr5", .without = &no_altivec },
        .ppc970, .cell => .{ .model = "970", .without = &no_altivec },
        .ppc440, .ppc440fp2 => .{ .model = "440" },
        .ppcg4 => .{ .model = "g4", .without = &no_altivec },
    };
}

// ---- riscv64 (`Makefile.riscv64`; `RISCV64_OPT` is `rv64imafdcv`) ------------

pub fn riscv(core: coretypes.RiscvCoreType) Spec {
    return switch (core) {
        // `Makefile.riscv64`: the T-Head core has the pre-ratification vector
        // extension, which this toolchain cannot assemble (it has no
        // `xtheadvector`), SiFive's U74 none, and SiFive's X280 the ratified one
        .riscv64_generic, .c910v, .u74 => .{ .features = &.{ "m", "a", "f", "d", "c" } },
        .riscv64_zvl128b => .{ .features = &.{ "m", "a", "f", "d", "c", "v", "zvl128b" } },
        .riscv64_zvl256b => .{ .features = &.{ "m", "a", "f", "d", "c", "v", "zvl256b" } },
        .x280 => .{ .features = &.{ "m", "a", "f", "d", "c", "v", "zvl512b" } },
    };
}

// ---- mips (`Makefile.mips64`: `MSA_FLAGS = -mmsa -mfp64`) --------------------

/// The FPU on the named cores' models: upstream lists only what every
/// configuration of a core has (Zig 2c5f1044fc, 2026-10-28), so `fpu_d` and
/// `msa` are opt-in there; older Zig has no `fpu_d` and implies both.
const mips_fpu = if (@hasField(std.Target.mips.Feature, "fpu_d")) [_][]const u8{"fpu_d"} else [_][]const u8{};
const mips_simd = [_][]const u8{"msa"} ++ mips_fpu;

pub fn mips(core: coretypes.MipsCoreType) Spec {
    return switch (core) {
        // MSA needs the revision-5 feature bit explicitly; GCC's `-mmsa` implies
        // it, LLVM's integrated assembler does not. Zig adds `-mfp64` for MSA.
        // `-march=loongson3a`: GCC fuses `madd.fmt` there, as GS464 and GS464E
        // do (`fma.fused_madd`)
        .loongson3r4 => .{ .features = &.{ "msa", "mips64r5" }, .fuse_fma = true },
        // `Makefile.mips`, `Makefile.mips64`: release 5 and release 6 cores with
        // the SIMD architecture; the 24K and 1004K have none
        .mips24k, .mips1004k => .{ .model = "mips32r2" },
        .p5600 => .{ .model = "p5600", .features = &mips_simd },
        .i6400, .i6500, .p6600 => .{ .model = "mips64r6", .features = &.{"msa"} },
        // `Makefile.system`'s `-march=loongson3a`, release 2 to gcc; the
        // 4-operand `madd.s` needs it
        .loongson3r3 => .{ .features = &.{"mips64r2"}, .fuse_fma = true },
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
        .x86 => &.{ .katmai, .coppermine, .northwood, .prescott, .banias, .atom, .core2, .penryn, .dunnington, .nehalem, .athlon, .opteron, .opteron_sse3, .barcelona, .nano, .bobcat, .p5, .p6, .viac3, .yonah },
        .aarch64 => &.{ .armv8, .cortexa53, .cortexa57, .cortexa72, .cortexa73, .falkor, .thunderx, .thunderx2t99, .tsv110, .emag8180, .neoversen1, .neoversev1, .neoversev2, .neoversen2, .thunderx3t110, .cortexa55, .armv8sve, .a64fx, .armv9sme, .vortex, .vortexm4, .cortexa76, .cortexa510, .cortexa710, .cortexx1, .cortexx2, .ampere1, .ft2000 },
        .arm => &.{ .armv6, .armv7, .cortexa9, .cortexa15 },
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => &.{ .power6, .power8, .power9, .power10, .power3, .power4, .power5, .ppc970, .cell, .ppc440, .ppc440fp2, .ppcg4 },
        .mips, .mipsel, .mips64, .mips64el => &.{ .loongson3r3, .loongson3r4, .p5600, .i6400, .i6500, .p6600 },
        .loongarch64 => &.{ .la64_generic, .la264, .la464 },
        .riscv64 => &.{ .riscv64_zvl128b, .riscv64_zvl256b, .x280 },
        .s390x => &.{ .z13, .z14 },
        .sparc, .sparc64 => &.{ .sparcv7, .sparc },
        .wasm32, .wasm64 => &.{.wasm128_generic},
        else => &.{},
    };
}

/// The core whose object serves `core` on `arch`. On 32-bit x86 the AVX cores
/// and Zen take `Makefile.system`'s 32-bit fallbacks (Nehalem; Barcelona for
/// the Bulldozer family and Zen), whose kernels and blocking they share there,
/// so they get no object of their own.
pub fn objectCore(comptime arch: std.Target.Cpu.Arch, core: CoreType(arch)) CoreType(arch) {
    return switch (arch) {
        .x86 => switch (core) {
            .sandybridge, .haswell, .skylakex, .cooperlake, .sapphirerapids => .nehalem,
            .bulldozer, .piledriver, .steamroller, .excavator, .zen => .barcelona,
            else => core,
        },
        else => core,
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
            .ppc440, .ppc440fp2, .ppcg4 => !cpu.arch.isPowerpc64(),
            else => true,
        },
        // a module's features are fixed when it is built: no SIMD128, no SIMD core
        coretypes.WasmCoreType => core == .wasm_generic or cpu.has(.wasm, .simd128),
        else => true,
    };
}

fn runsOnMips(cpu: std.Target.Cpu, core: coretypes.MipsCoreType) bool {
    const wide = cpu.arch.isMips64();
    return switch (core) {
        .mips24k, .mips1004k, .p5600 => !wide and (core != .p5600 or cpu.has(.mips, .mips32r5)),
        .i6400, .i6500, .p6600 => wide and cpu.has(.mips, .mips64r6),
        .loongson3r3, .loongson3r4 => wide and !cpu.has(.mips, .mips64r6),
        .mips64_generic, .unknown => wide,
    };
}

/// The query a core's kernel object is built with, for a target of `query`
/// resolved to `cpu` and `os`: the spec's model, or the baseline plus its
/// features, less those it goes without. A target short of some of the
/// baseline's ISA is the base instead: an x86 one without SSE2 (the Pentium 4
/// is the baseline), a macOS one older than Nehalem, a SPARC V9 without VIS.
/// So a Pentium III build's objects compute `f64` on the x87, as its cores and
/// its main library do.
pub fn objectQuery(query: std.Target.Query, os: std.Target.Os, cpu: std.Target.Cpu, spec: Spec) std.Target.Query {
    const arch = cpu.arch;
    var q = query;
    q.cpu_features_sub = featureSet(arch, spec.without);
    q.cpu_features_add = featureSet(arch, spec.features);
    if (spec.model) |name| {
        q.cpu_model = .{ .explicit = findModel(arch, name) };
    } else if (!isaOf(arch, cpu.features).isSuperSetOf(isaOf(arch, std.Target.Cpu.baseline(arch, os).features))) {
        q.cpu_model = .{ .explicit = cpu.model };
        q.cpu_features_add.addFeatureSet(cpu.features);
    } else {
        q.cpu_model = .baseline;
    }
    return q;
}

/// `set` without the features `notIsa` names.
fn isaOf(arch: std.Target.Cpu.Arch, set: std.Target.Cpu.Feature.Set) std.Target.Cpu.Feature.Set {
    var out = set;
    for (arch.allFeaturesList(), 0..) |f, i| {
        if (notIsa(f.name)) out.removeFeature(@intCast(i));
    }
    return out;
}

fn findModel(arch: std.Target.Cpu.Arch, name: []const u8) *const std.Target.Cpu.Model {
    for (arch.allCpuModels()) |m| {
        if (std.mem.eql(u8, m.name, name)) return m;
    }
    std.debug.panic("no CPU model '{s}' for {s}", .{ name, @tagName(arch) });
}

fn featureSet(arch: std.Target.Cpu.Arch, names: []const []const u8) std.Target.Cpu.Feature.Set {
    var set: std.Target.Cpu.Feature.Set = .empty;
    const all = arch.allFeaturesList();
    for (names) |name| {
        const i = for (all, 0..) |feat, k| {
            if (std.mem.eql(u8, feat.name, name)) break k;
        } else std.debug.panic("no CPU feature '{s}' for {s}", .{ name, @tagName(arch) });
        set.addFeature(@intCast(i));
    }
    return set;
}

// ---- what the cores' CPUs have ----------------------------------------------
// The specs say how a core is built; `machine` says what the CPUs that run it
// have, from the hardware rather than the reference's flags. `excess` holds
// every kernel object to it: the build fails on a violation
// (`build/blaskern.zig`), and a test tries each architecture's targets.

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

/// Features that say nothing about what code can execute: scheduling hints
/// (`-mtune`), which models carry along, and supervisor-mode protections.
pub fn notIsa(name: []const u8) bool {
    const prefixes = [_][]const u8{ "false_deps", "fast", "slow", "prefer", "idiv", "lea_", "no_bypass", "macrofusion", "tuning", "use_", "fuse", "balance" };
    for (prefixes) |p| if (std.mem.startsWith(u8, name, p)) return true;
    for ([_][]const u8{ "smap", "smep" }) |s| if (std.mem.eql(u8, name, s)) return true;
    return false;
}

/// `spec`'s features on `arch` as Zig resolves them: its model's (else the
/// generic one's) with what they imply, then those it adds, less those it goes
/// without.
pub fn specFeatures(arch: std.Target.Cpu.Arch, spec: Spec) std.Target.Cpu.Feature.Set {
    const model = if (spec.model) |name| findModel(arch, name) else std.Target.Cpu.Model.generic(arch);
    return resolve(arch, model, featureSet(arch, spec.features), featureSet(arch, spec.without));
}

/// A CPU's features from its model and a query's additions and removals, in
/// `std.zig.system.resolveTargetQuery`'s order.
fn resolve(
    arch: std.Target.Cpu.Arch,
    model: *const std.Target.Cpu.Model,
    add: std.Target.Cpu.Feature.Set,
    sub: std.Target.Cpu.Feature.Set,
) std.Target.Cpu.Feature.Set {
    const all = arch.allFeaturesList();
    var set = model.toCpu(arch).features;
    set.removeFeatureSet(sub);
    set.addFeatureSet(add);
    set.populateDependencies(all);
    set.removeFeatureSet(sub);
    return set;
}

/// A core's feature that its spec adds and some of its CPUs lack.
pub const Waiver = struct { core: []const u8, feature: []const u8 };

/// The features a spec adds that `machine` says some CPUs of the core lack,
/// each awaiting a decision: drop it from the spec, or have `runnable` require
/// it and add it to `machine`. `excess` passes them; keep the list empty.
pub const undecided = [_]Waiver{};

/// The features an object built with `object` uses that a CPU running `core`
/// (`machine`) may lack and the target, which the whole program needs, does
/// not promise, `undecided` aside. Empty where every such CPU runs the object.
pub fn excess(core: anytype, arch: std.Target.Cpu.Arch, object: std.Target.Cpu.Feature.Set, target: std.Target.Cpu.Feature.Set) std.Target.Cpu.Feature.Set {
    var set = isaOf(arch, object);
    set.removeFeatureSet(specFeatures(arch, machine(core)));
    set.removeFeatureSet(target);
    for (undecided) |w| {
        if (std.mem.eql(u8, w.core, @tagName(core))) set.removeFeatureSet(featureSet(arch, &.{w.feature}));
    }
    return set;
}

// ---- tests -----------------------------------------------------------------

/// Every model and feature named for `Core` exists in `Target` (a
/// `std.Target.<arch>` namespace).
fn check(comptime Target: type, comptime Core: type) !void {
    inline for (@typeInfo(Core).@"enum".field_names) |core_name| {
        inline for (.{ comptime of(@field(Core, core_name)), comptime machine(@field(Core, core_name)) }) |s| {
            if (s.model) |m| try std.testing.expect(@hasDecl(Target.cpu, m));
            for (s.features ++ s.without) |name| try std.testing.expect(std.meta.stringToEnum(Target.Feature, name) != null);
        }
    }
}

test "every model and feature named in a spec or machine exists for its architecture" {
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

/// How many of `arch`'s kernel objects, on `os`, need a feature a CPU they
/// serve may lack (`excess`), each printed. The targets are the baseline and
/// every kernel core's own CPU (`machine`), the weakest a build may name.
fn countExcess(comptime arch: std.Target.Cpu.Arch, os_tag: std.Target.Os.Tag) usize {
    const Model = std.Target.Cpu.Model;
    const os: std.Target.Os = .{ .tag = os_tag, .version_range = .{ .none = {} } };
    const cores = kernelCores(arch);
    var bad: usize = 0;
    for (0..cores.len + 1) |t| {
        const m: Spec = if (t == 0) .{ .model = Model.baseline(arch, os).name } else machine(cores[t - 1]);
        const model = if (m.model) |name| findModel(arch, name) else Model.generic(arch);
        var add = featureSet(arch, m.features);
        // the 32-bit cores' CPUs as 64-bit targets
        if (arch == .x86_64) add.addFeatureSet(featureSet(arch, &.{"64bit"}));
        const query: std.Target.Query = .{
            .cpu_arch = arch,
            .os_tag = os_tag,
            .cpu_model = .{ .explicit = model },
            .cpu_features_add = add,
            .cpu_features_sub = featureSet(arch, m.without),
        };
        const cpu: std.Target.Cpu = .{ .arch = arch, .model = model, .features = resolve(arch, model, add, query.cpu_features_sub) };
        for (cores) |core| {
            if (!runsOn(cpu, core)) continue;
            const q = objectQuery(query, os, cpu, of(core));
            const object_model = switch (q.cpu_model) {
                .explicit => |e| e,
                else => Model.baseline(arch, os),
            };
            const ex = excess(core, arch, resolve(arch, object_model, q.cpu_features_add, q.cpu_features_sub), cpu.features);
            if (ex.isEmpty()) continue;
            bad += 1;
            std.debug.print("{t}-{t} for {s}: {t}'s object needs", .{ arch, os_tag, model.name, core });
            for (arch.allFeaturesList(), 0..) |f, i| {
                if (ex.isEnabled(@intCast(i))) std.debug.print(" {s}", .{f.name});
            }
            std.debug.print("\n", .{});
        }
    }
    return bad;
}

test "no kernel object needs a feature the CPUs of its cores may lack" {
    const Os = std.Target.Os.Tag;
    var bad: usize = 0;
    inline for (.{
        .{ .x86, &[_]Os{ .linux, .freestanding } },
        .{ .x86_64, &[_]Os{ .linux, .macos } },
        .{ .aarch64, &[_]Os{ .linux, .macos } },
        .{ .arm, &[_]Os{ .linux, .freestanding } },
        .{ .powerpc, &[_]Os{ .linux, .freestanding } },
        .{ .powerpc64, &[_]Os{ .linux, .freestanding } },
        .{ .powerpc64le, &[_]Os{.linux} },
        .{ .mips, &[_]Os{.linux} },
        .{ .mipsel, &[_]Os{.linux} },
        .{ .mips64, &[_]Os{ .linux, .openbsd } },
        .{ .mips64el, &[_]Os{.linux} },
        .{ .riscv64, &[_]Os{.linux} },
        .{ .loongarch64, &[_]Os{.linux} },
        .{ .s390x, &[_]Os{.linux} },
        .{ .sparc, &[_]Os{ .linux, .freestanding } },
        .{ .sparc64, &[_]Os{.linux} },
        .{ .wasm32, &[_]Os{.wasi} },
    }) |a| {
        for (a[1]) |os| bad += countExcess(a[0], os);
    }
    try std.testing.expectEqual(@as(usize, 0), bad);
}

test "every undecided feature is one its core's spec adds and its machine lacks" {
    const os: std.Target.Os = .{ .tag = .freestanding, .version_range = .{ .none = {} } };
    for (undecided) |w| {
        var found = false;
        inline for (.{ .x86_64, .aarch64, .arm, .powerpc64, .riscv64, .mips64, .loongarch64, .s390x, .sparc64, .wasm32 }) |arch| {
            if (std.meta.stringToEnum(CoreType(arch), w.core)) |core| {
                found = true;
                const f = featureSet(arch, &.{w.feature});
                const base = std.Target.Cpu.Model.baseline(arch, os).name;
                try std.testing.expect(specFeatures(arch, .{ .model = of(core).model orelse base, .features = of(core).features }).isSuperSetOf(f));
                try std.testing.expect(!specFeatures(arch, machine(core)).isSuperSetOf(f));
            }
        }
        try std.testing.expect(found);
    }
}
