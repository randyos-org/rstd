//! What each core type is compiled with, from OpenBLAS's `DYNAMIC_ARCH`
//! flags (`Makefile.<arch>`'s primary-compiler branch, `getarch.c`'s `HAVE_*`
//! tables, `kernel/Makefile`'s overrides): plain data the build
//! (`build/blaskern.zig`) and the runtime share. Names are `std.Target`'s
//! (a test checks); only the ISA, not `-mtune=`.

const std = @import("std");
const coretypes = @import("coretypes.zig");
const machines = @import("machine.zig");

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
        // without NEON, which 0.17.0's LLVM has `v7a` imply (as `vfpv4`)
        .armv7 => .{ .features = &.{ "v7a", "vfp3" }, .without = &.{"neon"} },
        // Cortex-A9 and A15: `-mfpu=neon` for their `HAVE_NEON`
        .neon, .cortexa9, .cortexa15 => .{ .features = &.{ "v7a", "neon" } },
        .vfpv4 => .{ .features = &.{ "v7a", "vfp4" }, .without = &.{"neon"} },
    };
}

// ---- power (`Makefile.power`) -----------------------------------------------

pub const no_altivec = [_][]const u8{"altivec"};

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
pub const mips_simd = [_][]const u8{"msa"} ++ mips_fpu;

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
    return coretypes.CoreType(coretypes.familyOf(arch));
}

/// The detector of `arch` (`select.detect`): light, unlike `cpu.Arch`,
/// which adds the kernels.
pub fn Detector(comptime arch: std.Target.Cpu.Arch) type {
    return DetectorOf(coretypes.familyOf(arch));
}

/// The detector of `family`.
pub fn DetectorOf(comptime family: coretypes.Family) type {
    return switch (family) {
        .x86_64, .x86 => @import("detect/x86.zig"),
        .arm64 => @import("detect/arm64.zig"),
        .arm => @import("detect/arm.zig"),
        .power => @import("detect/power.zig"),
        .riscv => @import("detect/riscv.zig"),
        .s390x => @import("detect/s390x.zig"),
        .mips => @import("detect/mips.zig"),
        .loongarch => @import("detect/loongarch.zig"),
        .sparc => @import("detect/sparc.zig"),
        .wasm => @import("detect/wasm.zig"),
        .generic => @import("detect/generic.zig"),
    };
}

/// The cores of `arch` with kernels: each registry's `hasKernels`, as compiled
/// for that architecture. Derived from the registries; re-derive rather than
/// edit (`cpu/__root__.zig`'s test of it prints the right list).
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

/// Whether a target's features can tell that its CPU runs `core`'s kernels:
/// not the Loongson 3A ones, whose MMI and EXT instructions no Zig feature
/// names. A build for the target never picks those (`-Dblas-cores=target`);
/// detection reaches them where they are linked.
pub fn targetPickable(core: anytype) bool {
    return switch (@TypeOf(core)) {
        coretypes.MipsCoreType => core != .loongson3r3 and core != .loongson3r4,
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
/// Panics on a model or feature name this Zig lacks: a broken core table
/// (`check` tests every name), met at build time.
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

/// Panics when `arch` has no model `name` (a broken core table).
pub fn findModel(arch: std.Target.Cpu.Arch, name: []const u8) *const std.Target.Cpu.Model {
    for (arch.allCpuModels()) |m| {
        if (std.mem.eql(u8, m.name, name)) return m;
    }
    std.debug.panic("no CPU model '{s}' for {s}", .{ name, @tagName(arch) });
}

/// Panics when `arch` lacks a feature in `names` (a broken core table).
pub fn featureSet(arch: std.Target.Cpu.Arch, names: []const []const u8) std.Target.Cpu.Feature.Set {
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

pub const machine = machines.machine;

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
/// Panics as `objectQuery` does.
pub fn specFeatures(arch: std.Target.Cpu.Arch, spec: Spec) std.Target.Cpu.Feature.Set {
    const model = if (spec.model) |name| findModel(arch, name) else std.Target.Cpu.Model.generic(arch);
    return resolve(arch, model, featureSet(arch, spec.features), featureSet(arch, spec.without));
}

/// A CPU's features from its model and a query's additions and removals, in
/// `std.zig.system.resolveTargetQuery`'s order.
pub fn resolve(
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

/// The features an object built with `object` uses that a CPU running `core`
/// (`machine`) may lack and the target, which the whole program needs, does
/// not promise. Empty where every such CPU runs the object.
pub fn excess(core: anytype, arch: std.Target.Cpu.Arch, object: std.Target.Cpu.Feature.Set, target: std.Target.Cpu.Feature.Set) std.Target.Cpu.Feature.Set {
    var set = isaOf(arch, object);
    set.removeFeatureSet(specFeatures(arch, machine(core)));
    set.removeFeatureSet(target);
    return set;
}
