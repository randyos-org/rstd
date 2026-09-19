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
const zigstd = @import("../../../zigstd.zig");
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
        .piledriver, .steamroller, .excavator => .{ .features = &.{ "sse3", "sse4_1", "avx", "fma" }, .fuse_fma = false },
        .haswell, .zen => .{ .features = &(x86_sse41 ++ [_][]const u8{ "avx", "avx2", "fma" }), .fuse_fma = true },
        .skylakex => .{ .model = "skylake_avx512", .fuse_fma = true },
        .cooperlake => .{ .model = "cooperlake", .fuse_fma = true },
        .sapphirerapids => .{ .model = "sapphirerapids", .fuse_fma = true },
        // Older cores: what `getarch.c` lists as `HAVE_*`.
        .prescott, .opteron_sse3, .barcelona => .{ .features = &.{"sse3"}, .fuse_fma = false },
        .core2, .atom, .nano, .bobcat => .{ .features = &.{ "sse3", "ssse3" }, .fuse_fma = false },
        .penryn, .dunnington => .{ .features = &x86_sse41, .fuse_fma = false },
        .opteron, .katmai, .coppermine, .northwood, .banias, .athlon, .unknown => .{ .fuse_fma = false },
    };
}

// ---- aarch64 (`Makefile.arm64`, `kernel/Makefile`) --------------------------

pub fn arm64(core: coretypes.Arm64CoreType) Spec {
    return switch (core) {
        .armv8, .cortexa53, .cortexa57, .cortexa72, .cortexa73, .falkor, .thunderx, .emag8180, .unknown => .{},
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
        .armv7 => .{ .features = &.{ "v7a", "neon" } },
        .neon => .{ .features = &.{ "v7a", "neon" } },
        .vfpv4 => .{ .features = &.{ "v7a", "vfp4" } },
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
    };
}

// ---- riscv64 (`Makefile.riscv64`; `RISCV64_OPT` is `rv64imafdcv`) ------------

pub fn riscv(core: coretypes.RiscvCoreType) Spec {
    return switch (core) {
        .riscv64_generic => .{ .features = &.{ "m", "a", "f", "d", "c" } },
        .riscv64_zvl128b => .{ .features = &.{ "m", "a", "f", "d", "c", "v", "zvl128b" } },
        .riscv64_zvl256b => .{ .features = &.{ "m", "a", "f", "d", "c", "v", "zvl256b" } },
    };
}

// ---- mips (`Makefile.mips64`: `MSA_FLAGS = -mmsa -mfp64`) --------------------

pub fn mips(core: coretypes.MipsCoreType) Spec {
    return switch (core) {
        .loongson3r4 => .{ .features = &.{ "msa", "fp64" } },
        .mips64_generic, .loongson3r3, .unknown => .{},
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
        else => .{},
    };
}

// ---- tests -----------------------------------------------------------------

/// Every model and feature named for `Core` exists in `Target` (a
/// `std.Target.<arch>` namespace).
fn check(comptime Target: type, comptime Core: type) !void {
    inline for (@typeInfo(Core).@"enum".field_names) |core_name| {
        const s = comptime of(@field(Core, core_name));
        if (s.model) |m| try zigstd.testExpectTrue(@hasDecl(Target.cpu, m));
        for (s.features) |name| try zigstd.testExpectTrue(std.meta.stringToEnum(Target.Feature, name) != null);
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
}

test "cores that share flags compare equal, others do not" {
    try zigstd.testExpectTrue(x86(.haswell).eql(x86(.zen)));
    try zigstd.testExpectTrue(x86(.piledriver).eql(x86(.excavator)));
    try zigstd.testExpectTrue(!x86(.haswell).eql(x86(.sandybridge)));
    try zigstd.testExpectTrue(!x86(.skylakex).eql(x86(.cooperlake)));
}
