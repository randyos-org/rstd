//! Which build variant of the converted assembly a target needs, from its C
//! calling convention (`Target.cCallingConvention`), not its OS: x86_64 UEFI
//! calls like Windows. The names are the `Variant`s of
//! `tools/blas_asm/arch/<arch>.zig`; `registry.Kernels` fails to compile
//! for a core that lacks the variant, so the two cannot drift apart.
//!
//! Import-free, since the build script reads it too (`build/blaskern.zig`).

const std = @import("std");
const sysinfo = @import("builtin");

pub const Result = union(enum) {
    /// The variant's name (empty where the architecture has only one).
    variant: []const u8,
    /// No kernels are converted for this architecture or word size; the
    /// portable ones serve it.
    none,
    /// Kernels exist for the architecture but not this ABI.
    unsupported: []const u8,
};

pub fn variantOf(t: *const std.Target) Result {
    const cc = t.cCallingConvention() orelse return .none;
    // only arm has a soft-float build; the others read floating-point
    // arguments from FP registers
    if (!t.cpu.arch.isArm() and !t.cpu.arch.isThumb() and !hardFloat(t))
        return .{ .unsupported = "a soft-float ABI passes floats in integer registers; the kernels read FP registers" };
    return switch (t.cpu.arch) {
        .x86_64 => switch (cc) {
            .x86_64_sysv => .{ .variant = "sysv" },
            .x86_64_win => .{ .variant = "win" },
            .x86_64_x32 => .{ .unsupported = "x32 passes 32-bit pointers and longs; the kernels take 64-bit ones" },
            else => unreachable,
        },
        .x86 => .{ .variant = "" },
        .arm, .armeb, .thumb, .thumbeb => switch (cc) {
            .arm_aapcs_vfp => .{ .variant = "hf" },
            else => .{ .variant = "sf" },
        },
        .aarch64, .aarch64_be, .s390x, .loongarch64, .riscv64 => .{ .variant = "" },
        .wasm32, .wasm64 => .{ .variant = "" },
        .mips64, .mips64el => switch (cc) {
            .mips64_n64 => .{ .variant = "" },
            else => .{ .unsupported = "the kernels are n64; this is n32" },
        },
        // only the `@Vector` kernels of the MSA core (p5600); no assembly
        .mips, .mipsel => .{ .variant = "o32" },
        .sparc64 => .{ .variant = "64" },
        .sparc => .{ .variant = "32" },
        .powerpc64le => .{ .variant = "64" },
        // Zig emits ELFv2 (`.abiversion 2`) here even where `cc` says ELFv1
        .powerpc64 => .{ .variant = "64be" },
        .powerpc => .{ .variant = "32" },
        .powerpcle => .{ .variant = "32le" },
        else => .none,
    };
}

/// Whether floating-point arguments travel in FP registers.
fn hardFloat(t: *const std.Target) bool {
    if (t.abi.float() == .soft) return false;
    const f = t.cpu.features;
    return switch (t.cpu.arch) {
        .riscv64, .riscv64be => std.Target.riscv.featureSetHas(f, .d),
        .loongarch64 => std.Target.loongarch.featureSetHas(f, .d),
        .aarch64, .aarch64_be => std.Target.aarch64.featureSetHas(f, .fp_armv8),
        else => true,
    };
}

/// The variant for the compiled target; empty when there is none (nothing is
/// linked then, or the build stopped).
pub fn variantName(t: *const std.Target) []const u8 {
    return switch (variantOf(t)) {
        .variant => |v| v,
        else => "",
    };
}

test "variants follow the calling convention" {
    // resolving a query needs an `Io`, which freestanding builds lack
    if (comptime sysinfo.target.os.tag == .freestanding) return error.SkipZigTest;
    const cases = [_]struct { []const u8, []const u8 }{
        .{ "x86_64-linux-gnu", "sysv" },
        .{ "x86_64-linux-musl", "sysv" },
        .{ "x86_64-macos", "sysv" },
        .{ "x86_64-windows-gnu", "win" },
        .{ "x86_64-windows-msvc", "win" },
        .{ "x86_64-uefi", "win" },
        .{ "arm-linux-gnueabihf", "hf" },
        .{ "arm-linux-gnueabi", "sf" },
        .{ "thumb-linux-musleabihf", "hf" },
        .{ "aarch64-macos", "" },
        .{ "sparc64-linux-gnu", "64" },
        .{ "sparc-linux-gnu", "32" },
        .{ "powerpc64le-linux-gnu", "64" },
        .{ "powerpc64-linux-gnu", "64be" },
        .{ "powerpc64-linux-musl", "64be" },
        .{ "powerpc-linux-gnu", "32" },
        .{ "mips64el-linux-gnuabi64", "" },
        .{ "mips-linux-gnu", "o32" },
        .{ "mipsel-linux-gnueabihf", "o32" },
    };
    for (cases) |c| {
        const q = try std.Target.Query.parse(.{ .arch_os_abi = c[0] });
        const t = try std.zig.system.resolveTargetQuery(std.testing.io, q);
        try std.testing.expectEqualStrings(c[1], variantName(&t));
    }
    const unsupported = [_]struct { []const u8, []const u8 }{
        .{ "x86_64-linux-gnux32", "baseline" },
        .{ "mips64-linux-gnuabin32", "baseline" },
        // soft float: by the ABI, or no FP registers
        .{ "powerpc-linux-gnueabi", "baseline" },
        .{ "loongarch64-linux-gnusf", "baseline" },
        .{ "mipsel-linux-gnueabi", "baseline" },
        .{ "riscv64-linux-gnu", "baseline-d" },
        .{ "aarch64-linux-gnu", "baseline-fp_armv8" },
    };
    for (unsupported) |c| {
        const q = try std.Target.Query.parse(.{ .arch_os_abi = c[0], .cpu_features = c[1] });
        const t = try std.zig.system.resolveTargetQuery(std.testing.io, q);
        try std.testing.expect(variantOf(&t) == .unsupported);
    }
    for ([_][]const u8{ "riscv32-linux-gnu", "x86_16-freestanding" }) |s| {
        const t = try std.zig.system.resolveTargetQuery(std.testing.io, try std.Target.Query.parse(.{ .arch_os_abi = s }));
        try std.testing.expect(variantOf(&t) == .none);
    }
}
