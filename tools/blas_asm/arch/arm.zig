//! 32-bit ARM with VFP. Each source has a path for the hard-float procedure
//! call standard (`__ARM_PCS_VFP`) and one for the soft: the `hf` and `sf`
//! variants. It assumes ARM state.

const std = @import("std");
const arch = @import("__root__.zig");
const map = @import("../kernel_map.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    // what getarch's `config.h` defines, where it knows the core
    if (try map.getarch(core)) |d| return d;
    var out: std.ArrayList([]const u8) = .empty;
    try out.append(a, try std.fmt.allocPrint(a, "-D{s}", .{core}));
    try out.appendSlice(a, &.{ "-DARMV7", "-DHAVE_VFP", "-DHAVE_VFPV3" });
    return out.items;
}

pub const arm: arch.Arch = .{
    .name = "arm",
    .triple = "arm-linux-gnueabihf",
    .config_h = "#define ARCH_ARM 1\n#define BINARY32 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_ARM", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .core_defines = coreDefines,
    .patches = &patches,
    .variants = &.{
        .{ .name = "hf", .defines = &.{"-D__ARM_PCS_VFP"} },
        .{ .name = "sf", .defines = &.{"-U__ARM_PCS_VFP"}, .patches = &zig_complex_return },
    },
};

/// Fixes for bugs in upstream sources (OpenBLAS v0.3.34-317).
const patches = [_]arch.Patch{
    // the minima moved on `lt`, which an unordered compare also sets, so a
    // NaN replaced the pick; `mi` is less than only
    minOnMi("amax_vfp.S"),
    minOnMi("iamax_vfp.S"),
};

fn minOnMi(comptime source: []const u8) arch.Patch {
    return .{
        .source = source,
        .find = "#define\tMOVCOND\t\tmovlt\n\n#if\tdefined(DOUBLE)\n\n#define\tVMOVCOND\tvmovlt.f64\n\n#else\n\n#define\tVMOVCOND\tvmovlt.f32\n",
        .replace = "#define\tMOVCOND\t\tmovmi\n\n#if\tdefined(DOUBLE)\n\n#define\tVMOVCOND\tvmovmi.f64\n\n#else\n\n#define\tVMOVCOND\tvmovmi.f32\n",
    };
}

/// Not an upstream bug: under the soft-float convention Zig returns a complex
/// value's struct in r0 to r3 where C passes a hidden pointer in r0
/// (issues.md), so the `sf` complex dots take no pointer and return there,
/// as the power `zdot.S` patch does.
const zig_complex_return = [_]arch.Patch{
    .{ .source = "cdot_vfp.S", .find = "#if !defined(__ARM_PCS_VFP)\n#define\tOLD_RETURN_ADDR", .replace = "#if 0\n#define\tOLD_RETURN_ADDR" },
    .{ .source = "cdot_vfp.S", .find = "#if !defined(__ARM_PCS_VFP)\n\tmov\tRETURN_ADDR, OLD_RETURN_ADDR\n", .replace = "#if 0\n\tmov\tRETURN_ADDR, OLD_RETURN_ADDR\n" },
    .{ .source = "cdot_vfp.S", .find = "\tvstm\tRETURN_ADDR, {s0 - s1}\n", .replace = "\tvmov\tr0, s0\n\tvmov\tr1, s1\n" },
    .{ .source = "zdot_vfp.S", .find = "#if !defined(__ARM_PCS_VFP)\n#define\tOLD_RETURN_ADDR", .replace = "#if 0\n#define\tOLD_RETURN_ADDR" },
    .{ .source = "zdot_vfp.S", .find = "#if !defined(__ARM_PCS_VFP)\n\tmov\tRETURN_ADDR, OLD_RETURN_ADDR\n", .replace = "#if 0\n\tmov\tRETURN_ADDR, OLD_RETURN_ADDR\n" },
    .{ .source = "zdot_vfp.S", .find = "\tvstm\tRETURN_ADDR, {d0 - d1}\n", .replace = "\tvmov\tr0, r1, d0\n\tvmov\tr2, r3, d1\n" },
};
