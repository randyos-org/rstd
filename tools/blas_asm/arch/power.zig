//! PowerPC and POWER, little-endian ELFv2. Each source names its argument
//! registers per pointer width (`__64BIT__`): the `32` and `64` variants.

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    const out = try a.alloc([]const u8, 1);
    out[0] = try std.fmt.allocPrint(a, "-D{s}", .{core});
    return out;
}

pub const power: arch.Arch = .{
    .name = "power",
    .triple = "powerpc64le-linux-gnu",
    .config_h = "#define ARCH_POWER 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_POWER", "-DOS_LINUX", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-D__BYTE_ORDER__=__ORDER_LITTLE_ENDIAN__" },
    .core_defines = coreDefines,
    // the Blue Gene double-precision FPU instructions are not in LLVM
    .unsupported = &.{ "axpy_hummer.S", "zaxpy_hummer.S" },
    .variants = &.{
        .{ .name = "32", .defines = &.{} },
        .{ .name = "64", .defines = &.{"-D__64BIT__"}, .patches = &zdot_gpr },
    },
};

/// The complex dots return in f1/f2 as ELFv2 says, but Zig (0.18.0-dev.249)
/// reads a two-float struct from r3/r4. Each source already has the copy into
/// r3/r4 for gfortran; turn it on, keeping f1/f2 for callers that follow the
/// ABI. Drop this once Zig does.
const zdot_gpr = [_]arch.Patch{
    .{ .source = "zdot.S", .find = "#if defined(F_INTERFACE) && defined(F_INTERFACE_GFORT)", .replace = "#if 1" },
    .{ .source = "zdot_cell.S", .find = "#if defined(F_INTERFACE) && defined(F_INTERFACE_GFORT)", .replace = "#if 1" },
};
