//! MIPS64 (n64): the sources read their arguments from the n64 registers
//! and schedule their own delay slots (`.set noreorder`).

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    const out = try a.alloc([]const u8, 1);
    out[0] = try std.fmt.allocPrint(a, "-D{s}", .{core});
    return out;
}

pub const mips64: arch.Arch = .{
    .name = "mips64",
    .triple = "mips64el-linux-gnuabi64",
    .config_h = "#define ARCH_MIPS64 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_MIPS64", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .core_defines = coreDefines,
};
