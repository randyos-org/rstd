//! LoongArch64: the sources read their arguments from the LP64D registers.

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    const out = try a.alloc([]const u8, 1);
    out[0] = try std.fmt.allocPrint(a, "-D{s}", .{core});
    return out;
}

pub const loongarch64: arch.Arch = .{
    .name = "loongarch64",
    .triple = "loongarch64-linux-gnu",
    .config_h = "#define ARCH_LOONGARCH64 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_LOONGARCH64", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .core_defines = coreDefines,
};
