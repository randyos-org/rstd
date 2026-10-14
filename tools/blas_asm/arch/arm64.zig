//! AArch64: the sources read their arguments from the AAPCS64 registers.

const std = @import("std");
const arch = @import("__root__.zig");
const map = @import("../kernel_map.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    // what getarch's `config.h` defines, where it knows the core
    if (try map.getarch(core)) |d| return d;
    const out = try a.alloc([]const u8, 1);
    out[0] = try std.fmt.allocPrint(a, "-D{s}", .{core});
    return out;
}

pub const arm64: arch.Arch = .{
    .name = "arm64",
    .triple = "aarch64-linux-gnu",
    .config_h = "#define ARCH_ARM64 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_ARM64", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .core_defines = coreDefines,
};
