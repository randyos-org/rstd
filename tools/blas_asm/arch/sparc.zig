//! SPARC V9: each source `save`s its own register window and returns with
//! `return %i7 + 8`; register names differ per pointer width (`__64BIT__`):
//! the `32` and `64` variants.

const std = @import("std");
const arch = @import("__root__.zig");
const map = @import("../kernel_map.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    // what getarch's `config.h` defines, where it knows the core
    if (try map.getarch(core)) |d| return d;
    const out = try a.alloc([]const u8, 2);
    out[0] = try std.fmt.allocPrint(a, "-D{s}", .{core});
    out[1] = "-DSPARC";
    return out;
}

pub const sparc: arch.Arch = .{
    .name = "sparc",
    .triple = "sparc64-linux-gnu",
    .config_h = "#define ARCH_SPARC 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_SPARC", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .core_defines = coreDefines,
    .variants = &.{
        .{ .name = "32", .defines = &.{} },
        .{ .name = "64", .defines = &.{"-D__64BIT__"} },
    },
};
