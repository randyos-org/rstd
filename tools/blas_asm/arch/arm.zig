//! 32-bit ARM with VFP. Each source has a path for the hard-float procedure
//! call standard (`__ARM_PCS_VFP`) and one for the soft: the `hf` and `sf`
//! variants. It assumes ARM state.

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
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
    .variants = &.{
        .{ .name = "hf", .defines = &.{"-D__ARM_PCS_VFP"} },
        .{ .name = "sf", .defines = &.{"-U__ARM_PCS_VFP"} },
    },
};
