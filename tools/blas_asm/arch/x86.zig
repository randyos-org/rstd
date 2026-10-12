//! 32-bit x86: the x87 kernels. Every argument is on the stack, the same
//! on every OS.

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    const out = try a.alloc([]const u8, 1);
    out[0] = try std.fmt.allocPrint(a, "-D{s}", .{core});
    return out;
}

pub const x86: arch.Arch = .{
    .name = "x86",
    .triple = "x86-linux-gnu",
    .config_h = "#define ARCH_X86 1\n#define BINARY32 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_X86", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .core_defines = coreDefines,
    .align_bytes = true,
};
