//! x86_64: the SSE, SSE2 and Atom kernels. The sources read their arguments
//! where the calling convention puts them, with a path each for SysV and
//! Windows (`WINDOWS_ABI`, which upstream's `common.h` sets for `OS_WINNT`):
//! the two are the `sysv` and `win` variants.

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    // `OPTERON_SSE3` is `OPTERON` with SSE3; the headers know one macro
    const macro = if (std.mem.eql(u8, core, "OPTERON_SSE3")) "OPTERON" else core;
    var out: std.ArrayList([]const u8) = .empty;
    try out.append(a, try std.fmt.allocPrint(a, "-D{s}", .{macro}));
    try out.append(a, try std.fmt.allocPrint(a, "-DCORE_{s}", .{macro}));
    try out.appendSlice(a, &.{ "-DHAVE_SSE", "-DHAVE_SSE2" });
    if (!std.mem.eql(u8, core, "OPTERON")) try out.append(a, "-DHAVE_SSE3");
    return out.items;
}

pub const x86_64: arch.Arch = .{
    .name = "x86_64",
    .triple = "x86_64-linux-gnu",
    .config_h = "#define ARCH_X86_64 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_X86_64", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k" },
    .core_defines = coreDefines,
    .align_bytes = true,
    .variants = &.{
        .{ .name = "sysv", .defines = &.{"-DOS_LINUX"} },
        .{ .name = "win", .defines = &.{"-DOS_WINNT"} },
    },
};
