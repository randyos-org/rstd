//! SPARC V9: each source `save`s its own register window and returns with
//! `return %i7 + 8`; register names differ per pointer width (`__64BIT__`):
//! the `32` and `64` variants. The `64` variant shifts with `sllx`: the
//! sources scale 64-bit strides with `sll`, and under qemu (sun4u) `sll`
//! keeps only the low 32 bits, so a negative stride became a huge positive
//! one (`issues.md`). `sllx` is right either way.

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

/// `sll` as `sllx`, for 64-bit values.
fn wideShifts(a: std.mem.Allocator, statement: []const u8) Error![]const u8 {
    const t = std.mem.trimStart(u8, statement, " ");
    // tabs are spaces by now (`expandTabs`)
    if (!std.mem.startsWith(u8, t, "sll ")) return statement;
    return std.fmt.allocPrint(a, "{s}sllx{s}", .{ statement[0 .. statement.len - t.len], t[3..] });
}

pub const sparc: arch.Arch = .{
    .name = "sparc",
    .triple = "sparc64-linux-gnu",
    .config_h = "#define ARCH_SPARC 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_SPARC", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .core_defines = coreDefines,
    .variants = &.{
        .{ .name = "32", .defines = &.{}, .triple = "sparc-linux-gnu" },
        .{ .name = "64", .defines = &.{"-D__64BIT__"}, .rewrite = wideShifts },
    },
};
