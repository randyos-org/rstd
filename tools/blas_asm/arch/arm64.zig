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
    out[0] = try a.print("-D{s}", .{core});
    return out;
}

pub const arm64: arch.Arch = .{
    .name = "arm64",
    .triple = "aarch64-linux-gnu",
    .config_h = "#define ARCH_ARM64 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_ARM64", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .core_defines = coreDefines,
    .patches = &patches,
    .check = platformRegister,
};

/// x18 is the platform register on Windows (the TEB) and Apple (reserved; the
/// kernel may zero it on a context switch), and the sources are built once for
/// every OS. The gemvs only save and restore it, as callee-saved, which
/// writes back the value they found; any other use stops the generator.
fn platformRegister(source: []const u8, statement: []const u8) void {
    if (!usesX18(statement)) return;
    map.die("{s} uses x18, the Windows and Apple platform register: {s}", .{ source, statement });
}

/// Whether `statement` uses x18 other than to save or restore it on the stack.
pub fn usesX18(statement: []const u8) bool {
    const t = std.mem.trim(u8, statement, " ");
    var it = std.mem.tokenizeAny(u8, t, " ,[]!#");
    const op = it.next() orelse return false;
    var uses = false;
    while (it.next()) |tok| uses = uses or std.mem.eql(u8, tok, "x18") or std.mem.eql(u8, tok, "w18");
    if (!uses) return false;
    const saves = for ([_][]const u8{ "stp", "ldp", "str", "ldr" }) |m| {
        if (std.mem.eql(u8, op, m)) break true;
    } else false;
    return !(saves and std.mem.find(u8, t, "[sp") != null);
}

/// Fixes for upstream sources (OpenBLAS v0.3.34-317).
const patches = [_]arch.Patch{
    // a zero alpha returned at once, so `0 * x` never reached `y` and NaN and
    // infinities in `x` were lost, unlike the other architectures' kernels
    .{ .source = "axpy.S", .find = "\tfcmp\tDA, #0.0\n\tbeq\t.Laxpy_kernel_L999\n", .replace = "" },
    .{ .source = "zaxpy.S", .find = "\tfcmp\tDA_R, #0.0\n\tbne\t.L1\n\tfcmp\tDA_I, #0.0\n\tbeq\t.Lzaxpy_kernel_L999\n", .replace = "" },
    .{ .source = "daxpy_thunderx2t99.S", .find = "\tfcmp\tDA, #0.0\n\tbeq\t.Ldaxpy_kernel_L999\n", .replace = "" },
};
