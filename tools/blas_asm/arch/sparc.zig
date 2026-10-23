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

/// Fixes for bugs in upstream sources (OpenBLAS v0.3.34-317).
const patches = [_]arch.Patch{
    // the index kernels keep a pick per lane, and of equal picks the
    // reduction kept the first lane's, which may be the later index; take
    // the smaller index of equal picks before each step
    .{ .source = "iamax.S", .find = reduce, .replace = reduce_first, .count = 2 },
    .{ .source = "imax.S", .find = reduce, .replace = reduce_first, .count = 2 },
    .{ .source = "izamax.S", .find = reduce, .replace = reduce_first, .count = 2 },
    // with a zero stride, the strided loop loads a block before storing it,
    // so `n` operations on one element act as fewer; zero strides take the
    // one-at-a-time loop instead (utest's `swap_inc_0`, per-core `rot`)
    zeroStrides("swap.S", "\tsra\tN, 3, I\n"),
    zeroStrides("zswap.S", "\tsra\tN, 2, I\n"),
    zeroStrides("rot.S", "\tmov\tX, XX\n\tmov\tY, YY\n\n\tsra\tN, 3, I\n"),
    zeroStrides("zrot.S", "\tmov\tX, XX\n\tmov\tY, YY\n\n\tsra\tN, 2, I\n"),
    // and with a zero `incy`, `axpy` (per-core `axpy` checks)
    zeroStrides("axpy.S", "\tsra\tN, 3, I\n"),
    zeroStrides("zaxpy.S", "\tsra\tN, 2, I\n"),
};

const reduce =
    "\tFCMP\t%fcc0, c2, c1\n\tadd\tv2, 1, v2\n\tFCMP\t%fcc1, c4, c3\n\tadd\tv3, 2, v3\n\tadd\tv4, 3, v4\n\n" ++
    "\tFCMOV\t%fcc0, c2, c1\n\tCMOV\t%fcc0, v2, v1\n\tFCMOV\t%fcc1, c4, c3\n\tCMOV\t%fcc1, v4, v3\n" ++
    "\tFCMP\t%fcc0, c3, c1\n\tCMOV\t%fcc0, v3, v1\n";

const reduce_first =
    "\tFCMP\t%fcc0, c2, c1\n\tadd\tv2, 1, v2\n\tFCMP\t%fcc1, c4, c3\n\tadd\tv3, 2, v3\n\tadd\tv4, 3, v4\n\n" ++
    tie("%fcc0", "v1", "v2") ++ tie("%fcc1", "v3", "v4") ++
    "\tFCMOV\t%fcc0, c2, c1\n\tCMOV\t%fcc0, v2, v1\n\tFCMOV\t%fcc1, c4, c3\n\tCMOV\t%fcc1, v4, v3\n" ++
    "\tFCMP\t%fcc0, c3, c1\n" ++ tie("%fcc0", "v1", "v3") ++ "\tCMOV\t%fcc0, v3, v1\n";

/// `va` becomes `vb` where `fcc` says equal and `vb < va`.
fn tie(comptime fcc: []const u8, comptime va: []const u8, comptime vb: []const u8) []const u8 {
    return "\tmov\t" ++ va ++ ", %o5\n\tcmp\t" ++ vb ++ ", " ++ va ++ "\n" ++
        "#ifdef __64BIT__\n\tmovl\t%xcc, " ++ vb ++ ", %o5\n#else\n\tmovl\t%icc, " ++ vb ++ ", %o5\n#endif\n" ++
        "\tmove\t" ++ fcc ++ ", %o5, " ++ va ++ "\n";
}

/// `sll` as `sllx`, for 64-bit values.
fn wideShifts(a: std.mem.Allocator, statement: []const u8) Error![]const u8 {
    const t = std.mem.trimStart(u8, statement, " ");
    // tabs are spaces by now (`expandTabs`)
    if (!std.mem.startsWith(u8, t, "sll ")) return statement;
    return std.fmt.allocPrint(a, "{s}sllx{s}", .{ statement[0 .. statement.len - t.len], t[3..] });
}

/// Zero strides branch from the strided entry `.LL50`, whose first
/// statements are `head`, to the one-at-a-time loop `.LL56`; it runs once
/// before testing the count, so `n` = 0 returns (`.LL59`).
fn zeroStrides(comptime source: []const u8, comptime head: []const u8) arch.Patch {
    return .{
        .source = source,
        .find = "\n.LL50:\n" ++ head,
        .replace = "\n.LL57:\n\tcmp\tN, 0\n\tble\t.LL59\n\tmov\tN, I\n\tb\t.LL56\n\tnop\n\n.LL50:\n\tcmp\tINCX, 0\n\tbe\t.LL57\n\tnop\n" ++
            "\tcmp\tINCY, 0\n\tbe\t.LL57\n\tnop\n" ++ head,
    };
}

pub const sparc: arch.Arch = .{
    .name = "sparc",
    .triple = "sparc64-linux-gnu",
    .config_h = "#define ARCH_SPARC 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_SPARC", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .core_defines = coreDefines,
    .patches = &patches,
    .variants = &.{
        .{ .name = "32", .defines = &.{}, .triple = "sparc-linux-gnu" },
        .{ .name = "64", .defines = &.{"-D__64BIT__"}, .rewrite = wideShifts },
    },
};
