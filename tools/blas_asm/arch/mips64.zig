//! MIPS64 (n64): the sources read their arguments from the n64 registers
//! and schedule their own delay slots (`.set noreorder`).

const arch = @import("__root__.zig");
const arch_table = @import("../arches.zig");

pub const mips64: arch.Arch = .{
    .arch = arch_table.mips64,
    .triple = "mips64el-linux-gnuabi64",
    // getarch's `config.h` defines `__64BIT__` for a 64-bit build; without it
    // `LDARG` is `lw`, which reads a stack argument's high half big-endian
    .config_h = "#define ARCH_MIPS64 1\n#define BINARY64 1\n#define __64BIT__ 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_MIPS64", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .patches = &patches,
    .retired = &retired,
};

/// Fixes for bugs in upstream sources (OpenBLAS v0.3.34-317).
const patches = [_]arch.Patch{
    // with a zero stride, the strided loop loads a block before storing it,
    // so `n` swaps of one element act as fewer; zero strides take the
    // one-at-a-time loop (`I` set in each delay slot)
    .{
        .source = "swap.S",
        .find = ".L20:\n\tdsra\tI, N, 3\n",
        .replace = ".L20:\n\tbeqz\tINCX, .L26\n\tmove\tI, N\n\tbeqz\tINCY, .L26\n\tmove\tI, N\n\tdsra\tI, N, 3\n",
    },
    // Loongson's prefetch, `ld $0`, is still a load: it takes a bus error on
    // `f32` data not 8-aligned and faults past the array's end; `pref` cannot
    .{ .source = "axpy_loongson3a.S", .find = "#define PREFETCH_DISTANCE 48\n", .replace = "#define PREFETCH_DISTANCE 48\n" ++ pref },
    .{ .source = "daxpy_loongson3a_simd.S", .find = "#define PREFETCH_DISTANCE 2016\n", .replace = "#define PREFETCH_DISTANCE 2016\n" ++ pref },
};

const pref = "#undef PREFETCHD\n#define PREFETCHD(x) pref 0, x\n";

/// Fixes to sources no longer converted: the extremes.
const retired = firstOnTies("iamax.S", false) ++ firstOnTies("imax.S", false) ++ firstOnTies("izamax.S", false) ++
    firstOnTies("iamin.S", true) ++ firstOnTies("imin.S", true) ++ firstOnTies("izamin.S", true);

/// The index kernels keep a pick per lane, and of equal picks the reduction
/// kept the first lane's, which may be the later index; before each step,
/// take the smaller index of equal picks.
fn firstOnTies(comptime source: []const u8, comptime min: bool) [2]arch.Patch {
    const pairs = if (min) "\tCMPLT\t$fcc0, s2, s1\n\tCMPLT\t$fcc1, s4, s3\n" else "\tCMPLT\t$fcc0, s1, s2\n\tCMPLT\t$fcc1, s3, s4\n";
    const last = if (min) "\tCMPLT\t$fcc0, s3, s1\n" else "\tCMPLT\t$fcc0, s1, s3\n";
    return .{
        .{ .source = source, .find = ".L998:\n" ++ pairs, .replace = ".L998:\n" ++ tie("s1", "s2", "x1", "x2") ++ tie("s3", "s4", "x3", "x4") ++ pairs },
        .{ .source = source, .find = last, .replace = tie("s1", "s3", "x1", "x3") ++ last },
    };
}

/// `xa` becomes `xb` where `sa == sb` and `xb < xa`.
fn tie(comptime sa: []const u8, comptime sb: []const u8, comptime xa: []const u8, comptime xb: []const u8) []const u8 {
    return "\tCMPEQ\t$fcc0, " ++ sa ++ ", " ++ sb ++ "\n" ++
        "\tslt\tTEMP, " ++ xb ++ ", " ++ xa ++ "\n" ++
        "\tmovf\tTEMP, $0, $fcc0\n" ++
        "\tmovn\t" ++ xa ++ ", " ++ xb ++ ", TEMP\n";
}
