//! PowerPC and POWER. Each source names its argument registers per pointer
//! width (`__64BIT__`) and the prologue and frame follow the ABI: `32` and
//! `32le` (SysV), `64` (little-endian ELFv2) and `64be` (big-endian ELFv2;
//! Zig emits ELFv2 for every 64-bit target, glibc's included).

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

pub const power: arch.Arch = .{
    .name = "power",
    .triple = "powerpc64le-linux-gnu",
    .config_h = "#define ARCH_POWER 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_POWER", "-DOS_LINUX", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k" },
    .core_defines = coreDefines,
    .patches = &(flag_slot ++ zero_strides),
    // the Blue Gene double-precision FPU instructions are not in LLVM
    .unsupported = &.{ "axpy_hummer.S", "zaxpy_hummer.S" },
    .variants = &.{
        .{ .name = "32", .defines = &.{}, .triple = "powerpc-linux-gnu" },
        .{ .name = "32le", .defines = &.{}, .triple = "powerpcle-linux-gnu" },
        .{ .name = "64", .defines = &.{"-D__64BIT__"}, .patches = &zdot_gpr },
        .{ .name = "64be", .defines = &.{"-D__64BIT__"}, .patches = &zdot_gpr, .triple = "powerpc64-linux-musl" },
    },
};

/// `scal.S` and `zscal.S` read the flag argument (keep NaN for a zero alpha)
/// at Linux ELFv1's stack offset, so ELFv2 and 32-bit read past the caller's
/// arguments. `FRAMESLOT` has every ABI's offset; FreeBSD's ELFv2 one matches.
const flag_slot = [_]arch.Patch{
    .{
        .source = "scal.S",
        .find = "#if defined(__FreeBSD__)\n    LDLONG  FLAG, 104(SP)\n#else\n    LDLONG  FLAG, 120(SP)\n#endif\n",
        .replace = "#ifdef __64BIT__\n    LDLONG  FLAG, FRAMESLOT(1)(SP)\n#else\n    LDLONG  FLAG, FRAMESLOT(0)(SP)\n#endif\n",
    },
    .{
        .source = "zscal.S",
        .find = "#if defined(__FreeBSD__)\n    LDLONG  FLAG, 112(SP)\n#else\n\tLDLONG\tFLAG, 128(SP)\n#endif\n",
        .replace = "#ifdef __64BIT__\n    LDLONG  FLAG, FRAMESLOT(2)(SP)\n#else\n    LDLONG  FLAG, FRAMESLOT(0)(SP)\n#endif\n",
    },
};

/// With a zero stride, the strided loop (`LL(110)`) loads a block before
/// storing it, so `n` operations on one element act as fewer (utest's
/// `zdrot_inc_0`); zero strides take the one-at-a-time loop, `LL(160)`.
const zero_strides = [_]arch.Patch{
    zeroStrides("rot.S", "3"),
    zeroStrides("zrot.S", "2"),
    zeroStrides("swap.S", "4"),
    zeroStrides("zswap.S", "3"),
    zeroStrides("rot_ppc440.S", "3"),
    zeroStrides("zrot_ppc440.S", "2"),
};

fn zeroStrides(comptime source: []const u8, comptime shift: []const u8) arch.Patch {
    const blocks = "\tsrawi.\tr0, N, " ++ shift ++ "\n\tmtspr\tCTR,  r0\n\tbeq-\tLL(150)\n";
    return .{
        .source = source,
        .find = blocks,
        .replace = "\tmtspr\tCTR, N\n\tcmpwi\tcr0, INCX, 0\n\tbeq-\tcr0, LL(160)\n" ++
            "\tcmpwi\tcr0, INCY, 0\n\tbeq-\tcr0, LL(160)\n" ++ blocks,
    };
}

/// The complex dots return in f1/f2 as ELFv2 says, but Zig (0.18.0-dev.249)
/// reads a two-float struct from r3/r4. Each source already has the copy into
/// r3/r4 for gfortran; turn it on, keeping f1/f2 for callers that follow the
/// ABI. Drop this once Zig does.
const zdot_gpr = [_]arch.Patch{
    .{ .source = "zdot.S", .find = "#if defined(F_INTERFACE) && defined(F_INTERFACE_GFORT)", .replace = "#if 1" },
    .{ .source = "zdot_cell.S", .find = "#if defined(F_INTERFACE) && defined(F_INTERFACE_GFORT)", .replace = "#if 1" },
};
