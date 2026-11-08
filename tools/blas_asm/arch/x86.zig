//! 32-bit x86: the x87 kernels. Every argument is on the stack, the same
//! on every OS.

const std = @import("std");
const arch = @import("__root__.zig");
const map = @import("../openblas.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    // what getarch's `config.h` defines, where it knows the core; the AVX cores
    // take their 32-bit fallback's (their `KERNEL.<core>` already matches it)
    if (try map.getarch(map.target32(core))) |d| return d;
    const out = try a.alloc([]const u8, 1);
    out[0] = try a.print("-D{s}", .{core});
    return out;
}

pub const x86: arch.Arch = .{
    .name = "x86",
    .triple = "x86-linux-gnu",
    .config_h = "#define ARCH_X86 1\n#define BINARY32 1\n" ++ @import("x86_64.zig").dtb_default,
    .defines = &.{ "-DASSEMBLER", "-DARCH_X86", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .core_defines = coreDefines,
    .align_bytes = true,
    .patches = &patches,
    .retired = &retired,
};

/// Fixes for bugs in upstream sources (OpenBLAS v0.3.34-317).
const patches = [_]arch.Patch{
    // with both strides zero, `swap.S` returned at once, right only for an
    // even `n`; with either zero, both load two or four elements before
    // storing any. Zero strides take one element at a time.
    .{
        .source = "swap.S",
        .find = "\tcmpl\t$0, %ebx\n\tjne\t.L141\n\tcmpl\t$0, %ecx\n\tjne\t.L141\n/* INCX == 0 and INCY == 0 */\t\n\tjmp\t.L27\n",
        .replace = "\tmovl\t%edx, %eax\n\tcmpl\t$0, %ebx\n\tje\t.L35\n\tcmpl\t$0, %ecx\n\tje\t.L35\n",
    },
    .{ .source = "zswap.S", .find = zswap_strided, .replace = zswap_one ++ zswap_strided ++ "\tcmpl\t$0, %ebx\n\tje\t.L36\n\tcmpl\t$0, %ecx\n\tje\t.L36\n" },
};

/// Fixes to sources no longer converted: the extremes.
const retired = [_]arch.Patch{
    // the minima took the new element unless it was `>` the pick, so a NaN
    // (unordered) replaced it; compare the pick with the new element, and
    // the maxima's `fcmovbe` and `cmovnbe` keep NaN out
    minFirst("amax.S", 18),
    minFirst("iamax.S", 18),
    minFirst("zamax.S", 10),
    minFirst("izamax.S", 10),
    .{ .source = "amax.S", .find = "#define FMOV\tfcmovnbe\n", .replace = "#define FMOV\tfcmovbe\n" },
    .{ .source = "zamax.S", .find = "#define FMOV\tfcmovnbe\n", .replace = "#define FMOV\tfcmovbe\n" },
    .{ .source = "iamax.S", .find = "#define FMOV\tfcmovnbe\n#define IMOV\tcmovb\n", .replace = "#define FMOV\tfcmovbe\n#define IMOV\tcmovnbe\n" },
    .{ .source = "izamax.S", .find = "#define FMOV\tfcmovnb\n#define IMOV\tcmovb\n", .replace = "#define FMOV\tfcmovbe\n#define IMOV\tcmovnbe\n" },
};

const zswap_strided = "/* INCX != 1 or INCY != 1 */\n\n.L14:\n\tmovl\t%edx, %eax\n";

/// `.L35`'s swap of one element, as a loop.
const zswap_one = ".L36:\n#if defined(DOUBLE) || defined(XDOUBLE)\n" ++
    "\tFLD\t1 * SIZE(%esi)\n\tFLD\t0 * SIZE(%esi)\n\tFLD\t1 * SIZE(%edi)\n\tFLD\t0 * SIZE(%edi)\n" ++
    "\tFST\t0 * SIZE(%esi)\n\tFST\t1 * SIZE(%esi)\n\tFST\t0 * SIZE(%edi)\n\tFST\t1 * SIZE(%edi)\n" ++
    "#else\n\tfldl\t0 * SIZE(%esi)\n\tfldl\t0 * SIZE(%edi)\n\tfstpl\t0 * SIZE(%esi)\n\tfstpl\t0 * SIZE(%edi)\n#endif\n" ++
    "\taddl\t%ebx, %esi\n\taddl\t%ecx, %edi\n\tdecl\t%eax\n\tjg\t.L36\n\tjmp\t.L27\n\n";

/// For a minimum, `fcomi` with the pick in `st(0)`, so the flags say how it
/// compares with the new element.
fn minFirst(comptime source: []const u8, comptime count: usize) arch.Patch {
    return .{
        .source = source,
        .find = "\tfcomi\t%st(1), %st\n",
        .replace = "#ifdef USE_MIN\n\tfxch\n#endif\n\tfcomi\t%st(1), %st\n#ifdef USE_MIN\n\tfxch\n#endif\n",
        .count = count,
    };
}
