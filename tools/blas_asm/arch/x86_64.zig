//! x86_64: the SSE, SSE2 and Atom kernels. The sources read their arguments
//! where the calling convention puts them, with a path each for SysV and
//! Windows (`WINDOWS_ABI`, which upstream's `common.h` sets for `OS_WINNT`):
//! the two are the `sysv` and `win` variants.

const std = @import("std");
const arch = @import("__root__.zig");
const map = @import("../kernel_map.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    // what getarch's `config.h` defines, where it knows the core: `PRESCOTT`
    // builds as `PENTIUM4`, and the gemvs read `DTB_DEFAULT_ENTRIES`
    if (try map.getarch(core)) |d| return d;
    // `OPTERON_SSE3` is `OPTERON` with SSE3; the headers know one macro
    const macro = if (std.mem.eql(u8, core, "OPTERON_SSE3")) "OPTERON" else core;
    var out: std.ArrayList([]const u8) = .empty;
    try out.append(a, try std.fmt.allocPrint(a, "-D{s}", .{macro}));
    try out.append(a, try std.fmt.allocPrint(a, "-DCORE_{s}", .{macro}));
    try out.appendSlice(a, &.{ "-DHAVE_SSE", "-DHAVE_SSE2" });
    if (!std.mem.eql(u8, core, "OPTERON")) try out.append(a, "-DHAVE_SSE3");
    return out.items;
}

/// getarch's usual value, for a core it has no block for.
pub const dtb_default = "#ifndef DTB_DEFAULT_ENTRIES\n#define DTB_DEFAULT_ENTRIES 64\n#endif\n";

pub const x86_64: arch.Arch = .{
    .name = "x86_64",
    .triple = "x86_64-linux-gnu",
    .config_h = "#define ARCH_X86_64 1\n#define BINARY64 1\n" ++ dtb_default,
    .defines = &.{ "-DASSEMBLER", "-DARCH_X86_64", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k" },
    .core_defines = coreDefines,
    .align_bytes = true,
    .patches = &patches,
    .variants = &.{
        .{ .name = "sysv", .defines = &.{"-DOS_LINUX"} },
        .{ .name = "win", .defines = &.{"-DOS_WINNT"} },
    },
};

/// Fixes for bugs in upstream sources (OpenBLAS v0.3.34-317).
const patches = [_]arch.Patch{
    // Windows passes the double `alpha_i` on the stack; this read it with a
    // single-precision load, as every sibling `zgemv` kernel does not
    .{ .source = "zgemv_t_atom.S", .find = "\tmovss\tOLD_ALPHA_I, %xmm1\n", .replace = "\tmovsd\tOLD_ALPHA_I, %xmm1\n" },
    // Windows passes the address of a 16-byte result first and shifts the
    // arguments along; `zdot_atom.S` read that address as `n`. Taken from
    // `zdot_sse2.S`, which gets it right.
    .{
        .source = "zdot_atom.S",
        .find = "#define N\tARG1\t/* rdi */\n#define X\tARG2\t/* rsi */\n#define INCX\tARG3\t/* rdx */\n#define Y\tARG4\t/* rcx */\n#ifndef WINDOWS_ABI\n#define INCY\tARG5\t/* r8  */\n#else\n#define INCY\t%r10\n#endif\n",
        .replace = "#ifndef WINDOWS_ABI\n#define N\tARG1\n#define X\tARG2\n#define INCX\tARG3\n#define Y\tARG4\n#define INCY\tARG5\n#else\n#define RESULT_ADDRESS\tARG1\n#define N\tARG2\n#define X\tARG3\n#define INCX\tARG4\n#define Y\t%r10\n#define INCY\t%r11\n#endif\n",
    },
    .{
        .source = "zdot_atom.S",
        .find = "#ifdef WINDOWS_ABI\n\tmovq\t40(%rsp), INCY\n#endif\n",
        .replace = "#ifdef WINDOWS_ABI\n\tmovq\t40(%rsp), Y\n\tmovq\t48(%rsp), INCY\n#endif\n",
    },
    .{
        .source = "zdot_atom.S",
        .find = "#endif\n\n\tRESTOREREGISTERS\n\n\tret\n",
        .replace = "#endif\n\n#ifdef WINDOWS_ABI\n\tmovq\tRESULT_ADDRESS, %rax\n\tmovsd\t%xmm0, (%rax)\n\tmovsd\t%xmm1, 8(%rax)\n#endif\n\n\tRESTOREREGISTERS\n\n\tret\n",
    },
};
