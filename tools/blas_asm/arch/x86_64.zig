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
    try out.append(a, try a.print("-D{s}", .{macro}));
    try out.append(a, try a.print("-DCORE_{s}", .{macro}));
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
    // Windows: read the flag from incy's slot, 64(%rsp); the tenth argument
    // is at 80(%rsp), so a zero alpha stored zeros where flag 1 keeps NaN
    .{ .source = "scal_sse.S", .find = "\tmovq\t64(%rsp), %r9\n", .replace = "\tmovq\t80(%rsp), %r9\n" },
    .{ .source = "scal_sse2.S", .find = "\tmovq\t64(%rsp), FLAG\n", .replace = "\tmovq\t80(%rsp), FLAG\n" },
    .{ .source = "scal_atom.S", .find = "\tmovq    64(%rsp), %r9\n", .replace = "\tmovq\t80(%rsp), %r9\n" },
    // `comiss` sets ZF for NaN too, so a NaN alpha took the zero path;
    // `scal_sse2.S` has the `jp`
    .{ .source = "scal_sse.S", .find = "\tjne\t.L100\t\t# Alpha != ZERO\n", .replace = "\tjne\t.L100\t\t# Alpha != ZERO\n\tjp\t.L100\n" },
    // `zscal_atom.S` never read the flag (whether a zero alpha keeps NaN), so
    // a zero alpha always stored zeros; read and test it as `scal_atom.S` does.
    // The `jp`: a NaN imaginary part beside a zero real one counted as zero
    .{
        .source = "zscal_atom.S",
        .find = "\tmovq\t56(%rsp), INCX\n#endif\n",
        .replace = "\tmovq\t56(%rsp), INCX\n\tmovq\t88(%rsp), %r9\n#else\n\tmovq\t24(%rsp), %r9\n#endif\n",
    },
    .{
        .source = "zscal_atom.S",
        .find = "\tjne\t.L30\t\t# Alpha_i != ZERO\n",
        .replace = "\tjne\t.L30\t\t# Alpha_i != ZERO\n\tjp\t.L30\n\tcmpq\t$1, %r9\n\tje\t.L30\n",
    },
    // the index search's last step counts the last element as the match
    // unread; with an even count it was read already, so a NaN extreme
    // (which no element equals) returned n + 1
    inRange("iamax_sse.S", ".L27"),
    inRange("iamax_sse.S", ".L47"),
    inRange("iamax_sse.S", ".L97"),
    inRange("iamax_sse2.S", ".L28"),
    inRange("iamax_sse2.S", ".L68"),
    inRange("iamax_sse2.S", ".L98"),
    inRange("izamax_sse.S", ".L47"),
    inRange("izamax_sse.S", ".L87"),
    inRange("izamax_sse2.S", ".L36"),
    inRange("izamax_sse2.S", ".L76"),
    // with a zero `incy`, the strided loops load a block of `y` before
    // storing it, so `n` updates of one element act as fewer;
    // `zaxpy_sse2.S` took its one-at-a-time loop (which never stepped) only
    // with both strides zero. Zero strides take one element at a time.
    .{ .source = "zaxpy_sse2.S", .find = "\tcmpq\t$0, INCX\n\tjne  \t.L59\n", .replace = "\tcmpq\t$0, INCX\n\tje\t.L58\n" },
    .{
        .source = "zaxpy_sse2.S",
        .find = "\tmovhpd\t%xmm8,   1 * SIZE(YY)\n\n\tdecq %rax\n\tjg\t.L58\n",
        .replace = "\tmovhpd\t%xmm8,   1 * SIZE(YY)\n\taddq\tINCX, X\n\taddq\tINCY, Y\n\taddq\tINCY, YY\n\n\tdecq %rax\n\tjg\t.L58\n",
    },
    atomZeroIncy("axpy_atom.S", axpy_atom_one, "3"),
    atomZeroIncy("zaxpy_atom.S", zaxpy_atom_one, "2"),
    // OPTERON loaded by adding to zero, so -0 came out +0 and a signaling
    // NaN quiet; load as the other cores do
    .{ .source = "copy_sse.S", .find = opteronLoad("addps"), .replace = "#ifdef OPTERON\n#define LOAD(OFFSET, ADDR, REG)\t\tmovaps\tOFFSET(ADDR), REG\n" },
    .{ .source = "copy_sse2.S", .find = opteronLoad("addpd"), .replace = "#ifdef OPTERON\n#define LOAD(OFFSET, ADDR, REG)\t\tmovups\tOFFSET(ADDR), REG\n" },
    .{ .source = "zcopy_sse.S", .find = opteronLoad("addps"), .replace = "#ifdef OPTERON\n#define LOAD(OFFSET, ADDR, REG)\t\tmovaps\tOFFSET(ADDR), REG\n" },
    .{ .source = "zcopy_sse2.S", .find = opteronLoad("addpd"), .replace = "#ifdef OPTERON\n#define LOAD(OFFSET, ADDR, REG)\t\tmovaps\tOFFSET(ADDR), REG\n" },
    // the last two rows (`.L17`) and two columns (`.L20`) zeroed the upper
    // half of `A`'s operand but not of `x`'s, which holds the next elements,
    // so an infinity there met a 0 and the sum got NaN
    .{ .source = "symv_L_sse.S", .find = "\tpxor\t xtemp2, xtemp2\n", .replace = "\tpxor\t xtemp2, xtemp2\n\tmovlhps  xtemp2, xtemp1\n" },
    .{
        .source = "symv_L_sse.S",
        .find = "\tmovaps\t 0 * SIZE(NEW_X, IS, SIZE), atemp4\n",
        .replace = "\tmovaps\t 0 * SIZE(NEW_X, IS, SIZE), atemp4\n\txorps\t xt1, xt1\n\tmovlhps  xt1, atemp4\n",
    },
};

/// A zero `incy` branches from the strided entry `.L20` to `one`, a loop
/// `.L28` put before it; `shift` is the block's.
fn atomZeroIncy(comptime source: []const u8, comptime one: []const u8, comptime shift: []const u8) arch.Patch {
    const head = ".L20:\n\tmovq\tY, YY\n\n\tmovq\tM,  %rax\n";
    const blocks = "\tsarq\t$" ++ shift ++ ", %rax\n";
    return .{ .source = source, .find = head ++ blocks, .replace = one ++ head ++ "\tcmpq\t$0, INCY\n\tje\t.L28\n" ++ blocks };
}

/// `y += alpha*x` for `M` elements, one at a time, as the tails do one.
const axpy_atom_one = ".L28:\n\tmovsd\t(X), %xmm0\n\tmulsd\tALPHA, %xmm0\n\taddsd\t(Y), %xmm0\n\tmovsd\t%xmm0, (Y)\n" ++
    "\taddq\tINCX, X\n\taddq\tINCY, Y\n\tdecq\t%rax\n\tjg\t.L28\n\tjmp\t.L29\n\n";
const zaxpy_atom_one = ".L28:\n\tmovsd\t 0 * SIZE(X), %xmm0\n\tmovsd\t 1 * SIZE(X), %xmm1\n\tmovsd\t 0 * SIZE(Y), %xmm8\n\tmovsd\t 1 * SIZE(Y), %xmm9\n" ++
    "\tmovaps\t %xmm0,  %xmm2\n\tmulsd\tALPHA_R, %xmm0\n\tmovaps\t %xmm1,  %xmm3\n\tmulsd\tALPHA_R, %xmm1\n\tmulsd\tALPHA_I, %xmm3\n\tmulsd\tALPHA_I, %xmm2\n" ++
    "\taddsd\t%xmm0, %xmm8\n\tADD2\t%xmm1, %xmm9\n\tADD1\t%xmm3, %xmm8\n\taddsd\t%xmm2, %xmm9\n" ++
    "\tmovsd\t%xmm8,   0 * SIZE(Y)\n\tmovsd\t%xmm9,   1 * SIZE(Y)\n" ++
    "\taddq\t INCX, X\n\taddq\t INCY, Y\n\tdecq\t%rax\n\tjg\t.L28\n\tjmp\t.L999\n\n";

fn opteronLoad(comptime add: []const u8) []const u8 {
    return "#ifdef OPTERON\n#define LOAD(OFFSET, ADDR, REG)\t\txorps\tREG, REG; " ++ add ++ "\tOFFSET(ADDR), REG\n";
}

/// The index search's fallback at `label`: count the last element only when
/// it is left over (`M` odd).
fn inRange(comptime source: []const u8, comptime label: []const u8) arch.Patch {
    return .{ .source = source, .find = "\n" ++ label ++ ":\n\tincq\tRET\n", .replace = "\n" ++ label ++ ":\n\ttestq\t$1, M\n\tje\t.L999\n\tincq\tRET\n" };
}
