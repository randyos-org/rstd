//! LoongArch64: the sources read their arguments from the LP64D registers.

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

/// Fixes for bugs in the LSX/LASX sources, still in OpenBLAS v0.3.34-317;
/// the qemu rig's per-core tests found each one.
const patches = [_]arch.Patch{
    // the double result returned through `fmov.s`, leaving its upper half
    // undefined; `MOV` is `fmov.d` for double
    .{ .source = "sum_lsx.S", .find = "fmov.s $f0, $f16", .replace = "MOV $f0, $f16" },
    .{ .source = "sum_lasx.S", .find = "fmov.s $f0, $f16", .replace = "MOV $f0, $f16" },
    .{ .source = "csum_lsx.S", .find = "fmov.s $f0, $f16", .replace = "MOV $f0, $f16" },
    .{ .source = "csum_lasx.S", .find = "fmov.s $f0, $f16", .replace = "MOV $f0, $f16" },
    // double's contiguous loop added with `xvfadd.s`, as if eight floats
    .{
        .source = "sum_lasx.S",
        .find = "    xvfadd.s res1, res1, VX0\n#ifdef DOUBLE\n    xvld VX1, X, 32\n    xvfadd.s res1, res1, VX1\n#endif\n",
        .replace = "#ifdef DOUBLE\n    xvfadd.d res1, res1, VX0\n    xvld VX1, X, 32\n    xvfadd.d res1, res1, VX1\n#else\n    xvfadd.s res1, res1, VX0\n#endif\n",
    },
    // single's final reduction added the last block (still in `res2`) again,
    // and casum added lane 6 where lane 7 belongs
    .{ .source = "casum_lasx.S", .find = "xvfadd.s res2, res1, res2", .replace = "xvor.v res2, res1, res1", .count = 2 },
    .{ .source = "csum_lasx.S", .find = "xvfadd.s res2, res1, res2", .replace = "xvor.v res2, res1, res1", .count = 2 },
    .{
        .source = "casum_lasx.S",
        .find = "    xvfadd.s res1, VX2, res1\n    xvfadd.s res1, VX2, res1\n",
        .replace = "    xvfadd.s res1, VX2, res1\n    xvfadd.s res1, VX3, res1\n",
        .count = 2,
    },
    // with beta zero and y strided, the loops advance only `YY` but the tail
    // resumes from `Y`, overwriting y from its start
    .{ .source = "axpby_lsx.S", .find = "    blt $r0, I, .L122\n    b .L997\n", .replace = "    blt $r0, I, .L122\n    move Y, YY\n    b .L997\n" },
    .{ .source = "axpby_lsx.S", .find = "    blt $r0, I, .L124\n    b .L997\n", .replace = "    blt $r0, I, .L124\n    move Y, YY\n    b .L997\n" },
    .{ .source = "axpby_lsx.S", .find = "    blt $r0, I, .L224\n    b .L997\n", .replace = "    blt $r0, I, .L224\n    move Y, YY\n    b .L997\n" },
    // double's `.L121` counted each block twice, skipping every other one,
    // and its `.L122` never counted at all
    .{
        .source = "axpby_lasx.S",
        .find = "    addi.d  I, I, -1\n    blt $r0, I, .L121\n",
        .replace = "#ifndef DOUBLE\n    addi.d  I, I, -1\n#endif\n    blt $r0, I, .L121\n",
    },
    .{
        .source = "axpby_lasx.S",
        .find = ".L122: // ALPHA!=0 BETA==0\n    xvld VX0, X, 0 * SIZE\n#ifdef DOUBLE\n",
        .replace = ".L122: // ALPHA!=0 BETA==0\n    xvld VX0, X, 0 * SIZE\n#ifdef DOUBLE\n    addi.d  I, I, -1\n",
    },
    // as in axpby_lsx, with beta zero and x and y strided
    .{ .source = "caxpby_lsx.S", .find = "    blt $r0, I, .L222\n    b .L997\n", .replace = "    blt $r0, I, .L222\n    move Y, YY\n    b .L997\n", .count = 2 },
    // with alpha and beta zero and y strided, stored each zero twice at the
    // real part, never the imaginary; lasx single also stored one to `YY`,
    // which still points at y's start
    .{ .source = "caxpby_lsx.S", .find = "    vstelm.d VXZ, Y, 0, 0\n    vstelm.d VXZ, Y, 0, 0\n", .replace = "    vstelm.d VXZ, Y, 0 * SIZE, 0\n    vstelm.d VXZ, Y, 1 * SIZE, 0\n", .count = 4 },
    .{ .source = "caxpby_lsx.S", .find = "    vstelm.w VXZ, Y, 0, 0\n    vstelm.w VXZ, Y, 0, 0\n", .replace = "    vstelm.w VXZ, Y, 0 * SIZE, 0\n    vstelm.w VXZ, Y, 1 * SIZE, 0\n", .count = 4 },
    .{ .source = "caxpby_lasx.S", .find = "    xvstelm.w VXZ, YY, 0, 0\n    xvstelm.w VXZ, YY, 0, 0\n", .replace = "    xvstelm.w VXZ, Y, 0, 0\n    xvstelm.w VXZ, Y, 0, 0\n" },
    .{ .source = "caxpby_lasx.S", .find = "    xvstelm.d VXZ, Y, 0, 0\n    xvstelm.d VXZ, Y, 0, 0\n", .replace = "    xvstelm.d VXZ, Y, 0 * SIZE, 0\n    xvstelm.d VXZ, Y, 1 * SIZE, 0\n", .count = 4 },
    .{ .source = "caxpby_lasx.S", .find = "    xvstelm.w VXZ, Y, 0, 0\n    xvstelm.w VXZ, Y, 0, 0\n", .replace = "    xvstelm.w VXZ, Y, 0 * SIZE, 0\n    xvstelm.w VXZ, Y, 1 * SIZE, 0\n", .count = 8 },
    // strided x reduced with `vfmaxa`, by magnitude
    .{ .source = "max_lsx.S", .find = "vfmaxa.", .replace = "vfmax.", .count = 6 },
    .{ .source = "min_lsx.S", .find = "vfmaxa.", .replace = "vfmin.", .count = 6 },
    // strided x started the minimum at `|re(x0)|` rather than `|re|+|im|`
    .{ .source = "camin_lsx.S", .find = "    ADD s1, a1, a0\n", .replace = "    ADD s1, a1, a0\n    MOV a0, s1\n" },
    .{ .source = "camin_lasx.S", .find = "    ADD s1, a1, a0\n", .replace = "    ADD s1, a1, a0\n    MOV a0, s1\n" },
    // double, strided x: added 8 to the indices twice per block of eight
    .{ .source = "imax_lasx.S", .find = "    xvinsgr2vr.d VX0, t4, 3\n    xvadd.d VI1, VI1, VINC8\n", .replace = "    xvinsgr2vr.d VX0, t4, 3\n" },
    .{ .source = "imin_lasx.S", .find = "    xvinsgr2vr.d VX0, t4, 3\n    xvadd.d VI1, VI1, VINC8\n", .replace = "    xvinsgr2vr.d VX0, t4, 3\n" },
    // double, strided x: kept the new minimum by the previous comparison's
    // mask, before making this one
    .{
        .source = "iamin_lsx.S",
        .find = "    vfmina.d VM1, VM0, x3\n    vbitsel.v VM0, VM1, VM0, VT0\n    vfcmp.ceq.d VT0, VM0, VM1\n",
        .replace = "    vfmina.d VM1, VM0, x3\n    vfcmp.ceq.d VT0, VM0, VM1\n    vbitsel.v VM0, VM1, VM0, VT0\n",
    },
    // strided x: after the loop, skipped the reduction across lanes (`.L15`)
    .{ .source = "iamin_lasx.S", .find = "    blt $r0, I, .L24\n    .align 3\n\n.L26:", .replace = "    blt $r0, I, .L24\n    b .L15\n    .align 3\n\n.L26:" },
    // below one vector's worth, the tail started from the index the vector
    // set-up had counted to (4) rather than 1, so a tie with x0 returned it
    .{ .source = "icamin_lsx.S", .find = "bge $r0, I, .L21", .replace = "bge $r0, I, .L21_first", .count = 4 },
    .{ .source = "icamin_lsx.S", .find = ".L21: //N<4\n", .replace = "    b .L21\n.L21_first:\n    li.d i0, 1\n.L21: //N<4\n" },
    // single, strided x: started the minimum from `|x|` taken by selecting
    // the negations `x3`/`x4`, never computed here
    .{
        .source = "icamin_lsx.S",
        .find = "    vfcmp.clt.s VT0, x1, VI3\n    vfcmp.clt.s VINC8, x2, VI3\n",
        .replace = "    vfmul.s x3, VI4, x1\n    vfmul.s x4, VI4, x2\n    vfcmp.clt.s VT0, x1, VI3\n    vfcmp.clt.s VINC8, x2, VI3\n",
    },
    // strided x: the loop used `VINC8` as scratch, but single steps the
    // indices by it; `VM1` is free until the reduction
    .{
        .source = "icamin_lasx.S",
        .find = "    XVCMPLT VINC8, x2, VI3\n    xvbitsel.v x1, x1, x3, VT0\n    xvbitsel.v x2, x2, x4, VINC8\n",
        .replace = "    XVCMPLT VM1, x2, VI3\n    xvbitsel.v x1, x1, x3, VT0\n    xvbitsel.v x2, x2, x4, VM1\n",
    },
    // single, lanes 4 to 7: kept their pairwise minimum in `x1`, so the tie
    // check after read it as lane 4 and took that lane's index; keep it in
    // `VX1`, free by then (`VM1`, as `.L15` uses, holds the first half's index)
    .{
        .source = "icamin_lasx.S",
        .find = "    xvbitsel.v x1, x1, x2, VT0\n    xvbitsel.v VINC4, VI1, VI2, VT0\n    xvfcmp.clt.s VT0, x4, x3\n    xvbitsel.v VM0, x3, x4, VT0\n    xvbitsel.v VINC8, VI3, VI4, VT0\n    xvfcmp.clt.s VT0, x1, VM0\n    xvbitsel.v VM0, VM0, x1, VT0\n",
        .replace = "    xvbitsel.v VX1, x1, x2, VT0\n    xvbitsel.v VINC4, VI1, VI2, VT0\n    xvfcmp.clt.s VT0, x4, x3\n    xvbitsel.v VM0, x3, x4, VT0\n    xvbitsel.v VINC8, VI3, VI4, VT0\n    xvfcmp.clt.s VT0, VX1, VM0\n    xvbitsel.v VM0, VM0, VX1, VT0\n",
    },
    // single, strided x below 8 compared against a first element never loaded
    .{
        .source = "imax_lasx.S",
        .find = "    addi.w i0, i0, 1\n    ld.w t1, TEMP, 0 * SIZE\n    add.d TEMP, TEMP, INCX\n    srai.d I, N, 3\n",
        .replace = "    addi.w i0, i0, 1\n    ld.w t1, TEMP, 0 * SIZE\n    add.d TEMP, TEMP, INCX\n    xvinsgr2vr.w VM0, t1, 0\n    srai.d I, N, 3\n",
    },
    .{
        .source = "imin_lasx.S",
        .find = "    addi.w i0, i0, 1\n    ld.w t1, TEMP, 0 * SIZE\n    add.d TEMP, TEMP, INCX\n    srai.d I, N, 3\n",
        .replace = "    addi.w i0, i0, 1\n    ld.w t1, TEMP, 0 * SIZE\n    add.d TEMP, TEMP, INCX\n    xvinsgr2vr.w VM0, t1, 0\n    srai.d I, N, 3\n",
    },
};

pub const loongarch64: arch.Arch = .{
    .name = "loongarch64",
    .triple = "loongarch64-linux-gnu",
    .config_h = "#define ARCH_LOONGARCH64 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_LOONGARCH64", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .core_defines = coreDefines,
    .patches = &patches,
};
