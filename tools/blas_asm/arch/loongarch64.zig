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
    out[0] = try a.print("-D{s}", .{core});
    return out;
}

/// Fixes for bugs in the LSX/LASX sources, still in OpenBLAS v0.3.34-317;
/// the qemu rig's per-core tests found each one.
const patches = [_]arch.Patch{
    // `axpy` returned at once for a zero alpha, so `0 * x` never reached `y`
    // and NaN and infinities in `x` were lost, unlike on every other core
    .{ .source = "axpy_lsx.S", .find = "    CMPEQ   $fcc0, ALPHA, a1\n    bcnez $fcc0, .L999\n", .replace = "" },
    .{ .source = "axpy_lasx.S", .find = "    CMPEQ   $fcc0, ALPHA, a1\n    bcnez $fcc0, .L999\n", .replace = "" },
    .{ .source = "caxpy_lsx.S", .find = "    bceqz $fcc0, .L10\n    bcnez $fcc1, .L999\n", .replace = "" },
    .{ .source = "caxpy_lasx.S", .find = "    bceqz $fcc0, .L10\n    bcnez $fcc1, .L999\n", .replace = "" },
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
    // double, x strided and y not: the second pair of each block of four
    // stored its upper half over both halves
    .{ .source = "caxpby_lsx.S", .find = "    vst VX3, Y, 4 * SIZE\n    vst VX3, Y, 6 * SIZE\n", .replace = "    vst VX2, Y, 4 * SIZE\n    vst VX3, Y, 6 * SIZE\n" },
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
    // the tail after the blocks read both operands whatever the scalars, so
    // a NaN or infinity met a zero scalar's 0 (`kernel/arm/axpby.c` reads
    // neither), from elements the zero-scalar blocks left their pointer at
    .{ .source = "axpby_lsx.S", .find = axpby_tail, .replace = axpby_tail_fixed },
    .{ .source = "axpby_lasx.S", .find = axpby_tail, .replace = axpby_tail_fixed },
    .{ .source = "caxpby_lsx.S", .find = caxpby_lsx_tail, .replace = caxpby_tail_fixed },
    .{ .source = "caxpby_lasx.S", .find = caxpby_lasx_tail, .replace = caxpby_tail_fixed },
    // double compared the scalars with a single-precision 0 (`ffint.s.l`),
    // whose upper half is undefined, so it never took the zero-scalar blocks
    .{ .source = "axpby_lsx.S", .find = axpby_zero, .replace = axpby_zero_fixed },
    .{ .source = "axpby_lasx.S", .find = axpby_zero, .replace = axpby_zero_fixed },
    // with one stride 1 and the other not, there are no zero-scalar blocks:
    // send a zero scalar to the strided ones (`.L22`), which have them
    .{ .source = "caxpby_lsx.S", .find = caxpby_l12, .replace = caxpby_l12 ++ zeroToStrided("L12") },
    .{ .source = "caxpby_lasx.S", .find = caxpby_l12, .replace = caxpby_l12 ++ zeroToStrided("L12") },
    .{ .source = "caxpby_lsx.S", .find = caxpby_l21, .replace = caxpby_l21 ++ zeroToStrided("L21") },
    .{ .source = "caxpby_lasx.S", .find = caxpby_l21, .replace = caxpby_l21 ++ zeroToStrided("L21") },
    // negative stride: `?amax_k` returns 0 for `incx <= 0`, so the scale
    // came out 0 and so did the norm; walk x forwards from its lowest element
    .{ .source = "snrm2_lsx.S", .find = "beq $r0, INCX, .L999\n", .replace = nrm2Forwards("BASE_SHIFT") },
    .{ .source = "snrm2_lasx.S", .find = "beq $r0, INCX, .L999\n", .replace = nrm2Forwards("BASE_SHIFT") },
    .{ .source = "cnrm2_lsx.S", .find = "beq $r0, INCX, .L999\n", .replace = nrm2Forwards("ZBASE_SHIFT") },
    .{ .source = "cnrm2_lasx.S", .find = "beq $r0, INCX, .L999\n", .replace = nrm2Forwards("ZBASE_SHIFT") },
    // with a zero stride, the strided loops load a block before storing it,
    // so `n` operations on one element act as fewer; zero strides take the
    // one-at-a-time loop, as `crot_lsx.S` already does
    .{ .source = "rot_lsx.S", .find = rot_blocks, .replace = rot_zero ++ rot_blocks },
    .{ .source = "rot_lasx.S", .find = rot_blocks, .replace = rot_zero ++ rot_blocks },
    .{ .source = "swap.S", .find = swap_strided, .replace = ".L20:\n   move I, N\n   beq INCX, $r0, .L26\n   beq INCY, $r0, .L26\n" ++ swap_strided[6..] },
    .{ .source = "cswap_lsx.S", .find = cswap_strided, .replace = cswap_zero ++ cswap_strided[6..] },
    .{ .source = "cswap_lasx.S", .find = cswap_strided, .replace = cswap_zero ++ cswap_strided[6..] },
    .{ .source = "cswap_lsx.S", .find = cswap_tail, .replace = cswap_tail ++ cswap_one },
    .{ .source = "cswap_lasx.S", .find = cswap_tail, .replace = cswap_tail ++ cswap_one },
    // and with a zero `incy`, `axpy` (`y` loaded a block at a time from `.L12`
    // and `.L22`)
    zeroIncy("axpy_lsx.S", ".L224"),
    zeroIncy("axpy_lasx.S", ".L224"),
    zeroIncy("caxpy_lsx.S", ".L998"),
    zeroIncy("caxpy_lasx.S", ".L998"),
};

const cswap_strided = ".L22:\n    bge $r0, I, .L223\n";
const cswap_zero = ".L22:\n    beq INCX, $r0, .L22_one\n    beq INCY, $r0, .L22_one\n";
const cswap_tail = ".L223:\n    andi I, N, 3\n    bge $r0, I, .L999\n";
const cswap_one = "    b .L224\n.L22_one:\n    move I, N\n";

/// A zero `incy` takes the one-at-a-time loop `tail`, before the dispatch on
/// unit strides.
fn zeroIncy(comptime source: []const u8, comptime tail: []const u8) arch.Patch {
    const dispatch = "    bne INCX, TEMP, .L20\n    bne INCY, TEMP, .L12 // INCX==1 and INCY!=1\n";
    return .{ .source = source, .find = dispatch, .replace = "    bne INCY, $r0, .L_incy\n    move I, N\n    b " ++ tail ++ "\n.L_incy:\n" ++ dispatch };
}

const rot_blocks = "    srai.d I, N, 3\n    bge $r0, I, .L997\n    bne INCX, TEMP, .L20\n";
const rot_zero = "    move I, N\n    beq INCX, $r0, .L998\n    beq INCY, $r0, .L998\n";
const swap_strided = ".L20:\n   srai.d  I, N, 3\n   move    XX, X\n";

const axpby_zero = "    movgr2fr.d a1, $r0\n    ffint.s.l a1, a1\n";
const axpby_zero_fixed = "    movgr2fr.d a1, $r0\n#ifdef DOUBLE\n    ffint.d.l a1, a1\n#else\n    ffint.s.l a1, a1\n#endif\n";

const caxpby_l12 = ".L12: // INCX==1 and INCY!=1\n    bge $r0, I, .L997\n";
const caxpby_l21 = ".L21:// INCX!=1 and INCY==1\n    bge $r0, I, .L997\n";

/// To `.L22` if `alpha` or `beta` is zero.
fn zeroToStrided(comptime tag: []const u8) []const u8 {
    return "    CMPEQ $fcc0, BETAR, a1\n    CMPEQ $fcc1, BETAI, a1\n    CMPEQ $fcc2, ALPHAR, a1\n    CMPEQ $fcc3, ALPHAI, a1\n" ++
        "    bceqz $fcc0, ." ++ tag ++ "_a\n    bcnez $fcc1, .L22\n." ++ tag ++ "_a:\n" ++
        "    bceqz $fcc2, ." ++ tag ++ "_go\n    bcnez $fcc3, .L22\n." ++ tag ++ "_go:\n";
}

const axpby_tail =
    \\.L998:
    \\    LD   $f12, X, 0 * SIZE
    \\    LD   $f13, Y, 0 * SIZE
    \\    addi.d I, I, -1
    \\    MUL  $f12, $f12, ALPHA
    \\    MADD $f13, $f13, BETA, $f12
    \\    ST   $f13, Y, 0 * SIZE
    \\
;

/// Reads `x` only for a nonzero alpha and `y` for a nonzero beta; a zero
/// alpha adds -0, keeping `beta*y` exact, and both zero store +0.
const axpby_tail_fixed =
    \\.L998:
    \\    movgr2fr.d $f12, $r0
    \\    CMPEQ $fcc0, ALPHA, $f12
    \\    CMPEQ $fcc1, BETA, $f12
    \\    addi.d I, I, -1
    \\    bceqz $fcc0, .L998_ax
    \\    bcnez $fcc1, .L998_st
    \\    NEG  $f12, $f12
    \\    b .L998_by
    \\.L998_ax:
    \\    LD   $f12, X, 0 * SIZE
    \\    MUL  $f12, $f12, ALPHA
    \\    bcnez $fcc1, .L998_st
    \\.L998_by:
    \\    LD   $f13, Y, 0 * SIZE
    \\    MADD $f12, $f13, BETA, $f12
    \\.L998_st:
    \\    ST   $f12, Y, 0 * SIZE
    \\
;

const caxpby_lsx_tail =
    \\.L998:
    \\#ifdef DOUBLE
    \\    fld.d a1, X, 0 * SIZE
    \\    fld.d a2, X, 1 * SIZE
    \\    fld.d a3, Y, 0 * SIZE
    \\    fld.d a4, Y, 1 * SIZE
    \\    addi.d I, I, -1
    \\    fmul.d s1, ALPHAI, a2
    \\    fmul.d s2, ALPHAI, a1
    \\    fmul.d s3, BETAI, a4
    \\    fmul.d s4, BETAI, a3
    \\    fmsub.d s1, ALPHAR, a1, s1
    \\    fmadd.d s2, a2, ALPHAR, s2
    \\    fmsub.d s3, BETAR, a3, s3
    \\    fmadd.d s4, a4, BETAR, s4
    \\    fadd.d s3, s3, s1
    \\    fadd.d s4, s4, s2
    \\    fst.d s3, Y, 0 * SIZE
    \\    fst.d s4, Y, 1 * SIZE
    \\    add.d X, X, INCX
    \\    add.d Y, Y, INCY
    \\    blt $r0, I, .L998
    \\    .align 3
    \\#else
    \\    fld.s a1, X, 0 * SIZE
    \\    fld.s a2, X, 1 * SIZE
    \\    fld.s a3, Y, 0 * SIZE
    \\    fld.s a4, Y, 1 * SIZE
    \\    addi.d I, I, -1
    \\    fmul.s s1, ALPHAI, a2
    \\    fmul.s s2, ALPHAI, a1
    \\    fmul.s s3, BETAI, a4
    \\    fmul.s s4, BETAI, a3
    \\    fmsub.s s1, ALPHAR, a1, s1
    \\    fmadd.s s2, a2, ALPHAR, s2
    \\    fmsub.s s3, BETAR, a3, s3
    \\    fmadd.s s4, a4, BETAR, s4
    \\    fadd.s s3, s3, s1
    \\    fadd.s s4, s4, s2
    \\    fst.s s3, Y, 0 * SIZE
    \\    fst.s s4, Y, 1 * SIZE
    \\    add.d X, X, INCX
    \\    add.d Y, Y, INCY
    \\    blt $r0, I, .L998
    \\    .align 3
    \\#endif
    \\
;

const caxpby_lasx_tail =
    \\.L998:
    \\    LD a1, X, 0 * SIZE
    \\    LD a2, X, 1 * SIZE
    \\    LD a3, Y, 0 * SIZE
    \\    LD a4, Y, 1 * SIZE
    \\    addi.d I, I, -1
    \\    MUL s1, ALPHAI, a2
    \\    MUL s2, ALPHAI, a1
    \\    MUL s3, BETAI, a4
    \\    MUL s4, BETAI, a3
    \\    MSUB s1, ALPHAR, a1, s1
    \\    MADD s2, a2, ALPHAR, s2
    \\    MSUB s3, BETAR, a3, s3
    \\    MADD s4, a4, BETAR, s4
    \\    ADD s3, s3, s1
    \\    ADD s4, s4, s2
    \\    ST s3, Y, 0 * SIZE
    \\    ST s4, Y, 1 * SIZE
    \\    add.d X, X, INCX
    \\    add.d Y, Y, INCY
    \\    blt $r0, I, .L998
    \\    .align 3
    \\
;

/// As `axpby_tail_fixed`, complex.
const caxpby_tail_fixed =
    \\.L998:
    \\    movgr2fr.d s1, $r0
    \\    CMPEQ $fcc0, ALPHAR, s1
    \\    CMPEQ $fcc1, ALPHAI, s1
    \\    CMPEQ $fcc2, BETAR, s1
    \\    CMPEQ $fcc3, BETAI, s1
    \\    addi.d I, I, -1
    \\    bceqz $fcc0, .L998_ax
    \\    bceqz $fcc1, .L998_ax
    \\    MOV s2, s1
    \\    bceqz $fcc2, .L998_neg
    \\    bcnez $fcc3, .L998_st
    \\.L998_neg:
    \\    NEG s1, s1
    \\    NEG s2, s2
    \\    b .L998_by
    \\.L998_ax:
    \\    LD a1, X, 0 * SIZE
    \\    LD a2, X, 1 * SIZE
    \\    MUL s1, ALPHAI, a2
    \\    MUL s2, ALPHAI, a1
    \\    MSUB s1, ALPHAR, a1, s1
    \\    MADD s2, a2, ALPHAR, s2
    \\    bceqz $fcc2, .L998_by
    \\    bcnez $fcc3, .L998_st
    \\.L998_by:
    \\    LD a3, Y, 0 * SIZE
    \\    LD a4, Y, 1 * SIZE
    \\    MUL s3, BETAI, a4
    \\    MUL s4, BETAI, a3
    \\    MSUB s3, BETAR, a3, s3
    \\    MADD s4, a4, BETAR, s4
    \\    ADD s1, s3, s1
    \\    ADD s2, s4, s2
    \\.L998_st:
    \\    ST s1, Y, 0 * SIZE
    \\    ST s2, Y, 1 * SIZE
    \\    add.d X, X, INCX
    \\    add.d Y, Y, INCY
    \\    blt $r0, I, .L998
    \\    .align 3
    \\
;

/// After the zero-stride exit: for `incx < 0`, X to the lowest element and
/// `incx` negated.
fn nrm2Forwards(comptime shift: []const u8) []const u8 {
    return "beq $r0, INCX, .L999\n" ++
        "   bge INCX, $r0, .L_fwd\n" ++
        "   addi.d TEMP, N, -1\n" ++
        "   mul.d TEMP, TEMP, INCX\n" ++
        "   slli.d TEMP, TEMP, " ++ shift ++ "\n" ++
        "   add.d X, X, TEMP\n" ++
        "   sub.d INCX, $r0, INCX\n" ++
        ".L_fwd:\n";
}

pub const loongarch64: arch.Arch = .{
    .name = "loongarch64",
    .triple = "loongarch64-linux-gnu",
    // getarch's `config.h` defines `__64BIT__` for a 64-bit build; without it
    // `LDARG` is `ld.w`, which truncates a pointer passed on the stack
    .config_h = "#define ARCH_LOONGARCH64 1\n#define BINARY64 1\n#define __64BIT__ 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_LOONGARCH64", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-DOS_LINUX" },
    .core_defines = coreDefines,
    .patches = &patches,
};
