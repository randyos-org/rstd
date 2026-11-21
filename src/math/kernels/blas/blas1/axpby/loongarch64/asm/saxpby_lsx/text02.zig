//! Section 2 of 2 of the text of `saxpby`, which `__root__.zig` joins, from OpenBLAS `kernel/loongarch64/axpby_lsx.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .align 3
    \\ 20016:
    \\ ld.w $r14, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r16, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r15, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r17, $r5, 0 * 4
    \\ vinsgr2vr.w $vr8, $r14, 0
    \\ vinsgr2vr.w $vr8, $r16, 1
    \\ vinsgr2vr.w $vr8, $r15, 2
    \\ vinsgr2vr.w $vr8, $r17, 3
    \\ add.d $r5, $r5, $r6
    \\ vfmul.s $vr8, $vr23, $vr8
    \\ ld.w $r14, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r16, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r15, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r17, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.w $vr20, $r14, 0
    \\ vinsgr2vr.w $vr20, $r16, 1
    \\ vinsgr2vr.w $vr20, $r15, 2
    \\ vinsgr2vr.w $vr20, $r17, 3
    \\ vst $vr8, $r7, 0 * 4
    \\ vfmul.s $vr20, $vr20, $vr23
    \\ addi.d $r12, $r12, -1
    \\ vst $vr20, $r7, 4 * 4
    \\ addi.d $r7, $r7, 8 * 4
    \\ blt $r0, $r12, 20016b
    \\ b 20025f
    \\ .align 3
    \\ 20017:
    \\ vld $vr21, $r7, 0 * 4
    \\ vld $vr22, $r7, 4 * 4
    \\ vfmul.s $vr21, $vr21, $vr9
    \\ vfmul.s $vr22, $vr22, $vr9
    \\ vst $vr21, $r7, 0 * 4
    \\ vst $vr22, $r7, 4 * 4
    \\ addi.d $r7, $r7, 8 * 4
    \\ addi.d $r12, $r12, -1
    \\ blt $r0, $r12, 20017b
    \\ b 20025f
    \\ .align 3
    \\ 20018:
    \\ vst $vr19, $r7, 0 * 4
    \\ vst $vr19, $r7, 4 * 4
    \\ addi.d $r7, $r7, 8 * 4
    \\ addi.d $r12, $r12, -1
    \\ blt $r0, $r12, 20018b
    \\ b 20025f
    \\ .align 3
    \\ 20019:
    \\ bge $r0, $r12, 20025f
    \\ move $r19, $r7
    \\ fcmp.ceq.s $fcc0, $f0, $f12
    \\ bcnez $fcc0, 20020f
    \\ fcmp.ceq.s $fcc0, $f1, $f12
    \\ bcnez $fcc0, 20022f
    \\ b 20021f
    \\ .align 3
    \\ 20020:
    \\ fcmp.ceq.s $fcc0, $f1, $f12
    \\ bcnez $fcc0, 20024f
    \\ b 20023f
    \\ .align 3
    \\ 20021:
    \\ ld.w $r14, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r16, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r15, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r17, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.w $vr8, $r14, 0
    \\ vinsgr2vr.w $vr8, $r16, 1
    \\ vinsgr2vr.w $vr8, $r15, 2
    \\ vinsgr2vr.w $vr8, $r17, 3
    \\ ld.w $r14, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.w $r16, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.w $r15, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.w $r17, $r7, 0 * 4
    \\ vinsgr2vr.w $vr21, $r14, 0
    \\ vinsgr2vr.w $vr21, $r16, 1
    \\ vinsgr2vr.w $vr21, $r15, 2
    \\ vinsgr2vr.w $vr21, $r17, 3
    \\ add.d $r7, $r7, $r8
    \\ vfmul.s $vr8, $vr8, $vr23
    \\ ld.w $r14, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r16, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r15, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r17, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ vfmadd.s $vr21, $vr21, $vr9, $vr8
    \\ vinsgr2vr.w $vr20, $r14, 0
    \\ vinsgr2vr.w $vr20, $r16, 1
    \\ vinsgr2vr.w $vr20, $r15, 2
    \\ vinsgr2vr.w $vr20, $r17, 3
    \\ vstelm.w $vr21, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr21, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr21, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr21, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ ld.w $r14, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.w $r16, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.w $r15, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.w $r17, $r7, 0 * 4
    \\ vinsgr2vr.w $vr22, $r14, 0
    \\ vinsgr2vr.w $vr22, $r16, 1
    \\ vinsgr2vr.w $vr22, $r15, 2
    \\ vinsgr2vr.w $vr22, $r17, 3
    \\ add.d $r7, $r7, $r8
    \\ vfmul.s $vr20, $vr20, $vr23
    \\ addi.d $r12, $r12, -1
    \\ vfmadd.s $vr22, $vr22, $vr9, $vr20
    \\ vstelm.w $vr22, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr22, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr22, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr22, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ blt $r0, $r12, 20021b
    \\ b 20025f
    \\ .align 3
    \\ 20022:
    \\ ld.w $r14, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r16, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r15, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r17, $r5, 0 * 4
    \\ vinsgr2vr.w $vr8, $r14, 0
    \\ vinsgr2vr.w $vr8, $r16, 1
    \\ vinsgr2vr.w $vr8, $r15, 2
    \\ vinsgr2vr.w $vr8, $r17, 3
    \\ add.d $r5, $r5, $r6
    \\ vfmul.s $vr8, $vr8, $vr23
    \\ ld.w $r14, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r16, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r15, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.w $r17, $r5, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.w $vr20, $r14, 0
    \\ vinsgr2vr.w $vr20, $r16, 1
    \\ vinsgr2vr.w $vr20, $r15, 2
    \\ vinsgr2vr.w $vr20, $r17, 3
    \\ vstelm.w $vr8, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr8, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr8, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr8, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ vfmul.s $vr20, $vr20, $vr23
    \\ addi.d $r12, $r12, -1
    \\ vstelm.w $vr20, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr20, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr20, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr20, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ blt $r0, $r12, 20022b
    \\ move $r7, $r19
    \\ b 20025f
    \\ .align 3
    \\ 20023:
    \\ ld.w $r14, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.w $r16, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.w $r15, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.w $r17, $r7, 0 * 4
    \\ vinsgr2vr.w $vr21, $r14, 0
    \\ vinsgr2vr.w $vr21, $r16, 1
    \\ vinsgr2vr.w $vr21, $r15, 2
    \\ vinsgr2vr.w $vr21, $r17, 3
    \\ add.d $r7, $r7, $r8
    \\ vfmul.s $vr21, $vr21, $vr9
    \\ ld.w $r14, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.w $r16, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.w $r15, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.w $r17, $r7, 0 * 4
    \\ add.d $r7, $r7, $r8
    \\ vinsgr2vr.w $vr22, $r14, 0
    \\ vinsgr2vr.w $vr22, $r16, 1
    \\ vinsgr2vr.w $vr22, $r15, 2
    \\ vinsgr2vr.w $vr22, $r17, 3
    \\ vstelm.w $vr21, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr21, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr21, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr21, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ vfmul.s $vr22, $vr22, $vr9
    \\ addi.d $r12, $r12, -1
    \\ vstelm.w $vr22, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr22, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr22, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr22, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ blt $r0, $r12, 20023b
    \\ b 20025f
    \\ .align 3
    \\ 20024:
    \\ vstelm.w $vr19, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr19, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr19, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr19, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr19, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr19, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr19, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ vstelm.w $vr19, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ addi.d $r12, $r12, -1
    \\ blt $r0, $r12, 20024b
    \\ move $r7, $r19
    \\ b 20025f
    \\ .align 3
    \\ 20025:
    \\ andi $r12, $r4, 7
    \\ bge $r0, $r12, 20030f
    \\ .align 3
    \\ 20026:
    \\ movgr2fr.d $f12, $r0
    \\ fcmp.ceq.s $fcc0, $f0, $f12
    \\ fcmp.ceq.s $fcc1, $f1, $f12
    \\ addi.d $r12, $r12, -1
    \\ bceqz $fcc0, 20027f
    \\ bcnez $fcc1, 20029f
    \\ fneg.s $f12, $f12
    \\ b 20028f
    \\ 20027:
    \\ fld.s $f12, $r5, 0 * 4
    \\ fmul.s $f12, $f12, $f0
    \\ bcnez $fcc1, 20029f
    \\ 20028:
    \\ fld.s $f13, $r7, 0 * 4
    \\ fmadd.s $f12, $f13, $f1, $f12
    \\ 20029:
    \\ fst.s $f12, $r7, 0 * 4
    \\ add.d $r5, $r5, $r6
    \\ add.d $r7, $r7, $r8
    \\ blt $r0, $r12, 20026b
    \\ .align 3
    \\ 20030:
    \\ move $r4, $r12
    \\ jirl $r0, $r1, 0x0
    \\ .align 3
    \\
;
