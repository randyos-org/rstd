//! Section 2 of 2 of the text of `daxpby`, which `__root__.zig` joins, from OpenBLAS `kernel/loongarch64/axpby_lsx.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .align 3
    \\ 20016:
    \\ ld.d $r14, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r16, $r5, 0 * 8
    \\ vinsgr2vr.d $vr8, $r14, 0
    \\ vinsgr2vr.d $vr8, $r16, 1
    \\ add.d $r5, $r5, $r6
    \\ vfmul.d $vr8, $vr23, $vr8
    \\ ld.d $r15, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r17, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.d $vr20, $r15, 0
    \\ vinsgr2vr.d $vr20, $r17, 1
    \\ vst $vr8, $r7, 0 * 8
    \\ vfmul.d $vr20, $vr23, $vr20
    \\ ld.d $r14, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r16, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.d $vr8, $r14, 0
    \\ vinsgr2vr.d $vr8, $r16, 1
    \\ vst $vr20, $r7, 2 * 8
    \\ vfmul.d $vr8, $vr8, $vr23
    \\ ld.d $r15, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r17, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.d $vr20, $r15, 0
    \\ vinsgr2vr.d $vr20, $r17, 1
    \\ vst $vr8, $r7, 4 * 8
    \\ vfmul.d $vr20, $vr20, $vr23
    \\ addi.d $r12, $r12, -1
    \\ vst $vr20, $r7, 6 * 8
    \\ addi.d $r7, $r7, 8 * 8
    \\ blt $r0, $r12, 20016b
    \\ b 20025f
    \\ .align 3
    \\ 20017:
    \\ vld $vr21, $r7, 0 * 8
    \\ vld $vr22, $r7, 2 * 8
    \\ vfmul.d $vr21, $vr21, $vr9
    \\ vfmul.d $vr22, $vr22, $vr9
    \\ vst $vr21, $r7, 0 * 8
    \\ vst $vr22, $r7, 2 * 8
    \\ vld $vr21, $r7, 4 * 8
    \\ vld $vr22, $r7, 6 * 8
    \\ vfmul.d $vr21, $vr21, $vr9
    \\ vfmul.d $vr22, $vr22, $vr9
    \\ vst $vr21, $r7, 4 * 8
    \\ vst $vr22, $r7, 6 * 8
    \\ addi.d $r7, $r7, 8 * 8
    \\ addi.d $r12, $r12, -1
    \\ blt $r0, $r12, 20017b
    \\ b 20025f
    \\ .align 3
    \\ 20018:
    \\ vst $vr19, $r7, 0 * 8
    \\ vst $vr19, $r7, 2 * 8
    \\ vst $vr19, $r7, 4 * 8
    \\ vst $vr19, $r7, 6 * 8
    \\ addi.d $r7, $r7, 8 * 8
    \\ addi.d $r12, $r12, -1
    \\ blt $r0, $r12, 20018b
    \\ b 20025f
    \\ .align 3
    \\ 20019:
    \\ bge $r0, $r12, 20025f
    \\ move $r19, $r7
    \\ fcmp.ceq.d $fcc0, $f0, $f12
    \\ bcnez $fcc0, 20020f
    \\ fcmp.ceq.d $fcc0, $f1, $f12
    \\ bcnez $fcc0, 20022f
    \\ b 20021f
    \\ .align 3
    \\ 20020:
    \\ fcmp.ceq.d $fcc0, $f1, $f12
    \\ bcnez $fcc0, 20024f
    \\ b 20023f
    \\ .align 3
    \\ 20021:
    \\ ld.d $r14, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r16, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.d $vr8, $r14, 0
    \\ vinsgr2vr.d $vr8, $r16, 1
    \\ ld.d $r15, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r17, $r7, 0 * 8
    \\ vinsgr2vr.d $vr21, $r15, 0
    \\ vinsgr2vr.d $vr21, $r17, 1
    \\ add.d $r7, $r7, $r8
    \\ vfmul.d $vr8, $vr8, $vr23
    \\ vfmadd.d $vr21, $vr21, $vr9, $vr8
    \\ ld.d $r14, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r16, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.d $vr20, $r14, 0
    \\ vinsgr2vr.d $vr20, $r16, 1
    \\ vstelm.d $vr21, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr21, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ ld.d $r15, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r17, $r7, 0 * 8
    \\ vinsgr2vr.d $vr22, $r15, 0
    \\ vinsgr2vr.d $vr22, $r17, 1
    \\ add.d $r7, $r7, $r8
    \\ vfmul.d $vr20, $vr20, $vr23
    \\ vfmadd.d $vr22, $vr22, $vr9, $vr20
    \\ ld.d $r14, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r16, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.d $vr8, $r14, 0
    \\ vinsgr2vr.d $vr8, $r16, 1
    \\ vstelm.d $vr22, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr22, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ ld.d $r15, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r17, $r7, 0 * 8
    \\ vinsgr2vr.d $vr21, $r15, 0
    \\ vinsgr2vr.d $vr21, $r17, 1
    \\ add.d $r7, $r7, $r8
    \\ vfmul.d $vr8, $vr8, $vr23
    \\ vfmadd.d $vr21, $vr21, $vr9, $vr8
    \\ ld.d $r15, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r17, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.d $vr20, $r15, 0
    \\ vinsgr2vr.d $vr20, $r17, 1
    \\ vstelm.d $vr21, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr21, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ ld.d $r14, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r16, $r7, 0 * 8
    \\ vinsgr2vr.d $vr22, $r14, 0
    \\ vinsgr2vr.d $vr22, $r16, 1
    \\ add.d $r7, $r7, $r8
    \\ vfmul.d $vr20, $vr20, $vr23
    \\ vfmadd.d $vr22, $vr22, $vr9, $vr20
    \\ addi.d $r12, $r12, -1
    \\ vstelm.d $vr22, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr22, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ blt $r0, $r12, 20021b
    \\ b 20025f
    \\ .align 3
    \\ 20022:
    \\ ld.d $r14, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r16, $r5, 0 * 8
    \\ vinsgr2vr.d $vr8, $r14, 0
    \\ vinsgr2vr.d $vr8, $r16, 1
    \\ add.d $r5, $r5, $r6
    \\ vfmul.d $vr8, $vr8, $vr23
    \\ ld.d $r15, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r17, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.d $vr20, $r15, 0
    \\ vinsgr2vr.d $vr20, $r17, 1
    \\ vstelm.d $vr8, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr8, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vfmul.d $vr20, $vr20, $vr23
    \\ ld.d $r14, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r16, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.d $vr8, $r14, 0
    \\ vinsgr2vr.d $vr8, $r16, 1
    \\ vstelm.d $vr20, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr20, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vfmul.d $vr8, $vr8, $vr23
    \\ ld.d $r15, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r17, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ vinsgr2vr.d $vr20, $r15, 0
    \\ vinsgr2vr.d $vr20, $r17, 1
    \\ vstelm.d $vr8, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr8, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vfmul.d $vr20, $vr20, $vr23
    \\ addi.d $r12, $r12, -1
    \\ vstelm.d $vr20, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr20, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ blt $r0, $r12, 20022b
    \\ move $r7, $r19
    \\ b 20025f
    \\ .align 3
    \\ 20023:
    \\ ld.d $r14, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r16, $r7, 0 * 8
    \\ vinsgr2vr.d $vr21, $r14, 0
    \\ vinsgr2vr.d $vr21, $r16, 1
    \\ add.d $r7, $r7, $r8
    \\ vfmul.d $vr21, $vr21, $vr9
    \\ ld.d $r15, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r17, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ vinsgr2vr.d $vr22, $r15, 0
    \\ vinsgr2vr.d $vr22, $r17, 1
    \\ vstelm.d $vr21, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr21, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vfmul.d $vr22, $vr22, $vr9
    \\ ld.d $r14, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r16, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ vinsgr2vr.d $vr21, $r14, 0
    \\ vinsgr2vr.d $vr21, $r16, 1
    \\ vstelm.d $vr22, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr22, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vfmul.d $vr21, $vr21, $vr9
    \\ ld.d $r15, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r17, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ vinsgr2vr.d $vr22, $r15, 0
    \\ vinsgr2vr.d $vr22, $r17, 1
    \\ vstelm.d $vr21, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr21, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vfmul.d $vr22, $vr22, $vr9
    \\ addi.d $r12, $r12, -1
    \\ vstelm.d $vr22, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr22, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ blt $r0, $r12, 20023b
    \\ b 20025f
    \\ .align 3
    \\ 20024:
    \\ vstelm.d $vr19, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr19, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr19, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr19, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr19, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr19, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr19, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ vstelm.d $vr19, $r19, 0, 1
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
    \\ fcmp.ceq.d $fcc0, $f0, $f12
    \\ fcmp.ceq.d $fcc1, $f1, $f12
    \\ addi.d $r12, $r12, -1
    \\ bceqz $fcc0, 20027f
    \\ bcnez $fcc1, 20029f
    \\ fneg.d $f12, $f12
    \\ b 20028f
    \\ 20027:
    \\ fld.d $f12, $r5, 0 * 8
    \\ fmul.d $f12, $f12, $f0
    \\ bcnez $fcc1, 20029f
    \\ 20028:
    \\ fld.d $f13, $r7, 0 * 8
    \\ fmadd.d $f12, $f13, $f1, $f12
    \\ 20029:
    \\ fst.d $f12, $r7, 0 * 8
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
