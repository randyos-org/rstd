//! Section 2 of 2 of the text of `caxpby`, which `__root__.zig` joins, from OpenBLAS `kernel/loongarch64/caxpby_lasx.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .align 3
    \\ 20020:
    \\ ld.d $r14, $r5, 0 * 4
    \\ ld.d $r16, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r15, $r5, 0 * 4
    \\ ld.d $r17, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ xvinsgr2vr.w $xr18, $r14, 0
    \\ xvinsgr2vr.w $xr17, $r16, 0
    \\ xvinsgr2vr.w $xr18, $r15, 1
    \\ xvinsgr2vr.w $xr17, $r17, 1
    \\ ld.d $r14, $r5, 0 * 4
    \\ ld.d $r16, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r15, $r5, 0 * 4
    \\ ld.d $r17, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ xvinsgr2vr.w $xr18, $r14, 2
    \\ xvinsgr2vr.w $xr17, $r16, 2
    \\ xvinsgr2vr.w $xr18, $r15, 3
    \\ xvinsgr2vr.w $xr17, $r17, 3
    \\ ld.d $r14, $r5, 0 * 4
    \\ ld.d $r16, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r15, $r5, 0 * 4
    \\ ld.d $r17, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ xvinsgr2vr.w $xr18, $r14, 4
    \\ xvinsgr2vr.w $xr17, $r16, 4
    \\ xvinsgr2vr.w $xr18, $r15, 5
    \\ xvinsgr2vr.w $xr17, $r17, 5
    \\ ld.d $r14, $r5, 0 * 4
    \\ ld.d $r16, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r15, $r5, 0 * 4
    \\ ld.d $r17, $r5, 1 * 4
    \\ xvinsgr2vr.w $xr18, $r14, 6
    \\ xvinsgr2vr.w $xr17, $r16, 6
    \\ xvinsgr2vr.w $xr18, $r15, 7
    \\ xvinsgr2vr.w $xr17, $r17, 7
    \\ add.d $r5, $r5, $r6
    \\ xvfmul.s $xr16, $xr19, $xr17
    \\ xvfmul.s $xr15, $xr19, $xr18
    \\ xvfmsub.s $xr16, $xr23, $xr18, $xr16
    \\ xvfmadd.s $xr15, $xr23, $xr17, $xr15
    \\ addi.d $r12, $r12, -1
    \\ xvstelm.w $xr16, $r19, 0 * 4, 0
    \\ xvstelm.w $xr15, $r19, 1 * 4, 0
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 1
    \\ xvstelm.w $xr15, $r19, 1 * 4, 1
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 2
    \\ xvstelm.w $xr15, $r19, 1 * 4, 2
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 3
    \\ xvstelm.w $xr15, $r19, 1 * 4, 3
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 4
    \\ xvstelm.w $xr15, $r19, 1 * 4, 4
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 5
    \\ xvstelm.w $xr15, $r19, 1 * 4, 5
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 6
    \\ xvstelm.w $xr15, $r19, 1 * 4, 6
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 7
    \\ xvstelm.w $xr15, $r19, 1 * 4, 7
    \\ add.d $r19, $r19, $r8
    \\ blt $r0, $r12, 20020b
    \\ move $r7, $r19
    \\ b 20023f
    \\ .align 3
    \\ 20021:
    \\ ld.d $r14, $r7, 0 * 4
    \\ ld.d $r16, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r15, $r7, 0 * 4
    \\ ld.d $r17, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ xvinsgr2vr.w $xr18, $r14, 0
    \\ xvinsgr2vr.w $xr17, $r16, 0
    \\ xvinsgr2vr.w $xr18, $r15, 1
    \\ xvinsgr2vr.w $xr17, $r17, 1
    \\ ld.d $r14, $r7, 0 * 4
    \\ ld.d $r16, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r15, $r7, 0 * 4
    \\ ld.d $r17, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ xvinsgr2vr.w $xr18, $r14, 2
    \\ xvinsgr2vr.w $xr17, $r16, 2
    \\ xvinsgr2vr.w $xr18, $r15, 3
    \\ xvinsgr2vr.w $xr17, $r17, 3
    \\ ld.d $r14, $r7, 0 * 4
    \\ ld.d $r16, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r15, $r7, 0 * 4
    \\ ld.d $r17, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ xvinsgr2vr.w $xr18, $r14, 4
    \\ xvinsgr2vr.w $xr17, $r16, 4
    \\ xvinsgr2vr.w $xr18, $r15, 5
    \\ xvinsgr2vr.w $xr17, $r17, 5
    \\ ld.d $r14, $r7, 0 * 4
    \\ ld.d $r16, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r15, $r7, 0 * 4
    \\ ld.d $r17, $r7, 1 * 4
    \\ xvinsgr2vr.w $xr18, $r14, 6
    \\ xvinsgr2vr.w $xr17, $r16, 6
    \\ xvinsgr2vr.w $xr18, $r15, 7
    \\ xvinsgr2vr.w $xr17, $r17, 7
    \\ add.d $r7, $r7, $r8
    \\ xvfmul.s $xr16, $xr13, $xr17
    \\ xvfmul.s $xr15, $xr13, $xr18
    \\ xvfmsub.s $xr16, $xr14, $xr18, $xr16
    \\ xvfmadd.s $xr15, $xr14, $xr17, $xr15
    \\ addi.d $r12, $r12, -1
    \\ xvstelm.w $xr16, $r19, 0 * 4, 0
    \\ xvstelm.w $xr15, $r19, 1 * 4, 0
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 1
    \\ xvstelm.w $xr15, $r19, 1 * 4, 1
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 2
    \\ xvstelm.w $xr15, $r19, 1 * 4, 2
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 3
    \\ xvstelm.w $xr15, $r19, 1 * 4, 3
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 4
    \\ xvstelm.w $xr15, $r19, 1 * 4, 4
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 5
    \\ xvstelm.w $xr15, $r19, 1 * 4, 5
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 6
    \\ xvstelm.w $xr15, $r19, 1 * 4, 6
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 7
    \\ xvstelm.w $xr15, $r19, 1 * 4, 7
    \\ add.d $r19, $r19, $r8
    \\ blt $r0, $r12, 20021b
    \\ b 20023f
    \\ .align 3
    \\ 20022:
    \\ ld.d $r14, $r5, 0 * 4
    \\ ld.d $r16, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r15, $r5, 0 * 4
    \\ ld.d $r17, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ xvinsgr2vr.w $xr18, $r14, 0
    \\ xvinsgr2vr.w $xr17, $r16, 0
    \\ xvinsgr2vr.w $xr18, $r15, 1
    \\ xvinsgr2vr.w $xr17, $r17, 1
    \\ ld.d $r14, $r5, 0 * 4
    \\ ld.d $r16, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r15, $r5, 0 * 4
    \\ ld.d $r17, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ xvinsgr2vr.w $xr18, $r14, 2
    \\ xvinsgr2vr.w $xr17, $r16, 2
    \\ xvinsgr2vr.w $xr18, $r15, 3
    \\ xvinsgr2vr.w $xr17, $r17, 3
    \\ ld.d $r14, $r5, 0 * 4
    \\ ld.d $r16, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r15, $r5, 0 * 4
    \\ ld.d $r17, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ xvinsgr2vr.w $xr18, $r14, 4
    \\ xvinsgr2vr.w $xr17, $r16, 4
    \\ xvinsgr2vr.w $xr18, $r15, 5
    \\ xvinsgr2vr.w $xr17, $r17, 5
    \\ ld.d $r14, $r5, 0 * 4
    \\ ld.d $r16, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r15, $r5, 0 * 4
    \\ ld.d $r17, $r5, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ xvinsgr2vr.w $xr18, $r14, 6
    \\ xvinsgr2vr.w $xr17, $r16, 6
    \\ xvinsgr2vr.w $xr18, $r15, 7
    \\ xvinsgr2vr.w $xr17, $r17, 7
    \\ ld.d $r14, $r7, 0 * 4
    \\ ld.d $r16, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r15, $r7, 0 * 4
    \\ ld.d $r17, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ xvinsgr2vr.w $xr16, $r14, 0
    \\ xvinsgr2vr.w $xr15, $r16, 0
    \\ xvinsgr2vr.w $xr16, $r15, 1
    \\ xvinsgr2vr.w $xr15, $r17, 1
    \\ ld.d $r14, $r7, 0 * 4
    \\ ld.d $r16, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r15, $r7, 0 * 4
    \\ ld.d $r17, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ xvinsgr2vr.w $xr16, $r14, 2
    \\ xvinsgr2vr.w $xr15, $r16, 2
    \\ xvinsgr2vr.w $xr16, $r15, 3
    \\ xvinsgr2vr.w $xr15, $r17, 3
    \\ ld.d $r14, $r7, 0 * 4
    \\ ld.d $r16, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r15, $r7, 0 * 4
    \\ ld.d $r17, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ xvinsgr2vr.w $xr16, $r14, 4
    \\ xvinsgr2vr.w $xr15, $r16, 4
    \\ xvinsgr2vr.w $xr16, $r15, 5
    \\ xvinsgr2vr.w $xr15, $r17, 5
    \\ ld.d $r14, $r7, 0 * 4
    \\ ld.d $r16, $r7, 1 * 4
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r15, $r7, 0 * 4
    \\ ld.d $r17, $r7, 1 * 4
    \\ xvinsgr2vr.w $xr16, $r14, 6
    \\ xvinsgr2vr.w $xr15, $r16, 6
    \\ xvinsgr2vr.w $xr16, $r15, 7
    \\ xvinsgr2vr.w $xr15, $r17, 7
    \\ add.d $r7, $r7, $r8
    \\ xvfmul.s $xr8, $xr19, $xr17
    \\ xvfmul.s $xr20, $xr19, $xr18
    \\ xvfmul.s $xr21, $xr13, $xr15
    \\ xvfmul.s $xr22, $xr13, $xr16
    \\ xvfmsub.s $xr8, $xr23, $xr18, $xr8
    \\ xvfmadd.s $xr20, $xr23, $xr17, $xr20
    \\ xvfmsub.s $xr21, $xr14, $xr16, $xr21
    \\ xvfmadd.s $xr22, $xr14, $xr15, $xr22
    \\ xvfadd.s $xr16, $xr8, $xr21
    \\ xvfadd.s $xr15, $xr20, $xr22
    \\ addi.d $r12, $r12, -1
    \\ xvstelm.w $xr16, $r19, 0 * 4, 0
    \\ xvstelm.w $xr15, $r19, 1 * 4, 0
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 1
    \\ xvstelm.w $xr15, $r19, 1 * 4, 1
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 2
    \\ xvstelm.w $xr15, $r19, 1 * 4, 2
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 3
    \\ xvstelm.w $xr15, $r19, 1 * 4, 3
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 4
    \\ xvstelm.w $xr15, $r19, 1 * 4, 4
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 5
    \\ xvstelm.w $xr15, $r19, 1 * 4, 5
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 6
    \\ xvstelm.w $xr15, $r19, 1 * 4, 6
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.w $xr16, $r19, 0 * 4, 7
    \\ xvstelm.w $xr15, $r19, 1 * 4, 7
    \\ add.d $r19, $r19, $r8
    \\ blt $r0, $r12, 20022b
    \\ b 20023f
    \\ .align 3
    \\ 20023:
    \\ andi $r12, $r4, 7
    \\ bge $r0, $r12, 20029f
    \\ .align 3
    \\ 20024:
    \\ movgr2fr.d $f16, $r0
    \\ fcmp.ceq.s $fcc0, $f0, $f16
    \\ fcmp.ceq.s $fcc1, $f1, $f16
    \\ fcmp.ceq.s $fcc2, $f2, $f16
    \\ fcmp.ceq.s $fcc3, $f3, $f16
    \\ addi.d $r12, $r12, -1
    \\ bceqz $fcc0, 20026f
    \\ bceqz $fcc1, 20026f
    \\ fmov.s $f17, $f16
    \\ bceqz $fcc2, 20025f
    \\ bcnez $fcc3, 20028f
    \\ 20025:
    \\ fneg.s $f16, $f16
    \\ fneg.s $f17, $f17
    \\ b 20027f
    \\ 20026:
    \\ fld.s $f12, $r5, 0 * 4
    \\ fld.s $f13, $r5, 1 * 4
    \\ fmul.s $f16, $f1, $f13
    \\ fmul.s $f17, $f1, $f12
    \\ fmsub.s $f16, $f0, $f12, $f16
    \\ fmadd.s $f17, $f13, $f0, $f17
    \\ bceqz $fcc2, 20027f
    \\ bcnez $fcc3, 20028f
    \\ 20027:
    \\ fld.s $f14, $r7, 0 * 4
    \\ fld.s $f15, $r7, 1 * 4
    \\ fmul.s $f18, $f3, $f15
    \\ fmul.s $f19, $f3, $f14
    \\ fmsub.s $f18, $f2, $f14, $f18
    \\ fmadd.s $f19, $f15, $f2, $f19
    \\ fadd.s $f16, $f18, $f16
    \\ fadd.s $f17, $f19, $f17
    \\ 20028:
    \\ fst.s $f16, $r7, 0 * 4
    \\ fst.s $f17, $r7, 1 * 4
    \\ add.d $r5, $r5, $r6
    \\ add.d $r7, $r7, $r8
    \\ blt $r0, $r12, 20024b
    \\ .align 3
    \\ 20029:
    \\ move $r4, $r12
    \\ jirl $r0, $r1, 0x0
    \\ .align 3
    \\
;
