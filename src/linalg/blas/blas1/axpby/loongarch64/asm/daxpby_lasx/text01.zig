//! Section 1 of 2 of the text of `daxpby`, which `__root__.zig` joins, from OpenBLAS `kernel/loongarch64/axpby_lasx.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .align 5
    \\ bge $r0, $r4, 20030f
    \\ movgr2fr.d $f12, $r0
    \\ ffint.d.l $f12, $f12
    \\ slli.d $r6, $r6, 3
    \\ slli.d $r8, $r8, 3
    \\ movfr2gr.d $r14, $f0
    \\ movfr2gr.d $r16, $f1
    \\ movfr2gr.d $r15, $f12
    \\ xvreplgr2vr.d $xr23, $r14
    \\ xvreplgr2vr.d $xr9, $r16
    \\ xvreplgr2vr.d $xr19, $r15
    \\ and $r13, $r6, $r8
    \\ or $r12, $r4, $r4
    \\ beqz $r13, 20026f
    \\ li.d $r13, 1
    \\ slli.d $r13, $r13, 3
    \\ srai.d $r12, $r4, 3
    \\ bne $r6, $r13, 20000f
    \\ bne $r8, $r13, 20007f
    \\ b 20001f
    \\ 20000:
    \\ bne $r8, $r13, 20019f
    \\ b 20013f
    \\ 20001:
    \\ bge $r0, $r12, 20025f
    \\ fcmp.ceq.d $fcc0, $f0, $f12
    \\ bcnez $fcc0, 20002f
    \\ fcmp.ceq.d $fcc0, $f1, $f12
    \\ bcnez $fcc0, 20004f
    \\ b 20003f
    \\ .align 3
    \\ 20002:
    \\ fcmp.ceq.d $fcc0, $f1, $f12
    \\ bcnez $fcc0, 20006f
    \\ b 20005f
    \\ .align 3
    \\ 20003:
    \\ xvld $xr8, $r5, 0 * 8
    \\ xvld $xr21, $r7, 0 * 8
    \\ xvld $xr20, $r5, 4 * 8
    \\ xvld $xr22, $r7, 4 * 8
    \\ xvfmul.d $xr8, $xr8, $xr23
    \\ xvfmul.d $xr20, $xr20, $xr23
    \\ xvfmadd.d $xr21, $xr21, $xr9, $xr8
    \\ xvfmadd.d $xr22, $xr22, $xr9, $xr20
    \\ addi.d $r12, $r12, -1
    \\ xvst $xr21, $r7, 0 * 8
    \\ xvst $xr22, $r7, 4 * 8
    \\ addi.d $r5, $r5, 8 * 8
    \\ addi.d $r7, $r7, 8 * 8
    \\ blt $r0, $r12, 20003b
    \\ b 20025f
    \\ .align 3
    \\ 20004:
    \\ xvld $xr8, $r5, 0 * 8
    \\ xvld $xr20, $r5, 4 * 8
    \\ xvfmul.d $xr8, $xr8, $xr23
    \\ xvfmul.d $xr20, $xr20, $xr23
    \\ xvst $xr8, $r7, 0 * 8
    \\ xvst $xr20, $r7, 4 * 8
    \\ addi.d $r12, $r12, -1
    \\ addi.d $r5, $r5, 8 * 8
    \\ addi.d $r7, $r7, 8 * 8
    \\ blt $r0, $r12, 20004b
    \\ b 20025f
    \\ .align 3
    \\ 20005:
    \\ xvld $xr21, $r7, 0 * 8
    \\ xvld $xr22, $r7, 4 * 8
    \\ xvfmul.d $xr21, $xr21, $xr9
    \\ xvfmul.d $xr22, $xr22, $xr9
    \\ xvst $xr21, $r7, 0 * 8
    \\ xvst $xr22, $r7, 4 * 8
    \\ addi.d $r12, $r12, -1
    \\ addi.d $r7, $r7, 8 * 8
    \\ blt $r0, $r12, 20005b
    \\ b 20025f
    \\ .align 3
    \\ 20006:
    \\ xvst $xr19, $r7, 0 * 8
    \\ xvst $xr19, $r7, 4 * 8
    \\ addi.d $r7, $r7, 8 * 8
    \\ addi.d $r12, $r12, -1
    \\ blt $r0, $r12, 20006b
    \\ b 20025f
    \\ .align 3
    \\ 20007:
    \\ bge $r0, $r12, 20025f
    \\ move $r19, $r7
    \\ fcmp.ceq.d $fcc0, $f0, $f12
    \\ bcnez $fcc0, 20008f
    \\ fcmp.ceq.d $fcc0, $f1, $f12
    \\ bcnez $fcc0, 20010f
    \\ b 20009f
    \\ .align 3
    \\ 20008:
    \\ fcmp.ceq.d $fcc0, $f1, $f12
    \\ bcnez $fcc0, 20012f
    \\ b 20011f
    \\ .align 3
    \\ 20009:
    \\ xvld $xr8, $r5, 0 * 8
    \\ ld.d $r14, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r16, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r15, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r17, $r7, 0 * 8
    \\ xvinsgr2vr.d $xr21, $r14, 0
    \\ xvinsgr2vr.d $xr21, $r16, 1
    \\ xvinsgr2vr.d $xr21, $r15, 2
    \\ xvinsgr2vr.d $xr21, $r17, 3
    \\ add.d $r7, $r7, $r8
    \\ xvfmul.d $xr8, $xr8, $xr23
    \\ xvld $xr20, $r5, 4 * 8
    \\ xvfmadd.d $xr21, $xr21, $xr9, $xr8
    \\ ld.d $r14, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r16, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r15, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r17, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ xvinsgr2vr.d $xr22, $r14, 0
    \\ xvinsgr2vr.d $xr22, $r16, 1
    \\ xvinsgr2vr.d $xr22, $r15, 2
    \\ xvinsgr2vr.d $xr22, $r17, 3
    \\ xvstelm.d $xr21, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr21, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr21, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr21, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ xvfmul.d $xr20, $xr20, $xr23
    \\ xvfmadd.d $xr22, $xr22, $xr9, $xr20
    \\ addi.d $r12, $r12, -1
    \\ xvstelm.d $xr22, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr22, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr22, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr22, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ addi.d $r5, $r5, 8 * 8
    \\ blt $r0, $r12, 20009b
    \\ move $r7, $r19
    \\ b 20025f
    \\ .align 3
    \\ 20010:
    \\ xvld $xr8, $r5, 0 * 8
    \\ addi.d $r12, $r12, -1
    \\ xvld $xr20, $r5, 4 * 8
    \\ xvfmul.d $xr8, $xr8, $xr23
    \\ xvfmul.d $xr20, $xr20, $xr23
    \\ xvstelm.d $xr8, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr8, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr8, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr8, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr20, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr20, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr20, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr20, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ addi.d $r5, $r5, 8 * 8
    \\ blt $r0, $r12, 20010b
    \\ move $r7, $r19
    \\ b 20025f
    \\ .align 3
    \\ 20011:
    \\ ld.d $r14, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r16, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r15, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r17, $r7, 0 * 8
    \\ xvinsgr2vr.d $xr21, $r14, 0
    \\ xvinsgr2vr.d $xr21, $r16, 1
    \\ xvinsgr2vr.d $xr21, $r15, 2
    \\ xvinsgr2vr.d $xr21, $r17, 3
    \\ add.d $r7, $r7, $r8
    \\ xvfmul.d $xr21, $xr21, $xr9
    \\ ld.d $r14, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r16, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r15, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ ld.d $r17, $r7, 0 * 8
    \\ add.d $r7, $r7, $r8
    \\ xvinsgr2vr.d $xr22, $r14, 0
    \\ xvinsgr2vr.d $xr22, $r16, 1
    \\ xvinsgr2vr.d $xr22, $r15, 2
    \\ xvinsgr2vr.d $xr22, $r17, 3
    \\ xvstelm.d $xr21, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr21, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr21, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr21, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ xvfmul.d $xr22, $xr22, $xr9
    \\ xvstelm.d $xr22, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr22, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr22, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr22, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ addi.d $r12, $r12, -1
    \\ blt $r0, $r12, 20011b
    \\ move $r7, $r19
    \\ b 20025f
    \\ .align 3
    \\ 20012:
    \\ xvstelm.d $xr19, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr19, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr19, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr19, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr19, $r19, 0, 0
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr19, $r19, 0, 1
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr19, $r19, 0, 2
    \\ add.d $r19, $r19, $r8
    \\ xvstelm.d $xr19, $r19, 0, 3
    \\ add.d $r19, $r19, $r8
    \\ addi.d $r12, $r12, -1
    \\ blt $r0, $r12, 20012b
    \\ move $r7, $r19
    \\ b 20025f
    \\ .align 3
    \\ 20013:
    \\ bge $r0, $r12, 20025f
    \\ fcmp.ceq.d $fcc0, $f0, $f12
    \\ bcnez $fcc0, 20014f
    \\ fcmp.ceq.d $fcc0, $f1, $f12
    \\ bcnez $fcc0, 20016f
    \\ b 20015f
    \\ .align 3
    \\ 20014:
    \\ fcmp.ceq.d $fcc0, $f1, $f12
    \\ bcnez $fcc0, 20018f
    \\ b 20017f
    \\ .align 3
    \\ 20015:
    \\ xvld $xr21, $r7, 0 * 8
    \\ ld.d $r14, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r16, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r15, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r17, $r5, 0 * 8
    \\ xvinsgr2vr.d $xr8, $r14, 0
    \\ xvinsgr2vr.d $xr8, $r16, 1
    \\ xvinsgr2vr.d $xr8, $r15, 2
    \\ xvinsgr2vr.d $xr8, $r17, 3
    \\ add.d $r5, $r5, $r6
    \\ xvfmul.d $xr8, $xr23, $xr8
    \\ xvfmadd.d $xr21, $xr21, $xr9, $xr8
    \\ xvld $xr22, $r7, 4 * 8
    \\ xvst $xr21, $r7, 0 * 8
    \\ ld.d $r14, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r16, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r15, $r5, 0 * 8
    \\ add.d $r5, $r5, $r6
    \\ ld.d $r17, $r5, 0 * 8
    \\ xvinsgr2vr.d $xr20, $r14, 0
    \\ xvinsgr2vr.d $xr20, $r16, 1
    \\ xvinsgr2vr.d $xr20, $r15, 2
    \\ xvinsgr2vr.d $xr20, $r17, 3
    \\ add.d $r5, $r5, $r6
    \\ xvfmul.d $xr20, $xr20, $xr23
    \\ xvfmadd.d $xr22, $xr22, $xr9, $xr20
    \\ addi.d $r12, $r12, -1
    \\ xvst $xr22, $r7, 4 * 8
    \\ addi.d $r7, $r7, 8 * 8
    \\ blt $r0, $r12, 20015b
    \\ b 20025f
    \\
;
