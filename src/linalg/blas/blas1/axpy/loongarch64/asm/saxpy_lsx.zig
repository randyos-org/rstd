//! `saxpy` for la264, from OpenBLAS `kernel/loongarch64/axpy_lsx.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

/// Each function's core and variant.
pub const keys = .{
    .{ .core = .la264, .variant = "" },
};

pub const Body = struct {
    pub fn entry() callconv(.naked) void {
        asm volatile (
            \\ .align 5
            \\ bge $r0, $r4, 20024f
            \\ li.d $r13, 1
            \\ movgr2fr.d $f12, $r0
            \\ ffint.s.l $f12, $f12
            \\ movgr2fr.d $f13, $r13
            \\ ffint.s.l $f13, $f13
            \\ slli.d $r13, $r13, 2
            \\ slli.d $r8, $r8, 2
            \\ slli.d $r10, $r10, 2
            \\ movfr2gr.s $r14, $f0
            \\ vreplgr2vr.w $vr23, $r14
            \\ beq $r0, $r10, 20000f
            \\ srai.d $r12, $r4, 3
            \\ b 20001f
            \\ 20000:
            \\ srai.d $r12, $r4, 0
            \\ 20001:
            \\ bne $r8, $r13, 20002f
            \\ bne $r10, $r13, 20008f
            \\ b 20003f
            \\ 20002:
            \\ bne $r10, $r13, 20018f
            \\ b 20014f
            \\ 20003:
            \\ bge $r0, $r12, 20006f
            \\ fcmp.ceq.s $fcc0, $f0, $f13
            \\ bceqz $fcc0, 20005f
            \\ .align 3
            \\ 20004:
            \\ vld $vr8, $r7, 0 * 4
            \\ vld $vr21, $r9, 0 * 4
            \\ vld $vr20, $r7, 4 * 4
            \\ vld $vr22, $r9, 4 * 4
            \\ vfadd.s $vr21, $vr8, $vr21
            \\ vfadd.s $vr22, $vr20, $vr22
            \\ vst $vr21, $r9, 0 * 4
            \\ vst $vr22, $r9, 4 * 4
            \\ addi.d $r7, $r7, 8 * 4
            \\ addi.d $r9, $r9, 8 * 4
            \\ addi.d $r12, $r12, -1
            \\ blt $r0, $r12, 20004b
            \\ b 20006f
            \\ .align 3
            \\ 20005:
            \\ vld $vr8, $r7, 0 * 4
            \\ vld $vr21, $r9, 0 * 4
            \\ vld $vr20, $r7, 4 * 4
            \\ vld $vr22, $r9, 4 * 4
            \\ vfmadd.s $vr21, $vr8, $vr23, $vr21
            \\ vfmadd.s $vr22, $vr20, $vr23, $vr22
            \\ vst $vr21, $r9, 0 * 4
            \\ vst $vr22, $r9, 4 * 4
            \\ addi.d $r7, $r7, 8 * 4
            \\ addi.d $r9, $r9, 8 * 4
            \\ addi.d $r12, $r12, -1
            \\ blt $r0, $r12, 20005b
            \\ .align 3
            \\ 20006:
            \\ andi $r12, $r4, 7
            \\ bge $r0, $r12, 20024f
            \\ .align 3
            \\ 20007:
            \\ fld.s $f12, $r7, 0 * 4
            \\ fld.s $f14, $r9, 0 * 4
            \\ addi.d $r12, $r12, -1
            \\ fmadd.s $f14, $f12, $f0, $f14
            \\ fst.s $f14, $r9, 0 * 4
            \\ addi.d $r7, $r7, 4
            \\ addi.d $r9, $r9, 4
            \\ blt $r0, $r12, 20007b
            \\ b 20024f
            \\ .align 3
            \\ 20008:
            \\ bge $r0, $r12, 20010f
            \\ move $r6, $r9
            \\ .align 3
            \\ 20009:
            \\ beq $r0, $r10, 20012f
            \\ vld $vr8, $r7, 0 * 4
            \\ ld.w $r14, $r9, 0 * 4
            \\ add.d $r9, $r9, $r10
            \\ ld.w $r16, $r9, 0 * 4
            \\ add.d $r9, $r9, $r10
            \\ ld.w $r15, $r9, 0 * 4
            \\ add.d $r9, $r9, $r10
            \\ ld.w $r17, $r9, 0 * 4
            \\ vinsgr2vr.w $vr21, $r14, 0
            \\ vinsgr2vr.w $vr21, $r16, 1
            \\ vinsgr2vr.w $vr21, $r15, 2
            \\ vinsgr2vr.w $vr21, $r17, 3
            \\ add.d $r9, $r9, $r10
            \\ vfmadd.s $vr21, $vr8, $vr23, $vr21
            \\ vld $vr20, $r7, 4 * 4
            \\ vstelm.w $vr21, $r6, 0, 0
            \\ add.d $r6, $r6, $r10
            \\ vstelm.w $vr21, $r6, 0, 1
            \\ add.d $r6, $r6, $r10
            \\ vstelm.w $vr21, $r6, 0, 2
            \\ add.d $r6, $r6, $r10
            \\ vstelm.w $vr21, $r6, 0, 3
            \\ add.d $r6, $r6, $r10
            \\ ld.w $r14, $r9, 0 * 4
            \\ add.d $r9, $r9, $r10
            \\ ld.w $r16, $r9, 0 * 4
            \\ add.d $r9, $r9, $r10
            \\ ld.w $r15, $r9, 0 * 4
            \\ add.d $r9, $r9, $r10
            \\ ld.w $r17, $r9, 0 * 4
            \\ vinsgr2vr.w $vr22, $r14, 0
            \\ vinsgr2vr.w $vr22, $r16, 1
            \\ vinsgr2vr.w $vr22, $r15, 2
            \\ vinsgr2vr.w $vr22, $r17, 3
            \\ add.d $r9, $r9, $r10
            \\ vfmadd.s $vr22, $vr20, $vr23, $vr22
            \\ addi.d $r12, $r12, -1
            \\ vstelm.w $vr22, $r6, 0, 0
            \\ add.d $r6, $r6, $r10
            \\ vstelm.w $vr22, $r6, 0, 1
            \\ add.d $r6, $r6, $r10
            \\ vstelm.w $vr22, $r6, 0, 2
            \\ add.d $r6, $r6, $r10
            \\ vstelm.w $vr22, $r6, 0, 3
            \\ add.d $r6, $r6, $r10
            \\ addi.d $r7, $r7, 8 * 4
            \\ blt $r0, $r12, 20009b
            \\ .align 3
            \\ 20010:
            \\ andi $r12, $r4, 7
            \\ bge $r0, $r12, 20024f
            \\ .align 3
            \\ 20011:
            \\ beq $r0, $r10, 20012f
            \\ fld.s $f12, $r7, 0 * 4
            \\ fld.s $f14, $r9, 0 * 4
            \\ addi.d $r12, $r12, -1
            \\ fmadd.s $f14, $f12, $f0, $f14
            \\ fst.s $f14, $r9, 0 * 4
            \\ addi.d $r7, $r7, 4
            \\ add.d $r9, $r9, $r10
            \\ blt $r0, $r12, 20011b
            \\ b 20024f
            \\ .align 3
            \\ 20012:
            \\ fld.s $f14, $r9, 0 * 4
            \\ 20013:
            \\ fld.s $f12, $r7, 0 * 4
            \\ addi.d $r12, $r12, -1
            \\ fmadd.s $f14, $f12, $f0, $f14
            \\ addi.d $r7, $r7, 4
            \\ blt $r0, $r12, 20013b
            \\ fst.s $f14, $r9, 0 * 4
            \\ b 20024f
            \\ .align 3
            \\ 20014:
            \\ bge $r0, $r12, 20016f
            \\ .align 3
            \\ 20015:
            \\ vld $vr21, $r9, 0 * 4
            \\ ld.w $r14, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r16, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r15, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r17, $r7, 0 * 4
            \\ vinsgr2vr.w $vr8, $r14, 0
            \\ vinsgr2vr.w $vr8, $r16, 1
            \\ vinsgr2vr.w $vr8, $r15, 2
            \\ vinsgr2vr.w $vr8, $r17, 3
            \\ add.d $r7, $r7, $r8
            \\ vfmadd.s $vr21, $vr8, $vr23, $vr21
            \\ vld $vr22, $r9, 4 * 4
            \\ vst $vr21, $r9, 0 * 4
            \\ ld.w $r14, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r16, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r15, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r17, $r7, 0 * 4
            \\ vinsgr2vr.w $vr20, $r14, 0
            \\ vinsgr2vr.w $vr20, $r16, 1
            \\ vinsgr2vr.w $vr20, $r15, 2
            \\ vinsgr2vr.w $vr20, $r17, 3
            \\ add.d $r7, $r7, $r8
            \\ vfmadd.s $vr22, $vr20, $vr23, $vr22
            \\ addi.d $r12, $r12, -1
            \\ vst $vr22, $r9, 4 * 4
            \\ addi.d $r9, $r9, 8 * 4
            \\ blt $r0, $r12, 20015b
            \\ .align 3
            \\ 20016:
            \\ andi $r12, $r4, 7
            \\ bge $r0, $r12, 20024f
            \\ .align 3
            \\ 20017:
            \\ fld.s $f12, $r7, 0 * 4
            \\ fld.s $f14, $r9, 0 * 4
            \\ addi.d $r12, $r12, -1
            \\ fmadd.s $f14, $f12, $f0, $f14
            \\ fst.s $f14, $r9, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ addi.d $r9, $r9, 4
            \\ blt $r0, $r12, 20017b
            \\ b 20024f
            \\ .align 3
            \\ 20018:
            \\ bge $r0, $r12, 20020f
            \\ move $r6, $r9
            \\ .align 3
            \\ 20019:
            \\ beq $r0, $r10, 20022f
            \\ ld.w $r14, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r16, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r15, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r17, $r7, 0 * 4
            \\ vinsgr2vr.w $vr8, $r14, 0
            \\ vinsgr2vr.w $vr8, $r16, 1
            \\ vinsgr2vr.w $vr8, $r15, 2
            \\ vinsgr2vr.w $vr8, $r17, 3
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r14, $r9, 0 * 4
            \\ add.d $r9, $r9, $r10
            \\ ld.w $r16, $r9, 0 * 4
            \\ add.d $r9, $r9, $r10
            \\ ld.w $r15, $r9, 0 * 4
            \\ add.d $r9, $r9, $r10
            \\ ld.w $r17, $r9, 0 * 4
            \\ vinsgr2vr.w $vr21, $r14, 0
            \\ vinsgr2vr.w $vr21, $r16, 1
            \\ vinsgr2vr.w $vr21, $r15, 2
            \\ vinsgr2vr.w $vr21, $r17, 3
            \\ add.d $r9, $r9, $r10
            \\ vfmadd.s $vr21, $vr8, $vr23, $vr21
            \\ ld.w $r14, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r16, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r15, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ ld.w $r17, $r7, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ vinsgr2vr.w $vr20, $r14, 0
            \\ vinsgr2vr.w $vr20, $r16, 1
            \\ vinsgr2vr.w $vr20, $r15, 2
            \\ vinsgr2vr.w $vr20, $r17, 3
            \\ vstelm.w $vr21, $r6, 0, 0
            \\ add.d $r6, $r6, $r10
            \\ vstelm.w $vr21, $r6, 0, 1
            \\ add.d $r6, $r6, $r10
            \\ vstelm.w $vr21, $r6, 0, 2
            \\ add.d $r6, $r6, $r10
            \\ vstelm.w $vr21, $r6, 0, 3
            \\ add.d $r6, $r6, $r10
            \\ ld.w $r14, $r9, 0 * 4
            \\ add.d $r9, $r9, $r10
            \\ ld.w $r16, $r9, 0 * 4
            \\ add.d $r9, $r9, $r10
            \\ ld.w $r15, $r9, 0 * 4
            \\ add.d $r9, $r9, $r10
            \\ ld.w $r17, $r9, 0 * 4
            \\ vinsgr2vr.w $vr22, $r14, 0
            \\ vinsgr2vr.w $vr22, $r16, 1
            \\ vinsgr2vr.w $vr22, $r15, 2
            \\ vinsgr2vr.w $vr22, $r17, 3
            \\ add.d $r9, $r9, $r10
            \\ vfmadd.s $vr22, $vr20, $vr23, $vr22
            \\ addi.d $r12, $r12, -1
            \\ vstelm.w $vr22, $r6, 0, 0
            \\ add.d $r6, $r6, $r10
            \\ vstelm.w $vr22, $r6, 0, 1
            \\ add.d $r6, $r6, $r10
            \\ vstelm.w $vr22, $r6, 0, 2
            \\ add.d $r6, $r6, $r10
            \\ vstelm.w $vr22, $r6, 0, 3
            \\ add.d $r6, $r6, $r10
            \\ blt $r0, $r12, 20019b
            \\ .align 3
            \\ 20020:
            \\ andi $r12, $r4, 7
            \\ bge $r0, $r12, 20024f
            \\ .align 3
            \\ 20021:
            \\ beq $r0, $r10, 20022f
            \\ fld.s $f12, $r7, 0 * 4
            \\ fld.s $f14, $r9, 0 * 4
            \\ addi.d $r12, $r12, -1
            \\ fmadd.s $f14, $f12, $f0, $f14
            \\ fst.s $f14, $r9, 0 * 4
            \\ add.d $r7, $r7, $r8
            \\ add.d $r9, $r9, $r10
            \\ blt $r0, $r12, 20021b
            \\ .align 3
            \\ b 20024f
            \\ 20022:
            \\ fld.s $f14, $r9, 0 * 4
            \\ 20023:
            \\ fld.s $f12, $r7, 0 * 4
            \\ addi.d $r12, $r12, -1
            \\ fmadd.s $f14, $f12, $f0, $f14
            \\ add.d $r7, $r7, $r8
            \\ blt $r0, $r12, 20023b
            \\ fst.s $f14, $r9, 0 * 4
            \\ .align 3
            \\ 20024:
            \\ move $r4, $r12
            \\ jirl $r0, $r1, 0x0
            \\ .align 3
            \\
        );
    }
};
