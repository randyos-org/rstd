//! `daxpy` for la264, from OpenBLAS `kernel/loongarch64/axpy_lsx.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

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
            \\ ffint.d.l $f12, $f12
            \\ movgr2fr.d $f13, $r13
            \\ ffint.d.l $f13, $f13
            \\ slli.d $r13, $r13, 3
            \\ slli.d $r8, $r8, 3
            \\ slli.d $r10, $r10, 3
            \\ movfr2gr.d $r14, $f0
            \\ vreplgr2vr.d $vr23, $r14
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
            \\ fcmp.ceq.d $fcc0, $f0, $f13
            \\ bceqz $fcc0, 20005f
            \\ .align 3
            \\ 20004:
            \\ vld $vr8, $r7, 0 * 8
            \\ vld $vr21, $r9, 0 * 8
            \\ vld $vr20, $r7, 2 * 8
            \\ vld $vr22, $r9, 2 * 8
            \\ vfadd.d $vr21, $vr8, $vr21
            \\ vfadd.d $vr22, $vr20, $vr22
            \\ vst $vr21, $r9, 0 * 8
            \\ vst $vr22, $r9, 2 * 8
            \\ vld $vr8, $r7, 4 * 8
            \\ vld $vr21, $r9, 4 * 8
            \\ vld $vr20, $r7, 6 * 8
            \\ vld $vr22, $r9, 6 * 8
            \\ vfadd.d $vr21, $vr8, $vr21
            \\ vfadd.d $vr22, $vr20, $vr22
            \\ vst $vr21, $r9, 4 * 8
            \\ vst $vr22, $r9, 6 * 8
            \\ addi.d $r7, $r7, 8 * 8
            \\ addi.d $r9, $r9, 8 * 8
            \\ addi.d $r12, $r12, -1
            \\ blt $r0, $r12, 20004b
            \\ b 20006f
            \\ .align 3
            \\ 20005:
            \\ vld $vr8, $r7, 0 * 8
            \\ vld $vr21, $r9, 0 * 8
            \\ vld $vr20, $r7, 2 * 8
            \\ vld $vr22, $r9, 2 * 8
            \\ vfmadd.d $vr21, $vr8, $vr23, $vr21
            \\ vfmadd.d $vr22, $vr20, $vr23, $vr22
            \\ addi.d $r12, $r12, -1
            \\ vst $vr21, $r9, 0 * 8
            \\ vst $vr22, $r9, 2 * 8
            \\ vld $vr8, $r7, 4 * 8
            \\ vld $vr21, $r9, 4 * 8
            \\ vld $vr20, $r7, 6 * 8
            \\ vld $vr22, $r9, 6 * 8
            \\ addi.d $r7, $r7, 8 * 8
            \\ vfmadd.d $vr21, $vr8, $vr23, $vr21
            \\ vfmadd.d $vr22, $vr20, $vr23, $vr22
            \\ vst $vr21, $r9, 4 * 8
            \\ vst $vr22, $r9, 6 * 8
            \\ addi.d $r9, $r9, 8 * 8
            \\ blt $r0, $r12, 20005b
            \\ .align 3
            \\ 20006:
            \\ andi $r12, $r4, 7
            \\ bge $r0, $r12, 20024f
            \\ .align 3
            \\ 20007:
            \\ fld.d $f12, $r7, 0 * 8
            \\ fld.d $f14, $r9, 0 * 8
            \\ addi.d $r12, $r12, -1
            \\ fmadd.d $f14, $f12, $f0, $f14
            \\ fst.d $f14, $r9, 0 * 8
            \\ addi.d $r7, $r7, 8
            \\ addi.d $r9, $r9, 8
            \\ blt $r0, $r12, 20007b
            \\ b 20024f
            \\ .align 3
            \\ 20008:
            \\ bge $r0, $r12, 20010f
            \\ move $r6, $r9
            \\ .align 3
            \\ 20009:
            \\ beq $r0, $r10, 20012f
            \\ vld $vr8, $r7, 0 * 8
            \\ ld.d $r14, $r9, 0 * 8
            \\ add.d $r9, $r9, $r10
            \\ ld.d $r16, $r9, 0 * 8
            \\ vinsgr2vr.d $vr21, $r14, 0
            \\ vinsgr2vr.d $vr21, $r16, 1
            \\ add.d $r9, $r9, $r10
            \\ vfmadd.d $vr21, $vr8, $vr23, $vr21
            \\ vld $vr20, $r7, 2 * 8
            \\ vstelm.d $vr21, $r6, 0, 0
            \\ add.d $r6, $r6, $r10
            \\ vstelm.d $vr21, $r6, 0, 1
            \\ add.d $r6, $r6, $r10
            \\ ld.d $r15, $r9, 0 * 8
            \\ add.d $r9, $r9, $r10
            \\ ld.d $r17, $r9, 0 * 8
            \\ vinsgr2vr.d $vr22, $r15, 0
            \\ vinsgr2vr.d $vr22, $r17, 1
            \\ add.d $r9, $r9, $r10
            \\ vfmadd.d $vr22, $vr20, $vr23, $vr22
            \\ vld $vr8, $r7, 4 * 8
            \\ vstelm.d $vr22, $r6, 0, 0
            \\ add.d $r6, $r6, $r10
            \\ vstelm.d $vr22, $r6, 0, 1
            \\ add.d $r6, $r6, $r10
            \\ ld.d $r14, $r9, 0 * 8
            \\ add.d $r9, $r9, $r10
            \\ ld.d $r16, $r9, 0 * 8
            \\ vinsgr2vr.d $vr21, $r14, 0
            \\ vinsgr2vr.d $vr21, $r16, 1
            \\ add.d $r9, $r9, $r10
            \\ vfmadd.d $vr21, $vr8, $vr23, $vr21
            \\ vld $vr20, $r7, 6 * 8
            \\ vstelm.d $vr21, $r6, 0, 0
            \\ add.d $r6, $r6, $r10
            \\ vstelm.d $vr21, $r6, 0, 1
            \\ add.d $r6, $r6, $r10
            \\ ld.d $r15, $r9, 0 * 8
            \\ add.d $r9, $r9, $r10
            \\ ld.d $r17, $r9, 0 * 8
            \\ vinsgr2vr.d $vr22, $r15, 0
            \\ vinsgr2vr.d $vr22, $r17, 1
            \\ add.d $r9, $r9, $r10
            \\ vfmadd.d $vr22, $vr20, $vr23, $vr22
            \\ vstelm.d $vr22, $r6, 0, 0
            \\ add.d $r6, $r6, $r10
            \\ vstelm.d $vr22, $r6, 0, 1
            \\ add.d $r6, $r6, $r10
            \\ addi.d $r7, $r7, 8 * 8
            \\ addi.d $r12, $r12, -1
            \\ blt $r0, $r12, 20009b
            \\ .align 3
            \\ 20010:
            \\ andi $r12, $r4, 7
            \\ bge $r0, $r12, 20024f
            \\ .align 3
            \\ 20011:
            \\ beq $r0, $r10, 20012f
            \\ fld.d $f12, $r7, 0 * 8
            \\ fld.d $f14, $r9, 0 * 8
            \\ addi.d $r12, $r12, -1
            \\ fmadd.d $f14, $f12, $f0, $f14
            \\ fst.d $f14, $r9, 0 * 8
            \\ addi.d $r7, $r7, 8
            \\ add.d $r9, $r9, $r10
            \\ blt $r0, $r12, 20011b
            \\ b 20024f
            \\ .align 3
            \\ 20012:
            \\ fld.d $f14, $r9, 0 * 8
            \\ 20013:
            \\ fld.d $f12, $r7, 0 * 8
            \\ addi.d $r12, $r12, -1
            \\ fmadd.d $f14, $f12, $f0, $f14
            \\ addi.d $r7, $r7, 8
            \\ blt $r0, $r12, 20013b
            \\ fst.d $f14, $r9, 0 * 8
            \\ b 20024f
            \\ .align 3
            \\ 20014:
            \\ bge $r0, $r12, 20016f
            \\ .align 3
            \\ 20015:
            \\ vld $vr21, $r9, 0 * 8
            \\ ld.d $r14, $r7, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ ld.d $r16, $r7, 0 * 8
            \\ vinsgr2vr.d $vr8, $r14, 0
            \\ vinsgr2vr.d $vr8, $r16, 1
            \\ add.d $r7, $r7, $r8
            \\ vfmadd.d $vr21, $vr8, $vr23, $vr21
            \\ vld $vr22, $r9, 2 * 8
            \\ vst $vr21, $r9, 0 * 8
            \\ ld.d $r15, $r7, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ ld.d $r17, $r7, 0 * 8
            \\ vinsgr2vr.d $vr20, $r15, 0
            \\ vinsgr2vr.d $vr20, $r17, 1
            \\ add.d $r7, $r7, $r8
            \\ vfmadd.d $vr22, $vr20, $vr23, $vr22
            \\ vld $vr21, $r9, 4 * 8
            \\ vst $vr22, $r9, 2 * 8
            \\ ld.d $r14, $r7, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ ld.d $r16, $r7, 0 * 8
            \\ vinsgr2vr.d $vr8, $r14, 0
            \\ vinsgr2vr.d $vr8, $r16, 1
            \\ add.d $r7, $r7, $r8
            \\ vfmadd.d $vr21, $vr8, $vr23, $vr21
            \\ vld $vr22, $r9, 6 * 8
            \\ vst $vr21, $r9, 4 * 8
            \\ ld.d $r15, $r7, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ ld.d $r17, $r7, 0 * 8
            \\ vinsgr2vr.d $vr20, $r15, 0
            \\ vinsgr2vr.d $vr20, $r17, 1
            \\ add.d $r7, $r7, $r8
            \\ vfmadd.d $vr22, $vr20, $vr23, $vr22
            \\ addi.d $r12, $r12, -1
            \\ vst $vr22, $r9, 6 * 8
            \\ addi.d $r9, $r9, 8 * 8
            \\ blt $r0, $r12, 20015b
            \\ .align 3
            \\ 20016:
            \\ andi $r12, $r4, 7
            \\ bge $r0, $r12, 20024f
            \\ .align 3
            \\ 20017:
            \\ fld.d $f12, $r7, 0 * 8
            \\ fld.d $f14, $r9, 0 * 8
            \\ addi.d $r12, $r12, -1
            \\ fmadd.d $f14, $f12, $f0, $f14
            \\ fst.d $f14, $r9, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ addi.d $r9, $r9, 8
            \\ blt $r0, $r12, 20017b
            \\ b 20024f
            \\ .align 3
            \\ 20018:
            \\ bge $r0, $r12, 20020f
            \\ move $r6, $r9
            \\ .align 3
            \\ 20019:
            \\ beq $r0, $r10, 20022f
            \\ ld.d $r14, $r7, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ ld.d $r16, $r7, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ vinsgr2vr.d $vr8, $r14, 0
            \\ vinsgr2vr.d $vr8, $r16, 1
            \\ ld.d $r14, $r9, 0 * 8
            \\ add.d $r9, $r9, $r10
            \\ ld.d $r16, $r9, 0 * 8
            \\ vinsgr2vr.d $vr21, $r14, 0
            \\ vinsgr2vr.d $vr21, $r16, 1
            \\ add.d $r9, $r9, $r10
            \\ vfmadd.d $vr21, $vr8, $vr23, $vr21
            \\ ld.d $r15, $r7, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ ld.d $r17, $r7, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ vinsgr2vr.d $vr20, $r15, 0
            \\ vinsgr2vr.d $vr20, $r17, 1
            \\ vstelm.d $vr21, $r6, 0, 0
            \\ add.d $r6, $r6, $r10
            \\ vstelm.d $vr21, $r6, 0, 1
            \\ add.d $r6, $r6, $r10
            \\ ld.d $r15, $r9, 0 * 8
            \\ add.d $r9, $r9, $r10
            \\ ld.d $r17, $r9, 0 * 8
            \\ vinsgr2vr.d $vr22, $r15, 0
            \\ vinsgr2vr.d $vr22, $r17, 1
            \\ add.d $r9, $r9, $r10
            \\ vfmadd.d $vr22, $vr20, $vr23, $vr22
            \\ ld.d $r14, $r7, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ ld.d $r16, $r7, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ vinsgr2vr.d $vr8, $r14, 0
            \\ vinsgr2vr.d $vr8, $r16, 1
            \\ vstelm.d $vr22, $r6, 0, 0
            \\ add.d $r6, $r6, $r10
            \\ vstelm.d $vr22, $r6, 0, 1
            \\ add.d $r6, $r6, $r10
            \\ ld.d $r14, $r9, 0 * 8
            \\ add.d $r9, $r9, $r10
            \\ ld.d $r16, $r9, 0 * 8
            \\ vinsgr2vr.d $vr21, $r14, 0
            \\ vinsgr2vr.d $vr21, $r16, 1
            \\ add.d $r9, $r9, $r10
            \\ vfmadd.d $vr21, $vr8, $vr23, $vr21
            \\ ld.d $r15, $r7, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ ld.d $r17, $r7, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ vinsgr2vr.d $vr20, $r15, 0
            \\ vinsgr2vr.d $vr20, $r17, 1
            \\ vstelm.d $vr21, $r6, 0, 0
            \\ add.d $r6, $r6, $r10
            \\ vstelm.d $vr21, $r6, 0, 1
            \\ add.d $r6, $r6, $r10
            \\ ld.d $r14, $r9, 0 * 8
            \\ add.d $r9, $r9, $r10
            \\ ld.d $r16, $r9, 0 * 8
            \\ vinsgr2vr.d $vr22, $r14, 0
            \\ vinsgr2vr.d $vr22, $r16, 1
            \\ add.d $r9, $r9, $r10
            \\ vfmadd.d $vr22, $vr20, $vr23, $vr22
            \\ addi.d $r12, $r12, -1
            \\ vstelm.d $vr22, $r6, 0, 0
            \\ add.d $r6, $r6, $r10
            \\ vstelm.d $vr22, $r6, 0, 1
            \\ add.d $r6, $r6, $r10
            \\ blt $r0, $r12, 20019b
            \\ .align 3
            \\ 20020:
            \\ andi $r12, $r4, 7
            \\ bge $r0, $r12, 20024f
            \\ .align 3
            \\ 20021:
            \\ beq $r0, $r10, 20022f
            \\ fld.d $f12, $r7, 0 * 8
            \\ fld.d $f14, $r9, 0 * 8
            \\ addi.d $r12, $r12, -1
            \\ fmadd.d $f14, $f12, $f0, $f14
            \\ fst.d $f14, $r9, 0 * 8
            \\ add.d $r7, $r7, $r8
            \\ add.d $r9, $r9, $r10
            \\ blt $r0, $r12, 20021b
            \\ .align 3
            \\ b 20024f
            \\ 20022:
            \\ fld.d $f14, $r9, 0 * 8
            \\ 20023:
            \\ fld.d $f12, $r7, 0 * 8
            \\ addi.d $r12, $r12, -1
            \\ fmadd.d $f14, $f12, $f0, $f14
            \\ add.d $r7, $r7, $r8
            \\ blt $r0, $r12, 20023b
            \\ fst.d $f14, $r9, 0 * 8
            \\ .align 3
            \\ 20024:
            \\ move $r4, $r12
            \\ jirl $r0, $r1, 0x0
            \\ .align 3
            \\
        );
    }
};
