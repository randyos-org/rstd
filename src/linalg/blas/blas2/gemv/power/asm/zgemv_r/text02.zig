//! Section 2 of 12 of the text of `zgemv_r`, which `__root__.zig` joins, from OpenBLAS `kernel/power/zgemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .align 4
    \\ 20001:
    \\ fmadd 0, 16, 24, 0
    \\ fmadd 1, 17, 24, 1
    \\ fmadd 2, 16, 26, 2
    \\ fmadd 3, 17, 26, 3
    \\ fmadd 4, 16, 28, 4
    \\ fmadd 5, 17, 28, 5
    \\ fmadd 6, 16, 30, 6
    \\ fmadd 7, 17, 30, 7
    \\ lfd 24, 8 * 8(14)
    \\ lfd 26, 10 * 8(14)
    \\ lfd 28, 12 * 8(14)
    \\ lfd 30, 14 * 8(14)
    \\ fmadd 0, 17, 25, 0
    \\ fnmsub 1, 16, 25, 1
    \\ fmadd 2, 17, 27, 2
    \\ fnmsub 3, 16, 27, 3
    \\ fmadd 4, 17, 29, 4
    \\ fnmsub 5, 16, 29, 5
    \\ fmadd 6, 17, 31, 6
    \\ fnmsub 7, 16, 31, 7
    \\ lfd 25, 9 * 8(14)
    \\ lfd 27, 11 * 8(14)
    \\ lfd 29, 13 * 8(14)
    \\ lfd 31, 15 * 8(14)
    \\ addi 14, 14, 16 * 8
    \\ nop
    \\ .if \bp8
    \\ dcbt 21, 14
    \\ .endif
    \\ .if \bp9
    \\ dcbt 8, 21, 14
    \\ .endif
    \\ .if \bp6
    \\ dcbtst 21, 14
    \\ .endif
    \\ nop
    \\ fmadd 8, 16, 24, 8
    \\ fmadd 9, 17, 24, 9
    \\ fmadd 10, 16, 26, 10
    \\ fmadd 11, 17, 26, 11
    \\ fmadd 12, 16, 28, 12
    \\ fmadd 13, 17, 28, 13
    \\ fmadd 14, 16, 30, 14
    \\ fmadd 15, 17, 30, 15
    \\ lfd 24, 0 * 8(15)
    \\ lfd 26, 2 * 8(15)
    \\ lfd 28, 4 * 8(15)
    \\ lfd 30, 6 * 8(15)
    \\ fmadd 8, 17, 25, 8
    \\ fnmsub 9, 16, 25, 9
    \\ fmadd 10, 17, 27, 10
    \\ fnmsub 11, 16, 27, 11
    \\ fmadd 12, 17, 29, 12
    \\ fnmsub 13, 16, 29, 13
    \\ fmadd 14, 17, 31, 14
    \\ fnmsub 15, 16, 31, 15
    \\ lfd 25, 1 * 8(15)
    \\ lfd 27, 3 * 8(15)
    \\ lfd 29, 5 * 8(15)
    \\ lfd 31, 7 * 8(15)
    \\ fmadd 0, 18, 24, 0
    \\ fmadd 1, 19, 24, 1
    \\ fmadd 2, 18, 26, 2
    \\ fmadd 3, 19, 26, 3
    \\ fmadd 4, 18, 28, 4
    \\ fmadd 5, 19, 28, 5
    \\ fmadd 6, 18, 30, 6
    \\ fmadd 7, 19, 30, 7
    \\ lfd 24, 8 * 8(15)
    \\ lfd 26, 10 * 8(15)
    \\ lfd 28, 12 * 8(15)
    \\ lfd 30, 14 * 8(15)
    \\ fmadd 0, 19, 25, 0
    \\ fnmsub 1, 18, 25, 1
    \\ fmadd 2, 19, 27, 2
    \\ fnmsub 3, 18, 27, 3
    \\ fmadd 4, 19, 29, 4
    \\ fnmsub 5, 18, 29, 5
    \\ fmadd 6, 19, 31, 6
    \\ fnmsub 7, 18, 31, 7
    \\ lfd 25, 9 * 8(15)
    \\ lfd 27, 11 * 8(15)
    \\ lfd 29, 13 * 8(15)
    \\ lfd 31, 15 * 8(15)
    \\ addi 15, 15, 16 * 8
    \\ nop
    \\ .if \bp8
    \\ dcbt 21, 15
    \\ .endif
    \\ .if \bp9
    \\ dcbt 8, 21, 15
    \\ .endif
    \\ .if \bp6
    \\ dcbtst 21, 15
    \\ .endif
    \\ nop
    \\ fmadd 8, 18, 24, 8
    \\ fmadd 9, 19, 24, 9
    \\ fmadd 10, 18, 26, 10
    \\ fmadd 11, 19, 26, 11
    \\ fmadd 12, 18, 28, 12
    \\ fmadd 13, 19, 28, 13
    \\ fmadd 14, 18, 30, 14
    \\ fmadd 15, 19, 30, 15
    \\ lfd 24, 0 * 8(16)
    \\ lfd 26, 2 * 8(16)
    \\ lfd 28, 4 * 8(16)
    \\ lfd 30, 6 * 8(16)
    \\ fmadd 8, 19, 25, 8
    \\ fnmsub 9, 18, 25, 9
    \\ fmadd 10, 19, 27, 10
    \\ fnmsub 11, 18, 27, 11
    \\ fmadd 12, 19, 29, 12
    \\ fnmsub 13, 18, 29, 13
    \\ fmadd 14, 19, 31, 14
    \\ fnmsub 15, 18, 31, 15
    \\ lfd 25, 1 * 8(16)
    \\ lfd 27, 3 * 8(16)
    \\ lfd 29, 5 * 8(16)
    \\ lfd 31, 7 * 8(16)
    \\ fmadd 0, 20, 24, 0
    \\ fmadd 1, 21, 24, 1
    \\ fmadd 2, 20, 26, 2
    \\ fmadd 3, 21, 26, 3
    \\ fmadd 4, 20, 28, 4
    \\ fmadd 5, 21, 28, 5
    \\ fmadd 6, 20, 30, 6
    \\ fmadd 7, 21, 30, 7
    \\ lfd 24, 8 * 8(16)
    \\ lfd 26, 10 * 8(16)
    \\ lfd 28, 12 * 8(16)
    \\ lfd 30, 14 * 8(16)
    \\ fmadd 0, 21, 25, 0
    \\ fnmsub 1, 20, 25, 1
    \\ fmadd 2, 21, 27, 2
    \\ fnmsub 3, 20, 27, 3
    \\ fmadd 4, 21, 29, 4
    \\ fnmsub 5, 20, 29, 5
    \\ fmadd 6, 21, 31, 6
    \\ fnmsub 7, 20, 31, 7
    \\ lfd 25, 9 * 8(16)
    \\ lfd 27, 11 * 8(16)
    \\ lfd 29, 13 * 8(16)
    \\ lfd 31, 15 * 8(16)
    \\ addi 16, 16, 16 * 8
    \\ nop
    \\ .if \bp8
    \\ dcbt 21, 16
    \\ .endif
    \\ .if \bp9
    \\ dcbt 8, 21, 16
    \\ .endif
    \\ .if \bp6
    \\ dcbtst 21, 16
    \\ .endif
    \\ nop
    \\ fmadd 8, 20, 24, 8
    \\ fmadd 9, 21, 24, 9
    \\ fmadd 10, 20, 26, 10
    \\ fmadd 11, 21, 26, 11
    \\ fmadd 12, 20, 28, 12
    \\ fmadd 13, 21, 28, 13
    \\ fmadd 14, 20, 30, 14
    \\ fmadd 15, 21, 30, 15
    \\ lfd 24, 0 * 8(17)
    \\ lfd 26, 2 * 8(17)
    \\ lfd 28, 4 * 8(17)
    \\ lfd 30, 6 * 8(17)
    \\ fmadd 8, 21, 25, 8
    \\ fnmsub 9, 20, 25, 9
    \\ fmadd 10, 21, 27, 10
    \\ fnmsub 11, 20, 27, 11
    \\ fmadd 12, 21, 29, 12
    \\ fnmsub 13, 20, 29, 13
    \\ fmadd 14, 21, 31, 14
    \\ fnmsub 15, 20, 31, 15
    \\ lfd 25, 1 * 8(17)
    \\ lfd 27, 3 * 8(17)
    \\ lfd 29, 5 * 8(17)
    \\ lfd 31, 7 * 8(17)
    \\ fmadd 0, 22, 24, 0
    \\ fmadd 1, 23, 24, 1
    \\ fmadd 2, 22, 26, 2
    \\ fmadd 3, 23, 26, 3
    \\ fmadd 4, 22, 28, 4
    \\ fmadd 5, 23, 28, 5
    \\ fmadd 6, 22, 30, 6
    \\ fmadd 7, 23, 30, 7
    \\ lfd 24, 8 * 8(17)
    \\ lfd 26, 10 * 8(17)
    \\ lfd 28, 12 * 8(17)
    \\ lfd 30, 14 * 8(17)
    \\ fmadd 0, 23, 25, 0
    \\ fnmsub 1, 22, 25, 1
    \\ fmadd 2, 23, 27, 2
    \\ fnmsub 3, 22, 27, 3
    \\ stfd 0, 0 * 8(20)
    \\ stfd 1, 1 * 8(20)
    \\ stfd 2, 2 * 8(20)
    \\ stfd 3, 3 * 8(20)
    \\ lfd 0, 0 * 8(19)
    \\ lfd 1, 1 * 8(19)
    \\ lfd 2, 2 * 8(19)
    \\ lfd 3, 3 * 8(19)
    \\ fmadd 4, 23, 29, 4
    \\ fnmsub 5, 22, 29, 5
    \\ fmadd 6, 23, 31, 6
    \\ fnmsub 7, 22, 31, 7
    \\ lfd 25, 9 * 8(17)
    \\ lfd 27, 11 * 8(17)
    \\ lfd 29, 13 * 8(17)
    \\ lfd 31, 15 * 8(17)
    \\ addi 17, 17, 16 * 8
    \\ nop
    \\ .if \bp8
    \\ dcbt 21, 17
    \\ .endif
    \\ .if \bp9
    \\ dcbt 8, 21, 17
    \\ .endif
    \\ .if \bp6
    \\ dcbtst 21, 17
    \\ .endif
    \\ nop
    \\ stfd 4, 4 * 8(20)
    \\ stfd 5, 5 * 8(20)
    \\ stfd 6, 6 * 8(20)
    \\ stfd 7, 7 * 8(20)
    \\ lfd 4, 4 * 8(19)
    \\ lfd 5, 5 * 8(19)
    \\ lfd 6, 6 * 8(19)
    \\ lfd 7, 7 * 8(19)
    \\ fmadd 8, 22, 24, 8
    \\ fmadd 9, 23, 24, 9
    \\ fmadd 10, 22, 26, 10
    \\ fmadd 11, 23, 26, 11
    \\ fmadd 12, 22, 28, 12
    \\ fmadd 13, 23, 28, 13
    \\ fmadd 14, 22, 30, 14
    \\ fmadd 15, 23, 30, 15
    \\ lfd 24, 0 * 8(14)
    \\ lfd 26, 2 * 8(14)
    \\ lfd 28, 4 * 8(14)
    \\ lfd 30, 6 * 8(14)
    \\ fmadd 8, 23, 25, 8
    \\ fnmsub 9, 22, 25, 9
    \\ fmadd 10, 23, 27, 10
    \\ fnmsub 11, 22, 27, 11
    \\ stfd 8, 8 * 8(20)
    \\ stfd 9, 9 * 8(20)
    \\ stfd 10, 10 * 8(20)
    \\ stfd 11, 11 * 8(20)
    \\ lfd 8, 8 * 8(19)
    \\ lfd 9, 9 * 8(19)
    \\ lfd 10, 10 * 8(19)
    \\ lfd 11, 11 * 8(19)
    \\ fmadd 12, 23, 29, 12
    \\ fnmsub 13, 22, 29, 13
    \\ fmadd 14, 23, 31, 14
    \\ fnmsub 15, 22, 31, 15
    \\ lfd 25, 1 * 8(14)
    \\ lfd 27, 3 * 8(14)
    \\ lfd 29, 5 * 8(14)
    \\ lfd 31, 7 * 8(14)
    \\ stfd 12, 12 * 8(20)
    \\ stfd 13, 13 * 8(20)
    \\ stfd 14, 14 * 8(20)
    \\ stfd 15, 15 * 8(20)
    \\ lfd 12, 12 * 8(19)
    \\ lfd 13, 13 * 8(19)
    \\ lfd 14, 14 * 8(19)
    \\ lfd 15, 15 * 8(19)
    \\ addi 20, 20, 16 * 8
    \\ addi 19, 19, 16 * 8
    \\ .if \bp8
    \\ dcbt 22, 19
    \\ .endif
    \\ .if \bp9
    \\ dcbt 8, 22, 19
    \\ .endif
    \\ .if \bp6
    \\ dcbtst 22, 19
    \\ .endif
    \\ bdnz 20001b
    \\
;
