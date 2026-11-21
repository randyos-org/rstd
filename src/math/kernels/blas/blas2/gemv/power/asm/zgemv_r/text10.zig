//! Section 10 of 12 of the text of `zgemv_r`, which `__root__.zig` joins, from OpenBLAS `kernel/power/zgemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .align 4
    \\ 20032:
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
    \\ stfd 0, 0 * 8(20)
    \\ addi 14, 14, 16 * 8
    \\ stfd 1, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfd 2, 0 * 8(20)
    \\ nop
    \\ stfd 3, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ fmadd 4, 19, 29, 4
    \\ fnmsub 5, 18, 29, 5
    \\ fmadd 6, 19, 31, 6
    \\ fnmsub 7, 18, 31, 7
    \\ lfd 25, 9 * 8(15)
    \\ lfd 27, 11 * 8(15)
    \\ lfd 29, 13 * 8(15)
    \\ lfd 31, 15 * 8(15)
    \\ stfd 4, 0 * 8(20)
    \\ addi 15, 15, 16 * 8
    \\ stfd 5, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfd 6, 0 * 8(20)
    \\ nop
    \\ stfd 7, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ fmadd 8, 18, 24, 8
    \\ fmadd 9, 19, 24, 9
    \\ fmadd 10, 18, 26, 10
    \\ fmadd 11, 19, 26, 11
    \\ fmadd 12, 18, 28, 12
    \\ fmadd 13, 19, 28, 13
    \\ fmadd 14, 18, 30, 14
    \\ fmadd 15, 19, 30, 15
    \\ fmadd 8, 19, 25, 8
    \\ fnmsub 9, 18, 25, 9
    \\ fmadd 10, 19, 27, 10
    \\ fnmsub 11, 18, 27, 11
    \\ stfd 8, 0 * 8(20)
    \\ nop
    \\ stfd 9, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfd 10, 0 * 8(20)
    \\ nop
    \\ stfd 11, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ fmadd 12, 19, 29, 12
    \\ fnmsub 13, 18, 29, 13
    \\ fmadd 14, 19, 31, 14
    \\ fnmsub 15, 18, 31, 15
    \\ stfd 12, 0 * 8(20)
    \\ nop
    \\ stfd 13, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfd 14, 0 * 8(20)
    \\ nop
    \\ stfd 15, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ .align 4
    \\ 20033:
    \\ andi. 0, 3, 7
    \\ ble 20036f
    \\ andi. 0, 3, 4
    \\ ble 20034f
    \\ lfd 0, 0 * 8(19)
    \\ lfd 1, 1 * 8(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfd 2, 0 * 8(19)
    \\ lfd 3, 1 * 8(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfd 24, 0 * 8(14)
    \\ lfd 26, 2 * 8(14)
    \\ lfd 28, 4 * 8(14)
    \\ lfd 30, 6 * 8(14)
    \\ lfd 4, 0 * 8(19)
    \\ lfd 5, 1 * 8(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfd 6, 0 * 8(19)
    \\ lfd 7, 1 * 8(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfd 25, 1 * 8(14)
    \\ lfd 27, 3 * 8(14)
    \\ lfd 29, 5 * 8(14)
    \\ lfd 31, 7 * 8(14)
    \\ fmadd 0, 16, 24, 0
    \\ fmadd 1, 17, 24, 1
    \\ fmadd 2, 16, 26, 2
    \\ fmadd 3, 17, 26, 3
    \\ fmadd 4, 16, 28, 4
    \\ fmadd 5, 17, 28, 5
    \\ fmadd 6, 16, 30, 6
    \\ fmadd 7, 17, 30, 7
    \\ lfd 24, 0 * 8(15)
    \\ lfd 26, 2 * 8(15)
    \\ lfd 28, 4 * 8(15)
    \\ lfd 30, 6 * 8(15)
    \\ fmadd 0, 17, 25, 0
    \\ fnmsub 1, 16, 25, 1
    \\ fmadd 2, 17, 27, 2
    \\ fnmsub 3, 16, 27, 3
    \\ fmadd 4, 17, 29, 4
    \\ fnmsub 5, 16, 29, 5
    \\ fmadd 6, 17, 31, 6
    \\ fnmsub 7, 16, 31, 7
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
    \\ fmadd 0, 19, 25, 0
    \\ fnmsub 1, 18, 25, 1
    \\ fmadd 2, 19, 27, 2
    \\ fnmsub 3, 18, 27, 3
    \\ stfd 0, 0 * 8(20)
    \\ addi 14, 14, 8 * 8
    \\ stfd 1, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfd 2, 0 * 8(20)
    \\ addi 15, 15, 8 * 8
    \\ stfd 3, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ fmadd 4, 19, 29, 4
    \\ fnmsub 5, 18, 29, 5
    \\ fmadd 6, 19, 31, 6
    \\ fnmsub 7, 18, 31, 7
    \\ stfd 4, 0 * 8(20)
    \\ nop
    \\ stfd 5, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfd 6, 0 * 8(20)
    \\ nop
    \\ stfd 7, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ .align 4
    \\ 20034:
    \\ andi. 0, 3, 2
    \\ ble 20035f
    \\ lfd 24, 0 * 8(14)
    \\ lfd 25, 1 * 8(14)
    \\ lfd 26, 2 * 8(14)
    \\ lfd 27, 3 * 8(14)
    \\ lfd 0, 0 * 8(19)
    \\ lfd 1, 1 * 8(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfd 2, 0 * 8(19)
    \\ lfd 3, 1 * 8(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfd 28, 0 * 8(15)
    \\ lfd 29, 1 * 8(15)
    \\ lfd 30, 2 * 8(15)
    \\ lfd 31, 3 * 8(15)
    \\ fmadd 0, 16, 24, 0
    \\ fmadd 1, 17, 24, 1
    \\ fmadd 2, 16, 26, 2
    \\ fmadd 3, 17, 26, 3
    \\ fmadd 0, 17, 25, 0
    \\ fnmsub 1, 16, 25, 1
    \\ fmadd 2, 17, 27, 2
    \\ fnmsub 3, 16, 27, 3
    \\ fmadd 0, 18, 28, 0
    \\ fmadd 1, 19, 28, 1
    \\ fmadd 2, 18, 30, 2
    \\ fmadd 3, 19, 30, 3
    \\ fmadd 0, 19, 29, 0
    \\ fnmsub 1, 18, 29, 1
    \\ fmadd 2, 19, 31, 2
    \\ fnmsub 3, 18, 31, 3
    \\ stfd 0, 0 * 8(20)
    \\ addi 14, 14, 4 * 8
    \\ stfd 1, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfd 2, 0 * 8(20)
    \\ addi 15, 15, 4 * 8
    \\ stfd 3, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\
;
