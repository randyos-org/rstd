//! Section 8 of 12 of the text of `zgemv_n`, which `__root__.zig` joins, from OpenBLAS `kernel/power/zgemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .align 4
    \\ 20025:
    \\ andi. 0, 3, 7
    \\ ble 20028f
    \\ andi. 0, 3, 4
    \\ ble 20026f
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
    \\ fnmsub 0, 17, 25, 0
    \\ fmadd 1, 16, 25, 1
    \\ fnmsub 2, 17, 27, 2
    \\ fmadd 3, 16, 27, 3
    \\ fnmsub 4, 17, 29, 4
    \\ fmadd 5, 16, 29, 5
    \\ fnmsub 6, 17, 31, 6
    \\ fmadd 7, 16, 31, 7
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
    \\ lfd 24, 0 * 8(16)
    \\ lfd 26, 2 * 8(16)
    \\ lfd 28, 4 * 8(16)
    \\ lfd 30, 6 * 8(16)
    \\ fnmsub 0, 19, 25, 0
    \\ fmadd 1, 18, 25, 1
    \\ fnmsub 2, 19, 27, 2
    \\ fmadd 3, 18, 27, 3
    \\ fnmsub 4, 19, 29, 4
    \\ fmadd 5, 18, 29, 5
    \\ fnmsub 6, 19, 31, 6
    \\ fmadd 7, 18, 31, 7
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
    \\ lfd 24, 0 * 8(17)
    \\ lfd 26, 2 * 8(17)
    \\ lfd 28, 4 * 8(17)
    \\ lfd 30, 6 * 8(17)
    \\ fnmsub 0, 21, 25, 0
    \\ fmadd 1, 20, 25, 1
    \\ fnmsub 2, 21, 27, 2
    \\ fmadd 3, 20, 27, 3
    \\ fnmsub 4, 21, 29, 4
    \\ fmadd 5, 20, 29, 5
    \\ fnmsub 6, 21, 31, 6
    \\ fmadd 7, 20, 31, 7
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
    \\ fnmsub 0, 23, 25, 0
    \\ fmadd 1, 22, 25, 1
    \\ fnmsub 2, 23, 27, 2
    \\ fmadd 3, 22, 27, 3
    \\ fnmsub 4, 23, 29, 4
    \\ fmadd 5, 22, 29, 5
    \\ fnmsub 6, 23, 31, 6
    \\ fmadd 7, 22, 31, 7
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
    \\ stfd 4, 0 * 8(20)
    \\ addi 16, 16, 8 * 8
    \\ stfd 5, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfd 6, 0 * 8(20)
    \\ addi 17, 17, 8 * 8
    \\ stfd 7, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ .align 4
    \\ 20026:
    \\ andi. 0, 3, 2
    \\ ble 20027f
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
    \\ fnmsub 0, 17, 25, 0
    \\ fmadd 1, 16, 25, 1
    \\ fnmsub 2, 17, 27, 2
    \\ fmadd 3, 16, 27, 3
    \\ lfd 24, 0 * 8(16)
    \\ lfd 25, 1 * 8(16)
    \\ lfd 26, 2 * 8(16)
    \\ lfd 27, 3 * 8(16)
    \\ fmadd 0, 18, 28, 0
    \\ fmadd 1, 19, 28, 1
    \\ fmadd 2, 18, 30, 2
    \\ fmadd 3, 19, 30, 3
    \\ fnmsub 0, 19, 29, 0
    \\ fmadd 1, 18, 29, 1
    \\ fnmsub 2, 19, 31, 2
    \\ fmadd 3, 18, 31, 3
    \\ lfd 28, 0 * 8(17)
    \\ lfd 29, 1 * 8(17)
    \\ lfd 30, 2 * 8(17)
    \\ lfd 31, 3 * 8(17)
    \\ fmadd 0, 20, 24, 0
    \\ fmadd 1, 21, 24, 1
    \\ fmadd 2, 20, 26, 2
    \\ fmadd 3, 21, 26, 3
    \\ fnmsub 0, 21, 25, 0
    \\ fmadd 1, 20, 25, 1
    \\ fnmsub 2, 21, 27, 2
    \\ fmadd 3, 20, 27, 3
    \\ fmadd 0, 22, 28, 0
    \\ fmadd 1, 23, 28, 1
    \\ fmadd 2, 22, 30, 2
    \\ fmadd 3, 23, 30, 3
    \\ fnmsub 0, 23, 29, 0
    \\ fmadd 1, 22, 29, 1
    \\ fnmsub 2, 23, 31, 2
    \\ fmadd 3, 22, 31, 3
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
    \\ addi 16, 16, 4 * 8
    \\ addi 17, 17, 4 * 8
    \\ .align 4
    \\ 20027:
    \\ andi. 0, 3, 1
    \\ ble 20028f
    \\ lfd 0, 0 * 8(19)
    \\ lfd 1, 1 * 8(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfd 24, 0 * 8(14)
    \\ lfd 25, 1 * 8(14)
    \\ lfd 26, 0 * 8(15)
    \\ lfd 27, 1 * 8(15)
    \\ lfd 28, 0 * 8(16)
    \\ lfd 29, 1 * 8(16)
    \\ lfd 30, 0 * 8(17)
    \\ lfd 31, 1 * 8(17)
    \\ fmadd 0, 16, 24, 0
    \\ fmadd 1, 17, 24, 1
    \\ fnmsub 0, 17, 25, 0
    \\ fmadd 1, 16, 25, 1
    \\ fmadd 0, 18, 26, 0
    \\ fmadd 1, 19, 26, 1
    \\ fnmsub 0, 19, 27, 0
    \\ fmadd 1, 18, 27, 1
    \\ fmadd 0, 20, 28, 0
    \\ fmadd 1, 21, 28, 1
    \\ fnmsub 0, 21, 29, 0
    \\ fmadd 1, 20, 29, 1
    \\ fmadd 0, 22, 30, 0
    \\ fmadd 1, 23, 30, 1
    \\ fnmsub 0, 23, 31, 0
    \\ fmadd 1, 22, 31, 1
    \\ stfd 0, 0 * 8(20)
    \\ stfd 1, 1 * 8(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ .align 4
    \\ 20028:
    \\ addi 12, 12, -1
    \\ cmpi 0, 0, 12, 0
    \\ bgt 20022b
    \\ .align 4
    \\ 20029:
    \\ andi. 12, 4, 2
    \\ ble 20036f
    \\
;
