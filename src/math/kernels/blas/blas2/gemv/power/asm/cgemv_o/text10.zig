//! Section 10 of 12 of the text of `cgemv_o`, which `__root__.zig` joins, from OpenBLAS `kernel/power/zgemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .align 4
    \\ 20032:
    \\ fmadds 0, 16, 24, 0
    \\ fmadds 1, 17, 24, 1
    \\ fmadds 2, 16, 26, 2
    \\ fmadds 3, 17, 26, 3
    \\ fmadds 4, 16, 28, 4
    \\ fmadds 5, 17, 28, 5
    \\ fmadds 6, 16, 30, 6
    \\ fmadds 7, 17, 30, 7
    \\ lfs 24, 8 * 4(14)
    \\ lfs 26, 10 * 4(14)
    \\ lfs 28, 12 * 4(14)
    \\ lfs 30, 14 * 4(14)
    \\ fnmsubs 0, 17, 25, 0
    \\ fmadds 1, 16, 25, 1
    \\ fnmsubs 2, 17, 27, 2
    \\ fmadds 3, 16, 27, 3
    \\ fnmsubs 4, 17, 29, 4
    \\ fmadds 5, 16, 29, 5
    \\ fnmsubs 6, 17, 31, 6
    \\ fmadds 7, 16, 31, 7
    \\ lfs 25, 9 * 4(14)
    \\ lfs 27, 11 * 4(14)
    \\ lfs 29, 13 * 4(14)
    \\ lfs 31, 15 * 4(14)
    \\ fmadds 8, 16, 24, 8
    \\ fmadds 9, 17, 24, 9
    \\ fmadds 10, 16, 26, 10
    \\ fmadds 11, 17, 26, 11
    \\ fmadds 12, 16, 28, 12
    \\ fmadds 13, 17, 28, 13
    \\ fmadds 14, 16, 30, 14
    \\ fmadds 15, 17, 30, 15
    \\ lfs 24, 0 * 4(15)
    \\ lfs 26, 2 * 4(15)
    \\ lfs 28, 4 * 4(15)
    \\ lfs 30, 6 * 4(15)
    \\ fnmsubs 8, 17, 25, 8
    \\ fmadds 9, 16, 25, 9
    \\ fnmsubs 10, 17, 27, 10
    \\ fmadds 11, 16, 27, 11
    \\ fnmsubs 12, 17, 29, 12
    \\ fmadds 13, 16, 29, 13
    \\ fnmsubs 14, 17, 31, 14
    \\ fmadds 15, 16, 31, 15
    \\ lfs 25, 1 * 4(15)
    \\ lfs 27, 3 * 4(15)
    \\ lfs 29, 5 * 4(15)
    \\ lfs 31, 7 * 4(15)
    \\ fmadds 0, 18, 24, 0
    \\ fmadds 1, 19, 24, 1
    \\ fmadds 2, 18, 26, 2
    \\ fmadds 3, 19, 26, 3
    \\ fmadds 4, 18, 28, 4
    \\ fmadds 5, 19, 28, 5
    \\ fmadds 6, 18, 30, 6
    \\ fmadds 7, 19, 30, 7
    \\ lfs 24, 8 * 4(15)
    \\ lfs 26, 10 * 4(15)
    \\ lfs 28, 12 * 4(15)
    \\ lfs 30, 14 * 4(15)
    \\ fnmsubs 0, 19, 25, 0
    \\ fmadds 1, 18, 25, 1
    \\ fnmsubs 2, 19, 27, 2
    \\ fmadds 3, 18, 27, 3
    \\ stfs 0, 0 * 4(20)
    \\ addi 14, 14, 16 * 4
    \\ stfs 1, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfs 2, 0 * 4(20)
    \\ nop
    \\ stfs 3, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ fnmsubs 4, 19, 29, 4
    \\ fmadds 5, 18, 29, 5
    \\ fnmsubs 6, 19, 31, 6
    \\ fmadds 7, 18, 31, 7
    \\ lfs 25, 9 * 4(15)
    \\ lfs 27, 11 * 4(15)
    \\ lfs 29, 13 * 4(15)
    \\ lfs 31, 15 * 4(15)
    \\ stfs 4, 0 * 4(20)
    \\ addi 15, 15, 16 * 4
    \\ stfs 5, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfs 6, 0 * 4(20)
    \\ nop
    \\ stfs 7, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ fmadds 8, 18, 24, 8
    \\ fmadds 9, 19, 24, 9
    \\ fmadds 10, 18, 26, 10
    \\ fmadds 11, 19, 26, 11
    \\ fmadds 12, 18, 28, 12
    \\ fmadds 13, 19, 28, 13
    \\ fmadds 14, 18, 30, 14
    \\ fmadds 15, 19, 30, 15
    \\ fnmsubs 8, 19, 25, 8
    \\ fmadds 9, 18, 25, 9
    \\ fnmsubs 10, 19, 27, 10
    \\ fmadds 11, 18, 27, 11
    \\ stfs 8, 0 * 4(20)
    \\ nop
    \\ stfs 9, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfs 10, 0 * 4(20)
    \\ nop
    \\ stfs 11, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ fnmsubs 12, 19, 29, 12
    \\ fmadds 13, 18, 29, 13
    \\ fnmsubs 14, 19, 31, 14
    \\ fmadds 15, 18, 31, 15
    \\ stfs 12, 0 * 4(20)
    \\ nop
    \\ stfs 13, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfs 14, 0 * 4(20)
    \\ nop
    \\ stfs 15, 1 * 4(20)
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
    \\ lfs 0, 0 * 4(19)
    \\ lfs 1, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfs 2, 0 * 4(19)
    \\ lfs 3, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfs 24, 0 * 4(14)
    \\ lfs 26, 2 * 4(14)
    \\ lfs 28, 4 * 4(14)
    \\ lfs 30, 6 * 4(14)
    \\ lfs 4, 0 * 4(19)
    \\ lfs 5, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfs 6, 0 * 4(19)
    \\ lfs 7, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfs 25, 1 * 4(14)
    \\ lfs 27, 3 * 4(14)
    \\ lfs 29, 5 * 4(14)
    \\ lfs 31, 7 * 4(14)
    \\ fmadds 0, 16, 24, 0
    \\ fmadds 1, 17, 24, 1
    \\ fmadds 2, 16, 26, 2
    \\ fmadds 3, 17, 26, 3
    \\ fmadds 4, 16, 28, 4
    \\ fmadds 5, 17, 28, 5
    \\ fmadds 6, 16, 30, 6
    \\ fmadds 7, 17, 30, 7
    \\ lfs 24, 0 * 4(15)
    \\ lfs 26, 2 * 4(15)
    \\ lfs 28, 4 * 4(15)
    \\ lfs 30, 6 * 4(15)
    \\ fnmsubs 0, 17, 25, 0
    \\ fmadds 1, 16, 25, 1
    \\ fnmsubs 2, 17, 27, 2
    \\ fmadds 3, 16, 27, 3
    \\ fnmsubs 4, 17, 29, 4
    \\ fmadds 5, 16, 29, 5
    \\ fnmsubs 6, 17, 31, 6
    \\ fmadds 7, 16, 31, 7
    \\ lfs 25, 1 * 4(15)
    \\ lfs 27, 3 * 4(15)
    \\ lfs 29, 5 * 4(15)
    \\ lfs 31, 7 * 4(15)
    \\ fmadds 0, 18, 24, 0
    \\ fmadds 1, 19, 24, 1
    \\ fmadds 2, 18, 26, 2
    \\ fmadds 3, 19, 26, 3
    \\ fmadds 4, 18, 28, 4
    \\ fmadds 5, 19, 28, 5
    \\ fmadds 6, 18, 30, 6
    \\ fmadds 7, 19, 30, 7
    \\ fnmsubs 0, 19, 25, 0
    \\ fmadds 1, 18, 25, 1
    \\ fnmsubs 2, 19, 27, 2
    \\ fmadds 3, 18, 27, 3
    \\ stfs 0, 0 * 4(20)
    \\ addi 14, 14, 8 * 4
    \\ stfs 1, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfs 2, 0 * 4(20)
    \\ addi 15, 15, 8 * 4
    \\ stfs 3, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ fnmsubs 4, 19, 29, 4
    \\ fmadds 5, 18, 29, 5
    \\ fnmsubs 6, 19, 31, 6
    \\ fmadds 7, 18, 31, 7
    \\ stfs 4, 0 * 4(20)
    \\ nop
    \\ stfs 5, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfs 6, 0 * 4(20)
    \\ nop
    \\ stfs 7, 1 * 4(20)
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
    \\ lfs 24, 0 * 4(14)
    \\ lfs 25, 1 * 4(14)
    \\ lfs 26, 2 * 4(14)
    \\ lfs 27, 3 * 4(14)
    \\ lfs 0, 0 * 4(19)
    \\ lfs 1, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfs 2, 0 * 4(19)
    \\ lfs 3, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfs 28, 0 * 4(15)
    \\ lfs 29, 1 * 4(15)
    \\ lfs 30, 2 * 4(15)
    \\ lfs 31, 3 * 4(15)
    \\ fmadds 0, 16, 24, 0
    \\ fmadds 1, 17, 24, 1
    \\ fmadds 2, 16, 26, 2
    \\ fmadds 3, 17, 26, 3
    \\ fnmsubs 0, 17, 25, 0
    \\ fmadds 1, 16, 25, 1
    \\ fnmsubs 2, 17, 27, 2
    \\ fmadds 3, 16, 27, 3
    \\ fmadds 0, 18, 28, 0
    \\ fmadds 1, 19, 28, 1
    \\ fmadds 2, 18, 30, 2
    \\ fmadds 3, 19, 30, 3
    \\ fnmsubs 0, 19, 29, 0
    \\ fmadds 1, 18, 29, 1
    \\ fnmsubs 2, 19, 31, 2
    \\ fmadds 3, 18, 31, 3
    \\ stfs 0, 0 * 4(20)
    \\ addi 14, 14, 4 * 4
    \\ stfs 1, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ stfs 2, 0 * 4(20)
    \\ addi 15, 15, 4 * 4
    \\ stfs 3, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\
;
