//! Section 9 of 12 of the text of `cgemv_r`, which `__root__.zig` joins, from OpenBLAS `kernel/power/zgemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .align 4
    \\ 20030:
    \\ .if \bp0
    \\ lfd 14, 208(1)
    \\ lfd 15, 216(1)
    \\ lfs 24, 0 * 4(8)
    \\ lfs 25, 1 * 4(8)
    \\ add 8, 8, 9
    \\ lfs 26, 0 * 4(8)
    \\ lfs 27, 1 * 4(8)
    \\ add 8, 8, 9
    \\ .endif
    \\ .if \bp1
    \\ lfd 14, 256(1)
    \\ lfd 15, 264(1)
    \\ lfs 24, 0 * 4(10)
    \\ lfs 25, 1 * 4(10)
    \\ add 10, 10, 5
    \\ lfs 26, 0 * 4(10)
    \\ lfs 27, 1 * 4(10)
    \\ add 10, 10, 5
    \\ .endif
    \\ fmuls 16, 14, 24
    \\ fmuls 17, 15, 24
    \\ fmuls 18, 14, 26
    \\ fmuls 19, 15, 26
    \\ fnmsubs 16, 15, 25, 16
    \\ fmadds 17, 14, 25, 17
    \\ fnmsubs 18, 15, 27, 18
    \\ fmadds 19, 14, 27, 19
    \\ .if \bp0
    \\ mr 14, 6
    \\ add 15, 6, 7
    \\ add 6, 15, 7
    \\ mr 19, 10
    \\ mr 20, 10
    \\ srawi. 0, 3, 3
    \\ .endif
    \\ .if \bp1
    \\ mr 14, 8
    \\ add 15, 8, 9
    \\ add 8, 15, 9
    \\ mr 19, 6
    \\ mr 20, 6
    \\ sradi. 0, 3, 3
    \\ .endif
    \\ mtspr 9, 0
    \\ ble 20033f
    \\ .align 4
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
    \\ lfs 8, 0 * 4(19)
    \\ lfs 9, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfs 10, 0 * 4(19)
    \\ lfs 11, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfs 12, 0 * 4(19)
    \\ lfs 13, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfs 14, 0 * 4(19)
    \\ lfs 15, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ bdz 20032f
    \\ .align 4
    \\ 20031:
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
    \\ fmadds 0, 17, 25, 0
    \\ fnmsubs 1, 16, 25, 1
    \\ fmadds 2, 17, 27, 2
    \\ fnmsubs 3, 16, 27, 3
    \\ fmadds 4, 17, 29, 4
    \\ fnmsubs 5, 16, 29, 5
    \\ fmadds 6, 17, 31, 6
    \\ fnmsubs 7, 16, 31, 7
    \\ lfs 25, 9 * 4(14)
    \\ lfs 27, 11 * 4(14)
    \\ lfs 29, 13 * 4(14)
    \\ lfs 31, 15 * 4(14)
    \\ addi 14, 14, 16 * 4
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
    \\ fmadds 8, 17, 25, 8
    \\ fnmsubs 9, 16, 25, 9
    \\ fmadds 10, 17, 27, 10
    \\ fnmsubs 11, 16, 27, 11
    \\ fmadds 12, 17, 29, 12
    \\ fnmsubs 13, 16, 29, 13
    \\ fmadds 14, 17, 31, 14
    \\ fnmsubs 15, 16, 31, 15
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
    \\ fmadds 0, 19, 25, 0
    \\ fnmsubs 1, 18, 25, 1
    \\ fmadds 2, 19, 27, 2
    \\ fnmsubs 3, 18, 27, 3
    \\ stfs 0, 0 * 4(20)
    \\ nop
    \\ stfs 1, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ lfs 0, 0 * 4(19)
    \\ nop
    \\ lfs 1, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
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
    \\ lfs 2, 0 * 4(19)
    \\ nop
    \\ lfs 3, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ fmadds 4, 19, 29, 4
    \\ fnmsubs 5, 18, 29, 5
    \\ fmadds 6, 19, 31, 6
    \\ fnmsubs 7, 18, 31, 7
    \\ lfs 25, 9 * 4(15)
    \\ lfs 27, 11 * 4(15)
    \\ lfs 29, 13 * 4(15)
    \\ lfs 31, 15 * 4(15)
    \\ addi 15, 15, 16 * 4
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
    \\ stfs 4, 0 * 4(20)
    \\ nop
    \\ stfs 5, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ lfs 4, 0 * 4(19)
    \\ nop
    \\ lfs 5, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
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
    \\ lfs 6, 0 * 4(19)
    \\ nop
    \\ lfs 7, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ fmadds 8, 18, 24, 8
    \\ fmadds 9, 19, 24, 9
    \\ fmadds 10, 18, 26, 10
    \\ fmadds 11, 19, 26, 11
    \\ fmadds 12, 18, 28, 12
    \\ fmadds 13, 19, 28, 13
    \\ fmadds 14, 18, 30, 14
    \\ fmadds 15, 19, 30, 15
    \\ lfs 24, 0 * 4(14)
    \\ lfs 26, 2 * 4(14)
    \\ lfs 28, 4 * 4(14)
    \\ lfs 30, 6 * 4(14)
    \\ fmadds 8, 19, 25, 8
    \\ fnmsubs 9, 18, 25, 9
    \\ fmadds 10, 19, 27, 10
    \\ fnmsubs 11, 18, 27, 11
    \\ stfs 8, 0 * 4(20)
    \\ nop
    \\ stfs 9, 1 * 4(20)
    \\ .if \bp0
    \\ add 20, 20, 5
    \\ .endif
    \\ .if \bp1
    \\ add 20, 20, 7
    \\ .endif
    \\ lfs 8, 0 * 4(19)
    \\ nop
    \\ lfs 9, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
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
    \\ lfs 10, 0 * 4(19)
    \\ nop
    \\ lfs 11, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ fmadds 12, 19, 29, 12
    \\ fnmsubs 13, 18, 29, 13
    \\ fmadds 14, 19, 31, 14
    \\ fnmsubs 15, 18, 31, 15
    \\ lfs 25, 1 * 4(14)
    \\ lfs 27, 3 * 4(14)
    \\ lfs 29, 5 * 4(14)
    \\ lfs 31, 7 * 4(14)
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
    \\ lfs 12, 0 * 4(19)
    \\ nop
    \\ lfs 13, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ lfs 14, 0 * 4(19)
    \\ nop
    \\ lfs 15, 1 * 4(19)
    \\ .if \bp0
    \\ add 19, 19, 5
    \\ .endif
    \\ .if \bp1
    \\ add 19, 19, 7
    \\ .endif
    \\ .if \bp8
    \\ dcbt 22, 19
    \\ .endif
    \\ .if \bp9
    \\ dcbt 8, 22, 19
    \\ .endif
    \\ .if \bp6
    \\ dcbtst 22, 19
    \\ .endif
    \\ bdnz 20031b
    \\
;
