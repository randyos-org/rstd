//! Section 1 of 12 of the text of `zgemv_n`, which `__root__.zig` joins, from OpenBLAS `kernel/power/zgemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1, bp2, bp3, bp4, bp5, bp6, bp7, bp8, bp9
    \\ .align 6
    \\ .if \bp0
    \\ addi 1, 1, -224
    \\ .endif
    \\ .if \bp1
    \\ addi 1, 1, -280
    \\ .endif
    \\ li 0, 0
    \\ stfd 14, 0(1)
    \\ stfd 15, 8(1)
    \\ stfd 16, 16(1)
    \\ stfd 17, 24(1)
    \\ stfd 18, 32(1)
    \\ stfd 19, 40(1)
    \\ stfd 20, 48(1)
    \\ stfd 21, 56(1)
    \\ stfd 22, 64(1)
    \\ stfd 23, 72(1)
    \\ stfd 24, 80(1)
    \\ stfd 25, 88(1)
    \\ stfd 26, 96(1)
    \\ stfd 27, 104(1)
    \\ stfd 28, 112(1)
    \\ stfd 29, 120(1)
    \\ stfd 30, 128(1)
    \\ stfd 31, 136(1)
    \\ .if \bp0
    \\ stw 14, 144(1)
    \\ stw 15, 148(1)
    \\ stw 16, 152(1)
    \\ stw 17, 156(1)
    \\ stw 18, 160(1)
    \\ stw 19, 164(1)
    \\ stw 20, 168(1)
    \\ stw 21, 172(1)
    \\ stw 22, 176(1)
    \\ lwz 5, (((0) * 4) + 8) + 224(1)
    \\ stfd 1, 208(1)
    \\ stfd 2, 216(1)
    \\ slwi 18, 7, 4 + 2
    \\ slwi 7, 7, 4
    \\ slwi 9, 9, 4
    \\ slwi 5, 5, 4
    \\ .endif
    \\ .if \bp1
    \\ std 14, 144(1)
    \\ std 15, 152(1)
    \\ std 16, 160(1)
    \\ std 17, 168(1)
    \\ std 18, 176(1)
    \\ std 19, 184(1)
    \\ std 20, 192(1)
    \\ std 21, 200(1)
    \\ std 22, 208(1)
    \\ ld 5, (((0) * 8) + 96) + 280(1)
    \\ ld 6, (((1) * 8) + 96) + 280(1)
    \\ ld 7, (((2) * 8) + 96) + 280(1)
    \\ stfd 1, 256(1)
    \\ stfd 2, 264(1)
    \\ sldi 18, 9, 4 + 2
    \\ sldi 9, 9, 4
    \\ sldi 5, 5, 4
    \\ sldi 7, 7, 4
    \\ .endif
    \\ .if \bp2
    \\ li 21, 56 * 8
    \\ .endif
    \\ .if \bp3
    \\ li 21, 34 * 8
    \\ .endif
    \\ .if \bp4
    \\ li 22, 16 * 8
    \\ .endif
    \\ .if \bp5
    \\ li 21, 40 * 8
    \\ .endif
    \\ .if \bp6
    \\ li 21, 24 * 8
    \\ .endif
    \\ .if \bp7
    \\ li 22, 24 * 8
    \\ .endif
    \\ .if \bp0
    \\ cmpwi 0, 3, 0
    \\ .endif
    \\ .if \bp1
    \\ cmpdi 0, 3, 0
    \\ .endif
    \\ ble- 20043f
    \\ .if \bp0
    \\ cmpwi 0, 4, 0
    \\ .endif
    \\ .if \bp1
    \\ cmpdi 0, 4, 0
    \\ .endif
    \\ ble- 20043f
    \\ .if \bp0
    \\ cmpi 0, 0, 5, 2 * 8
    \\ .endif
    \\ .if \bp1
    \\ cmpi 0, 0, 7, 2 * 8
    \\ .endif
    \\ bne 20021f
    \\ .if \bp0
    \\ srawi. 12, 4, 2
    \\ .endif
    \\ .if \bp1
    \\ sradi. 12, 4, 2
    \\ .endif
    \\ ble 20007f
    \\ .align 4
    \\ 20000:
    \\ .if \bp0
    \\ lfd 14, 208(1)
    \\ lfd 15, 216(1)
    \\ lfd 24, 0 * 8(8)
    \\ lfd 25, 1 * 8(8)
    \\ add 8, 8, 9
    \\ lfd 26, 0 * 8(8)
    \\ lfd 27, 1 * 8(8)
    \\ add 8, 8, 9
    \\ lfd 28, 0 * 8(8)
    \\ lfd 29, 1 * 8(8)
    \\ add 8, 8, 9
    \\ lfd 30, 0 * 8(8)
    \\ lfd 31, 1 * 8(8)
    \\ add 8, 8, 9
    \\ .endif
    \\ .if \bp1
    \\ lfd 14, 256(1)
    \\ lfd 15, 264(1)
    \\ lfd 24, 0 * 8(10)
    \\ lfd 25, 1 * 8(10)
    \\ add 10, 10, 5
    \\ lfd 26, 0 * 8(10)
    \\ lfd 27, 1 * 8(10)
    \\ add 10, 10, 5
    \\ lfd 28, 0 * 8(10)
    \\ lfd 29, 1 * 8(10)
    \\ add 10, 10, 5
    \\ lfd 30, 0 * 8(10)
    \\ lfd 31, 1 * 8(10)
    \\ add 10, 10, 5
    \\ .endif
    \\ fmul 16, 14, 24
    \\ fmul 17, 15, 24
    \\ fmul 18, 14, 26
    \\ fmul 19, 15, 26
    \\ fmul 20, 14, 28
    \\ fmul 21, 15, 28
    \\ fmul 22, 14, 30
    \\ fmul 23, 15, 30
    \\ fnmsub 16, 15, 25, 16
    \\ fmadd 17, 14, 25, 17
    \\ fnmsub 18, 15, 27, 18
    \\ fmadd 19, 14, 27, 19
    \\ fnmsub 20, 15, 29, 20
    \\ fmadd 21, 14, 29, 21
    \\ fnmsub 22, 15, 31, 22
    \\ fmadd 23, 14, 31, 23
    \\ .if \bp0
    \\ mr 14, 6
    \\ add 15, 6, 7
    \\ add 16, 15, 7
    \\ add 17, 16, 7
    \\ add 6, 17, 7
    \\ mr 19, 10
    \\ mr 20, 10
    \\ srawi. 0, 3, 3
    \\ .endif
    \\ .if \bp1
    \\ mr 14, 8
    \\ add 15, 8, 9
    \\ add 16, 15, 9
    \\ add 17, 16, 9
    \\ add 8, 17, 9
    \\ mr 19, 6
    \\ mr 20, 6
    \\ sradi. 0, 3, 3
    \\ .endif
    \\ mtspr 9, 0
    \\ ble 20003f
    \\ .align 4
    \\ lfd 24, 0 * 8(14)
    \\ lfd 25, 1 * 8(14)
    \\ lfd 26, 2 * 8(14)
    \\ lfd 27, 3 * 8(14)
    \\ lfd 0, 0 * 8(19)
    \\ lfd 1, 1 * 8(19)
    \\ lfd 2, 2 * 8(19)
    \\ lfd 3, 3 * 8(19)
    \\ lfd 28, 4 * 8(14)
    \\ lfd 29, 5 * 8(14)
    \\ lfd 30, 6 * 8(14)
    \\ lfd 31, 7 * 8(14)
    \\ lfd 4, 4 * 8(19)
    \\ lfd 5, 5 * 8(19)
    \\ lfd 6, 6 * 8(19)
    \\ lfd 7, 7 * 8(19)
    \\ lfd 8, 8 * 8(19)
    \\ lfd 9, 9 * 8(19)
    \\ lfd 10, 10 * 8(19)
    \\ lfd 11, 11 * 8(19)
    \\ lfd 12, 12 * 8(19)
    \\ lfd 13, 13 * 8(19)
    \\ lfd 14, 14 * 8(19)
    \\ lfd 15, 15 * 8(19)
    \\ addi 19, 19, 16 * 8
    \\ bdz 20002f
    \\
;
