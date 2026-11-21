//! Section 1 of 8 of the text of `dgemv_n`, which `__root__.zig` joins, from OpenBLAS `kernel/power/gemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1, bp2, bp3, bp4, bp5, bp6, bp7, bp8
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
    \\ stw 0, 0 + 208(1)
    \\ stw 0, 4 + 208(1)
    \\ stw 14, 144(1)
    \\ stw 15, 148(1)
    \\ stw 16, 152(1)
    \\ stw 17, 156(1)
    \\ stw 18, 160(1)
    \\ stw 19, 164(1)
    \\ stw 20, 168(1)
    \\ stw 21, 172(1)
    \\ stw 22, 176(1)
    \\ stw 23, 180(1)
    \\ stw 24, 184(1)
    \\ stw 25, 188(1)
    \\ stw 26, 192(1)
    \\ stw 27, 196(1)
    \\ lwz 5, (((0) * 4) + 8) + 224(1)
    \\ lwz 27, (((1) * 4) + 8) + 224(1)
    \\ stfd 1, 200(1)
    \\ .endif
    \\ .if \bp1
    \\ std 0, 264(1)
    \\ std 14, 144(1)
    \\ std 15, 152(1)
    \\ std 16, 160(1)
    \\ std 17, 168(1)
    \\ std 18, 176(1)
    \\ std 19, 184(1)
    \\ std 20, 192(1)
    \\ std 21, 200(1)
    \\ std 22, 208(1)
    \\ std 23, 216(1)
    \\ std 24, 224(1)
    \\ std 25, 232(1)
    \\ std 26, 240(1)
    \\ std 27, 248(1)
    \\ ld 5, (((0) * 8) + 96) + 280(1)
    \\ ld 6, (((1) * 8) + 96) + 280(1)
    \\ ld 27, (((2) * 8) + 96) + 280(1)
    \\ stfd 1, 256(1)
    \\ .endif
    \\ fmr 31, 1
    \\ .if \bp0
    \\ slwi 22, 7, 3 + 3
    \\ slwi 7, 7, 3
    \\ slwi 9, 9, 3
    \\ slwi 5, 5, 3
    \\ .endif
    \\ .if \bp1
    \\ sldi 22, 8, 3 + 3
    \\ sldi 8, 8, 3
    \\ sldi 10, 10, 3
    \\ sldi 6, 6, 3
    \\ .endif
    \\ .if \bp2
    \\ li 24, 16 * 8
    \\ .endif
    \\ .if \bp3
    \\ li 24, 24 * 8
    \\ .endif
    \\ .if \bp4
    \\ li 25, 16 * 8
    \\ .endif
    \\ .if \bp5
    \\ li 24, 40 * 8
    \\ li 25, 24 * 8
    \\ .endif
    \\ .if \bp6
    \\ li 24, 96 * 8
    \\ li 25, 40 * 8
    \\ .endif
    \\ .if \bp0
    \\ cmpwi 0, 3, 0
    \\ .endif
    \\ .if \bp1
    \\ cmpdi 0, 3, 0
    \\ .endif
    \\ ble- 20036f
    \\ .if \bp0
    \\ cmpwi 0, 4, 0
    \\ .endif
    \\ .if \bp1
    \\ cmpdi 0, 4, 0
    \\ .endif
    \\ ble- 20036f
    \\ .if \bp0
    \\ mr 26, 10
    \\ lfd 0, 208(1)
    \\ cmpi 0, 0, 5, 8
    \\ .endif
    \\ .if \bp1
    \\ mr 26, 5
    \\ lfd 0, 264(1)
    \\ cmpi 0, 0, 6, 8
    \\ .endif
    \\ beq 20001f
    \\ mr 26, 27
    \\ mr 23, 27
    \\ addi 0, 3, 7
    \\ .if \bp0
    \\ srawi. 0, 0, 3
    \\ .endif
    \\ .if \bp1
    \\ sradi. 0, 0, 3
    \\ .endif
    \\ mtspr 9, 0
    \\ .align 4
    \\ 20000:
    \\ stfd 0, 0 * 8(23)
    \\ stfd 0, 1 * 8(23)
    \\ stfd 0, 2 * 8(23)
    \\ stfd 0, 3 * 8(23)
    \\ stfd 0, 4 * 8(23)
    \\ stfd 0, 5 * 8(23)
    \\ stfd 0, 6 * 8(23)
    \\ stfd 0, 7 * 8(23)
    \\ addi 23, 23, 8 * 8
    \\ bdnz 20000b
    \\ .align 4
    \\ 20001:
    \\ .if \bp0
    \\ srawi. 12, 4, 3
    \\ .endif
    \\ .if \bp1
    \\ sradi. 12, 4, 3
    \\ .endif
    \\ ble 20010f
    \\
;
