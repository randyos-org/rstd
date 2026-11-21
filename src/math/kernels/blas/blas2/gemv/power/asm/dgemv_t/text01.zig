//! Section 1 of 7 of the text of `dgemv_t`, which `__root__.zig` joins, from OpenBLAS `kernel/power/gemv_t.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1, bp2, bp3, bp4, bp5, bp6, bp7, bp8, bp9, bp10
    \\ .align 6
    \\ .if \bp0
    \\ addi 1, 1, -224
    \\ .endif
    \\ .if \bp1
    \\ addi 1, 1, -288
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
    \\ stw 0, 0 + 144(1)
    \\ stw 0, 4 + 144(1)
    \\ .endif
    \\ .if \bp1
    \\ std 0, 144(1)
    \\ .endif
    \\ stfd 1, 152(1)
    \\ .if \bp0
    \\ stw 14, 160(1)
    \\ stw 15, 164(1)
    \\ stw 16, 168(1)
    \\ stw 17, 172(1)
    \\ stw 18, 176(1)
    \\ stw 19, 180(1)
    \\ stw 20, 184(1)
    \\ stw 21, 188(1)
    \\ stw 22, 192(1)
    \\ stw 23, 196(1)
    \\ stw 24, 200(1)
    \\ stw 25, 204(1)
    \\ stw 26, 208(1)
    \\ stw 27, 212(1)
    \\ stw 28, 216(1)
    \\ stw 29, 220(1)
    \\ lwz 5, (((0) * 4) + 8) + 224(1)
    \\ lwz 11, (((1) * 4) + 8) + 224(1)
    \\ mullw 28, 7, 4
    \\ .endif
    \\ .if \bp1
    \\ std 14, 160(1)
    \\ std 15, 168(1)
    \\ std 16, 176(1)
    \\ std 17, 184(1)
    \\ std 18, 192(1)
    \\ std 19, 200(1)
    \\ std 20, 208(1)
    \\ std 21, 216(1)
    \\ std 22, 224(1)
    \\ std 23, 232(1)
    \\ std 24, 240(1)
    \\ std 25, 248(1)
    \\ std 26, 256(1)
    \\ std 27, 264(1)
    \\ std 28, 272(1)
    \\ std 29, 280(1)
    \\ ld 5, (((0) * 8) + 96) + 288(1)
    \\ ld 6, (((1) * 8) + 96) + 288(1)
    \\ ld 11, (((2) * 8) + 96) + 288(1)
    \\ mulld 28, 8, 4
    \\ .endif
    \\ li 12, 2048
    \\ subf 28, 12, 28
    \\ .if \bp0
    \\ slwi 28, 28, 3
    \\ slwi 7, 7, 3
    \\ slwi 9, 9, 3
    \\ slwi 5, 5, 3
    \\ subf 10, 5, 10
    \\ .endif
    \\ .if \bp1
    \\ sldi 28, 28, 3
    \\ sldi 8, 8, 3
    \\ sldi 10, 10, 3
    \\ sldi 6, 6, 3
    \\ subf 5, 6, 5
    \\ .endif
    \\ li 29, 0
    \\ .if \bp0
    \\ addi 6, 6, -8
    \\ .endif
    \\ .if \bp1
    \\ addi 7, 7, -8
    \\ .endif
    \\ .if \bp2
    \\ li 25, 42 * 8
    \\ .endif
    \\ .if \bp3
    \\ li 25, 16 * 8
    \\ .endif
    \\ .if \bp4
    \\ li 25, 48 * 8
    \\ .endif
    \\ .if \bp5
    \\ li 26, 16 * 8
    \\ .endif
    \\ .if \bp6
    \\ li 25, 40 * 8
    \\ .endif
    \\ .if \bp7
    \\ li 25, 96 * 8
    \\ .endif
    \\ .if \bp8
    \\ li 26, 8 * 8
    \\ .endif
    \\ cmpi 0, 0, 3, 0
    \\ ble 20043f
    \\ cmpi 0, 0, 4, 0
    \\ ble 20043f
    \\ .align 4
    \\ 20000:
    \\ subf 22, 29, 3
    \\ .if \bp0
    \\ slwi 0, 29, 3
    \\ .endif
    \\ .if \bp1
    \\ sldi 0, 29, 3
    \\ .endif
    \\ cmpi 0, 0, 22, 2048
    \\ ble+ 20001f
    \\ li 22, 2048
    \\ 20001:
    \\ .if \bp0
    \\ add 12, 8, 0
    \\ cmpi 0, 0, 9, 8
    \\ .endif
    \\ .if \bp1
    \\ add 12, 9, 0
    \\ cmpi 0, 0, 10, 8
    \\ .endif
    \\ beq 20005f
    \\ mr 12, 11
    \\ addi 24, 11, -8
    \\ .if \bp0
    \\ srawi. 0, 22, 3
    \\ .endif
    \\ .if \bp1
    \\ sradi. 0, 22, 3
    \\ .endif
    \\ mtspr 9, 0
    \\ ble 20003f
    \\ .align 4
    \\ 20002:
    \\ .if \bp0
    \\ lfd 0, 0 * 8(8)
    \\ add 8, 8, 9
    \\ lfd 1, 0 * 8(8)
    \\ add 8, 8, 9
    \\ lfd 2, 0 * 8(8)
    \\ add 8, 8, 9
    \\ lfd 3, 0 * 8(8)
    \\ add 8, 8, 9
    \\ lfd 4, 0 * 8(8)
    \\ add 8, 8, 9
    \\ lfd 5, 0 * 8(8)
    \\ add 8, 8, 9
    \\ lfd 6, 0 * 8(8)
    \\ add 8, 8, 9
    \\ lfd 7, 0 * 8(8)
    \\ add 8, 8, 9
    \\ .endif
    \\ .if \bp1
    \\ lfd 0, 0 * 8(9)
    \\ add 9, 9, 10
    \\ lfd 1, 0 * 8(9)
    \\ add 9, 9, 10
    \\ lfd 2, 0 * 8(9)
    \\ add 9, 9, 10
    \\ lfd 3, 0 * 8(9)
    \\ add 9, 9, 10
    \\ lfd 4, 0 * 8(9)
    \\ add 9, 9, 10
    \\ lfd 5, 0 * 8(9)
    \\ add 9, 9, 10
    \\ lfd 6, 0 * 8(9)
    \\ add 9, 9, 10
    \\ lfd 7, 0 * 8(9)
    \\ add 9, 9, 10
    \\ .endif
    \\ stfd 0, 1 * 8(24)
    \\ stfd 1, 2 * 8(24)
    \\ stfd 2, 3 * 8(24)
    \\ stfd 3, 4 * 8(24)
    \\ stfd 4, 5 * 8(24)
    \\ stfd 5, 6 * 8(24)
    \\ stfd 6, 7 * 8(24)
    \\ stfdu 7, 8 * 8(24)
    \\ bdnz 20002b
    \\ .align 4
    \\ 20003:
    \\ andi. 0, 22, 7
    \\ mtspr 9, 0
    \\ ble 20005f
    \\ .align 4
    \\ 20004:
    \\ .if \bp0
    \\ lfd 0, 0 * 8(8)
    \\ add 8, 8, 9
    \\ .endif
    \\ .if \bp1
    \\ lfd 0, 0 * 8(9)
    \\ add 9, 9, 10
    \\ .endif
    \\ stfdu 0, 1 * 8(24)
    \\ bdnz 20004b
    \\ .align 4
    \\ 20005:
    \\ .if \bp0
    \\ mr 24, 10
    \\ .endif
    \\ .if \bp1
    \\ mr 24, 5
    \\ .endif
    \\ addi 12, 12, -8
    \\ .if \bp0
    \\ srawi. 23, 4, 3
    \\ .endif
    \\ .if \bp1
    \\ sradi. 23, 4, 3
    \\ .endif
    \\ ble 20015f
    \\ .align 4
    \\ 20006:
    \\ .if \bp0
    \\ mr 14, 6
    \\ add 15, 6, 7
    \\ add 16, 15, 7
    \\ add 17, 16, 7
    \\ add 18, 17, 7
    \\ add 19, 18, 7
    \\ add 20, 19, 7
    \\ add 21, 20, 7
    \\ add 6, 21, 7
    \\ .endif
    \\ .if \bp1
    \\ mr 14, 7
    \\ add 15, 7, 8
    \\ add 16, 15, 8
    \\ add 17, 16, 8
    \\ add 18, 17, 8
    \\ add 19, 18, 8
    \\ add 20, 19, 8
    \\ add 21, 20, 8
    \\ add 7, 21, 8
    \\ .endif
    \\ mr 27, 12
    \\ lfd 0, 144(1)
    \\ fmr 1, 0
    \\ fmr 2, 0
    \\ fmr 3, 0
    \\ fmr 4, 0
    \\ fmr 5, 0
    \\ fmr 6, 0
    \\ fmr 7, 0
    \\ .if \bp9
    \\ dcbt 26, 24
    \\ .endif
    \\ .if \bp10
    \\ dcbt 8, 26, 24
    \\ .endif
    \\ .if \bp7
    \\ dcbtst 26, 24
    \\ .endif
    \\ .if \bp0
    \\ srawi. 0, 22, 4
    \\ .endif
    \\ .if \bp1
    \\ sradi. 0, 22, 4
    \\ .endif
    \\ mtspr 9, 0
    \\ ble 20009f
    \\ lfd 16, 1 * 8(14)
    \\ lfd 17, 1 * 8(15)
    \\ lfd 18, 1 * 8(16)
    \\ lfd 19, 1 * 8(17)
    \\ lfd 20, 1 * 8(18)
    \\ lfd 21, 1 * 8(19)
    \\ lfd 22, 1 * 8(20)
    \\ lfd 23, 1 * 8(21)
    \\ lfd 24, 1 * 8(27)
    \\ lfd 25, 2 * 8(27)
    \\ lfd 26, 3 * 8(27)
    \\ lfd 27, 4 * 8(27)
    \\ lfd 28, 5 * 8(27)
    \\ lfd 29, 6 * 8(27)
    \\ lfd 30, 7 * 8(27)
    \\ lfd 31, 8 * 8(27)
    \\ bdz 20008f
    \\
;
