//! Section 1 of 2 of the text of `zgemv_n`, which `__root__.zig` joins, from OpenBLAS `kernel/sparc/zgemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1
    \\ .align 32
    \\ .proc 07
    \\ .if \bp0
    \\ save %sp, -64, %sp
    \\ st %i3, [%sp + 128 + 16]
    \\ st %i4, [%sp + 128 + 20]
    \\ st %i5, [%sp + 128 + 24]
    \\ ld [%sp + 128 + 32], %i5
    \\ ld [%sp + 128 + 36], %i2
    \\ ld [%sp + 128 + 40], %i3
    \\ ld [%sp + 128 + 44], %i4
    \\ ld [%sp + 128 + 48], %l0
    \\ ld [%sp + 128 + 52], %l1
    \\ ld [%sp + 128 + 56], %l2
    \\ sll %i2, 4, %i2
    \\ .endif
    \\ .if \bp1
    \\ save %sp, -256, %sp
    \\ ldx [%sp + 2423 + 56], %i2
    \\ ldx [%sp + 2423 + 64], %i3
    \\ ldx [%sp + 2423 + 72], %i4
    \\ ldx [%sp + 2423 + 80], %l0
    \\ ldx [%sp + 2423 + 88], %l1
    \\ ldx [%sp + 2423 + 96], %l2
    \\ std %f6, [%sp + 2423 + 32]
    \\ std %f8, [%sp + 2423 + 40]
    \\ sllx %i2, 4, %i2
    \\ .endif
    \\ cmp %i0, 0
    \\ ble %icc, 20019f
    \\ .if \bp0
    \\ sll %i4, 4, %i4
    \\ .endif
    \\ .if \bp1
    \\ sllx %i4, 4, %i4
    \\ .endif
    \\ cmp %i1, 0
    \\ ble %icc, 20019f
    \\ .if \bp0
    \\ sll %l1, 4, %l1
    \\ .endif
    \\ .if \bp1
    \\ sllx %l1, 4, %l1
    \\ .endif
    \\ cmp %l1, 2 * 8
    \\ be %icc, 20001f
    \\ mov %l0, %l6
    \\ .word ((2 << 30) | ((19) << 25) | ( 0x36 << 19) | ( 0x60 << 5))
    \\ add %i0, 3, %l5
    \\ sra %l5, 2, %l5
    \\ mov %l2, %l6
    \\ mov %l2, %l4
    \\ 20000:
    \\ std %f50, [%l4 + 0 * 8]
    \\ nop
    \\ std %f50, [%l4 + 1 * 8]
    \\ std %f50, [%l4 + 2 * 8]
    \\ std %f50, [%l4 + 3 * 8]
    \\ std %f50, [%l4 + 4 * 8]
    \\ nop
    \\ std %f50, [%l4 + 5 * 8]
    \\ deccc %l5
    \\ std %f50, [%l4 + 6 * 8]
    \\ nop
    \\ std %f50, [%l4 + 7 * 8]
    \\ bg,pn %icc, 20000b
    \\ add %l4, 8 * 8, %l4
    \\ 20001:
    \\ sra %i1, 1, %l5
    \\ cmp %l5, 0
    \\ ble,pn %icc, 20009f
    \\ nop
    \\ 20002:
    \\ mov %l6, %l4
    \\ mov %i5, %o0
    \\ .if \bp0
    \\ ldd [%sp + 128 + 16], %f52
    \\ ldd [%sp + 128 + 24], %f54
    \\ .endif
    \\ .if \bp1
    \\ ldd [%sp + 2423 + 32], %f52
    \\ ldd [%sp + 2423 + 40], %f54
    \\ .endif
    \\ add %i5, %i2, %o1
    \\ add %o1, %i2, %i5
    \\ ldd [%i3 + 0 * 8], %f56
    \\ ldd [%i3 + 1 * 8], %f58
    \\ add %i3, %i4, %i3
    \\ ldd [%i3 + 0 * 8], %f60
    \\ ldd [%i3 + 1 * 8], %f62
    \\ add %i3, %i4, %i3
    \\ fmuld %f52, %f56, %f24
    \\ fmuld %f54, %f58, %f30
    \\ fmuld %f54, %f56, %f26
    \\ fmuld %f52, %f58, %f28
    \\ fmuld %f52, %f60, %f32
    \\ fmuld %f54, %f62, %f38
    \\ fmuld %f54, %f60, %f34
    \\ fmuld %f52, %f62, %f36
    \\ fsubd %f24, %f30, %f56
    \\ faddd %f26, %f28, %f58
    \\ fsubd %f32, %f38, %f60
    \\ faddd %f34, %f36, %f62
    \\ sra %i0, 2, %l3
    \\ cmp %l3, 0
    \\ ble,pn %icc, 20006f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o0 + 1 * 8], %f26
    \\ ldd [%o0 + 2 * 8], %f28
    \\ ldd [%o0 + 3 * 8], %f30
    \\ ldd [%o0 + 4 * 8], %f40
    \\ ldd [%o0 + 5 * 8], %f42
    \\ ldd [%o0 + 6 * 8], %f44
    \\ ldd [%o0 + 7 * 8], %f46
    \\ ldd [%o1 + 0 * 8], %f32
    \\ ldd [%o1 + 1 * 8], %f34
    \\ ldd [%o1 + 2 * 8], %f36
    \\ ldd [%o1 + 3 * 8], %f38
    \\ ldd [%o1 + 4 * 8], %f48
    \\ ldd [%o1 + 5 * 8], %f50
    \\ ldd [%o1 + 6 * 8], %f52
    \\ ldd [%o1 + 7 * 8], %f54
    \\ ldd [%l4 + 0 * 8], %f8
    \\ ldd [%l4 + 1 * 8], %f10
    \\ ldd [%l4 + 2 * 8], %f12
    \\ fmuld %f24, %f56, %f0
    \\ deccc %l3
    \\ fmuld %f24, %f58, %f2
    \\ ldd [%o0 + 8 * 8], %f24
    \\ fmuld %f28, %f56, %f4
    \\ fmuld %f28, %f58, %f6
    \\ ble,pn %icc, 20005f
    \\ ldd [%o0 + 10 * 8], %f28
    \\ faddd %f8, %f0, %f8
    \\ ldd [%l4 + 3 * 8], %f14
    \\ fmuld %f26, %f58, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f26, %f56, %f2
    \\ ldd [%o0 + 9 * 8], %f26
    \\ faddd %f12, %f4, %f12
    \\ ldd [%l4 + 4 * 8], %f16
    \\ fmuld %f30, %f58, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f30, %f56, %f6
    \\ ldd [%o0 + 11 * 8], %f30
    \\ fsubd %f8, %f0, %f8
    \\ ldd [%l4 + 5 * 8], %f18
    \\ fmuld %f32, %f60, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f32, %f62, %f2
    \\ ldd [%o1 + 8 * 8], %f32
    \\ fsubd %f12, %f4, %f12
    \\ ldd [%l4 + 6 * 8], %f20
    \\ fmuld %f36, %f60, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f36, %f62, %f6
    \\ ldd [%o1 + 10 * 8], %f36
    \\ faddd %f8, %f0, %f8
    \\ ldd [%l4 + 7 * 8], %f22
    \\ fmuld %f34, %f62, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f34, %f60, %f2
    \\ ldd [%o1 + 9 * 8], %f34
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f38, %f62, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f38, %f60, %f6
    \\ ldd [%o1 + 11 * 8], %f38
    \\ fsubd %f8, %f0, %f8
    \\ fmuld %f40, %f56, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f40, %f58, %f2
    \\ ldd [%o0 + 12 * 8], %f40
    \\ fsubd %f12, %f4, %f12
    \\ deccc %l3
    \\ fmuld %f44, %f56, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f44, %f58, %f6
    \\ ble,pn %icc, 20004f
    \\ ldd [%o0 + 14 * 8], %f44
    \\ 20003:
    \\ faddd %f16, %f0, %f16
    \\ prefetch [%o0 + 44 * 8], 1
    \\ fmuld %f42, %f58, %f0
    \\ ldd [%l4 + 7 * 8], %f22
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f42, %f56, %f2
    \\ ldd [%o0 + 13 * 8], %f42
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f46, %f58, %f4
    \\ std %f8, [%l4 + 0 * 8]
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f46, %f56, %f6
    \\ ldd [%o0 + 15 * 8], %f46
    \\ fsubd %f16, %f0, %f16
    \\ fmuld %f48, %f60, %f0
    \\ std %f10, [%l4 + 1 * 8]
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f48, %f62, %f2
    \\ ldd [%o1 + 12 * 8], %f48
    \\ fsubd %f20, %f4, %f20
    \\ fmuld %f52, %f60, %f4
    \\ std %f12, [%l4 + 2 * 8]
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f52, %f62, %f6
    \\ ldd [%o1 + 14 * 8], %f52
    \\ faddd %f16, %f0, %f16
    \\ fmuld %f50, %f62, %f0
    \\ std %f14, [%l4 + 3 * 8]
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f50, %f60, %f2
    \\ ldd [%o1 + 13 * 8], %f50
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f54, %f62, %f4
    \\ ldd [%l4 + 8 * 8], %f8
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f54, %f60, %f6
    \\ ldd [%o1 + 15 * 8], %f54
    \\ fsubd %f16, %f0, %f16
    \\ fmuld %f24, %f56, %f0
    \\ ldd [%l4 + 9 * 8], %f10
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f24, %f58, %f2
    \\ ldd [%o0 + 16 * 8], %f24
    \\ fsubd %f20, %f4, %f20
    \\ fmuld %f28, %f56, %f4
    \\ ldd [%l4 + 10 * 8], %f12
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f28, %f58, %f6
    \\ ldd [%o0 + 18 * 8], %f28
    \\ faddd %f8, %f0, %f8
    \\ prefetch [%o1 + 44 * 8], 1
    \\ fmuld %f26, %f58, %f0
    \\ ldd [%l4 + 11 * 8], %f14
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f26, %f56, %f2
    \\ ldd [%o0 + 17 * 8], %f26
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f30, %f58, %f4
    \\ std %f16, [%l4 + 4 * 8]
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f30, %f56, %f6
    \\ ldd [%o0 + 19 * 8], %f30
    \\ fsubd %f8, %f0, %f8
    \\ fmuld %f32, %f60, %f0
    \\ std %f18, [%l4 + 5 * 8]
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f32, %f62, %f2
    \\ ldd [%o1 + 16 * 8], %f32
    \\ fsubd %f12, %f4, %f12
    \\ fmuld %f36, %f60, %f4
    \\ std %f20, [%l4 + 6 * 8]
    \\ faddd %f14, %f6, %f14
    \\ deccc %l3
    \\ fmuld %f36, %f62, %f6
    \\ ldd [%o1 + 18 * 8], %f36
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f34, %f62, %f0
    \\ std %f22, [%l4 + 7 * 8]
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f34, %f60, %f2
    \\ ldd [%o1 + 17 * 8], %f34
    \\ faddd %f12, %f4, %f12
    \\ add %o0, 8 * 8, %o0
    \\ fmuld %f38, %f62, %f4
    \\ ldd [%l4 + 12 * 8], %f16
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f38, %f60, %f6
    \\ ldd [%o1 + 19 * 8], %f38
    \\ fsubd %f8, %f0, %f8
    \\ add %o1, 8 * 8, %o1
    \\ fmuld %f40, %f56, %f0
    \\ ldd [%l4 + 13 * 8], %f18
    \\ faddd %f10, %f2, %f10
    \\ add %l4, 8 * 8, %l4
    \\ fmuld %f40, %f58, %f2
    \\ ldd [%o0 + 12 * 8], %f40
    \\ fsubd %f12, %f4, %f12
    \\ fmuld %f44, %f56, %f4
    \\ ldd [%l4 + 6 * 8], %f20
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f44, %f58, %f6
    \\ bg,pn %icc, 20003b
    \\ ldd [%o0 + 14 * 8], %f44
    \\ 20004:
    \\ faddd %f16, %f0, %f16
    \\ fmuld %f42, %f58, %f0
    \\ ldd [%l4 + 7 * 8], %f22
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f42, %f56, %f2
    \\ ldd [%o0 + 13 * 8], %f42
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f46, %f58, %f4
    \\ std %f8, [%l4 + 0 * 8]
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f46, %f56, %f6
    \\ ldd [%o0 + 15 * 8], %f46
    \\ fsubd %f16, %f0, %f16
    \\ fmuld %f48, %f60, %f0
    \\ std %f10, [%l4 + 1 * 8]
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f48, %f62, %f2
    \\ ldd [%o1 + 12 * 8], %f48
    \\ fsubd %f20, %f4, %f20
    \\ fmuld %f52, %f60, %f4
    \\ std %f12, [%l4 + 2 * 8]
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f52, %f62, %f6
    \\ ldd [%o1 + 14 * 8], %f52
    \\ faddd %f16, %f0, %f16
    \\ fmuld %f50, %f62, %f0
    \\ std %f14, [%l4 + 3 * 8]
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f50, %f60, %f2
    \\ ldd [%o1 + 13 * 8], %f50
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f54, %f62, %f4
    \\ ldd [%l4 + 8 * 8], %f8
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f54, %f60, %f6
    \\ ldd [%o1 + 15 * 8], %f54
    \\ fsubd %f16, %f0, %f16
    \\ add %o0, 8 * 8, %o0
    \\ fmuld %f24, %f56, %f0
    \\ ldd [%l4 + 9 * 8], %f10
    \\ faddd %f18, %f2, %f18
    \\ add %o1, 8 * 8, %o1
    \\ fmuld %f24, %f58, %f2
    \\ ldd [%o0 + 8 * 8], %f24
    \\ fsubd %f20, %f4, %f20
    \\ fmuld %f28, %f56, %f4
    \\ ldd [%l4 + 10 * 8], %f12
    \\ faddd %f22, %f6, %f22
    \\ add %l4, 8 * 8, %l4
    \\ fmuld %f28, %f58, %f6
    \\ ldd [%o0 + 10 * 8], %f28
    \\ std %f16, [%l4 - 4 * 8]
    \\ std %f18, [%l4 - 3 * 8]
    \\ std %f20, [%l4 - 2 * 8]
    \\ std %f22, [%l4 - 1 * 8]
    \\ 20005:
    \\ faddd %f8, %f0, %f8
    \\ ldd [%l4 + 3 * 8], %f14
    \\ fmuld %f26, %f58, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f26, %f56, %f2
    \\ faddd %f12, %f4, %f12
    \\ ldd [%l4 + 4 * 8], %f16
    \\ fmuld %f30, %f58, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f30, %f56, %f6
    \\ fsubd %f8, %f0, %f8
    \\ ldd [%l4 + 5 * 8], %f18
    \\ fmuld %f32, %f60, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f32, %f62, %f2
    \\ fsubd %f12, %f4, %f12
    \\ ldd [%l4 + 6 * 8], %f20
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f36, %f62, %f6
    \\ faddd %f8, %f0, %f8
    \\ ldd [%l4 + 7 * 8], %f22
    \\ fmuld %f36, %f60, %f4
    \\ fmuld %f34, %f62, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f34, %f60, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f38, %f62, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f38, %f60, %f6
    \\ fsubd %f8, %f0, %f8
    \\ fmuld %f40, %f56, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f40, %f58, %f2
    \\ fsubd %f12, %f4, %f12
    \\ fmuld %f44, %f56, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f44, %f58, %f6
    \\ faddd %f16, %f0, %f16
    \\ fmuld %f42, %f58, %f0
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f42, %f56, %f2
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f46, %f58, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f46, %f56, %f6
    \\ fsubd %f16, %f0, %f16
    \\ fmuld %f48, %f60, %f0
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f48, %f62, %f2
    \\ fsubd %f20, %f4, %f20
    \\ fmuld %f52, %f60, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f52, %f62, %f6
    \\ faddd %f16, %f0, %f16
    \\ fmuld %f50, %f62, %f0
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f50, %f60, %f2
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f54, %f62, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f54, %f60, %f6
    \\ std %f8, [%l4 + 0 * 8]
    \\ fsubd %f16, %f0, %f16
    \\ std %f10, [%l4 + 1 * 8]
    \\ faddd %f18, %f2, %f18
    \\ std %f12, [%l4 + 2 * 8]
    \\ fsubd %f20, %f4, %f20
    \\ std %f14, [%l4 + 3 * 8]
    \\ faddd %f22, %f6, %f22
    \\ std %f16, [%l4 + 4 * 8]
    \\ add %o0, 8 * 8, %o0
    \\ std %f18, [%l4 + 5 * 8]
    \\ add %o1, 8 * 8, %o1
    \\ std %f20, [%l4 + 6 * 8]
    \\ std %f22, [%l4 + 7 * 8]
    \\ add %l4, 8 * 8, %l4
    \\
;
