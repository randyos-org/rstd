//! Section 3 of 3 of the text of `dgemv_n`, which `__root__.zig` joins, from OpenBLAS `kernel/sparc/gemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ 20014:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20015f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o1 + 0 * 8], %f26
    \\ ldd [%l4 + 0 * 8], %f8
    \\ ldd [%o0 + 1 * 8], %f32
    \\ ldd [%o1 + 1 * 8], %f34
    \\ add %o0, 2 * 8, %o0
    \\ ldd [%l4 + 1 * 8], %f10
    \\ add %o1, 2 * 8, %o1
    \\ fmuld %f24, %f56, %f0
    \\ fmuld %f26, %f58, %f2
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f32, %f56, %f0
    \\ faddd %f8, %f2, %f8
    \\ fmuld %f34, %f58, %f2
    \\ faddd %f10, %f0, %f10
    \\ faddd %f10, %f2, %f10
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ add %l4, 2 * 8, %l4
    \\ 20015:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20016f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o1 + 0 * 8], %f26
    \\ ldd [%l4 + 0 * 8], %f8
    \\ fmuld %f24, %f56, %f0
    \\ fmuld %f26, %f58, %f2
    \\ faddd %f8, %f0, %f8
    \\ faddd %f8, %f2, %f8
    \\ std %f8, [%l4]
    \\ 20016:
    \\ andcc %i1, 1, %l5
    \\ ble,pn %icc, 20023f
    \\ nop
    \\ 20017:
    \\ mov %l6, %l4
    \\ .if \bp0
    \\ mov %i5, %o0
    \\ ldd [%sp + 128 + 16], %f54
    \\ ldd [%i3], %f56
    \\ add %i3, %i4, %i3
    \\ .endif
    \\ .if \bp1
    \\ mov %i4, %o0
    \\ ldd [%sp + 2423 + 40], %f54
    \\ ldd [%i2], %f56
    \\ add %i2, %i3, %i2
    \\ .endif
    \\ fmuld %f54, %f56, %f56
    \\ sra %i0, 3, %l3
    \\ cmp %l3, 0
    \\ ble,pn %icc, 20020f
    \\ nop
    \\ ldd [%l4 + 0 * 8], %f8
    \\ ldd [%l4 + 1 * 8], %f10
    \\ ldd [%l4 + 2 * 8], %f12
    \\ ldd [%l4 + 3 * 8], %f14
    \\ ldd [%l4 + 4 * 8], %f16
    \\ ldd [%l4 + 5 * 8], %f18
    \\ ldd [%l4 + 6 * 8], %f20
    \\ ldd [%l4 + 7 * 8], %f22
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o0 + 1 * 8], %f26
    \\ ldd [%o0 + 2 * 8], %f28
    \\ ldd [%o0 + 3 * 8], %f30
    \\ ldd [%o0 + 4 * 8], %f32
    \\ ldd [%o0 + 5 * 8], %f34
    \\ ldd [%o0 + 6 * 8], %f36
    \\ ldd [%o0 + 7 * 8], %f38
    \\ fmuld %f24, %f56, %f0
    \\ deccc %l3
    \\ ldd [%o0 + 8 * 8], %f24
    \\ fmuld %f26, %f56, %f2
    \\ ldd [%o0 + 9 * 8], %f26
    \\ fmuld %f28, %f56, %f4
    \\ ldd [%o0 + 10 * 8], %f28
    \\ fmuld %f30, %f56, %f6
    \\ ble,pn %icc, 20019f
    \\ ldd [%o0 + 11 * 8], %f30
    \\ 20018:
    \\ faddd %f8, %f0, %f8
    \\ prefetch [%o0 + 20 * 8], 1
    \\ fmuld %f32, %f56, %f0
    \\ ldd [%o0 + 12 * 8], %f32
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f34, %f56, %f2
    \\ ldd [%o0 + 13 * 8], %f34
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f36, %f56, %f4
    \\ ldd [%o0 + 14 * 8], %f36
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f38, %f56, %f6
    \\ ldd [%o0 + 15 * 8], %f38
    \\ faddd %f16, %f0, %f16
    \\ fmuld %f24, %f56, %f0
    \\ ldd [%o0 + 16 * 8], %f24
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f26, %f56, %f2
    \\ ldd [%o0 + 17 * 8], %f26
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f28, %f56, %f4
    \\ ldd [%o0 + 18 * 8], %f28
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f30, %f56, %f6
    \\ ldd [%o0 + 19 * 8], %f30
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ std %f12, [%l4 + 2 * 8]
    \\ std %f14, [%l4 + 3 * 8]
    \\ std %f16, [%l4 + 4 * 8]
    \\ std %f18, [%l4 + 5 * 8]
    \\ std %f20, [%l4 + 6 * 8]
    \\ std %f22, [%l4 + 7 * 8]
    \\ ldd [%l4 + 8 * 8], %f8
    \\ ldd [%l4 + 9 * 8], %f10
    \\ ldd [%l4 + 10 * 8], %f12
    \\ ldd [%l4 + 11 * 8], %f14
    \\ ldd [%l4 + 12 * 8], %f16
    \\ deccc %l3
    \\ ldd [%l4 + 13 * 8], %f18
    \\ add %o0, 8 * 8, %o0
    \\ ldd [%l4 + 14 * 8], %f20
    \\ add %l4, 8 * 8, %l4
    \\ bg,pn %icc, 20018b
    \\ ldd [%l4 + 7 * 8], %f22
    \\ 20019:
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f32, %f56, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f34, %f56, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f36, %f56, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f38, %f56, %f6
    \\ std %f8, [%l4 + 0 * 8]
    \\ faddd %f16, %f0, %f16
    \\ std %f10, [%l4 + 1 * 8]
    \\ faddd %f18, %f2, %f18
    \\ std %f12, [%l4 + 2 * 8]
    \\ faddd %f20, %f4, %f20
    \\ std %f14, [%l4 + 3 * 8]
    \\ faddd %f22, %f6, %f22
    \\ std %f16, [%l4 + 4 * 8]
    \\ std %f18, [%l4 + 5 * 8]
    \\ std %f20, [%l4 + 6 * 8]
    \\ add %o0, 8 * 8, %o0
    \\ std %f22, [%l4 + 7 * 8]
    \\ add %l4, 8 * 8, %l4
    \\ 20020:
    \\ andcc %i0, 4, %l3
    \\ ble,pn %icc, 20021f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o0 + 1 * 8], %f26
    \\ ldd [%o0 + 2 * 8], %f28
    \\ ldd [%o0 + 3 * 8], %f30
    \\ ldd [%l4 + 0 * 8], %f8
    \\ add %o0, 4 * 8, %o0
    \\ ldd [%l4 + 1 * 8], %f10
    \\ ldd [%l4 + 2 * 8], %f12
    \\ ldd [%l4 + 3 * 8], %f14
    \\ fmuld %f24, %f56, %f0
    \\ fmuld %f26, %f56, %f2
    \\ fmuld %f28, %f56, %f4
    \\ fmuld %f30, %f56, %f6
    \\ faddd %f8, %f0, %f8
    \\ faddd %f10, %f2, %f10
    \\ faddd %f12, %f4, %f12
    \\ faddd %f14, %f6, %f14
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ std %f12, [%l4 + 2 * 8]
    \\ std %f14, [%l4 + 3 * 8]
    \\ add %l4, 4 * 8, %l4
    \\ 20021:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20022f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%l4 + 0 * 8], %f8
    \\ ldd [%o0 + 1 * 8], %f32
    \\ ldd [%l4 + 1 * 8], %f10
    \\ add %o0, 2 * 8, %o0
    \\ fmuld %f24, %f56, %f0
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f32, %f56, %f0
    \\ faddd %f10, %f0, %f10
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ add %l4, 2 * 8, %l4
    \\ 20022:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20023f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%l4 + 0 * 8], %f8
    \\ fmuld %f24, %f56, %f0
    \\ faddd %f8, %f0, %f8
    \\ std %f8, [%l4]
    \\ 20023:
    \\ cmp %l1, 8
    \\ be %icc, 20028f
    \\ mov %l0, %l4
    \\ sra %i0, 3, %l3
    \\ cmp %l3, 0
    \\ ble,pn %icc, 20025f
    \\ nop
    \\ 20024:
    \\ ldd [%l2 + 0 * 8], %f24
    \\ ldd [%l0], %f8
    \\ add %l0, %l1, %l0
    \\ ldd [%l2 + 1 * 8], %f26
    \\ ldd [%l0], %f10
    \\ add %l0, %l1, %l0
    \\ ldd [%l2 + 2 * 8], %f28
    \\ ldd [%l0], %f12
    \\ add %l0, %l1, %l0
    \\ ldd [%l2 + 3 * 8], %f30
    \\ ldd [%l0], %f14
    \\ add %l0, %l1, %l0
    \\ ldd [%l2 + 4 * 8], %f32
    \\ faddd %f8, %f24, %f8
    \\ ldd [%l0], %f16
    \\ add %l0, %l1, %l0
    \\ ldd [%l2 + 5 * 8], %f34
    \\ faddd %f10, %f26, %f10
    \\ ldd [%l0], %f18
    \\ add %l0, %l1, %l0
    \\ ldd [%l2 + 6 * 8], %f36
    \\ faddd %f12, %f28, %f12
    \\ ldd [%l0], %f20
    \\ add %l0, %l1, %l0
    \\ ldd [%l2 + 7 * 8], %f38
    \\ faddd %f14, %f30, %f14
    \\ ldd [%l0], %f22
    \\ add %l0, %l1, %l0
    \\ std %f8, [%l4]
    \\ faddd %f16, %f32, %f16
    \\ add %l4, %l1, %l4
    \\ std %f10, [%l4]
    \\ faddd %f18, %f34, %f18
    \\ add %l4, %l1, %l4
    \\ std %f12, [%l4]
    \\ faddd %f20, %f36, %f20
    \\ add %l4, %l1, %l4
    \\ std %f14, [%l4]
    \\ faddd %f22, %f38, %f22
    \\ add %l4, %l1, %l4
    \\ std %f16, [%l4]
    \\ add %l4, %l1, %l4
    \\ std %f18, [%l4]
    \\ add %l4, %l1, %l4
    \\ std %f20, [%l4]
    \\ add %l4, %l1, %l4
    \\ std %f22, [%l4]
    \\ add %l4, %l1, %l4
    \\ deccc %l3
    \\ bg,pn %icc, 20024b
    \\ add %l2, 8 * 8, %l2
    \\ 20025:
    \\ andcc %i0, 7, %l3
    \\ ble,pn %icc, 20028f
    \\ nop
    \\ andcc %i0, 4, %l3
    \\ ble,pn %icc, 20026f
    \\ nop
    \\ ldd [%l2 + 0 * 8], %f24
    \\ ldd [%l2 + 1 * 8], %f26
    \\ ldd [%l2 + 2 * 8], %f28
    \\ ldd [%l2 + 3 * 8], %f30
    \\ add %l2, 4 * 8, %l2
    \\ ldd [%l0], %f8
    \\ add %l0, %l1, %l0
    \\ ldd [%l0], %f10
    \\ add %l0, %l1, %l0
    \\ ldd [%l0], %f12
    \\ add %l0, %l1, %l0
    \\ ldd [%l0], %f14
    \\ add %l0, %l1, %l0
    \\ faddd %f8, %f24, %f8
    \\ faddd %f10, %f26, %f10
    \\ faddd %f12, %f28, %f12
    \\ faddd %f14, %f30, %f14
    \\ std %f8, [%l4]
    \\ add %l4, %l1, %l4
    \\ std %f10, [%l4]
    \\ add %l4, %l1, %l4
    \\ std %f12, [%l4]
    \\ add %l4, %l1, %l4
    \\ std %f14, [%l4]
    \\ add %l4, %l1, %l4
    \\ 20026:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20027f
    \\ nop
    \\ ldd [%l2 + 0 * 8], %f24
    \\ ldd [%l2 + 1 * 8], %f26
    \\ add %l2, 2 * 8, %l2
    \\ ldd [%l0], %f8
    \\ add %l0, %l1, %l0
    \\ ldd [%l0], %f10
    \\ add %l0, %l1, %l0
    \\ faddd %f8, %f24, %f8
    \\ faddd %f10, %f26, %f10
    \\ std %f8, [%l4]
    \\ add %l4, %l1, %l4
    \\ std %f10, [%l4]
    \\ add %l4, %l1, %l4
    \\ 20027:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20028f
    \\ nop
    \\ ldd [%l2 + 0 * 8], %f24
    \\ ldd [%l0], %f8
    \\ faddd %f8, %f24, %f8
    \\ std %f8, [%l4]
    \\ 20028:
    \\ return %i7 + 8
    \\ clr %o0
    \\ .endm
    \\
;
