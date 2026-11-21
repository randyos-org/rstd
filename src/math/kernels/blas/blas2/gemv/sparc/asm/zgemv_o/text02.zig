//! Section 2 of 2 of the text of `zgemv_o`, which `__root__.zig` joins, from OpenBLAS `kernel/sparc/zgemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ 20006:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20007f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o0 + 1 * 8], %f26
    \\ ldd [%o0 + 2 * 8], %f28
    \\ ldd [%o0 + 3 * 8], %f30
    \\ ldd [%l4 + 0 * 8], %f8
    \\ ldd [%l4 + 1 * 8], %f10
    \\ ldd [%l4 + 2 * 8], %f12
    \\ ldd [%l4 + 3 * 8], %f14
    \\ fmuld %f24, %f56, %f0
    \\ ldd [%o1 + 0 * 8], %f32
    \\ fmuld %f24, %f58, %f2
    \\ ldd [%o1 + 1 * 8], %f34
    \\ fmuld %f28, %f56, %f4
    \\ ldd [%o1 + 2 * 8], %f36
    \\ fmuld %f28, %f58, %f6
    \\ ldd [%o1 + 3 * 8], %f38
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f26, %f58, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f26, %f56, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f30, %f58, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f30, %f56, %f6
    \\ fsubd %f8, %f0, %f8
    \\ fmuld %f32, %f60, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f32, %f62, %f2
    \\ fsubd %f12, %f4, %f12
    \\ fmuld %f36, %f60, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f36, %f62, %f6
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f34, %f62, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f34, %f60, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f38, %f62, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f38, %f60, %f6
    \\ fsubd %f8, %f0, %f8
    \\ faddd %f10, %f2, %f10
    \\ fsubd %f12, %f4, %f12
    \\ faddd %f14, %f6, %f14
    \\ std %f8, [%l4 + 0 * 8]
    \\ add %o0, 4 * 8, %o0
    \\ std %f10, [%l4 + 1 * 8]
    \\ add %o1, 4 * 8, %o1
    \\ std %f12, [%l4 + 2 * 8]
    \\ nop
    \\ std %f14, [%l4 + 3 * 8]
    \\ add %l4, 4 * 8, %l4
    \\ 20007:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20008f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o0 + 1 * 8], %f26
    \\ ldd [%o1 + 0 * 8], %f28
    \\ ldd [%o1 + 1 * 8], %f30
    \\ ldd [%l4 + 0 * 8], %f8
    \\ ldd [%l4 + 1 * 8], %f10
    \\ fmuld %f24, %f56, %f0
    \\ fmuld %f24, %f58, %f2
    \\ fmuld %f26, %f58, %f4
    \\ fmuld %f26, %f56, %f6
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f28, %f60, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f28, %f62, %f2
    \\ fsubd %f8, %f4, %f8
    \\ fmuld %f30, %f62, %f4
    \\ faddd %f10, %f6, %f10
    \\ fmuld %f30, %f60, %f6
    \\ faddd %f8, %f0, %f8
    \\ faddd %f10, %f2, %f10
    \\ fsubd %f8, %f4, %f8
    \\ faddd %f10, %f6, %f10
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ 20008:
    \\ deccc %l5
    \\ bg %icc, 20002b
    \\ nop
    \\ 20009:
    \\ andcc %i1, 1, %l5
    \\ ble,pn %icc, 20015f
    \\ nop
    \\ 20010:
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
    \\ ldd [%i3 + 0 * 8], %f56
    \\ ldd [%i3 + 1 * 8], %f58
    \\ fmuld %f52, %f56, %f24
    \\ fmuld %f54, %f56, %f26
    \\ fmuld %f52, %f58, %f28
    \\ fmuld %f54, %f58, %f30
    \\ faddd %f24, %f30, %f56
    \\ fsubd %f26, %f28, %f58
    \\ sra %i0, 2, %l3
    \\ cmp %l3, 0
    \\ ble,pn %icc, 20013f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o0 + 1 * 8], %f26
    \\ ldd [%o0 + 2 * 8], %f28
    \\ ldd [%o0 + 3 * 8], %f30
    \\ ldd [%o0 + 4 * 8], %f40
    \\ ldd [%o0 + 5 * 8], %f42
    \\ ldd [%o0 + 6 * 8], %f44
    \\ ldd [%o0 + 7 * 8], %f46
    \\ ldd [%l4 + 0 * 8], %f8
    \\ ldd [%l4 + 1 * 8], %f10
    \\ ldd [%l4 + 2 * 8], %f12
    \\ ldd [%l4 + 3 * 8], %f14
    \\ ldd [%l4 + 4 * 8], %f16
    \\ ldd [%l4 + 5 * 8], %f18
    \\ ldd [%l4 + 6 * 8], %f20
    \\ ldd [%l4 + 7 * 8], %f22
    \\ fmuld %f24, %f56, %f0
    \\ deccc %l3
    \\ fmuld %f24, %f58, %f2
    \\ ldd [%o0 + 8 * 8], %f24
    \\ fmuld %f28, %f56, %f4
    \\ fmuld %f28, %f58, %f6
    \\ ble,pn %icc, 20012f
    \\ ldd [%o0 + 10 * 8], %f28
    \\ 20011:
    \\ faddd %f8, %f0, %f8
    \\ prefetch [%o0 + 44 * 8], 1
    \\ fmuld %f26, %f58, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f26, %f56, %f2
    \\ ldd [%o0 + 9 * 8], %f26
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f30, %f58, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f30, %f56, %f6
    \\ ldd [%o0 + 11 * 8], %f30
    \\ fsubd %f8, %f0, %f8
    \\ fmuld %f40, %f56, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f40, %f58, %f2
    \\ ldd [%o0 + 12 * 8], %f40
    \\ fsubd %f12, %f4, %f12
    \\ fmuld %f44, %f56, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f44, %f58, %f6
    \\ ldd [%o0 + 14 * 8], %f44
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ std %f12, [%l4 + 2 * 8]
    \\ std %f14, [%l4 + 3 * 8]
    \\ faddd %f16, %f0, %f16
    \\ fmuld %f42, %f58, %f0
    \\ ldd [%l4 + 8 * 8], %f8
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f42, %f56, %f2
    \\ ldd [%o0 + 13 * 8], %f42
    \\ faddd %f20, %f4, %f20
    \\ deccc %l3
    \\ fmuld %f46, %f58, %f4
    \\ ldd [%l4 + 9 * 8], %f10
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f46, %f56, %f6
    \\ ldd [%o0 + 15 * 8], %f46
    \\ fsubd %f16, %f0, %f16
    \\ add %o0, 8 * 8, %o0
    \\ fmuld %f24, %f56, %f0
    \\ ldd [%l4 + 10 * 8], %f12
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f24, %f58, %f2
    \\ ldd [%o0 + 8 * 8], %f24
    \\ fsubd %f20, %f4, %f20
    \\ fmuld %f28, %f56, %f4
    \\ ldd [%l4 + 11 * 8], %f14
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f28, %f58, %f6
    \\ ldd [%o0 + 10 * 8], %f28
    \\ std %f16, [%l4 + 4 * 8]
    \\ std %f18, [%l4 + 5 * 8]
    \\ std %f20, [%l4 + 6 * 8]
    \\ std %f22, [%l4 + 7 * 8]
    \\ ldd [%l4 + 12 * 8], %f16
    \\ ldd [%l4 + 13 * 8], %f18
    \\ ldd [%l4 + 14 * 8], %f20
    \\ add %l4, 8 * 8, %l4
    \\ bg,pn %icc, 20011b
    \\ ldd [%l4 + 7 * 8], %f22
    \\ 20012:
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f26, %f58, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f26, %f56, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f30, %f58, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f30, %f56, %f6
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
    \\ faddd %f18, %f2, %f18
    \\ fsubd %f20, %f4, %f20
    \\ faddd %f22, %f6, %f22
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ std %f12, [%l4 + 2 * 8]
    \\ std %f14, [%l4 + 3 * 8]
    \\ std %f16, [%l4 + 4 * 8]
    \\ std %f18, [%l4 + 5 * 8]
    \\ std %f20, [%l4 + 6 * 8]
    \\ std %f22, [%l4 + 7 * 8]
    \\ add %o0, 8 * 8, %o0
    \\ add %l4, 8 * 8, %l4
    \\ 20013:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20014f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o0 + 1 * 8], %f26
    \\ ldd [%o0 + 2 * 8], %f28
    \\ ldd [%o0 + 3 * 8], %f30
    \\ ldd [%l4 + 0 * 8], %f8
    \\ fmuld %f24, %f56, %f0
    \\ ldd [%l4 + 1 * 8], %f10
    \\ fmuld %f24, %f58, %f2
    \\ ldd [%l4 + 2 * 8], %f12
    \\ fmuld %f28, %f56, %f4
    \\ ldd [%l4 + 3 * 8], %f14
    \\ fmuld %f28, %f58, %f6
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f26, %f58, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f26, %f56, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f30, %f58, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f30, %f56, %f6
    \\ fsubd %f8, %f0, %f8
    \\ faddd %f10, %f2, %f10
    \\ fsubd %f12, %f4, %f12
    \\ faddd %f14, %f6, %f14
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ std %f12, [%l4 + 2 * 8]
    \\ std %f14, [%l4 + 3 * 8]
    \\ add %o0, 4 * 8, %o0
    \\ add %l4, 4 * 8, %l4
    \\ 20014:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20015f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o0 + 1 * 8], %f26
    \\ ldd [%l4 + 0 * 8], %f8
    \\ ldd [%l4 + 1 * 8], %f10
    \\ fmuld %f24, %f56, %f0
    \\ fmuld %f24, %f58, %f2
    \\ fmuld %f26, %f58, %f4
    \\ fmuld %f26, %f56, %f6
    \\ faddd %f8, %f0, %f8
    \\ faddd %f10, %f2, %f10
    \\ fsubd %f8, %f4, %f8
    \\ faddd %f10, %f6, %f10
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ 20015:
    \\ cmp %l1, 2 * 8
    \\ be %icc, 20019f
    \\ mov %l0, %l4
    \\ sra %i0, 2, %l3
    \\ cmp %l3, 0
    \\ ble,pn %icc, 20017f
    \\ nop
    \\ 20016:
    \\ ldd [%l2 + 0 * 8], %f24
    \\ ldd [%l2 + 1 * 8], %f26
    \\ ldd [%l0 + 0 * 8], %f8
    \\ ldd [%l0 + 1 * 8], %f10
    \\ add %l0, %l1, %l0
    \\ ldd [%l2 + 2 * 8], %f28
    \\ ldd [%l2 + 3 * 8], %f30
    \\ ldd [%l0 + 0 * 8], %f12
    \\ ldd [%l0 + 1 * 8], %f14
    \\ add %l0, %l1, %l0
    \\ ldd [%l2 + 4 * 8], %f32
    \\ ldd [%l2 + 5 * 8], %f34
    \\ ldd [%l0 + 0 * 8], %f16
    \\ ldd [%l0 + 1 * 8], %f18
    \\ add %l0, %l1, %l0
    \\ ldd [%l2 + 6 * 8], %f36
    \\ ldd [%l2 + 7 * 8], %f38
    \\ ldd [%l0 + 0 * 8], %f20
    \\ ldd [%l0 + 1 * 8], %f22
    \\ add %l0, %l1, %l0
    \\ faddd %f8, %f24, %f8
    \\ faddd %f10, %f26, %f10
    \\ faddd %f12, %f28, %f12
    \\ faddd %f14, %f30, %f14
    \\ faddd %f16, %f32, %f16
    \\ faddd %f18, %f34, %f18
    \\ faddd %f20, %f36, %f20
    \\ faddd %f22, %f38, %f22
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ add %l4, %l1, %l4
    \\ std %f12, [%l4 + 0 * 8]
    \\ std %f14, [%l4 + 1 * 8]
    \\ add %l4, %l1, %l4
    \\ std %f16, [%l4 + 0 * 8]
    \\ std %f18, [%l4 + 1 * 8]
    \\ add %l4, %l1, %l4
    \\ std %f20, [%l4 + 0 * 8]
    \\ std %f22, [%l4 + 1 * 8]
    \\ add %l4, %l1, %l4
    \\ deccc %l3
    \\ bg,pn %icc, 20016b
    \\ add %l2, 8 * 8, %l2
    \\ 20017:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20018f
    \\ nop
    \\ ldd [%l2 + 0 * 8], %f24
    \\ ldd [%l2 + 1 * 8], %f26
    \\ ldd [%l0 + 0 * 8], %f8
    \\ ldd [%l0 + 1 * 8], %f10
    \\ add %l0, %l1, %l0
    \\ ldd [%l2 + 2 * 8], %f28
    \\ ldd [%l2 + 3 * 8], %f30
    \\ ldd [%l0 + 0 * 8], %f12
    \\ ldd [%l0 + 1 * 8], %f14
    \\ add %l0, %l1, %l0
    \\ faddd %f8, %f24, %f8
    \\ faddd %f10, %f26, %f10
    \\ faddd %f12, %f28, %f12
    \\ faddd %f14, %f30, %f14
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ add %l4, %l1, %l4
    \\ std %f12, [%l4 + 0 * 8]
    \\ std %f14, [%l4 + 1 * 8]
    \\ add %l4, %l1, %l4
    \\ add %l2, 4 * 8, %l2
    \\ 20018:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20019f
    \\ nop
    \\ ldd [%l2 + 0 * 8], %f24
    \\ ldd [%l2 + 1 * 8], %f26
    \\ ldd [%l0 + 0 * 8], %f8
    \\ ldd [%l0 + 1 * 8], %f10
    \\ faddd %f8, %f24, %f8
    \\ faddd %f10, %f26, %f10
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ 20019:
    \\ return %i7 + 8
    \\ clr %o0
    \\ .endm
    \\
;
