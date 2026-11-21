//! Section 2 of 3 of the text of `dgemv_n`, which `__root__.zig` joins, from OpenBLAS `kernel/sparc/gemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ 20005:
    \\ andcc %i0, 4, %l3
    \\ ble,pn %icc, 20006f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o0 + 1 * 8], %f26
    \\ ldd [%o0 + 2 * 8], %f28
    \\ ldd [%o0 + 3 * 8], %f30
    \\ ldd [%o1 + 0 * 8], %f32
    \\ ldd [%o1 + 1 * 8], %f34
    \\ ldd [%o1 + 2 * 8], %f36
    \\ ldd [%o1 + 3 * 8], %f38
    \\ ldd [%o2 + 0 * 8], %f40
    \\ ldd [%o2 + 1 * 8], %f42
    \\ ldd [%o2 + 2 * 8], %f44
    \\ ldd [%o2 + 3 * 8], %f46
    \\ ldd [%o3 + 0 * 8], %f48
    \\ ldd [%o3 + 1 * 8], %f50
    \\ ldd [%o3 + 2 * 8], %f52
    \\ ldd [%o3 + 3 * 8], %f54
    \\ ldd [%l4 + 0 * 8], %f8
    \\ add %o0, 4 * 8, %o0
    \\ ldd [%l4 + 1 * 8], %f10
    \\ add %o1, 4 * 8, %o1
    \\ ldd [%l4 + 2 * 8], %f12
    \\ add %o2, 4 * 8, %o2
    \\ ldd [%l4 + 3 * 8], %f14
    \\ add %o3, 4 * 8, %o3
    \\ fmuld %f24, %f56, %f0
    \\ fmuld %f26, %f56, %f2
    \\ fmuld %f28, %f56, %f4
    \\ fmuld %f30, %f56, %f6
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f32, %f58, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f34, %f58, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f36, %f58, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f38, %f58, %f6
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f40, %f60, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f42, %f60, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f44, %f60, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f46, %f60, %f6
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f48, %f62, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f50, %f62, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f52, %f62, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f54, %f62, %f6
    \\ faddd %f8, %f0, %f8
    \\ faddd %f10, %f2, %f10
    \\ faddd %f12, %f4, %f12
    \\ faddd %f14, %f6, %f14
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ std %f12, [%l4 + 2 * 8]
    \\ std %f14, [%l4 + 3 * 8]
    \\ add %l4, 4 * 8, %l4
    \\ 20006:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20007f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o1 + 0 * 8], %f26
    \\ ldd [%o2 + 0 * 8], %f28
    \\ ldd [%o3 + 0 * 8], %f30
    \\ ldd [%l4 + 0 * 8], %f8
    \\ ldd [%o0 + 1 * 8], %f32
    \\ ldd [%o1 + 1 * 8], %f34
    \\ ldd [%o2 + 1 * 8], %f36
    \\ ldd [%o3 + 1 * 8], %f38
    \\ ldd [%l4 + 1 * 8], %f10
    \\ add %o0, 2 * 8, %o0
    \\ add %o1, 2 * 8, %o1
    \\ add %o2, 2 * 8, %o2
    \\ add %o3, 2 * 8, %o3
    \\ fmuld %f24, %f56, %f0
    \\ fmuld %f26, %f58, %f2
    \\ fmuld %f28, %f60, %f4
    \\ fmuld %f30, %f62, %f6
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f32, %f56, %f0
    \\ faddd %f8, %f2, %f8
    \\ fmuld %f34, %f58, %f2
    \\ faddd %f8, %f4, %f8
    \\ fmuld %f36, %f60, %f4
    \\ faddd %f8, %f6, %f8
    \\ fmuld %f38, %f62, %f6
    \\ faddd %f10, %f0, %f10
    \\ faddd %f10, %f2, %f10
    \\ faddd %f10, %f4, %f10
    \\ faddd %f10, %f6, %f10
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ add %l4, 2 * 8, %l4
    \\ 20007:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20008f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o1 + 0 * 8], %f26
    \\ ldd [%o2 + 0 * 8], %f28
    \\ ldd [%o3 + 0 * 8], %f30
    \\ ldd [%l4 + 0 * 8], %f8
    \\ fmuld %f24, %f56, %f0
    \\ fmuld %f26, %f58, %f2
    \\ fmuld %f28, %f60, %f4
    \\ fmuld %f30, %f62, %f6
    \\ faddd %f8, %f0, %f8
    \\ faddd %f8, %f2, %f8
    \\ faddd %f8, %f4, %f8
    \\ faddd %f8, %f6, %f8
    \\ std %f8, [%l4]
    \\ 20008:
    \\ deccc %l5
    \\ bg %icc, 20002b
    \\ nop
    \\ 20009:
    \\ andcc %i1, 2, %l5
    \\ ble,pn %icc, 20016f
    \\ nop
    \\ 20010:
    \\ mov %l6, %l4
    \\ .if \bp0
    \\ mov %i5, %o0
    \\ add %i5, %i2, %o1
    \\ add %o1, %i2, %i5
    \\ ldd [%sp + 128 + 16], %f54
    \\ ldd [%i3], %f56
    \\ add %i3, %i4, %i3
    \\ ldd [%i3], %f58
    \\ add %i3, %i4, %i3
    \\ .endif
    \\ .if \bp1
    \\ mov %i4, %o0
    \\ add %i4, %i5, %o1
    \\ add %o1, %i5, %i4
    \\ ldd [%sp + 2423 + 40], %f54
    \\ ldd [%i2], %f56
    \\ add %i2, %i3, %i2
    \\ ldd [%i2], %f58
    \\ add %i2, %i3, %i2
    \\ .endif
    \\ fmuld %f54, %f56, %f56
    \\ fmuld %f54, %f58, %f58
    \\ sra %i0, 3, %l3
    \\ cmp %l3, 0
    \\ ble,pn %icc, 20013f
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
    \\ ldd [%o1 + 0 * 8], %f40
    \\ ldd [%o1 + 1 * 8], %f42
    \\ ldd [%o1 + 2 * 8], %f44
    \\ ldd [%o1 + 3 * 8], %f46
    \\ ldd [%o1 + 4 * 8], %f48
    \\ ldd [%o1 + 5 * 8], %f50
    \\ ldd [%o1 + 6 * 8], %f52
    \\ ldd [%o1 + 7 * 8], %f54
    \\ fmuld %f24, %f56, %f0
    \\ deccc %l3
    \\ ldd [%o0 + 8 * 8], %f24
    \\ fmuld %f26, %f56, %f2
    \\ ldd [%o0 + 9 * 8], %f26
    \\ fmuld %f28, %f56, %f4
    \\ ldd [%o0 + 10 * 8], %f28
    \\ fmuld %f30, %f56, %f6
    \\ ble,pn %icc, 20012f
    \\ ldd [%o0 + 11 * 8], %f30
    \\ 20011:
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
    \\ fmuld %f40, %f58, %f0
    \\ ldd [%o1 + 8 * 8], %f40
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f42, %f58, %f2
    \\ ldd [%o1 + 9 * 8], %f42
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f44, %f58, %f4
    \\ ldd [%o1 + 10 * 8], %f44
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f46, %f58, %f6
    \\ ldd [%o1 + 11 * 8], %f46
    \\ faddd %f8, %f0, %f8
    \\ prefetch [%o1 + 20 * 8], 1
    \\ fmuld %f48, %f58, %f0
    \\ ldd [%o1 + 12 * 8], %f48
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f50, %f58, %f2
    \\ ldd [%o1 + 13 * 8], %f50
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f52, %f58, %f4
    \\ ldd [%o1 + 14 * 8], %f52
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f54, %f58, %f6
    \\ ldd [%o1 + 15 * 8], %f54
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
    \\ add %o0, 8 * 8, %o0
    \\ ldd [%l4 + 9 * 8], %f10
    \\ add %o1, 8 * 8, %o1
    \\ ldd [%l4 + 10 * 8], %f12
    \\ deccc %l3
    \\ ldd [%l4 + 11 * 8], %f14
    \\ ldd [%l4 + 12 * 8], %f16
    \\ ldd [%l4 + 13 * 8], %f18
    \\ ldd [%l4 + 14 * 8], %f20
    \\ ldd [%l4 + 15 * 8], %f22
    \\ bg,pn %icc, 20011b
    \\ add %l4, 8 * 8, %l4
    \\ 20012:
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f32, %f56, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f34, %f56, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f36, %f56, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f38, %f56, %f6
    \\ faddd %f16, %f0, %f16
    \\ fmuld %f40, %f58, %f0
    \\ faddd %f18, %f2, %f18
    \\ fmuld %f42, %f58, %f2
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f44, %f58, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f46, %f58, %f6
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f48, %f58, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f50, %f58, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f52, %f58, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f54, %f58, %f6
    \\ std %f8, [%l4 + 0 * 8]
    \\ faddd %f16, %f0, %f16
    \\ std %f10, [%l4 + 1 * 8]
    \\ faddd %f18, %f2, %f18
    \\ std %f12, [%l4 + 2 * 8]
    \\ faddd %f20, %f4, %f20
    \\ std %f14, [%l4 + 3 * 8]
    \\ faddd %f22, %f6, %f22
    \\ std %f16, [%l4 + 4 * 8]
    \\ add %o0, 8 * 8, %o0
    \\ std %f18, [%l4 + 5 * 8]
    \\ add %o1, 8 * 8, %o1
    \\ std %f20, [%l4 + 6 * 8]
    \\ nop
    \\ std %f22, [%l4 + 7 * 8]
    \\ add %l4, 8 * 8, %l4
    \\ 20013:
    \\ andcc %i0, 4, %l3
    \\ ble,pn %icc, 20014f
    \\ nop
    \\ ldd [%o0 + 0 * 8], %f24
    \\ ldd [%o0 + 1 * 8], %f26
    \\ ldd [%o0 + 2 * 8], %f28
    \\ ldd [%o0 + 3 * 8], %f30
    \\ ldd [%o1 + 0 * 8], %f32
    \\ ldd [%o1 + 1 * 8], %f34
    \\ ldd [%o1 + 2 * 8], %f36
    \\ ldd [%o1 + 3 * 8], %f38
    \\ ldd [%l4 + 0 * 8], %f8
    \\ add %o0, 4 * 8, %o0
    \\ ldd [%l4 + 1 * 8], %f10
    \\ add %o1, 4 * 8, %o1
    \\ ldd [%l4 + 2 * 8], %f12
    \\ ldd [%l4 + 3 * 8], %f14
    \\ fmuld %f24, %f56, %f0
    \\ fmuld %f26, %f56, %f2
    \\ fmuld %f28, %f56, %f4
    \\ fmuld %f30, %f56, %f6
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f32, %f58, %f0
    \\ faddd %f10, %f2, %f10
    \\ fmuld %f34, %f58, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f36, %f58, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f38, %f58, %f6
    \\ faddd %f8, %f0, %f8
    \\ faddd %f10, %f2, %f10
    \\ faddd %f12, %f4, %f12
    \\ faddd %f14, %f6, %f14
    \\ std %f8, [%l4 + 0 * 8]
    \\ std %f10, [%l4 + 1 * 8]
    \\ std %f12, [%l4 + 2 * 8]
    \\ std %f14, [%l4 + 3 * 8]
    \\ add %l4, 4 * 8, %l4
    \\
;
