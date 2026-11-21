//! Section 3 of 4 of the text of `zgemv_d`, which `__root__.zig` joins, from OpenBLAS `kernel/sparc/zgemv_t.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ 20010:
    \\ faddd %f32, %f0, %f32
    \\ ldd [%o5 + 0 * 8], %f40
    \\ fsubd %f34, %f2, %f34
    \\ ldd [%o5 + 1 * 8] ,%f42
    \\ add %o5, %l1, %o5
    \\ faddd %f36, %f4, %f36
    \\ ldd [%o5 + 0 * 8], %f44
    \\ faddd %f38, %f6, %f38
    \\ ldd [%o5 + 1 * 8] ,%f46
    \\ add %o5, %l1, %o5
    \\ fsubd %f8, %f14, %f8
    \\ ldd [%o5 + 0 * 8], %f48
    \\ fsubd %f16, %f22, %f16
    \\ ldd [%o5 + 1 * 8] ,%f50
    \\ add %o5, %l1, %o5
    \\ fsubd %f24, %f30, %f24
    \\ ldd [%o5 + 0 * 8], %f52
    \\ fsubd %f32, %f38, %f32
    \\ ldd [%o5 + 1 * 8] ,%f54
    \\ add %o5, %l1, %o5
    \\ fsubd %f10, %f12, %f10
    \\ .word ((2 << 30) | ((0) << 25) | ( 0x36 << 19) | ( 0x60 << 5))
    \\ fsubd %f18, %f20, %f18
    \\ fsubd %f26, %f28, %f26
    \\ fsubd %f34, %f36, %f34
    \\ fmuld %f60, %f8, %f12
    \\ fmovd %f0, %f2
    \\ fmuld %f62, %f8, %f14
    \\ fmovd %f0, %f4
    \\ fmuld %f62, %f10, %f8
    \\ fmovd %f0, %f6
    \\ fmuld %f60, %f10, %f10
    \\ fmuld %f60, %f16, %f20
    \\ fmuld %f62, %f16, %f22
    \\ fmuld %f62, %f18, %f16
    \\ fmuld %f60, %f18, %f18
    \\ fmuld %f60, %f24, %f28
    \\ fmuld %f62, %f24, %f30
    \\ fmuld %f62, %f26, %f24
    \\ fmuld %f60, %f26, %f26
    \\ fmuld %f60, %f32, %f36
    \\ fsubd %f12, %f8, %f8
    \\ fmuld %f62, %f32, %f38
    \\ faddd %f10, %f14, %f10
    \\ fmuld %f62, %f34, %f32
    \\ fsubd %f20, %f16, %f16
    \\ fmuld %f60, %f34, %f34
    \\ faddd %f18, %f22, %f18
    \\ fsubd %f28, %f24, %f24
    \\ faddd %f26, %f30, %f26
    \\ fsubd %f36, %f32, %f32
    \\ faddd %f34, %f38, %f34
    \\ faddd %f40, %f8, %f40
    \\ faddd %f42, %f10, %f42
    \\ faddd %f44, %f16, %f44
    \\ faddd %f46, %f18, %f46
    \\ std %f40, [%o7 + 0 * 8]
    \\ faddd %f48, %f24, %f48
    \\ std %f42, [%o7 + 1 * 8]
    \\ faddd %f50, %f26, %f50
    \\ add %o7, %l1, %o7
    \\ std %f44, [%o7 + 0 * 8]
    \\ faddd %f52, %f32, %f52
    \\ std %f46, [%o7 + 1 * 8]
    \\ faddd %f54, %f34, %f54
    \\ add %o7, %l1, %o7
    \\ std %f48, [%o7 + 0 * 8]
    \\ fmovd %f0, %f8
    \\ add %l5, -1, %l5
    \\ std %f50, [%o7 + 1 * 8]
    \\ fmovd %f0, %f10
    \\ cmp %l5, 0
    \\ add %o7, %l1, %o7
    \\ std %f52, [%o7 + 0 * 8]
    \\ fmovd %f0, %f12
    \\ std %f54, [%o7 + 1 * 8]
    \\ fmovd %f0, %f14
    \\ add %o7, %l1, %o7
    \\ fmovd %f0, %f16
    \\ bg %icc, 20005b
    \\ fmovd %f0, %f18
    \\ 20011:
    \\ .word ((2 << 30) | ((0) << 25) | ( 0x36 << 19) | ( 0x60 << 5))
    \\ and %i1, 2, %l5
    \\ cmp %l5, 0
    \\ fmovd %f0, %f8
    \\ ble %icc, 20018f
    \\ fmovd %f0, %f10
    \\ sra %l6, 2, %l3
    \\ fmovd %f0, %f2
    \\ add %i5, %i2, %o1
    \\ fmovd %f0, %f12
    \\ mov %i5, %o0
    \\ fmovd %f0, %f4
    \\ cmp %l3, 0
    \\ fmovd %f0, %f14
    \\ fmovd %f0, %f16
    \\ fmovd %f0, %f18
    \\ fmovd %f0, %f20
    \\ fmovd %f0, %f22
    \\ add %o1, %i2, %i5
    \\ fmovd %f0, %f6
    \\ ble %icc, 20014f
    \\ mov %l7, %o4
    \\ ldd [%o0 + 0 * 8], %f40
    \\ ldd [%o0 + 1 * 8], %f42
    \\ ldd [%o0 + 2 * 8], %f48
    \\ ldd [%o0 + 3 * 8], %f50
    \\ add %o0, 4 * 8, %o0
    \\ ldd [%o1 + 0 * 8], %f44
    \\ ldd [%o1 + 1 * 8], %f46
    \\ ldd [%o1 + 2 * 8], %f52
    \\ ldd [%o1 + 3 * 8], %f54
    \\ add %o1, 4 * 8, %o1
    \\ ldd [%o4 + 0 * 8], %f56
    \\ add %l3, -1, %l3
    \\ ldd [%o4 + 1 * 8], %f58
    \\ cmp %l3, 0
    \\ ldd [%o4 + 2 * 8], %f60
    \\ ldd [%o4 + 3 * 8], %f62
    \\ ble %icc, 20013f
    \\ add %o4, 4 * 8, %o4
    \\ 20012:
    \\ prefetch [%o0 + 18 * 8], 1
    \\ faddd %f16, %f0, %f16
    \\ fmuld %f40, %f56, %f0
    \\ fsubd %f18, %f2, %f18
    \\ fmuld %f40, %f58, %f2
    \\ ldd [%o0 + 0 * 8], %f40
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f42, %f56, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f42, %f58, %f6
    \\ ldd [%o0 + 1 * 8], %f42
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f44, %f56, %f0
    \\ fsubd %f10, %f2, %f10
    \\ fmuld %f44, %f58, %f2
    \\ ldd [%o1 + 0 * 8], %f44
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f46, %f56, %f4
    \\ ldd [%o4 + 0 * 8], %f56
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f46, %f58, %f6
    \\ ldd [%o1 + 1 * 8], %f46
    \\ faddd %f16, %f0, %f16
    \\ ldd [%o4 + 1 * 8], %f58
    \\ fmuld %f48, %f60, %f0
    \\ fsubd %f18, %f2, %f18
    \\ fmuld %f48, %f62, %f2
    \\ ldd [%o0 + 2 * 8], %f48
    \\ faddd %f20, %f4, %f20
    \\ add %l3, -1, %l3
    \\ fmuld %f50, %f60, %f4
    \\ faddd %f22, %f6, %f22
    \\ cmp %l3, 0
    \\ fmuld %f50, %f62, %f6
    \\ ldd [%o0 + 3 * 8], %f50
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f52, %f60, %f0
    \\ fsubd %f10, %f2, %f10
    \\ fmuld %f52, %f62, %f2
    \\ ldd [%o1 + 2 * 8], %f52
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f54, %f60, %f4
    \\ ldd [%o4 + 2 * 8], %f60
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f54, %f62, %f6
    \\ ldd [%o1 + 3 * 8], %f54
    \\ prefetch [%o1 + 18 * 8], 1
    \\ faddd %f16, %f0, %f16
    \\ ldd [%o4 + 3 * 8], %f62
    \\ fmuld %f40, %f56, %f0
    \\ fsubd %f18, %f2, %f18
    \\ fmuld %f40, %f58, %f2
    \\ ldd [%o0 + 4 * 8], %f40
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f42, %f56, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f42, %f58, %f6
    \\ ldd [%o0 + 5 * 8], %f42
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f44, %f56, %f0
    \\ fsubd %f10, %f2, %f10
    \\ fmuld %f44, %f58, %f2
    \\ ldd [%o1 + 4 * 8], %f44
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f46, %f56, %f4
    \\ ldd [%o4 + 4 * 8], %f56
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f46, %f58, %f6
    \\ ldd [%o1 + 5 * 8], %f46
    \\ faddd %f16, %f0, %f16
    \\ ldd [%o4 + 5 * 8], %f58
    \\ fmuld %f48, %f60, %f0
    \\ fsubd %f18, %f2, %f18
    \\ fmuld %f48, %f62, %f2
    \\ ldd [%o0 + 6 * 8], %f48
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f50, %f60, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f50, %f62, %f6
    \\ ldd [%o0 + 7 * 8], %f50
    \\ add %o0, 8 * 8, %o0
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f52, %f60, %f0
    \\ fsubd %f10, %f2, %f10
    \\ fmuld %f52, %f62, %f2
    \\ ldd [%o1 + 6 * 8], %f52
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f54, %f60, %f4
    \\ ldd [%o4 + 6 * 8], %f60
    \\ faddd %f14, %f6, %f14
    \\ add %o4, 8 * 8, %o4
    \\ fmuld %f54, %f62, %f6
    \\ ldd [%o1 + 7 * 8], %f54
    \\ add %o1, 8 * 8, %o1
    \\ bg,pn %icc, 20012b
    \\ ldd [%o4 - 1 * 8], %f62
    \\ 20013:
    \\ faddd %f16, %f0, %f16
    \\ fmuld %f40, %f56, %f0
    \\ fsubd %f18, %f2, %f18
    \\ fmuld %f40, %f58, %f2
    \\ ldd [%o0 + 0 * 8], %f40
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f42, %f56, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f42, %f58, %f6
    \\ ldd [%o0 + 1 * 8], %f42
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f44, %f56, %f0
    \\ fsubd %f10, %f2, %f10
    \\ fmuld %f44, %f58, %f2
    \\ ldd [%o1 + 0 * 8], %f44
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f46, %f56, %f4
    \\ ldd [%o4 + 0 * 8], %f56
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f46, %f58, %f6
    \\ ldd [%o1 + 1 * 8], %f46
    \\ faddd %f16, %f0, %f16
    \\ ldd [%o4 + 1 * 8], %f58
    \\ fmuld %f48, %f60, %f0
    \\ fsubd %f18, %f2, %f18
    \\ fmuld %f48, %f62, %f2
    \\ ldd [%o0 + 2 * 8], %f48
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f50, %f60, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f50, %f62, %f6
    \\ ldd [%o0 + 3 * 8], %f50
    \\ add %o0, 4 * 8, %o0
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f52, %f60, %f0
    \\ fsubd %f10, %f2, %f10
    \\ fmuld %f52, %f62, %f2
    \\ ldd [%o1 + 2 * 8], %f52
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f54, %f60, %f4
    \\ ldd [%o4 + 2 * 8], %f60
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f54, %f62, %f6
    \\ ldd [%o1 + 3 * 8], %f54
    \\ add %o1, 4 * 8, %o1
    \\ faddd %f16, %f0, %f16
    \\ ldd [%o4 + 3 * 8], %f62
    \\ add %o4, 4 * 8, %o4
    \\ fmuld %f40, %f56, %f0
    \\ fsubd %f18, %f2, %f18
    \\ fmuld %f40, %f58, %f2
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f42, %f56, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f42, %f58, %f6
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f44, %f56, %f0
    \\ fsubd %f10, %f2, %f10
    \\ fmuld %f44, %f58, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f46, %f56, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f46, %f58, %f6
    \\ faddd %f16, %f0, %f16
    \\ fmuld %f48, %f60, %f0
    \\ fsubd %f18, %f2, %f18
    \\ fmuld %f48, %f62, %f2
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f50, %f60, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f50, %f62, %f6
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f52, %f60, %f0
    \\ fsubd %f10, %f2, %f10
    \\ fmuld %f52, %f62, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f54, %f60, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f54, %f62, %f6
    \\
;
