//! Section 2 of 4 of the text of `zgemv_d`, which `__root__.zig` joins, from OpenBLAS `kernel/sparc/zgemv_t.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ 20007:
    \\ faddd %f32, %f0, %f32
    \\ nop
    \\ fmuld %f40, %f56, %f0
    \\ ldd [%o4 - 1 * 8], %f62
    \\ fsubd %f34, %f2, %f34
    \\ nop
    \\ fmuld %f40, %f58, %f2
    \\ ldd [%o0 + 0 * 8], %f40
    \\ faddd %f36, %f4, %f36
    \\ nop
    \\ fmuld %f42, %f56, %f4
    \\ ldd [%o4 - 1 * 8], %f62
    \\ faddd %f38, %f6, %f38
    \\ nop
    \\ fmuld %f42, %f58, %f6
    \\ ldd [%o0 + 1 * 8], %f42
    \\ faddd %f8, %f0, %f8
    \\ nop
    \\ fmuld %f44, %f56, %f0
    \\ nop
    \\ fsubd %f10, %f2, %f10
    \\ nop
    \\ fmuld %f44, %f58, %f2
    \\ ldd [%o1 + 0 * 8], %f44
    \\ faddd %f12, %f4, %f12
    \\ nop
    \\ fmuld %f46, %f56, %f4
    \\ nop
    \\ faddd %f14, %f6, %f14
    \\ nop
    \\ fmuld %f46, %f58, %f6
    \\ ldd [%o1 + 1 * 8], %f46
    \\ faddd %f16, %f0, %f16
    \\ nop
    \\ fmuld %f48, %f56, %f0
    \\ nop
    \\ fsubd %f18, %f2, %f18
    \\ nop
    \\ fmuld %f48, %f58, %f2
    \\ ldd [%o2 + 0 * 8], %f48
    \\ faddd %f20, %f4, %f20
    \\ nop
    \\ fmuld %f50, %f56, %f4
    \\ nop
    \\ faddd %f22, %f6, %f22
    \\ nop
    \\ fmuld %f50, %f58, %f6
    \\ ldd [%o2 + 1 * 8], %f50
    \\ faddd %f24, %f0, %f24
    \\ nop
    \\ fmuld %f52, %f56, %f0
    \\ nop
    \\ fsubd %f26, %f2, %f26
    \\ nop
    \\ fmuld %f52, %f58, %f2
    \\ ldd [%o3 + 0 * 8], %f52
    \\ faddd %f28, %f4, %f28
    \\ nop
    \\ fmuld %f54, %f56, %f4
    \\ ldd [%o4 + 0 * 8], %f56
    \\ faddd %f30, %f6, %f30
    \\ nop
    \\ fmuld %f54, %f58, %f6
    \\ ldd [%o3 + 1 * 8], %f54
    \\ faddd %f32, %f0, %f32
    \\ nop
    \\ fmuld %f40, %f60, %f0
    \\ ldd [%o4 + 1 * 8], %f58
    \\ fsubd %f34, %f2, %f34
    \\ nop
    \\ fmuld %f40, %f62, %f2
    \\ ldd [%o0 + 2 * 8], %f40
    \\ faddd %f36, %f4, %f36
    \\ nop
    \\ fmuld %f42, %f60, %f4
    \\ nop
    \\ faddd %f38, %f6, %f38
    \\ nop
    \\ fmuld %f42, %f62, %f6
    \\ ldd [%o0 + 3 * 8], %f42
    \\ faddd %f8, %f0, %f8
    \\ nop
    \\ fmuld %f44, %f60, %f0
    \\ nop
    \\ fsubd %f10, %f2, %f10
    \\ nop
    \\ fmuld %f44, %f62, %f2
    \\ ldd [%o1 + 2 * 8], %f44
    \\ faddd %f12, %f4, %f12
    \\ nop
    \\ fmuld %f46, %f60, %f4
    \\ nop
    \\ faddd %f14, %f6, %f14
    \\ nop
    \\ fmuld %f46, %f62, %f6
    \\ ldd [%o1 + 3 * 8], %f46
    \\ faddd %f16, %f0, %f16
    \\ nop
    \\ fmuld %f48, %f60, %f0
    \\ nop
    \\ fsubd %f18, %f2, %f18
    \\ nop
    \\ fmuld %f48, %f62, %f2
    \\ ldd [%o2 + 2 * 8], %f48
    \\ faddd %f20, %f4, %f20
    \\ nop
    \\ fmuld %f50, %f60, %f4
    \\ nop
    \\ faddd %f22, %f6, %f22
    \\ nop
    \\ fmuld %f50, %f62, %f6
    \\ ldd [%o2 + 3 * 8], %f50
    \\ faddd %f24, %f0, %f24
    \\ nop
    \\ fmuld %f52, %f60, %f0
    \\ nop
    \\ fsubd %f26, %f2, %f26
    \\ nop
    \\ fmuld %f52, %f62, %f2
    \\ ldd [%o3 + 2 * 8], %f52
    \\ faddd %f28, %f4, %f28
    \\ nop
    \\ fmuld %f54, %f60, %f4
    \\ ldd [%o4 + 2 * 8], %f60
    \\ faddd %f30, %f6, %f30
    \\ nop
    \\ fmuld %f54, %f62, %f6
    \\ ldd [%o3 + 3 * 8], %f54
    \\ faddd %f32, %f0, %f32
    \\ nop
    \\ fmuld %f40, %f56, %f0
    \\ ldd [%o4 + 3 * 8], %f62
    \\ fsubd %f34, %f2, %f34
    \\ add %o4, 4 * 8, %o4
    \\ fmuld %f40, %f58, %f2
    \\ ldd [%o0 + 4 * 8], %f40
    \\ faddd %f36, %f4, %f36
    \\ nop
    \\ fmuld %f42, %f56, %f4
    \\ nop
    \\ faddd %f38, %f6, %f38
    \\ nop
    \\ fmuld %f42, %f58, %f6
    \\ ldd [%o0 + 5 * 8], %f42
    \\ faddd %f8, %f0, %f8
    \\ add %o0, 6 * 8, %o0
    \\ fmuld %f44, %f56, %f0
    \\ nop
    \\ fsubd %f10, %f2, %f10
    \\ nop
    \\ fmuld %f44, %f58, %f2
    \\ ldd [%o1 + 4 * 8], %f44
    \\ faddd %f12, %f4, %f12
    \\ nop
    \\ fmuld %f46, %f56, %f4
    \\ nop
    \\ faddd %f14, %f6, %f14
    \\ nop
    \\ fmuld %f46, %f58, %f6
    \\ ldd [%o1 + 5 * 8], %f46
    \\ faddd %f16, %f0, %f16
    \\ add %o1, 6 * 8, %o1
    \\ fmuld %f48, %f56, %f0
    \\ nop
    \\ fsubd %f18, %f2, %f18
    \\ nop
    \\ fmuld %f48, %f58, %f2
    \\ ldd [%o2 + 4 * 8], %f48
    \\ faddd %f20, %f4, %f20
    \\ nop
    \\ fmuld %f50, %f56, %f4
    \\ nop
    \\ faddd %f22, %f6, %f22
    \\ nop
    \\ fmuld %f50, %f58, %f6
    \\ ldd [%o2 + 5 * 8], %f50
    \\ faddd %f24, %f0, %f24
    \\ add %o2, 6 * 8, %o2
    \\ fmuld %f52, %f56, %f0
    \\ nop
    \\ fsubd %f26, %f2, %f26
    \\ nop
    \\ fmuld %f52, %f58, %f2
    \\ ldd [%o3 + 4 * 8], %f52
    \\ faddd %f28, %f4, %f28
    \\ nop
    \\ fmuld %f54, %f56, %f4
    \\ nop
    \\ faddd %f30, %f6, %f30
    \\ nop
    \\ fmuld %f54, %f58, %f6
    \\ ldd [%o3 + 5 * 8], %f54
    \\ faddd %f32, %f0, %f32
    \\ add %o3, 6 * 8, %o3
    \\ fmuld %f40, %f60, %f0
    \\ nop
    \\ fsubd %f34, %f2, %f34
    \\ nop
    \\ fmuld %f40, %f62, %f2
    \\ nop
    \\ faddd %f36, %f4, %f36
    \\ fmuld %f42, %f60, %f4
    \\ faddd %f38, %f6, %f38
    \\ fmuld %f42, %f62, %f6
    \\ faddd %f8, %f0, %f8
    \\ fmuld %f44, %f60, %f0
    \\ fsubd %f10, %f2, %f10
    \\ fmuld %f44, %f62, %f2
    \\ faddd %f12, %f4, %f12
    \\ fmuld %f46, %f60, %f4
    \\ faddd %f14, %f6, %f14
    \\ fmuld %f46, %f62, %f6
    \\ faddd %f16, %f0, %f16
    \\ fmuld %f48, %f60, %f0
    \\ fsubd %f18, %f2, %f18
    \\ fmuld %f48, %f62, %f2
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f50, %f60, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f50, %f62, %f6
    \\ faddd %f24, %f0, %f24
    \\ fmuld %f52, %f60, %f0
    \\ fsubd %f26, %f2, %f26
    \\ fmuld %f52, %f62, %f2
    \\ faddd %f28, %f4, %f28
    \\ fmuld %f54, %f60, %f4
    \\ faddd %f30, %f6, %f30
    \\ fmuld %f54, %f62, %f6
    \\ 20008:
    \\ andcc %l6, 3, %l3
    \\ .if \bp0
    \\ ldd [%sp + 128 + 16], %f60
    \\ .endif
    \\ .if \bp1
    \\ ldd [%sp + 2423 + 32], %f60
    \\ .endif
    \\ mov %o5, %o7
    \\ ble,pn %icc, 20010f
    \\ .if \bp0
    \\ ldd [%sp + 128 + 24], %f62
    \\ .endif
    \\ .if \bp1
    \\ ldd [%sp + 2423 + 40], %f62
    \\ .endif
    \\ 20009:
    \\ ldd [%o0 + 0 * 8], %f40
    \\ ldd [%o0 + 1 * 8], %f42
    \\ add %o0, 2 * 8, %o0
    \\ ldd [%o4 + 0 * 8], %f56
    \\ ldd [%o4 + 1 * 8], %f58
    \\ add %o4, 2 * 8, %o4
    \\ ldd [%o1 + 0 * 8], %f44
    \\ ldd [%o1 + 1 * 8], %f46
    \\ add %o1, 2 * 8, %o1
    \\ ldd [%o2 + 0 * 8], %f48
    \\ ldd [%o2 + 1 * 8], %f50
    \\ add %o2, 2 * 8, %o2
    \\ ldd [%o3 + 0 * 8], %f52
    \\ ldd [%o3 + 1 * 8], %f54
    \\ add %o3, 2 * 8, %o3
    \\ faddd %f32, %f0, %f32
    \\ fmuld %f40, %f56, %f0
    \\ fsubd %f34, %f2, %f34
    \\ fmuld %f40, %f58, %f2
    \\ faddd %f36, %f4, %f36
    \\ fmuld %f42, %f56, %f4
    \\ faddd %f38, %f6, %f38
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
    \\ fmuld %f48, %f56, %f0
    \\ fsubd %f18, %f2, %f18
    \\ fmuld %f48, %f58, %f2
    \\ faddd %f20, %f4, %f20
    \\ fmuld %f50, %f56, %f4
    \\ faddd %f22, %f6, %f22
    \\ fmuld %f50, %f58, %f6
    \\ faddd %f24, %f0, %f24
    \\ fmuld %f52, %f56, %f0
    \\ fsubd %f26, %f2, %f26
    \\ fmuld %f52, %f58, %f2
    \\ faddd %f28, %f4, %f28
    \\ fmuld %f54, %f56, %f4
    \\ faddd %f30, %f6, %f30
    \\ fmuld %f54, %f58, %f6
    \\ deccc %l3
    \\ bg %icc, 20009b
    \\ nop
    \\
;
