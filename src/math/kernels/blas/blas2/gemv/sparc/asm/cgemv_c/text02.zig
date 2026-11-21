//! Section 2 of 4 of the text of `cgemv_c`, which `__root__.zig` joins, from OpenBLAS `kernel/sparc/zgemv_t.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ 20007:
    \\ fadds %f16, %f0, %f16
    \\ nop
    \\ fmuls %f20, %f28, %f0
    \\ ld [%o4 - 1 * 4], %f31
    \\ fadds %f17, %f1, %f17
    \\ nop
    \\ fmuls %f20, %f29, %f1
    \\ ld [%o0 + 0 * 4], %f20
    \\ fadds %f18, %f2, %f18
    \\ nop
    \\ fmuls %f21, %f28, %f2
    \\ ld [%o4 - 1 * 4], %f31
    \\ fadds %f19, %f3, %f19
    \\ nop
    \\ fmuls %f21, %f29, %f3
    \\ ld [%o0 + 1 * 4], %f21
    \\ fadds %f4, %f0, %f4
    \\ nop
    \\ fmuls %f22, %f28, %f0
    \\ nop
    \\ fadds %f5, %f1, %f5
    \\ nop
    \\ fmuls %f22, %f29, %f1
    \\ ld [%o1 + 0 * 4], %f22
    \\ fadds %f6, %f2, %f6
    \\ nop
    \\ fmuls %f23, %f28, %f2
    \\ nop
    \\ fadds %f7, %f3, %f7
    \\ nop
    \\ fmuls %f23, %f29, %f3
    \\ ld [%o1 + 1 * 4], %f23
    \\ fadds %f8, %f0, %f8
    \\ nop
    \\ fmuls %f24, %f28, %f0
    \\ nop
    \\ fadds %f9, %f1, %f9
    \\ nop
    \\ fmuls %f24, %f29, %f1
    \\ ld [%o2 + 0 * 4], %f24
    \\ fadds %f10, %f2, %f10
    \\ nop
    \\ fmuls %f25, %f28, %f2
    \\ nop
    \\ fadds %f11, %f3, %f11
    \\ nop
    \\ fmuls %f25, %f29, %f3
    \\ ld [%o2 + 1 * 4], %f25
    \\ fadds %f12, %f0, %f12
    \\ nop
    \\ fmuls %f26, %f28, %f0
    \\ nop
    \\ fadds %f13, %f1, %f13
    \\ nop
    \\ fmuls %f26, %f29, %f1
    \\ ld [%o3 + 0 * 4], %f26
    \\ fadds %f14, %f2, %f14
    \\ nop
    \\ fmuls %f27, %f28, %f2
    \\ ld [%o4 + 0 * 4], %f28
    \\ fadds %f15, %f3, %f15
    \\ nop
    \\ fmuls %f27, %f29, %f3
    \\ ld [%o3 + 1 * 4], %f27
    \\ fadds %f16, %f0, %f16
    \\ nop
    \\ fmuls %f20, %f30, %f0
    \\ ld [%o4 + 1 * 4], %f29
    \\ fadds %f17, %f1, %f17
    \\ nop
    \\ fmuls %f20, %f31, %f1
    \\ ld [%o0 + 2 * 4], %f20
    \\ fadds %f18, %f2, %f18
    \\ nop
    \\ fmuls %f21, %f30, %f2
    \\ nop
    \\ fadds %f19, %f3, %f19
    \\ nop
    \\ fmuls %f21, %f31, %f3
    \\ ld [%o0 + 3 * 4], %f21
    \\ fadds %f4, %f0, %f4
    \\ nop
    \\ fmuls %f22, %f30, %f0
    \\ nop
    \\ fadds %f5, %f1, %f5
    \\ nop
    \\ fmuls %f22, %f31, %f1
    \\ ld [%o1 + 2 * 4], %f22
    \\ fadds %f6, %f2, %f6
    \\ nop
    \\ fmuls %f23, %f30, %f2
    \\ nop
    \\ fadds %f7, %f3, %f7
    \\ nop
    \\ fmuls %f23, %f31, %f3
    \\ ld [%o1 + 3 * 4], %f23
    \\ fadds %f8, %f0, %f8
    \\ nop
    \\ fmuls %f24, %f30, %f0
    \\ nop
    \\ fadds %f9, %f1, %f9
    \\ nop
    \\ fmuls %f24, %f31, %f1
    \\ ld [%o2 + 2 * 4], %f24
    \\ fadds %f10, %f2, %f10
    \\ nop
    \\ fmuls %f25, %f30, %f2
    \\ nop
    \\ fadds %f11, %f3, %f11
    \\ nop
    \\ fmuls %f25, %f31, %f3
    \\ ld [%o2 + 3 * 4], %f25
    \\ fadds %f12, %f0, %f12
    \\ nop
    \\ fmuls %f26, %f30, %f0
    \\ nop
    \\ fadds %f13, %f1, %f13
    \\ nop
    \\ fmuls %f26, %f31, %f1
    \\ ld [%o3 + 2 * 4], %f26
    \\ fadds %f14, %f2, %f14
    \\ nop
    \\ fmuls %f27, %f30, %f2
    \\ ld [%o4 + 2 * 4], %f30
    \\ fadds %f15, %f3, %f15
    \\ nop
    \\ fmuls %f27, %f31, %f3
    \\ ld [%o3 + 3 * 4], %f27
    \\ fadds %f16, %f0, %f16
    \\ nop
    \\ fmuls %f20, %f28, %f0
    \\ ld [%o4 + 3 * 4], %f31
    \\ fadds %f17, %f1, %f17
    \\ add %o4, 4 * 4, %o4
    \\ fmuls %f20, %f29, %f1
    \\ ld [%o0 + 4 * 4], %f20
    \\ fadds %f18, %f2, %f18
    \\ nop
    \\ fmuls %f21, %f28, %f2
    \\ nop
    \\ fadds %f19, %f3, %f19
    \\ nop
    \\ fmuls %f21, %f29, %f3
    \\ ld [%o0 + 5 * 4], %f21
    \\ fadds %f4, %f0, %f4
    \\ add %o0, 6 * 4, %o0
    \\ fmuls %f22, %f28, %f0
    \\ nop
    \\ fadds %f5, %f1, %f5
    \\ nop
    \\ fmuls %f22, %f29, %f1
    \\ ld [%o1 + 4 * 4], %f22
    \\ fadds %f6, %f2, %f6
    \\ nop
    \\ fmuls %f23, %f28, %f2
    \\ nop
    \\ fadds %f7, %f3, %f7
    \\ nop
    \\ fmuls %f23, %f29, %f3
    \\ ld [%o1 + 5 * 4], %f23
    \\ fadds %f8, %f0, %f8
    \\ add %o1, 6 * 4, %o1
    \\ fmuls %f24, %f28, %f0
    \\ nop
    \\ fadds %f9, %f1, %f9
    \\ nop
    \\ fmuls %f24, %f29, %f1
    \\ ld [%o2 + 4 * 4], %f24
    \\ fadds %f10, %f2, %f10
    \\ nop
    \\ fmuls %f25, %f28, %f2
    \\ nop
    \\ fadds %f11, %f3, %f11
    \\ nop
    \\ fmuls %f25, %f29, %f3
    \\ ld [%o2 + 5 * 4], %f25
    \\ fadds %f12, %f0, %f12
    \\ add %o2, 6 * 4, %o2
    \\ fmuls %f26, %f28, %f0
    \\ nop
    \\ fadds %f13, %f1, %f13
    \\ nop
    \\ fmuls %f26, %f29, %f1
    \\ ld [%o3 + 4 * 4], %f26
    \\ fadds %f14, %f2, %f14
    \\ nop
    \\ fmuls %f27, %f28, %f2
    \\ nop
    \\ fadds %f15, %f3, %f15
    \\ nop
    \\ fmuls %f27, %f29, %f3
    \\ ld [%o3 + 5 * 4], %f27
    \\ fadds %f16, %f0, %f16
    \\ add %o3, 6 * 4, %o3
    \\ fmuls %f20, %f30, %f0
    \\ nop
    \\ fadds %f17, %f1, %f17
    \\ nop
    \\ fmuls %f20, %f31, %f1
    \\ nop
    \\ fadds %f18, %f2, %f18
    \\ fmuls %f21, %f30, %f2
    \\ fadds %f19, %f3, %f19
    \\ fmuls %f21, %f31, %f3
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f22, %f30, %f0
    \\ fadds %f5, %f1, %f5
    \\ fmuls %f22, %f31, %f1
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f23, %f30, %f2
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f23, %f31, %f3
    \\ fadds %f8, %f0, %f8
    \\ fmuls %f24, %f30, %f0
    \\ fadds %f9, %f1, %f9
    \\ fmuls %f24, %f31, %f1
    \\ fadds %f10, %f2, %f10
    \\ fmuls %f25, %f30, %f2
    \\ fadds %f11, %f3, %f11
    \\ fmuls %f25, %f31, %f3
    \\ fadds %f12, %f0, %f12
    \\ fmuls %f26, %f30, %f0
    \\ fadds %f13, %f1, %f13
    \\ fmuls %f26, %f31, %f1
    \\ fadds %f14, %f2, %f14
    \\ fmuls %f27, %f30, %f2
    \\ fadds %f15, %f3, %f15
    \\ fmuls %f27, %f31, %f3
    \\ 20008:
    \\ andcc %l6, 3, %l3
    \\ .if \bp0
    \\ ld [%sp + 128 + 16], %f30
    \\ .endif
    \\ .if \bp1
    \\ ld [%sp + 2423 + 32], %f30
    \\ .endif
    \\ mov %o5, %o7
    \\ ble,pn %icc, 20010f
    \\ .if \bp0
    \\ ld [%sp + 128 + 20], %f31
    \\ .endif
    \\ .if \bp1
    \\ ld [%sp + 2423 + 40], %f31
    \\ .endif
    \\ 20009:
    \\ ld [%o0 + 0 * 4], %f20
    \\ ld [%o0 + 1 * 4], %f21
    \\ add %o0, 2 * 4, %o0
    \\ ld [%o4 + 0 * 4], %f28
    \\ ld [%o4 + 1 * 4], %f29
    \\ add %o4, 2 * 4, %o4
    \\ ld [%o1 + 0 * 4], %f22
    \\ ld [%o1 + 1 * 4], %f23
    \\ add %o1, 2 * 4, %o1
    \\ ld [%o2 + 0 * 4], %f24
    \\ ld [%o2 + 1 * 4], %f25
    \\ add %o2, 2 * 4, %o2
    \\ ld [%o3 + 0 * 4], %f26
    \\ ld [%o3 + 1 * 4], %f27
    \\ add %o3, 2 * 4, %o3
    \\ fadds %f16, %f0, %f16
    \\ fmuls %f20, %f28, %f0
    \\ fadds %f17, %f1, %f17
    \\ fmuls %f20, %f29, %f1
    \\ fadds %f18, %f2, %f18
    \\ fmuls %f21, %f28, %f2
    \\ fadds %f19, %f3, %f19
    \\ fmuls %f21, %f29, %f3
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f22, %f28, %f0
    \\ fadds %f5, %f1, %f5
    \\ fmuls %f22, %f29, %f1
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f23, %f28, %f2
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f23, %f29, %f3
    \\ fadds %f8, %f0, %f8
    \\ fmuls %f24, %f28, %f0
    \\ fadds %f9, %f1, %f9
    \\ fmuls %f24, %f29, %f1
    \\ fadds %f10, %f2, %f10
    \\ fmuls %f25, %f28, %f2
    \\ fadds %f11, %f3, %f11
    \\ fmuls %f25, %f29, %f3
    \\ fadds %f12, %f0, %f12
    \\ fmuls %f26, %f28, %f0
    \\ fadds %f13, %f1, %f13
    \\ fmuls %f26, %f29, %f1
    \\ fadds %f14, %f2, %f14
    \\ fmuls %f27, %f28, %f2
    \\ fadds %f15, %f3, %f15
    \\ fmuls %f27, %f29, %f3
    \\ deccc %l3
    \\ bg %icc, 20009b
    \\ nop
    \\
;
