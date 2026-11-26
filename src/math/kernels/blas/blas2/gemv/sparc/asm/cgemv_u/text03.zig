//! Section 3 of 4 of the text of `cgemv_u`, which `__root__.zig` joins, from OpenBLAS `kernel/sparc/zgemv_t.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ 20010:
    \\ fadds %f16, %f0, %f16
    \\ ld [%o5 + 0 * 4], %f20
    \\ fsubs %f17, %f1, %f17
    \\ ld [%o5 + 1 * 4] ,%f21
    \\ add %o5, %l1, %o5
    \\ fadds %f18, %f2, %f18
    \\ ld [%o5 + 0 * 4], %f22
    \\ fadds %f19, %f3, %f19
    \\ ld [%o5 + 1 * 4] ,%f23
    \\ add %o5, %l1, %o5
    \\ fadds %f4, %f7, %f4
    \\ ld [%o5 + 0 * 4], %f24
    \\ fadds %f8, %f11, %f8
    \\ ld [%o5 + 1 * 4] ,%f25
    \\ add %o5, %l1, %o5
    \\ fadds %f12, %f15, %f12
    \\ ld [%o5 + 0 * 4], %f26
    \\ fadds %f16, %f19, %f16
    \\ ld [%o5 + 1 * 4] ,%f27
    \\ add %o5, %l1, %o5
    \\ fadds %f5, %f6, %f5
    \\ .word ((2 << 30) | ((0) << 25) | ( 0x36 << 19) | ( 0x61 << 5))
    \\ fadds %f9, %f10, %f9
    \\ fadds %f13, %f14, %f13
    \\ fadds %f17, %f18, %f17
    \\ fmuls %f30, %f4, %f6
    \\ fmovs %f0, %f1
    \\ fmuls %f31, %f4, %f7
    \\ fmovs %f0, %f2
    \\ fmuls %f31, %f5, %f4
    \\ fmovs %f0, %f3
    \\ fmuls %f30, %f5, %f5
    \\ fmuls %f30, %f8, %f10
    \\ fmuls %f31, %f8, %f11
    \\ fmuls %f31, %f9, %f8
    \\ fmuls %f30, %f9, %f9
    \\ fmuls %f30, %f12, %f14
    \\ fmuls %f31, %f12, %f15
    \\ fmuls %f31, %f13, %f12
    \\ fmuls %f30, %f13, %f13
    \\ fmuls %f30, %f16, %f18
    \\ fsubs %f6, %f4, %f4
    \\ fmuls %f31, %f16, %f19
    \\ fadds %f5, %f7, %f5
    \\ fmuls %f31, %f17, %f16
    \\ fsubs %f10, %f8, %f8
    \\ fmuls %f30, %f17, %f17
    \\ fadds %f9, %f11, %f9
    \\ fsubs %f14, %f12, %f12
    \\ fadds %f13, %f15, %f13
    \\ fsubs %f18, %f16, %f16
    \\ fadds %f17, %f19, %f17
    \\ fadds %f20, %f4, %f20
    \\ fadds %f21, %f5, %f21
    \\ fadds %f22, %f8, %f22
    \\ fadds %f23, %f9, %f23
    \\ st %f20, [%o7 + 0 * 4]
    \\ fadds %f24, %f12, %f24
    \\ st %f21, [%o7 + 1 * 4]
    \\ fadds %f25, %f13, %f25
    \\ add %o7, %l1, %o7
    \\ st %f22, [%o7 + 0 * 4]
    \\ fadds %f26, %f16, %f26
    \\ st %f23, [%o7 + 1 * 4]
    \\ fadds %f27, %f17, %f27
    \\ add %o7, %l1, %o7
    \\ st %f24, [%o7 + 0 * 4]
    \\ fmovs %f0, %f4
    \\ add %l5, -1, %l5
    \\ st %f25, [%o7 + 1 * 4]
    \\ fmovs %f0, %f5
    \\ cmp %l5, 0
    \\ add %o7, %l1, %o7
    \\ st %f26, [%o7 + 0 * 4]
    \\ fmovs %f0, %f6
    \\ st %f27, [%o7 + 1 * 4]
    \\ fmovs %f0, %f7
    \\ add %o7, %l1, %o7
    \\ fmovs %f0, %f8
    \\ bg %icc, 20005b
    \\ fmovs %f0, %f9
    \\ 20011:
    \\ .word ((2 << 30) | ((0) << 25) | ( 0x36 << 19) | ( 0x61 << 5))
    \\ and %i1, 2, %l5
    \\ cmp %l5, 0
    \\ fmovs %f0, %f4
    \\ ble %icc, 20018f
    \\ fmovs %f0, %f5
    \\ sra %l6, 2, %l3
    \\ fmovs %f0, %f1
    \\ add %i5, %i2, %o1
    \\ fmovs %f0, %f6
    \\ mov %i5, %o0
    \\ fmovs %f0, %f2
    \\ cmp %l3, 0
    \\ fmovs %f0, %f7
    \\ fmovs %f0, %f8
    \\ fmovs %f0, %f9
    \\ fmovs %f0, %f10
    \\ fmovs %f0, %f11
    \\ add %o1, %i2, %i5
    \\ fmovs %f0, %f3
    \\ ble %icc, 20014f
    \\ mov %l7, %o4
    \\ ld [%o0 + 0 * 4], %f20
    \\ ld [%o0 + 1 * 4], %f21
    \\ ld [%o0 + 2 * 4], %f24
    \\ ld [%o0 + 3 * 4], %f25
    \\ add %o0, 4 * 4, %o0
    \\ ld [%o1 + 0 * 4], %f22
    \\ ld [%o1 + 1 * 4], %f23
    \\ ld [%o1 + 2 * 4], %f26
    \\ ld [%o1 + 3 * 4], %f27
    \\ add %o1, 4 * 4, %o1
    \\ ld [%o4 + 0 * 4], %f28
    \\ add %l3, -1, %l3
    \\ ld [%o4 + 1 * 4], %f29
    \\ cmp %l3, 0
    \\ ld [%o4 + 2 * 4], %f30
    \\ ld [%o4 + 3 * 4], %f31
    \\ ble %icc, 20013f
    \\ add %o4, 4 * 4, %o4
    \\ 20012:
    \\ prefetch [%o0 + 36 * 4], 1
    \\ fadds %f8, %f0, %f8
    \\ fmuls %f20, %f28, %f0
    \\ fsubs %f9, %f1, %f9
    \\ fmuls %f20, %f29, %f1
    \\ ld [%o0 + 0 * 4], %f20
    \\ fadds %f10, %f2, %f10
    \\ fmuls %f21, %f28, %f2
    \\ fadds %f11, %f3, %f11
    \\ fmuls %f21, %f29, %f3
    \\ ld [%o0 + 1 * 4], %f21
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f22, %f28, %f0
    \\ fsubs %f5, %f1, %f5
    \\ fmuls %f22, %f29, %f1
    \\ ld [%o1 + 0 * 4], %f22
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f23, %f28, %f2
    \\ ld [%o4 + 0 * 4], %f28
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f23, %f29, %f3
    \\ ld [%o1 + 1 * 4], %f23
    \\ fadds %f8, %f0, %f8
    \\ ld [%o4 + 1 * 4], %f29
    \\ fmuls %f24, %f30, %f0
    \\ fsubs %f9, %f1, %f9
    \\ fmuls %f24, %f31, %f1
    \\ ld [%o0 + 2 * 4], %f24
    \\ fadds %f10, %f2, %f10
    \\ add %l3, -1, %l3
    \\ fmuls %f25, %f30, %f2
    \\ fadds %f11, %f3, %f11
    \\ cmp %l3, 0
    \\ fmuls %f25, %f31, %f3
    \\ ld [%o0 + 3 * 4], %f25
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f26, %f30, %f0
    \\ fsubs %f5, %f1, %f5
    \\ fmuls %f26, %f31, %f1
    \\ ld [%o1 + 2 * 4], %f26
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f27, %f30, %f2
    \\ ld [%o4 + 2 * 4], %f30
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f27, %f31, %f3
    \\ ld [%o1 + 3 * 4], %f27
    \\ prefetch [%o1 + 36 * 4], 1
    \\ fadds %f8, %f0, %f8
    \\ ld [%o4 + 3 * 4], %f31
    \\ fmuls %f20, %f28, %f0
    \\ fsubs %f9, %f1, %f9
    \\ fmuls %f20, %f29, %f1
    \\ ld [%o0 + 4 * 4], %f20
    \\ fadds %f10, %f2, %f10
    \\ fmuls %f21, %f28, %f2
    \\ fadds %f11, %f3, %f11
    \\ fmuls %f21, %f29, %f3
    \\ ld [%o0 + 5 * 4], %f21
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f22, %f28, %f0
    \\ fsubs %f5, %f1, %f5
    \\ fmuls %f22, %f29, %f1
    \\ ld [%o1 + 4 * 4], %f22
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f23, %f28, %f2
    \\ ld [%o4 + 4 * 4], %f28
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f23, %f29, %f3
    \\ ld [%o1 + 5 * 4], %f23
    \\ fadds %f8, %f0, %f8
    \\ ld [%o4 + 5 * 4], %f29
    \\ fmuls %f24, %f30, %f0
    \\ fsubs %f9, %f1, %f9
    \\ fmuls %f24, %f31, %f1
    \\ ld [%o0 + 6 * 4], %f24
    \\ fadds %f10, %f2, %f10
    \\ fmuls %f25, %f30, %f2
    \\ fadds %f11, %f3, %f11
    \\ fmuls %f25, %f31, %f3
    \\ ld [%o0 + 7 * 4], %f25
    \\ add %o0, 8 * 4, %o0
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f26, %f30, %f0
    \\ fsubs %f5, %f1, %f5
    \\ fmuls %f26, %f31, %f1
    \\ ld [%o1 + 6 * 4], %f26
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f27, %f30, %f2
    \\ ld [%o4 + 6 * 4], %f30
    \\ fadds %f7, %f3, %f7
    \\ add %o4, 8 * 4, %o4
    \\ fmuls %f27, %f31, %f3
    \\ ld [%o1 + 7 * 4], %f27
    \\ add %o1, 8 * 4, %o1
    \\ bg,pn %icc, 20012b
    \\ ld [%o4 - 1 * 4], %f31
    \\ 20013:
    \\ fadds %f8, %f0, %f8
    \\ fmuls %f20, %f28, %f0
    \\ fsubs %f9, %f1, %f9
    \\ fmuls %f20, %f29, %f1
    \\ ld [%o0 + 0 * 4], %f20
    \\ fadds %f10, %f2, %f10
    \\ fmuls %f21, %f28, %f2
    \\ fadds %f11, %f3, %f11
    \\ fmuls %f21, %f29, %f3
    \\ ld [%o0 + 1 * 4], %f21
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f22, %f28, %f0
    \\ fsubs %f5, %f1, %f5
    \\ fmuls %f22, %f29, %f1
    \\ ld [%o1 + 0 * 4], %f22
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f23, %f28, %f2
    \\ ld [%o4 + 0 * 4], %f28
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f23, %f29, %f3
    \\ ld [%o1 + 1 * 4], %f23
    \\ fadds %f8, %f0, %f8
    \\ ld [%o4 + 1 * 4], %f29
    \\ fmuls %f24, %f30, %f0
    \\ fsubs %f9, %f1, %f9
    \\ fmuls %f24, %f31, %f1
    \\ ld [%o0 + 2 * 4], %f24
    \\ fadds %f10, %f2, %f10
    \\ fmuls %f25, %f30, %f2
    \\ fadds %f11, %f3, %f11
    \\ fmuls %f25, %f31, %f3
    \\ ld [%o0 + 3 * 4], %f25
    \\ add %o0, 4 * 4, %o0
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f26, %f30, %f0
    \\ fsubs %f5, %f1, %f5
    \\ fmuls %f26, %f31, %f1
    \\ ld [%o1 + 2 * 4], %f26
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f27, %f30, %f2
    \\ ld [%o4 + 2 * 4], %f30
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f27, %f31, %f3
    \\ ld [%o1 + 3 * 4], %f27
    \\ add %o1, 4 * 4, %o1
    \\ fadds %f8, %f0, %f8
    \\ ld [%o4 + 3 * 4], %f31
    \\ add %o4, 4 * 4, %o4
    \\ fmuls %f20, %f28, %f0
    \\ fsubs %f9, %f1, %f9
    \\ fmuls %f20, %f29, %f1
    \\ fadds %f10, %f2, %f10
    \\ fmuls %f21, %f28, %f2
    \\ fadds %f11, %f3, %f11
    \\ fmuls %f21, %f29, %f3
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f22, %f28, %f0
    \\ fsubs %f5, %f1, %f5
    \\ fmuls %f22, %f29, %f1
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f23, %f28, %f2
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f23, %f29, %f3
    \\ fadds %f8, %f0, %f8
    \\ fmuls %f24, %f30, %f0
    \\ fsubs %f9, %f1, %f9
    \\ fmuls %f24, %f31, %f1
    \\ fadds %f10, %f2, %f10
    \\ fmuls %f25, %f30, %f2
    \\ fadds %f11, %f3, %f11
    \\ fmuls %f25, %f31, %f3
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f26, %f30, %f0
    \\ fsubs %f5, %f1, %f5
    \\ fmuls %f26, %f31, %f1
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f27, %f30, %f2
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f27, %f31, %f3
    \\
;
