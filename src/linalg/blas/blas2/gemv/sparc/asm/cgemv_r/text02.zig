//! Section 2 of 2 of the text of `cgemv_r`, which `__root__.zig` joins, from OpenBLAS `kernel/sparc/zgemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ 20006:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20007f
    \\ nop
    \\ ld [%o0 + 0 * 4], %f12
    \\ ld [%o0 + 1 * 4], %f13
    \\ ld [%o0 + 2 * 4], %f14
    \\ ld [%o0 + 3 * 4], %f15
    \\ ld [%l4 + 0 * 4], %f4
    \\ ld [%l4 + 1 * 4], %f5
    \\ ld [%l4 + 2 * 4], %f6
    \\ ld [%l4 + 3 * 4], %f7
    \\ fmuls %f12, %f28, %f0
    \\ ld [%o1 + 0 * 4], %f16
    \\ fmuls %f12, %f29, %f1
    \\ ld [%o1 + 1 * 4], %f17
    \\ fmuls %f14, %f28, %f2
    \\ ld [%o1 + 2 * 4], %f18
    \\ fmuls %f14, %f29, %f3
    \\ ld [%o1 + 3 * 4], %f19
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f13, %f29, %f0
    \\ fadds %f5, %f1, %f5
    \\ fmuls %f13, %f28, %f1
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f15, %f29, %f2
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f15, %f28, %f3
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f16, %f30, %f0
    \\ fsubs %f5, %f1, %f5
    \\ fmuls %f16, %f31, %f1
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f18, %f30, %f2
    \\ fsubs %f7, %f3, %f7
    \\ fmuls %f18, %f31, %f3
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f17, %f31, %f0
    \\ fadds %f5, %f1, %f5
    \\ fmuls %f17, %f30, %f1
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f19, %f31, %f2
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f19, %f30, %f3
    \\ fadds %f4, %f0, %f4
    \\ fsubs %f5, %f1, %f5
    \\ fadds %f6, %f2, %f6
    \\ fsubs %f7, %f3, %f7
    \\ st %f4, [%l4 + 0 * 4]
    \\ add %o0, 4 * 4, %o0
    \\ st %f5, [%l4 + 1 * 4]
    \\ add %o1, 4 * 4, %o1
    \\ st %f6, [%l4 + 2 * 4]
    \\ nop
    \\ st %f7, [%l4 + 3 * 4]
    \\ add %l4, 4 * 4, %l4
    \\ 20007:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20008f
    \\ nop
    \\ ld [%o0 + 0 * 4], %f12
    \\ ld [%o0 + 1 * 4], %f13
    \\ ld [%o1 + 0 * 4], %f14
    \\ ld [%o1 + 1 * 4], %f15
    \\ ld [%l4 + 0 * 4], %f4
    \\ ld [%l4 + 1 * 4], %f5
    \\ fmuls %f12, %f28, %f0
    \\ fmuls %f12, %f29, %f1
    \\ fmuls %f13, %f29, %f2
    \\ fmuls %f13, %f28, %f3
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f14, %f30, %f0
    \\ fadds %f5, %f1, %f5
    \\ fmuls %f14, %f31, %f1
    \\ fadds %f4, %f2, %f4
    \\ fmuls %f15, %f31, %f2
    \\ fsubs %f5, %f3, %f5
    \\ fmuls %f15, %f30, %f3
    \\ fadds %f4, %f0, %f4
    \\ fadds %f5, %f1, %f5
    \\ fadds %f4, %f2, %f4
    \\ fsubs %f5, %f3, %f5
    \\ st %f4, [%l4 + 0 * 4]
    \\ st %f5, [%l4 + 1 * 4]
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
    \\ ld [%sp + 128 + 16], %f26
    \\ ld [%sp + 128 + 20], %f27
    \\ .endif
    \\ .if \bp1
    \\ ld [%sp + 2423 + 32], %f26
    \\ ld [%sp + 2423 + 40], %f27
    \\ .endif
    \\ ld [%i3 + 0 * 4], %f28
    \\ ld [%i3 + 1 * 4], %f29
    \\ fmuls %f26, %f28, %f12
    \\ fmuls %f27, %f28, %f13
    \\ fmuls %f26, %f29, %f14
    \\ fmuls %f27, %f29, %f15
    \\ fsubs %f12, %f15, %f28
    \\ fadds %f13, %f14, %f29
    \\ sra %i0, 2, %l3
    \\ cmp %l3, 0
    \\ ble,pn %icc, 20013f
    \\ nop
    \\ ld [%o0 + 0 * 4], %f12
    \\ ld [%o0 + 1 * 4], %f13
    \\ ld [%o0 + 2 * 4], %f14
    \\ ld [%o0 + 3 * 4], %f15
    \\ ld [%o0 + 4 * 4], %f20
    \\ ld [%o0 + 5 * 4], %f21
    \\ ld [%o0 + 6 * 4], %f22
    \\ ld [%o0 + 7 * 4], %f23
    \\ ld [%l4 + 0 * 4], %f4
    \\ ld [%l4 + 1 * 4], %f5
    \\ ld [%l4 + 2 * 4], %f6
    \\ ld [%l4 + 3 * 4], %f7
    \\ ld [%l4 + 4 * 4], %f8
    \\ ld [%l4 + 5 * 4], %f9
    \\ ld [%l4 + 6 * 4], %f10
    \\ ld [%l4 + 7 * 4], %f11
    \\ fmuls %f12, %f28, %f0
    \\ deccc %l3
    \\ fmuls %f12, %f29, %f1
    \\ ld [%o0 + 8 * 4], %f12
    \\ fmuls %f14, %f28, %f2
    \\ fmuls %f14, %f29, %f3
    \\ ble,pn %icc, 20012f
    \\ ld [%o0 + 10 * 4], %f14
    \\ 20011:
    \\ fadds %f4, %f0, %f4
    \\ prefetch [%o0 + 88 * 4], 1
    \\ fmuls %f13, %f29, %f0
    \\ fadds %f5, %f1, %f5
    \\ fmuls %f13, %f28, %f1
    \\ ld [%o0 + 9 * 4], %f13
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f15, %f29, %f2
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f15, %f28, %f3
    \\ ld [%o0 + 11 * 4], %f15
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f20, %f28, %f0
    \\ fsubs %f5, %f1, %f5
    \\ fmuls %f20, %f29, %f1
    \\ ld [%o0 + 12 * 4], %f20
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f22, %f28, %f2
    \\ fsubs %f7, %f3, %f7
    \\ fmuls %f22, %f29, %f3
    \\ ld [%o0 + 14 * 4], %f22
    \\ st %f4, [%l4 + 0 * 4]
    \\ st %f5, [%l4 + 1 * 4]
    \\ st %f6, [%l4 + 2 * 4]
    \\ st %f7, [%l4 + 3 * 4]
    \\ fadds %f8, %f0, %f8
    \\ fmuls %f21, %f29, %f0
    \\ ld [%l4 + 8 * 4], %f4
    \\ fadds %f9, %f1, %f9
    \\ fmuls %f21, %f28, %f1
    \\ ld [%o0 + 13 * 4], %f21
    \\ fadds %f10, %f2, %f10
    \\ deccc %l3
    \\ fmuls %f23, %f29, %f2
    \\ ld [%l4 + 9 * 4], %f5
    \\ fadds %f11, %f3, %f11
    \\ fmuls %f23, %f28, %f3
    \\ ld [%o0 + 15 * 4], %f23
    \\ fadds %f8, %f0, %f8
    \\ add %o0, 8 * 4, %o0
    \\ fmuls %f12, %f28, %f0
    \\ ld [%l4 + 10 * 4], %f6
    \\ fsubs %f9, %f1, %f9
    \\ fmuls %f12, %f29, %f1
    \\ ld [%o0 + 8 * 4], %f12
    \\ fadds %f10, %f2, %f10
    \\ fmuls %f14, %f28, %f2
    \\ ld [%l4 + 11 * 4], %f7
    \\ fsubs %f11, %f3, %f11
    \\ fmuls %f14, %f29, %f3
    \\ ld [%o0 + 10 * 4], %f14
    \\ st %f8, [%l4 + 4 * 4]
    \\ st %f9, [%l4 + 5 * 4]
    \\ st %f10, [%l4 + 6 * 4]
    \\ st %f11, [%l4 + 7 * 4]
    \\ ld [%l4 + 12 * 4], %f8
    \\ ld [%l4 + 13 * 4], %f9
    \\ ld [%l4 + 14 * 4], %f10
    \\ add %l4, 8 * 4, %l4
    \\ bg,pn %icc, 20011b
    \\ ld [%l4 + 7 * 4], %f11
    \\ 20012:
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f13, %f29, %f0
    \\ fadds %f5, %f1, %f5
    \\ fmuls %f13, %f28, %f1
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f15, %f29, %f2
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f15, %f28, %f3
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f20, %f28, %f0
    \\ fsubs %f5, %f1, %f5
    \\ fmuls %f20, %f29, %f1
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f22, %f28, %f2
    \\ fsubs %f7, %f3, %f7
    \\ fmuls %f22, %f29, %f3
    \\ fadds %f8, %f0, %f8
    \\ fmuls %f21, %f29, %f0
    \\ fadds %f9, %f1, %f9
    \\ fmuls %f21, %f28, %f1
    \\ fadds %f10, %f2, %f10
    \\ fmuls %f23, %f29, %f2
    \\ fadds %f11, %f3, %f11
    \\ fmuls %f23, %f28, %f3
    \\ fadds %f8, %f0, %f8
    \\ fsubs %f9, %f1, %f9
    \\ fadds %f10, %f2, %f10
    \\ fsubs %f11, %f3, %f11
    \\ st %f4, [%l4 + 0 * 4]
    \\ st %f5, [%l4 + 1 * 4]
    \\ st %f6, [%l4 + 2 * 4]
    \\ st %f7, [%l4 + 3 * 4]
    \\ st %f8, [%l4 + 4 * 4]
    \\ st %f9, [%l4 + 5 * 4]
    \\ st %f10, [%l4 + 6 * 4]
    \\ st %f11, [%l4 + 7 * 4]
    \\ add %o0, 8 * 4, %o0
    \\ add %l4, 8 * 4, %l4
    \\ 20013:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20014f
    \\ nop
    \\ ld [%o0 + 0 * 4], %f12
    \\ ld [%o0 + 1 * 4], %f13
    \\ ld [%o0 + 2 * 4], %f14
    \\ ld [%o0 + 3 * 4], %f15
    \\ ld [%l4 + 0 * 4], %f4
    \\ fmuls %f12, %f28, %f0
    \\ ld [%l4 + 1 * 4], %f5
    \\ fmuls %f12, %f29, %f1
    \\ ld [%l4 + 2 * 4], %f6
    \\ fmuls %f14, %f28, %f2
    \\ ld [%l4 + 3 * 4], %f7
    \\ fmuls %f14, %f29, %f3
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f13, %f29, %f0
    \\ fadds %f5, %f1, %f5
    \\ fmuls %f13, %f28, %f1
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f15, %f29, %f2
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f15, %f28, %f3
    \\ fadds %f4, %f0, %f4
    \\ fsubs %f5, %f1, %f5
    \\ fadds %f6, %f2, %f6
    \\ fsubs %f7, %f3, %f7
    \\ st %f4, [%l4 + 0 * 4]
    \\ st %f5, [%l4 + 1 * 4]
    \\ st %f6, [%l4 + 2 * 4]
    \\ st %f7, [%l4 + 3 * 4]
    \\ add %o0, 4 * 4, %o0
    \\ add %l4, 4 * 4, %l4
    \\ 20014:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20015f
    \\ nop
    \\ ld [%o0 + 0 * 4], %f12
    \\ ld [%o0 + 1 * 4], %f13
    \\ ld [%l4 + 0 * 4], %f4
    \\ ld [%l4 + 1 * 4], %f5
    \\ fmuls %f12, %f28, %f0
    \\ fmuls %f12, %f29, %f1
    \\ fmuls %f13, %f29, %f2
    \\ fmuls %f13, %f28, %f3
    \\ fadds %f4, %f0, %f4
    \\ fadds %f5, %f1, %f5
    \\ fadds %f4, %f2, %f4
    \\ fsubs %f5, %f3, %f5
    \\ st %f4, [%l4 + 0 * 4]
    \\ st %f5, [%l4 + 1 * 4]
    \\ 20015:
    \\ cmp %l1, 2 * 4
    \\ be %icc, 20019f
    \\ mov %l0, %l4
    \\ sra %i0, 2, %l3
    \\ cmp %l3, 0
    \\ ble,pn %icc, 20017f
    \\ nop
    \\ 20016:
    \\ ld [%l2 + 0 * 4], %f12
    \\ ld [%l2 + 1 * 4], %f13
    \\ ld [%l0 + 0 * 4], %f4
    \\ ld [%l0 + 1 * 4], %f5
    \\ add %l0, %l1, %l0
    \\ ld [%l2 + 2 * 4], %f14
    \\ ld [%l2 + 3 * 4], %f15
    \\ ld [%l0 + 0 * 4], %f6
    \\ ld [%l0 + 1 * 4], %f7
    \\ add %l0, %l1, %l0
    \\ ld [%l2 + 4 * 4], %f16
    \\ ld [%l2 + 5 * 4], %f17
    \\ ld [%l0 + 0 * 4], %f8
    \\ ld [%l0 + 1 * 4], %f9
    \\ add %l0, %l1, %l0
    \\ ld [%l2 + 6 * 4], %f18
    \\ ld [%l2 + 7 * 4], %f19
    \\ ld [%l0 + 0 * 4], %f10
    \\ ld [%l0 + 1 * 4], %f11
    \\ add %l0, %l1, %l0
    \\ fadds %f4, %f12, %f4
    \\ fadds %f5, %f13, %f5
    \\ fadds %f6, %f14, %f6
    \\ fadds %f7, %f15, %f7
    \\ fadds %f8, %f16, %f8
    \\ fadds %f9, %f17, %f9
    \\ fadds %f10, %f18, %f10
    \\ fadds %f11, %f19, %f11
    \\ st %f4, [%l4 + 0 * 4]
    \\ st %f5, [%l4 + 1 * 4]
    \\ add %l4, %l1, %l4
    \\ st %f6, [%l4 + 0 * 4]
    \\ st %f7, [%l4 + 1 * 4]
    \\ add %l4, %l1, %l4
    \\ st %f8, [%l4 + 0 * 4]
    \\ st %f9, [%l4 + 1 * 4]
    \\ add %l4, %l1, %l4
    \\ st %f10, [%l4 + 0 * 4]
    \\ st %f11, [%l4 + 1 * 4]
    \\ add %l4, %l1, %l4
    \\ deccc %l3
    \\ bg,pn %icc, 20016b
    \\ add %l2, 8 * 4, %l2
    \\ 20017:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20018f
    \\ nop
    \\ ld [%l2 + 0 * 4], %f12
    \\ ld [%l2 + 1 * 4], %f13
    \\ ld [%l0 + 0 * 4], %f4
    \\ ld [%l0 + 1 * 4], %f5
    \\ add %l0, %l1, %l0
    \\ ld [%l2 + 2 * 4], %f14
    \\ ld [%l2 + 3 * 4], %f15
    \\ ld [%l0 + 0 * 4], %f6
    \\ ld [%l0 + 1 * 4], %f7
    \\ add %l0, %l1, %l0
    \\ fadds %f4, %f12, %f4
    \\ fadds %f5, %f13, %f5
    \\ fadds %f6, %f14, %f6
    \\ fadds %f7, %f15, %f7
    \\ st %f4, [%l4 + 0 * 4]
    \\ st %f5, [%l4 + 1 * 4]
    \\ add %l4, %l1, %l4
    \\ st %f6, [%l4 + 0 * 4]
    \\ st %f7, [%l4 + 1 * 4]
    \\ add %l4, %l1, %l4
    \\ add %l2, 4 * 4, %l2
    \\ 20018:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20019f
    \\ nop
    \\ ld [%l2 + 0 * 4], %f12
    \\ ld [%l2 + 1 * 4], %f13
    \\ ld [%l0 + 0 * 4], %f4
    \\ ld [%l0 + 1 * 4], %f5
    \\ fadds %f4, %f12, %f4
    \\ fadds %f5, %f13, %f5
    \\ st %f4, [%l4 + 0 * 4]
    \\ st %f5, [%l4 + 1 * 4]
    \\ 20019:
    \\ return %i7 + 8
    \\ clr %o0
    \\ .endm
    \\
;
