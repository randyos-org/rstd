//! Section 3 of 3 of the text of `sgemv_n`, which `__root__.zig` joins, from OpenBLAS `kernel/sparc/gemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ 20014:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20015f
    \\ nop
    \\ ld [%o0 + 0 * 4], %f12
    \\ ld [%o1 + 0 * 4], %f13
    \\ ld [%l4 + 0 * 4], %f4
    \\ ld [%o0 + 1 * 4], %f16
    \\ ld [%o1 + 1 * 4], %f17
    \\ add %o0, 2 * 4, %o0
    \\ ld [%l4 + 1 * 4], %f5
    \\ add %o1, 2 * 4, %o1
    \\ fmuls %f12, %f28, %f0
    \\ fmuls %f13, %f29, %f1
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f16, %f28, %f0
    \\ fadds %f4, %f1, %f4
    \\ fmuls %f17, %f29, %f1
    \\ fadds %f5, %f0, %f5
    \\ fadds %f5, %f1, %f5
    \\ st %f4, [%l4 + 0 * 4]
    \\ st %f5, [%l4 + 1 * 4]
    \\ add %l4, 2 * 4, %l4
    \\ 20015:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20016f
    \\ nop
    \\ ld [%o0 + 0 * 4], %f12
    \\ ld [%o1 + 0 * 4], %f13
    \\ ld [%l4 + 0 * 4], %f4
    \\ fmuls %f12, %f28, %f0
    \\ fmuls %f13, %f29, %f1
    \\ fadds %f4, %f0, %f4
    \\ fadds %f4, %f1, %f4
    \\ st %f4, [%l4]
    \\ 20016:
    \\ andcc %i1, 1, %l5
    \\ ble,pn %icc, 20023f
    \\ nop
    \\ 20017:
    \\ mov %l6, %l4
    \\ mov %i4, %o0
    \\ .if \bp0
    \\ ld [%sp + 128 + 16], %f27
    \\ .endif
    \\ .if \bp1
    \\ ld [%sp + 2423 + 40], %f27
    \\ .endif
    \\ ld [%i2], %f28
    \\ add %i2, %i3, %i2
    \\ fmuls %f27, %f28, %f28
    \\ sra %i0, 3, %l3
    \\ cmp %l3, 0
    \\ ble,pn %icc, 20020f
    \\ nop
    \\ ld [%l4 + 0 * 4], %f4
    \\ ld [%l4 + 1 * 4], %f5
    \\ ld [%l4 + 2 * 4], %f6
    \\ ld [%l4 + 3 * 4], %f7
    \\ ld [%l4 + 4 * 4], %f8
    \\ ld [%l4 + 5 * 4], %f9
    \\ ld [%l4 + 6 * 4], %f10
    \\ ld [%l4 + 7 * 4], %f11
    \\ ld [%o0 + 0 * 4], %f12
    \\ ld [%o0 + 1 * 4], %f13
    \\ ld [%o0 + 2 * 4], %f14
    \\ ld [%o0 + 3 * 4], %f15
    \\ ld [%o0 + 4 * 4], %f16
    \\ ld [%o0 + 5 * 4], %f17
    \\ ld [%o0 + 6 * 4], %f18
    \\ ld [%o0 + 7 * 4], %f19
    \\ fmuls %f12, %f28, %f0
    \\ deccc %l3
    \\ ld [%o0 + 8 * 4], %f12
    \\ fmuls %f13, %f28, %f1
    \\ ld [%o0 + 9 * 4], %f13
    \\ fmuls %f14, %f28, %f2
    \\ ld [%o0 + 10 * 4], %f14
    \\ fmuls %f15, %f28, %f3
    \\ ble,pn %icc, 20019f
    \\ ld [%o0 + 11 * 4], %f15
    \\ 20018:
    \\ fadds %f4, %f0, %f4
    \\ prefetch [%o0 + 40 * 4], 1
    \\ fmuls %f16, %f28, %f0
    \\ ld [%o0 + 12 * 4], %f16
    \\ fadds %f5, %f1, %f5
    \\ fmuls %f17, %f28, %f1
    \\ ld [%o0 + 13 * 4], %f17
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f18, %f28, %f2
    \\ ld [%o0 + 14 * 4], %f18
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f19, %f28, %f3
    \\ ld [%o0 + 15 * 4], %f19
    \\ fadds %f8, %f0, %f8
    \\ fmuls %f12, %f28, %f0
    \\ ld [%o0 + 16 * 4], %f12
    \\ fadds %f9, %f1, %f9
    \\ fmuls %f13, %f28, %f1
    \\ ld [%o0 + 17 * 4], %f13
    \\ fadds %f10, %f2, %f10
    \\ fmuls %f14, %f28, %f2
    \\ ld [%o0 + 18 * 4], %f14
    \\ fadds %f11, %f3, %f11
    \\ fmuls %f15, %f28, %f3
    \\ ld [%o0 + 19 * 4], %f15
    \\ st %f4, [%l4 + 0 * 4]
    \\ st %f5, [%l4 + 1 * 4]
    \\ st %f6, [%l4 + 2 * 4]
    \\ st %f7, [%l4 + 3 * 4]
    \\ st %f8, [%l4 + 4 * 4]
    \\ st %f9, [%l4 + 5 * 4]
    \\ st %f10, [%l4 + 6 * 4]
    \\ st %f11, [%l4 + 7 * 4]
    \\ ld [%l4 + 8 * 4], %f4
    \\ ld [%l4 + 9 * 4], %f5
    \\ ld [%l4 + 10 * 4], %f6
    \\ ld [%l4 + 11 * 4], %f7
    \\ ld [%l4 + 12 * 4], %f8
    \\ deccc %l3
    \\ ld [%l4 + 13 * 4], %f9
    \\ add %o0, 8 * 4, %o0
    \\ ld [%l4 + 14 * 4], %f10
    \\ add %l4, 8 * 4, %l4
    \\ bg,pn %icc, 20018b
    \\ ld [%l4 + 7 * 4], %f11
    \\ 20019:
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f16, %f28, %f0
    \\ fadds %f5, %f1, %f5
    \\ fmuls %f17, %f28, %f1
    \\ fadds %f6, %f2, %f6
    \\ fmuls %f18, %f28, %f2
    \\ fadds %f7, %f3, %f7
    \\ fmuls %f19, %f28, %f3
    \\ st %f4, [%l4 + 0 * 4]
    \\ fadds %f8, %f0, %f8
    \\ st %f5, [%l4 + 1 * 4]
    \\ fadds %f9, %f1, %f9
    \\ st %f6, [%l4 + 2 * 4]
    \\ fadds %f10, %f2, %f10
    \\ st %f7, [%l4 + 3 * 4]
    \\ fadds %f11, %f3, %f11
    \\ st %f8, [%l4 + 4 * 4]
    \\ st %f9, [%l4 + 5 * 4]
    \\ st %f10, [%l4 + 6 * 4]
    \\ add %o0, 8 * 4, %o0
    \\ st %f11, [%l4 + 7 * 4]
    \\ add %l4, 8 * 4, %l4
    \\ 20020:
    \\ andcc %i0, 4, %l3
    \\ ble,pn %icc, 20021f
    \\ nop
    \\ ld [%o0 + 0 * 4], %f12
    \\ ld [%o0 + 1 * 4], %f13
    \\ ld [%o0 + 2 * 4], %f14
    \\ ld [%o0 + 3 * 4], %f15
    \\ ld [%l4 + 0 * 4], %f4
    \\ add %o0, 4 * 4, %o0
    \\ ld [%l4 + 1 * 4], %f5
    \\ ld [%l4 + 2 * 4], %f6
    \\ ld [%l4 + 3 * 4], %f7
    \\ fmuls %f12, %f28, %f0
    \\ fmuls %f13, %f28, %f1
    \\ fmuls %f14, %f28, %f2
    \\ fmuls %f15, %f28, %f3
    \\ fadds %f4, %f0, %f4
    \\ fadds %f5, %f1, %f5
    \\ fadds %f6, %f2, %f6
    \\ fadds %f7, %f3, %f7
    \\ st %f4, [%l4 + 0 * 4]
    \\ st %f5, [%l4 + 1 * 4]
    \\ st %f6, [%l4 + 2 * 4]
    \\ st %f7, [%l4 + 3 * 4]
    \\ add %l4, 4 * 4, %l4
    \\ 20021:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20022f
    \\ nop
    \\ ld [%o0 + 0 * 4], %f12
    \\ ld [%l4 + 0 * 4], %f4
    \\ ld [%o0 + 1 * 4], %f16
    \\ ld [%l4 + 1 * 4], %f5
    \\ add %o0, 2 * 4, %o0
    \\ fmuls %f12, %f28, %f0
    \\ fadds %f4, %f0, %f4
    \\ fmuls %f16, %f28, %f0
    \\ fadds %f5, %f0, %f5
    \\ st %f4, [%l4 + 0 * 4]
    \\ st %f5, [%l4 + 1 * 4]
    \\ add %l4, 2 * 4, %l4
    \\ 20022:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20023f
    \\ nop
    \\ ld [%o0 + 0 * 4], %f12
    \\ ld [%l4 + 0 * 4], %f4
    \\ fmuls %f12, %f28, %f0
    \\ fadds %f4, %f0, %f4
    \\ st %f4, [%l4]
    \\ 20023:
    \\ cmp %l1, 4
    \\ be %icc, 20028f
    \\ mov %l0, %l4
    \\ sra %i0, 3, %l3
    \\ cmp %l3, 0
    \\ ble,pn %icc, 20025f
    \\ nop
    \\ 20024:
    \\ ld [%l2 + 0 * 4], %f12
    \\ ld [%l0], %f4
    \\ add %l0, %l1, %l0
    \\ ld [%l2 + 1 * 4], %f13
    \\ ld [%l0], %f5
    \\ add %l0, %l1, %l0
    \\ ld [%l2 + 2 * 4], %f14
    \\ ld [%l0], %f6
    \\ add %l0, %l1, %l0
    \\ ld [%l2 + 3 * 4], %f15
    \\ ld [%l0], %f7
    \\ add %l0, %l1, %l0
    \\ ld [%l2 + 4 * 4], %f16
    \\ fadds %f4, %f12, %f4
    \\ ld [%l0], %f8
    \\ add %l0, %l1, %l0
    \\ ld [%l2 + 5 * 4], %f17
    \\ fadds %f5, %f13, %f5
    \\ ld [%l0], %f9
    \\ add %l0, %l1, %l0
    \\ ld [%l2 + 6 * 4], %f18
    \\ fadds %f6, %f14, %f6
    \\ ld [%l0], %f10
    \\ add %l0, %l1, %l0
    \\ ld [%l2 + 7 * 4], %f19
    \\ fadds %f7, %f15, %f7
    \\ ld [%l0], %f11
    \\ add %l0, %l1, %l0
    \\ st %f4, [%l4]
    \\ fadds %f8, %f16, %f8
    \\ add %l4, %l1, %l4
    \\ st %f5, [%l4]
    \\ fadds %f9, %f17, %f9
    \\ add %l4, %l1, %l4
    \\ st %f6, [%l4]
    \\ fadds %f10, %f18, %f10
    \\ add %l4, %l1, %l4
    \\ st %f7, [%l4]
    \\ fadds %f11, %f19, %f11
    \\ add %l4, %l1, %l4
    \\ st %f8, [%l4]
    \\ add %l4, %l1, %l4
    \\ st %f9, [%l4]
    \\ add %l4, %l1, %l4
    \\ st %f10, [%l4]
    \\ add %l4, %l1, %l4
    \\ st %f11, [%l4]
    \\ add %l4, %l1, %l4
    \\ deccc %l3
    \\ bg,pn %icc, 20024b
    \\ add %l2, 8 * 4, %l2
    \\ 20025:
    \\ andcc %i0, 7, %l3
    \\ ble,pn %icc, 20028f
    \\ nop
    \\ andcc %i0, 4, %l3
    \\ ble,pn %icc, 20026f
    \\ nop
    \\ ld [%l2 + 0 * 4], %f12
    \\ ld [%l2 + 1 * 4], %f13
    \\ ld [%l2 + 2 * 4], %f14
    \\ ld [%l2 + 3 * 4], %f15
    \\ add %l2, 4 * 4, %l2
    \\ ld [%l0], %f4
    \\ add %l0, %l1, %l0
    \\ ld [%l0], %f5
    \\ add %l0, %l1, %l0
    \\ ld [%l0], %f6
    \\ add %l0, %l1, %l0
    \\ ld [%l0], %f7
    \\ add %l0, %l1, %l0
    \\ fadds %f4, %f12, %f4
    \\ fadds %f5, %f13, %f5
    \\ fadds %f6, %f14, %f6
    \\ fadds %f7, %f15, %f7
    \\ st %f4, [%l4]
    \\ add %l4, %l1, %l4
    \\ st %f5, [%l4]
    \\ add %l4, %l1, %l4
    \\ st %f6, [%l4]
    \\ add %l4, %l1, %l4
    \\ st %f7, [%l4]
    \\ add %l4, %l1, %l4
    \\ 20026:
    \\ andcc %i0, 2, %l3
    \\ ble,pn %icc, 20027f
    \\ nop
    \\ ld [%l2 + 0 * 4], %f12
    \\ ld [%l2 + 1 * 4], %f13
    \\ add %l2, 2 * 4, %l2
    \\ ld [%l0], %f4
    \\ add %l0, %l1, %l0
    \\ ld [%l0], %f5
    \\ add %l0, %l1, %l0
    \\ fadds %f4, %f12, %f4
    \\ fadds %f5, %f13, %f5
    \\ st %f4, [%l4]
    \\ add %l4, %l1, %l4
    \\ st %f5, [%l4]
    \\ add %l4, %l1, %l4
    \\ 20027:
    \\ andcc %i0, 1, %l3
    \\ ble,pn %icc, 20028f
    \\ nop
    \\ ld [%l2 + 0 * 4], %f12
    \\ ld [%l0], %f4
    \\ fadds %f4, %f12, %f4
    \\ st %f4, [%l4]
    \\ 20028:
    \\ return %i7 + 8
    \\ clr %o0
    \\ .endm
    \\
;
