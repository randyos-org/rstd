//! Section 1 of 4 of the text of `cgemv_d`, which `__root__.zig` joins, from OpenBLAS `kernel/sparc/zgemv_t.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1
    \\ .align 32
    \\ .proc 07
    \\ .if \bp0
    \\ save %sp, -64, %sp
    \\ .endif
    \\ .if \bp1
    \\ save %sp, -256, %sp
    \\ .endif
    \\ nop
    \\ .if \bp0
    \\ st %i3, [%sp + 128 + 16]
    \\ st %i4, [%sp + 128 + 20]
    \\ ld [%sp + 128 + 28], %i2
    \\ ld [%sp + 128 + 32], %i3
    \\ ld [%sp + 128 + 36], %i4
    \\ ld [%sp + 128 + 40], %l0
    \\ ld [%sp + 128 + 44], %l1
    \\ ld [%sp + 128 + 48], %l2
    \\ .endif
    \\ .if \bp1
    \\ ldx [%sp + 2423 + 56], %i2
    \\ ldx [%sp + 2423 + 64], %i3
    \\ ldx [%sp + 2423 + 72], %i4
    \\ ldx [%sp + 2423 + 80], %l0
    \\ ldx [%sp + 2423 + 88], %l1
    \\ ldx [%sp + 2423 + 96], %l2
    \\ st %f7, [%sp + 2423 + 32]
    \\ st %f9, [%sp + 2423 + 40]
    \\ .endif
    \\ clr %l4
    \\ mov 4000, %l3
    \\ .if \bp0
    \\ sll %i2, 3, %i2
    \\ sll %l3, 3, %l3
    \\ .endif
    \\ .if \bp1
    \\ sllx %i2, 3, %i2
    \\ sllx %l3, 3, %l3
    \\ .endif
    \\ smul %i2, %i1, %g1
    \\ .if \bp0
    \\ sll %i4, 3, %i4
    \\ sll %l1, 3, %l1
    \\ .endif
    \\ .if \bp1
    \\ sllx %i4, 3, %i4
    \\ sllx %l1, 3, %l1
    \\ .endif
    \\ sub %l3, %g1, %g1
    \\ 20000:
    \\ .if \bp0
    \\ sll %l4, 3, %l3
    \\ .endif
    \\ .if \bp1
    \\ sllx %l4, 3, %l3
    \\ .endif
    \\ sub %i0, %l4, %l6
    \\ mov 4000, %l5
    \\ cmp %l6, %l5
    \\ nop
    \\ movg %icc, %l5, %l6
    \\ nop
    \\ cmp %i4, 2 * 4
    \\ beq 20004f
    \\ add %i3, %l3, %l7
    \\ sra %l6, 2, %l3
    \\ mov %l2, %l7
    \\ cmp %l3, 0
    \\ ble,pn %icc, 20002f
    \\ mov %l2, %o5
    \\ 20001:
    \\ ld [%i3 + 0 * 4], %f20
    \\ ld [%i3 + 1 * 4], %f21
    \\ add %i3, %i4, %i3
    \\ ld [%i3 + 0 * 4], %f22
    \\ ld [%i3 + 1 * 4], %f23
    \\ add %i3, %i4, %i3
    \\ ld [%i3 + 0 * 4], %f24
    \\ ld [%i3 + 1 * 4], %f25
    \\ add %i3, %i4, %i3
    \\ ld [%i3 + 0 * 4], %f26
    \\ ld [%i3 + 1 * 4], %f27
    \\ add %i3, %i4, %i3
    \\ st %f20, [%o5 + 0 * 4]
    \\ add %l3, -1, %l3
    \\ st %f21, [%o5 + 1 * 4]
    \\ cmp %l3, 0
    \\ st %f22, [%o5 + 2 * 4]
    \\ st %f23, [%o5 + 3 * 4]
    \\ st %f24, [%o5 + 4 * 4]
    \\ st %f25, [%o5 + 5 * 4]
    \\ st %f26, [%o5 + 6 * 4]
    \\ st %f27, [%o5 + 7 * 4]
    \\ bg,pn %icc, 20001b
    \\ add %o5, 8 * 4, %o5
    \\ 20002:
    \\ and %l6, 3, %l3
    \\ cmp %l3, 0
    \\ ble,pn %icc, 20004f
    \\ nop
    \\ 20003:
    \\ ld [%i3 + 0 * 4], %f20
    \\ ld [%i3 + 1 * 4], %f21
    \\ add %i3, %i4, %i3
    \\ add %l3, -1, %l3
    \\ cmp %l3, 0
    \\ nop
    \\ st %f20, [%o5 + 0 * 4]
    \\ st %f21, [%o5 + 1 * 4]
    \\ bg,pn %icc, 20003b
    \\ add %o5, 2 * 4, %o5
    \\ 20004:
    \\ sra %i1, 2, %l5
    \\ cmp %l5, 0
    \\ ble %icc, 20011f
    \\ mov %l0, %o5
    \\ 20005:
    \\ .word ((2 << 30) | ((0) << 25) | ( 0x36 << 19) | ( 0x61 << 5))
    \\ fmovs %f0, %f4
    \\ sra %l6, 2, %l3
    \\ fmovs %f0, %f5
    \\ add %i5, %i2, %o1
    \\ fmovs %f0, %f6
    \\ mov %i5, %o0
    \\ fmovs %f0, %f7
    \\ add %o1, %i2, %o2
    \\ fmovs %f0, %f8
    \\ fmovs %f0, %f9
    \\ fmovs %f0, %f10
    \\ fmovs %f0, %f11
    \\ fmovs %f0, %f12
    \\ fmovs %f0, %f13
    \\ fmovs %f0, %f14
    \\ fmovs %f0, %f15
    \\ fmovs %f0, %f16
    \\ fmovs %f0, %f17
    \\ fmovs %f0, %f18
    \\ fmovs %f0, %f19
    \\ add %o2, %i2, %o3
    \\ fmovs %f0, %f1
    \\ mov %l7, %o4
    \\ fmovs %f0, %f2
    \\ add %o3, %i2, %i5
    \\ cmp %l3, 0
    \\ ble %icc, 20008f
    \\ fmovs %f0, %f3
    \\ ld [%o0 + 0 * 4], %f20
    \\ nop
    \\ ld [%o0 + 1 * 4], %f21
    \\ add %o0, 2 * 4, %o0
    \\ ld [%o1 + 0 * 4], %f22
    \\ ld [%o1 + 1 * 4], %f23
    \\ add %o1, 2 * 4, %o1
    \\ ld [%o2 + 0 * 4], %f24
    \\ ld [%o2 + 1 * 4], %f25
    \\ add %o2, 2 * 4, %o2
    \\ ld [%o3 + 0 * 4], %f26
    \\ ld [%o3 + 1 * 4], %f27
    \\ add %o3, 2 * 4, %o3
    \\ ld [%o4 + 0 * 4], %f28
    \\ nop
    \\ ld [%o4 + 1 * 4], %f29
    \\ nop
    \\ ld [%o4 + 2 * 4], %f30
    \\ add %o4, 4 * 4, %o4
    \\ deccc %l3
    \\ ble 20007f
    \\ prefetch [%o5 + 7 * 4], 2
    \\ 20006:
    \\ fadds %f16, %f0, %f16
    \\ prefetch [%o0 + 36 * 4], 1
    \\ fmuls %f20, %f28, %f0
    \\ nop
    \\ fsubs %f17, %f1, %f17
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
    \\ fsubs %f5, %f1, %f5
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
    \\ fsubs %f9, %f1, %f9
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
    \\ fsubs %f13, %f1, %f13
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
    \\ prefetch [%o1 + 36 * 4], 1
    \\ fsubs %f17, %f1, %f17
    \\ nop
    \\ fmuls %f20, %f31, %f1
    \\ ld [%o0 + 2 * 4], %f20
    \\ fadds %f18, %f2, %f18
    \\ nop
    \\ fmuls %f21, %f30, %f2
    \\ ld [%o4 + 1 * 4], %f29
    \\ fadds %f19, %f3, %f19
    \\ nop
    \\ fmuls %f21, %f31, %f3
    \\ ld [%o0 + 3 * 4], %f21
    \\ fadds %f4, %f0, %f4
    \\ nop
    \\ fmuls %f22, %f30, %f0
    \\ nop
    \\ fsubs %f5, %f1, %f5
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
    \\ fsubs %f9, %f1, %f9
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
    \\ fsubs %f13, %f1, %f13
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
    \\ prefetch [%o2 + 36 * 4], 1
    \\ fmuls %f20, %f28, %f0
    \\ nop
    \\ fsubs %f17, %f1, %f17
    \\ nop
    \\ fmuls %f20, %f29, %f1
    \\ ld [%o0 + 4 * 4], %f20
    \\ fadds %f18, %f2, %f18
    \\ nop
    \\ fmuls %f21, %f28, %f2
    \\ ld [%o4 + 3 * 4], %f31
    \\ fadds %f19, %f3, %f19
    \\ nop
    \\ fmuls %f21, %f29, %f3
    \\ ld [%o0 + 5 * 4], %f21
    \\ fadds %f4, %f0, %f4
    \\ nop
    \\ fmuls %f22, %f28, %f0
    \\ nop
    \\ fsubs %f5, %f1, %f5
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
    \\ nop
    \\ fmuls %f24, %f28, %f0
    \\ nop
    \\ fsubs %f9, %f1, %f9
    \\ nop
    \\ fmuls %f24, %f29, %f1
    \\ ld [%o2 + 4 * 4], %f24
    \\ fadds %f10, %f2, %f10
    \\ deccc %l3
    \\ fmuls %f25, %f28, %f2
    \\ nop
    \\ fadds %f11, %f3, %f11
    \\ nop
    \\ fmuls %f25, %f29, %f3
    \\ ld [%o2 + 5 * 4], %f25
    \\ fadds %f12, %f0, %f12
    \\ nop
    \\ fmuls %f26, %f28, %f0
    \\ nop
    \\ fsubs %f13, %f1, %f13
    \\ nop
    \\ fmuls %f26, %f29, %f1
    \\ ld [%o3 + 4 * 4], %f26
    \\ fadds %f14, %f2, %f14
    \\ nop
    \\ fmuls %f27, %f28, %f2
    \\ ld [%o4 + 4 * 4], %f28
    \\ fadds %f15, %f3, %f15
    \\ nop
    \\ fmuls %f27, %f29, %f3
    \\ ld [%o3 + 5 * 4], %f27
    \\ fadds %f16, %f0, %f16
    \\ prefetch [%o3 + 36 * 4], 1
    \\ fmuls %f20, %f30, %f0
    \\ nop
    \\ fsubs %f17, %f1, %f17
    \\ nop
    \\ fmuls %f20, %f31, %f1
    \\ ld [%o0 + 6 * 4], %f20
    \\ fadds %f18, %f2, %f18
    \\ nop
    \\ fmuls %f21, %f30, %f2
    \\ ld [%o4 + 5 * 4], %f29
    \\ fadds %f19, %f3, %f19
    \\ nop
    \\ fmuls %f21, %f31, %f3
    \\ ld [%o0 + 7 * 4], %f21
    \\ fadds %f4, %f0, %f4
    \\ add %o0, 8 * 4, %o0
    \\ fmuls %f22, %f30, %f0
    \\ nop
    \\ fsubs %f5, %f1, %f5
    \\ nop
    \\ fmuls %f22, %f31, %f1
    \\ ld [%o1 + 6 * 4], %f22
    \\ fadds %f6, %f2, %f6
    \\ nop
    \\ fmuls %f23, %f30, %f2
    \\ nop
    \\ fadds %f7, %f3, %f7
    \\ nop
    \\ fmuls %f23, %f31, %f3
    \\ ld [%o1 + 7 * 4], %f23
    \\ fadds %f8, %f0, %f8
    \\ add %o1, 8 * 4, %o1
    \\ fmuls %f24, %f30, %f0
    \\ nop
    \\ fsubs %f9, %f1, %f9
    \\ nop
    \\ fmuls %f24, %f31, %f1
    \\ ld [%o2 + 6 * 4], %f24
    \\ fadds %f10, %f2, %f10
    \\ add %o3, 8 * 4, %o3
    \\ fmuls %f25, %f30, %f2
    \\ nop
    \\ fadds %f11, %f3, %f11
    \\ nop
    \\ fmuls %f25, %f31, %f3
    \\ ld [%o2 + 7 * 4], %f25
    \\ fadds %f12, %f0, %f12
    \\ add %o2, 8 * 4, %o2
    \\ fmuls %f26, %f30, %f0
    \\ nop
    \\ fsubs %f13, %f1, %f13
    \\ add %o4, 8 * 4, %o4
    \\ fmuls %f26, %f31, %f1
    \\ ld [%o3 - 2 * 4], %f26
    \\ fadds %f14, %f2, %f14
    \\ nop
    \\ fmuls %f27, %f30, %f2
    \\ ld [%o4 - 2 * 4], %f30
    \\ fadds %f15, %f3, %f15
    \\ fmuls %f27, %f31, %f3
    \\ bg,pn %icc, 20006b
    \\ ld [%o3 - 1 * 4], %f27
    \\
;
