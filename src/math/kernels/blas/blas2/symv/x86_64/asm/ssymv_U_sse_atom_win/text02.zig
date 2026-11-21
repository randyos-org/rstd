//! Section 2 of 3 of the text of `ssymv_U`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/symv_U_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20009:
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm12, %xmm11
    \\ mulps %xmm0, %xmm12
    \\ addps %xmm11, %xmm4
    \\ addps %xmm12, %xmm10
    \\ .if \bp0
    \\ movsd 4 * 4(%rbx), %xmm12
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rbx), %xmm12
    \\ .endif
    \\ movhps 6 * 4(%rbx), %xmm12
    \\ .if \bp2
    \\ prefetcht0 (16 * 12)(%rbx)
    \\ .endif
    \\ .if \bp3
    \\ prefetch (16 * 16)(%rbx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 (16 * 24)(%rbx)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 8)(%rbx)
    \\ .endif
    \\ .if \bp5
    \\ prefetcht0 (16 * 20)(%rbx)
    \\ .endif
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm13, %xmm11
    \\ mulps %xmm1, %xmm13
    \\ addps %xmm11, %xmm5
    \\ addps %xmm13, %xmm10
    \\ .if \bp0
    \\ movsd 4 * 4(%rbx, %r8, 1), %xmm13
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rbx, %r8, 1), %xmm13
    \\ .endif
    \\ movhps 6 * 4(%rbx, %r8, 1), %xmm13
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm14, %xmm11
    \\ mulps %xmm2, %xmm14
    \\ addps %xmm11, %xmm6
    \\ addps %xmm14, %xmm10
    \\ .if \bp0
    \\ movsd 4 * 4(%rbp), %xmm14
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rbp), %xmm14
    \\ .endif
    \\ movhps 6 * 4(%rbp), %xmm14
    \\ .if \bp6
    \\ prefetcht0 (16 * 12)(%r13)
    \\ .endif
    \\ .if \bp3
    \\ prefetch (16 * 16)(%r13)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 (16 * 24)(%r13)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 8)(%r13)
    \\ .endif
    \\ .if \bp5
    \\ prefetcht0 (16 * 20)(%r13)
    \\ .endif
    \\ movaps %xmm8, %xmm11
    \\ movaps 8 * 4(%r13), %xmm8
    \\ mulps %xmm15, %xmm11
    \\ mulps %xmm3, %xmm15
    \\ addps %xmm11, %xmm7
    \\ addps %xmm15, %xmm10
    \\ .if \bp0
    \\ movsd 4 * 4(%rbp, %r8, 1), %xmm15
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rbp, %r8, 1), %xmm15
    \\ .endif
    \\ movhps 6 * 4(%rbp, %r8, 1), %xmm15
    \\ movlps %xmm10, 0 * 4(%r14)
    \\ movhps %xmm10, 2 * 4(%r14)
    \\ .if \bp0
    \\ movsd 4 * 4(%r14), %xmm10
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%r14), %xmm10
    \\ .endif
    \\ movhps 6 * 4(%r14), %xmm10
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm12, %xmm11
    \\ mulps %xmm0, %xmm12
    \\ addps %xmm11, %xmm4
    \\ addps %xmm12, %xmm10
    \\ .if \bp0
    \\ movsd 8 * 4(%rbx), %xmm12
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rbx), %xmm12
    \\ .endif
    \\ movhps 10 * 4(%rbx), %xmm12
    \\ .if \bp2
    \\ prefetcht0 (16 * 12)(%rbx, %r8, 1)
    \\ .endif
    \\ .if \bp3
    \\ prefetch (16 * 16)(%rbx, %r8, 1)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 (16 * 24)(%rbx, %r8, 1)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 8)(%rbx, %r8, 1)
    \\ .endif
    \\ .if \bp5
    \\ prefetcht0 (16 * 20)(%rbx, %r8, 1)
    \\ .endif
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm13, %xmm11
    \\ mulps %xmm1, %xmm13
    \\ addps %xmm11, %xmm5
    \\ addps %xmm13, %xmm10
    \\ .if \bp0
    \\ movsd 8 * 4(%rbx, %r8, 1), %xmm13
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rbx, %r8, 1), %xmm13
    \\ .endif
    \\ movhps 10 * 4(%rbx, %r8, 1), %xmm13
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm14, %xmm11
    \\ mulps %xmm2, %xmm14
    \\ addps %xmm11, %xmm6
    \\ addps %xmm14, %xmm10
    \\ .if \bp0
    \\ movsd 8 * 4(%rbp), %xmm14
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rbp), %xmm14
    \\ .endif
    \\ movhps 10 * 4(%rbp), %xmm14
    \\ movaps %xmm9, %xmm11
    \\ movaps 12 * 4(%r13), %xmm9
    \\ mulps %xmm15, %xmm11
    \\ mulps %xmm3, %xmm15
    \\ addps %xmm11, %xmm7
    \\ addps %xmm15, %xmm10
    \\ .if \bp0
    \\ movsd 8 * 4(%rbp, %r8, 1), %xmm15
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rbp, %r8, 1), %xmm15
    \\ .endif
    \\ movhps 10 * 4(%rbp, %r8, 1), %xmm15
    \\ movlps %xmm10, 4 * 4(%r14)
    \\ movhps %xmm10, 6 * 4(%r14)
    \\ .if \bp0
    \\ movsd 8 * 4(%r14), %xmm10
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%r14), %xmm10
    \\ .endif
    \\ movhps 10 * 4(%r14), %xmm10
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm12, %xmm11
    \\ mulps %xmm0, %xmm12
    \\ addps %xmm11, %xmm4
    \\ addps %xmm12, %xmm10
    \\ .if \bp0
    \\ movsd 12 * 4(%rbx), %xmm12
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%rbx), %xmm12
    \\ .endif
    \\ movhps 14 * 4(%rbx), %xmm12
    \\ .if \bp2
    \\ prefetcht0 (16 * 12)(%rbp)
    \\ .endif
    \\ .if \bp3
    \\ prefetch (16 * 16)(%rbp)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 (16 * 24)(%rbp)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 8)(%rbp)
    \\ .endif
    \\ .if \bp5
    \\ prefetcht0 (16 * 20)(%rbp)
    \\ .endif
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm13, %xmm11
    \\ mulps %xmm1, %xmm13
    \\ addps %xmm11, %xmm5
    \\ addps %xmm13, %xmm10
    \\ .if \bp0
    \\ movsd 12 * 4(%rbx, %r8, 1), %xmm13
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%rbx, %r8, 1), %xmm13
    \\ .endif
    \\ movhps 14 * 4(%rbx, %r8, 1), %xmm13
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm14, %xmm11
    \\ mulps %xmm2, %xmm14
    \\ addps %xmm11, %xmm6
    \\ addps %xmm14, %xmm10
    \\ .if \bp0
    \\ movsd 12 * 4(%rbp), %xmm14
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%rbp), %xmm14
    \\ .endif
    \\ movhps 14 * 4(%rbp), %xmm14
    \\ .if \bp6
    \\ prefetcht0 (16 * 12)(%r14)
    \\ .endif
    \\ .if \bp3
    \\ prefetchw (16 * 16)(%r14)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 (16 * 24)(%r14)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw (16 * 8)(%r14)
    \\ .endif
    \\ .if \bp5
    \\ prefetcht0 (16 * 20)(%r14)
    \\ .endif
    \\ movaps %xmm8, %xmm11
    \\ movaps 16 * 4(%r13), %xmm8
    \\ mulps %xmm15, %xmm11
    \\ mulps %xmm3, %xmm15
    \\ addps %xmm11, %xmm7
    \\ addps %xmm15, %xmm10
    \\ .if \bp0
    \\ movsd 12 * 4(%rbp, %r8, 1), %xmm15
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%rbp, %r8, 1), %xmm15
    \\ .endif
    \\ movhps 14 * 4(%rbp, %r8, 1), %xmm15
    \\ movlps %xmm10, 8 * 4(%r14)
    \\ movhps %xmm10, 10 * 4(%r14)
    \\ .if \bp0
    \\ movsd 12 * 4(%r14), %xmm10
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%r14), %xmm10
    \\ .endif
    \\ movhps 14 * 4(%r14), %xmm10
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm12, %xmm11
    \\ mulps %xmm0, %xmm12
    \\ addps %xmm11, %xmm4
    \\ addps %xmm12, %xmm10
    \\ .if \bp0
    \\ movsd 16 * 4(%rbx), %xmm12
    \\ .endif
    \\ .if \bp1
    \\ movlps 16 * 4(%rbx), %xmm12
    \\ .endif
    \\ movhps 18 * 4(%rbx), %xmm12
    \\ .if \bp2
    \\ prefetcht0 (16 * 12)(%rbp, %r8, 1)
    \\ .endif
    \\ .if \bp3
    \\ prefetch (16 * 16)(%rbp, %r8, 1)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 (16 * 24)(%rbp, %r8, 1)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 8)(%rbp, %r8, 1)
    \\ .endif
    \\ .if \bp5
    \\ prefetcht0 (16 * 20)(%rbp, %r8, 1)
    \\ .endif
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm13, %xmm11
    \\ mulps %xmm1, %xmm13
    \\ addps %xmm11, %xmm5
    \\ addps %xmm13, %xmm10
    \\ .if \bp0
    \\ movsd 16 * 4(%rbx, %r8, 1), %xmm13
    \\ .endif
    \\ .if \bp1
    \\ movlps 16 * 4(%rbx, %r8, 1), %xmm13
    \\ .endif
    \\ movhps 18 * 4(%rbx, %r8, 1), %xmm13
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm14, %xmm11
    \\ mulps %xmm2, %xmm14
    \\ addps %xmm11, %xmm6
    \\ addps %xmm14, %xmm10
    \\ .if \bp0
    \\ movsd 16 * 4(%rbp), %xmm14
    \\ .endif
    \\ .if \bp1
    \\ movlps 16 * 4(%rbp), %xmm14
    \\ .endif
    \\ movhps 18 * 4(%rbp), %xmm14
    \\ movaps %xmm9, %xmm11
    \\ movaps 20 * 4(%r13), %xmm9
    \\ mulps %xmm15, %xmm11
    \\ mulps %xmm3, %xmm15
    \\ addps %xmm11, %xmm7
    \\ addps %xmm15, %xmm10
    \\ .if \bp0
    \\ movsd 16 * 4(%rbp, %r8, 1), %xmm15
    \\ .endif
    \\ .if \bp1
    \\ movlps 16 * 4(%rbp, %r8, 1), %xmm15
    \\ .endif
    \\ movhps 18 * 4(%rbp, %r8, 1), %xmm15
    \\ movlps %xmm10, 12 * 4(%r14)
    \\ movhps %xmm10, 14 * 4(%r14)
    \\ .if \bp0
    \\ movsd 16 * 4(%r14), %xmm10
    \\ .endif
    \\ .if \bp1
    \\ movlps 16 * 4(%r14), %xmm10
    \\ .endif
    \\ movhps 18 * 4(%r14), %xmm10
    \\ addq $16 * 4, %r13
    \\ addq $16 * 4, %r14
    \\ addq $16 * 4, %rbx
    \\ addq $16 * 4, %rbp
    \\ decq %rax
    \\ jg 20009b
    \\ .balign 8
    \\ 20010:
    \\ testq $8, %rdx
    \\ jle 20011f
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm12, %xmm11
    \\ mulps %xmm0, %xmm12
    \\ addps %xmm11, %xmm4
    \\ addps %xmm12, %xmm10
    \\ .if \bp0
    \\ movsd 4 * 4(%rbx), %xmm12
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rbx), %xmm12
    \\ .endif
    \\ movhps 6 * 4(%rbx), %xmm12
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm13, %xmm11
    \\ mulps %xmm1, %xmm13
    \\ addps %xmm11, %xmm5
    \\ addps %xmm13, %xmm10
    \\ .if \bp0
    \\ movsd 4 * 4(%rbx, %r8, 1), %xmm13
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rbx, %r8, 1), %xmm13
    \\ .endif
    \\ movhps 6 * 4(%rbx, %r8, 1), %xmm13
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm14, %xmm11
    \\ mulps %xmm2, %xmm14
    \\ addps %xmm11, %xmm6
    \\ addps %xmm14, %xmm10
    \\ .if \bp0
    \\ movsd 4 * 4(%rbp), %xmm14
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rbp), %xmm14
    \\ .endif
    \\ movhps 6 * 4(%rbp), %xmm14
    \\ movaps %xmm8, %xmm11
    \\ movaps 8 * 4(%r13), %xmm8
    \\ mulps %xmm15, %xmm11
    \\ mulps %xmm3, %xmm15
    \\ addps %xmm11, %xmm7
    \\ addps %xmm15, %xmm10
    \\ .if \bp0
    \\ movsd 4 * 4(%rbp, %r8, 1), %xmm15
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rbp, %r8, 1), %xmm15
    \\ .endif
    \\ movhps 6 * 4(%rbp, %r8, 1), %xmm15
    \\ movlps %xmm10, 0 * 4(%r14)
    \\ movhps %xmm10, 2 * 4(%r14)
    \\ .if \bp0
    \\ movsd 4 * 4(%r14), %xmm10
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%r14), %xmm10
    \\ .endif
    \\ movhps 6 * 4(%r14), %xmm10
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm12, %xmm11
    \\ mulps %xmm0, %xmm12
    \\ addps %xmm11, %xmm4
    \\ addps %xmm12, %xmm10
    \\ .if \bp0
    \\ movsd 8 * 4(%rbx), %xmm12
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rbx), %xmm12
    \\ .endif
    \\ movhps 10 * 4(%rbx), %xmm12
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm13, %xmm11
    \\ mulps %xmm1, %xmm13
    \\ addps %xmm11, %xmm5
    \\ addps %xmm13, %xmm10
    \\ .if \bp0
    \\ movsd 8 * 4(%rbx, %r8, 1), %xmm13
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rbx, %r8, 1), %xmm13
    \\ .endif
    \\ movhps 10 * 4(%rbx, %r8, 1), %xmm13
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm14, %xmm11
    \\ mulps %xmm2, %xmm14
    \\ addps %xmm11, %xmm6
    \\ addps %xmm14, %xmm10
    \\ .if \bp0
    \\ movsd 8 * 4(%rbp), %xmm14
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rbp), %xmm14
    \\ .endif
    \\ movhps 10 * 4(%rbp), %xmm14
    \\ movaps %xmm9, %xmm11
    \\ movaps 12 * 4(%r13), %xmm9
    \\ mulps %xmm15, %xmm11
    \\ mulps %xmm3, %xmm15
    \\ addps %xmm11, %xmm7
    \\ addps %xmm15, %xmm10
    \\ .if \bp0
    \\ movsd 8 * 4(%rbp, %r8, 1), %xmm15
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rbp, %r8, 1), %xmm15
    \\ .endif
    \\ movhps 10 * 4(%rbp, %r8, 1), %xmm15
    \\ movlps %xmm10, 4 * 4(%r14)
    \\ movhps %xmm10, 6 * 4(%r14)
    \\ .if \bp0
    \\ movsd 8 * 4(%r14), %xmm10
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%r14), %xmm10
    \\ .endif
    \\ movhps 10 * 4(%r14), %xmm10
    \\ addq $8 * 4, %r13
    \\ addq $8 * 4, %r14
    \\ addq $8 * 4, %rbx
    \\ addq $8 * 4, %rbp
    \\
;
