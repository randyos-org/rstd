//! Section 2 of 3 of the text of `drot`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/rot_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20003:
    \\ testq $15, %rcx
    \\ jle 20017f
    \\ testq $8, %rcx
    \\ jle 20004f
    \\ movaps 0 * 8(%r9), %xmm1
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps 2 * 8(%r9), %xmm3
    \\ movaps 2 * 8(%rdx), %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 0 * 8(%rdx)
    \\ movaps %xmm2, 2 * 8(%rdx)
    \\ movaps %xmm4, 0 * 8(%r9)
    \\ movaps %xmm6, 2 * 8(%r9)
    \\ movaps 4 * 8(%r9), %xmm1
    \\ movaps 4 * 8(%rdx), %xmm0
    \\ movaps 6 * 8(%r9), %xmm3
    \\ movaps 6 * 8(%rdx), %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 4 * 8(%rdx)
    \\ movaps %xmm2, 6 * 8(%rdx)
    \\ movaps %xmm4, 4 * 8(%r9)
    \\ movaps %xmm6, 6 * 8(%r9)
    \\ addq $8 * 8, %rdx
    \\ addq $8 * 8, %r9
    \\ .balign 8
    \\ 20004:
    \\ testq $4, %rcx
    \\ jle 20005f
    \\ movaps 0 * 8(%r9), %xmm1
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps 2 * 8(%r9), %xmm3
    \\ movaps 2 * 8(%rdx), %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 0 * 8(%rdx)
    \\ movaps %xmm2, 2 * 8(%rdx)
    \\ movaps %xmm4, 0 * 8(%r9)
    \\ movaps %xmm6, 2 * 8(%r9)
    \\ addq $4 * 8, %rdx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20005:
    \\ testq $2, %rcx
    \\ jle 20006f
    \\ movaps 0 * 8(%r9), %xmm1
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm1, %xmm0
    \\ subpd %xmm3, %xmm2
    \\ movaps %xmm0, 0 * 8(%rdx)
    \\ movaps %xmm2, 0 * 8(%r9)
    \\ addq $2 * 8, %rdx
    \\ addq $2 * 8, %r9
    \\ .balign 8
    \\ 20006:
    \\ testq $1, %rcx
    \\ jle 20017f
    \\ .if \bp0
    \\ movsd 0 * 8(%r9), %xmm1
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 8(%r9), %xmm1
    \\ movlps 0 * 8(%rdx), %xmm0
    \\ .endif
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulsd %xmm14, %xmm0
    \\ mulsd %xmm15, %xmm1
    \\ mulsd %xmm14, %xmm2
    \\ mulsd %xmm15, %xmm3
    \\ addsd %xmm1, %xmm0
    \\ subsd %xmm3, %xmm2
    \\ .if \bp0
    \\ movsd %xmm0, 0 * 8(%rdx)
    \\ movsd %xmm2, 0 * 8(%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm0, 0 * 8(%rdx)
    \\ movlps %xmm2, 0 * 8(%r9)
    \\ .endif
    \\ jmp 20017f
    \\ .balign 8
    \\ 20007:
    \\ movaps -1 * 8(%r9), %xmm1
    \\ movq %rcx, %rax
    \\ sarq $4, %rax
    \\ jle 20009f
    \\ .balign 8
    \\ 20008:
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%rdx)
    \\ .endif
    \\ movaps 1 * 8(%r9), %xmm3
    \\ movaps 3 * 8(%r9), %xmm8
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps 2 * 8(%rdx), %xmm2
    \\ shufpd $1, %xmm3, %xmm1
    \\ shufpd $1, %xmm8, %xmm3
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 0 * 8(%rdx)
    \\ movaps %xmm2, 2 * 8(%rdx)
    \\ .if \bp5
    \\ movlpd %xmm4, 0 * 8(%r9)
    \\ .endif
    \\ .if \bp6
    \\ movlps %xmm4, 0 * 8(%r9)
    \\ .endif
    \\ movhps %xmm4, 1 * 8(%r9)
    \\ .if \bp5
    \\ movlpd %xmm6, 2 * 8(%r9)
    \\ .endif
    \\ .if \bp6
    \\ movlps %xmm6, 2 * 8(%r9)
    \\ .endif
    \\ movhps %xmm6, 3 * 8(%r9)
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%r9)
    \\ .endif
    \\ movaps 5 * 8(%r9), %xmm9
    \\ movaps 7 * 8(%r9), %xmm1
    \\ movaps 4 * 8(%rdx), %xmm0
    \\ movaps 6 * 8(%rdx), %xmm2
    \\ shufpd $1, %xmm9, %xmm8
    \\ shufpd $1, %xmm1, %xmm9
    \\ movaps %xmm8, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm9, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm9
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm8, %xmm0
    \\ addpd %xmm9, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 4 * 8(%rdx)
    \\ movaps %xmm2, 6 * 8(%rdx)
    \\ .if \bp5
    \\ movlpd %xmm4, 4 * 8(%r9)
    \\ .endif
    \\ .if \bp6
    \\ movlps %xmm4, 4 * 8(%r9)
    \\ .endif
    \\ movhps %xmm4, 5 * 8(%r9)
    \\ .if \bp5
    \\ movlpd %xmm6, 6 * 8(%r9)
    \\ .endif
    \\ .if \bp6
    \\ movlps %xmm6, 6 * 8(%r9)
    \\ .endif
    \\ movhps %xmm6, 7 * 8(%r9)
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%rdx)
    \\ .endif
    \\ movaps 9 * 8(%r9), %xmm3
    \\ movaps 11 * 8(%r9), %xmm8
    \\ movaps 8 * 8(%rdx), %xmm0
    \\ movaps 10 * 8(%rdx), %xmm2
    \\ shufpd $1, %xmm3, %xmm1
    \\ shufpd $1, %xmm8, %xmm3
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 8 * 8(%rdx)
    \\ movaps %xmm2, 10 * 8(%rdx)
    \\ .if \bp5
    \\ movlpd %xmm4, 8 * 8(%r9)
    \\ .endif
    \\ .if \bp6
    \\ movlps %xmm4, 8 * 8(%r9)
    \\ .endif
    \\ movhps %xmm4, 9 * 8(%r9)
    \\ .if \bp5
    \\ movlpd %xmm6, 10 * 8(%r9)
    \\ .endif
    \\ .if \bp6
    \\ movlps %xmm6, 10 * 8(%r9)
    \\ .endif
    \\ movhps %xmm6, 11 * 8(%r9)
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%r9)
    \\ .endif
    \\ movaps 13 * 8(%r9), %xmm9
    \\ movaps 15 * 8(%r9), %xmm1
    \\ movaps 12 * 8(%rdx), %xmm0
    \\ movaps 14 * 8(%rdx), %xmm2
    \\ shufpd $1, %xmm9, %xmm8
    \\ shufpd $1, %xmm1, %xmm9
    \\ movaps %xmm8, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm9, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm9
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm8, %xmm0
    \\ addpd %xmm9, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 12 * 8(%rdx)
    \\ movaps %xmm2, 14 * 8(%rdx)
    \\ .if \bp5
    \\ movlpd %xmm4, 12 * 8(%r9)
    \\ .endif
    \\ .if \bp6
    \\ movlps %xmm4, 12 * 8(%r9)
    \\ .endif
    \\ movhps %xmm4, 13 * 8(%r9)
    \\ .if \bp5
    \\ movlpd %xmm6, 14 * 8(%r9)
    \\ .endif
    \\ .if \bp6
    \\ movlps %xmm6, 14 * 8(%r9)
    \\ .endif
    \\ movhps %xmm6, 15 * 8(%r9)
    \\ addq $16 * 8, %rdx
    \\ addq $16 * 8, %r9
    \\ decq %rax
    \\ jg 20008b
    \\
;
