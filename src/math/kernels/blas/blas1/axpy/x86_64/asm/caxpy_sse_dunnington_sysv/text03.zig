//! Section 3 of 4 of the text of `caxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zaxpy_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 4
    \\ 20018:
    \\ testq $4, %rdi
    \\ jle 20019f
    \\ movaps -32 * 4(%rcx), %xmm1
    \\ movaps -28 * 4(%rcx), %xmm2
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ pshufd $0xb1, %xmm2, %xmm8
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm2
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps %xmm2, %xmm0
    \\ addq $8 * 4, %rcx
    \\ addq $8 * 4, %r9
    \\ .balign 4
    \\ 20019:
    \\ testq $2, %rdi
    \\ jle 20020f
    \\ movaps -32 * 4(%rcx), %xmm1
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, %xmm0
    \\ addq $4 * 4, %rcx
    \\ addq $4 * 4, %r9
    \\ .balign 4
    \\ 20020:
    \\ testq $1, %rdi
    \\ jle 20021f
    \\ .if \bp0
    \\ movsd -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ xorps %xmm1, %xmm1
    \\ movlps -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ jmp 20037f
    \\ .balign 8
    \\ 20021:
    \\ shufps $0x93, %xmm0, %xmm0
    \\ addss -32 * 4(%r9), %xmm0
    \\ movss %xmm0, -32 * 4(%r9)
    \\ jmp 20037f
    \\ .balign 8
    \\ 20022:
    \\ movq %rdi, %rax
    \\ sarq $4, %rax
    \\ jle 20025f
    \\ .if \bp0
    \\ movsd -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps -30 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd -28 * 4(%rcx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps -28 * 4(%rcx), %xmm2
    \\ .endif
    \\ movhps -26 * 4(%rcx), %xmm2
    \\ .if \bp0
    \\ movsd -24 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps -24 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps -22 * 4(%rcx), %xmm3
    \\ .if \bp0
    \\ movsd -20 * 4(%rcx), %xmm4
    \\ .endif
    \\ .if \bp1
    \\ movlps -20 * 4(%rcx), %xmm4
    \\ .endif
    \\ movhps -18 * 4(%rcx), %xmm4
    \\ decq %rax
    \\ jle 20024f
    \\ .balign 8
    \\ 20023:
    \\ .if \bp0
    \\ movsd -16 * 4(%rcx), %xmm5
    \\ .endif
    \\ .if \bp1
    \\ movlps -16 * 4(%rcx), %xmm5
    \\ .endif
    \\ movhps -14 * 4(%rcx), %xmm5
    \\ .if \bp0
    \\ movsd -12 * 4(%rcx), %xmm6
    \\ .endif
    \\ .if \bp1
    \\ movlps -12 * 4(%rcx), %xmm6
    \\ .endif
    \\ movhps -10 * 4(%rcx), %xmm6
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%r9)
    \\ .endif
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ pshufd $0xb1, %xmm2, %xmm8
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm2
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ .if \bp0
    \\ movsd -8 * 4(%rcx), %xmm7
    \\ .endif
    \\ .if \bp1
    \\ movlps -8 * 4(%rcx), %xmm7
    \\ .endif
    \\ movhps -6 * 4(%rcx), %xmm7
    \\ .if \bp0
    \\ movsd -4 * 4(%rcx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps -4 * 4(%rcx), %xmm0
    \\ .endif
    \\ movhps -2 * 4(%rcx), %xmm0
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 0) - 128(%rcx)
    \\ .endif
    \\ pshufd $0xb1, %xmm3, %xmm8
    \\ mulps %xmm14, %xmm3
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm3
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ pshufd $0xb1, %xmm4, %xmm8
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm4
    \\ movss %xmm4, %xmm3
    \\ shufps $0x93, %xmm4, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ .if \bp0
    \\ movsd 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd 4 * 4(%rcx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rcx), %xmm2
    \\ .endif
    \\ movhps 6 * 4(%rcx), %xmm2
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%r9)
    \\ .endif
    \\ pshufd $0xb1, %xmm5, %xmm8
    \\ mulps %xmm14, %xmm5
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm5
    \\ movss %xmm5, %xmm4
    \\ shufps $0x93, %xmm5, %xmm4
    \\ addps -16 * 4(%r9), %xmm4
    \\ movaps %xmm4, -16 * 4(%r9)
    \\ pshufd $0xb1, %xmm6, %xmm8
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm6
    \\ movss %xmm6, %xmm5
    \\ shufps $0x93, %xmm6, %xmm5
    \\ addps -12 * 4(%r9), %xmm5
    \\ movaps %xmm5, -12 * 4(%r9)
    \\ .if \bp0
    \\ movsd 8 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 10 * 4(%rcx), %xmm3
    \\ .if \bp0
    \\ movsd 12 * 4(%rcx), %xmm4
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%rcx), %xmm4
    \\ .endif
    \\ movhps 14 * 4(%rcx), %xmm4
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 64) - 128(%rcx)
    \\ .endif
    \\ pshufd $0xb1, %xmm7, %xmm8
    \\ mulps %xmm14, %xmm7
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm7
    \\ movss %xmm7, %xmm6
    \\ shufps $0x93, %xmm7, %xmm6
    \\ addps -8 * 4(%r9), %xmm6
    \\ movaps %xmm6, -8 * 4(%r9)
    \\ pshufd $0xb1, %xmm0, %xmm8
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm0
    \\ movss %xmm0, %xmm7
    \\ shufps $0x93, %xmm0, %xmm7
    \\ addps -4 * 4(%r9), %xmm7
    \\ movaps %xmm7, -4 * 4(%r9)
    \\ subq $-32 * 4, %rcx
    \\ subq $-32 * 4, %r9
    \\ decq %rax
    \\ jg 20023b
    \\ .balign 8
    \\ 20024:
    \\ .if \bp0
    \\ movsd -16 * 4(%rcx), %xmm5
    \\ .endif
    \\ .if \bp1
    \\ movlps -16 * 4(%rcx), %xmm5
    \\ .endif
    \\ movhps -14 * 4(%rcx), %xmm5
    \\ .if \bp0
    \\ movsd -12 * 4(%rcx), %xmm6
    \\ .endif
    \\ .if \bp1
    \\ movlps -12 * 4(%rcx), %xmm6
    \\ .endif
    \\ movhps -10 * 4(%rcx), %xmm6
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ pshufd $0xb1, %xmm2, %xmm8
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm2
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ .if \bp0
    \\ movsd -8 * 4(%rcx), %xmm7
    \\ .endif
    \\ .if \bp1
    \\ movlps -8 * 4(%rcx), %xmm7
    \\ .endif
    \\ movhps -6 * 4(%rcx), %xmm7
    \\ .if \bp0
    \\ movsd -4 * 4(%rcx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps -4 * 4(%rcx), %xmm0
    \\ .endif
    \\ movhps -2 * 4(%rcx), %xmm0
    \\ pshufd $0xb1, %xmm3, %xmm8
    \\ mulps %xmm14, %xmm3
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm3
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ pshufd $0xb1, %xmm4, %xmm8
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm4
    \\ movss %xmm4, %xmm3
    \\ shufps $0x93, %xmm4, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ pshufd $0xb1, %xmm5, %xmm8
    \\ mulps %xmm14, %xmm5
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm5
    \\ movss %xmm5, %xmm4
    \\ shufps $0x93, %xmm5, %xmm4
    \\ addps -16 * 4(%r9), %xmm4
    \\ movaps %xmm4, -16 * 4(%r9)
    \\ pshufd $0xb1, %xmm6, %xmm8
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm6
    \\ movss %xmm6, %xmm5
    \\ shufps $0x93, %xmm6, %xmm5
    \\ addps -12 * 4(%r9), %xmm5
    \\ movaps %xmm5, -12 * 4(%r9)
    \\ pshufd $0xb1, %xmm7, %xmm8
    \\ mulps %xmm14, %xmm7
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm7
    \\ movss %xmm7, %xmm6
    \\ shufps $0x93, %xmm7, %xmm6
    \\ addps -8 * 4(%r9), %xmm6
    \\ movaps %xmm6, -8 * 4(%r9)
    \\ pshufd $0xb1, %xmm0, %xmm8
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm0
    \\ movss %xmm0, %xmm7
    \\ shufps $0x93, %xmm0, %xmm7
    \\ addps -4 * 4(%r9), %xmm7
    \\ movaps %xmm7, -4 * 4(%r9)
    \\ subq $-32 * 4, %rcx
    \\ subq $-32 * 4, %r9
    \\ .balign 8
    \\ 20025:
    \\ testq $8, %rdi
    \\ jle 20026f
    \\ .if \bp0
    \\ movsd -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps -30 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd -28 * 4(%rcx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps -28 * 4(%rcx), %xmm2
    \\ .endif
    \\ movhps -26 * 4(%rcx), %xmm2
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ pshufd $0xb1, %xmm2, %xmm8
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm2
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ .if \bp0
    \\ movsd -24 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps -24 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps -22 * 4(%rcx), %xmm3
    \\ .if \bp0
    \\ movsd -20 * 4(%rcx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps -20 * 4(%rcx), %xmm0
    \\ .endif
    \\ movhps -18 * 4(%rcx), %xmm0
    \\ pshufd $0xb1, %xmm3, %xmm8
    \\ mulps %xmm14, %xmm3
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm3
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ pshufd $0xb1, %xmm0, %xmm8
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm0
    \\ movss %xmm0, %xmm3
    \\ shufps $0x93, %xmm0, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ addq $16 * 4, %rcx
    \\ addq $16 * 4, %r9
    \\ .balign 4
    \\ 20026:
    \\ testq $4, %rdi
    \\ jle 20027f
    \\ .if \bp0
    \\ movsd -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps -30 * 4(%rcx), %xmm1
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ .if \bp0
    \\ movsd -28 * 4(%rcx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps -28 * 4(%rcx), %xmm2
    \\ .endif
    \\ movhps -26 * 4(%rcx), %xmm2
    \\ pshufd $0xb1, %xmm2, %xmm8
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm2
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps %xmm2, %xmm0
    \\ addq $8 * 4, %rcx
    \\ addq $8 * 4, %r9
    \\
;
