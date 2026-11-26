//! Section 4 of 7 of the text of `caxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zaxpy_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20029:
    \\ movaps -16 * 4(%rcx), %xmm5
    \\ movaps -12 * 4(%rcx), %xmm6
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%r9)
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
    \\ movaps -8 * 4(%rcx), %xmm7
    \\ movaps -4 * 4(%rcx), %xmm0
    \\ .if \bp5
    \\ prefetcht0 (( 64 * 12 + 32) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 ((128 * 10) + 0) - 128(%rcx)
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
    \\ movaps 0 * 4(%rcx), %xmm1
    \\ movaps 4 * 4(%rcx), %xmm2
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%r9)
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
    \\ movaps 8 * 4(%rcx), %xmm3
    \\ movaps 12 * 4(%rcx), %xmm4
    \\ .if \bp5
    \\ prefetcht0 (( 64 * 12 + 32) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rcx)
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
    \\ jg 20029b
    \\ .balign 8
    \\ 20030:
    \\ movaps -16 * 4(%rcx), %xmm5
    \\ movaps -12 * 4(%rcx), %xmm6
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
    \\ movaps -8 * 4(%rcx), %xmm7
    \\ movaps -4 * 4(%rcx), %xmm0
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
    \\ 20031:
    \\ testq $8, %rdi
    \\ jle 20032f
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
    \\ movaps -24 * 4(%rcx), %xmm3
    \\ movaps -20 * 4(%rcx), %xmm0
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
    \\ 20032:
    \\ testq $4, %rdi
    \\ jle 20033f
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
    \\ 20033:
    \\ testq $2, %rdi
    \\ jle 20034f
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
    \\ 20034:
    \\ testq $1, %rdi
    \\ jle 20035f
    \\ movsd -32 * 4(%rcx), %xmm1
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ jmp 20067f
    \\ .balign 8
    \\ 20035:
    \\ shufps $0x93, %xmm0, %xmm0
    \\ addss -32 * 4(%r9), %xmm0
    \\ movss %xmm0, -32 * 4(%r9)
    \\ jmp 20067f
    \\ .balign 8
    \\ 20036:
    \\ testq $2 * 4, %rcx
    \\ jne 20044f
    \\ subq $1 * 4, %rcx
    \\ movaps -32 * 4(%rcx), %xmm1
    \\ movq %rdi, %rax
    \\ sarq $4, %rax
    \\ jle 20039f
    \\ movaps -28 * 4(%rcx), %xmm2
    \\ movaps -24 * 4(%rcx), %xmm3
    \\ decq %rax
    \\ jle 20038f
    \\ .balign 8
    \\ 20037:
    \\ movaps -20 * 4(%rcx), %xmm4
    \\ movaps -16 * 4(%rcx), %xmm5
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%r9)
    \\ .endif
    \\ movss %xmm2, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm3, %xmm2
    \\ shufps $0x39, %xmm2, %xmm2
    \\ pshufd $0xb1, %xmm2, %xmm8
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm2
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps -12 * 4(%rcx), %xmm6
    \\ movaps -8 * 4(%rcx), %xmm7
    \\ .if \bp5
    \\ prefetcht0 (( 64 * 12 + 32) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 ((128 * 10) + 0) - 128(%rcx)
    \\ .endif
    \\ movss %xmm4, %xmm3
    \\ shufps $0x39, %xmm3, %xmm3
    \\ pshufd $0xb1, %xmm3, %xmm8
    \\ mulps %xmm14, %xmm3
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm3
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ movss %xmm5, %xmm4
    \\ shufps $0x39, %xmm4, %xmm4
    \\ pshufd $0xb1, %xmm4, %xmm8
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm4
    \\ movss %xmm4, %xmm3
    \\ shufps $0x93, %xmm4, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%r9)
    \\ .endif
    \\ movaps -4 * 4(%rcx), %xmm0
    \\ movaps 0 * 4(%rcx), %xmm1
    \\ movss %xmm6, %xmm5
    \\ shufps $0x39, %xmm5, %xmm5
    \\ pshufd $0xb1, %xmm5, %xmm8
    \\ mulps %xmm14, %xmm5
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm5
    \\ movss %xmm5, %xmm4
    \\ shufps $0x93, %xmm5, %xmm4
    \\ addps -16 * 4(%r9), %xmm4
    \\ movaps %xmm4, -16 * 4(%r9)
    \\ movss %xmm7, %xmm6
    \\ shufps $0x39, %xmm6, %xmm6
    \\ pshufd $0xb1, %xmm6, %xmm8
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm6
    \\ movss %xmm6, %xmm5
    \\ shufps $0x93, %xmm6, %xmm5
    \\ addps -12 * 4(%r9), %xmm5
    \\ movaps %xmm5, -12 * 4(%r9)
    \\ movaps 4 * 4(%rcx), %xmm2
    \\ movaps 8 * 4(%rcx), %xmm3
    \\ .if \bp5
    \\ prefetcht0 (( 64 * 12 + 32) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rcx)
    \\ .endif
    \\ movss %xmm0, %xmm7
    \\ shufps $0x39, %xmm7, %xmm7
    \\ pshufd $0xb1, %xmm7, %xmm8
    \\ mulps %xmm14, %xmm7
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm7
    \\ movss %xmm7, %xmm6
    \\ shufps $0x93, %xmm7, %xmm6
    \\ addps -8 * 4(%r9), %xmm6
    \\ movaps %xmm6, -8 * 4(%r9)
    \\ movss %xmm1, %xmm0
    \\ shufps $0x39, %xmm0, %xmm0
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
    \\ jg 20037b
    \\
;
