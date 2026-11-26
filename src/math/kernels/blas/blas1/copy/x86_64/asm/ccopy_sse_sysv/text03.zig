//! Section 3 of 3 of the text of `ccopy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zcopy_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 16
    \\ 20029:
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%rcx)
    \\ .endif
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ movaps -3 * 4(%rsi), %xmm0
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ movaps %xmm1, -28 * 4(%rcx)
    \\ movaps 1 * 4(%rsi), %xmm1
    \\ .if \bp5
    \\ prefetcht0 (( 64 * 12 + 32) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp6
    \\ prefetcht0 ((128 * 10) + 0) - 128(%rsi)
    \\ .endif
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ movaps %xmm2, -24 * 4(%rcx)
    \\ movaps 5 * 4(%rsi), %xmm2
    \\ movss %xmm4, %xmm3
    \\ shufps $0x93, %xmm4, %xmm3
    \\ movaps %xmm3, -20 * 4(%rcx)
    \\ movaps 9 * 4(%rsi), %xmm3
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%rcx)
    \\ .endif
    \\ movss %xmm5, %xmm4
    \\ shufps $0x93, %xmm5, %xmm4
    \\ movaps %xmm4, -16 * 4(%rcx)
    \\ movaps 13 * 4(%rsi), %xmm4
    \\ movss %xmm6, %xmm5
    \\ shufps $0x93, %xmm6, %xmm5
    \\ movaps %xmm5, -12 * 4(%rcx)
    \\ movaps 17 * 4(%rsi), %xmm5
    \\ .if \bp5
    \\ prefetcht0 (( 64 * 12 + 32) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 64) - 128(%rsi)
    \\ .endif
    \\ movss %xmm7, %xmm6
    \\ shufps $0x93, %xmm7, %xmm6
    \\ movaps %xmm6, -8 * 4(%rcx)
    \\ movaps 21 * 4(%rsi), %xmm6
    \\ movss %xmm0, %xmm7
    \\ shufps $0x93, %xmm0, %xmm7
    \\ movaps %xmm7, -4 * 4(%rcx)
    \\ movaps 25 * 4(%rsi), %xmm7
    \\ subq $-32 * 4, %rsi
    \\ subq $-32 * 4, %rcx
    \\ decq %rax
    \\ jg 20029b
    \\ .balign 8
    \\ 20030:
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ movaps -3 * 4(%rsi), %xmm0
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ movaps %xmm1, -28 * 4(%rcx)
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ movaps %xmm2, -24 * 4(%rcx)
    \\ movss %xmm4, %xmm3
    \\ shufps $0x93, %xmm4, %xmm3
    \\ movaps %xmm3, -20 * 4(%rcx)
    \\ movss %xmm5, %xmm4
    \\ shufps $0x93, %xmm5, %xmm4
    \\ movaps %xmm4, -16 * 4(%rcx)
    \\ movss %xmm6, %xmm5
    \\ shufps $0x93, %xmm6, %xmm5
    \\ movaps %xmm5, -12 * 4(%rcx)
    \\ movss %xmm7, %xmm6
    \\ shufps $0x93, %xmm7, %xmm6
    \\ movaps %xmm6, -8 * 4(%rcx)
    \\ movss %xmm0, %xmm7
    \\ shufps $0x93, %xmm0, %xmm7
    \\ movaps %xmm7, -4 * 4(%rcx)
    \\ subq $-32 * 4, %rsi
    \\ subq $-32 * 4, %rcx
    \\ .balign 8
    \\ 20031:
    \\ testq $16, %rdi
    \\ jle 20032f
    \\ .balign 8
    \\ movaps -31 * 4(%rsi), %xmm1
    \\ movaps -27 * 4(%rsi), %xmm2
    \\ movaps -23 * 4(%rsi), %xmm3
    \\ movaps -19 * 4(%rsi), %xmm4
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ movaps %xmm1, -28 * 4(%rcx)
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ movaps %xmm2, -24 * 4(%rcx)
    \\ movss %xmm4, %xmm3
    \\ shufps $0x93, %xmm4, %xmm3
    \\ movaps %xmm3, -20 * 4(%rcx)
    \\ movaps %xmm4, %xmm0
    \\ addq $16 * 4, %rsi
    \\ addq $16 * 4, %rcx
    \\ .balign 8
    \\ 20032:
    \\ testq $8, %rdi
    \\ jle 20033f
    \\ .balign 8
    \\ movaps -31 * 4(%rsi), %xmm1
    \\ movaps -27 * 4(%rsi), %xmm2
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ movaps %xmm1, -28 * 4(%rcx)
    \\ movaps %xmm2, %xmm0
    \\ addq $8 * 4, %rsi
    \\ addq $8 * 4, %rcx
    \\ .balign 8
    \\ 20033:
    \\ testq $4, %rdi
    \\ jle 20034f
    \\ .balign 8
    \\ movaps -31 * 4(%rsi), %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ addq $4 * 4, %rsi
    \\ addq $4 * 4, %rcx
    \\ .balign 8
    \\ 20034:
    \\ testq $2, %rdi
    \\ jle 20035f
    \\ .balign 8
    \\ .if \bp0
    \\ movsd -32 * 4(%rsi), %xmm0
    \\ movsd %xmm0, -32 * 4(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rsi), %xmm0
    \\ movlps %xmm0, -32 * 4(%rcx)
    \\ .endif
    \\ addq $2 * 4, %rsi
    \\ addq $2 * 4, %rcx
    \\ .balign 8
    \\ 20035:
    \\ testq $1, %rdi
    \\ jle 20036f
    \\ .balign 8
    \\ movss -32 * 4(%rsi), %xmm0
    \\ movss %xmm0, -32 * 4(%rcx)
    \\ addq $4, %rcx
    \\ .balign 8
    \\ 20036:
    \\ xorq %rax,%rax
    \\ ret
    \\ .balign 16
    \\ 20037:
    \\ movq %rdi, %rax
    \\ sarq $3, %rax
    \\ jle 20039f
    \\ .balign 8
    \\ 20038:
    \\ .if \bp0
    \\ movsd (%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rsi), %xmm0
    \\ .endif
    \\ addq %rdx, %rsi
    \\ movhps (%rsi), %xmm0
    \\ addq %rdx, %rsi
    \\ .if \bp0
    \\ movsd (%rsi), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rsi), %xmm1
    \\ .endif
    \\ addq %rdx, %rsi
    \\ movhps (%rsi), %xmm1
    \\ addq %rdx, %rsi
    \\ .if \bp0
    \\ movsd (%rsi), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rsi), %xmm2
    \\ .endif
    \\ addq %rdx, %rsi
    \\ movhps (%rsi), %xmm2
    \\ addq %rdx, %rsi
    \\ .if \bp0
    \\ movsd (%rsi), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rsi), %xmm3
    \\ .endif
    \\ addq %rdx, %rsi
    \\ movhps (%rsi), %xmm3
    \\ addq %rdx, %rsi
    \\ .if \bp0
    \\ movsd %xmm0, (%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm0, (%rcx)
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps %xmm0, (%rcx)
    \\ addq %r8, %rcx
    \\ .if \bp0
    \\ movsd %xmm1, (%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm1, (%rcx)
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ .if \bp0
    \\ movsd %xmm2, (%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm2, (%rcx)
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps %xmm2, (%rcx)
    \\ addq %r8, %rcx
    \\ .if \bp0
    \\ movsd %xmm3, (%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm3, (%rcx)
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps %xmm3, (%rcx)
    \\ addq %r8, %rcx
    \\ decq %rax
    \\ jg 20038b
    \\ .balign 8
    \\ 20039:
    \\ testq $4, %rdi
    \\ jle 20040f
    \\ .if \bp0
    \\ movsd (%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rsi), %xmm0
    \\ .endif
    \\ addq %rdx, %rsi
    \\ movhps (%rsi), %xmm0
    \\ addq %rdx, %rsi
    \\ .if \bp0
    \\ movsd (%rsi), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rsi), %xmm1
    \\ .endif
    \\ addq %rdx, %rsi
    \\ movhps (%rsi), %xmm1
    \\ addq %rdx, %rsi
    \\ .if \bp0
    \\ movsd %xmm0, (%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm0, (%rcx)
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps %xmm0, (%rcx)
    \\ addq %r8, %rcx
    \\ .if \bp0
    \\ movsd %xmm1, (%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm1, (%rcx)
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ .balign 8
    \\ 20040:
    \\ testq $2, %rdi
    \\ jle 20041f
    \\ .if \bp0
    \\ movsd (%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rsi), %xmm0
    \\ .endif
    \\ addq %rdx, %rsi
    \\ movhps (%rsi), %xmm0
    \\ addq %rdx, %rsi
    \\ .if \bp0
    \\ movsd %xmm0, (%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm0, (%rcx)
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps %xmm0, (%rcx)
    \\ addq %r8, %rcx
    \\ .balign 8
    \\ 20041:
    \\ testq $1, %rdi
    \\ jle 20042f
    \\ .if \bp0
    \\ movsd (%rsi), %xmm0
    \\ movsd %xmm0, (%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rsi), %xmm0
    \\ movlps %xmm0, (%rcx)
    \\ .endif
    \\ .balign 8
    \\ 20042:
    \\ xorq %rax, %rax
    \\ ret
    \\ .endm
    \\
;
