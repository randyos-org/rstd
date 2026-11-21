//! Section 3 of 3 of the text of `sswap`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/swap_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20025:
    \\ movaps -35 * 4(%rcx), %xmm0
    \\ movaps -32 * 4(%r9), %xmm1
    \\ movss %xmm1, -32 * 4(%rcx)
    \\ subq $3, %rdi
    \\ movq %rdi, %rax
    \\ sarq $5, %rax
    \\ jle 20027f
    \\ .balign 16
    \\ 20026:
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
    \\ movaps -31 * 4(%rcx), %xmm2
    \\ movaps -28 * 4(%r9), %xmm3
    \\ movss %xmm2, %xmm0
    \\ shufps $0x93, %xmm2, %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm3, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ movaps %xmm1, -31 * 4(%rcx)
    \\ movaps -27 * 4(%rcx), %xmm0
    \\ movaps -24 * 4(%r9), %xmm1
    \\ movss %xmm0, %xmm2
    \\ shufps $0x93, %xmm0, %xmm2
    \\ movaps %xmm2, -28 * 4(%r9)
    \\ movss %xmm1, %xmm3
    \\ shufps $0x39, %xmm3, %xmm3
    \\ movaps %xmm3, -27 * 4(%rcx)
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
    \\ movaps -23 * 4(%rcx), %xmm2
    \\ movaps -20 * 4(%r9), %xmm3
    \\ movss %xmm2, %xmm0
    \\ shufps $0x93, %xmm2, %xmm0
    \\ movaps %xmm0, -24 * 4(%r9)
    \\ movss %xmm3, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ movaps %xmm1, -23 * 4(%rcx)
    \\ movaps -19 * 4(%rcx), %xmm0
    \\ movaps -16 * 4(%r9), %xmm1
    \\ movss %xmm0, %xmm2
    \\ shufps $0x93, %xmm0, %xmm2
    \\ movaps %xmm2, -20 * 4(%r9)
    \\ movss %xmm1, %xmm3
    \\ shufps $0x39, %xmm3, %xmm3
    \\ movaps %xmm3, -19 * 4(%rcx)
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
    \\ movaps -15 * 4(%rcx), %xmm2
    \\ movaps -12 * 4(%r9), %xmm3
    \\ movss %xmm2, %xmm0
    \\ shufps $0x93, %xmm2, %xmm0
    \\ movaps %xmm0, -16 * 4(%r9)
    \\ movss %xmm3, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ movaps %xmm1, -15 * 4(%rcx)
    \\ movaps -11 * 4(%rcx), %xmm0
    \\ movaps -8 * 4(%r9), %xmm1
    \\ movss %xmm0, %xmm2
    \\ shufps $0x93, %xmm0, %xmm2
    \\ movaps %xmm2, -12 * 4(%r9)
    \\ movss %xmm1, %xmm3
    \\ shufps $0x39, %xmm3, %xmm3
    \\ movaps %xmm3, -11 * 4(%rcx)
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
    \\ movaps -7 * 4(%rcx), %xmm2
    \\ movaps -4 * 4(%r9), %xmm3
    \\ movss %xmm2, %xmm0
    \\ shufps $0x93, %xmm2, %xmm0
    \\ movaps %xmm0, -8 * 4(%r9)
    \\ movss %xmm3, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ movaps %xmm1, -7 * 4(%rcx)
    \\ movaps -3 * 4(%rcx), %xmm0
    \\ movaps 0 * 4(%r9), %xmm1
    \\ movss %xmm0, %xmm2
    \\ shufps $0x93, %xmm0, %xmm2
    \\ movaps %xmm2, -4 * 4(%r9)
    \\ movss %xmm1, %xmm3
    \\ shufps $0x39, %xmm3, %xmm3
    \\ movaps %xmm3, -3 * 4(%rcx)
    \\ subq $-32 * 4, %rcx
    \\ subq $-32 * 4, %r9
    \\ decq %rax
    \\ jg 20026b
    \\ .balign 8
    \\ 20027:
    \\ testq $16, %rdi
    \\ jle 20028f
    \\ movaps -31 * 4(%rcx), %xmm2
    \\ movaps -28 * 4(%r9), %xmm3
    \\ movss %xmm2, %xmm0
    \\ shufps $0x93, %xmm2, %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm3, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ movaps %xmm1, -31 * 4(%rcx)
    \\ movaps -27 * 4(%rcx), %xmm0
    \\ movaps -24 * 4(%r9), %xmm1
    \\ movss %xmm0, %xmm2
    \\ shufps $0x93, %xmm0, %xmm2
    \\ movaps %xmm2, -28 * 4(%r9)
    \\ movss %xmm1, %xmm3
    \\ shufps $0x39, %xmm3, %xmm3
    \\ movaps %xmm3, -27 * 4(%rcx)
    \\ movaps -23 * 4(%rcx), %xmm2
    \\ movaps -20 * 4(%r9), %xmm3
    \\ movss %xmm2, %xmm0
    \\ shufps $0x93, %xmm2, %xmm0
    \\ movaps %xmm0, -24 * 4(%r9)
    \\ movss %xmm3, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ movaps %xmm1, -23 * 4(%rcx)
    \\ movaps -19 * 4(%rcx), %xmm0
    \\ movaps -16 * 4(%r9), %xmm1
    \\ movss %xmm0, %xmm2
    \\ shufps $0x93, %xmm0, %xmm2
    \\ movaps %xmm2, -20 * 4(%r9)
    \\ movss %xmm1, %xmm3
    \\ shufps $0x39, %xmm3, %xmm3
    \\ movaps %xmm3, -19 * 4(%rcx)
    \\ addq $16 * 4, %rcx
    \\ addq $16 * 4, %r9
    \\ .balign 8
    \\ 20028:
    \\ testq $8, %rdi
    \\ jle 20029f
    \\ movaps -31 * 4(%rcx), %xmm2
    \\ movaps -28 * 4(%r9), %xmm3
    \\ movss %xmm2, %xmm0
    \\ shufps $0x93, %xmm2, %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm3, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ movaps %xmm1, -31 * 4(%rcx)
    \\ movaps -27 * 4(%rcx), %xmm0
    \\ movaps -24 * 4(%r9), %xmm1
    \\ movss %xmm0, %xmm2
    \\ shufps $0x93, %xmm0, %xmm2
    \\ movaps %xmm2, -28 * 4(%r9)
    \\ movss %xmm1, %xmm3
    \\ shufps $0x39, %xmm3, %xmm3
    \\ movaps %xmm3, -27 * 4(%rcx)
    \\ addq $8 * 4, %rcx
    \\ addq $8 * 4, %r9
    \\ .balign 8
    \\ 20029:
    \\ testq $4, %rdi
    \\ jle 20030f
    \\ movaps -31 * 4(%rcx), %xmm2
    \\ movaps -28 * 4(%r9), %xmm3
    \\ movss %xmm2, %xmm0
    \\ shufps $0x93, %xmm2, %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm3, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ movaps %xmm1, -31 * 4(%rcx)
    \\ movaps %xmm2, %xmm0
    \\ movaps %xmm3, %xmm1
    \\ addq $4 * 4, %rcx
    \\ addq $4 * 4, %r9
    \\ .balign 8
    \\ 20030:
    \\ .if \bp0
    \\ movsd -31 * 4(%rcx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps -31 * 4(%rcx), %xmm2
    \\ .endif
    \\ pshufd $0x39, %xmm1, %xmm1
    \\ movlps %xmm1, -31 * 4(%rcx)
    \\ pshufd $0xff, %xmm0, %xmm0
    \\ movss %xmm0, -32 * 4(%r9)
    \\ movlps %xmm2, -31 * 4(%r9)
    \\ addq $3 * 4, %rcx
    \\ addq $3 * 4, %r9
    \\ testq $2, %rdi
    \\ jle 20031f
    \\ .if \bp0
    \\ movsd -32 * 4(%rcx), %xmm0
    \\ movsd -32 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rcx), %xmm0
    \\ movlps -32 * 4(%r9), %xmm1
    \\ .endif
    \\ movlps %xmm0, -32 * 4(%r9)
    \\ movlps %xmm1, -32 * 4(%rcx)
    \\ addq $2 * 4, %rcx
    \\ addq $2 * 4, %r9
    \\ .balign 8
    \\ 20031:
    \\ testq $1, %rdi
    \\ jle 20032f
    \\ movss -32 * 4(%rcx), %xmm0
    \\ movss -32 * 4(%r9), %xmm1
    \\ movss %xmm0, -32 * 4(%r9)
    \\ movss %xmm1, -32 * 4(%rcx)
    \\ .balign 8
    \\ 20032:
    \\ xorq %rax,%rax
    \\ ret
    \\ .balign 8
    \\ 20033:
    \\ movq %rdi, %rax
    \\ sarq $3, %rax
    \\ jle 20035f
    \\ .balign 8
    \\ 20034:
    \\ movss (%rcx), %xmm0
    \\ movss (%r9), %xmm1
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ movss (%rcx), %xmm0
    \\ movss (%r9), %xmm1
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ movss (%rcx), %xmm0
    \\ movss (%r9), %xmm1
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ movss (%rcx), %xmm0
    \\ movss (%r9), %xmm1
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ movss (%rcx), %xmm0
    \\ movss (%r9), %xmm1
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ movss (%rcx), %xmm0
    \\ movss (%r9), %xmm1
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ movss (%rcx), %xmm0
    \\ movss (%r9), %xmm1
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ movss (%rcx), %xmm0
    \\ movss (%r9), %xmm1
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ decq %rax
    \\ jg 20034b
    \\ .balign 8
    \\ 20035:
    \\ movq %rdi, %rax
    \\ andq $7, %rax
    \\ jle 20037f
    \\ .balign 8
    \\ 20036:
    \\ movss (%rcx), %xmm0
    \\ movss (%r9), %xmm1
    \\ movss %xmm1, (%rcx)
    \\ movss %xmm0, (%r9)
    \\ addq %r8, %rcx
    \\ addq %rsi, %r9
    \\ decq %rax
    \\ jg 20036b
    \\ .balign 8
    \\ 20037:
    \\ xorq %rax, %rax
    \\ ret
    \\ .endm
    \\
;
