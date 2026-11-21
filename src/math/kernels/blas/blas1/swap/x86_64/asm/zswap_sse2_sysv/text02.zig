//! Section 2 of 2 of the text of `zswap`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zswap_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20013:
    \\ testq $4, %rdi
    \\ jle 20014f
    \\ movaps -15 * 8(%r9), %xmm2
    \\ movaps -14 * 8(%rcx), %xmm3
    \\ shufpd $1, %xmm2, %xmm0
    \\ movaps %xmm0, -16 * 8(%rcx)
    \\ shufpd $1, %xmm3, %xmm1
    \\ movaps %xmm1, -15 * 8(%r9)
    \\ movaps -13 * 8(%r9), %xmm0
    \\ movaps -12 * 8(%rcx), %xmm1
    \\ shufpd $1, %xmm0, %xmm2
    \\ movaps %xmm2, -14 * 8(%rcx)
    \\ shufpd $1, %xmm1, %xmm3
    \\ movaps %xmm3, -13 * 8(%r9)
    \\ movaps -11 * 8(%r9), %xmm2
    \\ movaps -10 * 8(%rcx), %xmm3
    \\ shufpd $1, %xmm2, %xmm0
    \\ movaps %xmm0, -12 * 8(%rcx)
    \\ shufpd $1, %xmm3, %xmm1
    \\ movaps %xmm1, -11 * 8(%r9)
    \\ movaps -9 * 8(%r9), %xmm0
    \\ movaps -8 * 8(%rcx), %xmm1
    \\ shufpd $1, %xmm0, %xmm2
    \\ movaps %xmm2, -10 * 8(%rcx)
    \\ shufpd $1, %xmm1, %xmm3
    \\ movaps %xmm3, -9 * 8(%r9)
    \\ addq $8 * 8, %rcx
    \\ addq $8 * 8, %r9
    \\ .balign 8
    \\ 20014:
    \\ testq $2, %rdi
    \\ jle 20015f
    \\ movaps -15 * 8(%r9), %xmm2
    \\ movaps -14 * 8(%rcx), %xmm3
    \\ shufpd $1, %xmm2, %xmm0
    \\ movaps %xmm0, -16 * 8(%rcx)
    \\ shufpd $1, %xmm3, %xmm1
    \\ movaps %xmm1, -15 * 8(%r9)
    \\ movaps -13 * 8(%r9), %xmm0
    \\ movaps -12 * 8(%rcx), %xmm1
    \\ shufpd $1, %xmm0, %xmm2
    \\ movaps %xmm2, -14 * 8(%rcx)
    \\ shufpd $1, %xmm1, %xmm3
    \\ movaps %xmm3, -13 * 8(%r9)
    \\ addq $4 * 8, %rcx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20015:
    \\ testq $1, %rdi
    \\ jle 20016f
    \\ movaps -15 * 8(%r9), %xmm2
    \\ movaps -14 * 8(%rcx), %xmm3
    \\ shufpd $1, %xmm3, %xmm1
    \\ movaps %xmm1, -15 * 8(%r9)
    \\ shufpd $1, %xmm2, %xmm0
    \\ movaps %xmm0, -16 * 8(%rcx)
    \\ movaps %xmm2, %xmm0
    \\ movaps %xmm3, %xmm1
    \\ addq $2 * 8, %rcx
    \\ addq $2 * 8, %r9
    \\ .balign 8
    \\ 20016:
    \\ movaps -15 * 8(%r9), %xmm2
    \\ movhps %xmm1, -15 * 8(%r9)
    \\ shufpd $1, %xmm2, %xmm0
    \\ movaps %xmm0, -16 * 8(%rcx)
    \\ xorq %rax,%rax
    \\ ret
    \\ .balign 8
    \\ 20017:
    \\ .if \bp4
    \\ movsd -16 * 8(%rcx), %xmm0
    \\ movsd -16 * 8(%r9), %xmm1
    \\ .endif
    \\ .if \bp3
    \\ movlps -16 * 8(%rcx), %xmm0
    \\ movlps -16 * 8(%r9), %xmm1
    \\ .endif
    \\ movlps %xmm0, -16 * 8(%r9)
    \\ movlps %xmm1, -16 * 8(%rcx)
    \\ addq $8, %rcx
    \\ addq $8, %r9
    \\ decq %rdi
    \\ jle 20022f
    \\ movq %rdi, %rax
    \\ sarq $3, %rax
    \\ jle 20019f
    \\ .balign 8
    \\ 20018:
    \\ .if \bp0
    \\ prefetchw ((128 * 5) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp3
    \\ prefetchw ((128 * 3) + 0) - 128(%rcx)
    \\ .endif
    \\ movaps -16 * 8(%rcx), %xmm0
    \\ movaps -16 * 8(%r9), %xmm1
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps %xmm1, -16 * 8(%rcx)
    \\ movaps -14 * 8(%rcx), %xmm0
    \\ movaps -14 * 8(%r9), %xmm1
    \\ movaps %xmm0, -14 * 8(%r9)
    \\ movaps %xmm1, -14 * 8(%rcx)
    \\ .if \bp0
    \\ prefetchw ((128 * 5) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 ((128 * 4) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 ((128 * 12) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp3
    \\ prefetchw ((128 * 3) + 0) - 128(%r9)
    \\ .endif
    \\ movaps -12 * 8(%rcx), %xmm0
    \\ movaps -12 * 8(%r9), %xmm1
    \\ movaps %xmm0, -12 * 8(%r9)
    \\ movaps %xmm1, -12 * 8(%rcx)
    \\ movaps -10 * 8(%rcx), %xmm0
    \\ movaps -10 * 8(%r9), %xmm1
    \\ movaps %xmm0, -10 * 8(%r9)
    \\ movaps %xmm1, -10 * 8(%rcx)
    \\ .if \bp0
    \\ prefetchw ((128 * 5) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp3
    \\ prefetchw ((128 * 3) + 64) - 128(%rcx)
    \\ .endif
    \\ movaps -8 * 8(%rcx), %xmm0
    \\ movaps -8 * 8(%r9), %xmm1
    \\ movaps %xmm0, -8 * 8(%r9)
    \\ movaps %xmm1, -8 * 8(%rcx)
    \\ movaps -6 * 8(%rcx), %xmm0
    \\ movaps -6 * 8(%r9), %xmm1
    \\ movaps %xmm0, -6 * 8(%r9)
    \\ movaps %xmm1, -6 * 8(%rcx)
    \\ .if \bp0
    \\ prefetchw ((128 * 5) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 ((128 * 4) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 ((128 * 12) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp3
    \\ prefetchw ((128 * 3) + 64) - 128(%r9)
    \\ .endif
    \\ movaps -4 * 8(%rcx), %xmm0
    \\ movaps -4 * 8(%r9), %xmm1
    \\ movaps %xmm0, -4 * 8(%r9)
    \\ movaps %xmm1, -4 * 8(%rcx)
    \\ movaps -2 * 8(%rcx), %xmm0
    \\ movaps -2 * 8(%r9), %xmm1
    \\ movaps %xmm0, -2 * 8(%r9)
    \\ movaps %xmm1, -2 * 8(%rcx)
    \\ subq $-16 * 8, %r9
    \\ subq $-16 * 8, %rcx
    \\ decq %rax
    \\ jg 20018b
    \\ .balign 8
    \\ 20019:
    \\ testq $4, %rdi
    \\ jle 20020f
    \\ movaps -16 * 8(%rcx), %xmm0
    \\ movaps -16 * 8(%r9), %xmm1
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps %xmm1, -16 * 8(%rcx)
    \\ movaps -14 * 8(%rcx), %xmm0
    \\ movaps -14 * 8(%r9), %xmm1
    \\ movaps %xmm0, -14 * 8(%r9)
    \\ movaps %xmm1, -14 * 8(%rcx)
    \\ movaps -12 * 8(%rcx), %xmm0
    \\ movaps -12 * 8(%r9), %xmm1
    \\ movaps %xmm0, -12 * 8(%r9)
    \\ movaps %xmm1, -12 * 8(%rcx)
    \\ movaps -10 * 8(%rcx), %xmm0
    \\ movaps -10 * 8(%r9), %xmm1
    \\ movaps %xmm0, -10 * 8(%r9)
    \\ movaps %xmm1, -10 * 8(%rcx)
    \\ addq $8 * 8, %rcx
    \\ addq $8 * 8, %r9
    \\ .balign 8
    \\ 20020:
    \\ testq $2, %rdi
    \\ jle 20021f
    \\ movaps -16 * 8(%rcx), %xmm0
    \\ movaps -16 * 8(%r9), %xmm1
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps %xmm1, -16 * 8(%rcx)
    \\ movaps -14 * 8(%rcx), %xmm0
    \\ movaps -14 * 8(%r9), %xmm1
    \\ movaps %xmm0, -14 * 8(%r9)
    \\ movaps %xmm1, -14 * 8(%rcx)
    \\ addq $4 * 8, %rcx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20021:
    \\ testq $1, %rdi
    \\ jle 20022f
    \\ movaps -16 * 8(%rcx), %xmm0
    \\ movaps -16 * 8(%r9), %xmm1
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps %xmm1, -16 * 8(%rcx)
    \\ addq $2 * 8, %rcx
    \\ addq $2 * 8, %r9
    \\ .balign 8
    \\ 20022:
    \\ .if \bp4
    \\ movsd -16 * 8(%rcx), %xmm0
    \\ movsd -16 * 8(%r9), %xmm1
    \\ .endif
    \\ .if \bp3
    \\ movlps -16 * 8(%rcx), %xmm0
    \\ movlps -16 * 8(%r9), %xmm1
    \\ .endif
    \\ movlps %xmm0, -16 * 8(%r9)
    \\ movlps %xmm1, -16 * 8(%rcx)
    \\ xorq %rax,%rax
    \\ ret
    \\ .balign 8
    \\ 20023:
    \\ testq $8, %rcx
    \\ jne 20028f
    \\ testq $8, %r9
    \\ jne 20028f
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20025f
    \\ .balign 8
    \\ 20024:
    \\ movaps (%rcx), %xmm0
    \\ movaps (%r9), %xmm1
    \\ movaps %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movaps %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ movaps (%rcx), %xmm0
    \\ movaps (%r9), %xmm1
    \\ movaps %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movaps %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ movaps (%rcx), %xmm0
    \\ movaps (%r9), %xmm1
    \\ movaps %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movaps %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ movaps (%rcx), %xmm0
    \\ movaps (%r9), %xmm1
    \\ movaps %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movaps %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ decq %rax
    \\ jg 20024b
    \\ .balign 8
    \\ 20025:
    \\ movq %rdi, %rax
    \\ andq $3, %rax
    \\ jle 20027f
    \\ .balign 8
    \\ 20026:
    \\ movaps (%rcx), %xmm0
    \\ movaps (%r9), %xmm1
    \\ movaps %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movaps %xmm0, (%r9)
    \\ addq %rsi, %r9
    \\ decq %rax
    \\ jg 20026b
    \\ .balign 8
    \\ 20027:
    \\ xorq %rax, %rax
    \\ ret
    \\ .balign 8
    \\ 20028:
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20030f
    \\ .balign 8
    \\ 20029:
    \\ .if \bp4
    \\ movsd 0 * 8(%rcx), %xmm0
    \\ .endif
    \\ .if \bp3
    \\ movlps 0 * 8(%rcx), %xmm0
    \\ .endif
    \\ movhps 1 * 8(%rcx), %xmm0
    \\ .if \bp4
    \\ movsd 0 * 8(%r9), %xmm1
    \\ .endif
    \\ .if \bp3
    \\ movlps 0 * 8(%r9), %xmm1
    \\ .endif
    \\ movhps 1 * 8(%r9), %xmm1
    \\ movlps %xmm1, 0 * 8(%rcx)
    \\ movhps %xmm1, 1 * 8(%rcx)
    \\ addq %r8, %rcx
    \\ movlps %xmm0, 0 * 8(%r9)
    \\ movhps %xmm0, 1 * 8(%r9)
    \\ addq %rsi, %r9
    \\ .if \bp4
    \\ movsd 0 * 8(%rcx), %xmm0
    \\ .endif
    \\ .if \bp3
    \\ movlps 0 * 8(%rcx), %xmm0
    \\ .endif
    \\ movhps 1 * 8(%rcx), %xmm0
    \\ .if \bp4
    \\ movsd 0 * 8(%r9), %xmm1
    \\ .endif
    \\ .if \bp3
    \\ movlps 0 * 8(%r9), %xmm1
    \\ .endif
    \\ movhps 1 * 8(%r9), %xmm1
    \\ movlps %xmm1, 0 * 8(%rcx)
    \\ movhps %xmm1, 1 * 8(%rcx)
    \\ addq %r8, %rcx
    \\ movlps %xmm0, 0 * 8(%r9)
    \\ movhps %xmm0, 1 * 8(%r9)
    \\ addq %rsi, %r9
    \\ .if \bp4
    \\ movsd 0 * 8(%rcx), %xmm0
    \\ .endif
    \\ .if \bp3
    \\ movlps 0 * 8(%rcx), %xmm0
    \\ .endif
    \\ movhps 1 * 8(%rcx), %xmm0
    \\ .if \bp4
    \\ movsd 0 * 8(%r9), %xmm1
    \\ .endif
    \\ .if \bp3
    \\ movlps 0 * 8(%r9), %xmm1
    \\ .endif
    \\ movhps 1 * 8(%r9), %xmm1
    \\ movlps %xmm1, 0 * 8(%rcx)
    \\ movhps %xmm1, 1 * 8(%rcx)
    \\ addq %r8, %rcx
    \\ movlps %xmm0, 0 * 8(%r9)
    \\ movhps %xmm0, 1 * 8(%r9)
    \\ addq %rsi, %r9
    \\ .if \bp4
    \\ movsd 0 * 8(%rcx), %xmm0
    \\ .endif
    \\ .if \bp3
    \\ movlps 0 * 8(%rcx), %xmm0
    \\ .endif
    \\ movhps 1 * 8(%rcx), %xmm0
    \\ .if \bp4
    \\ movsd 0 * 8(%r9), %xmm1
    \\ .endif
    \\ .if \bp3
    \\ movlps 0 * 8(%r9), %xmm1
    \\ .endif
    \\ movhps 1 * 8(%r9), %xmm1
    \\ movlps %xmm1, 0 * 8(%rcx)
    \\ movhps %xmm1, 1 * 8(%rcx)
    \\ addq %r8, %rcx
    \\ movlps %xmm0, 0 * 8(%r9)
    \\ movhps %xmm0, 1 * 8(%r9)
    \\ addq %rsi, %r9
    \\ decq %rax
    \\ jg 20029b
    \\ .balign 8
    \\ 20030:
    \\ movq %rdi, %rax
    \\ andq $3, %rax
    \\ jle 20032f
    \\ .balign 8
    \\ 20031:
    \\ .if \bp4
    \\ movsd 0 * 8(%rcx), %xmm0
    \\ .endif
    \\ .if \bp3
    \\ movlps 0 * 8(%rcx), %xmm0
    \\ .endif
    \\ movhps 1 * 8(%rcx), %xmm0
    \\ .if \bp4
    \\ movsd 0 * 8(%r9), %xmm1
    \\ .endif
    \\ .if \bp3
    \\ movlps 0 * 8(%r9), %xmm1
    \\ .endif
    \\ movhps 1 * 8(%r9), %xmm1
    \\ movlps %xmm1, 0 * 8(%rcx)
    \\ movhps %xmm1, 1 * 8(%rcx)
    \\ addq %r8, %rcx
    \\ movlps %xmm0, 0 * 8(%r9)
    \\ movhps %xmm0, 1 * 8(%r9)
    \\ addq %rsi, %r9
    \\ decq %rax
    \\ jg 20031b
    \\ .balign 8
    \\ 20032:
    \\ xorq %rax, %rax
    \\ ret
    \\ .endm
    \\
;
