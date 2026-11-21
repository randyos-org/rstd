//! Section 2 of 3 of the text of `ccopy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zcopy_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20013:
    \\ testq $16, %rcx
    \\ jle 20014f
    \\ .balign 8
    \\ movaps -30 * 4(%rdx), %xmm1
    \\ movaps -26 * 4(%rdx), %xmm2
    \\ movaps -22 * 4(%rdx), %xmm3
    \\ movaps -18 * 4(%rdx), %xmm4
    \\ shufps $0x4e, %xmm1, %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ shufps $0x4e, %xmm2, %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ shufps $0x4e, %xmm3, %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ shufps $0x4e, %xmm4, %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ movaps %xmm4, %xmm0
    \\ addq $16 * 4, %rdx
    \\ addq $16 * 4, %r9
    \\ .balign 8
    \\ 20014:
    \\ testq $8, %rcx
    \\ jle 20015f
    \\ .balign 8
    \\ movaps -30 * 4(%rdx), %xmm1
    \\ movaps -26 * 4(%rdx), %xmm2
    \\ shufps $0x4e, %xmm1, %xmm0
    \\ shufps $0x4e, %xmm2, %xmm1
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps %xmm2, %xmm0
    \\ addq $8 * 4, %rdx
    \\ addq $8 * 4, %r9
    \\ .balign 8
    \\ 20015:
    \\ testq $4, %rcx
    \\ jle 20016f
    \\ .balign 8
    \\ movaps -30 * 4(%rdx), %xmm1
    \\ shufps $0x4e, %xmm1, %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ addq $4 * 4, %rdx
    \\ addq $4 * 4, %r9
    \\ .balign 8
    \\ 20016:
    \\ testq $2, %rcx
    \\ jle 20017f
    \\ .balign 8
    \\ .if \bp0
    \\ movsd -32 * 4(%rdx), %xmm0
    \\ movsd %xmm0, -32 * 4(%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rdx), %xmm0
    \\ movlps %xmm0, -32 * 4(%r9)
    \\ .endif
    \\ addq $2 * 4, %rdx
    \\ addq $2 * 4, %r9
    \\ .balign 8
    \\ 20017:
    \\ testq $1, %rcx
    \\ jle 20018f
    \\ .balign 8
    \\ movss -32 * 4(%rdx), %xmm0
    \\ movss %xmm0, -32 * 4(%r9)
    \\ addq $4, %r9
    \\ .balign 8
    \\ 20018:
    \\ xorq %rax,%rax
    \\ movups 0(%rsp), %xmm6
    \\ movups 16(%rsp), %xmm7
    \\ movups 32(%rsp), %xmm8
    \\ movups 48(%rsp), %xmm9
    \\ movups 64(%rsp), %xmm10
    \\ movups 80(%rsp), %xmm11
    \\ movups 96(%rsp), %xmm12
    \\ movups 112(%rsp), %xmm13
    \\ movups 128(%rsp), %xmm14
    \\ movups 144(%rsp), %xmm15
    \\ addq $256, %rsp
    \\ ret
    \\ .balign 8
    \\ 20019:
    \\ testq $2 * 4, %rdx
    \\ jne 20028f
    \\ movaps -33 * 4(%rdx), %xmm0
    \\ movq %rcx, %rax
    \\ sarq $5, %rax
    \\ jle 20022f
    \\ movaps -29 * 4(%rdx), %xmm1
    \\ movaps -25 * 4(%rdx), %xmm2
    \\ movaps -21 * 4(%rdx), %xmm3
    \\ movaps -17 * 4(%rdx), %xmm4
    \\ movaps -13 * 4(%rdx), %xmm5
    \\ movaps -9 * 4(%rdx), %xmm6
    \\ movaps -5 * 4(%rdx), %xmm7
    \\ decq %rax
    \\ jle 20021f
    \\ .balign 16
    \\ 20020:
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
    \\ movss %xmm1, %xmm0
    \\ shufps $0x39, %xmm0, %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps -1 * 4(%rdx), %xmm0
    \\ movss %xmm2, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps 3 * 4(%rdx), %xmm1
    \\ .if \bp5
    \\ prefetcht0 (( 64 * 12 + 32) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp6
    \\ prefetcht0 ((128 * 10) + 0) - 128(%rdx)
    \\ .endif
    \\ movss %xmm3, %xmm2
    \\ shufps $0x39, %xmm2, %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ movaps 7 * 4(%rdx), %xmm2
    \\ movss %xmm4, %xmm3
    \\ shufps $0x39, %xmm3, %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ movaps 11 * 4(%rdx), %xmm3
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
    \\ movss %xmm5, %xmm4
    \\ shufps $0x39, %xmm4, %xmm4
    \\ movaps %xmm4, -16 * 4(%r9)
    \\ movaps 15 * 4(%rdx), %xmm4
    \\ movss %xmm6, %xmm5
    \\ shufps $0x39, %xmm5, %xmm5
    \\ movaps %xmm5, -12 * 4(%r9)
    \\ movaps 19 * 4(%rdx), %xmm5
    \\ .if \bp5
    \\ prefetcht0 (( 64 * 12 + 32) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 64) - 128(%rdx)
    \\ .endif
    \\ movss %xmm7, %xmm6
    \\ shufps $0x39, %xmm6, %xmm6
    \\ movaps %xmm6, -8 * 4(%r9)
    \\ movaps 23 * 4(%rdx), %xmm6
    \\ movss %xmm0, %xmm7
    \\ shufps $0x39, %xmm7, %xmm7
    \\ movaps %xmm7, -4 * 4(%r9)
    \\ movaps 27 * 4(%rdx), %xmm7
    \\ subq $-32 * 4, %rdx
    \\ subq $-32 * 4, %r9
    \\ decq %rax
    \\ jg 20020b
    \\ .balign 8
    \\ 20021:
    \\ movss %xmm1, %xmm0
    \\ shufps $0x39, %xmm0, %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps -1 * 4(%rdx), %xmm0
    \\ movss %xmm2, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movss %xmm3, %xmm2
    \\ shufps $0x39, %xmm2, %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ movss %xmm4, %xmm3
    \\ shufps $0x39, %xmm3, %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ movss %xmm5, %xmm4
    \\ shufps $0x39, %xmm4, %xmm4
    \\ movaps %xmm4, -16 * 4(%r9)
    \\ movss %xmm6, %xmm5
    \\ shufps $0x39, %xmm5, %xmm5
    \\ movaps %xmm5, -12 * 4(%r9)
    \\ movss %xmm7, %xmm6
    \\ shufps $0x39, %xmm6, %xmm6
    \\ movaps %xmm6, -8 * 4(%r9)
    \\ movss %xmm0, %xmm7
    \\ shufps $0x39, %xmm7, %xmm7
    \\ movaps %xmm7, -4 * 4(%r9)
    \\ subq $-32 * 4, %rdx
    \\ subq $-32 * 4, %r9
    \\ .balign 8
    \\ 20022:
    \\ testq $16, %rcx
    \\ jle 20023f
    \\ .balign 8
    \\ movaps -29 * 4(%rdx), %xmm1
    \\ movaps -25 * 4(%rdx), %xmm2
    \\ movaps -21 * 4(%rdx), %xmm3
    \\ movaps -17 * 4(%rdx), %xmm4
    \\ movss %xmm1, %xmm0
    \\ shufps $0x39, %xmm0, %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm2, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movss %xmm3, %xmm2
    \\ shufps $0x39, %xmm2, %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ movss %xmm4, %xmm3
    \\ shufps $0x39, %xmm3, %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ movaps %xmm4, %xmm0
    \\ addq $16 * 4, %rdx
    \\ addq $16 * 4, %r9
    \\ .balign 8
    \\ 20023:
    \\ testq $8, %rcx
    \\ jle 20024f
    \\ .balign 8
    \\ movaps -29 * 4(%rdx), %xmm1
    \\ movaps -25 * 4(%rdx), %xmm2
    \\ movss %xmm1, %xmm0
    \\ shufps $0x39, %xmm0, %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm2, %xmm1
    \\ shufps $0x39, %xmm1, %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps %xmm2, %xmm0
    \\ addq $8 * 4, %rdx
    \\ addq $8 * 4, %r9
    \\ .balign 8
    \\ 20024:
    \\ testq $4, %rcx
    \\ jle 20025f
    \\ .balign 8
    \\ movaps -29 * 4(%rdx), %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x39, %xmm0, %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ addq $4 * 4, %rdx
    \\ addq $4 * 4, %r9
    \\ .balign 8
    \\ 20025:
    \\ testq $2, %rcx
    \\ jle 20026f
    \\ .balign 8
    \\ .if \bp0
    \\ movsd -32 * 4(%rdx), %xmm0
    \\ movsd %xmm0, -32 * 4(%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rdx), %xmm0
    \\ movlps %xmm0, -32 * 4(%r9)
    \\ .endif
    \\ addq $2 * 4, %rdx
    \\ addq $2 * 4, %r9
    \\ .balign 8
    \\ 20026:
    \\ testq $1, %rcx
    \\ jle 20027f
    \\ .balign 8
    \\ movss -32 * 4(%rdx), %xmm0
    \\ movss %xmm0, -32 * 4(%r9)
    \\ addq $4, %r9
    \\ .balign 8
    \\ 20027:
    \\ xorq %rax,%rax
    \\ movups 0(%rsp), %xmm6
    \\ movups 16(%rsp), %xmm7
    \\ movups 32(%rsp), %xmm8
    \\ movups 48(%rsp), %xmm9
    \\ movups 64(%rsp), %xmm10
    \\ movups 80(%rsp), %xmm11
    \\ movups 96(%rsp), %xmm12
    \\ movups 112(%rsp), %xmm13
    \\ movups 128(%rsp), %xmm14
    \\ movups 144(%rsp), %xmm15
    \\ addq $256, %rsp
    \\ ret
    \\ .balign 8
    \\ 20028:
    \\ movaps -35 * 4(%rdx), %xmm0
    \\ movq %rcx, %rax
    \\ sarq $5, %rax
    \\ jle 20031f
    \\ movaps -31 * 4(%rdx), %xmm1
    \\ movaps -27 * 4(%rdx), %xmm2
    \\ movaps -23 * 4(%rdx), %xmm3
    \\ movaps -19 * 4(%rdx), %xmm4
    \\ movaps -15 * 4(%rdx), %xmm5
    \\ movaps -11 * 4(%rdx), %xmm6
    \\ movaps -7 * 4(%rdx), %xmm7
    \\ decq %rax
    \\ jle 20030f
    \\
;
