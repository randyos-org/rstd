//! Section 1 of 3 of the text of `sswap`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/swap_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1, bp2, bp3, bp4
    \\ .balign 512
    \\ pushq %rbx
    \\ movq 48(%rsp), %rdx
    \\ movq 56(%rsp), %r8
    \\ movq 64(%rsp), %r9
    \\ movq 72(%rsp), %rbx
    \\ subq $256, %rsp
    \\ movups %xmm6, 0(%rsp)
    \\ movups %xmm7, 16(%rsp)
    \\ movups %xmm8, 32(%rsp)
    \\ movups %xmm9, 48(%rsp)
    \\ movups %xmm10, 64(%rsp)
    \\ movups %xmm11, 80(%rsp)
    \\ movups %xmm12, 96(%rsp)
    \\ movups %xmm13, 112(%rsp)
    \\ movups %xmm14, 128(%rsp)
    \\ movups %xmm15, 144(%rsp)
    \\ leaq (, %r8, 4), %r8
    \\ leaq (, %rbx, 4), %rbx
    \\ cmpq $4, %r8
    \\ jne 20033f
    \\ cmpq $4, %rbx
    \\ jne 20033f
    \\ subq $-32 * 4, %rdx
    \\ subq $-32 * 4, %r9
    \\ cmpq $3, %rcx
    \\ jle 20006f
    \\ testq $4, %r9
    \\ je 20000f
    \\ movss -32 * 4(%rdx), %xmm0
    \\ movss -32 * 4(%r9), %xmm1
    \\ movss %xmm1, -32 * 4(%rdx)
    \\ movss %xmm0, -32 * 4(%r9)
    \\ addq $1 * 4, %rdx
    \\ addq $1 * 4, %r9
    \\ decq %rcx
    \\ .balign 8
    \\ 20000:
    \\ testq $2 * 4, %r9
    \\ je 20001f
    \\ .if \bp0
    \\ movsd -32 * 4(%rdx), %xmm0
    \\ movsd -32 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rdx), %xmm0
    \\ movlps -32 * 4(%r9), %xmm1
    \\ .endif
    \\ movlps %xmm1, -32 * 4(%rdx)
    \\ movlps %xmm0, -32 * 4(%r9)
    \\ addq $2 * 4, %rdx
    \\ addq $2 * 4, %r9
    \\ subq $2, %rcx
    \\ jle 20008f
    \\ .balign 8
    \\ 20001:
    \\ cmpq $3, %rcx
    \\ jle 20006f
    \\ testq $2 * 4, %rdx
    \\ jne 20017f
    \\ testq $1 * 4, %rdx
    \\ jne 20009f
    \\ movq %rcx, %rax
    \\ sarq $5, %rax
    \\ jle 20003f
    \\ .balign 8
    \\ 20002:
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
    \\ movaps -32 * 4(%rdx), %xmm0
    \\ movaps -32 * 4(%r9), %xmm1
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, -32 * 4(%rdx)
    \\ movaps -28 * 4(%rdx), %xmm0
    \\ movaps -28 * 4(%r9), %xmm1
    \\ movaps %xmm0, -28 * 4(%r9)
    \\ movaps %xmm1, -28 * 4(%rdx)
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
    \\ movaps -24 * 4(%rdx), %xmm0
    \\ movaps -24 * 4(%r9), %xmm1
    \\ movaps %xmm0, -24 * 4(%r9)
    \\ movaps %xmm1, -24 * 4(%rdx)
    \\ movaps -20 * 4(%rdx), %xmm0
    \\ movaps -20 * 4(%r9), %xmm1
    \\ movaps %xmm0, -20 * 4(%r9)
    \\ movaps %xmm1, -20 * 4(%rdx)
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
    \\ movaps -16 * 4(%rdx), %xmm0
    \\ movaps -16 * 4(%r9), %xmm1
    \\ movaps %xmm0, -16 * 4(%r9)
    \\ movaps %xmm1, -16 * 4(%rdx)
    \\ movaps -12 * 4(%rdx), %xmm0
    \\ movaps -12 * 4(%r9), %xmm1
    \\ movaps %xmm0, -12 * 4(%r9)
    \\ movaps %xmm1, -12 * 4(%rdx)
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
    \\ movaps -8 * 4(%rdx), %xmm0
    \\ movaps -8 * 4(%r9), %xmm1
    \\ movaps %xmm0, -8 * 4(%r9)
    \\ movaps %xmm1, -8 * 4(%rdx)
    \\ movaps -4 * 4(%rdx), %xmm0
    \\ movaps -4 * 4(%r9), %xmm1
    \\ movaps %xmm0, -4 * 4(%r9)
    \\ movaps %xmm1, -4 * 4(%rdx)
    \\ subq $-32 * 4, %r9
    \\ subq $-32 * 4, %rdx
    \\ decq %rax
    \\ jg 20002b
    \\ .balign 8
    \\ 20003:
    \\ testq $16, %rcx
    \\ jle 20004f
    \\ movaps -32 * 4(%rdx), %xmm0
    \\ movaps -32 * 4(%r9), %xmm1
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, -32 * 4(%rdx)
    \\ movaps -28 * 4(%rdx), %xmm0
    \\ movaps -28 * 4(%r9), %xmm1
    \\ movaps %xmm0, -28 * 4(%r9)
    \\ movaps %xmm1, -28 * 4(%rdx)
    \\ movaps -24 * 4(%rdx), %xmm0
    \\ movaps -24 * 4(%r9), %xmm1
    \\ movaps %xmm0, -24 * 4(%r9)
    \\ movaps %xmm1, -24 * 4(%rdx)
    \\ movaps -20 * 4(%rdx), %xmm0
    \\ movaps -20 * 4(%r9), %xmm1
    \\ movaps %xmm0, -20 * 4(%r9)
    \\ movaps %xmm1, -20 * 4(%rdx)
    \\ addq $16 * 4, %rdx
    \\ addq $16 * 4, %r9
    \\ .balign 8
    \\ 20004:
    \\ testq $8, %rcx
    \\ jle 20005f
    \\ movaps -32 * 4(%rdx), %xmm0
    \\ movaps -32 * 4(%r9), %xmm1
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, -32 * 4(%rdx)
    \\ movaps -28 * 4(%rdx), %xmm0
    \\ movaps -28 * 4(%r9), %xmm1
    \\ movaps %xmm0, -28 * 4(%r9)
    \\ movaps %xmm1, -28 * 4(%rdx)
    \\ addq $8 * 4, %rdx
    \\ addq $8 * 4, %r9
    \\ .balign 8
    \\ 20005:
    \\ testq $4, %rcx
    \\ jle 20006f
    \\ movaps -32 * 4(%rdx), %xmm0
    \\ movaps -32 * 4(%r9), %xmm1
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, -32 * 4(%rdx)
    \\ addq $4 * 4, %rdx
    \\ addq $4 * 4, %r9
    \\ .balign 8
    \\ 20006:
    \\ testq $2, %rcx
    \\ jle 20007f
    \\ .if \bp0
    \\ movsd -32 * 4(%rdx), %xmm0
    \\ movsd -32 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rdx), %xmm0
    \\ movlps -32 * 4(%r9), %xmm1
    \\ .endif
    \\ movlps %xmm1, -32 * 4(%rdx)
    \\ addq $2 * 4, %rdx
    \\ movlps %xmm0, -32 * 4(%r9)
    \\ addq $2 * 4, %r9
    \\ .balign 8
    \\ 20007:
    \\ testq $1, %rcx
    \\ jle 20008f
    \\ movss -32 * 4(%rdx), %xmm0
    \\ movss -32 * 4(%r9), %xmm1
    \\ movss %xmm1, -32 * 4(%rdx)
    \\ movss %xmm0, -32 * 4(%r9)
    \\ .balign 8
    \\ 20008:
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
    \\ popq %rbx
    \\ ret
    \\ .balign 8
    \\ 20009:
    \\ movaps -33 * 4(%rdx), %xmm0
    \\ movaps -32 * 4(%r9), %xmm1
    \\ movss %xmm1, -32 * 4(%rdx)
    \\ pshufd $0x39, %xmm1, %xmm3
    \\ movlps %xmm3, -31 * 4(%rdx)
    \\ subq $3, %rcx
    \\ movq %rcx, %rax
    \\ sarq $5, %rax
    \\ jle 20011f
    \\ .balign 16
    \\ 20010:
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
    \\ movaps -29 * 4(%rdx), %xmm2
    \\ movaps -28 * 4(%r9), %xmm3
    \\ movss %xmm2, %xmm0
    \\ shufps $0x39, %xmm0, %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm3, %xmm1
    \\ shufps $0x93, %xmm3, %xmm1
    \\ movaps %xmm1, -29 * 4(%rdx)
    \\ movaps -25 * 4(%rdx), %xmm0
    \\ movaps -24 * 4(%r9), %xmm1
    \\ movss %xmm0, %xmm2
    \\ shufps $0x39, %xmm2, %xmm2
    \\ movaps %xmm2, -28 * 4(%r9)
    \\ movss %xmm1, %xmm3
    \\ shufps $0x93, %xmm1, %xmm3
    \\ movaps %xmm3, -25 * 4(%rdx)
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
    \\ movaps -21 * 4(%rdx), %xmm2
    \\ movaps -20 * 4(%r9), %xmm3
    \\ movss %xmm2, %xmm0
    \\ shufps $0x39, %xmm0, %xmm0
    \\ movaps %xmm0, -24 * 4(%r9)
    \\ movss %xmm3, %xmm1
    \\ shufps $0x93, %xmm3, %xmm1
    \\ movaps %xmm1, -21 * 4(%rdx)
    \\ movaps -17 * 4(%rdx), %xmm0
    \\ movaps -16 * 4(%r9), %xmm1
    \\ movss %xmm0, %xmm2
    \\ shufps $0x39, %xmm2, %xmm2
    \\ movaps %xmm2, -20 * 4(%r9)
    \\ movss %xmm1, %xmm3
    \\ shufps $0x93, %xmm1, %xmm3
    \\ movaps %xmm3, -17 * 4(%rdx)
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
    \\ movaps -13 * 4(%rdx), %xmm2
    \\ movaps -12 * 4(%r9), %xmm3
    \\ movss %xmm2, %xmm0
    \\ shufps $0x39, %xmm0, %xmm0
    \\ movaps %xmm0, -16 * 4(%r9)
    \\ movss %xmm3, %xmm1
    \\ shufps $0x93, %xmm3, %xmm1
    \\ movaps %xmm1, -13 * 4(%rdx)
    \\ movaps -9 * 4(%rdx), %xmm0
    \\ movaps -8 * 4(%r9), %xmm1
    \\ movss %xmm0, %xmm2
    \\ shufps $0x39, %xmm2, %xmm2
    \\ movaps %xmm2, -12 * 4(%r9)
    \\ movss %xmm1, %xmm3
    \\ shufps $0x93, %xmm1, %xmm3
    \\ movaps %xmm3, -9 * 4(%rdx)
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
    \\ movaps -5 * 4(%rdx), %xmm2
    \\ movaps -4 * 4(%r9), %xmm3
    \\ movss %xmm2, %xmm0
    \\ shufps $0x39, %xmm0, %xmm0
    \\ movaps %xmm0, -8 * 4(%r9)
    \\ movss %xmm3, %xmm1
    \\ shufps $0x93, %xmm3, %xmm1
    \\ movaps %xmm1, -5 * 4(%rdx)
    \\ movaps -1 * 4(%rdx), %xmm0
    \\ movaps 0 * 4(%r9), %xmm1
    \\ movss %xmm0, %xmm2
    \\ shufps $0x39, %xmm2, %xmm2
    \\ movaps %xmm2, -4 * 4(%r9)
    \\ movss %xmm1, %xmm3
    \\ shufps $0x93, %xmm1, %xmm3
    \\ movaps %xmm3, -1 * 4(%rdx)
    \\ subq $-32 * 4, %rdx
    \\ subq $-32 * 4, %r9
    \\ decq %rax
    \\ jg 20010b
    \\
;
