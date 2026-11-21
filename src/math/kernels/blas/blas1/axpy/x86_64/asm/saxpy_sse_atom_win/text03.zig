//! Section 3 of 3 of the text of `saxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/axpy_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
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
    \\ decq %rax
    \\ jle 20030f
    \\ .balign 16
    \\ 20029:
    \\ movaps -15 * 4(%rdx), %xmm5
    \\ movaps -11 * 4(%rdx), %xmm6
    \\ .if \bp0
    \\ prefetchw ((128 * 5) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 ((128 * 4) + 0) - 128(%r9)
    \\ .endif
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ mulps %xmm15, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps -7 * 4(%rdx), %xmm7
    \\ movaps -3 * 4(%rdx), %xmm0
    \\ .if \bp2
    \\ prefetcht0 (( 64 * 12 + 32) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp0
    \\ prefetch ((128 * 5) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 10) + 0) - 128(%rdx)
    \\ .endif
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ mulps %xmm15, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ movss %xmm4, %xmm3
    \\ shufps $0x93, %xmm4, %xmm3
    \\ mulps %xmm15, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ movaps 1 * 4(%rdx), %xmm1
    \\ movaps 5 * 4(%rdx), %xmm2
    \\ .if \bp0
    \\ prefetchw ((128 * 5) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 ((128 * 4) + 64) - 128(%r9)
    \\ .endif
    \\ movss %xmm5, %xmm4
    \\ shufps $0x93, %xmm5, %xmm4
    \\ mulps %xmm15, %xmm4
    \\ addps -16 * 4(%r9), %xmm4
    \\ movaps %xmm4, -16 * 4(%r9)
    \\ movss %xmm6, %xmm5
    \\ shufps $0x93, %xmm6, %xmm5
    \\ mulps %xmm15, %xmm5
    \\ addps -12 * 4(%r9), %xmm5
    \\ movaps %xmm5, -12 * 4(%r9)
    \\ movaps 9 * 4(%rdx), %xmm3
    \\ movaps 13 * 4(%rdx), %xmm4
    \\ .if \bp2
    \\ prefetcht0 (( 64 * 12 + 32) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp0
    \\ prefetch ((128 * 5) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rdx)
    \\ .endif
    \\ movss %xmm7, %xmm6
    \\ shufps $0x93, %xmm7, %xmm6
    \\ mulps %xmm15, %xmm6
    \\ addps -8 * 4(%r9), %xmm6
    \\ movaps %xmm6, -8 * 4(%r9)
    \\ movss %xmm0, %xmm7
    \\ shufps $0x93, %xmm0, %xmm7
    \\ mulps %xmm15, %xmm7
    \\ addps -4 * 4(%r9), %xmm7
    \\ movaps %xmm7, -4 * 4(%r9)
    \\ subq $-32 * 4, %rdx
    \\ subq $-32 * 4, %r9
    \\ decq %rax
    \\ jg 20029b
    \\ .balign 8
    \\ 20030:
    \\ movaps -15 * 4(%rdx), %xmm5
    \\ movaps -11 * 4(%rdx), %xmm6
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ mulps %xmm15, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps -7 * 4(%rdx), %xmm7
    \\ movaps -3 * 4(%rdx), %xmm0
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ mulps %xmm15, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ movss %xmm4, %xmm3
    \\ shufps $0x93, %xmm4, %xmm3
    \\ mulps %xmm15, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ movss %xmm5, %xmm4
    \\ shufps $0x93, %xmm5, %xmm4
    \\ mulps %xmm15, %xmm4
    \\ addps -16 * 4(%r9), %xmm4
    \\ movaps %xmm4, -16 * 4(%r9)
    \\ movss %xmm6, %xmm5
    \\ shufps $0x93, %xmm6, %xmm5
    \\ mulps %xmm15, %xmm5
    \\ addps -12 * 4(%r9), %xmm5
    \\ movaps %xmm5, -12 * 4(%r9)
    \\ movss %xmm7, %xmm6
    \\ shufps $0x93, %xmm7, %xmm6
    \\ mulps %xmm15, %xmm6
    \\ addps -8 * 4(%r9), %xmm6
    \\ movaps %xmm6, -8 * 4(%r9)
    \\ movss %xmm0, %xmm7
    \\ shufps $0x93, %xmm0, %xmm7
    \\ mulps %xmm15, %xmm7
    \\ addps -4 * 4(%r9), %xmm7
    \\ movaps %xmm7, -4 * 4(%r9)
    \\ subq $-32 * 4, %rdx
    \\ subq $-32 * 4, %r9
    \\ .balign 8
    \\ 20031:
    \\ movq %rcx, %rax
    \\ andq $16, %rax
    \\ jle 20032f
    \\ .balign 8
    \\ movaps -31 * 4(%rdx), %xmm1
    \\ movaps -27 * 4(%rdx), %xmm2
    \\ movaps -23 * 4(%rdx), %xmm3
    \\ movaps -19 * 4(%rdx), %xmm4
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ mulps %xmm15, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ mulps %xmm15, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movss %xmm4, %xmm3
    \\ shufps $0x93, %xmm4, %xmm3
    \\ mulps %xmm15, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ movaps %xmm4, %xmm0
    \\ addq $16 * 4, %rdx
    \\ addq $16 * 4, %r9
    \\ .balign 8
    \\ 20032:
    \\ movq %rcx, %rax
    \\ andq $8, %rax
    \\ jle 20033f
    \\ .balign 8
    \\ movaps -31 * 4(%rdx), %xmm1
    \\ movaps -27 * 4(%rdx), %xmm2
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ mulps %xmm15, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps %xmm2, %xmm0
    \\ addq $8 * 4, %rdx
    \\ addq $8 * 4, %r9
    \\ .balign 8
    \\ 20033:
    \\ movq %rcx, %rax
    \\ andq $4, %rax
    \\ jle 20034f
    \\ .balign 8
    \\ movaps -31 * 4(%rdx), %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ addq $4 * 4, %rdx
    \\ addq $4 * 4, %r9
    \\ .balign 8
    \\ 20034:
    \\ movq %rcx, %rax
    \\ andq $2, %rax
    \\ jle 20035f
    \\ .balign 8
    \\ movsd -32 * 4(%rdx), %xmm0
    \\ movsd -32 * 4(%r9), %xmm4
    \\ mulps %xmm15, %xmm0
    \\ addps %xmm4, %xmm0
    \\ movsd %xmm0, -32 * 4(%r9)
    \\ addq $2 * 4, %rdx
    \\ addq $2 * 4, %r9
    \\ .balign 8
    \\ 20035:
    \\ movq %rcx, %rax
    \\ andq $1, %rax
    \\ jle 20036f
    \\ .balign 8
    \\ movss -32 * 4(%rdx), %xmm0
    \\ mulss %xmm15, %xmm0
    \\ addss -32 * 4(%r9), %xmm0
    \\ movss %xmm0, -32 * 4(%r9)
    \\ addq $4, %r9
    \\ .balign 8
    \\ 20036:
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
    \\ 20037:
    \\ movq %rcx, %rax
    \\ movq %r9, %r11
    \\ cmpq $0, %r8
    \\ je 20040f
    \\ cmpq $0, %r10
    \\ je 20040f
    \\ sarq $3, %rax
    \\ jle 20039f
    \\ .balign 8
    \\ 20038:
    \\ movss (%rdx), %xmm0
    \\ addq %r8, %rdx
    \\ mulss %xmm15, %xmm0
    \\ movss (%r11), %xmm6
    \\ addq %r10, %r11
    \\ addss %xmm6, %xmm0
    \\ movss (%rdx), %xmm1
    \\ addq %r8, %rdx
    \\ mulss %xmm15, %xmm1
    \\ movss (%r11), %xmm6
    \\ addq %r10, %r11
    \\ addss %xmm6, %xmm1
    \\ movss (%rdx), %xmm2
    \\ addq %r8, %rdx
    \\ mulss %xmm15, %xmm2
    \\ movss (%r11), %xmm6
    \\ addq %r10, %r11
    \\ addss %xmm6, %xmm2
    \\ movss (%rdx), %xmm3
    \\ addq %r8, %rdx
    \\ mulss %xmm15, %xmm3
    \\ movss (%r11), %xmm6
    \\ addq %r10, %r11
    \\ addss %xmm6, %xmm3
    \\ movss %xmm0, (%r9)
    \\ addq %r10, %r9
    \\ movss %xmm1, (%r9)
    \\ addq %r10, %r9
    \\ movss %xmm2, (%r9)
    \\ addq %r10, %r9
    \\ movss %xmm3, (%r9)
    \\ addq %r10, %r9
    \\ movss (%rdx), %xmm0
    \\ addq %r8, %rdx
    \\ mulss %xmm15, %xmm0
    \\ movss (%r11), %xmm6
    \\ addq %r10, %r11
    \\ addss %xmm6, %xmm0
    \\ movss (%rdx), %xmm1
    \\ addq %r8, %rdx
    \\ mulss %xmm15, %xmm1
    \\ movss (%r11), %xmm6
    \\ addq %r10, %r11
    \\ addss %xmm6, %xmm1
    \\ movss (%rdx), %xmm2
    \\ addq %r8, %rdx
    \\ mulss %xmm15, %xmm2
    \\ movss (%r11), %xmm6
    \\ addq %r10, %r11
    \\ addss %xmm6, %xmm2
    \\ movss (%rdx), %xmm3
    \\ addq %r8, %rdx
    \\ mulss %xmm15, %xmm3
    \\ movss (%r11), %xmm6
    \\ addq %r10, %r11
    \\ addss %xmm6, %xmm3
    \\ movss %xmm0, (%r9)
    \\ addq %r10, %r9
    \\ movss %xmm1, (%r9)
    \\ addq %r10, %r9
    \\ movss %xmm2, (%r9)
    \\ addq %r10, %r9
    \\ movss %xmm3, (%r9)
    \\ addq %r10, %r9
    \\ decq %rax
    \\ jg 20038b
    \\ .balign 8
    \\ 20039:
    \\ movq %rcx, %rax
    \\ andq $7, %rax
    \\ jle 20041f
    \\ .balign 8
    \\ 20040:
    \\ movss (%rdx), %xmm0
    \\ addq %r8, %rdx
    \\ mulss %xmm15, %xmm0
    \\ movss (%r9), %xmm6
    \\ addss %xmm6, %xmm0
    \\ movss %xmm0, (%r9)
    \\ addq %r10, %r9
    \\ decq %rax
    \\ jg 20040b
    \\ .balign 8
    \\ 20041:
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
    \\ .endm
    \\
;
