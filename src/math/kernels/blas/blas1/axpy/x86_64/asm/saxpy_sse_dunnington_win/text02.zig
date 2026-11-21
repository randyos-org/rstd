//! Section 2 of 2 of the text of `saxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/axpy_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20012:
    \\ .if \bp0
    \\ movsd -16 * 4(%rdx), %xmm4
    \\ .endif
    \\ .if \bp1
    \\ movlps -16 * 4(%rdx), %xmm4
    \\ .endif
    \\ movhps -14 * 4(%rdx), %xmm4
    \\ .if \bp0
    \\ movsd -12 * 4(%rdx), %xmm5
    \\ .endif
    \\ .if \bp1
    \\ movlps -12 * 4(%rdx), %xmm5
    \\ .endif
    \\ movhps -10 * 4(%rdx), %xmm5
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ mulps %xmm15, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ .if \bp0
    \\ movsd -8 * 4(%rdx), %xmm6
    \\ .endif
    \\ .if \bp1
    \\ movlps -8 * 4(%rdx), %xmm6
    \\ .endif
    \\ movhps -6 * 4(%rdx), %xmm6
    \\ .if \bp0
    \\ movsd -4 * 4(%rdx), %xmm7
    \\ .endif
    \\ .if \bp1
    \\ movlps -4 * 4(%rdx), %xmm7
    \\ .endif
    \\ movhps -2 * 4(%rdx), %xmm7
    \\ mulps %xmm15, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ mulps %xmm15, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ mulps %xmm15, %xmm4
    \\ addps -16 * 4(%r9), %xmm4
    \\ movaps %xmm4, -16 * 4(%r9)
    \\ mulps %xmm15, %xmm5
    \\ addps -12 * 4(%r9), %xmm5
    \\ movaps %xmm5, -12 * 4(%r9)
    \\ mulps %xmm15, %xmm6
    \\ addps -8 * 4(%r9), %xmm6
    \\ movaps %xmm6, -8 * 4(%r9)
    \\ mulps %xmm15, %xmm7
    \\ addps -4 * 4(%r9), %xmm7
    \\ movaps %xmm7, -4 * 4(%r9)
    \\ subq $-32 * 4, %rdx
    \\ subq $-32 * 4, %r9
    \\ .balign 8
    \\ 20013:
    \\ movq %rcx, %rax
    \\ andq $16, %rax
    \\ jle 20014f
    \\ .balign 8
    \\ .if \bp0
    \\ movsd -32 * 4(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rdx), %xmm0
    \\ .endif
    \\ movhps -30 * 4(%rdx), %xmm0
    \\ .if \bp0
    \\ movsd -28 * 4(%rdx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps -28 * 4(%rdx), %xmm1
    \\ .endif
    \\ movhps -26 * 4(%rdx), %xmm1
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ mulps %xmm15, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ .if \bp0
    \\ movsd -24 * 4(%rdx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps -24 * 4(%rdx), %xmm2
    \\ .endif
    \\ movhps -22 * 4(%rdx), %xmm2
    \\ .if \bp0
    \\ movsd -20 * 4(%rdx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps -20 * 4(%rdx), %xmm3
    \\ .endif
    \\ movhps -18 * 4(%rdx), %xmm3
    \\ mulps %xmm15, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ mulps %xmm15, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ addq $16 * 4, %rdx
    \\ addq $16 * 4, %r9
    \\ .balign 8
    \\ 20014:
    \\ movq %rcx, %rax
    \\ andq $8, %rax
    \\ jle 20015f
    \\ .balign 8
    \\ .if \bp0
    \\ movsd -32 * 4(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rdx), %xmm0
    \\ .endif
    \\ movhps -30 * 4(%rdx), %xmm0
    \\ .if \bp0
    \\ movsd -28 * 4(%rdx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps -28 * 4(%rdx), %xmm1
    \\ .endif
    \\ movhps -26 * 4(%rdx), %xmm1
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ mulps %xmm15, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ addq $8 * 4, %rdx
    \\ addq $8 * 4, %r9
    \\ .balign 8
    \\ 20015:
    \\ movq %rcx, %rax
    \\ andq $4, %rax
    \\ jle 20016f
    \\ .balign 8
    \\ .if \bp0
    \\ movsd -32 * 4(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rdx), %xmm0
    \\ .endif
    \\ movhps -30 * 4(%rdx), %xmm0
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ addq $4 * 4, %rdx
    \\ addq $4 * 4, %r9
    \\ .balign 8
    \\ 20016:
    \\ movq %rcx, %rax
    \\ andq $2, %rax
    \\ jle 20017f
    \\ .balign 8
    \\ .if \bp0
    \\ movsd -32 * 4(%rdx), %xmm0
    \\ movsd -32 * 4(%r9), %xmm4
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rdx), %xmm0
    \\ movlps -32 * 4(%r9), %xmm4
    \\ .endif
    \\ mulps %xmm15, %xmm0
    \\ addps %xmm4, %xmm0
    \\ .if \bp0
    \\ movsd %xmm0, -32 * 4(%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm0, -32 * 4(%r9)
    \\ .endif
    \\ addq $2 * 4, %rdx
    \\ addq $2 * 4, %r9
    \\ .balign 8
    \\ 20017:
    \\ movq %rcx, %rax
    \\ andq $1, %rax
    \\ jle 20018f
    \\ .balign 8
    \\ movss -32 * 4(%rdx), %xmm0
    \\ mulss %xmm15, %xmm0
    \\ addss -32 * 4(%r9), %xmm0
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
    \\ movq %rcx, %rax
    \\ movq %r9, %r11
    \\ cmpq $0, %r8
    \\ je 20022f
    \\ cmpq $0, %r10
    \\ je 20022f
    \\ sarq $3, %rax
    \\ jle 20021f
    \\ .balign 8
    \\ 20020:
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
    \\ jg 20020b
    \\ .balign 8
    \\ 20021:
    \\ movq %rcx, %rax
    \\ andq $7, %rax
    \\ jle 20023f
    \\ .balign 8
    \\ 20022:
    \\ movss (%rdx), %xmm0
    \\ addq %r8, %rdx
    \\ mulss %xmm15, %xmm0
    \\ movss (%r9), %xmm6
    \\ addss %xmm6, %xmm0
    \\ movss %xmm0, (%r9)
    \\ addq %r10, %r9
    \\ decq %rax
    \\ jg 20022b
    \\ .balign 8
    \\ 20023:
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
