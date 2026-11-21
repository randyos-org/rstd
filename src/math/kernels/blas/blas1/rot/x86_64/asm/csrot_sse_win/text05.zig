//! Section 5 of 5 of the text of `csrot`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zrot_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20017:
    \\ testq $15, %rcx
    \\ jle 20025f
    \\ testq $8, %rcx
    \\ jle 20018f
    \\ .if \bp0
    \\ movsd 0 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%r9), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%r9), %xmm1
    \\ .if \bp0
    \\ movsd 0 * 4(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rdx), %xmm0
    \\ .endif
    \\ movhps 2 * 4(%rdx), %xmm0
    \\ .if \bp0
    \\ movsd 4 * 4(%r9), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%r9), %xmm3
    \\ .endif
    \\ movhps 6 * 4(%r9), %xmm3
    \\ .if \bp0
    \\ movsd 4 * 4(%rdx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rdx), %xmm2
    \\ .endif
    \\ movhps 6 * 4(%rdx), %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm5
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm7
    \\ addps %xmm1, %xmm0
    \\ addps %xmm3, %xmm2
    \\ subps %xmm5, %xmm4
    \\ subps %xmm7, %xmm6
    \\ movlps %xmm0, 0 * 4(%rdx)
    \\ movhps %xmm0, 2 * 4(%rdx)
    \\ movlps %xmm2, 4 * 4(%rdx)
    \\ movhps %xmm2, 6 * 4(%rdx)
    \\ movlps %xmm4, 0 * 4(%r9)
    \\ movhps %xmm4, 2 * 4(%r9)
    \\ movlps %xmm6, 4 * 4(%r9)
    \\ movhps %xmm6, 6 * 4(%r9)
    \\ .if \bp0
    \\ movsd 8 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%r9), %xmm1
    \\ .endif
    \\ movhps 10 * 4(%r9), %xmm1
    \\ .if \bp0
    \\ movsd 8 * 4(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rdx), %xmm0
    \\ .endif
    \\ movhps 10 * 4(%rdx), %xmm0
    \\ .if \bp0
    \\ movsd 12 * 4(%r9), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%r9), %xmm3
    \\ .endif
    \\ movhps 14 * 4(%r9), %xmm3
    \\ .if \bp0
    \\ movsd 12 * 4(%rdx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%rdx), %xmm2
    \\ .endif
    \\ movhps 14 * 4(%rdx), %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm5
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm7
    \\ addps %xmm1, %xmm0
    \\ addps %xmm3, %xmm2
    \\ subps %xmm5, %xmm4
    \\ subps %xmm7, %xmm6
    \\ movlps %xmm0, 8 * 4(%rdx)
    \\ movhps %xmm0, 10 * 4(%rdx)
    \\ movlps %xmm2, 12 * 4(%rdx)
    \\ movhps %xmm2, 14 * 4(%rdx)
    \\ movlps %xmm4, 8 * 4(%r9)
    \\ movhps %xmm4, 10 * 4(%r9)
    \\ movlps %xmm6, 12 * 4(%r9)
    \\ movhps %xmm6, 14 * 4(%r9)
    \\ addq $16 * 4, %rdx
    \\ addq $16 * 4, %r9
    \\ .balign 8
    \\ 20018:
    \\ testq $4, %rcx
    \\ jle 20019f
    \\ .if \bp0
    \\ movsd 0 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%r9), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%r9), %xmm1
    \\ .if \bp0
    \\ movsd 0 * 4(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rdx), %xmm0
    \\ .endif
    \\ movhps 2 * 4(%rdx), %xmm0
    \\ .if \bp0
    \\ movsd 4 * 4(%r9), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%r9), %xmm3
    \\ .endif
    \\ movhps 6 * 4(%r9), %xmm3
    \\ .if \bp0
    \\ movsd 4 * 4(%rdx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rdx), %xmm2
    \\ .endif
    \\ movhps 6 * 4(%rdx), %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm5
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm7
    \\ addps %xmm1, %xmm0
    \\ addps %xmm3, %xmm2
    \\ subps %xmm5, %xmm4
    \\ subps %xmm7, %xmm6
    \\ movlps %xmm0, 0 * 4(%rdx)
    \\ movhps %xmm0, 2 * 4(%rdx)
    \\ movlps %xmm2, 4 * 4(%rdx)
    \\ movhps %xmm2, 6 * 4(%rdx)
    \\ movlps %xmm4, 0 * 4(%r9)
    \\ movhps %xmm4, 2 * 4(%r9)
    \\ movlps %xmm6, 4 * 4(%r9)
    \\ movhps %xmm6, 6 * 4(%r9)
    \\ addq $8 * 4, %rdx
    \\ addq $8 * 4, %r9
    \\ .balign 8
    \\ 20019:
    \\ testq $2, %rcx
    \\ jle 20020f
    \\ .if \bp0
    \\ movsd 0 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%r9), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%r9), %xmm1
    \\ .if \bp0
    \\ movsd 0 * 4(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rdx), %xmm0
    \\ .endif
    \\ movhps 2 * 4(%rdx), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movlps %xmm0, 0 * 4(%rdx)
    \\ movhps %xmm0, 2 * 4(%rdx)
    \\ movlps %xmm2, 0 * 4(%r9)
    \\ movhps %xmm2, 2 * 4(%r9)
    \\ addq $4 * 4, %rdx
    \\ addq $4 * 4, %r9
    \\ .balign 8
    \\ 20020:
    \\ testq $1, %rcx
    \\ jle 20025f
    \\ .if \bp0
    \\ movsd 0 * 4(%r9), %xmm1
    \\ movsd 0 * 4(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%r9), %xmm1
    \\ movlps 0 * 4(%rdx), %xmm0
    \\ .endif
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movlps %xmm0, 0 * 4(%rdx)
    \\ movlps %xmm2, 0 * 4(%r9)
    \\ jmp 20025f
    \\ .balign 8
    \\ .balign 8
    \\ 20021:
    \\ movq %rcx, %rax
    \\ cmpq $0, %r8
    \\ je 20024f
    \\ cmpq $0, %r10
    \\ je 20024f
    \\ sarq $2, %rax
    \\ jle 20023f
    \\ .balign 8
    \\ 20022:
    \\ .if \bp0
    \\ movsd (%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps (%r9), %xmm1
    \\ .endif
    \\ movhps (%r9, %r10), %xmm1
    \\ .if \bp0
    \\ movsd (%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rdx), %xmm0
    \\ .endif
    \\ movhps (%rdx, %r8), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movlps %xmm0, (%rdx)
    \\ movhps %xmm0, (%rdx, %r8)
    \\ movlps %xmm2, (%r9)
    \\ movhps %xmm2, (%r9, %r10)
    \\ leaq (%rdx, %r8, 2), %rdx
    \\ leaq (%r9, %r10, 2), %r9
    \\ .if \bp0
    \\ movsd (%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps (%r9), %xmm1
    \\ .endif
    \\ movhps (%r9, %r10), %xmm1
    \\ .if \bp0
    \\ movsd (%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rdx), %xmm0
    \\ .endif
    \\ movhps (%rdx, %r8), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movlps %xmm0, (%rdx)
    \\ movhps %xmm0, (%rdx, %r8)
    \\ movlps %xmm2, (%r9)
    \\ movhps %xmm2, (%r9, %r10)
    \\ leaq (%rdx, %r8, 2), %rdx
    \\ leaq (%r9, %r10, 2), %r9
    \\ decq %rax
    \\ jg 20022b
    \\ .balign 8
    \\ 20023:
    \\ movq %rcx, %rax
    \\ andq $3, %rax
    \\ jle 20025f
    \\ .balign 8
    \\ 20024:
    \\ .if \bp0
    \\ movsd (%r9), %xmm1
    \\ movsd (%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%r9), %xmm1
    \\ movlps (%rdx), %xmm0
    \\ .endif
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movlps %xmm0, (%rdx)
    \\ movlps %xmm2, (%r9)
    \\ addq %r8, %rdx
    \\ addq %r10, %r9
    \\ decq %rax
    \\ jg 20024b
    \\ .balign 8
    \\ 20025:
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
    \\ .endm
    \\
;
