//! Section 1 of 2 of the text of `daxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/axpy_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1
    \\ .balign 512
    \\ movq 40(%rsp), %rdx
    \\ movq 48(%rsp), %r8
    \\ movq 56(%rsp), %r9
    \\ movq 64(%rsp), %r10
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
    \\ movaps %xmm3, %xmm15
    \\ unpcklpd %xmm15, %xmm15
    \\ leaq (, %r8, 8), %r8
    \\ leaq (, %r10, 8), %r10
    \\ testq %rcx, %rcx
    \\ jle 20020f
    \\ cmpq $8, %r8
    \\ jne 20016f
    \\ cmpq $8, %r10
    \\ jne 20016f
    \\ testq $8, %r9
    \\ je 20000f
    \\ .if \bp0
    \\ movsd (%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rdx), %xmm0
    \\ .endif
    \\ mulsd %xmm15, %xmm0
    \\ addsd (%r9), %xmm0
    \\ .if \bp0
    \\ movsd %xmm0, (%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm0, (%r9)
    \\ .endif
    \\ addq $1 * 8, %rdx
    \\ addq $1 * 8, %r9
    \\ decq %rcx
    \\ jle 20007f
    \\ .balign 16
    \\ 20000:
    \\ subq $-16 * 8, %rdx
    \\ subq $-16 * 8, %r9
    \\ testq $8, %rdx
    \\ jne 20008f
    \\ movq %rcx, %rax
    \\ sarq $4, %rax
    \\ jle 20003f
    \\ movaps -16 * 8(%rdx), %xmm0
    \\ movaps -14 * 8(%rdx), %xmm1
    \\ movaps -12 * 8(%rdx), %xmm2
    \\ movaps -10 * 8(%rdx), %xmm3
    \\ decq %rax
    \\ jle 20002f
    \\ .balign 8
    \\ 20001:
    \\ movaps -8 * 8(%rdx), %xmm4
    \\ movaps -6 * 8(%rdx), %xmm5
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%r9)
    \\ .endif
    \\ mulpd %xmm15, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ mulpd %xmm15, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movaps -4 * 8(%rdx), %xmm6
    \\ movaps -2 * 8(%rdx), %xmm7
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 0) - 128(%rdx)
    \\ .endif
    \\ mulpd %xmm15, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ mulpd %xmm15, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps 2 * 8(%rdx), %xmm1
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%r9)
    \\ .endif
    \\ mulpd %xmm15, %xmm4
    \\ addpd -8 * 8(%r9), %xmm4
    \\ movaps %xmm4, -8 * 8(%r9)
    \\ mulpd %xmm15, %xmm5
    \\ addpd -6 * 8(%r9), %xmm5
    \\ movaps %xmm5, -6 * 8(%r9)
    \\ movaps 4 * 8(%rdx), %xmm2
    \\ movaps 6 * 8(%rdx), %xmm3
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 64) - 128(%rdx)
    \\ .endif
    \\ mulpd %xmm15, %xmm6
    \\ addpd -4 * 8(%r9), %xmm6
    \\ movaps %xmm6, -4 * 8(%r9)
    \\ mulpd %xmm15, %xmm7
    \\ addpd -2 * 8(%r9), %xmm7
    \\ movaps %xmm7, -2 * 8(%r9)
    \\ subq $-16 * 8, %r9
    \\ subq $-16 * 8, %rdx
    \\ decq %rax
    \\ jg 20001b
    \\ .balign 8
    \\ 20002:
    \\ movaps -8 * 8(%rdx), %xmm4
    \\ movaps -6 * 8(%rdx), %xmm5
    \\ mulpd %xmm15, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ mulpd %xmm15, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movaps -4 * 8(%rdx), %xmm6
    \\ movaps -2 * 8(%rdx), %xmm7
    \\ mulpd %xmm15, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ mulpd %xmm15, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ mulpd %xmm15, %xmm4
    \\ addpd -8 * 8(%r9), %xmm4
    \\ movaps %xmm4, -8 * 8(%r9)
    \\ mulpd %xmm15, %xmm5
    \\ addpd -6 * 8(%r9), %xmm5
    \\ movaps %xmm5, -6 * 8(%r9)
    \\ mulpd %xmm15, %xmm6
    \\ addpd -4 * 8(%r9), %xmm6
    \\ movaps %xmm6, -4 * 8(%r9)
    \\ mulpd %xmm15, %xmm7
    \\ addpd -2 * 8(%r9), %xmm7
    \\ movaps %xmm7, -2 * 8(%r9)
    \\ subq $-16 * 8, %r9
    \\ subq $-16 * 8, %rdx
    \\ .balign 8
    \\ 20003:
    \\ movq %rcx, %rax
    \\ andq $8, %rax
    \\ jle 20004f
    \\ .balign 8
    \\ movaps -16 * 8(%rdx), %xmm0
    \\ movaps -14 * 8(%rdx), %xmm1
    \\ movaps -12 * 8(%rdx), %xmm2
    \\ movaps -10 * 8(%rdx), %xmm3
    \\ mulpd %xmm15, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ mulpd %xmm15, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ addq $8 * 8, %rdx
    \\ addq $8 * 8, %r9
    \\ .balign 8
    \\ 20004:
    \\ movq %rcx, %rax
    \\ andq $4, %rax
    \\ jle 20005f
    \\ .balign 8
    \\ movaps -16 * 8(%rdx), %xmm0
    \\ movaps -14 * 8(%rdx), %xmm1
    \\ mulpd %xmm15, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ addpd -16 * 8(%r9), %xmm0
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ addq $4 * 8, %rdx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20005:
    \\ movq %rcx, %rax
    \\ andq $2, %rax
    \\ jle 20006f
    \\ .balign 8
    \\ movaps -16 * 8(%rdx), %xmm0
    \\ mulpd %xmm15, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ addq $2 * 8, %rdx
    \\ addq $2 * 8, %r9
    \\ .balign 8
    \\ 20006:
    \\ movq %rcx, %rax
    \\ andq $1, %rax
    \\ jle 20007f
    \\ .balign 8
    \\ .if \bp0
    \\ movsd -16 * 8(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps -16 * 8(%rdx), %xmm0
    \\ .endif
    \\ mulsd %xmm15, %xmm0
    \\ addsd -16 * 8(%r9), %xmm0
    \\ .if \bp0
    \\ movsd %xmm0, -16 * 8(%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm0, -16 * 8(%r9)
    \\ .endif
    \\ .balign 8
    \\ 20007:
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
    \\ 20008:
    \\ movq %rcx, %rax
    \\ sarq $4, %rax
    \\ jle 20011f
    \\ .if \bp0
    \\ movsd -16 * 8(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps -16 * 8(%rdx), %xmm0
    \\ .endif
    \\ movhps -15 * 8(%rdx), %xmm0
    \\ .if \bp0
    \\ movsd -14 * 8(%rdx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps -14 * 8(%rdx), %xmm1
    \\ .endif
    \\ movhps -13 * 8(%rdx), %xmm1
    \\ .if \bp0
    \\ movsd -12 * 8(%rdx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps -12 * 8(%rdx), %xmm2
    \\ .endif
    \\ movhps -11 * 8(%rdx), %xmm2
    \\ .if \bp0
    \\ movsd -10 * 8(%rdx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps -10 * 8(%rdx), %xmm3
    \\ .endif
    \\ movhps -9 * 8(%rdx), %xmm3
    \\ decq %rax
    \\ jle 20010f
    \\ .balign 8
    \\ 20009:
    \\ .if \bp0
    \\ movsd -8 * 8(%rdx), %xmm4
    \\ .endif
    \\ .if \bp1
    \\ movlps -8 * 8(%rdx), %xmm4
    \\ .endif
    \\ movhps -7 * 8(%rdx), %xmm4
    \\ .if \bp0
    \\ movsd -6 * 8(%rdx), %xmm5
    \\ .endif
    \\ .if \bp1
    \\ movlps -6 * 8(%rdx), %xmm5
    \\ .endif
    \\ movhps -5 * 8(%rdx), %xmm5
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%r9)
    \\ .endif
    \\ mulpd %xmm15, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ mulpd %xmm15, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ .if \bp0
    \\ movsd -4 * 8(%rdx), %xmm6
    \\ .endif
    \\ .if \bp1
    \\ movlps -4 * 8(%rdx), %xmm6
    \\ .endif
    \\ movhps -3 * 8(%rdx), %xmm6
    \\ .if \bp0
    \\ movsd -2 * 8(%rdx), %xmm7
    \\ .endif
    \\ .if \bp1
    \\ movlps -2 * 8(%rdx), %xmm7
    \\ .endif
    \\ movhps -1 * 8(%rdx), %xmm7
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 0) - 128(%rdx)
    \\ .endif
    \\ mulpd %xmm15, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ mulpd %xmm15, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ .if \bp0
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 8(%rdx), %xmm0
    \\ .endif
    \\ movhps 1 * 8(%rdx), %xmm0
    \\ .if \bp0
    \\ movsd 2 * 8(%rdx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 2 * 8(%rdx), %xmm1
    \\ .endif
    \\ movhps 3 * 8(%rdx), %xmm1
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%r9)
    \\ .endif
    \\ mulpd %xmm15, %xmm4
    \\ addpd -8 * 8(%r9), %xmm4
    \\ movaps %xmm4, -8 * 8(%r9)
    \\ mulpd %xmm15, %xmm5
    \\ addpd -6 * 8(%r9), %xmm5
    \\ movaps %xmm5, -6 * 8(%r9)
    \\ .if \bp0
    \\ movsd 4 * 8(%rdx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 8(%rdx), %xmm2
    \\ .endif
    \\ movhps 5 * 8(%rdx), %xmm2
    \\ .if \bp0
    \\ movsd 6 * 8(%rdx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 6 * 8(%rdx), %xmm3
    \\ .endif
    \\ movhps 7 * 8(%rdx), %xmm3
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 64) - 128(%rdx)
    \\ .endif
    \\ mulpd %xmm15, %xmm6
    \\ addpd -4 * 8(%r9), %xmm6
    \\ movaps %xmm6, -4 * 8(%r9)
    \\ mulpd %xmm15, %xmm7
    \\ addpd -2 * 8(%r9), %xmm7
    \\ movaps %xmm7, -2 * 8(%r9)
    \\ subq $-16 * 8, %r9
    \\ subq $-16 * 8, %rdx
    \\ decq %rax
    \\ jg 20009b
    \\
;
