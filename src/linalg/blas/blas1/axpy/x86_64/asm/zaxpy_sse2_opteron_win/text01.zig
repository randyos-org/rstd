//! Section 1 of 3 of the text of `zaxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zaxpy_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1
    \\ .balign 512
    \\ movaps %xmm3, %xmm0
    \\ movlps 40(%rsp), %xmm1
    \\ movq 48(%rsp), %rdx
    \\ movq 56(%rsp), %r8
    \\ movq 64(%rsp), %r9
    \\ movq 72(%rsp), %r10
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
    \\ salq $4, %r8
    \\ salq $4, %r10
    \\ testq %rcx, %rcx
    \\ jle 20033f
    \\ cmpq $2 * 8, %r8
    \\ jne 20025f
    \\ cmpq $2 * 8, %r10
    \\ jne 20025f
    \\ subq $-16 * 8, %rdx
    \\ subq $-16 * 8, %r9
    \\ pcmpeqb %xmm7, %xmm7
    \\ psllq $63, %xmm7
    \\ .if \bp0
    \\ pshufd $0x44, %xmm0, %xmm14
    \\ pshufd $0x44, %xmm1, %xmm15
    \\ .endif
    \\ .if \bp1
    \\ movddup %xmm0, %xmm14
    \\ movddup %xmm1, %xmm15
    \\ .endif
    \\ shufps $0x0c, %xmm7, %xmm7
    \\ xorpd %xmm7, %xmm15
    \\ testq $8, %r9
    \\ jne 20011f
    \\ testq $8, %rdx
    \\ jne 20005f
    \\ movq %rcx, %rax
    \\ sarq $3, %rax
    \\ jle 20002f
    \\ movaps -16 * 8(%rdx), %xmm0
    \\ movaps -14 * 8(%rdx), %xmm1
    \\ movaps -12 * 8(%rdx), %xmm2
    \\ movaps -10 * 8(%rdx), %xmm3
    \\ decq %rax
    \\ jle 20001f
    \\ .balign 8
    \\ 20000:
    \\ movaps -8 * 8(%rdx), %xmm4
    \\ movaps -6 * 8(%rdx), %xmm5
    \\ prefetchw ((128 * 3) + 0) - 128(%r9)
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movaps -4 * 8(%rdx), %xmm6
    \\ movaps -2 * 8(%rdx), %xmm7
    \\ prefetch ((128 * 3) + 0) - 128(%rdx)
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps 2 * 8(%rdx), %xmm1
    \\ prefetchw ((128 * 3) + 64) - 128(%r9)
    \\ pshufd $0x4e, %xmm4, %xmm8
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ addpd -8 * 8(%r9), %xmm4
    \\ movaps %xmm4, -8 * 8(%r9)
    \\ pshufd $0x4e, %xmm5, %xmm8
    \\ mulpd %xmm14, %xmm5
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm5
    \\ addpd -6 * 8(%r9), %xmm5
    \\ movaps %xmm5, -6 * 8(%r9)
    \\ movaps 4 * 8(%rdx), %xmm2
    \\ movaps 6 * 8(%rdx), %xmm3
    \\ prefetch ((128 * 3) + 64) - 128(%rdx)
    \\ pshufd $0x4e, %xmm6, %xmm8
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm6
    \\ addpd -4 * 8(%r9), %xmm6
    \\ movaps %xmm6, -4 * 8(%r9)
    \\ pshufd $0x4e, %xmm7, %xmm8
    \\ mulpd %xmm14, %xmm7
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm7
    \\ addpd -2 * 8(%r9), %xmm7
    \\ movaps %xmm7, -2 * 8(%r9)
    \\ subq $-16 * 8, %rdx
    \\ subq $-16 * 8, %r9
    \\ decq %rax
    \\ jg 20000b
    \\ .balign 8
    \\ 20001:
    \\ movaps -8 * 8(%rdx), %xmm4
    \\ movaps -6 * 8(%rdx), %xmm5
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movaps -4 * 8(%rdx), %xmm6
    \\ movaps -2 * 8(%rdx), %xmm7
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ pshufd $0x4e, %xmm4, %xmm8
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ addpd -8 * 8(%r9), %xmm4
    \\ movaps %xmm4, -8 * 8(%r9)
    \\ pshufd $0x4e, %xmm5, %xmm8
    \\ mulpd %xmm14, %xmm5
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm5
    \\ addpd -6 * 8(%r9), %xmm5
    \\ movaps %xmm5, -6 * 8(%r9)
    \\ pshufd $0x4e, %xmm6, %xmm8
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm6
    \\ addpd -4 * 8(%r9), %xmm6
    \\ movaps %xmm6, -4 * 8(%r9)
    \\ pshufd $0x4e, %xmm7, %xmm8
    \\ mulpd %xmm14, %xmm7
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm7
    \\ addpd -2 * 8(%r9), %xmm7
    \\ movaps %xmm7, -2 * 8(%r9)
    \\ subq $-16 * 8, %rdx
    \\ subq $-16 * 8, %r9
    \\ .balign 8
    \\ 20002:
    \\ movq %rcx, %rax
    \\ andq $4, %rax
    \\ jle 20003f
    \\ movaps -16 * 8(%rdx), %xmm0
    \\ movaps -14 * 8(%rdx), %xmm1
    \\ movaps -12 * 8(%rdx), %xmm2
    \\ movaps -10 * 8(%rdx), %xmm3
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ addq $8 * 8, %rdx
    \\ addq $8 * 8, %r9
    \\ .balign 8
    \\ 20003:
    \\ movq %rcx, %rax
    \\ andq $2, %rax
    \\ jle 20004f
    \\ movaps -16 * 8(%rdx), %xmm0
    \\ movaps -14 * 8(%rdx), %xmm1
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ addq $4 * 8, %rdx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20004:
    \\ movq %rcx, %rax
    \\ andq $1, %rax
    \\ jle 20033f
    \\ movaps -16 * 8(%rdx), %xmm0
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ jmp 20033f
    \\ .balign 8
    \\ 20005:
    \\ movq %rcx, %rax
    \\ sarq $3, %rax
    \\ jle 20008f
    \\ movlps -16 * 8(%rdx), %xmm0
    \\ movhps -15 * 8(%rdx), %xmm0
    \\ movlps -14 * 8(%rdx), %xmm1
    \\ movhps -13 * 8(%rdx), %xmm1
    \\ movlps -12 * 8(%rdx), %xmm2
    \\ movhps -11 * 8(%rdx), %xmm2
    \\ movlps -10 * 8(%rdx), %xmm3
    \\ movhps -9 * 8(%rdx), %xmm3
    \\ decq %rax
    \\ jle 20007f
    \\ .balign 8
    \\ 20006:
    \\ movlps -8 * 8(%rdx), %xmm4
    \\ movhps -7 * 8(%rdx), %xmm4
    \\ prefetchw ((128 * 3) + 0) - 128(%r9)
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movlps -6 * 8(%rdx), %xmm5
    \\ movhps -5 * 8(%rdx), %xmm5
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movlps -4 * 8(%rdx), %xmm6
    \\ movhps -3 * 8(%rdx), %xmm6
    \\ prefetch ((128 * 3) + 0) - 128(%rdx)
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ movlps -2 * 8(%rdx), %xmm7
    \\ movhps -1 * 8(%rdx), %xmm7
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ movlps 0 * 8(%rdx), %xmm0
    \\ movhps 1 * 8(%rdx), %xmm0
    \\ prefetchw ((128 * 3) + 64) - 128(%r9)
    \\ pshufd $0x4e, %xmm4, %xmm8
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ addpd -8 * 8(%r9), %xmm4
    \\ movaps %xmm4, -8 * 8(%r9)
    \\ movlps 2 * 8(%rdx), %xmm1
    \\ movhps 3 * 8(%rdx), %xmm1
    \\ pshufd $0x4e, %xmm5, %xmm8
    \\ mulpd %xmm14, %xmm5
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm5
    \\ addpd -6 * 8(%r9), %xmm5
    \\ movaps %xmm5, -6 * 8(%r9)
    \\ movlps 4 * 8(%rdx), %xmm2
    \\ movhps 5 * 8(%rdx), %xmm2
    \\ prefetch ((128 * 3) + 64) - 128(%rdx)
    \\ pshufd $0x4e, %xmm6, %xmm8
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm6
    \\ addpd -4 * 8(%r9), %xmm6
    \\ movaps %xmm6, -4 * 8(%r9)
    \\ movlps 6 * 8(%rdx), %xmm3
    \\ movhps 7 * 8(%rdx), %xmm3
    \\ pshufd $0x4e, %xmm7, %xmm8
    \\ mulpd %xmm14, %xmm7
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm7
    \\ addpd -2 * 8(%r9), %xmm7
    \\ movaps %xmm7, -2 * 8(%r9)
    \\ subq $-16 * 8, %rdx
    \\ subq $-16 * 8, %r9
    \\ decq %rax
    \\ jg 20006b
    \\ .balign 8
    \\ 20007:
    \\ movlps -8 * 8(%rdx), %xmm4
    \\ movhps -7 * 8(%rdx), %xmm4
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movlps -6 * 8(%rdx), %xmm5
    \\ movhps -5 * 8(%rdx), %xmm5
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movlps -4 * 8(%rdx), %xmm6
    \\ movhps -3 * 8(%rdx), %xmm6
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ movlps -2 * 8(%rdx), %xmm7
    \\ movhps -1 * 8(%rdx), %xmm7
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ pshufd $0x4e, %xmm4, %xmm8
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ addpd -8 * 8(%r9), %xmm4
    \\ movaps %xmm4, -8 * 8(%r9)
    \\ pshufd $0x4e, %xmm5, %xmm8
    \\ mulpd %xmm14, %xmm5
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm5
    \\ addpd -6 * 8(%r9), %xmm5
    \\ movaps %xmm5, -6 * 8(%r9)
    \\ pshufd $0x4e, %xmm6, %xmm8
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm6
    \\ addpd -4 * 8(%r9), %xmm6
    \\ movaps %xmm6, -4 * 8(%r9)
    \\ pshufd $0x4e, %xmm7, %xmm8
    \\ mulpd %xmm14, %xmm7
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm7
    \\ addpd -2 * 8(%r9), %xmm7
    \\ movaps %xmm7, -2 * 8(%r9)
    \\ subq $-16 * 8, %rdx
    \\ subq $-16 * 8, %r9
    \\ .balign 8
    \\ 20008:
    \\ movq %rcx, %rax
    \\ andq $4, %rax
    \\ jle 20009f
    \\ movlps -16 * 8(%rdx), %xmm0
    \\ movhps -15 * 8(%rdx), %xmm0
    \\ movlps -14 * 8(%rdx), %xmm1
    \\ movhps -13 * 8(%rdx), %xmm1
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movlps -12 * 8(%rdx), %xmm2
    \\ movhps -11 * 8(%rdx), %xmm2
    \\ movlps -10 * 8(%rdx), %xmm3
    \\ movhps -9 * 8(%rdx), %xmm3
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ addq $8 * 8, %rdx
    \\ addq $8 * 8, %r9
    \\ .balign 8
    \\ 20009:
    \\ movq %rcx, %rax
    \\ andq $2, %rax
    \\ jle 20010f
    \\ movlps -16 * 8(%rdx), %xmm0
    \\ movhps -15 * 8(%rdx), %xmm0
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movlps -14 * 8(%rdx), %xmm1
    \\ movhps -13 * 8(%rdx), %xmm1
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ addq $4 * 8, %rdx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20010:
    \\ movq %rcx, %rax
    \\ andq $1, %rax
    \\ jle 20033f
    \\ movlps -16 * 8(%rdx), %xmm0
    \\ movhps -15 * 8(%rdx), %xmm0
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ jmp 20033f
    \\
;
