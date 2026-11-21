//! Section 2 of 3 of the text of `zgemv_r`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zgemv_n_dup.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20005:
    \\ testq $2, %rdi
    \\ je 20006f
    \\ movddup -16 * 8(%r12), %xmm4
    \\ movddup -15 * 8(%r12), %xmm5
    \\ movddup -14 * 8(%r12), %xmm6
    \\ movddup -13 * 8(%r12), %xmm7
    \\ mulpd %xmm8, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -16 * 8(%r12, %r8, 1), %xmm4
    \\ mulpd %xmm8, %xmm6
    \\ addpd %xmm6, %xmm1
    \\ movddup -14 * 8(%r12, %r8, 1), %xmm6
    \\ mulpd %xmm9, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ movddup -15 * 8(%r12, %r8, 1), %xmm5
    \\ mulpd %xmm9, %xmm7
    \\ addpd %xmm7, %xmm1
    \\ movddup -13 * 8(%r12, %r8, 1), %xmm7
    \\ mulpd %xmm10, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -16 * 8(%r13), %xmm4
    \\ mulpd %xmm10, %xmm6
    \\ addpd %xmm6, %xmm1
    \\ movddup -14 * 8(%r13), %xmm6
    \\ mulpd %xmm11, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ movddup -15 * 8(%r13), %xmm5
    \\ mulpd %xmm11, %xmm7
    \\ addpd %xmm7, %xmm1
    \\ movddup -13 * 8(%r13), %xmm7
    \\ mulpd %xmm12, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -16 * 8(%r13, %r8, 1), %xmm4
    \\ mulpd %xmm12, %xmm6
    \\ addpd %xmm6, %xmm1
    \\ movddup -14 * 8(%r13, %r8, 1), %xmm6
    \\ mulpd %xmm13, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ movddup -15 * 8(%r13, %r8, 1), %xmm5
    \\ mulpd %xmm13, %xmm7
    \\ addpd %xmm7, %xmm1
    \\ movddup -13 * 8(%r13, %r8, 1), %xmm7
    \\ mulpd %xmm14, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ mulpd %xmm14, %xmm6
    \\ addpd %xmm6, %xmm1
    \\ mulpd %xmm15, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm7, %xmm1
    \\ movlps %xmm0, -16 * 8(%r14)
    \\ movhps %xmm0, -16 * 8 + 8(%r14)
    \\ movlps %xmm1, -14 * 8(%r14)
    \\ movhps %xmm1, -14 * 8 + 8(%r14)
    \\ movaps %xmm2, %xmm0
    \\ addq $4 * 8, %r12
    \\ addq $4 * 8, %r13
    \\ addq $4 * 8, %r14
    \\ .balign 8
    \\ 20006:
    \\ testq $1, %rdi
    \\ je 20007f
    \\ movddup -16 * 8(%r12), %xmm4
    \\ movddup -15 * 8(%r12), %xmm5
    \\ movddup -16 * 8(%r12, %r8, 1), %xmm6
    \\ movddup -15 * 8(%r12, %r8, 1), %xmm7
    \\ mulpd %xmm8, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -16 * 8(%r13), %xmm4
    \\ mulpd %xmm9, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ movddup -15 * 8(%r13), %xmm5
    \\ mulpd %xmm10, %xmm6
    \\ addpd %xmm6, %xmm0
    \\ movddup -16 * 8(%r13, %r8, 1), %xmm6
    \\ mulpd %xmm11, %xmm7
    \\ addpd %xmm7, %xmm0
    \\ movddup -15 * 8(%r13, %r8, 1), %xmm7
    \\ mulpd %xmm12, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ mulpd %xmm13, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ mulpd %xmm14, %xmm6
    \\ addpd %xmm6, %xmm0
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm7, %xmm0
    \\ movlps %xmm0, -16 * 8(%r14)
    \\ movhps %xmm0, -16 * 8 + 8(%r14)
    \\ .balign 8
    \\ 20007:
    \\ cmpq $4, %rsi
    \\ jge 20002b
    \\ .balign 8
    \\ 20008:
    \\ cmpq $2, %rsi
    \\ jl 20013f
    \\ subq $2, %rsi
    \\ leaq 16 * 8(%r15), %r14
    \\ movq %rcx, %r12
    \\ leaq (%rcx, %r8, 1), %r13
    \\ leaq (%rcx, %r8, 2), %rcx
    \\ movddup 0 * 8(%r9), %xmm8
    \\ movddup 1 * 8(%r9), %xmm9
    \\ addq %rdx, %r9
    \\ movddup 0 * 8(%r9), %xmm10
    \\ movddup 1 * 8(%r9), %xmm11
    \\ addq %rdx, %r9
    \\ pcmpeqb %xmm5, %xmm5
    \\ psllq $63, %xmm5
    \\ shufps $0x40, %xmm5, %xmm5
    \\ movlps 48 (%rsp), %xmm6
    \\ movhps 56 (%rsp), %xmm6
    \\ pshufd $0x4e, %xmm6, %xmm7
    \\ xorps %xmm5, %xmm7
    \\ mulpd %xmm6, %xmm8
    \\ mulpd %xmm7, %xmm9
    \\ mulpd %xmm6, %xmm10
    \\ mulpd %xmm7, %xmm11
    \\ subpd %xmm9, %xmm8
    \\ subpd %xmm11, %xmm10
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ xorps %xmm5, %xmm9
    \\ xorps %xmm5, %xmm11
    \\ movlps -16 * 8(%r14), %xmm0
    \\ movhps -16 * 8 + 8(%r14), %xmm0
    \\ movlps -14 * 8(%r14), %xmm1
    \\ movhps -14 * 8 + 8(%r14), %xmm1
    \\ movlps -12 * 8(%r14), %xmm2
    \\ movhps -12 * 8 + 8(%r14), %xmm2
    \\ movlps -10 * 8(%r14), %xmm3
    \\ movhps -10 * 8 + 8(%r14), %xmm3
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20011f
    \\ movddup -16 * 8(%r12), %xmm4
    \\ movddup -14 * 8(%r12), %xmm5
    \\ movddup -12 * 8(%r12), %xmm6
    \\ movddup -10 * 8(%r12), %xmm7
    \\ decq %rax
    \\ jle 20010f
    \\ .balign 8
    \\ 20009:
    \\ prefetch (64 * 1) * 2 - 128 + (64 >> 1)(%r13)
    \\ mulpd %xmm8, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -15 * 8(%r12), %xmm4
    \\ mulpd %xmm8, %xmm5
    \\ addpd %xmm5, %xmm1
    \\ movddup -13 * 8(%r12), %xmm5
    \\ mulpd %xmm8, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ movddup -11 * 8(%r12), %xmm6
    \\ mulpd %xmm8, %xmm7
    \\ addpd %xmm7, %xmm3
    \\ movddup -9 * 8(%r12), %xmm7
    \\ mulpd %xmm9, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -16 * 8(%r13), %xmm4
    \\ mulpd %xmm9, %xmm5
    \\ addpd %xmm5, %xmm1
    \\ movddup -14 * 8(%r13), %xmm5
    \\ mulpd %xmm9, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ movddup -12 * 8(%r13), %xmm6
    \\ mulpd %xmm9, %xmm7
    \\ addpd %xmm7, %xmm3
    \\ movddup -10 * 8(%r13), %xmm7
    \\ prefetch (64 * 1) * 2 - 128 + (64 >> 1)(%r12)
    \\ mulpd %xmm10, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -15 * 8(%r13), %xmm4
    \\ mulpd %xmm10, %xmm5
    \\ addpd %xmm5, %xmm1
    \\ movddup -13 * 8(%r13), %xmm5
    \\ mulpd %xmm10, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ movddup -11 * 8(%r13), %xmm6
    \\ mulpd %xmm10, %xmm7
    \\ addpd %xmm7, %xmm3
    \\ movddup -9 * 8(%r13), %xmm7
    \\ mulpd %xmm11, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -8 * 8(%r12), %xmm4
    \\ mulpd %xmm11, %xmm5
    \\ addpd %xmm5, %xmm1
    \\ movddup -6 * 8(%r12), %xmm5
    \\ mulpd %xmm11, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ movddup -4 * 8(%r12), %xmm6
    \\ mulpd %xmm11, %xmm7
    \\ addpd %xmm7, %xmm3
    \\ movddup -2 * 8(%r12), %xmm7
    \\ prefetchw (64 * 1) * 2 - 128 + (64 >> 1)(%r14)
    \\ movlps %xmm0, -16 * 8(%r14)
    \\ movhps %xmm0, -16 * 8 + 8(%r14)
    \\ movlps %xmm1, -14 * 8(%r14)
    \\ movhps %xmm1, -14 * 8 + 8(%r14)
    \\ movlps %xmm2, -12 * 8(%r14)
    \\ movhps %xmm2, -12 * 8 + 8(%r14)
    \\ movlps %xmm3, -10 * 8(%r14)
    \\ movhps %xmm3, -10 * 8 + 8(%r14)
    \\ movlps -8 * 8(%r14), %xmm0
    \\ movhps -8 * 8 + 8(%r14), %xmm0
    \\ movlps -6 * 8(%r14), %xmm1
    \\ movhps -6 * 8 + 8(%r14), %xmm1
    \\ movlps -4 * 8(%r14), %xmm2
    \\ movhps -4 * 8 + 8(%r14), %xmm2
    \\ movlps -2 * 8(%r14), %xmm3
    \\ movhps -2 * 8 + 8(%r14), %xmm3
    \\ subq $-8 * 8, %r12
    \\ subq $-8 * 8, %r13
    \\ subq $-8 * 8, %r14
    \\ subq $1, %rax
    \\ .byte 0x3e
    \\ jg 20009b
    \\ .balign 8
    \\ 20010:
    \\ mulpd %xmm8, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -15 * 8(%r12), %xmm4
    \\ mulpd %xmm8, %xmm5
    \\ addpd %xmm5, %xmm1
    \\ movddup -13 * 8(%r12), %xmm5
    \\ mulpd %xmm8, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ movddup -11 * 8(%r12), %xmm6
    \\ mulpd %xmm8, %xmm7
    \\ addpd %xmm7, %xmm3
    \\ movddup -9 * 8(%r12), %xmm7
    \\ mulpd %xmm9, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -16 * 8(%r13), %xmm4
    \\ mulpd %xmm9, %xmm5
    \\ addpd %xmm5, %xmm1
    \\ movddup -14 * 8(%r13), %xmm5
    \\ mulpd %xmm9, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ movddup -12 * 8(%r13), %xmm6
    \\ mulpd %xmm9, %xmm7
    \\ addpd %xmm7, %xmm3
    \\ movddup -10 * 8(%r13), %xmm7
    \\ mulpd %xmm10, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -15 * 8(%r13), %xmm4
    \\ mulpd %xmm10, %xmm5
    \\ addpd %xmm5, %xmm1
    \\ movddup -13 * 8(%r13), %xmm5
    \\ mulpd %xmm10, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ movddup -11 * 8(%r13), %xmm6
    \\ mulpd %xmm10, %xmm7
    \\ addpd %xmm7, %xmm3
    \\ movddup -9 * 8(%r13), %xmm7
    \\ mulpd %xmm11, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ mulpd %xmm11, %xmm5
    \\ addpd %xmm5, %xmm1
    \\ mulpd %xmm11, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ mulpd %xmm11, %xmm7
    \\ addpd %xmm7, %xmm3
    \\ movlps %xmm0, -16 * 8(%r14)
    \\ movhps %xmm0, -16 * 8 + 8(%r14)
    \\ movlps %xmm1, -14 * 8(%r14)
    \\ movhps %xmm1, -14 * 8 + 8(%r14)
    \\ movlps %xmm2, -12 * 8(%r14)
    \\ movhps %xmm2, -12 * 8 + 8(%r14)
    \\ movlps %xmm3, -10 * 8(%r14)
    \\ movhps %xmm3, -10 * 8 + 8(%r14)
    \\ movlps -8 * 8(%r14), %xmm0
    \\ movhps -8 * 8 + 8(%r14), %xmm0
    \\ movlps -6 * 8(%r14), %xmm1
    \\ movhps -6 * 8 + 8(%r14), %xmm1
    \\ movlps -4 * 8(%r14), %xmm2
    \\ movhps -4 * 8 + 8(%r14), %xmm2
    \\ movlps -2 * 8(%r14), %xmm3
    \\ movhps -2 * 8 + 8(%r14), %xmm3
    \\ subq $-8 * 8, %r12
    \\ subq $-8 * 8, %r13
    \\ subq $-8 * 8, %r14
    \\ .balign 8
    \\ 20011:
    \\ testq $2, %rdi
    \\ je 20012f
    \\ movddup -16 * 8(%r12), %xmm4
    \\ movddup -15 * 8(%r12), %xmm5
    \\ movddup -14 * 8(%r12), %xmm6
    \\ movddup -13 * 8(%r12), %xmm7
    \\ mulpd %xmm8, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -16 * 8(%r13), %xmm4
    \\ mulpd %xmm8, %xmm6
    \\ addpd %xmm6, %xmm1
    \\ movddup -14 * 8(%r13), %xmm6
    \\ mulpd %xmm9, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ movddup -15 * 8(%r13), %xmm5
    \\ mulpd %xmm9, %xmm7
    \\ addpd %xmm7, %xmm1
    \\ movddup -13 * 8(%r13), %xmm7
    \\ mulpd %xmm10, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ mulpd %xmm10, %xmm6
    \\ addpd %xmm6, %xmm1
    \\ mulpd %xmm11, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ mulpd %xmm11, %xmm7
    \\ addpd %xmm7, %xmm1
    \\ movlps %xmm0, -16 * 8(%r14)
    \\ movhps %xmm0, -16 * 8 + 8(%r14)
    \\ movlps %xmm1, -14 * 8(%r14)
    \\ movhps %xmm1, -14 * 8 + 8(%r14)
    \\ movaps %xmm2, %xmm0
    \\ addq $4 * 8, %r12
    \\ addq $4 * 8, %r13
    \\ addq $4 * 8, %r14
    \\ .balign 8
    \\ 20012:
    \\ testq $1, %rdi
    \\ je 20013f
    \\ movddup -16 * 8(%r12), %xmm4
    \\ movddup -15 * 8(%r12), %xmm5
    \\ movddup -16 * 8(%r13), %xmm6
    \\ movddup -15 * 8(%r13), %xmm7
    \\ mulpd %xmm8, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ mulpd %xmm9, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ mulpd %xmm10, %xmm6
    \\ addpd %xmm6, %xmm0
    \\ mulpd %xmm11, %xmm7
    \\ addpd %xmm7, %xmm0
    \\ movlps %xmm0, -16 * 8(%r14)
    \\ movhps %xmm0, -16 * 8 + 8(%r14)
    \\ .balign 8
    \\ 20013:
    \\ cmpq $1, %rsi
    \\ jl 20018f
    \\ leaq 16 * 8(%r15), %r14
    \\ movq %rcx, %r12
    \\ movddup 0 * 8(%r9), %xmm8
    \\ movddup 1 * 8(%r9), %xmm9
    \\ addq %rdx, %r9
    \\ pcmpeqb %xmm5, %xmm5
    \\ psllq $63, %xmm5
    \\ shufps $0x40, %xmm5, %xmm5
    \\ movlps 48 (%rsp), %xmm6
    \\ movhps 56 (%rsp), %xmm6
    \\ pshufd $0x4e, %xmm6, %xmm7
    \\ xorps %xmm5, %xmm7
    \\ mulpd %xmm6, %xmm8
    \\ mulpd %xmm7, %xmm9
    \\ subpd %xmm9, %xmm8
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ xorps %xmm5, %xmm9
    \\ movlps -16 * 8(%r14), %xmm0
    \\ movhps -16 * 8 + 8(%r14), %xmm0
    \\ movlps -14 * 8(%r14), %xmm1
    \\ movhps -14 * 8 + 8(%r14), %xmm1
    \\ movlps -12 * 8(%r14), %xmm2
    \\ movhps -12 * 8 + 8(%r14), %xmm2
    \\ movlps -10 * 8(%r14), %xmm3
    \\ movhps -10 * 8 + 8(%r14), %xmm3
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20016f
    \\ movddup -16 * 8(%r12), %xmm4
    \\ movddup -14 * 8(%r12), %xmm5
    \\ movddup -12 * 8(%r12), %xmm6
    \\ movddup -10 * 8(%r12), %xmm7
    \\ decq %rax
    \\ jle 20015f
    \\
;
