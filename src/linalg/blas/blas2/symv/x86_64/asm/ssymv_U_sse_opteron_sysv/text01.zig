//! Section 1 of 2 of the text of `ssymv_U`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/symv_U_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 512
    \\ subq $80, %rsp
    \\ movq %rbx, 0(%rsp)
    \\ movq %rbp, 8(%rsp)
    \\ movq %r12, 16(%rsp)
    \\ movq %r13, 24(%rsp)
    \\ movq %r14, 32(%rsp)
    \\ movq %r15, 40(%rsp)
    \\ movq 8 + 80(%rsp), %r10
    \\ movq 16 + 80(%rsp), %r11
    \\ movq 24 + 80(%rsp), %r12
    \\ leaq (,%r9, 4), %r9
    \\ leaq (,%r11, 4), %r11
    \\ leaq (,%rcx, 4), %rcx
    \\ testq %rdi, %rdi
    \\ jle 20023f
    \\ negq %rsi
    \\ addq %rdi, %rsi
    \\ movq %rsi, %rax
    \\ imulq %rcx, %rax
    \\ addq %rax, %rdx
    \\ shufps $0, %xmm0, %xmm0
    \\ movq %r12, %r13
    \\ movq %rdi, %rax
    \\ sarq $3, %rax
    \\ jle 20001f
    \\ .balign 8
    \\ 20000:
    \\ movss 0 * 4(%r8), %xmm1
    \\ addq %r9, %r8
    \\ movss 0 * 4(%r8), %xmm2
    \\ addq %r9, %r8
    \\ movss 0 * 4(%r8), %xmm3
    \\ addq %r9, %r8
    \\ movss 0 * 4(%r8), %xmm4
    \\ addq %r9, %r8
    \\ movss 0 * 4(%r8), %xmm5
    \\ addq %r9, %r8
    \\ movss 0 * 4(%r8), %xmm6
    \\ addq %r9, %r8
    \\ movss 0 * 4(%r8), %xmm7
    \\ addq %r9, %r8
    \\ movss 0 * 4(%r8), %xmm8
    \\ addq %r9, %r8
    \\ mulss %xmm0, %xmm1
    \\ mulss %xmm0, %xmm2
    \\ mulss %xmm0, %xmm3
    \\ mulss %xmm0, %xmm4
    \\ mulss %xmm0, %xmm5
    \\ mulss %xmm0, %xmm6
    \\ mulss %xmm0, %xmm7
    \\ mulss %xmm0, %xmm8
    \\ movss %xmm1, 0 * 4(%r13)
    \\ movss %xmm2, 1 * 4(%r13)
    \\ movss %xmm3, 2 * 4(%r13)
    \\ movss %xmm4, 3 * 4(%r13)
    \\ movss %xmm5, 4 * 4(%r13)
    \\ movss %xmm6, 5 * 4(%r13)
    \\ movss %xmm7, 6 * 4(%r13)
    \\ movss %xmm8, 7 * 4(%r13)
    \\ addq $8 * 4, %r13
    \\ decq %rax
    \\ jg 20000b
    \\ .balign 8
    \\ 20001:
    \\ movq %rdi, %rax
    \\ andq $7, %rax
    \\ jle 20003f
    \\ .balign 8
    \\ 20002:
    \\ movss 0 * 4(%r8), %xmm1
    \\ addq %r9, %r8
    \\ mulss %xmm0, %xmm1
    \\ movss %xmm1, 0 * 4(%r13)
    \\ addq $1 * 4, %r13
    \\ decq %rax
    \\ jg 20002b
    \\ .balign 8
    \\ 20003:
    \\ movq %r10, %r8
    \\ addq $512, %r13
    \\ andq $-512, %r13
    \\ cmpq $4, %r11
    \\ je 20007f
    \\ movq %r10, %r14
    \\ movq %r13, %r8
    \\ movq %rdi, %rax
    \\ sarq $3, %rax
    \\ jle 20005f
    \\ .balign 8
    \\ 20004:
    \\ movss 0 * 4(%r14), %xmm0
    \\ addq %r11, %r14
    \\ movss 0 * 4(%r14), %xmm1
    \\ addq %r11, %r14
    \\ movss 0 * 4(%r14), %xmm2
    \\ addq %r11, %r14
    \\ movss 0 * 4(%r14), %xmm3
    \\ addq %r11, %r14
    \\ movss 0 * 4(%r14), %xmm4
    \\ addq %r11, %r14
    \\ movss 0 * 4(%r14), %xmm5
    \\ addq %r11, %r14
    \\ movss 0 * 4(%r14), %xmm6
    \\ addq %r11, %r14
    \\ movss 0 * 4(%r14), %xmm7
    \\ addq %r11, %r14
    \\ movss %xmm0, 0 * 4(%r13)
    \\ movss %xmm1, 1 * 4(%r13)
    \\ movss %xmm2, 2 * 4(%r13)
    \\ movss %xmm3, 3 * 4(%r13)
    \\ movss %xmm4, 4 * 4(%r13)
    \\ movss %xmm5, 5 * 4(%r13)
    \\ movss %xmm6, 6 * 4(%r13)
    \\ movss %xmm7, 7 * 4(%r13)
    \\ addq $8 * 4, %r13
    \\ decq %rax
    \\ jg 20004b
    \\ .balign 8
    \\ 20005:
    \\ movq %rdi, %rax
    \\ andq $7, %rax
    \\ jle 20007f
    \\ .balign 8
    \\ 20006:
    \\ movss 0 * 4(%r14), %xmm0
    \\ addq %r11, %r14
    \\ movss %xmm0, 0 * 4(%r13)
    \\ addq $1 * 4, %r13
    \\ decq %rax
    \\ jg 20006b
    \\ .balign 8
    \\ 20007:
    \\ movq %rsi, %rax
    \\ addq $4, %rax
    \\ cmpq %rdi, %rax
    \\ jg 20013f
    \\ .balign 8
    \\ 20008:
    \\ movq %rdx, %rbx
    \\ leaq (%rdx, %rcx, 2), %rbp
    \\ leaq (%rdx, %rcx, 4), %rdx
    \\ movaps 0 * 4(%r12, %rsi, 4), %xmm3
    \\ pshufd $0x00, %xmm3, %xmm0
    \\ pshufd $0x55, %xmm3, %xmm1
    \\ pshufd $0xaa, %xmm3, %xmm2
    \\ pshufd $0xff, %xmm3, %xmm3
    \\ pxor %xmm4, %xmm4
    \\ pxor %xmm5, %xmm5
    \\ pxor %xmm6, %xmm6
    \\ pxor %xmm7, %xmm7
    \\ movaps 0 * 4(%r12), %xmm8
    \\ movaps 4 * 4(%r12), %xmm9
    \\ movlps 0 * 4(%rbx), %xmm12
    \\ movhps 2 * 4(%rbx), %xmm12
    \\ movlps 0 * 4(%rbx, %rcx, 1), %xmm13
    \\ movhps 2 * 4(%rbx, %rcx, 1), %xmm13
    \\ movlps 0 * 4(%rbp), %xmm14
    \\ movhps 2 * 4(%rbp), %xmm14
    \\ movlps 0 * 4(%rbp, %rcx, 1), %xmm15
    \\ movhps 2 * 4(%rbp, %rcx, 1), %xmm15
    \\ movlps 0 * 4(%r8), %xmm10
    \\ movhps 2 * 4(%r8), %xmm10
    \\ movq %r12, %r13
    \\ movq %r8, %r14
    \\ movq %rsi, %rax
    \\ sarq $4, %rax
    \\ jle 20010f
    \\ .balign 8
    \\ 20009:
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm12, %xmm11
    \\ mulps %xmm0, %xmm12
    \\ addps %xmm11, %xmm4
    \\ addps %xmm12, %xmm10
    \\ movlps 4 * 4(%rbx), %xmm12
    \\ movhps 6 * 4(%rbx), %xmm12
    \\ prefetch (16 * 8)(%rbx)
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm13, %xmm11
    \\ mulps %xmm1, %xmm13
    \\ addps %xmm11, %xmm5
    \\ addps %xmm13, %xmm10
    \\ movlps 4 * 4(%rbx, %rcx, 1), %xmm13
    \\ movhps 6 * 4(%rbx, %rcx, 1), %xmm13
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm14, %xmm11
    \\ mulps %xmm2, %xmm14
    \\ addps %xmm11, %xmm6
    \\ addps %xmm14, %xmm10
    \\ movlps 4 * 4(%rbp), %xmm14
    \\ movhps 6 * 4(%rbp), %xmm14
    \\ prefetch (16 * 8)(%r13)
    \\ movaps %xmm8, %xmm11
    \\ movaps 8 * 4(%r13), %xmm8
    \\ mulps %xmm15, %xmm11
    \\ mulps %xmm3, %xmm15
    \\ addps %xmm11, %xmm7
    \\ addps %xmm15, %xmm10
    \\ movlps 4 * 4(%rbp, %rcx, 1), %xmm15
    \\ movhps 6 * 4(%rbp, %rcx, 1), %xmm15
    \\ movlps %xmm10, 0 * 4(%r14)
    \\ movhps %xmm10, 2 * 4(%r14)
    \\ movlps 4 * 4(%r14), %xmm10
    \\ movhps 6 * 4(%r14), %xmm10
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm12, %xmm11
    \\ mulps %xmm0, %xmm12
    \\ addps %xmm11, %xmm4
    \\ addps %xmm12, %xmm10
    \\ movlps 8 * 4(%rbx), %xmm12
    \\ movhps 10 * 4(%rbx), %xmm12
    \\ prefetch (16 * 8)(%rbx, %rcx, 1)
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm13, %xmm11
    \\ mulps %xmm1, %xmm13
    \\ addps %xmm11, %xmm5
    \\ addps %xmm13, %xmm10
    \\ movlps 8 * 4(%rbx, %rcx, 1), %xmm13
    \\ movhps 10 * 4(%rbx, %rcx, 1), %xmm13
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm14, %xmm11
    \\ mulps %xmm2, %xmm14
    \\ addps %xmm11, %xmm6
    \\ addps %xmm14, %xmm10
    \\ movlps 8 * 4(%rbp), %xmm14
    \\ movhps 10 * 4(%rbp), %xmm14
    \\ movaps %xmm9, %xmm11
    \\ movaps 12 * 4(%r13), %xmm9
    \\ mulps %xmm15, %xmm11
    \\ mulps %xmm3, %xmm15
    \\ addps %xmm11, %xmm7
    \\ addps %xmm15, %xmm10
    \\ movlps 8 * 4(%rbp, %rcx, 1), %xmm15
    \\ movhps 10 * 4(%rbp, %rcx, 1), %xmm15
    \\ movlps %xmm10, 4 * 4(%r14)
    \\ movhps %xmm10, 6 * 4(%r14)
    \\ movlps 8 * 4(%r14), %xmm10
    \\ movhps 10 * 4(%r14), %xmm10
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm12, %xmm11
    \\ mulps %xmm0, %xmm12
    \\ addps %xmm11, %xmm4
    \\ addps %xmm12, %xmm10
    \\ movlps 12 * 4(%rbx), %xmm12
    \\ movhps 14 * 4(%rbx), %xmm12
    \\ prefetch (16 * 8)(%rbp)
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm13, %xmm11
    \\ mulps %xmm1, %xmm13
    \\ addps %xmm11, %xmm5
    \\ addps %xmm13, %xmm10
    \\ movlps 12 * 4(%rbx, %rcx, 1), %xmm13
    \\ movhps 14 * 4(%rbx, %rcx, 1), %xmm13
    \\ movaps %xmm8, %xmm11
    \\ mulps %xmm14, %xmm11
    \\ mulps %xmm2, %xmm14
    \\ addps %xmm11, %xmm6
    \\ addps %xmm14, %xmm10
    \\ movlps 12 * 4(%rbp), %xmm14
    \\ movhps 14 * 4(%rbp), %xmm14
    \\ prefetchw (16 * 8)(%r14)
    \\ movaps %xmm8, %xmm11
    \\ movaps 16 * 4(%r13), %xmm8
    \\ mulps %xmm15, %xmm11
    \\ mulps %xmm3, %xmm15
    \\ addps %xmm11, %xmm7
    \\ addps %xmm15, %xmm10
    \\ movlps 12 * 4(%rbp, %rcx, 1), %xmm15
    \\ movhps 14 * 4(%rbp, %rcx, 1), %xmm15
    \\ movlps %xmm10, 8 * 4(%r14)
    \\ movhps %xmm10, 10 * 4(%r14)
    \\ movlps 12 * 4(%r14), %xmm10
    \\ movhps 14 * 4(%r14), %xmm10
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm12, %xmm11
    \\ mulps %xmm0, %xmm12
    \\ addps %xmm11, %xmm4
    \\ addps %xmm12, %xmm10
    \\ movlps 16 * 4(%rbx), %xmm12
    \\ movhps 18 * 4(%rbx), %xmm12
    \\ prefetch (16 * 8)(%rbp, %rcx, 1)
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm13, %xmm11
    \\ mulps %xmm1, %xmm13
    \\ addps %xmm11, %xmm5
    \\ addps %xmm13, %xmm10
    \\ movlps 16 * 4(%rbx, %rcx, 1), %xmm13
    \\ movhps 18 * 4(%rbx, %rcx, 1), %xmm13
    \\ movaps %xmm9, %xmm11
    \\ mulps %xmm14, %xmm11
    \\ mulps %xmm2, %xmm14
    \\ addps %xmm11, %xmm6
    \\ addps %xmm14, %xmm10
    \\ movlps 16 * 4(%rbp), %xmm14
    \\ movhps 18 * 4(%rbp), %xmm14
    \\ movaps %xmm9, %xmm11
    \\ movaps 20 * 4(%r13), %xmm9
    \\ mulps %xmm15, %xmm11
    \\ mulps %xmm3, %xmm15
    \\ addps %xmm11, %xmm7
    \\ addps %xmm15, %xmm10
    \\ movlps 16 * 4(%rbp, %rcx, 1), %xmm15
    \\ movhps 18 * 4(%rbp, %rcx, 1), %xmm15
    \\ movlps %xmm10, 12 * 4(%r14)
    \\ movhps %xmm10, 14 * 4(%r14)
    \\ movlps 16 * 4(%r14), %xmm10
    \\ movhps 18 * 4(%r14), %xmm10
    \\ addq $16 * 4, %r13
    \\ addq $16 * 4, %r14
    \\ addq $16 * 4, %rbx
    \\ addq $16 * 4, %rbp
    \\ decq %rax
    \\ jg 20009b
    \\
;
