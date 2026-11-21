//! Section 1 of 2 of the text of `dsymv_U`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/symv_U_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

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
    \\ leaq (,%r9, 8), %r9
    \\ leaq (,%r11, 8), %r11
    \\ leaq (,%rcx, 8), %rcx
    \\ testq %rdi, %rdi
    \\ jle 20024f
    \\ negq %rsi
    \\ addq %rdi, %rsi
    \\ movq %rsi, %rax
    \\ imulq %rcx, %rax
    \\ addq %rax, %rdx
    \\ unpcklpd %xmm0, %xmm0
    \\ movq %r12, %r13
    \\ movq %rdi, %rax
    \\ sarq $3, %rax
    \\ jle 20001f
    \\ .balign 8
    \\ 20000:
    \\ movsd 0 * 8(%r8), %xmm1
    \\ addq %r9, %r8
    \\ movhps 0 * 8(%r8), %xmm1
    \\ addq %r9, %r8
    \\ movsd 0 * 8(%r8), %xmm2
    \\ addq %r9, %r8
    \\ movhps 0 * 8(%r8), %xmm2
    \\ addq %r9, %r8
    \\ movsd 0 * 8(%r8), %xmm3
    \\ addq %r9, %r8
    \\ movhps 0 * 8(%r8), %xmm3
    \\ addq %r9, %r8
    \\ movsd 0 * 8(%r8), %xmm4
    \\ addq %r9, %r8
    \\ movhps 0 * 8(%r8), %xmm4
    \\ addq %r9, %r8
    \\ mulpd %xmm0, %xmm1
    \\ mulpd %xmm0, %xmm2
    \\ mulpd %xmm0, %xmm3
    \\ mulpd %xmm0, %xmm4
    \\ movaps %xmm1, 0 * 8(%r13)
    \\ movaps %xmm2, 2 * 8(%r13)
    \\ movaps %xmm3, 4 * 8(%r13)
    \\ movaps %xmm4, 6 * 8(%r13)
    \\ addq $8 * 8, %r13
    \\ decq %rax
    \\ jg 20000b
    \\ .balign 8
    \\ 20001:
    \\ movq %rdi, %rax
    \\ andq $7, %rax
    \\ jle 20003f
    \\ .balign 8
    \\ 20002:
    \\ movsd 0 * 8(%r8), %xmm1
    \\ addq %r9, %r8
    \\ mulsd %xmm0, %xmm1
    \\ movlps %xmm1, 0 * 8(%r13)
    \\ addq $1 * 8, %r13
    \\ decq %rax
    \\ jg 20002b
    \\ .balign 8
    \\ 20003:
    \\ movq %r10, %r8
    \\ addq $512, %r13
    \\ andq $-512, %r13
    \\ cmpq $8, %r11
    \\ je 20007f
    \\ movq %r10, %r14
    \\ movq %r13, %r8
    \\ movq %rdi, %rax
    \\ sarq $3, %rax
    \\ jle 20005f
    \\ .balign 8
    \\ 20004:
    \\ movsd 0 * 8(%r14), %xmm0
    \\ addq %r11, %r14
    \\ movhps 0 * 8(%r14), %xmm0
    \\ addq %r11, %r14
    \\ movsd 0 * 8(%r14), %xmm1
    \\ addq %r11, %r14
    \\ movhps 0 * 8(%r14), %xmm1
    \\ addq %r11, %r14
    \\ movsd 0 * 8(%r14), %xmm2
    \\ addq %r11, %r14
    \\ movhps 0 * 8(%r14), %xmm2
    \\ addq %r11, %r14
    \\ movsd 0 * 8(%r14), %xmm3
    \\ addq %r11, %r14
    \\ movhps 0 * 8(%r14), %xmm3
    \\ addq %r11, %r14
    \\ movaps %xmm0, 0 * 8(%r13)
    \\ movaps %xmm1, 2 * 8(%r13)
    \\ movaps %xmm2, 4 * 8(%r13)
    \\ movaps %xmm3, 6 * 8(%r13)
    \\ addq $8 * 8, %r13
    \\ decq %rax
    \\ jg 20004b
    \\ .balign 8
    \\ 20005:
    \\ movq %rdi, %rax
    \\ andq $7, %rax
    \\ jle 20007f
    \\ .balign 8
    \\ 20006:
    \\ movsd 0 * 8(%r14), %xmm0
    \\ addq %r11, %r14
    \\ movsd %xmm0, 0 * 8(%r13)
    \\ addq $1 * 8, %r13
    \\ decq %rax
    \\ jg 20006b
    \\ .balign 8
    \\ 20007:
    \\ movq %rsi, %rax
    \\ addq $4, %rax
    \\ cmpq %rdi, %rax
    \\ jg 20012f
    \\ .balign 8
    \\ 20008:
    \\ movq %rdx, %rbx
    \\ leaq (%rdx, %rcx, 2), %rbp
    \\ leaq (%rdx, %rcx, 4), %rdx
    \\ movddup 0 * 8(%r12, %rsi, 8), %xmm4
    \\ movddup 1 * 8(%r12, %rsi, 8), %xmm5
    \\ movddup 2 * 8(%r12, %rsi, 8), %xmm6
    \\ movddup 3 * 8(%r12, %rsi, 8), %xmm7
    \\ pxor %xmm8, %xmm8
    \\ pxor %xmm9, %xmm9
    \\ pxor %xmm10, %xmm10
    \\ pxor %xmm11, %xmm11
    \\ movaps 0 * 8(%r12), %xmm0
    \\ movaps 2 * 8(%r12), %xmm1
    \\ movsd 0 * 8(%rbx), %xmm12
    \\ movhps 1 * 8(%rbx), %xmm12
    \\ movsd 2 * 8(%rbx), %xmm13
    \\ movhps 3 * 8(%rbx), %xmm13
    \\ movsd 0 * 8(%rbx, %rcx, 1), %xmm14
    \\ movhps 1 * 8(%rbx, %rcx, 1), %xmm14
    \\ movsd 0 * 8(%r8), %xmm2
    \\ movhps 1 * 8(%r8), %xmm2
    \\ movsd 2 * 8(%r8), %xmm3
    \\ movhps 3 * 8(%r8), %xmm3
    \\ movq %r12, %r13
    \\ movq %r8, %r14
    \\ movq %rsi, %rax
    \\ sarq $3, %rax
    \\ jle 20010f
    \\ .balign 8
    \\ 20009:
    \\ movaps %xmm0, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm4, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm2
    \\ movsd 2 * 8(%rbx, %rcx, 1), %xmm12
    \\ movhps 3 * 8(%rbx, %rcx, 1), %xmm12
    \\ prefetcht0 (16 * 12)(%rbx)
    \\ movaps %xmm1, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm4, %xmm13
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm13, %xmm3
    \\ movsd 0 * 8(%rbp), %xmm13
    \\ movhps 1 * 8(%rbp), %xmm13
    \\ movaps %xmm0, %xmm15
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm5, %xmm14
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm14, %xmm2
    \\ movsd 2 * 8(%rbp), %xmm14
    \\ movhps 3 * 8(%rbp), %xmm14
    \\ movaps %xmm1, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm12, %xmm3
    \\ movsd 0 * 8(%rbp, %rcx, 1), %xmm12
    \\ movhps 1 * 8(%rbp, %rcx, 1), %xmm12
    \\ movaps %xmm0, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm6, %xmm13
    \\ addpd %xmm15, %xmm10
    \\ addpd %xmm13, %xmm2
    \\ movsd 2 * 8(%rbp, %rcx, 1), %xmm13
    \\ movhps 3 * 8(%rbp, %rcx, 1), %xmm13
    \\ prefetcht0 (16 * 12)(%rbx, %rcx, 1)
    \\ movaps %xmm1, %xmm15
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm6, %xmm14
    \\ addpd %xmm15, %xmm10
    \\ addpd %xmm14, %xmm3
    \\ movsd 4 * 8(%rbx), %xmm14
    \\ movhps 5 * 8(%rbx), %xmm14
    \\ movaps %xmm0, %xmm15
    \\ movaps 4 * 8(%r13), %xmm0
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm7, %xmm12
    \\ addpd %xmm15, %xmm11
    \\ addpd %xmm12, %xmm2
    \\ movsd 6 * 8(%rbx), %xmm12
    \\ movhps 7 * 8(%rbx), %xmm12
    \\ movaps %xmm1, %xmm15
    \\ movaps 6 * 8(%r13), %xmm1
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm7, %xmm13
    \\ addpd %xmm15, %xmm11
    \\ addpd %xmm13, %xmm3
    \\ movsd 4 * 8(%rbx, %rcx, 1), %xmm13
    \\ movhps 5 * 8(%rbx, %rcx, 1), %xmm13
    \\ movsd %xmm2, 0 * 8(%r14)
    \\ movhps %xmm2, 1 * 8(%r14)
    \\ movsd 4 * 8(%r14), %xmm2
    \\ movhps 5 * 8(%r14), %xmm2
    \\ movsd %xmm3, 2 * 8(%r14)
    \\ movhps %xmm3, 3 * 8(%r14)
    \\ movsd 6 * 8(%r14), %xmm3
    \\ movhps 7 * 8(%r14), %xmm3
    \\ movaps %xmm0, %xmm15
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm4, %xmm14
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm14, %xmm2
    \\ movsd 6 * 8(%rbx, %rcx, 1), %xmm14
    \\ movhps 7 * 8(%rbx, %rcx, 1), %xmm14
    \\ prefetcht0 (16 * 12)(%rbp)
    \\ movaps %xmm1, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm4, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm3
    \\ movsd 4 * 8(%rbp), %xmm12
    \\ movhps 5 * 8(%rbp), %xmm12
    \\ movaps %xmm0, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm5, %xmm13
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm13, %xmm2
    \\ movsd 6 * 8(%rbp), %xmm13
    \\ movhps 7 * 8(%rbp), %xmm13
    \\ movaps %xmm1, %xmm15
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm5, %xmm14
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm14, %xmm3
    \\ movsd 4 * 8(%rbp, %rcx, 1), %xmm14
    \\ movhps 5 * 8(%rbp, %rcx, 1), %xmm14
    \\ movaps %xmm0, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm6, %xmm12
    \\ addpd %xmm15, %xmm10
    \\ addpd %xmm12, %xmm2
    \\ movsd 6 * 8(%rbp, %rcx, 1), %xmm12
    \\ movhps 7 * 8(%rbp, %rcx, 1), %xmm12
    \\ prefetcht0 (16 * 12)(%rbp, %rcx, 1)
    \\ movaps %xmm1, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm6, %xmm13
    \\ addpd %xmm15, %xmm10
    \\ addpd %xmm13, %xmm3
    \\ movsd 10 * 8(%rbx), %xmm13
    \\ movhps 11 * 8(%rbx), %xmm13
    \\ movaps %xmm0, %xmm15
    \\ movaps 8 * 8(%r13), %xmm0
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm7, %xmm14
    \\ addpd %xmm15, %xmm11
    \\ addpd %xmm14, %xmm2
    \\ movsd 8 * 8(%rbx, %rcx, 1), %xmm14
    \\ movhps 9 * 8(%rbx, %rcx, 1), %xmm14
    \\ movaps %xmm1, %xmm15
    \\ movaps 10 * 8(%r13), %xmm1
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm7, %xmm12
    \\ addpd %xmm15, %xmm11
    \\ addpd %xmm12, %xmm3
    \\ movsd 8 * 8(%rbx), %xmm12
    \\ movhps 9 * 8(%rbx), %xmm12
    \\ movsd %xmm2, 4 * 8(%r14)
    \\ movhps %xmm2, 5 * 8(%r14)
    \\ movsd 8 * 8(%r14), %xmm2
    \\ movhps 9 * 8(%r14), %xmm2
    \\ movsd %xmm3, 6 * 8(%r14)
    \\ movhps %xmm3, 7 * 8(%r14)
    \\ movsd 10 * 8(%r14), %xmm3
    \\ movhps 11 * 8(%r14), %xmm3
    \\ addq $8 * 8, %r13
    \\ addq $8 * 8, %r14
    \\ addq $8 * 8, %rbx
    \\ addq $8 * 8, %rbp
    \\ decq %rax
    \\ jg 20009b
    \\
;
