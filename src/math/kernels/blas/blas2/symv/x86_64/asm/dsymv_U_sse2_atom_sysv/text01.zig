//! Section 1 of 2 of the text of `dsymv_U`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/symv_U_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1, bp2, bp3, bp4
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
    \\ movhpd 0 * 8(%r8), %xmm1
    \\ addq %r9, %r8
    \\ movsd 0 * 8(%r8), %xmm2
    \\ addq %r9, %r8
    \\ movhpd 0 * 8(%r8), %xmm2
    \\ addq %r9, %r8
    \\ movsd 0 * 8(%r8), %xmm3
    \\ addq %r9, %r8
    \\ movhpd 0 * 8(%r8), %xmm3
    \\ addq %r9, %r8
    \\ movsd 0 * 8(%r8), %xmm4
    \\ addq %r9, %r8
    \\ movhpd 0 * 8(%r8), %xmm4
    \\ addq %r9, %r8
    \\ mulpd %xmm0, %xmm1
    \\ mulpd %xmm0, %xmm2
    \\ mulpd %xmm0, %xmm3
    \\ mulpd %xmm0, %xmm4
    \\ movapd %xmm1, 0 * 8(%r13)
    \\ movapd %xmm2, 2 * 8(%r13)
    \\ movapd %xmm3, 4 * 8(%r13)
    \\ movapd %xmm4, 6 * 8(%r13)
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
    \\ movlpd %xmm1, 0 * 8(%r13)
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
    \\ movhpd 0 * 8(%r14), %xmm0
    \\ addq %r11, %r14
    \\ movsd 0 * 8(%r14), %xmm1
    \\ addq %r11, %r14
    \\ movhpd 0 * 8(%r14), %xmm1
    \\ addq %r11, %r14
    \\ movsd 0 * 8(%r14), %xmm2
    \\ addq %r11, %r14
    \\ movhpd 0 * 8(%r14), %xmm2
    \\ addq %r11, %r14
    \\ movsd 0 * 8(%r14), %xmm3
    \\ addq %r11, %r14
    \\ movhpd 0 * 8(%r14), %xmm3
    \\ addq %r11, %r14
    \\ movapd %xmm0, 0 * 8(%r13)
    \\ movapd %xmm1, 2 * 8(%r13)
    \\ movapd %xmm2, 4 * 8(%r13)
    \\ movapd %xmm3, 6 * 8(%r13)
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
    \\ movapd 0 * 8(%r12), %xmm0
    \\ movapd 2 * 8(%r12), %xmm1
    \\ movsd 0 * 8(%rbx), %xmm12
    \\ movhpd 1 * 8(%rbx), %xmm12
    \\ movsd 2 * 8(%rbx), %xmm13
    \\ movhpd 3 * 8(%rbx), %xmm13
    \\ movsd 0 * 8(%rbx, %rcx, 1), %xmm14
    \\ movhpd 1 * 8(%rbx, %rcx, 1), %xmm14
    \\ movsd 0 * 8(%r8), %xmm2
    \\ movhpd 1 * 8(%r8), %xmm2
    \\ movsd 2 * 8(%r8), %xmm3
    \\ movhpd 3 * 8(%r8), %xmm3
    \\ movq %r12, %r13
    \\ movq %r8, %r14
    \\ movq %rsi, %rax
    \\ sarq $3, %rax
    \\ jle 20010f
    \\ .balign 8
    \\ 20009:
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm4, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm2
    \\ movsd 2 * 8(%rbx, %rcx, 1), %xmm12
    \\ movhpd 3 * 8(%rbx, %rcx, 1), %xmm12
    \\ .if \bp0
    \\ prefetcht0 (16 * 12)(%rbx)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 16)(%rbx)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 (8 * 24)(%rbx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 (16 * 20)(%rbx)
    \\ .endif
    \\ movapd %xmm1, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm4, %xmm13
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm13, %xmm3
    \\ movsd 0 * 8(%rbp), %xmm13
    \\ movhpd 1 * 8(%rbp), %xmm13
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm5, %xmm14
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm14, %xmm2
    \\ movsd 2 * 8(%rbp), %xmm14
    \\ movhpd 3 * 8(%rbp), %xmm14
    \\ .if \bp4
    \\ prefetcht0 (16 * 12)(%r13)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 16)(%r13)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 (8 * 24)(%r13)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 (16 * 20)(%r13)
    \\ .endif
    \\ movapd %xmm1, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm12, %xmm3
    \\ movsd 0 * 8(%rbp, %rcx, 1), %xmm12
    \\ movhpd 1 * 8(%rbp, %rcx, 1), %xmm12
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm6, %xmm13
    \\ addpd %xmm15, %xmm10
    \\ addpd %xmm13, %xmm2
    \\ movsd 2 * 8(%rbp, %rcx, 1), %xmm13
    \\ movhpd 3 * 8(%rbp, %rcx, 1), %xmm13
    \\ .if \bp0
    \\ prefetcht0 (16 * 12)(%rbx, %rcx, 1)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 16)(%rbx, %rcx, 1)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 (8 * 24)(%rbx, %rcx, 1)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 (16 * 20)(%rbx, %rcx, 1)
    \\ .endif
    \\ movapd %xmm1, %xmm15
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm6, %xmm14
    \\ addpd %xmm15, %xmm10
    \\ addpd %xmm14, %xmm3
    \\ movsd 4 * 8(%rbx), %xmm14
    \\ movhpd 5 * 8(%rbx), %xmm14
    \\ movapd %xmm0, %xmm15
    \\ movapd 4 * 8(%r13), %xmm0
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm7, %xmm12
    \\ addpd %xmm15, %xmm11
    \\ addpd %xmm12, %xmm2
    \\ movsd 6 * 8(%rbx), %xmm12
    \\ movhpd 7 * 8(%rbx), %xmm12
    \\ movapd %xmm1, %xmm15
    \\ movapd 6 * 8(%r13), %xmm1
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm7, %xmm13
    \\ addpd %xmm15, %xmm11
    \\ addpd %xmm13, %xmm3
    \\ movsd 4 * 8(%rbx, %rcx, 1), %xmm13
    \\ movhpd 5 * 8(%rbx, %rcx, 1), %xmm13
    \\ movsd %xmm2, 0 * 8(%r14)
    \\ movhpd %xmm2, 1 * 8(%r14)
    \\ movsd 4 * 8(%r14), %xmm2
    \\ movhpd 5 * 8(%r14), %xmm2
    \\ movsd %xmm3, 2 * 8(%r14)
    \\ movhpd %xmm3, 3 * 8(%r14)
    \\ movsd 6 * 8(%r14), %xmm3
    \\ movhpd 7 * 8(%r14), %xmm3
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm4, %xmm14
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm14, %xmm2
    \\ movsd 6 * 8(%rbx, %rcx, 1), %xmm14
    \\ movhpd 7 * 8(%rbx, %rcx, 1), %xmm14
    \\ .if \bp0
    \\ prefetcht0 (16 * 12)(%rbp)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 16)(%rbp)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 (8 * 24)(%rbp)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 (16 * 20)(%rbp)
    \\ .endif
    \\ movapd %xmm1, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm4, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm3
    \\ movsd 4 * 8(%rbp), %xmm12
    \\ movhpd 5 * 8(%rbp), %xmm12
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm5, %xmm13
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm13, %xmm2
    \\ movsd 6 * 8(%rbp), %xmm13
    \\ movhpd 7 * 8(%rbp), %xmm13
    \\ .if \bp4
    \\ prefetcht0 (16 * 12)(%r14)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw (16 * 16)(%r14)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 (8 * 24)(%r14)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 (16 * 20)(%r14)
    \\ .endif
    \\ movapd %xmm1, %xmm15
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm5, %xmm14
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm14, %xmm3
    \\ movsd 4 * 8(%rbp, %rcx, 1), %xmm14
    \\ movhpd 5 * 8(%rbp, %rcx, 1), %xmm14
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm6, %xmm12
    \\ addpd %xmm15, %xmm10
    \\ addpd %xmm12, %xmm2
    \\ movsd 6 * 8(%rbp, %rcx, 1), %xmm12
    \\ movhpd 7 * 8(%rbp, %rcx, 1), %xmm12
    \\ .if \bp0
    \\ prefetcht0 (16 * 12)(%rbp, %rcx, 1)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 16)(%rbp, %rcx, 1)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 (8 * 24)(%rbp, %rcx, 1)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 (16 * 20)(%rbp, %rcx, 1)
    \\ .endif
    \\ movapd %xmm1, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm6, %xmm13
    \\ addpd %xmm15, %xmm10
    \\ addpd %xmm13, %xmm3
    \\ movsd 10 * 8(%rbx), %xmm13
    \\ movhpd 11 * 8(%rbx), %xmm13
    \\ movapd %xmm0, %xmm15
    \\ movapd 8 * 8(%r13), %xmm0
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm7, %xmm14
    \\ addpd %xmm15, %xmm11
    \\ addpd %xmm14, %xmm2
    \\ movsd 8 * 8(%rbx, %rcx, 1), %xmm14
    \\ movhpd 9 * 8(%rbx, %rcx, 1), %xmm14
    \\ movapd %xmm1, %xmm15
    \\ movapd 10 * 8(%r13), %xmm1
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm7, %xmm12
    \\ addpd %xmm15, %xmm11
    \\ addpd %xmm12, %xmm3
    \\ movsd 8 * 8(%rbx), %xmm12
    \\ movhpd 9 * 8(%rbx), %xmm12
    \\ movsd %xmm2, 4 * 8(%r14)
    \\ movhpd %xmm2, 5 * 8(%r14)
    \\ movsd 8 * 8(%r14), %xmm2
    \\ movhpd 9 * 8(%r14), %xmm2
    \\ movsd %xmm3, 6 * 8(%r14)
    \\ movhpd %xmm3, 7 * 8(%r14)
    \\ movsd 10 * 8(%r14), %xmm3
    \\ movhpd 11 * 8(%r14), %xmm3
    \\ addq $8 * 8, %r13
    \\ addq $8 * 8, %r14
    \\ addq $8 * 8, %rbx
    \\ addq $8 * 8, %rbp
    \\ decq %rax
    \\ jg 20009b
    \\
;
