//! Section 1 of 2 of the text of `zsymv_U`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zsymv_U_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1, bp2, bp3
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
    \\ salq $4, %r9
    \\ salq $4, %r11
    \\ salq $4, %rcx
    \\ testq %rdi, %rdi
    \\ jle 20019f
    \\ negq %rsi
    \\ addq %rdi, %rsi
    \\ movq %rsi, %rax
    \\ imulq %rcx, %rax
    \\ addq %rax, %rdx
    \\ pcmpeqb %xmm2, %xmm2
    \\ xorpd %xmm3, %xmm3
    \\ psllq $63, %xmm2
    \\ unpcklpd %xmm3, %xmm2
    \\ unpcklpd %xmm1, %xmm0
    \\ unpcklpd %xmm0, %xmm1
    \\ xorpd %xmm2, %xmm1
    \\ movq %r12, %r13
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20001f
    \\ .balign 8
    \\ 20000:
    \\ movddup 0 * 8(%r8), %xmm3
    \\ movddup 1 * 8(%r8), %xmm4
    \\ addq %r9, %r8
    \\ movddup 0 * 8(%r8), %xmm5
    \\ movddup 1 * 8(%r8), %xmm6
    \\ addq %r9, %r8
    \\ mulpd %xmm0, %xmm3
    \\ mulpd %xmm1, %xmm4
    \\ mulpd %xmm0, %xmm5
    \\ mulpd %xmm1, %xmm6
    \\ addpd %xmm4, %xmm3
    \\ addpd %xmm6, %xmm5
    \\ movapd %xmm3, 0 * 8(%r13)
    \\ shufpd $1, %xmm3, %xmm3
    \\ pxor %xmm2, %xmm3
    \\ movapd %xmm3, 2 * 8(%r13)
    \\ movapd %xmm5, 4 * 8(%r13)
    \\ shufpd $1, %xmm5, %xmm5
    \\ pxor %xmm2, %xmm5
    \\ movapd %xmm5, 6 * 8(%r13)
    \\ movddup 0 * 8(%r8), %xmm3
    \\ movddup 1 * 8(%r8), %xmm4
    \\ addq %r9, %r8
    \\ movddup 0 * 8(%r8), %xmm5
    \\ movddup 1 * 8(%r8), %xmm6
    \\ addq %r9, %r8
    \\ mulpd %xmm0, %xmm3
    \\ mulpd %xmm1, %xmm4
    \\ mulpd %xmm0, %xmm5
    \\ mulpd %xmm1, %xmm6
    \\ addpd %xmm4, %xmm3
    \\ addpd %xmm6, %xmm5
    \\ movapd %xmm3, 8 * 8(%r13)
    \\ shufpd $1, %xmm3, %xmm3
    \\ pxor %xmm2, %xmm3
    \\ movapd %xmm3, 10 * 8(%r13)
    \\ movapd %xmm5, 12 * 8(%r13)
    \\ shufpd $1, %xmm5, %xmm5
    \\ pxor %xmm2, %xmm5
    \\ movapd %xmm5, 14 * 8(%r13)
    \\ subq $-16 * 8, %r13
    \\ decq %rax
    \\ jg 20000b
    \\ .balign 8
    \\ 20001:
    \\ movq %rdi, %rax
    \\ andq $3, %rax
    \\ jle 20003f
    \\ .balign 8
    \\ 20002:
    \\ movddup 0 * 8(%r8), %xmm3
    \\ movddup 1 * 8(%r8), %xmm4
    \\ addq %r9, %r8
    \\ mulpd %xmm0, %xmm3
    \\ mulpd %xmm1, %xmm4
    \\ addpd %xmm4, %xmm3
    \\ movapd %xmm3, 0 * 8(%r13)
    \\ shufpd $1, %xmm3, %xmm3
    \\ pxor %xmm2, %xmm3
    \\ movapd %xmm3, 2 * 8(%r13)
    \\ addq $4 * 8, %r13
    \\ decq %rax
    \\ jg 20002b
    \\ .balign 8
    \\ 20003:
    \\ movq %r10, %r8
    \\ addq $512, %r13
    \\ andq $-512, %r13
    \\ cmpq $2 * 8, %r11
    \\ je 20007f
    \\ movq %r10, %r14
    \\ movq %r13, %r8
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20005f
    \\ .balign 8
    \\ 20004:
    \\ movsd 0 * 8(%r14), %xmm0
    \\ movhpd 1 * 8(%r14), %xmm0
    \\ addq %r11, %r14
    \\ movsd 0 * 8(%r14), %xmm1
    \\ movhpd 1 * 8(%r14), %xmm1
    \\ addq %r11, %r14
    \\ movsd 0 * 8(%r14), %xmm2
    \\ movhpd 1 * 8(%r14), %xmm2
    \\ addq %r11, %r14
    \\ movsd 0 * 8(%r14), %xmm3
    \\ movhpd 1 * 8(%r14), %xmm3
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
    \\ andq $3, %rax
    \\ jle 20007f
    \\ .balign 8
    \\ 20006:
    \\ movsd 0 * 8(%r14), %xmm0
    \\ movhpd 1 * 8(%r14), %xmm0
    \\ addq %r11, %r14
    \\ movapd %xmm0, 0 * 8(%r13)
    \\ addq $2 * 8, %r13
    \\ decq %rax
    \\ jg 20006b
    \\ .balign 8
    \\ 20007:
    \\ movq %rsi, %rax
    \\ addq $2, %rax
    \\ cmpq %rdi, %rax
    \\ jg 20012f
    \\ .balign 8
    \\ 20008:
    \\ movq %rdx, %rbx
    \\ leaq (%rdx, %rcx, 1), %rbp
    \\ leaq (%rdx, %rcx, 2), %rdx
    \\ leaq (, %rsi, 4), %rax
    \\ movapd 0 * 8(%r12, %rax, 8), %xmm4
    \\ movapd 2 * 8(%r12, %rax, 8), %xmm5
    \\ movapd 4 * 8(%r12, %rax, 8), %xmm6
    \\ movapd 6 * 8(%r12, %rax, 8), %xmm7
    \\ pxor %xmm8, %xmm8
    \\ pxor %xmm9, %xmm9
    \\ movsd 0 * 8(%r8), %xmm10
    \\ movhpd 1 * 8(%r8), %xmm10
    \\ movsd 2 * 8(%r8), %xmm11
    \\ movhpd 3 * 8(%r8), %xmm11
    \\ movapd 0 * 8(%r12), %xmm0
    \\ movapd 2 * 8(%r12), %xmm1
    \\ movapd 4 * 8(%r12), %xmm2
    \\ movapd 6 * 8(%r12), %xmm3
    \\ movddup 0 * 8(%rbx), %xmm12
    \\ movddup 2 * 8(%rbp), %xmm13
    \\ movddup 1 * 8(%rbx), %xmm14
    \\ movq %r12, %r13
    \\ movq %r8, %r14
    \\ movq %rsi, %rax
    \\ sarq $2, %rax
    \\ jle 20010f
    \\ .balign 8
    \\ 20009:
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm4, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm10
    \\ movddup 3 * 8(%rbp), %xmm12
    \\ .if \bp0
    \\ prefetcht0 (16 * 24)(%rbx)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 16)(%rbx)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 (8 * 24)(%rbx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 (16 * 28)(%rbx)
    \\ .endif
    \\ movapd %xmm2, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm6, %xmm13
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm13, %xmm11
    \\ movddup 2 * 8(%rbx), %xmm13
    \\ movapd %xmm1, %xmm15
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm5, %xmm14
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm14, %xmm10
    \\ movddup 0 * 8(%rbp), %xmm14
    \\ movapd %xmm3, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm7, %xmm12
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm12, %xmm11
    \\ movddup 3 * 8(%rbx), %xmm12
    \\ .if \bp0
    \\ prefetcht0 (16 * 24)(%r13)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 16)(%r13)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 (8 * 24)(%r13)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 (16 * 28)(%r13)
    \\ .endif
    \\ movapd %xmm2, %xmm15
    \\ movapd 12 * 8(%r13), %xmm2
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm4, %xmm13
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm13, %xmm11
    \\ movddup 1 * 8(%rbp), %xmm13
    \\ movapd %xmm0, %xmm15
    \\ movapd 8 * 8(%r13), %xmm0
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm6, %xmm14
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm14, %xmm10
    \\ movddup 4 * 8(%rbx), %xmm14
    \\ movapd %xmm3, %xmm15
    \\ movapd 14 * 8(%r13), %xmm3
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm11
    \\ movddup 6 * 8(%rbp), %xmm12
    \\ movlpd %xmm11, 2 * 8(%r14)
    \\ movhpd %xmm11, 3 * 8(%r14)
    \\ movsd 6 * 8(%r14), %xmm11
    \\ movhpd 7 * 8(%r14), %xmm11
    \\ movapd %xmm1, %xmm15
    \\ movapd 10 * 8(%r13), %xmm1
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm7, %xmm13
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm13, %xmm10
    \\ movddup 5 * 8(%rbx), %xmm13
    \\ .if \bp0
    \\ prefetcht0 (16 * 24)(%rbp)
    \\ .endif
    \\ .if \bp1
    \\ prefetch (16 * 16)(%rbp)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 (8 * 24)(%rbp)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 (16 * 28)(%rbp)
    \\ .endif
    \\ movlpd %xmm10, 0 * 8(%r14)
    \\ movhpd %xmm10, 1 * 8(%r14)
    \\ movsd 4 * 8(%r14), %xmm10
    \\ movhpd 5 * 8(%r14), %xmm10
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm4, %xmm14
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm14, %xmm10
    \\ movddup 7 * 8(%rbp), %xmm14
    \\ movapd %xmm2, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm6, %xmm12
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm12, %xmm11
    \\ movddup 6 * 8(%rbx), %xmm12
    \\ movapd %xmm1, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm5, %xmm13
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm13, %xmm10
    \\ movddup 4 * 8(%rbp), %xmm13
    \\ .if \bp0
    \\ prefetcht0 (16 * 24)(%r14)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw (16 * 16)(%r14)
    \\ .endif
    \\ .if \bp2
    \\ prefetcht0 (8 * 24)(%r14)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 (16 * 28)(%r14)
    \\ .endif
    \\ movapd %xmm3, %xmm15
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm7, %xmm14
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm14, %xmm11
    \\ movddup 7 * 8(%rbx), %xmm14
    \\ movapd %xmm2, %xmm15
    \\ movapd 20 * 8(%r13), %xmm2
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm4, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm11
    \\ movddup 5 * 8(%rbp), %xmm12
    \\ movapd %xmm0, %xmm15
    \\ movapd 16 * 8(%r13), %xmm0
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm6, %xmm13
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm13, %xmm10
    \\ movddup 10 * 8(%rbp), %xmm13
    \\ movapd %xmm3, %xmm15
    \\ movapd 22 * 8(%r13), %xmm3
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm5, %xmm14
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm14, %xmm11
    \\ movddup 9 * 8(%rbx), %xmm14
    \\ movlpd %xmm11, 6 * 8(%r14)
    \\ movhpd %xmm11, 7 * 8(%r14)
    \\ movsd 10 * 8(%r14), %xmm11
    \\ movhpd 11 * 8(%r14), %xmm11
    \\ movapd %xmm1, %xmm15
    \\ movapd 18 * 8(%r13), %xmm1
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm7, %xmm12
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm12, %xmm10
    \\ movddup 8 * 8(%rbx), %xmm12
    \\ movlpd %xmm10, 4 * 8(%r14)
    \\ movhpd %xmm10, 5 * 8(%r14)
    \\ movsd 8 * 8(%r14), %xmm10
    \\ movhpd 9 * 8(%r14), %xmm10
    \\ subq $-16 * 8, %r13
    \\ addq $ 8 * 8, %r14
    \\ addq $ 8 * 8, %rbx
    \\ addq $ 8 * 8, %rbp
    \\ decq %rax
    \\ jg 20009b
    \\
;
