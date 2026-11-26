//! Section 2 of 4 of the text of `zgemv_d`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zgemv_t.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20011:
    \\ pcmpeqb %xmm13, %xmm13
    \\ psllq $63, %xmm13
    \\ shufps $0xc0, %xmm13, %xmm13
    \\ xorpd %xmm13, %xmm0
    \\ xorpd %xmm13, %xmm2
    \\ xorpd %xmm13, %xmm4
    \\ xorpd %xmm13, %xmm6
    \\ haddpd %xmm1, %xmm0
    \\ haddpd %xmm3, %xmm2
    \\ haddpd %xmm5, %xmm4
    \\ haddpd %xmm7, %xmm6
    \\ pshufd $0x4e, %xmm0, %xmm1
    \\ pshufd $0x4e, %xmm2, %xmm3
    \\ pshufd $0x4e, %xmm4, %xmm5
    \\ pshufd $0x4e, %xmm6, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ xorpd %xmm13, %xmm1
    \\ xorpd %xmm13, %xmm3
    \\ xorpd %xmm13, %xmm5
    \\ xorpd %xmm13, %xmm7
    \\ subpd %xmm1, %xmm0
    \\ subpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movsd 0 * 8(%rbp), %xmm1
    \\ movhpd 1 * 8(%rbp), %xmm1
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm3
    \\ movhpd 1 * 8(%rbp), %xmm3
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm5
    \\ movhpd 1 * 8(%rbp), %xmm5
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm7
    \\ movhpd 1 * 8(%rbp), %xmm7
    \\ addq %r10, %rbp
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ addpd %xmm5, %xmm4
    \\ addpd %xmm7, %xmm6
    \\ movlpd %xmm0, 0 * 8(%r14)
    \\ movhpd %xmm0, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm2, 0 * 8(%r14)
    \\ movhpd %xmm2, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm4, 0 * 8(%r14)
    \\ movhpd %xmm4, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm6, 0 * 8(%r14)
    \\ movhpd %xmm6, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ cmpq $4, %rsi
    \\ jge 20006b
    \\ .balign 8
    \\ 20012:
    \\ cmpq $2, %rsi
    \\ jl 20018f
    \\ subq $2, %rsi
    \\ leaq 16 * 8(%r15), %r13
    \\ movq %rcx, %r11
    \\ leaq (%r11, %r8), %r12
    \\ leaq (%r11, %r8, 2), %rcx
    \\ xorpd %xmm0, %xmm0
    \\ xorpd %xmm1, %xmm1
    \\ xorpd %xmm2, %xmm2
    \\ xorpd %xmm3, %xmm3
    \\ movaps -16 * 8(%r13), %xmm4
    \\ movaps -14 * 8(%r13), %xmm5
    \\ .if \bp0
    \\ prefetchw 3 * 8(%r14)
    \\ .endif
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20015f
    \\ movaps -16 * 8(%r11), %xmm8
    \\ movaps -16 * 8(%r12), %xmm10
    \\ movaps -14 * 8(%r11), %xmm12
    \\ movaps -14 * 8(%r12), %xmm6
    \\ decq %rax
    \\ jle 20014f
    \\ .balign 8
    \\ 20013:
    \\ .if \bp0
    \\ prefetch (64 * 4) * 2 - 128 + (64 >> 1)(%r11)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 (64 * 2) * 2 - 128 + (64 >> 1)(%r11)
    \\ .endif
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm4, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movaps -12 * 8(%r11), %xmm8
    \\ mulpd %xmm4, %xmm9
    \\ subpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm4, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movaps -12 * 8(%r12), %xmm10
    \\ mulpd %xmm4, %xmm11
    \\ subpd %xmm11, %xmm3
    \\ movaps -12 * 8(%r13), %xmm4
    \\ pshufd $0x4e, %xmm12, %xmm13
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm12, %xmm0
    \\ movaps -10 * 8(%r11), %xmm12
    \\ mulpd %xmm5, %xmm13
    \\ subpd %xmm13, %xmm1
    \\ pshufd $0x4e, %xmm6, %xmm7
    \\ mulpd %xmm5, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ movaps -10 * 8(%r12), %xmm6
    \\ mulpd %xmm5, %xmm7
    \\ subpd %xmm7, %xmm3
    \\ movaps -10 * 8(%r13), %xmm5
    \\ .if \bp0
    \\ prefetch (64 * 4) * 2 - 128 + (64 >> 1)(%r12)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 (64 * 2) * 2 - 128 + (64 >> 1)(%r12)
    \\ .endif
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm4, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movaps -8 * 8(%r11), %xmm8
    \\ mulpd %xmm4, %xmm9
    \\ subpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm4, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movaps -8 * 8(%r12), %xmm10
    \\ mulpd %xmm4, %xmm11
    \\ subpd %xmm11, %xmm3
    \\ movaps -8 * 8(%r13), %xmm4
    \\ .if \bp0
    \\ prefetch (64 * 4) * 2 - 128 + (64 >> 1)(%r13)
    \\ .endif
    \\ pshufd $0x4e, %xmm12, %xmm13
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm12, %xmm0
    \\ movaps -6 * 8(%r11), %xmm12
    \\ mulpd %xmm5, %xmm13
    \\ subpd %xmm13, %xmm1
    \\ pshufd $0x4e, %xmm6, %xmm7
    \\ mulpd %xmm5, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ movaps -6 * 8(%r12), %xmm6
    \\ mulpd %xmm5, %xmm7
    \\ subpd %xmm7, %xmm3
    \\ movaps -6 * 8(%r13), %xmm5
    \\ subq $-8 * 8, %r11
    \\ subq $-8 * 8, %r12
    \\ subq $-8 * 8, %r13
    \\ subq $1, %rax
    \\ .byte 0x3e
    \\ jg 20013b
    \\ .balign 8
    \\ 20014:
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm4, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movaps -12 * 8(%r11), %xmm8
    \\ mulpd %xmm4, %xmm9
    \\ subpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm4, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movaps -12 * 8(%r12), %xmm10
    \\ mulpd %xmm4, %xmm11
    \\ subpd %xmm11, %xmm3
    \\ movaps -12 * 8(%r13), %xmm4
    \\ pshufd $0x4e, %xmm12, %xmm13
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm12, %xmm0
    \\ movaps -10 * 8(%r11), %xmm12
    \\ mulpd %xmm5, %xmm13
    \\ subpd %xmm13, %xmm1
    \\ pshufd $0x4e, %xmm6, %xmm7
    \\ mulpd %xmm5, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ movaps -10 * 8(%r12), %xmm6
    \\ mulpd %xmm5, %xmm7
    \\ subpd %xmm7, %xmm3
    \\ movaps -10 * 8(%r13), %xmm5
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm4, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ mulpd %xmm4, %xmm9
    \\ subpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm4, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ mulpd %xmm4, %xmm11
    \\ subpd %xmm11, %xmm3
    \\ movaps -8 * 8(%r13), %xmm4
    \\ pshufd $0x4e, %xmm12, %xmm13
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm12, %xmm0
    \\ mulpd %xmm5, %xmm13
    \\ subpd %xmm13, %xmm1
    \\ pshufd $0x4e, %xmm6, %xmm7
    \\ mulpd %xmm5, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ mulpd %xmm5, %xmm7
    \\ subpd %xmm7, %xmm3
    \\ movaps -6 * 8(%r13), %xmm5
    \\ subq $-8 * 8, %r11
    \\ subq $-8 * 8, %r12
    \\ subq $-8 * 8, %r13
    \\ .balign 8
    \\ 20015:
    \\ testq $2, %rdi
    \\ je 20016f
    \\ movaps -16 * 8(%r11), %xmm8
    \\ movaps -16 * 8(%r12), %xmm10
    \\ movaps -14 * 8(%r11), %xmm12
    \\ movaps -14 * 8(%r12), %xmm6
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm4, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ mulpd %xmm4, %xmm9
    \\ subpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm4, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ mulpd %xmm4, %xmm11
    \\ subpd %xmm11, %xmm3
    \\ movaps -12 * 8(%r13), %xmm4
    \\ pshufd $0x4e, %xmm12, %xmm13
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm12, %xmm0
    \\ mulpd %xmm5, %xmm13
    \\ subpd %xmm13, %xmm1
    \\ pshufd $0x4e, %xmm6, %xmm7
    \\ mulpd %xmm5, %xmm6
    \\ addpd %xmm6, %xmm2
    \\ mulpd %xmm5, %xmm7
    \\ subpd %xmm7, %xmm3
    \\ addq $4 * 8, %r11
    \\ addq $4 * 8, %r12
    \\ .balign 8
    \\ 20016:
    \\ testq $1, %rdi
    \\ je 20017f
    \\ movaps -16 * 8(%r11), %xmm8
    \\ movaps -16 * 8(%r12), %xmm10
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm4, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ mulpd %xmm4, %xmm9
    \\ subpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm4, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ mulpd %xmm4, %xmm11
    \\ subpd %xmm11, %xmm3
    \\ .balign 8
    \\ 20017:
    \\ pcmpeqb %xmm11, %xmm11
    \\ psllq $63, %xmm11
    \\ shufps $0xc0, %xmm11, %xmm11
    \\ xorpd %xmm11, %xmm0
    \\ xorpd %xmm11, %xmm2
    \\ haddpd %xmm1, %xmm0
    \\ haddpd %xmm3, %xmm2
    \\ pshufd $0x4e, %xmm0, %xmm1
    \\ pshufd $0x4e, %xmm2, %xmm3
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ xorpd %xmm11, %xmm1
    \\ xorpd %xmm11, %xmm3
    \\ subpd %xmm1, %xmm0
    \\ subpd %xmm3, %xmm2
    \\ movsd 0 * 8(%rbp), %xmm4
    \\ movhpd 1 * 8(%rbp), %xmm4
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm5
    \\ movhpd 1 * 8(%rbp), %xmm5
    \\ addq %r10, %rbp
    \\ addpd %xmm4, %xmm0
    \\ addpd %xmm5, %xmm2
    \\ movlpd %xmm0, 0 * 8(%r14)
    \\ movhpd %xmm0, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm2, 0 * 8(%r14)
    \\ movhpd %xmm2, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ .balign 8
    \\ 20018:
    \\ cmpq $1, %rsi
    \\ jl 20043f
    \\ leaq 16 * 8(%r15), %r13
    \\ movq %rcx, %r11
    \\ xorpd %xmm0, %xmm0
    \\ xorpd %xmm1, %xmm1
    \\ movaps -16 * 8(%r13), %xmm4
    \\ movaps -14 * 8(%r13), %xmm5
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20021f
    \\ movaps -16 * 8(%r11), %xmm8
    \\ movaps -14 * 8(%r11), %xmm12
    \\ decq %rax
    \\ jle 20020f
    \\ .balign 8
    \\ 20019:
    \\ .if \bp0
    \\ prefetch (64 * 4) * 4 - 128 + (64 >> 1)(%r11)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 (64 * 2) * 4 - 128 + (64 >> 1)(%r11)
    \\ .endif
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm4, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movaps -12 * 8(%r11), %xmm8
    \\ mulpd %xmm4, %xmm9
    \\ subpd %xmm9, %xmm1
    \\ movaps -12 * 8(%r13), %xmm4
    \\ pshufd $0x4e, %xmm12, %xmm13
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm12, %xmm0
    \\ movaps -10 * 8(%r11), %xmm12
    \\ mulpd %xmm5, %xmm13
    \\ subpd %xmm13, %xmm1
    \\ movaps -10 * 8(%r13), %xmm5
    \\ .if \bp0
    \\ prefetch (64 * 4) * 4 - 128 + (64 >> 1)(%r13)
    \\ .endif
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm4, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movaps -8 * 8(%r11), %xmm8
    \\ mulpd %xmm4, %xmm9
    \\ subpd %xmm9, %xmm1
    \\ movaps -8 * 8(%r13), %xmm4
    \\ pshufd $0x4e, %xmm12, %xmm13
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm12, %xmm0
    \\ movaps -6 * 8(%r11), %xmm12
    \\ mulpd %xmm5, %xmm13
    \\ subpd %xmm13, %xmm1
    \\ movaps -6 * 8(%r13), %xmm5
    \\ subq $-8 * 8, %r11
    \\ subq $-8 * 8, %r13
    \\ subq $1, %rax
    \\ .byte 0x3e
    \\ jg 20019b
    \\ .balign 8
    \\ 20020:
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm4, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movaps -12 * 8(%r11), %xmm8
    \\ mulpd %xmm4, %xmm9
    \\ subpd %xmm9, %xmm1
    \\ movaps -12 * 8(%r13), %xmm4
    \\ pshufd $0x4e, %xmm12, %xmm13
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm12, %xmm0
    \\ movaps -10 * 8(%r11), %xmm12
    \\ mulpd %xmm5, %xmm13
    \\ subpd %xmm13, %xmm1
    \\ movaps -10 * 8(%r13), %xmm5
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm4, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ mulpd %xmm4, %xmm9
    \\ subpd %xmm9, %xmm1
    \\ movaps -8 * 8(%r13), %xmm4
    \\ pshufd $0x4e, %xmm12, %xmm13
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm12, %xmm0
    \\ mulpd %xmm5, %xmm13
    \\ subpd %xmm13, %xmm1
    \\ movaps -6 * 8(%r13), %xmm5
    \\ subq $-8 * 8, %r11
    \\ subq $-8 * 8, %r13
    \\ .balign 8
    \\ 20021:
    \\ testq $2, %rdi
    \\ je 20022f
    \\ movaps -16 * 8(%r11), %xmm8
    \\ movaps -14 * 8(%r11), %xmm12
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm4, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ mulpd %xmm4, %xmm9
    \\ subpd %xmm9, %xmm1
    \\ movaps -12 * 8(%r13), %xmm4
    \\ pshufd $0x4e, %xmm12, %xmm13
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm12, %xmm0
    \\ mulpd %xmm5, %xmm13
    \\ subpd %xmm13, %xmm1
    \\ addq $4 * 8, %r11
    \\ .balign 8
    \\ 20022:
    \\ testq $1, %rdi
    \\ je 20023f
    \\ movaps -16 * 8(%r11), %xmm8
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm4, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ mulpd %xmm4, %xmm9
    \\ subpd %xmm9, %xmm1
    \\ .balign 8
    \\ 20023:
    \\ pcmpeqb %xmm11, %xmm11
    \\ psllq $63, %xmm11
    \\ shufps $0xc0, %xmm11, %xmm11
    \\ xorpd %xmm11, %xmm0
    \\ haddpd %xmm1, %xmm0
    \\ pshufd $0x4e, %xmm0, %xmm1
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ xorpd %xmm11, %xmm1
    \\ subpd %xmm1, %xmm0
    \\ movsd 0 * 8(%rbp), %xmm4
    \\ movhpd 1 * 8(%rbp), %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movlpd %xmm0, 0 * 8(%r14)
    \\ movhpd %xmm0, 1 * 8(%r14)
    \\ jmp 20043f
    \\ .balign 8
    \\ 20024:
    \\ cmpq $4, %rsi
    \\ jl 20031f
    \\
;
