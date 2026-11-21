//! Section 3 of 4 of the text of `zgemv_t`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zgemv_t.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20025:
    \\ subq $4, %rsi
    \\ leaq 16 * 8(%r15), %r13
    \\ movq %rcx, %r11
    \\ leaq (%r11, %r8, 2), %r12
    \\ leaq (%r11, %r8, 4), %rcx
    \\ movaps -16 * 8(%r13), %xmm12
    \\ xorpd %xmm0, %xmm0
    \\ xorpd %xmm1, %xmm1
    \\ xorpd %xmm2, %xmm2
    \\ xorpd %xmm3, %xmm3
    \\ movaps -14 * 8(%r13), %xmm13
    \\ xorpd %xmm4, %xmm4
    \\ xorpd %xmm5, %xmm5
    \\ xorpd %xmm6, %xmm6
    \\ xorpd %xmm7, %xmm7
    \\ .if \bp0
    \\ prefetchw 3 * 8(%r14)
    \\ .endif
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20028f
    \\ movsd -16 * 8(%r11), %xmm8
    \\ movhpd -15 * 8(%r11), %xmm8
    \\ movsd -16 * 8(%r11, %r8), %xmm10
    \\ movhpd -15 * 8(%r11, %r8), %xmm10
    \\ decq %rax
    \\ jle 20027f
    \\ .balign 8
    \\ 20026:
    \\ .if \bp0
    \\ prefetch (64 * 4) - 128 + (64 >> 1)(%r11)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 (64 * 2) - 128 + (64 >> 1)(%r11)
    \\ .endif
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movsd -16 * 8(%r12), %xmm8
    \\ movhpd -15 * 8(%r12), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movsd -16 * 8(%r12, %r8), %xmm10
    \\ movhpd -15 * 8(%r12, %r8), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movsd -14 * 8(%r11), %xmm8
    \\ movhpd -13 * 8(%r11), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movsd -14 * 8(%r11, %r8), %xmm10
    \\ movhpd -13 * 8(%r11, %r8), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ movaps -12 * 8(%r13), %xmm12
    \\ addpd %xmm11, %xmm7
    \\ .if \bp0
    \\ prefetch (64 * 4) - 128 + (64 >> 1)(%r11, %r8)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 (64 * 2) - 128 + (64 >> 1)(%r11, %r8)
    \\ .endif
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movsd -14 * 8(%r12), %xmm8
    \\ movhpd -13 * 8(%r12), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movsd -14 * 8(%r12, %r8), %xmm10
    \\ movhpd -13 * 8(%r12, %r8), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movsd -12 * 8(%r11), %xmm8
    \\ movhpd -11 * 8(%r11), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movsd -12 * 8(%r11, %r8), %xmm10
    \\ movhpd -11 * 8(%r11, %r8), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ movaps -10 * 8(%r13), %xmm13
    \\ addpd %xmm11, %xmm7
    \\ .if \bp0
    \\ prefetch (64 * 4) - 128 + (64 >> 1)(%r12)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 (64 * 2) - 128 + (64 >> 1)(%r12)
    \\ .endif
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movsd -12 * 8(%r12), %xmm8
    \\ movhpd -11 * 8(%r12), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movsd -12 * 8(%r12, %r8), %xmm10
    \\ movhpd -11 * 8(%r12, %r8), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movsd -10 * 8(%r11), %xmm8
    \\ movhpd -9 * 8(%r11), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movsd -10 * 8(%r11, %r8), %xmm10
    \\ movhpd -9 * 8(%r11, %r8), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ movaps -8 * 8(%r13), %xmm12
    \\ addpd %xmm11, %xmm7
    \\ .if \bp0
    \\ prefetch (64 * 4) - 128 + (64 >> 1)(%r12, %r8)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 (64 * 2) - 128 + (64 >> 1)(%r12, %r8)
    \\ .endif
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movsd -10 * 8(%r12), %xmm8
    \\ movhpd -9 * 8(%r12), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movsd -10 * 8(%r12, %r8), %xmm10
    \\ movhpd -9 * 8(%r12, %r8), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ .if \bp0
    \\ prefetch (64 * 4) - 128 + (64 >> 1)(%r13)
    \\ .endif
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movsd -8 * 8(%r11), %xmm8
    \\ movhpd -7 * 8(%r11), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movsd -8 * 8(%r11, %r8), %xmm10
    \\ movhpd -7 * 8(%r11, %r8), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ movaps -6 * 8(%r13), %xmm13
    \\ addpd %xmm11, %xmm7
    \\ subq $-8 * 8, %r11
    \\ subq $-8 * 8, %r12
    \\ subq $-8 * 8, %r13
    \\ subq $1, %rax
    \\ .byte 0x3e
    \\ jg 20026b
    \\ .balign 8
    \\ 20027:
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movsd -16 * 8(%r12), %xmm8
    \\ movhpd -15 * 8(%r12), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movsd -16 * 8(%r12, %r8), %xmm10
    \\ movhpd -15 * 8(%r12, %r8), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movsd -14 * 8(%r11), %xmm8
    \\ movhpd -13 * 8(%r11), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movsd -14 * 8(%r11, %r8), %xmm10
    \\ movhpd -13 * 8(%r11, %r8), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ movaps -12 * 8(%r13), %xmm12
    \\ addpd %xmm11, %xmm7
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movsd -14 * 8(%r12), %xmm8
    \\ movhpd -13 * 8(%r12), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movsd -14 * 8(%r12, %r8), %xmm10
    \\ movhpd -13 * 8(%r12, %r8), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movsd -12 * 8(%r11), %xmm8
    \\ movhpd -11 * 8(%r11), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movsd -12 * 8(%r11, %r8), %xmm10
    \\ movhpd -11 * 8(%r11, %r8), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ movaps -10 * 8(%r13), %xmm13
    \\ addpd %xmm11, %xmm7
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movsd -12 * 8(%r12), %xmm8
    \\ movhpd -11 * 8(%r12), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movsd -12 * 8(%r12, %r8), %xmm10
    \\ movhpd -11 * 8(%r12, %r8), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movsd -10 * 8(%r11), %xmm8
    \\ movhpd -9 * 8(%r11), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movsd -10 * 8(%r11, %r8), %xmm10
    \\ movhpd -9 * 8(%r11, %r8), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ movaps -8 * 8(%r13), %xmm12
    \\ addpd %xmm11, %xmm7
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movsd -10 * 8(%r12), %xmm8
    \\ movhpd -9 * 8(%r12), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movsd -10 * 8(%r12, %r8), %xmm10
    \\ movhpd -9 * 8(%r12, %r8), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ mulpd %xmm13, %xmm11
    \\ movaps -6 * 8(%r13), %xmm13
    \\ addpd %xmm11, %xmm7
    \\ subq $-8 * 8, %r11
    \\ subq $-8 * 8, %r12
    \\ subq $-8 * 8, %r13
    \\ .balign 8
    \\ 20028:
    \\ testq $2, %rdi
    \\ je 20029f
    \\ movsd -16 * 8(%r11), %xmm8
    \\ movhpd -15 * 8(%r11), %xmm8
    \\ movsd -16 * 8(%r11, %r8), %xmm10
    \\ movhpd -15 * 8(%r11, %r8), %xmm10
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movsd -16 * 8(%r12), %xmm8
    \\ movhpd -15 * 8(%r12), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movsd -16 * 8(%r12, %r8), %xmm10
    \\ movhpd -15 * 8(%r12, %r8), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movsd -14 * 8(%r11), %xmm8
    \\ movhpd -13 * 8(%r11), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movsd -14 * 8(%r11, %r8), %xmm10
    \\ movhpd -13 * 8(%r11, %r8), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ movaps -12 * 8(%r13), %xmm12
    \\ addpd %xmm11, %xmm7
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movsd -14 * 8(%r12), %xmm8
    \\ movhpd -13 * 8(%r12), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movsd -14 * 8(%r12, %r8), %xmm10
    \\ movhpd -13 * 8(%r12, %r8), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ mulpd %xmm13, %xmm11
    \\ addpd %xmm11, %xmm7
    \\ addq $4 * 8, %r11
    \\ addq $4 * 8, %r12
    \\ .balign 8
    \\ 20029:
    \\ testq $1, %rdi
    \\ je 20030f
    \\ movsd -16 * 8(%r11), %xmm8
    \\ movhpd -15 * 8(%r11), %xmm8
    \\ movsd -16 * 8(%r11, %r8), %xmm10
    \\ movhpd -15 * 8(%r11, %r8), %xmm10
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movsd -16 * 8(%r12), %xmm8
    \\ movhpd -15 * 8(%r12), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movsd -16 * 8(%r12, %r8), %xmm10
    \\ movhpd -15 * 8(%r12, %r8), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ mulpd %xmm12, %xmm11
    \\ addpd %xmm11, %xmm7
    \\ .balign 8
    \\ 20030:
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
    \\ jge 20025b
    \\
;
