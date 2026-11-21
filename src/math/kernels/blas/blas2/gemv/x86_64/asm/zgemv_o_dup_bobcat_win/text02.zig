//! Section 2 of 3 of the text of `zgemv_o`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zgemv_n_dup.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20005:
    \\ testq $2, %rcx
    \\ je 20006f
    \\ movddup -16 * 8(%r12), %xmm4
    \\ movddup -15 * 8(%r12), %xmm5
    \\ movddup -14 * 8(%r12), %xmm6
    \\ movddup -13 * 8(%r12), %xmm7
    \\ mulpd %xmm8, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -16 * 8(%r12, %r9, 1), %xmm4
    \\ mulpd %xmm8, %xmm6
    \\ addpd %xmm6, %xmm1
    \\ movddup -14 * 8(%r12, %r9, 1), %xmm6
    \\ mulpd %xmm9, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ movddup -15 * 8(%r12, %r9, 1), %xmm5
    \\ mulpd %xmm9, %xmm7
    \\ addpd %xmm7, %xmm1
    \\ movddup -13 * 8(%r12, %r9, 1), %xmm7
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
    \\ movddup -16 * 8(%r13, %r9, 1), %xmm4
    \\ mulpd %xmm12, %xmm6
    \\ addpd %xmm6, %xmm1
    \\ movddup -14 * 8(%r13, %r9, 1), %xmm6
    \\ mulpd %xmm13, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ movddup -15 * 8(%r13, %r9, 1), %xmm5
    \\ mulpd %xmm13, %xmm7
    \\ addpd %xmm7, %xmm1
    \\ movddup -13 * 8(%r13, %r9, 1), %xmm7
    \\ mulpd %xmm14, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ mulpd %xmm14, %xmm6
    \\ addpd %xmm6, %xmm1
    \\ mulpd %xmm15, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm7, %xmm1
    \\ movaps %xmm0, -16 * 8(%r14)
    \\ movaps %xmm1, -14 * 8(%r14)
    \\ movaps %xmm2, %xmm0
    \\ addq $4 * 8, %r12
    \\ addq $4 * 8, %r13
    \\ addq $4 * 8, %r14
    \\ .balign 8
    \\ 20006:
    \\ testq $1, %rcx
    \\ je 20007f
    \\ movddup -16 * 8(%r12), %xmm4
    \\ movddup -15 * 8(%r12), %xmm5
    \\ movddup -16 * 8(%r12, %r9, 1), %xmm6
    \\ movddup -15 * 8(%r12, %r9, 1), %xmm7
    \\ mulpd %xmm8, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ movddup -16 * 8(%r13), %xmm4
    \\ mulpd %xmm9, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ movddup -15 * 8(%r13), %xmm5
    \\ mulpd %xmm10, %xmm6
    \\ addpd %xmm6, %xmm0
    \\ movddup -16 * 8(%r13, %r9, 1), %xmm6
    \\ mulpd %xmm11, %xmm7
    \\ addpd %xmm7, %xmm0
    \\ movddup -15 * 8(%r13, %r9, 1), %xmm7
    \\ mulpd %xmm12, %xmm4
    \\ addpd %xmm4, %xmm0
    \\ mulpd %xmm13, %xmm5
    \\ addpd %xmm5, %xmm0
    \\ mulpd %xmm14, %xmm6
    \\ addpd %xmm6, %xmm0
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm7, %xmm0
    \\ movaps %xmm0, -16 * 8(%r14)
    \\ .balign 8
    \\ 20007:
    \\ cmpq $4, %rdx
    \\ jge 20002b
    \\ .balign 8
    \\ 20008:
    \\ cmpq $2, %rdx
    \\ jl 20013f
    \\ subq $2, %rdx
    \\ leaq 16 * 8(%r15), %r14
    \\ movq %r8, %r12
    \\ leaq (%r8, %r9, 1), %r13
    \\ leaq (%r8, %r9, 2), %r8
    \\ movddup 0 * 8(%rdi), %xmm8
    \\ movddup 1 * 8(%rdi), %xmm9
    \\ addq %rsi, %rdi
    \\ movddup 0 * 8(%rdi), %xmm10
    \\ movddup 1 * 8(%rdi), %xmm11
    \\ addq %rsi, %rdi
    \\ pcmpeqb %xmm5, %xmm5
    \\ psllq $63, %xmm5
    \\ shufps $0x40, %xmm5, %xmm5
    \\ movsd 224 (%rsp), %xmm6
    \\ movhps 232 (%rsp), %xmm6
    \\ pshufd $0x4e, %xmm6, %xmm7
    \\ xorps %xmm5, %xmm6
    \\ mulpd %xmm6, %xmm8
    \\ mulpd %xmm7, %xmm9
    \\ mulpd %xmm6, %xmm10
    \\ mulpd %xmm7, %xmm11
    \\ addpd %xmm9, %xmm8
    \\ addpd %xmm11, %xmm10
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ xorps %xmm5, %xmm8
    \\ xorps %xmm5, %xmm10
    \\ movaps -16 * 8(%r14), %xmm0
    \\ movaps -14 * 8(%r14), %xmm1
    \\ movaps -12 * 8(%r14), %xmm2
    \\ movaps -10 * 8(%r14), %xmm3
    \\ movq %rcx, %rax
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
    \\ .if \bp0
    \\ prefetch (64 * 4) * 2 - 128 + (64 >> 1)(%r13)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 (64 * 2) * 2 - 128 + (64 >> 1)(%r13)
    \\ .endif
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
    \\ .if \bp0
    \\ prefetch (64 * 4) * 2 - 128 + (64 >> 1)(%r12)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 (64 * 2) * 2 - 128 + (64 >> 1)(%r12)
    \\ .endif
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
    \\ .if \bp0
    \\ prefetchw (64 * 4) * 2 - 128 + (64 >> 1)(%r14)
    \\ .endif
    \\ movaps %xmm0, -16 * 8(%r14)
    \\ movaps %xmm1, -14 * 8(%r14)
    \\ movaps %xmm2, -12 * 8(%r14)
    \\ movaps %xmm3, -10 * 8(%r14)
    \\ movaps -8 * 8(%r14), %xmm0
    \\ movaps -6 * 8(%r14), %xmm1
    \\ movaps -4 * 8(%r14), %xmm2
    \\ movaps -2 * 8(%r14), %xmm3
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
    \\ movaps %xmm0, -16 * 8(%r14)
    \\ movaps %xmm1, -14 * 8(%r14)
    \\ movaps %xmm2, -12 * 8(%r14)
    \\ movaps %xmm3, -10 * 8(%r14)
    \\ movaps -8 * 8(%r14), %xmm0
    \\ movaps -6 * 8(%r14), %xmm1
    \\ movaps -4 * 8(%r14), %xmm2
    \\ movaps -2 * 8(%r14), %xmm3
    \\ subq $-8 * 8, %r12
    \\ subq $-8 * 8, %r13
    \\ subq $-8 * 8, %r14
    \\ .balign 8
    \\ 20011:
    \\ testq $2, %rcx
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
    \\ movaps %xmm0, -16 * 8(%r14)
    \\ movaps %xmm1, -14 * 8(%r14)
    \\ movaps %xmm2, %xmm0
    \\ addq $4 * 8, %r12
    \\ addq $4 * 8, %r13
    \\ addq $4 * 8, %r14
    \\ .balign 8
    \\ 20012:
    \\ testq $1, %rcx
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
    \\ movaps %xmm0, -16 * 8(%r14)
    \\ .balign 8
    \\ 20013:
    \\ cmpq $1, %rdx
    \\ jl 20018f
    \\ leaq 16 * 8(%r15), %r14
    \\ movq %r8, %r12
    \\ movddup 0 * 8(%rdi), %xmm8
    \\ movddup 1 * 8(%rdi), %xmm9
    \\ addq %rsi, %rdi
    \\ pcmpeqb %xmm5, %xmm5
    \\ psllq $63, %xmm5
    \\ shufps $0x40, %xmm5, %xmm5
    \\ movsd 224 (%rsp), %xmm6
    \\ movhps 232 (%rsp), %xmm6
    \\ pshufd $0x4e, %xmm6, %xmm7
    \\ xorps %xmm5, %xmm6
    \\ mulpd %xmm6, %xmm8
    \\ mulpd %xmm7, %xmm9
    \\ addpd %xmm9, %xmm8
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ xorps %xmm5, %xmm8
    \\ movaps -16 * 8(%r14), %xmm0
    \\ movaps -14 * 8(%r14), %xmm1
    \\ movaps -12 * 8(%r14), %xmm2
    \\ movaps -10 * 8(%r14), %xmm3
    \\ movq %rcx, %rax
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
