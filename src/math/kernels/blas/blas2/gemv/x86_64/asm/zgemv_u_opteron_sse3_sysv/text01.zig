//! Section 1 of 2 of the text of `zgemv_u`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zgemv_t.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 512
    \\ subq $128, %rsp
    \\ movq %rbx, 0(%rsp)
    \\ movq %rbp, 8(%rsp)
    \\ movq %r12, 16(%rsp)
    \\ movq %r13, 24(%rsp)
    \\ movq %r14, 32(%rsp)
    \\ movq %r15, 40(%rsp)
    \\ movq %rcx, 80(%rsp)
    \\ movq %rsi, 72(%rsp)
    \\ movq %rdi, 64(%rsp)
    \\ movq %r8, 88(%rsp)
    \\ movlps %xmm0,96(%rsp)
    \\ movlps %xmm1,104(%rsp)
    \\ 20000:
    \\ xorq %rax,%rax
    \\ addq $1,%rax
    \\ salq $19,%rax
    \\ subq %rax,64(%rsp)
    \\ movq %rax,%rdi
    \\ movlps 96(%rsp),%xmm0
    \\ movlps 104(%rsp),%xmm1
    \\ jge 20001f
    \\ movq 64(%rsp),%rdi
    \\ addq %rax,%rdi
    \\ jle 20025f
    \\ 20001:
    \\ movq 80(%rsp), %rcx
    \\ movq 72(%rsp), %rsi
    \\ movq 88(%rsp), %r8
    \\ movq 8 + 128(%rsp), %rdx
    \\ movq 16 + 128(%rsp), %rbp
    \\ movq 24 + 128(%rsp), %r10
    \\ movq 32 + 128(%rsp), %r15
    \\ salq $4, %r8
    \\ salq $4, %rdx
    \\ salq $4, %r10
    \\ movddup %xmm0, %xmm14
    \\ movddup %xmm1, %xmm15
    \\ subq $-16 * 8, %rcx
    \\ testq %rdi, %rdi
    \\ jle 20024f
    \\ testq %rsi, %rsi
    \\ jle 20024f
    \\ .balign 8
    \\ movq %r15, %r13
    \\ movq %rbp, %r14
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20003f
    \\ .balign 16
    \\ 20002:
    \\ movlps 0 * 8(%r9), %xmm0
    \\ movhpd 1 * 8(%r9), %xmm0
    \\ addq %rdx, %r9
    \\ movlps 0 * 8(%r9), %xmm1
    \\ movhpd 1 * 8(%r9), %xmm1
    \\ addq %rdx, %r9
    \\ movlps 0 * 8(%r9), %xmm2
    \\ movhpd 1 * 8(%r9), %xmm2
    \\ addq %rdx, %r9
    \\ movlps 0 * 8(%r9), %xmm3
    \\ movhpd 1 * 8(%r9), %xmm3
    \\ addq %rdx, %r9
    \\ movaps %xmm0, 0 * 8(%r13)
    \\ movaps %xmm1, 2 * 8(%r13)
    \\ movaps %xmm2, 4 * 8(%r13)
    \\ movaps %xmm3, 6 * 8(%r13)
    \\ addq $8 * 8, %r13
    \\ decq %rax
    \\ jg 20002b
    \\ .balign 16
    \\ 20003:
    \\ movq %rdi, %rax
    \\ andq $3, %rax
    \\ jle 20005f
    \\ .balign 4
    \\ 20004:
    \\ movlps 0 * 8(%r9), %xmm0
    \\ movhpd 1 * 8(%r9), %xmm0
    \\ addq %rdx, %r9
    \\ movaps %xmm0, 0 * 8(%r13)
    \\ addq $2 * 8, %r13
    \\ decq %rax
    \\ jg 20004b
    \\ .balign 16
    \\ 20005:
    \\ cmpq $4, %rsi
    \\ jl 20012f
    \\ .balign 8
    \\ 20006:
    \\ subq $4, %rsi
    \\ leaq 16 * 8(%r15), %r13
    \\ movq %rcx, %r11
    \\ leaq (%r11, %r8, 2), %r12
    \\ leaq (%r11, %r8, 4), %rcx
    \\ movlps -16 * 8(%r13), %xmm12
    \\ movhps -16 * 8 + 8(%r13), %xmm12
    \\ xorpd %xmm0, %xmm0
    \\ xorpd %xmm1, %xmm1
    \\ xorpd %xmm2, %xmm2
    \\ xorpd %xmm3, %xmm3
    \\ movlps -14 * 8(%r13), %xmm13
    \\ movhps -14 * 8 + 8(%r13), %xmm13
    \\ xorpd %xmm4, %xmm4
    \\ xorpd %xmm5, %xmm5
    \\ xorpd %xmm6, %xmm6
    \\ xorpd %xmm7, %xmm7
    \\ prefetchw 3 * 8(%r14)
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20009f
    \\ movlps -16 * 8(%r11), %xmm8
    \\ movhps -16 * 8 + 8(%r11), %xmm8
    \\ movlps -16 * 8(%r11, %r8, 1), %xmm10
    \\ movhps -16 * 8 + 8(%r11, %r8, 1), %xmm10
    \\ decq %rax
    \\ jle 20008f
    \\ .balign 8
    \\ 20007:
    \\ prefetch (64 * 1) - 128 + (64 >> 1)(%r11)
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movlps -16 * 8(%r12), %xmm8
    \\ movhps -16 * 8 + 8(%r12), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movlps -16 * 8(%r12, %r8, 1), %xmm10
    \\ movhps -16 * 8 + 8(%r12, %r8, 1), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movlps -14 * 8(%r11), %xmm8
    \\ movhps -14 * 8 + 8(%r11), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movlps -14 * 8(%r11, %r8, 1), %xmm10
    \\ movhps -14 * 8 + 8(%r11, %r8, 1), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ movlps -12 * 8(%r13), %xmm12
    \\ movhps -12 * 8 + 8(%r13), %xmm12
    \\ addpd %xmm11, %xmm7
    \\ prefetch (64 * 1) - 128 + (64 >> 1)(%r11, %r8)
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movlps -14 * 8(%r12), %xmm8
    \\ movhps -14 * 8 + 8(%r12), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movlps -14 * 8(%r12, %r8, 1), %xmm10
    \\ movhps -14 * 8 + 8(%r12, %r8, 1), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movlps -12 * 8(%r11), %xmm8
    \\ movhps -12 * 8 + 8(%r11), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movlps -12 * 8(%r11, %r8, 1), %xmm10
    \\ movhps -12 * 8 + 8(%r11, %r8, 1), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ movlps -10 * 8(%r13), %xmm13
    \\ movhps -10 * 8 + 8(%r13), %xmm13
    \\ addpd %xmm11, %xmm7
    \\ prefetch (64 * 1) - 128 + (64 >> 1)(%r12)
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movlps -12 * 8(%r12), %xmm8
    \\ movhps -12 * 8 + 8(%r12), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movlps -12 * 8(%r12, %r8, 1), %xmm10
    \\ movhps -12 * 8 + 8(%r12, %r8, 1), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movlps -10 * 8(%r11), %xmm8
    \\ movhps -10 * 8 + 8(%r11), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movlps -10 * 8(%r11, %r8, 1), %xmm10
    \\ movhps -10 * 8 + 8(%r11, %r8, 1), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ movlps -8 * 8(%r13), %xmm12
    \\ movhps -8 * 8 + 8(%r13), %xmm12
    \\ addpd %xmm11, %xmm7
    \\ prefetch (64 * 1) - 128 + (64 >> 1)(%r12, %r8)
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movlps -10 * 8(%r12), %xmm8
    \\ movhps -10 * 8 + 8(%r12), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movlps -10 * 8(%r12, %r8, 1), %xmm10
    \\ movhps -10 * 8 + 8(%r12, %r8, 1), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ prefetch (64 * 1) - 128 + (64 >> 1)(%r13)
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movlps -8 * 8(%r11), %xmm8
    \\ movhps -8 * 8 + 8(%r11), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movlps -8 * 8(%r11, %r8, 1), %xmm10
    \\ movhps -8 * 8 + 8(%r11, %r8, 1), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ movlps -6 * 8(%r13), %xmm13
    \\ movhps -6 * 8 + 8(%r13), %xmm13
    \\ addpd %xmm11, %xmm7
    \\ subq $-8 * 8, %r11
    \\ subq $-8 * 8, %r12
    \\ subq $-8 * 8, %r13
    \\ subq $1, %rax
    \\ .byte 0x3e
    \\ jg 20007b
    \\ .balign 8
    \\ 20008:
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movlps -16 * 8(%r12), %xmm8
    \\ movhps -16 * 8 + 8(%r12), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movlps -16 * 8(%r12, %r8, 1), %xmm10
    \\ movhps -16 * 8 + 8(%r12, %r8, 1), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movlps -14 * 8(%r11), %xmm8
    \\ movhps -14 * 8 + 8(%r11), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movlps -14 * 8(%r11, %r8, 1), %xmm10
    \\ movhps -14 * 8 + 8(%r11, %r8, 1), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ movlps -12 * 8(%r13), %xmm12
    \\ movhps -12 * 8 + 8(%r13), %xmm12
    \\ addpd %xmm11, %xmm7
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movlps -14 * 8(%r12), %xmm8
    \\ movhps -14 * 8 + 8(%r12), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movlps -14 * 8(%r12, %r8, 1), %xmm10
    \\ movhps -14 * 8 + 8(%r12, %r8, 1), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movlps -12 * 8(%r11), %xmm8
    \\ movhps -12 * 8 + 8(%r11), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movlps -12 * 8(%r11, %r8, 1), %xmm10
    \\ movhps -12 * 8 + 8(%r11, %r8, 1), %xmm10
    \\ mulpd %xmm13, %xmm11
    \\ movlps -10 * 8(%r13), %xmm13
    \\ movhps -10 * 8 + 8(%r13), %xmm13
    \\ addpd %xmm11, %xmm7
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movlps -12 * 8(%r12), %xmm8
    \\ movhps -12 * 8 + 8(%r12), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movlps -12 * 8(%r12, %r8, 1), %xmm10
    \\ movhps -12 * 8 + 8(%r12, %r8, 1), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movlps -10 * 8(%r11), %xmm8
    \\ movhps -10 * 8 + 8(%r11), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movlps -10 * 8(%r11, %r8, 1), %xmm10
    \\ movhps -10 * 8 + 8(%r11, %r8, 1), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ movlps -8 * 8(%r13), %xmm12
    \\ movhps -8 * 8 + 8(%r13), %xmm12
    \\ addpd %xmm11, %xmm7
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movlps -10 * 8(%r12), %xmm8
    \\ movhps -10 * 8 + 8(%r12), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movlps -10 * 8(%r12, %r8, 1), %xmm10
    \\ movhps -10 * 8 + 8(%r12, %r8, 1), %xmm10
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
    \\ movlps -6 * 8(%r13), %xmm13
    \\ movhps -6 * 8 + 8(%r13), %xmm13
    \\ addpd %xmm11, %xmm7
    \\ subq $-8 * 8, %r11
    \\ subq $-8 * 8, %r12
    \\ subq $-8 * 8, %r13
    \\ .balign 8
    \\ 20009:
    \\ testq $2, %rdi
    \\ je 20010f
    \\ movlps -16 * 8(%r11), %xmm8
    \\ movhps -16 * 8 + 8(%r11), %xmm8
    \\ movlps -16 * 8(%r11, %r8, 1), %xmm10
    \\ movhps -16 * 8 + 8(%r11, %r8, 1), %xmm10
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movlps -16 * 8(%r12), %xmm8
    \\ movhps -16 * 8 + 8(%r12), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movlps -16 * 8(%r12, %r8, 1), %xmm10
    \\ movhps -16 * 8 + 8(%r12, %r8, 1), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ addpd %xmm11, %xmm3
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ movlps -14 * 8(%r11), %xmm8
    \\ movhps -14 * 8 + 8(%r11), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm5
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm6
    \\ movlps -14 * 8(%r11, %r8, 1), %xmm10
    \\ movhps -14 * 8 + 8(%r11, %r8, 1), %xmm10
    \\ mulpd %xmm12, %xmm11
    \\ movlps -12 * 8(%r13), %xmm12
    \\ movhps -12 * 8 + 8(%r13), %xmm12
    \\ addpd %xmm11, %xmm7
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm13, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movlps -14 * 8(%r12), %xmm8
    \\ movhps -14 * 8 + 8(%r12), %xmm8
    \\ mulpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm13, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movlps -14 * 8(%r12, %r8, 1), %xmm10
    \\ movhps -14 * 8 + 8(%r12, %r8, 1), %xmm10
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
    \\ 20010:
    \\ testq $1, %rdi
    \\ je 20011f
    \\ movlps -16 * 8(%r11), %xmm8
    \\ movhps -16 * 8 + 8(%r11), %xmm8
    \\ movlps -16 * 8(%r11, %r8, 1), %xmm10
    \\ movhps -16 * 8 + 8(%r11, %r8, 1), %xmm10
    \\ pshufd $0x4e, %xmm8, %xmm9
    \\ mulpd %xmm12, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ movlps -16 * 8(%r12), %xmm8
    \\ movhps -16 * 8 + 8(%r12), %xmm8
    \\ mulpd %xmm12, %xmm9
    \\ addpd %xmm9, %xmm1
    \\ pshufd $0x4e, %xmm10, %xmm11
    \\ mulpd %xmm12, %xmm10
    \\ addpd %xmm10, %xmm2
    \\ movlps -16 * 8(%r12, %r8, 1), %xmm10
    \\ movhps -16 * 8 + 8(%r12, %r8, 1), %xmm10
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
    \\
;
