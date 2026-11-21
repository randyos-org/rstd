//! Section 2 of 2 of the text of `zgemv_s`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zgemv_n_atom.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20006:
    \\ .if \bp0
    \\ testq $1, %rdi
    \\ .endif
    \\ .if \bp1
    \\ testq $1, %rcx
    \\ .endif
    \\ je 20007f
    \\ movsd -16 * 8(%r12), %xmm8
    \\ movsd -15 * 8(%r12), %xmm9
    \\ movsd -16 * 8(%r13), %xmm10
    \\ movsd -15 * 8(%r13), %xmm11
    \\ movaps %xmm8, %xmm12
    \\ mulsd %xmm4, %xmm8
    \\ mulsd %xmm5, %xmm12
    \\ movaps %xmm9, %xmm13
    \\ mulsd %xmm5, %xmm9
    \\ addsd %xmm8, %xmm0
    \\ mulsd %xmm4, %xmm13
    \\ subsd %xmm12, %xmm1
    \\ movaps %xmm10, %xmm12
    \\ mulsd %xmm6, %xmm10
    \\ subsd %xmm9, %xmm0
    \\ mulsd %xmm7, %xmm12
    \\ subsd %xmm13, %xmm1
    \\ movaps %xmm11, %xmm13
    \\ mulsd %xmm7, %xmm11
    \\ addsd %xmm10, %xmm0
    \\ mulsd %xmm6, %xmm13
    \\ subsd %xmm12, %xmm1
    \\ subsd %xmm11, %xmm0
    \\ subsd %xmm13, %xmm1
    \\ movlpd %xmm0, -16 * 8(%r14)
    \\ movlpd %xmm1, -15 * 8(%r14)
    \\ .balign 8
    \\ 20007:
    \\ decq %r11
    \\ jg 20002b
    \\ .balign 8
    \\ 20008:
    \\ .if \bp0
    \\ testq $1, %rsi
    \\ .endif
    \\ .if \bp1
    \\ testq $1, %rdx
    \\ .endif
    \\ jle 20013f
    \\ leaq 16 * 8(%r15), %r14
    \\ .if \bp0
    \\ movq %rcx, %r12
    \\ movsd 0 * 8(%r9), %xmm4
    \\ movsd 1 * 8(%r9), %xmm5
    \\ .endif
    \\ .if \bp1
    \\ movq %r8, %r12
    \\ movsd 0 * 8(%rdi), %xmm4
    \\ movsd 1 * 8(%rdi), %xmm5
    \\ .endif
    \\ movaps %xmm4, %xmm8
    \\ mulsd %xmm14, %xmm4
    \\ mulsd %xmm15, %xmm8
    \\ movaps %xmm5, %xmm9
    \\ mulsd %xmm15, %xmm9
    \\ mulsd %xmm14, %xmm5
    \\ addsd %xmm9, %xmm4
    \\ subsd %xmm8, %xmm5
    \\ movsd -16 * 8(%r14), %xmm0
    \\ movsd -15 * 8(%r14), %xmm1
    \\ movsd -14 * 8(%r14), %xmm2
    \\ movsd -13 * 8(%r14), %xmm3
    \\ .balign 8
    \\ .if \bp0
    \\ movq %rdi, %rax
    \\ .endif
    \\ .if \bp1
    \\ movq %rcx, %rax
    \\ .endif
    \\ sarq $2, %rax
    \\ jle 20011f
    \\ movsd -16 * 8(%r12), %xmm8
    \\ movsd -15 * 8(%r12), %xmm9
    \\ movsd -14 * 8(%r12), %xmm10
    \\ movsd -13 * 8(%r12), %xmm11
    \\ movaps %xmm8, %xmm12
    \\ mulsd %xmm4, %xmm8
    \\ mulsd %xmm5, %xmm12
    \\ movaps %xmm10, %xmm13
    \\ mulsd %xmm4, %xmm10
    \\ addsd %xmm8, %xmm0
    \\ movsd -12 * 8(%r12), %xmm8
    \\ mulsd %xmm5, %xmm13
    \\ subsd %xmm12, %xmm1
    \\ decq %rax
    \\ jle 20010f
    \\ .balign 8
    \\ 20009:
    \\ prefetchnta ((8 * 6) + 0) * 8(%r13)
    \\ movaps %xmm9, %xmm12
    \\ mulsd %xmm5, %xmm9
    \\ addsd %xmm10, %xmm2
    \\ movsd -10 * 8(%r12), %xmm10
    \\ mulsd %xmm4, %xmm12
    \\ subsd %xmm13, %xmm3
    \\ movaps %xmm11, %xmm13
    \\ mulsd %xmm5, %xmm11
    \\ subsd %xmm9, %xmm0
    \\ movsd -11 * 8(%r12), %xmm9
    \\ mulsd %xmm4, %xmm13
    \\ subsd %xmm12, %xmm1
    \\ movaps %xmm8, %xmm12
    \\ movlpd %xmm0, -16 * 8(%r14)
    \\ mulsd %xmm4, %xmm8
    \\ movsd -12 * 8(%r14), %xmm0
    \\ subsd %xmm11, %xmm2
    \\ movsd -9 * 8(%r12), %xmm11
    \\ mulsd %xmm5, %xmm12
    \\ movlpd %xmm1, -15 * 8(%r14)
    \\ subsd %xmm13, %xmm3
    \\ movsd -11 * 8(%r14), %xmm1
    \\ movaps %xmm10, %xmm13
    \\ movlpd %xmm2, -14 * 8(%r14)
    \\ mulsd %xmm4, %xmm10
    \\ movsd -10 * 8(%r14), %xmm2
    \\ addsd %xmm8, %xmm0
    \\ movsd -8 * 8(%r12), %xmm8
    \\ mulsd %xmm5, %xmm13
    \\ movlpd %xmm3, -13 * 8(%r14)
    \\ subsd %xmm12, %xmm1
    \\ movsd -9 * 8(%r14), %xmm3
    \\ movaps %xmm9, %xmm12
    \\ mulsd %xmm5, %xmm9
    \\ addsd %xmm10, %xmm2
    \\ movsd -6 * 8(%r12), %xmm10
    \\ mulsd %xmm4, %xmm12
    \\ subsd %xmm13, %xmm3
    \\ movaps %xmm11, %xmm13
    \\ mulsd %xmm5, %xmm11
    \\ subsd %xmm9, %xmm0
    \\ movsd -7 * 8(%r12), %xmm9
    \\ mulsd %xmm4, %xmm13
    \\ subq $-8 * 8, %r12
    \\ subsd %xmm12, %xmm1
    \\ movaps %xmm8, %xmm12
    \\ movlpd %xmm0, -12 * 8(%r14)
    \\ mulsd %xmm4, %xmm8
    \\ movsd -8 * 8(%r14), %xmm0
    \\ subsd %xmm11, %xmm2
    \\ movsd -13 * 8(%r12), %xmm11
    \\ mulsd %xmm5, %xmm12
    \\ movlpd %xmm1, -11 * 8(%r14)
    \\ subsd %xmm13, %xmm3
    \\ movsd -7 * 8(%r14), %xmm1
    \\ movaps %xmm10, %xmm13
    \\ movlpd %xmm2, -10 * 8(%r14)
    \\ mulsd %xmm4, %xmm10
    \\ movsd -6 * 8(%r14), %xmm2
    \\ addsd %xmm8, %xmm0
    \\ movsd -12 * 8(%r12), %xmm8
    \\ mulsd %xmm5, %xmm13
    \\ movlpd %xmm3, -9 * 8(%r14)
    \\ subsd %xmm12, %xmm1
    \\ movsd -5 * 8(%r14), %xmm3
    \\ subq $-8 * 8, %r14
    \\ subq $1, %rax
    \\ .byte 0x3e
    \\ jg 20009b
    \\ .balign 8
    \\ 20010:
    \\ movaps %xmm9, %xmm12
    \\ mulsd %xmm5, %xmm9
    \\ addsd %xmm10, %xmm2
    \\ movsd -10 * 8(%r12), %xmm10
    \\ mulsd %xmm4, %xmm12
    \\ subsd %xmm13, %xmm3
    \\ movaps %xmm11, %xmm13
    \\ mulsd %xmm5, %xmm11
    \\ subsd %xmm9, %xmm0
    \\ movsd -11 * 8(%r12), %xmm9
    \\ mulsd %xmm4, %xmm13
    \\ subsd %xmm12, %xmm1
    \\ movaps %xmm8, %xmm12
    \\ movlpd %xmm0, -16 * 8(%r14)
    \\ mulsd %xmm4, %xmm8
    \\ movsd -12 * 8(%r14), %xmm0
    \\ subsd %xmm11, %xmm2
    \\ movsd -9 * 8(%r12), %xmm11
    \\ mulsd %xmm5, %xmm12
    \\ movlpd %xmm1, -15 * 8(%r14)
    \\ subsd %xmm13, %xmm3
    \\ movsd -11 * 8(%r14), %xmm1
    \\ movaps %xmm10, %xmm13
    \\ mulsd %xmm4, %xmm10
    \\ movlpd %xmm2, -14 * 8(%r14)
    \\ addsd %xmm8, %xmm0
    \\ movsd -10 * 8(%r14), %xmm2
    \\ mulsd %xmm5, %xmm13
    \\ movlpd %xmm3, -13 * 8(%r14)
    \\ subsd %xmm12, %xmm1
    \\ movsd -9 * 8(%r14), %xmm3
    \\ movaps %xmm9, %xmm12
    \\ mulsd %xmm5, %xmm9
    \\ addsd %xmm10, %xmm2
    \\ mulsd %xmm4, %xmm12
    \\ subsd %xmm13, %xmm3
    \\ movaps %xmm11, %xmm13
    \\ mulsd %xmm5, %xmm11
    \\ subsd %xmm9, %xmm0
    \\ mulsd %xmm4, %xmm13
    \\ subsd %xmm12, %xmm1
    \\ subsd %xmm11, %xmm2
    \\ movlpd %xmm0, -12 * 8(%r14)
    \\ movsd -8 * 8(%r14), %xmm0
    \\ subsd %xmm13, %xmm3
    \\ movlpd %xmm1, -11 * 8(%r14)
    \\ movsd -7 * 8(%r14), %xmm1
    \\ movlpd %xmm2, -10 * 8(%r14)
    \\ movlpd %xmm3, -9 * 8(%r14)
    \\ movsd -6 * 8(%r14), %xmm2
    \\ movsd -5 * 8(%r14), %xmm3
    \\ subq $-8 * 8, %r12
    \\ subq $-8 * 8, %r14
    \\ .balign 8
    \\ 20011:
    \\ .if \bp0
    \\ testq $2, %rdi
    \\ .endif
    \\ .if \bp1
    \\ testq $2, %rcx
    \\ .endif
    \\ je 20012f
    \\ movsd -16 * 8(%r12), %xmm8
    \\ movsd -15 * 8(%r12), %xmm9
    \\ movsd -14 * 8(%r12), %xmm10
    \\ movsd -13 * 8(%r12), %xmm11
    \\ movaps %xmm8, %xmm12
    \\ mulsd %xmm4, %xmm8
    \\ mulsd %xmm5, %xmm12
    \\ movaps %xmm10, %xmm13
    \\ mulsd %xmm4, %xmm10
    \\ addsd %xmm8, %xmm0
    \\ mulsd %xmm5, %xmm13
    \\ subsd %xmm12, %xmm1
    \\ movaps %xmm9, %xmm12
    \\ mulsd %xmm5, %xmm9
    \\ addsd %xmm10, %xmm2
    \\ mulsd %xmm4, %xmm12
    \\ subsd %xmm13, %xmm3
    \\ movaps %xmm11, %xmm13
    \\ mulsd %xmm5, %xmm11
    \\ subsd %xmm9, %xmm0
    \\ mulsd %xmm4, %xmm13
    \\ subsd %xmm12, %xmm1
    \\ subsd %xmm11, %xmm2
    \\ movlpd %xmm0, -16 * 8(%r14)
    \\ movsd -12 * 8(%r14), %xmm0
    \\ subsd %xmm13, %xmm3
    \\ movlpd %xmm1, -15 * 8(%r14)
    \\ movsd -11 * 8(%r14), %xmm1
    \\ movlpd %xmm2, -14 * 8(%r14)
    \\ movlpd %xmm3, -13 * 8(%r14)
    \\ addq $4 * 8, %r12
    \\ addq $4 * 8, %r14
    \\ .balign 8
    \\ 20012:
    \\ .if \bp0
    \\ testq $1, %rdi
    \\ .endif
    \\ .if \bp1
    \\ testq $1, %rcx
    \\ .endif
    \\ je 20013f
    \\ movsd -16 * 8(%r12), %xmm8
    \\ movsd -15 * 8(%r12), %xmm9
    \\ movaps %xmm8, %xmm12
    \\ mulsd %xmm4, %xmm8
    \\ mulsd %xmm5, %xmm12
    \\ movaps %xmm9, %xmm13
    \\ mulsd %xmm5, %xmm9
    \\ addsd %xmm8, %xmm0
    \\ mulsd %xmm4, %xmm13
    \\ subsd %xmm12, %xmm1
    \\ subsd %xmm9, %xmm0
    \\ subsd %xmm13, %xmm1
    \\ movlpd %xmm0, -16 * 8(%r14)
    \\ movlpd %xmm1, -15 * 8(%r14)
    \\ .balign 8
    \\ 20013:
    \\ movq %rbp, %r14
    \\ .if \bp0
    \\ movq %rdi, %rax
    \\ .endif
    \\ .if \bp1
    \\ movq %rcx, %rax
    \\ .endif
    \\ sarq $3, %rax
    \\ jle 20015f
    \\ .balign 8
    \\ 20014:
    \\ movsd 0 * 8(%rbp), %xmm0
    \\ movhpd 1 * 8(%rbp), %xmm0
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm1
    \\ movhpd 1 * 8(%rbp), %xmm1
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm2
    \\ movhpd 1 * 8(%rbp), %xmm2
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm3
    \\ movhpd 1 * 8(%rbp), %xmm3
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm4
    \\ movhpd 1 * 8(%rbp), %xmm4
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm5
    \\ movhpd 1 * 8(%rbp), %xmm5
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm6
    \\ movhpd 1 * 8(%rbp), %xmm6
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm7
    \\ movhpd 1 * 8(%rbp), %xmm7
    \\ addq %r10, %rbp
    \\ addpd 0 * 8(%r15), %xmm0
    \\ addpd 2 * 8(%r15), %xmm1
    \\ addpd 4 * 8(%r15), %xmm2
    \\ addpd 6 * 8(%r15), %xmm3
    \\ addpd 8 * 8(%r15), %xmm4
    \\ addpd 10 * 8(%r15), %xmm5
    \\ addpd 12 * 8(%r15), %xmm6
    \\ addpd 14 * 8(%r15), %xmm7
    \\ movlpd %xmm0, 0 * 8(%r14)
    \\ movhpd %xmm0, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm1, 0 * 8(%r14)
    \\ movhpd %xmm1, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm2, 0 * 8(%r14)
    \\ movhpd %xmm2, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm3, 0 * 8(%r14)
    \\ movhpd %xmm3, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm4, 0 * 8(%r14)
    \\ movhpd %xmm4, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm5, 0 * 8(%r14)
    \\ movhpd %xmm5, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm6, 0 * 8(%r14)
    \\ movhpd %xmm6, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm7, 0 * 8(%r14)
    \\ movhpd %xmm7, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ subq $-16 * 8, %r15
    \\ decq %rax
    \\ jg 20014b
    \\ .balign 8
    \\ 20015:
    \\ .if \bp0
    \\ testq $7, %rdi
    \\ .endif
    \\ .if \bp1
    \\ testq $7, %rcx
    \\ .endif
    \\ jle 20018f
    \\ .if \bp0
    \\ testq $4, %rdi
    \\ .endif
    \\ .if \bp1
    \\ testq $4, %rcx
    \\ .endif
    \\ jle 20016f
    \\ movsd 0 * 8(%rbp), %xmm0
    \\ movhpd 1 * 8(%rbp), %xmm0
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm1
    \\ movhpd 1 * 8(%rbp), %xmm1
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm2
    \\ movhpd 1 * 8(%rbp), %xmm2
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm3
    \\ movhpd 1 * 8(%rbp), %xmm3
    \\ addq %r10, %rbp
    \\ addpd 0 * 8(%r15), %xmm0
    \\ addpd 2 * 8(%r15), %xmm1
    \\ addpd 4 * 8(%r15), %xmm2
    \\ addpd 6 * 8(%r15), %xmm3
    \\ movlpd %xmm0, 0 * 8(%r14)
    \\ movhpd %xmm0, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm1, 0 * 8(%r14)
    \\ movhpd %xmm1, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm2, 0 * 8(%r14)
    \\ movhpd %xmm2, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm3, 0 * 8(%r14)
    \\ movhpd %xmm3, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ addq $8 * 8, %r15
    \\ .balign 8
    \\ 20016:
    \\ .if \bp0
    \\ testq $2, %rdi
    \\ .endif
    \\ .if \bp1
    \\ testq $2, %rcx
    \\ .endif
    \\ jle 20017f
    \\ movsd 0 * 8(%rbp), %xmm0
    \\ movhpd 1 * 8(%rbp), %xmm0
    \\ addq %r10, %rbp
    \\ movsd 0 * 8(%rbp), %xmm1
    \\ movhpd 1 * 8(%rbp), %xmm1
    \\ addq %r10, %rbp
    \\ addpd 0 * 8(%r15), %xmm0
    \\ addpd 2 * 8(%r15), %xmm1
    \\ movlpd %xmm0, 0 * 8(%r14)
    \\ movhpd %xmm0, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ movlpd %xmm1, 0 * 8(%r14)
    \\ movhpd %xmm1, 1 * 8(%r14)
    \\ addq %r10, %r14
    \\ addq $4 * 8, %r15
    \\ .balign 8
    \\ 20017:
    \\ .if \bp0
    \\ testq $1, %rdi
    \\ .endif
    \\ .if \bp1
    \\ testq $1, %rcx
    \\ .endif
    \\ jle 20018f
    \\ movsd 0 * 8(%rbp), %xmm0
    \\ movhpd 1 * 8(%rbp), %xmm0
    \\ addpd 0 * 8(%r15), %xmm0
    \\ movlpd %xmm0, 0 * 8(%r14)
    \\ movhpd %xmm0, 1 * 8(%r14)
    \\ .balign 8
    \\ 20018:
    \\ movq 0(%rsp), %rbx
    \\ movq 8(%rsp), %rbp
    \\ movq 16(%rsp), %r12
    \\ movq 24(%rsp), %r13
    \\ movq 32(%rsp), %r14
    \\ movq 40(%rsp), %r15
    \\ .if \bp0
    \\ addq $64, %rsp
    \\ .endif
    \\ .if \bp1
    \\ movq 48(%rsp), %rdi
    \\ movq 56(%rsp), %rsi
    \\ movups 64(%rsp), %xmm6
    \\ movups 80(%rsp), %xmm7
    \\ movups 96(%rsp), %xmm8
    \\ movups 112(%rsp), %xmm9
    \\ movups 128(%rsp), %xmm10
    \\ movups 144(%rsp), %xmm11
    \\ movups 160(%rsp), %xmm12
    \\ movups 176(%rsp), %xmm13
    \\ movups 192(%rsp), %xmm14
    \\ movups 208(%rsp), %xmm15
    \\ addq $256, %rsp
    \\ .endif
    \\ ret
    \\ .endm
    \\
;
