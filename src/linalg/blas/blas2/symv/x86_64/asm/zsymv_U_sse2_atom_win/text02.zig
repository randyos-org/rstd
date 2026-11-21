//! Section 2 of 2 of the text of `zsymv_U`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zsymv_U_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20010:
    \\ testq $2, %rdx
    \\ jle 20011f
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm4, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm10
    \\ movddup 1 * 8(%rbx), %xmm12
    \\ movapd %xmm2, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm6, %xmm13
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm13, %xmm11
    \\ movddup 3 * 8(%rbp), %xmm13
    \\ movapd %xmm1, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm10
    \\ movddup 2 * 8(%rbx), %xmm12
    \\ movapd %xmm3, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm7, %xmm13
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm13, %xmm11
    \\ movddup 0 * 8(%rbp), %xmm13
    \\ movapd %xmm2, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm4, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm11
    \\ movddup 3 * 8(%rbx), %xmm12
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm6, %xmm13
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm13, %xmm10
    \\ movddup 1 * 8(%rbp), %xmm13
    \\ movapd %xmm3, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm11
    \\ movlpd %xmm11, 2 * 8(%r14)
    \\ movhpd %xmm11, 3 * 8(%r14)
    \\ movsd 6 * 8(%r14), %xmm11
    \\ movhpd 7 * 8(%r14), %xmm11
    \\ movapd %xmm1, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm7, %xmm13
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm13, %xmm10
    \\ movlpd %xmm10, 0 * 8(%r14)
    \\ movhpd %xmm10, 1 * 8(%r14)
    \\ movsd 4 * 8(%r14), %xmm10
    \\ movhpd 5 * 8(%r14), %xmm10
    \\ addq $4 * 8, %r14
    \\ addq $4 * 8, %rbx
    \\ addq $4 * 8, %rbp
    \\ .balign 8
    \\ 20011:
    \\ movddup 0 * 8(%rbx), %xmm12
    \\ movddup 0 * 8(%rbp), %xmm13
    \\ mulpd %xmm4, %xmm12
    \\ mulpd %xmm4, %xmm13
    \\ addpd %xmm12, %xmm8
    \\ addpd %xmm13, %xmm9
    \\ movddup 1 * 8(%rbx), %xmm12
    \\ movddup 1 * 8(%rbp), %xmm13
    \\ mulpd %xmm5, %xmm12
    \\ mulpd %xmm5, %xmm13
    \\ addpd %xmm12, %xmm8
    \\ addpd %xmm13, %xmm9
    \\ movddup 0 * 8(%rbp), %xmm12
    \\ movddup 2 * 8(%rbp), %xmm13
    \\ mulpd %xmm6, %xmm12
    \\ mulpd %xmm6, %xmm13
    \\ addpd %xmm12, %xmm8
    \\ addpd %xmm13, %xmm9
    \\ movddup 1 * 8(%rbp), %xmm12
    \\ movddup 3 * 8(%rbp), %xmm13
    \\ mulpd %xmm7, %xmm12
    \\ mulpd %xmm7, %xmm13
    \\ addpd %xmm12, %xmm8
    \\ addpd %xmm13, %xmm9
    \\ addpd %xmm8, %xmm10
    \\ addpd %xmm9, %xmm11
    \\ movlpd %xmm10, 0 * 8(%r14)
    \\ movhpd %xmm10, 1 * 8(%r14)
    \\ movlpd %xmm11, 2 * 8(%r14)
    \\ movhpd %xmm11, 3 * 8(%r14)
    \\ addq $2, %rdx
    \\ movq %rdx, %rax
    \\ addq $2, %rax
    \\ cmpq %rcx, %rax
    \\ jle 20008b
    \\ .balign 8
    \\ 20012:
    \\ testq $1, %rcx
    \\ jle 20015f
    \\ movq %r9, %rbx
    \\ leaq (, %rdx, 4), %rax
    \\ movapd 0 * 8(%r12, %rax, 8), %xmm4
    \\ movapd 2 * 8(%r12, %rax, 8), %xmm5
    \\ pxor %xmm8, %xmm8
    \\ pxor %xmm9, %xmm9
    \\ movddup 0 * 8(%rbx), %xmm12
    \\ movddup 1 * 8(%rbx), %xmm13
    \\ movapd 0 * 8(%r12), %xmm0
    \\ movapd 2 * 8(%r12), %xmm1
    \\ movapd 4 * 8(%r12), %xmm2
    \\ movapd 6 * 8(%r12), %xmm3
    \\ movsd 0 * 8(%rdi), %xmm10
    \\ movhpd 1 * 8(%rdi), %xmm10
    \\ movsd 2 * 8(%rdi), %xmm11
    \\ movhpd 3 * 8(%rdi), %xmm11
    \\ movq %r12, %r13
    \\ movq %rdi, %r14
    \\ movq %rdx, %rax
    \\ sarq $1, %rax
    \\ jle 20014f
    \\ .balign 8
    \\ 20013:
    \\ movapd %xmm0, %xmm15
    \\ movapd 8 * 8(%r13), %xmm0
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm4, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm10
    \\ movddup 2 * 8(%rbx), %xmm12
    \\ movapd %xmm1, %xmm15
    \\ movapd 10 * 8(%r13), %xmm1
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm5, %xmm13
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm13, %xmm10
    \\ movddup 3 * 8(%rbx), %xmm13
    \\ movlpd %xmm10, 0 * 8(%r14)
    \\ movhpd %xmm10, 1 * 8(%r14)
    \\ movsd 4 * 8(%r14), %xmm10
    \\ movhpd 5 * 8(%r14), %xmm10
    \\ movapd %xmm2, %xmm15
    \\ movapd 12 * 8(%r13), %xmm2
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm4, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm11
    \\ movddup 4 * 8(%rbx), %xmm12
    \\ movapd %xmm3, %xmm15
    \\ movapd 14 * 8(%r13), %xmm3
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm5, %xmm13
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm13, %xmm11
    \\ movddup 5 * 8(%rbx), %xmm13
    \\ movlpd %xmm11, 2 * 8(%r14)
    \\ movhpd %xmm11, 3 * 8(%r14)
    \\ movsd 6 * 8(%r14), %xmm11
    \\ movhpd 7 * 8(%r14), %xmm11
    \\ addq $8 * 8, %r13
    \\ addq $4 * 8, %r14
    \\ addq $4 * 8, %rbx
    \\ decq %rax
    \\ jg 20013b
    \\ .balign 8
    \\ 20014:
    \\ movddup 0 * 8(%rbx), %xmm12
    \\ movddup 1 * 8(%rbx), %xmm13
    \\ mulpd %xmm4, %xmm12
    \\ mulpd %xmm5, %xmm13
    \\ addpd %xmm12, %xmm8
    \\ addpd %xmm13, %xmm9
    \\ addpd %xmm9, %xmm8
    \\ addpd %xmm8, %xmm10
    \\ movlpd %xmm10, 0 * 8(%r14)
    \\ movhpd %xmm10, 1 * 8(%r14)
    \\ .balign 8
    \\ 20015:
    \\ cmpq $2 * 8, %r11
    \\ je 20019f
    \\ movq %rcx, %rax
    \\ sarq $2, %rax
    \\ jle 20017f
    \\ .balign 8
    \\ 20016:
    \\ movapd 0 * 8(%rdi), %xmm0
    \\ movapd 2 * 8(%rdi), %xmm1
    \\ movapd 4 * 8(%rdi), %xmm2
    \\ movapd 6 * 8(%rdi), %xmm3
    \\ movsd %xmm0, 0 * 8(%r10)
    \\ movhpd %xmm0, 1 * 8(%r10)
    \\ addq %r11, %r10
    \\ movsd %xmm1, 0 * 8(%r10)
    \\ movhpd %xmm1, 1 * 8(%r10)
    \\ addq %r11, %r10
    \\ movsd %xmm2, 0 * 8(%r10)
    \\ movhpd %xmm2, 1 * 8(%r10)
    \\ addq %r11, %r10
    \\ movsd %xmm3, 0 * 8(%r10)
    \\ movhpd %xmm3, 1 * 8(%r10)
    \\ addq %r11, %r10
    \\ addq $8 * 8, %rdi
    \\ decq %rax
    \\ jg 20016b
    \\ .balign 8
    \\ 20017:
    \\ movq %rcx, %rax
    \\ andq $3, %rax
    \\ jle 20019f
    \\ .balign 8
    \\ 20018:
    \\ movapd 0 * 8(%rdi), %xmm0
    \\ movsd %xmm0, 0 * 8(%r10)
    \\ movhpd %xmm0, 1 * 8(%r10)
    \\ addq %r11, %r10
    \\ addq $2 * 8, %rdi
    \\ decq %rax
    \\ jg 20018b
    \\ .balign 8
    \\ 20019:
    \\ movq 0(%rsp), %rbx
    \\ movq 8(%rsp), %rbp
    \\ movq 16(%rsp), %r12
    \\ movq 24(%rsp), %r13
    \\ movq 32(%rsp), %r14
    \\ movq 40(%rsp), %r15
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
    \\ ret
    \\ .endm
    \\
;
