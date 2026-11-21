//! Section 2 of 2 of the text of `dsymv_L`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/symv_L_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20010:
    \\ movq %rcx, %rax
    \\ subq %r15, %rax
    \\ subq $4, %rax
    \\ test $4, %rax
    \\ jle 20011f
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm4, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm2
    \\ movsd 2 * 8(%rbx, %r8, 1), %xmm12
    \\ movhpd 3 * 8(%rbx, %r8, 1), %xmm12
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
    \\ movapd %xmm1, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm12, %xmm3
    \\ movsd 0 * 8(%rbp, %r8, 1), %xmm12
    \\ movhpd 1 * 8(%rbp, %r8, 1), %xmm12
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm6, %xmm13
    \\ addpd %xmm15, %xmm10
    \\ addpd %xmm13, %xmm2
    \\ movsd 2 * 8(%rbp, %r8, 1), %xmm13
    \\ movhpd 3 * 8(%rbp, %r8, 1), %xmm13
    \\ movapd %xmm1, %xmm15
    \\ mulpd %xmm14, %xmm15
    \\ mulpd %xmm6, %xmm14
    \\ addpd %xmm15, %xmm10
    \\ addpd %xmm14, %xmm3
    \\ movsd 4 * 8(%rbx, %r8, 1), %xmm14
    \\ movhpd 5 * 8(%rbx, %r8, 1), %xmm14
    \\ movapd %xmm0, %xmm15
    \\ movapd 4 * 8(%r13), %xmm0
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm7, %xmm12
    \\ addpd %xmm15, %xmm11
    \\ addpd %xmm12, %xmm2
    \\ movsd 4 * 8(%rbx), %xmm12
    \\ movhpd 5 * 8(%rbx), %xmm12
    \\ movapd %xmm1, %xmm15
    \\ movapd 6 * 8(%r13), %xmm1
    \\ mulpd %xmm13, %xmm15
    \\ mulpd %xmm7, %xmm13
    \\ addpd %xmm15, %xmm11
    \\ addpd %xmm13, %xmm3
    \\ movsd 6 * 8(%rbx), %xmm13
    \\ movhpd 7 * 8(%rbx), %xmm13
    \\ movsd %xmm2, 0 * 8(%r14)
    \\ movhpd %xmm2, 1 * 8(%r14)
    \\ movsd 4 * 8(%r14), %xmm2
    \\ movhpd 5 * 8(%r14), %xmm2
    \\ movsd %xmm3, 2 * 8(%r14)
    \\ movhpd %xmm3, 3 * 8(%r14)
    \\ movsd 6 * 8(%r14), %xmm3
    \\ movhpd 7 * 8(%r14), %xmm3
    \\ addq $4 * 8, %r13
    \\ addq $4 * 8, %r14
    \\ addq $4 * 8, %rbx
    \\ addq $4 * 8, %rbp
    \\ .balign 8
    \\ 20011:
    \\ testq $2, %rcx
    \\ jle 20012f
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm4, %xmm12
    \\ addpd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm2
    \\ movsd 0 * 8(%rbx, %r8, 1), %xmm12
    \\ movhpd 1 * 8(%rbx, %r8, 1), %xmm12
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm5, %xmm12
    \\ addpd %xmm15, %xmm9
    \\ addpd %xmm12, %xmm2
    \\ movsd 0 * 8(%rbp), %xmm12
    \\ movhpd 1 * 8(%rbp), %xmm12
    \\ movapd %xmm0, %xmm15
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm6, %xmm12
    \\ addpd %xmm15, %xmm10
    \\ addpd %xmm12, %xmm2
    \\ movsd 0 * 8(%rbp, %r8, 1), %xmm12
    \\ movhpd 1 * 8(%rbp, %r8, 1), %xmm12
    \\ movapd %xmm0, %xmm15
    \\ movapd 2 * 8(%r13), %xmm0
    \\ mulpd %xmm12, %xmm15
    \\ mulpd %xmm7, %xmm12
    \\ addpd %xmm15, %xmm11
    \\ addpd %xmm12, %xmm2
    \\ movsd 2 * 8(%rbx), %xmm12
    \\ movsd %xmm2, 0 * 8(%r14)
    \\ movhpd %xmm2, 1 * 8(%r14)
    \\ movsd 2 * 8(%r14), %xmm2
    \\ addq $2 * 8, %r13
    \\ addq $2 * 8, %r14
    \\ addq $2 * 8, %rbx
    \\ addq $2 * 8, %rbp
    \\ .balign 8
    \\ 20012:
    \\ testq $1, %rcx
    \\ jle 20013f
    \\ movapd %xmm0, %xmm15
    \\ mulsd %xmm12, %xmm15
    \\ mulsd %xmm4, %xmm12
    \\ addsd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm2
    \\ movsd 0 * 8(%rbx, %r8, 1), %xmm12
    \\ movapd %xmm0, %xmm15
    \\ mulsd %xmm12, %xmm15
    \\ mulsd %xmm5, %xmm12
    \\ addsd %xmm15, %xmm9
    \\ addsd %xmm12, %xmm2
    \\ movsd 0 * 8(%rbp), %xmm12
    \\ movapd %xmm0, %xmm15
    \\ mulsd %xmm12, %xmm15
    \\ mulsd %xmm6, %xmm12
    \\ addsd %xmm15, %xmm10
    \\ addsd %xmm12, %xmm2
    \\ movsd 0 * 8(%rbp, %r8, 1), %xmm12
    \\ movapd %xmm0, %xmm15
    \\ mulsd %xmm12, %xmm15
    \\ mulsd %xmm7, %xmm12
    \\ addsd %xmm15, %xmm11
    \\ addsd %xmm12, %xmm2
    \\ movsd %xmm2, 0 * 8(%r14)
    \\ .balign 8
    \\ 20013:
    \\ haddpd %xmm9, %xmm8
    \\ haddpd %xmm11, %xmm10
    \\ movsd 0 * 8(%rdi, %r15, 8), %xmm2
    \\ movhpd 1 * 8(%rdi, %r15, 8), %xmm2
    \\ movsd 2 * 8(%rdi, %r15, 8), %xmm3
    \\ movhpd 3 * 8(%rdi, %r15, 8), %xmm3
    \\ addpd %xmm8, %xmm2
    \\ addpd %xmm10, %xmm3
    \\ movsd %xmm2, 0 * 8(%rdi, %r15, 8)
    \\ movhpd %xmm2, 1 * 8(%rdi, %r15, 8)
    \\ movsd %xmm3, 2 * 8(%rdi, %r15, 8)
    \\ movhpd %xmm3, 3 * 8(%rdi, %r15, 8)
    \\ addq $4, %r15
    \\ movq %r15, %rax
    \\ addq $4, %rax
    \\ cmpq %rdx, %rax
    \\ jle 20008b
    \\ .balign 8
    \\ 20014:
    \\ testq $2, %rdx
    \\ jle 20016f
    \\ movq %r9, %rbx
    \\ leaq 2 * 8(%r9, %r8, 2), %r9
    \\ movapd 0 * 8(%r12, %r15, 8), %xmm5
    \\ movsd 0 * 8(%rbx), %xmm8
    \\ movhpd 1 * 8(%rbx), %xmm8
    \\ mulpd %xmm5, %xmm8
    \\ movsd 1 * 8(%rbx), %xmm9
    \\ movhpd 1 * 8(%rbx, %r8, 1), %xmm9
    \\ mulpd %xmm5, %xmm9
    \\ movddup %xmm5, %xmm4
    \\ unpckhpd %xmm5, %xmm5
    \\ testq $1, %rcx
    \\ jle 20015f
    \\ movsd 2 * 8(%rbx), %xmm12
    \\ movsd 2 * 8(%rbx, %r8, 1), %xmm13
    \\ movsd 2 * 8(%r12, %r15, 8), %xmm0
    \\ movsd 2 * 8(%rdi, %r15, 8), %xmm2
    \\ movapd %xmm0, %xmm15
    \\ mulsd %xmm12, %xmm15
    \\ mulsd %xmm4, %xmm12
    \\ addsd %xmm15, %xmm8
    \\ addpd %xmm12, %xmm2
    \\ movapd %xmm0, %xmm15
    \\ mulsd %xmm13, %xmm15
    \\ mulsd %xmm5, %xmm13
    \\ addsd %xmm15, %xmm9
    \\ addsd %xmm13, %xmm2
    \\ movsd %xmm2, 2 * 8(%rdi, %r15, 8)
    \\ .balign 8
    \\ 20015:
    \\ haddpd %xmm9, %xmm8
    \\ movsd 0 * 8(%rdi, %r15, 8), %xmm2
    \\ movhpd 1 * 8(%rdi, %r15, 8), %xmm2
    \\ addpd %xmm8, %xmm2
    \\ movsd %xmm2, 0 * 8(%rdi, %r15, 8)
    \\ movhpd %xmm2, 1 * 8(%rdi, %r15, 8)
    \\ addq $2, %r15
    \\ .balign 8
    \\ 20016:
    \\ testq $1, %rdx
    \\ jle 20017f
    \\ movsd 0 * 8(%r9), %xmm8
    \\ movsd 0 * 8(%r12, %r15, 8), %xmm4
    \\ movsd 0 * 8(%rdi, %r15, 8), %xmm2
    \\ mulsd %xmm4, %xmm8
    \\ addsd %xmm8, %xmm2
    \\ movsd %xmm2, 0 * 8(%rdi, %r15, 8)
    \\ .balign 8
    \\ 20017:
    \\ cmpq $8, %r11
    \\ je 20021f
    \\ movq %rcx, %rax
    \\ sarq $3, %rax
    \\ jle 20019f
    \\ .balign 8
    \\ 20018:
    \\ movapd 0 * 8(%rdi), %xmm0
    \\ movapd 2 * 8(%rdi), %xmm1
    \\ movapd 4 * 8(%rdi), %xmm2
    \\ movapd 6 * 8(%rdi), %xmm3
    \\ movsd %xmm0, 0 * 8(%r10)
    \\ addq %r11, %r10
    \\ movhpd %xmm0, 0 * 8(%r10)
    \\ addq %r11, %r10
    \\ movsd %xmm1, 0 * 8(%r10)
    \\ addq %r11, %r10
    \\ movhpd %xmm1, 0 * 8(%r10)
    \\ addq %r11, %r10
    \\ movsd %xmm2, 0 * 8(%r10)
    \\ addq %r11, %r10
    \\ movhpd %xmm2, 0 * 8(%r10)
    \\ addq %r11, %r10
    \\ movsd %xmm3, 0 * 8(%r10)
    \\ addq %r11, %r10
    \\ movhpd %xmm3, 0 * 8(%r10)
    \\ addq %r11, %r10
    \\ addq $8 * 8, %rdi
    \\ decq %rax
    \\ jg 20018b
    \\ .balign 8
    \\ 20019:
    \\ movq %rcx, %rax
    \\ andq $7, %rax
    \\ jle 20021f
    \\ .balign 8
    \\ 20020:
    \\ movsd 0 * 8(%rdi), %xmm0
    \\ movsd %xmm0, 0 * 8(%r10)
    \\ addq %r11, %r10
    \\ addq $1 * 8, %rdi
    \\ decq %rax
    \\ jg 20020b
    \\ .balign 8
    \\ 20021:
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
