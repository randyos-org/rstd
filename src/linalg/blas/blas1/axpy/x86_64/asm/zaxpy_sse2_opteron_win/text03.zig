//! Section 3 of 3 of the text of `zaxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zaxpy_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20022:
    \\ movq %rcx, %rax
    \\ andq $2, %rax
    \\ jle 20023f
    \\ movlps -16 * 8(%rdx), %xmm1
    \\ movhps -15 * 8(%rdx), %xmm1
    \\ movlps -14 * 8(%rdx), %xmm2
    \\ movhps -13 * 8(%rdx), %xmm2
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ shufpd $1, %xmm2, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movaps %xmm2, %xmm0
    \\ addq $4 * 8, %rdx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20023:
    \\ movq %rcx, %rax
    \\ andq $1, %rax
    \\ jle 20024f
    \\ movlps -16 * 8(%rdx), %xmm1
    \\ movhps -15 * 8(%rdx), %xmm1
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps %xmm1, %xmm0
    \\ addq $2 * 8, %rdx
    \\ addq $2 * 8, %r9
    \\ .balign 8
    \\ 20024:
    \\ shufpd $1, %xmm0, %xmm0
    \\ addsd -16 * 8(%r9), %xmm0
    \\ movlps %xmm0, -16 * 8(%r9)
    \\ jmp 20033f
    \\ .balign 8
    \\ 20025:
    \\ movaps %xmm0, %xmm14 # a 0
    \\ pxor %xmm15, %xmm15 # 0 0
    \\ subsd %xmm1, %xmm15 # -b 0
    \\ unpcklpd %xmm14, %xmm15 # -b a
    \\ unpcklpd %xmm1, %xmm14 # a b
    \\ movq %r9, %r11
    \\ movq %rcx, %rax
    \\ cmpq $0, %r8
    \\ je 20032f
    \\ cmpq $0, %r10
    \\ je 20032f
    \\ 20026:
    \\ sarq $3, %rax
    \\ jle 20029f
    \\ movlpd 0 * 8(%rdx), %xmm0
    \\ movhpd 0 * 8(%rdx), %xmm0
    \\ movlpd 1 * 8(%rdx), %xmm1
    \\ movhpd 1 * 8(%rdx), %xmm1
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm2
    \\ movhpd 0 * 8(%rdx), %xmm2
    \\ movlpd 1 * 8(%rdx), %xmm3
    \\ movhpd 1 * 8(%rdx), %xmm3
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm4
    \\ movhpd 0 * 8(%rdx), %xmm4
    \\ movlpd 1 * 8(%rdx), %xmm5
    \\ movhpd 1 * 8(%rdx), %xmm5
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm6
    \\ movhpd 0 * 8(%rdx), %xmm6
    \\ movlpd 1 * 8(%rdx), %xmm7
    \\ movhpd 1 * 8(%rdx), %xmm7
    \\ addq %r8, %rdx
    \\ movlps 0 * 8(%r9), %xmm8
    \\ movhpd 1 * 8(%r9), %xmm8
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm9
    \\ movhpd 1 * 8(%r9), %xmm9
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm10
    \\ movhpd 1 * 8(%r9), %xmm10
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm11
    \\ movhpd 1 * 8(%r9), %xmm11
    \\ addq %r10, %r9
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm14, %xmm6
    \\ decq %rax
    \\ jle 20028f
    \\ .balign 8
    \\ 20027:
    \\ addpd %xmm0, %xmm8
    \\ mulpd %xmm15, %xmm1
    \\ addpd %xmm2, %xmm9
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm4, %xmm10
    \\ mulpd %xmm15, %xmm5
    \\ addpd %xmm6, %xmm11
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm8
    \\ addpd %xmm3, %xmm9
    \\ addpd %xmm5, %xmm10
    \\ addpd %xmm7, %xmm11
    \\ movlpd 0 * 8(%rdx), %xmm0
    \\ movhpd 0 * 8(%rdx), %xmm0
    \\ movlpd 1 * 8(%rdx), %xmm1
    \\ movhpd 1 * 8(%rdx), %xmm1
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm2
    \\ movhpd 0 * 8(%rdx), %xmm2
    \\ movlpd 1 * 8(%rdx), %xmm3
    \\ movhpd 1 * 8(%rdx), %xmm3
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm4
    \\ movhpd 0 * 8(%rdx), %xmm4
    \\ movlpd 1 * 8(%rdx), %xmm5
    \\ movhpd 1 * 8(%rdx), %xmm5
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm6
    \\ movhpd 0 * 8(%rdx), %xmm6
    \\ movlpd 1 * 8(%rdx), %xmm7
    \\ movhpd 1 * 8(%rdx), %xmm7
    \\ addq %r8, %rdx
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm14, %xmm6
    \\ movlpd %xmm8, 0 * 8(%r11)
    \\ movhpd %xmm8, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm9, 0 * 8(%r11)
    \\ movhpd %xmm9, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm10, 0 * 8(%r11)
    \\ movhpd %xmm10, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm11, 0 * 8(%r11)
    \\ movhpd %xmm11, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlps 0 * 8(%r9), %xmm8
    \\ movhpd 1 * 8(%r9), %xmm8
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm9
    \\ movhpd 1 * 8(%r9), %xmm9
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm10
    \\ movhpd 1 * 8(%r9), %xmm10
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm11
    \\ movhpd 1 * 8(%r9), %xmm11
    \\ addq %r10, %r9
    \\ addpd %xmm0, %xmm8
    \\ mulpd %xmm15, %xmm1
    \\ addpd %xmm2, %xmm9
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm4, %xmm10
    \\ mulpd %xmm15, %xmm5
    \\ addpd %xmm6, %xmm11
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm8
    \\ addpd %xmm3, %xmm9
    \\ addpd %xmm5, %xmm10
    \\ addpd %xmm7, %xmm11
    \\ movlpd 0 * 8(%rdx), %xmm0
    \\ movhpd 0 * 8(%rdx), %xmm0
    \\ movlpd 1 * 8(%rdx), %xmm1
    \\ movhpd 1 * 8(%rdx), %xmm1
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm2
    \\ movhpd 0 * 8(%rdx), %xmm2
    \\ movlpd 1 * 8(%rdx), %xmm3
    \\ movhpd 1 * 8(%rdx), %xmm3
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm4
    \\ movhpd 0 * 8(%rdx), %xmm4
    \\ movlpd 1 * 8(%rdx), %xmm5
    \\ movhpd 1 * 8(%rdx), %xmm5
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm6
    \\ movhpd 0 * 8(%rdx), %xmm6
    \\ movlpd 1 * 8(%rdx), %xmm7
    \\ movhpd 1 * 8(%rdx), %xmm7
    \\ addq %r8, %rdx
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm14, %xmm6
    \\ movlpd %xmm8, 0 * 8(%r11)
    \\ movhpd %xmm8, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm9, 0 * 8(%r11)
    \\ movhpd %xmm9, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm10, 0 * 8(%r11)
    \\ movhpd %xmm10, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm11, 0 * 8(%r11)
    \\ movhpd %xmm11, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlps 0 * 8(%r9), %xmm8
    \\ movhpd 1 * 8(%r9), %xmm8
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm9
    \\ movhpd 1 * 8(%r9), %xmm9
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm10
    \\ movhpd 1 * 8(%r9), %xmm10
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm11
    \\ movhpd 1 * 8(%r9), %xmm11
    \\ addq %r10, %r9
    \\ decq %rax
    \\ jg 20027b
    \\ .balign 8
    \\ 20028:
    \\ addpd %xmm0, %xmm8
    \\ mulpd %xmm15, %xmm1
    \\ addpd %xmm2, %xmm9
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm4, %xmm10
    \\ mulpd %xmm15, %xmm5
    \\ addpd %xmm6, %xmm11
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm8
    \\ addpd %xmm3, %xmm9
    \\ addpd %xmm5, %xmm10
    \\ addpd %xmm7, %xmm11
    \\ movlpd 0 * 8(%rdx), %xmm0
    \\ movhpd 0 * 8(%rdx), %xmm0
    \\ movlpd 1 * 8(%rdx), %xmm1
    \\ movhpd 1 * 8(%rdx), %xmm1
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm2
    \\ movhpd 0 * 8(%rdx), %xmm2
    \\ movlpd 1 * 8(%rdx), %xmm3
    \\ movhpd 1 * 8(%rdx), %xmm3
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm4
    \\ movhpd 0 * 8(%rdx), %xmm4
    \\ movlpd 1 * 8(%rdx), %xmm5
    \\ movhpd 1 * 8(%rdx), %xmm5
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm6
    \\ movhpd 0 * 8(%rdx), %xmm6
    \\ movlpd 1 * 8(%rdx), %xmm7
    \\ movhpd 1 * 8(%rdx), %xmm7
    \\ addq %r8, %rdx
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm14, %xmm6
    \\ movlpd %xmm8, 0 * 8(%r11)
    \\ movhpd %xmm8, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm9, 0 * 8(%r11)
    \\ movhpd %xmm9, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm10, 0 * 8(%r11)
    \\ movhpd %xmm10, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm11, 0 * 8(%r11)
    \\ movhpd %xmm11, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlps 0 * 8(%r9), %xmm8
    \\ movhpd 1 * 8(%r9), %xmm8
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm9
    \\ movhpd 1 * 8(%r9), %xmm9
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm10
    \\ movhpd 1 * 8(%r9), %xmm10
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm11
    \\ movhpd 1 * 8(%r9), %xmm11
    \\ addq %r10, %r9
    \\ addpd %xmm0, %xmm8
    \\ mulpd %xmm15, %xmm1
    \\ addpd %xmm2, %xmm9
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm4, %xmm10
    \\ mulpd %xmm15, %xmm5
    \\ addpd %xmm6, %xmm11
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm8
    \\ addpd %xmm3, %xmm9
    \\ addpd %xmm5, %xmm10
    \\ addpd %xmm7, %xmm11
    \\ movlpd %xmm8, 0 * 8(%r11)
    \\ movhpd %xmm8, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm9, 0 * 8(%r11)
    \\ movhpd %xmm9, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm10, 0 * 8(%r11)
    \\ movhpd %xmm10, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm11, 0 * 8(%r11)
    \\ movhpd %xmm11, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ .balign 8
    \\ 20029:
    \\ movq %rcx, %rax
    \\ andq $4, %rax
    \\ jle 20030f
    \\ movlpd 0 * 8(%rdx), %xmm0
    \\ movhpd 0 * 8(%rdx), %xmm0
    \\ movlpd 1 * 8(%rdx), %xmm1
    \\ movhpd 1 * 8(%rdx), %xmm1
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm2
    \\ movhpd 0 * 8(%rdx), %xmm2
    \\ movlpd 1 * 8(%rdx), %xmm3
    \\ movhpd 1 * 8(%rdx), %xmm3
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm4
    \\ movhpd 0 * 8(%rdx), %xmm4
    \\ movlpd 1 * 8(%rdx), %xmm5
    \\ movhpd 1 * 8(%rdx), %xmm5
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm6
    \\ movhpd 0 * 8(%rdx), %xmm6
    \\ movlpd 1 * 8(%rdx), %xmm7
    \\ movhpd 1 * 8(%rdx), %xmm7
    \\ addq %r8, %rdx
    \\ movlps 0 * 8(%r9), %xmm8
    \\ movhpd 1 * 8(%r9), %xmm8
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm9
    \\ movhpd 1 * 8(%r9), %xmm9
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm10
    \\ movhpd 1 * 8(%r9), %xmm10
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm11
    \\ movhpd 1 * 8(%r9), %xmm11
    \\ addq %r10, %r9
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm14, %xmm6
    \\ addpd %xmm0, %xmm8
    \\ mulpd %xmm15, %xmm1
    \\ addpd %xmm2, %xmm9
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm4, %xmm10
    \\ mulpd %xmm15, %xmm5
    \\ addpd %xmm6, %xmm11
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm8
    \\ addpd %xmm3, %xmm9
    \\ addpd %xmm5, %xmm10
    \\ addpd %xmm7, %xmm11
    \\ movlpd %xmm8, 0 * 8(%r11)
    \\ movhpd %xmm8, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm9, 0 * 8(%r11)
    \\ movhpd %xmm9, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm10, 0 * 8(%r11)
    \\ movhpd %xmm10, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm11, 0 * 8(%r11)
    \\ movhpd %xmm11, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ .balign 8
    \\ 20030:
    \\ movq %rcx, %rax
    \\ andq $2, %rax
    \\ jle 20031f
    \\ movlpd 0 * 8(%rdx), %xmm0
    \\ movhpd 0 * 8(%rdx), %xmm0
    \\ movlpd 1 * 8(%rdx), %xmm1
    \\ movhpd 1 * 8(%rdx), %xmm1
    \\ addq %r8, %rdx
    \\ movlpd 0 * 8(%rdx), %xmm2
    \\ movhpd 0 * 8(%rdx), %xmm2
    \\ movlpd 1 * 8(%rdx), %xmm3
    \\ movhpd 1 * 8(%rdx), %xmm3
    \\ addq %r8, %rdx
    \\ movlps 0 * 8(%r9), %xmm8
    \\ movhpd 1 * 8(%r9), %xmm8
    \\ addq %r10, %r9
    \\ movlps 0 * 8(%r9), %xmm9
    \\ movhpd 1 * 8(%r9), %xmm9
    \\ addq %r10, %r9
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm0, %xmm8
    \\ addpd %xmm2, %xmm9
    \\ addpd %xmm1, %xmm8
    \\ addpd %xmm3, %xmm9
    \\ movlpd %xmm8, 0 * 8(%r11)
    \\ movhpd %xmm8, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ movlpd %xmm9, 0 * 8(%r11)
    \\ movhpd %xmm9, 1 * 8(%r11)
    \\ addq %r10, %r11
    \\ .balign 8
    \\ 20031:
    \\ movq %rcx, %rax
    \\ andq $1, %rax
    \\ jle 20033f
    \\ 20032:
    \\ movlpd 0 * 8(%rdx), %xmm0
    \\ movhpd 0 * 8(%rdx), %xmm0
    \\ movlpd 1 * 8(%rdx), %xmm1
    \\ movhpd 1 * 8(%rdx), %xmm1
    \\ movlps 0 * 8(%r9), %xmm8
    \\ movhpd 1 * 8(%r9), %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ addpd %xmm0, %xmm8
    \\ addpd %xmm1, %xmm8
    \\ movlpd %xmm8, 0 * 8(%r11)
    \\ movhpd %xmm8, 1 * 8(%r11)
    \\ addq %r8, %rdx
    \\ addq %r10, %r9
    \\ addq %r10, %r11
    \\ decq %rax
    \\ jg 20032b
    \\ .balign 8
    \\ 20033:
    \\ xorq %rax, %rax
    \\ movups 0(%rsp), %xmm6
    \\ movups 16(%rsp), %xmm7
    \\ movups 32(%rsp), %xmm8
    \\ movups 48(%rsp), %xmm9
    \\ movups 64(%rsp), %xmm10
    \\ movups 80(%rsp), %xmm11
    \\ movups 96(%rsp), %xmm12
    \\ movups 112(%rsp), %xmm13
    \\ movups 128(%rsp), %xmm14
    \\ movups 144(%rsp), %xmm15
    \\ addq $256, %rsp
    \\ ret
    \\ .endm
    \\
;
