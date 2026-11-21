//! Section 2 of 2 of the text of `daxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/axpy_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20010:
    \\ movaps -9 * 8(%rdx), %xmm4
    \\ movaps -7 * 8(%rdx), %xmm5
    \\ shufpd $1, %xmm1, %xmm0
    \\ mulpd %xmm15, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps -1 * 8(%rdx), %xmm0
    \\ shufpd $1, %xmm2, %xmm1
    \\ mulpd %xmm15, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movaps -5 * 8(%rdx), %xmm6
    \\ movaps -3 * 8(%rdx), %xmm7
    \\ shufpd $1, %xmm3, %xmm2
    \\ mulpd %xmm15, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ shufpd $1, %xmm4, %xmm3
    \\ mulpd %xmm15, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ shufpd $1, %xmm5, %xmm4
    \\ mulpd %xmm15, %xmm4
    \\ addpd -8 * 8(%r9), %xmm4
    \\ movaps %xmm4, -8 * 8(%r9)
    \\ shufpd $1, %xmm6, %xmm5
    \\ mulpd %xmm15, %xmm5
    \\ addpd -6 * 8(%r9), %xmm5
    \\ movaps %xmm5, -6 * 8(%r9)
    \\ shufpd $1, %xmm7, %xmm6
    \\ mulpd %xmm15, %xmm6
    \\ addpd -4 * 8(%r9), %xmm6
    \\ movaps %xmm6, -4 * 8(%r9)
    \\ shufpd $1, %xmm0, %xmm7
    \\ mulpd %xmm15, %xmm7
    \\ addpd -2 * 8(%r9), %xmm7
    \\ movaps %xmm7, -2 * 8(%r9)
    \\ subq $-16 * 8, %rdx
    \\ subq $-16 * 8, %r9
    \\ .balign 8
    \\ 20011:
    \\ movq %rcx, %rax
    \\ andq $8, %rax
    \\ jle 20012f
    \\ .balign 8
    \\ movaps -15 * 8(%rdx), %xmm1
    \\ movaps -13 * 8(%rdx), %xmm2
    \\ movaps -11 * 8(%rdx), %xmm3
    \\ movaps -9 * 8(%rdx), %xmm8
    \\ shufpd $1, %xmm1, %xmm0
    \\ mulpd %xmm15, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ shufpd $1, %xmm2, %xmm1
    \\ mulpd %xmm15, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ shufpd $1, %xmm3, %xmm2
    \\ mulpd %xmm15, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ shufpd $1, %xmm8, %xmm3
    \\ mulpd %xmm15, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ movaps %xmm8, %xmm0
    \\ addq $8 * 8, %rdx
    \\ addq $8 * 8, %r9
    \\ .balign 8
    \\ 20012:
    \\ movq %rcx, %rax
    \\ andq $4, %rax
    \\ jle 20013f
    \\ .balign 8
    \\ movaps -15 * 8(%rdx), %xmm1
    \\ movaps -13 * 8(%rdx), %xmm2
    \\ shufpd $1, %xmm1, %xmm0
    \\ shufpd $1, %xmm2, %xmm1
    \\ mulpd %xmm15, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ addpd -16 * 8(%r9), %xmm0
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movaps %xmm2, %xmm0
    \\ addq $4 * 8, %rdx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20013:
    \\ movq %rcx, %rax
    \\ andq $2, %rax
    \\ jle 20014f
    \\ .balign 8
    \\ movaps -15 * 8(%rdx), %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ mulpd %xmm15, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ addq $2 * 8, %rdx
    \\ addq $2 * 8, %r9
    \\ .balign 8
    \\ 20014:
    \\ movq %rcx, %rax
    \\ andq $1, %rax
    \\ jle 20015f
    \\ .balign 8
    \\ movsd -16 * 8(%rdx), %xmm0
    \\ mulsd %xmm15, %xmm0
    \\ addsd -16 * 8(%r9), %xmm0
    \\ movsd %xmm0, -16 * 8(%r9)
    \\ .balign 8
    \\ 20015:
    \\ xorq %rax,%rax
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
    \\ .balign 8
    \\ 20016:
    \\ movq %r9, %r11
    \\ movq %rcx, %rax
    \\ cmpq $0, %r8
    \\ je 20019f
    \\ cmpq $0, %r10
    \\ je 20019f
    \\ sarq $3, %rax
    \\ jle 20018f
    \\ .balign 8
    \\ 20017:
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ addq %r8, %rdx
    \\ .if \bp3
    \\ movhpd 0 * 8(%rdx), %xmm0
    \\ .endif
    \\ .if \bp4
    \\ movhps 0 * 8(%rdx), %xmm0
    \\ .endif
    \\ addq %r8, %rdx
    \\ mulpd %xmm15, %xmm0
    \\ movsd 0 * 8(%r11), %xmm6
    \\ addq %r10, %r11
    \\ .if \bp3
    \\ movhpd 0 * 8(%r11), %xmm6
    \\ .endif
    \\ .if \bp4
    \\ movhps 0 * 8(%r11), %xmm6
    \\ .endif
    \\ addq %r10, %r11
    \\ addpd %xmm6, %xmm0
    \\ movsd 0 * 8(%rdx), %xmm1
    \\ addq %r8, %rdx
    \\ .if \bp3
    \\ movhpd 0 * 8(%rdx), %xmm1
    \\ .endif
    \\ .if \bp4
    \\ movhps 0 * 8(%rdx), %xmm1
    \\ .endif
    \\ addq %r8, %rdx
    \\ mulpd %xmm15, %xmm1
    \\ movsd 0 * 8(%r11), %xmm6
    \\ addq %r10, %r11
    \\ .if \bp3
    \\ movhpd 0 * 8(%r11), %xmm6
    \\ .endif
    \\ .if \bp4
    \\ movhps 0 * 8(%r11), %xmm6
    \\ .endif
    \\ addq %r10, %r11
    \\ addpd %xmm6, %xmm1
    \\ movsd 0 * 8(%rdx), %xmm2
    \\ addq %r8, %rdx
    \\ .if \bp3
    \\ movhpd 0 * 8(%rdx), %xmm2
    \\ .endif
    \\ .if \bp4
    \\ movhps 0 * 8(%rdx), %xmm2
    \\ .endif
    \\ addq %r8, %rdx
    \\ mulpd %xmm15, %xmm2
    \\ movsd 0 * 8(%r11), %xmm6
    \\ addq %r10, %r11
    \\ .if \bp3
    \\ movhpd 0 * 8(%r11), %xmm6
    \\ .endif
    \\ .if \bp4
    \\ movhps 0 * 8(%r11), %xmm6
    \\ .endif
    \\ addq %r10, %r11
    \\ addpd %xmm6, %xmm2
    \\ movsd 0 * 8(%rdx), %xmm3
    \\ addq %r8, %rdx
    \\ .if \bp3
    \\ movhpd 0 * 8(%rdx), %xmm3
    \\ .endif
    \\ .if \bp4
    \\ movhps 0 * 8(%rdx), %xmm3
    \\ .endif
    \\ addq %r8, %rdx
    \\ mulpd %xmm15, %xmm3
    \\ movsd 0 * 8(%r11), %xmm6
    \\ addq %r10, %r11
    \\ .if \bp3
    \\ movhpd 0 * 8(%r11), %xmm6
    \\ .endif
    \\ .if \bp4
    \\ movhps 0 * 8(%r11), %xmm6
    \\ .endif
    \\ addq %r10, %r11
    \\ addpd %xmm6, %xmm3
    \\ movsd %xmm0, 0 * 8(%r9)
    \\ addq %r10, %r9
    \\ .if \bp3
    \\ movhpd %xmm0, 0 * 8(%r9)
    \\ .endif
    \\ .if \bp4
    \\ movhps %xmm0, 0 * 8(%r9)
    \\ .endif
    \\ addq %r10, %r9
    \\ movsd %xmm1, 0 * 8(%r9)
    \\ addq %r10, %r9
    \\ .if \bp3
    \\ movhpd %xmm1, 0 * 8(%r9)
    \\ .endif
    \\ .if \bp4
    \\ movhps %xmm1, 0 * 8(%r9)
    \\ .endif
    \\ addq %r10, %r9
    \\ movsd %xmm2, 0 * 8(%r9)
    \\ addq %r10, %r9
    \\ .if \bp3
    \\ movhpd %xmm2, 0 * 8(%r9)
    \\ .endif
    \\ .if \bp4
    \\ movhps %xmm2, 0 * 8(%r9)
    \\ .endif
    \\ addq %r10, %r9
    \\ movsd %xmm3, 0 * 8(%r9)
    \\ addq %r10, %r9
    \\ .if \bp3
    \\ movhpd %xmm3, 0 * 8(%r9)
    \\ .endif
    \\ .if \bp4
    \\ movhps %xmm3, 0 * 8(%r9)
    \\ .endif
    \\ addq %r10, %r9
    \\ decq %rax
    \\ jg 20017b
    \\ .balign 8
    \\ 20018:
    \\ movq %rcx, %rax
    \\ andq $7, %rax
    \\ jle 20020f
    \\ .balign 8
    \\ 20019:
    \\ movsd (%rdx), %xmm0
    \\ addq %r8, %rdx
    \\ mulsd %xmm15, %xmm0
    \\ addsd (%r9), %xmm0
    \\ movsd %xmm0, (%r9)
    \\ addq %r10, %r9
    \\ decq %rax
    \\ jg 20019b
    \\ .balign 8
    \\ 20020:
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
