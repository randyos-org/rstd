//! Section 2 of 3 of the text of `zdrot`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zrot_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20008:
    \\ testq $7, %rcx
    \\ jle 20027f
    \\ testq $4, %rcx
    \\ jle 20009f
    \\ movaps 1 * 8(%r9), %xmm3
    \\ movaps 3 * 8(%r9), %xmm8
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps 2 * 8(%rdx), %xmm2
    \\ shufpd $1, %xmm3, %xmm1
    \\ shufpd $1, %xmm8, %xmm3
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 0 * 8(%rdx)
    \\ movaps %xmm2, 2 * 8(%rdx)
    \\ movlps %xmm4, 0 * 8(%r9)
    \\ movhps %xmm4, 1 * 8(%r9)
    \\ movlps %xmm6, 2 * 8(%r9)
    \\ movhps %xmm6, 3 * 8(%r9)
    \\ movaps 5 * 8(%r9), %xmm9
    \\ movaps 7 * 8(%r9), %xmm1
    \\ movaps 4 * 8(%rdx), %xmm0
    \\ movaps 6 * 8(%rdx), %xmm2
    \\ shufpd $1, %xmm9, %xmm8
    \\ shufpd $1, %xmm1, %xmm9
    \\ movaps %xmm8, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm9, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm9
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm8, %xmm0
    \\ addpd %xmm9, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 4 * 8(%rdx)
    \\ movaps %xmm2, 6 * 8(%rdx)
    \\ movlps %xmm4, 4 * 8(%r9)
    \\ movhps %xmm4, 5 * 8(%r9)
    \\ movlps %xmm6, 6 * 8(%r9)
    \\ movhps %xmm6, 7 * 8(%r9)
    \\ addq $8 * 8, %rdx
    \\ addq $8 * 8, %r9
    \\ .balign 8
    \\ 20009:
    \\ testq $2, %rcx
    \\ jle 20010f
    \\ movaps 1 * 8(%r9), %xmm3
    \\ movaps 3 * 8(%r9), %xmm8
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps 2 * 8(%rdx), %xmm2
    \\ shufpd $1, %xmm3, %xmm1
    \\ shufpd $1, %xmm8, %xmm3
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 0 * 8(%rdx)
    \\ movaps %xmm2, 2 * 8(%rdx)
    \\ movlps %xmm4, 0 * 8(%r9)
    \\ movhps %xmm4, 1 * 8(%r9)
    \\ movlps %xmm6, 2 * 8(%r9)
    \\ movhps %xmm6, 3 * 8(%r9)
    \\ movaps %xmm8, %xmm1
    \\ addq $4 * 8, %rdx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20010:
    \\ testq $1, %rcx
    \\ jle 20027f
    \\ movaps 1 * 8(%r9), %xmm4
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ shufpd $1, %xmm4, %xmm1
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm1, %xmm0
    \\ subpd %xmm3, %xmm2
    \\ movaps %xmm0, 0 * 8(%rdx)
    \\ movlps %xmm2, 0 * 8(%r9)
    \\ movhps %xmm2, 1 * 8(%r9)
    \\ jmp 20027f
    \\ .balign 8
    \\ 20011:
    \\ testq $8, %r9
    \\ jne 20016f
    \\ movaps -1 * 8(%rdx), %xmm0
    \\ movq %rcx, %rax
    \\ sarq $3, %rax
    \\ jle 20013f
    \\ .balign 8
    \\ 20012:
    \\ movaps 1 * 8(%rdx), %xmm2
    \\ movaps 3 * 8(%rdx), %xmm8
    \\ movaps 0 * 8(%r9), %xmm1
    \\ movaps 2 * 8(%r9), %xmm3
    \\ shufpd $1, %xmm2, %xmm0
    \\ shufpd $1, %xmm8, %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movlps %xmm0, 0 * 8(%rdx)
    \\ movhps %xmm0, 1 * 8(%rdx)
    \\ movlps %xmm2, 2 * 8(%rdx)
    \\ movhps %xmm2, 3 * 8(%rdx)
    \\ movaps %xmm4, 0 * 8(%r9)
    \\ movaps %xmm6, 2 * 8(%r9)
    \\ movaps 5 * 8(%rdx), %xmm2
    \\ movaps 7 * 8(%rdx), %xmm0
    \\ movaps 4 * 8(%r9), %xmm1
    \\ movaps 6 * 8(%r9), %xmm3
    \\ shufpd $1, %xmm2, %xmm8
    \\ shufpd $1, %xmm0, %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm8, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm8
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm8
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movlps %xmm8, 4 * 8(%rdx)
    \\ movhps %xmm8, 5 * 8(%rdx)
    \\ movlps %xmm2, 6 * 8(%rdx)
    \\ movhps %xmm2, 7 * 8(%rdx)
    \\ movaps %xmm4, 4 * 8(%r9)
    \\ movaps %xmm6, 6 * 8(%r9)
    \\ movaps 9 * 8(%rdx), %xmm2
    \\ movaps 11 * 8(%rdx), %xmm8
    \\ movaps 8 * 8(%r9), %xmm1
    \\ movaps 10 * 8(%r9), %xmm3
    \\ shufpd $1, %xmm2, %xmm0
    \\ shufpd $1, %xmm8, %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movlps %xmm0, 8 * 8(%rdx)
    \\ movhps %xmm0, 9 * 8(%rdx)
    \\ movlps %xmm2, 10 * 8(%rdx)
    \\ movhps %xmm2, 11 * 8(%rdx)
    \\ movaps %xmm4, 8 * 8(%r9)
    \\ movaps %xmm6, 10 * 8(%r9)
    \\ movaps 13 * 8(%rdx), %xmm2
    \\ movaps 15 * 8(%rdx), %xmm0
    \\ movaps 12 * 8(%r9), %xmm1
    \\ movaps 14 * 8(%r9), %xmm3
    \\ shufpd $1, %xmm2, %xmm8
    \\ shufpd $1, %xmm0, %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm8, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm8
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm8
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movlps %xmm8, 12 * 8(%rdx)
    \\ movhps %xmm8, 13 * 8(%rdx)
    \\ movlps %xmm2, 14 * 8(%rdx)
    \\ movhps %xmm2, 15 * 8(%rdx)
    \\ movaps %xmm4, 12 * 8(%r9)
    \\ movaps %xmm6, 14 * 8(%r9)
    \\ addq $16 * 8, %r9
    \\ addq $16 * 8, %rdx
    \\ decq %rax
    \\ jg 20012b
    \\ .balign 8
    \\ 20013:
    \\ testq $7, %rcx
    \\ jle 20027f
    \\ testq $4, %rcx
    \\ jle 20014f
    \\ movaps 1 * 8(%rdx), %xmm2
    \\ movaps 3 * 8(%rdx), %xmm8
    \\ movaps 0 * 8(%r9), %xmm1
    \\ movaps 2 * 8(%r9), %xmm3
    \\ shufpd $1, %xmm2, %xmm0
    \\ shufpd $1, %xmm8, %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movlps %xmm0, 0 * 8(%rdx)
    \\ movhps %xmm0, 1 * 8(%rdx)
    \\ movlps %xmm2, 2 * 8(%rdx)
    \\ movhps %xmm2, 3 * 8(%rdx)
    \\ movaps %xmm4, 0 * 8(%r9)
    \\ movaps %xmm6, 2 * 8(%r9)
    \\ movaps 5 * 8(%rdx), %xmm2
    \\ movaps 7 * 8(%rdx), %xmm0
    \\ movaps 4 * 8(%r9), %xmm1
    \\ movaps 6 * 8(%r9), %xmm3
    \\ shufpd $1, %xmm2, %xmm8
    \\ shufpd $1, %xmm0, %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm8, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm8
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm8
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movlps %xmm8, 4 * 8(%rdx)
    \\ movhps %xmm8, 5 * 8(%rdx)
    \\ movlps %xmm2, 6 * 8(%rdx)
    \\ movhps %xmm2, 7 * 8(%rdx)
    \\ movaps %xmm4, 4 * 8(%r9)
    \\ movaps %xmm6, 6 * 8(%r9)
    \\ addq $8 * 8, %r9
    \\ addq $8 * 8, %rdx
    \\ .balign 8
    \\ 20014:
    \\ testq $2, %rcx
    \\ jle 20015f
    \\ movaps 1 * 8(%rdx), %xmm2
    \\ movaps 3 * 8(%rdx), %xmm8
    \\ movaps 0 * 8(%r9), %xmm1
    \\ movaps 2 * 8(%r9), %xmm3
    \\ shufpd $1, %xmm2, %xmm0
    \\ shufpd $1, %xmm8, %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movlps %xmm0, 0 * 8(%rdx)
    \\ movhps %xmm0, 1 * 8(%rdx)
    \\ movlps %xmm2, 2 * 8(%rdx)
    \\ movhps %xmm2, 3 * 8(%rdx)
    \\ movaps %xmm4, 0 * 8(%r9)
    \\ movaps %xmm6, 2 * 8(%r9)
    \\ movaps %xmm8, %xmm0
    \\ addq $4 * 8, %r9
    \\ addq $4 * 8, %rdx
    \\ .balign 8
    \\ 20015:
    \\ testq $1, %rcx
    \\ jle 20027f
    \\ movaps 1 * 8(%rdx), %xmm4
    \\ movaps 0 * 8(%r9), %xmm1
    \\ shufpd $1, %xmm4, %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm1, %xmm0
    \\ subpd %xmm3, %xmm2
    \\ movlps %xmm0, 0 * 8(%rdx)
    \\ movhps %xmm0, 1 * 8(%rdx)
    \\ movaps %xmm2, 0 * 8(%r9)
    \\ jmp 20027f
    \\ .balign 8
    \\ 20016:
    \\ movsd 0 * 8(%r9), %xmm1
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulsd %xmm14, %xmm0
    \\ mulsd %xmm15, %xmm1
    \\ mulsd %xmm14, %xmm2
    \\ mulsd %xmm15, %xmm3
    \\ addsd %xmm1, %xmm0
    \\ subsd %xmm3, %xmm2
    \\ movsd %xmm0, 0 * 8(%rdx)
    \\ movsd %xmm2, 0 * 8(%r9)
    \\ addq $1 * 8, %r9
    \\ addq $1 * 8, %rdx
    \\ decq %rcx
    \\ jle 20022f
    \\ movq %rcx, %rax
    \\ sarq $3, %rax
    \\ jle 20019f
    \\ movaps 0 * 8(%r9), %xmm1
    \\ movaps 2 * 8(%r9), %xmm3
    \\ movaps 4 * 8(%r9), %xmm9
    \\ movaps 6 * 8(%r9), %xmm11
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps 2 * 8(%rdx), %xmm2
    \\ movaps 4 * 8(%rdx), %xmm8
    \\ movaps 6 * 8(%rdx), %xmm10
    \\ decq %rax
    \\ jle 20018f
    \\
;
