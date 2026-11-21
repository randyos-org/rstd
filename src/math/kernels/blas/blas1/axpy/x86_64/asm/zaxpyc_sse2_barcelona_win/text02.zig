//! Section 2 of 5 of the text of `zaxpyc`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zaxpy_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20007:
    \\ movsd -8 * 8(%rdx), %xmm4
    \\ movhps -7 * 8(%rdx), %xmm4
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movsd -6 * 8(%rdx), %xmm5
    \\ movhps -5 * 8(%rdx), %xmm5
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movsd -4 * 8(%rdx), %xmm6
    \\ movhps -3 * 8(%rdx), %xmm6
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ movsd -2 * 8(%rdx), %xmm7
    \\ movhps -1 * 8(%rdx), %xmm7
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ pshufd $0x4e, %xmm4, %xmm8
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ addpd -8 * 8(%r9), %xmm4
    \\ movaps %xmm4, -8 * 8(%r9)
    \\ pshufd $0x4e, %xmm5, %xmm8
    \\ mulpd %xmm14, %xmm5
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm5
    \\ addpd -6 * 8(%r9), %xmm5
    \\ movaps %xmm5, -6 * 8(%r9)
    \\ pshufd $0x4e, %xmm6, %xmm8
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm6
    \\ addpd -4 * 8(%r9), %xmm6
    \\ movaps %xmm6, -4 * 8(%r9)
    \\ pshufd $0x4e, %xmm7, %xmm8
    \\ mulpd %xmm14, %xmm7
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm7
    \\ addpd -2 * 8(%r9), %xmm7
    \\ movaps %xmm7, -2 * 8(%r9)
    \\ subq $-16 * 8, %rdx
    \\ subq $-16 * 8, %r9
    \\ .balign 8
    \\ 20008:
    \\ movq %rcx, %rax
    \\ andq $4, %rax
    \\ jle 20009f
    \\ movsd -16 * 8(%rdx), %xmm0
    \\ movhps -15 * 8(%rdx), %xmm0
    \\ movsd -14 * 8(%rdx), %xmm1
    \\ movhps -13 * 8(%rdx), %xmm1
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movsd -12 * 8(%rdx), %xmm2
    \\ movhps -11 * 8(%rdx), %xmm2
    \\ movsd -10 * 8(%rdx), %xmm3
    \\ movhps -9 * 8(%rdx), %xmm3
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ addq $8 * 8, %rdx
    \\ addq $8 * 8, %r9
    \\ .balign 8
    \\ 20009:
    \\ movq %rcx, %rax
    \\ andq $2, %rax
    \\ jle 20010f
    \\ movsd -16 * 8(%rdx), %xmm0
    \\ movhps -15 * 8(%rdx), %xmm0
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movsd -14 * 8(%rdx), %xmm1
    \\ movhps -13 * 8(%rdx), %xmm1
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ addq $4 * 8, %rdx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20010:
    \\ movq %rcx, %rax
    \\ andq $1, %rax
    \\ jle 20033f
    \\ movsd -16 * 8(%rdx), %xmm0
    \\ movhps -15 * 8(%rdx), %xmm0
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ jmp 20033f
    \\ .balign 8
    \\ 20011:
    \\ testq $8, %rdx
    \\ jne 20018f
    \\ movaps -16 * 8(%rdx), %xmm1
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ .if \bp0
    \\ xorps %xmm0, %xmm0
    \\ .endif
    \\ .if \bp1
    \\ pxor %xmm0, %xmm0
    \\ .endif
    \\ shufpd $1, %xmm1, %xmm0
    \\ .if \bp0
    \\ xorps %xmm4, %xmm4
    \\ .endif
    \\ .if \bp1
    \\ pxor %xmm4, %xmm4
    \\ .endif
    \\ movhps -16 * 8(%r9), %xmm4
    \\ addpd %xmm0, %xmm4
    \\ movhps %xmm4, -16 * 8(%r9)
    \\ movaps %xmm1, %xmm0
    \\ addq $2 * 8, %rdx
    \\ addq $1 * 8, %r9
    \\ decq %rcx
    \\ jle 20017f
    \\ movq %rcx, %rax
    \\ sarq $3, %rax
    \\ jle 20014f
    \\ movaps -16 * 8(%rdx), %xmm1
    \\ movaps -14 * 8(%rdx), %xmm2
    \\ movaps -12 * 8(%rdx), %xmm3
    \\ decq %rax
    \\ jle 20013f
    \\ .balign 8
    \\ 20012:
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%r9)
    \\ .endif
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps -10 * 8(%rdx), %xmm0
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ shufpd $1, %xmm2, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movaps -8 * 8(%rdx), %xmm1
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp1
    \\ prefetcht0 ((128 * 10) + 0) - 128(%rdx)
    \\ .endif
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ shufpd $1, %xmm3, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ movaps -6 * 8(%rdx), %xmm2
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ shufpd $1, %xmm0, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ movaps -4 * 8(%rdx), %xmm3
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%r9)
    \\ .endif
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -8 * 8(%r9), %xmm0
    \\ movaps %xmm0, -8 * 8(%r9)
    \\ movaps -2 * 8(%rdx), %xmm0
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ shufpd $1, %xmm2, %xmm1
    \\ addpd -6 * 8(%r9), %xmm1
    \\ movaps %xmm1, -6 * 8(%r9)
    \\ movaps 0 * 8(%rdx), %xmm1
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rdx)
    \\ .endif
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ shufpd $1, %xmm3, %xmm2
    \\ addpd -4 * 8(%r9), %xmm2
    \\ movaps %xmm2, -4 * 8(%r9)
    \\ movaps 2 * 8(%rdx), %xmm2
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ shufpd $1, %xmm0, %xmm3
    \\ addpd -2 * 8(%r9), %xmm3
    \\ movaps %xmm3, -2 * 8(%r9)
    \\ movaps 4 * 8(%rdx), %xmm3
    \\ subq $-16 * 8, %rdx
    \\ subq $-16 * 8, %r9
    \\ decq %rax
    \\ jg 20012b
    \\ .balign 8
    \\ 20013:
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps -10 * 8(%rdx), %xmm0
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ shufpd $1, %xmm2, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movaps -8 * 8(%rdx), %xmm1
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ shufpd $1, %xmm3, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ movaps -6 * 8(%rdx), %xmm2
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ shufpd $1, %xmm0, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ movaps -4 * 8(%rdx), %xmm3
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -8 * 8(%r9), %xmm0
    \\ movaps %xmm0, -8 * 8(%r9)
    \\ movaps -2 * 8(%rdx), %xmm0
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ shufpd $1, %xmm2, %xmm1
    \\ addpd -6 * 8(%r9), %xmm1
    \\ movaps %xmm1, -6 * 8(%r9)
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ shufpd $1, %xmm3, %xmm2
    \\ addpd -4 * 8(%r9), %xmm2
    \\ movaps %xmm2, -4 * 8(%r9)
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ shufpd $1, %xmm0, %xmm3
    \\ addpd -2 * 8(%r9), %xmm3
    \\ movaps %xmm3, -2 * 8(%r9)
    \\ subq $-16 * 8, %rdx
    \\ subq $-16 * 8, %r9
    \\
;
