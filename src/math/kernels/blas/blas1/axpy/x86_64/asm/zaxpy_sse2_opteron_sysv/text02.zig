//! Section 2 of 3 of the text of `zaxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zaxpy_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20011:
    \\ testq $8, %rcx
    \\ jne 20018f
    \\ movaps -16 * 8(%rcx), %xmm1
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ xorps %xmm0, %xmm0
    \\ shufpd $1, %xmm1, %xmm0
    \\ xorps %xmm4, %xmm4
    \\ movhps -16 * 8(%r9), %xmm4
    \\ addpd %xmm0, %xmm4
    \\ movhps %xmm4, -16 * 8(%r9)
    \\ movaps %xmm1, %xmm0
    \\ addq $2 * 8, %rcx
    \\ addq $1 * 8, %r9
    \\ decq %rdi
    \\ jle 20017f
    \\ movq %rdi, %rax
    \\ sarq $3, %rax
    \\ jle 20014f
    \\ movaps -16 * 8(%rcx), %xmm1
    \\ movaps -14 * 8(%rcx), %xmm2
    \\ movaps -12 * 8(%rcx), %xmm3
    \\ decq %rax
    \\ jle 20013f
    \\ .balign 8
    \\ 20012:
    \\ prefetchw ((128 * 3) + 0) - 128(%r9)
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps -10 * 8(%rcx), %xmm0
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ shufpd $1, %xmm2, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movaps -8 * 8(%rcx), %xmm1
    \\ prefetch ((128 * 3) + 0) - 128(%rcx)
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ shufpd $1, %xmm3, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ movaps -6 * 8(%rcx), %xmm2
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ shufpd $1, %xmm0, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ movaps -4 * 8(%rcx), %xmm3
    \\ prefetchw ((128 * 3) + 64) - 128(%r9)
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -8 * 8(%r9), %xmm0
    \\ movaps %xmm0, -8 * 8(%r9)
    \\ movaps -2 * 8(%rcx), %xmm0
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ shufpd $1, %xmm2, %xmm1
    \\ addpd -6 * 8(%r9), %xmm1
    \\ movaps %xmm1, -6 * 8(%r9)
    \\ movaps 0 * 8(%rcx), %xmm1
    \\ prefetch ((128 * 3) + 64) - 128(%rcx)
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ shufpd $1, %xmm3, %xmm2
    \\ addpd -4 * 8(%r9), %xmm2
    \\ movaps %xmm2, -4 * 8(%r9)
    \\ movaps 2 * 8(%rcx), %xmm2
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ shufpd $1, %xmm0, %xmm3
    \\ addpd -2 * 8(%r9), %xmm3
    \\ movaps %xmm3, -2 * 8(%r9)
    \\ movaps 4 * 8(%rcx), %xmm3
    \\ subq $-16 * 8, %rcx
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
    \\ movaps -10 * 8(%rcx), %xmm0
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ shufpd $1, %xmm2, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movaps -8 * 8(%rcx), %xmm1
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ shufpd $1, %xmm3, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ movaps -6 * 8(%rcx), %xmm2
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ shufpd $1, %xmm0, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ movaps -4 * 8(%rcx), %xmm3
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -8 * 8(%r9), %xmm0
    \\ movaps %xmm0, -8 * 8(%r9)
    \\ movaps -2 * 8(%rcx), %xmm0
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
    \\ subq $-16 * 8, %rcx
    \\ subq $-16 * 8, %r9
    \\ .balign 8
    \\ 20014:
    \\ movq %rdi, %rax
    \\ andq $4, %rax
    \\ jle 20015f
    \\ movaps -16 * 8(%rcx), %xmm1
    \\ movaps -14 * 8(%rcx), %xmm2
    \\ movaps -12 * 8(%rcx), %xmm3
    \\ movaps -10 * 8(%rcx), %xmm4
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
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ shufpd $1, %xmm3, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ pshufd $0x4e, %xmm4, %xmm8
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ shufpd $1, %xmm4, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ movaps %xmm4, %xmm0
    \\ addq $8 * 8, %rcx
    \\ addq $8 * 8, %r9
    \\ .balign 8
    \\ 20015:
    \\ movq %rdi, %rax
    \\ andq $2, %rax
    \\ jle 20016f
    \\ movaps -16 * 8(%rcx), %xmm1
    \\ movaps -14 * 8(%rcx), %xmm2
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
    \\ addq $4 * 8, %rcx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20016:
    \\ movq %rdi, %rax
    \\ andq $1, %rax
    \\ jle 20017f
    \\ movaps -16 * 8(%rcx), %xmm1
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movaps %xmm1, %xmm0
    \\ addq $2 * 8, %rcx
    \\ addq $2 * 8, %r9
    \\ .balign 8
    \\ 20017:
    \\ shufpd $1, %xmm0, %xmm0
    \\ addsd -16 * 8(%r9), %xmm0
    \\ movlps %xmm0, -16 * 8(%r9)
    \\ jmp 20033f
    \\ .balign 8
    \\ 20018:
    \\ movlps -16 * 8(%rcx), %xmm1
    \\ movhps -15 * 8(%rcx), %xmm1
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ xorps %xmm0, %xmm0
    \\ shufpd $1, %xmm1, %xmm0
    \\ xorps %xmm4, %xmm4
    \\ movhps -16 * 8(%r9), %xmm4
    \\ addpd %xmm0, %xmm4
    \\ movhps %xmm4, -16 * 8(%r9)
    \\ movaps %xmm1, %xmm0
    \\ addq $2 * 8, %rcx
    \\ addq $1 * 8, %r9
    \\ decq %rdi
    \\ jle 20024f
    \\ movq %rdi, %rax
    \\ sarq $3, %rax
    \\ jle 20021f
    \\ movlps -16 * 8(%rcx), %xmm1
    \\ movhps -15 * 8(%rcx), %xmm1
    \\ movlps -14 * 8(%rcx), %xmm2
    \\ movhps -13 * 8(%rcx), %xmm2
    \\ movlps -12 * 8(%rcx), %xmm3
    \\ movhps -11 * 8(%rcx), %xmm3
    \\ decq %rax
    \\ jle 20020f
    \\ .balign 8
    \\ 20019:
    \\ prefetchw ((128 * 3) + 0) - 128(%r9)
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movlps -10 * 8(%rcx), %xmm0
    \\ movhps -9 * 8(%rcx), %xmm0
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ shufpd $1, %xmm2, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movlps -8 * 8(%rcx), %xmm1
    \\ movhps -7 * 8(%rcx), %xmm1
    \\ prefetch ((128 * 3) + 0) - 128(%rcx)
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ shufpd $1, %xmm3, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ movlps -6 * 8(%rcx), %xmm2
    \\ movhps -5 * 8(%rcx), %xmm2
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ shufpd $1, %xmm0, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ movlps -4 * 8(%rcx), %xmm3
    \\ movhps -3 * 8(%rcx), %xmm3
    \\ prefetchw ((128 * 3) + 64) - 128(%r9)
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -8 * 8(%r9), %xmm0
    \\ movaps %xmm0, -8 * 8(%r9)
    \\ movlps -2 * 8(%rcx), %xmm0
    \\ movhps -1 * 8(%rcx), %xmm0
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ shufpd $1, %xmm2, %xmm1
    \\ addpd -6 * 8(%r9), %xmm1
    \\ movaps %xmm1, -6 * 8(%r9)
    \\ movlps 0 * 8(%rcx), %xmm1
    \\ movhps 1 * 8(%rcx), %xmm1
    \\ prefetch ((128 * 3) + 64) - 128(%rcx)
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ shufpd $1, %xmm3, %xmm2
    \\ addpd -4 * 8(%r9), %xmm2
    \\ movaps %xmm2, -4 * 8(%r9)
    \\ movlps 2 * 8(%rcx), %xmm2
    \\ movhps 3 * 8(%rcx), %xmm2
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ shufpd $1, %xmm0, %xmm3
    \\ addpd -2 * 8(%r9), %xmm3
    \\ movaps %xmm3, -2 * 8(%r9)
    \\ movlps 4 * 8(%rcx), %xmm3
    \\ movhps 5 * 8(%rcx), %xmm3
    \\ subq $-16 * 8, %rcx
    \\ subq $-16 * 8, %r9
    \\ decq %rax
    \\ jg 20019b
    \\ .balign 8
    \\ 20020:
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -16 * 8(%r9), %xmm0
    \\ movaps %xmm0, -16 * 8(%r9)
    \\ movlps -10 * 8(%rcx), %xmm0
    \\ movhps -9 * 8(%rcx), %xmm0
    \\ pshufd $0x4e, %xmm2, %xmm8
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm2
    \\ shufpd $1, %xmm2, %xmm1
    \\ addpd -14 * 8(%r9), %xmm1
    \\ movaps %xmm1, -14 * 8(%r9)
    \\ movlps -8 * 8(%rcx), %xmm1
    \\ movhps -7 * 8(%rcx), %xmm1
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ shufpd $1, %xmm3, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ movlps -6 * 8(%rcx), %xmm2
    \\ movhps -5 * 8(%rcx), %xmm2
    \\ pshufd $0x4e, %xmm0, %xmm8
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm0
    \\ shufpd $1, %xmm0, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ movlps -4 * 8(%rcx), %xmm3
    \\ movhps -3 * 8(%rcx), %xmm3
    \\ pshufd $0x4e, %xmm1, %xmm8
    \\ mulpd %xmm14, %xmm1
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm1
    \\ shufpd $1, %xmm1, %xmm0
    \\ addpd -8 * 8(%r9), %xmm0
    \\ movaps %xmm0, -8 * 8(%r9)
    \\ movlps -2 * 8(%rcx), %xmm0
    \\ movhps -1 * 8(%rcx), %xmm0
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
    \\ subq $-16 * 8, %rcx
    \\ subq $-16 * 8, %r9
    \\ .balign 8
    \\ 20021:
    \\ movq %rdi, %rax
    \\ andq $4, %rax
    \\ jle 20022f
    \\ movlps -16 * 8(%rcx), %xmm1
    \\ movhps -15 * 8(%rcx), %xmm1
    \\ movlps -14 * 8(%rcx), %xmm2
    \\ movhps -13 * 8(%rcx), %xmm2
    \\ movlps -12 * 8(%rcx), %xmm3
    \\ movhps -11 * 8(%rcx), %xmm3
    \\ movlps -10 * 8(%rcx), %xmm4
    \\ movhps -9 * 8(%rcx), %xmm4
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
    \\ pshufd $0x4e, %xmm3, %xmm8
    \\ mulpd %xmm14, %xmm3
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm3
    \\ shufpd $1, %xmm3, %xmm2
    \\ addpd -12 * 8(%r9), %xmm2
    \\ movaps %xmm2, -12 * 8(%r9)
    \\ pshufd $0x4e, %xmm4, %xmm8
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm8
    \\ addpd %xmm8, %xmm4
    \\ shufpd $1, %xmm4, %xmm3
    \\ addpd -10 * 8(%r9), %xmm3
    \\ movaps %xmm3, -10 * 8(%r9)
    \\ movaps %xmm4, %xmm0
    \\ addq $8 * 8, %rcx
    \\ addq $8 * 8, %r9
    \\
;
