//! Section 1 of 3 of the text of `drot`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/rot_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1, bp2, bp3, bp4, bp5, bp6
    \\ .balign 512
    \\ movq 40(%rsp), %r10
    \\ .if \bp0
    \\ movsd 48(%rsp), %xmm0
    \\ movsd 56(%rsp), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 48(%rsp), %xmm0
    \\ movlps 56(%rsp), %xmm1
    \\ .endif
    \\ subq $256, %rsp
    \\ movups %xmm6, 0(%rsp)
    \\ movups %xmm7, 16(%rsp)
    \\ movups %xmm8, 32(%rsp)
    \\ movups %xmm9, 48(%rsp)
    \\ movups %xmm10, 64(%rsp)
    \\ movups %xmm11, 80(%rsp)
    \\ movups %xmm12, 96(%rsp)
    \\ movups %xmm13, 112(%rsp)
    \\ movups %xmm14, 128(%rsp)
    \\ movups %xmm15, 144(%rsp)
    \\ leaq (, %r8, 8), %r8
    \\ leaq (, %r10, 8), %r10
    \\ pshufd $0x44, %xmm0, %xmm14
    \\ pshufd $0x44, %xmm1, %xmm15
    \\ cmpq $0, %rcx
    \\ jle 20017f
    \\ cmpq $8, %r8
    \\ jne 20013f
    \\ cmpq $8, %r10
    \\ jne 20013f
    \\ testq $8, %rdx
    \\ je 20000f
    \\ .if \bp0
    \\ movsd 0 * 8(%r9), %xmm1
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 8(%r9), %xmm1
    \\ movlps 0 * 8(%rdx), %xmm0
    \\ .endif
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulsd %xmm14, %xmm0
    \\ mulsd %xmm15, %xmm1
    \\ mulsd %xmm14, %xmm2
    \\ mulsd %xmm15, %xmm3
    \\ addsd %xmm1, %xmm0
    \\ subsd %xmm3, %xmm2
    \\ .if \bp0
    \\ movsd %xmm0, 0 * 8(%rdx)
    \\ movsd %xmm2, 0 * 8(%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm0, 0 * 8(%rdx)
    \\ movlps %xmm2, 0 * 8(%r9)
    \\ .endif
    \\ addq $1 * 8, %rdx
    \\ addq $1 * 8, %r9
    \\ decq %rcx
    \\ jle 20017f
    \\ .balign 4
    \\ 20000:
    \\ testq $8, %r9
    \\ jne 20007f
    \\ movq %rcx, %rax
    \\ sarq $4, %rax
    \\ jle 20003f
    \\ movaps 0 * 8(%r9), %xmm1
    \\ movaps 2 * 8(%r9), %xmm3
    \\ movaps 4 * 8(%r9), %xmm9
    \\ movaps 6 * 8(%r9), %xmm11
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps 2 * 8(%rdx), %xmm2
    \\ movaps 4 * 8(%rdx), %xmm8
    \\ movaps 6 * 8(%rdx), %xmm10
    \\ decq %rax
    \\ jle 20002f
    \\ .balign 8
    \\ 20001:
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rdx)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%rdx)
    \\ .endif
    \\ movaps %xmm1, %xmm4
    \\ mulpd %xmm15, %xmm1
    \\ movaps %xmm3, %xmm6
    \\ mulpd %xmm15, %xmm3
    \\ movaps %xmm0, %xmm5
    \\ mulpd %xmm14, %xmm0
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ movaps 8 * 8(%r9), %xmm1
    \\ addpd %xmm3, %xmm2
    \\ movaps 10 * 8(%r9), %xmm3
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%r9)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%r9)
    \\ .endif
    \\ movaps %xmm0, 0 * 8(%rdx)
    \\ movaps 8 * 8(%rdx), %xmm0
    \\ movaps %xmm2, 2 * 8(%rdx)
    \\ movaps 10 * 8(%rdx), %xmm2
    \\ movaps %xmm4, 0 * 8(%r9)
    \\ movaps %xmm6, 2 * 8(%r9)
    \\ movaps %xmm9, %xmm4
    \\ mulpd %xmm15, %xmm9
    \\ movaps %xmm8, %xmm5
    \\ mulpd %xmm14, %xmm8
    \\ movaps %xmm11, %xmm6
    \\ mulpd %xmm15, %xmm11
    \\ movaps %xmm10, %xmm7
    \\ mulpd %xmm14, %xmm10
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm9, %xmm8
    \\ movaps 12 * 8(%r9), %xmm9
    \\ addpd %xmm11, %xmm10
    \\ movaps 14 * 8(%r9), %xmm11
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm8, 4 * 8(%rdx)
    \\ movaps 12 * 8(%rdx), %xmm8
    \\ movaps %xmm10,6 * 8(%rdx)
    \\ movaps 14 * 8(%rdx), %xmm10
    \\ movaps %xmm4, 4 * 8(%r9)
    \\ movaps %xmm6, 6 * 8(%r9)
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rdx)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%rdx)
    \\ .endif
    \\ movaps %xmm1, %xmm4
    \\ mulpd %xmm15, %xmm1
    \\ movaps %xmm3, %xmm6
    \\ mulpd %xmm15, %xmm3
    \\ movaps %xmm0, %xmm5
    \\ mulpd %xmm14, %xmm0
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ movaps 16 * 8(%r9), %xmm1
    \\ addpd %xmm3, %xmm2
    \\ movaps 18 * 8(%r9), %xmm3
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 8 * 8(%rdx)
    \\ movaps 16 * 8(%rdx), %xmm0
    \\ movaps %xmm2, 10 * 8(%rdx)
    \\ movaps 18 * 8(%rdx), %xmm2
    \\ movaps %xmm4, 8 * 8(%r9)
    \\ movaps %xmm6, 10 * 8(%r9)
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%r9)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%r9)
    \\ .endif
    \\ movaps %xmm9, %xmm4
    \\ mulpd %xmm15, %xmm9
    \\ movaps %xmm8, %xmm5
    \\ mulpd %xmm14, %xmm8
    \\ movaps %xmm11, %xmm6
    \\ mulpd %xmm15, %xmm11
    \\ movaps %xmm10, %xmm7
    \\ mulpd %xmm14, %xmm10
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm9, %xmm8
    \\ movaps 20 * 8(%r9), %xmm9
    \\ addpd %xmm11, %xmm10
    \\ movaps 22 * 8(%r9), %xmm11
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm8, 12 * 8(%rdx)
    \\ movaps 20 * 8(%rdx), %xmm8
    \\ movaps %xmm10, 14 * 8(%rdx)
    \\ movaps 22 * 8(%rdx), %xmm10
    \\ movaps %xmm4, 12 * 8(%r9)
    \\ movaps %xmm6, 14 * 8(%r9)
    \\ addq $16 * 8, %rdx
    \\ addq $16 * 8, %r9
    \\ decq %rax
    \\ jg 20001b
    \\ .balign 8
    \\ 20002:
    \\ movaps %xmm1, %xmm4
    \\ mulpd %xmm15, %xmm1
    \\ movaps %xmm3, %xmm6
    \\ mulpd %xmm15, %xmm3
    \\ movaps %xmm0, %xmm5
    \\ mulpd %xmm14, %xmm0
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ movaps 8 * 8(%r9), %xmm1
    \\ addpd %xmm3, %xmm2
    \\ movaps 10 * 8(%r9), %xmm3
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 0 * 8(%rdx)
    \\ movaps 8 * 8(%rdx), %xmm0
    \\ movaps %xmm2, 2 * 8(%rdx)
    \\ movaps 10 * 8(%rdx), %xmm2
    \\ movaps %xmm4, 0 * 8(%r9)
    \\ movaps %xmm6, 2 * 8(%r9)
    \\ movaps %xmm9, %xmm4
    \\ mulpd %xmm15, %xmm9
    \\ movaps %xmm8, %xmm5
    \\ mulpd %xmm14, %xmm8
    \\ movaps %xmm11, %xmm6
    \\ mulpd %xmm15, %xmm11
    \\ movaps %xmm10, %xmm7
    \\ mulpd %xmm14, %xmm10
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm9, %xmm8
    \\ movaps 12 * 8(%r9), %xmm9
    \\ addpd %xmm11, %xmm10
    \\ movaps 14 * 8(%r9), %xmm11
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm8, 4 * 8(%rdx)
    \\ movaps 12 * 8(%rdx), %xmm8
    \\ movaps %xmm10,6 * 8(%rdx)
    \\ movaps 14 * 8(%rdx), %xmm10
    \\ movaps %xmm4, 4 * 8(%r9)
    \\ movaps %xmm6, 6 * 8(%r9)
    \\ movaps %xmm1, %xmm4
    \\ mulpd %xmm15, %xmm1
    \\ movaps %xmm3, %xmm6
    \\ mulpd %xmm15, %xmm3
    \\ movaps %xmm0, %xmm5
    \\ mulpd %xmm14, %xmm0
    \\ movaps %xmm2, %xmm7
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm1, %xmm0
    \\ addpd %xmm3, %xmm2
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 8 * 8(%rdx)
    \\ movaps %xmm2, 10 * 8(%rdx)
    \\ movaps %xmm4, 8 * 8(%r9)
    \\ movaps %xmm6, 10 * 8(%r9)
    \\ movaps %xmm9, %xmm4
    \\ mulpd %xmm15, %xmm9
    \\ movaps %xmm8, %xmm5
    \\ mulpd %xmm14, %xmm8
    \\ movaps %xmm11, %xmm6
    \\ mulpd %xmm15, %xmm11
    \\ movaps %xmm10, %xmm7
    \\ mulpd %xmm14, %xmm10
    \\ mulpd %xmm14, %xmm4
    \\ mulpd %xmm15, %xmm5
    \\ mulpd %xmm14, %xmm6
    \\ mulpd %xmm15, %xmm7
    \\ addpd %xmm9, %xmm8
    \\ addpd %xmm11, %xmm10
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm8, 12 * 8(%rdx)
    \\ movaps %xmm10, 14 * 8(%rdx)
    \\ movaps %xmm4, 12 * 8(%r9)
    \\ movaps %xmm6, 14 * 8(%r9)
    \\ addq $16 * 8, %rdx
    \\ addq $16 * 8, %r9
    \\
;
