//! Section 1 of 3 of the text of `drot`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/rot_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1, bp2, bp3, bp4, bp5, bp6
    \\ .balign 512
    \\ leaq (, %rdx, 8), %rdx
    \\ leaq (, %r8, 8), %r8
    \\ pshufd $0x44, %xmm0, %xmm14
    \\ pshufd $0x44, %xmm1, %xmm15
    \\ cmpq $0, %rdi
    \\ jle 20017f
    \\ cmpq $8, %rdx
    \\ jne 20013f
    \\ cmpq $8, %r8
    \\ jne 20013f
    \\ testq $8, %rsi
    \\ je 20000f
    \\ .if \bp0
    \\ movsd 0 * 8(%rcx), %xmm1
    \\ movsd 0 * 8(%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 8(%rcx), %xmm1
    \\ movlps 0 * 8(%rsi), %xmm0
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
    \\ movsd %xmm0, 0 * 8(%rsi)
    \\ movsd %xmm2, 0 * 8(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm0, 0 * 8(%rsi)
    \\ movlps %xmm2, 0 * 8(%rcx)
    \\ .endif
    \\ addq $1 * 8, %rsi
    \\ addq $1 * 8, %rcx
    \\ decq %rdi
    \\ jle 20017f
    \\ .balign 4
    \\ 20000:
    \\ testq $8, %rcx
    \\ jne 20007f
    \\ movq %rdi, %rax
    \\ sarq $4, %rax
    \\ jle 20003f
    \\ movaps 0 * 8(%rcx), %xmm1
    \\ movaps 2 * 8(%rcx), %xmm3
    \\ movaps 4 * 8(%rcx), %xmm9
    \\ movaps 6 * 8(%rcx), %xmm11
    \\ movaps 0 * 8(%rsi), %xmm0
    \\ movaps 2 * 8(%rsi), %xmm2
    \\ movaps 4 * 8(%rsi), %xmm8
    \\ movaps 6 * 8(%rsi), %xmm10
    \\ decq %rax
    \\ jle 20002f
    \\ .balign 8
    \\ 20001:
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%rsi)
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
    \\ movaps 8 * 8(%rcx), %xmm1
    \\ addpd %xmm3, %xmm2
    \\ movaps 10 * 8(%rcx), %xmm3
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%rcx)
    \\ .endif
    \\ movaps %xmm0, 0 * 8(%rsi)
    \\ movaps 8 * 8(%rsi), %xmm0
    \\ movaps %xmm2, 2 * 8(%rsi)
    \\ movaps 10 * 8(%rsi), %xmm2
    \\ movaps %xmm4, 0 * 8(%rcx)
    \\ movaps %xmm6, 2 * 8(%rcx)
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
    \\ movaps 12 * 8(%rcx), %xmm9
    \\ addpd %xmm11, %xmm10
    \\ movaps 14 * 8(%rcx), %xmm11
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm8, 4 * 8(%rsi)
    \\ movaps 12 * 8(%rsi), %xmm8
    \\ movaps %xmm10,6 * 8(%rsi)
    \\ movaps 14 * 8(%rsi), %xmm10
    \\ movaps %xmm4, 4 * 8(%rcx)
    \\ movaps %xmm6, 6 * 8(%rcx)
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%rsi)
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
    \\ movaps 16 * 8(%rcx), %xmm1
    \\ addpd %xmm3, %xmm2
    \\ movaps 18 * 8(%rcx), %xmm3
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 8 * 8(%rsi)
    \\ movaps 16 * 8(%rsi), %xmm0
    \\ movaps %xmm2, 10 * 8(%rsi)
    \\ movaps 18 * 8(%rsi), %xmm2
    \\ movaps %xmm4, 8 * 8(%rcx)
    \\ movaps %xmm6, 10 * 8(%rcx)
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%rcx)
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
    \\ movaps 20 * 8(%rcx), %xmm9
    \\ addpd %xmm11, %xmm10
    \\ movaps 22 * 8(%rcx), %xmm11
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm8, 12 * 8(%rsi)
    \\ movaps 20 * 8(%rsi), %xmm8
    \\ movaps %xmm10, 14 * 8(%rsi)
    \\ movaps 22 * 8(%rsi), %xmm10
    \\ movaps %xmm4, 12 * 8(%rcx)
    \\ movaps %xmm6, 14 * 8(%rcx)
    \\ addq $16 * 8, %rsi
    \\ addq $16 * 8, %rcx
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
    \\ movaps 8 * 8(%rcx), %xmm1
    \\ addpd %xmm3, %xmm2
    \\ movaps 10 * 8(%rcx), %xmm3
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm0, 0 * 8(%rsi)
    \\ movaps 8 * 8(%rsi), %xmm0
    \\ movaps %xmm2, 2 * 8(%rsi)
    \\ movaps 10 * 8(%rsi), %xmm2
    \\ movaps %xmm4, 0 * 8(%rcx)
    \\ movaps %xmm6, 2 * 8(%rcx)
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
    \\ movaps 12 * 8(%rcx), %xmm9
    \\ addpd %xmm11, %xmm10
    \\ movaps 14 * 8(%rcx), %xmm11
    \\ subpd %xmm5, %xmm4
    \\ subpd %xmm7, %xmm6
    \\ movaps %xmm8, 4 * 8(%rsi)
    \\ movaps 12 * 8(%rsi), %xmm8
    \\ movaps %xmm10,6 * 8(%rsi)
    \\ movaps 14 * 8(%rsi), %xmm10
    \\ movaps %xmm4, 4 * 8(%rcx)
    \\ movaps %xmm6, 6 * 8(%rcx)
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
    \\ movaps %xmm0, 8 * 8(%rsi)
    \\ movaps %xmm2, 10 * 8(%rsi)
    \\ movaps %xmm4, 8 * 8(%rcx)
    \\ movaps %xmm6, 10 * 8(%rcx)
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
    \\ movaps %xmm8, 12 * 8(%rsi)
    \\ movaps %xmm10, 14 * 8(%rsi)
    \\ movaps %xmm4, 12 * 8(%rcx)
    \\ movaps %xmm6, 14 * 8(%rcx)
    \\ addq $16 * 8, %rsi
    \\ addq $16 * 8, %rcx
    \\
;
