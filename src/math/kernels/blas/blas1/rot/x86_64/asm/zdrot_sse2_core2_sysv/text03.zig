//! Section 3 of 3 of the text of `zdrot`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zrot_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20017:
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
    \\ jg 20017b
    \\ .balign 8
    \\ 20018:
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
    \\ .balign 8
    \\ 20019:
    \\ testq $4, %rdi
    \\ jle 20020f
    \\ movaps 0 * 8(%rcx), %xmm1
    \\ movaps 0 * 8(%rsi), %xmm0
    \\ movaps 2 * 8(%rcx), %xmm3
    \\ movaps 2 * 8(%rsi), %xmm2
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
    \\ movaps %xmm0, 0 * 8(%rsi)
    \\ movaps %xmm2, 2 * 8(%rsi)
    \\ movaps %xmm4, 0 * 8(%rcx)
    \\ movaps %xmm6, 2 * 8(%rcx)
    \\ movaps 4 * 8(%rcx), %xmm1
    \\ movaps 4 * 8(%rsi), %xmm0
    \\ movaps 6 * 8(%rcx), %xmm3
    \\ movaps 6 * 8(%rsi), %xmm2
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
    \\ movaps %xmm0, 4 * 8(%rsi)
    \\ movaps %xmm2, 6 * 8(%rsi)
    \\ movaps %xmm4, 4 * 8(%rcx)
    \\ movaps %xmm6, 6 * 8(%rcx)
    \\ addq $8 * 8, %rsi
    \\ addq $8 * 8, %rcx
    \\ .balign 8
    \\ 20020:
    \\ testq $2, %rdi
    \\ jle 20021f
    \\ movaps 0 * 8(%rcx), %xmm1
    \\ movaps 0 * 8(%rsi), %xmm0
    \\ movaps 2 * 8(%rcx), %xmm3
    \\ movaps 2 * 8(%rsi), %xmm2
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
    \\ movaps %xmm0, 0 * 8(%rsi)
    \\ movaps %xmm2, 2 * 8(%rsi)
    \\ movaps %xmm4, 0 * 8(%rcx)
    \\ movaps %xmm6, 2 * 8(%rcx)
    \\ addq $4 * 8, %rsi
    \\ addq $4 * 8, %rcx
    \\ .balign 8
    \\ 20021:
    \\ testq $1, %rdi
    \\ jle 20022f
    \\ movaps 0 * 8(%rcx), %xmm1
    \\ movaps 0 * 8(%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm1, %xmm0
    \\ subpd %xmm3, %xmm2
    \\ movaps %xmm0, 0 * 8(%rsi)
    \\ movaps %xmm2, 0 * 8(%rcx)
    \\ addq $2 * 8, %rcx
    \\ addq $2 * 8, %rsi
    \\ .balign 8
    \\ 20022:
    \\ movsd 0 * 8(%rcx), %xmm1
    \\ movsd 0 * 8(%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulsd %xmm14, %xmm0
    \\ mulsd %xmm15, %xmm1
    \\ mulsd %xmm14, %xmm2
    \\ mulsd %xmm15, %xmm3
    \\ addsd %xmm1, %xmm0
    \\ subsd %xmm3, %xmm2
    \\ movsd %xmm0, 0 * 8(%rsi)
    \\ movsd %xmm2, 0 * 8(%rcx)
    \\ jmp 20027f
    \\ .balign 8
    \\ 20023:
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20025f
    \\ .balign 8
    \\ 20024:
    \\ movsd 0 * 8(%rcx), %xmm1
    \\ movhps 1 * 8(%rcx), %xmm1
    \\ movsd 0 * 8(%rsi), %xmm0
    \\ movhps 1 * 8(%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm1, %xmm0
    \\ subpd %xmm3, %xmm2
    \\ movlps %xmm0, 0 * 8(%rsi)
    \\ movhps %xmm0, 1 * 8(%rsi)
    \\ movlps %xmm2, 0 * 8(%rcx)
    \\ movhps %xmm2, 1 * 8(%rcx)
    \\ addq %rdx, %rsi
    \\ addq %r8, %rcx
    \\ movsd 0 * 8(%rcx), %xmm1
    \\ movhps 1 * 8(%rcx), %xmm1
    \\ movsd 0 * 8(%rsi), %xmm0
    \\ movhps 1 * 8(%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm1, %xmm0
    \\ subpd %xmm3, %xmm2
    \\ movlps %xmm0, 0 * 8(%rsi)
    \\ movhps %xmm0, 1 * 8(%rsi)
    \\ movlps %xmm2, 0 * 8(%rcx)
    \\ movhps %xmm2, 1 * 8(%rcx)
    \\ addq %rdx, %rsi
    \\ addq %r8, %rcx
    \\ movsd 0 * 8(%rcx), %xmm1
    \\ movhps 1 * 8(%rcx), %xmm1
    \\ movsd 0 * 8(%rsi), %xmm0
    \\ movhps 1 * 8(%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm1, %xmm0
    \\ subpd %xmm3, %xmm2
    \\ movlps %xmm0, 0 * 8(%rsi)
    \\ movhps %xmm0, 1 * 8(%rsi)
    \\ movlps %xmm2, 0 * 8(%rcx)
    \\ movhps %xmm2, 1 * 8(%rcx)
    \\ addq %rdx, %rsi
    \\ addq %r8, %rcx
    \\ movsd 0 * 8(%rcx), %xmm1
    \\ movhps 1 * 8(%rcx), %xmm1
    \\ movsd 0 * 8(%rsi), %xmm0
    \\ movhps 1 * 8(%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm1, %xmm0
    \\ subpd %xmm3, %xmm2
    \\ movlps %xmm0, 0 * 8(%rsi)
    \\ movhps %xmm0, 1 * 8(%rsi)
    \\ movlps %xmm2, 0 * 8(%rcx)
    \\ movhps %xmm2, 1 * 8(%rcx)
    \\ addq %rdx, %rsi
    \\ addq %r8, %rcx
    \\ decq %rax
    \\ jg 20024b
    \\ .balign 8
    \\ 20025:
    \\ movq %rdi, %rax
    \\ andq $3, %rax
    \\ jle 20027f
    \\ .balign 8
    \\ 20026:
    \\ movsd 0 * 8(%rcx), %xmm1
    \\ movhps 1 * 8(%rcx), %xmm1
    \\ movsd 0 * 8(%rsi), %xmm0
    \\ movhps 1 * 8(%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm1, %xmm0
    \\ subpd %xmm3, %xmm2
    \\ movlps %xmm0, 0 * 8(%rsi)
    \\ movhps %xmm0, 1 * 8(%rsi)
    \\ movlps %xmm2, 0 * 8(%rcx)
    \\ movhps %xmm2, 1 * 8(%rcx)
    \\ addq %rdx, %rsi
    \\ addq %r8, %rcx
    \\ decq %rax
    \\ jg 20026b
    \\ .balign 8
    \\ 20027:
    \\ ret
    \\
;
