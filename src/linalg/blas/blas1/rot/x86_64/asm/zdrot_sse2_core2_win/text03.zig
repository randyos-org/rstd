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
    \\ .balign 8
    \\ 20019:
    \\ testq $4, %rcx
    \\ jle 20020f
    \\ movaps 0 * 8(%r9), %xmm1
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps 2 * 8(%r9), %xmm3
    \\ movaps 2 * 8(%rdx), %xmm2
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
    \\ movaps %xmm4, 0 * 8(%r9)
    \\ movaps %xmm6, 2 * 8(%r9)
    \\ movaps 4 * 8(%r9), %xmm1
    \\ movaps 4 * 8(%rdx), %xmm0
    \\ movaps 6 * 8(%r9), %xmm3
    \\ movaps 6 * 8(%rdx), %xmm2
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
    \\ movaps %xmm0, 4 * 8(%rdx)
    \\ movaps %xmm2, 6 * 8(%rdx)
    \\ movaps %xmm4, 4 * 8(%r9)
    \\ movaps %xmm6, 6 * 8(%r9)
    \\ addq $8 * 8, %rdx
    \\ addq $8 * 8, %r9
    \\ .balign 8
    \\ 20020:
    \\ testq $2, %rcx
    \\ jle 20021f
    \\ movaps 0 * 8(%r9), %xmm1
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps 2 * 8(%r9), %xmm3
    \\ movaps 2 * 8(%rdx), %xmm2
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
    \\ movaps %xmm4, 0 * 8(%r9)
    \\ movaps %xmm6, 2 * 8(%r9)
    \\ addq $4 * 8, %rdx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20021:
    \\ testq $1, %rcx
    \\ jle 20022f
    \\ movaps 0 * 8(%r9), %xmm1
    \\ movaps 0 * 8(%rdx), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulpd %xmm14, %xmm0
    \\ mulpd %xmm15, %xmm1
    \\ mulpd %xmm14, %xmm2
    \\ mulpd %xmm15, %xmm3
    \\ addpd %xmm1, %xmm0
    \\ subpd %xmm3, %xmm2
    \\ movaps %xmm0, 0 * 8(%rdx)
    \\ movaps %xmm2, 0 * 8(%r9)
    \\ addq $2 * 8, %r9
    \\ addq $2 * 8, %rdx
    \\ .balign 8
    \\ 20022:
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
    \\ jmp 20027f
    \\ .balign 8
    \\ 20023:
    \\ movq %rcx, %rax
    \\ sarq $2, %rax
    \\ jle 20025f
    \\ .balign 8
    \\ 20024:
    \\ movsd 0 * 8(%r9), %xmm1
    \\ movhps 1 * 8(%r9), %xmm1
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ movhps 1 * 8(%rdx), %xmm0
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
    \\ movlps %xmm2, 0 * 8(%r9)
    \\ movhps %xmm2, 1 * 8(%r9)
    \\ addq %r8, %rdx
    \\ addq %r10, %r9
    \\ movsd 0 * 8(%r9), %xmm1
    \\ movhps 1 * 8(%r9), %xmm1
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ movhps 1 * 8(%rdx), %xmm0
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
    \\ movlps %xmm2, 0 * 8(%r9)
    \\ movhps %xmm2, 1 * 8(%r9)
    \\ addq %r8, %rdx
    \\ addq %r10, %r9
    \\ movsd 0 * 8(%r9), %xmm1
    \\ movhps 1 * 8(%r9), %xmm1
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ movhps 1 * 8(%rdx), %xmm0
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
    \\ movlps %xmm2, 0 * 8(%r9)
    \\ movhps %xmm2, 1 * 8(%r9)
    \\ addq %r8, %rdx
    \\ addq %r10, %r9
    \\ movsd 0 * 8(%r9), %xmm1
    \\ movhps 1 * 8(%r9), %xmm1
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ movhps 1 * 8(%rdx), %xmm0
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
    \\ movlps %xmm2, 0 * 8(%r9)
    \\ movhps %xmm2, 1 * 8(%r9)
    \\ addq %r8, %rdx
    \\ addq %r10, %r9
    \\ decq %rax
    \\ jg 20024b
    \\ .balign 8
    \\ 20025:
    \\ movq %rcx, %rax
    \\ andq $3, %rax
    \\ jle 20027f
    \\ .balign 8
    \\ 20026:
    \\ movsd 0 * 8(%r9), %xmm1
    \\ movhps 1 * 8(%r9), %xmm1
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ movhps 1 * 8(%rdx), %xmm0
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
    \\ movlps %xmm2, 0 * 8(%r9)
    \\ movhps %xmm2, 1 * 8(%r9)
    \\ addq %r8, %rdx
    \\ addq %r10, %r9
    \\ decq %rax
    \\ jg 20026b
    \\ .balign 8
    \\ 20027:
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
    \\
;
