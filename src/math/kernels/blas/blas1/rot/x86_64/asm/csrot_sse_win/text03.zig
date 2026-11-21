//! Section 3 of 5 of the text of `csrot`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zrot_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20009:
    \\ movaps %xmm1, %xmm4
    \\ mulps %xmm15, %xmm1
    \\ movaps %xmm3, %xmm6
    \\ mulps %xmm15, %xmm3
    \\ movaps %xmm0, %xmm5
    \\ mulps %xmm14, %xmm0
    \\ movaps %xmm2, %xmm7
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm5
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm7
    \\ addps %xmm1, %xmm0
    \\ .if \bp0
    \\ movsd 16 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 16 * 4(%r9), %xmm1
    \\ .endif
    \\ movhps 18 * 4(%r9), %xmm1
    \\ addps %xmm3, %xmm2
    \\ .if \bp0
    \\ movsd 20 * 4(%r9), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 20 * 4(%r9), %xmm3
    \\ .endif
    \\ movhps 22 * 4(%r9), %xmm3
    \\ subps %xmm5, %xmm4
    \\ subps %xmm7, %xmm6
    \\ movaps %xmm0, 0 * 4(%rdx)
    \\ movaps 16 * 4(%rdx), %xmm0
    \\ movaps %xmm2, 4 * 4(%rdx)
    \\ movaps 20 * 4(%rdx), %xmm2
    \\ .if \bp0
    \\ movsd %xmm4, 0 * 4(%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm4, 0 * 4(%r9)
    \\ .endif
    \\ movhps %xmm4, 2 * 4(%r9)
    \\ .if \bp0
    \\ movsd %xmm6, 4 * 4(%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm6, 4 * 4(%r9)
    \\ .endif
    \\ movhps %xmm6, 6 * 4(%r9)
    \\ movaps %xmm9, %xmm4
    \\ mulps %xmm15, %xmm9
    \\ movaps %xmm8, %xmm5
    \\ mulps %xmm14, %xmm8
    \\ movaps %xmm11, %xmm6
    \\ mulps %xmm15, %xmm11
    \\ movaps %xmm10, %xmm7
    \\ mulps %xmm14, %xmm10
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm5
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm7
    \\ addps %xmm9, %xmm8
    \\ .if \bp0
    \\ movsd 24 * 4(%r9), %xmm9
    \\ .endif
    \\ .if \bp1
    \\ movlps 24 * 4(%r9), %xmm9
    \\ .endif
    \\ movhps 26 * 4(%r9), %xmm9
    \\ addps %xmm11, %xmm10
    \\ .if \bp0
    \\ movsd 28 * 4(%r9), %xmm11
    \\ .endif
    \\ .if \bp1
    \\ movlps 28 * 4(%r9), %xmm11
    \\ .endif
    \\ movhps 30 * 4(%r9), %xmm11
    \\ subps %xmm5, %xmm4
    \\ subps %xmm7, %xmm6
    \\ movaps %xmm8, 8 * 4(%rdx)
    \\ movaps 24 * 4(%rdx), %xmm8
    \\ movaps %xmm10,12 * 4(%rdx)
    \\ movaps 28 * 4(%rdx), %xmm10
    \\ movlps %xmm4, 8 * 4(%r9)
    \\ movhps %xmm4, 10 * 4(%r9)
    \\ movlps %xmm6, 12 * 4(%r9)
    \\ movhps %xmm6, 14 * 4(%r9)
    \\ movaps %xmm1, %xmm4
    \\ mulps %xmm15, %xmm1
    \\ movaps %xmm3, %xmm6
    \\ mulps %xmm15, %xmm3
    \\ movaps %xmm0, %xmm5
    \\ mulps %xmm14, %xmm0
    \\ movaps %xmm2, %xmm7
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm5
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm7
    \\ addps %xmm1, %xmm0
    \\ addps %xmm3, %xmm2
    \\ subps %xmm5, %xmm4
    \\ subps %xmm7, %xmm6
    \\ movaps %xmm0, 16 * 4(%rdx)
    \\ movaps %xmm2, 20 * 4(%rdx)
    \\ movlps %xmm4, 16 * 4(%r9)
    \\ movhps %xmm4, 18 * 4(%r9)
    \\ movlps %xmm6, 20 * 4(%r9)
    \\ movhps %xmm6, 22 * 4(%r9)
    \\ movaps %xmm9, %xmm4
    \\ mulps %xmm15, %xmm9
    \\ movaps %xmm8, %xmm5
    \\ mulps %xmm14, %xmm8
    \\ movaps %xmm11, %xmm6
    \\ mulps %xmm15, %xmm11
    \\ movaps %xmm10, %xmm7
    \\ mulps %xmm14, %xmm10
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm5
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm7
    \\ addps %xmm9, %xmm8
    \\ addps %xmm11, %xmm10
    \\ subps %xmm5, %xmm4
    \\ subps %xmm7, %xmm6
    \\ movaps %xmm8, 24 * 4(%rdx)
    \\ movaps %xmm10, 28 * 4(%rdx)
    \\ movlps %xmm4, 24 * 4(%r9)
    \\ movhps %xmm4, 26 * 4(%r9)
    \\ movlps %xmm6, 28 * 4(%r9)
    \\ movhps %xmm6, 30 * 4(%r9)
    \\ addq $32 * 4, %rdx
    \\ addq $32 * 4, %r9
    \\ .balign 8
    \\ 20010:
    \\ testq $15, %rcx
    \\ jle 20025f
    \\ testq $8, %rcx
    \\ jle 20011f
    \\ .if \bp0
    \\ movsd 0 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%r9), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%r9), %xmm1
    \\ movaps 0 * 4(%rdx), %xmm0
    \\ .if \bp0
    \\ movsd 4 * 4(%r9), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%r9), %xmm3
    \\ .endif
    \\ movhps 6 * 4(%r9), %xmm3
    \\ movaps 4 * 4(%rdx), %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm5
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm7
    \\ addps %xmm1, %xmm0
    \\ addps %xmm3, %xmm2
    \\ subps %xmm5, %xmm4
    \\ subps %xmm7, %xmm6
    \\ movaps %xmm0, 0 * 4(%rdx)
    \\ movaps %xmm2, 4 * 4(%rdx)
    \\ movlps %xmm4, 0 * 4(%r9)
    \\ movhps %xmm4, 2 * 4(%r9)
    \\ movlps %xmm6, 4 * 4(%r9)
    \\ movhps %xmm6, 6 * 4(%r9)
    \\ .if \bp0
    \\ movsd 8 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%r9), %xmm1
    \\ .endif
    \\ movhps 10 * 4(%r9), %xmm1
    \\ movaps 8 * 4(%rdx), %xmm0
    \\ .if \bp0
    \\ movsd 12 * 4(%r9), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%r9), %xmm3
    \\ .endif
    \\ movhps 14 * 4(%r9), %xmm3
    \\ movaps 12 * 4(%rdx), %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm5
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm7
    \\ addps %xmm1, %xmm0
    \\ addps %xmm3, %xmm2
    \\ subps %xmm5, %xmm4
    \\ subps %xmm7, %xmm6
    \\ movaps %xmm0, 8 * 4(%rdx)
    \\ movaps %xmm2, 12 * 4(%rdx)
    \\ movlps %xmm4, 8 * 4(%r9)
    \\ movhps %xmm4, 10 * 4(%r9)
    \\ movlps %xmm6, 12 * 4(%r9)
    \\ movhps %xmm6, 14 * 4(%r9)
    \\ addq $16 * 4, %rdx
    \\ addq $16 * 4, %r9
    \\ .balign 8
    \\ 20011:
    \\ testq $4, %rcx
    \\ jle 20012f
    \\ .if \bp0
    \\ movsd 0 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%r9), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%r9), %xmm1
    \\ movaps 0 * 4(%rdx), %xmm0
    \\ .if \bp0
    \\ movsd 4 * 4(%r9), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%r9), %xmm3
    \\ .endif
    \\ movhps 6 * 4(%r9), %xmm3
    \\ movaps 4 * 4(%rdx), %xmm2
    \\ movaps %xmm1, %xmm4
    \\ movaps %xmm0, %xmm5
    \\ movaps %xmm3, %xmm6
    \\ movaps %xmm2, %xmm7
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm5
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm7
    \\ addps %xmm1, %xmm0
    \\ addps %xmm3, %xmm2
    \\ subps %xmm5, %xmm4
    \\ subps %xmm7, %xmm6
    \\ movaps %xmm0, 0 * 4(%rdx)
    \\ movaps %xmm2, 4 * 4(%rdx)
    \\ movlps %xmm4, 0 * 4(%r9)
    \\ movhps %xmm4, 2 * 4(%r9)
    \\ movlps %xmm6, 4 * 4(%r9)
    \\ movhps %xmm6, 6 * 4(%r9)
    \\ addq $8 * 4, %rdx
    \\ addq $8 * 4, %r9
    \\ .balign 8
    \\ 20012:
    \\ testq $2, %rcx
    \\ jle 20013f
    \\ .if \bp0
    \\ movsd 0 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%r9), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%r9), %xmm1
    \\ movaps 0 * 4(%rdx), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movaps %xmm0, 0 * 4(%rdx)
    \\ movlps %xmm2, 0 * 4(%r9)
    \\ movhps %xmm2, 2 * 4(%r9)
    \\ addq $4 * 4, %rdx
    \\ addq $4 * 4, %r9
    \\ .balign 8
    \\ 20013:
    \\ testq $1, %rcx
    \\ jle 20025f
    \\ .if \bp0
    \\ movsd 0 * 4(%r9), %xmm1
    \\ movsd 0 * 4(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%r9), %xmm1
    \\ movlps 0 * 4(%rdx), %xmm0
    \\ .endif
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movlps %xmm0, 0 * 4(%rdx)
    \\ movlps %xmm2, 0 * 4(%r9)
    \\ jmp 20025f
    \\ .balign 8
    \\ 20014:
    \\ movq %rcx, %rax
    \\ sarq $4, %rax
    \\ jle 20017f
    \\ .if \bp0
    \\ movsd 0 * 4(%r9), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%r9), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%r9), %xmm1
    \\ .if \bp0
    \\ movsd 4 * 4(%r9), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%r9), %xmm3
    \\ .endif
    \\ movhps 6 * 4(%r9), %xmm3
    \\ .if \bp0
    \\ movsd 8 * 4(%r9), %xmm9
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%r9), %xmm9
    \\ .endif
    \\ movhps 10 * 4(%r9), %xmm9
    \\ .if \bp0
    \\ movsd 12 * 4(%r9), %xmm11
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%r9), %xmm11
    \\ .endif
    \\ movhps 14 * 4(%r9), %xmm11
    \\ .if \bp0
    \\ movsd 0 * 4(%rdx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rdx), %xmm0
    \\ .endif
    \\ movhps 2 * 4(%rdx), %xmm0
    \\ .if \bp0
    \\ movsd 4 * 4(%rdx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rdx), %xmm2
    \\ .endif
    \\ movhps 6 * 4(%rdx), %xmm2
    \\ .if \bp0
    \\ movsd 8 * 4(%rdx), %xmm8
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rdx), %xmm8
    \\ .endif
    \\ movhps 10 * 4(%rdx), %xmm8
    \\ .if \bp0
    \\ movsd 12 * 4(%rdx), %xmm10
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%rdx), %xmm10
    \\ .endif
    \\ movhps 14 * 4(%rdx), %xmm10
    \\ decq %rax
    \\ jle 20016f
    \\
;
