//! Section 4 of 4 of the text of `csrot`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zrot_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20016:
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
    \\ movsd 16 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 16 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 18 * 4(%rcx), %xmm1
    \\ addps %xmm3, %xmm2
    \\ .if \bp0
    \\ movsd 20 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 20 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 22 * 4(%rcx), %xmm3
    \\ subps %xmm5, %xmm4
    \\ subps %xmm7, %xmm6
    \\ movlps %xmm0, 0 * 4(%rsi)
    \\ movhps %xmm0, 2 * 4(%rsi)
    \\ .if \bp0
    \\ movsd 16 * 4(%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 16 * 4(%rsi), %xmm0
    \\ .endif
    \\ movhps 18 * 4(%rsi), %xmm0
    \\ movlps %xmm2, 4 * 4(%rsi)
    \\ movhps %xmm2, 6 * 4(%rsi)
    \\ .if \bp0
    \\ movsd 20 * 4(%rsi), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps 20 * 4(%rsi), %xmm2
    \\ .endif
    \\ movhps 22 * 4(%rsi), %xmm2
    \\ .if \bp0
    \\ movsd %xmm4, 0 * 4(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm4, 0 * 4(%rcx)
    \\ .endif
    \\ movhps %xmm4, 2 * 4(%rcx)
    \\ .if \bp0
    \\ movsd %xmm6, 4 * 4(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm6, 4 * 4(%rcx)
    \\ .endif
    \\ movhps %xmm6, 6 * 4(%rcx)
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
    \\ movsd 24 * 4(%rcx), %xmm9
    \\ .endif
    \\ .if \bp1
    \\ movlps 24 * 4(%rcx), %xmm9
    \\ .endif
    \\ movhps 26 * 4(%rcx), %xmm9
    \\ addps %xmm11, %xmm10
    \\ .if \bp0
    \\ movsd 28 * 4(%rcx), %xmm11
    \\ .endif
    \\ .if \bp1
    \\ movlps 28 * 4(%rcx), %xmm11
    \\ .endif
    \\ movhps 30 * 4(%rcx), %xmm11
    \\ subps %xmm5, %xmm4
    \\ subps %xmm7, %xmm6
    \\ movlps %xmm8, 8 * 4(%rsi)
    \\ movhps %xmm8, 10 * 4(%rsi)
    \\ .if \bp0
    \\ movsd 24 * 4(%rsi), %xmm8
    \\ .endif
    \\ .if \bp1
    \\ movlps 24 * 4(%rsi), %xmm8
    \\ .endif
    \\ movhps 26 * 4(%rsi), %xmm8
    \\ movlps %xmm10, 12 * 4(%rsi)
    \\ movhps %xmm10, 14 * 4(%rsi)
    \\ .if \bp0
    \\ movsd 28 * 4(%rsi), %xmm10
    \\ .endif
    \\ .if \bp1
    \\ movlps 28 * 4(%rsi), %xmm10
    \\ .endif
    \\ movhps 30 * 4(%rsi), %xmm10
    \\ movlps %xmm4, 8 * 4(%rcx)
    \\ movhps %xmm4, 10 * 4(%rcx)
    \\ movlps %xmm6, 12 * 4(%rcx)
    \\ movhps %xmm6, 14 * 4(%rcx)
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
    \\ movlps %xmm0, 16 * 4(%rsi)
    \\ movhps %xmm0, 18 * 4(%rsi)
    \\ movlps %xmm2, 20 * 4(%rsi)
    \\ movhps %xmm2, 22 * 4(%rsi)
    \\ movlps %xmm4, 16 * 4(%rcx)
    \\ movhps %xmm4, 18 * 4(%rcx)
    \\ movlps %xmm6, 20 * 4(%rcx)
    \\ movhps %xmm6, 22 * 4(%rcx)
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
    \\ movlps %xmm8, 24 * 4(%rsi)
    \\ movhps %xmm8, 26 * 4(%rsi)
    \\ movlps %xmm10, 28 * 4(%rsi)
    \\ movhps %xmm10, 30 * 4(%rsi)
    \\ movlps %xmm4, 24 * 4(%rcx)
    \\ movhps %xmm4, 26 * 4(%rcx)
    \\ movlps %xmm6, 28 * 4(%rcx)
    \\ movhps %xmm6, 30 * 4(%rcx)
    \\ addq $32 * 4, %rsi
    \\ addq $32 * 4, %rcx
    \\ .balign 8
    \\ 20017:
    \\ testq $15, %rdi
    \\ jle 20025f
    \\ testq $8, %rdi
    \\ jle 20018f
    \\ .if \bp0
    \\ movsd 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd 0 * 4(%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rsi), %xmm0
    \\ .endif
    \\ movhps 2 * 4(%rsi), %xmm0
    \\ .if \bp0
    \\ movsd 4 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 6 * 4(%rcx), %xmm3
    \\ .if \bp0
    \\ movsd 4 * 4(%rsi), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rsi), %xmm2
    \\ .endif
    \\ movhps 6 * 4(%rsi), %xmm2
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
    \\ movlps %xmm0, 0 * 4(%rsi)
    \\ movhps %xmm0, 2 * 4(%rsi)
    \\ movlps %xmm2, 4 * 4(%rsi)
    \\ movhps %xmm2, 6 * 4(%rsi)
    \\ movlps %xmm4, 0 * 4(%rcx)
    \\ movhps %xmm4, 2 * 4(%rcx)
    \\ movlps %xmm6, 4 * 4(%rcx)
    \\ movhps %xmm6, 6 * 4(%rcx)
    \\ .if \bp0
    \\ movsd 8 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 10 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd 8 * 4(%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rsi), %xmm0
    \\ .endif
    \\ movhps 10 * 4(%rsi), %xmm0
    \\ .if \bp0
    \\ movsd 12 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 14 * 4(%rcx), %xmm3
    \\ .if \bp0
    \\ movsd 12 * 4(%rsi), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%rsi), %xmm2
    \\ .endif
    \\ movhps 14 * 4(%rsi), %xmm2
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
    \\ movlps %xmm0, 8 * 4(%rsi)
    \\ movhps %xmm0, 10 * 4(%rsi)
    \\ movlps %xmm2, 12 * 4(%rsi)
    \\ movhps %xmm2, 14 * 4(%rsi)
    \\ movlps %xmm4, 8 * 4(%rcx)
    \\ movhps %xmm4, 10 * 4(%rcx)
    \\ movlps %xmm6, 12 * 4(%rcx)
    \\ movhps %xmm6, 14 * 4(%rcx)
    \\ addq $16 * 4, %rsi
    \\ addq $16 * 4, %rcx
    \\ .balign 8
    \\ 20018:
    \\ testq $4, %rdi
    \\ jle 20019f
    \\ .if \bp0
    \\ movsd 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd 0 * 4(%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rsi), %xmm0
    \\ .endif
    \\ movhps 2 * 4(%rsi), %xmm0
    \\ .if \bp0
    \\ movsd 4 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 6 * 4(%rcx), %xmm3
    \\ .if \bp0
    \\ movsd 4 * 4(%rsi), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rsi), %xmm2
    \\ .endif
    \\ movhps 6 * 4(%rsi), %xmm2
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
    \\ movlps %xmm0, 0 * 4(%rsi)
    \\ movhps %xmm0, 2 * 4(%rsi)
    \\ movlps %xmm2, 4 * 4(%rsi)
    \\ movhps %xmm2, 6 * 4(%rsi)
    \\ movlps %xmm4, 0 * 4(%rcx)
    \\ movhps %xmm4, 2 * 4(%rcx)
    \\ movlps %xmm6, 4 * 4(%rcx)
    \\ movhps %xmm6, 6 * 4(%rcx)
    \\ addq $8 * 4, %rsi
    \\ addq $8 * 4, %rcx
    \\ .balign 8
    \\ 20019:
    \\ testq $2, %rdi
    \\ jle 20020f
    \\ .if \bp0
    \\ movsd 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd 0 * 4(%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rsi), %xmm0
    \\ .endif
    \\ movhps 2 * 4(%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movlps %xmm0, 0 * 4(%rsi)
    \\ movhps %xmm0, 2 * 4(%rsi)
    \\ movlps %xmm2, 0 * 4(%rcx)
    \\ movhps %xmm2, 2 * 4(%rcx)
    \\ addq $4 * 4, %rsi
    \\ addq $4 * 4, %rcx
    \\ .balign 8
    \\ 20020:
    \\ testq $1, %rdi
    \\ jle 20025f
    \\ .if \bp0
    \\ movsd 0 * 4(%rcx), %xmm1
    \\ movsd 0 * 4(%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rcx), %xmm1
    \\ movlps 0 * 4(%rsi), %xmm0
    \\ .endif
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movlps %xmm0, 0 * 4(%rsi)
    \\ movlps %xmm2, 0 * 4(%rcx)
    \\ jmp 20025f
    \\ .balign 8
    \\ .balign 8
    \\ 20021:
    \\ movq %rdi, %rax
    \\ cmpq $0, %rdx
    \\ je 20024f
    \\ cmpq $0, %r8
    \\ je 20024f
    \\ sarq $2, %rax
    \\ jle 20023f
    \\ .balign 8
    \\ 20022:
    \\ .if \bp0
    \\ movsd (%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rcx), %xmm1
    \\ .endif
    \\ movhps (%rcx, %r8), %xmm1
    \\ .if \bp0
    \\ movsd (%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rsi), %xmm0
    \\ .endif
    \\ movhps (%rsi, %rdx), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movlps %xmm0, (%rsi)
    \\ movhps %xmm0, (%rsi, %rdx)
    \\ movlps %xmm2, (%rcx)
    \\ movhps %xmm2, (%rcx, %r8)
    \\ leaq (%rsi, %rdx, 2), %rsi
    \\ leaq (%rcx, %r8, 2), %rcx
    \\ .if \bp0
    \\ movsd (%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rcx), %xmm1
    \\ .endif
    \\ movhps (%rcx, %r8), %xmm1
    \\ .if \bp0
    \\ movsd (%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rsi), %xmm0
    \\ .endif
    \\ movhps (%rsi, %rdx), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movlps %xmm0, (%rsi)
    \\ movhps %xmm0, (%rsi, %rdx)
    \\ movlps %xmm2, (%rcx)
    \\ movhps %xmm2, (%rcx, %r8)
    \\ leaq (%rsi, %rdx, 2), %rsi
    \\ leaq (%rcx, %r8, 2), %rcx
    \\ decq %rax
    \\ jg 20022b
    \\ .balign 8
    \\ 20023:
    \\ movq %rdi, %rax
    \\ andq $3, %rax
    \\ jle 20025f
    \\ .balign 8
    \\ 20024:
    \\ .if \bp0
    \\ movsd (%rcx), %xmm1
    \\ movsd (%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rcx), %xmm1
    \\ movlps (%rsi), %xmm0
    \\ .endif
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movlps %xmm0, (%rsi)
    \\ movlps %xmm2, (%rcx)
    \\ addq %rdx, %rsi
    \\ addq %r8, %rcx
    \\ decq %rax
    \\ jg 20024b
    \\ .balign 8
    \\ 20025:
    \\ ret
    \\ .endm
    \\
;
