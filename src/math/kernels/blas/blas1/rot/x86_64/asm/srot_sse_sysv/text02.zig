//! Section 2 of 2 of the text of `srot`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/rot_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20010:
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp0
    \\ movsd 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%rsi)
    \\ movlps 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd 4 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 6 * 4(%rcx), %xmm3
    \\ movaps 0 * 4(%rsi), %xmm0
    \\ movaps 4 * 4(%rsi), %xmm2
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
    \\ movaps %xmm0, 0 * 4(%rsi)
    \\ movaps %xmm2, 4 * 4(%rsi)
    \\ movlps %xmm4, 0 * 4(%rcx)
    \\ movhps %xmm4, 2 * 4(%rcx)
    \\ movlps %xmm6, 4 * 4(%rcx)
    \\ movhps %xmm6, 6 * 4(%rcx)
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp0
    \\ movsd 8 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%rcx)
    \\ movlps 8 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 10 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd 12 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 14 * 4(%rcx), %xmm3
    \\ movaps 8 * 4(%rsi), %xmm0
    \\ movaps 12 * 4(%rsi), %xmm2
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
    \\ movaps %xmm0, 8 * 4(%rsi)
    \\ movaps %xmm2, 12 * 4(%rsi)
    \\ movlps %xmm4, 8 * 4(%rcx)
    \\ movhps %xmm4, 10 * 4(%rcx)
    \\ movlps %xmm6, 12 * 4(%rcx)
    \\ movhps %xmm6, 14 * 4(%rcx)
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp0
    \\ movsd 16 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%rsi)
    \\ movlps 16 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 18 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd 20 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 20 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 22 * 4(%rcx), %xmm3
    \\ movaps 16 * 4(%rsi), %xmm0
    \\ movaps 20 * 4(%rsi), %xmm2
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
    \\ movaps %xmm0, 16 * 4(%rsi)
    \\ movaps %xmm2, 20 * 4(%rsi)
    \\ movlps %xmm4, 16 * 4(%rcx)
    \\ movhps %xmm4, 18 * 4(%rcx)
    \\ movlps %xmm6, 20 * 4(%rcx)
    \\ movhps %xmm6, 22 * 4(%rcx)
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp0
    \\ movsd 24 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%rcx)
    \\ movlps 24 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 26 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd 28 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 28 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 30 * 4(%rcx), %xmm3
    \\ movaps 24 * 4(%rsi), %xmm0
    \\ movaps 28 * 4(%rsi), %xmm2
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
    \\ movaps %xmm0, 24 * 4(%rsi)
    \\ movaps %xmm2, 28 * 4(%rsi)
    \\ movlps %xmm4, 24 * 4(%rcx)
    \\ movhps %xmm4, 26 * 4(%rcx)
    \\ movlps %xmm6, 28 * 4(%rcx)
    \\ movhps %xmm6, 30 * 4(%rcx)
    \\ addq $32 * 4, %rsi
    \\ addq $32 * 4, %rcx
    \\ decq %rax
    \\ jg 20010b
    \\ .balign 8
    \\ 20011:
    \\ testq $31, %rdi
    \\ jle 20020f
    \\ testq $16, %rdi
    \\ jle 20012f
    \\ .if \bp0
    \\ movsd 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd 4 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 6 * 4(%rcx), %xmm3
    \\ movaps 0 * 4(%rsi), %xmm0
    \\ movaps 4 * 4(%rsi), %xmm2
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
    \\ movaps %xmm0, 0 * 4(%rsi)
    \\ movaps %xmm2, 4 * 4(%rsi)
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
    \\ movsd 12 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 14 * 4(%rcx), %xmm3
    \\ movaps 8 * 4(%rsi), %xmm0
    \\ movaps 12 * 4(%rsi), %xmm2
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
    \\ movaps %xmm0, 8 * 4(%rsi)
    \\ movaps %xmm2, 12 * 4(%rsi)
    \\ movlps %xmm4, 8 * 4(%rcx)
    \\ movhps %xmm4, 10 * 4(%rcx)
    \\ movlps %xmm6, 12 * 4(%rcx)
    \\ movhps %xmm6, 14 * 4(%rcx)
    \\ addq $16 * 4, %rsi
    \\ addq $16 * 4, %rcx
    \\ .balign 8
    \\ 20012:
    \\ testq $8, %rdi
    \\ jle 20013f
    \\ .if \bp0
    \\ movsd 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd 4 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 6 * 4(%rcx), %xmm3
    \\ movaps 0 * 4(%rsi), %xmm0
    \\ movaps 4 * 4(%rsi), %xmm2
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
    \\ movaps %xmm0, 0 * 4(%rsi)
    \\ movaps %xmm2, 4 * 4(%rsi)
    \\ movlps %xmm4, 0 * 4(%rcx)
    \\ movhps %xmm4, 2 * 4(%rcx)
    \\ movlps %xmm6, 4 * 4(%rcx)
    \\ movhps %xmm6, 6 * 4(%rcx)
    \\ addq $8 * 4, %rsi
    \\ addq $8 * 4, %rcx
    \\ .balign 8
    \\ 20013:
    \\ testq $4, %rdi
    \\ jle 20014f
    \\ .if \bp0
    \\ movsd 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 2 * 4(%rcx), %xmm1
    \\ movaps 0 * 4(%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps %xmm1, %xmm0
    \\ subps %xmm3, %xmm2
    \\ movaps %xmm0, 0 * 4(%rsi)
    \\ movlps %xmm2, 0 * 4(%rcx)
    \\ movhps %xmm2, 2 * 4(%rcx)
    \\ addq $4 * 4, %rsi
    \\ addq $4 * 4, %rcx
    \\ .balign 8
    \\ 20014:
    \\ testq $2, %rdi
    \\ jle 20015f
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
    \\ addq $2 * 4, %rsi
    \\ addq $2 * 4, %rcx
    \\ .balign 8
    \\ 20015:
    \\ testq $1, %rdi
    \\ jle 20020f
    \\ movss 0 * 4(%rcx), %xmm1
    \\ movss 0 * 4(%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulss %xmm14, %xmm0
    \\ mulss %xmm15, %xmm1
    \\ mulss %xmm14, %xmm2
    \\ mulss %xmm15, %xmm3
    \\ addss %xmm1, %xmm0
    \\ subss %xmm3, %xmm2
    \\ movss %xmm0, 0 * 4(%rsi)
    \\ movss %xmm2, 0 * 4(%rcx)
    \\ jmp 20020f
    \\ .balign 8
    \\ 20016:
    \\ movq %rdi, %rax
    \\ sarq $2, %rax
    \\ jle 20018f
    \\ .balign 8
    \\ 20017:
    \\ movss (%rcx), %xmm1
    \\ movss (%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulss %xmm14, %xmm0
    \\ mulss %xmm15, %xmm1
    \\ mulss %xmm14, %xmm2
    \\ mulss %xmm15, %xmm3
    \\ addss %xmm1, %xmm0
    \\ subss %xmm3, %xmm2
    \\ movss %xmm0, (%rsi)
    \\ movss %xmm2, (%rcx)
    \\ addq %rdx, %rsi
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm1
    \\ movss (%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulss %xmm14, %xmm0
    \\ mulss %xmm15, %xmm1
    \\ mulss %xmm14, %xmm2
    \\ mulss %xmm15, %xmm3
    \\ addss %xmm1, %xmm0
    \\ subss %xmm3, %xmm2
    \\ movss %xmm0, (%rsi)
    \\ movss %xmm2, (%rcx)
    \\ addq %rdx, %rsi
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm1
    \\ movss (%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulss %xmm14, %xmm0
    \\ mulss %xmm15, %xmm1
    \\ mulss %xmm14, %xmm2
    \\ mulss %xmm15, %xmm3
    \\ addss %xmm1, %xmm0
    \\ subss %xmm3, %xmm2
    \\ movss %xmm0, (%rsi)
    \\ movss %xmm2, (%rcx)
    \\ addq %rdx, %rsi
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm1
    \\ movss (%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulss %xmm14, %xmm0
    \\ mulss %xmm15, %xmm1
    \\ mulss %xmm14, %xmm2
    \\ mulss %xmm15, %xmm3
    \\ addss %xmm1, %xmm0
    \\ subss %xmm3, %xmm2
    \\ movss %xmm0, (%rsi)
    \\ movss %xmm2, (%rcx)
    \\ addq %rdx, %rsi
    \\ addq %r8, %rcx
    \\ decq %rax
    \\ jg 20017b
    \\ .balign 8
    \\ 20018:
    \\ movq %rdi, %rax
    \\ andq $3, %rax
    \\ jle 20020f
    \\ .balign 8
    \\ 20019:
    \\ movss (%rcx), %xmm1
    \\ movss (%rsi), %xmm0
    \\ movaps %xmm1, %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulss %xmm14, %xmm0
    \\ mulss %xmm15, %xmm1
    \\ mulss %xmm14, %xmm2
    \\ mulss %xmm15, %xmm3
    \\ addss %xmm1, %xmm0
    \\ subss %xmm3, %xmm2
    \\ movss %xmm0, (%rsi)
    \\ movss %xmm2, (%rcx)
    \\ addq %rdx, %rsi
    \\ addq %r8, %rcx
    \\ decq %rax
    \\ jg 20019b
    \\ .balign 8
    \\ 20020:
    \\ ret
    \\ .endm
    \\
;
