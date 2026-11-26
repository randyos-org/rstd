//! Section 2 of 2 of the text of `sscal`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/scal_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 16
    \\ 20015:
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp7
    \\ mulps %xmm0, %xmm1
    \\ .endif
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 0) - 128(%rcx)
    \\ .endif
    \\ .if \bp8
    \\ mulps %xmm0, %xmm1
    \\ .endif
    \\ movaps %xmm1, -32 * 4(%rcx)
    \\ .if \bp5
    \\ movaps 0 * 4(%rcx), %xmm1
    \\ mulps %xmm0, %xmm2
    \\ .endif
    \\ movaps %xmm2, -28 * 4(%rcx)
    \\ .if \bp5
    \\ movaps 4 * 4(%rcx), %xmm2
    \\ mulps %xmm0, %xmm3
    \\ .endif
    \\ movaps %xmm3, -24 * 4(%rcx)
    \\ .if \bp5
    \\ movaps 8 * 4(%rcx), %xmm3
    \\ mulps %xmm0, %xmm4
    \\ .endif
    \\ movaps %xmm4, -20 * 4(%rcx)
    \\ .if \bp5
    \\ movaps 12 * 4(%rcx), %xmm4
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
    \\ .if \bp7
    \\ mulps %xmm0, %xmm5
    \\ .endif
    \\ .if \bp6
    \\ movaps %xmm0, %xmm1
    \\ mulps 0 * 4(%rcx), %xmm1
    \\ movaps %xmm0, %xmm2
    \\ mulps 4 * 4(%rcx), %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps 8 * 4(%rcx), %xmm3
    \\ movaps %xmm0, %xmm4
    \\ mulps 12 * 4(%rcx), %xmm4
    \\ .endif
    \\ .if \bp2
    \\ prefetchw ((128 * 5) + 64) - 128(%rcx)
    \\ .endif
    \\ .if \bp8
    \\ mulps %xmm0, %xmm5
    \\ .endif
    \\ movaps %xmm5, -16 * 4(%rcx)
    \\ .if \bp5
    \\ movaps 16 * 4(%rcx), %xmm5
    \\ mulps %xmm0, %xmm6
    \\ .endif
    \\ movaps %xmm6, -12 * 4(%rcx)
    \\ .if \bp5
    \\ movaps 20 * 4(%rcx), %xmm6
    \\ mulps %xmm0, %xmm7
    \\ .endif
    \\ movaps %xmm7, -8 * 4(%rcx)
    \\ .if \bp5
    \\ movaps 24 * 4(%rcx), %xmm7
    \\ mulps %xmm0, %xmm8
    \\ .endif
    \\ movaps %xmm8, -4 * 4(%rcx)
    \\ .if \bp5
    \\ movaps 28 * 4(%rcx), %xmm8
    \\ .endif
    \\ .if \bp6
    \\ movaps %xmm0, %xmm5
    \\ mulps 16 * 4(%rcx), %xmm5
    \\ movaps %xmm0, %xmm6
    \\ mulps 20 * 4(%rcx), %xmm6
    \\ movaps %xmm0, %xmm7
    \\ mulps 24 * 4(%rcx), %xmm7
    \\ movaps %xmm0, %xmm8
    \\ mulps 28 * 4(%rcx), %xmm8
    \\ .endif
    \\ subq $-32 * 4, %rcx
    \\ decq %rax
    \\ jg 20015b
    \\ .balign 16
    \\ 20016:
    \\ .if \bp5
    \\ mulps %xmm0, %xmm1
    \\ .endif
    \\ movaps %xmm1, -32 * 4(%rcx)
    \\ .if \bp5
    \\ mulps %xmm0, %xmm2
    \\ .endif
    \\ movaps %xmm2, -28 * 4(%rcx)
    \\ .if \bp5
    \\ mulps %xmm0, %xmm3
    \\ .endif
    \\ movaps %xmm3, -24 * 4(%rcx)
    \\ .if \bp5
    \\ mulps %xmm0, %xmm4
    \\ .endif
    \\ movaps %xmm4, -20 * 4(%rcx)
    \\ .if \bp5
    \\ mulps %xmm0, %xmm5
    \\ .endif
    \\ movaps %xmm5, -16 * 4(%rcx)
    \\ .if \bp5
    \\ mulps %xmm0, %xmm6
    \\ .endif
    \\ movaps %xmm6, -12 * 4(%rcx)
    \\ .if \bp5
    \\ mulps %xmm0, %xmm7
    \\ .endif
    \\ movaps %xmm7, -8 * 4(%rcx)
    \\ .if \bp5
    \\ mulps %xmm0, %xmm8
    \\ .endif
    \\ movaps %xmm8, -4 * 4(%rcx)
    \\ subq $-32 * 4, %rcx
    \\ .balign 8
    \\ 20017:
    \\ testq $31, %rdi
    \\ je 20027f
    \\ testq $16, %rdi
    \\ je 20018f
    \\ movaps -32 * 4(%rcx), %xmm1
    \\ movaps -28 * 4(%rcx), %xmm3
    \\ movaps -24 * 4(%rcx), %xmm5
    \\ movaps -20 * 4(%rcx), %xmm7
    \\ mulps %xmm0, %xmm1
    \\ movaps %xmm1, -32 * 4(%rcx)
    \\ mulps %xmm0, %xmm3
    \\ movaps %xmm3, -28 * 4(%rcx)
    \\ mulps %xmm0, %xmm5
    \\ movaps %xmm5, -24 * 4(%rcx)
    \\ mulps %xmm0, %xmm7
    \\ movaps %xmm7, -20 * 4(%rcx)
    \\ addq $16 * 4, %rcx
    \\ .balign 8
    \\ 20018:
    \\ testq $8, %rdi
    \\ je 20019f
    \\ movaps -32 * 4(%rcx), %xmm1
    \\ movaps -28 * 4(%rcx), %xmm3
    \\ mulps %xmm0, %xmm1
    \\ movaps %xmm1, -32 * 4(%rcx)
    \\ mulps %xmm0, %xmm3
    \\ movaps %xmm3, -28 * 4(%rcx)
    \\ addq $8 * 4, %rcx
    \\ .balign 8
    \\ 20019:
    \\ testq $4, %rdi
    \\ je 20020f
    \\ movaps -32 * 4(%rcx), %xmm1
    \\ mulps %xmm0, %xmm1
    \\ movaps %xmm1, -32 * 4(%rcx)
    \\ addq $4 * 4, %rcx
    \\ .balign 8
    \\ 20020:
    \\ testq $2, %rdi
    \\ je 20021f
    \\ .if \bp0
    \\ movsd -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ mulps %xmm0, %xmm1
    \\ .if \bp0
    \\ movsd %xmm1, -32 * 4(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm1, -32 * 4(%rcx)
    \\ .endif
    \\ addq $2 * 4, %rcx
    \\ .balign 8
    \\ 20021:
    \\ testq $1, %rdi
    \\ je 20027f
    \\ movss -32 * 4(%rcx), %xmm1
    \\ mulss %xmm0, %xmm1
    \\ movss %xmm1, -32 * 4(%rcx)
    \\ jmp 20027f
    \\ .balign 8
    \\ 20022:
    \\ movq %rcx, %r10
    \\ movq %rdi, %rax # rcx = n
    \\ sarq $3, %rax # (n >> 3)
    \\ jle 20024f
    \\ .balign 16
    \\ 20023:
    \\ movss (%rcx), %xmm1
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm2
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm3
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm4
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm5
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm6
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm7
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm8
    \\ addq %r8, %rcx
    \\ mulss %xmm0, %xmm1
    \\ mulss %xmm0, %xmm2
    \\ mulss %xmm0, %xmm3
    \\ mulss %xmm0, %xmm4
    \\ mulss %xmm0, %xmm5
    \\ mulss %xmm0, %xmm6
    \\ mulss %xmm0, %xmm7
    \\ mulss %xmm0, %xmm8
    \\ movss %xmm1, (%r10)
    \\ addq %r8, %r10
    \\ movss %xmm2, (%r10)
    \\ addq %r8, %r10
    \\ movss %xmm3, (%r10)
    \\ addq %r8, %r10
    \\ movss %xmm4, (%r10)
    \\ addq %r8, %r10
    \\ movss %xmm5, (%r10)
    \\ addq %r8, %r10
    \\ movss %xmm6, (%r10)
    \\ addq %r8, %r10
    \\ movss %xmm7, (%r10)
    \\ addq %r8, %r10
    \\ movss %xmm8, (%r10)
    \\ addq %r8, %r10
    \\ decq %rax
    \\ jg 20023b
    \\ .balign 16
    \\ 20024:
    \\ testq $7, %rdi
    \\ je 20027f
    \\ testq $4, %rdi
    \\ je 20025f
    \\ movss (%rcx), %xmm1
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm2
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm3
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm4
    \\ addq %r8, %rcx
    \\ mulss %xmm0, %xmm1
    \\ mulss %xmm0, %xmm2
    \\ mulss %xmm0, %xmm3
    \\ mulss %xmm0, %xmm4
    \\ movss %xmm1, (%r10)
    \\ addq %r8, %r10
    \\ movss %xmm2, (%r10)
    \\ addq %r8, %r10
    \\ movss %xmm3, (%r10)
    \\ addq %r8, %r10
    \\ movss %xmm4, (%r10)
    \\ addq %r8, %r10
    \\ .balign 8
    \\ 20025:
    \\ testq $2, %rdi
    \\ je 20026f
    \\ movss (%rcx), %xmm1
    \\ addq %r8, %rcx
    \\ movss (%rcx), %xmm2
    \\ addq %r8, %rcx
    \\ mulss %xmm0, %xmm1
    \\ mulss %xmm0, %xmm2
    \\ movss %xmm1, (%r10)
    \\ addq %r8, %r10
    \\ movss %xmm2, (%r10)
    \\ addq %r8, %r10
    \\ .balign 8
    \\ 20026:
    \\ testq $1, %rdi
    \\ je 20027f
    \\ movss (%rcx), %xmm1
    \\ mulss %xmm0, %xmm1
    \\ movss %xmm1, (%rcx)
    \\ .balign 16
    \\ 20027:
    \\ xorq %rax, %rax
    \\ ret
    \\ .endm
    \\
;
