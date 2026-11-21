//! Section 1 of 2 of the text of `saxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/axpy_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1
    \\ .balign 512
    \\ movq 8(%rsp), %rsi
    \\ movaps %xmm0, %xmm15
    \\ shufps $0, %xmm15, %xmm15
    \\ leaq (, %r8, 4), %r8
    \\ leaq (, %rsi, 4), %rsi
    \\ testq %rdi, %rdi
    \\ jle 20009f
    \\ cmpq $4, %r8
    \\ jne 20019f
    \\ cmpq $4, %rsi
    \\ jne 20019f
    \\ subq $-32 * 4, %rcx
    \\ subq $-32 * 4, %r9
    \\ cmpq $3, %rdi
    \\ jle 20007f
    \\ testq $4, %r9
    \\ je 20000f
    \\ movss -32 * 4(%rcx), %xmm0
    \\ mulss %xmm15, %xmm0
    \\ addss -32 * 4(%r9), %xmm0
    \\ movss %xmm0, -32 * 4(%r9)
    \\ addq $1 * 4, %rcx
    \\ addq $1 * 4, %r9
    \\ decq %rdi
    \\ jle 20009f
    \\ .balign 8
    \\ 20000:
    \\ testq $4 * 2, %r9
    \\ je 20001f
    \\ .if \bp0
    \\ movsd -32 * 4(%rcx), %xmm0
    \\ movsd -32 * 4(%r9), %xmm4
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rcx), %xmm0
    \\ movlps -32 * 4(%r9), %xmm4
    \\ .endif
    \\ mulps %xmm15, %xmm0
    \\ addps %xmm4, %xmm0
    \\ .if \bp0
    \\ movsd %xmm0, -32 * 4(%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm0, -32 * 4(%r9)
    \\ .endif
    \\ addq $2 * 4, %rcx
    \\ addq $2 * 4, %r9
    \\ subq $2, %rdi
    \\ jle 20009f
    \\ .balign 8
    \\ 20001:
    \\ testq $4 * 3, %rcx
    \\ jne 20010f
    \\ movq %rdi, %rax
    \\ sarq $5, %rax
    \\ jle 20004f
    \\ movaps -32 * 4(%rcx), %xmm0
    \\ movaps -28 * 4(%rcx), %xmm1
    \\ movaps -24 * 4(%rcx), %xmm2
    \\ movaps -20 * 4(%rcx), %xmm3
    \\ decq %rax
    \\ jle 20003f
    \\ .balign 16
    \\ 20002:
    \\ movaps -16 * 4(%rcx), %xmm4
    \\ movaps -12 * 4(%rcx), %xmm5
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%r9)
    \\ .endif
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ mulps %xmm15, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps -8 * 4(%rcx), %xmm6
    \\ movaps -4 * 4(%rcx), %xmm7
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 0) - 128(%rcx)
    \\ .endif
    \\ mulps %xmm15, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ mulps %xmm15, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ movaps 0 * 4(%rcx), %xmm0
    \\ movaps 4 * 4(%rcx), %xmm1
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%r9)
    \\ .endif
    \\ mulps %xmm15, %xmm4
    \\ addps -16 * 4(%r9), %xmm4
    \\ movaps %xmm4, -16 * 4(%r9)
    \\ mulps %xmm15, %xmm5
    \\ addps -12 * 4(%r9), %xmm5
    \\ movaps %xmm5, -12 * 4(%r9)
    \\ movaps 8 * 4(%rcx), %xmm2
    \\ movaps 12 * 4(%rcx), %xmm3
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 64) - 128(%rcx)
    \\ .endif
    \\ mulps %xmm15, %xmm6
    \\ addps -8 * 4(%r9), %xmm6
    \\ movaps %xmm6, -8 * 4(%r9)
    \\ mulps %xmm15, %xmm7
    \\ addps -4 * 4(%r9), %xmm7
    \\ movaps %xmm7, -4 * 4(%r9)
    \\ subq $-32 * 4, %rcx
    \\ subq $-32 * 4, %r9
    \\ decq %rax
    \\ jg 20002b
    \\ .balign 8
    \\ 20003:
    \\ movaps -16 * 4(%rcx), %xmm4
    \\ movaps -12 * 4(%rcx), %xmm5
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ mulps %xmm15, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps -8 * 4(%rcx), %xmm6
    \\ movaps -4 * 4(%rcx), %xmm7
    \\ mulps %xmm15, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ mulps %xmm15, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ mulps %xmm15, %xmm4
    \\ addps -16 * 4(%r9), %xmm4
    \\ movaps %xmm4, -16 * 4(%r9)
    \\ mulps %xmm15, %xmm5
    \\ addps -12 * 4(%r9), %xmm5
    \\ movaps %xmm5, -12 * 4(%r9)
    \\ mulps %xmm15, %xmm6
    \\ addps -8 * 4(%r9), %xmm6
    \\ movaps %xmm6, -8 * 4(%r9)
    \\ mulps %xmm15, %xmm7
    \\ addps -4 * 4(%r9), %xmm7
    \\ movaps %xmm7, -4 * 4(%r9)
    \\ subq $-32 * 4, %rcx
    \\ subq $-32 * 4, %r9
    \\ .balign 8
    \\ 20004:
    \\ movq %rdi, %rax
    \\ andq $16, %rax
    \\ jle 20005f
    \\ .balign 8
    \\ movaps -32 * 4(%rcx), %xmm0
    \\ movaps -28 * 4(%rcx), %xmm1
    \\ movaps -24 * 4(%rcx), %xmm2
    \\ movaps -20 * 4(%rcx), %xmm3
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ mulps %xmm15, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ mulps %xmm15, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ mulps %xmm15, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ addq $16 * 4, %rcx
    \\ addq $16 * 4, %r9
    \\ .balign 8
    \\ 20005:
    \\ movq %rdi, %rax
    \\ andq $8, %rax
    \\ jle 20006f
    \\ .balign 8
    \\ movaps -32 * 4(%rcx), %xmm0
    \\ movaps -28 * 4(%rcx), %xmm1
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ mulps %xmm15, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ addq $8 * 4, %rcx
    \\ addq $8 * 4, %r9
    \\ .balign 8
    \\ 20006:
    \\ movq %rdi, %rax
    \\ andq $4, %rax
    \\ jle 20007f
    \\ .balign 8
    \\ movaps -32 * 4(%rcx), %xmm0
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ addq $4 * 4, %rcx
    \\ addq $4 * 4, %r9
    \\ .balign 8
    \\ 20007:
    \\ movq %rdi, %rax
    \\ andq $2, %rax
    \\ jle 20008f
    \\ .balign 8
    \\ .if \bp0
    \\ movsd -32 * 4(%rcx), %xmm0
    \\ movsd -32 * 4(%r9), %xmm4
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rcx), %xmm0
    \\ movlps -32 * 4(%r9), %xmm4
    \\ .endif
    \\ mulps %xmm15, %xmm0
    \\ addps %xmm4, %xmm0
    \\ .if \bp0
    \\ movsd %xmm0, -32 * 4(%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm0, -32 * 4(%r9)
    \\ .endif
    \\ addq $2 * 4, %rcx
    \\ addq $2 * 4, %r9
    \\ .balign 8
    \\ 20008:
    \\ movq %rdi, %rax
    \\ andq $1, %rax
    \\ jle 20009f
    \\ .balign 8
    \\ movss -32 * 4(%rcx), %xmm0
    \\ mulss %xmm15, %xmm0
    \\ addss -32 * 4(%r9), %xmm0
    \\ movss %xmm0, -32 * 4(%r9)
    \\ .balign 8
    \\ 20009:
    \\ xorq %rax,%rax
    \\ ret
    \\ .balign 8
    \\ 20010:
    \\ movq %rdi, %rax
    \\ sarq $5, %rax
    \\ jle 20013f
    \\ .if \bp0
    \\ movsd -32 * 4(%rcx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rcx), %xmm0
    \\ .endif
    \\ movhps -30 * 4(%rcx), %xmm0
    \\ .if \bp0
    \\ movsd -28 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps -28 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps -26 * 4(%rcx), %xmm1
    \\ .if \bp0
    \\ movsd -24 * 4(%rcx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps -24 * 4(%rcx), %xmm2
    \\ .endif
    \\ movhps -22 * 4(%rcx), %xmm2
    \\ .if \bp0
    \\ movsd -20 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps -20 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps -18 * 4(%rcx), %xmm3
    \\ decq %rax
    \\ jle 20012f
    \\ .balign 16
    \\ 20011:
    \\ .if \bp0
    \\ movsd -16 * 4(%rcx), %xmm4
    \\ .endif
    \\ .if \bp1
    \\ movlps -16 * 4(%rcx), %xmm4
    \\ .endif
    \\ movhps -14 * 4(%rcx), %xmm4
    \\ .if \bp0
    \\ movsd -12 * 4(%rcx), %xmm5
    \\ .endif
    \\ .if \bp1
    \\ movlps -12 * 4(%rcx), %xmm5
    \\ .endif
    \\ movhps -10 * 4(%rcx), %xmm5
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 0) - 128(%r9)
    \\ .endif
    \\ mulps %xmm15, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ mulps %xmm15, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ .if \bp0
    \\ movsd -8 * 4(%rcx), %xmm6
    \\ .endif
    \\ .if \bp1
    \\ movlps -8 * 4(%rcx), %xmm6
    \\ .endif
    \\ movhps -6 * 4(%rcx), %xmm6
    \\ .if \bp0
    \\ movsd -4 * 4(%rcx), %xmm7
    \\ .endif
    \\ .if \bp1
    \\ movlps -4 * 4(%rcx), %xmm7
    \\ .endif
    \\ movhps -2 * 4(%rcx), %xmm7
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 0) - 128(%rcx)
    \\ .endif
    \\ mulps %xmm15, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ mulps %xmm15, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ .if \bp0
    \\ movsd 0 * 4(%rcx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps 0 * 4(%rcx), %xmm0
    \\ .endif
    \\ movhps 2 * 4(%rcx), %xmm0
    \\ .if \bp0
    \\ movsd 4 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps 4 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps 6 * 4(%rcx), %xmm1
    \\ .if \bp1
    \\ prefetchw ((128 * 3) + 64) - 128(%r9)
    \\ .endif
    \\ mulps %xmm15, %xmm4
    \\ addps -16 * 4(%r9), %xmm4
    \\ movaps %xmm4, -16 * 4(%r9)
    \\ mulps %xmm15, %xmm5
    \\ addps -12 * 4(%r9), %xmm5
    \\ movaps %xmm5, -12 * 4(%r9)
    \\ .if \bp0
    \\ movsd 8 * 4(%rcx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps 8 * 4(%rcx), %xmm2
    \\ .endif
    \\ movhps 10 * 4(%rcx), %xmm2
    \\ .if \bp0
    \\ movsd 12 * 4(%rcx), %xmm3
    \\ .endif
    \\ .if \bp1
    \\ movlps 12 * 4(%rcx), %xmm3
    \\ .endif
    \\ movhps 14 * 4(%rcx), %xmm3
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 64) - 128(%rcx)
    \\ .endif
    \\ mulps %xmm15, %xmm6
    \\ addps -8 * 4(%r9), %xmm6
    \\ movaps %xmm6, -8 * 4(%r9)
    \\ mulps %xmm15, %xmm7
    \\ addps -4 * 4(%r9), %xmm7
    \\ movaps %xmm7, -4 * 4(%r9)
    \\ subq $-32 * 4, %rcx
    \\ subq $-32 * 4, %r9
    \\ decq %rax
    \\ jg 20011b
    \\
;
