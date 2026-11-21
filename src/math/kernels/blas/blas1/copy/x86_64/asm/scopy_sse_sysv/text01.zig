//! Section 1 of 2 of the text of `scopy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/copy_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1, bp2, bp3, bp4, bp5, bp6
    \\ .balign 512
    \\ leaq (, %rdx, 4), %rdx
    \\ leaq (, %r8, 4), %r8
    \\ cmpq $4, %rdx
    \\ jne 20037f
    \\ cmpq $4, %r8
    \\ jne 20037f
    \\ cmpq $3, %rdi
    \\ jle 20039f
    \\ subq $-32 * 4, %rsi
    \\ subq $-32 * 4, %rcx
    \\ testq $4, %rcx
    \\ je 20000f
    \\ movss -32 * 4(%rsi), %xmm0
    \\ movss %xmm0, -32 * 4(%rcx)
    \\ addq $1 * 4, %rsi
    \\ addq $1 * 4, %rcx
    \\ decq %rdi
    \\ .balign 16
    \\ 20000:
    \\ testq $2 * 4, %rcx
    \\ je 20001f
    \\ .if \bp0
    \\ movsd -32 * 4(%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rsi), %xmm0
    \\ .endif
    \\ movlps %xmm0, -32 * 4(%rcx)
    \\ addq $2 * 4, %rsi
    \\ addq $2 * 4, %rcx
    \\ subq $2, %rdi
    \\ jle 20009f
    \\ .balign 16
    \\ 20001:
    \\ testq $3 * 4, %rsi
    \\ jne 20010f
    \\ movq %rdi, %rax
    \\ sarq $5, %rax
    \\ jle 20004f
    \\ movaps -32 * 4(%rsi), %xmm0
    \\ movaps -28 * 4(%rsi), %xmm1
    \\ movaps -24 * 4(%rsi), %xmm2
    \\ movaps -20 * 4(%rsi), %xmm3
    \\ movaps -16 * 4(%rsi), %xmm4
    \\ movaps -12 * 4(%rsi), %xmm5
    \\ movaps -8 * 4(%rsi), %xmm6
    \\ movaps -4 * 4(%rsi), %xmm7
    \\ decq %rax
    \\ jle 20003f
    \\ .balign 8
    \\ 20002:
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
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ movaps 0 * 4(%rsi), %xmm0
    \\ movaps %xmm1, -28 * 4(%rcx)
    \\ movaps 4 * 4(%rsi), %xmm1
    \\ .if \bp5
    \\ prefetcht0 (( 64 * 12 + 32) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp6
    \\ prefetcht0 ((128 * 10) + 0) - 128(%rsi)
    \\ .endif
    \\ movaps %xmm2, -24 * 4(%rcx)
    \\ movaps 8 * 4(%rsi), %xmm2
    \\ movaps %xmm3, -20 * 4(%rcx)
    \\ movaps 12 * 4(%rsi), %xmm3
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
    \\ movaps %xmm4,-16 * 4(%rcx)
    \\ movaps 16 * 4(%rsi), %xmm4
    \\ movaps %xmm5,-12 * 4(%rcx)
    \\ movaps 20 * 4(%rsi), %xmm5
    \\ .if \bp5
    \\ prefetcht0 (( 64 * 12 + 32) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 64) - 128(%rsi)
    \\ .endif
    \\ movaps %xmm6, -8 * 4(%rcx)
    \\ movaps 24 * 4(%rsi), %xmm6
    \\ movaps %xmm7, -4 * 4(%rcx)
    \\ movaps 28 * 4(%rsi), %xmm7
    \\ subq $-32 * 4, %rcx
    \\ subq $-32 * 4, %rsi
    \\ decq %rax
    \\ jg 20002b
    \\ .balign 8
    \\ 20003:
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ movaps %xmm1, -28 * 4(%rcx)
    \\ movaps %xmm2, -24 * 4(%rcx)
    \\ movaps %xmm3, -20 * 4(%rcx)
    \\ movaps %xmm4, -16 * 4(%rcx)
    \\ movaps %xmm5, -12 * 4(%rcx)
    \\ movaps %xmm6, -8 * 4(%rcx)
    \\ movaps %xmm7, -4 * 4(%rcx)
    \\ subq $-32 * 4, %rcx
    \\ subq $-32 * 4, %rsi
    \\ .balign 8
    \\ 20004:
    \\ testq $16, %rdi
    \\ jle 20005f
    \\ movaps -32 * 4(%rsi), %xmm0
    \\ movaps -28 * 4(%rsi), %xmm1
    \\ movaps -24 * 4(%rsi), %xmm2
    \\ movaps -20 * 4(%rsi), %xmm3
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ movaps %xmm1, -28 * 4(%rcx)
    \\ movaps %xmm2, -24 * 4(%rcx)
    \\ movaps %xmm3, -20 * 4(%rcx)
    \\ addq $16 * 4, %rsi
    \\ addq $16 * 4, %rcx
    \\ .balign 8
    \\ 20005:
    \\ testq $8, %rdi
    \\ jle 20006f
    \\ movaps -32 * 4(%rsi), %xmm0
    \\ movaps -28 * 4(%rsi), %xmm1
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ movaps %xmm1, -28 * 4(%rcx)
    \\ addq $8 * 4, %rsi
    \\ addq $8 * 4, %rcx
    \\ .balign 8
    \\ 20006:
    \\ testq $4, %rdi
    \\ jle 20007f
    \\ movaps -32 * 4(%rsi), %xmm0
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ addq $4 * 4, %rsi
    \\ addq $4 * 4, %rcx
    \\ .balign 8
    \\ 20007:
    \\ testq $2, %rdi
    \\ jle 20008f
    \\ .if \bp0
    \\ movsd -32 * 4(%rsi), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rsi), %xmm0
    \\ .endif
    \\ movlps %xmm0, -32 * 4(%rcx)
    \\ addq $2 * 4, %rsi
    \\ addq $2 * 4, %rcx
    \\ .balign 8
    \\ 20008:
    \\ testq $1, %rdi
    \\ jle 20009f
    \\ movss -32 * 4(%rsi), %xmm0
    \\ movss %xmm0, -32 * 4(%rcx)
    \\ .balign 8
    \\ 20009:
    \\ xorq %rax,%rax
    \\ ret
    \\ .balign 8
    \\ 20010:
    \\ testq $4, %rsi
    \\ jne 20019f
    \\ movhps -32 * 4(%rsi), %xmm0
    \\ movq %rdi, %rax
    \\ sarq $5, %rax
    \\ jle 20013f
    \\ movaps -30 * 4(%rsi), %xmm1
    \\ movaps -26 * 4(%rsi), %xmm2
    \\ movaps -22 * 4(%rsi), %xmm3
    \\ movaps -18 * 4(%rsi), %xmm4
    \\ movaps -14 * 4(%rsi), %xmm5
    \\ movaps -10 * 4(%rsi), %xmm6
    \\ movaps -6 * 4(%rsi), %xmm7
    \\ decq %rax
    \\ jle 20012f
    \\ .balign 16
    \\ 20011:
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
    \\ shufps $0x4e, %xmm1, %xmm0
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ movaps -2 * 4(%rsi), %xmm0
    \\ shufps $0x4e, %xmm2, %xmm1
    \\ movaps %xmm1, -28 * 4(%rcx)
    \\ movaps 2 * 4(%rsi), %xmm1
    \\ .if \bp5
    \\ prefetcht0 (( 64 * 12 + 32) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 0) - 128(%rsi)
    \\ .endif
    \\ .if \bp6
    \\ prefetcht0 ((128 * 10) + 0) - 128(%rsi)
    \\ .endif
    \\ shufps $0x4e, %xmm3, %xmm2
    \\ movaps %xmm2, -24 * 4(%rcx)
    \\ movaps 6 * 4(%rsi), %xmm2
    \\ shufps $0x4e, %xmm4, %xmm3
    \\ movaps %xmm3, -20 * 4(%rcx)
    \\ movaps 10 * 4(%rsi), %xmm3
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
    \\ shufps $0x4e, %xmm5, %xmm4
    \\ movaps %xmm4, -16 * 4(%rcx)
    \\ movaps 14 * 4(%rsi), %xmm4
    \\ shufps $0x4e, %xmm6, %xmm5
    \\ movaps %xmm5, -12 * 4(%rcx)
    \\ movaps 18 * 4(%rsi), %xmm5
    \\ .if \bp5
    \\ prefetcht0 (( 64 * 12 + 32) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp2
    \\ prefetch ((128 * 5) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp3
    \\ prefetcht0 ((128 * 4) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp4
    \\ prefetcht0 ((128 * 12) + 64) - 128(%rsi)
    \\ .endif
    \\ .if \bp1
    \\ prefetch ((128 * 3) + 64) - 128(%rsi)
    \\ .endif
    \\ shufps $0x4e, %xmm7, %xmm6
    \\ movaps %xmm6, -8 * 4(%rcx)
    \\ movaps 22 * 4(%rsi), %xmm6
    \\ shufps $0x4e, %xmm0, %xmm7
    \\ movaps %xmm7, -4 * 4(%rcx)
    \\ movaps 26 * 4(%rsi), %xmm7
    \\ subq $-32 * 4, %rsi
    \\ subq $-32 * 4, %rcx
    \\ decq %rax
    \\ jg 20011b
    \\ .balign 8
    \\ 20012:
    \\ shufps $0x4e, %xmm1, %xmm0
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ movaps -2 * 4(%rsi), %xmm0
    \\ shufps $0x4e, %xmm2, %xmm1
    \\ movaps %xmm1, -28 * 4(%rcx)
    \\ shufps $0x4e, %xmm3, %xmm2
    \\ movaps %xmm2, -24 * 4(%rcx)
    \\ shufps $0x4e, %xmm4, %xmm3
    \\ movaps %xmm3, -20 * 4(%rcx)
    \\ shufps $0x4e, %xmm5, %xmm4
    \\ movaps %xmm4, -16 * 4(%rcx)
    \\ shufps $0x4e, %xmm6, %xmm5
    \\ movaps %xmm5, -12 * 4(%rcx)
    \\ shufps $0x4e, %xmm7, %xmm6
    \\ movaps %xmm6, -8 * 4(%rcx)
    \\ shufps $0x4e, %xmm0, %xmm7
    \\ movaps %xmm7, -4 * 4(%rcx)
    \\ subq $-32 * 4, %rsi
    \\ subq $-32 * 4, %rcx
    \\ .balign 8
    \\ 20013:
    \\ testq $16, %rdi
    \\ jle 20014f
    \\ .balign 8
    \\ movaps -30 * 4(%rsi), %xmm1
    \\ movaps -26 * 4(%rsi), %xmm2
    \\ movaps -22 * 4(%rsi), %xmm3
    \\ movaps -18 * 4(%rsi), %xmm4
    \\ shufps $0x4e, %xmm1, %xmm0
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ shufps $0x4e, %xmm2, %xmm1
    \\ movaps %xmm1, -28 * 4(%rcx)
    \\ shufps $0x4e, %xmm3, %xmm2
    \\ movaps %xmm2, -24 * 4(%rcx)
    \\ shufps $0x4e, %xmm4, %xmm3
    \\ movaps %xmm3, -20 * 4(%rcx)
    \\ movaps %xmm4, %xmm0
    \\ addq $16 * 4, %rsi
    \\ addq $16 * 4, %rcx
    \\ .balign 8
    \\ 20014:
    \\ testq $8, %rdi
    \\ jle 20015f
    \\ .balign 8
    \\ movaps -30 * 4(%rsi), %xmm1
    \\ movaps -26 * 4(%rsi), %xmm2
    \\ shufps $0x4e, %xmm1, %xmm0
    \\ shufps $0x4e, %xmm2, %xmm1
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ movaps %xmm1, -28 * 4(%rcx)
    \\ movaps %xmm2, %xmm0
    \\ addq $8 * 4, %rsi
    \\ addq $8 * 4, %rcx
    \\ .balign 8
    \\ 20015:
    \\ testq $4, %rdi
    \\ jle 20016f
    \\ .balign 8
    \\ movaps -30 * 4(%rsi), %xmm1
    \\ shufps $0x4e, %xmm1, %xmm0
    \\ movaps %xmm0, -32 * 4(%rcx)
    \\ addq $4 * 4, %rsi
    \\ addq $4 * 4, %rcx
    \\ .balign 8
    \\ 20016:
    \\ testq $2, %rdi
    \\ jle 20017f
    \\ .balign 8
    \\ .if \bp0
    \\ movsd -32 * 4(%rsi), %xmm0
    \\ movsd %xmm0, -32 * 4(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rsi), %xmm0
    \\ movlps %xmm0, -32 * 4(%rcx)
    \\ .endif
    \\ addq $2 * 4, %rsi
    \\ addq $2 * 4, %rcx
    \\ .balign 8
    \\ 20017:
    \\ testq $1, %rdi
    \\ jle 20018f
    \\ .balign 8
    \\ movss -32 * 4(%rsi), %xmm0
    \\ movss %xmm0, -32 * 4(%rcx)
    \\ addq $4, %rcx
    \\ .balign 8
    \\ 20018:
    \\ xorq %rax,%rax
    \\ ret
    \\ .balign 8
    \\ 20019:
    \\ testq $2 * 4, %rsi
    \\ jne 20028f
    \\ movaps -33 * 4(%rsi), %xmm0
    \\ movq %rdi, %rax
    \\ sarq $5, %rax
    \\ jle 20022f
    \\ movaps -29 * 4(%rsi), %xmm1
    \\ movaps -25 * 4(%rsi), %xmm2
    \\ movaps -21 * 4(%rsi), %xmm3
    \\ movaps -17 * 4(%rsi), %xmm4
    \\ movaps -13 * 4(%rsi), %xmm5
    \\ movaps -9 * 4(%rsi), %xmm6
    \\ movaps -5 * 4(%rsi), %xmm7
    \\ decq %rax
    \\ jle 20021f
    \\
;
