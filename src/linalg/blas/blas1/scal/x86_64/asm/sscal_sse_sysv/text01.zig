//! Section 1 of 2 of the text of `sscal`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/scal_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .macro __bk bp0, bp1, bp2, bp3, bp4, bp5, bp6, bp7, bp8
    \\ .balign 512
    \\ movq 24(%rsp), %r9
    \\ testq %rdi, %rdi
    \\ jle 20027f
    \\ lea (, %r8, 4), %r8
    \\ pxor %xmm1, %xmm1
    \\ comiss %xmm0, %xmm1
    \\ shufps $0, %xmm0, %xmm0
    \\ jne 20012f # Alpha != 0.e0f
    \\ jp 20012f
    \\ cmpq $1, %r9
    \\ je 20012f
    \\ cmpq $4, %r8
    \\ jne 20007f
    \\ cmpq $3, %rdi
    \\ jle 20005f
    \\ testq $4, %rcx # aligned for double word?
    \\ je 20000f
    \\ movss %xmm1, 0 * 4(%rcx)
    \\ addq $4, %rcx
    \\ decq %rdi
    \\ jle 20027f
    \\ .balign 8
    \\ 20000:
    \\ testq $8, %rcx # aligned for quad word?
    \\ je 20001f
    \\ .if \bp0
    \\ movsd %xmm1, 0 * 4(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm1, 0 * 4(%rcx)
    \\ .endif
    \\ addq $2 * 4, %rcx
    \\ subq $2, %rdi
    \\ jle 20027f
    \\ .balign 8
    \\ 20001:
    \\ movq %rdi, %rax
    \\ sarq $4, %rax
    \\ jle 20003f
    \\ .balign 16
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
    \\ movaps %xmm1, 0 * 4(%rcx)
    \\ movaps %xmm1, 4 * 4(%rcx)
    \\ movaps %xmm1, 8 * 4(%rcx)
    \\ movaps %xmm1, 12 * 4(%rcx)
    \\ addq $16 * 4, %rcx
    \\ decq %rax
    \\ jg 20002b
    \\ .balign 16
    \\ 20003:
    \\ testq $15, %rdi
    \\ je 20027f
    \\ testq $8, %rdi
    \\ je 20004f
    \\ movaps %xmm1, 0 * 4(%rcx)
    \\ movaps %xmm1, 4 * 4(%rcx)
    \\ addq $8 * 4, %rcx
    \\ .balign 8
    \\ 20004:
    \\ testq $4, %rdi
    \\ je 20005f
    \\ movaps %xmm1, 0 * 4(%rcx)
    \\ addq $4 * 4, %rcx
    \\ .balign 8
    \\ 20005:
    \\ testq $2, %rdi
    \\ je 20006f
    \\ .if \bp0
    \\ movsd %xmm1, 0 * 4(%rcx)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm1, 0 * 4(%rcx)
    \\ .endif
    \\ addq $2 * 4, %rcx
    \\ .balign 8
    \\ 20006:
    \\ testq $1, %rdi
    \\ je 20027f
    \\ movss %xmm1, 0 * 4(%rcx)
    \\ jmp 20027f
    \\ .balign 16
    \\ 20007:
    \\ movq %rdi, %rax # rcx = n
    \\ sarq $3, %rax # (n >> 3)
    \\ jle 20009f
    \\ .balign 16
    \\ 20008:
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
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ decq %rax
    \\ jg 20008b
    \\ .balign 16
    \\ 20009:
    \\ testq $7, %rdi
    \\ je 20027f
    \\ testq $4, %rdi
    \\ je 20010f
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ .balign 8
    \\ 20010:
    \\ testq $2, %rdi
    \\ je 20011f
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ movss %xmm1, (%rcx)
    \\ addq %r8, %rcx
    \\ .balign 8
    \\ 20011:
    \\ testq $1, %rdi
    \\ je 20027f
    \\ movss %xmm1, (%rcx)
    \\ jmp 20027f
    \\ .balign 16
    \\ 20012:
    \\ cmpq $4, %r8
    \\ jne 20022f
    \\ subq $-32 * 4, %rcx
    \\ cmpq $3, %rdi
    \\ jle 20020f
    \\ testq $4, %rcx
    \\ je 20013f
    \\ movss -32 * 4(%rcx), %xmm1
    \\ mulss %xmm0, %xmm1
    \\ movss %xmm1, -32 * 4(%rcx)
    \\ addq $4, %rcx
    \\ decq %rdi
    \\ jle 20027f
    \\ .balign 8
    \\ 20013:
    \\ testq $2 * 4, %rcx
    \\ je 20014f
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
    \\ subq $2, %rdi
    \\ jle 20027f
    \\ .balign 8
    \\ 20014:
    \\ movq %rdi, %rax
    \\ sarq $5, %rax
    \\ jle 20017f
    \\ .if \bp5
    \\ movaps -32 * 4(%rcx), %xmm1
    \\ movaps -28 * 4(%rcx), %xmm2
    \\ movaps -24 * 4(%rcx), %xmm3
    \\ movaps -20 * 4(%rcx), %xmm4
    \\ movaps -16 * 4(%rcx), %xmm5
    \\ movaps -12 * 4(%rcx), %xmm6
    \\ movaps -8 * 4(%rcx), %xmm7
    \\ movaps -4 * 4(%rcx), %xmm8
    \\ .endif
    \\ .if \bp6
    \\ movaps %xmm0, %xmm1
    \\ mulps -32 * 4(%rcx), %xmm1
    \\ movaps %xmm0, %xmm2
    \\ mulps -28 * 4(%rcx), %xmm2
    \\ movaps %xmm0, %xmm3
    \\ mulps -24 * 4(%rcx), %xmm3
    \\ movaps %xmm0, %xmm4
    \\ mulps -20 * 4(%rcx), %xmm4
    \\ movaps %xmm0, %xmm5
    \\ mulps -16 * 4(%rcx), %xmm5
    \\ movaps %xmm0, %xmm6
    \\ mulps -12 * 4(%rcx), %xmm6
    \\ movaps %xmm0, %xmm7
    \\ mulps -8 * 4(%rcx), %xmm7
    \\ movaps %xmm0, %xmm8
    \\ mulps -4 * 4(%rcx), %xmm8
    \\ .endif
    \\ decq %rax
    \\ jle 20016f
    \\
;
