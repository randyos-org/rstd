//! Section 7 of 7 of the text of `caxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zaxpy_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 8
    \\ 20054:
    \\ movaps -20 * 4(%rcx), %xmm4
    \\ movaps -16 * 4(%rcx), %xmm5
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ pshufd $0xb1, %xmm2, %xmm8
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm2
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps -12 * 4(%rcx), %xmm6
    \\ movaps -8 * 4(%rcx), %xmm7
    \\ movss %xmm4, %xmm3
    \\ shufps $0x93, %xmm4, %xmm3
    \\ pshufd $0xb1, %xmm3, %xmm8
    \\ mulps %xmm14, %xmm3
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm3
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ movss %xmm5, %xmm4
    \\ shufps $0x93, %xmm5, %xmm4
    \\ pshufd $0xb1, %xmm4, %xmm8
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm4
    \\ movss %xmm4, %xmm3
    \\ shufps $0x93, %xmm4, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ movss %xmm6, %xmm5
    \\ shufps $0x93, %xmm6, %xmm5
    \\ pshufd $0xb1, %xmm5, %xmm8
    \\ mulps %xmm14, %xmm5
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm5
    \\ movss %xmm5, %xmm4
    \\ shufps $0x93, %xmm5, %xmm4
    \\ addps -16 * 4(%r9), %xmm4
    \\ movaps %xmm4, -16 * 4(%r9)
    \\ movaps -4 * 4(%rcx), %xmm0
    \\ movaps 0 * 4(%rcx), %xmm1
    \\ movss %xmm7, %xmm6
    \\ shufps $0x93, %xmm7, %xmm6
    \\ pshufd $0xb1, %xmm6, %xmm8
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm6
    \\ movss %xmm6, %xmm5
    \\ shufps $0x93, %xmm6, %xmm5
    \\ addps -12 * 4(%r9), %xmm5
    \\ movaps %xmm5, -12 * 4(%r9)
    \\ movss %xmm0, %xmm7
    \\ shufps $0x93, %xmm0, %xmm7
    \\ pshufd $0xb1, %xmm7, %xmm8
    \\ mulps %xmm14, %xmm7
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm7
    \\ movss %xmm7, %xmm6
    \\ shufps $0x93, %xmm7, %xmm6
    \\ addps -8 * 4(%r9), %xmm6
    \\ movaps %xmm6, -8 * 4(%r9)
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ pshufd $0xb1, %xmm0, %xmm8
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm0
    \\ movss %xmm0, %xmm7
    \\ shufps $0x93, %xmm0, %xmm7
    \\ addps -4 * 4(%r9), %xmm7
    \\ movaps %xmm7, -4 * 4(%r9)
    \\ subq $-32 * 4, %rcx
    \\ subq $-32 * 4, %r9
    \\ .balign 8
    \\ 20055:
    \\ testq $8, %rdi
    \\ jle 20056f
    \\ movaps -28 * 4(%rcx), %xmm2
    \\ movaps -24 * 4(%rcx), %xmm3
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ pshufd $0xb1, %xmm2, %xmm8
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm2
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps -20 * 4(%rcx), %xmm0
    \\ movaps -16 * 4(%rcx), %xmm1
    \\ movss %xmm0, %xmm3
    \\ shufps $0x93, %xmm0, %xmm3
    \\ pshufd $0xb1, %xmm3, %xmm8
    \\ mulps %xmm14, %xmm3
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm3
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ addps -24 * 4(%r9), %xmm2
    \\ movaps %xmm2, -24 * 4(%r9)
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ pshufd $0xb1, %xmm0, %xmm8
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm0
    \\ movss %xmm0, %xmm3
    \\ shufps $0x93, %xmm0, %xmm3
    \\ addps -20 * 4(%r9), %xmm3
    \\ movaps %xmm3, -20 * 4(%r9)
    \\ addq $16 * 4, %rcx
    \\ addq $16 * 4, %r9
    \\ .balign 4
    \\ 20056:
    \\ testq $4, %rdi
    \\ jle 20057f
    \\ movaps -28 * 4(%rcx), %xmm2
    \\ movaps -24 * 4(%rcx), %xmm3
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movss %xmm3, %xmm2
    \\ shufps $0x93, %xmm3, %xmm2
    \\ pshufd $0xb1, %xmm2, %xmm8
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm2
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ addps -28 * 4(%r9), %xmm1
    \\ movaps %xmm1, -28 * 4(%r9)
    \\ movaps %xmm2, %xmm0
    \\ movaps %xmm3, %xmm1
    \\ addq $8 * 4, %rcx
    \\ addq $8 * 4, %r9
    \\ .balign 4
    \\ 20057:
    \\ testq $2, %rdi
    \\ jle 20058f
    \\ movaps -28 * 4(%rcx), %xmm2
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, %xmm0
    \\ movaps %xmm2, %xmm1
    \\ addq $4 * 4, %rcx
    \\ addq $4 * 4, %r9
    \\ .balign 4
    \\ 20058:
    \\ testq $1, %rdi
    \\ jle 20059f
    \\ movaps -28 * 4(%rcx), %xmm2
    \\ movss %xmm2, %xmm1
    \\ shufps $0x93, %xmm2, %xmm1
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movlps %xmm0, -32 * 4(%r9)
    \\ movhlps %xmm0, %xmm0
    \\ movss %xmm0, -30 * 4(%r9)
    \\ jmp 20067f
    \\ 20059:
    \\ shufps $0x93, %xmm0, %xmm0
    \\ addss -32 * 4(%r9), %xmm0
    \\ movss %xmm0, -32 * 4(%r9)
    \\ jmp 20067f
    \\ .balign 8
    \\ 20060:
    \\ pshufd $0, %xmm0, %xmm14
    \\ pshufd $0, %xmm1, %xmm15
    \\ pxor %xmm13, %xmm13
    \\ subps %xmm15, %xmm13
    \\ unpcklps %xmm14, %xmm13
    \\ unpcklps %xmm15, %xmm14
    \\ movaps %xmm13, %xmm15
    \\ cmpq $0, %r8
    \\ je 20065f
    \\ cmpq $0, %rsi
    \\ je 20065f
    \\ movq %r9, %r11
    \\ movq %rdi, %rax
    \\ sarq $3, %rax
    \\ jle 20062f
    \\ .balign 8
    \\ 20061:
    \\ movsd (%rcx), %xmm0
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm0
    \\ addq %r8, %rcx
    \\ movsd (%rcx), %xmm2
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm2
    \\ addq %r8, %rcx
    \\ movsd (%rcx), %xmm4
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm4
    \\ addq %r8, %rcx
    \\ movsd (%rcx), %xmm6
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm6
    \\ addq %r8, %rcx
    \\ movshdup %xmm0, %xmm1
    \\ movsldup %xmm0, %xmm0
    \\ movshdup %xmm2, %xmm3
    \\ movsldup %xmm2, %xmm2
    \\ movshdup %xmm4, %xmm5
    \\ movsldup %xmm4, %xmm4
    \\ movshdup %xmm6, %xmm7
    \\ movsldup %xmm6, %xmm6
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm5
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm7
    \\ movsd (%r9), %xmm8
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm8
    \\ addq %rsi, %r9
    \\ movsd (%r9), %xmm9
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm9
    \\ addq %rsi, %r9
    \\ movsd (%r9), %xmm10
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm10
    \\ addq %rsi, %r9
    \\ movsd (%r9), %xmm11
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm11
    \\ addq %rsi, %r9
    \\ addps %xmm0, %xmm8
    \\ addps %xmm1, %xmm8
    \\ addps %xmm2, %xmm9
    \\ addps %xmm3, %xmm9
    \\ addps %xmm4, %xmm10
    \\ addps %xmm5, %xmm10
    \\ addps %xmm6, %xmm11
    \\ addps %xmm7, %xmm11
    \\ movsd %xmm8, (%r11)
    \\ addq %rsi, %r11
    \\ movhps %xmm8, (%r11)
    \\ addq %rsi, %r11
    \\ movsd %xmm9, (%r11)
    \\ addq %rsi, %r11
    \\ movhps %xmm9, (%r11)
    \\ addq %rsi, %r11
    \\ movsd %xmm10, (%r11)
    \\ addq %rsi, %r11
    \\ movhps %xmm10, (%r11)
    \\ addq %rsi, %r11
    \\ movsd %xmm11, (%r11)
    \\ addq %rsi, %r11
    \\ movhps %xmm11, (%r11)
    \\ addq %rsi, %r11
    \\ decq %rax
    \\ jg 20061b
    \\ .balign 8
    \\ 20062:
    \\ testq $4, %rdi
    \\ jle 20063f
    \\ movsd (%rcx), %xmm0
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm0
    \\ addq %r8, %rcx
    \\ movsd (%rcx), %xmm2
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm2
    \\ addq %r8, %rcx
    \\ movshdup %xmm0, %xmm1
    \\ movsldup %xmm0, %xmm0
    \\ movshdup %xmm2, %xmm3
    \\ movsldup %xmm2, %xmm2
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ movsd (%r9), %xmm8
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm8
    \\ addq %rsi, %r9
    \\ movsd (%r9), %xmm9
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm9
    \\ addq %rsi, %r9
    \\ addps %xmm0, %xmm8
    \\ addps %xmm1, %xmm8
    \\ addps %xmm2, %xmm9
    \\ addps %xmm3, %xmm9
    \\ movsd %xmm8, (%r11)
    \\ addq %rsi, %r11
    \\ movhps %xmm8, (%r11)
    \\ addq %rsi, %r11
    \\ movsd %xmm9, (%r11)
    \\ addq %rsi, %r11
    \\ movhps %xmm9, (%r11)
    \\ addq %rsi, %r11
    \\ .balign 8
    \\ 20063:
    \\ testq $2, %rdi
    \\ jle 20064f
    \\ movsd (%rcx), %xmm0
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm0
    \\ addq %r8, %rcx
    \\ movshdup %xmm0, %xmm1
    \\ movsldup %xmm0, %xmm0
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ movsd (%r9), %xmm8
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm8
    \\ addq %rsi, %r9
    \\ addps %xmm0, %xmm8
    \\ addps %xmm1, %xmm8
    \\ movsd %xmm8, (%r11)
    \\ addq %rsi, %r11
    \\ movhps %xmm8, (%r11)
    \\ addq %rsi, %r11
    \\ .balign 8
    \\ 20064:
    \\ testq $1, %rdi
    \\ jle 20067f
    \\ movsd (%rcx), %xmm0
    \\ movshdup %xmm0, %xmm1
    \\ movsldup %xmm0, %xmm0
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ movsd (%r9), %xmm8
    \\ addps %xmm0, %xmm8
    \\ addps %xmm1, %xmm8
    \\ movsd %xmm8, (%r9)
    \\ jmp 20067f
    \\ .balign 8
    \\ 20065:
    \\ movq %rdi, %rax
    \\ cmpq $0, %rax
    \\ jle 20067f
    \\ .balign 8
    \\ 20066:
    \\ movsd (%rcx), %xmm0
    \\ addq %r8, %rcx
    \\ movshdup %xmm0, %xmm1
    \\ movsldup %xmm0, %xmm0
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ movsd (%r9), %xmm8
    \\ addps %xmm0, %xmm8
    \\ addps %xmm1, %xmm8
    \\ movsd %xmm8, (%r9)
    \\ addq %rsi, %r9
    \\ decq %rax
    \\ jg 20066b
    \\ .balign 8
    \\ 20067:
    \\ xorq %rax, %rax
    \\ ret
    \\ .endm
    \\
;
