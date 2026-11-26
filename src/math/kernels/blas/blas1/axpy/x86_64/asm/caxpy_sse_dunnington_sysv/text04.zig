//! Section 4 of 4 of the text of `caxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zaxpy_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 4
    \\ 20027:
    \\ testq $2, %rdi
    \\ jle 20028f
    \\ .if \bp0
    \\ movsd -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ movlps -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ movhps -30 * 4(%rcx), %xmm1
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ movaps %xmm1, %xmm0
    \\ addq $4 * 4, %rcx
    \\ addq $4 * 4, %r9
    \\ .balign 4
    \\ 20028:
    \\ testq $1, %rdi
    \\ jle 20029f
    \\ .if \bp0
    \\ movsd -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ .if \bp1
    \\ xorps %xmm1, %xmm1
    \\ movlps -32 * 4(%rcx), %xmm1
    \\ .endif
    \\ pshufd $0xb1, %xmm1, %xmm8
    \\ mulps %xmm14, %xmm1
    \\ mulps %xmm15, %xmm8
    \\ addps %xmm8, %xmm1
    \\ movss %xmm1, %xmm0
    \\ shufps $0x93, %xmm1, %xmm0
    \\ addps -32 * 4(%r9), %xmm0
    \\ movaps %xmm0, -32 * 4(%r9)
    \\ jmp 20037f
    \\ .balign 8
    \\ 20029:
    \\ shufps $0x93, %xmm0, %xmm0
    \\ addss -32 * 4(%r9), %xmm0
    \\ movss %xmm0, -32 * 4(%r9)
    \\ jmp 20037f
    \\ .balign 8
    \\ 20030:
    \\ pshufd $0, %xmm0, %xmm14
    \\ pshufd $0, %xmm1, %xmm15
    \\ pxor %xmm13, %xmm13
    \\ subps %xmm15, %xmm13
    \\ unpcklps %xmm14, %xmm13
    \\ unpcklps %xmm15, %xmm14
    \\ movaps %xmm13, %xmm15
    \\ cmpq $0, %r8
    \\ je 20035f
    \\ cmpq $0, %rsi
    \\ je 20035f
    \\ movq %r9, %r11
    \\ movq %rdi, %rax
    \\ sarq $3, %rax
    \\ jle 20032f
    \\ .balign 8
    \\ 20031:
    \\ .if \bp0
    \\ movsd (%rcx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rcx), %xmm0
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm0
    \\ addq %r8, %rcx
    \\ .if \bp0
    \\ movsd (%rcx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rcx), %xmm2
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm2
    \\ addq %r8, %rcx
    \\ .if \bp0
    \\ movsd (%rcx), %xmm4
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rcx), %xmm4
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm4
    \\ addq %r8, %rcx
    \\ .if \bp0
    \\ movsd (%rcx), %xmm6
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rcx), %xmm6
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm6
    \\ addq %r8, %rcx
    \\ .if \bp2
    \\ movshdup %xmm0, %xmm1
    \\ movsldup %xmm0, %xmm0
    \\ movshdup %xmm2, %xmm3
    \\ movsldup %xmm2, %xmm2
    \\ movshdup %xmm4, %xmm5
    \\ movsldup %xmm4, %xmm4
    \\ movshdup %xmm6, %xmm7
    \\ movsldup %xmm6, %xmm6
    \\ .endif
    \\ .if \bp3
    \\ pshufd $0xf5, %xmm0, %xmm1
    \\ shufps $0xa0, %xmm0, %xmm0
    \\ pshufd $0xf5, %xmm2, %xmm3
    \\ shufps $0xa0, %xmm2, %xmm2
    \\ pshufd $0xf5, %xmm4, %xmm5
    \\ shufps $0xa0, %xmm4, %xmm4
    \\ pshufd $0xf5, %xmm6, %xmm7
    \\ shufps $0xa0, %xmm6, %xmm6
    \\ .endif
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ mulps %xmm14, %xmm4
    \\ mulps %xmm15, %xmm5
    \\ mulps %xmm14, %xmm6
    \\ mulps %xmm15, %xmm7
    \\ .if \bp0
    \\ movsd (%r9), %xmm8
    \\ .endif
    \\ .if \bp1
    \\ movlps (%r9), %xmm8
    \\ .endif
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm8
    \\ addq %rsi, %r9
    \\ .if \bp0
    \\ movsd (%r9), %xmm9
    \\ .endif
    \\ .if \bp1
    \\ movlps (%r9), %xmm9
    \\ .endif
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm9
    \\ addq %rsi, %r9
    \\ .if \bp0
    \\ movsd (%r9), %xmm10
    \\ .endif
    \\ .if \bp1
    \\ movlps (%r9), %xmm10
    \\ .endif
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm10
    \\ addq %rsi, %r9
    \\ .if \bp0
    \\ movsd (%r9), %xmm11
    \\ .endif
    \\ .if \bp1
    \\ movlps (%r9), %xmm11
    \\ .endif
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
    \\ .if \bp0
    \\ movsd %xmm8, (%r11)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm8, (%r11)
    \\ .endif
    \\ addq %rsi, %r11
    \\ movhps %xmm8, (%r11)
    \\ addq %rsi, %r11
    \\ .if \bp0
    \\ movsd %xmm9, (%r11)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm9, (%r11)
    \\ .endif
    \\ addq %rsi, %r11
    \\ movhps %xmm9, (%r11)
    \\ addq %rsi, %r11
    \\ .if \bp0
    \\ movsd %xmm10, (%r11)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm10, (%r11)
    \\ .endif
    \\ addq %rsi, %r11
    \\ movhps %xmm10, (%r11)
    \\ addq %rsi, %r11
    \\ .if \bp0
    \\ movsd %xmm11, (%r11)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm11, (%r11)
    \\ .endif
    \\ addq %rsi, %r11
    \\ movhps %xmm11, (%r11)
    \\ addq %rsi, %r11
    \\ decq %rax
    \\ jg 20031b
    \\ .balign 8
    \\ 20032:
    \\ testq $4, %rdi
    \\ jle 20033f
    \\ .if \bp0
    \\ movsd (%rcx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rcx), %xmm0
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm0
    \\ addq %r8, %rcx
    \\ .if \bp0
    \\ movsd (%rcx), %xmm2
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rcx), %xmm2
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm2
    \\ addq %r8, %rcx
    \\ .if \bp2
    \\ movshdup %xmm0, %xmm1
    \\ movsldup %xmm0, %xmm0
    \\ movshdup %xmm2, %xmm3
    \\ movsldup %xmm2, %xmm2
    \\ .endif
    \\ .if \bp3
    \\ pshufd $0xf5, %xmm0, %xmm1
    \\ shufps $0xa0, %xmm0, %xmm0
    \\ pshufd $0xf5, %xmm2, %xmm3
    \\ shufps $0xa0, %xmm2, %xmm2
    \\ .endif
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ mulps %xmm14, %xmm2
    \\ mulps %xmm15, %xmm3
    \\ .if \bp0
    \\ movsd (%r9), %xmm8
    \\ .endif
    \\ .if \bp1
    \\ movlps (%r9), %xmm8
    \\ .endif
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm8
    \\ addq %rsi, %r9
    \\ .if \bp0
    \\ movsd (%r9), %xmm9
    \\ .endif
    \\ .if \bp1
    \\ movlps (%r9), %xmm9
    \\ .endif
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm9
    \\ addq %rsi, %r9
    \\ addps %xmm0, %xmm8
    \\ addps %xmm1, %xmm8
    \\ addps %xmm2, %xmm9
    \\ addps %xmm3, %xmm9
    \\ .if \bp0
    \\ movsd %xmm8, (%r11)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm8, (%r11)
    \\ .endif
    \\ addq %rsi, %r11
    \\ movhps %xmm8, (%r11)
    \\ addq %rsi, %r11
    \\ .if \bp0
    \\ movsd %xmm9, (%r11)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm9, (%r11)
    \\ .endif
    \\ addq %rsi, %r11
    \\ movhps %xmm9, (%r11)
    \\ addq %rsi, %r11
    \\ .balign 8
    \\ 20033:
    \\ testq $2, %rdi
    \\ jle 20034f
    \\ .if \bp0
    \\ movsd (%rcx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rcx), %xmm0
    \\ .endif
    \\ addq %r8, %rcx
    \\ movhps (%rcx), %xmm0
    \\ addq %r8, %rcx
    \\ .if \bp2
    \\ movshdup %xmm0, %xmm1
    \\ movsldup %xmm0, %xmm0
    \\ .endif
    \\ .if \bp3
    \\ pshufd $0xf5, %xmm0, %xmm1
    \\ shufps $0xa0, %xmm0, %xmm0
    \\ .endif
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ .if \bp0
    \\ movsd (%r9), %xmm8
    \\ .endif
    \\ .if \bp1
    \\ movlps (%r9), %xmm8
    \\ .endif
    \\ addq %rsi, %r9
    \\ movhps (%r9), %xmm8
    \\ addq %rsi, %r9
    \\ addps %xmm0, %xmm8
    \\ addps %xmm1, %xmm8
    \\ .if \bp0
    \\ movsd %xmm8, (%r11)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm8, (%r11)
    \\ .endif
    \\ addq %rsi, %r11
    \\ movhps %xmm8, (%r11)
    \\ addq %rsi, %r11
    \\ .balign 8
    \\ 20034:
    \\ testq $1, %rdi
    \\ jle 20037f
    \\ .if \bp0
    \\ movsd (%rcx), %xmm0
    \\ .endif
    \\ .if \bp4
    \\ movlps (%rcx), %xmm0
    \\ .endif
    \\ .if \bp2
    \\ movshdup %xmm0, %xmm1
    \\ movsldup %xmm0, %xmm0
    \\ .endif
    \\ .if \bp3
    \\ movlps (%rcx), %xmm0
    \\ pshufd $0xf5, %xmm0, %xmm1
    \\ shufps $0xa0, %xmm0, %xmm0
    \\ .endif
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ .if \bp0
    \\ movsd (%r9), %xmm8
    \\ .endif
    \\ .if \bp1
    \\ movlps (%r9), %xmm8
    \\ .endif
    \\ addps %xmm0, %xmm8
    \\ addps %xmm1, %xmm8
    \\ .if \bp0
    \\ movsd %xmm8, (%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm8, (%r9)
    \\ .endif
    \\ jmp 20037f
    \\ .balign 8
    \\ 20035:
    \\ movq %rdi, %rax
    \\ cmpq $0, %rax
    \\ jle 20037f
    \\ .balign 8
    \\ 20036:
    \\ .if \bp0
    \\ movsd (%rcx), %xmm0
    \\ .endif
    \\ .if \bp1
    \\ movlps (%rcx), %xmm0
    \\ .endif
    \\ addq %r8, %rcx
    \\ .if \bp2
    \\ movshdup %xmm0, %xmm1
    \\ movsldup %xmm0, %xmm0
    \\ .endif
    \\ .if \bp3
    \\ pshufd $0xf5, %xmm0, %xmm1
    \\ shufps $0xa0, %xmm0, %xmm0
    \\ .endif
    \\ mulps %xmm14, %xmm0
    \\ mulps %xmm15, %xmm1
    \\ .if \bp0
    \\ movsd (%r9), %xmm8
    \\ .endif
    \\ .if \bp1
    \\ movlps (%r9), %xmm8
    \\ .endif
    \\ addps %xmm0, %xmm8
    \\ addps %xmm1, %xmm8
    \\ .if \bp0
    \\ movsd %xmm8, (%r9)
    \\ .endif
    \\ .if \bp1
    \\ movlps %xmm8, (%r9)
    \\ .endif
    \\ addq %rsi, %r9
    \\ decq %rax
    \\ jg 20036b
    \\ .balign 8
    \\ 20037:
    \\ xorq %rax, %rax
    \\ ret
    \\ .endm
    \\
;
