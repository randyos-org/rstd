//! Section 2 of 2 of the text of `dgemv_t`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/dgemv_t_bulldozer.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 16
    \\ 20014:
    \\ cmpq $4, %r11
    \\ jl 20022f
    \\ subq $4, %r11
    \\ leaq 16 * 8(%r13), %r15
    \\ movq %r12, %rbx
    \\ leaq (%rbx, %r8, 2), %rcx
    \\ leaq (%rbx, %r8, 4), %r12
    \\ vxorps %xmm0 , %xmm0, %xmm0
    \\ vxorps %xmm1 , %xmm1, %xmm1
    \\ vxorps %xmm2 , %xmm2, %xmm2
    \\ vxorps %xmm3 , %xmm3, %xmm3
    \\ movq %r10, %rax
    \\ sarq $3, %rax
    \\ jle 20017f
    \\ vmovups -16 * 8(%r15), %xmm12
    \\ vmovups -14 * 8(%r15), %xmm13
    \\ decq %rax
    \\ jle 20016f
    \\ .balign 16
    \\ 20015:
    \\ prefetchnta 256(%rbx)
    \\ prefetchnta 256(%rcx)
    \\ vfmaddpd %xmm0 , -16 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm1 , -16 * 8(%rbx , %r8 , 1) , %xmm12 , %xmm1
    \\ vfmaddpd %xmm2 , -16 * 8(%rcx) , %xmm12 , %xmm2
    \\ vfmaddpd %xmm3 , -16 * 8(%rcx , %r8 , 1) , %xmm12 , %xmm3
    \\ vmovups -12 * 8(%r15), %xmm12
    \\ prefetchnta 256(%rbx,%r8,1)
    \\ prefetchnta 256(%rcx,%r8,1)
    \\ vfmaddpd %xmm0 , -14 * 8(%rbx) , %xmm13 , %xmm0
    \\ vfmaddpd %xmm1 , -14 * 8(%rbx , %r8 , 1) , %xmm13 , %xmm1
    \\ vfmaddpd %xmm2 , -14 * 8(%rcx) , %xmm13 , %xmm2
    \\ vfmaddpd %xmm3 , -14 * 8(%rcx , %r8 , 1) , %xmm13 , %xmm3
    \\ prefetchnta 256(%r15)
    \\ vmovups -10 * 8(%r15), %xmm13
    \\ vfmaddpd %xmm0 , -12 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm1 , -12 * 8(%rbx , %r8 , 1) , %xmm12 , %xmm1
    \\ vfmaddpd %xmm2 , -12 * 8(%rcx) , %xmm12 , %xmm2
    \\ vfmaddpd %xmm3 , -12 * 8(%rcx , %r8 , 1) , %xmm12 , %xmm3
    \\ vmovups -8 * 8(%r15), %xmm12
    \\ vfmaddpd %xmm0 , -10 * 8(%rbx) , %xmm13 , %xmm0
    \\ vfmaddpd %xmm1 , -10 * 8(%rbx , %r8 , 1) , %xmm13 , %xmm1
    \\ vfmaddpd %xmm2 , -10 * 8(%rcx) , %xmm13 , %xmm2
    \\ vfmaddpd %xmm3 , -10 * 8(%rcx , %r8 , 1) , %xmm13 , %xmm3
    \\ vmovups -6 * 8(%r15), %xmm13
    \\ addq $8 * 8, %rbx
    \\ addq $8 * 8, %rcx
    \\ addq $8 * 8, %r15
    \\ decq %rax
    \\ jg 20015b
    \\ .balign 16
    \\ 20016:
    \\ vfmaddpd %xmm0 , -16 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm1 , -16 * 8(%rbx , %r8 , 1) , %xmm12 , %xmm1
    \\ vfmaddpd %xmm2 , -16 * 8(%rcx) , %xmm12 , %xmm2
    \\ vfmaddpd %xmm3 , -16 * 8(%rcx , %r8 , 1) , %xmm12 , %xmm3
    \\ vmovups -12 * 8(%r15), %xmm12
    \\ vfmaddpd %xmm0 , -14 * 8(%rbx) , %xmm13 , %xmm0
    \\ vfmaddpd %xmm1 , -14 * 8(%rbx , %r8 , 1) , %xmm13 , %xmm1
    \\ vfmaddpd %xmm2 , -14 * 8(%rcx) , %xmm13 , %xmm2
    \\ vfmaddpd %xmm3 , -14 * 8(%rcx , %r8 , 1) , %xmm13 , %xmm3
    \\ vmovups -10 * 8(%r15), %xmm13
    \\ vfmaddpd %xmm0 , -12 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm1 , -12 * 8(%rbx , %r8 , 1) , %xmm12 , %xmm1
    \\ vfmaddpd %xmm2 , -12 * 8(%rcx) , %xmm12 , %xmm2
    \\ vfmaddpd %xmm3 , -12 * 8(%rcx , %r8 , 1) , %xmm12 , %xmm3
    \\ vfmaddpd %xmm0 , -10 * 8(%rbx) , %xmm13 , %xmm0
    \\ vfmaddpd %xmm1 , -10 * 8(%rbx , %r8 , 1) , %xmm13 , %xmm1
    \\ vfmaddpd %xmm2 , -10 * 8(%rcx) , %xmm13 , %xmm2
    \\ vfmaddpd %xmm3 , -10 * 8(%rcx , %r8 , 1) , %xmm13 , %xmm3
    \\ addq $8 * 8, %rbx
    \\ addq $8 * 8, %rcx
    \\ addq $8 * 8, %r15
    \\ .balign 16
    \\ 20017:
    \\ testq $4, %r10
    \\ jle 20018f
    \\ vmovups -16 * 8(%r15), %xmm12
    \\ vmovups -14 * 8(%r15), %xmm13
    \\ vfmaddpd %xmm0 , -16 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm1 , -16 * 8(%rbx , %r8 , 1) , %xmm12 , %xmm1
    \\ vfmaddpd %xmm2 , -16 * 8(%rcx) , %xmm12 , %xmm2
    \\ vfmaddpd %xmm3 , -16 * 8(%rcx , %r8 , 1) , %xmm12 , %xmm3
    \\ vfmaddpd %xmm0 , -14 * 8(%rbx) , %xmm13 , %xmm0
    \\ vfmaddpd %xmm1 , -14 * 8(%rbx , %r8 , 1) , %xmm13 , %xmm1
    \\ vfmaddpd %xmm2 , -14 * 8(%rcx) , %xmm13 , %xmm2
    \\ vfmaddpd %xmm3 , -14 * 8(%rcx , %r8 , 1) , %xmm13 , %xmm3
    \\ addq $4 * 8, %rbx
    \\ addq $4 * 8, %rcx
    \\ addq $4 * 8, %r15
    \\ .balign 16
    \\ 20018:
    \\ testq $2, %r10
    \\ jle 20019f
    \\ vmovups -16 * 8(%r15), %xmm12
    \\ vfmaddpd %xmm0 , -16 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm1 , -16 * 8(%rbx , %r8 , 1) , %xmm12 , %xmm1
    \\ vfmaddpd %xmm2 , -16 * 8(%rcx) , %xmm12 , %xmm2
    \\ vfmaddpd %xmm3 , -16 * 8(%rcx , %r8 , 1) , %xmm12 , %xmm3
    \\ addq $2 * 8, %rbx
    \\ addq $2 * 8, %rcx
    \\ addq $2 * 8, %r15
    \\ .balign 16
    \\ 20019:
    \\ testq $1, %r10
    \\ je 20020f
    \\ vmovsd -16 * 8(%r15), %xmm12
    \\ vmovsd -16 * 8(%rbx), %xmm8
    \\ vmovsd -16 * 8(%rbx, %r8), %xmm9
    \\ vmovsd -16 * 8(%rcx), %xmm10
    \\ vmovsd -16 * 8(%rcx, %r8), %xmm11
    \\ vfmaddpd %xmm0, %xmm8 , %xmm12, %xmm0
    \\ vfmaddpd %xmm1, %xmm9 , %xmm12, %xmm1
    \\ vfmaddpd %xmm2, %xmm10, %xmm12, %xmm2
    \\ vfmaddpd %xmm3, %xmm11, %xmm12, %xmm3
    \\ .balign 16
    \\ 20020:
    \\ vhaddpd %xmm1, %xmm0 , %xmm0
    \\ vhaddpd %xmm3, %xmm2 , %xmm2
    \\ vmulpd %xmm15, %xmm0 , %xmm0
    \\ vmulpd %xmm15, %xmm2 , %xmm2
    \\ cmpq $8, %rdi
    \\ jne 20021f
    \\ vmovups 0 * 8(%r14), %xmm4
    \\ vmovups 2 * 8(%r14), %xmm5
    \\ addq $4 * 8, %r14
    \\ vaddpd %xmm4, %xmm0 , %xmm0
    \\ vaddpd %xmm5, %xmm2 , %xmm2
    \\ vmovups %xmm0, 0 * 8(%rbp)
    \\ vmovups %xmm2, 2 * 8(%rbp)
    \\ addq $4 * 8, %rbp
    \\ jmp 20022f
    \\ .balign 16
    \\ 20021:
    \\ vmovsd (%r14), %xmm4
    \\ addq %rdi, %r14
    \\ vmovhpd (%r14), %xmm4 , %xmm4
    \\ addq %rdi, %r14
    \\ vmovsd (%r14), %xmm5
    \\ addq %rdi, %r14
    \\ vmovhpd (%r14), %xmm5 , %xmm5
    \\ addq %rdi, %r14
    \\ vaddpd %xmm4, %xmm0 , %xmm0
    \\ vaddpd %xmm5, %xmm2 , %xmm2
    \\ vmovlpd %xmm0, (%rbp)
    \\ addq %rdi, %rbp
    \\ vmovhpd %xmm0, (%rbp)
    \\ addq %rdi, %rbp
    \\ vmovlpd %xmm2, (%rbp)
    \\ addq %rdi, %rbp
    \\ vmovhpd %xmm2, (%rbp)
    \\ addq %rdi, %rbp
    \\ .balign 16
    \\ 20022:
    \\ cmpq $2, %r11
    \\ jl 20029f
    \\ subq $2, %r11
    \\ leaq 16 * 8(%r13), %r15
    \\ movq %r12, %rbx
    \\ leaq (%rbx, %r8), %rcx
    \\ leaq (%rbx, %r8, 2), %r12
    \\ vxorps %xmm0 , %xmm0, %xmm0
    \\ vxorps %xmm1 , %xmm1, %xmm1
    \\ vxorps %xmm2 , %xmm2, %xmm2
    \\ vxorps %xmm3 , %xmm3, %xmm3
    \\ movq %r10, %rax
    \\ sarq $3, %rax
    \\ jle 20025f
    \\ vmovups -16 * 8(%r15), %xmm12
    \\ vmovups -14 * 8(%r15), %xmm13
    \\ decq %rax
    \\ jle 20024f
    \\ .balign 16
    \\ 20023:
    \\ prefetchnta 256(%rbx)
    \\ prefetchnta 256(%rcx)
    \\ vfmaddpd %xmm0 , -16 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm1 , -16 * 8(%rcx) , %xmm12 , %xmm1
    \\ prefetchnta 256(%r15)
    \\ vmovups -12 * 8(%r15), %xmm12
    \\ vfmaddpd %xmm2 , -14 * 8(%rbx) , %xmm13 , %xmm2
    \\ vfmaddpd %xmm3 , -14 * 8(%rcx) , %xmm13 , %xmm3
    \\ vmovups -10 * 8(%r15), %xmm13
    \\ vfmaddpd %xmm0 , -12 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm1 , -12 * 8(%rcx) , %xmm12 , %xmm1
    \\ vmovups -8 * 8(%r15), %xmm12
    \\ vfmaddpd %xmm2 , -10 * 8(%rbx) , %xmm13 , %xmm2
    \\ vfmaddpd %xmm3 , -10 * 8(%rcx) , %xmm13 , %xmm3
    \\ vmovups -6 * 8(%r15), %xmm13
    \\ addq $8 * 8, %rbx
    \\ addq $8 * 8, %rcx
    \\ addq $8 * 8, %r15
    \\ decq %rax
    \\ jg 20023b
    \\ .balign 16
    \\ 20024:
    \\ vfmaddpd %xmm0 , -16 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm1 , -16 * 8(%rcx) , %xmm12 , %xmm1
    \\ vmovups -12 * 8(%r15), %xmm12
    \\ vfmaddpd %xmm2 , -14 * 8(%rbx) , %xmm13 , %xmm2
    \\ vfmaddpd %xmm3 , -14 * 8(%rcx) , %xmm13 , %xmm3
    \\ vmovups -10 * 8(%r15), %xmm13
    \\ vfmaddpd %xmm0 , -12 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm1 , -12 * 8(%rcx) , %xmm12 , %xmm1
    \\ vfmaddpd %xmm2 , -10 * 8(%rbx) , %xmm13 , %xmm2
    \\ vfmaddpd %xmm3 , -10 * 8(%rcx) , %xmm13 , %xmm3
    \\ addq $8 * 8, %rbx
    \\ addq $8 * 8, %rcx
    \\ addq $8 * 8, %r15
    \\ .balign 16
    \\ 20025:
    \\ testq $4, %r10
    \\ jle 20026f
    \\ vmovups -16 * 8(%r15), %xmm12
    \\ vmovups -14 * 8(%r15), %xmm13
    \\ vfmaddpd %xmm0 , -16 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm1 , -16 * 8(%rcx) , %xmm12 , %xmm1
    \\ vfmaddpd %xmm2 , -14 * 8(%rbx) , %xmm13 , %xmm2
    \\ vfmaddpd %xmm3 , -14 * 8(%rcx) , %xmm13 , %xmm3
    \\ addq $4 * 8, %rbx
    \\ addq $4 * 8, %rcx
    \\ addq $4 * 8, %r15
    \\ .balign 16
    \\ 20026:
    \\ testq $2, %r10
    \\ jle 20027f
    \\ vmovups -16 * 8(%r15), %xmm12
    \\ vfmaddpd %xmm0 , -16 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm1 , -16 * 8(%rcx) , %xmm12 , %xmm1
    \\ addq $2 * 8, %rbx
    \\ addq $2 * 8, %rcx
    \\ addq $2 * 8, %r15
    \\ .balign 16
    \\ 20027:
    \\ testq $1, %r10
    \\ je 20028f
    \\ vmovsd -16 * 8(%r15), %xmm12
    \\ vmovsd -16 * 8(%rbx), %xmm8
    \\ vmovsd -16 * 8(%rcx), %xmm9
    \\ vfmaddpd %xmm0, %xmm8 , %xmm12, %xmm0
    \\ vfmaddpd %xmm1, %xmm9 , %xmm12, %xmm1
    \\ .balign 16
    \\ 20028:
    \\ vaddpd %xmm2, %xmm0 , %xmm0
    \\ vaddpd %xmm3, %xmm1 , %xmm1
    \\ vhaddpd %xmm1, %xmm0 , %xmm0
    \\ mulpd %xmm15, %xmm0
    \\ vmovsd (%r14), %xmm4
    \\ addq %rdi, %r14
    \\ vmovhpd (%r14), %xmm4 , %xmm4
    \\ addq %rdi, %r14
    \\ vaddpd %xmm4, %xmm0 , %xmm0
    \\ vmovlpd %xmm0, (%rbp)
    \\ addq %rdi, %rbp
    \\ vmovhpd %xmm0, (%rbp)
    \\ addq %rdi, %rbp
    \\ .balign 16
    \\ 20029:
    \\ cmpq $1, %r11
    \\ jl 20036f
    \\ leaq 16 * 8(%r13), %r15
    \\ movq %r12, %rbx
    \\ vxorps %xmm0 , %xmm0, %xmm0
    \\ vxorps %xmm1 , %xmm1, %xmm1
    \\ vxorps %xmm2 , %xmm2, %xmm2
    \\ vxorps %xmm3 , %xmm3, %xmm3
    \\ movq %r10, %rax
    \\ sarq $3, %rax
    \\ jle 20032f
    \\ vmovups -16 * 8(%r15), %xmm12
    \\ vmovups -14 * 8(%r15), %xmm13
    \\ decq %rax
    \\ jle 20031f
    \\ .balign 16
    \\ 20030:
    \\ prefetchnta 256(%rbx)
    \\ vfmaddpd %xmm0 , -16 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm2 , -14 * 8(%rbx) , %xmm13 , %xmm2
    \\ prefetchnta 256(%r15)
    \\ vmovups -12 * 8(%r15), %xmm12
    \\ vmovups -10 * 8(%r15), %xmm13
    \\ vfmaddpd %xmm0 , -12 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm2 , -10 * 8(%rbx) , %xmm13 , %xmm2
    \\ vmovups -8 * 8(%r15), %xmm12
    \\ vmovups -6 * 8(%r15), %xmm13
    \\ addq $8 * 8, %rbx
    \\ addq $8 * 8, %r15
    \\ decq %rax
    \\ jg 20030b
    \\ .balign 16
    \\ 20031:
    \\ vfmaddpd %xmm0 , -16 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm2 , -14 * 8(%rbx) , %xmm13 , %xmm2
    \\ vmovups -12 * 8(%r15), %xmm12
    \\ vmovups -10 * 8(%r15), %xmm13
    \\ vfmaddpd %xmm0 , -12 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm2 , -10 * 8(%rbx) , %xmm13 , %xmm2
    \\ addq $8 * 8, %rbx
    \\ addq $8 * 8, %r15
    \\ .balign 16
    \\ 20032:
    \\ testq $4, %r10
    \\ jle 20033f
    \\ vmovups -16 * 8(%r15), %xmm12
    \\ vmovups -14 * 8(%r15), %xmm13
    \\ vfmaddpd %xmm0 , -16 * 8(%rbx) , %xmm12 , %xmm0
    \\ vfmaddpd %xmm2 , -14 * 8(%rbx) , %xmm13 , %xmm2
    \\ addq $4 * 8, %rbx
    \\ addq $4 * 8, %r15
    \\ .balign 16
    \\ 20033:
    \\ testq $2, %r10
    \\ jle 20034f
    \\ vmovups -16 * 8(%r15), %xmm12
    \\ vfmaddpd %xmm0 , -16 * 8(%rbx) , %xmm12 , %xmm0
    \\ addq $2 * 8, %rbx
    \\ addq $2 * 8, %r15
    \\ .balign 16
    \\ 20034:
    \\ testq $1, %r10
    \\ je 20035f
    \\ vmovsd -16 * 8(%r15), %xmm12
    \\ vmovsd -16 * 8(%rbx), %xmm8
    \\ vfmaddpd %xmm0, %xmm8 , %xmm12, %xmm0
    \\ .balign 16
    \\ 20035:
    \\ vaddpd %xmm2, %xmm0 , %xmm0
    \\ vaddpd %xmm3, %xmm1 , %xmm1
    \\ vaddpd %xmm1, %xmm0 , %xmm0
    \\ vhaddpd %xmm1, %xmm0 , %xmm0
    \\ vmulsd %xmm15, %xmm0 , %xmm0
    \\ vmovsd (%r14), %xmm4
    \\ addq %rdi, %r14
    \\ vaddsd %xmm4, %xmm0 , %xmm0
    \\ vmovlpd %xmm0, (%rbp)
    \\ addq %rdi, %rbp
    \\ .balign 16
    \\ 20036:
    \\ leaq (, %r10, 8), %rax
    \\ .if \bp0
    \\ addq %rax,72(%rsp)
    \\ .endif
    \\ .if \bp1
    \\ addq %rax,240(%rsp)
    \\ .endif
    \\ jmp 20000b
    \\ .balign 16
    \\ 20037:
    \\ movq 0(%rsp), %rbx
    \\ movq 8(%rsp), %rbp
    \\ movq 16(%rsp), %r12
    \\ movq 24(%rsp), %r13
    \\ movq 32(%rsp), %r14
    \\ movq 40(%rsp), %r15
    \\ .if \bp0
    \\ addq $128, %rsp
    \\ .endif
    \\ .if \bp1
    \\ movq 48(%rsp), %rdi
    \\ movq 56(%rsp), %rsi
    \\ vmovups 64(%rsp), %xmm6
    \\ vmovups 80(%rsp), %xmm7
    \\ vmovups 96(%rsp), %xmm8
    \\ vmovups 112(%rsp), %xmm9
    \\ vmovups 128(%rsp), %xmm10
    \\ vmovups 144(%rsp), %xmm11
    \\ vmovups 160(%rsp), %xmm12
    \\ vmovups 176(%rsp), %xmm13
    \\ vmovups 192(%rsp), %xmm14
    \\ vmovups 208(%rsp), %xmm15
    \\ addq $256, %rsp
    \\ .endif
    \\ ret
    \\ .balign 16
    \\ .endm
    \\
;
