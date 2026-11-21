//! Section 1 of 2 of the text of `zaxpy`, which `__root__.zig` joins, from OpenBLAS `kernel/x86_64/zaxpy_atom.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

pub const text =
    \\ .balign 512
    \\ movaps %xmm3, %xmm0
    \\ movsd 40(%rsp), %xmm1
    \\ movq 48(%rsp), %rdx
    \\ movq 56(%rsp), %r8
    \\ movq 64(%rsp), %r9
    \\ movq 72(%rsp), %r10
    \\ subq $256, %rsp
    \\ movups %xmm6, 0(%rsp)
    \\ movups %xmm7, 16(%rsp)
    \\ movups %xmm8, 32(%rsp)
    \\ movups %xmm9, 48(%rsp)
    \\ movups %xmm10, 64(%rsp)
    \\ movups %xmm11, 80(%rsp)
    \\ movups %xmm12, 96(%rsp)
    \\ movups %xmm13, 112(%rsp)
    \\ movups %xmm14, 128(%rsp)
    \\ movups %xmm15, 144(%rsp)
    \\ salq $4, %r8
    \\ movaps %xmm0, %xmm14
    \\ salq $4, %r10
    \\ movaps %xmm1, %xmm15
    \\ testq %rcx, %rcx
    \\ jle 20010f
    \\ cmpq $2 * 8, %r8
    \\ jne 20005f
    \\ cmpq $2 * 8, %r10
    \\ jne 20005f
    \\ movq %rcx, %rax
    \\ sarq $2, %rax
    \\ jle 20002f
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ movsd 1 * 8(%rdx), %xmm1
    \\ movsd 0 * 8(%r9), %xmm8
    \\ movsd 1 * 8(%r9), %xmm9
    \\ movsd 2 * 8(%rdx), %xmm4
    \\ movsd 3 * 8(%rdx), %xmm5
    \\ movsd 2 * 8(%r9), %xmm10
    \\ movsd 3 * 8(%r9), %xmm11
    \\ movaps %xmm0, %xmm2
    \\ mulsd %xmm14, %xmm0
    \\ movaps %xmm1, %xmm3
    \\ mulsd %xmm14, %xmm1
    \\ mulsd %xmm15, %xmm3
    \\ mulsd %xmm15, %xmm2
    \\ movaps %xmm4, %xmm6
    \\ mulsd %xmm14, %xmm4
    \\ addsd %xmm0, %xmm8
    \\ movsd 4 * 8(%rdx), %xmm0
    \\ movaps %xmm5, %xmm7
    \\ mulsd %xmm14, %xmm5
    \\ addsd %xmm1, %xmm9
    \\ movsd 5 * 8(%rdx), %xmm1
    \\ decq %rax
    \\ jle 20001f
    \\ .balign 8
    \\ 20000:
    \\ prefetcht0 (( 64 * 12 + 32) + 0) - 128(%rdx)
    \\ mulsd %xmm15, %xmm7
    \\ movsd 4 * 8(%r9), %xmm12
    \\ subsd %xmm3, %xmm8
    \\ mulsd %xmm15, %xmm6
    \\ movsd 5 * 8(%r9), %xmm13
    \\ addsd %xmm2, %xmm9
    \\ addsd %xmm4, %xmm10
    \\ movsd 6 * 8(%rdx), %xmm4
    \\ movaps %xmm0, %xmm2
    \\ mulsd %xmm14, %xmm0
    \\ addsd %xmm5, %xmm11
    \\ movsd 7 * 8(%rdx), %xmm5
    \\ movaps %xmm1, %xmm3
    \\ mulsd %xmm14, %xmm1
    \\ subsd %xmm7, %xmm10
    \\ movsd %xmm8, 0 * 8(%r9)
    \\ mulsd %xmm15, %xmm3
    \\ addsd %xmm6, %xmm11
    \\ movsd %xmm9, 1 * 8(%r9)
    \\ mulsd %xmm15, %xmm2
    \\ movaps %xmm4, %xmm6
    \\ movsd %xmm10, 2 * 8(%r9)
    \\ mulsd %xmm14, %xmm4
    \\ movsd 6 * 8(%r9), %xmm10
    \\ addsd %xmm0, %xmm12
    \\ movsd 8 * 8(%rdx), %xmm0
    \\ movaps %xmm5, %xmm7
    \\ movsd %xmm11, 3 * 8(%r9)
    \\ mulsd %xmm14, %xmm5
    \\ movsd 7 * 8(%r9), %xmm11
    \\ addsd %xmm1, %xmm13
    \\ movsd 9 * 8(%rdx), %xmm1
    \\ mulsd %xmm15, %xmm7
    \\ movsd 8 * 8(%r9), %xmm8
    \\ subsd %xmm3, %xmm12
    \\ mulsd %xmm15, %xmm6
    \\ movsd 9 * 8(%r9), %xmm9
    \\ addsd %xmm2, %xmm13
    \\ movaps %xmm0, %xmm2
    \\ mulsd %xmm14, %xmm0
    \\ addsd %xmm4, %xmm10
    \\ movsd 10 * 8(%rdx), %xmm4
    \\ movaps %xmm1, %xmm3
    \\ mulsd %xmm14, %xmm1
    \\ addsd %xmm5, %xmm11
    \\ movsd 11 * 8(%rdx), %xmm5
    \\ mulsd %xmm15, %xmm3
    \\ movsd %xmm12, 4 * 8(%r9)
    \\ subsd %xmm7, %xmm10
    \\ mulsd %xmm15, %xmm2
    \\ movsd %xmm13, 5 * 8(%r9)
    \\ addsd %xmm6, %xmm11
    \\ movaps %xmm4, %xmm6
    \\ movsd %xmm10, 6 * 8(%r9)
    \\ mulsd %xmm14, %xmm4
    \\ addsd %xmm0, %xmm8
    \\ movsd 10 * 8(%r9), %xmm10
    \\ movsd 12 * 8(%rdx), %xmm0
    \\ movaps %xmm5, %xmm7
    \\ movsd %xmm11, 7 * 8(%r9)
    \\ mulsd %xmm14, %xmm5
    \\ movsd 11 * 8(%r9), %xmm11
    \\ addsd %xmm1, %xmm9
    \\ movsd 13 * 8(%rdx), %xmm1
    \\ addq $8 * 8, %rdx
    \\ addq $8 * 8, %r9
    \\ decq %rax
    \\ jg 20000b
    \\ .balign 8
    \\ 20001:
    \\ mulsd %xmm15, %xmm7
    \\ movsd 4 * 8(%r9), %xmm12
    \\ subsd %xmm3, %xmm8
    \\ mulsd %xmm15, %xmm6
    \\ movsd 5 * 8(%r9), %xmm13
    \\ addsd %xmm2, %xmm9
    \\ addsd %xmm4, %xmm10
    \\ movsd 6 * 8(%rdx), %xmm4
    \\ movaps %xmm0, %xmm2
    \\ mulsd %xmm14, %xmm0
    \\ addsd %xmm5, %xmm11
    \\ movsd 7 * 8(%rdx), %xmm5
    \\ movaps %xmm1, %xmm3
    \\ mulsd %xmm14, %xmm1
    \\ subsd %xmm7, %xmm10
    \\ movsd %xmm8, 0 * 8(%r9)
    \\ mulsd %xmm15, %xmm3
    \\ addsd %xmm6, %xmm11
    \\ movsd %xmm9, 1 * 8(%r9)
    \\ mulsd %xmm15, %xmm2
    \\ movaps %xmm4, %xmm6
    \\ movsd %xmm10, 2 * 8(%r9)
    \\ mulsd %xmm14, %xmm4
    \\ movsd 6 * 8(%r9), %xmm10
    \\ addsd %xmm0, %xmm12
    \\ movaps %xmm5, %xmm7
    \\ movsd %xmm11, 3 * 8(%r9)
    \\ mulsd %xmm14, %xmm5
    \\ addsd %xmm1, %xmm13
    \\ movsd 7 * 8(%r9), %xmm11
    \\ mulsd %xmm15, %xmm7
    \\ subsd %xmm3, %xmm12
    \\ mulsd %xmm15, %xmm6
    \\ addsd %xmm2, %xmm13
    \\ movaps %xmm0, %xmm2
    \\ mulsd %xmm14, %xmm0
    \\ addsd %xmm4, %xmm10
    \\ movaps %xmm1, %xmm3
    \\ mulsd %xmm14, %xmm1
    \\ addsd %xmm5, %xmm11
    \\ mulsd %xmm15, %xmm3
    \\ subsd %xmm7, %xmm10
    \\ addsd %xmm6, %xmm11
    \\ mulsd %xmm15, %xmm2
    \\ movsd %xmm12, 4 * 8(%r9)
    \\ movsd %xmm13, 5 * 8(%r9)
    \\ movsd %xmm10, 6 * 8(%r9)
    \\ movsd %xmm11, 7 * 8(%r9)
    \\ addq $8 * 8, %rdx
    \\ addq $8 * 8, %r9
    \\ .balign 8
    \\ 20002:
    \\ movq %rcx, %rax
    \\ andq $2, %rax
    \\ jle 20003f
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ movsd 1 * 8(%rdx), %xmm1
    \\ movsd 2 * 8(%rdx), %xmm4
    \\ movsd 3 * 8(%rdx), %xmm5
    \\ movaps %xmm0, %xmm2
    \\ movsd 0 * 8(%r9), %xmm8
    \\ mulsd %xmm14, %xmm0
    \\ movaps %xmm1, %xmm3
    \\ movsd 1 * 8(%r9), %xmm9
    \\ mulsd %xmm14, %xmm1
    \\ movsd 2 * 8(%r9), %xmm10
    \\ mulsd %xmm15, %xmm3
    \\ movsd 3 * 8(%r9), %xmm11
    \\ mulsd %xmm15, %xmm2
    \\ movaps %xmm4, %xmm6
    \\ mulsd %xmm14, %xmm4
    \\ addsd %xmm0, %xmm8
    \\ movaps %xmm5, %xmm7
    \\ mulsd %xmm14, %xmm5
    \\ addsd %xmm1, %xmm9
    \\ mulsd %xmm15, %xmm7
    \\ subsd %xmm3, %xmm8
    \\ mulsd %xmm15, %xmm6
    \\ addsd %xmm2, %xmm9
    \\ addsd %xmm4, %xmm10
    \\ movsd %xmm8, 0 * 8(%r9)
    \\ addsd %xmm5, %xmm11
    \\ movsd %xmm9, 1 * 8(%r9)
    \\ subsd %xmm7, %xmm10
    \\ addsd %xmm6, %xmm11
    \\ movsd %xmm10, 2 * 8(%r9)
    \\ movsd %xmm11, 3 * 8(%r9)
    \\ addq $4 * 8, %rdx
    \\ addq $4 * 8, %r9
    \\ .balign 8
    \\ 20003:
    \\ movq %rcx, %rax
    \\ andq $1, %rax
    \\ jle 20010f
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ movsd 1 * 8(%rdx), %xmm1
    \\ movsd 0 * 8(%r9), %xmm8
    \\ movsd 1 * 8(%r9), %xmm9
    \\ movaps %xmm0, %xmm2
    \\ mulsd %xmm14, %xmm0
    \\ movaps %xmm1, %xmm3
    \\ mulsd %xmm14, %xmm1
    \\ mulsd %xmm15, %xmm3
    \\ mulsd %xmm15, %xmm2
    \\ addsd %xmm0, %xmm8
    \\ addsd %xmm1, %xmm9
    \\ subsd %xmm3, %xmm8
    \\ addsd %xmm2, %xmm9
    \\ movsd %xmm8, 0 * 8(%r9)
    \\ movsd %xmm9, 1 * 8(%r9)
    \\ jmp 20010f
    \\ .balign 8
    \\ 20004:
    \\ movsd 0 * 8(%rdx), %xmm0
    \\ movsd 1 * 8(%rdx), %xmm1
    \\ movsd 0 * 8(%r9), %xmm8
    \\ movsd 1 * 8(%r9), %xmm9
    \\ movaps %xmm0, %xmm2
    \\ mulsd %xmm14, %xmm0
    \\ movaps %xmm1, %xmm3
    \\ mulsd %xmm14, %xmm1
    \\ mulsd %xmm15, %xmm3
    \\ mulsd %xmm15, %xmm2
    \\ addsd %xmm0, %xmm8
    \\ addsd %xmm1, %xmm9
    \\ subsd %xmm3, %xmm8
    \\ addsd %xmm2, %xmm9
    \\ movsd %xmm8, 0 * 8(%r9)
    \\ movsd %xmm9, 1 * 8(%r9)
    \\ addq %r8, %rdx
    \\ addq %r10, %r9
    \\ decq %rax
    \\ jg 20004b
    \\ jmp 20010f
    \\
;
