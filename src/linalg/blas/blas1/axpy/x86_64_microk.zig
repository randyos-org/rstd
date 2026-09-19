//! Copyright (c) 2014, The OpenBLAS Project, under the BSD 3-clause license
//! reproduced in the OpenBLAS section of the repository `LICENSE` file.
//!
//! Inline-assembly bodies of the x86_64 axpy microkernels, mechanically
//! transcribed from `kernel/x86_64/{s,d,c,z}axpy_microk_*-2.c` in
//! https://github.com/OpenMathLib/OpenBLAS at commit
//! 9724481b599bb75fd260a9e7e991b445ef02a0bd (v0.3.34).
//! Re-derive from that source rather than hand-editing if it goes stale.
//!
//! The C operands %0..%5 became named operands: `i`, `n` (in/out), `x`, `y`,
//! `a` (pointer to alpha) and `m` (pointer to the sign vector `mvec`), and
//! `clobbers` is the C clobber list. The C comments are kept as assembler `#` comments; `#if !defined(CONJ)` regions
//! are resolved into `text` (plain) and `text_conj`. `threshold` is the
//! smallest `n` that takes the second block (the source's `if (n < N)`, or
//! `n <= N` as N+1), where it has two.

pub const daxpy_nehalem = struct {
    pub const threshold: ?comptime_int = null;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\movsd (%[a]), %%xmm0 # alpha
        \\shufpd $0, %%xmm0, %%xmm0
        \\.p2align 4
        \\1:
        \\# "prefetcht0 192(%%2,%%0,8) "
        \\# "prefetcht0 192(%%3,%%0,8) "
        \\movups (%[x],%[i],8), %%xmm12 # 2 * x
        \\movups 16(%[x],%[i],8), %%xmm13 # 2 * x
        \\movups 32(%[x],%[i],8), %%xmm14 # 2 * x
        \\movups 48(%[x],%[i],8), %%xmm15 # 2 * x
        \\movups (%[y],%[i],8), %%xmm8 # 2 * y
        \\movups 16(%[y],%[i],8), %%xmm9 # 2 * y
        \\movups 32(%[y],%[i],8), %%xmm10 # 2 * y
        \\movups 48(%[y],%[i],8), %%xmm11 # 2 * y
        \\mulpd %%xmm0 , %%xmm12 # alpha * x
        \\mulpd %%xmm0 , %%xmm13
        \\mulpd %%xmm0 , %%xmm14
        \\mulpd %%xmm0 , %%xmm15
        \\addpd %%xmm12, %%xmm8 # y += alpha *x
        \\addpd %%xmm13, %%xmm9
        \\addpd %%xmm14, %%xmm10
        \\addpd %%xmm15, %%xmm11
        \\movups %%xmm8 , (%[y],%[i],8)
        \\movups %%xmm9 , 16(%[y],%[i],8)
        \\movups %%xmm10, 32(%[y],%[i],8)
        \\movups %%xmm11, 48(%[y],%[i],8)
        \\addq $8 , %[i]
        \\subq $8 , %[n]
        \\jnz 1b
        ,
    };
    pub const text_conj = [_][]const u8{
        \\movsd (%[a]), %%xmm0 # alpha
        \\shufpd $0, %%xmm0, %%xmm0
        \\.p2align 4
        \\1:
        \\# "prefetcht0 192(%%2,%%0,8) "
        \\# "prefetcht0 192(%%3,%%0,8) "
        \\movups (%[x],%[i],8), %%xmm12 # 2 * x
        \\movups 16(%[x],%[i],8), %%xmm13 # 2 * x
        \\movups 32(%[x],%[i],8), %%xmm14 # 2 * x
        \\movups 48(%[x],%[i],8), %%xmm15 # 2 * x
        \\movups (%[y],%[i],8), %%xmm8 # 2 * y
        \\movups 16(%[y],%[i],8), %%xmm9 # 2 * y
        \\movups 32(%[y],%[i],8), %%xmm10 # 2 * y
        \\movups 48(%[y],%[i],8), %%xmm11 # 2 * y
        \\mulpd %%xmm0 , %%xmm12 # alpha * x
        \\mulpd %%xmm0 , %%xmm13
        \\mulpd %%xmm0 , %%xmm14
        \\mulpd %%xmm0 , %%xmm15
        \\addpd %%xmm12, %%xmm8 # y += alpha *x
        \\addpd %%xmm13, %%xmm9
        \\addpd %%xmm14, %%xmm10
        \\addpd %%xmm15, %%xmm11
        \\movups %%xmm8 , (%[y],%[i],8)
        \\movups %%xmm9 , 16(%[y],%[i],8)
        \\movups %%xmm10, 32(%[y],%[i],8)
        \\movups %%xmm11, 48(%[y],%[i],8)
        \\addq $8 , %[i]
        \\subq $8 , %[n]
        \\jnz 1b
        ,
    };
};

pub const daxpy_sandy = struct {
    pub const threshold: ?comptime_int = null;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vbroadcastsd (%[a]), %%ymm0 # alpha
        \\vmovups (%[y],%[i],8), %%ymm8
        \\vmovups 32(%[y],%[i],8), %%ymm9
        \\vmovups 64(%[y],%[i],8), %%ymm10
        \\vmovups 96(%[y],%[i],8), %%ymm11
        \\vmovups (%[x],%[i],8), %%ymm4
        \\vmovups 32(%[x],%[i],8), %%ymm5
        \\vmovups 64(%[x],%[i],8), %%ymm6
        \\vmovups 96(%[x],%[i],8), %%ymm7
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jz 2f
        \\.p2align 4
        \\1:
        \\vmulpd %%ymm4, %%ymm0, %%ymm4
        \\vaddpd %%ymm8 , %%ymm4, %%ymm12
        \\vmulpd %%ymm5, %%ymm0, %%ymm5
        \\vaddpd %%ymm9 , %%ymm5, %%ymm13
        \\vmulpd %%ymm6, %%ymm0, %%ymm6
        \\vaddpd %%ymm10, %%ymm6, %%ymm14
        \\vmulpd %%ymm7, %%ymm0, %%ymm7
        \\vaddpd %%ymm11, %%ymm7, %%ymm15
        \\vmovups (%[y],%[i],8), %%ymm8
        \\vmovups 32(%[y],%[i],8), %%ymm9
        \\vmovups 64(%[y],%[i],8), %%ymm10
        \\vmovups 96(%[y],%[i],8), %%ymm11
        \\vmovups (%[x],%[i],8), %%ymm4
        \\vmovups 32(%[x],%[i],8), %%ymm5
        \\vmovups 64(%[x],%[i],8), %%ymm6
        \\vmovups 96(%[x],%[i],8), %%ymm7
        \\vmovups %%ymm12, -128(%[y],%[i],8)
        \\vmovups %%ymm13, -96(%[y],%[i],8)
        \\vmovups %%ymm14, -64(%[y],%[i],8)
        \\vmovups %%ymm15, -32(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        \\2:
        \\vmulpd %%ymm4, %%ymm0, %%ymm4
        \\vmulpd %%ymm5, %%ymm0, %%ymm5
        \\vmulpd %%ymm6, %%ymm0, %%ymm6
        \\vmulpd %%ymm7, %%ymm0, %%ymm7
        \\vaddpd %%ymm8 , %%ymm4, %%ymm12
        \\vaddpd %%ymm9 , %%ymm5, %%ymm13
        \\vaddpd %%ymm10, %%ymm6, %%ymm14
        \\vaddpd %%ymm11, %%ymm7, %%ymm15
        \\vmovups %%ymm12, -128(%[y],%[i],8)
        \\vmovups %%ymm13, -96(%[y],%[i],8)
        \\vmovups %%ymm14, -64(%[y],%[i],8)
        \\vmovups %%ymm15, -32(%[y],%[i],8)
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vbroadcastsd (%[a]), %%ymm0 # alpha
        \\vmovups (%[y],%[i],8), %%ymm8
        \\vmovups 32(%[y],%[i],8), %%ymm9
        \\vmovups 64(%[y],%[i],8), %%ymm10
        \\vmovups 96(%[y],%[i],8), %%ymm11
        \\vmovups (%[x],%[i],8), %%ymm4
        \\vmovups 32(%[x],%[i],8), %%ymm5
        \\vmovups 64(%[x],%[i],8), %%ymm6
        \\vmovups 96(%[x],%[i],8), %%ymm7
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jz 2f
        \\.p2align 4
        \\1:
        \\vmulpd %%ymm4, %%ymm0, %%ymm4
        \\vaddpd %%ymm8 , %%ymm4, %%ymm12
        \\vmulpd %%ymm5, %%ymm0, %%ymm5
        \\vaddpd %%ymm9 , %%ymm5, %%ymm13
        \\vmulpd %%ymm6, %%ymm0, %%ymm6
        \\vaddpd %%ymm10, %%ymm6, %%ymm14
        \\vmulpd %%ymm7, %%ymm0, %%ymm7
        \\vaddpd %%ymm11, %%ymm7, %%ymm15
        \\vmovups (%[y],%[i],8), %%ymm8
        \\vmovups 32(%[y],%[i],8), %%ymm9
        \\vmovups 64(%[y],%[i],8), %%ymm10
        \\vmovups 96(%[y],%[i],8), %%ymm11
        \\vmovups (%[x],%[i],8), %%ymm4
        \\vmovups 32(%[x],%[i],8), %%ymm5
        \\vmovups 64(%[x],%[i],8), %%ymm6
        \\vmovups 96(%[x],%[i],8), %%ymm7
        \\vmovups %%ymm12, -128(%[y],%[i],8)
        \\vmovups %%ymm13, -96(%[y],%[i],8)
        \\vmovups %%ymm14, -64(%[y],%[i],8)
        \\vmovups %%ymm15, -32(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        \\2:
        \\vmulpd %%ymm4, %%ymm0, %%ymm4
        \\vmulpd %%ymm5, %%ymm0, %%ymm5
        \\vmulpd %%ymm6, %%ymm0, %%ymm6
        \\vmulpd %%ymm7, %%ymm0, %%ymm7
        \\vaddpd %%ymm8 , %%ymm4, %%ymm12
        \\vaddpd %%ymm9 , %%ymm5, %%ymm13
        \\vaddpd %%ymm10, %%ymm6, %%ymm14
        \\vaddpd %%ymm11, %%ymm7, %%ymm15
        \\vmovups %%ymm12, -128(%[y],%[i],8)
        \\vmovups %%ymm13, -96(%[y],%[i],8)
        \\vmovups %%ymm14, -64(%[y],%[i],8)
        \\vmovups %%ymm15, -32(%[y],%[i],8)
        \\vzeroupper
        ,
    };
};

pub const daxpy_haswell = struct {
    pub const threshold: ?comptime_int = null;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vbroadcastsd (%[a]), %%ymm0 # alpha
        \\.p2align 4
        \\1:
        \\vmovups (%[y],%[i],8), %%ymm12 # 4 * y
        \\vmovups 32(%[y],%[i],8), %%ymm13 # 4 * y
        \\vmovups 64(%[y],%[i],8), %%ymm14 # 4 * y
        \\vmovups 96(%[y],%[i],8), %%ymm15 # 4 * y
        \\vfmadd231pd (%[x],%[i],8), %%ymm0 , %%ymm12 # y += alpha * x
        \\vfmadd231pd 32(%[x],%[i],8), %%ymm0 , %%ymm13 # y += alpha * x
        \\vfmadd231pd 64(%[x],%[i],8), %%ymm0 , %%ymm14 # y += alpha * x
        \\vfmadd231pd 96(%[x],%[i],8), %%ymm0 , %%ymm15 # y += alpha * x
        \\vmovups %%ymm12, (%[y],%[i],8)
        \\vmovups %%ymm13, 32(%[y],%[i],8)
        \\vmovups %%ymm14, 64(%[y],%[i],8)
        \\vmovups %%ymm15, 96(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vbroadcastsd (%[a]), %%ymm0 # alpha
        \\.p2align 4
        \\1:
        \\vmovups (%[y],%[i],8), %%ymm12 # 4 * y
        \\vmovups 32(%[y],%[i],8), %%ymm13 # 4 * y
        \\vmovups 64(%[y],%[i],8), %%ymm14 # 4 * y
        \\vmovups 96(%[y],%[i],8), %%ymm15 # 4 * y
        \\vfmadd231pd (%[x],%[i],8), %%ymm0 , %%ymm12 # y += alpha * x
        \\vfmadd231pd 32(%[x],%[i],8), %%ymm0 , %%ymm13 # y += alpha * x
        \\vfmadd231pd 64(%[x],%[i],8), %%ymm0 , %%ymm14 # y += alpha * x
        \\vfmadd231pd 96(%[x],%[i],8), %%ymm0 , %%ymm15 # y += alpha * x
        \\vmovups %%ymm12, (%[y],%[i],8)
        \\vmovups %%ymm13, 32(%[y],%[i],8)
        \\vmovups %%ymm14, 64(%[y],%[i],8)
        \\vmovups %%ymm15, 96(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
};

pub const daxpy_bulldozer = struct {
    pub const threshold: ?comptime_int = null;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vmovddup (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\prefetcht0 768(%[y],%[i],8)
        \\vmovups (%[x],%[i],8), %%xmm12 # 2 * x
        \\vfmaddpd (%[y],%[i],8), %%xmm0 , %%xmm12, %%xmm8 # y += alpha * x
        \\vmovups 16(%[x],%[i],8), %%xmm13 # 2 * x
        \\.align 2
        \\vmovups %%xmm8 , (%[y],%[i],8)
        \\vfmaddpd 16(%[y],%[i],8), %%xmm0 , %%xmm13, %%xmm9 # y += alpha * x
        \\.align 2
        \\vmovups 32(%[x],%[i],8), %%xmm14 # 2 * x
        \\vmovups %%xmm9 , 16(%[y],%[i],8)
        \\prefetcht0 768(%[x],%[i],8)
        \\.align 2
        \\vfmaddpd 32(%[y],%[i],8), %%xmm0 , %%xmm14, %%xmm10 # y += alpha * x
        \\vmovups 48(%[x],%[i],8), %%xmm15 # 2 * x
        \\vmovups %%xmm10, 32(%[y],%[i],8)
        \\vfmaddpd 48(%[y],%[i],8), %%xmm0 , %%xmm15, %%xmm11 # y += alpha * x
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\addq $8 , %[i]
        \\subq $8 , %[n]
        \\jnz 1b
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vmovddup (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\prefetcht0 768(%[y],%[i],8)
        \\vmovups (%[x],%[i],8), %%xmm12 # 2 * x
        \\vfmaddpd (%[y],%[i],8), %%xmm0 , %%xmm12, %%xmm8 # y += alpha * x
        \\vmovups 16(%[x],%[i],8), %%xmm13 # 2 * x
        \\.align 2
        \\vmovups %%xmm8 , (%[y],%[i],8)
        \\vfmaddpd 16(%[y],%[i],8), %%xmm0 , %%xmm13, %%xmm9 # y += alpha * x
        \\.align 2
        \\vmovups 32(%[x],%[i],8), %%xmm14 # 2 * x
        \\vmovups %%xmm9 , 16(%[y],%[i],8)
        \\prefetcht0 768(%[x],%[i],8)
        \\.align 2
        \\vfmaddpd 32(%[y],%[i],8), %%xmm0 , %%xmm14, %%xmm10 # y += alpha * x
        \\vmovups 48(%[x],%[i],8), %%xmm15 # 2 * x
        \\vmovups %%xmm10, 32(%[y],%[i],8)
        \\vfmaddpd 48(%[y],%[i],8), %%xmm0 , %%xmm15, %%xmm11 # y += alpha * x
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\addq $8 , %[i]
        \\subq $8 , %[n]
        \\jnz 1b
        ,
    };
};

pub const daxpy_piledriver = struct {
    pub const threshold: ?comptime_int = 640;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vmovddup (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\vmovups (%[y],%[i],8), %%xmm8 # 2 y
        \\vmovups 16(%[y],%[i],8), %%xmm9 # 2 y
        \\vmovups 32(%[y],%[i],8), %%xmm10 # 2 y
        \\vmovups 48(%[y],%[i],8), %%xmm11 # 2 y
        \\vmovups 64(%[y],%[i],8), %%xmm12 # 2 y
        \\vmovups 80(%[y],%[i],8), %%xmm13 # 2 y
        \\vmovups 96(%[y],%[i],8), %%xmm14 # 2 y
        \\vmovups 112(%[y],%[i],8), %%xmm15 # 2 y
        \\vfmadd231pd (%[x],%[i],8), %%xmm0 , %%xmm8 # y += alpha * x
        \\vfmadd231pd 16(%[x],%[i],8), %%xmm0 , %%xmm9 # y += alpha * x
        \\vfmadd231pd 32(%[x],%[i],8), %%xmm0 , %%xmm10 # y += alpha * x
        \\vfmadd231pd 48(%[x],%[i],8), %%xmm0 , %%xmm11 # y += alpha * x
        \\vfmadd231pd 64(%[x],%[i],8), %%xmm0 , %%xmm12 # y += alpha * x
        \\vfmadd231pd 80(%[x],%[i],8), %%xmm0 , %%xmm13 # y += alpha * x
        \\vfmadd231pd 96(%[x],%[i],8), %%xmm0 , %%xmm14 # y += alpha * x
        \\vfmadd231pd 112(%[x],%[i],8), %%xmm0 , %%xmm15 # y += alpha * x
        \\vmovups %%xmm8 , (%[y],%[i],8)
        \\vmovups %%xmm9 , 16(%[y],%[i],8)
        \\vmovups %%xmm10, 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\vmovups %%xmm12, 64(%[y],%[i],8)
        \\vmovups %%xmm13, 80(%[y],%[i],8)
        \\vmovups %%xmm14, 96(%[y],%[i],8)
        \\vmovups %%xmm15,112(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        ,
        \\vmovddup (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\prefetcht0 512(%[y],%[i],8)
        \\vmovups (%[y],%[i],8), %%xmm8 # 2 y
        \\vmovups 16(%[y],%[i],8), %%xmm9 # 2 y
        \\vmovups 32(%[y],%[i],8), %%xmm10 # 2 y
        \\vmovups 48(%[y],%[i],8), %%xmm11 # 2 y
        \\prefetcht0 576(%[y],%[i],8)
        \\vmovups 64(%[y],%[i],8), %%xmm12 # 2 y
        \\vmovups 80(%[y],%[i],8), %%xmm13 # 2 y
        \\vmovups 96(%[y],%[i],8), %%xmm14 # 2 y
        \\vmovups 112(%[y],%[i],8), %%xmm15 # 2 y
        \\prefetcht0 512(%[x],%[i],8)
        \\vfmadd231pd (%[x],%[i],8), %%xmm0 , %%xmm8 # y += alpha * x
        \\vfmadd231pd 16(%[x],%[i],8), %%xmm0 , %%xmm9 # y += alpha * x
        \\vfmadd231pd 32(%[x],%[i],8), %%xmm0 , %%xmm10 # y += alpha * x
        \\vfmadd231pd 48(%[x],%[i],8), %%xmm0 , %%xmm11 # y += alpha * x
        \\prefetcht0 576(%[x],%[i],8)
        \\vfmadd231pd 64(%[x],%[i],8), %%xmm0 , %%xmm12 # y += alpha * x
        \\vfmadd231pd 80(%[x],%[i],8), %%xmm0 , %%xmm13 # y += alpha * x
        \\vfmadd231pd 96(%[x],%[i],8), %%xmm0 , %%xmm14 # y += alpha * x
        \\vfmadd231pd 112(%[x],%[i],8), %%xmm0 , %%xmm15 # y += alpha * x
        \\vmovups %%xmm8 , (%[y],%[i],8)
        \\vmovups %%xmm9 , 16(%[y],%[i],8)
        \\vmovups %%xmm10, 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\vmovups %%xmm12, 64(%[y],%[i],8)
        \\vmovups %%xmm13, 80(%[y],%[i],8)
        \\vmovups %%xmm14, 96(%[y],%[i],8)
        \\vmovups %%xmm15,112(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vmovddup (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\vmovups (%[y],%[i],8), %%xmm8 # 2 y
        \\vmovups 16(%[y],%[i],8), %%xmm9 # 2 y
        \\vmovups 32(%[y],%[i],8), %%xmm10 # 2 y
        \\vmovups 48(%[y],%[i],8), %%xmm11 # 2 y
        \\vmovups 64(%[y],%[i],8), %%xmm12 # 2 y
        \\vmovups 80(%[y],%[i],8), %%xmm13 # 2 y
        \\vmovups 96(%[y],%[i],8), %%xmm14 # 2 y
        \\vmovups 112(%[y],%[i],8), %%xmm15 # 2 y
        \\vfmadd231pd (%[x],%[i],8), %%xmm0 , %%xmm8 # y += alpha * x
        \\vfmadd231pd 16(%[x],%[i],8), %%xmm0 , %%xmm9 # y += alpha * x
        \\vfmadd231pd 32(%[x],%[i],8), %%xmm0 , %%xmm10 # y += alpha * x
        \\vfmadd231pd 48(%[x],%[i],8), %%xmm0 , %%xmm11 # y += alpha * x
        \\vfmadd231pd 64(%[x],%[i],8), %%xmm0 , %%xmm12 # y += alpha * x
        \\vfmadd231pd 80(%[x],%[i],8), %%xmm0 , %%xmm13 # y += alpha * x
        \\vfmadd231pd 96(%[x],%[i],8), %%xmm0 , %%xmm14 # y += alpha * x
        \\vfmadd231pd 112(%[x],%[i],8), %%xmm0 , %%xmm15 # y += alpha * x
        \\vmovups %%xmm8 , (%[y],%[i],8)
        \\vmovups %%xmm9 , 16(%[y],%[i],8)
        \\vmovups %%xmm10, 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\vmovups %%xmm12, 64(%[y],%[i],8)
        \\vmovups %%xmm13, 80(%[y],%[i],8)
        \\vmovups %%xmm14, 96(%[y],%[i],8)
        \\vmovups %%xmm15,112(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        ,
        \\vmovddup (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\prefetcht0 512(%[y],%[i],8)
        \\vmovups (%[y],%[i],8), %%xmm8 # 2 y
        \\vmovups 16(%[y],%[i],8), %%xmm9 # 2 y
        \\vmovups 32(%[y],%[i],8), %%xmm10 # 2 y
        \\vmovups 48(%[y],%[i],8), %%xmm11 # 2 y
        \\prefetcht0 576(%[y],%[i],8)
        \\vmovups 64(%[y],%[i],8), %%xmm12 # 2 y
        \\vmovups 80(%[y],%[i],8), %%xmm13 # 2 y
        \\vmovups 96(%[y],%[i],8), %%xmm14 # 2 y
        \\vmovups 112(%[y],%[i],8), %%xmm15 # 2 y
        \\prefetcht0 512(%[x],%[i],8)
        \\vfmadd231pd (%[x],%[i],8), %%xmm0 , %%xmm8 # y += alpha * x
        \\vfmadd231pd 16(%[x],%[i],8), %%xmm0 , %%xmm9 # y += alpha * x
        \\vfmadd231pd 32(%[x],%[i],8), %%xmm0 , %%xmm10 # y += alpha * x
        \\vfmadd231pd 48(%[x],%[i],8), %%xmm0 , %%xmm11 # y += alpha * x
        \\prefetcht0 576(%[x],%[i],8)
        \\vfmadd231pd 64(%[x],%[i],8), %%xmm0 , %%xmm12 # y += alpha * x
        \\vfmadd231pd 80(%[x],%[i],8), %%xmm0 , %%xmm13 # y += alpha * x
        \\vfmadd231pd 96(%[x],%[i],8), %%xmm0 , %%xmm14 # y += alpha * x
        \\vfmadd231pd 112(%[x],%[i],8), %%xmm0 , %%xmm15 # y += alpha * x
        \\vmovups %%xmm8 , (%[y],%[i],8)
        \\vmovups %%xmm9 , 16(%[y],%[i],8)
        \\vmovups %%xmm10, 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\vmovups %%xmm12, 64(%[y],%[i],8)
        \\vmovups %%xmm13, 80(%[y],%[i],8)
        \\vmovups %%xmm14, 96(%[y],%[i],8)
        \\vmovups %%xmm15,112(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        ,
    };
};

pub const daxpy_steamroller = struct {
    pub const threshold: ?comptime_int = 2048;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vmovddup (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\vmovups (%[y],%[i],8), %%xmm8 # 2 y
        \\vmovups 16(%[y],%[i],8), %%xmm9 # 2 y
        \\vmovups 32(%[y],%[i],8), %%xmm10 # 2 y
        \\vmovups 48(%[y],%[i],8), %%xmm11 # 2 y
        \\vmovups 64(%[y],%[i],8), %%xmm12 # 2 y
        \\vmovups 80(%[y],%[i],8), %%xmm13 # 2 y
        \\vmovups 96(%[y],%[i],8), %%xmm14 # 2 y
        \\vmovups 112(%[y],%[i],8), %%xmm15 # 2 y
        \\vfmadd231pd (%[x],%[i],8), %%xmm0 , %%xmm8 # y += alpha * x
        \\vfmadd231pd 16(%[x],%[i],8), %%xmm0 , %%xmm9 # y += alpha * x
        \\vfmadd231pd 32(%[x],%[i],8), %%xmm0 , %%xmm10 # y += alpha * x
        \\vfmadd231pd 48(%[x],%[i],8), %%xmm0 , %%xmm11 # y += alpha * x
        \\vfmadd231pd 64(%[x],%[i],8), %%xmm0 , %%xmm12 # y += alpha * x
        \\vfmadd231pd 80(%[x],%[i],8), %%xmm0 , %%xmm13 # y += alpha * x
        \\vfmadd231pd 96(%[x],%[i],8), %%xmm0 , %%xmm14 # y += alpha * x
        \\vfmadd231pd 112(%[x],%[i],8), %%xmm0 , %%xmm15 # y += alpha * x
        \\vmovups %%xmm8 , (%[y],%[i],8)
        \\vmovups %%xmm9 , 16(%[y],%[i],8)
        \\vmovups %%xmm10, 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\vmovups %%xmm12, 64(%[y],%[i],8)
        \\vmovups %%xmm13, 80(%[y],%[i],8)
        \\vmovups %%xmm14, 96(%[y],%[i],8)
        \\vmovups %%xmm15,112(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        ,
        \\vmovddup (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\prefetcht0 512(%[y],%[i],8)
        \\vmovups (%[y],%[i],8), %%xmm8 # 2 y
        \\vmovups 16(%[y],%[i],8), %%xmm9 # 2 y
        \\vmovups 32(%[y],%[i],8), %%xmm10 # 2 y
        \\vmovups 48(%[y],%[i],8), %%xmm11 # 2 y
        \\prefetcht0 576(%[y],%[i],8)
        \\vmovups 64(%[y],%[i],8), %%xmm12 # 2 y
        \\vmovups 80(%[y],%[i],8), %%xmm13 # 2 y
        \\vmovups 96(%[y],%[i],8), %%xmm14 # 2 y
        \\vmovups 112(%[y],%[i],8), %%xmm15 # 2 y
        \\prefetcht0 512(%[x],%[i],8)
        \\vfmadd231pd (%[x],%[i],8), %%xmm0 , %%xmm8 # y += alpha * x
        \\vfmadd231pd 16(%[x],%[i],8), %%xmm0 , %%xmm9 # y += alpha * x
        \\vfmadd231pd 32(%[x],%[i],8), %%xmm0 , %%xmm10 # y += alpha * x
        \\vfmadd231pd 48(%[x],%[i],8), %%xmm0 , %%xmm11 # y += alpha * x
        \\prefetcht0 576(%[x],%[i],8)
        \\vfmadd231pd 64(%[x],%[i],8), %%xmm0 , %%xmm12 # y += alpha * x
        \\vfmadd231pd 80(%[x],%[i],8), %%xmm0 , %%xmm13 # y += alpha * x
        \\vfmadd231pd 96(%[x],%[i],8), %%xmm0 , %%xmm14 # y += alpha * x
        \\vfmadd231pd 112(%[x],%[i],8), %%xmm0 , %%xmm15 # y += alpha * x
        \\vmovups %%xmm8 , (%[y],%[i],8)
        \\vmovups %%xmm9 , 16(%[y],%[i],8)
        \\vmovups %%xmm10, 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\vmovups %%xmm12, 64(%[y],%[i],8)
        \\vmovups %%xmm13, 80(%[y],%[i],8)
        \\vmovups %%xmm14, 96(%[y],%[i],8)
        \\vmovups %%xmm15,112(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vmovddup (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\vmovups (%[y],%[i],8), %%xmm8 # 2 y
        \\vmovups 16(%[y],%[i],8), %%xmm9 # 2 y
        \\vmovups 32(%[y],%[i],8), %%xmm10 # 2 y
        \\vmovups 48(%[y],%[i],8), %%xmm11 # 2 y
        \\vmovups 64(%[y],%[i],8), %%xmm12 # 2 y
        \\vmovups 80(%[y],%[i],8), %%xmm13 # 2 y
        \\vmovups 96(%[y],%[i],8), %%xmm14 # 2 y
        \\vmovups 112(%[y],%[i],8), %%xmm15 # 2 y
        \\vfmadd231pd (%[x],%[i],8), %%xmm0 , %%xmm8 # y += alpha * x
        \\vfmadd231pd 16(%[x],%[i],8), %%xmm0 , %%xmm9 # y += alpha * x
        \\vfmadd231pd 32(%[x],%[i],8), %%xmm0 , %%xmm10 # y += alpha * x
        \\vfmadd231pd 48(%[x],%[i],8), %%xmm0 , %%xmm11 # y += alpha * x
        \\vfmadd231pd 64(%[x],%[i],8), %%xmm0 , %%xmm12 # y += alpha * x
        \\vfmadd231pd 80(%[x],%[i],8), %%xmm0 , %%xmm13 # y += alpha * x
        \\vfmadd231pd 96(%[x],%[i],8), %%xmm0 , %%xmm14 # y += alpha * x
        \\vfmadd231pd 112(%[x],%[i],8), %%xmm0 , %%xmm15 # y += alpha * x
        \\vmovups %%xmm8 , (%[y],%[i],8)
        \\vmovups %%xmm9 , 16(%[y],%[i],8)
        \\vmovups %%xmm10, 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\vmovups %%xmm12, 64(%[y],%[i],8)
        \\vmovups %%xmm13, 80(%[y],%[i],8)
        \\vmovups %%xmm14, 96(%[y],%[i],8)
        \\vmovups %%xmm15,112(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        ,
        \\vmovddup (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\prefetcht0 512(%[y],%[i],8)
        \\vmovups (%[y],%[i],8), %%xmm8 # 2 y
        \\vmovups 16(%[y],%[i],8), %%xmm9 # 2 y
        \\vmovups 32(%[y],%[i],8), %%xmm10 # 2 y
        \\vmovups 48(%[y],%[i],8), %%xmm11 # 2 y
        \\prefetcht0 576(%[y],%[i],8)
        \\vmovups 64(%[y],%[i],8), %%xmm12 # 2 y
        \\vmovups 80(%[y],%[i],8), %%xmm13 # 2 y
        \\vmovups 96(%[y],%[i],8), %%xmm14 # 2 y
        \\vmovups 112(%[y],%[i],8), %%xmm15 # 2 y
        \\prefetcht0 512(%[x],%[i],8)
        \\vfmadd231pd (%[x],%[i],8), %%xmm0 , %%xmm8 # y += alpha * x
        \\vfmadd231pd 16(%[x],%[i],8), %%xmm0 , %%xmm9 # y += alpha * x
        \\vfmadd231pd 32(%[x],%[i],8), %%xmm0 , %%xmm10 # y += alpha * x
        \\vfmadd231pd 48(%[x],%[i],8), %%xmm0 , %%xmm11 # y += alpha * x
        \\prefetcht0 576(%[x],%[i],8)
        \\vfmadd231pd 64(%[x],%[i],8), %%xmm0 , %%xmm12 # y += alpha * x
        \\vfmadd231pd 80(%[x],%[i],8), %%xmm0 , %%xmm13 # y += alpha * x
        \\vfmadd231pd 96(%[x],%[i],8), %%xmm0 , %%xmm14 # y += alpha * x
        \\vfmadd231pd 112(%[x],%[i],8), %%xmm0 , %%xmm15 # y += alpha * x
        \\vmovups %%xmm8 , (%[y],%[i],8)
        \\vmovups %%xmm9 , 16(%[y],%[i],8)
        \\vmovups %%xmm10, 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\vmovups %%xmm12, 64(%[y],%[i],8)
        \\vmovups %%xmm13, 80(%[y],%[i],8)
        \\vmovups %%xmm14, 96(%[y],%[i],8)
        \\vmovups %%xmm15,112(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        ,
    };
};

pub const saxpy_nehalem = struct {
    pub const threshold: ?comptime_int = null;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\movss (%[a]), %%xmm0 # alpha
        \\shufps $0, %%xmm0, %%xmm0
        \\.p2align 4
        \\1:
        \\# "prefetcht0 192(%%2,%%0,4) "
        \\# "prefetcht0 192(%%3,%%0,4) "
        \\movups (%[x],%[i],4), %%xmm12 # 4 * x
        \\movups 16(%[x],%[i],4), %%xmm13 # 4 * x
        \\movups 32(%[x],%[i],4), %%xmm14 # 4 * x
        \\movups 48(%[x],%[i],4), %%xmm15 # 4 * x
        \\movups (%[y],%[i],4), %%xmm8 # 4 * y
        \\movups 16(%[y],%[i],4), %%xmm9 # 4 * y
        \\movups 32(%[y],%[i],4), %%xmm10 # 4 * y
        \\movups 48(%[y],%[i],4), %%xmm11 # 4 * y
        \\mulps %%xmm0 , %%xmm12 # alpha * x
        \\mulps %%xmm0 , %%xmm13
        \\mulps %%xmm0 , %%xmm14
        \\mulps %%xmm0 , %%xmm15
        \\addps %%xmm12, %%xmm8 # y += alpha *x
        \\addps %%xmm13, %%xmm9
        \\addps %%xmm14, %%xmm10
        \\addps %%xmm15, %%xmm11
        \\movups %%xmm8 , (%[y],%[i],4)
        \\movups %%xmm9 , 16(%[y],%[i],4)
        \\movups %%xmm10, 32(%[y],%[i],4)
        \\movups %%xmm11, 48(%[y],%[i],4)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        ,
    };
    pub const text_conj = [_][]const u8{
        \\movss (%[a]), %%xmm0 # alpha
        \\shufps $0, %%xmm0, %%xmm0
        \\.p2align 4
        \\1:
        \\# "prefetcht0 192(%%2,%%0,4) "
        \\# "prefetcht0 192(%%3,%%0,4) "
        \\movups (%[x],%[i],4), %%xmm12 # 4 * x
        \\movups 16(%[x],%[i],4), %%xmm13 # 4 * x
        \\movups 32(%[x],%[i],4), %%xmm14 # 4 * x
        \\movups 48(%[x],%[i],4), %%xmm15 # 4 * x
        \\movups (%[y],%[i],4), %%xmm8 # 4 * y
        \\movups 16(%[y],%[i],4), %%xmm9 # 4 * y
        \\movups 32(%[y],%[i],4), %%xmm10 # 4 * y
        \\movups 48(%[y],%[i],4), %%xmm11 # 4 * y
        \\mulps %%xmm0 , %%xmm12 # alpha * x
        \\mulps %%xmm0 , %%xmm13
        \\mulps %%xmm0 , %%xmm14
        \\mulps %%xmm0 , %%xmm15
        \\addps %%xmm12, %%xmm8 # y += alpha *x
        \\addps %%xmm13, %%xmm9
        \\addps %%xmm14, %%xmm10
        \\addps %%xmm15, %%xmm11
        \\movups %%xmm8 , (%[y],%[i],4)
        \\movups %%xmm9 , 16(%[y],%[i],4)
        \\movups %%xmm10, 32(%[y],%[i],4)
        \\movups %%xmm11, 48(%[y],%[i],4)
        \\addq $16, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        ,
    };
};

pub const saxpy_sandy = struct {
    pub const threshold: ?comptime_int = null;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vbroadcastss (%[a]), %%ymm0 # alpha
        \\vmovups (%[y],%[i],4), %%ymm8
        \\vmovups 32(%[y],%[i],4), %%ymm9
        \\vmovups 64(%[y],%[i],4), %%ymm10
        \\vmovups 96(%[y],%[i],4), %%ymm11
        \\vmovups (%[x],%[i],4), %%ymm4
        \\vmovups 32(%[x],%[i],4), %%ymm5
        \\vmovups 64(%[x],%[i],4), %%ymm6
        \\vmovups 96(%[x],%[i],4), %%ymm7
        \\addq $32, %[i]
        \\subq $32, %[n]
        \\jz 2f
        \\.p2align 4
        \\1:
        \\vmulps %%ymm4, %%ymm0, %%ymm4
        \\vaddps %%ymm8 , %%ymm4, %%ymm12
        \\vmulps %%ymm5, %%ymm0, %%ymm5
        \\vaddps %%ymm9 , %%ymm5, %%ymm13
        \\vmulps %%ymm6, %%ymm0, %%ymm6
        \\vaddps %%ymm10, %%ymm6, %%ymm14
        \\vmulps %%ymm7, %%ymm0, %%ymm7
        \\vaddps %%ymm11, %%ymm7, %%ymm15
        \\vmovups (%[y],%[i],4), %%ymm8
        \\vmovups 32(%[y],%[i],4), %%ymm9
        \\vmovups 64(%[y],%[i],4), %%ymm10
        \\vmovups 96(%[y],%[i],4), %%ymm11
        \\vmovups (%[x],%[i],4), %%ymm4
        \\vmovups 32(%[x],%[i],4), %%ymm5
        \\vmovups 64(%[x],%[i],4), %%ymm6
        \\vmovups 96(%[x],%[i],4), %%ymm7
        \\vmovups %%ymm12, -128(%[y],%[i],4)
        \\vmovups %%ymm13, -96(%[y],%[i],4)
        \\vmovups %%ymm14, -64(%[y],%[i],4)
        \\vmovups %%ymm15, -32(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $32, %[n]
        \\jnz 1b
        \\2:
        \\vmulps %%ymm4, %%ymm0, %%ymm4
        \\vmulps %%ymm5, %%ymm0, %%ymm5
        \\vmulps %%ymm6, %%ymm0, %%ymm6
        \\vmulps %%ymm7, %%ymm0, %%ymm7
        \\vaddps %%ymm8 , %%ymm4, %%ymm12
        \\vaddps %%ymm9 , %%ymm5, %%ymm13
        \\vaddps %%ymm10, %%ymm6, %%ymm14
        \\vaddps %%ymm11, %%ymm7, %%ymm15
        \\vmovups %%ymm12, -128(%[y],%[i],4)
        \\vmovups %%ymm13, -96(%[y],%[i],4)
        \\vmovups %%ymm14, -64(%[y],%[i],4)
        \\vmovups %%ymm15, -32(%[y],%[i],4)
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vbroadcastss (%[a]), %%ymm0 # alpha
        \\vmovups (%[y],%[i],4), %%ymm8
        \\vmovups 32(%[y],%[i],4), %%ymm9
        \\vmovups 64(%[y],%[i],4), %%ymm10
        \\vmovups 96(%[y],%[i],4), %%ymm11
        \\vmovups (%[x],%[i],4), %%ymm4
        \\vmovups 32(%[x],%[i],4), %%ymm5
        \\vmovups 64(%[x],%[i],4), %%ymm6
        \\vmovups 96(%[x],%[i],4), %%ymm7
        \\addq $32, %[i]
        \\subq $32, %[n]
        \\jz 2f
        \\.p2align 4
        \\1:
        \\vmulps %%ymm4, %%ymm0, %%ymm4
        \\vaddps %%ymm8 , %%ymm4, %%ymm12
        \\vmulps %%ymm5, %%ymm0, %%ymm5
        \\vaddps %%ymm9 , %%ymm5, %%ymm13
        \\vmulps %%ymm6, %%ymm0, %%ymm6
        \\vaddps %%ymm10, %%ymm6, %%ymm14
        \\vmulps %%ymm7, %%ymm0, %%ymm7
        \\vaddps %%ymm11, %%ymm7, %%ymm15
        \\vmovups (%[y],%[i],4), %%ymm8
        \\vmovups 32(%[y],%[i],4), %%ymm9
        \\vmovups 64(%[y],%[i],4), %%ymm10
        \\vmovups 96(%[y],%[i],4), %%ymm11
        \\vmovups (%[x],%[i],4), %%ymm4
        \\vmovups 32(%[x],%[i],4), %%ymm5
        \\vmovups 64(%[x],%[i],4), %%ymm6
        \\vmovups 96(%[x],%[i],4), %%ymm7
        \\vmovups %%ymm12, -128(%[y],%[i],4)
        \\vmovups %%ymm13, -96(%[y],%[i],4)
        \\vmovups %%ymm14, -64(%[y],%[i],4)
        \\vmovups %%ymm15, -32(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $32, %[n]
        \\jnz 1b
        \\2:
        \\vmulps %%ymm4, %%ymm0, %%ymm4
        \\vmulps %%ymm5, %%ymm0, %%ymm5
        \\vmulps %%ymm6, %%ymm0, %%ymm6
        \\vmulps %%ymm7, %%ymm0, %%ymm7
        \\vaddps %%ymm8 , %%ymm4, %%ymm12
        \\vaddps %%ymm9 , %%ymm5, %%ymm13
        \\vaddps %%ymm10, %%ymm6, %%ymm14
        \\vaddps %%ymm11, %%ymm7, %%ymm15
        \\vmovups %%ymm12, -128(%[y],%[i],4)
        \\vmovups %%ymm13, -96(%[y],%[i],4)
        \\vmovups %%ymm14, -64(%[y],%[i],4)
        \\vmovups %%ymm15, -32(%[y],%[i],4)
        \\vzeroupper
        ,
    };
};

pub const saxpy_haswell = struct {
    pub const threshold: ?comptime_int = null;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vbroadcastss (%[a]), %%ymm0 # alpha
        \\.p2align 4
        \\1:
        \\vmovups (%[y],%[i],4), %%ymm12 # 8 * y
        \\vmovups 32(%[y],%[i],4), %%ymm13 # 8 * y
        \\vmovups 64(%[y],%[i],4), %%ymm14 # 8 * y
        \\vmovups 96(%[y],%[i],4), %%ymm15 # 8 * y
        \\vfmadd231ps (%[x],%[i],4), %%ymm0 , %%ymm12 # y += alpha * x
        \\vfmadd231ps 32(%[x],%[i],4), %%ymm0 , %%ymm13 # y += alpha * x
        \\vfmadd231ps 64(%[x],%[i],4), %%ymm0 , %%ymm14 # y += alpha * x
        \\vfmadd231ps 96(%[x],%[i],4), %%ymm0 , %%ymm15 # y += alpha * x
        \\vmovups %%ymm12, (%[y],%[i],4)
        \\vmovups %%ymm13, 32(%[y],%[i],4)
        \\vmovups %%ymm14, 64(%[y],%[i],4)
        \\vmovups %%ymm15, 96(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $32, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vbroadcastss (%[a]), %%ymm0 # alpha
        \\.p2align 4
        \\1:
        \\vmovups (%[y],%[i],4), %%ymm12 # 8 * y
        \\vmovups 32(%[y],%[i],4), %%ymm13 # 8 * y
        \\vmovups 64(%[y],%[i],4), %%ymm14 # 8 * y
        \\vmovups 96(%[y],%[i],4), %%ymm15 # 8 * y
        \\vfmadd231ps (%[x],%[i],4), %%ymm0 , %%ymm12 # y += alpha * x
        \\vfmadd231ps 32(%[x],%[i],4), %%ymm0 , %%ymm13 # y += alpha * x
        \\vfmadd231ps 64(%[x],%[i],4), %%ymm0 , %%ymm14 # y += alpha * x
        \\vfmadd231ps 96(%[x],%[i],4), %%ymm0 , %%ymm15 # y += alpha * x
        \\vmovups %%ymm12, (%[y],%[i],4)
        \\vmovups %%ymm13, 32(%[y],%[i],4)
        \\vmovups %%ymm14, 64(%[y],%[i],4)
        \\vmovups %%ymm15, 96(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $32, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
};

pub const saxpy_piledriver = struct {
    pub const threshold: ?comptime_int = 1024;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vzeroupper
        \\vbroadcastss (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\vmovups (%[y],%[i],4), %%xmm8 # 4 * y
        \\vmovups 16(%[y],%[i],4), %%xmm9 # 4 * y
        \\vmovups 32(%[y],%[i],4), %%xmm10 # 4 * y
        \\vmovups 48(%[y],%[i],4), %%xmm11 # 4 * y
        \\vmovups 64(%[y],%[i],4), %%xmm12 # 4 * y
        \\vmovups 80(%[y],%[i],4), %%xmm13 # 4 * y
        \\vmovups 96(%[y],%[i],4), %%xmm14 # 4 * y
        \\vmovups 112(%[y],%[i],4), %%xmm15 # 4 * y
        \\vfmadd231ps (%[x],%[i],4), %%xmm0 , %%xmm8 # y += alpha * x
        \\vfmadd231ps 16(%[x],%[i],4), %%xmm0 , %%xmm9 # y += alpha * x
        \\vfmadd231ps 32(%[x],%[i],4), %%xmm0 , %%xmm10 # y += alpha * x
        \\vfmadd231ps 48(%[x],%[i],4), %%xmm0 , %%xmm11 # y += alpha * x
        \\vfmadd231ps 64(%[x],%[i],4), %%xmm0 , %%xmm12 # y += alpha * x
        \\vfmadd231ps 80(%[x],%[i],4), %%xmm0 , %%xmm13 # y += alpha * x
        \\vfmadd231ps 96(%[x],%[i],4), %%xmm0 , %%xmm14 # y += alpha * x
        \\vfmadd231ps 112(%[x],%[i],4), %%xmm0 , %%xmm15 # y += alpha * x
        \\vmovups %%xmm8 , (%[y],%[i],4)
        \\vmovups %%xmm9 , 16(%[y],%[i],4)
        \\vmovups %%xmm10, 32(%[y],%[i],4)
        \\vmovups %%xmm11, 48(%[y],%[i],4)
        \\vmovups %%xmm12, 64(%[y],%[i],4)
        \\vmovups %%xmm13, 80(%[y],%[i],4)
        \\vmovups %%xmm14, 96(%[y],%[i],4)
        \\vmovups %%xmm15,112(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $32, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
        \\vzeroupper
        \\vbroadcastss (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\prefetcht0 512(%[y],%[i],4)
        \\prefetcht0 576(%[y],%[i],4)
        \\vmovups (%[y],%[i],4), %%xmm8 # 4 * y
        \\vmovups 16(%[y],%[i],4), %%xmm9 # 4 * y
        \\vmovups 32(%[y],%[i],4), %%xmm10 # 4 * y
        \\vmovups 48(%[y],%[i],4), %%xmm11 # 4 * y
        \\vmovups 64(%[y],%[i],4), %%xmm12 # 4 * y
        \\vmovups 80(%[y],%[i],4), %%xmm13 # 4 * y
        \\vmovups 96(%[y],%[i],4), %%xmm14 # 4 * y
        \\vmovups 112(%[y],%[i],4), %%xmm15 # 4 * y
        \\prefetcht0 512(%[x],%[i],4)
        \\prefetcht0 576(%[x],%[i],4)
        \\vfmadd231ps (%[x],%[i],4), %%xmm0 , %%xmm8 # y += alpha * x
        \\vfmadd231ps 16(%[x],%[i],4), %%xmm0 , %%xmm9 # y += alpha * x
        \\vfmadd231ps 32(%[x],%[i],4), %%xmm0 , %%xmm10 # y += alpha * x
        \\vfmadd231ps 48(%[x],%[i],4), %%xmm0 , %%xmm11 # y += alpha * x
        \\vfmadd231ps 64(%[x],%[i],4), %%xmm0 , %%xmm12 # y += alpha * x
        \\vfmadd231ps 80(%[x],%[i],4), %%xmm0 , %%xmm13 # y += alpha * x
        \\vfmadd231ps 96(%[x],%[i],4), %%xmm0 , %%xmm14 # y += alpha * x
        \\vfmadd231ps 112(%[x],%[i],4), %%xmm0 , %%xmm15 # y += alpha * x
        \\vmovups %%xmm8 , (%[y],%[i],4)
        \\vmovups %%xmm9 , 16(%[y],%[i],4)
        \\vmovups %%xmm10, 32(%[y],%[i],4)
        \\vmovups %%xmm11, 48(%[y],%[i],4)
        \\vmovups %%xmm12, 64(%[y],%[i],4)
        \\vmovups %%xmm13, 80(%[y],%[i],4)
        \\vmovups %%xmm14, 96(%[y],%[i],4)
        \\vmovups %%xmm15,112(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $32, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vzeroupper
        \\vbroadcastss (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\vmovups (%[y],%[i],4), %%xmm8 # 4 * y
        \\vmovups 16(%[y],%[i],4), %%xmm9 # 4 * y
        \\vmovups 32(%[y],%[i],4), %%xmm10 # 4 * y
        \\vmovups 48(%[y],%[i],4), %%xmm11 # 4 * y
        \\vmovups 64(%[y],%[i],4), %%xmm12 # 4 * y
        \\vmovups 80(%[y],%[i],4), %%xmm13 # 4 * y
        \\vmovups 96(%[y],%[i],4), %%xmm14 # 4 * y
        \\vmovups 112(%[y],%[i],4), %%xmm15 # 4 * y
        \\vfmadd231ps (%[x],%[i],4), %%xmm0 , %%xmm8 # y += alpha * x
        \\vfmadd231ps 16(%[x],%[i],4), %%xmm0 , %%xmm9 # y += alpha * x
        \\vfmadd231ps 32(%[x],%[i],4), %%xmm0 , %%xmm10 # y += alpha * x
        \\vfmadd231ps 48(%[x],%[i],4), %%xmm0 , %%xmm11 # y += alpha * x
        \\vfmadd231ps 64(%[x],%[i],4), %%xmm0 , %%xmm12 # y += alpha * x
        \\vfmadd231ps 80(%[x],%[i],4), %%xmm0 , %%xmm13 # y += alpha * x
        \\vfmadd231ps 96(%[x],%[i],4), %%xmm0 , %%xmm14 # y += alpha * x
        \\vfmadd231ps 112(%[x],%[i],4), %%xmm0 , %%xmm15 # y += alpha * x
        \\vmovups %%xmm8 , (%[y],%[i],4)
        \\vmovups %%xmm9 , 16(%[y],%[i],4)
        \\vmovups %%xmm10, 32(%[y],%[i],4)
        \\vmovups %%xmm11, 48(%[y],%[i],4)
        \\vmovups %%xmm12, 64(%[y],%[i],4)
        \\vmovups %%xmm13, 80(%[y],%[i],4)
        \\vmovups %%xmm14, 96(%[y],%[i],4)
        \\vmovups %%xmm15,112(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $32, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
        \\vzeroupper
        \\vbroadcastss (%[a]), %%xmm0 # alpha
        \\.align 16
        \\1:
        \\prefetcht0 512(%[y],%[i],4)
        \\prefetcht0 576(%[y],%[i],4)
        \\vmovups (%[y],%[i],4), %%xmm8 # 4 * y
        \\vmovups 16(%[y],%[i],4), %%xmm9 # 4 * y
        \\vmovups 32(%[y],%[i],4), %%xmm10 # 4 * y
        \\vmovups 48(%[y],%[i],4), %%xmm11 # 4 * y
        \\vmovups 64(%[y],%[i],4), %%xmm12 # 4 * y
        \\vmovups 80(%[y],%[i],4), %%xmm13 # 4 * y
        \\vmovups 96(%[y],%[i],4), %%xmm14 # 4 * y
        \\vmovups 112(%[y],%[i],4), %%xmm15 # 4 * y
        \\prefetcht0 512(%[x],%[i],4)
        \\prefetcht0 576(%[x],%[i],4)
        \\vfmadd231ps (%[x],%[i],4), %%xmm0 , %%xmm8 # y += alpha * x
        \\vfmadd231ps 16(%[x],%[i],4), %%xmm0 , %%xmm9 # y += alpha * x
        \\vfmadd231ps 32(%[x],%[i],4), %%xmm0 , %%xmm10 # y += alpha * x
        \\vfmadd231ps 48(%[x],%[i],4), %%xmm0 , %%xmm11 # y += alpha * x
        \\vfmadd231ps 64(%[x],%[i],4), %%xmm0 , %%xmm12 # y += alpha * x
        \\vfmadd231ps 80(%[x],%[i],4), %%xmm0 , %%xmm13 # y += alpha * x
        \\vfmadd231ps 96(%[x],%[i],4), %%xmm0 , %%xmm14 # y += alpha * x
        \\vfmadd231ps 112(%[x],%[i],4), %%xmm0 , %%xmm15 # y += alpha * x
        \\vmovups %%xmm8 , (%[y],%[i],4)
        \\vmovups %%xmm9 , 16(%[y],%[i],4)
        \\vmovups %%xmm10, 32(%[y],%[i],4)
        \\vmovups %%xmm11, 48(%[y],%[i],4)
        \\vmovups %%xmm12, 64(%[y],%[i],4)
        \\vmovups %%xmm13, 80(%[y],%[i],4)
        \\vmovups %%xmm14, 96(%[y],%[i],4)
        \\vmovups %%xmm15,112(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $32, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
};

pub const caxpy_sandy = struct {
    pub const threshold: ?comptime_int = null;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vzeroupper
        \\vbroadcastss (%[a]), %%ymm0 # real part of alpha
        \\vbroadcastss 4(%[a]), %%ymm1 # imag part of alpha
        \\vmulps (%[m]), %%ymm1 , %%ymm1
        \\.p2align 4
        \\1:
        \\vmovups (%[x],%[i],4), %%ymm5 # 4 complex values from x
        \\.p2align 1
        \\vmovups 32(%[x],%[i],4), %%ymm7 # 4 complex values from x
        \\vmovups 64(%[x],%[i],4), %%ymm9 # 4 complex values from x
        \\vmovups 96(%[x],%[i],4), %%ymm11 # 4 complex values from x
        \\vpermilps $0xb1 , %%ymm5 , %%ymm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm7 , %%ymm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm9 , %%ymm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm11, %%ymm10 # exchange real and imag part
        \\vmulps %%ymm5 , %%ymm0 , %%ymm5
        \\vmulps %%ymm7 , %%ymm0 , %%ymm7
        \\vmulps %%ymm9 , %%ymm0 , %%ymm9
        \\vmulps %%ymm11, %%ymm0 , %%ymm11
        \\vaddps (%[y],%[i],4), %%ymm5 , %%ymm5
        \\vaddps 32(%[y],%[i],4), %%ymm7 , %%ymm7
        \\vaddps 64(%[y],%[i],4), %%ymm9 , %%ymm9
        \\vaddps 96(%[y],%[i],4), %%ymm11, %%ymm11
        \\vmulps %%ymm4 , %%ymm1 , %%ymm4
        \\vmulps %%ymm6 , %%ymm1 , %%ymm6
        \\vmulps %%ymm8 , %%ymm1 , %%ymm8
        \\vmulps %%ymm10, %%ymm1 , %%ymm10
        \\vaddps %%ymm4 , %%ymm5 , %%ymm5
        \\vaddps %%ymm6 , %%ymm7 , %%ymm7
        \\vaddps %%ymm8 , %%ymm9 , %%ymm9
        \\vaddps %%ymm10, %%ymm11, %%ymm11
        \\vmovups %%ymm5 , (%[y],%[i],4)
        \\.p2align 1
        \\vmovups %%ymm7 , 32(%[y],%[i],4)
        \\vmovups %%ymm9 , 64(%[y],%[i],4)
        \\vmovups %%ymm11, 96(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vzeroupper
        \\vbroadcastss (%[a]), %%ymm0 # real part of alpha
        \\vbroadcastss 4(%[a]), %%ymm1 # imag part of alpha
        \\vmulps (%[m]), %%ymm0 , %%ymm0
        \\.p2align 4
        \\1:
        \\vmovups (%[x],%[i],4), %%ymm5 # 4 complex values from x
        \\.p2align 1
        \\vmovups 32(%[x],%[i],4), %%ymm7 # 4 complex values from x
        \\vmovups 64(%[x],%[i],4), %%ymm9 # 4 complex values from x
        \\vmovups 96(%[x],%[i],4), %%ymm11 # 4 complex values from x
        \\vpermilps $0xb1 , %%ymm5 , %%ymm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm7 , %%ymm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm9 , %%ymm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm11, %%ymm10 # exchange real and imag part
        \\vmulps %%ymm5 , %%ymm0 , %%ymm5
        \\vmulps %%ymm7 , %%ymm0 , %%ymm7
        \\vmulps %%ymm9 , %%ymm0 , %%ymm9
        \\vmulps %%ymm11, %%ymm0 , %%ymm11
        \\vaddps (%[y],%[i],4), %%ymm5 , %%ymm5
        \\vaddps 32(%[y],%[i],4), %%ymm7 , %%ymm7
        \\vaddps 64(%[y],%[i],4), %%ymm9 , %%ymm9
        \\vaddps 96(%[y],%[i],4), %%ymm11, %%ymm11
        \\vmulps %%ymm4 , %%ymm1 , %%ymm4
        \\vmulps %%ymm6 , %%ymm1 , %%ymm6
        \\vmulps %%ymm8 , %%ymm1 , %%ymm8
        \\vmulps %%ymm10, %%ymm1 , %%ymm10
        \\vaddps %%ymm4 , %%ymm5 , %%ymm5
        \\vaddps %%ymm6 , %%ymm7 , %%ymm7
        \\vaddps %%ymm8 , %%ymm9 , %%ymm9
        \\vaddps %%ymm10, %%ymm11, %%ymm11
        \\vmovups %%ymm5 , (%[y],%[i],4)
        \\.p2align 1
        \\vmovups %%ymm7 , 32(%[y],%[i],4)
        \\vmovups %%ymm9 , 64(%[y],%[i],4)
        \\vmovups %%ymm11, 96(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const mvec = [_]f64{ -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0 };
    pub const mvec_conj = [_]f64{ 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0 };
};

pub const caxpy_haswell = struct {
    pub const threshold: ?comptime_int = null;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vzeroupper
        \\vbroadcastss (%[a]), %%ymm0 # real part of alpha
        \\vbroadcastss 4(%[a]), %%ymm1 # imag part of alpha
        \\vmulps (%[m]), %%ymm1 , %%ymm1
        \\.p2align 4
        \\1:
        \\vmovups (%[x],%[i],4), %%ymm5 # 4 complex values from x
        \\.p2align 1
        \\vmovups 32(%[x],%[i],4), %%ymm7 # 4 complex values from x
        \\vmovups 64(%[x],%[i],4), %%ymm9 # 4 complex values from x
        \\vmovups 96(%[x],%[i],4), %%ymm11 # 4 complex values from x
        \\vmovups 128(%[x],%[i],4), %%ymm12 # 4 complex values from x
        \\vmovups 160(%[x],%[i],4), %%ymm13 # 4 complex values from x
        \\vmovups 192(%[x],%[i],4), %%ymm14 # 4 complex values from x
        \\vmovups 224(%[x],%[i],4), %%ymm15 # 4 complex values from x
        \\vpermilps $0xb1 , %%ymm5 , %%ymm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm7 , %%ymm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm9 , %%ymm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm11, %%ymm10 # exchange real and imag part
        \\vfmadd213ps (%[y],%[i],4), %%ymm0 , %%ymm5
        \\.p2align 1
        \\vfmadd213ps 32(%[y],%[i],4), %%ymm0 , %%ymm7
        \\vfmadd213ps 64(%[y],%[i],4), %%ymm0 , %%ymm9
        \\vfmadd213ps 96(%[y],%[i],4), %%ymm0 , %%ymm11
        \\vfmadd231ps %%ymm1 , %%ymm4 , %%ymm5
        \\vfmadd231ps %%ymm1 , %%ymm6 , %%ymm7
        \\vfmadd231ps %%ymm1 , %%ymm8 , %%ymm9
        \\vfmadd231ps %%ymm1 , %%ymm10, %%ymm11
        \\vpermilps $0xb1 , %%ymm12, %%ymm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm13, %%ymm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm14, %%ymm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm15, %%ymm10 # exchange real and imag part
        \\vfmadd213ps 128(%[y],%[i],4), %%ymm0 , %%ymm12
        \\vfmadd213ps 160(%[y],%[i],4), %%ymm0 , %%ymm13
        \\vfmadd213ps 192(%[y],%[i],4), %%ymm0 , %%ymm14
        \\vfmadd213ps 224(%[y],%[i],4), %%ymm0 , %%ymm15
        \\vfmadd231ps %%ymm1 , %%ymm4 , %%ymm12
        \\vfmadd231ps %%ymm1 , %%ymm6 , %%ymm13
        \\vfmadd231ps %%ymm1 , %%ymm8 , %%ymm14
        \\vfmadd231ps %%ymm1 , %%ymm10, %%ymm15
        \\vmovups %%ymm5 , (%[y],%[i],4)
        \\.p2align 1
        \\vmovups %%ymm7 , 32(%[y],%[i],4)
        \\vmovups %%ymm9 , 64(%[y],%[i],4)
        \\vmovups %%ymm11, 96(%[y],%[i],4)
        \\vmovups %%ymm12,128(%[y],%[i],4)
        \\vmovups %%ymm13,160(%[y],%[i],4)
        \\vmovups %%ymm14,192(%[y],%[i],4)
        \\vmovups %%ymm15,224(%[y],%[i],4)
        \\addq $64, %[i]
        \\subq $32, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vzeroupper
        \\vbroadcastss (%[a]), %%ymm0 # real part of alpha
        \\vbroadcastss 4(%[a]), %%ymm1 # imag part of alpha
        \\vmulps (%[m]), %%ymm0 , %%ymm0
        \\.p2align 4
        \\1:
        \\vmovups (%[x],%[i],4), %%ymm5 # 4 complex values from x
        \\.p2align 1
        \\vmovups 32(%[x],%[i],4), %%ymm7 # 4 complex values from x
        \\vmovups 64(%[x],%[i],4), %%ymm9 # 4 complex values from x
        \\vmovups 96(%[x],%[i],4), %%ymm11 # 4 complex values from x
        \\vmovups 128(%[x],%[i],4), %%ymm12 # 4 complex values from x
        \\vmovups 160(%[x],%[i],4), %%ymm13 # 4 complex values from x
        \\vmovups 192(%[x],%[i],4), %%ymm14 # 4 complex values from x
        \\vmovups 224(%[x],%[i],4), %%ymm15 # 4 complex values from x
        \\vpermilps $0xb1 , %%ymm5 , %%ymm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm7 , %%ymm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm9 , %%ymm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm11, %%ymm10 # exchange real and imag part
        \\vfmadd213ps (%[y],%[i],4), %%ymm0 , %%ymm5
        \\.p2align 1
        \\vfmadd213ps 32(%[y],%[i],4), %%ymm0 , %%ymm7
        \\vfmadd213ps 64(%[y],%[i],4), %%ymm0 , %%ymm9
        \\vfmadd213ps 96(%[y],%[i],4), %%ymm0 , %%ymm11
        \\vfmadd231ps %%ymm1 , %%ymm4 , %%ymm5
        \\vfmadd231ps %%ymm1 , %%ymm6 , %%ymm7
        \\vfmadd231ps %%ymm1 , %%ymm8 , %%ymm9
        \\vfmadd231ps %%ymm1 , %%ymm10, %%ymm11
        \\vpermilps $0xb1 , %%ymm12, %%ymm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm13, %%ymm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm14, %%ymm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%ymm15, %%ymm10 # exchange real and imag part
        \\vfmadd213ps 128(%[y],%[i],4), %%ymm0 , %%ymm12
        \\vfmadd213ps 160(%[y],%[i],4), %%ymm0 , %%ymm13
        \\vfmadd213ps 192(%[y],%[i],4), %%ymm0 , %%ymm14
        \\vfmadd213ps 224(%[y],%[i],4), %%ymm0 , %%ymm15
        \\vfmadd231ps %%ymm1 , %%ymm4 , %%ymm12
        \\vfmadd231ps %%ymm1 , %%ymm6 , %%ymm13
        \\vfmadd231ps %%ymm1 , %%ymm8 , %%ymm14
        \\vfmadd231ps %%ymm1 , %%ymm10, %%ymm15
        \\vmovups %%ymm5 , (%[y],%[i],4)
        \\.p2align 1
        \\vmovups %%ymm7 , 32(%[y],%[i],4)
        \\vmovups %%ymm9 , 64(%[y],%[i],4)
        \\vmovups %%ymm11, 96(%[y],%[i],4)
        \\vmovups %%ymm12,128(%[y],%[i],4)
        \\vmovups %%ymm13,160(%[y],%[i],4)
        \\vmovups %%ymm14,192(%[y],%[i],4)
        \\vmovups %%ymm15,224(%[y],%[i],4)
        \\addq $64, %[i]
        \\subq $32, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const mvec = [_]f64{ -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0 };
    pub const mvec_conj = [_]f64{ 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0 };
};

pub const caxpy_bulldozer = struct {
    pub const threshold: ?comptime_int = 640;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vzeroupper
        \\vbroadcastss (%[a]), %%xmm0 # real part of alpha
        \\vbroadcastss 4(%[a]), %%xmm1 # imag part of alpha
        \\vmulps (%[m]), %%xmm1 , %%xmm1
        \\.align 16
        \\1:
        \\vmovups (%[x],%[i],4), %%xmm5 # 2 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],4), %%xmm7 # 2 complex values from x
        \\vmovups 32(%[x],%[i],4), %%xmm9 # 2 complex values from x
        \\vmovups 48(%[x],%[i],4), %%xmm11 # 2 complex values from x
        \\vmovups 64(%[x],%[i],4), %%xmm12 # 2 complex values from x
        \\vmovups 80(%[x],%[i],4), %%xmm13 # 2 complex values from x
        \\vmovups 96(%[x],%[i],4), %%xmm14 # 2 complex values from x
        \\vmovups 112(%[x],%[i],4), %%xmm15 # 2 complex values from x
        \\vpermilps $0xb1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\vfmaddps (%[y],%[i],4), %%xmm0 , %%xmm5, %%xmm5
        \\.align 2
        \\vfmaddps 16(%[y],%[i],4), %%xmm0 , %%xmm7, %%xmm7
        \\vfmaddps 32(%[y],%[i],4), %%xmm0 , %%xmm9, %%xmm9
        \\vfmaddps 48(%[y],%[i],4), %%xmm0 , %%xmm11,%%xmm11
        \\vfmaddps %%xmm5 , %%xmm1 , %%xmm4 , %%xmm5
        \\vfmaddps %%xmm7 , %%xmm1 , %%xmm6 , %%xmm7
        \\vfmaddps %%xmm9 , %%xmm1 , %%xmm8 , %%xmm9
        \\vfmaddps %%xmm11, %%xmm1 , %%xmm10, %%xmm11
        \\vpermilps $0xb1 , %%xmm12, %%xmm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm13, %%xmm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm14, %%xmm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm15, %%xmm10 # exchange real and imag part
        \\vfmaddps 64(%[y],%[i],4), %%xmm0 , %%xmm12, %%xmm12
        \\vfmaddps 80(%[y],%[i],4), %%xmm0 , %%xmm13, %%xmm13
        \\vfmaddps 96(%[y],%[i],4), %%xmm0 , %%xmm14, %%xmm14
        \\vfmaddps 112(%[y],%[i],4), %%xmm0 , %%xmm15, %%xmm15
        \\vfmaddps %%xmm12, %%xmm1 , %%xmm4 , %%xmm12
        \\vfmaddps %%xmm13, %%xmm1 , %%xmm6 , %%xmm13
        \\vfmaddps %%xmm14, %%xmm1 , %%xmm8 , %%xmm14
        \\vfmaddps %%xmm15, %%xmm1 , %%xmm10, %%xmm15
        \\vmovups %%xmm5 , (%[y],%[i],4)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],4)
        \\vmovups %%xmm9 , 32(%[y],%[i],4)
        \\vmovups %%xmm11, 48(%[y],%[i],4)
        \\vmovups %%xmm12, 64(%[y],%[i],4)
        \\vmovups %%xmm13, 80(%[y],%[i],4)
        \\vmovups %%xmm14, 96(%[y],%[i],4)
        \\vmovups %%xmm15,112(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
        \\vzeroupper
        \\vbroadcastss (%[a]), %%xmm0 # real part of alpha
        \\vbroadcastss 4(%[a]), %%xmm1 # imag part of alpha
        \\vmulps (%[m]), %%xmm1 , %%xmm1
        \\.align 16
        \\1:
        \\prefetcht0 512(%[x],%[i],4)
        \\vmovups (%[x],%[i],4), %%xmm5 # 2 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],4), %%xmm7 # 2 complex values from x
        \\vmovups 32(%[x],%[i],4), %%xmm9 # 2 complex values from x
        \\vmovups 48(%[x],%[i],4), %%xmm11 # 2 complex values from x
        \\vpermilps $0xb1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\prefetcht0 512(%[y],%[i],4)
        \\vfmaddps (%[y],%[i],4), %%xmm0 , %%xmm5, %%xmm5
        \\.align 2
        \\vfmaddps 16(%[y],%[i],4), %%xmm0 , %%xmm7, %%xmm7
        \\vfmaddps 32(%[y],%[i],4), %%xmm0 , %%xmm9, %%xmm9
        \\vfmaddps 48(%[y],%[i],4), %%xmm0 , %%xmm11,%%xmm11
        \\vfmaddps %%xmm5 , %%xmm1 , %%xmm4 , %%xmm5
        \\vfmaddps %%xmm7 , %%xmm1 , %%xmm6 , %%xmm7
        \\vfmaddps %%xmm9 , %%xmm1 , %%xmm8 , %%xmm9
        \\vfmaddps %%xmm11, %%xmm1 , %%xmm10, %%xmm11
        \\vmovups %%xmm5 , (%[y],%[i],4)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],4)
        \\vmovups %%xmm9 , 32(%[y],%[i],4)
        \\vmovups %%xmm11, 48(%[y],%[i],4)
        \\addq $16, %[i]
        \\subq $8, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vzeroupper
        \\vbroadcastss (%[a]), %%xmm0 # real part of alpha
        \\vbroadcastss 4(%[a]), %%xmm1 # imag part of alpha
        \\vmulps (%[m]), %%xmm0 , %%xmm0
        \\.align 16
        \\1:
        \\vmovups (%[x],%[i],4), %%xmm5 # 2 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],4), %%xmm7 # 2 complex values from x
        \\vmovups 32(%[x],%[i],4), %%xmm9 # 2 complex values from x
        \\vmovups 48(%[x],%[i],4), %%xmm11 # 2 complex values from x
        \\vmovups 64(%[x],%[i],4), %%xmm12 # 2 complex values from x
        \\vmovups 80(%[x],%[i],4), %%xmm13 # 2 complex values from x
        \\vmovups 96(%[x],%[i],4), %%xmm14 # 2 complex values from x
        \\vmovups 112(%[x],%[i],4), %%xmm15 # 2 complex values from x
        \\vpermilps $0xb1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\vfmaddps (%[y],%[i],4), %%xmm0 , %%xmm5, %%xmm5
        \\.align 2
        \\vfmaddps 16(%[y],%[i],4), %%xmm0 , %%xmm7, %%xmm7
        \\vfmaddps 32(%[y],%[i],4), %%xmm0 , %%xmm9, %%xmm9
        \\vfmaddps 48(%[y],%[i],4), %%xmm0 , %%xmm11,%%xmm11
        \\vfmaddps %%xmm5 , %%xmm1 , %%xmm4 , %%xmm5
        \\vfmaddps %%xmm7 , %%xmm1 , %%xmm6 , %%xmm7
        \\vfmaddps %%xmm9 , %%xmm1 , %%xmm8 , %%xmm9
        \\vfmaddps %%xmm11, %%xmm1 , %%xmm10, %%xmm11
        \\vpermilps $0xb1 , %%xmm12, %%xmm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm13, %%xmm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm14, %%xmm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm15, %%xmm10 # exchange real and imag part
        \\vfmaddps 64(%[y],%[i],4), %%xmm0 , %%xmm12, %%xmm12
        \\vfmaddps 80(%[y],%[i],4), %%xmm0 , %%xmm13, %%xmm13
        \\vfmaddps 96(%[y],%[i],4), %%xmm0 , %%xmm14, %%xmm14
        \\vfmaddps 112(%[y],%[i],4), %%xmm0 , %%xmm15, %%xmm15
        \\vfmaddps %%xmm12, %%xmm1 , %%xmm4 , %%xmm12
        \\vfmaddps %%xmm13, %%xmm1 , %%xmm6 , %%xmm13
        \\vfmaddps %%xmm14, %%xmm1 , %%xmm8 , %%xmm14
        \\vfmaddps %%xmm15, %%xmm1 , %%xmm10, %%xmm15
        \\vmovups %%xmm5 , (%[y],%[i],4)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],4)
        \\vmovups %%xmm9 , 32(%[y],%[i],4)
        \\vmovups %%xmm11, 48(%[y],%[i],4)
        \\vmovups %%xmm12, 64(%[y],%[i],4)
        \\vmovups %%xmm13, 80(%[y],%[i],4)
        \\vmovups %%xmm14, 96(%[y],%[i],4)
        \\vmovups %%xmm15,112(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
        \\vzeroupper
        \\vbroadcastss (%[a]), %%xmm0 # real part of alpha
        \\vbroadcastss 4(%[a]), %%xmm1 # imag part of alpha
        \\vmulps (%[m]), %%xmm0 , %%xmm0
        \\.align 16
        \\1:
        \\prefetcht0 512(%[x],%[i],4)
        \\vmovups (%[x],%[i],4), %%xmm5 # 2 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],4), %%xmm7 # 2 complex values from x
        \\vmovups 32(%[x],%[i],4), %%xmm9 # 2 complex values from x
        \\vmovups 48(%[x],%[i],4), %%xmm11 # 2 complex values from x
        \\vpermilps $0xb1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\prefetcht0 512(%[y],%[i],4)
        \\vfmaddps (%[y],%[i],4), %%xmm0 , %%xmm5, %%xmm5
        \\.align 2
        \\vfmaddps 16(%[y],%[i],4), %%xmm0 , %%xmm7, %%xmm7
        \\vfmaddps 32(%[y],%[i],4), %%xmm0 , %%xmm9, %%xmm9
        \\vfmaddps 48(%[y],%[i],4), %%xmm0 , %%xmm11,%%xmm11
        \\vfmaddps %%xmm5 , %%xmm1 , %%xmm4 , %%xmm5
        \\vfmaddps %%xmm7 , %%xmm1 , %%xmm6 , %%xmm7
        \\vfmaddps %%xmm9 , %%xmm1 , %%xmm8 , %%xmm9
        \\vfmaddps %%xmm11, %%xmm1 , %%xmm10, %%xmm11
        \\vmovups %%xmm5 , (%[y],%[i],4)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],4)
        \\vmovups %%xmm9 , 32(%[y],%[i],4)
        \\vmovups %%xmm11, 48(%[y],%[i],4)
        \\addq $16, %[i]
        \\subq $8, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const mvec = [_]f64{ -1.0, 1.0, -1.0, 1.0 };
    pub const mvec_conj = [_]f64{ 1.0, -1.0, 1.0, -1.0 };
};

pub const caxpy_steamroller = struct {
    pub const threshold: ?comptime_int = 2049;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vzeroupper
        \\vbroadcastss (%[a]), %%xmm0 # real part of alpha
        \\vbroadcastss 4(%[a]), %%xmm1 # imag part of alpha
        \\vmulps (%[m]), %%xmm1 , %%xmm1
        \\.align 16
        \\1:
        \\vmovups (%[x],%[i],4), %%xmm5 # 2 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],4), %%xmm7 # 2 complex values from x
        \\vmovups 32(%[x],%[i],4), %%xmm9 # 2 complex values from x
        \\vmovups 48(%[x],%[i],4), %%xmm11 # 2 complex values from x
        \\vmovups 64(%[x],%[i],4), %%xmm12 # 2 complex values from x
        \\vmovups 80(%[x],%[i],4), %%xmm13 # 2 complex values from x
        \\vmovups 96(%[x],%[i],4), %%xmm14 # 2 complex values from x
        \\vmovups 112(%[x],%[i],4), %%xmm15 # 2 complex values from x
        \\vpermilps $0xb1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\vfmadd213ps (%[y],%[i],4), %%xmm0 , %%xmm5
        \\.align 2
        \\vfmadd213ps 16(%[y],%[i],4), %%xmm0 , %%xmm7
        \\vfmadd213ps 32(%[y],%[i],4), %%xmm0 , %%xmm9
        \\vfmadd213ps 48(%[y],%[i],4), %%xmm0 , %%xmm11
        \\vfmadd231ps %%xmm1 , %%xmm4 , %%xmm5
        \\vfmadd231ps %%xmm1 , %%xmm6 , %%xmm7
        \\vfmadd231ps %%xmm1 , %%xmm8 , %%xmm9
        \\vfmadd231ps %%xmm1 , %%xmm10, %%xmm11
        \\vpermilps $0xb1 , %%xmm12, %%xmm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm13, %%xmm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm14, %%xmm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm15, %%xmm10 # exchange real and imag part
        \\vfmadd213ps 64(%[y],%[i],4), %%xmm0 , %%xmm12
        \\vfmadd213ps 80(%[y],%[i],4), %%xmm0 , %%xmm13
        \\vfmadd213ps 96(%[y],%[i],4), %%xmm0 , %%xmm14
        \\vfmadd213ps 112(%[y],%[i],4), %%xmm0 , %%xmm15
        \\vfmadd231ps %%xmm1 , %%xmm4 , %%xmm12
        \\vfmadd231ps %%xmm1 , %%xmm6 , %%xmm13
        \\vfmadd231ps %%xmm1 , %%xmm8 , %%xmm14
        \\vfmadd231ps %%xmm1 , %%xmm10, %%xmm15
        \\vmovups %%xmm5 , (%[y],%[i],4)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],4)
        \\vmovups %%xmm9 , 32(%[y],%[i],4)
        \\vmovups %%xmm11, 48(%[y],%[i],4)
        \\vmovups %%xmm12, 64(%[y],%[i],4)
        \\vmovups %%xmm13, 80(%[y],%[i],4)
        \\vmovups %%xmm14, 96(%[y],%[i],4)
        \\vmovups %%xmm15,112(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
        \\vzeroupper
        \\vbroadcastss (%[a]), %%xmm0 # real part of alpha
        \\vbroadcastss 4(%[a]), %%xmm1 # imag part of alpha
        \\vmulps (%[m]), %%xmm1 , %%xmm1
        \\.align 16
        \\1:
        \\prefetcht0 512(%[x],%[i],4)
        \\vmovups (%[x],%[i],4), %%xmm5 # 2 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],4), %%xmm7 # 2 complex values from x
        \\vmovups 32(%[x],%[i],4), %%xmm9 # 2 complex values from x
        \\vmovups 48(%[x],%[i],4), %%xmm11 # 2 complex values from x
        \\vpermilps $0xb1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\prefetcht0 512(%[y],%[i],4)
        \\vfmadd213ps (%[y],%[i],4), %%xmm0 , %%xmm5
        \\.align 2
        \\vfmadd213ps 16(%[y],%[i],4), %%xmm0 , %%xmm7
        \\vfmadd213ps 32(%[y],%[i],4), %%xmm0 , %%xmm9
        \\vfmadd213ps 48(%[y],%[i],4), %%xmm0 , %%xmm11
        \\vfmadd231ps %%xmm1 , %%xmm4 , %%xmm5
        \\vfmadd231ps %%xmm1 , %%xmm6 , %%xmm7
        \\vfmadd231ps %%xmm1 , %%xmm8 , %%xmm9
        \\vfmadd231ps %%xmm1 , %%xmm10, %%xmm11
        \\vmovups %%xmm5 , (%[y],%[i],4)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],4)
        \\vmovups %%xmm9 , 32(%[y],%[i],4)
        \\vmovups %%xmm11, 48(%[y],%[i],4)
        \\addq $16, %[i]
        \\subq $8, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vzeroupper
        \\vbroadcastss (%[a]), %%xmm0 # real part of alpha
        \\vbroadcastss 4(%[a]), %%xmm1 # imag part of alpha
        \\vmulps (%[m]), %%xmm0 , %%xmm0
        \\.align 16
        \\1:
        \\vmovups (%[x],%[i],4), %%xmm5 # 2 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],4), %%xmm7 # 2 complex values from x
        \\vmovups 32(%[x],%[i],4), %%xmm9 # 2 complex values from x
        \\vmovups 48(%[x],%[i],4), %%xmm11 # 2 complex values from x
        \\vmovups 64(%[x],%[i],4), %%xmm12 # 2 complex values from x
        \\vmovups 80(%[x],%[i],4), %%xmm13 # 2 complex values from x
        \\vmovups 96(%[x],%[i],4), %%xmm14 # 2 complex values from x
        \\vmovups 112(%[x],%[i],4), %%xmm15 # 2 complex values from x
        \\vpermilps $0xb1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\vfmadd213ps (%[y],%[i],4), %%xmm0 , %%xmm5
        \\.align 2
        \\vfmadd213ps 16(%[y],%[i],4), %%xmm0 , %%xmm7
        \\vfmadd213ps 32(%[y],%[i],4), %%xmm0 , %%xmm9
        \\vfmadd213ps 48(%[y],%[i],4), %%xmm0 , %%xmm11
        \\vfmadd231ps %%xmm1 , %%xmm4 , %%xmm5
        \\vfmadd231ps %%xmm1 , %%xmm6 , %%xmm7
        \\vfmadd231ps %%xmm1 , %%xmm8 , %%xmm9
        \\vfmadd231ps %%xmm1 , %%xmm10, %%xmm11
        \\vpermilps $0xb1 , %%xmm12, %%xmm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm13, %%xmm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm14, %%xmm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm15, %%xmm10 # exchange real and imag part
        \\vfmadd213ps 64(%[y],%[i],4), %%xmm0 , %%xmm12
        \\vfmadd213ps 80(%[y],%[i],4), %%xmm0 , %%xmm13
        \\vfmadd213ps 96(%[y],%[i],4), %%xmm0 , %%xmm14
        \\vfmadd213ps 112(%[y],%[i],4), %%xmm0 , %%xmm15
        \\vfmadd231ps %%xmm1 , %%xmm4 , %%xmm12
        \\vfmadd231ps %%xmm1 , %%xmm6 , %%xmm13
        \\vfmadd231ps %%xmm1 , %%xmm8 , %%xmm14
        \\vfmadd231ps %%xmm1 , %%xmm10, %%xmm15
        \\vmovups %%xmm5 , (%[y],%[i],4)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],4)
        \\vmovups %%xmm9 , 32(%[y],%[i],4)
        \\vmovups %%xmm11, 48(%[y],%[i],4)
        \\vmovups %%xmm12, 64(%[y],%[i],4)
        \\vmovups %%xmm13, 80(%[y],%[i],4)
        \\vmovups %%xmm14, 96(%[y],%[i],4)
        \\vmovups %%xmm15,112(%[y],%[i],4)
        \\addq $32, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
        \\vzeroupper
        \\vbroadcastss (%[a]), %%xmm0 # real part of alpha
        \\vbroadcastss 4(%[a]), %%xmm1 # imag part of alpha
        \\vmulps (%[m]), %%xmm0 , %%xmm0
        \\.align 16
        \\1:
        \\prefetcht0 512(%[x],%[i],4)
        \\vmovups (%[x],%[i],4), %%xmm5 # 2 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],4), %%xmm7 # 2 complex values from x
        \\vmovups 32(%[x],%[i],4), %%xmm9 # 2 complex values from x
        \\vmovups 48(%[x],%[i],4), %%xmm11 # 2 complex values from x
        \\vpermilps $0xb1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilps $0xb1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\prefetcht0 512(%[y],%[i],4)
        \\vfmadd213ps (%[y],%[i],4), %%xmm0 , %%xmm5
        \\.align 2
        \\vfmadd213ps 16(%[y],%[i],4), %%xmm0 , %%xmm7
        \\vfmadd213ps 32(%[y],%[i],4), %%xmm0 , %%xmm9
        \\vfmadd213ps 48(%[y],%[i],4), %%xmm0 , %%xmm11
        \\vfmadd231ps %%xmm1 , %%xmm4 , %%xmm5
        \\vfmadd231ps %%xmm1 , %%xmm6 , %%xmm7
        \\vfmadd231ps %%xmm1 , %%xmm8 , %%xmm9
        \\vfmadd231ps %%xmm1 , %%xmm10, %%xmm11
        \\vmovups %%xmm5 , (%[y],%[i],4)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],4)
        \\vmovups %%xmm9 , 32(%[y],%[i],4)
        \\vmovups %%xmm11, 48(%[y],%[i],4)
        \\addq $16, %[i]
        \\subq $8, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const mvec = [_]f64{ -1.0, 1.0, -1.0, 1.0 };
    pub const mvec_conj = [_]f64{ 1.0, -1.0, 1.0, -1.0 };
};

pub const zaxpy_sandy = struct {
    pub const threshold: ?comptime_int = 1280;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vzeroupper
        \\vbroadcastsd (%[a]), %%ymm0 # real part of alpha
        \\vbroadcastsd 8(%[a]), %%ymm1 # imag part of alpha
        \\vmulpd (%[m]), %%ymm1 , %%ymm1
        \\.p2align 4
        \\1:
        \\vmovups (%[x],%[i],8), %%ymm5 # 4 complex values from x
        \\.p2align 1
        \\vmovups 32(%[x],%[i],8), %%ymm7 # 4 complex values from x
        \\vmovups 64(%[x],%[i],8), %%ymm9 # 4 complex values from x
        \\vmovups 96(%[x],%[i],8), %%ymm11 # 4 complex values from x
        \\vpermilpd $0x5 , %%ymm5 , %%ymm4 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm7 , %%ymm6 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm9 , %%ymm8 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm11, %%ymm10 # exchange real and imag part
        \\vmulpd %%ymm5 , %%ymm0 , %%ymm5
        \\vmulpd %%ymm7 , %%ymm0 , %%ymm7
        \\vmulpd %%ymm9 , %%ymm0 , %%ymm9
        \\vmulpd %%ymm11, %%ymm0 , %%ymm11
        \\vaddpd (%[y],%[i],8), %%ymm5 , %%ymm5
        \\vaddpd 32(%[y],%[i],8), %%ymm7 , %%ymm7
        \\vaddpd 64(%[y],%[i],8), %%ymm9 , %%ymm9
        \\vaddpd 96(%[y],%[i],8), %%ymm11, %%ymm11
        \\vmulpd %%ymm4 , %%ymm1 , %%ymm4
        \\vmulpd %%ymm6 , %%ymm1 , %%ymm6
        \\vmulpd %%ymm8 , %%ymm1 , %%ymm8
        \\vmulpd %%ymm10, %%ymm1 , %%ymm10
        \\vaddpd %%ymm4 , %%ymm5 , %%ymm5
        \\vaddpd %%ymm6 , %%ymm7 , %%ymm7
        \\vaddpd %%ymm8 , %%ymm9 , %%ymm9
        \\vaddpd %%ymm10, %%ymm11, %%ymm11
        \\vmovups %%ymm5 , (%[y],%[i],8)
        \\.p2align 1
        \\vmovups %%ymm7 , 32(%[y],%[i],8)
        \\vmovups %%ymm9 , 64(%[y],%[i],8)
        \\vmovups %%ymm11, 96(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $8 , %[n]
        \\jnz 1b
        \\vzeroupper
        ,
        \\vzeroupper
        \\vbroadcastsd (%[a]), %%ymm0 # real part of alpha
        \\vbroadcastsd 8(%[a]), %%ymm1 # imag part of alpha
        \\vmulpd (%[m]), %%ymm1 , %%ymm1
        \\.p2align 4
        \\1:
        \\prefetcht0 512(%[x],%[i],8)
        \\prefetcht0 576(%[x],%[i],8)
        \\vmovups (%[x],%[i],8), %%ymm5 # 4 complex values from x
        \\.p2align 1
        \\vmovups 32(%[x],%[i],8), %%ymm7 # 4 complex values from x
        \\vmovups 64(%[x],%[i],8), %%ymm9 # 4 complex values from x
        \\vmovups 96(%[x],%[i],8), %%ymm11 # 4 complex values from x
        \\vpermilpd $0x5 , %%ymm5 , %%ymm4 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm7 , %%ymm6 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm9 , %%ymm8 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm11, %%ymm10 # exchange real and imag part
        \\vmulpd %%ymm5 , %%ymm0 , %%ymm5
        \\vmulpd %%ymm7 , %%ymm0 , %%ymm7
        \\vmulpd %%ymm9 , %%ymm0 , %%ymm9
        \\vmulpd %%ymm11, %%ymm0 , %%ymm11
        \\prefetcht0 512(%[y],%[i],8)
        \\prefetcht0 576(%[y],%[i],8)
        \\vaddpd (%[y],%[i],8), %%ymm5 , %%ymm5
        \\vaddpd 32(%[y],%[i],8), %%ymm7 , %%ymm7
        \\vaddpd 64(%[y],%[i],8), %%ymm9 , %%ymm9
        \\vaddpd 96(%[y],%[i],8), %%ymm11, %%ymm11
        \\vmulpd %%ymm4 , %%ymm1 , %%ymm4
        \\vmulpd %%ymm6 , %%ymm1 , %%ymm6
        \\vmulpd %%ymm8 , %%ymm1 , %%ymm8
        \\vmulpd %%ymm10, %%ymm1 , %%ymm10
        \\vaddpd %%ymm4 , %%ymm5 , %%ymm5
        \\vaddpd %%ymm6 , %%ymm7 , %%ymm7
        \\vaddpd %%ymm8 , %%ymm9 , %%ymm9
        \\vaddpd %%ymm10, %%ymm11, %%ymm11
        \\vmovups %%ymm5 , (%[y],%[i],8)
        \\.p2align 1
        \\vmovups %%ymm7 , 32(%[y],%[i],8)
        \\vmovups %%ymm9 , 64(%[y],%[i],8)
        \\vmovups %%ymm11, 96(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $8 , %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vzeroupper
        \\vbroadcastsd (%[a]), %%ymm0 # real part of alpha
        \\vbroadcastsd 8(%[a]), %%ymm1 # imag part of alpha
        \\vmulpd (%[m]), %%ymm0 , %%ymm0
        \\.p2align 4
        \\1:
        \\vmovups (%[x],%[i],8), %%ymm5 # 4 complex values from x
        \\.p2align 1
        \\vmovups 32(%[x],%[i],8), %%ymm7 # 4 complex values from x
        \\vmovups 64(%[x],%[i],8), %%ymm9 # 4 complex values from x
        \\vmovups 96(%[x],%[i],8), %%ymm11 # 4 complex values from x
        \\vpermilpd $0x5 , %%ymm5 , %%ymm4 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm7 , %%ymm6 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm9 , %%ymm8 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm11, %%ymm10 # exchange real and imag part
        \\vmulpd %%ymm5 , %%ymm0 , %%ymm5
        \\vmulpd %%ymm7 , %%ymm0 , %%ymm7
        \\vmulpd %%ymm9 , %%ymm0 , %%ymm9
        \\vmulpd %%ymm11, %%ymm0 , %%ymm11
        \\vaddpd (%[y],%[i],8), %%ymm5 , %%ymm5
        \\vaddpd 32(%[y],%[i],8), %%ymm7 , %%ymm7
        \\vaddpd 64(%[y],%[i],8), %%ymm9 , %%ymm9
        \\vaddpd 96(%[y],%[i],8), %%ymm11, %%ymm11
        \\vmulpd %%ymm4 , %%ymm1 , %%ymm4
        \\vmulpd %%ymm6 , %%ymm1 , %%ymm6
        \\vmulpd %%ymm8 , %%ymm1 , %%ymm8
        \\vmulpd %%ymm10, %%ymm1 , %%ymm10
        \\vaddpd %%ymm4 , %%ymm5 , %%ymm5
        \\vaddpd %%ymm6 , %%ymm7 , %%ymm7
        \\vaddpd %%ymm8 , %%ymm9 , %%ymm9
        \\vaddpd %%ymm10, %%ymm11, %%ymm11
        \\vmovups %%ymm5 , (%[y],%[i],8)
        \\.p2align 1
        \\vmovups %%ymm7 , 32(%[y],%[i],8)
        \\vmovups %%ymm9 , 64(%[y],%[i],8)
        \\vmovups %%ymm11, 96(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $8 , %[n]
        \\jnz 1b
        \\vzeroupper
        ,
        \\vzeroupper
        \\vbroadcastsd (%[a]), %%ymm0 # real part of alpha
        \\vbroadcastsd 8(%[a]), %%ymm1 # imag part of alpha
        \\vmulpd (%[m]), %%ymm0 , %%ymm0
        \\.p2align 4
        \\1:
        \\prefetcht0 512(%[x],%[i],8)
        \\prefetcht0 576(%[x],%[i],8)
        \\vmovups (%[x],%[i],8), %%ymm5 # 4 complex values from x
        \\.p2align 1
        \\vmovups 32(%[x],%[i],8), %%ymm7 # 4 complex values from x
        \\vmovups 64(%[x],%[i],8), %%ymm9 # 4 complex values from x
        \\vmovups 96(%[x],%[i],8), %%ymm11 # 4 complex values from x
        \\vpermilpd $0x5 , %%ymm5 , %%ymm4 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm7 , %%ymm6 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm9 , %%ymm8 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm11, %%ymm10 # exchange real and imag part
        \\vmulpd %%ymm5 , %%ymm0 , %%ymm5
        \\vmulpd %%ymm7 , %%ymm0 , %%ymm7
        \\vmulpd %%ymm9 , %%ymm0 , %%ymm9
        \\vmulpd %%ymm11, %%ymm0 , %%ymm11
        \\prefetcht0 512(%[y],%[i],8)
        \\prefetcht0 576(%[y],%[i],8)
        \\vaddpd (%[y],%[i],8), %%ymm5 , %%ymm5
        \\vaddpd 32(%[y],%[i],8), %%ymm7 , %%ymm7
        \\vaddpd 64(%[y],%[i],8), %%ymm9 , %%ymm9
        \\vaddpd 96(%[y],%[i],8), %%ymm11, %%ymm11
        \\vmulpd %%ymm4 , %%ymm1 , %%ymm4
        \\vmulpd %%ymm6 , %%ymm1 , %%ymm6
        \\vmulpd %%ymm8 , %%ymm1 , %%ymm8
        \\vmulpd %%ymm10, %%ymm1 , %%ymm10
        \\vaddpd %%ymm4 , %%ymm5 , %%ymm5
        \\vaddpd %%ymm6 , %%ymm7 , %%ymm7
        \\vaddpd %%ymm8 , %%ymm9 , %%ymm9
        \\vaddpd %%ymm10, %%ymm11, %%ymm11
        \\vmovups %%ymm5 , (%[y],%[i],8)
        \\.p2align 1
        \\vmovups %%ymm7 , 32(%[y],%[i],8)
        \\vmovups %%ymm9 , 64(%[y],%[i],8)
        \\vmovups %%ymm11, 96(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $8 , %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const mvec = [_]f64{ -1.0, 1.0, -1.0, 1.0 };
    pub const mvec_conj = [_]f64{ 1.0, -1.0, 1.0, -1.0 };
};

pub const zaxpy_haswell = struct {
    pub const threshold: ?comptime_int = null;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vzeroupper
        \\vbroadcastsd (%[a]), %%ymm0 # real part of alpha
        \\vbroadcastsd 8(%[a]), %%ymm1 # imag part of alpha
        \\vmulpd (%[m]), %%ymm1 , %%ymm1
        \\.p2align 4
        \\1:
        \\vmovups (%[x],%[i],8), %%ymm5 # 2 complex values from x
        \\.p2align 1
        \\vmovups 32(%[x],%[i],8), %%ymm7 # 2 complex values from x
        \\vmovups 64(%[x],%[i],8), %%ymm9 # 2 complex values from x
        \\vmovups 96(%[x],%[i],8), %%ymm11 # 2 complex values from x
        \\vmovups 128(%[x],%[i],8), %%ymm12 # 2 complex values from x
        \\vmovups 160(%[x],%[i],8), %%ymm13 # 2 complex values from x
        \\vmovups 192(%[x],%[i],8), %%ymm14 # 2 complex values from x
        \\vmovups 224(%[x],%[i],8), %%ymm15 # 2 complex values from x
        \\vpermilpd $0x5 , %%ymm5 , %%ymm4 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm7 , %%ymm6 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm9 , %%ymm8 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm11, %%ymm10 # exchange real and imag part
        \\vfmadd213pd (%[y],%[i],8), %%ymm0 , %%ymm5
        \\.p2align 1
        \\vfmadd213pd 32(%[y],%[i],8), %%ymm0 , %%ymm7
        \\vfmadd213pd 64(%[y],%[i],8), %%ymm0 , %%ymm9
        \\vfmadd213pd 96(%[y],%[i],8), %%ymm0 , %%ymm11
        \\vfmadd231pd %%ymm1 , %%ymm4 , %%ymm5
        \\vfmadd231pd %%ymm1 , %%ymm6 , %%ymm7
        \\vfmadd231pd %%ymm1 , %%ymm8 , %%ymm9
        \\vfmadd231pd %%ymm1 , %%ymm10, %%ymm11
        \\vpermilpd $0x5 , %%ymm12, %%ymm4 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm13, %%ymm6 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm14, %%ymm8 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm15, %%ymm10 # exchange real and imag part
        \\vfmadd213pd 128(%[y],%[i],8), %%ymm0 , %%ymm12
        \\vfmadd213pd 160(%[y],%[i],8), %%ymm0 , %%ymm13
        \\vfmadd213pd 192(%[y],%[i],8), %%ymm0 , %%ymm14
        \\vfmadd213pd 224(%[y],%[i],8), %%ymm0 , %%ymm15
        \\vfmadd231pd %%ymm1 , %%ymm4 , %%ymm12
        \\vfmadd231pd %%ymm1 , %%ymm6 , %%ymm13
        \\vfmadd231pd %%ymm1 , %%ymm8 , %%ymm14
        \\vfmadd231pd %%ymm1 , %%ymm10, %%ymm15
        \\vmovups %%ymm5 , (%[y],%[i],8)
        \\.p2align 1
        \\vmovups %%ymm7 , 32(%[y],%[i],8)
        \\vmovups %%ymm9 , 64(%[y],%[i],8)
        \\vmovups %%ymm11, 96(%[y],%[i],8)
        \\vmovups %%ymm12,128(%[y],%[i],8)
        \\vmovups %%ymm13,160(%[y],%[i],8)
        \\vmovups %%ymm14,192(%[y],%[i],8)
        \\vmovups %%ymm15,224(%[y],%[i],8)
        \\addq $32, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vzeroupper
        \\vbroadcastsd (%[a]), %%ymm0 # real part of alpha
        \\vbroadcastsd 8(%[a]), %%ymm1 # imag part of alpha
        \\vmulpd (%[m]), %%ymm0 , %%ymm0
        \\.p2align 4
        \\1:
        \\vmovups (%[x],%[i],8), %%ymm5 # 2 complex values from x
        \\.p2align 1
        \\vmovups 32(%[x],%[i],8), %%ymm7 # 2 complex values from x
        \\vmovups 64(%[x],%[i],8), %%ymm9 # 2 complex values from x
        \\vmovups 96(%[x],%[i],8), %%ymm11 # 2 complex values from x
        \\vmovups 128(%[x],%[i],8), %%ymm12 # 2 complex values from x
        \\vmovups 160(%[x],%[i],8), %%ymm13 # 2 complex values from x
        \\vmovups 192(%[x],%[i],8), %%ymm14 # 2 complex values from x
        \\vmovups 224(%[x],%[i],8), %%ymm15 # 2 complex values from x
        \\vpermilpd $0x5 , %%ymm5 , %%ymm4 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm7 , %%ymm6 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm9 , %%ymm8 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm11, %%ymm10 # exchange real and imag part
        \\vfmadd213pd (%[y],%[i],8), %%ymm0 , %%ymm5
        \\.p2align 1
        \\vfmadd213pd 32(%[y],%[i],8), %%ymm0 , %%ymm7
        \\vfmadd213pd 64(%[y],%[i],8), %%ymm0 , %%ymm9
        \\vfmadd213pd 96(%[y],%[i],8), %%ymm0 , %%ymm11
        \\vfmadd231pd %%ymm1 , %%ymm4 , %%ymm5
        \\vfmadd231pd %%ymm1 , %%ymm6 , %%ymm7
        \\vfmadd231pd %%ymm1 , %%ymm8 , %%ymm9
        \\vfmadd231pd %%ymm1 , %%ymm10, %%ymm11
        \\vpermilpd $0x5 , %%ymm12, %%ymm4 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm13, %%ymm6 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm14, %%ymm8 # exchange real and imag part
        \\vpermilpd $0x5 , %%ymm15, %%ymm10 # exchange real and imag part
        \\vfmadd213pd 128(%[y],%[i],8), %%ymm0 , %%ymm12
        \\vfmadd213pd 160(%[y],%[i],8), %%ymm0 , %%ymm13
        \\vfmadd213pd 192(%[y],%[i],8), %%ymm0 , %%ymm14
        \\vfmadd213pd 224(%[y],%[i],8), %%ymm0 , %%ymm15
        \\vfmadd231pd %%ymm1 , %%ymm4 , %%ymm12
        \\vfmadd231pd %%ymm1 , %%ymm6 , %%ymm13
        \\vfmadd231pd %%ymm1 , %%ymm8 , %%ymm14
        \\vfmadd231pd %%ymm1 , %%ymm10, %%ymm15
        \\vmovups %%ymm5 , (%[y],%[i],8)
        \\.p2align 1
        \\vmovups %%ymm7 , 32(%[y],%[i],8)
        \\vmovups %%ymm9 , 64(%[y],%[i],8)
        \\vmovups %%ymm11, 96(%[y],%[i],8)
        \\vmovups %%ymm12,128(%[y],%[i],8)
        \\vmovups %%ymm13,160(%[y],%[i],8)
        \\vmovups %%ymm14,192(%[y],%[i],8)
        \\vmovups %%ymm15,224(%[y],%[i],8)
        \\addq $32, %[i]
        \\subq $16, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const mvec = [_]f64{ -1.0, 1.0, -1.0, 1.0 };
    pub const mvec_conj = [_]f64{ 1.0, -1.0, 1.0, -1.0 };
};

pub const zaxpy_bulldozer = struct {
    pub const threshold: ?comptime_int = 384;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vzeroupper
        \\vmovddup (%[a]), %%xmm0 # real part of alpha
        \\vmovddup 8(%[a]), %%xmm1 # imag part of alpha
        \\vmulpd (%[m]), %%xmm1 , %%xmm1
        \\.align 16
        \\1:
        \\vmovups (%[x],%[i],8), %%xmm5 # 1 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],8), %%xmm7 # 1 complex values from x
        \\vmovups 32(%[x],%[i],8), %%xmm9 # 1 complex values from x
        \\vmovups 48(%[x],%[i],8), %%xmm11 # 1 complex values from x
        \\vmovups 64(%[x],%[i],8), %%xmm12 # 1 complex values from x
        \\vmovups 80(%[x],%[i],8), %%xmm13 # 1 complex values from x
        \\vmovups 96(%[x],%[i],8), %%xmm14 # 1 complex values from x
        \\vmovups 112(%[x],%[i],8), %%xmm15 # 1 complex values from x
        \\vpermilpd $0x1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\vfmaddpd (%[y],%[i],8), %%xmm0 , %%xmm5, %%xmm5
        \\.align 2
        \\vfmaddpd 16(%[y],%[i],8), %%xmm0 , %%xmm7, %%xmm7
        \\vfmaddpd 32(%[y],%[i],8), %%xmm0 , %%xmm9, %%xmm9
        \\vfmaddpd 48(%[y],%[i],8), %%xmm0 , %%xmm11,%%xmm11
        \\vfmaddpd %%xmm5 , %%xmm1 , %%xmm4 , %%xmm5
        \\vfmaddpd %%xmm7 , %%xmm1 , %%xmm6 , %%xmm7
        \\vfmaddpd %%xmm9 , %%xmm1 , %%xmm8 , %%xmm9
        \\vfmaddpd %%xmm11, %%xmm1 , %%xmm10, %%xmm11
        \\vpermilpd $0x1 , %%xmm12, %%xmm4 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm13, %%xmm6 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm14, %%xmm8 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm15, %%xmm10 # exchange real and imag part
        \\vfmaddpd 64(%[y],%[i],8), %%xmm0 , %%xmm12, %%xmm12
        \\vfmaddpd 80(%[y],%[i],8), %%xmm0 , %%xmm13, %%xmm13
        \\vfmaddpd 96(%[y],%[i],8), %%xmm0 , %%xmm14, %%xmm14
        \\vfmaddpd 112(%[y],%[i],8), %%xmm0 , %%xmm15, %%xmm15
        \\vfmaddpd %%xmm12, %%xmm1 , %%xmm4 , %%xmm12
        \\vfmaddpd %%xmm13, %%xmm1 , %%xmm6 , %%xmm13
        \\vfmaddpd %%xmm14, %%xmm1 , %%xmm8 , %%xmm14
        \\vfmaddpd %%xmm15, %%xmm1 , %%xmm10, %%xmm15
        \\vmovups %%xmm5 , (%[y],%[i],8)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],8)
        \\vmovups %%xmm9 , 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\vmovups %%xmm12, 64(%[y],%[i],8)
        \\vmovups %%xmm13, 80(%[y],%[i],8)
        \\vmovups %%xmm14, 96(%[y],%[i],8)
        \\vmovups %%xmm15,112(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $8 , %[n]
        \\jnz 1b
        \\vzeroupper
        ,
        \\vzeroupper
        \\vmovddup (%[a]), %%xmm0 # real part of alpha
        \\vmovddup 8(%[a]), %%xmm1 # imag part of alpha
        \\vmulpd (%[m]), %%xmm1 , %%xmm1
        \\.align 16
        \\1:
        \\prefetcht0 512(%[x],%[i],8)
        \\vmovups (%[x],%[i],8), %%xmm5 # 1 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],8), %%xmm7 # 1 complex values from x
        \\vmovups 32(%[x],%[i],8), %%xmm9 # 1 complex values from x
        \\vmovups 48(%[x],%[i],8), %%xmm11 # 1 complex values from x
        \\vpermilpd $0x1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\prefetcht0 512(%[y],%[i],8)
        \\vfmaddpd (%[y],%[i],8), %%xmm0 , %%xmm5, %%xmm5
        \\.align 2
        \\vfmaddpd 16(%[y],%[i],8), %%xmm0 , %%xmm7, %%xmm7
        \\vfmaddpd 32(%[y],%[i],8), %%xmm0 , %%xmm9, %%xmm9
        \\vfmaddpd 48(%[y],%[i],8), %%xmm0 , %%xmm11,%%xmm11
        \\vfmaddpd %%xmm5 , %%xmm1 , %%xmm4 , %%xmm5
        \\vfmaddpd %%xmm7 , %%xmm1 , %%xmm6 , %%xmm7
        \\vfmaddpd %%xmm9 , %%xmm1 , %%xmm8 , %%xmm9
        \\vfmaddpd %%xmm11, %%xmm1 , %%xmm10, %%xmm11
        \\vmovups %%xmm5 , (%[y],%[i],8)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],8)
        \\vmovups %%xmm9 , 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\addq $8 , %[i]
        \\subq $4, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vzeroupper
        \\vmovddup (%[a]), %%xmm0 # real part of alpha
        \\vmovddup 8(%[a]), %%xmm1 # imag part of alpha
        \\vmulpd (%[m]), %%xmm0 , %%xmm0
        \\.align 16
        \\1:
        \\vmovups (%[x],%[i],8), %%xmm5 # 1 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],8), %%xmm7 # 1 complex values from x
        \\vmovups 32(%[x],%[i],8), %%xmm9 # 1 complex values from x
        \\vmovups 48(%[x],%[i],8), %%xmm11 # 1 complex values from x
        \\vmovups 64(%[x],%[i],8), %%xmm12 # 1 complex values from x
        \\vmovups 80(%[x],%[i],8), %%xmm13 # 1 complex values from x
        \\vmovups 96(%[x],%[i],8), %%xmm14 # 1 complex values from x
        \\vmovups 112(%[x],%[i],8), %%xmm15 # 1 complex values from x
        \\vpermilpd $0x1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\vfmaddpd (%[y],%[i],8), %%xmm0 , %%xmm5, %%xmm5
        \\.align 2
        \\vfmaddpd 16(%[y],%[i],8), %%xmm0 , %%xmm7, %%xmm7
        \\vfmaddpd 32(%[y],%[i],8), %%xmm0 , %%xmm9, %%xmm9
        \\vfmaddpd 48(%[y],%[i],8), %%xmm0 , %%xmm11,%%xmm11
        \\vfmaddpd %%xmm5 , %%xmm1 , %%xmm4 , %%xmm5
        \\vfmaddpd %%xmm7 , %%xmm1 , %%xmm6 , %%xmm7
        \\vfmaddpd %%xmm9 , %%xmm1 , %%xmm8 , %%xmm9
        \\vfmaddpd %%xmm11, %%xmm1 , %%xmm10, %%xmm11
        \\vpermilpd $0x1 , %%xmm12, %%xmm4 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm13, %%xmm6 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm14, %%xmm8 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm15, %%xmm10 # exchange real and imag part
        \\vfmaddpd 64(%[y],%[i],8), %%xmm0 , %%xmm12, %%xmm12
        \\vfmaddpd 80(%[y],%[i],8), %%xmm0 , %%xmm13, %%xmm13
        \\vfmaddpd 96(%[y],%[i],8), %%xmm0 , %%xmm14, %%xmm14
        \\vfmaddpd 112(%[y],%[i],8), %%xmm0 , %%xmm15, %%xmm15
        \\vfmaddpd %%xmm12, %%xmm1 , %%xmm4 , %%xmm12
        \\vfmaddpd %%xmm13, %%xmm1 , %%xmm6 , %%xmm13
        \\vfmaddpd %%xmm14, %%xmm1 , %%xmm8 , %%xmm14
        \\vfmaddpd %%xmm15, %%xmm1 , %%xmm10, %%xmm15
        \\vmovups %%xmm5 , (%[y],%[i],8)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],8)
        \\vmovups %%xmm9 , 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\vmovups %%xmm12, 64(%[y],%[i],8)
        \\vmovups %%xmm13, 80(%[y],%[i],8)
        \\vmovups %%xmm14, 96(%[y],%[i],8)
        \\vmovups %%xmm15,112(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $8 , %[n]
        \\jnz 1b
        \\vzeroupper
        ,
        \\vzeroupper
        \\vmovddup (%[a]), %%xmm0 # real part of alpha
        \\vmovddup 8(%[a]), %%xmm1 # imag part of alpha
        \\vmulpd (%[m]), %%xmm0 , %%xmm0
        \\.align 16
        \\1:
        \\prefetcht0 512(%[x],%[i],8)
        \\vmovups (%[x],%[i],8), %%xmm5 # 1 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],8), %%xmm7 # 1 complex values from x
        \\vmovups 32(%[x],%[i],8), %%xmm9 # 1 complex values from x
        \\vmovups 48(%[x],%[i],8), %%xmm11 # 1 complex values from x
        \\vpermilpd $0x1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\prefetcht0 512(%[y],%[i],8)
        \\vfmaddpd (%[y],%[i],8), %%xmm0 , %%xmm5, %%xmm5
        \\.align 2
        \\vfmaddpd 16(%[y],%[i],8), %%xmm0 , %%xmm7, %%xmm7
        \\vfmaddpd 32(%[y],%[i],8), %%xmm0 , %%xmm9, %%xmm9
        \\vfmaddpd 48(%[y],%[i],8), %%xmm0 , %%xmm11,%%xmm11
        \\vfmaddpd %%xmm5 , %%xmm1 , %%xmm4 , %%xmm5
        \\vfmaddpd %%xmm7 , %%xmm1 , %%xmm6 , %%xmm7
        \\vfmaddpd %%xmm9 , %%xmm1 , %%xmm8 , %%xmm9
        \\vfmaddpd %%xmm11, %%xmm1 , %%xmm10, %%xmm11
        \\vmovups %%xmm5 , (%[y],%[i],8)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],8)
        \\vmovups %%xmm9 , 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\addq $8 , %[i]
        \\subq $4, %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const mvec = [_]f64{ -1.0, 1.0 };
    pub const mvec_conj = [_]f64{ 1.0, -1.0 };
};

pub const zaxpy_steamroller = struct {
    pub const threshold: ?comptime_int = 640;
    pub const clobbers = [_][]const u8{ "cc", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15", "memory" };
    pub const text = [_][]const u8{
        \\vzeroupper
        \\vmovddup (%[a]), %%xmm0 # real part of alpha
        \\vmovddup 8(%[a]), %%xmm1 # imag part of alpha
        \\vmulpd (%[m]), %%xmm1 , %%xmm1
        \\.align 16
        \\1:
        \\vmovups (%[x],%[i],8), %%xmm5 # 2 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],8), %%xmm7 # 2 complex values from x
        \\vmovups 32(%[x],%[i],8), %%xmm9 # 2 complex values from x
        \\vmovups 48(%[x],%[i],8), %%xmm11 # 2 complex values from x
        \\vmovups 64(%[x],%[i],8), %%xmm12 # 2 complex values from x
        \\vmovups 80(%[x],%[i],8), %%xmm13 # 2 complex values from x
        \\vmovups 96(%[x],%[i],8), %%xmm14 # 2 complex values from x
        \\vmovups 112(%[x],%[i],8), %%xmm15 # 2 complex values from x
        \\vpermilpd $0x1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\vfmadd213pd (%[y],%[i],8), %%xmm0 , %%xmm5
        \\.align 2
        \\vfmadd213pd 16(%[y],%[i],8), %%xmm0 , %%xmm7
        \\vfmadd213pd 32(%[y],%[i],8), %%xmm0 , %%xmm9
        \\vfmadd213pd 48(%[y],%[i],8), %%xmm0 , %%xmm11
        \\vfmadd231pd %%xmm1 , %%xmm4 , %%xmm5
        \\vfmadd231pd %%xmm1 , %%xmm6 , %%xmm7
        \\vfmadd231pd %%xmm1 , %%xmm8 , %%xmm9
        \\vfmadd231pd %%xmm1 , %%xmm10, %%xmm11
        \\vpermilpd $0x1 , %%xmm12, %%xmm4 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm13, %%xmm6 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm14, %%xmm8 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm15, %%xmm10 # exchange real and imag part
        \\vfmadd213pd 64(%[y],%[i],8), %%xmm0 , %%xmm12
        \\vfmadd213pd 80(%[y],%[i],8), %%xmm0 , %%xmm13
        \\vfmadd213pd 96(%[y],%[i],8), %%xmm0 , %%xmm14
        \\vfmadd213pd 112(%[y],%[i],8), %%xmm0 , %%xmm15
        \\vfmadd231pd %%xmm1 , %%xmm4 , %%xmm12
        \\vfmadd231pd %%xmm1 , %%xmm6 , %%xmm13
        \\vfmadd231pd %%xmm1 , %%xmm8 , %%xmm14
        \\vfmadd231pd %%xmm1 , %%xmm10, %%xmm15
        \\vmovups %%xmm5 , (%[y],%[i],8)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],8)
        \\vmovups %%xmm9 , 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\vmovups %%xmm12, 64(%[y],%[i],8)
        \\vmovups %%xmm13, 80(%[y],%[i],8)
        \\vmovups %%xmm14, 96(%[y],%[i],8)
        \\vmovups %%xmm15,112(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $8 , %[n]
        \\jnz 1b
        \\vzeroupper
        ,
        \\vzeroupper
        \\vmovddup (%[a]), %%xmm0 # real part of alpha
        \\vmovddup 8(%[a]), %%xmm1 # imag part of alpha
        \\vmulpd (%[m]), %%xmm1 , %%xmm1
        \\.align 16
        \\1:
        \\prefetcht0 512(%[x],%[i],8)
        \\vmovups (%[x],%[i],8), %%xmm5 # 2 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],8), %%xmm7 # 2 complex values from x
        \\vmovups 32(%[x],%[i],8), %%xmm9 # 2 complex values from x
        \\vmovups 48(%[x],%[i],8), %%xmm11 # 2 complex values from x
        \\vpermilpd $0x1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\prefetcht0 512(%[y],%[i],8)
        \\vfmadd213pd (%[y],%[i],8), %%xmm0 , %%xmm5
        \\.align 2
        \\vfmadd213pd 16(%[y],%[i],8), %%xmm0 , %%xmm7
        \\vfmadd213pd 32(%[y],%[i],8), %%xmm0 , %%xmm9
        \\vfmadd213pd 48(%[y],%[i],8), %%xmm0 , %%xmm11
        \\vfmadd231pd %%xmm1 , %%xmm4 , %%xmm5
        \\vfmadd231pd %%xmm1 , %%xmm6 , %%xmm7
        \\vfmadd231pd %%xmm1 , %%xmm8 , %%xmm9
        \\vfmadd231pd %%xmm1 , %%xmm10, %%xmm11
        \\vmovups %%xmm5 , (%[y],%[i],8)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],8)
        \\vmovups %%xmm9 , 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\addq $8 , %[i]
        \\subq $4 , %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const text_conj = [_][]const u8{
        \\vzeroupper
        \\vmovddup (%[a]), %%xmm0 # real part of alpha
        \\vmovddup 8(%[a]), %%xmm1 # imag part of alpha
        \\vmulpd (%[m]), %%xmm0 , %%xmm0
        \\.align 16
        \\1:
        \\vmovups (%[x],%[i],8), %%xmm5 # 2 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],8), %%xmm7 # 2 complex values from x
        \\vmovups 32(%[x],%[i],8), %%xmm9 # 2 complex values from x
        \\vmovups 48(%[x],%[i],8), %%xmm11 # 2 complex values from x
        \\vmovups 64(%[x],%[i],8), %%xmm12 # 2 complex values from x
        \\vmovups 80(%[x],%[i],8), %%xmm13 # 2 complex values from x
        \\vmovups 96(%[x],%[i],8), %%xmm14 # 2 complex values from x
        \\vmovups 112(%[x],%[i],8), %%xmm15 # 2 complex values from x
        \\vpermilpd $0x1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\vfmadd213pd (%[y],%[i],8), %%xmm0 , %%xmm5
        \\.align 2
        \\vfmadd213pd 16(%[y],%[i],8), %%xmm0 , %%xmm7
        \\vfmadd213pd 32(%[y],%[i],8), %%xmm0 , %%xmm9
        \\vfmadd213pd 48(%[y],%[i],8), %%xmm0 , %%xmm11
        \\vfmadd231pd %%xmm1 , %%xmm4 , %%xmm5
        \\vfmadd231pd %%xmm1 , %%xmm6 , %%xmm7
        \\vfmadd231pd %%xmm1 , %%xmm8 , %%xmm9
        \\vfmadd231pd %%xmm1 , %%xmm10, %%xmm11
        \\vpermilpd $0x1 , %%xmm12, %%xmm4 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm13, %%xmm6 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm14, %%xmm8 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm15, %%xmm10 # exchange real and imag part
        \\vfmadd213pd 64(%[y],%[i],8), %%xmm0 , %%xmm12
        \\vfmadd213pd 80(%[y],%[i],8), %%xmm0 , %%xmm13
        \\vfmadd213pd 96(%[y],%[i],8), %%xmm0 , %%xmm14
        \\vfmadd213pd 112(%[y],%[i],8), %%xmm0 , %%xmm15
        \\vfmadd231pd %%xmm1 , %%xmm4 , %%xmm12
        \\vfmadd231pd %%xmm1 , %%xmm6 , %%xmm13
        \\vfmadd231pd %%xmm1 , %%xmm8 , %%xmm14
        \\vfmadd231pd %%xmm1 , %%xmm10, %%xmm15
        \\vmovups %%xmm5 , (%[y],%[i],8)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],8)
        \\vmovups %%xmm9 , 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\vmovups %%xmm12, 64(%[y],%[i],8)
        \\vmovups %%xmm13, 80(%[y],%[i],8)
        \\vmovups %%xmm14, 96(%[y],%[i],8)
        \\vmovups %%xmm15,112(%[y],%[i],8)
        \\addq $16, %[i]
        \\subq $8 , %[n]
        \\jnz 1b
        \\vzeroupper
        ,
        \\vzeroupper
        \\vmovddup (%[a]), %%xmm0 # real part of alpha
        \\vmovddup 8(%[a]), %%xmm1 # imag part of alpha
        \\vmulpd (%[m]), %%xmm0 , %%xmm0
        \\.align 16
        \\1:
        \\prefetcht0 512(%[x],%[i],8)
        \\vmovups (%[x],%[i],8), %%xmm5 # 2 complex values from x
        \\.align 2
        \\vmovups 16(%[x],%[i],8), %%xmm7 # 2 complex values from x
        \\vmovups 32(%[x],%[i],8), %%xmm9 # 2 complex values from x
        \\vmovups 48(%[x],%[i],8), %%xmm11 # 2 complex values from x
        \\vpermilpd $0x1 , %%xmm5 , %%xmm4 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm7 , %%xmm6 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm9 , %%xmm8 # exchange real and imag part
        \\vpermilpd $0x1 , %%xmm11, %%xmm10 # exchange real and imag part
        \\prefetcht0 512(%[y],%[i],8)
        \\vfmadd213pd (%[y],%[i],8), %%xmm0 , %%xmm5
        \\.align 2
        \\vfmadd213pd 16(%[y],%[i],8), %%xmm0 , %%xmm7
        \\vfmadd213pd 32(%[y],%[i],8), %%xmm0 , %%xmm9
        \\vfmadd213pd 48(%[y],%[i],8), %%xmm0 , %%xmm11
        \\vfmadd231pd %%xmm1 , %%xmm4 , %%xmm5
        \\vfmadd231pd %%xmm1 , %%xmm6 , %%xmm7
        \\vfmadd231pd %%xmm1 , %%xmm8 , %%xmm9
        \\vfmadd231pd %%xmm1 , %%xmm10, %%xmm11
        \\vmovups %%xmm5 , (%[y],%[i],8)
        \\.align 2
        \\vmovups %%xmm7 , 16(%[y],%[i],8)
        \\vmovups %%xmm9 , 32(%[y],%[i],8)
        \\vmovups %%xmm11, 48(%[y],%[i],8)
        \\addq $8 , %[i]
        \\subq $4 , %[n]
        \\jnz 1b
        \\vzeroupper
        ,
    };
    pub const mvec = [_]f64{ -1.0, 1.0 };
    pub const mvec_conj = [_]f64{ 1.0, -1.0 };
};
