//! `cgemv_n` block kernel, the assembly as written, the source preprocessed with `XCONJ`, from OpenBLAS `kernel/x86_64/cgemv_n_microk_haswell-4.c` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

/// `add_y` past its C branch `if (inc_dest != 2)`, which is the driver's.
pub noinline fn add_y_x_asm(n: isize, src: [*]f32, dest: [*]f32, inc_dest: isize, alpha_r: f32, alpha_i: f32) void {
    _ = inc_dest;
    var i: isize = undefined;
    i = 0;
    var n1: isize = n & -8;
    const n2: isize = n & 4;
    asm volatile (
        \\vzeroupper
        \\        vbroadcastss      (%[op4]), %%ymm0
        \\        vbroadcastss      (%[op5]), %%ymm1
        \\        cmpq            $0 , %[n1]
        \\        je              2f
        \\        1:
        \\        vmovups (%[src],%[i],4), %%ymm8
        \\        vmovups      32(%[src],%[i],4), %%ymm9
        \\        vmulps      %%ymm8 , %%ymm0, %%ymm12
        \\        vmulps      %%ymm8 , %%ymm1, %%ymm13
        \\        vmulps      %%ymm9 , %%ymm0, %%ymm14
        \\        vmulps      %%ymm9 , %%ymm1, %%ymm15
        \\        vmovups   (%[dest],%[i],4),  %%ymm10
        \\        vmovups 32(%[dest],%[i],4),  %%ymm11
        \\        vpermilps      $0xb1 , %%ymm12, %%ymm12
        \\        vpermilps      $0xb1 , %%ymm14, %%ymm14
        \\        vaddsubps      %%ymm12, %%ymm13, %%ymm8
        \\        vaddsubps      %%ymm14, %%ymm15, %%ymm9
        \\        vpermilps      $0xb1 , %%ymm8 , %%ymm8
        \\        vpermilps      $0xb1 , %%ymm9 , %%ymm9
        \\        addq            $16, %[i]
        \\        vaddps         %%ymm8, %%ymm10, %%ymm12
        \\        vaddps         %%ymm9, %%ymm11, %%ymm13
        \\        subq            $8 , %[n1]
        \\        vmovups  %%ymm12,-64(%[dest],%[i],4)
        \\        vmovups  %%ymm13,-32(%[dest],%[i],4)
        \\        jnz             1b
        \\        2:
        \\        cmpq            $4, %[n2]
        \\        jne             3f
        \\        vmovups (%[src],%[i],4), %%ymm8
        \\        vmulps      %%ymm8 , %%ymm0, %%ymm12
        \\        vmulps      %%ymm8 , %%ymm1, %%ymm13
        \\        vmovups   (%[dest],%[i],4),  %%ymm10
        \\        vpermilps      $0xb1 , %%ymm12, %%ymm12
        \\        vaddsubps      %%ymm12, %%ymm13, %%ymm8
        \\        vpermilps      $0xb1 , %%ymm8 , %%ymm8
        \\        vaddps         %%ymm8, %%ymm10, %%ymm12
        \\        vmovups  %%ymm12,   (%[dest],%[i],4)
        \\        3:
        \\        vzeroupper
        : [i] "+r" (i),
          [n1] "+r" (n1),
        : [src] "r" (src),
          [dest] "r" (dest),
          [op4] "r" (&alpha_r),
          [op5] "r" (&alpha_i),
          [n2] "r" (n2),
        : .{ .cc = true, .ymm0 = true, .ymm1 = true, .xmm2 = true, .xmm3 = true, .xmm4 = true, .xmm5 = true, .xmm6 = true, .xmm7 = true, .ymm8 = true, .ymm9 = true, .ymm10 = true, .ymm11 = true, .ymm12 = true, .ymm13 = true, .ymm14 = true, .ymm15 = true, .memory = true });
}
