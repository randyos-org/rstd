# gemm kernel notes

How the Zig kernels of `kernel.zig` compare with the hardware, and where converted assembly would pay. Measured 2026-10-15. The kernels' design is in [`design.md`](../../../design.md#level-3).

## Setup

- **Host:** AMD Ryzen 7 7730U (Zen 3), a laptop, one thread, Windows. Other builds were running, so every figure is the median of 5 to 9 rounds with the configurations interleaved inside each round; ranges are min to max.
- **Build:** a scratch benchmark over `kernel.zig` and `driver.zig`, `-OReleaseFast` with the flags the per-core object of this host gets (`baseline+avx+avx2+fma+sse3+sse4_1+ssse3`, the Haswell and Zen group).
- **Blocking:** `param.h`'s for Zen outside the Windows ABI (whose `dgemm` Q is 128): `sgemm` P 320, Q 320; `dgemm` P 512, Q 256, R 13824; `sgemm` R as x86 derives it.
- **Not available:** no OpenBLAS binary on this machine, so no direct comparison with the reference's kernels.

## Peak

Zen 3 has two 256-bit FMA pipes (AMD, "Software Optimization Guide for AMD EPYC 7003 Processors", pub. 56665 rev. 3.00, section 2.11), so 32 single or 16 double flops per cycle: 64 / 32 GFLOPS at the 2.0 GHz base clock, 144 / 72 at the 4.5 GHz boost (AMD's 7730U specifications, <https://www.amd.com/en/products/processors/laptop/ryzen/7000-series/amd-ryzen-7-7730u.html>). An inline-assembly loop of 12 independent `vfmadd231ps` (or `pd`) chains measures what this laptop sustains: 84 to 89 median single, 44 to 45 double (ranges 80 to 107 and 37 to 55), about 2.7 GHz. The percentages below are of the measured peak.

## Results

| What                              | GFLOPS (median) | Range     | Of peak  |
|-----------------------------------|-----------------|-----------|----------|
| kernel f32, 16 x 6 (chosen)       | 84 to 86        | 57 to 104 | ~100%    |
| kernel f32, 24 x 4                | 38 to 40        | 31 to 48  | 45%      |
| kernel f32, 8 x 12                | 76 to 77        | 51 to 94  | 88%      |
| kernel f64, 8 x 6 (chosen)        | 39 to 46        | 30 to 51  | 90-100%  |
| kernel f64, 12 x 4                | 19 to 21        | 18 to 24  | 45%      |
| `sgemm` NN, n = 1024 to 2048      | 70 to 83        | 54 to 92  | 80-93%   |
| `dgemm` NN and TN, n = 1024, 2048 | 32 to 39        | 23 to 44  | 75-88%   |
| `sgemm` NN, n = 256, 512          | 55 to 74        | 36 to 80  | 65-85%   |

The kernel figures are the kernel alone on panels resident in L1 and L2 (`mr * 4` rows, `nr * 2` columns, depth 256). Packing alone (the driver with a kernel that does nothing) costs 5 to 11% of a 1024 gemm's time. P and Q swept from 96 to 512 moved `sgemm` and `dgemm` by less than the noise, so the reference's Zen blocking stays.

## Code generation

The hot loop of the full tile, from `-femit-asm`:

- **AVX2 (Haswell, Zen; generic tuning, as the object is built):** the chosen 16 x 6 single and 8 x 6 double tiles are 24 instructions: 2 loads of A, 6 broadcasts of B, 12 FMAs, pointer updates; no spills. The 24 x 4 tile fills all 16 registers, and LLVM's generic tuning then keeps one accumulator on the stack (a load and a store every step, plus three register moves), which halves its speed; with `-mcpu=znver3` the same source is clean and runs at 83, hence the spare register in `geometry`.
- **AVX-512 (Skylake-X, compile only):** `skylake_avx512` gets 512-bit registers (no `prefer-256-bit` split): 32 x 14 single, 16 x 14 double, 28 FMAs, 14 broadcasts, 2 loads, no spills. The broadcasts are separate instructions, not embedded `{1to16}` operands.
- **AArch64 NEON (compile only):** LLVM multiplies by a lane of B (`fmla v.4s, v.4s, v.s[i]`), so B takes `nr / vl` registers. Before that was counted, 28 accumulators (8 x 14) collapsed into 50 register moves and 8 stack accesses per step; the chosen 16 x 6 single and 8 x 6 double are 24 `fmla`, 4 loads, no spills.
- **SSE2 x86_64 baseline and Cortex-A15 NEON (compile only):** 8 x 6 and 8 x 5 single, multiply and add or `vfma`, no spills. ARMv7 NEON and pre-VSX AltiVec `f32` tiles are scalar since: those lanes flush denormals (`util.vectorIeee`).

## Tails

Row tails are fully unrolled tiles, one per pair of widths, as OpenBLAS's edge tiles are: on Skylake-X `f32`, 30 tiles. One `thin` function per column width for the tails narrower than a vector moved `sgemm` and `dgemm` at n = 37, 100, 250 and 1000 within the noise (for example `sgemm` n = 100: 36.8 against 36.7 median), so the simpler shape stays. Stripped (`-OReleaseFast`, COFF), `thin` against unrolled: the `s` and `d` `gemm_kernel`s 30 KB against 35 KB on Skylake-X, 20 KB against 18 KB on the AVX2 group; the whole per-core objects 405 KB against 410 KB, and 390 KB against 387 KB. The tiles move `A` and `C` as arrays: the self-hosted x86_64 backend stores a two-lane `f32` vector as sixteen bytes.

## Recommendation

Keep the generic Zig kernel for the AVX2 cores (Haswell, Zen and the cores they serve): it runs at the measured peak, so converting `sgemm`/`dgemm` `haswell` assembly cannot win in the kernel. What the full gemm loses (7 to 25%) is packing, the `C` update and edge tiles; improve the Zig copies first (whole-vector transposes in `oncopy` and `incopy`), and measure again.

Convert assembly first where Zig cannot express what the core's kernel uses, and measure there before keeping it:

1. The arm64 SVE kernels (A64FX, Neoverse V1, N2, V2): `@Vector` has no scalable vectors, so Zig falls back to NEON;
2. POWER10's MMA kernels, and the AMX and AVX-512 BF16 ones for `sbgemm`, which `@Vector` cannot reach at all;
3. Skylake-X and its successors, only if a machine with AVX-512 shows the Zig tile under 85% of the converted `skylakex` kernel (our own threshold, chosen in 2cf8757, not from a source); its loop looks clean but is unmeasured.

## Design notes

The headers of the level-3 files point here; the facts below were cut from them.

**Driver loops** (`driver.zig`). `C` is split into `R`-wide column blocks, the depth into `Q`-deep slices and the rows into `P`-tall blocks; each slice of `op(B)` is packed once (`sb`), each block of `op(A)` once per slice (`sa`). The same loops serve `gemm`, `symm` and `hemm` (the symmetric operand packed whole from its stored triangle, `symm_k.c`, `zhemm_k.c`) and the rank-k updates (`level3_syrk.c`, `level3_syr2k.c`). Unlike the reference, which hands the kernel a slice's whole `sb` for the later row blocks, the kernel takes the `3 * unroll_n` chunks it was packed in: the same access order, right for any unroll. `beta` (`trmm`'s and `trsm`'s `alpha`) scales each block of `C` just before its first depth slice, as BLIS does on its first `KC` slice: the same multiplies as a pass first, so the same bits, and a zero `beta` never reads `C`. More than one worker goes to `threaded.zig`.

**Packed layout and tile** (`kernel.zig`). Panels: `n / u` of width `u`, then one of each power of two below `u` set in `n % u`; panel `p` at column `c0` of width `w` sits at `m * c0` as `m` groups of `w` (OpenBLAS's `kernel/generic` layout for a power-of-two `u`); a complex element is its two parts, real first. The tile is `mr x nr`: `mv` vectors of `vl` lanes by `nr` columns, the most accumulators the registers hold beside the loads of `A` and a broadcast of `B` (`geometry`); complex lanes hold whole elements with two accumulators each, one per part of `B` (`zgemm_kernel_4x2_haswell.S`'s scheme).

**Blocking rules** (`params.zig`). The table has the unrolls and `p`, `q`, `r`; what `setparam-ref.c` computes is ported: `p` of the older x86 cores follows the L2 size (`Affine`); on x86 `r` comes from the packing buffer and `p` of `s`, `d`, `c`, `z` and 3M rounds up to the unroll, and 3M takes the real `p` and `q` unless the header says; mips64's `dgemm` `r` is 640; half and bfloat16 are set per family; LA464 and the AMD AVX-512 parts override by cache size and threads. `OPENBLAS_BLOCK_FACTOR` (`Config.block_factor`) scales a computed `p` as `driver/others/parameter.c` does, here in dynamic builds too.

**Half gemms** (`half.zig`). `driver.zig`'s loops for one input and one output type; `BGEMM` rounds `C` to bfloat16 after every slice, as upstream's kernel stores it, so its bits follow `Q`. Workers take ranges of columns (rows for a tall `C`) as problems of their own. The scalars are tested by value, where upstream compares bfloat16 bits (`level3.c`, `gemm_beta.c`: a `-0` multiplies, a `beta` of 1 still scales) and its `bgemm` beta truncates ([#236](https://github.com/randyos-org/rstd/issues/236)).

**Small kernels** (`small_x86_64.zig`). Each block sums an element over the depth in order, fused, then stores `alpha*sum` fused with `beta*C` as `vfmadd231p[sd]` does: the generic kernels' bits, but for the dot blocks (`tn`, and `nn`'s last rows), whose lanes sum every `L`-th product and reduce as `REDUCE_4` or `_mm512_reduce_add_p[sd]`. Upstream's `nn` copies those rows to a buffer; here they are read in place. The tests run the same code split into narrower vectors. A zero `beta` takes the `B0` forms, which never read `C`; the generic permit allows nothing (`gemm_small_matrix_permit.c`), as complex does upstream.

**3M** (`driver3m.zig`, `gemm3m/kernel.zig`). Each slice takes three passes, packing a real form of `op(A)` and of `alpha * op(B)` (the `o` copies scale first) and adding their product to both parts of `C` with a sign each. One part comes from a difference of products, so it is only normwise stable (Higham, SIAM J. Matrix Anal. Appl. 13(3), 681-687, 1992, <https://doi.org/10.1137/0613043>; `comparisons.md`, 7). Workers take parts on whole panels.

**Triangles** (`trmm/`, `trsm/`). `trmm`'s copies pack `m` deep from `posX`, `n` lines from `posY`, one on a unit diagonal; on the left a row panel takes depths from (`LN`) or up to (`LT`) its row, on the right the column panel up to (`RN`) or from (`RT`) its column. `trsm`'s panels take the gemm update by what is solved, then `solve` their diagonal block, writing `C` and the packed right side; `LT` and `RN` go forward, `LN` and `RT` backward, and `LR`, `LC`, `RR`, `RC` conjugate. A unit diagonal is flagged by the driver's dummy `alpha` of 1.

**Threads** (`threaded.zig`). The `ic` workers of a `jc` group pack each depth slice of their columns of `op(B)` together into one buffer, between a barrier that frees it and one that fills it (`bli_packm_int.c`); each packs its own blocks of `op(A)`, and row ranges start on row panels.
