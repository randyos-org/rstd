# gemm kernel notes

How the Zig kernels of `kernel.zig` compare with the hardware, and where
converted assembly would pay. Measured 2026-10-15.

## Setup

- **Host:** AMD Ryzen 7 7730U (Zen 3), a laptop, one thread, Windows. Other
  builds were running, so every figure is the median of 5 to 9 rounds with
  the configurations interleaved inside each round; ranges are min to max.
- **Build:** a scratch benchmark over `kernel.zig` and `driver.zig`,
  `-OReleaseFast` with the flags the per-core object of this host gets
  (`baseline+avx+avx2+fma+sse3+sse4_1+ssse3`, the Haswell and Zen group).
- **Blocking:** `param.h`'s for Zen: `sgemm` P 320, Q 320; `dgemm` P 512,
  Q 256; R as x86 derives it.
- **Not available:** no OpenBLAS binary on this machine, so no direct
  comparison with the reference's kernels.

## Peak

Assumption: Zen 3 has two 256-bit FMA pipes, so 32 single or 16 double
flops per cycle: 64 / 32 GFLOPS at the 2.0 GHz base clock, 144 / 72 at the
4.5 GHz boost. An inline-assembly loop of 12 independent `vfmadd231ps` (or
`pd`) chains measures what this laptop sustains: 84 to 89 median single,
44 to 45 double (ranges 80 to 107 and 37 to 55), about 2.7 GHz. The
percentages below are of the measured peak.

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

The kernel figures are the kernel alone on panels resident in L1 and L2
(`mr * 4` rows, `nr * 2` columns, depth 256). Packing alone (the driver with
a kernel that does nothing) costs 5 to 11% of a 1024 gemm's time. P and Q
swept from 96 to 512 moved `sgemm` and `dgemm` by less than the noise, so
the reference's Zen blocking stays.

## Code generation

The hot loop of the full tile, from `-femit-asm`:

- **AVX2 (Haswell, Zen; generic tuning, as the object is built):** the
  chosen 16 x 6 single and 8 x 6 double tiles are 24 instructions: 2 loads
  of A, 6 broadcasts of B, 12 FMAs, pointer updates; no spills. The 24 x 4
  tile fills all 16 registers, and LLVM's generic tuning then keeps one
  accumulator on the stack (a load and a store every step, plus three
  register moves), which halves its speed; with `-mcpu=znver3` the same
  source is clean and runs at 83, hence the spare register in `geometry`.
- **AVX-512 (Skylake-X, compile only):** `skylake_avx512` gets 512-bit
  registers (no `prefer-256-bit` split): 32 x 14 single, 16 x 14 double,
  28 FMAs, 14 broadcasts, 2 loads, no spills. The broadcasts are separate
  instructions, not embedded `{1to16}` operands.
- **AArch64 NEON (compile only):** LLVM multiplies by a lane of B
  (`fmla v.4s, v.4s, v.s[i]`), so B takes `nr / vl` registers. Before that
  was counted, 28 accumulators (8 x 14) collapsed into 50 register moves
  and 8 stack accesses per step; the chosen 16 x 6 single and 8 x 6 double
  are 24 `fmla`, 4 loads, no spills.
- **SSE2 x86_64 baseline and Cortex-A15 NEON (compile only):** 8 x 6 and
  8 x 5 single, multiply and add or `vfma`, no spills. ARMv7 NEON and
  pre-VSX AltiVec `f32` tiles are scalar since: those lanes flush
  denormals (`util.vectorIeee`).

## Tails

Row tails narrower than a vector go through one function per column width
(`thin`: each row against the columns as one vector) rather than a fully
unrolled tile per pair of widths: on Skylake-X `f32`, 10 unrolled tiles
and 5 `thin` instead of 30 tiles. Against the earlier all-unrolled tails,
interleaved on this host, `sgemm` and `dgemm` at n = 37, 100, 250 and 1000
moved within the noise (for example `sgemm` n = 100: 36.8 against 36.7
median). The objects did not shrink (Skylake-X 34 KB against 30 KB, the
AVX2 group 22 KB against 15 KB, stripped): LLVM unrolls the `thin` loops.

## Recommendation

Keep the generic Zig kernel for the AVX2 cores (Haswell, Zen and the cores
they serve): it runs at the measured peak, so converting `sgemm`/`dgemm`
`haswell` assembly cannot win in the kernel. What the full gemm loses (7 to
25%) is packing, the `C` update and edge tiles; improve the Zig copies first
(whole-vector transposes in `oncopy` and `incopy`), and measure again.

Convert assembly first where Zig cannot express what the core's kernel uses,
and measure there before keeping it:

1. The arm64 SVE kernels (A64FX, Neoverse V1, N2, V2): `@Vector` has no
   scalable vectors, so Zig falls back to NEON;
2. POWER10's MMA kernels, and the AMX and AVX-512 BF16 ones for `sbgemm`,
   which `@Vector` cannot reach at all;
3. Skylake-X and its successors, only if a machine with AVX-512 shows the
   Zig tile under 85% of the converted `skylakex` kernel (the decision rule
   of the research); its loop looks clean but is unmeasured.
