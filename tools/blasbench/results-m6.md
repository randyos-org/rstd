# BLAS kernels on the Apple M6 Mac mini

Whether a generic Zig `@Vector` kernel is as fast as the converted OpenBLAS
assembly on the Apple cores (`vortex`, `vortexm4`), per routine and
precision. Measured 2026-10-19 with `blas_bench.zig`:

```sh
zig build ssh-bench-aarch64-macos -Dopenblas=<OpenBLAS checkout> -j3
```

## Setup

- **Machine:** Apple M6 Mac mini, macOS 27.0.1 (26A434), 16 GB. Two "Super"
  cores (perflevel 0, 128 KiB L1d) and four "Performance" cores (perflevel 1,
  64 KiB L1d); `sysctl` gives 20 MiB of L2 for each. rstd detects `vortexm4` (SME, no SVE).
- **Build:** Zig 0.18.0-dev.267+de94864b9, ReleaseFast, `aarch64-macos`
  (CPU model `apple_m1`, Zig's baseline there), `-Dblas-cores=all`.
- **Reference:** OpenBLAS v0.3.34-317-ga69062524 for the gemm `.S`.
- **Peak:** 16 independent `fmla` chains measure 149.5 GFLOP/s single and
  74.7 double (up to 153.2 and 76.6 in other runs): four 128-bit FMA
  pipes at 4.78 GHz.

## Method

- **Harness:** `bench.zig`'s `race`. Every contestant is called through the
  vtable slot's kernel ABI (`daxpy_k`, `dgemv_n`, ...), from a function
  pointer, the same way for all.
- **Correctness first:** before timing, each contestant runs once on fresh
  data at every size and is compared with an `f64` reference: element-wise
  results within 8 eps of their terms, reductions within `(n + 8)` eps of
  the sum of magnitudes, `i?amax` exactly the first maximum, and every
  float outside the vectors (stride gaps, a guard past the end) unchanged.
  gemm checks 256 sampled entries of `C` (all of it for the kernel alone).
  No contestant failed.
- **Statistic:** the minimum time per call. Each batch runs at least 10 ms
  (20 ms for gemm), its call count doubled until it does; the best of 15
  batches (7 for gemm) counts. Batches go round the contestants in turn, so
  a disturbance falls on all of them alike.
- **Core placement:** the bench sets its thread's QoS to user-interactive,
  and reports the clock from a chain of dependent adds before each level
  and at the end: 4.78 to 4.79 GHz throughout, the Super cores. For
  contrast, the probe alone: default 4.78 GHz and 149.5 GFLOP/s,
  `taskpolicy -c utility` the same, `taskpolicy -c background` 2.90 GHz and
  69.6 GFLOP/s (the Performance cluster). Nothing was pinned.
- **Repeatability:** three full runs. Over the columns no change touched,
  the median spread between runs is 0.3%, and under 4% but for the cells in
  [Fairness](#fairness). The tables are the third run.
- **Sizes:** level 1, bytes per vector: 4 KiB and 32 KiB (L1), 1 MiB (L2),
  64 MiB (memory), and one strided case, `inc = 2` over a 32 KiB span.
  Level 2, square, bytes of `A`: 32 KiB, 4 MiB, 64 MiB. Rates count useful
  bytes (strided too) or flops (`2mn` real, `8mn` complex gemv; `2n^3`
  gemm).

### Contestants

- **asm:** the vtable slot of `cpu.vtableFromCoreType(.vortex)` where vortex
  has a converted kernel of its own. Where it has none (upstream's vortex
  kernel there is C: `dot.c`, `*_thunderx2t99.c`, `sgemv_n_neon.c`, not
  ported, so rstd runs the portable kernel), `armv8`'s converted kernel
  stands in; `none` where neither has one.
- **portable:** rstd's portable kernel (`cpu.generic`'s table): scalar
  loops.
- **vec, vec x4:** `vector_kernels.zig`, plain `@Vector` kernels one native
  vector (`vec`) or four (`vec x4`) wide per step; reductions keep two such
  accumulators. `i?amax` keeps the first block of `4 * W` with the largest
  maximum and searches only that block. `nrm2` sums unscaled squares and
  falls back to the scaled loop only where that could have overflowed or
  underflowed. Strided paths use four scalar accumulators (real) or one
  2-lane vector per element (complex).
- **gemm:** the Zig `@Vector` kernel of the `vortex`, `vortexm4` and
  portable builds, each with its own copies, against OpenBLAS's vortex
  kernels (`KERNEL.VORTEX` includes `KERNEL.NEOVERSEN1`):
  `sgemm_kernel_16x4.S` with `sgemm_tcopy_16.S` and `sgemm_ncopy_4.S`,
  `dgemm_kernel_8x4.S` with `dgemm_tcopy_8.S` and `dgemm_ncopy_4.S`,
  assembled from the checkout for Mach-O. `C += A*B`, `NN`, through a copy
  of `blas3/driver.zig`'s loops with vortex's `P`, `Q`, `R`.

`vortexm4` has no level-1 or level-2 kernel of its own: `KERNEL.VORTEXM4`
includes the same `KERNEL.NEOVERSEN1`, and the converted assembly is the
same text in both objects, so those tables hold for both cores. Its gemm is
raced separately. Upstream's one distinct vortexm4 path, the SME
`sgemm_direct` and `sme_sgemm_kernel` that `interface/gemm.c` takes for
clang builds, is not ported and was not measured.

## Level 1

GB/s; ratios are over the asm column.

### Single (`s`)

| routine | p | n | inc | asm | asm GB/s | portable | vec | vec x4 | portable/asm | vec/asm | vec x4/asm |
| ------- | - | -: | --: | --- | -------: | -------: | --: | -----: | -----------: | ------: | ---------: |
| axpy | s | 1024 | 1 | vortex | 217.1 | 207.1 | 207.6 | 293.2 | 0.95 | 0.96 | 1.35 |
| axpy | s | 8192 | 1 | vortex | 217.1 | 224.3 | 221.9 | 293.8 | 1.03 | 1.02 | 1.35 |
| axpy | s | 262144 | 1 | vortex | 159.5 | 170.2 | 169.0 | 170.4 | 1.07 | 1.06 | 1.07 |
| axpy | s | 16777216 | 1 | vortex | 138.7 | 136.7 | 136.1 | 139.8 | 0.99 | 0.98 | 1.01 |
| axpy | s | 4096 | 2 | vortex | 28.7 | 69.6 | 55.6 | 55.5 | 2.43 | 1.94 | 1.94 |
| dot | s | 1024 | 1 | none | - | 18.0 | 131.1 | 218.1 | - | - | - |
| dot | s | 8192 | 1 | none | - | 16.3 | 88.4 | 222.9 | - | - | - |
| dot | s | 262144 | 1 | none | - | 16.1 | 83.4 | 150.1 | - | - | - |
| dot | s | 16777216 | 1 | none | - | 16.1 | 83.2 | 127.2 | - | - | - |
| dot | s | 4096 | 2 | none | - | 16.6 | 42.3 | 42.2 | - | - | - |
| scal | s | 1024 | 1 | vortex | 233.4 | 64.2 | 139.5 | 275.6 | 0.27 | 0.60 | 1.18 |
| scal | s | 8192 | 1 | vortex | 231.1 | 66.0 | 150.3 | 288.8 | 0.29 | 0.65 | 1.25 |
| scal | s | 262144 | 1 | vortex | 164.2 | 65.0 | 152.8 | 162.8 | 0.40 | 0.93 | 0.99 |
| scal | s | 16777216 | 1 | vortex | 135.5 | 63.1 | 132.2 | 136.1 | 0.47 | 0.98 | 1.00 |
| scal | s | 4096 | 2 | vortex | 36.4 | 66.0 | 37.9 | 37.8 | 1.81 | 1.04 | 1.04 |
| nrm2 | s | 1024 | 1 | vortex | 5.1 | 6.9 | 72.1 | 159.7 | 1.35 | 14.22 | 31.50 |
| nrm2 | s | 8192 | 1 | vortex | 4.8 | 6.8 | 48.0 | 161.8 | 1.40 | 9.98 | 33.60 |
| nrm2 | s | 262144 | 1 | vortex | 4.8 | 6.8 | 45.7 | 151.4 | 1.41 | 9.55 | 31.65 |
| nrm2 | s | 16777216 | 1 | vortex | 4.8 | 6.7 | 46.1 | 111.7 | 1.41 | 9.65 | 23.37 |
| nrm2 | s | 4096 | 2 | vortex | 4.8 | 6.8 | 21.9 | 21.9 | 1.40 | 4.53 | 4.53 |
| asum | s | 1024 | 1 | armv8 | 87.3 | 10.4 | 83.6 | 141.0 | 0.12 | 0.96 | 1.61 |
| asum | s | 8192 | 1 | armv8 | 61.5 | 8.8 | 58.5 | 148.3 | 0.14 | 0.95 | 2.41 |
| asum | s | 262144 | 1 | armv8 | 58.3 | 8.6 | 55.6 | 143.1 | 0.15 | 0.95 | 2.45 |
| asum | s | 16777216 | 1 | armv8 | 58.1 | 8.6 | 56.0 | 109.7 | 0.15 | 0.96 | 1.89 |
| asum | s | 4096 | 2 | armv8 | 4.9 | 9.0 | 29.4 | 29.0 | 1.85 | 6.05 | 5.98 |
| iamax | s | 1024 | 1 | armv8 | 69.9 | 10.1 | 101.9 | 114.9 | 0.14 | 1.46 | 1.64 |
| iamax | s | 8192 | 1 | armv8 | 74.5 | 9.6 | 79.4 | 134.8 | 0.13 | 1.06 | 1.81 |
| iamax | s | 262144 | 1 | armv8 | 74.4 | 9.6 | 76.2 | 141.9 | 0.13 | 1.02 | 1.91 |
| iamax | s | 16777216 | 1 | armv8 | 73.8 | 9.6 | 74.8 | 107.1 | 0.13 | 1.01 | 1.45 |
| iamax | s | 4096 | 2 | armv8 | 4.8 | 9.7 | 51.3 | 38.8 | 2.02 | 10.72 | 8.10 |
| copy | s | 1024 | 1 | armv8 | 140.1 | 67.5 | 140.3 | 277.9 | 0.48 | 1.00 | 1.98 |
| copy | s | 8192 | 1 | armv8 | 149.4 | 69.6 | 146.2 | 291.8 | 0.47 | 0.98 | 1.95 |
| copy | s | 262144 | 1 | armv8 | 143.9 | 67.9 | 140.2 | 151.2 | 0.47 | 0.97 | 1.05 |
| copy | s | 16777216 | 1 | armv8 | 119.8 | 66.1 | 119.5 | 120.2 | 0.55 | 1.00 | 1.00 |
| copy | s | 4096 | 2 | armv8 | 37.8 | 69.3 | 37.7 | 37.7 | 1.83 | 1.00 | 1.00 |
| swap | s | 1024 | 1 | vortex | 188.0 | 65.5 | 248.7 | 275.1 | 0.35 | 1.32 | 1.46 |
| swap | s | 8192 | 1 | vortex | 195.6 | 67.9 | 253.5 | 292.9 | 0.35 | 1.30 | 1.50 |
| swap | s | 262144 | 1 | vortex | 102.2 | 67.0 | 156.5 | 164.2 | 0.66 | 1.53 | 1.61 |
| swap | s | 16777216 | 1 | vortex | 135.4 | 65.0 | 141.0 | 143.1 | 0.48 | 1.04 | 1.06 |
| swap | s | 4096 | 2 | vortex | 67.7 | 67.6 | 67.4 | 67.4 | 1.00 | 1.00 | 1.00 |
| rot | s | 1024 | 1 | vortex | 241.2 | 48.2 | 219.9 | 228.0 | 0.20 | 0.91 | 0.94 |
| rot | s | 8192 | 1 | vortex | 250.2 | 48.3 | 224.5 | 227.1 | 0.19 | 0.90 | 0.91 |
| rot | s | 262144 | 1 | vortex | 153.0 | 48.4 | 156.9 | 157.6 | 0.32 | 1.03 | 1.03 |
| rot | s | 16777216 | 1 | vortex | 141.1 | 48.0 | 141.1 | 142.1 | 0.34 | 1.00 | 1.01 |
| rot | s | 4096 | 2 | vortex | 19.0 | 48.4 | 64.3 | 64.3 | 2.55 | 3.38 | 3.38 |

### Double (`d`)

| routine | p | n | inc | asm | asm GB/s | portable | vec | vec x4 | portable/asm | vec/asm | vec x4/asm |
| ------- | - | -: | --: | --- | -------: | -------: | --: | -----: | -----------: | ------: | ---------: |
| axpy | d | 512 | 1 | vortex | 287.7 | 201.1 | 202.3 | 288.4 | 0.70 | 0.70 | 1.00 |
| axpy | d | 4096 | 1 | vortex | 254.4 | 216.8 | 217.3 | 289.7 | 0.85 | 0.85 | 1.14 |
| axpy | d | 131072 | 1 | vortex | 165.8 | 167.2 | 167.1 | 166.5 | 1.01 | 1.01 | 1.00 |
| axpy | d | 8388608 | 1 | vortex | 147.5 | 144.5 | 145.6 | 146.7 | 0.98 | 0.99 | 0.99 |
| axpy | d | 2048 | 2 | vortex | 56.5 | 136.7 | 108.0 | 55.9 | 2.42 | 1.91 | 0.99 |
| dot | d | 512 | 1 | armv8 | 93.5 | 48.8 | 129.3 | 211.9 | 0.52 | 1.38 | 2.27 |
| dot | d | 4096 | 1 | armv8 | 65.2 | 34.9 | 88.3 | 215.9 | 0.54 | 1.35 | 3.31 |
| dot | d | 131072 | 1 | armv8 | 62.8 | 33.6 | 83.2 | 146.0 | 0.53 | 1.33 | 2.33 |
| dot | d | 8388608 | 1 | armv8 | 61.9 | 33.5 | 82.9 | 124.6 | 0.54 | 1.34 | 2.01 |
| dot | d | 2048 | 2 | armv8 | 23.9 | 36.4 | 86.0 | 85.6 | 1.52 | 3.61 | 3.59 |
| scal | d | 512 | 1 | vortex | 255.5 | 129.3 | 136.3 | 273.2 | 0.51 | 0.53 | 1.07 |
| scal | d | 4096 | 1 | vortex | 256.7 | 130.4 | 147.6 | 288.4 | 0.51 | 0.58 | 1.12 |
| scal | d | 131072 | 1 | vortex | 162.8 | 121.5 | 148.3 | 159.3 | 0.75 | 0.91 | 0.98 |
| scal | d | 8388608 | 1 | vortex | 138.2 | 119.9 | 133.5 | 138.0 | 0.87 | 0.97 | 1.00 |
| scal | d | 2048 | 2 | vortex | 71.0 | 129.7 | 74.9 | 74.8 | 1.83 | 1.06 | 1.05 |
| nrm2 | d | 512 | 1 | vortex | 10.7 | 13.9 | 72.2 | 169.5 | 1.29 | 6.72 | 15.77 |
| nrm2 | d | 4096 | 1 | vortex | 9.7 | 13.6 | 48.1 | 160.9 | 1.40 | 4.96 | 16.59 |
| nrm2 | d | 131072 | 1 | vortex | 9.6 | 13.5 | 45.7 | 148.2 | 1.41 | 4.78 | 15.49 |
| nrm2 | d | 8388608 | 1 | vortex | 9.6 | 13.5 | 46.4 | 108.5 | 1.41 | 4.85 | 11.35 |
| nrm2 | d | 2048 | 2 | vortex | 9.7 | 13.6 | 46.1 | 46.0 | 1.40 | 4.74 | 4.72 |
| asum | d | 512 | 1 | armv8 | 92.7 | 26.7 | 84.2 | 143.2 | 0.29 | 0.91 | 1.55 |
| asum | d | 4096 | 1 | armv8 | 67.2 | 18.0 | 58.2 | 146.9 | 0.27 | 0.87 | 2.19 |
| asum | d | 131072 | 1 | armv8 | 61.2 | 17.3 | 55.2 | 141.7 | 0.28 | 0.90 | 2.32 |
| asum | d | 8388608 | 1 | armv8 | 61.7 | 17.2 | 55.9 | 108.1 | 0.28 | 0.91 | 1.75 |
| asum | d | 2048 | 2 | armv8 | 9.8 | 18.8 | 57.2 | 56.6 | 1.91 | 5.81 | 5.75 |
| iamax | d | 512 | 1 | armv8 | 122.9 | 21.1 | 99.2 | 138.5 | 0.17 | 0.81 | 1.13 |
| iamax | d | 4096 | 1 | armv8 | 129.3 | 19.1 | 77.9 | 138.9 | 0.15 | 0.60 | 1.07 |
| iamax | d | 131072 | 1 | armv8 | 129.9 | 18.7 | 74.5 | 141.1 | 0.14 | 0.57 | 1.09 |
| iamax | d | 8388608 | 1 | armv8 | 99.7 | 18.7 | 73.4 | 106.3 | 0.19 | 0.74 | 1.07 |
| iamax | d | 2048 | 2 | armv8 | 9.6 | 19.6 | 77.3 | 93.9 | 2.05 | 8.08 | 9.81 |
| copy | d | 512 | 1 | armv8 | 158.2 | 137.0 | 138.2 | 275.4 | 0.87 | 0.87 | 1.74 |
| copy | d | 4096 | 1 | armv8 | 148.5 | 137.6 | 143.5 | 288.5 | 0.93 | 0.97 | 1.94 |
| copy | d | 131072 | 1 | armv8 | 143.7 | 109.5 | 139.0 | 147.7 | 0.76 | 0.97 | 1.03 |
| copy | d | 8388608 | 1 | armv8 | 122.1 | 93.0 | 120.4 | 121.8 | 0.76 | 0.99 | 1.00 |
| copy | d | 2048 | 2 | armv8 | 73.1 | 136.5 | 74.3 | 74.4 | 1.87 | 1.02 | 1.02 |
| swap | d | 512 | 1 | vortex | 187.0 | 124.6 | 249.3 | 272.9 | 0.67 | 1.33 | 1.46 |
| swap | d | 4096 | 1 | vortex | 194.8 | 134.0 | 253.5 | 291.1 | 0.69 | 1.30 | 1.49 |
| swap | d | 131072 | 1 | vortex | 100.3 | 114.9 | 154.3 | 161.8 | 1.15 | 1.54 | 1.61 |
| swap | d | 8388608 | 1 | vortex | 137.2 | 110.8 | 147.2 | 150.0 | 0.81 | 1.07 | 1.09 |
| swap | d | 2048 | 2 | vortex | 133.3 | 127.2 | 128.1 | 128.2 | 0.95 | 0.96 | 0.96 |
| rot | d | 512 | 1 | vortex | 232.3 | 96.2 | 219.7 | 227.9 | 0.41 | 0.95 | 0.98 |
| rot | d | 4096 | 1 | vortex | 233.0 | 96.6 | 224.7 | 227.8 | 0.41 | 0.96 | 0.98 |
| rot | d | 131072 | 1 | vortex | 159.0 | 96.7 | 156.3 | 157.6 | 0.61 | 0.98 | 0.99 |
| rot | d | 8388608 | 1 | vortex | 150.4 | 95.2 | 148.6 | 149.7 | 0.63 | 0.99 | 1.00 |
| rot | d | 2048 | 2 | vortex | 36.8 | 96.6 | 124.6 | 125.2 | 2.63 | 3.39 | 3.40 |

### Single complex (`c`)

| routine | p | n | inc | asm | asm GB/s | portable | vec | vec x4 | portable/asm | vec/asm | vec x4/asm |
| ------- | - | -: | --: | --- | -------: | -------: | --: | -----: | -----------: | ------: | ---------: |
| axpy | c | 512 | 1 | vortex | 144.4 | 124.1 | 105.7 | 283.4 | 0.86 | 0.73 | 1.96 |
| axpy | c | 4096 | 1 | vortex | 145.6 | 124.9 | 109.4 | 283.5 | 0.86 | 0.75 | 1.95 |
| axpy | c | 131072 | 1 | vortex | 146.5 | 125.7 | 110.0 | 165.8 | 0.86 | 0.75 | 1.13 |
| axpy | c | 8388608 | 1 | vortex | 140.4 | 123.5 | 109.5 | 144.9 | 0.88 | 0.78 | 1.03 |
| axpy | c | 2048 | 2 | vortex | 109.4 | 72.5 | 107.7 | 107.7 | 0.66 | 0.98 | 0.98 |
| dot | c | 512 | 1 | none | - | 28.1 | 100.5 | 149.3 | - | - | - |
| dot | c | 4096 | 1 | none | - | 26.3 | 79.6 | 184.8 | - | - | - |
| dot | c | 131072 | 1 | none | - | 26.0 | 77.0 | 141.4 | - | - | - |
| dot | c | 8388608 | 1 | none | - | 26.0 | 76.6 | 121.2 | - | - | - |
| dot | c | 2048 | 2 | none | - | 26.8 | 41.5 | 41.7 | - | - | - |
| scal | c | 512 | 1 | vortex | 96.0 | 48.6 | 132.7 | 189.9 | 0.51 | 1.38 | 1.98 |
| scal | c | 4096 | 1 | vortex | 97.0 | 48.6 | 143.4 | 191.2 | 0.50 | 1.48 | 1.97 |
| scal | c | 131072 | 1 | vortex | 97.7 | 48.6 | 143.5 | 157.2 | 0.50 | 1.47 | 1.61 |
| scal | c | 8388608 | 1 | vortex | 97.2 | 48.6 | 133.3 | 138.1 | 0.50 | 1.37 | 1.42 |
| scal | c | 2048 | 2 | vortex | 42.0 | 48.7 | 71.3 | 71.5 | 1.16 | 1.70 | 1.70 |
| nrm2 | c | 512 | 1 | vortex | 5.1 | 6.8 | 71.9 | 148.7 | 1.34 | 14.13 | 29.23 |
| nrm2 | c | 4096 | 1 | vortex | 4.8 | 6.7 | 48.0 | 152.3 | 1.38 | 9.97 | 31.63 |
| nrm2 | c | 131072 | 1 | vortex | 4.8 | 6.6 | 45.6 | 147.1 | 1.39 | 9.57 | 30.83 |
| nrm2 | c | 8388608 | 1 | vortex | 4.8 | 6.6 | 45.5 | 109.8 | 1.39 | 9.54 | 23.00 |
| nrm2 | c | 2048 | 2 | vortex | 4.8 | 6.7 | 23.6 | 23.6 | 1.38 | 4.88 | 4.87 |
| asum | c | 512 | 1 | armv8 | 92.3 | 9.7 | 80.9 | 137.7 | 0.10 | 0.88 | 1.49 |
| asum | c | 4096 | 1 | armv8 | 65.4 | 8.6 | 56.3 | 145.8 | 0.13 | 0.86 | 2.23 |
| asum | c | 131072 | 1 | armv8 | 60.6 | 8.2 | 54.2 | 141.8 | 0.14 | 0.90 | 2.34 |
| asum | c | 8388608 | 1 | armv8 | 61.1 | 8.3 | 54.2 | 107.4 | 0.14 | 0.89 | 1.76 |
| asum | c | 2048 | 2 | armv8 | 16.6 | 8.6 | 28.6 | 28.6 | 0.52 | 1.72 | 1.72 |
| iamax | c | 512 | 1 | armv8 | 120.7 | 20.3 | 122.6 | 86.6 | 0.17 | 1.02 | 0.72 |
| iamax | c | 4096 | 1 | armv8 | 129.1 | 18.7 | 130.4 | 145.5 | 0.15 | 1.01 | 1.13 |
| iamax | c | 131072 | 1 | armv8 | 130.1 | 18.4 | 130.1 | 144.4 | 0.14 | 1.00 | 1.11 |
| iamax | c | 8388608 | 1 | armv8 | 98.8 | 18.4 | 106.3 | 106.7 | 0.19 | 1.08 | 1.08 |
| iamax | c | 2048 | 2 | armv8 | 19.0 | 19.0 | 51.2 | 48.0 | 1.00 | 2.69 | 2.52 |
| copy | c | 512 | 1 | armv8 | 157.2 | 64.2 | 136.8 | 274.8 | 0.41 | 0.87 | 1.75 |
| copy | c | 4096 | 1 | armv8 | 146.6 | 69.6 | 144.3 | 288.6 | 0.47 | 0.98 | 1.97 |
| copy | c | 131072 | 1 | armv8 | 142.1 | 69.7 | 138.2 | 151.9 | 0.49 | 0.97 | 1.07 |
| copy | c | 8388608 | 1 | armv8 | 121.9 | 68.4 | 120.4 | 121.7 | 0.56 | 0.99 | 1.00 |
| copy | c | 2048 | 2 | armv8 | 73.3 | 68.7 | 135.3 | 133.6 | 0.94 | 1.85 | 1.82 |
| swap | c | 512 | 1 | vortex | 186.3 | 66.4 | 246.7 | 272.9 | 0.36 | 1.32 | 1.46 |
| swap | c | 4096 | 1 | vortex | 194.3 | 67.8 | 263.3 | 290.3 | 0.35 | 1.35 | 1.49 |
| swap | c | 131072 | 1 | vortex | 100.9 | 67.0 | 157.9 | 159.9 | 0.66 | 1.57 | 1.59 |
| swap | c | 8388608 | 1 | vortex | 136.9 | 65.0 | 146.8 | 147.9 | 0.47 | 1.07 | 1.08 |
| swap | c | 2048 | 2 | vortex | 131.7 | 67.6 | 127.4 | 127.0 | 0.51 | 0.97 | 0.96 |
| rot | c | 512 | 1 | vortex | 264.1 | 48.1 | 219.0 | 227.7 | 0.18 | 0.83 | 0.86 |
| rot | c | 4096 | 1 | vortex | 274.5 | 48.0 | 224.4 | 226.6 | 0.18 | 0.82 | 0.83 |
| rot | c | 131072 | 1 | vortex | 158.9 | 48.2 | 157.0 | 157.7 | 0.30 | 0.99 | 0.99 |
| rot | c | 8388608 | 1 | vortex | 148.0 | 47.8 | 146.7 | 147.5 | 0.32 | 0.99 | 1.00 |
| rot | c | 2048 | 2 | vortex | 126.2 | 48.1 | 113.5 | 113.4 | 0.38 | 0.90 | 0.90 |

### Double complex (`z`)

| routine | p | n | inc | asm | asm GB/s | portable | vec | vec x4 | portable/asm | vec/asm | vec x4/asm |
| ------- | - | -: | --: | --- | -------: | -------: | --: | -----: | -----------: | ------: | ---------: |
| axpy | z | 256 | 1 | vortex | 144.4 | 124.1 | 200.5 | 281.9 | 0.86 | 1.39 | 1.95 |
| axpy | z | 2048 | 1 | vortex | 145.6 | 124.9 | 215.4 | 281.2 | 0.86 | 1.48 | 1.93 |
| axpy | z | 65536 | 1 | vortex | 146.6 | 125.7 | 167.0 | 166.1 | 0.86 | 1.14 | 1.13 |
| axpy | z | 4194304 | 1 | vortex | 142.1 | 123.4 | 148.6 | 148.5 | 0.87 | 1.05 | 1.05 |
| axpy | z | 1024 | 2 | vortex | 217.8 | 145.1 | 210.8 | 209.3 | 0.67 | 0.97 | 0.96 |
| dot | z | 256 | 1 | none | - | 62.1 | 101.9 | 171.7 | - | - | - |
| dot | z | 2048 | 1 | none | - | 53.0 | 77.4 | 188.9 | - | - | - |
| dot | z | 65536 | 1 | none | - | 52.0 | 74.0 | 141.5 | - | - | - |
| dot | z | 4194304 | 1 | none | - | 51.9 | 74.0 | 121.5 | - | - | - |
| dot | z | 1024 | 2 | none | - | 53.8 | 84.7 | 84.8 | - | - | - |
| scal | z | 256 | 1 | vortex | 96.0 | 97.2 | 128.4 | 190.0 | 1.01 | 1.34 | 1.98 |
| scal | z | 2048 | 1 | vortex | 97.0 | 97.2 | 144.1 | 190.6 | 1.00 | 1.48 | 1.96 |
| scal | z | 65536 | 1 | vortex | 97.8 | 97.2 | 146.5 | 157.2 | 0.99 | 1.50 | 1.61 |
| scal | z | 4194304 | 1 | vortex | 97.2 | 97.0 | 133.5 | 138.5 | 1.00 | 1.37 | 1.42 |
| scal | z | 1024 | 2 | vortex | 83.5 | 97.2 | 139.7 | 139.8 | 1.16 | 1.67 | 1.67 |
| nrm2 | z | 256 | 1 | vortex | 10.8 | 13.8 | 71.5 | 168.7 | 1.28 | 6.65 | 15.68 |
| nrm2 | z | 2048 | 1 | vortex | 9.7 | 13.3 | 48.1 | 159.4 | 1.37 | 4.95 | 16.43 |
| nrm2 | z | 65536 | 1 | vortex | 9.6 | 13.2 | 45.5 | 147.8 | 1.38 | 4.76 | 15.47 |
| nrm2 | z | 4194304 | 1 | vortex | 9.5 | 13.2 | 45.8 | 108.1 | 1.38 | 4.79 | 11.32 |
| nrm2 | z | 1024 | 2 | vortex | 9.8 | 13.4 | 49.2 | 48.1 | 1.36 | 5.02 | 4.91 |
| asum | z | 256 | 1 | armv8 | 92.2 | 24.1 | 81.5 | 142.9 | 0.26 | 0.88 | 1.55 |
| asum | z | 2048 | 1 | armv8 | 64.6 | 17.3 | 56.3 | 146.3 | 0.27 | 0.87 | 2.27 |
| asum | z | 65536 | 1 | armv8 | 62.7 | 16.5 | 53.9 | 141.6 | 0.26 | 0.86 | 2.26 |
| asum | z | 4194304 | 1 | armv8 | 63.7 | 16.6 | 54.3 | 105.8 | 0.26 | 0.85 | 1.66 |
| asum | z | 1024 | 2 | armv8 | 37.4 | 17.8 | 63.4 | 62.8 | 0.48 | 1.69 | 1.68 |
| iamax | z | 256 | 1 | armv8 | 127.0 | 22.1 | 124.2 | 109.5 | 0.17 | 0.98 | 0.86 |
| iamax | z | 2048 | 1 | armv8 | 136.8 | 19.4 | 129.6 | 144.0 | 0.14 | 0.95 | 1.05 |
| iamax | z | 65536 | 1 | armv8 | 137.5 | 18.7 | 130.2 | 144.8 | 0.14 | 0.95 | 1.05 |
| iamax | z | 4194304 | 1 | armv8 | 101.2 | 18.7 | 104.9 | 105.5 | 0.18 | 1.04 | 1.04 |
| iamax | z | 1024 | 2 | armv8 | 39.0 | 19.9 | 90.0 | 89.4 | 0.51 | 2.31 | 2.29 |
| copy | z | 256 | 1 | armv8 | 154.9 | 130.1 | 138.6 | 274.5 | 0.84 | 0.89 | 1.77 |
| copy | z | 2048 | 1 | armv8 | 147.7 | 137.4 | 145.7 | 288.8 | 0.93 | 0.99 | 1.96 |
| copy | z | 65536 | 1 | armv8 | 141.2 | 135.7 | 139.0 | 146.1 | 0.96 | 0.98 | 1.03 |
| copy | z | 4194304 | 1 | armv8 | 121.6 | 117.4 | 121.1 | 121.8 | 0.96 | 1.00 | 1.00 |
| copy | z | 1024 | 2 | armv8 | 146.4 | 134.0 | 229.0 | 229.1 | 0.92 | 1.56 | 1.57 |
| swap | z | 256 | 1 | vortex | 186.3 | 131.8 | 247.1 | 272.6 | 0.71 | 1.33 | 1.46 |
| swap | z | 2048 | 1 | vortex | 194.3 | 135.2 | 263.9 | 290.1 | 0.70 | 1.36 | 1.49 |
| swap | z | 65536 | 1 | vortex | 99.8 | 115.0 | 157.3 | 159.5 | 1.15 | 1.58 | 1.60 |
| swap | z | 4194304 | 1 | vortex | 136.9 | 110.5 | 149.0 | 150.0 | 0.81 | 1.09 | 1.10 |
| swap | z | 1024 | 2 | vortex | 237.9 | 134.5 | 216.6 | 216.5 | 0.57 | 0.91 | 0.91 |
| rot | z | 256 | 1 | vortex | 265.4 | 96.2 | 140.7 | 227.8 | 0.36 | 0.53 | 0.86 |
| rot | z | 2048 | 1 | vortex | 274.5 | 96.1 | 145.4 | 227.8 | 0.35 | 0.53 | 0.83 |
| rot | z | 65536 | 1 | vortex | 159.0 | 96.4 | 146.6 | 157.7 | 0.61 | 0.92 | 0.99 |
| rot | z | 4194304 | 1 | vortex | 150.2 | 95.0 | 143.9 | 149.9 | 0.63 | 0.96 | 1.00 |
| rot | z | 1024 | 2 | vortex | 226.9 | 96.1 | 222.4 | 222.7 | 0.42 | 0.98 | 0.98 |

## Level 2

GFLOP/s, square, unit strides.

| routine | p | n | asm | asm GFLOP/s | portable | vec | vec x4 | portable/asm | vec/asm | vec x4/asm |
| ------- | - | -: | --- | ----------: | -------: | --: | -----: | -----------: | ------: | ---------: |
| gemv_n | s | 90 | armv8 | 36.4 | 4.5 | 42.1 | 61.3 | 0.12 | 1.16 | 1.69 |
| gemv_n | s | 1024 | armv8 | 37.8 | 4.5 | 47.6 | 52.5 | 0.12 | 1.26 | 1.39 |
| gemv_n | s | 4096 | armv8 | 38.1 | 4.6 | 35.9 | 38.5 | 0.12 | 0.94 | 1.01 |
| gemv_t | s | 90 | vortex | 38.3 | 7.2 | 60.3 | 41.5 | 0.19 | 1.58 | 1.08 |
| gemv_t | s | 1024 | vortex | 40.2 | 4.0 | 42.7 | 55.1 | 0.10 | 1.06 | 1.37 |
| gemv_t | s | 4096 | vortex | 38.9 | 3.8 | 33.7 | 40.3 | 0.10 | 0.87 | 1.04 |
| gemv_n | d | 64 | vortex | 20.5 | 8.6 | 23.5 | 36.6 | 0.42 | 1.15 | 1.78 |
| gemv_n | d | 724 | vortex | 19.0 | 8.7 | 25.4 | 29.0 | 0.46 | 1.34 | 1.53 |
| gemv_n | d | 2896 | vortex | 19.0 | 8.9 | 17.6 | 18.8 | 0.47 | 0.93 | 0.99 |
| gemv_t | d | 64 | vortex | 20.1 | 7.8 | 35.3 | 35.6 | 0.39 | 1.75 | 1.77 |
| gemv_t | d | 724 | vortex | 20.0 | 4.2 | 21.2 | 29.7 | 0.21 | 1.06 | 1.48 |
| gemv_t | d | 2896 | vortex | 19.5 | 3.8 | 16.2 | 19.2 | 0.20 | 0.83 | 0.98 |
| gemv_n | c | 64 | vortex | 34.8 | 18.0 | 85.1 | 85.5 | 0.52 | 2.44 | 2.45 |
| gemv_n | c | 724 | vortex | 35.8 | 18.2 | 92.9 | 91.3 | 0.51 | 2.59 | 2.55 |
| gemv_n | c | 2896 | vortex | 36.4 | 18.2 | 66.0 | 63.2 | 0.50 | 1.81 | 1.74 |
| gemv_t | c | 64 | vortex | 45.8 | 18.0 | 72.1 | 50.7 | 0.39 | 1.58 | 1.11 |
| gemv_t | c | 724 | vortex | 25.7 | 14.0 | 77.6 | 96.0 | 0.54 | 3.01 | 3.73 |
| gemv_t | c | 2896 | vortex | 22.7 | 13.3 | 60.5 | 69.7 | 0.59 | 2.66 | 3.06 |
| gemv_n | z | 45 | vortex | 23.4 | 17.9 | 42.8 | 42.3 | 0.77 | 1.83 | 1.81 |
| gemv_n | z | 512 | vortex | 24.3 | 18.2 | 46.1 | 45.9 | 0.75 | 1.89 | 1.89 |
| gemv_n | z | 2048 | vortex | 24.2 | 18.2 | 35.2 | 35.9 | 0.75 | 1.45 | 1.48 |
| gemv_t | z | 45 | vortex | 27.7 | 17.9 | 39.4 | 31.8 | 0.65 | 1.42 | 1.15 |
| gemv_t | z | 512 | vortex | 22.8 | 14.3 | 38.4 | 44.5 | 0.63 | 1.69 | 1.95 |
| gemv_t | z | 2048 | vortex | 21.2 | 13.4 | 33.3 | 35.4 | 0.63 | 1.57 | 1.67 |

## Level 3

`sgemm` kernel alone: 64 x 24, depth 256, on panels its own copies packed
(in L1 and L2).

| kernel | tile | GFLOP/s | over asm |
| ------ | ---- | ------: | -------: |
| asm 16x4 | 16 x 4 | 129.0 | 1.00 |
| zig vortex | 16 x 6 | 125.9 | 0.98 |
| zig vortexm4 | 16 x 6 | 122.7 | 0.95 |
| zig portable | 16 x 6 | 122.7 | 0.95 |

`sgemm` NN, square (vortex `P` 512, `Q` 1024, `R` 4096).

| n | asm 16x4 | zig vortex | zig vortexm4 | zig portable | zig vortex/asm | zig vortexm4/asm | zig portable/asm |
| -: | -: | -: | -: | -: | -: | -: | -: |
| 64 | 121.9 | 114.3 | 112.8 | 112.9 | 0.94 | 0.93 | 0.93 |
| 128 | 122.9 | 117.3 | 113.8 | 113.8 | 0.95 | 0.93 | 0.93 |
| 256 | 122.9 | 120.8 | 118.7 | 118.8 | 0.98 | 0.97 | 0.97 |
| 512 | 126.2 | 122.7 | 120.5 | 120.6 | 0.97 | 0.95 | 0.96 |
| 1024 | 127.5 | 124.7 | 122.2 | 122.2 | 0.98 | 0.96 | 0.96 |
| 2048 | 127.3 | 124.3 | 121.8 | 121.8 | 0.98 | 0.96 | 0.96 |

`dgemm` kernel alone, as above.

| kernel | tile | GFLOP/s | over asm |
| ------ | ---- | ------: | -------: |
| asm 8x4 | 8 x 4 | 64.7 | 1.00 |
| zig vortex | 8 x 6 | 62.0 | 0.96 |
| zig vortexm4 | 8 x 6 | 62.0 | 0.96 |
| zig portable | 8 x 6 | 61.9 | 0.96 |

`dgemm` NN, square (vortex `P` 256, `Q` 512, `R` 4096).

| n | asm 8x4 | zig vortex | zig vortexm4 | zig portable | zig vortex/asm | zig vortexm4/asm | zig portable/asm |
| -: | -: | -: | -: | -: | -: | -: | -: |
| 64 | 60.1 | 56.9 | 57.1 | 57.3 | 0.95 | 0.95 | 0.95 |
| 128 | 60.7 | 58.7 | 58.7 | 58.7 | 0.97 | 0.97 | 0.97 |
| 256 | 61.7 | 59.8 | 59.6 | 59.7 | 0.97 | 0.97 | 0.97 |
| 512 | 63.4 | 60.6 | 60.4 | 60.4 | 0.96 | 0.95 | 0.95 |
| 1024 | 63.8 | 61.0 | 60.9 | 61.0 | 0.96 | 0.96 | 0.96 |
| 2048 | 64.0 | 61.0 | 61.0 | 60.9 | 0.95 | 0.95 | 0.95 |

`cgemm` and `zgemm`: no Zig kernel yet (only real gemm is ported), so not
raced.

## Recommendations

By the design rule, tuning only where `@Vector` can't. "x4" is the 4-wide
`@Vector` kernel, plain Zig; the ratios are over the asm.

| Family | Recommendation | Evidence |
| ------ | -------------- | -------- |
| `axpy` | Drop the asm for `@Vector` (x4) | x4 0.99 to 1.35 real, 1.03 to 1.96 complex, unit stride; strided real Zig about 1.9 times the asm (once 0.99, see [Fairness](#fairness)), strided complex 0.96 to 0.98 |
| `dot` | `@Vector` (x4); there is no vortex asm | vortex's `dot.c` and `zdot_thunderx2t99.c` are not ported; x4 2.0 to 3.3 times armv8's `ddot` asm, 2.3 to 13.7 times portable elsewhere |
| `scal` | Drop the asm for `@Vector` (x4) | x4 0.98 to 1.25 real, 1.42 to 1.98 complex; strided 1.04 to 1.70 |
| `nrm2` | Drop the asm for `@Vector` | the asm is the scaled one-pass loop: x4 11 to 34 times faster at unit stride, 4.5 to 4.9 strided |
| `asum` | `@Vector` x4; no vortex asm | armv8 asm against x4: 1.49 to 2.45; plain `vec` 0.85 to 0.96, so the unroll matters |
| `i?amax` | Tuned `@Vector` (block maximum); no vortex asm | armv8 asm against x4: 1.07 to 1.91 real, 0.72 to 1.13 complex (below 1 only at 4 KiB); a per-lane index `@Vector` was 0.44 to 0.79 |
| `copy` | `@Vector` x4; no vortex asm | armv8 asm against x4: 1.00 to 1.98; strided 1.00 to 1.82 |
| `swap` | Drop the asm for `@Vector` (x4) | x4 1.06 to 1.61; strided 0.91 to 1.00 |
| `rot` | Needs a tuned `@Vector`; keep the asm until then | in L1 x4 is 0.91 to 0.98 real and 0.83 to 0.86 complex; L2 and memory 0.99 to 1.03; strided real 3.4 times faster in Zig, complex 0.90 to 0.98 |
| `gemv_n` | Drop the asm for `@Vector` (x4) | x4 1.39 to 1.78 real in L1 and L2, 0.99 to 1.01 from memory; complex 1.48 to 2.55 |
| `gemv_t` | Drop the asm for `@Vector` (x4) | x4 1.08 to 1.77 real in L1 and L2, 0.98 to 1.04 from memory; complex 1.11 to 3.73 |
| `gemm` (s, d) | Keep the Zig kernel; no asm | kernel alone 0.95 to 0.98 of the asm, full gemm 0.93 to 0.98 (lowest at 64 and 128); the asm kernel reaches 86% of peak, Zig 83 to 84% |
| `gemm` (c, z) | Not decided | no Zig complex kernel to race |
| `vortexm4` | Nothing of its own to keep | no level-1 or level-2 kernels of its own; SME is where `@Vector` cannot reach, unmeasured |

## Fairness

What could tilt the comparison:

- **Stand-in assembly:** for `dot`, `asum`, `i?amax`, `copy` and `sgemv_n`
  the asm column is armv8's generic arm64 `.S`, not vortex's upstream kernel,
  which is C with inline assembly or intrinsics and is not ported. A
  vortex-tuned asm could do better than the column shows.
- **gemm driver and assembly:** rstd's `driver.gemm` is not public, so the
  bench runs a copy of its loops; the rates are kernel plus copies through
  identical loops, not rstd's gemm entry. The upstream `.S` are assembled
  directly, not through `tools/blas_asm`, which converts no arm64 level-3
  kernel yet. The Zig tile is 16 x 6 and 8 x 6, the asm's 16 x 4 and 8 x 4.
- **Codegen swings:** three cells changed by up to 2 times between builds
  whose source for them did not change: portable `sgemv_n` (8.8 GFLOP/s in
  runs 1 and 2, 4.5 in run 3), the plain `vec` `zrot` at unit stride (0.83
  of asm in run 1, 0.53 after), and `vec x4` strided `daxpy` (1.89 in runs 1
  and 2, 0.99 in run 3). Plain Zig loops are sensitive to LLVM's choices; a
  kept `@Vector` kernel wants checking per build.
- **vortex and vortexm4 Zig gemm:** the same source runs 2 to 3% slower in
  the vortexm4 object than in the vortex one (`sgemm`, every run); the
  objects differ only in target features.
- **Semantics:** the candidates handle positive strides only. `vec` `sdot`
  sums in `f32`, as the asm does; rstd's portable `sdot` sums in `f64`.
- **Machine:** `audiomxd` and `configd` used about one core between them
  throughout; with six cores the bench kept a Super core (the clock checks),
  and the minimum over interleaved batches filters what is left.
