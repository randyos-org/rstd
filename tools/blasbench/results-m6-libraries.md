# BLAS libraries on the Apple M6 Mac mini

[`results-m6.md`](results-m6.md) found plain Zig `@Vector` x4 kernels
matching or beating rstd's converted OpenBLAS assembly on `vortex`. This
races rstd's public routines against the real libraries, to tell whether
that is our conversion, OpenBLAS's Apple kernels, or Zig. Measured
2026-10-22 with `blas_bench.zig`'s `lib1`, `lib2` and `lib3` sections:

```sh
zig build remote-bench-aarch64-macos -Dopenblas=<OpenBLAS checkout> \
  -Dopenblas-lib=<libopenblas.dylib on the Mac> -j3 -- probe lib1 lib2 lib3
```

## Answers

- **Our conversion costs nothing.** Where rstd and OpenBLAS run the same
  kernel (`axpy`, `scal`, `nrm2`, `rot`, `swap`, every `gemv` but `sgemv_n`),
  rstd's public routine runs at 0.98 to 1.04 of OpenBLAS's `cblas_*`, and
  the converted kernel alone at 0.98 to 1.06. OpenBLAS's own gemm equals
  its vortex `.S` driven by the bench's copy of `level3.c` (0.99 to 1.01).
- **OpenBLAS's real C kernels are better than the armv8 stand-ins**, and
  `@Vector` x4 still matches or beats them at unit stride. Real `ddot` is
  1.7 to 2.6 times armv8's, real `?asum` 1.5 to 2.5 times. x4 against the
  real kernels: `asum` 0.96 to 1.04, `i?amax` 0.86 to 1.40, real `dot` 1.25
  to 1.34 in L1 and 1.01 to 1.02 beyond, complex `dot` 1.14 to 1.79, `copy`
  1.95 to 2.05 in L1, `sgemv_n` 0.82 to 0.97. From memory they are within
  4% of each other, but for complex `dot` (1.14).
- **rstd itself is far behind wherever vortex has no converted kernel.**
  The real kernels are C, never ported, so rstd runs its portable loops:
  `?dot` at 0.09 to 0.61 of OpenBLAS, `?asum` 0.06 to 0.18, `i?amax` 0.07 to
  0.39, `copy` 0.46 to 0.98 and `sgemv_n` 0.12 to 0.22.
- **Accelerate is in another class wherever SME can work.** Single-threaded
  `sgemm` runs at up to 1,923 GFLOP/s, 15 times OpenBLAS and 13 times the
  NEON peak of 149.5; `dgemm` reaches 511, 8 times. No core's NEON pipes do
  that: this is the SME unit (`FEAT_SME`, `FEAT_SME2`), similar in design to the undocumented AMX of M1 to M3 (Remke and Breuer, "Hello SME! Generating Fast Matrix Multiplication Kernels Using the Scalable Matrix Extension", 2024, <https://arxiv.org/abs/2409.18779>).
  With 1 MiB vectors, level 1 runs at 2.6 to 8.3 times OpenBLAS (all but
  `asum`, `i?amax`, complex `dot` and `dznrm2`), and `gemv` in L2
  at 6.6 to 18 times (`sgemv_n` 476 GFLOP/s): rates past a NEON core's,
  which fit SME's streaming loads too. In L1 it wins less (1.2 to 3.9
  times), from memory it is 0.8 to 1.6 times, and its `i?amax` (0.26 to
  0.75) and `cdotu` (0.73 to 0.96) are slower than OpenBLAS. Its `dznrm2`
  barely beats OpenBLAS's scaled loop, where x4 is 11 to 16 times faster.
- **The `@Vector` win is the kernels' tuning, not Zig.** The hot loops
  (below) are what LLVM emits for any 4-wide unroll: eight independent
  accumulators, paired 128-bit loads, immediate offsets. The OpenBLAS
  kernels it beats were written for ThunderX2 and Neoverse N1: software
  prefetches every 128 bytes that the M6 does not need and that cost load
  slots, 16- or 32-byte loops (`copy`, `saxpy`), and the scaled one-pass
  `nrm2`. Where a loop is bound by the FP pipes instead (`asum`), the two
  tie. Zig adds nothing that clang on the same source would not.

## Setup

- **Machine, build and peak:** as [`results-m6.md`](results-m6.md#setup):
  Apple M6 Mac mini, macOS 27.0.1 (26A434); Zig 0.18.0-dev.267+de94864b9,
  ReleaseFast, CPU model `apple_m1`, `-Dblas-cores=all`; rstd detects
  `vortexm4`. Peak `fmla` 149.5 GFLOP/s single and 74.7 double, the clock
  4.78 to 4.79 GHz at every check of the single-threaded runs.
- **OpenBLAS:** 0.3.33 from nixpkgs' binary cache
  (`/nix/store/s2xzgzc6dbj1bgz5qf9331zha3g0w1w3-openblas-0.3.33`), built
  `TARGET=VORTEX`, `DYNAMIC_ARCH`, 32-bit integers. `openblas_get_config()`:
  `OpenBLAS 0.3.33 DYNAMIC_ARCH NO_AFFINITY neoversen1 MAX_THREADS=64`;
  `openblas_get_corename()`: `neoversen1`, and `openblas_get_num_threads()`
  1. A `DYNAMIC_ARCH` build defines `gotoblas_VORTEX` as
  `gotoblas_NEOVERSEN1`, so `neoversen1` is how it names the vortex core:
  the same kernels (`KERNEL.VORTEX` includes `KERNEL.NEOVERSEN1`), with
  N1's blocking. This binary maps `VORTEXM4` there too:
  `OPENBLAS_CORETYPE=vortexm4` still reports `neoversen1`, with the same
  gemm rates, so its SME gemm path never runs.
- **Accelerate:** the system framework, Accelerate 1.11 (vecLib 3.11) of
  macOS 27.0.1, through its LP64 `cblas_*$NEWLAPACK` entries.
  `BLASSetThreading(BLAS_THREADING_SINGLE_THREADED)` on the calling thread;
  `BLASGetThreading()` reports single-threaded. A `lib2` run took 29.33 s
  of user time in 29.56 s, so nothing ran on a second core.
- **Threads:** the bench sets `OPENBLAS_NUM_THREADS=1` and
  `VECLIB_MAXIMUM_THREADS=1` before it `dlopen`s the libraries, except for
  the one labeled default-threading gemm run. rstd's routines ran on
  `Io.failing`, one thread, in every run.

## Method

As [`results-m6.md`](results-m6.md#method): the same `race`, sizes, data,
correctness checks, timing (best of 15 batches of at least 10 ms; gemm 7 of
20 ms) and core placement. No contestant failed a check. What differs:

- **Libraries:** both are `dlopen`ed by absolute path (`std.DynLib`) and
  called through `cblas_*`, column-major, with the bench's strides. Each
  call is wrapped in the kernel ABI of the vtable slot it stands for, so
  the checks and calls are the ones the kernels get. `i?amax` adds one to
  cblas's 0-based index; complex `dot` is `cblas_?dotu_sub`. `gemv` and
  `gemm` pass `beta` 1, so every contestant does `y += alpha*op(A)*x` and
  `C += A*B`.
- **Repeatability:** two full runs, plus the gemm tables a third and fourth
  time (`OPENBLAS_CORETYPE=vortexm4`, and an earlier `lib3`). Between the
  full runs the median cell moved 0.2%, and 7 of 977 moved more than 4%:
  x4 `copy` with 1 MiB vectors (149 to 225 GB/s), x4 `dasum` from memory,
  and Accelerate's `copy` from memory (5%). The earlier `lib3` ran 3 to 5%
  faster for all three NEON contestants alike, their ratios within 0.02.
  The tables are the first full run.

### Contestants

- **openblas:** OpenBLAS's `cblas_*`. Ratios are over it.
- **rstd:** rstd's public routine (`blas.dot`, `blas.gemv`, `blas.gemm`
  and the rest), on the table rstd selects. `rstd runs` says whether that slot holds a
  converted kernel (`asm`) or the portable loop. `gemv` and `gemm` get
  malloc's allocator; `gemm` takes its packing space from it every call.
- **asm:** the converted kernel alone, as `results-m6.md`'s asm column:
  vortex's own, else armv8's stand-in (`from`), else none.
- **accelerate:** Accelerate's `cblas_*`.
- **vec x4:** `vector_kernels.zig`'s 4-wide `@Vector` kernel.
- **ob .S** (gemm only): OpenBLAS's vortex `sgemm_kernel_16x4.S` or
  `dgemm_kernel_8x4.S` with their copies, through the bench's copy of
  `level3.c`'s loops with vortex's `P`, `Q`, `R`.

## Level 1

GB/s; ratios are over openblas.

### Single (`s`)

| routine | p | n | inc | rstd runs | from | openblas | rstd | asm | accelerate | vec x4 | rstd/ob | asm/ob | accel/ob | x4/ob |
| ------- | - | -: | --: | --------- | ---- | -------: | ---: | --: | ---------: | -----: | ------: | -----: | -------: | ----: |
| axpy | s | 1024 | 1 | asm | vortex | 209.9 | 211.6 | 215.0 | 545.3 | 292.2 | 1.01 | 1.02 | 2.60 | 1.39 |
| axpy | s | 8192 | 1 | asm | vortex | 213.9 | 213.9 | 215.2 | 740.3 | 293.3 | 1.00 | 1.01 | 3.46 | 1.37 |
| axpy | s | 262144 | 1 | asm | vortex | 158.5 | 158.6 | 158.4 | 742.5 | 170.2 | 1.00 | 1.00 | 4.68 | 1.07 |
| axpy | s | 16777216 | 1 | asm | vortex | 140.2 | 140.6 | 140.4 | 144.6 | 137.3 | 1.00 | 1.00 | 1.03 | 0.98 |
| axpy | s | 4096 | 2 | asm | vortex | 28.7 | 28.7 | 28.7 | 55.8 | 56.0 | 1.00 | 1.00 | 1.95 | 1.95 |
| dot | s | 1024 | 1 | portable | none | 163.0 | 17.9 | - | 202.4 | 215.8 | 0.11 | - | 1.24 | 1.32 |
| dot | s | 8192 | 1 | portable | none | 176.3 | 16.4 | - | 502.4 | 220.7 | 0.09 | - | 2.85 | 1.25 |
| dot | s | 262144 | 1 | portable | none | 148.5 | 16.2 | - | 1039.1 | 150.0 | 0.11 | - | 7.00 | 1.01 |
| dot | s | 16777216 | 1 | portable | none | 125.6 | 16.1 | - | 154.1 | 127.3 | 0.13 | - | 1.23 | 1.01 |
| dot | s | 4096 | 2 | portable | none | 11.9 | 16.6 | - | 42.4 | 42.0 | 1.40 | - | 3.57 | 3.54 |
| scal | s | 1024 | 1 | asm | vortex | 224.1 | 223.9 | 232.2 | 382.7 | 275.0 | 1.00 | 1.04 | 1.71 | 1.23 |
| scal | s | 8192 | 1 | asm | vortex | 229.0 | 229.0 | 230.0 | 571.8 | 288.5 | 1.00 | 1.00 | 2.50 | 1.26 |
| scal | s | 262144 | 1 | asm | vortex | 162.6 | 162.4 | 162.2 | 571.5 | 162.5 | 1.00 | 1.00 | 3.52 | 1.00 |
| scal | s | 16777216 | 1 | asm | vortex | 137.8 | 137.7 | 137.7 | 161.8 | 134.8 | 1.00 | 1.00 | 1.17 | 0.98 |
| scal | s | 4096 | 2 | asm | vortex | 35.8 | 35.8 | 35.9 | 37.8 | 37.9 | 1.00 | 1.00 | 1.05 | 1.06 |
| nrm2 | s | 1024 | 1 | asm | vortex | 5.1 | 5.1 | 5.1 | 53.1 | 157.6 | 1.00 | 1.00 | 10.46 | 31.05 |
| nrm2 | s | 8192 | 1 | asm | vortex | 4.8 | 4.8 | 4.8 | 41.5 | 160.1 | 1.00 | 1.00 | 8.63 | 33.27 |
| nrm2 | s | 262144 | 1 | asm | vortex | 4.8 | 4.8 | 4.8 | 39.3 | 149.8 | 1.00 | 1.00 | 8.21 | 31.34 |
| nrm2 | s | 16777216 | 1 | asm | vortex | 4.8 | 4.8 | 4.8 | 40.1 | 111.8 | 1.00 | 1.00 | 8.38 | 23.37 |
| nrm2 | s | 4096 | 2 | asm | vortex | 4.8 | 4.8 | 4.8 | 19.8 | 21.9 | 1.00 | 1.00 | 4.10 | 4.54 |
| asum | s | 1024 | 1 | portable | armv8 | 144.5 | 10.4 | 87.4 | 134.6 | 138.7 | 0.07 | 0.60 | 0.93 | 0.96 |
| asum | s | 8192 | 1 | portable | armv8 | 147.1 | 8.8 | 61.6 | 176.7 | 148.3 | 0.06 | 0.42 | 1.20 | 1.01 |
| asum | s | 262144 | 1 | portable | armv8 | 146.4 | 8.6 | 58.3 | 273.9 | 143.5 | 0.06 | 0.40 | 1.87 | 0.98 |
| asum | s | 16777216 | 1 | portable | armv8 | 107.4 | 8.6 | 58.1 | 163.3 | 111.9 | 0.08 | 0.54 | 1.52 | 1.04 |
| asum | s | 4096 | 2 | portable | armv8 | 8.5 | 9.0 | 4.9 | 30.5 | 28.6 | 1.06 | 0.57 | 3.60 | 3.37 |
| iamax | s | 1024 | 1 | portable | armv8 | 80.3 | 10.0 | 69.6 | 60.2 | 112.2 | 0.12 | 0.87 | 0.75 | 1.40 |
| iamax | s | 8192 | 1 | portable | armv8 | 127.3 | 9.6 | 73.9 | 62.5 | 133.5 | 0.08 | 0.58 | 0.49 | 1.05 |
| iamax | s | 262144 | 1 | portable | armv8 | 141.5 | 9.6 | 74.1 | 62.2 | 141.4 | 0.07 | 0.52 | 0.44 | 1.00 |
| iamax | s | 16777216 | 1 | portable | armv8 | 107.5 | 9.6 | 73.3 | 61.7 | 107.7 | 0.09 | 0.68 | 0.57 | 1.00 |
| iamax | s | 4096 | 2 | portable | armv8 | 9.7 | 9.7 | 4.8 | 4.8 | 38.7 | 1.00 | 0.49 | 0.50 | 3.99 |
| copy | s | 1024 | 1 | portable | armv8 | 136.6 | 65.9 | 139.9 | 433.3 | 274.4 | 0.48 | 1.02 | 3.17 | 2.01 |
| copy | s | 8192 | 1 | portable | armv8 | 149.4 | 69.1 | 150.2 | 571.9 | 291.9 | 0.46 | 1.01 | 3.83 | 1.95 |
| copy | s | 262144 | 1 | portable | armv8 | 144.5 | 67.7 | 144.0 | 571.8 | 167.5 | 0.47 | 1.00 | 3.96 | 1.16 |
| copy | s | 16777216 | 1 | portable | armv8 | 120.7 | 66.2 | 120.5 | 96.2 | 119.8 | 0.55 | 1.00 | 0.80 | 0.99 |
| copy | s | 4096 | 2 | portable | armv8 | 37.4 | 68.7 | 37.3 | 37.7 | 37.7 | 1.83 | 1.00 | 1.01 | 1.01 |
| swap | s | 1024 | 1 | asm | vortex | 183.9 | 183.7 | 187.3 | 569.5 | 275.0 | 1.00 | 1.02 | 3.10 | 1.50 |
| swap | s | 8192 | 1 | asm | vortex | 195.0 | 194.9 | 195.5 | 568.5 | 294.3 | 1.00 | 1.00 | 2.92 | 1.51 |
| swap | s | 262144 | 1 | asm | vortex | 103.6 | 101.4 | 101.1 | 566.9 | 163.9 | 0.98 | 0.98 | 5.47 | 1.58 |
| swap | s | 16777216 | 1 | asm | vortex | 134.2 | 134.2 | 134.8 | 145.8 | 137.7 | 1.00 | 1.00 | 1.09 | 1.03 |
| swap | s | 4096 | 2 | asm | vortex | 67.3 | 67.2 | 67.3 | 61.0 | 67.2 | 1.00 | 1.00 | 0.91 | 1.00 |
| rot | s | 1024 | 1 | asm | vortex | 237.4 | 236.1 | 241.6 | 548.4 | 228.3 | 0.99 | 1.02 | 2.31 | 0.96 |
| rot | s | 8192 | 1 | asm | vortex | 249.7 | 249.9 | 250.1 | 571.6 | 227.4 | 1.00 | 1.00 | 2.29 | 0.91 |
| rot | s | 262144 | 1 | asm | vortex | 153.1 | 154.5 | 153.0 | 571.2 | 158.3 | 1.01 | 1.00 | 3.73 | 1.03 |
| rot | s | 16777216 | 1 | asm | vortex | 140.5 | 141.9 | 143.3 | 145.2 | 137.7 | 1.01 | 1.02 | 1.03 | 0.98 |
| rot | s | 4096 | 2 | asm | vortex | 18.9 | 19.1 | 19.1 | 56.0 | 64.3 | 1.01 | 1.01 | 2.97 | 3.40 |

### Double (`d`)

| routine | p | n | inc | rstd runs | from | openblas | rstd | asm | accelerate | vec x4 | rstd/ob | asm/ob | accel/ob | x4/ob |
| ------- | - | -: | --: | --------- | ---- | -------: | ---: | --: | ---------: | -----: | ------: | -----: | -------: | ----: |
| axpy | d | 512 | 1 | asm | vortex | 280.9 | 279.1 | 287.4 | 545.5 | 294.6 | 0.99 | 1.02 | 1.94 | 1.05 |
| axpy | d | 4096 | 1 | asm | vortex | 251.9 | 250.4 | 253.9 | 740.8 | 291.0 | 0.99 | 1.01 | 2.94 | 1.16 |
| axpy | d | 131072 | 1 | asm | vortex | 165.5 | 165.4 | 165.3 | 743.7 | 167.1 | 1.00 | 1.00 | 4.49 | 1.01 |
| axpy | d | 8388608 | 1 | asm | vortex | 149.5 | 149.2 | 149.8 | 147.6 | 144.2 | 1.00 | 1.00 | 0.99 | 0.96 |
| axpy | d | 2048 | 2 | asm | vortex | 55.9 | 55.9 | 55.9 | 107.8 | 106.4 | 1.00 | 1.00 | 1.93 | 1.90 |
| dot | d | 512 | 1 | portable | armv8 | 158.9 | 48.0 | 93.8 | 204.4 | 212.2 | 0.30 | 0.59 | 1.29 | 1.34 |
| dot | d | 4096 | 1 | portable | armv8 | 172.4 | 34.9 | 65.2 | 506.9 | 218.2 | 0.20 | 0.38 | 2.94 | 1.27 |
| dot | d | 131072 | 1 | portable | armv8 | 143.7 | 33.5 | 62.8 | 1039.7 | 146.1 | 0.23 | 0.44 | 7.24 | 1.02 |
| dot | d | 8388608 | 1 | portable | armv8 | 123.7 | 33.6 | 61.1 | 156.5 | 126.4 | 0.27 | 0.49 | 1.27 | 1.02 |
| dot | d | 2048 | 2 | portable | armv8 | 24.4 | 36.3 | 23.9 | 85.2 | 85.6 | 1.49 | 0.98 | 3.49 | 3.51 |
| scal | d | 512 | 1 | asm | vortex | 247.0 | 247.0 | 255.4 | 382.7 | 274.0 | 1.00 | 1.03 | 1.55 | 1.11 |
| scal | d | 4096 | 1 | asm | vortex | 255.2 | 255.3 | 256.7 | 571.4 | 288.0 | 1.00 | 1.01 | 2.24 | 1.13 |
| scal | d | 131072 | 1 | asm | vortex | 161.9 | 161.8 | 161.8 | 571.5 | 160.1 | 1.00 | 1.00 | 3.53 | 0.99 |
| scal | d | 8388608 | 1 | asm | vortex | 138.6 | 138.7 | 138.6 | 159.2 | 135.2 | 1.00 | 1.00 | 1.15 | 0.98 |
| scal | d | 2048 | 2 | asm | vortex | 70.6 | 70.4 | 70.0 | 72.3 | 74.3 | 1.00 | 0.99 | 1.02 | 1.05 |
| nrm2 | d | 512 | 1 | asm | vortex | 10.8 | 10.7 | 10.8 | 48.2 | 171.8 | 0.99 | 1.00 | 4.48 | 15.96 |
| nrm2 | d | 4096 | 1 | asm | vortex | 9.7 | 9.7 | 9.7 | 41.5 | 159.8 | 1.00 | 1.00 | 4.28 | 16.47 |
| nrm2 | d | 131072 | 1 | asm | vortex | 9.6 | 9.6 | 9.6 | 40.6 | 148.3 | 1.00 | 1.00 | 4.25 | 15.51 |
| nrm2 | d | 8388608 | 1 | asm | vortex | 9.6 | 9.6 | 9.6 | 40.4 | 107.7 | 1.00 | 1.00 | 4.22 | 11.26 |
| nrm2 | d | 2048 | 2 | asm | vortex | 9.7 | 9.7 | 9.7 | 21.6 | 45.7 | 1.00 | 1.00 | 2.22 | 4.70 |
| asum | d | 512 | 1 | portable | armv8 | 144.1 | 26.6 | 93.2 | 133.3 | 143.0 | 0.18 | 0.65 | 0.92 | 0.99 |
| asum | d | 4096 | 1 | portable | armv8 | 146.7 | 18.0 | 68.0 | 177.7 | 147.6 | 0.12 | 0.46 | 1.21 | 1.01 |
| asum | d | 131072 | 1 | portable | armv8 | 145.3 | 17.2 | 62.4 | 274.5 | 146.4 | 0.12 | 0.43 | 1.89 | 1.01 |
| asum | d | 8388608 | 1 | portable | armv8 | 106.2 | 17.2 | 61.7 | 164.5 | 108.9 | 0.16 | 0.58 | 1.55 | 1.03 |
| asum | d | 2048 | 2 | portable | armv8 | 17.7 | 18.8 | 9.8 | 62.1 | 56.7 | 1.06 | 0.56 | 3.51 | 3.20 |
| iamax | d | 512 | 1 | portable | armv8 | 109.8 | 20.9 | 123.2 | 60.2 | 138.4 | 0.19 | 1.12 | 0.55 | 1.26 |
| iamax | d | 4096 | 1 | portable | armv8 | 131.7 | 19.1 | 129.5 | 62.0 | 139.2 | 0.14 | 0.98 | 0.47 | 1.06 |
| iamax | d | 131072 | 1 | portable | armv8 | 137.3 | 18.6 | 130.1 | 61.8 | 140.9 | 0.14 | 0.95 | 0.45 | 1.03 |
| iamax | d | 8388608 | 1 | portable | armv8 | 103.9 | 18.7 | 99.2 | 61.4 | 104.9 | 0.18 | 0.95 | 0.59 | 1.01 |
| iamax | d | 2048 | 2 | portable | armv8 | 15.7 | 19.6 | 9.6 | 9.8 | 94.4 | 1.25 | 0.61 | 0.62 | 6.02 |
| copy | d | 512 | 1 | portable | armv8 | 136.8 | 132.2 | 159.0 | 436.0 | 273.0 | 0.97 | 1.16 | 3.19 | 1.99 |
| copy | d | 4096 | 1 | portable | armv8 | 145.9 | 137.4 | 148.8 | 571.9 | 289.6 | 0.94 | 1.02 | 3.92 | 1.99 |
| copy | d | 131072 | 1 | portable | armv8 | 140.6 | 109.7 | 141.8 | 571.8 | 201.3 | 0.78 | 1.01 | 4.07 | 1.43 |
| copy | d | 8388608 | 1 | portable | armv8 | 121.4 | 93.1 | 121.8 | 98.2 | 121.2 | 0.77 | 1.00 | 0.81 | 1.00 |
| copy | d | 2048 | 2 | portable | armv8 | 73.1 | 135.6 | 72.7 | 72.4 | 73.9 | 1.85 | 0.99 | 0.99 | 1.01 |
| swap | d | 512 | 1 | asm | vortex | 183.3 | 182.5 | 187.0 | 569.5 | 273.7 | 1.00 | 1.02 | 3.11 | 1.49 |
| swap | d | 4096 | 1 | asm | vortex | 194.4 | 194.2 | 194.9 | 568.4 | 292.8 | 1.00 | 1.00 | 2.92 | 1.51 |
| swap | d | 131072 | 1 | asm | vortex | 100.9 | 102.3 | 101.6 | 567.8 | 163.9 | 1.01 | 1.01 | 5.63 | 1.62 |
| swap | d | 8388608 | 1 | asm | vortex | 135.9 | 136.0 | 135.9 | 147.7 | 141.6 | 1.00 | 1.00 | 1.09 | 1.04 |
| swap | d | 2048 | 2 | asm | vortex | 127.5 | 127.2 | 132.8 | 120.4 | 127.9 | 1.00 | 1.04 | 0.94 | 1.00 |
| rot | d | 512 | 1 | asm | vortex | 229.3 | 227.6 | 232.5 | 548.4 | 228.0 | 0.99 | 1.01 | 2.39 | 0.99 |
| rot | d | 4096 | 1 | asm | vortex | 232.5 | 232.3 | 233.0 | 571.6 | 228.1 | 1.00 | 1.00 | 2.46 | 0.98 |
| rot | d | 131072 | 1 | asm | vortex | 159.1 | 158.9 | 158.9 | 571.2 | 157.5 | 1.00 | 1.00 | 3.59 | 0.99 |
| rot | d | 8388608 | 1 | asm | vortex | 144.9 | 147.0 | 149.2 | 147.1 | 140.6 | 1.01 | 1.03 | 1.02 | 0.97 |
| rot | d | 2048 | 2 | asm | vortex | 36.8 | 36.9 | 37.0 | 119.0 | 126.8 | 1.00 | 1.00 | 3.23 | 3.44 |

### Single complex (`c`)

| routine | p | n | inc | rstd runs | from | openblas | rstd | asm | accelerate | vec x4 | rstd/ob | asm/ob | accel/ob | x4/ob |
| ------- | - | -: | --: | --------- | ---- | -------: | ---: | --: | ---------: | -----: | ------: | -----: | -------: | ----: |
| axpy | c | 512 | 1 | asm | vortex | 141.5 | 141.2 | 144.4 | 360.3 | 283.6 | 1.00 | 1.02 | 2.55 | 2.00 |
| axpy | c | 4096 | 1 | asm | vortex | 145.1 | 145.2 | 145.6 | 374.3 | 284.0 | 1.00 | 1.00 | 2.58 | 1.96 |
| axpy | c | 131072 | 1 | asm | vortex | 146.3 | 146.3 | 146.4 | 376.4 | 165.9 | 1.00 | 1.00 | 2.57 | 1.13 |
| axpy | c | 8388608 | 1 | asm | vortex | 140.0 | 140.7 | 141.2 | 145.0 | 141.1 | 1.01 | 1.01 | 1.04 | 1.01 |
| axpy | c | 2048 | 2 | asm | vortex | 108.6 | 108.9 | 109.4 | 89.7 | 107.8 | 1.00 | 1.01 | 0.83 | 0.99 |
| dot | c | 512 | 1 | portable | none | 101.6 | 28.0 | - | 97.7 | 149.5 | 0.28 | - | 0.96 | 1.47 |
| dot | c | 4096 | 1 | portable | none | 105.8 | 26.4 | - | 79.9 | 185.1 | 0.25 | - | 0.76 | 1.75 |
| dot | c | 131072 | 1 | portable | none | 106.6 | 26.0 | - | 77.4 | 141.8 | 0.24 | - | 0.73 | 1.33 |
| dot | c | 8388608 | 1 | portable | none | 106.0 | 26.0 | - | 77.0 | 121.3 | 0.25 | - | 0.73 | 1.14 |
| dot | c | 2048 | 2 | portable | none | 9.9 | 26.9 | - | 20.6 | 41.8 | 2.72 | - | 2.09 | 4.24 |
| scal | c | 512 | 1 | asm | vortex | 94.7 | 95.0 | 96.2 | 253.4 | 190.3 | 1.00 | 1.02 | 2.68 | 2.01 |
| scal | c | 4096 | 1 | asm | vortex | 97.1 | 97.2 | 97.3 | 263.6 | 191.8 | 1.00 | 1.00 | 2.72 | 1.98 |
| scal | c | 131072 | 1 | asm | vortex | 97.7 | 97.8 | 97.8 | 265.0 | 157.7 | 1.00 | 1.00 | 2.71 | 1.61 |
| scal | c | 8388608 | 1 | asm | vortex | 97.3 | 97.4 | 97.4 | 157.7 | 135.4 | 1.00 | 1.00 | 1.62 | 1.39 |
| scal | c | 2048 | 2 | asm | vortex | 41.9 | 42.0 | 42.0 | 71.0 | 71.3 | 1.00 | 1.00 | 1.69 | 1.70 |
| nrm2 | c | 512 | 1 | asm | vortex | 5.1 | 5.1 | 5.1 | 50.0 | 154.2 | 1.00 | 1.00 | 9.85 | 30.36 |
| nrm2 | c | 4096 | 1 | asm | vortex | 4.8 | 4.8 | 4.8 | 40.9 | 157.0 | 1.00 | 1.00 | 8.49 | 32.59 |
| nrm2 | c | 131072 | 1 | asm | vortex | 4.8 | 4.8 | 4.8 | 39.8 | 148.6 | 1.00 | 1.00 | 8.33 | 31.08 |
| nrm2 | c | 8388608 | 1 | asm | vortex | 4.8 | 4.8 | 4.8 | 39.5 | 108.8 | 1.00 | 1.00 | 8.26 | 22.78 |
| nrm2 | c | 2048 | 2 | asm | vortex | 4.8 | 4.8 | 4.8 | 11.5 | 23.6 | 1.00 | 1.00 | 2.36 | 4.87 |
| asum | c | 512 | 1 | portable | armv8 | 143.6 | 9.6 | 92.5 | 133.3 | 137.2 | 0.07 | 0.64 | 0.93 | 0.96 |
| asum | c | 4096 | 1 | portable | armv8 | 146.3 | 8.6 | 66.3 | 176.7 | 146.6 | 0.06 | 0.45 | 1.21 | 1.00 |
| asum | c | 131072 | 1 | portable | armv8 | 145.9 | 8.3 | 60.3 | 273.9 | 142.4 | 0.06 | 0.41 | 1.88 | 0.98 |
| asum | c | 8388608 | 1 | portable | armv8 | 106.8 | 8.2 | 61.5 | 164.3 | 109.8 | 0.08 | 0.58 | 1.54 | 1.03 |
| asum | c | 2048 | 2 | portable | armv8 | 15.8 | 8.6 | 16.6 | 16.0 | 28.9 | 0.55 | 1.05 | 1.02 | 1.83 |
| iamax | c | 512 | 1 | portable | armv8 | 100.7 | 20.4 | 119.6 | 36.6 | 86.7 | 0.20 | 1.19 | 0.36 | 0.86 |
| iamax | c | 4096 | 1 | portable | armv8 | 134.4 | 18.8 | 129.2 | 36.7 | 145.6 | 0.14 | 0.96 | 0.27 | 1.08 |
| iamax | c | 131072 | 1 | portable | armv8 | 141.8 | 18.5 | 130.0 | 36.7 | 144.7 | 0.13 | 0.92 | 0.26 | 1.02 |
| iamax | c | 8388608 | 1 | portable | armv8 | 104.6 | 18.4 | 99.3 | 36.6 | 107.2 | 0.18 | 0.95 | 0.35 | 1.02 |
| iamax | c | 2048 | 2 | portable | armv8 | 19.1 | 19.1 | 19.1 | 19.6 | 48.2 | 1.00 | 1.00 | 1.03 | 2.52 |
| copy | c | 512 | 1 | portable | armv8 | 134.6 | 63.2 | 157.4 | 436.0 | 272.3 | 0.47 | 1.17 | 3.24 | 2.02 |
| copy | c | 4096 | 1 | portable | armv8 | 145.4 | 69.4 | 147.3 | 571.9 | 289.2 | 0.48 | 1.01 | 3.93 | 1.99 |
| copy | c | 131072 | 1 | portable | armv8 | 140.4 | 69.8 | 141.5 | 571.8 | 174.5 | 0.50 | 1.01 | 4.07 | 1.24 |
| copy | c | 8388608 | 1 | portable | armv8 | 121.7 | 68.4 | 122.2 | 96.4 | 121.2 | 0.56 | 1.00 | 0.79 | 1.00 |
| copy | c | 2048 | 2 | portable | armv8 | 73.1 | 68.4 | 73.4 | 137.1 | 134.8 | 0.94 | 1.00 | 1.88 | 1.84 |
| swap | c | 512 | 1 | asm | vortex | 182.4 | 182.7 | 186.7 | 569.2 | 272.8 | 1.00 | 1.02 | 3.12 | 1.50 |
| swap | c | 4096 | 1 | asm | vortex | 194.0 | 194.2 | 194.7 | 568.4 | 292.0 | 1.00 | 1.00 | 2.93 | 1.50 |
| swap | c | 131072 | 1 | asm | vortex | 101.6 | 102.4 | 101.5 | 567.3 | 163.2 | 1.01 | 1.00 | 5.58 | 1.61 |
| swap | c | 8388608 | 1 | asm | vortex | 136.3 | 136.4 | 136.1 | 146.8 | 140.5 | 1.00 | 1.00 | 1.08 | 1.03 |
| swap | c | 2048 | 2 | asm | vortex | 127.2 | 127.3 | 131.5 | 125.3 | 126.9 | 1.00 | 1.03 | 0.99 | 1.00 |
| rot | c | 512 | 1 | asm | vortex | 263.1 | 262.3 | 264.7 | 548.4 | 228.1 | 1.00 | 1.01 | 2.08 | 0.87 |
| rot | c | 4096 | 1 | asm | vortex | 274.2 | 269.1 | 274.4 | 571.6 | 226.7 | 0.98 | 1.00 | 2.08 | 0.83 |
| rot | c | 131072 | 1 | asm | vortex | 158.9 | 158.9 | 158.9 | 571.2 | 157.8 | 1.00 | 1.00 | 3.60 | 0.99 |
| rot | c | 8388608 | 1 | asm | vortex | 145.8 | 147.4 | 148.2 | 146.4 | 139.4 | 1.01 | 1.02 | 1.00 | 0.96 |
| rot | c | 2048 | 2 | asm | vortex | 126.2 | 127.4 | 126.2 | 111.3 | 113.4 | 1.01 | 1.00 | 0.88 | 0.90 |

### Double complex (`z`)

| routine | p | n | inc | rstd runs | from | openblas | rstd | asm | accelerate | vec x4 | rstd/ob | asm/ob | accel/ob | x4/ob |
| ------- | - | -: | --: | --------- | ---- | -------: | ---: | --: | ---------: | -----: | ------: | -----: | -------: | ----: |
| axpy | z | 256 | 1 | asm | vortex | 141.9 | 141.2 | 144.4 | 360.5 | 279.6 | 1.00 | 1.02 | 2.54 | 1.97 |
| axpy | z | 2048 | 1 | asm | vortex | 145.1 | 145.2 | 145.6 | 374.5 | 280.2 | 1.00 | 1.00 | 2.58 | 1.93 |
| axpy | z | 65536 | 1 | asm | vortex | 146.6 | 146.6 | 146.6 | 376.5 | 165.9 | 1.00 | 1.00 | 2.57 | 1.13 |
| axpy | z | 4194304 | 1 | asm | vortex | 141.4 | 142.0 | 142.4 | 147.0 | 142.8 | 1.00 | 1.01 | 1.04 | 1.01 |
| axpy | z | 1024 | 2 | asm | vortex | 216.0 | 215.5 | 217.6 | 175.8 | 209.4 | 1.00 | 1.01 | 0.81 | 0.97 |
| dot | z | 256 | 1 | portable | none | 101.5 | 61.5 | - | 166.4 | 174.5 | 0.61 | - | 1.64 | 1.72 |
| dot | z | 2048 | 1 | portable | none | 105.4 | 53.0 | - | 157.4 | 189.0 | 0.50 | - | 1.49 | 1.79 |
| dot | z | 65536 | 1 | portable | none | 106.6 | 52.0 | - | 143.6 | 141.9 | 0.49 | - | 1.35 | 1.33 |
| dot | z | 4194304 | 1 | portable | none | 106.0 | 51.9 | - | 123.0 | 121.3 | 0.49 | - | 1.16 | 1.14 |
| dot | z | 1024 | 2 | portable | none | 20.3 | 53.9 | - | 42.7 | 84.7 | 2.66 | - | 2.10 | 4.18 |
| scal | z | 256 | 1 | asm | vortex | 94.6 | 94.8 | 96.0 | 253.6 | 190.9 | 1.00 | 1.02 | 2.68 | 2.02 |
| scal | z | 2048 | 1 | asm | vortex | 96.9 | 96.9 | 97.0 | 263.7 | 190.9 | 1.00 | 1.00 | 2.72 | 1.97 |
| scal | z | 65536 | 1 | asm | vortex | 97.7 | 97.7 | 97.8 | 265.1 | 157.6 | 1.00 | 1.00 | 2.71 | 1.61 |
| scal | z | 4194304 | 1 | asm | vortex | 97.3 | 97.4 | 97.4 | 153.9 | 135.0 | 1.00 | 1.00 | 1.58 | 1.39 |
| scal | z | 1024 | 2 | asm | vortex | 83.1 | 83.5 | 83.5 | 139.1 | 140.0 | 1.00 | 1.01 | 1.67 | 1.69 |
| nrm2 | z | 256 | 1 | asm | vortex | 10.8 | 10.7 | 10.8 | 11.8 | 169.5 | 0.99 | 1.00 | 1.09 | 15.65 |
| nrm2 | z | 2048 | 1 | asm | vortex | 9.7 | 9.7 | 9.7 | 10.7 | 159.4 | 1.00 | 1.00 | 1.10 | 16.44 |
| nrm2 | z | 65536 | 1 | asm | vortex | 9.6 | 9.6 | 9.6 | 10.6 | 148.0 | 1.00 | 1.00 | 1.10 | 15.49 |
| nrm2 | z | 4194304 | 1 | asm | vortex | 9.6 | 9.6 | 9.6 | 10.5 | 109.1 | 1.00 | 1.00 | 1.10 | 11.41 |
| nrm2 | z | 1024 | 2 | asm | vortex | 9.8 | 9.8 | 9.8 | 10.8 | 49.1 | 1.00 | 1.00 | 1.11 | 5.01 |
| asum | z | 256 | 1 | portable | armv8 | 144.1 | 23.5 | 92.2 | 133.2 | 142.7 | 0.16 | 0.64 | 0.92 | 0.99 |
| asum | z | 2048 | 1 | portable | armv8 | 146.3 | 17.2 | 64.9 | 177.7 | 146.8 | 0.12 | 0.44 | 1.21 | 1.00 |
| asum | z | 65536 | 1 | portable | armv8 | 146.5 | 16.4 | 62.9 | 274.5 | 141.6 | 0.11 | 0.43 | 1.87 | 0.97 |
| asum | z | 4194304 | 1 | portable | armv8 | 103.7 | 16.4 | 61.9 | 163.3 | 106.9 | 0.16 | 0.60 | 1.57 | 1.03 |
| asum | z | 1024 | 2 | portable | armv8 | 36.1 | 17.7 | 37.4 | 33.2 | 63.1 | 0.49 | 1.04 | 0.92 | 1.75 |
| iamax | z | 256 | 1 | portable | armv8 | 117.1 | 45.4 | 128.9 | 65.9 | 109.7 | 0.39 | 1.10 | 0.56 | 0.94 |
| iamax | z | 2048 | 1 | portable | armv8 | 132.7 | 37.9 | 137.1 | 50.8 | 144.4 | 0.29 | 1.03 | 0.38 | 1.09 |
| iamax | z | 65536 | 1 | portable | armv8 | 133.9 | 36.9 | 137.4 | 49.1 | 144.7 | 0.28 | 1.03 | 0.37 | 1.08 |
| iamax | z | 4194304 | 1 | portable | armv8 | 104.2 | 36.7 | 100.4 | 48.9 | 104.4 | 0.35 | 0.96 | 0.47 | 1.00 |
| iamax | z | 1024 | 2 | portable | armv8 | 38.9 | 38.9 | 39.1 | 40.1 | 89.9 | 1.00 | 1.01 | 1.03 | 2.31 |
| copy | z | 256 | 1 | portable | armv8 | 133.5 | 126.7 | 154.2 | 435.9 | 273.9 | 0.95 | 1.15 | 3.27 | 2.05 |
| copy | z | 2048 | 1 | portable | armv8 | 145.0 | 137.0 | 148.2 | 571.9 | 288.8 | 0.95 | 1.02 | 3.94 | 1.99 |
| copy | z | 65536 | 1 | portable | armv8 | 139.9 | 136.7 | 139.6 | 571.8 | 174.3 | 0.98 | 1.00 | 4.09 | 1.25 |
| copy | z | 4194304 | 1 | portable | armv8 | 121.1 | 117.1 | 121.6 | 98.3 | 121.1 | 0.97 | 1.00 | 0.81 | 1.00 |
| copy | z | 1024 | 2 | portable | armv8 | 145.3 | 133.5 | 146.5 | 232.2 | 229.0 | 0.92 | 1.01 | 1.60 | 1.58 |
| swap | z | 256 | 1 | asm | vortex | 182.7 | 182.5 | 186.3 | 569.5 | 272.6 | 1.00 | 1.02 | 3.12 | 1.49 |
| swap | z | 2048 | 1 | asm | vortex | 193.8 | 193.8 | 194.3 | 568.7 | 291.5 | 1.00 | 1.00 | 2.93 | 1.50 |
| swap | z | 65536 | 1 | asm | vortex | 100.6 | 101.6 | 101.0 | 567.8 | 161.3 | 1.01 | 1.00 | 5.64 | 1.60 |
| swap | z | 4194304 | 1 | asm | vortex | 136.7 | 136.5 | 136.0 | 147.6 | 141.4 | 1.00 | 0.99 | 1.08 | 1.03 |
| swap | z | 1024 | 2 | asm | vortex | 223.0 | 222.8 | 236.7 | 228.8 | 217.8 | 1.00 | 1.06 | 1.03 | 0.98 |
| rot | z | 256 | 1 | asm | vortex | 265.6 | 261.0 | 265.4 | 548.4 | 227.9 | 0.98 | 1.00 | 2.06 | 0.86 |
| rot | z | 2048 | 1 | asm | vortex | 270.9 | 274.0 | 273.5 | 571.6 | 227.8 | 1.01 | 1.01 | 2.11 | 0.84 |
| rot | z | 65536 | 1 | asm | vortex | 158.8 | 158.9 | 158.9 | 571.2 | 157.4 | 1.00 | 1.00 | 3.60 | 0.99 |
| rot | z | 4194304 | 1 | asm | vortex | 146.6 | 148.4 | 149.2 | 147.3 | 140.6 | 1.01 | 1.02 | 1.00 | 0.96 |
| rot | z | 1024 | 2 | asm | vortex | 228.7 | 238.7 | 226.5 | 219.5 | 220.6 | 1.04 | 0.99 | 0.96 | 0.96 |

## Level 2

GFLOP/s, square, unit strides, `beta` 1.

| routine | p | n | rstd runs | from | openblas | rstd | asm | accelerate | vec x4 | rstd/ob | asm/ob | accel/ob | x4/ob |
| ------- | - | -: | --------- | ---- | -------: | ---: | --: | ---------: | -----: | ------: | -----: | -------: | ----: |
| gemv_n | s | 90 | portable | armv8 | 74.5 | 8.7 | 36.3 | 56.3 | 61.3 | 0.12 | 0.49 | 0.76 | 0.82 |
| gemv_n | s | 1024 | portable | armv8 | 53.8 | 8.8 | 38.0 | 476.0 | 52.5 | 0.16 | 0.71 | 8.85 | 0.97 |
| gemv_n | s | 4096 | portable | armv8 | 40.3 | 9.0 | 38.1 | 68.7 | 38.7 | 0.22 | 0.95 | 1.71 | 0.96 |
| gemv_t | s | 90 | asm | vortex | 37.3 | 37.9 | 38.2 | 71.8 | 41.3 | 1.02 | 1.02 | 1.92 | 1.11 |
| gemv_t | s | 1024 | asm | vortex | 40.3 | 40.1 | 40.1 | 266.7 | 55.3 | 1.00 | 1.00 | 6.62 | 1.37 |
| gemv_t | s | 4096 | asm | vortex | 38.9 | 38.9 | 38.9 | 67.6 | 39.5 | 1.00 | 1.00 | 1.74 | 1.01 |
| gemv_n | d | 64 | asm | vortex | 19.9 | 20.3 | 20.5 | 46.2 | 36.5 | 1.02 | 1.03 | 2.32 | 1.83 |
| gemv_n | d | 724 | asm | vortex | 19.0 | 18.9 | 18.9 | 198.6 | 28.9 | 1.00 | 1.00 | 10.47 | 1.52 |
| gemv_n | d | 2896 | asm | vortex | 18.9 | 18.9 | 19.0 | 36.5 | 18.8 | 1.00 | 1.00 | 1.93 | 0.99 |
| gemv_t | d | 64 | asm | vortex | 19.6 | 20.0 | 20.1 | 41.5 | 35.7 | 1.02 | 1.03 | 2.12 | 1.82 |
| gemv_t | d | 724 | asm | vortex | 20.1 | 20.0 | 20.1 | 172.7 | 29.8 | 1.00 | 1.00 | 8.60 | 1.48 |
| gemv_t | d | 2896 | asm | vortex | 19.5 | 19.5 | 19.5 | 38.4 | 18.8 | 1.00 | 1.00 | 1.97 | 0.97 |
| gemv_n | c | 64 | asm | vortex | 34.7 | 34.7 | 34.8 | 127.2 | 85.8 | 1.00 | 1.01 | 3.67 | 2.48 |
| gemv_n | c | 724 | asm | vortex | 35.8 | 35.8 | 35.8 | 494.9 | 91.5 | 1.00 | 1.00 | 13.83 | 2.56 |
| gemv_n | c | 2896 | asm | vortex | 36.4 | 36.4 | 36.4 | 131.2 | 62.2 | 1.00 | 1.00 | 3.61 | 1.71 |
| gemv_t | c | 64 | asm | vortex | 45.5 | 45.5 | 45.7 | 81.8 | 51.1 | 1.00 | 1.01 | 1.80 | 1.12 |
| gemv_t | c | 724 | asm | vortex | 25.8 | 25.9 | 26.0 | 472.5 | 95.9 | 1.00 | 1.01 | 18.29 | 3.71 |
| gemv_t | c | 2896 | asm | vortex | 22.8 | 22.8 | 23.1 | 128.7 | 68.3 | 1.00 | 1.01 | 5.63 | 2.99 |
| gemv_n | z | 45 | asm | vortex | 23.2 | 23.2 | 23.4 | 44.5 | 42.3 | 1.00 | 1.01 | 1.91 | 1.82 |
| gemv_n | z | 512 | asm | vortex | 24.3 | 24.3 | 24.3 | 266.6 | 46.0 | 1.00 | 1.00 | 10.95 | 1.89 |
| gemv_n | z | 2048 | asm | vortex | 24.2 | 24.2 | 24.2 | 43.5 | 35.9 | 1.00 | 1.00 | 1.80 | 1.48 |
| gemv_t | z | 45 | asm | vortex | 27.5 | 27.5 | 27.7 | 39.2 | 32.1 | 1.00 | 1.01 | 1.43 | 1.17 |
| gemv_t | z | 512 | asm | vortex | 22.8 | 22.9 | 22.8 | 250.0 | 44.8 | 1.00 | 1.00 | 10.97 | 1.97 |
| gemv_t | z | 2048 | asm | vortex | 21.3 | 21.3 | 21.3 | 58.5 | 35.0 | 1.00 | 1.00 | 2.75 | 1.65 |

## Level 3

gemm `NN`, square, `C += A*B`, GFLOP/s. rstd has no complex gemm yet.

### One thread

| p | n | openblas | rstd | ob .S | accelerate | rstd/ob | ob .S/ob | accel/ob |
| - | -: | -------: | ---: | ----: | ---------: | ------: | -------: | -------: |
| s | 64 | 115.9 | 111.8 | 115.4 | 1035.8 | 0.96 | 1.00 | 8.94 |
| s | 128 | 123.3 | 113.5 | 122.8 | 1311.2 | 0.92 | 1.00 | 10.63 |
| s | 256 | 123.8 | 118.8 | 123.2 | 1684.1 | 0.96 | 1.00 | 13.60 |
| s | 512 | 126.3 | 120.7 | 126.0 | 1764.0 | 0.96 | 1.00 | 13.97 |
| s | 1024 | 127.7 | 122.2 | 127.4 | 1923.3 | 0.96 | 1.00 | 15.07 |
| s | 2048 | 127.6 | 121.6 | 127.3 | 1837.5 | 0.95 | 1.00 | 14.40 |
| d | 64 | 60.0 | 56.2 | 60.1 | 367.1 | 0.94 | 1.00 | 6.12 |
| d | 128 | 60.6 | 58.6 | 60.8 | 458.0 | 0.97 | 1.00 | 7.55 |
| d | 256 | 61.4 | 59.8 | 61.9 | 478.0 | 0.97 | 1.01 | 7.78 |
| d | 512 | 63.5 | 60.5 | 63.4 | 497.1 | 0.95 | 1.00 | 7.83 |
| d | 1024 | 64.0 | 60.9 | 63.8 | 510.6 | 0.95 | 1.00 | 7.98 |
| d | 2048 | 64.1 | 61.0 | 64.0 | 448.0 | 0.95 | 1.00 | 6.99 |
| c | 64 | 109.5 | - | - | 987.8 | - | - | 9.02 |
| c | 128 | 111.5 | - | - | 1414.5 | - | - | 12.69 |
| c | 256 | 113.3 | - | - | 1625.6 | - | - | 14.35 |
| c | 512 | 114.6 | - | - | 1842.9 | - | - | 16.08 |
| c | 1024 | 114.9 | - | - | 1944.0 | - | - | 16.92 |
| c | 2048 | 115.1 | - | - | 1801.5 | - | - | 15.65 |
| z | 64 | 54.2 | - | - | 356.3 | - | - | 6.57 |
| z | 128 | 54.9 | - | - | 405.9 | - | - | 7.40 |
| z | 256 | 56.2 | - | - | 463.3 | - | - | 8.25 |
| z | 512 | 56.8 | - | - | 510.9 | - | - | 9.00 |
| z | 1024 | 56.6 | - | - | 460.2 | - | - | 8.13 |
| z | 2048 | 56.9 | - | - | 460.2 | - | - | 8.09 |

### Default threading

Both libraries as they come: OpenBLAS with 12 threads, Accelerate
multi-threaded. rstd (`Io.failing`) and `ob .S` stay single-threaded. The clock read 4.47
GHz at the start and end of this run.

| p | n | openblas | rstd | ob .S | accelerate | rstd/ob | ob .S/ob | accel/ob |
| - | -: | -------: | ---: | ----: | ---------: | ------: | -------: | -------: |
| s | 64 | 115.9 | 111.8 | 115.0 | 1036.7 | 0.96 | 0.99 | 8.94 |
| s | 128 | 324.7 | 111.6 | 119.8 | 1294.1 | 0.34 | 0.37 | 3.99 |
| s | 256 | 556.5 | 116.7 | 120.0 | 1663.1 | 0.21 | 0.22 | 2.99 |
| s | 512 | 642.3 | 118.4 | 123.0 | 1732.1 | 0.18 | 0.19 | 2.70 |
| s | 1024 | 552.2 | 119.8 | 124.4 | 1982.3 | 0.22 | 0.23 | 3.59 |
| s | 2048 | 719.2 | 119.1 | 124.8 | 2108.2 | 0.17 | 0.17 | 2.93 |
| d | 64 | 60.0 | 56.1 | 60.2 | 366.5 | 0.93 | 1.00 | 6.11 |
| d | 128 | 209.6 | 57.2 | 58.6 | 450.3 | 0.27 | 0.28 | 2.15 |
| d | 256 | 250.8 | 58.1 | 59.4 | 455.9 | 0.23 | 0.24 | 1.82 |
| d | 512 | 317.7 | 58.2 | 60.6 | 490.5 | 0.18 | 0.19 | 1.54 |
| d | 1024 | 291.8 | 58.7 | 60.0 | 516.4 | 0.20 | 0.21 | 1.77 |
| d | 2048 | 354.6 | 59.5 | 64.0 | 534.4 | 0.17 | 0.18 | 1.51 |
| c | 64 | 270.3 | - | - | 975.2 | - | - | 3.61 |
| c | 128 | 426.0 | - | - | 1397.7 | - | - | 3.28 |
| c | 256 | 546.8 | - | - | 1589.5 | - | - | 2.91 |
| c | 512 | 603.7 | - | - | 1812.0 | - | - | 3.00 |
| c | 1024 | 615.7 | - | - | 2089.4 | - | - | 3.39 |
| c | 2048 | 652.1 | - | - | 2074.3 | - | - | 3.18 |
| z | 64 | 143.1 | - | - | 350.9 | - | - | 2.45 |
| z | 128 | 254.9 | - | - | 400.5 | - | - | 1.57 |
| z | 256 | 273.4 | - | - | 456.2 | - | - | 1.67 |
| z | 512 | 302.9 | - | - | 500.4 | - | - | 1.65 |
| z | 1024 | 313.7 | - | - | 537.9 | - | - | 1.71 |
| z | 2048 | 325.1 | - | - | 536.0 | - | - | 1.65 |

Accelerate gains at most 1.2 times from more threads: its single thread
already drives the SME unit.

## Hot loops

The x4 kernels' unit-stride loops as Zig emits them (`zig build-obj
-OReleaseFast -target aarch64-macos -femit-asm`), against the OpenBLAS
source of the kernel the library runs.

| routine | x4: per iteration | OpenBLAS: per iteration | x4/ob in L1 |
| ------- | ----------------- | ----------------------- | ----------: |
| `sdot` | 32 floats: 8 `ldp q`, 8 `fmla .4s` into 8 accumulators | `dot_kernel_asimd.c` (inline asm), per 32 floats: the same 8 `ldp q` and 8 `fmla` into 8 accumulators, plus 4 `prfm` | 1.25 to 1.32 |
| `saxpy` | 16 floats: 4 `ldp q`, 4 `fmla` by element, 2 `stp q` | `axpy.S`, 8 floats: 2 `ld1` of two `q`, 2 `fmla`, 1 `st1` of two `q`, 2 `prfm` | 1.37 to 1.39 |
| `scopy` | 16 floats: 2 `ldp q`, 2 `stp q` | `copy_thunderx2t99.c`, 4 floats: 1 `ldr q`, 1 `str q` | 1.95 to 2.01 |
| `sasum` | 32 floats: 4 `ldp q`, 8 `fabs`, 8 `fadd` into 8 accumulators | `sasum_thunderx2t99.c`, 64 floats: 16 `q` loads, 16 `fabs`, 16 `fadd`, 4 `prfm` | 0.96 to 1.01 |
| `sgemv_n` | 16 rows by 4 columns: 10 `ldp q` (2 of `y`), 16 `fmla` by element, 2 `stp q` | `sgemv_n_neon.c` (intrinsics), 8 rows by 4 columns: 10 `q` loads (2 of `y`), 8 FMAs, 2 stores | 0.82 (n 90) |

Every x4 loop addresses with immediate offsets from one pointer per
vector, bumped once per iteration; LLVM splits each `@Vector(16, f32)`
into four `q` registers, so the two accumulators become eight independent
chains. None prefetches. The numbers fit a load-slot budget:

- **`sdot`:** x4 reaches 221 GB/s, 46 bytes a cycle at 4.78 GHz, nearly
  three 16-byte loads. OpenBLAS's loop issues four `prfm` for every 16
  loads; if each takes a load slot, it gets 16/20 of that, 177 GB/s. It
  measures 176.
- **`saxpy`:** `axpy.S` spends 6 load-slot operations per 32 bytes of `x`
  where x4 spends 8 per 64, 1.5 times as many; x4 runs 1.37 to 1.39 times
  faster.
- **`scopy`:** one 16-byte load and store per loop iteration caps
  `copy_thunderx2t99.c` near one iteration a cycle; x4 moves 64 bytes.
- **`sasum`** is bound by the FP pipes, not loads: x4's two vectors a cycle
  need four FP operations (`fabs`, `fadd` each), the probe's four `fmla` a
  cycle. Both run at 31 bytes a cycle, prefetches and all.

The prefetch reading is an inference from these rates, not a measured
pipe count. `nrm2` and `i?amax` are different algorithms rather than
different code (`results-m6.md`'s contestants).

## Fairness

What could tilt these comparisons:

- **OpenBLAS's build:** one binary, nixpkgs' `DYNAMIC_ARCH` 0.3.33, on
  the core it selects. A build with the `vortexm4` SME gemm path that
  upstream's `interface/gemm.c` takes for clang would be much faster at
  gemm than the openblas column; it was not available to measure.
- **Call overhead:** every library call goes through a wrapper in the
  kernel ABI, one indirect call more than the asm column pays. asm/ob is
  0.98 to 1.06, highest at 4 KiB and strided `zswap`. OpenBLAS's `cblas_*`
  also does its own argument checks, as rstd's routines do.
- **rstd's allocations:** its `gemm` allocates its packing space with
  malloc on every call where OpenBLAS keeps a buffer pool. rstd/ob is
  lowest at `n` 128 (0.92); kernel against kernel was 0.95 to 0.98.
- **Accelerate's internals:** the SME attribution follows from rates
  past the NEON peak on one thread; the bench cannot see which units run.
- **Strided cases:** OpenBLAS's real `sdot` and `ddot` strided are slower
  than rstd's portable loop (11.9 and 24.4 GB/s against 16.6 and 36.3).
- **The rest:** as [`results-m6.md`](results-m6.md#fairness), notably
  the codegen swings of plain Zig loops between builds.
