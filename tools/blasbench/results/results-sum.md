# `sum`'s kernels on Zen 3 and the M6

Why the public `sum` ran 5 to 30 times slower per element than `asum` (T192, [#545](https://github.com/randyos-org/rstd/issues/545)), what replaced the kernels, and what is left for later tuning. Measured 2026-11-19 with `blas_bench`'s `l1` and `lib1` sections (`--routines=sum,asum`) and its `compensated` section, one thread.

## What ran before

- **x86_64, every core.** `ssum` and `dsum` are `kernel/arm/sum.c` upstream (`kernel/x86_64/KERNEL`), which the coverage page listed as plain C, so rstd ran its portable loop: one `f64` accumulator, 0.65 to 0.95 ns an element on Zen 3, bound by the add's latency. But `sum.c` is not plain C: at unit stride it runs `kernel/simd/intrin.h`'s vectors, four sums of one vector each in the element's type (`dsum` from AVX2 up). The survey's text search does not see those intrinsics, as the file names no intrinsics header itself. `csum` and `zsum` (`zsum_sse.S`, `zsum_sse2.S`, which peel) already ran `vector.zig`'s `sum`.
- **arm64 on Apple's systems.** `registry.zig` left every slot that is not routed to Accelerate there. Accelerate has no `sum` (it is OpenBLAS's extension), so `ssum`, `dsum` and `zsum`, C upstream on `vortex` and `vortexm4` (`?sum_thunderx2t99.c`), ran the portable loop, 0.3 to 0.5 ns an element on the M6, and `csum` the converted `csum.S`. Off Apple the same cores already ran `vector.zig`'s `sum` in those three slots.
- **`asum`** runs `vector.zig`'s kernel on Zen (its `.S` peels) and Accelerate or the converted kernels on the M6: hence the gap.
- **Elsewhere** (not measured, from [`coverage.md`](../../../src/math/kernels/blas/blas1/sum/coverage.md)): 32-bit x86 runs the converted x87 `sum.S`; arm64's other cores the converted `sum.S`, `csum.S` and `zsum.S`, or `vector.zig`'s where upstream is C; ARMv6 and ARMv7 `sum_vfp.S`; POWER, SPARC, LoongArch's LSX and LASX cores and Loongson their converted assembly; RISC-V's vector cores `vector.zig`'s. The portable loop is left where upstream is `sum.c` too: ARMv5, MIPS, the scalar RISC-V and LoongArch cores, z/Architecture's generic core and WebAssembly (there `sum.c` is `wasm_simd128` vectors upstream, the same gap as x86_64's was), and Z13 and Z14, whose kernels are not ported (the page's six gaps).

## What runs now

`vector.zig`'s `sum`, four native vectors a step (`Real(F, 4).sum`, `Complex(F, 4).sum`), in every `s`, `d`, `c` and `z` sum slot of x86_64's cores (`blas1/sum/x86_64/kernel.zig`) and of `vortex` and `vortexm4` on every system (`registry.unaccelerated`; `csum` left `registry.displacing`'s `keep`). The kernel itself did not change, so the cores that already ran it keep their bits. Its order is fixed by `n` and the stride alone (no peeling, unaligned loads), `f32` sums in `f64` lanes and rounds once, and the public routine folds its chunks in order as before: one core's bits at any worker count.

## Timings

GB/s of `x` read, best of 15 batches of at least 10 ms; a cell with two numbers is two runs. `n` is 4 KiB and 32 KiB of `x` (L1), 1 MiB (L2), 64 MiB (memory), and a strided case over 32 KiB. Neither machine was idle: other builds ran on both all day. On Zen 3 the bench's clock read 4.49 GHz, its idle rate, only in the first run's kernel section, and 2.6 to 4.1 GHz in the rest, so those rates are low by up to a third; the ratios, which agree between the runs, are what to read.

### Zen 3, the kernels

Ryzen 7 7730U (Zen 3, core `zen`), Windows 11, Zig 0.18.0-dev.1029, ReleaseFast. `OpenBLAS` is its own kernel for ZEN compiled from the pin (`-Dopenblas`, v0.3.34-794-gdbf2c24c9: `../arm/sum.c` with AVX2, `zsum_sse.S`, `zsum_sse2.S`; `asum_sse.S`, `asum_sse2.S`, `zasum_sse.S`, `zasum_sse2.S`), `new` the selected core's slot, `old` the portable loop that `ssum` and `dsum` ran (`csum` and `zsum` already ran `new`).

| routine | p | n | inc | OpenBLAS | new | old | new/OpenBLAS | new/old |
| ------- | - | -: | --: | -------: | --: | --: | -----------: | ------: |
| asum | s | 1024 | 1 | 84.6 / 65.4 | 107.1 / 82.6 | 5.7 / 4.6 | 1.27 / 1.26 | 18.8 / 18.0 |
| asum | s | 8192 | 1 | 82.6 / 68.3 | 205.3 / 167.9 | 5.6 / 4.6 | 2.49 / 2.46 | 36.7 / 36.5 |
| asum | s | 262144 | 1 | 68.6 / 62.0 | 93.6 / 78.1 | 5.2 / 4.4 | 1.36 / 1.26 | 18.0 / 17.7 |
| asum | s | 16777216 | 1 | 22.2 / 20.6 | 26.3 / 24.2 | 5.5 / 4.7 | 1.18 / 1.17 | 4.8 / 5.1 |
| asum | s | 4096 | 2 | 17.1 / 15.5 | 22.9 / 19.8 | 5.8 / 5.2 | 1.34 / 1.28 | 3.9 / 3.8 |
| sum | s | 1024 | 1 | 162.6 / 152.0 | 57.6 / 52.9 | 5.3 / 4.9 | 0.35 / 0.35 | 10.9 / 10.8 |
| sum | s | 8192 | 1 | 170.2 / 146.4 | 67.3 / 57.5 | 5.8 / 5.0 | 0.40 / 0.39 | 11.6 / 11.5 |
| sum | s | 262144 | 1 | 97.7 / 82.7 | 59.0 / 49.9 | 5.3 / 4.4 | 0.60 / 0.60 | 11.1 / 11.3 |
| sum | s | 16777216 | 1 | 26.6 / 23.2 | 22.9 / 20.1 | 5.5 / 4.1 | 0.86 / 0.87 | 4.2 / 4.9 |
| sum | s | 4096 | 2 | 5.9 / 4.7 | 23.2 / 19.1 | 5.9 / 5.0 | 3.93 / 4.06 | 3.9 / 3.8 |
| asum | d | 512 | 1 | 85.3 / 60.3 | 146.2 / 108.0 | 12.3 / 8.8 | 1.71 / 1.79 | 11.9 / 12.3 |
| asum | d | 4096 | 1 | 74.5 / 74.3 | 191.3 / 192.9 | 10.1 / 10.1 | 2.57 / 2.60 | 18.9 / 19.1 |
| asum | d | 131072 | 1 | 68.8 / 64.7 | 91.3 / 83.5 | 10.1 / 9.7 | 1.33 / 1.29 | 9.0 / 8.6 |
| asum | d | 8388608 | 1 | 21.3 / 21.3 | 26.1 / 26.1 | 9.5 / 8.9 | 1.23 / 1.23 | 2.7 / 2.9 |
| asum | d | 2048 | 2 | 33.1 / 30.3 | 44.1 / 40.1 | 11.6 / 10.5 | 1.33 / 1.32 | 3.8 / 3.8 |
| sum | d | 512 | 1 | 185.1 / 165.3 | 190.4 / 173.7 | 12.8 / 11.7 | 1.03 / 1.05 | 14.9 / 14.8 |
| sum | d | 4096 | 1 | 167.3 / 148.1 | 223.2 / 190.4 | 11.4 / 10.4 | 1.33 / 1.29 | 19.6 / 18.3 |
| sum | d | 131072 | 1 | 95.5 / 80.2 | 90.3 / 81.2 | 10.0 / 9.2 | 0.95 / 1.01 | 9.0 / 8.8 |
| sum | d | 8388608 | 1 | 21.6 / 24.2 | 22.8 / 24.2 | 8.4 / 8.5 | 1.06 / 1.00 | 2.7 / 2.8 |
| sum | d | 2048 | 2 | 12.1 / 10.7 | 47.4 / 42.5 | 12.3 / 10.9 | 3.92 / 3.97 | 3.9 / 3.9 |
| asum | c | 512 | 1 | 84.5 / 75.3 | 106.0 / 90.1 | 5.9 / 5.3 | 1.25 / 1.20 | 18.0 / 17.0 |
| asum | c | 4096 | 1 | 79.7 / 74.7 | 196.3 / 183.3 | 5.5 / 5.0 | 2.46 / 2.45 | 35.7 / 36.7 |
| asum | c | 131072 | 1 | 59.3 / 59.7 | 83.0 / 82.5 | 4.6 / 4.7 | 1.40 / 1.38 | 18.0 / 17.6 |
| asum | c | 8388608 | 1 | 20.5 / 21.1 | 25.0 / 24.6 | 4.9 / 4.7 | 1.22 / 1.17 | 5.1 / 5.2 |
| asum | c | 2048 | 2 | 20.1 / 29.7 | 11.3 / 19.9 | 3.4 / 5.0 | 0.56 / 0.67 | 3.3 / 4.0 |
| sum | c | 512 | 1 | 52.1 / 65.8 | 39.5 / 52.2 | 3.7 / 4.9 | 0.76 / 0.79 | 10.7 / 10.7 |
| sum | c | 4096 | 1 | 56.2 / 73.7 | 47.9 / 62.3 | 4.0 / 5.3 | 0.85 / 0.85 | 12.0 / 11.8 |
| sum | c | 131072 | 1 | 52.6 / 59.6 | 44.1 / 49.2 | 3.8 / 4.3 | 0.84 / 0.83 | 11.6 / 11.4 |
| sum | c | 8388608 | 1 | 20.1 / 20.5 | 20.1 / 19.9 | 3.8 / 4.2 | 1.00 / 0.97 | 5.3 / 4.7 |
| sum | c | 2048 | 2 | 23.9 / 27.5 | 16.2 / 18.4 | 4.2 / 4.7 | 0.68 / 0.67 | 3.9 / 3.9 |
| asum | z | 256 | 1 | 56.2 / 77.1 | 95.1 / 126.1 | 7.7 / 10.8 | 1.69 / 1.64 | 12.4 / 11.7 |
| asum | z | 2048 | 1 | 57.3 / 75.0 | 158.0 / 195.2 | 7.9 / 10.1 | 2.76 / 2.60 | 20.0 / 19.3 |
| asum | z | 65536 | 1 | 50.1 / 54.4 | 68.1 / 75.7 | 7.4 / 8.3 | 1.36 / 1.39 | 9.2 / 9.1 |
| asum | z | 4194304 | 1 | 18.4 / 18.6 | 21.7 / 22.1 | 6.5 / 7.7 | 1.18 / 1.19 | 3.3 / 2.9 |
| asum | z | 1024 | 2 | 38.7 / 54.2 | 27.5 / 40.4 | 7.3 / 10.5 | 0.71 / 0.75 | 3.8 / 3.8 |
| sum | z | 256 | 1 | 54.5 / 75.6 | 125.1 / 172.1 | 8.3 / 11.7 | 2.30 / 2.28 | 15.1 / 14.7 |
| sum | z | 2048 | 1 | 53.0 / 67.4 | 151.7 / 196.6 | 7.5 / 9.9 | 2.86 / 2.92 | 20.2 / 19.9 |
| sum | z | 65536 | 1 | 49.4 / 63.7 | 65.9 / 86.8 | 7.2 / 9.5 | 1.33 / 1.36 | 9.2 / 9.1 |
| sum | z | 4194304 | 1 | 19.9 / 21.7 | 23.1 / 24.2 | 6.7 / 9.0 | 1.16 / 1.12 | 3.4 / 2.7 |
| sum | z | 1024 | 2 | 42.3 / 57.0 | 32.6 / 43.3 | 8.2 / 10.9 | 0.77 / 0.76 | 4.0 / 4.0 |

### Zen 3, the public routine

`blas.sum` and `blas.asum` on `Io.Threaded` with one thread against the official OpenBLAS 0.3.34 DLL's `cblas_?sum` and `cblas_?asum`.

| routine | p | n | inc | OpenBLAS | rstd | rstd/OpenBLAS |
| ------- | - | -: | --: | -------: | ---: | ------------: |
| asum | s | 1024 | 1 | 53.9 / 74.1 | 62.5 / 84.9 | 1.16 / 1.15 |
| asum | s | 8192 | 1 | 55.7 / 75.4 | 130.1 / 178.6 | 2.34 / 2.37 |
| asum | s | 262144 | 1 | 46.2 / 60.9 | 62.2 / 81.1 | 1.35 / 1.33 |
| asum | s | 16777216 | 1 | 18.3 / 21.2 | 22.6 / 25.7 | 1.23 / 1.21 |
| asum | s | 4096 | 2 | 10.6 / 15.4 | 14.3 / 20.1 | 1.35 / 1.31 |
| sum | s | 1024 | 1 | 75.8 / 109.1 | 37.6 / 54.2 | 0.50 / 0.50 |
| sum | s | 8192 | 1 | 88.8 / 114.8 | 46.1 / 59.5 | 0.52 / 0.52 |
| sum | s | 262144 | 1 | 58.5 / 75.5 | 40.3 / 53.6 | 0.69 / 0.71 |
| sum | s | 16777216 | 1 | 19.3 / 20.3 | 19.0 / 19.9 | 0.98 / 0.98 |
| sum | s | 4096 | 2 | 3.8 / 4.1 | 15.0 / 16.2 | 3.95 / 3.95 |
| asum | d | 512 | 1 | 55.7 / 70.1 | 81.8 / 103.7 | 1.47 / 1.48 |
| asum | d | 4096 | 1 | 53.4 / 74.2 | 132.2 / 187.3 | 2.48 / 2.52 |
| asum | d | 131072 | 1 | 52.0 / 61.2 | 69.2 / 83.4 | 1.33 / 1.36 |
| asum | d | 8388608 | 1 | 18.1 / 21.1 | 21.8 / 24.8 | 1.20 / 1.18 |
| asum | d | 2048 | 2 | 22.1 / 30.6 | 29.3 / 40.1 | 1.33 / 1.31 |
| sum | d | 512 | 1 | 83.8 / 84.8 | 105.6 / 108.9 | 1.26 / 1.28 |
| sum | d | 4096 | 1 | 87.8 / 113.2 | 151.0 / 185.2 | 1.72 / 1.64 |
| sum | d | 131072 | 1 | 59.9 / 72.4 | 67.1 / 84.0 | 1.12 / 1.16 |
| sum | d | 8388608 | 1 | 19.9 / 21.2 | 22.5 / 26.7 | 1.13 / 1.26 |
| sum | d | 2048 | 2 | 7.8 / 10.4 | 30.4 / 40.4 | 3.90 / 3.88 |
| asum | c | 512 | 1 | 56.3 / 73.9 | 61.2 / 78.7 | 1.09 / 1.06 |
| asum | c | 4096 | 1 | 55.7 / 67.8 | 121.5 / 148.3 | 2.18 / 2.19 |
| asum | c | 131072 | 1 | 47.5 / 61.4 | 62.4 / 80.3 | 1.31 / 1.31 |
| asum | c | 8388608 | 1 | 18.9 / 20.7 | 23.0 / 24.5 | 1.22 / 1.18 |
| asum | c | 2048 | 2 | 23.5 / 29.9 | 15.3 / 20.0 | 0.65 / 0.67 |
| sum | c | 512 | 1 | 57.8 / 70.1 | 37.9 / 46.4 | 0.66 / 0.66 |
| sum | c | 4096 | 1 | 55.0 / 73.4 | 44.5 / 60.1 | 0.81 / 0.82 |
| sum | c | 131072 | 1 | 50.9 / 65.3 | 40.9 / 54.0 | 0.80 / 0.83 |
| sum | c | 8388608 | 1 | 20.1 / 20.2 | 19.5 / 20.1 | 0.97 / 1.00 |
| sum | c | 2048 | 2 | 22.3 / 31.1 | 15.1 / 21.2 | 0.68 / 0.68 |
| asum | z | 256 | 1 | 58.7 / 68.0 | 77.6 / 90.7 | 1.32 / 1.33 |
| asum | z | 2048 | 1 | 57.9 / 68.8 | 133.9 / 164.1 | 2.31 / 2.39 |
| asum | z | 65536 | 1 | 47.5 / 55.9 | 64.6 / 74.5 | 1.36 / 1.33 |
| asum | z | 4194304 | 1 | 17.7 / 19.4 | 21.1 / 22.4 | 1.19 / 1.15 |
| asum | z | 1024 | 2 | 35.1 / 50.0 | 26.8 / 36.1 | 0.76 / 0.72 |
| sum | z | 256 | 1 | 52.0 / 71.0 | 79.6 / 109.9 | 1.53 / 1.55 |
| sum | z | 2048 | 1 | 54.4 / 65.6 | 142.0 / 175.3 | 2.61 / 2.67 |
| sum | z | 65536 | 1 | 51.8 / 64.1 | 66.9 / 82.9 | 1.29 / 1.29 |
| sum | z | 4194304 | 1 | 18.9 / 21.9 | 22.0 / 26.1 | 1.16 / 1.19 |
| sum | z | 1024 | 2 | 37.4 / 58.9 | 26.6 / 41.6 | 0.71 / 0.71 |

The public `sum` before and after, ns an element (`compensated` section, the `sum` and `asum` columns; `before` is the tree without the change):

| type | n | `sum` before | `sum` after | `asum` | before/after |
| ---- | -: | -----------: | ----------: | -----: | -----------: |
| f32 | 100 | 0.647 / 0.665 | 0.236 / 0.214 | 0.473 / 0.428 | 2.7 / 3.1 |
| f32 | 1000 | 0.784 / 0.735 | 0.096 / 0.084 | 0.087 / 0.083 | 8.2 / 8.8 |
| f32 | 10000 | 1.047 / 0.882 | 0.071 / 0.077 | 0.036 / 0.038 | 14.7 / 11.5 |
| f32 | 100000 | 0.808 / 0.880 | 0.067 / 0.069 | 0.036 / 0.036 | 12.1 / 12.8 |
| f32 | 1000000 | 0.769 / 0.788 | 0.081 / 0.080 | 0.050 / 0.048 | 9.5 / 9.8 |
| f32 | 10000000 | 0.835 / 0.847 | 0.223 / 0.189 | 0.187 / 0.157 | 3.7 / 4.5 |
| f64 | 100 | 0.532 / 0.559 | 0.203 / 0.217 | 0.211 / 0.230 | 2.6 / 2.6 |
| f64 | 1000 | 0.822 / 0.725 | 0.111 / 0.086 | 0.103 / 0.083 | 7.4 / 8.4 |
| f64 | 10000 | 0.739 / 0.750 | 0.073 / 0.080 | 0.072 / 0.078 | 10.1 / 9.4 |
| f64 | 100000 | 0.835 / 0.917 | 0.103 / 0.101 | 0.105 / 0.097 | 8.1 / 9.1 |
| f64 | 1000000 | 0.797 / 0.812 | 0.109 / 0.122 | 0.108 / 0.119 | 7.3 / 6.7 |
| f64 | 10000000 | 0.867 / 0.885 | 0.360 / 0.434 | 0.377 / 0.422 | 2.4 / 2.0 |

### M6, the kernels

Apple M6 Mac mini (`vortexm4`), macOS, Zig 0.18.0-dev.859, ReleaseFast, rstd's own kernels (`--rstd-accelerate=0`); one run, of the tree before the change, the clock at 4.78 GHz. `converted` is armv8's converted `sum.S`, `csum.S` and `zsum.S` (`csum.S` is also what `vortex` ran for `csum`), `old` the portable loop that `ssum`, `dsum` and `zsum` ran, `new` `vector.zig`'s kernel four vectors a step, the bench's own instance of what the slots now hold (`vec x4`).

| routine | p | n | inc | converted | old | new | new/converted | new/old |
| ------- | - | -: | --: | --------: | --: | --: | ------------: | ------: |
| asum | s | 1024 | 1 | 84.1 | 9.9 | 135.7 | 1.61 | 13.7 |
| asum | s | 8192 | 1 | 59.7 | 8.3 | 144.8 | 2.43 | 17.4 |
| asum | s | 262144 | 1 | 56.2 | 8.1 | 140.5 | 2.50 | 17.3 |
| asum | s | 16777216 | 1 | 56.3 | 8.1 | 112.0 | 1.99 | 13.8 |
| asum | s | 4096 | 2 | 4.9 | 8.5 | 27.7 | 5.65 | 3.3 |
| sum | s | 1024 | 1 | 93.1 | 10.3 | 71.3 | 0.77 | 6.9 |
| sum | s | 8192 | 1 | 65.6 | 8.8 | 73.5 | 1.12 | 8.4 |
| sum | s | 262144 | 1 | 59.7 | 8.6 | 73.5 | 1.23 | 8.5 |
| sum | s | 16777216 | 1 | 59.4 | 8.6 | 73.1 | 1.23 | 8.5 |
| sum | s | 4096 | 2 | 8.5 | 9.0 | 26.6 | 3.13 | 3.0 |
| asum | d | 512 | 1 | 92.7 | 26.7 | 142.7 | 1.54 | 5.3 |
| asum | d | 4096 | 1 | 66.2 | 18.0 | 146.3 | 2.21 | 8.1 |
| asum | d | 131072 | 1 | 61.1 | 17.3 | 134.7 | 2.20 | 7.8 |
| asum | d | 8388608 | 1 | 61.4 | 17.2 | 113.4 | 1.85 | 6.6 |
| asum | d | 2048 | 2 | 9.8 | 18.8 | 56.4 | 5.76 | 3.0 |
| sum | d | 512 | 1 | 104.7 | 28.5 | 210.2 | 2.01 | 7.4 |
| sum | d | 4096 | 1 | 68.4 | 18.8 | 212.6 | 3.11 | 11.3 |
| sum | d | 131072 | 1 | 62.7 | 18.0 | 139.4 | 2.22 | 7.7 |
| sum | d | 8388608 | 1 | 54.2 | 14.9 | 99.7 | 1.84 | 6.7 |
| sum | d | 2048 | 2 | 17.5 | 19.9 | 58.1 | 3.32 | 2.9 |
| asum | c | 512 | 1 | 92.6 | 9.4 | 136.4 | 1.47 | 14.5 |
| asum | c | 4096 | 1 | 66.2 | 8.1 | 145.2 | 2.19 | 17.9 |
| asum | c | 131072 | 1 | 63.6 | 7.9 | 140.8 | 2.21 | 17.8 |
| asum | c | 8388608 | 1 | 60.4 | 7.7 | 112.4 | 1.86 | 14.6 |
| asum | c | 2048 | 2 | 16.6 | 8.2 | 29.9 | 1.80 | 3.6 |
| sum | c | 512 | 1 | 103.8 | 9.7 | 71.8 | 0.69 | 7.4 |
| sum | c | 4096 | 1 | 66.1 | 8.3 | 73.0 | 1.10 | 8.8 |
| sum | c | 131072 | 1 | 59.0 | 7.6 | 68.8 | 1.17 | 9.1 |
| sum | c | 8388608 | 1 | 51.4 | 6.3 | 58.8 | 1.14 | 9.3 |
| sum | c | 2048 | 2 | 17.4 | 8.6 | 29.4 | 1.69 | 3.4 |
| asum | z | 256 | 1 | 90.7 | 23.1 | 142.3 | 1.57 | 6.2 |
| asum | z | 2048 | 1 | 66.1 | 17.3 | 147.5 | 2.23 | 8.5 |
| asum | z | 65536 | 1 | 62.9 | 15.7 | 133.8 | 2.13 | 8.5 |
| asum | z | 4194304 | 1 | 53.5 | 14.0 | 100.3 | 1.87 | 7.2 |
| asum | z | 1024 | 2 | 37.4 | 17.8 | 63.6 | 1.70 | 3.6 |
| sum | z | 256 | 1 | 106.0 | 25.6 | 207.4 | 1.96 | 8.1 |
| sum | z | 2048 | 1 | 69.1 | 18.4 | 213.4 | 3.09 | 11.6 |
| sum | z | 65536 | 1 | 63.6 | 17.7 | 139.3 | 2.19 | 7.9 |
| sum | z | 4194304 | 1 | 52.4 | 15.2 | 99.2 | 1.89 | 6.5 |
| sum | z | 1024 | 2 | 39.3 | 19.0 | 70.8 | 1.80 | 3.7 |

Not timed on the M6: the public routine after the change, and OpenBLAS's VORTEX kernels (`?sum_thunderx2t99.c`, C with inline assembly, through its library's `cblas_?sum`). The Mac was held by other work; `blas_bench --rstd-accelerate=0 --routines=sum,asum l1 lib1` gives both.

## The price of `f64` sums for `f32`

`ssum` and `csum` widen each `f32` to `f64` before adding, so the result rounds once. Both machines convert about four floats a cycle (`vcvtps2pd` on Zen 3, `fcvtl` on the M6's four vector pipes), which bounds those two near 0.055 ns an element in cache on the M6 and 0.06 to 0.07 on Zen 3, whatever the accumulators' shape; `f32` lanes run two to three times faster there. Shapes raced alone (scratch code, not in the tree), ns an element; `f64 2x16` is two accumulators of 16 `f64` lanes, the shape of `vector.zig`'s `ssum` on both machines, `f32 4x8` upstream `sum.c`'s with AVX2, `abs f32 2x32` `vector.zig`'s `sasum`:

M6:

| n | scalar | f64 2x16 | f64 4x8 | f64 8x8 | f64 4x16 | f64 8x4 | f64 4x4 | f32 4x8 | f32 8x8 | f32 4x16 | f32 2x32 | abs f32 2x32 |
| -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: |
| 1000 | 0.3876 | 0.0571 | 0.0558 | 0.0578 | 0.0584 | 0.0556 | 0.0552 | 0.0206 | 0.0207 | 0.0264 | 0.0347 | 0.0414 |
| 10000 | 0.4605 | 0.0547 | 0.0547 | 0.0560 | 0.0549 | 0.0547 | 0.0546 | 0.0188 | 0.0184 | 0.0184 | 0.0211 | 0.0290 |
| 100000 | 0.4679 | 0.0546 | 0.0546 | 0.0558 | 0.0547 | 0.0546 | 0.0546 | 0.0272 | 0.0267 | 0.0271 | 0.0274 | 0.0276 |
| 1000000 | 0.4686 | 0.0546 | 0.0546 | 0.0557 | 0.0547 | 0.0546 | 0.0546 | 0.0305 | 0.0308 | 0.0309 | 0.0299 | 0.0319 |
| 10000000 | 0.4507 | 0.0543 | 0.0546 | 0.0558 | 0.0547 | 0.0546 | 0.0546 | 0.0344 | 0.0343 | 0.0342 | 0.0338 | 0.0353 |

| n | scalar | f64 2x16 | f64 4x4 | f64 8x4 | f64 4x8 | f64 8x2 | f64 4x2 | abs f64 2x16 |
| -: | -: | -: | -: | -: | -: | -: | -: | -: |
| 1000 | 0.3728 | 0.0479 | 0.0382 | 0.0371 | 0.0372 | 0.0380 | 0.0549 | 0.0593 |
| 10000 | 0.4438 | 0.0387 | 0.0368 | 0.0367 | 0.0367 | 0.0367 | 0.0740 | 0.0547 |
| 100000 | 0.4602 | 0.0548 | 0.0545 | 0.0538 | 0.0542 | 0.0544 | 0.0760 | 0.0548 |
| 1000000 | 0.4643 | 0.0610 | 0.0623 | 0.0631 | 0.0631 | 0.0624 | 0.0763 | 0.0654 |
| 10000000 | 0.4459 | 0.0735 | 0.0734 | 0.0731 | 0.0733 | 0.0728 | 0.0807 | 0.0763 |

Zen 3, the machine shared (about half its idle rates; the ratios between shapes are the point):

| n | scalar | f64 2x16 | f64 4x8 | f64 8x8 | f64 4x16 | f64 8x4 | f64 4x4 | f32 4x8 | f32 8x8 | f32 4x16 | f32 2x32 | abs f32 2x32 |
| -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: |
| 1000 | 1.1566 | 0.1157 | 0.1034 | 0.1133 | 0.1239 | 0.1046 | 0.1050 | 0.0393 | 0.0333 | 0.0432 | 0.0537 | 0.0789 |
| 10000 | 1.1960 | 0.1004 | 0.1024 | 0.1008 | 0.1011 | 0.1026 | 0.1057 | 0.0533 | 0.0531 | 0.0519 | 0.0506 | 0.0560 |
| 100000 | 1.2462 | 0.1033 | 0.1029 | 0.1054 | 0.1073 | 0.1067 | 0.1077 | 0.0558 | 0.0564 | 0.0548 | 0.0562 | 0.0554 |
| 1000000 | 1.2199 | 0.1028 | 0.1054 | 0.1082 | 0.1087 | 0.1026 | 0.1051 | 0.0637 | 0.0640 | 0.0650 | 0.0642 | 0.0657 |
| 10000000 | 1.2004 | 0.2827 | 0.3010 | 0.3065 | 0.3161 | 0.2963 | 0.2895 | 0.2806 | 0.2928 | 0.2853 | 0.2856 | 0.2844 |

| n | scalar | f64 2x16 | f64 4x4 | f64 8x4 | f64 4x8 | f64 8x2 | f64 4x2 | abs f64 2x16 |
| -: | -: | -: | -: | -: | -: | -: | -: | -: |
| 1000 | 1.1037 | 0.0637 | 0.0762 | 0.0583 | 0.0589 | 0.1141 | 0.1586 | 0.0804 |
| 10000 | 1.1766 | 0.1056 | 0.1034 | 0.1008 | 0.1034 | 0.1336 | 0.1597 | 0.1015 |
| 100000 | 1.2110 | 0.1311 | 0.1333 | 0.1310 | 0.1327 | 0.1543 | 0.1821 | 0.1334 |
| 1000000 | 1.2650 | 0.1548 | 0.1643 | 0.1678 | 0.1497 | 0.1834 | 0.2351 | 0.1604 |
| 10000000 | 1.2037 | 0.5206 | 0.4761 | 0.4663 | 0.4834 | 0.5201 | 0.4689 | 0.4580 |

What `f32` lanes would cost in accuracy: the error of each shape in ulps of the exact `f32` sum, on uniform data (Zen 3; the shapes' arithmetic is the same on any machine):

`f32`, uniform in [-1,1):

| n | scalar | f64 2x16 | f64 4x8 | f64 8x8 | f64 4x16 | f64 8x4 | f64 4x4 | f32 4x8 | f32 8x8 | f32 4x16 | f32 2x32 |
| -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: |
| 1000 | 0.18 | 0.18 | 0.18 | 0.18 | 0.18 | 0.18 | 0.18 | 1.82 | 2.18 | 1.18 | 3.18 |
| 10000 | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 | 3.17 | 4.83 | 4.83 | 1.83 |
| 100000 | 0.23 | 0.23 | 0.23 | 0.23 | 0.23 | 0.23 | 0.23 | 1.77 | 5.23 | 7.23 | 6.23 |
| 1000000 | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 | 0.15 | 24.85 | 55.85 | 55.85 | 55.85 |
| 10000000 | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 | 186.11 | 80.11 | 80.11 | 80.11 |

`f32`, uniform in [0,1):

| n | scalar | f64 2x16 | f64 4x8 | f64 8x8 | f64 4x16 | f64 8x4 | f64 4x4 | f32 4x8 | f32 8x8 | f32 4x16 | f32 2x32 |
| -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: |
| 1000 | 0.38 | 0.38 | 0.38 | 0.38 | 0.38 | 0.38 | 0.38 | 0.38 | 0.38 | 0.62 | 0.38 |
| 10000 | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 | 1.17 | 1.17 | 0.83 |
| 100000 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 3.05 | 2.95 | 3.95 | 3.95 |
| 1000000 | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 | 0.11 | 0.89 | 0.89 | 2.89 | 1.89 |
| 10000000 | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 | 14.07 | 0.07 | 1.93 | 0.93 |

`f64`, uniform in [-1,1):

| n | scalar | f64 2x16 | f64 4x4 | f64 8x4 | f64 4x8 | f64 8x2 | f64 4x2 |
| -: | -: | -: | -: | -: | -: | -: | -: |
| 1000 | 14.56 | 0.56 | 4.56 | 0.44 | 0.56 | 4.56 | 3.44 |
| 10000 | 152.22 | 9.78 | 29.22 | 9.78 | 9.78 | 30.22 | 7.78 |
| 100000 | 69.45 | 5.45 | 0.55 | 5.45 | 3.45 | 0.55 | 23.55 |
| 1000000 | 201.05 | 50.05 | 47.05 | 51.05 | 50.05 | 47.05 | 45.05 |
| 10000000 | 585.54 | 86.54 | 191.54 | 88.54 | 88.54 | 191.54 | 133.46 |

`f64`, uniform in [0,1):

| n | scalar | f64 2x16 | f64 4x4 | f64 8x4 | f64 4x8 | f64 8x2 | f64 4x2 |
| -: | -: | -: | -: | -: | -: | -: | -: |
| 1000 | 4.93 | 0.07 | 1.07 | 1.07 | 1.07 | 0.93 | 1.07 |
| 10000 | 10.57 | 0.43 | 0.57 | 0.57 | 0.43 | 0.57 | 0.43 |
| 100000 | 135.59 | 0.41 | 0.59 | 0.41 | 0.41 | 1.59 | 4.41 |
| 1000000 | 198.77 | 4.23 | 1.23 | 5.23 | 6.23 | 1.23 | 0.77 |
| 10000000 | 426.11 | 16.89 | 60.11 | 15.89 | 15.89 | 60.11 | 54.89 |

Kept in `f64`, the contract rstd documents ([`comparisons.md`](../../../src/math/kernels/blas/comparisons.md)): `sum` is at or above `asum`'s rate in `d` and `z`; in `s` and `c` it is within the factor above in cache and equal in memory, and behind OpenBLAS's `f32` lanes on Zen 3 by the same factor (0.35 to 0.6 of `ssum` in cache, 0.76 to 0.85 of `csum`). Whether `ssum` and `csum` move to `f32` lanes is the user's call (#545). Switching `ssum` and `csum` to `f32` lanes is small (the accumulators' type in `Real.plainSum` and the strided loops, and the two tests that state the order and the single rounding); it would match OpenBLAS's `ssum` on Zen 3, at the error shown.

## For later tuning

- **`ssum` and `csum` in cache** are bound by the widening, not by the kernel's shape: see above. `sdot` has the same bound ([`comparisons.md`](../../../src/math/kernels/blas/comparisons.md), "`f64` `sdot` cost").
- **Short vectors.** `Real.plainSum` adds what is left of a `2 * AW` block in one scalar chain: up to 31 floats in `d` with AVX2. A one-vector loop before the scalar tail took 20% off an `f64` sum at `n = 1000` on the M6 in the scratch race (`f64 4x4` against `f64 2x16`, AVX2's shape); it changes the stated order, and so the bits, on every core that runs the kernel.
- **Complex strided `sum` and `asum` on Zen 3** run at 0.55 to 0.8 of `zsum_sse*.S`'s and `zasum_sse*.S`'s rate: two 2-lane sums an element pair. Not touched here.
- **`csum` on the M6 at 4 KiB** is the one case the kernel `vortex` ran before won (`csum.S`, `f32` lanes: 1.4 times `new`); `new` is ahead from 32 KiB on and strided, and sums in `f64`. armv8's `sum.S`, which `vortex` never ran, is likewise 1.3 times `new` in `s` at 4 KiB and behind it from 32 KiB on.
- **Not measured:** the x86_64 cores but Zen 3 (the same kernel, on each core's own vector width; upstream's `sum.c` is one file for all of them), AVX-512 cores in particular, and the architectures listed above where `sum.c` leaves the portable loop.
- **The survey** cannot see `kernel/simd/intrin.h`'s vectors, so it lists `sum.c` (and `kernel/generic/dot.c`) as plain C where upstream vectorizes them.
