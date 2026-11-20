# `rstd.math` design

How rstd's math is organized: one low-level kernel layer, modeled on Intel's MKL and Apple's Accelerate, under a high-level API with no tuning options. Decided with the user on 2026-11-20 ([#597](https://github.com/randyos-org/rstd/issues/597), [#606](https://github.com/randyos-org/rstd/issues/606)); the moves that get there are [#607](https://github.com/randyos-org/rstd/issues/607) to [#616](https://github.com/randyos-org/rstd/issues/616) and [#618](https://github.com/randyos-org/rstd/issues/618). The details of each domain stay in its own documents: BLAS and LAPACK in [`linalg/design.md`](../linalg/design.md), threading in [`threading.md`](../linalg/blas/threading.md).

## Goals

- **The API is settled by v0.1.0.** Later work is tuning, benchmarking and smaller changes, not a new shape.
- **The best algorithm for each routine, collected natively.** rstd is an anthology of the best open-source methods (OpenBLAS, BLIS, ReLAPACK, faer, PocketFFT and so on), and its algorithms are the seed of a later port to the user's language Cesium. A clean, correct, well-sourced implementation matters more than the last few percent of speed.
- **Users never tune.** The plain functions give the best answer: accurate, deterministic, fast. Knobs exist for those who want them, one step away, never in the main signatures.
- **One library, every domain.** Linear algebra, FFT, signal processing and statistics share one kernel layer and its infrastructure, as MKL's and Accelerate's domains do, rather than each growing its own.

## Layout

```text
rstd.math                  scalar functions (today's math.zig)
rstd.math.kernels          the kernel layer: every optimized routine, every knob explicit
  blas                     BLAS levels 1 to 3 and rstd's extensions (sum, axpby, gemmt, batch, bf16)
  lapack                   LAPACK
  sparse                   sparse kernels (products, triangular solves, formats as kernels take them)
  fft                      transform kernels (pffft, PocketFFT's cfftp, rfftp, Bluestein, twiddles)
  elementwise              transcendental functions over arrays (MKL's VM, Accelerate's vForce)
  dsp                      vector and signal primitives (vDSP's): scalar ops on vectors, convolution, windows
  stats                    summary-statistics reductions (MKL's VS): moments, deviations, histograms
rstd.math.linalg           vectors, matrices, dense and sparse algorithms for users
rstd.math.fft              plans and transforms for users
rstd.math.signal           windows, STFT, convolvers for users
rstd.math.stats            mean, variance, quantiles, covariance for users
```

Each high-level module calls the kernel layer and nothing below it. Random numbers (MKL's VS RNG) would go in `kernels.stats` and `rstd.math.stats` when wanted; nothing plans them yet.

Where MKL's and Accelerate's domains fall:

| Domain | MKL | Accelerate | rstd |
| --- | --- | --- | --- |
| BLAS | BLAS | BLAS | `kernels.blas` |
| LAPACK | LAPACK | LAPACK | `kernels.lapack` |
| Sparse | Sparse BLAS, PARDISO | Sparse Solvers | `kernels.sparse` |
| FFT | DFTI | vDSP's FFTs | `kernels.fft` |
| Elementwise math | VM | vForce | `kernels.elementwise` |
| Vector and signal primitives | VM, VS (in part) | vDSP | `kernels.dsp` |
| Summary statistics | VS | vDSP (in part) | `kernels.stats` |

FFT and statistics are not built on BLAS, in MKL, Accelerate, FFTW or PocketFFT alike: a butterfly or a deviation sum is not a BLAS operation. The domains share infrastructure (below), not each other's kernels. [#607](https://github.com/randyos-org/rstd/issues/607) surveys MKL's and Accelerate's routines to name what each `kernels` domain offers for v0.1.0 and what waits in the backlog.

## The kernel layer

- **Every knob is an explicit argument:** the element type, strides, the summation mode (`Summation`), the `Io` that threads it, and the allocator first wherever it allocates. No options structs: the kernel layer is called by code that knows what it wants.
- **Zig-native and column-major:** no CBLAS or Fortran names, no row-major layer (`linalg/design.md`).
- **Each routine states its contract:** what it computes and its source, its accuracy, and that its bits are the same for any worker count on a core, reductions in a fixed, size-determined order.
- **Exact error sets** with plain names; panics only for a broken invariant, documented.
- **Per-core dispatch where it pays:** a routine reaches its kernels through slots (`blas.kernel(.slot)`, `portable.fetch`), filled when a core is selected; a domain without per-core kernels yet uses the portable `@Vector` code.
- **Backends below the signature:** a routine may run rstd's kernels or a vendor's (below), chosen per process, never per call.

## The high-level API

- **A plain function with the best defaults and no options:** `stats.variance(io, xs, normalize)`, `fft` transforms, `Matrix` operations.
- **A `*With` sibling for advanced users, with an options struct:** `stats.varianceWith(io, xs, normalize, .{ .summation = .compensated })`. Each field defaults to what the plain function uses, so a caller names only what they change; a field takes the kernel layer's own type (`blas.Summation`) and passes it through unchanged, so a user reaches the knobs without calling the kernel layer themselves.
- **No options on the plain function** (user, 2026-11-20): nobody passes `.{}` just to say "the defaults". The options struct belongs to the `*With` sibling only. `fft`'s `Options` (`.norm`) and `stats.quantile`'s `QuantileOptions` (`.method`) move from the plain functions to new `*With` siblings in their moves ([#611](https://github.com/randyos-org/rstd/issues/611), [#613](https://github.com/randyos-org/rstd/issues/613)); today's `meanWith`, `varianceWith` and `stdDevWith` take a struct in place of their `Summation` argument ([#605](https://github.com/randyos-org/rstd/issues/605)). Build-time settings a root declares once (`rstd_linalg_options`) are not call options and stay.
- **No switch on a kernel option:** a high-level function that would have to switch on the summation mode or a method itself is a sign the algorithm belongs in a kernel ([#597](https://github.com/randyos-org/rstd/issues/597)): variance's deviation pass became a `kernels.stats` reduction ([#605](https://github.com/randyos-org/rstd/issues/605)).

## Shared infrastructure

Today inside `blas`; it moves up to `rstd.math.kernels` when a second domain uses it:

- **CPU detection and core selection** (`blas/cpu/`): one detected core per process, its kernels per domain.
- **Threading** (`blas/threading/`, `threading.zig`): worker teams sized per call, chunks by size alone, results folded in order, so bits do not depend on the worker count.
- **Vector lanes** (`blas/vector/`, `blas1/fast.zig`): the `@Vector` accumulation patterns the reductions share, first needed outside `blas` by [#605](https://github.com/randyos-org/rstd/issues/605).
- **Configuration** (`Config.zig`): the environment's settings, read once.

## Backends

- **rstd's own kernels, everywhere.** RandyOS's targets (bare-metal x86_64 and aarch64, UEFI, 32-bit Arm, the PowerPC 750) and Linux have no vendor library, and the kernels are the port's source.
- **Accelerate, the default on macOS:** part of the OS, and the only full access to Apple's matrix hardware. Routed in `blas/apple_accel.zig` today; vDSP and vForce are not bound yet.
- **MKL and Arm Performance Libraries: backlog** ([#617](https://github.com/randyos-org/rstd/issues/617)). Each is one more entry in the wrappers that already switch between rstd's kernels and Accelerate, loaded at run time (`dll.zig`), never shipped.
- **A vendor backend keeps the reproducibility rule:** verified to give the same bits for any worker count, or run single-threaded.

## Tuning bar

- Slower than the reference (OpenBLAS, or the domain's) on a target core is milestone work only where no vendor backend covers the platform. On macOS Accelerate covers Apple Silicon, so the SME2 kernels moved to the backlog ([#586](https://github.com/randyos-org/rstd/issues/586) to [#589](https://github.com/randyos-org/rstd/issues/589)).
- Work is on tiers 1 and 2 only for now ([`support.md`](../../support.md), [#603](https://github.com/randyos-org/rstd/issues/603)).
- The converted kernels and derived tables are maintained by hand; nothing is regenerated from OpenBLAS ([#618](https://github.com/randyos-org/rstd/issues/618)).

## Getting there

In the order of the `dev` to `main` pull request's Open work list:

1. This document ([#606](https://github.com/randyos-org/rstd/issues/606)), then the API survey against MKL and Accelerate ([#607](https://github.com/randyos-org/rstd/issues/607)).
2. Retire `tools/blas_asm` ([#618](https://github.com/randyos-org/rstd/issues/618)).
3. `rstd.math` with the scalar functions ([#608](https://github.com/randyos-org/rstd/issues/608)); `blas` and `lapack` into `kernels` ([#609](https://github.com/randyos-org/rstd/issues/609)); `sparse`, `fft` and `signal` split into kernels and high-level parts ([#610](https://github.com/randyos-org/rstd/issues/610) to [#612](https://github.com/randyos-org/rstd/issues/612)).
4. The deviation reduction in `kernels.stats` ([#605](https://github.com/randyos-org/rstd/issues/605)); `stats` split ([#613](https://github.com/randyos-org/rstd/issues/613)); `kernels.elementwise` ([#614](https://github.com/randyos-org/rstd/issues/614)); the high-level `linalg` ([#615](https://github.com/randyos-org/rstd/issues/615)); rgui and RandyOS follow ([#616](https://github.com/randyos-org/rstd/issues/616)).

Each move changes paths, not code: test names change with them, so the rig filters in `build/qemu/machines.zig` follow in the same commit.
