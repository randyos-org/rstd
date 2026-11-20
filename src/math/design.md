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

## Kernel API

Every routine of `rstd.math.kernels`, by domain ([#607](https://github.com/randyos-org/rstd/issues/607)): what rstd implements today, what it maps to in Accelerate (the macOS backend) and MKL (a later backend, [#617](https://github.com/randyos-org/rstd/issues/617)), and whether v0.1.0 has it. The wrappers here are the only entry points: each runs rstd's kernels or, where the process routes it, the vendor's.

**Signatures.** `name(allocator, io, comptime T, scalars, inputs, outputs)`, as the BLAS already is (`axpy(io, T, alpha, x, y)`): the allocator only where the routine allocates, the `Io` where it may thread, vectors as `ConstStrided(T)` in and `Strided(T)` out, matrices as column-major views with a leading dimension. Knobs are trailing plain arguments (`Summation`, a norm, a method). An output may alias an input only where the routine says so.

**Status.** "yes": in rstd today, kept for v0.1.0 as it is or renamed in its move. "stub": its API is in place and returns `error.NotImplemented`. An issue number: work planned for v0.1.0. "backlog": named here so the API has room for it; each domain's backlog rows are one `backlog` issue, split when picked up: `lapack` [#623](https://github.com/randyos-org/rstd/issues/623), `sparse` [#624](https://github.com/randyos-org/rstd/issues/624), `fft` [#625](https://github.com/randyos-org/rstd/issues/625), `elementwise` [#626](https://github.com/randyos-org/rstd/issues/626), `dsp` [#627](https://github.com/randyos-org/rstd/issues/627), `stats` [#628](https://github.com/randyos-org/rstd/issues/628).

### `kernels.blas`

| Routines | What | rstd today | Accelerate | MKL | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `copy`, `swap`, `scal`, `scalReal`, `axpy`, `axpyc`, `axpby` | vector updates | `blas1`, per-core kernels | `cblas_?copy` and the rest, `catlas_?axpby`; `axpyc` not | `cblas_?copy` and the rest, `cblas_?axpby` | yes |
| `dot`, `dotc`, `dsdot`, `sdsdot`, `dotWith`, `dotcWith` | dot products; `*With` takes a `Summation` | `blas1` | `cblas_?dot`, `?dotc_sub`; `*With` not | `cblas_?dot` | yes |
| `nrm2`, `asum`, `sum`, `asumWith`, `sumWith` | norms and sums | `blas1` | `cblas_?nrm2`, `?asum`; `sum` none | `cblas_?nrm2`, `?asum`; `sum` none | yes |
| `index`, `extreme` | position and value of the largest or smallest (magnitude or value) | `blas1` | `cblas_i?amax`; not routed (rstd's NaN rule) | `cblas_i?amax`, `i?amin` | yes |
| `rot`, `rotm`, `rotg`, `zrotg`, `rotmg`, `cabs1` | plane rotations | `blas1` | `cblas_?rot`, `?rotm`; the generators not routed | `cblas_?rot` and the rest | yes |
| `sbdot`, `bscal`, `sbtobf16`, `bf16to` | bfloat16 level 1 | `blas1` | none | none | yes |
| `gemv`, `gbmv`, `ger`, `gerc`, `gemvHalf` | general level 2; `Half` in `f16` or bfloat16 | `blas2` | `cblas_?gemv` and the rest; `Half` none | `cblas_?gemv` and the rest | yes |
| `symv`, `hemv`, `sbmv`, `hbmv`, `spmv`, `hpmv`, `syr`, `her`, `spr`, `hpr`, `syr2`, `her2`, `spr2`, `hpr2` | symmetric and Hermitian level 2 | `blas2` | the CBLAS forms; complex symmetric none | the CBLAS forms | yes |
| `trmv`, `trsv`, `tbmv`, `tbsv`, `tpmv`, `tpsv` | triangular level 2 | `blas2` | the CBLAS forms | the CBLAS forms | yes |
| `gemm`, `symm`, `hemm`, `syrk`, `herk`, `syr2k`, `her2k`, `trmm`, `trsm` | level 3 | `blas3` | the CBLAS forms | the CBLAS forms | yes |
| `gemmt`, `gemm3m`, `symm3m`, `hemm3m` | one triangle of a product; complex in three real products | `blas3` | none | `cblas_?gemmt`, `cblas_?gemm3m` | yes |
| `gemmHalf`, `gemmtHalf` | `f16` and bfloat16 products | `blas3` | none | `cblas_hgemm`, `cblas_gemm_bf16bf16f32` | yes |
| `gemmBatch`, `gemmBatchStrided`, `gemmHalfBatch`, `gemmHalfBatchStrided` | batched products | `blas3` | none | `cblas_?gemm_batch`, `_batch_strided` | yes |
| `geadd`, `omatcopy`, `imatcopy`, `laswp` | matrix add, copies and transposes, row swaps | `ext` | none; `laswp` in LAPACK | `mkl_?omatadd`, `mkl_?omatcopy`, `mkl_?imatcopy`; `?laswp` | yes |

### `kernels.lapack`

Accelerate (macOS 13.3's LAPACK, `$NEWLAPACK$ILP64`) and MKL have every routine here; rstd routes only `potrf` to Accelerate so far.

| Routines | What | v0.1.0 | Backlog |
| --- | --- | --- | --- |
| `potrf`, `potf2`, `potrs`, `potri`, `pbtrf` | Cholesky: factor, solve, inverse, banded factor | yes | `posv`, `pbtrs`, `pptrf` (packed), `pstrf` (pivoted), `pocon` |
| `getrf`, `getrs`, `gesv`, `getri`, `gbtrf` | LU with partial pivoting | yes | `gbtrs`, `gbsv`, `gesvx` (expert) |
| `trtri`, `trti2`, `lauum`, `lauu2`, `trtrs` | triangular inverse, products, solves | yes | `trcon` |
| `sytrf`, `hetrf`, `sytrfRook`, `hetrfRook`, `sygst`, `hegst` | symmetric indefinite factor; generalized eigenproblem's reduction | yes | `sytrs`, `sysv`, `sytri` |
| `geqrf`, `orgqr`, `ungqr`, `ormqr`, `unmqr`, `gelqf`, `orglq`, `unglq`, `ormlq`, `unmlq`, `geqp3`, `gels`, `gelsy` | QR and LQ, pivoted QR, least squares | yes | `gelsd`, `gelss` (SVD least squares), `gglse`, `ggglm` |
| `syevd`, `heevd`, `syev`, `heev` | symmetric and Hermitian eigenproblems (with `sytrd`, `stebr` inside) | yes | `syevr`, `syevx` (subsets), `sygvd` (generalized) |
| `geev`, `hseqr`, `gees`, `trevc3`, `trevc`, `gebal`, `gebak`, `gehrd`, `orghr`, `unghr`, `trexc`, `trsen`, `trsyl`, `tgsyl` | nonsymmetric eigenproblems, Schur forms, Sylvester equations | yes | `geevx`, `ggev`, `gges` (generalized) |
| `gesvd`, `gesdd` | singular value decomposition | yes | `gesvdx`, `gejsv`, `ggsvd3` |
| `lange`, `gecon`, `laset`, `lapy2` | norms, condition estimate, fills | yes | `lansy`, `lantr` |

### `kernels.sparse`

| Routines | What | rstd today | Accelerate | MKL | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `Csr`, `Csc`, `Triplets`; `toCsc`, `toCsr`, `fromDense`, `toDense`, `transpose`, `ctranspose`, `selectCols`, `scal`, `sub` | formats, conversions, scaling, difference | `sparse/` | `SparseConvertToCoordinate` and the like | `mkl_sparse_convert_*`, `mkl_sparse_?_add` | yes |
| `mulVec`, `mulDense` | sparse times dense vector or matrix | `sparse/mul.zig` | `SparseMultiply`, `SparseMultiplyAdd` | `mkl_sparse_?_mv`, `mkl_sparse_?_mm` | yes |
| `spgemm` (`symbolic`, `numeric`) | sparse times sparse | `sparse/spgemm.zig` | none | `mkl_sparse_sp2m` | yes |
| `solveTri` | sparse triangular solve | `sparse/triangular.zig` | `SparseSolve` (triangular) | `mkl_sparse_?_trsv` | yes |
| `cg`, `gmres` (over an `Operator`) | iterative solvers | `sparse/iterative.zig` | `SparseConjugateGradient`, `SparseGMRES` | RCI ISS (`dcg`, `dfgmres`) | yes |
| `add`; Cholesky, LDLᵀ, LU, QR factorizations and solves; `lsmr`; preconditioners | direct solvers and the rest | none | `SparseFactor`, `SparseSolve`, `SparseLSMR`, `SparsePreconditioner` | PARDISO, `mkl_sparse_qr` | backlog |

### `kernels.fft`

`Layout` describes the data (count, stride, distance), not a knob; the norm is a trailing argument in the kernel layer, where the high-level `Options` struct moves to the `*With` siblings ([#611](https://github.com/randyos-org/rstd/issues/611)).

| Routines | What | rstd today | Accelerate | MKL | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `Plan(T)`: `transform`, `transformInto` | complex transforms of any length, batched by `Layout` | PocketFFT's `cfftp` and Bluestein; pffft for `f32` at its lengths | `vDSP_DFT_zop`, `vDSP_fft_zip` (powers of two) | DFTI, complex | yes |
| `RealPlan(T)`: `forward`, `backward`, `forwardUnordered`, `backwardUnordered`, `convolveAccumulate` | real transforms; the unordered forms and their products for fast convolution | `rfftp`, pffft | `vDSP_DFT_zrop`, `vDSP_fft_zrip` | DFTI, real | yes |
| `StaticPlan`, `StaticRealPlan` | lengths fixed at compile time | stub | none | none | stub |
| `TrigPlan(T)` (`dct1` to `dct4`, `dst1` to `dst4`) | trigonometric transforms | stub | `vDSP_DCT_CreateSetup` (II, III, IV) | trigonometric transforms (TT) | stub |
| `Mdct(T)` | modified DCT, for audio frames | stub | none | none | stub |
| 2-D and 3-D transforms | multi-dimensional | none | `vDSP_fft2d_zip` | DFTI | backlog |

### `kernels.elementwise`

`f(io, T, x, y)`: `y[i] = f(x[i])`, and `f(io, T, x, y, z)` for two inputs; `y` may be `x`. v0.1.0 has the set Zig's vector builtins compute ([#614](https://github.com/randyos-org/rstd/issues/614)), each stating its accuracy; the rest waits for the backlog, where faster vector polynomials may replace builtins.

| Routines | Accelerate (vForce, vDSP) | MKL (VM) | v0.1.0 |
| --- | --- | --- | --- |
| `add`, `sub`, `mul`, `div` | `vDSP_vadd`, `vDSP_vsub`, `vDSP_vmul`, `vvdiv` | `v?Add`, `v?Sub`, `v?Mul`, `v?Div` | #614 |
| `sqrt` | `vvsqrt` | `v?Sqrt` | #614 |
| `exp`, `exp2`, `log`, `log2`, `log10` | `vvexp`, `vvexp2`, `vvlog`, `vvlog2`, `vvlog10` | `v?Exp`, `v?Exp2`, `v?Ln`, `v?Log2`, `v?Log10` | #614 |
| `sin`, `cos`, `tan` | `vvsin`, `vvcos`, `vvtan` | `v?Sin`, `v?Cos`, `v?Tan` | #614 |
| `abs`, `floor`, `ceil`, `trunc`, `round` | `vvfabs`, `vvfloor`, `vvceil`, `vvint`, `vvnint` | `v?Abs`, `v?Floor`, `v?Ceil`, `v?Trunc`, `v?Round` | #614 |
| `rsqrt`, `recip`, `cbrt`, `hypot` | `vvrsqrt`, `vvrec`; `cbrt`, `hypot` none | `v?InvSqrt`, `v?Inv`, `v?Cbrt`, `v?Hypot` | backlog |
| `expm1`, `log1p`, `pow`, `powx` | `vvexpm1`, `vvlog1p`, `vvpow`, `vvpows` | `v?Expm1`, `v?Log1p`, `v?Pow`, `v?Powx` | backlog |
| `sincos`, `asin`, `acos`, `atan`, `atan2` | `vvsincos`, `vvasin`, `vvacos`, `vvatan`, `vvatan2` | `v?SinCos`, `v?Asin`, `v?Acos`, `v?Atan`, `v?Atan2` | backlog |
| `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh` | `vvsinh` and the rest | `v?Sinh` and the rest | backlog |
| `erf`, `erfc`, `lgamma`, `tgamma` | none | `v?Erf`, `v?Erfc`, `v?LGamma`, `v?TGamma` | backlog |
| complex `abs`, `arg`, `conj`, `exp` | `vDSP_zvabs`, `vDSP_zvphas`, `vDSP_zvconj` | `v?Abs`, `v?Arg`, `v?Conj`, `v?Exp` | backlog |

### `kernels.dsp`

| Routines | What | rstd today | Accelerate (vDSP) | MKL | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `window(T, out, kind, symmetry)` | Hann, Hamming, Blackman, Kaiser, sine windows | `signal/window.zig` | `vDSP_hann_window`, `vDSP_hamm_window`, `vDSP_blkman_window` | none | yes, moved by #612 |
| partitioned convolution (`RealPlan.convolveAccumulate`) | long FIR filters through FFT products | `fft/real_plan.zig`, used by `signal.Convolver` | none | VS convolution (`vslsConvExec`) | yes |
| `sadd`, `sdiv`, `ma`, `fill`, `ramp`, `neg`, `clip`, `reverse` | scalar operations on vectors, fills, clipping | none (`scal` is in `blas`) | `vDSP_vsadd`, `vDSP_vsdiv`, `vDSP_vma`, `vDSP_vfill`, `vDSP_vramp`, `vDSP_vneg`, `vDSP_vclip`, `vDSP_vrvrs` | VM in part (`v?LinearFrac`) | backlog |
| `convolve`, `correlate` (direct) | short FIR filters, correlation | none | `vDSP_conv` | VS (`vslsConvExec`, `vslsCorrExec`) | backlog |
| `biquad`, `decimate`, `runningSum`, `poly`, `interpolate` | IIR sections, resampling, running sums, polynomials, linear interpolation | none | `vDSP_biquad`, `vDSP_desamp`, `vDSP_vrsum`, `vDSP_vpoly`, `vDSP_vlint` | none | backlog |
| `magnitudes`, `phases`, `toDecibels`, split and interleaved complex | spectra | none | `vDSP_zvmags`, `vDSP_zvphas`, `vDSP_vdbcon`, `vDSP_ctoz`, `vDSP_ztoc` | VM in part | backlog |

### `kernels.stats`

| Routines | What | rstd today | Accelerate | MKL (VS) | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `sum`, `sumWith`, `index`, `extreme` | sums and extremes | `kernels.blas` (reused, not copied) | `vDSP_sve`, `vDSP_maxv`, `vDSP_minv` | VS | yes |
| `deviations(io, T, x, m, s)` | `Σ(x − m)` and `Σ(x − m)²` in one pass, every `Summation` | none | `vDSP_measqv`, `vDSP_mvessq` in part | VS central moments | #605 |
| running moments (`Running`'s update) | Welford's one-pass mean and variance | `stats/running.zig` | none | VS streaming | yes, moved by #613 |
| histogram counting | counts into bins | `stats/histogram.zig` | none | none | yes, moved by #613 |
| selection for quantiles | order statistics | `stats/quantile.zig` | none | VS quantiles, order statistics | yes, moved by #613 |
| covariance, correlation | through `kernels.blas.syrk` | `stats/cov.zig` | none | VS | yes |
| third and fourth moments (skewness, kurtosis), median and mean absolute deviation, robust estimates, random numbers | the rest of a summary-statistics library | none | none | VS, VS RNG | backlog |

### `rstd.math` (scalar)

`floor`, `ceil`, `round` (exact where the builtins are not), `x87_f64`, `Fraction` and `rationalApprox`: today's `math.zig`, moved by [#608](https://github.com/randyos-org/rstd/issues/608).

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
