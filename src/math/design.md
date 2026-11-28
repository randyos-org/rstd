# `rstd.math` design

How rstd's math is organized: one low-level kernel layer, modeled on Intel's MKL and Apple's Accelerate, under a high-level API with no tuning options. Decided with the user on 2026-11-20 ([#597](https://github.com/randyos-org/rstd/issues/597), [#606](https://github.com/randyos-org/rstd/issues/606)); the moves that get there are [#607](https://github.com/randyos-org/rstd/issues/607) to [#616](https://github.com/randyos-org/rstd/issues/616) and [#618](https://github.com/randyos-org/rstd/issues/618). The details of each domain stay in its own documents: BLAS and LAPACK in [`linalg/design.md`](linalg/design.md), threading in [`threading.md`](kernels/blas/threading.md).

## Goals

- **The API is settled by v0.1.0.** Later work is tuning, benchmarking and smaller changes, not a new shape.
- **The best algorithm for each routine, collected natively.** Rstd is an anthology of the best open-source methods (OpenBLAS, BLIS, ReLAPACK, faer, PocketFFT and so on), and its algorithms are the seed of a later port to the user's language Cesium. A clean, correct, well-sourced implementation matters more than the last few percent of speed.
- **Users never tune.** The plain functions give the best answer: accurate, deterministic, fast. Knobs exist for those who want them, one step away, never in the main signatures.
- **One library, every domain.** Linear algebra, FFT, signal processing and statistics share one kernel layer and its infrastructure, as MKL's and Accelerate's domains do, rather than each growing its own.

## Layout

```text
rstd.math                  scalar functions (math/scalar.zig)
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

Every routine and type of `rstd.math.kernels`, with its literal signature ([#607](https://github.com/randyos-org/rstd/issues/607)). Signatures of what exists are taken from the code; the rest are this design's, and the code follows them. Each wrapper here is the only entry point: it runs rstd's kernels or, where the process routes it, the vendor's.

**Shape.** `name(allocator, io, T, knobs, scalars, inputs, outputs)`, the three shared arguments written bare below: `allocator: Allocator`, `io: Io` and `comptime T: type` (or `F` for a real float). In order: the allocator only where the routine allocates, the `Io` where it may thread, the element type `T` (`f16` to `f128`, `std.math.Complex` of them where the routine takes complex), then the routine's choices (`Trans`, `Uplo`, a `Summation`), scalars, inputs, outputs. Vectors are `ConstStrided(T)` in and `Strided(T)` out; matrices are a column-major slice with its leading dimension (`a: []const T, lda: usize`) and their dimensions as `usize`. No options structs in this layer. An output aliases an input only where its row says so.

**Errors.** Each routine returns an exact error set named for it (`PotrfError`) or its level (`blas2.Error`): bad input (`SliceTooShort`, `CountMismatch`, `LeadingDimensionTooSmall`) is an error, `error.OutOfMemory` where it allocates, `error.Canceled` where it threads; it panics only for a broken invariant.

**Columns.** "Accelerate" and "MKL" name the vendor routine a wrapper maps to; "(routed)" marks the calls rstd sends to Accelerate today, the rest are mapped when that backend comes ([#617](https://github.com/randyos-org/rstd/issues/617) for MKL). "Source" names what rstd's implementation is ported from, or will be. "v0.1.0": "yes" exists; "stub" has its API and returns `error.NotImplemented`; an issue number is planned for v0.1.0; "backlog" is named here, so the API has room for it, each domain's backlog in one issue: `lapack` [#623](https://github.com/randyos-org/rstd/issues/623), `sparse` [#624](https://github.com/randyos-org/rstd/issues/624), `fft` [#625](https://github.com/randyos-org/rstd/issues/625), `elementwise` [#626](https://github.com/randyos-org/rstd/issues/626), `dsp` [#627](https://github.com/randyos-org/rstd/issues/627), `stats` [#628](https://github.com/randyos-org/rstd/issues/628).

**References for what is new:** NumPy and SciPy (names, modes, defaults), MATLAB (argument order), Julia's `LinearAlgebra`, `Statistics`, `StatsBase`, `DSP.jl` and `SpecialFunctions.jl` (algorithms, MIT), Fortran's intrinsics, Accelerate's vDSP and vForce, MKL's VM and VS, faer (sparse solvers, MIT), and muhammad-fiaz's num.zig ([github.com/muhammad-fiaz/num.zig](https://github.com/muhammad-fiaz/num.zig), MIT: a starting point for some backlog routines, its algorithms unoptimized).

### Types

Each type the kernel layer defines or takes, after the moves; the high-level option structs are at the end. BLAS and LAPACK keep their routines' historical names and the four BLAS argument types every BLAS user knows (`Trans`, `Uplo`, `Diag`, `Side`); every type rstd names itself has a descriptive name, and "Today" gives the name in the code now, renamed in its domain's move (user, 2026-11-20).

| Type | Definition | Today | Module | v0.1.0 |
| --- | --- | --- | --- | --- |
| `Strided(T)`, `ConstStrided(T)` | `struct { data: []T, n: usize, stride: isize }` (`[]const T` for `ConstStrided`): `n` elements, `stride` apart, negative backward | same | `kernels.blas` | yes |
| `BFloat16` | `enum(u16) { _ }`: a bfloat16's bits | same | `kernels.blas` | yes |
| `Trans` | `enum { n, t, r, c, o, u, s, d }`: none, transpose, conjugate, conjugate transpose; the rest OpenBLAS's conjugated-`x` letters | same | `kernels.blas` | yes |
| `Uplo` | `enum { u, l, v, m }`: upper, lower; `v`, `m` OpenBLAS's reversed forms | same | `kernels.blas` | yes |
| `Diag` | `enum { u, n }`: unit, non-unit | same | `kernels.blas` | yes |
| `Side` | `enum { l, r }`: left, right | same | `kernels.blas` | yes |
| `Extreme` | `enum { largest_magnitude, smallest_magnitude, largest, smallest }`: by magnitude (`\|re\| + \|im\|`) or by value | same | `kernels.blas` | yes |
| `Summation` | `enum { plain, compensated, fast }`: `f32` in `f64`; Neumaier's; `f32` in `f32` lanes | same | `kernels.blas` | yes |
| `HalfPrecision` | `enum { bf16_to_f32, bf16, f16_to_f32 }`: bfloat16 in, `f32` out; bfloat16 in and out; `f16` in, `f32` out | same | `kernels.blas` | yes |
| `GemmGroup(In, Out)` | `struct { transa: Trans = .n, transb: Trans = .n, m: usize, n: usize, k: usize, alpha: Out, a: []const []const In, lda: usize, b: []const []const In, ldb: usize, beta: Out, c: []const []Out, ldc: usize }`: products of one shape for `gemmBatch` | same | `kernels.blas` | yes |
| `Pivot` | `packed struct(usize) { row: Row, block: enum(u1) { one, two } }`, `Row` the unsigned integer one bit narrower than `usize` (`u63` on 64-bit): a Bunch-Kaufman pivot | same | `kernels.lapack` | yes |
| `GeneralizedProblem` | `enum { ax_equals_lambda_bx, abx_equals_lambda_x, bax_equals_lambda_x }`: `Ax = λBx`, `ABx = λx`, `BAx = λx` | same | `kernels.lapack` | yes |
| `VectorJob` | `enum { none, compute }`: eigenvalues or singular values only, or their vectors too | same | `kernels.lapack` | yes |
| `SingularVectorJob` | `enum { all, leading, overwrite, none }`: all of `U` or `Vᴴ`, the leading `min(m, n)`, those written over `A`, none | same | `kernels.lapack` | yes |
| `NormType` | `enum { max_abs, one, infinity, frobenius }` | same | `kernels.lapack` | yes |
| `Eigenvalues(T)` | `[]T` for complex `T`; `struct { re: []T, im: []T }` for real | same | `kernels.lapack` | yes |
| `EigenvalueSelector(T)` | `*const fn (T) bool` for complex; `*const fn (Float(T), Float(T)) bool` for real: which eigenvalues `gees` sorts first | same | `kernels.lapack` | yes |
| `PairSelector(T)` | `*const fn (T, T) bool` for complex; `*const fn (Float(T), Float(T), Float(T)) bool` for real (`α_re`, `α_im`, `β`): `gges`'s | none | `kernels.lapack` | backlog (#623) |
| `SchurJob` | `enum { eigenvalues, schur_form }` | same | `kernels.lapack` | yes |
| `SchurVectors` | `enum { none, from_identity, accumulate }`: Schur vectors not computed, computed from the identity, or applied to the given `Z` | same | `kernels.lapack` | yes |
| `ReorderVectors` | `enum { none, update }`: whether reordering a Schur form updates `Q` | same | `kernels.lapack` | yes |
| `EigenvectorSide` | `enum { right, left, both }` | same | `kernels.lapack` | yes |
| `EigenvectorSelection` | `enum { all, back_transformed, selected }` | same | `kernels.lapack` | yes |
| `Balancing` | `enum { none, permute, scale, both }` | same | `kernels.lapack` | yes |
| `BalancedRange` | `struct { lo: usize, hi: usize }`: the rows and columns balancing left to reduce | same | `kernels.lapack` | yes |
| `ConditionJob` | `enum { none, eigenvalues, subspace, both }`: which condition numbers reordering computes | same | `kernels.lapack` | yes |
| `ReorderedSchur(T)` | `struct { selected: usize, eigenvalue_condition: ?Float(T) = null, subspace_separation: ?Float(T) = null }` | same | `kernels.lapack` | yes |
| `SylvesterSolution(F)` | `struct { scale: F, perturbed: bool }`: `X` was scaled by `scale`; `perturbed` when eigenvalues were nudged apart | same | `kernels.lapack` | yes |
| `GeneralizedSylvesterSolution(F)` | `struct { scale: F, separation: ?F = null, perturbed: bool }` | same | `kernels.lapack` | yes |
| `SpectrumRange(F)` | `union(enum) { all, values: struct { lo: F, hi: F }, indices: struct { first: usize, last: usize } }`: which eigenvalues or singular values | same | `kernels.lapack` | yes |
| `SpectrumOrder` | `enum { ascending, by_block }`: `stebz`'s order, all ascending or by split block | same | `kernels.lapack` | yes |
| `BisectedEigenvalues` | `struct { count: usize, blocks: usize, unconverged: bool = false, incomplete: bool = false }`: what `stebz` found, the reference's `M`, `NSPLIT` and `INFO` 1 to 3 | same | `kernels.lapack` | yes |
| `FactorizationInput` | `enum { factor, factored, equilibrate_and_factor }`: `gesvx` factors `A`, takes its given factor, or equilibrates first | same | `kernels.lapack` | yes |
| `Equilibration` | `enum { none, rows, columns, both }`: what `laqge` and `gesvx` scaled | none | `kernels.lapack` | yes |
| `LuSolveReport(F)` | `struct { rcond: F, pivot_growth: F, singular: ?usize = null, ill_conditioned: bool = false }`: what `gesvx` reports besides `X`, the reference's `RCOND`, `RPVGRW` and `INFO > 0` | same | `kernels.lapack` | yes |
| `CompressedColumns(T)`, `CompressedRows(T)` | `struct { rows: u32, cols: u32, col_start: []u32, row: []u32, value: []T }` (rows: `row_start`, `col`): compressed sparse columns and rows | same | `kernels.sparse` | yes |
| `Triplets(T)` | `struct { rows: u32, cols: u32, i: ArrayList(u32), j: ArrayList(u32), x: ArrayList(T) }`: assembly | same | `kernels.sparse` | yes |
| `SparseVector(T)` | `struct { len: usize, index: []const u32, value: []const T }`: a sparse column or row | same | `kernels.sparse` | yes |
| `LinearOperator(T)` | `struct { ctx: *const anyopaque, apply: *const fn (ctx: *const anyopaque, io: Io, x: []const T, y: []T) Io.Cancelable!void }`: `y ← Ax` (SciPy's `LinearOperator`) | same | `kernels.sparse` | yes |
| `SolveReport` | `struct { iterations: usize, residual: f64 }`: the result of an iterative solve | same | `kernels.sparse` | yes |
| `SparseCholesky(T)`, `SparseLdlt(T)`, `SparseLu(T)`, `SparseQr(T)` | `struct { perm: []u32, parent: []u32, l: CompressedColumns(T) }` (fill-reducing order, elimination tree, factor); `SparseLdlt` adds `d: []T`; `SparseLu`: `struct { row_perm: []u32, col_perm: []u32, l: CompressedColumns(T), u: CompressedColumns(T) }`; `SparseQr`: `struct { col_perm: []u32, v: CompressedColumns(T), tau: []T, r: CompressedColumns(T) }`. Analyzed once from the pattern, refactored per values | none | `kernels.sparse` | backlog (#624) |
| `JacobiPreconditioner(T)`, `IncompleteLu(T)` | `struct { inv_diag: []T }`; `struct { lu: CompressedRows(T) }`: each with `operator(self) LinearOperator(T)` | none | `kernels.sparse` | backlog (#624) |
| `Direction` | `enum { forward, backward }` | same | `kernels.fft` | yes |
| `Normalization` | `enum { backward, ortho, forward }`: which direction scales by `1/n`; `ortho` scales both by `1/√n` | same | `kernels.fft` | yes |
| `BatchLayout`, `BatchLayouts` | `struct { count: usize = 1, stride: usize = 1, distance: ?usize = null }`; `struct { in: BatchLayout = .{}, out: BatchLayout = .{} }`: a batch's shape in memory | same | `kernels.fft` | yes |
| `Plan(T)`, `RealPlan(T)` | `struct { n: usize, scratch_len: usize, kind: Kind }`: a complex or real transform of length `n`, immutable, shared by threads | same | `kernels.fft` | yes |
| `StaticPlan(n, T)`, `StaticRealPlan(n, T)` | `struct {}`: `Plan` with `n` fixed at compile time, no allocation | same | `kernels.fft` | stub (#625) |
| `TrigTransformKind`, `TrigTransformPlan(T)` | `enum { dct1, dct2, dct3, dct4, dst1, dst2, dst3, dst4 }`; `struct { n: usize, kind: TrigTransformKind, scratch_len: usize }`: the discrete cosine and sine transforms of types I to IV | same | `kernels.fft` | stub (#625) |
| `ModifiedCosineTransform(T)` | `struct { n: usize, scratch_len: usize }`: the MDCT of audio codecs | same | `kernels.fft` | stub (#625) |
| `MultiDimPlan(T)` | `struct { shape: []const usize, plans: []Plan(T), scratch_len: usize }`: complex transforms over every axis | none | `kernels.fft` | backlog (#625) |
| `Window`, `Symmetry` | `union(enum) { hann, hamming, blackman, kaiser: f64, sine, kbd: f64 }`; `enum { symmetric, periodic }` | same | `kernels.dsp` | yes |
| `ConvolutionMode` | `enum { full, same, valid }`: NumPy's `convolve` modes | none | `kernels.dsp` | backlog (#627) |
| `BiquadSection(T)` | `struct { b0: T, b1: T, b2: T, a1: T, a2: T }`: one second-order section, `a0` 1 (SciPy's `sos` row) | none | `kernels.dsp` | backlog (#627) |
| `DecibelScale` | `enum { power, amplitude }`: `10·log10` or `20·log10` | none | `kernels.dsp` | backlog (#627) |
| `DeviationSums(T)` | `struct { sum: T, sum_sq: Float(T) }`: `Σ(x − m)` and `Σ\|x − m\|²` | same | `kernels.stats` | yes (#605) |
| `CentralMomentSums(F)` | `struct { m2: F, m3: F, m4: F }`: `Σ(x − m)²`, `Σ(x − m)³`, `Σ(x − m)⁴` | none | `kernels.stats` | backlog (#628) |
| `RunningMoments(F)` | `struct { count: u64 = 0, avg: F = 0, m2: F = 0 }`: Welford's state | same | `kernels.stats` | yes |
| `Bins(F)`, `Histogram(F)` | `union(enum) { count: usize, edges: []const F }`; `struct { counts: []u64, edges: []F }` | same | `kernels.stats` | yes |
| `RandomStream` | `struct { state: std.Random.Xoshiro256 }`: one stream; a chunk's stream is the seed's jumped by its index, so the bits do not depend on the worker count | none | `kernels.stats` | backlog (#628) |
| `Fraction(F)` | `struct { numerator: F, denominator: F }` | same | `rstd.math` | yes |
| `MomentOptions` | `struct { summation: Summation = .plain }`: `meanWith`, `varianceWith`, `stdDevWith` | same | `rstd.math.stats` | yes (#605) |
| `QuantileOptions` | `struct { method: QuantileMethod = .linear }`: `quantileWith`; `QuantileMethod` NumPy's nine | same | `rstd.math.stats` | yes (#613) |
| `FftOptions` | `struct { normalization: Normalization = .backward }`: the `*With` transforms | same | `rstd.math.fft` | yes (#611) |

### `kernels.blas`

The BLAS levels 1 to 3 and rstd's extensions, as `src/math/kernels/blas/` has them today. Every routine is generic over `T` (`f16` to `f128`, and complex); a CBLAS or MKL name with `?` stands for its letters (`s`, `d`, `c`, `z`).

#### Level 1

| Signature | Computes | Accelerate | MKL |
| --- | --- | --- | --- |
| `copy(io, T, x: ConstStrided(T), y: Strided(T)) PairError!void` | `y ← x` | `cblas_?copy` (routed) | `cblas_?copy` |
| `swap(io, T, x: Strided(T), y: Strided(T)) PairError!void` | `x ↔ y` | `cblas_?swap` (routed) | `cblas_?swap` |
| `scal(io, T, alpha: T, x: Strided(T)) Io.Cancelable!void` | `x ← αx` | `cblas_?scal` (routed) | `cblas_?scal` |
| `scalReal(io, T, alpha: Float(T), x: Strided(T)) Io.Cancelable!void` | `x ← αx`, real `α`, complex `x` | `cblas_csscal`, `cblas_zdscal` (routed) | `cblas_csscal`, `cblas_zdscal` |
| `axpy(io, T, alpha: T, x: ConstStrided(T), y: Strided(T)) PairError!void` | `y ← αx + y` | `cblas_?axpy` (routed) | `cblas_?axpy` |
| `axpyc(io, T, alpha: T, x: ConstStrided(T), y: Strided(T)) PairError!void` | `y ← α·conj(x) + y` | none | none |
| `axpby(io, T, alpha: T, x: ConstStrided(T), beta: T, y: Strided(T)) LengthError!void` | `y ← αx + βy` | `catlas_?axpby` (routed) | `cblas_?axpby` |
| `dot(io, T, x: ConstStrided(T), y: ConstStrided(T)) PairError!T` | `Σ xᵢ·yᵢ` | `cblas_?dot`, `cblas_?dotu_sub` (routed) | `cblas_?dot`, `cblas_?dotu_sub` |
| `dotc(io, T, x: ConstStrided(T), y: ConstStrided(T)) PairError!T` | `Σ conj(xᵢ)yᵢ` | `cblas_?dotc_sub` (routed) | `cblas_?dotc_sub` |
| `dotWith(io, T, x: ConstStrided(T), y: ConstStrided(T), s: Summation) PairError!T` | `dot`, summed as `s` says | none | none |
| `dotcWith(io, T, x: ConstStrided(T), y: ConstStrided(T), s: Summation) PairError!T` | `dotc`, summed as `s` says | none | none |
| `dsdot(io, x: ConstStrided(f32), y: ConstStrided(f32)) PairError!f64` | `Σ xᵢ·yᵢ` of `f32`, in `f64` | `cblas_dsdot` (routed) | `cblas_dsdot` |
| `sdsdot(io, alpha: f32, x: ConstStrided(f32), y: ConstStrided(f32)) PairError!f32` | `α + Σ xᵢ·yᵢ` in `f64`, rounded to `f32` | `cblas_sdsdot` (routed) | `cblas_sdsdot` |
| `nrm2(io, T, x: ConstStrided(T)) Float(T)` | `‖x‖₂`, without overflow or underflow | `cblas_?nrm2` (routed) | `cblas_?nrm2` |
| `asum(io, T, x: ConstStrided(T)) Io.Cancelable!Float(T)` | `Σ (\|re xᵢ\| + \|im xᵢ\|)` | `cblas_?asum` (routed) | `cblas_?asum` |
| `sum(io, T, x: ConstStrided(T)) Io.Cancelable!Float(T)` | `Σ xᵢ` (complex: `Σ (re xᵢ + im xᵢ)`, OpenBLAS's `?sum`) | none | none |
| `asumWith(io, T, x: ConstStrided(T), s: Summation) Io.Cancelable!Float(T)` | `asum`, summed as `s` says | none | none |
| `sumWith(io, T, x: ConstStrided(T), s: Summation) Io.Cancelable!Float(T)` | `sum`, summed as `s` says | none | none |
| `index(io, T, comptime which: Extreme, x: ConstStrided(T)) Io.Cancelable!?usize` | position of the `which` extreme, the first on ties; null when `x` is empty | `cblas_i?amax` (not routed: rstd's NaN rule) | `cblas_i?amax`, `cblas_i?amin` |
| `extreme(io, T, comptime which: Extreme, x: ConstStrided(T)) Io.Cancelable!?Float(T)` | value of the `which` extreme; null when `x` is empty | none | none |
| `rot(io, T, x: Strided(T), y: Strided(T), c: Float(T), s: Float(T)) PairError!void` | applies the plane rotation `(c, s)` to `(x, y)` | `cblas_?rot` (routed) | `cblas_?rot` |
| `rotm(io, T, x: Strided(T), y: Strided(T), param: *const [5]T) LengthError!void` | applies the modified Givens rotation `param` | `cblas_?rotm` (routed) | `cblas_?rotm` |
| `rotg(T, a: *T, b: *T, c: *T, s: *T) void` | the rotation zeroing `b`: `c`, `s`, and `r` in `a` | `cblas_?rotg` (not routed) | `cblas_?rotg` |
| `zrotg(T, a: *T, b: T, c: *Float(T), s: *T) void` | `rotg` for complex | `cblas_?rotg` (not routed) | `cblas_?rotg` |
| `rotmg(T, d1: *T, d2: *T, x1: *T, y1: T, param: *[5]T) void` | the modified Givens rotation zeroing `y1` | `cblas_?rotmg` (not routed) | `cblas_?rotmg` |
| `cabs1(T, z: T) Float(T)` | `\|re z\| + \|im z\|` | none | none |
| `sbdot(io, x: ConstStrided(BFloat16), y: ConstStrided(BFloat16)) LengthError!f32` | `Σ xᵢ·yᵢ` of bfloat16, in `f32` | none | none |
| `bscal(io, alpha: BFloat16, x: Strided(BFloat16)) void` | `x ← αx`, bfloat16 | none | none |
| `sbtobf16(io, F, in: ConstStrided(F), out: Strided(BFloat16)) LengthError!void` | `f32` or `f64` to bfloat16, rounded to nearest even | none | none |
| `bf16to(io, F, in: ConstStrided(BFloat16), out: Strided(F)) LengthError!void` | bfloat16 to `f32` or `f64` | none | none |

#### Level 2

| Signature | Computes | Accelerate | MKL |
| --- | --- | --- | --- |
| `gemv(allocator, io, T, trans: Trans, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `y ← α·op(A)x + βy` | `cblas_?gemv` (routed) | `cblas_?gemv` |
| `gemvHalf(io, comptime kind: HalfPrecision, trans: Trans, alpha: kind.Out(), a: []const kind.In(), lda: usize, x: ConstStrided(kind.In()), beta: kind.Out(), y: Strided(kind.Out())) HalfError!void` | `gemv` with `f16` or bfloat16 `A` and `x` (`kind`) | none | none |
| `gbmv(allocator, io, T, trans: Trans, kl: usize, ku: usize, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `gemv` for band `A` (`kl`, `ku`) | `cblas_?gbmv` (routed) | `cblas_?gbmv` |
| `ger(allocator, io, T, alpha: T, x: ConstStrided(T), y: ConstStrided(T), a: []T, lda: usize) Error!void` | `A ← α·x·yᵀ + A` | `cblas_?ger`, `cblas_?geru` (routed) | `cblas_?ger`, `cblas_?geru` |
| `gerc(allocator, io, T, alpha: T, x: ConstStrided(T), y: ConstStrided(T), a: []T, lda: usize) Error!void` | `A ← α·x·yᴴ + A` | `cblas_?gerc` (routed) | `cblas_?gerc` |
| `symv(allocator, io, T, uplo: Uplo, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `y ← αAx + βy`, symmetric `A` | `cblas_?symv`, real (routed) | `cblas_?symv`, real |
| `hemv(allocator, io, T, uplo: Uplo, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `y ← αAx + βy`, Hermitian `A` | `cblas_?hemv` (routed) | `cblas_?hemv` |
| `sbmv(allocator, io, T, uplo: Uplo, k: usize, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `symv` for symmetric band `A` | `cblas_?sbmv`, real (routed) | `cblas_?sbmv` |
| `hbmv(allocator, io, T, uplo: Uplo, k: usize, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `hemv` for Hermitian band `A` | `cblas_?hbmv` (routed) | `cblas_?hbmv` |
| `spmv(allocator, io, T, uplo: Uplo, alpha: T, ap: []const T, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `symv` for packed `A` | `cblas_?spmv`, real (routed) | `cblas_?spmv` |
| `hpmv(allocator, io, T, uplo: Uplo, alpha: T, ap: []const T, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `hemv` for packed `A` | `cblas_?hpmv` (routed) | `cblas_?hpmv` |
| `trmv(allocator, io, T, uplo: Uplo, trans: Trans, diag: Diag, a: []const T, lda: usize, x: Strided(T)) Error!void` | `x ← op(A)x`, triangular `A` | `cblas_?trmv` (routed) | `cblas_?trmv` |
| `trsv(allocator, io, T, uplo: Uplo, trans: Trans, diag: Diag, a: []const T, lda: usize, x: Strided(T)) Error!void` | `x ← op(A)⁻¹x` | `cblas_?trsv` (routed) | `cblas_?trsv` |
| `tbmv(allocator, io, T, uplo: Uplo, trans: Trans, diag: Diag, k: usize, a: []const T, lda: usize, x: Strided(T)) Error!void` | `trmv` for triangular band `A` | `cblas_?tbmv` (routed) | `cblas_?tbmv` |
| `tbsv(allocator, io, T, uplo: Uplo, trans: Trans, diag: Diag, k: usize, a: []const T, lda: usize, x: Strided(T)) Error!void` | `trsv` for triangular band `A` | `cblas_?tbsv` (routed) | `cblas_?tbsv` |
| `tpmv(allocator, io, T, uplo: Uplo, trans: Trans, diag: Diag, ap: []const T, x: Strided(T)) Error!void` | `trmv` for packed `A` | `cblas_?tpmv` (routed) | `cblas_?tpmv` |
| `tpsv(allocator, io, T, uplo: Uplo, trans: Trans, diag: Diag, ap: []const T, x: Strided(T)) Error!void` | `trsv` for packed `A` | `cblas_?tpsv` (routed) | `cblas_?tpsv` |
| `syr(allocator, io, T, uplo: Uplo, alpha: T, x: ConstStrided(T), a: []T, lda: usize) Error!void` | `A ← α·x·xᵀ + A`, symmetric | `cblas_?syr`, real (routed) | `cblas_?syr` |
| `her(allocator, io, T, uplo: Uplo, alpha: Float(T), x: ConstStrided(T), a: []T, lda: usize) Error!void` | `A ← α·x·xᴴ + A`, Hermitian, real `α` | `cblas_?her` (routed) | `cblas_?her` |
| `spr(allocator, io, T, uplo: Uplo, alpha: T, x: ConstStrided(T), ap: []T) Error!void` | `syr` for packed `A` | `cblas_?spr`, real (routed) | `cblas_?spr` |
| `hpr(allocator, io, T, uplo: Uplo, alpha: Float(T), x: ConstStrided(T), ap: []T) Error!void` | `her` for packed `A` | `cblas_?hpr` (routed) | `cblas_?hpr` |
| `syr2(allocator, io, T, uplo: Uplo, alpha: T, x: ConstStrided(T), y: ConstStrided(T), a: []T, lda: usize) Error!void` | `A ← α·x·yᵀ + α·y·xᵀ + A` | `cblas_?syr2`, real (routed) | `cblas_?syr2` |
| `her2(allocator, io, T, uplo: Uplo, alpha: T, x: ConstStrided(T), y: ConstStrided(T), a: []T, lda: usize) Error!void` | `A ← α·x·yᴴ + conj(α)yxᴴ + A` | `cblas_?her2` (routed) | `cblas_?her2` |
| `spr2(allocator, io, T, uplo: Uplo, alpha: T, x: ConstStrided(T), y: ConstStrided(T), ap: []T) Error!void` | `syr2` for packed `A` | `cblas_?spr2`, real (routed) | `cblas_?spr2` |
| `hpr2(allocator, io, T, uplo: Uplo, alpha: T, x: ConstStrided(T), y: ConstStrided(T), ap: []T) Error!void` | `her2` for packed `A` | `cblas_?hpr2` (routed) | `cblas_?hpr2` |

#### Level 3

| Signature | Computes | Accelerate | MKL |
| --- | --- | --- | --- |
| `gemm(allocator, io, T, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `C ← α·op(A)op(B) + βC` | `cblas_?gemm` (routed) | `cblas_?gemm` |
| `gemm3m(allocator, io, T, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `gemm` for complex, three real products | none | `cblas_?gemm3m` |
| `gemmt(allocator, io, T, uplo: Uplo, transa: Trans, transb: Trans, m: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `gemm`, only the `uplo` triangle of `C` | none | `cblas_?gemmt` |
| `gemmHalf(allocator, io, comptime kind: HalfPrecision, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: kind.Out(), a: []const kind.In(), lda: usize, b: []const kind.In(), ldb: usize, beta: kind.Out(), c: []kind.Out(), ldc: usize) Error!void` | `gemm` with `f16` or bfloat16 inputs (`kind`) | none | `cblas_hgemm`, `cblas_gemm_bf16bf16f32` |
| `gemmtHalf(io, comptime kind: HalfPrecision, uplo: Uplo, transa: Trans, transb: Trans, m: usize, k: usize, alpha: kind.Out(), a: []const kind.In(), lda: usize, b: []const kind.In(), ldb: usize, beta: kind.Out(), c: []kind.Out(), ldc: usize) HalfError!void` | `gemmt` with `f16` or bfloat16 inputs | none | none |
| `symm(allocator, io, T, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `C ← αAB + βC` or `αBA + βC`, symmetric `A` | `cblas_?symm` (routed) | `cblas_?symm` |
| `hemm(allocator, io, T, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `symm`, Hermitian `A` | `cblas_?hemm` (routed) | `cblas_?hemm` |
| `symm3m(allocator, io, T, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `symm` for complex, three real products | none | none |
| `hemm3m(allocator, io, T, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `hemm` for complex, three real products | none | none |
| `syrk(allocator, io, T, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: T, a: []const T, lda: usize, beta: T, c: []T, ldc: usize) Error!void` | `C ← α·op(A)op(A)ᵀ + βC`, one triangle | `cblas_?syrk` (routed) | `cblas_?syrk` |
| `herk(allocator, io, T, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: Float(T), a: []const T, lda: usize, beta: Float(T), c: []T, ldc: usize) Error!void` | `C ← α·op(A)op(A)ᴴ + βC`, real `α`, `β` | `cblas_?herk` (routed) | `cblas_?herk` |
| `syr2k(allocator, io, T, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `C ← α·op(A)op(B)ᵀ + α·op(B)op(A)ᵀ + βC` | `cblas_?syr2k` (routed) | `cblas_?syr2k` |
| `her2k(allocator, io, T, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: Float(T), c: []T, ldc: usize) Error!void` | `C ← α·op(A)op(B)ᴴ + conj(α)·op(B)op(A)ᴴ + βC` | `cblas_?her2k` (routed) | `cblas_?her2k` |
| `trmm(allocator, io, T, side: Side, uplo: Uplo, transa: Trans, diag: Diag, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []T, ldb: usize) Error!void` | `B ← α·op(A)B` or `α·B·op(A)`, triangular `A` | `cblas_?trmm` (routed) | `cblas_?trmm` |
| `trsm(allocator, io, T, side: Side, uplo: Uplo, transa: Trans, diag: Diag, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []T, ldb: usize) Error!void` | `B ← α·op(A)⁻¹B` or `α·B·op(A)⁻¹` | `cblas_?trsm` (routed) | `cblas_?trsm` |
| `gemmBatch(allocator, io, T, groups: []const GemmGroup(T, T)) BatchError!void` | `gemm` over every product of every `Group` | none | `cblas_?gemm_batch` |
| `gemmHalfBatch(allocator, io, comptime kind: HalfPrecision, groups: []const GemmGroup(kind.In(), kind.Out())) BatchError!void` | `gemmBatch` with `f16` or bfloat16 inputs | none | none |
| `gemmBatchStrided(allocator, io, T, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: T, a: []const T, lda: usize, stride_a: usize, b: []const T, ldb: usize, stride_b: usize, beta: T, c: []T, ldc: usize, stride_c: usize, count: usize) Error!void` | `gemm` over `count` products at fixed strides | none | `cblas_?gemm_batch_strided` |
| `gemmHalfBatchStrided(allocator, io, comptime kind: HalfPrecision, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: kind.Out(), a: []const kind.In(), lda: usize, stride_a: usize, b: []const kind.In(), ldb: usize, stride_b: usize, beta: kind.Out(), c: []kind.Out(), ldc: usize, stride_c: usize, count: usize) Error!void` | `gemmBatchStrided` with `f16` or bfloat16 inputs | none | none |

#### Extensions

| Signature | Computes | Accelerate | MKL |
| --- | --- | --- | --- |
| `geadd(io, T, m: usize, n: usize, alpha: T, a: []const T, lda: usize, beta: T, c: []T, ldc: usize, transa: Trans, transc: Trans) GeaddError!void` | `C ← α·op(A) + β·op(C)` | none | `mkl_?omatadd` |
| `omatcopy(io, T, trans: Trans, rows: usize, cols: usize, alpha: T, a: []const T, lda: usize, b: []T, ldb: usize) OmatcopyError!void` | `B ← α·op(A)`, out of place | none | `mkl_?omatcopy` |
| `imatcopy(allocator, io, T, trans: Trans, rows: usize, cols: usize, alpha: T, a: []T, lda: usize, ldb: usize) ImatcopyError!void` | `A ← α·op(A)`, in place | none | `mkl_?imatcopy` |
| `laswp(io, T, m: usize, n: usize, a: []T, lda: usize, first: usize, ipiv: ConstStrided(usize)) LaswpError!void` | swaps rows `first` to `first + len(ipiv)` of `A` as `ipiv` says | LAPACK `?laswp` | LAPACK `?laswp` |

### `kernels.lapack`

As `src/math/kernels/lapack/` has it today. Accelerate (macOS 13.3's `$NEWLAPACK$ILP64`) and MKL have every routine under its LAPACK name; rstd routes only `potrf` to Accelerate so far.

#### Cholesky

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `potrf(allocator, io, T, uplo: Uplo, n: usize, a: []T, lda: usize) PotrfError!void` | `A = UᴴU` or `LLᴴ`, blocked | `?potrf` (Accelerate: routed) | OpenBLAS `interface/lapack/potrf.c`; ReLAPACK with `-Drelapack` |
| `potf2(T, uplo: Uplo, n: usize, a: []T, lda: usize) Potf2Error!void` | `potrf`, unblocked | `?potf2` | OpenBLAS `potf2.c` |
| `potrs(allocator, io, T, uplo: Uplo, n: usize, nrhs: usize, a: []const T, lda: usize, b: []T, ldb: usize) PotrsError!void` | solves `AX = B` from `potrf`'s factor | `?potrs` | Reference-LAPACK `dpotrs.f` |
| `posv(allocator, io, T, uplo: Uplo, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize) PosvError!void` | `potrf`, then `potrs` | `?posv` | Reference-LAPACK `dposv.f` (v0.2.0) |
| `pocon(allocator, io, T, uplo: Uplo, n: usize, a: []const T, lda: usize, anorm: Float(T)) PoconError!Float(T)` | reciprocal condition number from `potrf`'s factor | `?pocon` | Reference-LAPACK `dpocon.f` (v0.2.0) |
| `potri(allocator, io, T, uplo: Uplo, n: usize, a: []T, lda: usize) PotriError!void` | `A⁻¹` from `potrf`'s factor | `?potri` | OpenBLAS `potri.c` |
| `pbtrf(allocator, io, T, uplo: Uplo, n: usize, kd: usize, ab: []T, ldab: usize) PbtrfError!void` | `potrf` for band `A` (`kd`) | `?pbtrf` | ReLAPACK's recursive `pbtrf` |
| `pbtrs(allocator, io, T, uplo: Uplo, n: usize, kd: usize, nrhs: usize, ab: []const T, ldab: usize, b: []T, ldb: usize) PbtrsError!void` | solves from `pbtrf`'s factor | `?pbtrs` | Reference-LAPACK `dpbtrs.f` (v0.2.0) |
| `pptrf(allocator, io, T, uplo: Uplo, n: usize, ap: []T) PptrfError!void` | `potrf` for packed `A` | `?pptrf` | Reference-LAPACK `dpptrf.f` (v0.2.0) |
| `pptrs(allocator, io, T, uplo: Uplo, n: usize, nrhs: usize, ap: []const T, b: []T, ldb: usize) PptrsError!void` | solves from `pptrf`'s factor | `?pptrs` | Reference-LAPACK `dpptrs.f` (v0.2.0) |
| `pstrf(allocator, io, T, uplo: Uplo, n: usize, a: []T, lda: usize, piv: []usize, tol: Float(T)) PstrfError!usize` | `Pᵀ·A·P = UᴴU` or `LLᴴ`, complete pivoting, blocked; returns the rank | `?pstrf` | Reference-LAPACK `dpstrf.f`, `?pstf2` as its single block |

#### LU

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `getrf(allocator, io, T, m: usize, n: usize, a: []T, lda: usize, ipiv: []usize) GetrfError!?usize` | `PA = LU`, partial pivoting; the first zero `U(i, i)`, if any | `?getrf` | OpenBLAS `getrf_single.c`, `getrf_parallel.c`; ReLAPACK with `-Drelapack` |
| `getrs(allocator, io, T, trans: Trans, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const usize, b: []T, ldb: usize) GetrsError!void` | solves `op(A)X = B` from `getrf`'s factor | `?getrs` | OpenBLAS `getrs` |
| `gesv(allocator, io, T, n: usize, nrhs: usize, a: []T, lda: usize, ipiv: []usize, b: []T, ldb: usize) GesvError!void` | `getrf`, then `getrs` | `?gesv` | OpenBLAS `gesv` |
| `gerfs(allocator, io, T, trans: Trans, n: usize, nrhs: usize, a: []const T, lda: usize, af: []const T, ldaf: usize, ipiv: []const usize, b: []const T, ldb: usize, x: []T, ldx: usize, ferr: []Float(T), berr: []Float(T)) GerfsError!void` | refines `getrs`'s `X` by `af`, `getrf`'s factor, with each column's backward error and forward bound (`lacn2`) | `?gerfs` | Reference-LAPACK `dgerfs.f`, `zgerfs.f` (v0.2.0) |
| `gesvx(allocator, io, T, fact: FactorizationInput, trans: Trans, n: usize, nrhs: usize, a: []T, lda: usize, af: []T, ldaf: usize, ipiv: []usize, equed: *Equilibration, r: []Float(T), c: []Float(T), b: []T, ldb: usize, x: []T, ldx: usize, ferr: []Float(T), berr: []Float(T)) GesvxError!LuSolveReport(Float(T))` | `gesv` with equilibration (`geequ`, `laqge`), refinement (`gerfs`) and error bounds; reports `rcond` (`gecon`), the pivot growth and a singular or ill-conditioned `A` | `?gesvx` | Reference-LAPACK `dgesvx.f`, `zgesvx.f` (v0.2.0) |
| `getri(allocator, io, T, n: usize, a: []T, lda: usize, ipiv: []const usize) GetriError!void` | `A⁻¹` from `getrf`'s factor | `?getri` | Reference-LAPACK `dgetri.f` |
| `gbtrf(T, m: usize, n: usize, kl: usize, ku: usize, ab: []T, ldab: usize, ipiv: []usize) GbtrfError!?usize` | `getrf` for band `A` (`kl`, `ku`) | `?gbtrf` | Reference-LAPACK `dgbtf2.f` |
| `gbtrs(allocator, io, T, trans: Trans, n: usize, kl: usize, ku: usize, nrhs: usize, ab: []const T, ldab: usize, ipiv: []const usize, b: []T, ldb: usize) GbtrsError!void` | solves from `gbtrf`'s factor | `?gbtrs` | Reference-LAPACK `dgbtrs.f` (v0.2.0) |
| `gbsv(allocator, io, T, n: usize, kl: usize, ku: usize, nrhs: usize, ab: []T, ldab: usize, ipiv: []usize, b: []T, ldb: usize) GbsvError!void` | `gbtrf`, then `gbtrs` | `?gbsv` | Reference-LAPACK `dgbsv.f` (v0.2.0) |
| `getrfComplete(io, T, m: usize, n: usize, a: []T, lda: usize, ipiv: []usize, jpiv: []usize) GetrfCompleteError!usize` | `PAQ = LU`, complete pivoting, unblocked, its trailing updates threaded; returns the rank | none (LAPACK's `?getc2` perturbs small pivots) | faer `lu/full_pivoting/factor.rs` (v0.2.0) |
| `getrsComplete(allocator, io, T, trans: Trans, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const usize, jpiv: []const usize, b: []T, ldb: usize) GetrsCompleteError!void` | solves `op(A)X = B` from `getrfComplete`'s factor | none | faer `lu/full_pivoting/solve.rs` (v0.2.0) |

#### Triangular

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `trtri(allocator, io, T, uplo: Uplo, diag: Diag, n: usize, a: []T, lda: usize) TrtriError!void` | `A⁻¹`, triangular, blocked | `?trtri` | OpenBLAS `trtri.c`; ReLAPACK with `-Drelapack` |
| `trti2(T, uplo: Uplo, diag: Diag, n: usize, a: []T, lda: usize) ArgError!void` | `trtri`, unblocked | `?trti2` | OpenBLAS `trti2.c` |
| `lauum(allocator, io, T, uplo: Uplo, n: usize, a: []T, lda: usize) LauumError!void` | `UUᴴ` or `LᴴL`, blocked | `?lauum` | OpenBLAS `lauum.c`; ReLAPACK with `-Drelapack` |
| `lauu2(T, uplo: Uplo, n: usize, a: []T, lda: usize) ArgError!void` | `lauum`, unblocked | `?lauu2` | OpenBLAS `lauu2.c` |
| `trtrs(allocator, io, T, uplo: Uplo, trans: Trans, diag: Diag, n: usize, nrhs: usize, a: []const T, lda: usize, b: []T, ldb: usize) TrtrsError!void` | solves `op(A)X = B`, triangular `A` | `?trtrs` | OpenBLAS `trtrs.c` |
| `trcon(allocator, io, T, norm: NormType, uplo: Uplo, diag: Diag, n: usize, a: []const T, lda: usize) TrconError!Float(T)` | reciprocal condition number of triangular `A` | `?trcon` | Reference-LAPACK `dtrcon.f` (v0.2.0) |

#### Symmetric indefinite

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `sytrf(allocator, io, T, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize` | `A = U·D·Uᵀ` or `L·D·Lᵀ`, Bunch-Kaufman; the first zero pivot, if any | `?sytrf` | Reference-LAPACK `dsytrf.f` |
| `hetrf(allocator, io, T, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize` | `sytrf` for Hermitian `A` | `?hetrf` | Reference-LAPACK `zhetrf.f` |
| `sytrfRook(allocator, io, T, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize` | `sytrf` with rook pivoting | `?sytrf_rook` | Reference-LAPACK `dsytrf_rook.f` |
| `hetrfRook(allocator, io, T, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize` | `hetrf` with rook pivoting | `?hetrf_rook` | Reference-LAPACK `zhetrf_rook.f` |
| `sytrs(T, uplo: Uplo, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const Pivot, b: []T, ldb: usize) SytrsError!void`, `hetrs`, `sytrsRook` | solves from `sytrf`'s (`sytrfRook`'s) factor | `?sytrs`, `?hetrs`, `?sytrs_rook` | Reference-LAPACK `dsytrs.f`, `dsytrs_rook.f` (v0.2.0) |
| `sysv(allocator, io, T, uplo: Uplo, n: usize, nrhs: usize, a: []T, lda: usize, ipiv: []Pivot, b: []T, ldb: usize) SysvError!void`, `hesv` | `sytrf` (`hetrf`), then `sytrs` (`hetrs`) | `?sysv`, `?hesv` | Reference-LAPACK `dsysv.f`, `zhesv.f`, over `?SYTRS` (v0.2.0) |
| `sytri(allocator, io, T, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []const Pivot) SytriError!void` | `A⁻¹` from `sytrf`'s factor | `?sytri`, `?hetri` | Reference-LAPACK `dsytri.f`, `zhetri.f` (v0.2.0) |
| `sygst(allocator, io, T, problem: GeneralizedProblem, uplo: Uplo, n: usize, a: []T, lda: usize, b: []const T, ldb: usize) SygstError!void` | reduces `Ax = λBx` (and the other `problem`s) to standard form, `B` from `potrf` | `?sygst` | ReLAPACK's recursive `sygst` |
| `hegst(allocator, io, T, problem: GeneralizedProblem, uplo: Uplo, n: usize, a: []T, lda: usize, b: []const T, ldb: usize) SygstError!void` | `sygst` for Hermitian `A` | `?hegst` | ReLAPACK's recursive `hegst` |

#### QR and LQ

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `geqrf(allocator, io, T, m: usize, n: usize, a: []T, lda: usize, tau: []T) GeqrfError!void` | `A = QR`, `Q` as reflectors and `tau` | `?geqrf` | faer v0.24.4's staircase QR |
| `orgqr(allocator, io, T, m: usize, n: usize, k: usize, a: []T, lda: usize, tau: []const T) OrgqrError!void` | forms `Q` from `geqrf`'s reflectors | `?orgqr` | faer v0.24.4 |
| `ungqr(allocator, io, T, m: usize, n: usize, k: usize, a: []T, lda: usize, tau: []const T) OrgqrError!void` | `orgqr` for complex | `?ungqr` | faer v0.24.4 |
| `ormqr(allocator, io, T, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmqrError!void` | `C ← op(Q)C` or `C·op(Q)` | `?ormqr` | faer v0.24.4 |
| `unmqr(allocator, io, T, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmqrError!void` | `ormqr` for complex | `?unmqr` | faer v0.24.4 |
| `gelqf(allocator, io, T, m: usize, n: usize, a: []T, lda: usize, tau: []T) GelqfError!void` | `A = LQ` | `?gelqf` | libflame 5.2.0 `FLA_LQ_UT` |
| `orglq(allocator, io, T, m: usize, n: usize, k: usize, a: []T, lda: usize, tau: []const T) OrglqError!void` | forms `Q` from `gelqf`'s reflectors | `?orglq` | libflame 5.2.0 |
| `unglq(allocator, io, T, m: usize, n: usize, k: usize, a: []T, lda: usize, tau: []const T) OrglqError!void` | `orglq` for complex | `?unglq` | libflame 5.2.0 |
| `ormlq(allocator, io, T, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmlqError!void` | `C ← op(Q)C` or `C·op(Q)`, LQ's `Q` | `?ormlq` | libflame 5.2.0 |
| `unmlq(allocator, io, T, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmlqError!void` | `ormlq` for complex | `?unmlq` | libflame 5.2.0 |
| `gerqf(allocator, io, T, m: usize, n: usize, a: []T, lda: usize, tau: []T) GerqfError!void` | `A = RQ` | `?gerqf` | `gelqf` on `A` reversed (v0.2.0) |
| `orgrq(allocator, io, T, m: usize, n: usize, k: usize, a: []T, lda: usize, tau: []const T) OrgrqError!void` | forms `Q` from `gerqf`'s reflectors | `?orgrq` | `orglq` reversed (v0.2.0) |
| `ungrq(allocator, io, T, m: usize, n: usize, k: usize, a: []T, lda: usize, tau: []const T) OrgrqError!void` | `orgrq` for complex | `?ungrq` | `orglq` reversed (v0.2.0) |
| `ormrq(allocator, io, T, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmrqError!void` | `C ← op(Q)C` or `C·op(Q)`, RQ's `Q` | `?ormrq` | `ormlq` reversed (v0.2.0) |
| `unmrq(allocator, io, T, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmrqError!void` | `ormrq` for complex | `?unmrq` | `ormlq` reversed (v0.2.0) |
| `ggqrf(allocator, io, T, n: usize, m: usize, p: usize, a: []T, lda: usize, taua: []T, b: []T, ldb: usize, taub: []T) GgqrfError!void` | `A = QR`, `B = QTZ` | `?ggqrf` | Reference-LAPACK `dggqrf.f` (v0.2.0) |
| `ggrqf(allocator, io, T, m: usize, p: usize, n: usize, a: []T, lda: usize, taua: []T, b: []T, ldb: usize, taub: []T) GgrqfError!void` | `A = RQ`, `B = ZTQ` | `?ggrqf` | Reference-LAPACK `dggrqf.f` (v0.2.0) |
| `geqp3(allocator, io, T, m: usize, n: usize, a: []T, lda: usize, jpvt: []usize, tau: []T) Geqp3Error!void` | `AP = QR`, column pivoting | `?geqp3` | faer v0.24.4's pivoted QR |

#### Least squares

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `gels(allocator, io, T, trans: Trans, m: usize, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize) GelsError!void` | the least-squares or minimum-norm solution, full rank, by QR or LQ | `?gels` | faer v0.24.4 `qr/no_pivoting/solve.rs` |
| `gelsy(allocator, io, T, m: usize, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize, jpvt: []usize, rcond: Float(T)) GelsyError!usize` | the least-squares solution, any rank, by pivoted QR; returns the rank | `?gelsy` | Reference-LAPACK `dgelsy.f` over faer's pivoted QR |
| `gelss(allocator, io, T, m: usize, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize, s: []Float(T), rcond: Float(T)) GelssError!usize` | `gelsd` by QR iteration | `?gelss` | Reference-LAPACK `dgelss.f` over rstd's SVD (v0.2.0) |
| `gelsd(allocator, io, T, m: usize, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize, s: []Float(T), rcond: Float(T)) GelsdError!usize` | the least-squares solution, any rank, by SVD (divide and conquer); returns the rank | `?gelsd` | Reference-LAPACK `dgelsd.f` over rstd's divide-and-conquer SVD (v0.2.0) |
| `gglse(allocator, io, T, m: usize, n: usize, p: usize, a: []T, lda: usize, b: []T, ldb: usize, c: []T, d: []T, x: []T) GglseError!void` | the least-squares solution with equality constraints `Bx = d` | `?gglse` | Reference-LAPACK `dgglse.f` over `ggrqf` (v0.2.0) |
| `ggglm(allocator, io, T, n: usize, m: usize, p: usize, a: []T, lda: usize, b: []T, ldb: usize, d: []T, x: []T, y: []T) GgglmError!void` | general Gauss-Markov linear model | `?ggglm` | Reference-LAPACK `dggglm.f` over `ggqrf` (v0.2.0) |

#### Symmetric eigenproblems

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `syevd(allocator, io, T, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, w: []Float(T)) SyevError!void` | eigenvalues, and vectors if `jobz`, of symmetric `A`, divide and conquer | `?syevd` | Reference-LAPACK `dsyevd.f`'s arguments over `sytrd`, `stebr` |
| `heevd(allocator, io, T, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, w: []Float(T)) SyevError!void` | `syevd` for Hermitian `A` | `?heevd` | Reference-LAPACK `zheevd.f` |
| `syev(allocator, io, T, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, w: []Float(T)) SyevError!void` | `syevd` by QR iteration | `?syev` | Reference-LAPACK `dsyev.f` |
| `heev(allocator, io, T, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, w: []Float(T)) SyevError!void` | `syev` for Hermitian `A` | `?heev` | Reference-LAPACK `zheev.f` |
| `sygv(allocator, io, T, problem: GeneralizedProblem, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, b: []T, ldb: usize, w: []Float(T)) SygvError!void`, `sygvd`, `hegv`, `hegvd` | the Hermitian definite pencil by `potrf`, `sygst` and `syev` (`syevd`) | `?sygv`, `?sygvd`, `?hegv`, `?hegvd` | Reference-LAPACK `dsygv.f`, `dsygvd.f` (v0.2.0) |
| `sytrd(allocator, io, T, uplo: Uplo, n: usize, a: []T, lda: usize, d: []Float(T), e: []Float(T), tau: []T) SytrdError!void`, `hetrd` | `Qᴴ·A·Q = T` tridiagonal, blocked | `?sytrd`, `?hetrd` | Reference-LAPACK `dsytrd.f` (v0.2.0) |
| `stebz(allocator, F, range: SpectrumRange(F), order: SpectrumOrder, d: []const F, e: []const F, abstol: F, w: []F, iblock: []usize, isplit: []usize) StebzError!BisectedEigenvalues` | selected eigenvalues of a symmetric tridiagonal by bisection, with their split blocks | `?stebz` | Reference-LAPACK `dstebz.f`, `dlaebz.f` (v0.2.0) |
| `stein(allocator, io, T, d: []const Float(T), e: []const Float(T), w: []const Float(T), iblock: []const usize, isplit: []const usize, z: []T, ldz: usize, ifail: []usize) SteinError!usize` | eigenvectors of a symmetric tridiagonal for `stebz`'s eigenvalues by inverse iteration; returns how many did not converge, listed in `ifail` | `?stein` | Reference-LAPACK `dstein.f`, `zstein.f` (v0.2.0) |

#### Nonsymmetric eigenproblems

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `geev(allocator, io, T, jobvl: VectorJob, jobvr: VectorJob, n: usize, a: []T, lda: usize, w: Eigenvalues(T), vl: []T, ldvl: usize, vr: []T, ldvr: usize, unconverged: ?*usize) GeevError!void` | eigenvalues, and left or right vectors, of general `A` | `?geev` | faer v0.24.4's Schur form, behind `dgeev.f`'s arguments |
| `hseqr(io, T, job: SchurJob, compz: SchurVectors, n: usize, lo: usize, hi: usize, h: []T, ldh: usize, w: Eigenvalues(T), z: []T, ldz: usize, unconverged: ?*usize) HseqrError!void` | Schur form of a Hessenberg matrix | `?hseqr` | faer v0.24.4's QR iterations |
| `gees(allocator, io, T, jobvs: VectorJob, select: ?EigenvalueSelector(T), n: usize, a: []T, lda: usize, w: Eigenvalues(T), vs: []T, ldvs: usize, unconverged: ?*usize) GeesError!usize` | Schur form `A = Z·T·Zᴴ`, sorted by `select`; returns the count selected | `?gees` | Reference-LAPACK `dgees.f` |
| `trevc3(allocator, io, T, side: EigenvectorSide, howmny: EigenvectorSelection, select: ?[]bool, n: usize, t: []T, ldt: usize, vl: []T, ldvl: usize, vr: []T, ldvr: usize, mm: usize) Trevc3Error!usize` | eigenvectors of a Schur form | `?trevc3` | Reference-LAPACK `dtrevc3.f` |
| `trevc(allocator, io, T, side: EigenvectorSide, howmny: EigenvectorSelection, select: ?[]bool, n: usize, t: []T, ldt: usize, vl: []T, ldvl: usize, vr: []T, ldvr: usize, mm: usize) Trevc3Error!usize` | `trevc3` | `?trevc` | Reference-LAPACK `dtrevc3.f` |
| `gebal(io, T, job: Balancing, n: usize, a: []T, lda: usize, scale: []Float(T)) GebalError!BalancedRange` | balances `A`; returns the block `lo..hi` | `?gebal` | Reference-LAPACK `dgebal.f` |
| `gebak(io, T, job: Balancing, side: Side, n: usize, lo: usize, hi: usize, scale: []const Float(T), m: usize, v: []T, ldv: usize) GebakError!void` | undoes `gebal` on eigenvectors | `?gebak` | Reference-LAPACK `dgebak.f` |
| `gehrd(allocator, io, T, n: usize, lo: usize, hi: usize, a: []T, lda: usize, tau: []T) GehrdError!void` | `A = Q·H·Qᴴ`, Hessenberg | `?gehrd` | faer v0.24.4 `evd/hessenberg.rs` |
| `orghr(allocator, io, T, n: usize, lo: usize, hi: usize, a: []T, lda: usize, tau: []const T) OrghrError!void` | forms `Q` from `gehrd` | `?orghr` | faer v0.24.4 |
| `unghr(allocator, io, T, n: usize, lo: usize, hi: usize, a: []T, lda: usize, tau: []const T) OrghrError!void` | `orghr` for complex | `?unghr` | faer v0.24.4 |
| `trexc(io, T, compq: ReorderVectors, n: usize, t: []T, ldt: usize, q: []T, ldq: usize, ifst: *usize, ilst: *usize) TrexcError!void` | moves one block of a Schur form | `?trexc` | Reference-LAPACK `dtrexc.f` |
| `trsen(allocator, io, T, job: ConditionJob, compq: ReorderVectors, select: []const bool, n: usize, t_in: []T, ldt: usize, q: []T, ldq: usize, w: Eigenvalues(T)) TrsenError!ReorderedSchur(T)` | reorders a Schur form by `select`, with condition numbers | `?trsen` | Reference-LAPACK `dtrsen.f` |
| `trsyl(allocator, io, T, trana: Trans, tranb: Trans, isgn: isize, m: usize, n: usize, a: []const T, lda: usize, b: []const T, ldb: usize, c: []T, ldc: usize) TrsylError!SylvesterSolution(Float(T))` | solves `op(A)X ± X·op(B) = scale·C` | `?trsyl` | ReLAPACK `dtrsyl.c` |
| `tgsyl(allocator, io, T, trans: Trans, ijob: usize, m: usize, n: usize, a: []const T, lda: usize, b: []const T, ldb: usize, c: []T, ldc: usize, d: []const T, ldd: usize, e: []const T, lde: usize, f: []T, ldf: usize) TgsylError!GeneralizedSylvesterSolution(Float(T))` | the generalized Sylvester equation, with its `dif` estimate | `?tgsyl` | ReLAPACK `dtgsyl.c` |
| `gghrd(T, compq: SchurVectors, compz: SchurVectors, n: usize, lo: usize, hi: usize, a: []T, lda: usize, b: []T, ldb: usize, q: []T, ldq: usize, z: []T, ldz: usize) GghrdError!void` | `(A, B)` to Hessenberg-triangular form | `?gghrd` | faer v0.24.4, unblocked (v0.2.0) |
| `hgeqzComplex(T, job: SchurJob, compq: SchurVectors, compz: SchurVectors, n: usize, lo: usize, hi: usize, h: []T, ldh: usize, t: []T, ldt: usize, alpha: []T, beta: []T, q: []T, ldq: usize, z: []T, ldz: usize) HgeqzError!void`, `hgeqzReal(io, T, ...)` (`alphar`, `alphai`) | the QZ iteration: the generalized Schur form | `?hgeqz` | faer v0.24.4 with `?hgeqz.f`'s ranges (v0.2.0) |
| `ggevComplex(allocator, io, T, jobvl: VectorJob, jobvr: VectorJob, n: usize, a: []T, lda: usize, b: []T, ldb: usize, alpha: []T, beta: []T, vl: []T, ldvl: usize, vr: []T, ldvr: usize) GgevError!void`, `ggevReal(...)` (`alphar`, `alphai`) | generalized eigenvalues `α/β`, and vectors | `?ggev` | faer v0.24.4's driver (v0.2.0) |
| `ggesComplex(allocator, io, T, jobvsl: VectorJob, jobvsr: VectorJob, n: usize, a: []T, lda: usize, b: []T, ldb: usize, alpha: []T, beta: []T, vsl: []T, ldvsl: usize, vsr: []T, ldvsr: usize) GgevError!void`, `ggesReal(...)` | the generalized Schur form, unordered | `?gges` | faer v0.24.4 (v0.2.0) |
| `tgexc(allocator, io, T, compq: ReorderVectors, compz: ReorderVectors, n: usize, a: []T, lda: usize, b: []T, ldb: usize, q: []T, ldq: usize, z: []T, ldz: usize, ifst: *usize, ilst: *usize) TgexcError!void` | moves one block of a generalized Schur form | `?tgexc` | Reference-LAPACK `dtgexc.f` (v0.2.0) |
| `tgsen(allocator, io, T, job: TgsenJob, compq: ReorderVectors, compz: ReorderVectors, select: []const bool, n: usize, a: []T, lda: usize, b: []T, ldb: usize, alpha: Eigenvalues(T), beta: []T, q: []T, ldq: usize, z: []T, ldz: usize) TgsenError!ReorderedPencil(Float(T))` | reorders a generalized Schur form by `select`, with condition estimates | `?tgsen` | Reference-LAPACK `dtgsen.f` (v0.2.0) |

#### SVD

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `gesvd(allocator, io, T, jobu: SingularVectorJob, jobvt: SingularVectorJob, m: usize, n: usize, a: []T, lda: usize, s: []Float(T), u: []T, ldu: usize, vt: []T, ldvt: usize, unconverged: ?*usize) GesvdError!void` | `A = U·Σ·Vᴴ`, by QR iteration's arguments | `?gesvd` | faer v0.24.4's SVD (bidiagonal divide and conquer, Reference-LAPACK's `bdsqr` as its QR iteration), behind `dgesvd.f`'s arguments, jobs and scaling |
| `gesdd(allocator, io, T, jobz: SingularVectorJob, m: usize, n: usize, a: []T, lda: usize, s: []Float(T), u: []T, ldu: usize, vt: []T, ldvt: usize) GesddError!void` | `gesvd` with the divide-and-conquer driver's arguments | `?gesdd` | faer v0.24.4's SVD, behind `dgesdd.f`'s arguments, jobs and scaling |

#### Auxiliary

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `lange(allocator, io, T, norm: NormType, m: usize, n: usize, a: []const T, lda: usize) LangeError!Float(T)` | `‖A‖` by `norm`: max, one, infinity or Frobenius | `?lange` | Reference-LAPACK `dlange.f` |
| `lansy(allocator, io, T, norm: NormType, uplo: Uplo, n: usize, a: []const T, lda: usize) LansyError!Float(T)`, `lanhe` | `lange` for symmetric (Hermitian) `A` from one triangle | `?lansy`, `?lanhe` | Reference-LAPACK `dlansy.f`, `zlanhe.f` (v0.2.0) |
| `lantr(allocator, io, T, norm: NormType, uplo: Uplo, diag: Diag, m: usize, n: usize, a: []const T, lda: usize) LantrError!Float(T)` | `lange` for trapezoidal `A` | `?lantr` | Reference-LAPACK `dlantr.f` (v0.2.0) |
| `gecon(allocator, io, T, norm: NormType, n: usize, a: []const T, lda: usize, anorm: Float(T)) GeconError!Float(T)` | reciprocal condition number from `getrf`'s factor | `?gecon` | Reference-LAPACK `dgecon.f` |
| `geequ(T, m: usize, n: usize, a: []const T, lda: usize, r: []Float(T), c: []Float(T)) GeequError!Scaling(Float(T))` | row and column factors that equilibrate `A`, with `rowcnd`, `colcnd`, `amax`; a zero row or column is an error | `?geequ` | Reference-LAPACK `dgeequ.f` |
| `laqge(T, m: usize, n: usize, a: []T, lda: usize, r: []const Float(T), c: []const Float(T), scaling: Scaling(Float(T))) LaqgeError!Equilibration` | applies `geequ`'s factors where they matter; returns what it scaled | `?laqge` | Reference-LAPACK `dlaqge.f` |
| `laset(T, uplo: ?Uplo, m: usize, n: usize, alpha: T, beta: T, a: []T, lda: usize) LasetError!void` | off-diagonal `alpha`, diagonal `beta` | `?laset` | Reference-LAPACK `dlaset.f` |
| `lapy2(F, x: F, y: F) F` | `sqrt(x² + y²)` without overflow | `?lapy2` | Reference-LAPACK `dlapy2.f` |

#### Backlog ([#623](https://github.com/randyos-org/rstd/issues/623))

Each in Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/`, as OpenBLAS ships it; modified BSD), the source of its port, and in Accelerate's and MKL's LAPACK under the name shown. Complex forms take the same signature (`he` for `sy` where LAPACK names them so).

| Signature | Computes | LAPACK name (Accelerate, MKL) |
| --- | --- | --- |
| `syevr(allocator, io, T, jobz: VectorJob, range: SpectrumRange(Float(T)), uplo: Uplo, n: usize, a: []T, lda: usize, abstol: Float(T), w: []Float(T), z: []T, ldz: usize, isuppz: []usize) SyevrError!usize` | selected eigenvalues and vectors, relatively robust representations; returns the count | `?syevr`, `?heevr` |
| `syevx(allocator, io, T, jobz: VectorJob, range: SpectrumRange(Float(T)), uplo: Uplo, n: usize, a: []T, lda: usize, abstol: Float(T), w: []Float(T), z: []T, ldz: usize, ifail: []usize) SyevxError!usize` | `syevr` by bisection and inverse iteration | `?syevx`, `?heevx` |
| `geevx(allocator, io, T, balance: Balancing, jobvl: VectorJob, jobvr: VectorJob, n: usize, a: []T, lda: usize, w: Eigenvalues(T), vl: []T, ldvl: usize, vr: []T, ldvr: usize, scale: []Float(T), rconde: []Float(T), rcondv: []Float(T)) GeevxError!BalancedRange` | `geev` with balancing and condition numbers | `?geevx` |
| `gges(allocator, io, T, jobvsl: VectorJob, jobvsr: VectorJob, select: ?PairSelector(T), n: usize, a: []T, lda: usize, b: []T, ldb: usize, alpha: Eigenvalues(T), beta: []T, vsl: []T, ldvsl: usize, vsr: []T, ldvsr: usize) GgesError!usize` | `gges` with `select`'s ordering (over `tgsen`), the shipped `ggesComplex`, `ggesReal` being unordered | `?gges3` |
| `gesvdx(allocator, io, T, jobu: VectorJob, jobvt: VectorJob, range: SpectrumRange(Float(T)), m: usize, n: usize, a: []T, lda: usize, s: []Float(T), u: []T, ldu: usize, vt: []T, ldvt: usize) GesvdxError!usize` | selected singular values and vectors; returns the count | `?gesvdx` |
| `gejsv(allocator, io, T, jobu: SingularVectorJob, jobv: SingularVectorJob, m: usize, n: usize, a: []T, lda: usize, s: []Float(T), u: []T, ldu: usize, v: []T, ldv: usize) GejsvError!void` | SVD by preconditioned Jacobi, high relative accuracy | `?gejsv` |
| `ggsvd3(allocator, io, T, jobu: VectorJob, jobv: VectorJob, jobq: VectorJob, m: usize, n: usize, p: usize, a: []T, lda: usize, b: []T, ldb: usize, alpha: []Float(T), beta: []Float(T), u: []T, ldu: usize, v: []T, ldv: usize, q: []T, ldq: usize) Ggsvd3Error!struct { k: usize, l: usize }` | generalized SVD of `(A, B)` | `?ggsvd3` |

`PairSelector(T)` is `EigenvalueSelector`'s form for a pair: `*const fn (T, T) bool` for complex, `*const fn (Float(T), Float(T), Float(T)) bool` for real (`α_re`, `α_im`, `β`).

### `kernels.sparse`

The routines are methods of the formats. Since [#610](https://github.com/randyos-org/rstd/issues/610) the iterative solvers take their options as arguments, in place of an `Options(T)` struct, and the methods rstd names itself have descriptive names: `toCompressedColumns`, `toCompressedRows`, `asColumnsOfTranspose`, `asRowsOfTranspose` and `solveTriangular` for the former `toCsc`, `toCsr`, `asCscOfTranspose`, `asCsrOfTranspose` and `solveTri`.

| Signature | Computes | Accelerate | MKL | Source | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `CompressedColumns(T)`: `alloc(allocator, rows: u32, cols: u32, nnz: usize) Allocator.Error!Self`, `deinit(self: *Self, allocator) void`, `nonzeros(self: Self) usize`, `col(self: Self, j: usize) SparseVector(T)` | storage | `SparseMatrix_Double` | `mkl_sparse_?_create_csc` | rstd | yes |
| `CompressedColumns(T)`: `fromDense(allocator, a: []const T, lda: usize, rows: u32, cols: u32, drop_below: Float(T)) FromDenseError!Self`, `toDense(self: Self, a: []T, lda: usize) ToDenseError!void` | conversion from and to dense, dropping `\|aᵢⱼ\| ≤ drop_below` | `SparseConvertFromCoordinate` | `mkl_sparse_convert_*` | rstd | yes; takes a slice and `lda` (#610) |
| `CompressedColumns(T)`: `toCompressedRows(self: Self, allocator) Allocator.Error!CompressedRows(T)`, `asRowsOfTranspose(self: Self) CompressedRows(T)`; `CompressedRows(T)`: `toCompressedColumns(self: Self, allocator) Allocator.Error!CompressedColumns(T)`, `asColumnsOfTranspose(self: Self) CompressedColumns(T)`, `row(self: Self, i: usize) SparseVector(T)` | format changes; the `as` forms reuse the arrays | none | `mkl_sparse_convert_csr` | rstd | yes |
| `CompressedColumns(T)`: `transpose(self: Self, allocator) Allocator.Error!Self`, `ctranspose(self: Self, allocator) Allocator.Error!Self` | `Aᵀ`, `Aᴴ` | `SparseGetTranspose` | `mkl_sparse_convert_csr` with `SPARSE_OPERATION_TRANSPOSE` | rstd | yes |
| `CompressedColumns(T)`: `scal(self: Self, alpha: T) void`, `sub(self: Self, allocator, i: usize, j: usize, rows: usize, cols: usize) SubError!Self`, `selectCols(self: Self, allocator, cols: []const u32) SelectError!Self` | scaling, a block, chosen columns | none | none | rstd | yes |
| `CompressedColumns(T)`, `CompressedRows(T)`: `mulVec(self: Self, allocator, io, trans: Trans, alpha: T, x: []const T, beta: T, y: []T) MulError!void` | `y ← α·op(A)x + βy` | `SparseMultiplyAdd` | `mkl_sparse_?_mv` | rstd | yes |
| `CompressedColumns(T)`: `mulDense(self: Self, allocator, io, b: []const T, ldb: usize, cols: usize, c: []T, ldc: usize) MulError!void` | `C ← AB`, dense `B` | `SparseMultiply` | `mkl_sparse_?_mm` | rstd | yes (#610) |
| `CompressedColumns(T)`: `spgemmSymbolic(self: Self, allocator, b: Self) SpgemmError!Self`, `spgemmNumeric(self: Self, allocator, b: Self, c: Self) NumericError!void`, `spgemm(self: Self, allocator, b: Self) SpgemmError!Self` | `C ← AB`, sparse `B`: the pattern, the values, both | none | `mkl_sparse_sp2m` | rstd (Gustavson) | yes |
| `CompressedColumns(T)`: `solveTriangular(self: Self, uplo: Uplo, trans: Trans, diag: Diag, x: []T) TriSolveError!void` | `x ← op(A)⁻¹x`, triangular `A` | `SparseSolve` (triangular) | `mkl_sparse_?_trsv` | rstd | yes |
| `Triplets(T)`: `init(rows: usize, cols: usize) error{TooLarge}!Self`, `deinit(self: *Self, allocator) void`, `append(self: *Self, allocator, row: usize, col: usize, value: T) (Allocator.Error \|\| error{IndexOutOfRange})!void`, `toCompressedColumns(self: Self, allocator) Error!CompressedColumns(T)` | assembly, duplicates summed | `SparseConvertFromCoordinate` | `mkl_sparse_?_create_coo` | rstd | yes |
| `SparseVector(T)`: `axpy(x: Self, alpha: T, y: []T) VectorError!void`, `dot(x: Self, y: []const T) VectorError!T`, `dotc(x: Self, y: []const T) VectorError!T` | sparse level 1 | none | `cblas_?axpyi`, `cblas_?doti`, `cblas_?dotci` | rstd | yes |
| `cg(allocator, io, T, a: LinearOperator(T), preconditioner: ?LinearOperator(T), tol: Float(T), max_iterations: usize, b: []const T, x: []T) SolveError!SolveReport` | conjugate gradients, symmetric positive definite `A` | `SparseConjugateGradient` | RCI ISS `dcg` | faer | yes; #161 (#610) |
| `gmres(allocator, io, T, a: LinearOperator(T), preconditioner: ?LinearOperator(T), tol: Float(T), max_iterations: usize, restart: usize, b: []const T, x: []T) SolveError!SolveReport` | restarted GMRES | `SparseGMRES` | RCI ISS `dfgmres` | SciPy | yes; #161 (#610) |
| `bicgstab(allocator, io, T, a: LinearOperator(T), preconditioner: ?LinearOperator(T), tol: Float(T), max_iterations: usize, b: []const T, x: []T) SolveError!SolveReport` | BiCGSTAB, right preconditioned | none | none | faer | yes; #162 |
| `minres(allocator, io, T, a: LinearOperator(T), preconditioner: ?LinearOperator(T), tol: Float(T), max_iterations: usize, b: []const T, x: []T) SolveError!SolveReport` | MINRES, Hermitian `A`, definite or not | none | none | SciPy | yes; #162 |
| `CompressedColumns(T)`: `add(self: Self, allocator, alpha: T, b: Self, beta: T) AddError!Self` | `αA + βB` | none | `mkl_sparse_?_add` | SciPy `csr_plus_csr`; rstd's `spgemm` merge | backlog |
| `SparseCholesky(T)`: `analyze(allocator, a: CompressedColumns(T), uplo: Uplo) AnalyzeError!Self`, `factor(self: *Self, allocator, io, a: CompressedColumns(T)) FactorError!void`, `solve(self: Self, io, b: []T) SolveError!void`, `deinit(self: *Self, allocator) void` | sparse `A = LLᴴ`: AMD ordering and elimination tree, then values, then solves | `SparseFactor(SparseFactorizationCholesky)` | PARDISO (`mtype` 2, 4) | faer v0.24.4 `sparse/linalg/cholesky.rs`, `amd.rs` (MIT) | backlog |
| `SparseLdlt(T)`: as `SparseCholesky(T)` | sparse `A = L·D·Lᴴ` | `SparseFactorizationLDLT` | PARDISO (`mtype` -2) | faer v0.24.4 (MIT) | backlog |
| `SparseLu(T)`: as `SparseCholesky(T)`, `analyze` without `uplo` | sparse `PAQ = LU` | none | PARDISO (`mtype` 11, 13) | faer v0.24.4 `sparse/linalg/lu.rs` (MIT) | backlog |
| `SparseQr(T)`: as `SparseLu(T)`, `solve` minimizing `‖Ax − b‖₂` | sparse `AP = QR` | `SparseFactorizationQR` | `mkl_sparse_?_qr` | faer v0.24.4 `sparse/linalg/qr.rs` (MIT) | backlog |
| `lsmr(allocator, io, T, a: LinearOperator(T), at: LinearOperator(T), damp: Float(T), atol: Float(T), btol: Float(T), max_iterations: usize, b: []const T, x: []T) SolveError!SolveReport` | the least-squares solution by LSMR | `SparseLSMR` | none | SciPy `sparse/linalg/_isolve/lsmr.py` (BSD-3) | backlog |
| `JacobiPreconditioner(T)`: `init(allocator, a: CompressedColumns(T)) JacobiError!Self`, `operator(self: *const Self) LinearOperator(T)`, `deinit(self: *Self, allocator) void` | diagonal preconditioner | `SparsePreconditionerDiagonal` | none | Saad, *Iterative Methods for Sparse Linear Systems*, 10.2 | backlog |
| `IncompleteLu(T)`: as `JacobiPreconditioner(T)` | incomplete LU without fill | none | `dcsrilu0` | Saad, algorithm 10.4 | backlog |

### `kernels.fft`

Since [#611](https://github.com/randyos-org/rstd/issues/611) the transforms take the normalization as their last argument in place of an `Options` struct; `rstd.math.fft`'s `*With` siblings take the struct.

| Signature | Computes | Accelerate | MKL | Source | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `Plan(T)`: `init(allocator, n: usize) InitError!Self`, `deinit(self: *Self, allocator) void`, `scratchLen(self: Self) usize` | plans a complex transform of length `n` | `vDSP_DFT_zop_CreateSetup` (routed, #155) | `DftiCreateDescriptor`, `DftiCommitDescriptor` | PocketFFT `cfftp`, Bluestein; pffft for `f32` | yes |
| `Plan(T)`: `transform(self: *const Self, io, data: []Complex(T), scratch: []Complex(T), direction: Direction, layout: BatchLayout, normalization: Normalization) TransformError!void` | in place, a batch by `layout` | `vDSP_DFT_Execute` (routed), `vDSP_fft_zip` | `DftiComputeForward`, `DftiComputeBackward` | as above | yes (#611) |
| `Plan(T)`: `transformInto(self: *const Self, io, out: []Complex(T), in: []const Complex(T), scratch: []Complex(T), direction: Direction, layouts: BatchLayouts, normalization: Normalization) TransformIntoError!void` | out of place | `vDSP_DFT_Execute` | `DftiComputeForward` | as above | yes (#611) |
| `RealPlan(T)`: `init(allocator, n: usize) InitError!Self`, `deinit(self: *Self, allocator) void`, `scratchLen(self: Self) usize` | plans a real transform | `vDSP_DFT_zrop_CreateSetup` (routed, #155) | `DftiCreateDescriptor` (real) | PocketFFT `rfftp`; pffft | yes |
| `RealPlan(T)`: `forward(self: *const Self, io, out: []Complex(T), in: []const T, scratch: []T, layouts: BatchLayouts, normalization: Normalization) TransformIntoError!void`, `backward(self: *const Self, io, out: []T, in: []const Complex(T), scratch: []T, layouts: BatchLayouts, normalization: Normalization) TransformIntoError!void` | real to half-spectrum and back | `vDSP_DFT_Execute` (routed), `vDSP_fft_zrip` | `DftiComputeForward` | as above | yes (#611) |
| `RealPlan(T)`: `forwardUnordered(self: *const Self, out: []T, in: []const T, scratch: []T) UnorderedError!void`, `backwardUnordered(self: *const Self, out: []T, in: []const T, scratch: []T) UnorderedError!void`, `convolveAccumulate(self: *const Self, acc: []T, a: []const T, b: []const T, scale: T) UnorderedError!void` | the spectrum in the transform's own order, and `acc += scale·a⊙b` in it: fast convolution | none | none | pffft `pffft_zconvolve_accumulate` | yes |
| `StaticPlan(n, T)`: `scratchLen(_: Self) usize`, `transform(self: Self, io, data: []Complex(T), scratch: []Complex(T), direction: Direction, layout: BatchLayout, normalization: Normalization) TransformError!void`, `transformInto(self: Self, io, out: []Complex(T), in: []const Complex(T), scratch: []Complex(T), direction: Direction, layouts: BatchLayouts, normalization: Normalization) TransformIntoError!void` | `Plan` with `n` at compile time | none | none | PocketFFT, unrolled | stub (#625) |
| `StaticRealPlan(n, T)`: `scratchLen(_: Self) usize`, `forward(self: Self, io, out: []Complex(T), in: []const T, scratch: []T, layouts: BatchLayouts, normalization: Normalization) TransformIntoError!void`, `backward(self: Self, io, out: []T, in: []const Complex(T), scratch: []T, layouts: BatchLayouts, normalization: Normalization) TransformIntoError!void` | `RealPlan` with `n` at compile time | none | none | PocketFFT, unrolled | stub (#625) |
| `TrigTransformPlan(T)`: `init(allocator, n: usize, kind: TrigTransformKind) InitError!Self`, `deinit(self: *Self, allocator) void`, `scratchLen(self: Self) usize`, `transform(self: *const Self, io, data: []T, scratch: []T, layout: BatchLayout, normalization: Normalization) TransformError!void` | DCT and DST types I to IV | `vDSP_DCT_CreateSetup` (II, III, IV) | trigonometric transforms (`d_init_trig_transform`) | PocketFFT C++ `pocketfft_hdronly.hpp` `T_dcst4`, `T_dct1` (BSD-3) | stub (#625) |
| `ModifiedCosineTransform(T)`: `init(allocator, n: usize) InitError!Self`, `deinit(self: *Self, allocator) void`, `scratchLen(self: Self) usize`, `forward(self: *const Self, out: []T, in: []const T, scratch: []T) FrameError!void`, `backward(self: *const Self, out: []T, in: []const T, scratch: []T) FrameError!void` | modified DCT of `n` samples to `n/2` and back, by a DCT-IV | none | none | Princen-Bradley via `TrigTransformPlan(.dct4)` | stub (#625) |
| `MultiDimPlan(T)`: `init(allocator, shape: []const usize) InitError!Self`, `deinit(self: *Self, allocator) void`, `scratchLen(self: Self) usize`, `transform(self: *const Self, io, data: []Complex(T), scratch: []Complex(T), direction: Direction, normalization: Normalization) TransformError!void` | complex transform over every axis, column-major | `vDSP_fft2d_zip` (powers of two) | DFTI, multidimensional | PocketFFT C++ `c2c` over axes (BSD-3) | backlog |
| `binFrequencies(T, n: usize, d: T, out: []T) LengthError!void`, `realBinFrequencies(T, n: usize, d: T, out: []T) LengthError!void` | the frequencies of the bins | none | none | NumPy `numpy.fft.fftfreq`, `rfftfreq` | backlog |
| `centerZeroFrequency(T, x: []T) void`, `uncenterZeroFrequency(T, x: []T) void` | moves the zero frequency to the middle and back | none | none | NumPy `numpy.fft.fftshift` | backlog |

### `kernels.elementwise`

`Error = LengthError || Io.Cancelable` for all. One input: `y[i] = f(x[i])`, `y` may be `x`; two: `z[i] = f(x[i], y[i])`, `z` may be either. Each states its accuracy: the v0.1.0 set ([#614](https://github.com/randyos-org/rstd/issues/614)) is Zig's vector builtins, the accuracy of Zig's `std.math` (musl's algorithms), or exact. Real floats throughout; complex for the arithmetic and `abs`, the other functions' complex forms with the rest in [#626](https://github.com/randyos-org/rstd/issues/626).

| Signature | Computes | Accelerate | MKL | Source | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `add(io, T, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void` | `x + y`, exact rounding | `vDSP_vadd`, `vDSP_zvadd` | `v?Add` | `@Vector` | yes (#614) |
| `sub(io, T, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void` | `x − y` | `vDSP_vsub` (operands swapped: vDSP subtracts its first from its second) | `v?Sub` | `@Vector` | yes (#614) |
| `mul(io, T, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void` | `x · y` | `vDSP_vmul`, `vDSP_zvmul` | `v?Mul` | `@Vector` | yes (#614) |
| `div(io, T, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void` | `x / y` | `vvdiv`, `vDSP_zvdiv` | `v?Div` | `@Vector` | yes (#614) |
| `sqrt(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `√x`, correctly rounded | `vvsqrt` | `v?Sqrt` | `@sqrt` | yes (#614) |
| `exp(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `eˣ` | `vvexp` | `v?Exp` | `@exp` | yes (#614) |
| `exp2(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `2ˣ` | `vvexp2` | `v?Exp2` | `@exp2` | yes (#614) |
| `log(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `ln x` | `vvlog` | `v?Ln` | `@log` | yes (#614) |
| `log2(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `log₂ x` | `vvlog2` | `v?Log2` | `@log2` | yes (#614) |
| `log10(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `log₁₀ x` | `vvlog10` | `v?Log10` | `@log10` | yes (#614) |
| `sin(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `sin x` | `vvsin` | `v?Sin` | `@sin` | yes (#614) |
| `cos(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `cos x` | `vvcos` | `v?Cos` | `@cos` | yes (#614) |
| `tan(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `tan x` | `vvtan` | `v?Tan` | `@tan` | yes (#614) |
| `abs(io, T, x: ConstStrided(T), y: Strided(Float(T))) Error!void` | `\|x\|`; complex: the modulus, without overflow | `vvfabs`, `vDSP_zvabs` | `v?Abs` | `@abs`; complex by `std.math.hypot` | yes (#614) |
| `floor(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `⌊x⌋` | `vvfloor` | `v?Floor` | `@floor` | yes (#614) |
| `ceil(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `⌈x⌉` | `vvceil` | `v?Ceil` | `@ceil` | yes (#614) |
| `trunc(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | toward zero | `vvint` | `v?Trunc` | `@trunc` | yes (#614) |
| `round(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | to nearest, ties away from zero (Fortran's `ANINT`) | none (`vvnint` ties to even) | `v?Round` | `@round` | yes (#614) |
| `inverseSqrt(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `1/√x` | `vvrsqrt` | `v?InvSqrt` | `1 / @sqrt` | backlog |
| `reciprocal(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `1/x` | `vvrec` | `v?Inv` | `@Vector` | backlog |
| `cbrt(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `∛x` | none | `v?Cbrt` | `std.math.cbrt`; Julia `Base.Math.cbrt` (MIT) | backlog |
| `hypot(io, T, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void` | `√(x² + y²)`, without overflow | none | `v?Hypot` | Julia `Base.Math.hypot` (MIT) | backlog |
| `expm1(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `eˣ − 1`, accurate near 0 | `vvexpm1` | `v?Expm1` | `std.math.expm1`; Julia `Base.Math.expm1` | backlog |
| `log1p(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `ln(1 + x)`, accurate near 0 | `vvlog1p` | `v?Log1p` | `std.math.log1p` | backlog |
| `pow(io, T, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void` | `xʸ` | `vvpow` | `v?Pow` | `std.math.pow`; Julia `Base.Math.pow` | backlog |
| `powScalar(io, T, b: T, x: ConstStrided(T), y: Strided(T)) Error!void` | `xᵇ`, one exponent | `vvpows` | `v?Powx` | `std.math.pow` | backlog |
| `sincos(io, T, x: ConstStrided(T), s: Strided(T), c: Strided(T)) Error!void` | `sin x` and `cos x` in one reduction | `vvsincos` | `v?SinCos` | Julia `Base.Math.sincos` | backlog |
| `asin(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `arcsin x` | `vvasin` | `v?Asin` | `std.math.asin` | backlog |
| `acos(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `arccos x` | `vvacos` | `v?Acos` | `std.math.acos` | backlog |
| `atan(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `arctan x` | `vvatan` | `v?Atan` | `std.math.atan` | backlog |
| `atan2(io, T, y: ConstStrided(T), x: ConstStrided(T), z: Strided(T)) Error!void` | the angle of `(x, y)`, Fortran's `ATAN2(Y, X)` order | `vvatan2` | `v?Atan2` | `std.math.atan2` | backlog |
| `sinh(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `sinh x` | `vvsinh` | `v?Sinh` | `std.math.sinh` | backlog |
| `cosh(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `cosh x` | `vvcosh` | `v?Cosh` | `std.math.cosh` | backlog |
| `tanh(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `tanh x` | `vvtanh` | `v?Tanh` | `std.math.tanh` | backlog |
| `asinh(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `arsinh x` | `vvasinh` | `v?Asinh` | `std.math.asinh` | backlog |
| `acosh(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `arcosh x` | `vvacosh` | `v?Acosh` | `std.math.acosh` | backlog |
| `atanh(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `artanh x` | `vvatanh` | `v?Atanh` | `std.math.atanh` | backlog |
| `erf(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | the error function | none | `v?Erf` | openlibm `s_erf.c` (BSD), as SpecialFunctions.jl calls it; num.zig `erf` | backlog |
| `erfc(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `1 − erf x`, accurate for large `x` | none | `v?Erfc` | openlibm `s_erf.c` | backlog |
| `lgamma(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `ln \|Γ(x)\|` | none | `v?LGamma` | openlibm `e_lgamma_r.c`; SpecialFunctions.jl (MIT) | backlog |
| `tgamma(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `Γ(x)` | none | `v?TGamma` | openlibm `e_tgamma.c` (Fortran's `GAMMA`) | backlog |
| `angle(io, T, x: ConstStrided(T), y: Strided(Float(T))) Error!void` | complex `x`'s angle | `vDSP_zvphas` | `v?Arg` | `std.math.atan2(im, re)` | backlog |
| `conjugate(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | complex conjugate | `vDSP_zvconj` | `v?Conj` | `@Vector` | backlog |

### `kernels.dsp`

`Error` as in `elementwise`.

| Signature | Computes | Accelerate | MKL | Source | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `window(T, out: []T, kind: Window, symmetry: Symmetry) WindowError!void` | a window of `out.len` points | `vDSP_hann_window`, `vDSP_hamm_window`, `vDSP_blkman_window` | none | MATLAB `hann`, `hamming`, `blackman`, `kaiser` | yes |
| `addScalar(io, T, alpha: T, x: ConstStrided(T), y: Strided(T)) Error!void` | `x + α` | `vDSP_vsadd` | none | `@Vector` | backlog |
| `divideByScalar(io, T, alpha: T, x: ConstStrided(T), y: Strided(T)) Error!void` | `x / α` | `vDSP_vsdiv` | none | `@Vector` | backlog |
| `multiplyAdd(io, T, a: ConstStrided(T), b: ConstStrided(T), c: ConstStrided(T), d: Strided(T)) Error!void` | `a·b + c`, fused where the core has an instruction | `vDSP_vma` | `v?Fma` | `@mulAdd` | backlog |
| `fill(io, T, alpha: T, y: Strided(T)) Io.Cancelable!void` | every element `α` | `vDSP_vfill` | none | NumPy `full` | backlog |
| `ramp(io, T, start: T, step: T, y: Strided(T)) Io.Cancelable!void` | `start + i·step`, each computed, not accumulated | `vDSP_vramp` | none | NumPy `arange` | backlog |
| `linspace(io, T, start: T, stop: T, y: Strided(T)) Io.Cancelable!void` | `y.n` points from `start` to `stop`, both included | `vDSP_vgen` | none | NumPy `linspace` | backlog |
| `negate(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `−x` | `vDSP_vneg` | none | `@Vector` | backlog |
| `clip(io, T, lo: T, hi: T, x: ConstStrided(T), y: Strided(T)) Error!void` | `x` held to `[lo, hi]`; NaN stays NaN | `vDSP_vclip` | none | NumPy `clip` | backlog |
| `reverse(io, T, x: ConstStrided(T), y: Strided(T)) Error!void` | `y[i] = x[n − 1 − i]` | `vDSP_vrvrs` | none | NumPy `flip` | backlog |
| `convolve(io, T, mode: ConvolutionMode, x: ConstStrided(T), h: ConstStrided(T), y: Strided(T)) Error!void` | `y = x ∗ h`, `y.n` by `mode`: `n + m − 1`, `max(n, m)`, `max(n, m) − min(n, m) + 1` | `vDSP_conv` (`h` reversed) | `vslsConvExec` | NumPy `convolve`; direct, `RealPlan` products past a size | backlog |
| `correlate(io, T, mode: ConvolutionMode, x: ConstStrided(T), h: ConstStrided(T), y: Strided(T)) Error!void` | cross-correlation, NumPy's conjugation of `h` | `vDSP_conv` | `vslsCorrExec` | NumPy `correlate` | backlog |
| `biquad(T, sections: []const BiquadSection(T), state: []T, x: []const T, y: []T) Error!void` | cascaded second-order sections, direct form II transposed; `state` two per section, kept across calls | `vDSP_biquad` | none | SciPy `sosfilt` | backlog |
| `decimate(io, T, factor: usize, filter: []const T, x: ConstStrided(T), y: Strided(T)) Error!void` | FIR filter, then every `factor`th sample | `vDSP_desamp` | none | SciPy `upfirdn` | backlog |
| `cumulativeSum(T, x: []const T, y: []T) LengthError!void` | running sums, in order | `vDSP_vrsum` | none | NumPy `cumsum`, MATLAB `cumsum` | backlog |
| `evaluatePolynomial(io, T, coeffs: []const T, x: ConstStrided(T), y: Strided(T)) Error!void` | polynomial at each `x`, coefficients highest first, by Horner | `vDSP_vpoly` | none | NumPy `polyval`; num.zig `poly.val` | backlog |
| `interpolateLinear(io, T, xp: []const T, fp: []const T, x: ConstStrided(T), y: Strided(T)) Error!void` | linear interpolation of `(xp, fp)` at `x`, ends held | none | `dfdInterpolate1D` | NumPy `interp` | backlog |
| `absSquared(io, T, x: ConstStrided(Complex(T)), y: Strided(T)) Error!void` | `\|x\|²`, a power spectrum | `vDSP_zvmags` | none | Julia `abs2` | backlog |
| `toDecibels(io, T, scale: DecibelScale, reference: T, x: ConstStrided(T), y: Strided(T)) Error!void` | `10·log10(x/ref)` or `20·log10(x/ref)` | `vDSP_vdbcon` | none | librosa `power_to_db` | backlog |
| `splitComplex(io, T, z: ConstStrided(Complex(T)), re: Strided(T), im: Strided(T)) Error!void`, `joinComplex(io, T, re: ConstStrided(T), im: ConstStrided(T), z: Strided(Complex(T))) Error!void` | interleaved complex to split and back | `vDSP_ctoz`, `vDSP_ztoc` | none | `@Vector` shuffles | backlog |

### `kernels.stats`

| Signature | Computes | Accelerate | MKL | Source | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `kernels.blas.sum`, `sumWith`, `index`, `extreme` | sums and extremes: the BLAS's, not copies | `vDSP_sve`, `vDSP_maxv`, `vDSP_minv` | VS | rstd | yes |
| `deviationSums(io, T, s: Summation, m: T, x: ConstStrided(T)) Io.Cancelable!DeviationSums(T)` | `Σ(x − m)` and `Σ\|x − m\|²` in one pass; the same bits at any worker count | none | `vsldSSCompute` (`VSL_SS_2C_MOM`) | Chan, Golub and LeVeque's corrected two-pass; rstd's `fast.zig` lanes | yes (#605) |
| `centralMoments(io, T, s: Summation, m: T, x: ConstStrided(T)) Io.Cancelable!CentralMomentSums(T)` | `Σ(x − m)²`, `³`, `⁴`, for skewness and kurtosis | none | `vsldSSCompute` (`VSL_SS_4C_MOM`) | SciPy `stats.moment`; StatsBase `moment` | backlog |
| `RunningMoments(F)`: `push(self: *Self, x: F) void`, `mean(self: Self) error{EmptyInput}!F`, `variance(self: Self, normalize: Normalize) error{EmptyInput}!F` | Welford's running mean and variance | none | VS streaming | Welford 1962 | yes (#613) |
| `RunningMoments(F)`: `merge(self: *Self, other: Self) void` | combines two runs (Chan's formula), for threads | none | none | Chan, Golub and LeVeque 1979 | backlog |
| `histogramInto(T, edges: []const T, x: ConstStrided(T), counts: []u64) HistogramIntoError!void` | counts into the bins `edges` bound, NumPy's half-open bins, the last closed | `vImageHistogramCalculation` in part | none | NumPy `histogram` | yes (#613: `ConstStrided` in place of `anytype`, arguments in the kernel order) |
| `selectNth(T, x: []T, k: usize) T` | the `k`th smallest, `x` partitioned around it, in O(n); NaNs last | none | `vsldSSCompute` (`VSL_SS_QUANTS`) | NumPy `partition` (introselect); Rust `select_nth_unstable` (MIT) | backlog; `quantile` sorts a copy today |
| `medianAbsDeviation(allocator, T, x: []const T) MadError!T` | `median(\|x − median(x)\|)` | none | `vsldSSCompute` (`VSL_SS_MDAD`) | SciPy `stats.median_abs_deviation` | backlog |
| `RandomStream`: `init(seed: u64) RandomStream`, `jump(self: *RandomStream, chunks: u64) void` | a stream and its chunk's stream | none | `vslNewStream`, `vslSkipAheadStream` | Zig `std.Random.Xoshiro256` (MIT) | backlog |
| `uniform(io, T, rng: *RandomStream, lo: T, hi: T, y: Strided(T)) Io.Cancelable!void` | uniform on `[lo, hi)` | none | `v?RngUniform` | NumPy `Generator.uniform`; num.zig `random.uniform` | backlog |
| `normal(io, T, rng: *RandomStream, mean: T, sd: T, y: Strided(T)) Io.Cancelable!void` | normal, by the ziggurat | none | `v?RngGaussian` | Zig `std.Random.floatNorm` (Marsaglia and Tsang) | backlog |

### `rstd.math` (scalar)

| Signature | Computes | Source | v0.1.0 |
| --- | --- | --- | --- |
| `floor(x: anytype) @TypeOf(x)` | `⌊x⌋`; `@floor`, but exact for `f64` on the x87 (`x87_f64`) | rstd | yes (#608) |
| `ceil(x: anytype) @TypeOf(x)` | `⌈x⌉`, as `floor` | rstd | yes (#608) |
| `round(x: anytype) @TypeOf(x)` | to nearest, ties away from zero, as `floor` | rstd | yes (#608) |
| `rationalApprox(x: anytype, tol: ?@TypeOf(x)) Fraction(@TypeOf(x))` | the continued fraction of `x` within `tol` (MATLAB `rat`) | MATLAB `rat` | yes (#608) |
| `x87_f64: bool` | whether `f64` math runs on the x87 (32-bit x86 without SSE2), where compiler-rt's `f64` `@floor`, `@ceil` and `@round` keep fraction bits ([#264](https://github.com/randyos-org/rstd/issues/264)) | rstd | yes (#608) |

## The high-level API

- **A plain function with the best defaults and no options:** `stats.variance(io, xs, normalize)`, `fft` transforms, `Matrix` operations.
- **For advanced users, a sibling ending in `With`, taking a struct of options:** `stats.varianceWith(io, xs, normalize, .{ .summation = .compensated })`. Each field defaults to what the plain function uses, so a caller names only what they change; a field takes the kernel layer's own type (`blas.Summation`) and passes it through unchanged, so a user reaches the knobs without calling the kernel layer themselves.
- **No options on the plain function** (user, 2026-11-20): nobody passes `.{}` just to ask for the default behavior. That struct belongs only to the sibling ending in `With`. `fft`'s `Options` (`.norm`) and `stats.quantile`'s `QuantileOptions` (`.method`) move from the plain functions to new `*With` siblings in their moves ([#611](https://github.com/randyos-org/rstd/issues/611), [#613](https://github.com/randyos-org/rstd/issues/613)); today's `meanWith`, `varianceWith` and `stdDevWith` take a struct in place of their `Summation` argument ([#605](https://github.com/randyos-org/rstd/issues/605)). Build-time settings a root declares once (`rstd_linalg_options`) are not call options and stay.
- **No switch on a kernel option:** a high-level function that would have to switch on the summation mode or a method itself is a sign the algorithm belongs in a kernel ([#597](https://github.com/randyos-org/rstd/issues/597)): variance's deviation pass became a `kernels.stats` reduction ([#605](https://github.com/randyos-org/rstd/issues/605)).

## High-level functions

What users call: plain functions with the best defaults, and `*With` siblings taking an options struct where a knob is worth exposing ([The high-level API](#the-high-level-api)). Existing functions are given by their literal signatures, the shared arguments bare as in the Kernel API (`allocator`, `io`; `a: anytype` takes a `Mat`, `View`, `ConstView`, `StaticMat` or slice of the element type); functions that will exist are named; capabilities not yet designed are bullets, spec to be decided. Each module moves under `rstd.math` in its issue ([#611](https://github.com/randyos-org/rstd/issues/611) to [#615](https://github.com/randyos-org/rstd/issues/615)).

### `rstd.math.linalg`

MATLAB's names and argument order where MATLAB has the routine; NumPy's where it does not. Today's `dense` helpers named `*With` (`detWith`, `mulWith`, `sumWith` and the rest) take an execution path, not options; they are internal, and [#615](https://github.com/randyos-org/rstd/issues/615) renames them so `*With` means only the options siblings.

#### The `linalg` types

| Type | What | v0.1.0 |
| --- | --- | --- |
| `Mat(T)` | an owned column-major matrix: `init`, `zeros`, `deinit`, `at`, `col`, `sub`, `view`, `constView` | yes |
| `View(T)`, `ConstView(T)` | a strided window on a matrix, writable or not | yes |
| `StaticMat(R, C, T)`, `Vec(N, T)` | compile-time shapes, on Zig's `@Matrix` where the compiler has it (`_matrix/`) | yes |
| `Array(T)`, `StaticArray(dims, T)`, `ArrayView(T)`, `ConstArrayView(T)` | N-d arrays, column-major: `init`, `deinit`, `view`, `slice`, `reshape`, `permute`, `broadcastTo`, `asMatrix`, `count` | stub: types and storage only |
| `Lu(T)`, `Qr(T)`, `Svd(T)`, `Eig(T)`, `Hess(T)`, `Schur(T)`, `Rref(T)` | factorizations' results, each with `deinit` (`Lu` with `l` and `u`) | yes |
| `QrMode`, `SvdMode`, `MatrixNorm`, `VectorNorm` | economy or full factors; which norm | yes |
| `Options` | build-time `rstd_linalg_options` (`inline_limit`: up to which size the inline kernels run) | yes |

#### The `linalg` functions

The last column is `yes` for what v0.1.0 has, or the release that added the row.

| Signature | Computes | v0.1.0 |
| --- | --- | --- |
| `dot(io, x: anytype, y: anytype) VectorError!ElementOf(@TypeOf(x))` | `xᴴy` | yes |
| `norm(io, x: anytype) VectorError!Real(ElementOf(@TypeOf(x)))` | `‖x‖₂` | yes |
| `vectorNorm(io, x: anytype, which: VectorNorm) VectorNormError!Real(ElementOf(@TypeOf(x)))` | `‖x‖ₚ`, MATLAB's `vecnorm` | yes |
| `sum(io, x: anytype) VectorError!ElementOf(@TypeOf(x))` | `Σ xᵢ` | yes |
| `max(io, x: anytype) VectorError!ElementOf(@TypeOf(x))`, `min(io, x: anytype) VectorError!ElementOf(@TypeOf(x))` | the largest, the smallest | yes |
| `axpy(io, a: anytype, x: anytype, y: anytype) VectorError!void`, `scale(a: anytype, x: anytype) void` | `y ← ax + y`; `x ← ax` | yes |
| `mul(allocator, io, a: anytype, b: anytype) MulError!Mat(ElementOf(@TypeOf(a)))`, `mulInto(allocator, io, out: anytype, a: anytype, b: anytype) MulError!void` | `AB` | yes |
| `transpose(allocator, io, a: anytype) (Allocator.Error \|\| Io.Cancelable)!Mat(ElementOf(@TypeOf(a)))`, `ctranspose(...)`, `transposeInto(io, out: anytype, a: anytype) TransposeError!void`, `ctransposeInto(...)` | `Aᵀ`, `Aᴴ` | yes |
| `add(allocator, a: anytype, b: anytype) ElementwiseError!Mat(ElementOf(@TypeOf(a)))`, `sub`, `times`, `rdivide`, and their `*Into(out: anytype, a: anytype, b: anytype) ElementwiseIntoError!void` | elementwise `+`, `−`, `.*`, `./` (MATLAB's names) | yes |
| `map(allocator, a: anytype, f: anytype) Allocator.Error!Mat(ElementOf(@TypeOf(a)))`, `mapInto(out: anytype, a: anytype, f: anytype) ElementwiseIntoError!void` | `f` of each element | yes |
| `prod(a: anytype) ElementOf(@TypeOf(a))` | `Π aᵢⱼ` | yes |
| `realPart(allocator, a: anytype) Allocator.Error!Mat(Real(ElementOf(@TypeOf(a))))`, `imagPart(...)`, `conj(allocator, a: anytype) Allocator.Error!Mat(ElementOf(@TypeOf(a)))` | `real`, `imag`, `conj` | yes |
| `det(allocator, io, a: anytype) DetError!ElementOf(@TypeOf(a))` | `det(A)` | yes |
| `inv(allocator, io, a: anytype) InvError!Mat(ElementOf(@TypeOf(a)))`, `invInto(allocator, io, out: anytype, a: anytype) InvIntoError!void` | `A⁻¹` | yes |
| `solve(allocator, io, a: anytype, b: anytype) SolveError!Mat(ElementOf(@TypeOf(a)))`, `solveInto(allocator, io, out: anytype, a: anytype, b: anytype) SolveError!void` | `A \ B`: square by its structure (triangular, Cholesky, LU), else least squares | yes |
| `lu(allocator, io, a: anytype) (Allocator.Error \|\| Io.Cancelable)!Lu(ElementOf(@TypeOf(a)))` | `PA = LU` | yes |
| `luComplete(allocator, io, a: anytype) (Allocator.Error \|\| Io.Cancelable)!LuComplete(ElementOf(@TypeOf(a)))` | `PAQ = LU` with complete pivoting and its rank (Eigen's `FullPivLU`) | v0.2.0 |
| `chol(allocator, io, a: anytype) CholError!Mat(ElementOf(@TypeOf(a)))` | `A = RᴴR` | yes |
| `qr(allocator, io, a: anytype, mode: QrMode) (Allocator.Error \|\| Io.Cancelable)!Qr(ElementOf(@TypeOf(a)))` | `A = QR` | yes |
| `svd(allocator, io, a: anytype, mode: SvdMode) SvdError!Svd(ElementOf(@TypeOf(a)))` | `A = U·Σ·Vᴴ` | yes |
| `eig(allocator, io, a: anytype, want_vectors: bool) EigError!Eig(ElementOf(@TypeOf(a)))` | eigenvalues, and vectors; Hermitian `A` by the symmetric path | yes |
| `hess(allocator, io, a: anytype) HessError!Hess(ElementOf(@TypeOf(a)))` | `A = Q·H·Qᴴ` | yes |
| `schur(allocator, io, a: anytype) SchurError!Schur(ElementOf(@TypeOf(a)))`, `schurComplex(...)` | real or complex Schur form | yes |
| `ldl(allocator, io, a: anytype) LdlError!Ldl(ElementOf(@TypeOf(a)))` | `A = L·D·Lᴴ` by Bunch and Kaufman (MATLAB's `ldl`) | v0.2.0 |
| `cholupdate(allocator, r: anytype, x: anytype, op: CholupdateOp) CholupdateError!Mat(ElementOf(@TypeOf(r)))` | `R1ᴴR1 = RᴴR ± xxᴴ` (MATLAB's `cholupdate`) | v0.2.0 |
| `lsqminnorm(allocator, io, a: anytype, b: anytype) LsqminnormError!Mat(ElementOf(@TypeOf(a)))`, `lsqminnormWith(..., options: LsqminnormOptions(...))` | the minimum-norm least-squares solution | v0.2.0 |
| `balance(allocator, io, a: anytype, mode: BalanceMode) BalanceError!Balance(ElementOf(@TypeOf(a)))` | `T⁻¹AT` balanced (MATLAB's `balance`) | v0.2.0 |
| `eigGeneralized(allocator, io, a: anytype, b: anytype, want_vectors: bool) EigGeneralizedError!Eig(ElementOf(@TypeOf(a)))` | `eig(A, B)`: a Hermitian definite pencil by `sygvd`, any other by `ggev` | v0.2.0 |
| `qz(allocator, io, a: anytype, b: anytype) QzError!Qz(ElementOf(@TypeOf(a)))`, `qzComplex(...)` | `QAZ = AA`, `QBZ = BB` (MATLAB's `qz`) | v0.2.0 |
| `ordschur(allocator, io, u: anytype, t: anytype, select: OrdschurSelect) OrdschurError!Schur(ElementOf(@TypeOf(t)))`, `ordqz(allocator, io, aa: anytype, bb: anytype, q: anytype, z: anytype, select: OrdschurSelect) OrdqzError!Qz(ElementOf(@TypeOf(aa)))` | a Schur or generalized Schur form reordered by a mask, a region (`lhp`, `rhp`, `udi`, `udo`) or clusters | v0.2.0 |
| `sylvester(allocator, io, a: anytype, b: anytype, c: anytype) SylvesterError!Mat(ElementOf(@TypeOf(a)))`, `lyap(allocator, io, a: anytype, q: anytype)`, `lyapSylvester(allocator, io, a, b, c)`, `lyapGeneralized(allocator, io, a, q, e)` | `AX + XB = C`; `AX + XAᴴ + Q = 0`, `AX + XB + C = 0`, `AXEᴴ + EXAᴴ + Q = 0` | v0.2.0 |
| `matrixNorm(allocator, io, a: anytype, which: MatrixNorm) MatrixNormError!Real(ElementOf(@TypeOf(a)))` | `norm(A, p)` | yes |
| `cond(allocator, io, a: anytype) SvdError!Real(ElementOf(@TypeOf(a)))`, `rcond(allocator, io, a: anytype) RcondError!Real(ElementOf(@TypeOf(a)))` | the 2-norm condition number; its one-norm reciprocal estimate | yes |
| `rank(allocator, io, a: anytype, tol: ?Real(ElementOf(@TypeOf(a)))) SvdError!usize` | numerical rank | yes |
| `pinv(allocator, io, a: anytype, tol: ?Real(ElementOf(@TypeOf(a)))) SvdError!Mat(ElementOf(@TypeOf(a)))` | the pseudoinverse | yes |
| `orth(allocator, io, a: anytype, tol: ?Real(...)) SvdError!Mat(ElementOf(@TypeOf(a)))`, `nullSpace(...)` | bases of the range and the null space (MATLAB's `orth`, `null`) | yes |
| `expm(allocator, io, a: anytype) ExpmError!Mat(ElementOf(@TypeOf(a)))`, `logm`, `sqrtm` | matrix exponential, logarithm, square root | yes |
| `powm(allocator, io, a: anytype, p: f64) PowmError!Mat(ElementOf(@TypeOf(a)))` | `A^p` | yes |
| `funm(allocator, io, a: anytype, f: Derivatives(...)) FunmError!Mat(...)`, `cosm`, `sinm`, `atanm` | a matrix function by Schur-Parlett | yes |
| `identity(allocator, T, rows: usize, cols: usize) Allocator.Error!Mat(T)`, `ones(...)` | `eye`, `ones` | yes |
| `tril(allocator, T, a: ConstView(T), k: isize) Allocator.Error!Mat(T)`, `triu(...)`, `diagMatrix(allocator, T, v: []const T, k: isize) Allocator.Error!Mat(T)`, `diagonal(allocator, T, a: ConstView(T), k: isize) Allocator.Error![]T`, `kron(allocator, T, a: ConstView(T), b: ConstView(T)) Allocator.Error!Mat(T)` | `tril`, `triu`, `diag` both ways, `kron` | yes |
| `rref(allocator, T, a: ConstView(T), tol: ?Real(T)) Allocator.Error!Rref(T)` | reduced row echelon form | yes |
| `roots(allocator, io, T, c: []const T) RootsError![]Complex(T)`, `poly(allocator, T, r: []const T) Allocator.Error![]T`, `charPoly(allocator, io, T, a: ConstView(T)) CharPolyError![]T` | polynomial roots, the polynomial of given roots, the characteristic polynomial | yes |
| `magic(allocator, T, n: usize) Allocator.Error!Mat(T)`, `hilbert`, `invHilbert`; `uniform(allocator, T, rows: usize, cols: usize, rnd: std.Random) Allocator.Error!Mat(T)`, `normal`, `fillUniform(T, v: View(T), rnd: std.Random) void`, `fillNormal` | MATLAB's gallery and random matrices | yes |
| `trace`, `logDet`, `cross`, `outer` | `trace`, `log(det)` without overflow, `cross`, `xyᴴ` | named; #615 |
| `sumWith`, `dotWith`, `normWith` | the reductions with `.{ .summation = ... }` | named; #615 |
| `svdGeneralized` | `gsvd` over the backlog's `ggsvd3` | named; backlog (#623) |
| `polyfit`, `polyval` | least-squares fit; evaluation over `kernels.dsp.evaluatePolynomial` | named; backlog (#627) |

- N-d arrays beyond storage: indexing, broadcasting, reductions along a dimension. Spec to be decided.
- A sparse matrix type for users over `CompressedColumns`, with `mul`, `solve` and the iterative solvers behind one call. Spec to be decided.
- Partial eigenvalues and singular values of large and sparse matrices (MATLAB's `eigs`, `svds`). Spec to be decided.
- Which routines take a `*With` sibling, and its options, beyond the reductions. Spec to be decided; algorithm choices stay out (#554).

### `rstd.math.fft`

The plans (`Plan`, `RealPlan` and the stubs) are in `kernels.fft`; this layer, new in [#611](https://github.com/randyos-org/rstd/issues/611), plans per call. Errors: `FftError` (`OutOfMemory`, `ZeroLength`, `Canceled`), and `IrfftError` adds `LengthMismatch`.

| Signature | Computes | v0.1.0 |
| --- | --- | --- |
| `fft(allocator, io, T, x: []const Complex(T)) FftError![]Complex(T)`, `ifft` the same | complex transform of a vector and its inverse, allocating the result (NumPy's `fft.fft`, MATLAB's `fft`) | yes |
| `rfft(allocator, io, T, x: []const T) FftError![]Complex(T)`, `irfft(allocator, io, T, n: usize, x: []const Complex(T)) IrfftError![]T` | real input to the `n / 2 + 1` values of the half spectrum and back | yes |
| `fftWith`, `ifftWith`, `rfftWith`, `irfftWith` | the same with `options: FftOptions` last (`.normalization`) | yes |
| `binFrequencies`, `realBinFrequencies`, `centerZeroFrequency`, `uncenterZeroFrequency` | `fftfreq`, `rfftfreq`, `fftshift`, `ifftshift` | named; backlog (#625) |
| `dct`, `idct`, `dst`, `idst`, `mdct`, `imdct` | the trigonometric transforms, over the stubs | named; backlog (#625) |

- Whether the plain functions plan each call or keep a plan per length (a cache, or a plan the caller holds) for speed. Spec to be decided.
- 2-D and N-d transforms (`fft2`, `fftn`) over `kernels.fft.MultiDimPlan`. Spec to be decided (#625).

### `rstd.math.signal`

| Signature | Computes | v0.1.0 |
| --- | --- | --- |
| `window(T, out: []T, kind: Window, symmetry: Symmetry) WindowError!void` | a window, from `kernels.dsp` | yes |
| `Stft(T)`: `init(allocator, window: []const T, hop: usize, nfft: usize) InitError!Self`, `deinit(self: *Self, allocator) void`, `bins(self: Self) usize`, `frames(self: Self, len: usize) usize`, `scratchLen(self: Self) usize`, `forward(self: *const Self, io, out: []Complex(T), in: []const T, scratch: []T) Error!void`, `inverse(self: *const Self, io, out: []T, in: []const Complex(T), frame_count: usize, scratch: []T) Error!void` | the short-time Fourier transform and its inverse (overlap-add) | yes |
| `Convolver(T)`: `init(allocator, filter: []const T, block: usize) InitError!Self`, `deinit(self: *Self, allocator) void`, `reset(self: *Self) void`, `process(self: *Self, out: []T, in: []const T) Error!void` | streaming convolution with a long filter, uniformly partitioned | yes |
| `convolve`, `correlate` | whole-signal convolution and correlation, allocating, by direct sums or FFT products by size (`kernels.dsp.convolve` with `ConvolutionMode`) | named; backlog (#627) |
| `filter` | IIR filtering by second-order sections (SciPy's `sosfilt`, over `kernels.dsp.biquad`) | named; backlog (#627) |
| `spectrogram` | power spectra of a short-time Fourier transform, in decibels | named; backlog (#627) |
| `resample` | rational resampling (SciPy's `resample_poly`) | named; backlog (#627) |

- Filter design (Butterworth, Chebyshev, FIR windows: SciPy's `butter`, `firwin`). Spec to be decided.

### `rstd.math.stats`

| Signature | Computes | v0.1.0 |
| --- | --- | --- |
| `mean(io, R: type, xs: anytype) MomentError!R` | the mean, of any integer, float or complex elements, as `R` | yes |
| `meanWith(io, R: type, xs: anytype, options: MomentOptions) MomentError!R` | the mean with a `Summation` | yes (#605) |
| `variance(io, R: type, xs: anytype, normalize: Normalize) MomentError!R`, `varianceWith(io, R: type, xs: anytype, normalize: Normalize, options: MomentOptions) MomentError!R` | the corrected two-pass variance (#605) | yes (#605) |
| `stdDev(io, R: type, xs: anytype, normalize: Normalize) MomentError!R`, `stdDevWith(...)` | its square root | yes (#605) |
| `RunningMoments(F)`: `push(self: *Self, x: F) void`, `mean(self: Self) error{EmptyInput}!F`, `variance(self: Self, normalize: Normalize) error{EmptyInput}!F` | a running mean and variance | yes (#613) |
| `median(allocator, R: type, xs: anytype) MedianError!R` | the median; NaN if any | yes |
| `quantile(allocator, R: type, xs: anytype, p: []const R) QuantileError![]R`, `quantileInto(allocator, R: type, out: []R, xs: anytype, p: []const R) QuantileIntoError!void` | quantiles at `p`, NumPy's `linear` | yes (#613: options moved to `quantileWith`) |
| `quantileWith(allocator, R: type, xs: anytype, p: []const R, options: QuantileOptions) QuantileError![]R`, `quantileIntoWith(allocator, R: type, out: []R, xs: anytype, p: []const R, options: QuantileOptions) QuantileIntoError!void` | quantiles by any of NumPy's nine methods | yes (#613) |
| `histogram(allocator, F: type, xs: anytype, bins: Bins(F)) HistogramError!Histogram(F)`, `histogramInto(F: type, counts: []u64, xs: anytype, edges: []const F) HistogramIntoError!void` | counts into bins | yes |
| `cov(allocator, io, x: anytype) CovError!Mat(...)`, `corrcoef(allocator, io, x: anytype) CovError!Mat(...)` | covariance and correlation of a matrix's columns (MATLAB's `cov`, `corrcoef`) | yes |
| `skewness`, `kurtosis`, `medianAbsDeviation`, `zscore` | over `kernels.stats` | named; backlog (#628) |
| `uniform`, `normal` | random samples into a slice, over `kernels.stats.RandomStream` | named; backlog (#628) |

- Weighted statistics, and statistics along a dimension of a matrix or N-d array. Spec to be decided.
- Distributions (probability densities, cumulative and inverse distribution functions: SciPy's `stats.norm` and the like). Spec to be decided.

### Scalar functions

`floor`, `ceil`, `round`, `rationalApprox`, `Fraction` and `x87_f64`, as in the Kernel API's last table; the scalar functions users call directly.

## Shared infrastructure

Today inside `blas`; it moves up to `rstd.math.kernels` when a second domain uses it:

- **CPU detection and core selection** (`blas/cpu/`): one detected core per process, its kernels per domain.
- **Threading** (`blas/threading/`, `threading.zig`): worker teams sized per call, chunks by size alone, results folded in order, so bits do not depend on the worker count.
- **Vector lanes** (`blas/vector/`, `blas1/fast.zig`): the `@Vector` accumulation patterns the reductions share, first needed outside `blas` by [#605](https://github.com/randyos-org/rstd/issues/605).
- **Configuration** (`Config.zig`): the environment's settings, read once.

## Backends

- **Rstd's own kernels, everywhere.** RandyOS's targets (bare-metal x86_64 and aarch64, UEFI, 32-bit Arm, the PowerPC 750) and Linux have no vendor library, and the kernels are the port's source.
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
