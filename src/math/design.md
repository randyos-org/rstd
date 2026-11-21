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

Every routine and type of `rstd.math.kernels`, with its literal signature ([#607](https://github.com/randyos-org/rstd/issues/607)). Signatures of what exists are taken from the code; the rest are this design's, and the code follows them. Each wrapper here is the only entry point: it runs rstd's kernels or, where the process routes it, the vendor's.

**Shape.** `name(allocator, io, comptime T, knobs, scalars, inputs, outputs)`: the allocator only where the routine allocates, the `Io` where it may thread, the element type `T` (`f16` to `f128`, `std.math.Complex` of them where the routine takes complex), then the routine's choices (`Trans`, `Uplo`, a `Summation`), scalars, inputs, outputs. Vectors are `ConstStrided(T)` in and `Strided(T)` out; matrices are a column-major slice with its leading dimension (`a: []const T, lda: usize`) and their dimensions as `usize`. No options structs in this layer. An output aliases an input only where its row says so.

**Errors.** Each routine returns an exact error set named for it (`PotrfError`) or its level (`blas2.Error`): bad input (`SliceTooShort`, `CountMismatch`, `LeadingDimensionTooSmall`) is an error, `error.OutOfMemory` where it allocates, `error.Canceled` where it threads; it panics only for a broken invariant.

**Columns.** "Accelerate" and "MKL" name the vendor routine a wrapper maps to; "(routed)" marks the calls rstd sends to Accelerate today, the rest are mapped when that backend comes ([#617](https://github.com/randyos-org/rstd/issues/617) for MKL). "Source" names what rstd's implementation is ported from, or will be. "v0.1.0": "yes" exists; "stub" has its API and returns `error.NotImplemented`; an issue number is planned for v0.1.0; "backlog" is named here so the API has room for it, each domain's backlog in one issue: `lapack` [#623](https://github.com/randyos-org/rstd/issues/623), `sparse` [#624](https://github.com/randyos-org/rstd/issues/624), `fft` [#625](https://github.com/randyos-org/rstd/issues/625), `elementwise` [#626](https://github.com/randyos-org/rstd/issues/626), `dsp` [#627](https://github.com/randyos-org/rstd/issues/627), `stats` [#628](https://github.com/randyos-org/rstd/issues/628).

**References for what is new:** NumPy and SciPy (names, modes, defaults), MATLAB (argument order), Julia's `LinearAlgebra`, `Statistics`, `StatsBase`, `DSP.jl` and `SpecialFunctions.jl` (algorithms, MIT), Fortran's intrinsics, Accelerate's vDSP and vForce, MKL's VM and VS, faer (sparse solvers, MIT), and muhammad-fiaz's num.zig ([github.com/muhammad-fiaz/num.zig](https://github.com/muhammad-fiaz/num.zig), MIT: a starting point for some backlog routines, its algorithms unoptimized).

### Types

Each type the kernel layer defines or takes, after the moves; the high-level option structs are at the end.

| Type | Definition | Module | v0.1.0 |
| --- | --- | --- | --- |
| `Strided(T)`, `ConstStrided(T)` | `struct { data: []T, n: usize, stride: isize }` (`[]const T` for `ConstStrided`): `n` elements, `stride` apart, negative backward | `kernels.blas` | yes |
| `BF16` | `enum(u16) { _ }`: a bfloat16's bits | `kernels.blas` | yes |
| `Trans` | `enum { n, t, r, c, o, u, s, d }`: none, transpose, conjugate, conjugate transpose; the rest OpenBLAS's conjugated-`x` letters | `kernels.blas` | yes |
| `Uplo` | `enum { u, l, v, m }`: upper, lower; `v`, `m` OpenBLAS's reversed forms | `kernels.blas` | yes |
| `Diag` | `enum { u, n }`: unit, non-unit | `kernels.blas` | yes |
| `Side` | `enum { l, r }` | `kernels.blas` | yes |
| `Extreme` | `enum { amax, amin, max, min }`: by magnitude (`\|re\| + \|im\|`) or by value | `kernels.blas` | yes |
| `Summation` | `enum { plain, compensated, fast }`: `f32` in `f64`; Neumaier's; `f32` in `f32` lanes | `kernels.blas` | yes |
| `HalfGemv` | `enum { sbgemv, bgemv, shgemv }`: bfloat16 in, `f32` out; bfloat16 in and out; `f16` in, `f32` out | `kernels.blas` | yes |
| `HalfGemm` | `enum { sbgemm, bgemm, shgemm }`: as `HalfGemv` | `kernels.blas` | yes |
| `Group(In, Out)` | `struct { transa: Trans = .n, transb: Trans = .n, m: usize, n: usize, k: usize, alpha: Out, a: []const []const In, lda: usize, b: []const []const In, ldb: usize, beta: Out, c: []const []Out, ldc: usize }`: products of one shape for `gemmBatch` | `kernels.blas` | yes |
| `Pivot` | `packed struct(usize) { row: Row, block: enum(u1) { one, two } }`, `Row` the unsigned integer one bit narrower than `usize` (`u63` on 64-bit): a Bunch-Kaufman pivot | `kernels.lapack` | yes |
| `Problem` | `enum { ax_lbx, abx_lx, bax_lx }`: `Ax = λBx`, `ABx = λx`, `BAx = λx` | `kernels.lapack` | yes |
| `Jobz` | `enum { n, v }`: values only, or vectors too | `kernels.lapack` | yes |
| `Job` | `enum { a, s, o, n }`: all, the leading, overwrite `A`, none (SVD vectors) | `kernels.lapack` | yes |
| `Norm` | `enum { m, o, i, f }`: max abs, one, infinity, Frobenius | `kernels.lapack` | yes |
| `Eigenvalues(T)` | `[]T` for complex `T`; `struct { re: []T, im: []T }` for real | `kernels.lapack` | yes |
| `GeesSelect(T)` | `*const fn (T) bool` for complex; `*const fn (Float(T), Float(T)) bool` for real | `kernels.lapack` | yes |
| `HseqrJob`, `Compz`, `Compq` | `enum { e, s }`; `enum { n, i, v }`; `enum { n, v }`: eigenvalues or Schur form; Schur vectors none, from identity, accumulated | `kernels.lapack` | yes |
| `TrevcSide`, `Howmny` | `enum { r, l, b }`; `enum { a, b, s }`: right, left, both; all, back-transformed, selected | `kernels.lapack` | yes |
| `BalanceJob`, `BalanceBlock` | `enum { n, p, s, b }`; `struct { lo: usize, hi: usize }` | `kernels.lapack` | yes |
| `TrsenJob`, `Reordered(T)` | `enum { n, e, v, b }`; `struct { m: usize, s: ?Float(T) = null, sep: ?Float(T) = null }` | `kernels.lapack` | yes |
| `Solved(F)`, `SolvedWithDif(F)` | `struct { scale: F, perturbed: bool }`; `struct { scale: F, dif: ?F = null, perturbed: bool }` | `kernels.lapack` | yes |
| `Range(F)` | `union(enum) { all, values: struct { lo: F, hi: F }, indices: struct { first: usize, last: usize } }`: which eigenvalues or singular values | `kernels.lapack` | backlog (#623) |
| `Fact`, `Equed` | `enum { n, f, e }`; `enum { n, r, c, b }`: `gesvx`'s factor and equilibration | `kernels.lapack` | backlog (#623) |
| `Csc(T)`, `Csr(T)` | `struct { rows: u32, cols: u32, col_start: []u32, row: []u32, value: []T }` (`Csr`: `row_start`, `col`) | `kernels.sparse` | yes |
| `Triplets(T)` | `struct { rows: u32, cols: u32, i: ArrayList(u32), j: ArrayList(u32), x: ArrayList(T) }`: assembly | `kernels.sparse` | yes |
| `Vector(T)` | `struct { len: usize, index: []const u32, value: []const T }`: a sparse column or row | `kernels.sparse` | yes |
| `Operator(T)` | `struct { ctx: *const anyopaque, apply: *const fn (ctx: *const anyopaque, io: Io, x: []const T, y: []T) Io.Cancelable!void }`: `y ← Ax` | `kernels.sparse` | yes |
| `Stats` | `struct { iterations: usize, residual: f64 }`: an iterative solve's | `kernels.sparse` | yes |
| `SparseCholesky(T)`, `SparseLdlt(T)`, `SparseLu(T)`, `SparseQr(T)` | `struct { perm: []u32, parent: []u32, l: Csc(T) }` (Cholesky: fill-reducing order, elimination tree, factor); `SparseLdlt` adds `d: []T`; `SparseLu`: `struct { row_perm: []u32, col_perm: []u32, l: Csc(T), u: Csc(T) }`; `SparseQr`: `struct { col_perm: []u32, v: Csc(T), tau: []T, r: Csc(T) }`. Analyzed once from the pattern, refactored per values | `kernels.sparse` | backlog (#624) |
| `Jacobi(T)`, `Ilu0(T)` | `struct { inv_diag: []T }`; `struct { lu: Csr(T) }`: preconditioners, each with `operator(self) Operator(T)` | `kernels.sparse` | backlog (#624) |
| `Direction`, `Norm` (fft) | `enum { forward, backward }`; `enum { backward, ortho, forward }`: which direction scales by `1/n` | `kernels.fft` | yes |
| `Layout`, `Layouts` | `struct { count: usize = 1, stride: usize = 1, distance: ?usize = null }`; `struct { in: Layout = .{}, out: Layout = .{} }`: a batch's shape in memory | `kernels.fft` | yes |
| `Plan(T)`, `RealPlan(T)` | `struct { n: usize, scratch_len: usize, kind: Kind }`: a complex or real transform of length `n`, immutable, shared by threads | `kernels.fft` | yes |
| `StaticPlan(n, T)`, `StaticRealPlan(n, T)` | `struct {}`: `Plan` with `n` fixed at compile time, no allocation | `kernels.fft` | stub (#625) |
| `TrigKind`, `TrigPlan(T)` | `enum { dct1, dct2, dct3, dct4, dst1, dst2, dst3, dst4 }`; `struct { n: usize, kind: TrigKind, scratch_len: usize }` | `kernels.fft` | stub (#625) |
| `Mdct(T)` | `struct { n: usize, scratch_len: usize }` | `kernels.fft` | stub (#625) |
| `PlanNd(T)` | `struct { shape: []const usize, plans: []Plan(T), scratch_len: usize }`: complex transforms over every axis | `kernels.fft` | backlog (#625) |
| `Window`, `Symmetry` | `union(enum) { hann, hamming, blackman, kaiser: f64, sine, kbd: f64 }`; `enum { symmetric, periodic }` | `kernels.dsp` | yes |
| `ConvMode` | `enum { full, same, valid }`: NumPy's `convolve` modes | `kernels.dsp` | backlog (#627) |
| `Biquad(T)` | `struct { b0: T, b1: T, b2: T, a1: T, a2: T }`: one second-order section, `a0` 1 (SciPy's `sos` row) | `kernels.dsp` | backlog (#627) |
| `DbScale` | `enum { power, amplitude }`: `10·log10` or `20·log10` | `kernels.dsp` | backlog (#627) |
| `Deviations(T)` | `struct { sum: T, sum_sq: Float(T) }`: `Σ(x − m)` and `Σ\|x − m\|²` | `kernels.stats` | #605 |
| `CentralMoments(F)` | `struct { m2: F, m3: F, m4: F }`: `Σ(x − m)²`, `Σ(x − m)³`, `Σ(x − m)⁴` | `kernels.stats` | backlog (#628) |
| `Running(F)` | `struct { count: u64 = 0, avg: F = 0, m2: F = 0 }`: Welford's state | `kernels.stats` | yes |
| `Bins(F)`, `Histogram(F)` | `union(enum) { count: usize, edges: []const F }`; `struct { counts: []u64, edges: []F }` | `kernels.stats` | yes |
| `Rng` | `struct { state: std.Random.Xoshiro256 }`: one stream; a chunk's stream is the seed's jumped by its index, so the bits do not depend on the worker count | `kernels.stats` | backlog (#628) |
| `Fraction(F)` | `struct { numerator: F, denominator: F }` | `rstd.math` | yes |
| `MomentOptions` | `struct { summation: Summation = .plain }`: `meanWith`, `varianceWith`, `stdDevWith` | `rstd.math.stats` | #605 |
| `QuantileOptions` | `struct { method: QuantileMethod = .linear }`: `quantileWith`; `QuantileMethod` NumPy's nine | `rstd.math.stats` | yes, moved by #613 |
| `FftOptions` | `struct { norm: Norm = .backward }`: the `*With` transforms (today's `fft.Options`) | `rstd.math.fft` | yes, moved by #611 |

### `kernels.blas`

The BLAS levels 1 to 3 and rstd's extensions, as `src/linalg/blas/` has them today. Every routine is generic over `T` (`f16` to `f128`, and complex); a CBLAS or MKL name with `?` stands for its letters (`s`, `d`, `c`, `z`).

#### Level 1

| Signature | Computes | Accelerate | MKL |
| --- | --- | --- | --- |
| `pub fn copy(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) PairError!void` | `y ← x` | `cblas_?copy` (routed) | `cblas_?copy` |
| `pub fn swap(io: Io, comptime T: type, x: Strided(T), y: Strided(T)) PairError!void` | `x ↔ y` | `cblas_?swap` (routed) | `cblas_?swap` |
| `pub fn scal(io: Io, comptime T: type, alpha: T, x: Strided(T)) Io.Cancelable!void` | `x ← αx` | `cblas_?scal` (routed) | `cblas_?scal` |
| `pub fn scalReal(io: Io, comptime T: type, alpha: Float(T), x: Strided(T)) Io.Cancelable!void` | `x ← αx`, real `α`, complex `x` | `cblas_csscal`, `cblas_zdscal` (routed) | `cblas_csscal`, `cblas_zdscal` |
| `pub fn axpy(io: Io, comptime T: type, alpha: T, x: ConstStrided(T), y: Strided(T)) PairError!void` | `y ← αx + y` | `cblas_?axpy` (routed) | `cblas_?axpy` |
| `pub fn axpyc(io: Io, comptime T: type, alpha: T, x: ConstStrided(T), y: Strided(T)) PairError!void` | `y ← α·conj(x) + y` | none | none |
| `pub fn axpby(io: Io, comptime T: type, alpha: T, x: ConstStrided(T), beta: T, y: Strided(T)) LengthError!void` | `y ← αx + βy` | `catlas_?axpby` (routed) | `cblas_?axpby` |
| `pub fn dot(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T)) PairError!T` | `Σ xᵢyᵢ` | `cblas_?dot`, `cblas_?dotu_sub` (routed) | `cblas_?dot`, `cblas_?dotu_sub` |
| `pub fn dotc(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T)) PairError!T` | `Σ conj(xᵢ)yᵢ` | `cblas_?dotc_sub` (routed) | `cblas_?dotc_sub` |
| `pub fn dotWith(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), s: Summation) PairError!T` | `dot`, summed as `s` says | none | none |
| `pub fn dotcWith(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), s: Summation) PairError!T` | `dotc`, summed as `s` says | none | none |
| `pub fn dsdot(io: Io, x: ConstStrided(f32), y: ConstStrided(f32)) PairError!f64` | `Σ xᵢyᵢ` of `f32`, in `f64` | `cblas_dsdot` (routed) | `cblas_dsdot` |
| `pub fn sdsdot(io: Io, alpha: f32, x: ConstStrided(f32), y: ConstStrided(f32)) PairError!f32` | `α + Σ xᵢyᵢ` in `f64`, rounded to `f32` | `cblas_sdsdot` (routed) | `cblas_sdsdot` |
| `pub fn nrm2(io: Io, comptime T: type, x: ConstStrided(T)) Float(T)` | `‖x‖₂`, without overflow or underflow | `cblas_?nrm2` (routed) | `cblas_?nrm2` |
| `pub fn asum(io: Io, comptime T: type, x: ConstStrided(T)) Io.Cancelable!Float(T)` | `Σ (\|re xᵢ\| + \|im xᵢ\|)` | `cblas_?asum` (routed) | `cblas_?asum` |
| `pub fn sum(io: Io, comptime T: type, x: ConstStrided(T)) Io.Cancelable!Float(T)` | `Σ xᵢ` (complex: `Σ (re xᵢ + im xᵢ)`, OpenBLAS's `?sum`) | none | none |
| `pub fn asumWith(io: Io, comptime T: type, x: ConstStrided(T), s: Summation) Io.Cancelable!Float(T)` | `asum`, summed as `s` says | none | none |
| `pub fn sumWith(io: Io, comptime T: type, x: ConstStrided(T), s: Summation) Io.Cancelable!Float(T)` | `sum`, summed as `s` says | none | none |
| `pub fn index(io: Io, comptime T: type, comptime which: Extreme, x: ConstStrided(T)) Io.Cancelable!?usize` | position of the `which` extreme, the first on ties; null when `x` is empty | `cblas_i?amax` (not routed: rstd's NaN rule) | `cblas_i?amax`, `cblas_i?amin` |
| `pub fn extreme(io: Io, comptime T: type, comptime which: Extreme, x: ConstStrided(T)) Io.Cancelable!?Float(T)` | value of the `which` extreme; null when `x` is empty | none | none |
| `pub fn rot(io: Io, comptime T: type, x: Strided(T), y: Strided(T), c: Float(T), s: Float(T)) PairError!void` | applies the plane rotation `(c, s)` to `(x, y)` | `cblas_?rot` (routed) | `cblas_?rot` |
| `pub fn rotm(io: Io, comptime T: type, x: Strided(T), y: Strided(T), param: *const [5]T) LengthError!void` | applies the modified Givens rotation `param` | `cblas_?rotm` (routed) | `cblas_?rotm` |
| `pub fn rotg(comptime T: type, a: *T, b: *T, c: *T, s: *T) void` | the rotation zeroing `b`: `c`, `s`, and `r` in `a` | `cblas_?rotg` (not routed) | `cblas_?rotg` |
| `pub fn zrotg(comptime T: type, a: *T, b: T, c: *Float(T), s: *T) void` | `rotg` for complex | `cblas_?rotg` (not routed) | `cblas_?rotg` |
| `pub fn rotmg(comptime T: type, d1: *T, d2: *T, x1: *T, y1: T, param: *[5]T) void` | the modified Givens rotation zeroing `y1` | `cblas_?rotmg` (not routed) | `cblas_?rotmg` |
| `pub fn cabs1(comptime T: type, z: T) Float(T)` | `\|re z\| + \|im z\|` | none | none |
| `pub fn sbdot(io: Io, x: ConstStrided(BF16), y: ConstStrided(BF16)) LengthError!f32` | `Σ xᵢyᵢ` of bfloat16, in `f32` | none | none |
| `pub fn bscal(io: Io, alpha: BF16, x: Strided(BF16)) void` | `x ← αx`, bfloat16 | none | none |
| `pub fn sbtobf16(io: Io, comptime F: type, in: ConstStrided(F), out: Strided(BF16)) LengthError!void` | `f32` or `f64` to bfloat16, rounded to nearest even | none | none |
| `pub fn bf16to(io: Io, comptime F: type, in: ConstStrided(BF16), out: Strided(F)) LengthError!void` | bfloat16 to `f32` or `f64` | none | none |

#### Level 2

| Signature | Computes | Accelerate | MKL |
| --- | --- | --- | --- |
| `pub fn gemv(allocator: Allocator, io: Io, comptime T: type, trans: Trans, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `y ← α·op(A)x + βy` | `cblas_?gemv` (routed) | `cblas_?gemv` |
| `pub fn gemvHalf(io: Io, comptime kind: HalfGemv, trans: Trans, alpha: kind.Out(), a: []const kind.In(), lda: usize, x: ConstStrided(kind.In()), beta: kind.Out(), y: Strided(kind.Out())) HalfError!void` | `gemv` with `f16` or bfloat16 `A` and `x` (`kind`) | none | none |
| `pub fn gbmv(allocator: Allocator, io: Io, comptime T: type, trans: Trans, kl: usize, ku: usize, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `gemv` for band `A` (`kl`, `ku`) | `cblas_?gbmv` (routed) | `cblas_?gbmv` |
| `pub fn ger(allocator: Allocator, io: Io, comptime T: type, alpha: T, x: ConstStrided(T), y: ConstStrided(T), a: []T, lda: usize) Error!void` | `A ← αxyᵀ + A` | `cblas_?ger`, `cblas_?geru` (routed) | `cblas_?ger`, `cblas_?geru` |
| `pub fn gerc(allocator: Allocator, io: Io, comptime T: type, alpha: T, x: ConstStrided(T), y: ConstStrided(T), a: []T, lda: usize) Error!void` | `A ← αxyᴴ + A` | `cblas_?gerc` (routed) | `cblas_?gerc` |
| `pub fn symv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `y ← αAx + βy`, symmetric `A` | `cblas_?symv`, real (routed) | `cblas_?symv`, real |
| `pub fn hemv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `y ← αAx + βy`, Hermitian `A` | `cblas_?hemv` (routed) | `cblas_?hemv` |
| `pub fn sbmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, k: usize, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `symv` for symmetric band `A` | `cblas_?sbmv`, real (routed) | `cblas_?sbmv` |
| `pub fn hbmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, k: usize, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `hemv` for Hermitian band `A` | `cblas_?hbmv` (routed) | `cblas_?hbmv` |
| `pub fn spmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, ap: []const T, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `symv` for packed `A` | `cblas_?spmv`, real (routed) | `cblas_?spmv` |
| `pub fn hpmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, ap: []const T, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void` | `hemv` for packed `A` | `cblas_?hpmv` (routed) | `cblas_?hpmv` |
| `pub fn trmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, a: []const T, lda: usize, x: Strided(T)) Error!void` | `x ← op(A)x`, triangular `A` | `cblas_?trmv` (routed) | `cblas_?trmv` |
| `pub fn trsv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, a: []const T, lda: usize, x: Strided(T)) Error!void` | `x ← op(A)⁻¹x` | `cblas_?trsv` (routed) | `cblas_?trsv` |
| `pub fn tbmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, k: usize, a: []const T, lda: usize, x: Strided(T)) Error!void` | `trmv` for triangular band `A` | `cblas_?tbmv` (routed) | `cblas_?tbmv` |
| `pub fn tbsv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, k: usize, a: []const T, lda: usize, x: Strided(T)) Error!void` | `trsv` for triangular band `A` | `cblas_?tbsv` (routed) | `cblas_?tbsv` |
| `pub fn tpmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, ap: []const T, x: Strided(T)) Error!void` | `trmv` for packed `A` | `cblas_?tpmv` (routed) | `cblas_?tpmv` |
| `pub fn tpsv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, ap: []const T, x: Strided(T)) Error!void` | `trsv` for packed `A` | `cblas_?tpsv` (routed) | `cblas_?tpsv` |
| `pub fn syr(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, x: ConstStrided(T), a: []T, lda: usize) Error!void` | `A ← αxxᵀ + A`, symmetric | `cblas_?syr`, real (routed) | `cblas_?syr` |
| `pub fn her(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: Float(T), x: ConstStrided(T), a: []T, lda: usize) Error!void` | `A ← αxxᴴ + A`, Hermitian, real `α` | `cblas_?her` (routed) | `cblas_?her` |
| `pub fn spr(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, x: ConstStrided(T), ap: []T) Error!void` | `syr` for packed `A` | `cblas_?spr`, real (routed) | `cblas_?spr` |
| `pub fn hpr(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: Float(T), x: ConstStrided(T), ap: []T) Error!void` | `her` for packed `A` | `cblas_?hpr` (routed) | `cblas_?hpr` |
| `pub fn syr2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, x: ConstStrided(T), y: ConstStrided(T), a: []T, lda: usize) Error!void` | `A ← αxyᵀ + αyxᵀ + A` | `cblas_?syr2`, real (routed) | `cblas_?syr2` |
| `pub fn her2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, x: ConstStrided(T), y: ConstStrided(T), a: []T, lda: usize) Error!void` | `A ← αxyᴴ + conj(α)yxᴴ + A` | `cblas_?her2` (routed) | `cblas_?her2` |
| `pub fn spr2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, x: ConstStrided(T), y: ConstStrided(T), ap: []T) Error!void` | `syr2` for packed `A` | `cblas_?spr2`, real (routed) | `cblas_?spr2` |
| `pub fn hpr2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, x: ConstStrided(T), y: ConstStrided(T), ap: []T) Error!void` | `her2` for packed `A` | `cblas_?hpr2` (routed) | `cblas_?hpr2` |

#### Level 3

| Signature | Computes | Accelerate | MKL |
| --- | --- | --- | --- |
| `pub fn gemm(allocator: Allocator, io: Io, comptime T: type, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `C ← α·op(A)op(B) + βC` | `cblas_?gemm` (routed) | `cblas_?gemm` |
| `pub fn gemm3m(allocator: Allocator, io: Io, comptime T: type, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `gemm` for complex, three real products | none | `cblas_?gemm3m` |
| `pub fn gemmt(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, transa: Trans, transb: Trans, m: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `gemm`, only the `uplo` triangle of `C` | none | `cblas_?gemmt` |
| `pub fn gemmHalf(allocator: Allocator, io: Io, comptime kind: HalfGemm, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: kind.Out(), a: []const kind.In(), lda: usize, b: []const kind.In(), ldb: usize, beta: kind.Out(), c: []kind.Out(), ldc: usize) Error!void` | `gemm` with `f16` or bfloat16 inputs (`kind`) | none | `cblas_hgemm`, `cblas_gemm_bf16bf16f32` |
| `pub fn gemmtHalf(io: Io, comptime kind: HalfGemm, uplo: Uplo, transa: Trans, transb: Trans, m: usize, k: usize, alpha: kind.Out(), a: []const kind.In(), lda: usize, b: []const kind.In(), ldb: usize, beta: kind.Out(), c: []kind.Out(), ldc: usize) HalfError!void` | `gemmt` with `f16` or bfloat16 inputs | none | none |
| `pub fn symm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `C ← αAB + βC` or `αBA + βC`, symmetric `A` | `cblas_?symm` (routed) | `cblas_?symm` |
| `pub fn hemm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `symm`, Hermitian `A` | `cblas_?hemm` (routed) | `cblas_?hemm` |
| `pub fn symm3m(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `symm` for complex, three real products | none | none |
| `pub fn hemm3m(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `hemm` for complex, three real products | none | none |
| `pub fn syrk(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: T, a: []const T, lda: usize, beta: T, c: []T, ldc: usize) Error!void` | `C ← α·op(A)op(A)ᵀ + βC`, one triangle | `cblas_?syrk` (routed) | `cblas_?syrk` |
| `pub fn herk(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: Float(T), a: []const T, lda: usize, beta: Float(T), c: []T, ldc: usize) Error!void` | `C ← α·op(A)op(A)ᴴ + βC`, real `α`, `β` | `cblas_?herk` (routed) | `cblas_?herk` |
| `pub fn syr2k(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void` | `C ← α·op(A)op(B)ᵀ + α·op(B)op(A)ᵀ + βC` | `cblas_?syr2k` (routed) | `cblas_?syr2k` |
| `pub fn her2k(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: Float(T), c: []T, ldc: usize) Error!void` | `C ← α·op(A)op(B)ᴴ + conj(α)·op(B)op(A)ᴴ + βC` | `cblas_?her2k` (routed) | `cblas_?her2k` |
| `pub fn trmm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, transa: Trans, diag: Diag, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []T, ldb: usize) Error!void` | `B ← α·op(A)B` or `α·B·op(A)`, triangular `A` | `cblas_?trmm` (routed) | `cblas_?trmm` |
| `pub fn trsm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, transa: Trans, diag: Diag, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []T, ldb: usize) Error!void` | `B ← α·op(A)⁻¹B` or `α·B·op(A)⁻¹` | `cblas_?trsm` (routed) | `cblas_?trsm` |
| `pub fn gemmBatch(allocator: Allocator, io: Io, comptime T: type, groups: []const Group(T, T)) BatchError!void` | `gemm` over every product of every `Group` | none | `cblas_?gemm_batch` |
| `pub fn gemmHalfBatch(allocator: Allocator, io: Io, comptime kind: HalfGemm, groups: []const Group(kind.In(), kind.Out())) BatchError!void` | `gemmBatch` with `f16` or bfloat16 inputs | none | none |
| `pub fn gemmBatchStrided(allocator: Allocator, io: Io, comptime T: type, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: T, a: []const T, lda: usize, stride_a: usize, b: []const T, ldb: usize, stride_b: usize, beta: T, c: []T, ldc: usize, stride_c: usize, count: usize) Error!void` | `gemm` over `count` products at fixed strides | none | `cblas_?gemm_batch_strided` |
| `pub fn gemmHalfBatchStrided(allocator: Allocator, io: Io, comptime kind: HalfGemm, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: kind.Out(), a: []const kind.In(), lda: usize, stride_a: usize, b: []const kind.In(), ldb: usize, stride_b: usize, beta: kind.Out(), c: []kind.Out(), ldc: usize, stride_c: usize, count: usize) Error!void` | `gemmBatchStrided` with `f16` or bfloat16 inputs | none | none |

#### Extensions

| Signature | Computes | Accelerate | MKL |
| --- | --- | --- | --- |
| `pub fn geadd(io: Io, comptime T: type, m: usize, n: usize, alpha: T, a: []const T, lda: usize, beta: T, c: []T, ldc: usize, transa: Trans, transc: Trans) GeaddError!void` | `C ← α·op(A) + β·op(C)` | none | `mkl_?omatadd` |
| `pub fn omatcopy(io: Io, comptime T: type, trans: Trans, rows: usize, cols: usize, alpha: T, a: []const T, lda: usize, b: []T, ldb: usize) OmatcopyError!void` | `B ← α·op(A)`, out of place | none | `mkl_?omatcopy` |
| `pub fn imatcopy(allocator: Allocator, io: Io, comptime T: type, trans: Trans, rows: usize, cols: usize, alpha: T, a: []T, lda: usize, ldb: usize) ImatcopyError!void` | `A ← α·op(A)`, in place | none | `mkl_?imatcopy` |
| `pub fn laswp(io: Io, comptime T: type, m: usize, n: usize, a: []T, lda: usize, first: usize, ipiv: ConstStrided(usize)) LaswpError!void` | swaps rows `first` to `first + len(ipiv)` of `A` as `ipiv` says | LAPACK `?laswp` | LAPACK `?laswp` |

### `kernels.lapack`

As `src/linalg/lapack/` has it today. Accelerate (macOS 13.3's `$NEWLAPACK$ILP64`) and MKL have every routine under its LAPACK name; rstd routes only `potrf` to Accelerate so far.

#### Cholesky

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `pub fn potrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) PotrfError!void` | `A = UᴴU` or `LLᴴ`, blocked | `?potrf` (Accelerate: routed) | OpenBLAS `interface/lapack/potrf.c`; ReLAPACK with `-Drelapack` |
| `pub fn potf2(comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) Potf2Error!void` | `potrf`, unblocked | `?potf2` | OpenBLAS `potf2.c` |
| `pub fn potrs(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, nrhs: usize, a: []const T, lda: usize, b: []T, ldb: usize) PotrsError!void` | solves `AX = B` from `potrf`'s factor | `?potrs` | Reference-LAPACK `dpotrs.f` |
| `pub fn potri(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) PotriError!void` | `A⁻¹` from `potrf`'s factor | `?potri` | OpenBLAS `potri.c` |
| `pub fn pbtrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, kd: usize, ab: []T, ldab: usize) PbtrfError!void` | `potrf` for band `A` (`kd`) | `?pbtrf` | ReLAPACK's recursive `pbtrf` |

#### LU

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `pub fn getrf(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, a: []T, lda: usize, ipiv: []usize) GetrfError!?usize` | `PA = LU`, partial pivoting; the first zero `U(i, i)`, if any | `?getrf` | OpenBLAS `getrf_single.c`, `getrf_parallel.c`; ReLAPACK with `-Drelapack` |
| `pub fn getrs(allocator: Allocator, io: Io, comptime T: type, trans: Trans, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const usize, b: []T, ldb: usize) GetrsError!void` | solves `op(A)X = B` from `getrf`'s factor | `?getrs` | OpenBLAS `getrs` |
| `pub fn gesv(allocator: Allocator, io: Io, comptime T: type, n: usize, nrhs: usize, a: []T, lda: usize, ipiv: []usize, b: []T, ldb: usize) GesvError!void` | `getrf`, then `getrs` | `?gesv` | OpenBLAS `gesv` |
| `pub fn getri(allocator: Allocator, io: Io, comptime T: type, n: usize, a: []T, lda: usize, ipiv: []const usize) GetriError!void` | `A⁻¹` from `getrf`'s factor | `?getri` | Reference-LAPACK `dgetri.f` |
| `pub fn gbtrf(comptime T: type, m: usize, n: usize, kl: usize, ku: usize, ab: []T, ldab: usize, ipiv: []usize) GbtrfError!?usize` | `getrf` for band `A` (`kl`, `ku`) | `?gbtrf` | Reference-LAPACK `dgbtf2.f` |

#### Triangular

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `pub fn trtri(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, diag: Diag, n: usize, a: []T, lda: usize) TrtriError!void` | `A⁻¹`, triangular, blocked | `?trtri` | OpenBLAS `trtri.c`; ReLAPACK with `-Drelapack` |
| `pub fn trti2(comptime T: type, uplo: Uplo, diag: Diag, n: usize, a: []T, lda: usize) ArgError!void` | `trtri`, unblocked | `?trti2` | OpenBLAS `trti2.c` |
| `pub fn lauum(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) LauumError!void` | `UUᴴ` or `LᴴL`, blocked | `?lauum` | OpenBLAS `lauum.c`; ReLAPACK with `-Drelapack` |
| `pub fn lauu2(comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) ArgError!void` | `lauum`, unblocked | `?lauu2` | OpenBLAS `lauu2.c` |
| `pub fn trtrs(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: usize, nrhs: usize, a: []const T, lda: usize, b: []T, ldb: usize) TrtrsError!void` | solves `op(A)X = B`, triangular `A` | `?trtrs` | OpenBLAS `trtrs.c` |

#### Symmetric indefinite

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `pub fn sytrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize` | `A = UDUᵀ` or `LDLᵀ`, Bunch-Kaufman; the first zero pivot, if any | `?sytrf` | Reference-LAPACK `dsytrf.f` |
| `pub fn hetrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize` | `sytrf` for Hermitian `A` | `?hetrf` | Reference-LAPACK `zhetrf.f` |
| `pub fn sytrfRook(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize` | `sytrf` with rook pivoting | `?sytrf_rook` | Reference-LAPACK `dsytrf_rook.f` |
| `pub fn hetrfRook(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize` | `hetrf` with rook pivoting | `?hetrf_rook` | Reference-LAPACK `zhetrf_rook.f` |
| `pub fn sygst(allocator: Allocator, io: Io, comptime T: type, problem: Problem, uplo: Uplo, n: usize, a: []T, lda: usize, b: []const T, ldb: usize) SygstError!void` | reduces `Ax = λBx` (and the other `problem`s) to standard form, `B` from `potrf` | `?sygst` | ReLAPACK's recursive `sygst` |
| `pub fn hegst(allocator: Allocator, io: Io, comptime T: type, problem: Problem, uplo: Uplo, n: usize, a: []T, lda: usize, b: []const T, ldb: usize) SygstError!void` | `sygst` for Hermitian `A` | `?hegst` | ReLAPACK's recursive `hegst` |

#### QR and LQ

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `pub fn geqrf(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, a: []T, lda: usize, tau: []T) GeqrfError!void` | `A = QR`, `Q` as reflectors and `tau` | `?geqrf` | faer v0.24.4's staircase QR |
| `pub fn orgqr(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, k: usize, a: []T, lda: usize, tau: []const T) OrgqrError!void` | forms `Q` from `geqrf`'s reflectors | `?orgqr` | faer v0.24.4 |
| `pub fn ungqr(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, k: usize, a: []T, lda: usize, tau: []const T) OrgqrError!void` | `orgqr` for complex | `?ungqr` | faer v0.24.4 |
| `pub fn ormqr(allocator: Allocator, io: Io, comptime T: type, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmqrError!void` | `C ← op(Q)C` or `C·op(Q)` | `?ormqr` | faer v0.24.4 |
| `pub fn unmqr(allocator: Allocator, io: Io, comptime T: type, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmqrError!void` | `ormqr` for complex | `?unmqr` | faer v0.24.4 |
| `pub fn gelqf(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, a: []T, lda: usize, tau: []T) GelqfError!void` | `A = LQ` | `?gelqf` | libflame 5.2.0 `FLA_LQ_UT` |
| `pub fn orglq(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, k: usize, a: []T, lda: usize, tau: []const T) OrglqError!void` | forms `Q` from `gelqf`'s reflectors | `?orglq` | libflame 5.2.0 |
| `pub fn unglq(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, k: usize, a: []T, lda: usize, tau: []const T) OrglqError!void` | `orglq` for complex | `?unglq` | libflame 5.2.0 |
| `pub fn ormlq(allocator: Allocator, io: Io, comptime T: type, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmlqError!void` | `C ← op(Q)C` or `C·op(Q)`, LQ's `Q` | `?ormlq` | libflame 5.2.0 |
| `pub fn unmlq(allocator: Allocator, io: Io, comptime T: type, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmlqError!void` | `ormlq` for complex | `?unmlq` | libflame 5.2.0 |
| `pub fn geqp3(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, a: []T, lda: usize, jpvt: []usize, tau: []T) Geqp3Error!void` | `AP = QR`, column pivoting | `?geqp3` | faer v0.24.4's pivoted QR |

#### Least squares

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `pub fn gels(allocator: Allocator, io: Io, comptime T: type, trans: Trans, m: usize, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize) GelsError!void` | least squares or minimum norm, full rank, by QR or LQ | `?gels` | faer v0.24.4 `qr/no_pivoting/solve.rs` |
| `pub fn gelsy(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize, jpvt: []usize, rcond: Float(T)) GelsyError!usize` | least squares, any rank, pivoted QR; returns the rank | `?gelsy` | Reference-LAPACK `dgelsy.f` over faer's pivoted QR |

#### Symmetric eigenproblems

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `pub fn syevd(allocator: Allocator, io: Io, comptime T: type, jobz: Jobz, uplo: Uplo, n: usize, a: []T, lda: usize, w: []Float(T)) SyevError!void` | eigenvalues, and vectors if `jobz`, of symmetric `A`, divide and conquer | `?syevd` | Reference-LAPACK `dsyevd.f`'s arguments over `sytrd`, `stebr` |
| `pub fn heevd(allocator: Allocator, io: Io, comptime T: type, jobz: Jobz, uplo: Uplo, n: usize, a: []T, lda: usize, w: []Float(T)) SyevError!void` | `syevd` for Hermitian `A` | `?heevd` | Reference-LAPACK `zheevd.f` |
| `pub fn syev(allocator: Allocator, io: Io, comptime T: type, jobz: Jobz, uplo: Uplo, n: usize, a: []T, lda: usize, w: []Float(T)) SyevError!void` | `syevd` by QR iteration | `?syev` | Reference-LAPACK `dsyev.f` |
| `pub fn heev(allocator: Allocator, io: Io, comptime T: type, jobz: Jobz, uplo: Uplo, n: usize, a: []T, lda: usize, w: []Float(T)) SyevError!void` | `syev` for Hermitian `A` | `?heev` | Reference-LAPACK `zheev.f` |

#### Nonsymmetric eigenproblems

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `pub fn geev(allocator: Allocator, io: Io, comptime T: type, jobvl: Jobz, jobvr: Jobz, n: usize, a: []T, lda: usize, w: Eigenvalues(T), vl: []T, ldvl: usize, vr: []T, ldvr: usize, unconverged: ?*usize) GeevError!void` | eigenvalues, and left or right vectors, of general `A` | `?geev` | Reference-LAPACK `dgeev.f` over faer's Schur |
| `pub fn hseqr(io: Io, comptime T: type, job: HseqrJob, compz: Compz, n: usize, lo: usize, hi: usize, h: []T, ldh: usize, w: Eigenvalues(T), z: []T, ldz: usize, unconverged: ?*usize) HseqrError!void` | Schur form of a Hessenberg matrix | `?hseqr` | faer v0.24.4's QR iterations |
| `pub fn gees(allocator: Allocator, io: Io, comptime T: type, jobvs: Jobz, select: ?GeesSelect(T), n: usize, a: []T, lda: usize, w: Eigenvalues(T), vs: []T, ldvs: usize, unconverged: ?*usize) GeesError!usize` | Schur form `A = ZTZᴴ`, sorted by `select`; returns the count selected | `?gees` | Reference-LAPACK `dgees.f` |
| `pub fn trevc3(allocator: Allocator, io: Io, comptime T: type, side: Side, howmny: Howmny, select: ?[]bool, n: usize, t: []T, ldt: usize, vl: []T, ldvl: usize, vr: []T, ldvr: usize, mm: usize) Trevc3Error!usize` | eigenvectors of a Schur form | `?trevc3` | Reference-LAPACK `dtrevc3.f` |
| `pub fn trevc(allocator: Allocator, io: Io, comptime T: type, side: Side, howmny: Howmny, select: ?[]bool, n: usize, t: []T, ldt: usize, vl: []T, ldvl: usize, vr: []T, ldvr: usize, mm: usize) Trevc3Error!usize` | `trevc3` | `?trevc` | Reference-LAPACK `dtrevc3.f` |
| `pub fn gebal(io: Io, comptime T: type, job: BalanceJob, n: usize, a: []T, lda: usize, scale: []Float(T)) GebalError!Block` | balances `A`; returns the block `lo..hi` | `?gebal` | Reference-LAPACK `dgebal.f` |
| `pub fn gebak(io: Io, comptime T: type, job: BalanceJob, side: Side, n: usize, lo: usize, hi: usize, scale: []const Float(T), m: usize, v: []T, ldv: usize) GebakError!void` | undoes `gebal` on eigenvectors | `?gebak` | Reference-LAPACK `dgebak.f` |
| `pub fn gehrd(allocator: Allocator, io: Io, comptime T: type, n: usize, lo: usize, hi: usize, a: []T, lda: usize, tau: []T) GehrdError!void` | `A = QHQᴴ`, Hessenberg | `?gehrd` | faer v0.24.4 `evd/hessenberg.rs` |
| `pub fn orghr(allocator: Allocator, io: Io, comptime T: type, n: usize, lo: usize, hi: usize, a: []T, lda: usize, tau: []const T) OrghrError!void` | forms `Q` from `gehrd` | `?orghr` | faer v0.24.4 |
| `pub fn unghr(allocator: Allocator, io: Io, comptime T: type, n: usize, lo: usize, hi: usize, a: []T, lda: usize, tau: []const T) OrghrError!void` | `orghr` for complex | `?unghr` | faer v0.24.4 |
| `pub fn trexc(io: Io, comptime T: type, compq: Compq, n: usize, t: []T, ldt: usize, q: []T, ldq: usize, ifst: *usize, ilst: *usize) TrexcError!void` | moves one block of a Schur form | `?trexc` | Reference-LAPACK `dtrexc.f` |
| `pub fn trsen(allocator: Allocator, io: Io, comptime T: type, job: TrsenJob, compq: Compq, select: []const bool, n: usize, t_in: []T, ldt: usize, q: []T, ldq: usize, w: Eigenvalues(T)) TrsenError!Reordered(T)` | reorders a Schur form by `select`, with condition numbers | `?trsen` | Reference-LAPACK `dtrsen.f` |
| `pub fn trsyl(allocator: Allocator, io: Io, comptime T: type, trana: Trans, tranb: Trans, isgn: isize, m: usize, n: usize, a: []const T, lda: usize, b: []const T, ldb: usize, c: []T, ldc: usize) TrsylError!Solved(Float(T))` | solves `op(A)X ± X·op(B) = scale·C` | `?trsyl` | ReLAPACK `dtrsyl.c` |
| `pub fn tgsyl(allocator: Allocator, io: Io, comptime T: type, trans: Trans, ijob: usize, m: usize, n: usize, a: []const T, lda: usize, b: []const T, ldb: usize, c: []T, ldc: usize, d: []const T, ldd: usize, e: []const T, lde: usize, f: []T, ldf: usize) TgsylError!SolvedWithDif(Float(T))` | the generalized Sylvester equation, with its `dif` estimate | `?tgsyl` | ReLAPACK `dtgsyl.c` |

#### SVD

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `pub fn gesvd(allocator: Allocator, io: Io, comptime T: type, jobu: Job, jobvt: Job, m: usize, n: usize, a: []T, lda: usize, s: []Float(T), u: []T, ldu: usize, vt: []T, ldvt: usize, unconverged: ?*usize) GesvdError!void` | `A = UΣVᴴ`, by QR iteration's arguments | `?gesvd` | Reference-LAPACK `dgesvd.f` over faer's SVD |
| `pub fn gesdd(allocator: Allocator, io: Io, comptime T: type, jobz: Job, m: usize, n: usize, a: []T, lda: usize, s: []Float(T), u: []T, ldu: usize, vt: []T, ldvt: usize) GesddError!void` | `gesvd`, divide and conquer's arguments | `?gesdd` | Reference-LAPACK `dgesdd.f` over faer's SVD |

#### Auxiliary

| Signature | Computes | LAPACK name (Accelerate, MKL) | rstd's source |
| --- | --- | --- | --- |
| `pub fn lange(allocator: Allocator, io: Io, comptime T: type, norm: Norm, m: usize, n: usize, a: []const T, lda: usize) LangeError!Float(T)` | `‖A‖` by `norm`: max, one, infinity or Frobenius | `?lange` | Reference-LAPACK `dlange.f` |
| `pub fn gecon(allocator: Allocator, io: Io, comptime T: type, norm: Norm, n: usize, a: []const T, lda: usize, anorm: Float(T)) GeconError!Float(T)` | reciprocal condition number from `getrf`'s factor | `?gecon` | Reference-LAPACK `dgecon.f` |
| `pub fn laset(comptime T: type, uplo: ?Uplo, m: usize, n: usize, alpha: T, beta: T, a: []T, lda: usize) LasetError!void` | off-diagonal `alpha`, diagonal `beta` | `?laset` | Reference-LAPACK `dlaset.f` |
| `pub fn lapy2(comptime F: type, x: F, y: F) F` | `sqrt(x² + y²)` without overflow | `?lapy2` | Reference-LAPACK `dlapy2.f` |

#### Backlog ([#623](https://github.com/randyos-org/rstd/issues/623))

Each in Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/`, as OpenBLAS ships it; modified BSD), the source of its port, and in Accelerate's and MKL's LAPACK under the name shown. Complex forms take the same signature (`he` for `sy` where LAPACK names them so).

| Signature | Computes | LAPACK name (Accelerate, MKL) |
| --- | --- | --- |
| `pub fn posv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize) PosvError!void` | `potrf`, then `potrs` | `?posv` |
| `pub fn pbtrs(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, kd: usize, nrhs: usize, ab: []const T, ldab: usize, b: []T, ldb: usize) PbtrsError!void` | solves from `pbtrf`'s factor | `?pbtrs` |
| `pub fn pptrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, ap: []T) PptrfError!void` | `potrf` for packed `A` | `?pptrf` |
| `pub fn pstrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize, piv: []usize, tol: Float(T)) PstrfError!usize` | `PᵀAP = UᴴU`, complete pivoting; returns the rank | `?pstrf` |
| `pub fn pocon(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []const T, lda: usize, anorm: Float(T)) PoconError!Float(T)` | reciprocal condition number from `potrf`'s factor | `?pocon` |
| `pub fn gbtrs(allocator: Allocator, io: Io, comptime T: type, trans: Trans, n: usize, kl: usize, ku: usize, nrhs: usize, ab: []const T, ldab: usize, ipiv: []const usize, b: []T, ldb: usize) GbtrsError!void` | solves from `gbtrf`'s factor | `?gbtrs` |
| `pub fn gbsv(allocator: Allocator, io: Io, comptime T: type, n: usize, kl: usize, ku: usize, nrhs: usize, ab: []T, ldab: usize, ipiv: []usize, b: []T, ldb: usize) GbsvError!void` | `gbtrf`, then `gbtrs` | `?gbsv` |
| `pub fn gesvx(allocator: Allocator, io: Io, comptime T: type, fact: Fact, trans: Trans, n: usize, nrhs: usize, a: []T, lda: usize, af: []T, ldaf: usize, ipiv: []usize, equed: *Equed, r: []Float(T), c: []Float(T), b: []T, ldb: usize, x: []T, ldx: usize, ferr: []Float(T), berr: []Float(T)) GesvxError!Float(T)` | `gesv` with equilibration, refinement and error bounds; returns the reciprocal condition number | `?gesvx` |
| `pub fn trcon(allocator: Allocator, io: Io, comptime T: type, norm: Norm, uplo: Uplo, diag: Diag, n: usize, a: []const T, lda: usize) TrconError!Float(T)` | reciprocal condition number of triangular `A` | `?trcon` |
| `pub fn sytrs(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const Pivot, b: []T, ldb: usize) SytrsError!void` | solves from `sytrf`'s factor | `?sytrs`, `?hetrs` |
| `pub fn sysv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, nrhs: usize, a: []T, lda: usize, ipiv: []Pivot, b: []T, ldb: usize) SysvError!void` | `sytrf`, then `sytrs` | `?sysv`, `?hesv` |
| `pub fn sytri(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []const Pivot) SytriError!void` | `A⁻¹` from `sytrf`'s factor | `?sytri`, `?hetri` |
| `pub fn gelsd(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize, s: []Float(T), rcond: Float(T)) GelsdError!usize` | least squares, any rank, by SVD (divide and conquer); returns the rank | `?gelsd` |
| `pub fn gelss(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize, s: []Float(T), rcond: Float(T)) GelssError!usize` | `gelsd` by QR iteration | `?gelss` |
| `pub fn gglse(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, p: usize, a: []T, lda: usize, b: []T, ldb: usize, c: []T, d: []T, x: []T) GglseError!void` | least squares with equality constraints `Bx = d` | `?gglse` |
| `pub fn ggglm(allocator: Allocator, io: Io, comptime T: type, n: usize, m: usize, p: usize, a: []T, lda: usize, b: []T, ldb: usize, d: []T, x: []T, y: []T) GgglmError!void` | general Gauss-Markov linear model | `?ggglm` |
| `pub fn syevr(allocator: Allocator, io: Io, comptime T: type, jobz: Jobz, range: Range(Float(T)), uplo: Uplo, n: usize, a: []T, lda: usize, abstol: Float(T), w: []Float(T), z: []T, ldz: usize, isuppz: []usize) SyevrError!usize` | selected eigenvalues and vectors, relatively robust representations; returns the count | `?syevr`, `?heevr` |
| `pub fn syevx(allocator: Allocator, io: Io, comptime T: type, jobz: Jobz, range: Range(Float(T)), uplo: Uplo, n: usize, a: []T, lda: usize, abstol: Float(T), w: []Float(T), z: []T, ldz: usize, ifail: []usize) SyevxError!usize` | `syevr` by bisection and inverse iteration | `?syevx`, `?heevx` |
| `pub fn sygvd(allocator: Allocator, io: Io, comptime T: type, problem: Problem, jobz: Jobz, uplo: Uplo, n: usize, a: []T, lda: usize, b: []T, ldb: usize, w: []Float(T)) SygvdError!void` | generalized symmetric-definite eigenproblem | `?sygvd`, `?hegvd` |
| `pub fn geevx(allocator: Allocator, io: Io, comptime T: type, balance: BalanceJob, jobvl: Jobz, jobvr: Jobz, n: usize, a: []T, lda: usize, w: Eigenvalues(T), vl: []T, ldvl: usize, vr: []T, ldvr: usize, scale: []Float(T), rconde: []Float(T), rcondv: []Float(T)) GeevxError!BalanceBlock` | `geev` with balancing and condition numbers | `?geevx` |
| `pub fn ggev(allocator: Allocator, io: Io, comptime T: type, jobvl: Jobz, jobvr: Jobz, n: usize, a: []T, lda: usize, b: []T, ldb: usize, alpha: Eigenvalues(T), beta: []T, vl: []T, ldvl: usize, vr: []T, ldvr: usize) GgevError!void` | generalized eigenproblem `Ax = λBx`, `λ = α/β` | `?ggev3` |
| `pub fn gges(allocator: Allocator, io: Io, comptime T: type, jobvsl: Jobz, jobvsr: Jobz, select: ?GgesSelect(T), n: usize, a: []T, lda: usize, b: []T, ldb: usize, alpha: Eigenvalues(T), beta: []T, vsl: []T, ldvsl: usize, vsr: []T, ldvsr: usize) GgesError!usize` | generalized Schur form, sorted; returns the count selected | `?gges3` |
| `pub fn gesvdx(allocator: Allocator, io: Io, comptime T: type, jobu: Jobz, jobvt: Jobz, range: Range(Float(T)), m: usize, n: usize, a: []T, lda: usize, s: []Float(T), u: []T, ldu: usize, vt: []T, ldvt: usize) GesvdxError!usize` | selected singular values and vectors; returns the count | `?gesvdx` |
| `pub fn gejsv(allocator: Allocator, io: Io, comptime T: type, jobu: Job, jobv: Job, m: usize, n: usize, a: []T, lda: usize, s: []Float(T), u: []T, ldu: usize, v: []T, ldv: usize) GejsvError!void` | SVD by preconditioned Jacobi, high relative accuracy | `?gejsv` |
| `pub fn ggsvd3(allocator: Allocator, io: Io, comptime T: type, jobu: Jobz, jobv: Jobz, jobq: Jobz, m: usize, n: usize, p: usize, a: []T, lda: usize, b: []T, ldb: usize, alpha: []Float(T), beta: []Float(T), u: []T, ldu: usize, v: []T, ldv: usize, q: []T, ldq: usize) Ggsvd3Error!struct { k: usize, l: usize }` | generalized SVD of `(A, B)` | `?ggsvd3` |
| `pub fn lansy(allocator: Allocator, io: Io, comptime T: type, norm: Norm, uplo: Uplo, n: usize, a: []const T, lda: usize) LansyError!Float(T)` | `lange` for symmetric `A` | `?lansy`, `?lanhe` |
| `pub fn lantr(allocator: Allocator, io: Io, comptime T: type, norm: Norm, uplo: Uplo, diag: Diag, m: usize, n: usize, a: []const T, lda: usize) LantrError!Float(T)` | `lange` for trapezoidal `A` | `?lantr` |

`GgesSelect(T)` is `GeesSelect`'s form for a pair: `*const fn (T, T) bool` for complex, `*const fn (Float(T), Float(T), Float(T)) bool` for real (`α_re`, `α_im`, `β`).

### `kernels.sparse`

Today the routines are methods of the formats, and stay so. The iterative solvers' `Options(T)` struct becomes their trailing arguments in the move ([#610](https://github.com/randyos-org/rstd/issues/610)).

| Signature | Computes | Accelerate | MKL | Source | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `Csc(T)`: `pub fn alloc(allocator: Allocator, rows: u32, cols: u32, nnz: usize) Allocator.Error!Self`, `pub fn deinit(self: *Self, allocator: Allocator) void`, `pub fn nonzeros(self: Self) usize`, `pub fn col(self: Self, j: usize) Vector(T)` | storage | `SparseMatrix_Double` | `mkl_sparse_?_create_csc` | rstd | yes |
| `Csc(T)`: `pub fn fromDense(allocator: Allocator, a: []const T, lda: usize, rows: u32, cols: u32, drop_below: Float(T)) Error!Self`, `pub fn toDense(self: Self, allocator: Allocator, a: []T, lda: usize) Allocator.Error!void` | conversion from and to dense, dropping `\|aᵢⱼ\| < drop_below` | `SparseConvertFromCoordinate` | `mkl_sparse_convert_*` | rstd | yes; takes a slice and `lda` in place of today's views (#610) |
| `Csc(T)`: `pub fn toCsr(self: Self, allocator: Allocator) Allocator.Error!Csr(T)`, `pub fn asCsrOfTranspose(self: Self) Csr(T)`; `Csr(T)`: `pub fn toCsc(self: Self, allocator: Allocator) Allocator.Error!Csc(T)`, `pub fn asCscOfTranspose(self: Self) Csc(T)`, `pub fn row(self: Self, i: usize) Vector(T)` | format changes; the `as` forms reuse the arrays | none | `mkl_sparse_convert_csr` | rstd | yes |
| `Csc(T)`: `pub fn transpose(self: Self, allocator: Allocator) Allocator.Error!Self`, `pub fn ctranspose(self: Self, allocator: Allocator) Allocator.Error!Self` | `Aᵀ`, `Aᴴ` | `SparseGetTranspose` | `mkl_sparse_convert_csr` with `SPARSE_OPERATION_TRANSPOSE` | rstd | yes |
| `Csc(T)`: `pub fn scal(self: Self, alpha: T) void`, `pub fn sub(self: Self, allocator: Allocator, i: usize, j: usize, rows: usize, cols: usize) SubError!Self`, `pub fn selectCols(self: Self, allocator: Allocator, cols: []const u32) SelectError!Self` | scaling, a block, chosen columns | none | none | rstd | yes |
| `Csc(T)`, `Csr(T)`: `pub fn mulVec(self: Self, allocator: Allocator, io: Io, trans: Trans, alpha: T, x: []const T, beta: T, y: []T) MulError!void` | `y ← α·op(A)x + βy` | `SparseMultiplyAdd` | `mkl_sparse_?_mv` | rstd | yes |
| `Csc(T)`: `pub fn mulDense(self: Self, allocator: Allocator, io: Io, b: []const T, ldb: usize, cols: usize, c: []T, ldc: usize) MulError!void` | `C ← AB`, dense `B` | `SparseMultiply` | `mkl_sparse_?_mm` | rstd | yes; writes into `c` in place of returning a matrix (#610) |
| `Csc(T)`: `pub fn spgemmSymbolic(self: Self, allocator: Allocator, b: Self) SpgemmError!Self`, `pub fn spgemmNumeric(self: Self, allocator: Allocator, b: Self, c: Self) NumericError!void`, `pub fn spgemm(self: Self, allocator: Allocator, b: Self) SpgemmError!Self` | `C ← AB`, sparse `B`: the pattern, the values, both | none | `mkl_sparse_sp2m` | rstd (Gustavson) | yes |
| `Csc(T)`: `pub fn solveTri(self: Self, uplo: Uplo, trans: Trans, diag: Diag, x: []T) TriSolveError!void` | `x ← op(A)⁻¹x`, triangular `A` | `SparseSolve` (triangular) | `mkl_sparse_?_trsv` | rstd | yes |
| `Triplets(T)`: `pub fn init(rows: usize, cols: usize) error{TooLarge}!Self`, `pub fn deinit(self: *Self, allocator: Allocator) void`, `pub fn append(self: *Self, allocator: Allocator, row: usize, col: usize, value: T) (Allocator.Error \|\| error{IndexOutOfRange})!void`, `pub fn toCsc(self: Self, allocator: Allocator) Error!Csc(T)` | assembly, duplicates summed | `SparseConvertFromCoordinate` | `mkl_sparse_?_create_coo` | rstd | yes |
| `Vector(T)`: `pub fn axpy(x: Self, alpha: T, y: []T) VectorError!void`, `pub fn dot(x: Self, y: []const T) VectorError!T`, `pub fn dotc(x: Self, y: []const T) VectorError!T` | sparse level 1 | none | `cblas_?axpyi`, `cblas_?doti`, `cblas_?dotci` | rstd | yes |
| `pub fn cg(allocator: Allocator, io: Io, comptime T: type, a: Operator(T), preconditioner: ?Operator(T), tol: Float(T), max_iterations: usize, b: []const T, x: []T) SolveError!Stats` | conjugate gradients, symmetric positive definite `A` | `SparseConjugateGradient` | RCI ISS `dcg` | rstd | yes; arguments in place of `Options(T)` (#610) |
| `pub fn gmres(allocator: Allocator, io: Io, comptime T: type, a: Operator(T), preconditioner: ?Operator(T), tol: Float(T), max_iterations: usize, restart: usize, b: []const T, x: []T) SolveError!Stats` | restarted GMRES | `SparseGMRES` | RCI ISS `dfgmres` | rstd | yes; as `cg` (#610) |
| `Csc(T)`: `pub fn add(self: Self, allocator: Allocator, alpha: T, b: Self, beta: T) AddError!Self` | `αA + βB` | none | `mkl_sparse_?_add` | SciPy `csr_plus_csr`; rstd's `spgemm` merge | backlog |
| `SparseCholesky(T)`: `pub fn analyze(allocator: Allocator, a: Csc(T), uplo: Uplo) AnalyzeError!Self`, `pub fn factor(self: *Self, allocator: Allocator, io: Io, a: Csc(T)) FactorError!void`, `pub fn solve(self: Self, io: Io, b: []T) SolveError!void`, `pub fn deinit(self: *Self, allocator: Allocator) void` | sparse `A = LLᴴ`: AMD ordering and elimination tree, then values, then solves | `SparseFactor(SparseFactorizationCholesky)` | PARDISO (`mtype` 2, 4) | faer v0.24.4 `sparse/linalg/cholesky.rs`, `amd.rs` (MIT) | backlog |
| `SparseLdlt(T)`: as `SparseCholesky(T)` | sparse `A = LDLᴴ` | `SparseFactorizationLDLT` | PARDISO (`mtype` -2) | faer v0.24.4 (MIT) | backlog |
| `SparseLu(T)`: as `SparseCholesky(T)`, `analyze` without `uplo` | sparse `PAQ = LU` | none | PARDISO (`mtype` 11, 13) | faer v0.24.4 `sparse/linalg/lu.rs` (MIT) | backlog |
| `SparseQr(T)`: as `SparseLu(T)`, `solve` minimizing `‖Ax − b‖₂` | sparse `AP = QR` | `SparseFactorizationQR` | `mkl_sparse_?_qr` | faer v0.24.4 `sparse/linalg/qr.rs` (MIT) | backlog |
| `pub fn lsmr(allocator: Allocator, io: Io, comptime T: type, a: Operator(T), at: Operator(T), damp: Float(T), atol: Float(T), btol: Float(T), max_iterations: usize, b: []const T, x: []T) SolveError!Stats` | least squares by LSMR | `SparseLSMR` | none | SciPy `sparse/linalg/_isolve/lsmr.py` (BSD-3) | backlog |
| `Jacobi(T)`: `pub fn init(allocator: Allocator, a: Csc(T)) JacobiError!Self`, `pub fn operator(self: *const Self) Operator(T)`, `pub fn deinit(self: *Self, allocator: Allocator) void` | diagonal preconditioner | `SparsePreconditionerDiagonal` | none | Saad, *Iterative Methods for Sparse Linear Systems*, 10.2 | backlog |
| `Ilu0(T)`: as `Jacobi(T)` | incomplete LU without fill | none | `dcsrilu0` | Saad, algorithm 10.4 | backlog |

### `kernels.fft`

In the move the transforms take the norm as a trailing argument in place of `Options` ([#611](https://github.com/randyos-org/rstd/issues/611)); `rstd.math.fft`'s `*With` siblings take the struct.

| Signature | Computes | Accelerate | MKL | Source | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `Plan(T)`: `pub fn init(allocator: Allocator, n: usize) InitError!Self`, `pub fn deinit(self: *Self, allocator: Allocator) void`, `pub fn scratchLen(self: Self) usize` | plans a complex transform of length `n` | `vDSP_DFT_zop_CreateSetup` | `DftiCreateDescriptor`, `DftiCommitDescriptor` | PocketFFT `cfftp`, Bluestein; pffft for `f32` | yes |
| `Plan(T)`: `pub fn transform(self: *const Self, io: Io, data: []Complex(T), scratch: []Complex(T), direction: Direction, layout: Layout, norm: Norm) TransformError!void` | in place, a batch by `layout` | `vDSP_DFT_Execute`, `vDSP_fft_zip` | `DftiComputeForward`, `DftiComputeBackward` | as above | yes (#611: `norm`) |
| `Plan(T)`: `pub fn transformInto(self: *const Self, io: Io, out: []Complex(T), in: []const Complex(T), scratch: []Complex(T), direction: Direction, layouts: Layouts, norm: Norm) TransformIntoError!void` | out of place | `vDSP_DFT_Execute` | `DftiComputeForward` | as above | yes (#611) |
| `RealPlan(T)`: `pub fn init(allocator: Allocator, n: usize) InitError!Self`, `pub fn deinit(self: *Self, allocator: Allocator) void`, `pub fn scratchLen(self: Self) usize` | plans a real transform | `vDSP_DFT_zrop_CreateSetup` | `DftiCreateDescriptor` (real) | PocketFFT `rfftp`; pffft | yes |
| `RealPlan(T)`: `pub fn forward(self: *const Self, io: Io, out: []Complex(T), in: []const T, scratch: []T, layouts: Layouts, norm: Norm) TransformIntoError!void`, `pub fn backward(self: *const Self, io: Io, out: []T, in: []const Complex(T), scratch: []T, layouts: Layouts, norm: Norm) TransformIntoError!void` | real to half-spectrum and back | `vDSP_DFT_Execute`, `vDSP_fft_zrip` | `DftiComputeForward` | as above | yes (#611) |
| `RealPlan(T)`: `pub fn forwardUnordered(self: *const Self, out: []T, in: []const T, scratch: []T) UnorderedError!void`, `pub fn backwardUnordered(self: *const Self, out: []T, in: []const T, scratch: []T) UnorderedError!void`, `pub fn convolveAccumulate(self: *const Self, acc: []T, a: []const T, b: []const T, scale: T) UnorderedError!void` | the spectrum in the transform's own order, and `acc += scale·a⊙b` in it: fast convolution | none | none | pffft `pffft_zconvolve_accumulate` | yes |
| `StaticPlan(n, T)`: `pub fn scratchLen(_: Self) usize`, `pub fn transform(self: Self, io: Io, data: []Complex(T), scratch: []Complex(T), direction: Direction, layout: Layout, norm: Norm) TransformError!void`, `pub fn transformInto(self: Self, io: Io, out: []Complex(T), in: []const Complex(T), scratch: []Complex(T), direction: Direction, layouts: Layouts, norm: Norm) TransformIntoError!void` | `Plan` with `n` at compile time | none | none | PocketFFT, unrolled | stub (#625) |
| `StaticRealPlan(n, T)`: `pub fn scratchLen(_: Self) usize`, `pub fn forward(self: Self, io: Io, out: []Complex(T), in: []const T, scratch: []T, layouts: Layouts, norm: Norm) TransformIntoError!void`, `pub fn backward(self: Self, io: Io, out: []T, in: []const Complex(T), scratch: []T, layouts: Layouts, norm: Norm) TransformIntoError!void` | `RealPlan` with `n` at compile time | none | none | PocketFFT, unrolled | stub (#625) |
| `TrigPlan(T)`: `pub fn init(allocator: Allocator, n: usize, kind: TrigKind) InitError!Self`, `pub fn deinit(self: *Self, allocator: Allocator) void`, `pub fn scratchLen(self: Self) usize`, `pub fn transform(self: *const Self, io: Io, data: []T, scratch: []T, layout: Layout, norm: Norm) TransformError!void` | DCT and DST types I to IV | `vDSP_DCT_CreateSetup` (II, III, IV) | trigonometric transforms (`d_init_trig_transform`) | PocketFFT C++ `pocketfft_hdronly.hpp` `T_dcst4`, `T_dct1` (BSD-3) | stub (#625) |
| `Mdct(T)`: `pub fn init(allocator: Allocator, n: usize) InitError!Self`, `pub fn deinit(self: *Self, allocator: Allocator) void`, `pub fn scratchLen(self: Self) usize`, `pub fn forward(self: *const Self, out: []T, in: []const T, scratch: []T) FrameError!void`, `pub fn backward(self: *const Self, out: []T, in: []const T, scratch: []T) FrameError!void` | modified DCT of `n` samples to `n/2` and back, by a DCT-IV | none | none | Princen-Bradley via `TrigPlan(.dct4)` | stub (#625) |
| `PlanNd(T)`: `pub fn init(allocator: Allocator, shape: []const usize) InitError!Self`, `pub fn deinit(self: *Self, allocator: Allocator) void`, `pub fn scratchLen(self: Self) usize`, `pub fn transform(self: *const Self, io: Io, data: []Complex(T), scratch: []Complex(T), direction: Direction, norm: Norm) TransformError!void` | complex transform over every axis, column-major | `vDSP_fft2d_zip` (powers of two) | DFTI, multi-dimensional | PocketFFT C++ `c2c` over axes (BSD-3) | backlog |
| `pub fn fftfreq(comptime T: type, n: usize, d: T, out: []T) LengthError!void`, `pub fn rfftfreq(comptime T: type, n: usize, d: T, out: []T) LengthError!void` | the frequencies of the bins | none | none | NumPy `numpy.fft.fftfreq`, `rfftfreq` | backlog |
| `pub fn fftshift(comptime T: type, x: []T) void`, `pub fn ifftshift(comptime T: type, x: []T) void` | moves the zero frequency to the middle and back | none | none | NumPy `numpy.fft.fftshift` | backlog |

### `kernels.elementwise`

`pub const Error = LengthError || Io.Cancelable;` for all. One input: `y[i] = f(x[i])`, `y` may be `x`; two: `z[i] = f(x[i], y[i])`, `z` may be either. Each states its accuracy: the v0.1.0 set ([#614](https://github.com/randyos-org/rstd/issues/614)) is Zig's vector builtins, the accuracy of Zig's `std.math` (musl's algorithms), or exact.

| Signature | Computes | Accelerate | MKL | Source | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `pub fn add(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void` | `x + y`, exact rounding | `vDSP_vadd`, `vDSP_zvadd` | `v?Add` | `@Vector` | #614 |
| `pub fn sub(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void` | `x − y` | `vDSP_vsub` (operands swapped: vDSP subtracts its first from its second) | `v?Sub` | `@Vector` | #614 |
| `pub fn mul(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void` | `x · y` | `vDSP_vmul`, `vDSP_zvmul` | `v?Mul` | `@Vector` | #614 |
| `pub fn div(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void` | `x / y` | `vvdiv`, `vDSP_zvdiv` | `v?Div` | `@Vector` | #614 |
| `pub fn sqrt(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `√x`, correctly rounded | `vvsqrt` | `v?Sqrt` | `@sqrt` | #614 |
| `pub fn exp(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `eˣ` | `vvexp` | `v?Exp` | `@exp` | #614 |
| `pub fn exp2(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `2ˣ` | `vvexp2` | `v?Exp2` | `@exp2` | #614 |
| `pub fn log(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `ln x` | `vvlog` | `v?Ln` | `@log` | #614 |
| `pub fn log2(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `log₂ x` | `vvlog2` | `v?Log2` | `@log2` | #614 |
| `pub fn log10(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `log₁₀ x` | `vvlog10` | `v?Log10` | `@log10` | #614 |
| `pub fn sin(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `sin x` | `vvsin` | `v?Sin` | `@sin` | #614 |
| `pub fn cos(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `cos x` | `vvcos` | `v?Cos` | `@cos` | #614 |
| `pub fn tan(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `tan x` | `vvtan` | `v?Tan` | `@tan` | #614 |
| `pub fn abs(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(Float(T))) Error!void` | `\|x\|`; complex: the modulus, without overflow | `vvfabs`, `vDSP_zvabs` | `v?Abs` | `@abs`; complex by `std.math.hypot` | #614 |
| `pub fn floor(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `⌊x⌋` | `vvfloor` | `v?Floor` | `@floor` | #614 |
| `pub fn ceil(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `⌈x⌉` | `vvceil` | `v?Ceil` | `@ceil` | #614 |
| `pub fn trunc(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | toward zero | `vvint` | `v?Trunc` | `@trunc` | #614 |
| `pub fn round(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | to nearest, ties away from zero (Fortran's `ANINT`) | none (`vvnint` ties to even) | `v?Round` | `@round` | #614 |
| `pub fn rsqrt(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `1/√x` | `vvrsqrt` | `v?InvSqrt` | `1 / @sqrt` | backlog |
| `pub fn recip(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `1/x` | `vvrec` | `v?Inv` | `@Vector` | backlog |
| `pub fn cbrt(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `∛x` | none | `v?Cbrt` | `std.math.cbrt`; Julia `Base.Math.cbrt` (MIT) | backlog |
| `pub fn hypot(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void` | `√(x² + y²)`, without overflow | none | `v?Hypot` | Julia `Base.Math.hypot` (MIT) | backlog |
| `pub fn expm1(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `eˣ − 1`, accurate near 0 | `vvexpm1` | `v?Expm1` | `std.math.expm1`; Julia `Base.Math.expm1` | backlog |
| `pub fn log1p(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `ln(1 + x)`, accurate near 0 | `vvlog1p` | `v?Log1p` | `std.math.log1p` | backlog |
| `pub fn pow(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void` | `xʸ` | `vvpow` | `v?Pow` | `std.math.pow`; Julia `Base.Math.pow` | backlog |
| `pub fn powx(io: Io, comptime T: type, b: T, x: ConstStrided(T), y: Strided(T)) Error!void` | `xᵇ`, one exponent | `vvpows` | `v?Powx` | `std.math.pow` | backlog |
| `pub fn sincos(io: Io, comptime T: type, x: ConstStrided(T), s: Strided(T), c: Strided(T)) Error!void` | `sin x` and `cos x` in one reduction | `vvsincos` | `v?SinCos` | Julia `Base.Math.sincos` | backlog |
| `pub fn asin(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `arcsin x` | `vvasin` | `v?Asin` | `std.math.asin` | backlog |
| `pub fn acos(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `arccos x` | `vvacos` | `v?Acos` | `std.math.acos` | backlog |
| `pub fn atan(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `arctan x` | `vvatan` | `v?Atan` | `std.math.atan` | backlog |
| `pub fn atan2(io: Io, comptime T: type, y: ConstStrided(T), x: ConstStrided(T), z: Strided(T)) Error!void` | the angle of `(x, y)`, Fortran's `ATAN2(Y, X)` order | `vvatan2` | `v?Atan2` | `std.math.atan2` | backlog |
| `pub fn sinh(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `sinh x` | `vvsinh` | `v?Sinh` | `std.math.sinh` | backlog |
| `pub fn cosh(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `cosh x` | `vvcosh` | `v?Cosh` | `std.math.cosh` | backlog |
| `pub fn tanh(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `tanh x` | `vvtanh` | `v?Tanh` | `std.math.tanh` | backlog |
| `pub fn asinh(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `arsinh x` | `vvasinh` | `v?Asinh` | `std.math.asinh` | backlog |
| `pub fn acosh(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `arcosh x` | `vvacosh` | `v?Acosh` | `std.math.acosh` | backlog |
| `pub fn atanh(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `artanh x` | `vvatanh` | `v?Atanh` | `std.math.atanh` | backlog |
| `pub fn erf(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | the error function | none | `v?Erf` | openlibm `s_erf.c` (BSD), as SpecialFunctions.jl calls it; num.zig `erf` | backlog |
| `pub fn erfc(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `1 − erf x`, accurate for large `x` | none | `v?Erfc` | openlibm `s_erf.c` | backlog |
| `pub fn lgamma(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `ln \|Γ(x)\|` | none | `v?LGamma` | openlibm `e_lgamma_r.c`; SpecialFunctions.jl (MIT) | backlog |
| `pub fn tgamma(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `Γ(x)` | none | `v?TGamma` | openlibm `e_tgamma.c` (Fortran's `GAMMA`) | backlog |
| `pub fn arg(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(Float(T))) Error!void` | complex `x`'s angle | `vDSP_zvphas` | `v?Arg` | `std.math.atan2(im, re)` | backlog |
| `pub fn conj(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | complex conjugate | `vDSP_zvconj` | `v?Conj` | `@Vector` | backlog |

### `kernels.dsp`

`Error` as in `elementwise`.

| Signature | Computes | Accelerate | MKL | Source | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `pub fn window(comptime T: type, out: []T, kind: Window, symmetry: Symmetry) WindowError!void` | a window of `out.len` points | `vDSP_hann_window`, `vDSP_hamm_window`, `vDSP_blkman_window` | none | MATLAB `hann`, `hamming`, `blackman`, `kaiser` | yes (#612) |
| `pub fn sadd(io: Io, comptime T: type, alpha: T, x: ConstStrided(T), y: Strided(T)) Error!void` | `x + α` | `vDSP_vsadd` | none | `@Vector` | backlog |
| `pub fn sdiv(io: Io, comptime T: type, alpha: T, x: ConstStrided(T), y: Strided(T)) Error!void` | `x / α` | `vDSP_vsdiv` | none | `@Vector` | backlog |
| `pub fn ma(io: Io, comptime T: type, a: ConstStrided(T), b: ConstStrided(T), c: ConstStrided(T), d: Strided(T)) Error!void` | `a·b + c`, fused where the core has an instruction | `vDSP_vma` | `v?Fma` | `@mulAdd` | backlog |
| `pub fn fill(io: Io, comptime T: type, alpha: T, y: Strided(T)) Io.Cancelable!void` | every element `α` | `vDSP_vfill` | none | NumPy `full` | backlog |
| `pub fn ramp(io: Io, comptime T: type, start: T, step: T, y: Strided(T)) Io.Cancelable!void` | `start + i·step`, each computed, not accumulated | `vDSP_vramp` | none | NumPy `arange` | backlog |
| `pub fn linspace(io: Io, comptime T: type, start: T, stop: T, y: Strided(T)) Io.Cancelable!void` | `y.n` points from `start` to `stop`, both included | `vDSP_vgen` | none | NumPy `linspace` | backlog |
| `pub fn neg(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `−x` | `vDSP_vneg` | none | `@Vector` | backlog |
| `pub fn clip(io: Io, comptime T: type, lo: T, hi: T, x: ConstStrided(T), y: Strided(T)) Error!void` | `x` held to `[lo, hi]`; NaN stays NaN | `vDSP_vclip` | none | NumPy `clip` | backlog |
| `pub fn reverse(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) Error!void` | `y[i] = x[n − 1 − i]` | `vDSP_vrvrs` | none | NumPy `flip` | backlog |
| `pub fn convolve(io: Io, comptime T: type, mode: ConvMode, x: ConstStrided(T), h: ConstStrided(T), y: Strided(T)) Error!void` | `y = x ∗ h`, `y.n` by `mode`: `n + m − 1`, `max(n, m)`, `max(n, m) − min(n, m) + 1` | `vDSP_conv` (`h` reversed) | `vslsConvExec` | NumPy `convolve`; direct, `RealPlan` products past a size | backlog |
| `pub fn correlate(io: Io, comptime T: type, mode: ConvMode, x: ConstStrided(T), h: ConstStrided(T), y: Strided(T)) Error!void` | cross-correlation, NumPy's conjugation of `h` | `vDSP_conv` | `vslsCorrExec` | NumPy `correlate` | backlog |
| `pub fn biquad(comptime T: type, sections: []const Biquad(T), state: []T, x: []const T, y: []T) Error!void` | cascaded second-order sections, direct form II transposed; `state` two per section, kept across calls | `vDSP_biquad` | none | SciPy `sosfilt` | backlog |
| `pub fn decimate(io: Io, comptime T: type, factor: usize, filter: []const T, x: ConstStrided(T), y: Strided(T)) Error!void` | FIR filter, then every `factor`th sample | `vDSP_desamp` | none | SciPy `upfirdn` | backlog |
| `pub fn cumsum(comptime T: type, x: []const T, y: []T) LengthError!void` | running sums, in order | `vDSP_vrsum` | none | NumPy `cumsum`, MATLAB `cumsum` | backlog |
| `pub fn polyval(io: Io, comptime T: type, coeffs: []const T, x: ConstStrided(T), y: Strided(T)) Error!void` | polynomial at each `x`, coefficients highest first, by Horner | `vDSP_vpoly` | none | NumPy `polyval`; num.zig `poly.val` | backlog |
| `pub fn interp(io: Io, comptime T: type, xp: []const T, fp: []const T, x: ConstStrided(T), y: Strided(T)) Error!void` | linear interpolation of `(xp, fp)` at `x`, ends held | none | `dfdInterpolate1D` | NumPy `interp` | backlog |
| `pub fn abs2(io: Io, comptime T: type, x: ConstStrided(Complex(T)), y: Strided(T)) Error!void` | `\|x\|²`, a power spectrum | `vDSP_zvmags` | none | Julia `abs2` | backlog |
| `pub fn toDecibels(io: Io, comptime T: type, scale: DbScale, reference: T, x: ConstStrided(T), y: Strided(T)) Error!void` | `10·log10(x/ref)` or `20·log10(x/ref)` | `vDSP_vdbcon` | none | librosa `power_to_db` | backlog |
| `pub fn deinterleave(io: Io, comptime T: type, z: ConstStrided(Complex(T)), re: Strided(T), im: Strided(T)) Error!void`, `pub fn interleave(io: Io, comptime T: type, re: ConstStrided(T), im: ConstStrided(T), z: Strided(Complex(T))) Error!void` | interleaved complex to split and back | `vDSP_ctoz`, `vDSP_ztoc` | none | `@Vector` shuffles | backlog |

### `kernels.stats`

| Signature | Computes | Accelerate | MKL | Source | v0.1.0 |
| --- | --- | --- | --- | --- | --- |
| `kernels.blas.sum`, `sumWith`, `index`, `extreme` | sums and extremes: the BLAS's, not copies | `vDSP_sve`, `vDSP_maxv`, `vDSP_minv` | VS | rstd | yes |
| `pub fn deviations(io: Io, comptime T: type, s: Summation, m: T, x: ConstStrided(T)) Io.Cancelable!Deviations(T)` | `Σ(x − m)` and `Σ\|x − m\|²` in one pass; the same bits at any worker count | none | `vsldSSCompute` (`VSL_SS_2C_MOM`) | Chan, Golub and LeVeque's corrected two-pass; rstd's `fast.zig` lanes | #605 |
| `pub fn centralMoments(io: Io, comptime T: type, s: Summation, m: T, x: ConstStrided(T)) Io.Cancelable!CentralMoments(T)` | `Σ(x − m)²`, `³`, `⁴`, for skewness and kurtosis | none | `vsldSSCompute` (`VSL_SS_4C_MOM`) | SciPy `stats.moment`; StatsBase `moment` | backlog |
| `Running(F)`: `pub fn push(self: *Self, x: F) void`, `pub fn mean(self: Self) error{EmptyInput}!F`, `pub fn variance(self: Self, normalize: Normalize) error{EmptyInput}!F` | Welford's running mean and variance | none | VS streaming | Welford 1962 | yes (#613) |
| `Running(F)`: `pub fn merge(self: *Self, other: Self) void` | combines two runs (Chan's formula), for threads | none | none | Chan, Golub and LeVeque 1979 | backlog |
| `pub fn histogramInto(comptime T: type, edges: []const T, x: ConstStrided(T), counts: []u64) HistogramIntoError!void` | counts into the bins `edges` bound, NumPy's half-open bins, the last closed | `vImageHistogramCalculation` in part | none | NumPy `histogram` | yes (#613: `ConstStrided` in place of today's `anytype`, arguments in the kernel order) |
| `pub fn select(comptime T: type, x: []T, k: usize) T` | the `k`th smallest, `x` partitioned around it, in O(n); NaNs last | none | `vsldSSCompute` (`VSL_SS_QUANTS`) | NumPy `partition` (introselect); Rust `select_nth_unstable` (MIT) | backlog; `quantile` sorts a copy today |
| `pub fn medianAbsDeviation(allocator: Allocator, comptime T: type, x: []const T) MadError!T` | `median(\|x − median(x)\|)` | none | `vsldSSCompute` (`VSL_SS_MDAD`) | SciPy `stats.median_abs_deviation` | backlog |
| `Rng`: `pub fn init(seed: u64) Rng`, `pub fn jump(self: *Rng, chunks: u64) void` | a stream and its chunk's stream | none | `vslNewStream`, `vslSkipAheadStream` | Zig `std.Random.Xoshiro256` (MIT) | backlog |
| `pub fn uniform(io: Io, comptime T: type, rng: *Rng, lo: T, hi: T, y: Strided(T)) Io.Cancelable!void` | uniform on `[lo, hi)` | none | `v?RngUniform` | NumPy `Generator.uniform`; num.zig `random.uniform` | backlog |
| `pub fn normal(io: Io, comptime T: type, rng: *Rng, mean: T, sd: T, y: Strided(T)) Io.Cancelable!void` | normal, by the ziggurat | none | `v?RngGaussian` | Zig `std.Random.floatNorm` (Marsaglia and Tsang) | backlog |

### `rstd.math` (scalar)

| Signature | Computes | Source | v0.1.0 |
| --- | --- | --- | --- |
| `pub fn floor(x: anytype) @TypeOf(x)` | `⌊x⌋`; `@floor`, but exact for `f64` on the x87 (`x87_f64`) | rstd | yes (#608) |
| `pub fn ceil(x: anytype) @TypeOf(x)` | `⌈x⌉`, as `floor` | rstd | yes (#608) |
| `pub fn round(x: anytype) @TypeOf(x)` | to nearest, ties away from zero, as `floor` | rstd | yes (#608) |
| `pub fn rationalApprox(x: anytype, tol: ?@TypeOf(x)) Fraction(@TypeOf(x))` | the continued fraction of `x` within `tol` (MATLAB `rat`) | MATLAB `rat` | yes (#608) |
| `pub const x87_f64: bool` | whether `f64` math runs on the x87 (32-bit x86 without SSE2), where compiler-rt's `f64` `@floor`, `@ceil` and `@round` keep fraction bits ([#264](https://github.com/randyos-org/rstd/issues/264)) | rstd | yes (#608) |

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
