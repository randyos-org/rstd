# rstd design

This document specifies rstd's modules beyond linear algebra's own (that is in
[`src/linalg/design.md`](src/linalg/design.md)): statistics, sorting and sets, the FFT, signal
processing, and the N-d arrays and sparse matrices that extend `rstd.linalg`. It describes the
code as it is, and says "Planned" where it is not. The work left is in [`todo.md`](todo.md).

## Requirements

Shared with the public `rstd.linalg` layer ([`design.md`](src/linalg/design.md), "Operations over
views"):

- `(allocator, io, ...)` where a function allocates, `(io, ...)` where it only calls the BLAS;
  allocator first. A slice's element type is an explicit `comptime` parameter after them, as in
  the BLAS (`stats.mean(io, f64, xs)`).
- Types are unmanaged and store no allocator. A new result comes back owned, and each such
  function has an `...Into` form that writes caller storage.
- Error sets are exact and named, with plain error names. Operands are slices, not pointers and
  lengths. Storage is column-major.
- Results are bit-identical for any thread count, as the BLAS's are
  ([`threading.md`](src/linalg/blas/threading.md)).
- Zig-native: no dependencies but Apple's Accelerate on macOS, which ships with the OS.

Until a body lands, it returns `error.NotImplemented`, which is in its error set until then, and
its tests skip (each file's `skip`). Signatures below are abbreviated: `v` is a view, `T` its
element type, `F` a float type.

## Background

### num.zig

A survey of num.zig (<https://github.com/muhammad-fiaz/num.zig>, MIT, v0.0.4, at a43e109) found
the areas it covers that rstd lacks. Only ideas and API shape come from it, no code: its arrays
carry their element type at run time, while rstd dispatches on comptime types. Its files that
informed each section are named there.

### FFT libraries

- **No usable Zig FFT exists:** none in Zig's std, and none on Codeberg or GitHub. fftiny.zig's
  author advises against it, cliFFTop's scope is undocumented, num.zig's is radix 2 with a DFT
  fallback, and zfft-orchard is abandoned.
- **FFTW is out:** it is GPLv2+, and its `FFTW_MEASURE` planner chooses by timing, so the bits
  can differ from run to run. Its speed comes from `genfft`, an OCaml code generator.
- **PocketFFT** (BSD-3, NumPy's and SciPy's FFT): its C++ header (`pocketfft_hdronly.h`, cpp
  branch, 3,888 lines) has complex and real transforms for any `n` (mixed radix 2, 3, 4, 5, 7,
  11, plus Bluestein), DCT and DST types 1 to 4, transforms over several axes, and threads.
- **pffft** (marton78's fork, active, single and double precision): FFTPACK's UCAR license, the
  same three conditions as BSD-3. Its SIMD runs within one transform, sizes `2^a 3^b 5^c`, and
  its unordered output speeds up convolution.
- **Accelerate:** vDSP's DFT takes sizes `f 2^n` with `f` in 1, 3, 5, 15; other sizes stay
  native. Its real transforms pack and scale differently, so routing converts to rstd's layout.
- **VkFFT** is a separate library, not part of Vulkan, so it is out.
- **Other references:** lifft (MIT); djbfft (public domain, power-of-two sizes only).

### Licenses

PocketFFT (BSD-3) and pffft (UCAR) are approved as OpenBLAS is: a `LICENSE` section (copyright
notice and license text) with the first ported code, and a provenance line in each ported file.
For sparse matrices, SuiteSparse's orderings (AMD, COLAMD) are BSD-3, but its factorizations
(CHOLMOD, UMFPACK, KLU) and CSparse are LGPL or GPL: direct solvers need a license decision first.

## Design

### Statistics: `rstd.stats`

[`src/stats.zig`](src/stats.zig), over slices. The forms along a dimension are planned with N-d
arrays; `cov` and `corrcoef` with the public layer's `ConstView` and `Mat`.

```zig
pub const Normalize = enum { sample, population }; // n - 1 (MATLAB's default) or n
pub fn mean(io, xs: []const F) error{ EmptyInput, Canceled }!F
pub fn variance(io, xs, normalize: Normalize) error{ EmptyInput, Canceled }!F
pub fn stdDev(io, xs, normalize: Normalize) error{ EmptyInput, Canceled }!F
pub fn Running(comptime F: type) type; // streaming mean and variance (Welford's update)
pub fn median(allocator, xs) error{ EmptyInput, OutOfMemory }!F // selects on a copy
pub const QuantileMethod = enum {
    inverted_cdf, averaged_inverted_cdf, closest_observation, interpolated_inverted_cdf,
    hazen, weibull, linear, median_unbiased, normal_unbiased,
};
pub const QuantileOptions = struct { method: QuantileMethod = .linear };
pub fn quantile(allocator, xs, p: []const F, options) error{ ... }![]F
pub fn quantileInto(allocator, out: []F, xs, p: []const F, options) error{ ... }!void
pub fn Bins(comptime F: type) type; // union(enum) { count: usize, edges: []const F }
pub fn histogram(allocator, xs, bins: Bins(F)) error{ ... }!Histogram(F) // counts, edges
pub fn histogramInto(counts: []u64, xs, edges: []const F) error{ ... }!void // adds to counts
// planned
pub fn cov(allocator, io, x: ConstView(F)) error{ ... }!Mat(F) // columns are variables
pub fn corrcoef(allocator, io, x: ConstView(F)) error{ ... }!Mat(F)
```

- **Quantile methods** are Hyndman and Fan's nine (The American Statistician 50(4), 1996), under
  NumPy's names; num.zig has five of them. They differ in where `p` falls between neighboring
  sorted values: for `[1, 2, 3, 4]` at 0.25, `linear` (type 7; NumPy, R, Excel's `.INC`) gives
  1.75, `hazen` (type 5; MATLAB) 1.5, `weibull` (type 6; SPSS, Excel's `.EXC`) 1.25. The gap
  shrinks below one spacing as `n` grows. The default is `.linear`, the modern common choice; a
  caller matching MATLAB passes `.{ .method = .hazen }`.
- **NaN** propagates, as in MATLAB. Planned: an option to omit NaNs.
- **Variance** takes two passes (mean, then squared deviations) through the BLAS `sum`, and a
  sample of one has variance 0 either way, as MATLAB's `var`. `Running` is the streaming form;
  its bits differ from `variance`'s. `cov` is `syrk` over the centered columns.
- **Histogram:** the last bin includes its right edge, as MATLAB's `histcounts` and NumPy do.
  Values outside the edges, and NaNs, are not counted. Errors include `EdgesNotSorted`,
  `ZeroBins` and `SliceTooShort`.
- **Status:** `Running` is real; the rest are stubs.
- num.zig's `stats/describe.zig` and `stats/correlate.zig` informed this.

### Sorting and sets: `rstd.sort`

[`src/sort.zig`](src/sort.zig). `std.sort` already has the algorithms (`pdq`, stable `block`,
`lowerBound`, `upperBound`, `partitionPoint`); this module adds what it lacks.

```zig
pub const Direction = enum { ascend, descend };
pub fn order(comptime T: type, comptime direction: Direction) fn (void, T, T) bool
pub fn sort(xs: []T, direction: Direction) void // stable
pub fn argsort(allocator, xs, direction) error{OutOfMemory}![]usize // MATLAB's [~, i] = sort
pub fn searchSorted(xs, x, direction, side: Side) usize // Side = enum { left, right }
pub const Keep = enum { sorted, stable }; // every set operation's output order
pub fn Unique(comptime T: type) type; // values, first, inverse: MATLAB's C, ia, ic
pub fn unique(allocator, xs, keep) error{OutOfMemory}!Unique(T)
pub fn intersect(allocator, a, b, keep) error{OutOfMemory}!Matched(T) // values, a_index, b_index
pub fn setUnion(allocator, a, b, keep) error{OutOfMemory}!Matched(T) // `union` is a keyword
pub fn setDiff(allocator, a, b, keep) error{OutOfMemory}!Selected(T) // values, index into a
pub fn isMember(allocator, a, set) error{OutOfMemory}![]?usize // ismember's; found: != null
```

- **Floats:** `std.sort.asc(f64)` is no strict weak order once a NaN appears, so the result order
  is unspecified. `order` puts NaNs last ascending and first descending, as MATLAB does; -0 and +0
  compare equal and keep their order. `unique` treats each NaN as distinct, as MATLAB does.
- `searchSorted` takes the direction the slice is sorted in. Every set operation takes `Keep`, as
  MATLAB's `'stable'`. `setDiff` has no `b` side, so it returns `Selected(T)`.
- **Status:** `order`, `sort`, `argsort` and `searchSorted` are real, over `std.sort`'s `block`
  and `partitionPoint`; the set operations are stubs.
- num.zig's `sort/ordering.zig` and `sort/search.zig` informed this.

### FFT: `rstd.fft`

[`src/fft/__root__.zig`](src/fft/__root__.zig). PocketFFT's algorithms in Zig as the baseline
(its C++ header), pffft's SIMD kernels ported for speed, Accelerate on macOS by default. The use
cases are audio, numerics and images.

```zig
pub const Direction = enum { forward, backward }; // forward is exp(-2πi jk/n), as NumPy and FFTW
pub const Norm = enum { backward, ortho, forward }; // which way divides by n; ortho: √n both ways
pub const Options = struct { norm: Norm = .backward }; // NumPy's default: round trips are exact
pub const Layout = struct { count: usize = 1, stride: usize = 1, distance: ?usize = null };
pub const Layouts = struct { in: Layout = .{}, out: Layout = .{} }; // out of place

// complex, any n; immutable once built, so one plan serves every thread
pub fn Plan(comptime T: type) type; // init(allocator, n), deinit(allocator), scratchLen()
fn transform(plan, io, data: []Complex(T), scratch: []Complex(T), dir, layout, options) !void
fn transformInto(plan, io, out, in: []const Complex(T), scratch, dir, layouts, options) !void

// real input: n reals to the n/2 + 1 complex values of the half spectrum (NumPy's rfft)
pub fn RealPlan(comptime T: type) type;
fn forward(plan, io, out: []Complex(T), in: []const T, scratch: []T, layouts, options) !void
fn backward(plan, io, out: []T, in: []const Complex(T), scratch: []T, layouts, options) !void

// fast convolution in the plan's own order (pffft's), skipping the reordering both ways
fn forwardUnordered(plan, out: []T, in: []const T, scratch) !void
fn convolveAccumulate(plan, acc: []T, a: []const T, b: []const T, scale: T) !void
fn backwardUnordered(plan, out: []T, in: []const T, scratch) !void

// comptime size: factors and twiddles at comptime, unrolled, no allocator; the same calls
pub fn StaticPlan(comptime n: usize, comptime T: type) type;
pub fn StaticRealPlan(comptime n: usize, comptime T: type) type;

// DCT and DST types 1 to 4; the MDCT (audio codecs) on DCT-IV
pub const TrigKind = enum { dct1, dct2, dct3, dct4, dst1, dst2, dst3, dst4 };
pub fn TrigPlan(comptime T: type) type; // init(allocator, n, kind)
pub fn Mdct(comptime T: type) type; // init(allocator, n), forward, backward

// planned, with N-d arrays: several axes of a view; plans per axis, scratch allocated
pub fn transformAxes(allocator, io, v: ArrayView(Complex(T)), axes: []const u8, dir, options) !void
pub fn realForwardAxes(allocator, io, out, in: ConstArrayView(T), axes, options) !void
```

- **Real time:** a transform allocates nothing; the plan and the scratch are made up front.
  `scratchLen()` is per transform, any layout: at least `n` (PocketFFT's ping-pong buffer), more
  for Bluestein's. Errors: `SliceTooShort`, `CountMismatch` (out of place), and `Canceled` from
  `io`. A plan's `init` gives `OutOfMemory` and `ZeroLength`.
- **Layouts:** transform `k`'s element `i` is at `k * distance + i * stride`. A side's `distance`
  defaults to its elements per transform times its stride, so `.{}` is contiguous on both sides
  even where `n` reals become `n / 2 + 1` complex values. Transforms must not overlap.
- **`io`** runs a batch's transforms concurrently, as many at once as the scratch holds
  (`scratch.len / scratchLen()`). They are independent, so the bits depend on neither. Planned: a
  single huge transform's threads, by a fixed split.
- **Frames:** the unordered forms and the `Mdct` take no `io` or layout, one frame at a time. The
  `Mdct` leaves windowing to the caller; its scaling makes windowed overlap-add with a
  Princen-Bradley window (sine, KBD) reconstruct the input (TDAC).
- **Element types:** every float type through the portable code, f16 to f128. Planned: f32 SIMD
  first (audio), then f64 (numerics).
- **Bits:** the plan depends on `n`, the type and the core only, never on timing. The SIMD
  kernels change the operation order, so bits may differ between cores, as with the BLAS.
  `StaticPlan` dispatches by size alone, so it gives `Plan`'s bits.
- **Tests:** a direct DFT in a wider type (as the BLAS's `Wide` reference), round trips,
  Parseval's identity, and every size from 1 to a few hundred (primes take Bluestein's path).
- **Status:** stubs; `init` already rejects `n = 0`.

### Planned: signal processing, `rstd.signal`

Built on `rstd.fft`, for audio first:

- **Window functions:** Hann, Hamming, Blackman, Kaiser, and the Princen-Bradley windows the
  `Mdct` needs (sine, KBD); symmetric or periodic (MATLAB's `'symmetric'`, `'periodic'`).
- **STFT** and its inverse by weighted overlap-add, over `RealPlan`; hop, window and FFT size
  free.
- **Partitioned convolution:** uniform-partitioned overlap-add over the unordered forms, for long
  filters (reverb) in real time: fixed latency, no allocation per block.

### Planned: N-d arrays

In `rstd.linalg` beside `Mat(T)`, in two kinds, as the public layer's `StaticMat` and
`Mat`:

- **`StaticArray`:** the shape, and so the rank, is comptime. It's a value with no allocator,
  loops are unrolled, rank-changing operations get their result types at comptime, and shape
  errors are compile errors.
- **`Array` and its views:** rank and shape are run time, up to `max_rank`, for sizes known only
  at run time, large data and files. Shape errors are errors.

Both give views, every operation's core takes views, and dispatch is by size alone (design.md's
Dispatch), so the same numbers give the same bits in either. Strides are column-major (the first
index varies fastest), as in MATLAB, LAPACK and `Mat(T)`, so a rank-2 view with a first stride
of 1 is a `View(T)` without a copy.

```zig
pub const max_rank = 8;
pub fn StaticArray(comptime dims: anytype, comptime T: type) type; // .{ 4, 4, 3 }; data: [n]T
pub const Shape = struct { dims: [max_rank]usize, rank: u8 };
pub fn Array(comptime T: type) type; // data: []T, shape; init(allocator, dims), deinit(allocator)
pub fn ArrayView(comptime T: type) type; // data, shape, strides: [max_rank]isize; also Const
pub const Range = struct { start: usize, end: usize, step: isize = 1 };
pub const Index = union(enum) { all, at: usize, range: Range };

// views, no copy
fn slice(v, indices: []const Index) error{ IndexOutOfRange, RankMismatch }!ArrayView(T)
fn permute(v, order: []const u8) error{InvalidPermutation}!ArrayView(T) // MATLAB's permute
fn reshape(v, dims: []const usize) error{ CountMismatch, NotContiguous }!ArrayView(T)
fn squeeze(v) ArrayView(T)
fn broadcastTo(v, shape: Shape) error{ RankMismatch, DimensionMismatch }!ArrayView(T) // stride 0
fn asMatrix(v) error{ RankMismatch, NotColumnMajor }!View(T) // rank 2, first stride 1
pub fn broadcastShape(a: Shape, b: Shape) error{ RankMismatch, DimensionMismatch }!Shape

// iteration in column-major order: element offsets into n views over one broadcast shape
pub fn Iterator(comptime n: usize) type;

// copies and element-wise
pub fn copy(allocator, v) error{OutOfMemory}!Array(T) // contiguous; reshape's copying form
pub fn map(out: ArrayView(T), a: anytype, comptime f: anytype) error{DimensionMismatch}!void
pub fn zip(out, a, b, comptime f) error{DimensionMismatch}!void // broadcasts a and b to out
pub fn cat(allocator, dim: u8, parts: []const ConstArrayView(T)) error{ ... }!Array(T)
pub fn repmat(allocator, v, counts: []const usize) error{OutOfMemory}!Array(T)
pub fn sumDim(allocator, io, v, dim: u8) error{ ... }!Array(T) // and prod, max, min, mean
```

- **Broadcasting** needs equal ranks: a dimension of 1 expands to match, and different ranks are
  `error.RankMismatch`, never padded. Reductions keep the reduced dimension as 1 (as MATLAB's
  `mean`), so their results still broadcast back. If padding is ever added, it follows MATLAB's
  implicit expansion (trailing 1s, free in column-major), not NumPy's leading 1s.
- **Reductions along a dimension** use the BLAS `sum` where that dimension is contiguous, and fixed
  chunks combined in order elsewhere.
- num.zig's `core/shape.zig`, `core/iterator.zig`, `ops/broadcast.zig` and `manip/` informed this.

### Planned: sparse matrices

In `rstd.linalg.sparse`, with 0-based indices.

```zig
pub fn Csc(comptime T: type) type; // rows, cols: u32, col_start: []u32, row: []u32, value: []T
pub fn Csr(comptime T: type) type; // the row-major twin: a Csc of A is a Csr of Aᵀ, same arrays
pub fn Triplets(comptime T: type) type; // builder: append(allocator, i, j, x), then toCsc
fn toCsc(allocator, t) error{OutOfMemory}!Csc(T) // sums duplicates in insertion order
fn fromDense(allocator, a: ConstView(T), drop_below: F) error{OutOfMemory}!Csc(T)
fn toDense(allocator, a) error{OutOfMemory}!Mat(T)
fn mulVec(io, a, trans: Trans, alpha: T, x: []const T, beta: T, y: []T) error{ ... }!void
fn mulDense(allocator, io, a, b: ConstView(T)) error{ ... }!Mat(T)

// iterative solvers, matrix-free: anything with apply
pub fn Operator(comptime T: type) type; // context plus apply(io, x: []const T, y: []T)
pub fn Options(comptime T: type) type; // tol (relative residual), max_iterations, preconditioner
pub const Stats = struct { iterations: usize, residual: f64 };
pub fn cg(allocator, io, a: Operator(T), b: []const T, x: []T, options) error{ ... }!Stats
pub fn gmres(allocator, io, a, b, x, options, restart: usize) error{ ... }!Stats
```

- **Indices:** `u32` throughout (shape, column starts, rows), so under 2³² rows, columns and
  nonzeros; widening is an easy change if ever needed. Reads need no cast (`u32` coerces to
  `usize`); only the constructors narrow, checked (`error.TooLarge`). With `f64` values a nonzero
  costs 12 bytes, not 16, and `A x` is bandwidth-bound.
- **CSC comes first**, as in MATLAB and SuiteSparse. A threaded `A x` on CSC combines per-chunk
  partial vectors in chunk order (the BLAS's ordered combine); on CSR each row's sum stays whole.
- **Errors** beyond `OutOfMemory` and `Canceled`: `DimensionMismatch`, `NoConvergence`,
  `NotPositiveDefinite` (CG's breakdown).
- num.zig's `core/sparse.zig` and `linalg/sparse_solve.zig` (CG, restarted GMRES, no
  preconditioner) informed this.

## Decisions

By the user, 2026-10-31:

- `rstd.stats`, `rstd.sort` and `rstd.fft` are top-level modules; N-d arrays and sparse matrices
  go in `rstd.linalg`.
- Broadcasting needs equal ranks; MATLAB's implicit expansion is the fallback if padding is ever
  added.
- Sparse indices are `u32` throughout.
- Quantiles take `QuantileOptions`, `.linear` by default.
- Arrays come in two kinds, `StaticArray` (comptime shape) and `Array` (run time), as the BLAS
  design's `StaticMat` and `Mat`.
- FFT: Zig-native, PocketFFT as the baseline and pffft to inform the SIMD; Accelerate on macOS
  only; no VkFFT or other dependency. Audio, numerics and images are all in scope: audio's SIMD
  comes first, every other feature gets a stub so it can follow.
- PocketFFT's and pffft's licenses are approved as BSD-3 is (Licenses, above).
- The API sketches are approved. N-d arrays and sparse matrices start after the public layer
  (linalg's M3), on its real `View(T)` and `Mat(T)`.
- The stubs' changes to the sketches stand for now (explicit `comptime` element types,
  `stats.quantile`, `sort.Keep`, `Selected(T)`, `[]?usize`, `fft.Layouts`, scratch-bounded
  batches); they may be revisited.
- The FFT comes before statistics and sets.
- PocketFFT's source is its C++ header (cpp branch), pinned to one commit for provenance.
- Accelerate's FFT routing is on by default on macOS, as the BLAS's, with a way to force native.
- Kernels follow the BLAS's strategy: port the upstream SIMD kernels first, faithfully; native
  rewrites come later and stay only where a benchmark backs them. Where an upstream kernel
  repeats one rstd already has (BLAS, LAPACK), call rstd's instead of porting it.
- pffft is C intrinsics behind a macro layer, so its faithful port maps each macro to the
  `@Vector` or `@shuffle` form that emits the same instruction, in pffft's order, without fused
  multiply-adds pffft lacks; checked by reading the code generated and by matching pffft's bits.
  The port is a chance to improve `rstd.linalg`'s APIs for such code: its
  vector primitives (interleave, transpose, half-swap) go in the `Vec` and `StaticMat` code that
  linalg's M7 redesigns around `@Matrix`, which also decides the N-d arrays' fixed-size kind.
- Kernels per platform: pffft's own first (SSE on x86_64, NEON on aarch64 and the Pi 3's 32-bit
  arm), then derived or new ones for other tier 1 and 2 targets as needed; where none fits (the
  PowerPC 750 has no AltiVec), generic `@Vector` or `@Matrix` code.
- The bit guarantee is the BLAS's: identical for any thread count on one machine, or one kind of
  machine. Results may differ between cores.
- `stats` sums through the BLAS and takes its accuracy from there. Compensated summation becomes
  a BLAS option (linalg's T138); benchmarks choose the defaults, there and per use case.
- A caller chooses a result type explicitly, as `@as` does for `@floatFromInt`: Zig gives a
  function no access to its result type. So `stats.mean(io, R, xs)` takes the result type `R`,
  and the element type comes from `xs`; for floats today's calls keep their form. The matlab
  port gets the same (linalg's M6).
- Signal processing gets its own module, `rstd.signal`: window functions, the STFT and
  partitioned convolution (Planned: signal processing, above).
- Front ends are type-agnostic, as the BLAS's: any element type, dispatched on it, integers and
  `Complex` included where the operation means something to them.
- Complex values sort in MATLAB's order (by magnitude, then angle), and everything that orders
  them (`sort`, `unique`, the set operations, `median`) uses it.
