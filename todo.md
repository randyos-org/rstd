# rstd to-do

Work on rstd beyond a single module's backlog, in priority order. IDs stay fixed when items move.
Linear algebra keeps its own backlog in [`src/linalg/todo.md`](src/linalg/todo.md) (the BLAS,
LAPACK, the public layer, the matlab port). Its MVP comes first, and nothing here starts before
it (user, 2026-10-30).

## New modules, after the MVP

Planned 2026-10-31 from a survey of num.zig (<https://github.com/muhammad-fiaz/num.zig>, MIT,
v0.0.4, at a43e109): the areas it covers that rstd lacks. Only ideas and API shape come from it,
no code: its arrays carry their element type at run time, while rstd dispatches on comptime types.
The sketches are approved (user, 2026-10-31). **On hold until the BLAS work settles** (user,
2026-10-31), stubs included. Then each item starts as stubs: the module with the sketched
signatures, bodies returning `error.NotImplemented`, and tests that skip until implemented (as
`linalg/optim`'s `solveOrSkip`). N3, N4 and N5 take slices and come first; N1 and N2 follow M3,
on its real `View(T)` and `Matrix(T)`.

- [ ] **N1** N-d arrays: strided views, broadcasting, an N-d iterator
  ([sketch](#n1-n-d-arrays)). **Needs:** M3's views and `Matrix(T)` (linalg's todo.md), which
  they extend.
- [ ] **N2** Sparse matrices and iterative solvers ([sketch](#n2-sparse-matrices)). A likely first
  user is `linalg/optim`'s sparse Jacobian.
- [ ] **N3** Statistics ([sketch](#n3-statistics)). The forms along a dimension need N1.
- [ ] **N4** Sorting and set operations ([sketch](#n4-sorting-and-sets)).
- [ ] **N5** FFT ([sketch](#n5-fft)): PocketFFT's algorithms in Zig, pffft's SIMD for audio
  first, stubs for the rest. The multi-axis forms need N1.

## API sketches

Shared conventions are those of the public `rstd.linalg` layer
([`design.md`](src/linalg/design.md), "Operations over views"):

- `(allocator, io, ...)` where a function allocates, `(io, ...)` where it only calls the BLAS;
  allocator first.
- Types are unmanaged and store no allocator. A new result comes back owned, and each such
  function has an `...Into` form that writes caller storage.
- Error sets are exact and named, with plain error names. Operands are slices, not pointers and
  lengths. Storage is column-major.
- Results are bit-identical for any thread count, as the BLAS's are
  ([`threading.md`](src/linalg/blas/threading.md)).

Signatures are abbreviated: `v` is a view, `T` its element type, `F` a float type.

### N1: N-d arrays

These go in `rstd.linalg` beside `Matrix(T)`, in two kinds (user, 2026-10-31), as the public
layer's `StaticMatrix` and `Matrix`:

- **`StaticArray`:** the shape, and so the rank, is comptime. It's a value with no allocator,
  loops are unrolled, rank-changing operations get their result types at comptime, and shape
  errors are compile errors.
- **`Array` and its views:** rank and shape are run time, up to `max_rank`, for sizes known only
  at run time, large data and files. Shape errors are errors.

Both give views, every operation's core takes views, and dispatch is by size alone (design.md's
Dispatch), so the same numbers give the same bits in either. Strides are column-major (the first
index varies fastest), as in MATLAB, LAPACK and `Matrix(T)`, so a rank-2 view with a first stride
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

- **Broadcasting** (user, 2026-10-31) needs equal ranks: a dimension of 1 expands to match, and
  different ranks are `error.RankMismatch`, never padded. Reductions keep the reduced dimension as
  1 (as MATLAB's `mean`), so their results still broadcast back. If padding is ever added, it
  follows MATLAB's implicit expansion (trailing 1s, free in column-major), not NumPy's leading 1s.
- **Reductions along a dimension** use the BLAS `sum` where that dimension is contiguous, and fixed
  chunks combined in order elsewhere.
- num.zig's `core/shape.zig`, `core/iterator.zig`, `ops/broadcast.zig` and `manip/` informed this.

### N2: Sparse matrices

These go in `rstd.linalg.sparse`, with 0-based indices. N3, N4 and N5 are top-level modules
(user, 2026-10-31).

```zig
pub fn Csc(comptime T: type) type; // rows, cols: u32, col_start: []u32, row: []u32, value: []T
pub fn Csr(comptime T: type) type; // the row-major twin: a Csc of A is a Csr of Aᵀ, same arrays
pub fn Triplets(comptime T: type) type; // builder: append(allocator, i, j, x), then toCsc
fn toCsc(allocator, t) error{OutOfMemory}!Csc(T) // sums duplicates in insertion order
fn fromDense(allocator, a: ConstView(T), drop_below: F) error{OutOfMemory}!Csc(T)
fn toDense(allocator, a) error{OutOfMemory}!Matrix(T)
fn mulVec(io, a, trans: Trans, alpha: T, x: []const T, beta: T, y: []T) error{ ... }!void
fn mulDense(allocator, io, a, b: ConstView(T)) error{ ... }!Matrix(T)

// iterative solvers, matrix-free: anything with apply
pub fn Operator(comptime T: type) type; // context plus apply(io, x: []const T, y: []T)
pub fn Options(comptime T: type) type; // tol (relative residual), max_iterations, preconditioner
pub const Stats = struct { iterations: usize, residual: f64 };
pub fn cg(allocator, io, a: Operator(T), b: []const T, x: []T, options) error{ ... }!Stats
pub fn gmres(allocator, io, a, b, x, options, restart: usize) error{ ... }!Stats
```

- **Indices** (user, 2026-10-31): `u32` throughout (shape, column starts, rows), so under 2³²
  rows, columns and nonzeros; widening is an easy change if ever needed. Reads need no cast (`u32`
  coerces to `usize`); only the constructors narrow, checked (`error.TooLarge`). With `f64`
  values a nonzero costs 12 bytes, not 16, and `A x` is bandwidth-bound.
- **CSC comes first**, as in MATLAB and SuiteSparse. A threaded `A x` on CSC combines per-chunk
  partial vectors in chunk order (the BLAS's ordered combine); on CSR each row's sum stays whole.
- **Errors** beyond `OutOfMemory` and `Canceled`: `DimensionMismatch`, `NoConvergence`,
  `NotPositiveDefinite` (CG's breakdown).
- **Later:** BiCGSTAB and MINRES. Direct factorizations need a license decision first:
  SuiteSparse's orderings (AMD, COLAMD) are BSD-3, but its factorizations (CHOLMOD, UMFPACK, KLU)
  and CSparse are LGPL or GPL.
- num.zig's `core/sparse.zig` and `linalg/sparse_solve.zig` (CG, restarted GMRES, no
  preconditioner) informed this.

### N3: Statistics

These go in `rstd.stats`, over slices first; the forms along a dimension come with N1.

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
pub fn quantileInto(allocator, out: []F, xs, p: []const F, options) error{ ... }!void
pub fn cov(allocator, io, x: ConstView(F)) error{ ... }!Matrix(F) // columns are variables
pub fn corrcoef(allocator, io, x: ConstView(F)) error{ ... }!Matrix(F)
pub const Bins = union(enum) { count: usize, edges: []const F };
pub fn histogram(allocator, xs, bins: Bins) error{ ... }!Histogram(F) // counts, edges
pub fn histogramInto(counts: []u64, xs, edges: []const F) error{ ... }!void
```

- **Quantile methods** are Hyndman and Fan's nine (The American Statistician 50(4), 1996), under
  NumPy's names; num.zig has five of them. They differ in where `p` falls between neighboring
  sorted values: for `[1, 2, 3, 4]` at 0.25, `linear` (type 7; NumPy, R, Excel's `.INC`) gives
  1.75, `hazen` (type 5; MATLAB) 1.5, `weibull` (type 6; SPSS, Excel's `.EXC`) 1.25. The gap
  shrinks below one spacing as `n` grows. The default is `.linear`, the modern common choice
  (user, 2026-10-31); a caller matching MATLAB passes `.{ .method = .hazen }`.
- **NaN** propagates by default, as in MATLAB; an option to omit NaNs comes later.
- **Variance** takes two passes (mean, then squared deviations) through the BLAS `sum`. `cov` is
  `syrk` over the centered columns.
- **Histogram:** the last bin includes its right edge, as MATLAB's `histcounts` and NumPy do.
  Errors include `EdgesNotSorted` and `SliceTooShort`.
- num.zig's `stats/describe.zig` and `stats/correlate.zig` informed this.

### N4: Sorting and sets

These go in `rstd.sort`. `std.sort` already has the algorithms (`pdq`, stable `block`,
`lowerBound`, `upperBound`, `equalRange`); this module adds what it lacks.

```zig
pub const Direction = enum { ascend, descend };
pub fn order(comptime T: type, direction: Direction) fn (void, T, T) bool // NaN-aware for floats
pub fn sort(xs: anytype, direction: Direction) void // stable
pub fn argsort(allocator, xs, direction) error{OutOfMemory}![]usize // MATLAB's [~, i] = sort
pub fn searchSorted(xs, x, side: enum { left, right }) usize
pub fn Unique(comptime T: type) type; // values, first, inverse: MATLAB's C, ia, ic
pub fn unique(allocator, xs, keep: enum { sorted, stable }) error{OutOfMemory}!Unique(T)
pub fn intersect(allocator, a, b) error{OutOfMemory}!Matched(T) // values, a_index, b_index
pub fn setUnion(allocator, a, b) error{OutOfMemory}!Matched(T) // `union` is a keyword
pub fn setDiff(allocator, a, b) error{OutOfMemory}!Matched(T)
pub fn isMember(allocator, a, set) error{OutOfMemory}!Membership // found: []bool, at: []?usize
```

- **Floats:** `std.sort.asc(f64)` is no strict weak order once a NaN appears, so the result order
  is unspecified. `order` puts NaNs last ascending and first descending, as MATLAB does; -0 and +0
  compare equal and keep their order. `unique` treats each NaN as distinct, as MATLAB does.
- num.zig's `sort/ordering.zig` and `sort/search.zig` informed this.

### N5: FFT

These go in `rstd.fft`. Decided (user, 2026-10-31): Zig-native, PocketFFT as the baseline, and
pffft (license approved like BSD-3) to inform the SIMD. Accelerate on macOS only, since it ships
with the OS; no other dependencies, so no VkFFT (a separate library, not part of Vulkan). All
three use cases (audio, numerics, images) are in scope. Audio's SIMD comes first, and every other
feature gets a stub so it can follow, f64 numerics among them. The sketch is approved.

```zig
pub const Direction = enum { forward, backward }; // forward is exp(-2πi jk/n), as NumPy and FFTW
pub const Norm = enum { backward, ortho, forward }; // which way divides by n; ortho: √n both ways
pub const Options = struct { norm: Norm = .backward }; // NumPy's default: round trips are exact
pub const Layout = struct { count: usize = 1, stride: usize = 1, distance: ?usize = null };

// complex, any n; immutable once built, so one plan serves every thread
pub fn Plan(comptime T: type) type; // init(allocator, n), deinit(allocator), scratchLen()
fn transform(plan, io, data: []Complex(T), scratch: []Complex(T), dir, layout, options) !void
fn transformInto(plan, io, out, in: []const Complex(T), scratch, dir, layout, options) !void

// real input: n reals to the n/2 + 1 complex values of the half spectrum (NumPy's rfft)
pub fn RealPlan(comptime T: type) type;
fn forward(plan, io, out: []Complex(T), in: []const T, scratch, layout, options) !void
fn backward(plan, io, out: []T, in: []const Complex(T), scratch, layout, options) !void

// fast convolution in the plan's own order (pffft's), skipping the reordering both ways
fn forwardUnordered(plan, out: []T, in: []const T, scratch) !void
fn convolveAccumulate(plan, acc: []T, a: []const T, b: []const T, scale: T) !void
fn backwardUnordered(plan, out: []T, in: []const T, scratch) !void

// comptime size: factors and twiddles at comptime, unrolled, no allocator; the same calls
pub fn StaticPlan(comptime n: usize, comptime T: type) type;
pub fn StaticRealPlan(comptime n: usize, comptime T: type) type;

// several axes of an N1 view; plans per axis, scratch allocated
pub fn transformAxes(allocator, io, v: ArrayView(Complex(T)), axes: []const u8, dir, options) !void
pub fn realForwardAxes(allocator, io, out, in: ConstArrayView(T), axes, options) !void

// DCT and DST types 1 to 4; the MDCT (audio codecs) on DCT-IV
pub const TrigKind = enum { dct1, dct2, dct3, dct4, dst1, dst2, dst3, dst4 };
pub fn TrigPlan(comptime T: type) type; // init(allocator, n, kind)
pub fn Mdct(comptime T: type) type; // init(allocator, n), forward, backward
```

- **Real time:** a transform allocates nothing. The plan and the scratch (`scratchLen()`, 0 for
  most sizes; Bluestein's needs it) are made up front. Errors: `SliceTooShort`, and `Canceled`
  from `io`. A plan's `init` gives `OutOfMemory` and `ZeroLength`.
- **Element types:** every float type through the portable code, f16 to f128. f32 SIMD comes
  first (audio), then f64 (numerics).
- **`io`** runs a batch's transforms concurrently. They are independent, so the bits don't depend
  on the thread count; a single huge transform's threads come later, by a fixed split.
- **Bits:** the plan depends on `n`, the type and the core only, never on timing. The SIMD
  kernels change the operation order, so bits may differ between cores, as with the BLAS.
- **Tests:** a direct DFT in a wider type (as the BLAS's `Wide` reference), round trips,
  Parseval's identity, and every size from 1 to a few hundred (primes take Bluestein's path).
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
  native. Its real transforms pack and scale differently, so the routing converts to rstd's
  layout.
- **Other references:** lifft (MIT); djbfft (public domain, power-of-two sizes only).

Phases. Every function above is stubbed at the start; these fill them in order:

1. `Plan`: complex 1-D, portable, every float type (PocketFFT's mixed radix and Bluestein). The
   correctness baseline every later phase is tested against.
2. `RealPlan` (PocketFFT's real transforms).
3. `Layout` batches, run concurrently through `io`.
4. Audio: pffft-style f32 kernels within one transform, sizes `2^a 3^b 5^c`; the unordered
   convolution.
5. Expansions, each is its own later phase:
   - f64 SIMD (pffft's double) and a cache-blocked split for `n` past L2, for numerics.
   - SIMD across a batch's transforms, for images.
   - `transformAxes` (needs N1).
   - `TrigPlan` and `Mdct`.
   - `StaticPlan`.
   - Accelerate's routing.
   - Threads within one large transform.

## Open questions for the user

Decided 2026-10-31: top-level `rstd.stats`, `rstd.sort` and `rstd.fft`; `u32` sparse indices
throughout; broadcasting needs equal ranks; quantiles take `QuantileOptions`, `.linear` by
default; arrays in two kinds, `StaticArray` (comptime shape) and `Array` (run time); FFT as in
N5; N1 to N5's sketches approved; everything, stubs included, on hold until the BLAS work
settles, then N1's and N2's stubs after M3.

Nothing open.
