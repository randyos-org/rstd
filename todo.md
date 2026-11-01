# rstd to-do

The backlog of rstd beyond a single module's own, in priority order: work goes top to bottom. IDs
stay fixed when items move. **Needs:** names items that must land first; **Blocked:** names an
outside cause. The design is in [`design.md`](design.md). Linear algebra keeps its own backlog
in [`src/linalg/todo.md`](src/linalg/todo.md); its MVP comes first (user, 2026-10-30).

## In progress

### Working on

Nothing here; linalg's T139 and M3 come first (its todo.md).

## Queue

### FFT

- [ ] **N7** `fft.Plan`: complex 1-D, portable, every float type (PocketFFT's mixed radix and
  Bluestein, from its C++ header at a pinned commit). The correctness baseline every later FFT
  item is tested against. Add PocketFFT's `LICENSE` section and per-file provenance with the
  first ported code. **Needs:** nothing (N6 done).
- [ ] **N8** `fft.RealPlan` (PocketFFT's real transforms). **Needs:** N7.
- [ ] **N9** `Layout` batches, run concurrently through `io`, as many at once as the scratch
  holds. **Needs:** N7.
- [ ] **N10** Audio: port pffft's f32 SIMD kernels faithfully (design.md, Decisions), within one
  transform, sizes `2^a 3^b 5^c`, and the unordered convolution (`forwardUnordered`,
  `convolveAccumulate`, `backwardUnordered`): SSE on x86_64, NEON on aarch64 and 32-bit arm,
  generic code on the PowerPC 750. Where a kernel repeats one rstd has (BLAS, LAPACK), call
  rstd's. Its vector primitives (interleave, transpose, half-swap) go in the `Vec` and
  `StaticMat` code linalg's M7 designs. Add pffft's `LICENSE` section with its first ported code.
  **Needs:** N8, M7.
- [ ] **N35** FFT kernels for tier 1 and 2 targets pffft lacks or serves poorly (AVX on x86_64,
  say), derived from its own or new. **Needs:** N10.

### Signal processing

- [ ] **N29** `rstd.signal` window functions: Hann, Hamming, Blackman, Kaiser, sine and KBD,
  symmetric or periodic. **Needs:** nothing (N6 done).
- [ ] **N30** `signal` STFT and its inverse (weighted overlap-add) over `RealPlan`. **Needs:**
  N8, N29.
- [ ] **N31** `signal` uniform-partitioned convolution over the unordered forms: fixed latency,
  no allocation per block. **Needs:** N10.

### Statistics

- [ ] **N11** `stats.mean`, `variance` and `stdDev` through the BLAS `sum`, with its accuracy
  options once linalg's T138 adds them. **Needs:** nothing (N6 done).
- [ ] **N12** `stats.median` and the quantiles, all nine methods. **Needs:** nothing (N6 done).
- [ ] **N13** `stats.histogram` and `histogramInto`. **Needs:** nothing (N6 done).
- [ ] **N33** Type-agnostic `rstd.stats`, as the BLAS: integer and `Complex` inputs beside
  floats, the result type chosen by the caller (`mean(io, R, xs)`, the element type from `xs`).
  **Needs:** N11.

### Sorting and sets

- [ ] **N14** `sort.unique`, `intersect`, `setUnion`, `setDiff` and `isMember`, each NaN distinct.
  **Needs:** nothing (N6 done).
- [ ] **N32** `Complex` in `rstd.sort`, in MATLAB's order (by magnitude, then angle), used by
  every ordering operation, `stats.median` included. **Needs:** nothing (N6 done).

### After the public layer

- [ ] **N1** N-d array stubs (design.md, "Planned: N-d arrays"). **Needs:** M3 (linalg's
  todo.md), whose `View(T)` and `Mat(T)` they extend.
- [ ] **N2** Sparse matrix stubs (design.md, "Planned: sparse matrices"). **Needs:** M3.
- [ ] **N15** `stats.cov` and `corrcoef`: `syrk` over the centered columns. **Needs:** M3.

## Later

Picked up after the Queue, in this order.

- [ ] **N16** FFT: f64 SIMD (pffft's double) and a cache-blocked split for `n` past L2, for
  numerics. **Needs:** N10.
- [ ] **N17** FFT: SIMD across a batch's transforms, for images. **Needs:** N9.
- [ ] **N18** `fft.TrigPlan` (DCT and DST types 1 to 4) and `Mdct`. **Needs:** N8.
- [ ] **N19** `fft.StaticPlan` and `StaticRealPlan`. **Needs:** N8.
- [ ] **N20** FFT: route vDSP's sizes (`f 2^n`, `f` in 1, 3, 5, 15) to Accelerate on macOS by
  default, as the BLAS does, converting its real layout and scaling; a way to force native.
  **Needs:** N8.
- [ ] **N21** FFT: threads within one large transform, by a fixed split. **Needs:** N9.
- [ ] **N22** N-d array bodies. **Needs:** N1.
- [ ] **N23** `fft.transformAxes` and `realForwardAxes`. **Needs:** N22.
- [ ] **N24** `rstd.stats` along a dimension. **Needs:** N22.
- [ ] **N25** `rstd.stats`: an option to omit NaNs.
- [ ] **N26** Sparse bodies: `Csc` first, `mulVec`, `cg` and `gmres`. A likely first user is
  `linalg/optim`'s sparse Jacobian. **Needs:** N2.
- [ ] **N27** Sparse: BiCGSTAB and MINRES. **Needs:** N26.
- [ ] **N28** Sparse direct factorizations. **Blocked:** a license decision (design.md,
  Licenses).
- [ ] **N34** FFT: native rewrites of pffft's kernels, kept only where a benchmark shows them
  faster, or as fast and cleaner. **Needs:** N10.

## Done

- [X] The survey and API design ([`design.md`](design.md)), approved by the user 2026-10-31.
- [X] **N3** to **N5** The `rstd.stats`, `rstd.sort` and `rstd.fft` stubs, 2026-10-31; the real
  parts are `Running`, `order`, `sort`, `argsort` and `searchSorted`.
- [X] **N6** Their first compile, 2026-11-01: each file's `zig test` alone on Zig 0.18 and 0.17,
  the real parts passing and the stubs skipping; rstd's own tests with T139.

## Waiting on the user

Nothing waiting.
