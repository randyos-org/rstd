# rstd to-do

The backlog of rstd beyond a single module's own, in priority order: work goes top to bottom. IDs
stay fixed when items move. **Needs:** names items that must land first; **Blocked:** names an
outside cause. The design is in [`design.md`](design.md). Linear algebra keeps its own backlog
in [`src/linalg/todo.md`](src/linalg/todo.md); its MVP comes first (user, 2026-10-30).

## In progress

### Working on

- [ ] **N3** `rstd.stats` stubs ([`src/stats.zig`](src/stats.zig)), written 2026-10-31, not yet
  compiled. `Running` is real. Left: N6.
- [ ] **N4** `rstd.sort` stubs ([`src/sort.zig`](src/sort.zig)), written 2026-10-31, not yet
  compiled. `order`, `sort`, `argsort` and `searchSorted` are real. Left: N6.
- [ ] **N5** `rstd.fft` stubs ([`src/fft/__root__.zig`](src/fft/__root__.zig)), written
  2026-10-31, not yet compiled. Left: N6.

## Queue

### Checks

- [ ] **N6** First compile and test of N3 to N5: `zig test` each file alone (they import only
  std), then rstd's tests; the real parts' tests pass and the stubs' skip. **Needs:** local
  builds again (stopped 2026-10-31).

### FFT

- [ ] **N7** `fft.Plan`: complex 1-D, portable, every float type (PocketFFT's mixed radix and
  Bluestein). The correctness baseline every later FFT item is tested against. Add PocketFFT's
  `LICENSE` section and per-file provenance with the first ported code. **Needs:** N6.
- [ ] **N8** `fft.RealPlan` (PocketFFT's real transforms). **Needs:** N7.
- [ ] **N9** `Layout` batches, run concurrently through `io`, as many at once as the scratch
  holds. **Needs:** N7.
- [ ] **N10** Audio: pffft-style f32 kernels within one transform, sizes `2^a 3^b 5^c`, and the
  unordered convolution (`forwardUnordered`, `convolveAccumulate`, `backwardUnordered`). Add
  pffft's `LICENSE` section with its first ported code. **Needs:** N8.

### Statistics

- [ ] **N11** `stats.mean`, `variance` and `stdDev` through the BLAS `sum`. **Needs:** N6.
- [ ] **N12** `stats.median` and the quantiles, all nine methods. **Needs:** N6.
- [ ] **N13** `stats.histogram` and `histogramInto`. **Needs:** N6.

### Sorting and sets

- [ ] **N14** `sort.unique`, `intersect`, `setUnion`, `setDiff` and `isMember`, each NaN distinct.
  **Needs:** N6.

### After the public layer

- [ ] **N1** N-d array stubs (design.md, "Planned: N-d arrays"). **Needs:** M3 (linalg's
  todo.md), whose `View(T)` and `Matrix(T)` they extend.
- [ ] **N2** Sparse matrix stubs (design.md, "Planned: sparse matrices"). **Needs:** M3.
- [ ] **N15** `stats.cov` and `corrcoef`: `syrk` over the centered columns. **Needs:** M3.

## Later

Picked up after the Queue, in this order.

- [ ] **N16** FFT: f64 SIMD (pffft's double) and a cache-blocked split for `n` past L2, for
  numerics. **Needs:** N10.
- [ ] **N17** FFT: SIMD across a batch's transforms, for images. **Needs:** N9.
- [ ] **N18** `fft.TrigPlan` (DCT and DST types 1 to 4) and `Mdct`. **Needs:** N8.
- [ ] **N19** `fft.StaticPlan` and `StaticRealPlan`. **Needs:** N8.
- [ ] **N20** FFT: route vDSP's sizes (`f 2^n`, `f` in 1, 3, 5, 15) to Accelerate on macOS,
  converting its real layout and scaling. **Needs:** N8.
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

## Done

- [X] The survey and API design ([`design.md`](design.md)), approved by the user 2026-10-31.

## Waiting on the user

- Review the stubs' changes to the approved sketches, now in design.md: explicit `comptime`
  element types; `stats.quantile`; `sort.Keep` on every set operation, `Selected(T)` and
  `[]?usize` results; `fft.Layouts` for out-of-place forms; the scratch sets a batch's
  concurrency.
- The Queue puts the FFT ahead of statistics and sets; reorder if that's wrong.
