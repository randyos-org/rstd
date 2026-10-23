# rstd BLAS design decisions against BLIS and industry references

Date: 2026-10-22. Research only; nothing in rstd or OpenBLAS was changed.

**Implemented** (2026-10-23, `wip/semantics`): 2 (`nrm2`: `f64` sums for
`s`/`c`, Blue's for `d`/`z`), 3 (first NaN wins in the extremes, NaN from
the value routines, the `nan_lanes` quirk gone), X1 (`sdot` in `f64`, one
`nrm2` and the extremes' rule on every core, and no zero-`alpha` return in
any `axpy` kernel) and pitfall 1 of 8 (x86_64's alignment-peeling
reductions routed to `@Vector` kernels). `readme.md`, "One contract on every
core", has the details. Also 10.2 and 10.3 (2026-10-23, `wip/denormals`):
no `f32` vector lanes where they flush (`util.vectorIeee`).

## Sources and versions

| Source | Where | Version |
| --- | --- | --- |
| rstd | this repository, branch `review` | `72b2efc` |
| rstd design doc | `design.md`, then on branch `wip/blas-docs` | `c33fca8` |
| rstd arm64 `@Vector` kernels | `vector.zig`, then on branch `wip/arm64-vector` | `4407e5a` |
| BLIS | <https://github.com/flame/blis> | tag `2.1`, `caf0db6be1202d9c83c79b51e75eceb96aa8b556`, 2026-06-25 |
| OpenBLAS | the checkout `tools/blas_asm` converts | `v0.3.34-317-ga69062524` |
| Reference BLAS (Reference-LAPACK) | <https://github.com/Reference-LAPACK/lapack> | last release v3.12.1 (2025-01-08); `master` as of 2026-10-22 |
| Exception-handling proposal | Demmel, Dongarra, Gates, Henry, Langou, Li, Luszczek, Pereira, Riedy, Rubio-Gonzalez, "Proposed Consistent Exception Handling for the BLAS and LAPACK", <https://arxiv.org/abs/2207.09281> (2022) | v1 |
| Intel oneMKL CNR | <https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2025-2/reproducibility-conditions.html> | 2025.2 docs |

BLIS and OpenBLAS paths below are relative to their source trees, rstd paths to
`src/linalg/blas/` of the branch named.
"Reference master" means unreleased commits on Reference-LAPACK `master`; they are not in 3.12.1.

Note on the review branch vs `wip/blas-docs`: `review`'s readme still says only real `gemm` is
ported; `wip/blas-docs` (`5fe7276`, `839628a`, `a6350dd`) has every level-3 routine, real and
complex. Level-3 statements below are about `wip/blas-docs`.

## Summary table

| # | Decision | Verdict | Priority |
| --- | --- | --- | --- |
| 2 | Portable `nrm2` algorithm (one-pass scaled loop) | **Change**: f64 sum for `s`/`c`, Blue's algorithm for `d`/`z` | High |
| 3 | `i?amax` NaN rule (NaN only if first) | **Change**: first NaN wins (reference master, BLIS, proposal) | High |
| X1 | Per-core differences in numerical contract (`sdot` precision, `nrm2`, SVE `axpy` zero skip) | **Change**: one contract per routine, enforced above the kernels | High |
| 10 | trsm pre-inverted diagonal; ARMv7 NEON / AltiVec flush to zero in `@Vector` tiles | **Investigate**, likely change | Medium-high |
| 5 | Rank updates skip zero columns (`ger`, `syr`, `her`, ...) | **Investigate**, lean to not skipping | Medium |
| 8 | Reproducibility rule | **Keep**, but close the pitfalls listed | Medium |
| 7 | Level 3: Goto/OpenBLAS driver, `level3_thread.c` protocol | **Keep driver; investigate** BLIS-style barrier packing for threading; fuse `beta` | Medium |
| 9 | `@Vector` kernels, unroll, accumulators, prefetch | **Keep**; try C-tile prefetch and Apple `KC` | Low-medium |
| 1 | `sdot` in f64 | **Keep** (more accurate than all references), measure cost | Low |
| 6 | alpha = 0 / beta = 0 semantics, `scal` flag | **Keep** | Low |
| 4 | Zero strides | **Keep** | Low |
| 11 | Named errors instead of `xerbla` | **Keep** | Low |
| X2 | trmm sums only the triangle (no packed zeros) | **Keep**; rstd is ahead of BLIS and OpenBLAS | n/a |

---

## 1. Reductions: `sdot`, `dsdot`, summation order

**rstd.** Portable `sdot` multiplies in f64 (exact for f32 operands) and sums sequentially in
f64 (`review:blas1/dot/__root__.zig` lines 69-120, `Dot(f32, f64, f32)`). The arm64 `@Vector`
`sdot` does the same with two vector accumulators of `4 * lanes` (`wip/arm64-vector:vector.zig`,
`Real.dot`, commit `4407e5a`). `dsdot` sums in f64 and returns f64. `sdsdot` adds `alpha` after
the f64 sum, then rounds.

**OpenBLAS.** Inconsistent across kernels:

- `kernel/arm/dot.c` (the generic C): `double dot`, but the product `y[iy] * x[ix]` is computed in
  `FLOAT` (f32) before it is added, unless `DSDOT`.
- `kernel/x86_64/sdot.c:51`: `FLOAT dot = 0.0` (f32 accumulation, with f32 vector microkernels).
- `kernel/arm64/dot_kernel_asimd.c:34-46`: `RETURN_TYPE` is `double` only for `DSDOT`; `sdot`
  accumulates in f32 NEON lanes.

**BLIS.** `ref_kernels/1/bli_dotv_ref.c` accumulates in `ctype` (f32 for `s`).
`kernels/zen/1/bli_dotv_zen_int10.c:77-105` keeps ten `__m256` f32 accumulators.
`dsdot`/`sdsdot` are plain scalar double loops in `frame/compat/bla_dot.c:154-211`.

**Reference.** `sdot.f` accumulates in `REAL`; `dsdot.f`/`sdsdot.f` in `DOUBLE PRECISION`
(`sdsdot` seeds the sum with `SB`, then rounds once).

**Assessment.** Accumulating `sdot` in f64 is more accurate than every reference: products are
exact and the f64 sum's error, `n * 2^-53 * sum|x*y|`, is far below one f32 ulp of the result
unless there is heavy cancellation, so the result is close to correctly rounded. It also makes
results nearly independent of summation order, which helps the reproducibility rule. No
reference requires f32 accumulation; no caller can be hurt by extra accuracy. The cost is real
in cache: widening halves lanes and adds `fcvtl`/`cvtps2pd` work. The M6 table in
`tools/blasbench/results-m6.md` (`dot`, `s`, 218 GB/s at 4 KiB for `vec x4`) predates `4407e5a`,
so the f64 `sdot` is unmeasured. Out of cache it should be close to bandwidth either way
(uncertain; not measured).

The real problem is consistency (see X1): x86_64, arm (VFP), POWER and others keep OpenBLAS's
f32 `sdot`, so the accuracy of `sdot` depends on the core, by far more than "the last bits".

**Recommendation: keep** f64 accumulation, and measure it on M6 and Zen. If the in-cache cost
matters, Ogita-Rump-Oishi `Dot2` (error-free `TwoProduct` with `@mulAdd`, compensated sum)
gives similar accuracy in f32 lanes; measure before choosing. Summation order (sequential in
the portable kernel, two vector accumulators in `@Vector`): keep; no reference prescribes an
order and none uses pairwise summation. `dsdot`/`sdsdot`: keep. Adding `alpha` after rather than
before the f64 sum differs from `sdsdot.f` only in f64 rounding; harmless.

## 2. `nrm2`

**rstd.**

- Portable kernel (`review:blas1/nrm2/__root__.zig:54-95`): the classic one-pass scaled loop
  (`scale`, `ssq`, two divisions per nonzero element), in `types.Acc(F)`, which is `F` for f32
  and f64. If the result is NaN it rescans for a NaN and otherwise returns `+Inf`. So the
  *semantics* match LAPACK 3.10 (Inf for any Inf, NaN for NaN), but the *algorithm* is the
  pre-3.10 one (old `dnrm2.f`, LAPACK's old `dlassq`), not Blue's.
- arm64 `@Vector` (`wip/arm64-vector:vector.zig`, `squares`, `root`, `scaledNrm2`): plain sum of
  squares in `F`, accepted when `floatMin/eps * n <= s <= floatMax`, else the scaled loop. The
  acceptance test is sound: squares lost to underflow total at most about `n * floatMin`, which
  the lower bound keeps below `eps * s`.
- Interface: `incx == 0` gives `sqrt(n) * |x[0]|` (`blas1/interface.zig`, `nrm2`), as
  `interface/nrm2.c` does.
- Converted kernels: x86_64 keeps OpenBLAS's `nrm2_sse.S` (f32 widened to f64, unscaled) and
  `nrm2.S` (x87 f80, unscaled). POWER, Loongson, loongarch64 and SPARC fall back to the portable
  loop because their kernels mishandle Inf, NaN and denormals (`issues.md`).

**Reference (LAPACK 3.10+).** `BLAS/SRC/dnrm2.f90`
(<https://github.com/Reference-LAPACK/lapack/blob/master/BLAS/SRC/dnrm2.f90>, Anderson,
"Algorithm 978: Safe Scaling in the Level 1 BLAS", ACM TOMS 44(1), 2017, citing Blue, ACM TOMS
4, 1978). One pass, three accumulators (`asml`, `amed`, `abig`) selected by thresholds
`tsml = 2^-511`, `tbig = 2^486` (double), scaled by constants `ssml = 2^537`, `sbig = 2^-538`;
no divisions in the loop. Inf lands in `abig` and gives Inf; NaN propagates; `incx = 0`
counts `x(1)` `n` times; `incx < 0` walks backwards.

**OpenBLAS.** Per-arch. x86_64 widens (`kernel/x86_64/nrm2_sse.S:72,99`, `cvtss2sd`/`cvtps2pd`;
`kernel/x86_64/KERNEL:257-274` selects it for `snrm2`, x87 `nrm2.S` for `dnrm2`). arm64
ThunderX2 `scnrm2_thunderx2t99.c` also widens with `fcvt`. The others are scaled loops, several
broken on special values (rstd's `issues.md`).

**BLIS.** `frame/util/bli_util_unb_var1.c`: `normfv` calls `sumsqv_unb_var1` (lines
1060-1173), a scalar one-pass scaled loop "mimicking LAPACK's `?lassq()`", with explicit NaN and
Inf handling (NaN sticks; Inf sets `sumsq = Inf, scale = 1`). A faster variant (unscaled `dotv`,
then fall back if `FE_OVERFLOW`/`FE_INVALID` were raised) exists but is disabled under `#if 0`
(lines 312-400, GitHub issue #332). No optimized `normfv` kernel exists in `kernels/`. BLIS's
`nrm2` is therefore accurate on special values but slow.

**Assessment.** rstd's special-value semantics are right. Its portable algorithm is the weakest
of the available choices:

- Speed: two divisions per element and a data-dependent branch; unvectorizable as written. Blue's
  algorithm has no division and vectorizes with three accumulators; the plain sum is faster
  still.
- Accuracy: each step rounds a division, a square and a rescale, so the error bound is worse
  than a plain or Blue's sum (qualitative; not measured here).
- For `s` and `c` there is a strictly better option that OpenBLAS already uses on x86_64: sum the
  squares in f64, unscaled. An f32 square is exact in f64; the largest (`2^256`) and the smallest
  denormal (`2^-298`) squared are inside f64's normal range, so no scaling is needed for any `n`
  below about `2^767`. One rounding at the end; the result is nearly correctly rounded.

**Recommendation: change.**

1. `snrm2`/`scnrm2` (portable and arm64 `@Vector`): f64 accumulation of exact squares, no
   scaling, no fallback (NaN and Inf propagate naturally: Inf squared is Inf, and two Infs give
   Inf, not NaN).
2. `dnrm2`/`dznrm2` portable: Blue's algorithm as in `dnrm2.f90`, which needs no NaN rescan.
3. Keep the arm64 fast path for `d`/`z`, with Blue's algorithm as its fallback instead of the
   scaled loop.
4. Fix the wording in `readme.md`/`design.md`: "LAPACK 3.10's" applies to the results, not the
   algorithm.

Threading note: Blue's partials are a triple (`asml`, `amed`, `abig`) and combine in chunk order
as easily as `(scale, ssq)`.

## 3. `i?amax` / `i?amin` with NaN

**rstd.** Portable (`review:blas1/extremum.zig`): start at element 0, replace on strict `>`; a
NaN is chosen only if it is first, and then is never replaced. The arm64 `@Vector` `firstMax`
copies that rule. The `?amax` *value* routines follow the same scan, so a NaN after the first
element is ignored. Several converted x86_64, POWER, loongarch64 and arm64 kernels are left as
`nan_lanes` quirks (`blas1/edge.zig`), so the rule is not uniform across cores. Complex keys
`|re| + |im|` can overflow to Inf for finite elements; nothing handles that.

**Reference 3.12.1 (released).** `idamax.f` replaced on `DABS(DX(I)).GT.DMAX`: rstd's rule.

**Reference master (unreleased).**

- `izamax.f90`/`icamax.f90` rewritten for the exception-handling proposal (`50dc1bf`,
  2025-03-27): "If the vector contains a NaN, the index of the first NaN is returned. Otherwise,
  if it contains an Inf or -Inf, the index of the first of those is returned." Finite elements
  whose `|re| + |im|` overflows switch the comparison to `0.25*|re| + 0.25*|im|`.
- `idamax.f`/`isamax.f` (`b358788`, 2026-09-14): "make i?amax return the first NaN"; the test is
  now `.NOT.(DABS(DX(I)).LE.DMAX)` plus `DMAX.NE.DMAX`.

**BLIS.** `ref_kernels/1/bli_amaxv_ref.c`: replace if `abs_chi1_max < abs_chi1 ||
(isnan(abs_chi1) && !isnan(abs_chi1_max))`, so the first NaN wins (the comment there says
"smaller than any previously seen", but the code selects the NaN). The Zen kernel
`kernels/zen/1/bli_amaxv_zen_int.c:70-92` implements the same rule in vector lanes. The value
form `normiv` (`frame/util/bli_util_unb_var1.c:449-486`) propagates NaN "to mimic LAPACK's
`?lange()`".

**Proposal (arXiv 2207.09281).** Real: first NaN, else first Inf, else first largest. Complex:
the overflow-safe version in its Appendix A.

**Assessment.** rstd's rule matches only the old reference and OpenBLAS's generic C. The
direction of every maintained source (reference master for all four precisions, BLIS, the
proposal) is "first NaN". It is also cheaper to make order-independent, which matters for
threading: "first NaN, else first max" combines exactly over chunks, as does rstd's rule, but it
removes the position dependence that the reference itself called a bug.

**Recommendation: change** to: first NaN; else (complex) first Inf part; else first largest key,
with `izamax.f90`'s quarter scaling when a finite key overflows. Make `?amax`/`?amin` values
return NaN if any element is NaN (BLIS `normiv`, LAPACK `lange`). For `i?amin`, the mirror rule
(first NaN) is what reference master's new docs imply for the max; `i?amin` is an OpenBLAS
extension with no reference, so choose and document. The `nan_lanes` quirk kernels would then
need replacing (the arm64 `@Vector` `i?amax` matched or beat the converted kernel on M6;
x86_64 unmeasured) or a NaN pre-check.

## 4. Zero strides (`incx`/`incy` = 0)

**rstd.**

- `axpy`: both zero gives `y[0] += n*alpha*x[0]` (in `Acc`), as `interface/axpy.c:80-81`;
  otherwise the kernels run sequentially (the converted `swap`/`rot` kernels that batched a block
  before storing were patched to go one element at a time; `readme.md`, "Fixes").
- `scal`, `asum`, `i?amax`: nothing (or 0) for `incx <= 0`, as the reference.
- `nrm2`: `sqrt(n) * |x[0]|`.
- Level 2 and rank updates: `incx == 0` is `InvalidIncx`, as the reference's `xerbla`.

**Reference.** Level-1 loops simply step by `incx`, so a zero input stride broadcasts `x(1)`
and a zero output stride applies `n` sequential updates to one element. `dscal`, `dasum`,
`idamax` return early for `incx <= 0`. Level 2 rejects `incx = 0`.

**BLIS.** The compat layer (`frame/include/bli_blas_macro_defs.h:50-75`) only flips negative
strides; zero passes through to kernels, which loop sequentially. `frame/compat/bla_scal.c` has no
`incx <= 0` guard, so BLIS `dscal` with a negative stride scales the vector and with zero
scales `x(1)` `n` times: a deviation from the reference.

**MKL** documents a zero increment on an input vector as all elements equal to `x[0]`
(<https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-2/vector-arguments-in-blas.html>).
A zero stride on an *output* vector is not specified by any standard I found; the reference's
sequential semantics is the de facto definition (uncertain: I could not locate the BLAS
Technical Forum text on this).

**Assessment.** rstd follows the reference's sequential meaning everywhere but `axpy` with both
strides zero, where it follows OpenBLAS's closed form. That closed form rounds once instead of `n`
times and can overflow (`n*alpha*x`) where the sequential sum would too, so it is a defensible,
documented deviation.

**Recommendation: keep.** Document the `axpy` closed form as an intentional OpenBLAS deviation.
Do not copy BLIS's `scal` behavior.

## 5. Rank updates skipping zero columns

**rstd.** `ger` splits `y` into runs of nonzeros and calls the kernel once per run
(`review:blas2/interface.zig`, `gerImpl`); the `syr`, `her`, `spr`, `hpr` and rank-2 drivers skip
`x(j) == 0` (`blas2/banded.zig:350,379`), as reference BLAS does. OpenBLAS's `zher_k.c`,
`zhpr_k.c` and rank-2 drivers do not skip (rstd `issues.md`, "Rank updates and a zero factor").

**Reference.** 3.12.1 and master `dger.f`/`dsyr.f` still test `Y(JY).NE.ZERO` / `X(J).NE.ZERO`
(no recent change to `dger.f`, `dsyr.f`). But master removed the analogous checks from the
triangular routines in July 2026: `dtrsv.f` (`6e7b86c`, "remove the IF (X(J).NE.ZERO) THEN that
saves FLOPS but might prevent NaN propagation"), `dtrmv.f`, `dtrmm.f`, and `dtrsm.f` (`38ccf80`,
"remove all IF (A(K,J).NE.ZERO)"). `dgemm.f` never skips zero elements of `B`.

**Proposal (arXiv 2207.09281).** "Would allow checking for alpha=0, but not checking for zeros
inside" `x` or `y` for `GER`; calls the `SYR` case a worse inconsistency (upper and lower can
differ).

**BLIS.** `ger`/`her`/`syr` are built on `axpyv`/`axpy2v`, whose reference kernel returns
for `alpha == 0` (`ref_kernels/1/bli_axpyv_ref.c:57`), so BLIS skips a column whose scaled factor
is exactly zero, by accident of composition rather than by contract.

**Assessment.** rstd matches the letter of today's reference for `ger`/`syr`, but the field is
moving the other way (reference master for the triangular family, the proposal), and skipping
has two costs in rstd:

- Performance pathology in `ger`: a `y` with many zeros turns one kernel call into many; with a
  strided `x` each call re-copies `x` into the scratch buffer (`blas2/buffers.md`), so the cost
  can approach doubling.
- It is NaN-inconsistent with `gemm` (a rank-1 `gemm` propagates the NaN that `ger` hides).

**Recommendation: investigate, leaning to change**: drop interior zero skipping in `ger` and the
rank-update drivers (keep the `alpha == 0` quick return), so that `ger` agrees with `gemm` and
with the OpenBLAS kernels it calls, which also removes the run splitting. If strict
reference-3.12 parity is wanted instead, keep it but document the `ger` pathology.

## 6. NaN/Inf propagation for alpha = 0 and beta = 0; the `scal` flag

**rstd.** `gemm`/`gemv`/level 3: `beta == 0` never reads `C`/`y`; `alpha == 0` (or `k == 0`)
returns after scaling by `beta` without touching `A`/`B` (`wip/blas-docs:blas3/level3.zig:379-381`).
Public `scal` passes flag 1: a zero `alpha` turns NaN and Inf into NaN, finite into 0
(`review:blas1/interface.zig`, `scal`); internal `beta` scaling passes flag 0 (store zeros).
`axpy` returns for `alpha == 0`.

**Reference.** `dgemm.f` now documents it (master, 2025-01-19): "If alpha is zero, NaNs or Infs
in A or B do not affect the result"; "If both alpha and beta are zero, then a zero matrix is
returned in C"; "If only beta is zero, alpha*op(A)*op(B) is returned". `dscal.f` multiplies
(`DA*DX(I)`), so `0*NaN = NaN`, and returns early only for `DA == 1`. `daxpy.f` returns for
`DA == 0`.

**OpenBLAS.** The flag argument (`interface/scal.c:98` passes 1; `kernel/arm/scal.c:46` tests
`dummy2 == 0`) implements exactly rstd's split.

**BLIS.** `scalv` with `alpha == 0` calls `setv` (`ref_kernels/1/bli_scalv_ref.c:57-66`), so BLIS
`dscal(0, x)` stores zeros over NaN: a deviation from the reference. `gemm` microkernel
overwrites `C` for `beta == 0` (`ref_kernels/3/bli_gemm_ref.c:124-126`). `axpyv` returns for
`alpha == 0`.

**Proposal.** Agrees with "alpha = 0 does not propagate from A/B" and "beta = 0 does not
propagate from C"; it notes `SCAL` does not test `alpha` and proposes no change.

**Recommendation: keep.** rstd matches the reference's documented contract in every case and is
more faithful than BLIS on `scal`.

## 7. Level 3

**rstd** (`wip/blas-docs`). `blas3/level3.zig` is `driver/level3/level3.c`: `R`-wide column blocks
(BLIS `NC`), `Q`-deep slices (`KC`), `P`-tall row blocks (`MC`), packed `B` once per slice,
packed `A` once per block, register tile (`MR x NR`) from `geometryFor` in
`blas3/gemm/kernel.zig`. One driver serves `gemm`, `symm`/`hemm` (pack the symmetric operand
densified), `syrk`/`herk`/`syr2k`/`her2k` (diagonal tiles through a scratch tile). `trsm`/`trmm`
follow `trsm_L.c`, `trsm_R.c`, `trmm_L.c`, `trmm_R.c` with Zig ports of the generic
`trsm_kernel_*.c` (gemm update, then `solve` per tile, pre-inverted diagonal). `beta` is a
separate pass over `C` before the loops (`level3.zig:380`). Complex: native `@Vector` kernel, two
accumulators per element (the `zgemm_kernel_4x2_haswell.S` scheme). Planned threading:
`level3_thread.c`'s shared packed `B` protocol on `Io.Group.concurrent` (`threading.md`).

**BLIS.**

- Five loops around the microkernel, the same Goto structure (`docs/Multithreading.md:29`;
  Smith et al., "Anatomy of High-Performance Many-Threaded Matrix Multiplication", IPDPS 2014,
  <https://doi.org/10.1109/IPDPS.2014.110>).
- Microkernel contract: `C := beta*C + alpha*A*B` on one `MR x NR` tile with general strides,
  `beta` fused (`ref_kernels/3/bli_gemm_ref.c`); `beta` is applied on the first `KC` slice only,
  so there is no separate pass over `C`.
- Threading by "ways" per loop: `JC` (n), `IC` (m), `JR`, `IR`; the `PC` (k) loop is "Unavailable;
  always 1" (`docs/Multithreading.md:160-172`). Threads in a group pack shared `A` and `B` blocks
  cooperatively and synchronize with barriers.
- Complex: native kernels where a configuration has them (Haswell `c`/`z`), else the **1m**
  method, which reuses the real microkernel by reorganizing packed panels (Van Zee, "Implementing
  High-Performance Complex Matrix Multiplication via the 1m Method", SIAM J. Sci. Comput. 42(5),
  C221-C244, 2020, <https://doi.org/10.1137/19M1282040>). 3m/4m were studied earlier (Van Zee and
  Smith, ACM TOMS 44(1), 2017, <https://doi.org/10.1145/3086466>); 1m superseded them in BLIS.
- `trsm`: fused `gemmtrsm` microkernel per tile (`kernels/haswell/3/bli_gemmtrsm_l_haswell_asm_d6x8.c`),
  pre-inverted diagonal by default with `--disable-trsm-preinversion` to divide instead
  (`configure:375-387`); only three of the five loops can be parallelized (`docs/Multithreading.md:29`).
- `herk`/`syrk`: `gemmt` with pruned off-triangle tiles; `her2k`/`syr2k` call `gemmt` twice, `beta`
  only the first time (`frame/3/bli_l3_oapi_ex.c:253-255`).
- Triangular packing writes explicit zeros in the unstored part (`frame/1m/packm/bli_packm_struc_cxk.c:180-190`).
- Small and skinny problems take the "sup" path (unpacked or partly packed kernels, thresholds
  `MT`/`NT`/`KT`, e.g. `config/zen3/bli_cntx_init_zen3.c:244-253`).

**Assessment, item by item.**

- *Driver structure.* rstd's driver and BLIS's five loops are the same algorithm; rstd's single
  driver for every packed routine mirrors BLIS's "everything through gemm/gemmt". **Keep.**
- *Separate `beta` pass.* BLIS fuses `beta` into the first slice's kernel call; OpenBLAS (and
  rstd) read and write all of `C` once more. `gemm/notes.md` attributes 7-25% of the full-gemm
  loss to packing, the `C` update and edges; for LAPACK-style updates with small `k` the extra
  pass is a larger share. **Change (performance)**: pass `beta` to the kernel for the first `Q`
  slice and 1 afterwards; keep the separate pass only for `alpha == 0`/`k == 0`. This also keeps
  `beta == 0` from reading `C`.
- *Threading protocol.* OpenBLAS's flag-based producer/consumer (`level3_thread.c`) has a
  history of races: `Changelog.txt:1703-1715` (0.3.9: "Fixed a race condition in the level3 blas
  code" on POWER, ARMV7, ARMV8), `:250-254` (0.3.31: reverted a 0.3.30 partitioning change "that
  could lead to race conditions and subsequent invalid results in GEMM"), `:7` (0.3.34),
  `:1955-1961` (`USE_SIMPLE_THREADED_LEVEL3` as the known-safe fallback). BLIS's cooperative
  packing plus barriers needs the same "all workers concurrent" guarantee but has a much smaller
  state space, and its `JR`/`IR` ways let threads that share an L2 (Apple clusters, SMT pairs)
  share the packed `A` block too. Neither splits `k`, so both meet rstd's reproducibility rule.
  **Investigate**: prototype BLIS-style `JC x IC` ways with cooperative packing and a futex
  barrier before committing to the flag protocol; rstd's start-barrier design already covers the
  "all or nothing" part.
- *Complex.* rstd's native `@Vector` complex kernel is generated for every core, so the main
  benefit of 1m (no hand-written complex kernels) is already had. 1m could help only where the
  complex tile is weak, and on the `split` (scalar) targets the real kernel is scalar too.
  **Keep native; no 1m.** For the planned `gemm3m`/`symm3m`/`hemm3m`: 3m's imaginary part is only
  normwise stable (Higham, SIAM J. Matrix Anal. Appl. 13(3), 1992), so keep it opt-in, as
  OpenBLAS does (separate entry points), never the default.
- *trsm.* rstd's per-tile "gemm update then solve" is the same as BLIS's fused `gemmtrsm`.
  **Keep**, except the pre-inversion (see 10).
- *herk/syrk.* Scratch tile for diagonal tiles equals BLIS `gemmt`'s approach. **Keep.**
- *trmm.* See X2: rstd is better than both.
- *Small/skinny.* rstd has no small-matrix path (OpenBLAS's `SMALL_MATRIX_OPT` kernels are not
  ported either). BLIS's sup path is its biggest practical win for small `m, n` or small `k`.
  **Investigate later**, after threading.

## 8. Threading and reproducibility

**rstd** (`threading.md`). Bit-identical results for any worker count, `Io`, scheduling, and
threaded vs single-threaded builds, for a given core selection. Levels 1 and 2: chunks fixed by
problem size (at most 64), claimed dynamically, partials combined in chunk order. Level 3: each
`C` element owned by one worker, `k` blocking independent of the split.

**References.**

- **MKL CNR**: run-to-run reproducibility requires a fixed thread count (`MKL_DYNAMIC=FALSE`) and a
  pinned code branch (`MKL_CBWR`); only "strict CNR" (`?gemm`, `?symm`, `?hemm`, `?trsm`, AVX2 or
  AVX-512, 64-bit) is bitwise reproducible "even when the number of threads changes", and
  "enabling strict CNR mode can reduce performance" (developer guide, reproducibility conditions,
  URL above). Older MKL versions also required aligned data for CNR (from memory; the current page
  says CNR "fully supports unaligned input and output data").
- **BLIS**: level 1 and 2 are not threaded; level 3 never splits `k` (`docs/Multithreading.md:161`),
  so thread-count invariance holds by construction, though the docs never state it as a guarantee
  (no "reproducib" in `docs/`).
- **OpenBLAS**: level-1 reductions split by thread count (`kernel/arm64/dot.c:161-199`,
  `kernel/x86_64/ddot.c:160-168`), so bits change with `OPENBLAS_NUM_THREADS`.
- **Literature**: two families. Fixed reduction trees (MKL CNR; oneTBB's
  `parallel_deterministic_reduce` with a fixed grain size, from memory) are cheap but reproducible
  only for the same tree. Order-independent summation (ReproBLAS: Ahrens, Demmel and Nguyen,
  "Algorithms for Efficient Reproducible Floating Point Summation", ACM TOMS 46(3), 2020, from
  memory) is reproducible across trees and machines at a higher cost.

**Assessment.** rstd's rule is stronger than MKL's default CNR and equal to strict CNR's
guarantee, without a performance mode switch, which is plausible precisely because the chunks
depend only on `n`. The approach (fixed tree) is the standard one; ReproBLAS is unnecessary for
the stated goal. Pitfalls found:

1. **Alignment-dependent converted kernels.** `kernel/x86_64/asum_sse2.S:75`, `dot_sse2.S:82,93`,
   `nrm2_sse.S:68`, `zdot_sse2.S:97-100` peel one element when `x` is 8- but not 16-byte aligned,
   which changes lane assignment and summation order. Same data in a differently aligned buffer
   gives different bits on cores that use them. This does not break thread-count invariance
   (chunks inherit the base alignment) but breaks the natural reading of "reproducible". State the
   condition, or make chunk sizes multiples of 64 bytes and document "for a given alignment", or
   replace these kernels with `@Vector` ones that do not peel.
2. **Single-threaded builds.** "Threading compiles out under `single_threaded`" must still chunk
   above the threshold, or single-threaded and threaded builds differ. Say so explicitly.
3. **Inputs to chunk sizes.** Chunk sizes and thresholds must not read the CPU count,
   `OPENBLAS_NUM_THREADS` or `Io` limits. Planned per-core tables are fine.
4. **`Q` must not depend on threads.** `blas3/gemm/params.zig` picks LA464's set by
   `env.threads == 1` (`la464_sets`, lines 212-224); today only `P` and `R` differ between the
   single- and multi-thread sets, and `Q` differs only by L3 size, so the rule holds. Add a test
   that `Q` is identical for every thread count, since a future set could break it silently.
5. **Edge tiles under repartitioning.** If thread partitions move tile boundaries, an element can
   move between the full tile and `thin`/tail paths. Rounding `M` ranges to `unroll_m` (planned)
   and `N` ranges to `unroll_n` avoids it; otherwise every path must do identical per-element
   operations (same `k` order, same fusion). The planned 1-9 worker bitwise tests will catch it if
   the shapes sweep remainders.
6. **FMA per core.** One core selection serves the whole process, so all threads run the same
   kernels; heterogeneous cores are not a problem unless dispatch ever becomes per thread.
   But `fma.zig`'s policy for the *portable* kernels follows the build target, so the same core
   selection can give different bits from builds with different `-mcpu` wherever a portable
   kernel serves. Phrase the guarantee as "for a given build and core selection" (or pin the
   portable policy).

**Recommendation: keep the rule**, and fix or document the six pitfalls. Also see X1: per-core
differences in `sdot` precision exceed "last bits".

## 9. Kernel strategy

**rstd.** One `@Vector` source per routine compiled per core group; level 3 tile from a register
budget (`geometryFor`); arm64 level 1/2 `@Vector` kernels with `U = 4` native vectors per step and
two accumulators for reductions (8 independent chains on NEON); no software prefetch.

**BLIS.** Reference kernels are portable C with `PRAGMA_SIMD`, compiled per sub-configuration
(`ref_kernels/`), and are not expected to be fast; performance comes from hand-written per-arch
inline-assembly microkernels (`kernels/<arch>/3/`). Level-1 Zen kernels unroll ten vectors with
unaligned loads and no prefetch (`kernels/zen/1/bli_axpyv_zen_int10.c:72-101`,
`bli_dotv_zen_int10.c`). gemm microkernels prefetch the `C` tile before the `k` loop and `A`
ahead (`kernels/haswell/3/bli_gemm_haswell_asm_d6x8.c:127-161, 175, 227`; ARMv8
`kernels/armv8a/3/bli_gemm_armv8a_asm_d6x8.c:41-47`, a prefetch fix for unaligned `C` columns).

Tile shapes agree with BLIS where comparable:

| Target | rstd (`notes.md`) | BLIS (`config/*/bli_cntx_init_*.c`) |
| --- | --- | --- |
| AVX2 `s` / `d` | 16 x 6 / 8 x 6 | Haswell 6 x 16 / 6 x 8, row-preferential (same registers) |
| AVX-512 `s` / `d` | 32 x 14 / 16 x 14 | SKX 32 x 12 / 16 x 14 |
| NEON `s` / `d` | 16 x 6 / 8 x 6 | Firestorm 12 x 8 / 8 x 6 |

**Assessment.** rstd's approach is BLIS's reference-kernel idea done with explicit vector types
and a register model, and the measurements (`gemm/notes.md`: about 100% of measured peak on
Zen 3; `results-m6.md`: level 1 at or above OpenBLAS's arm64 assembly) support it. Accumulator
counts are adequate: 8 chains cover Zen's 2 FMA pipes x 4-cycle latency, and on Apple cores dot
is load-bound before it is FMA-bound. BLIS's lack of level-1 prefetch agrees with rstd's finding
that OpenBLAS's prefetch-heavy arm64 kernels are slower on M6.

**Recommendation: keep.** Two experiments are worth their small cost:

1. `@prefetch` the `C` tile (one line per column) at tile entry, as BLIS does; `notes.md`
   attributes part of the gemm loss to the `C` update.
2. Blocking on Apple cores: BLIS Firestorm uses `KC` 4096 (`s`) and 3072 (`d`),
   `MC` 480/256 (`config/firestorm/bli_cntx_init_firestorm.c:96-100`), against OpenBLAS's VORTEX
   `Q` 1024/512, `P` 512/256 (`param.h`, VORTEX block). Sweep `Q` upward on M6.

Also note that OpenBLAS's `param.h` VORTEX block chooses `P`/`Q` by `NUM_CORES > 8` at build
time; check which branch `params_data.zig` captured.

## 10. Special values and denormals

**rstd.** The library never changes MXCSR/FPCR; tests require that denormals do not flush
(`blas1/edge.zig`; `check.Pattern.denormal` in level 3). Level-2 `trsv` divides by the diagonal for
real types and multiplies by a reciprocal for complex (`blas2/triangular.zig:153-163`); level-3
`trsm` packs the reciprocal of the diagonal (`blas3/trsm/kernel.zig:46-72`, OpenBLAS's
`compinv` for complex).

**References.**

- Neither OpenBLAS nor BLIS sets FTZ/DAZ; MKL also leaves it to the application (from memory;
  Intel compilers set FTZ in `main` at some optimization levels).
- BLIS `configure:375-387`: pre-inversion is the default, and disabling it improves "numerical
  robustness ... for certain cases involving denormal numbers that would otherwise result in
  overflow in the pre-inverted values".
- Reference master `dtrsm.f` (`e8785e4`, July 2026): "remove all TEMP = ONE/A(J,J) and force
  explicit division".

**Findings.**

1. **Pre-inverted `trsm` diagonal.** `1/d` overflows to Inf for `|d| < 2^-1024` (f64) or
   `2^-128` (f32), and loses accuracy where `1/d` is denormal, so `trsm` can give Inf where the
   reference gives a finite result, and `trsv` and `trsm` disagree on the same system. Each
   multiply-by-reciprocal also adds a rounding.
   **Recommendation: change, or at least guard**: divide in `solve` (the cost is one division per
   row of each diagonal block, `O(n^2)` against `O(n^3)` work), or keep the reciprocal only when
   it is finite and normal and fall back to division for that row. Measure first; BLIS's own note
   calls the penalty a performance cost, not a prohibitive one.
2. **ARMv7 NEON flushes denormals.** The ARMv7 Advanced SIMD unit always runs with the "standard
   FPSCR value" (flush-to-zero and default NaN), regardless of FPSCR (ARMv7-A Architecture Reference
   Manual, from memory but well established). Zig emits NEON for `@Vector(4, f32)` there: compiling
   a two-line test for `arm-linux-gnueabihf -mcpu cortex_a15 -OReleaseFast` gives `vmul.f32 q8, q9,
   q8` and `vfma.f32 q8, q10, q9`. `vectorBits` in `blas3/gemm/kernel.zig` returns 128 for `f32`
   with NEON on arm, so `sgemm`/`cgemm` and every level-3 routine on ARMv7 flush denormal inputs and
   outputs, unlike the scalar VFP level-1/2 kernels. The level-3 `denormal` pattern should fail
   there under the tolerance in `blas1/check.zig:317`; I did not run the arm rig, and whether QEMU
   applies FZ to AArch32 NEON is unverified. **Investigate**: run the level-3 denormal case on the
   `arm` rig and on hardware; if confirmed, either use scalar VFP for ARMv7 `f32` tiles (as
   OpenBLAS's ARMv7 `sgemm` VFP kernels do) or document flushing as an ARMv7 property.
   **Done**: confirmed on the `arm` and `linux-arm` rigs (the portable `saxpy`'s vectors too, not
   `arm-softfp`, whose vectors are soft-float calls); `f32` tiles and kernels are scalar VFP there.
3. **AltiVec.** Pre-VSX PowerPC (`ppcg4`, `ppc970`) vector float uses AltiVec, which flushes
   denormals in non-Java mode (`VSCR[NJ] = 1`); Linux starts processes with `NJ = 1` (from memory
   of `arch/powerpc/kernel/process.c`; unverified here). VSX (POWER7+) is IEEE. **Investigate**
   alongside item 2 on the `powerpc` rig.
   **Done**: confirmed on the `powerpc` and `powerpc-pseries` rigs; scalar FPU there, as item 2.
   Clearing `NJ` (Java mode, IEEE denormals) was not taken: it changes the caller's environment,
   per thread, and the 970 handles such denormals by a software assist.

**Policy recommendation.** Write down "rstd never changes the floating-point environment; results
assume IEEE default modes; denormals are honored except where the hardware's vector unit cannot
(ARMv7 NEON, AltiVec), listed here". That matches every reference and makes items 2 and 3 explicit.

## 11. Error handling

**rstd.** Named errors in the reference's `xerbla` argument order; no `xerbla`; out-of-memory as
an error; nothing panics.

**References.** Reference BLAS calls `XERBLA` (prints, stops). OpenBLAS calls `xerbla` and 0.3.34
added an extension to override it (`Changelog.txt:11`). BLIS's native API checks arguments when
the thread-local error-checking level is `BLIS_FULL_ERROR_CHECKING` (default) and then aborts
(`frame/base/bli_error.c:134-163`); its BLAS compat `xerbla` prints and returns
(`frame/compat/f2c/bla_xerbla.c:83`, abort commented out).

**Assessment.** Returning an error union is the idiomatic Zig equivalent and safer than all three
(no process abort, no silent return). The checks are O(1), so BLIS's "turn checking off" level
has no counterpart worth adding. **Keep.**

---

## X1. Cross-core consistency of the numerical contract

Found while comparing; not one of the listed decisions. rstd selects one core per process and
promises reproducibility per core, but some routines' *semantics* (not just rounding) depend on
which core's kernel is installed:

| Routine | Cores | Behavior |
| --- | --- | --- |
| `sdot` | portable, arm64 `@Vector` | f64 accumulation |
| `sdot` | x86_64 (`sdot.c`), arm (`sdot_vfp.S`), POWER (`sdot.c`), ... | f32 accumulation |
| `snrm2` | x86_64 (`nrm2_sse.S`) | f64 unscaled |
| `snrm2` | portable, arm64 `@Vector` | f32 scaled / plain f32 |
| `dnrm2` | x86_64 (`nrm2.S`) | x87 f80 unscaled |
| `dnrm2` | portable | f64 one-pass scaled |
| `axpy` with `alpha == 0` | SVE cores (`blas1/axpy/arm64/sve.zig:52`, from `kernel/arm64/axpy_sve.c:57`) | returns without reading `x` |
| `axpy` with `alpha == 0` | other cores | multiplies (`0 * Inf = NaN`) |
| `i?amax` with NaN | `nan_lanes` quirk cores | depends on lane |

The `axpy` row matters because the level-2 triangular and banded drivers call `axpy_k` with
`-x(j)` (OpenBLAS's `driver/level2` design), so `trsv`/`tbsv` NaN propagation from a zero
solution component differs between SVE and other arm64 cores.

**Recommendation: change.** Pick one contract per routine and enforce it above the kernels:
`sdot` and `snrm2` in f64 everywhere (replace or wrap the f32 kernels; the arm64 `@Vector`
kernels show the replacements can be fast), Blue's or wide accumulation for `dnrm2`, and a
driver-level rule for zero multipliers (call or skip, decided in the driver, not the kernel). The
`threading.md` sentence "across core types, results may differ in the last bits" is not true for
these routines today.

## X2. trmm without packed zeros (rstd ahead)

rstd's `trmm` kernels sum a diagonal block over its triangle only, so an Inf in `B` never meets a
structural zero (`design.md`, Level 3). BLIS writes zeros into the unstored part of a packed
triangular block (`frame/1m/packm/bli_packm_struc_cxk.c:180-190`) and OpenBLAS's tiles multiply
packed zeros, so both produce NaN there; reference master's `dtrmm.f` changes (July 2026) remove
zero *skipping*, but the reference never multiplies by the unstored triangle. **Keep.**

---

## Where BLIS is notably better, per kernel

Worth porting or emulating:

1. **`beta` fused into the gemm microkernel's first `KC` slice** (`ref_kernels/3/bli_gemm_ref.c`):
   removes a full pass over `C`; matters most for small `k`.
2. **Cooperative packing with barriers and per-loop "ways"** (`docs/Multithreading.md`, IPDPS'14):
   simpler and historically safer than `level3_thread.c`; `JR`/`IR` ways share packed `A` among
   threads on a shared L2.
3. **`i?amax` NaN rule** (`ref_kernels/1/bli_amaxv_ref.c`, vectorized in
   `kernels/zen/1/bli_amaxv_zen_int.c`): first NaN, order independent, cheap in lanes.
4. **`normiv` / max-norm NaN propagation** (`frame/util/bli_util_unb_var1.c:449-486`).
5. **trsm pre-inversion as an option** (`configure:375-387`) with division for robustness.
6. **The sup path for small and skinny gemm** (`frame/3/bli_l3_sup*.c`, thresholds in
   `config/zen3/bli_cntx_init_zen3.c:244-253`).
7. **C-tile and A-panel prefetch in gemm microkernels**
   (`kernels/haswell/3/bli_gemm_haswell_asm_d6x8.c:127-175`).
8. **Larger `KC` on Apple cores** (`config/firestorm/bli_cntx_init_firestorm.c:99`).

Where BLIS is not better and rstd should not follow it:

- `nrm2`: BLIS's `sumsqv` is a scalar scaled loop with divisions; LAPACK 3.10's Blue's algorithm
  (and wide accumulation for `s`/`c`) is better.
- `dscal`: BLIS zeroes NaN for `alpha == 0` and scales for `incx <= 0`, both unlike the
  reference.
- `trmm`: BLIS multiplies packed zeros.
- `dsdot`/`sdsdot`: BLIS's are unoptimized scalar loops.
- 1m: not needed with a native `@Vector` complex kernel on every core.

## Uncertain points

- f64 `sdot` throughput on M6 and Zen: unmeasured after `4407e5a`.
- ARMv7 NEON and AltiVec flushing: seen on the rigs, not on hardware; Linux's AltiVec `NJ` default
  and the 970's software assist from memory.
- MKL's historical alignment requirement for CNR and oneTBB / ReproBLAS citations: from memory.
- The BLAS Technical Forum standard's wording on zero increments: not located.
- BLIS level-3 thread-count invariance: true by construction (no `k` split), not documented.
