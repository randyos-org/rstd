# Audit of the matlab port against the BLAS and LAPACK (M1)

**Status:** done (2026-10-30). Verdicts from reading the code, then checked by experiments E1 to
E9 against an `f128` LAPACK reference (Results, below; the scratch experiments are not kept). The
default is the modern routine; the old one wins on E1, which M2 keeps (`lu.zig`'s test), and on E2
and E4, which the modern LAPACK is to fix (todo.md's T125 and T126) rather than M2 keeping the old
routines (user, 2026-10-30).

## What the port is and who uses it

- `src/linalg/matlab/` (about 6,200 lines) ports matlab88.f90's LINPACK and EISPACK routines over
  its own level-1 kernels (`blas1.zig`), with split complex storage (`CxMat`: `r`, `i`, `lda`, `m`,
  `n`; 1-based `ix`).
- **Nothing outside the directory calls it.** The only reference in rstd is the re-export in
  `src/linalg/__root__.zig`; nothing in RandyOS, rgui or zfls imports `linalg.matlab`. Changing its
  API or layout breaks no one.

## Mapping

"Public" means exported by `src/linalg/lapack/__root__.zig` or `src/linalg/blas/__root__.zig`
today. Verdicts are provisional.

| matlab88 | Called by | Modern counterpart | Public? | Provisional verdict |
| --- | --- | --- | --- | --- |
| `wnrm2` | `wqrdc`, norms | `blas.nrm2` (`Complex64`) | yes | replace; the old fast path is wrong at partial underflow (finding 3) |
| `wasum`, `wdotur/ui/cr/ci`, `waxpy`, `wscal`, `wrscal`, `wcopy`, `wswap`, `rswap`, `rrot` | everywhere | `asum`, `dot`, `dotc`, `axpy`, `scal`, `scalReal`, `copy`, `swap`, `rot` | yes | replace; same arithmetic class, rounding noise only (E8) |
| `iwamax` | `wgefa` | `blas.index(.amax)` | yes | replace; same `cabs1` metric, but `iwamax` skips NaN where `izamax` returns the first NaN |
| `rrotg` | `wsvdc` | `blas.rotg` | yes | replace; the modern one scales safely |
| `wset`, `rset` | fills | LAPACK `laset` | no (T123) | replace once T123 lands; `@memset` until then |
| `pythag` | everywhere | `std.math.hypot` already; `lapy2` (T123) | no | already modern |
| `wmul`, `wdiv`, `wsign`, `wlog`, `watan`, `cabs1` | scalars | none (`blas.cabs1` exists) | partly | keep as scalar helpers |
| `matmul` (STACK2's `*`) | `CxMat.matmul`, `matrixfun` | `blas.gemm` | yes | replace; the old one is an unblocked dot-per-entry loop |
| `wgefa` | `det`, `inv`, `mldivide`, `mrdivide` | `getrf` | yes | replace; same partial pivoting by `cabs1` |
| `wgesl` | `mldivide`, `mrdivide` | `getrs` (`.n`, `.c`) | yes | replace |
| `wgeco` | `inv`, divisions | `lange(.o)` + `getrf` + `gecon(.o)` | yes | replace; `gecon`'s estimate is tighter (E6) |
| `wgedi` inverse | `inv` | `getri` | no (T123) | replace once T123 lands |
| `wgedi` determinant | `det` | none (a composition) | n/a | **keep** the scaling loop over `getrf`'s diagonal (finding 1) |
| `wpofa` | `chol` | `potrf(.u)` | yes | replace (finding 4) |
| `wqrdc` unpivoted / pivoted | `orth`, `qr` | `geqrf` / `geqp3` | yes | replace; faer's `geqp3` recomputes column norms each step, so LINPACK's norm downdating can't fail there |
| `wqrsl` | least squares, `orth` | `ormqr`/`unmqr` + `trtrs`, or `gels` | yes | replace; `gels` returns `error.Singular` where `wqrsl` returns an index |
| `wsvdc` | `cond`, `rank`, `pinv`, `matNorm2` | `gesvd` (faer), `gesdd` | yes | replace once T126 keeps graded bidiagonals' small values (E4) |
| `htridi` + `imtql2` + `htribk` | Hermitian `eig` | `heevd`/`heev` (`syevd`/`syev` for real) | yes | replace once T125 reduces from the given triangle's end (E2) |
| `corth` | `hess`, `schur`, `eig` | `gehrd` | no (T123) | replace once T123 lands |
| `comqr3` jobs 0/1 (Schur) | `schur`, `roots` | `hseqr`, `gees` | no (T123) | replace once T123 lands |
| `comqr3` job 2 (eigenvectors) | general `eig` | `geev` | yes | replace; balancing helps (E3) |

## Findings so far

1. **Determinant range (keep).** `wgedi` holds `det` as `mantissa * 10^exponent`, so
   `det(diag(1e200, 1e200))` is `1 * 10^400`; a plain product of `getrf`'s diagonal overflows.
   LAPACK has no determinant, so this is no LAPACK flaw, but M2 keeps the scaling loop over
   `getrf`, with a test (E1).
2. **Balancing.** matlab88's `eig`, `schur` and `roots` never balance; `geev` always does (`gebal`
   job `B`). Watkins (2006) and MATLAB's `eig(A,'nobalance')` example give a 4x4 matrix that
   balancing ruins (E3). E3 did not show that harm (old and `geev` alike, 1e-16), and on a badly
   scaled `D B D^-1` the unbalanced old path is off by 2e1 where `geev` gets 2e-15: replace, with
   balancing (a no-balance option through T123's `gebal`, `gehrd`, `hseqr`, `trevc` stays possible
   for MATLAB's `'nobalance'`). `roots` should balance, as MATLAB's does.
3. **`wnrm2` loses digits at partial underflow (replace).** Its fast path squares unscaled and keeps
   any finite nonzero sum, so entries near 1e-160 give denormal squares and the fallback never runs.
   `blas.nrm2` is scaled. Caveat (`faer/qr_pivot.zig`): on cores whose `nrm2` kernel is an unscaled
   extended sum, a matrix near that type's underflow reads as zero; a per-core BLAS matter.
4. **Cholesky (replace; M2 picks the behavior).** `wpofa` takes a NaN diagonal as success (`s <= 0`
   is false for NaN); `potrf` reports it as not positive definite. `wpofa` rejects a nonzero
   imaginary part on the diagonal; `potrf` ignores it, as LAPACK and modern MATLAB do.
5. **Old edge-case bugs (replace).** `imtql2` forms `sqrt(g*g + 1)` without `pythag`, so it can
   overflow; `htridi` and `imtql2` index `[n-1]` for `n == 0`; `comqr3` and `imtql2` don't scale to
   a safe range (`geev` scales into `[SMLNUM, BIGNUM]`). Modern routines return `error.InvalidA` or
   `error.NoConvergence` for NaN or Inf where the old ones return an iteration-failure index.

## Results (E1 to E9, 2026-10-30)

Native (Zen 3, Debug), against an `f128` LAPACK reference; errors relative to the reference value
unless said.

| Exp. | Case | Old | Modern | Verdict |
| --- | --- | --- | --- | --- |
| E1 | `det(diag(1e200, 1e200, -3e-100))` | `-3 * 10^300` | `getrf` diagonal product `-inf` | **keep** the scaling (finding 1) |
| E2 | Hermitian graded large at the top left (`inc=false`) | 1e2 to 1e8 | `heevd`/`heev` 1e-15 | modern |
| E2 | graded large at the bottom right (`inc=true`) | 1e-15 to 1e-13 | 0.3 to 3, 3e8 at ratio 1e-4 | **old**: modern loses the small eigenvalues' relative accuracy |
| E3 | Watkins' balancing example, random | 1e-16 to 1e-15 | same | even: the balancing harm did not show |
| E3 | `D B D^-1` (badly scaled) | eigenvalue error 2e1 | 2e-15 | modern (balancing wins) |
| E4 | random, column- and row-graded | 1e-16 to 1e-15 | same | even |
| E4 | column-graded reversed | 7e-6 to 1e-4 | 2e-5 to 6e-5 | even (both lose) |
| E4 | graded bidiagonal (`smin` 9e-21) | 2e-16 | `gesvd`/`gesdd` 6e3 | **old**: modern loses the small singular values' relative accuracy |
| E5, E6 | rcond: random, Hilbert, Kahan, lower minus ones | 1.4x to 2.7x the true value | 1.0x to 1.8x | modern |
| E5 | solve error | same order | same order | even |
| E7 | least squares, random and column-graded | 1e-16 to 1e-15 | same | even |
| E8 | `nrm2` at 1e-160, mixed tiny | 6e-6 | 4e-18, 7e-17 | modern (finding 3) |
| E8 | `iamax([1, NaN, 5])` | 3 (skips NaN) | 2 (first NaN) | modern (LAPACK's rule) |
| E8 | `rotg` extremes; `dotc`, `asum` means | same | same | even |
| E9 | Cholesky accuracy | same order | same order | even |
| E9 | NaN diagonal; diagonal `4+i` | success; rejected | `NotPositiveDefinite`; accepted (`r11 = 2`) | modern (finding 4) |

- **E2:** `heevd` and `heev` give the same errors for `.l` and `.u`, so faer's tridiagonal
  reduction runs one way whichever triangle is given; matlab88's `htridi` reduces from the bottom.
  Reference LAPACK's `xSYTRD` reduces from the bottom for `.u` (`dsytrd.f`), which may be why the
  reference's own users rarely meet this.
- **E4:** reference LAPACK's `xBDSQR` computes a bidiagonal's singular values to high relative
  accuracy (Demmel and Kahan's zero-shift QR); LINPACK's `wsvdc` keeps it here, faer's
  bidiagonal SVD does not.

## M2's layout: split `xr`/`xi` against interleaved `Complex`

No external callers; every modern routine takes `(allocator, io, ...)` while most matlab kernels
take no allocator, so signatures change under M2 whatever the layout; LAPACK's real outputs are
already split (`geev`'s `w = .{ .re, .im }`); matlab88 users' matrices are mostly real.

- **A, convert at the boundary (recommended for M2).** A real input (`i` all zero) goes zero-copy to
  the `d` routines, in place on `r`, which also gives `dgeev`'s exact conjugate pairs. A complex
  input is copied once into interleaved scratch, run through the `z` routine and copied back. Cost:
  `toInterleaved`/`fromInterleaved`, one `16 * m * n`-byte allocation per complex call and an
  `O(mn)` zero scan: small against `O(n^3)` factorizations, about 2x for `O(mn)` operations. Open:
  whether switching real inputs to the `d` routines fits M2's "no behavior change beyond M1".
- **B, `CxMat` as `[]Complex(f64)`.** Matches the BLAS, LAPACK and M3 (MATLAB went interleaved in
  R2018a), but changes every file and fixture here, loses zero-copy real dispatch, and M3 is
  expected to redefine storage anyway.

Recommendation: A in M2, final storage left to M3. The split layout's only reason (`blas1.zig`'s
`@Vector` stride-1 paths) goes away with M2, so neither M2 nor M3 should keep it for that.

## Follow-ups

- T125 and T126 fix E2 and E4 in the modern LAPACK; their tests take E2's and E4's matrices.
- M2 keeps E1's scaling over `getrf` (`lu.zig`'s test "wgedi's determinant keeps its range").
- Not examined: `matrixfun.zig` (`expm` goes through `eig`; scaling and squaring is modern
  practice, outside M1); `svdapi.zig`'s `rank`/`pinv` tolerances; pivot-order ties between `wqrdc`
  and `geqp3`; `corth`/`comqr3` against `gehrd`/`hseqr` (needs T123).

For M2's calls: `getrs` takes `.n`, `.t` and `.c` but rejects conjugate-no-transpose
(`trans.conjX()`), and `ormqr` rejects `.t` for complex types; use `.c` for the `ctrans(A)` solves
that `wgesl`'s `job != 0` does.
