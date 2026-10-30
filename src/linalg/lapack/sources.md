# LAPACK and kernel sources

The anthology's source table (design.md, Settled decisions: "`rstd.linalg` is an anthology of
the best open-source methods", 2026-10-30; todo.md T113). For each LAPACK routine rstd has or
plans, and each BLAS kernel family with more than one source, it lists the candidates, which of
them implement it (read in their trees at the pins below, not guessed), the current choice, why,
and the evidence. A choice changes when a benchmark or an accuracy test favors another source.
Reference-LAPACK and Eigen are baselines only; Eigen is never a source.

## Candidates

| Candidate | Pin | License | What it offers |
| --- | --- | --- | --- |
| OpenBLAS | `v0.3.34-317-ga69062524` | BSD-3 | its own LAPACK (`interface/lapack/`, `lapack/`), the kernels rstd converts |
| ReLAPACK | OpenBLAS's `relapack/`, same pin | MIT (`relapack/LICENSE`) | recursive factorizations (E. Peise, P. Bientinesi, "Algorithm 979: Recursive Algorithms for Dense Linear Algebra: The ReLAPACK Collection", ACM TOMS 44(2), 2017, <https://doi.org/10.1145/3061664>); list in `relapack/coverage.md` |
| faer | tag `faer-v0.24.4` (`0539947ffb`), <https://codeberg.org/sarah-quinones/faer> | MIT | dense decompositions in Rust (`faer/src/linalg/`) |
| libflame | tag `5.2.0` (`53bd60c50b`, 2019-08-09), <https://github.com/flame/libflame> | BSD-3 (`LICENSE`); its f2c parts under AT&T, Lucent and Bellcore's permission notice (same file) | FLAME algorithms behind LAPACK names |
| AOCL-libFLAME | tag `5.3.2` (`3b3251adec`, 2026-08-20), <https://github.com/amd/libflame> | BSD-3 and MIT, per file (`LICENSE.txt`); third-party parts BSD-3 or MIT (`NOTICE.txt`) | libflame plus AMD's x86 tuning |
| BLIS | `2.1` (`caf0db6be1`) | BSD-3 | kernels and the level-3 framework; no LAPACK (no `getrf`, `potrf` or `geqrf` under `frame/` or `kernels/`) |
| Apple Accelerate | macOS's, bound at run time | Apple's; called, never ported | the routines `apple_accel.zig` routes ([BLAS design](../blas/design.md#apple-accelerate)) |
| Reference-LAPACK | OpenBLAS's `lapack-netlib/` | baseline | the semantics, and a benchmark baseline |
| Eigen | upstream, built through Nix (todo.md T73) | MPL 2.0 | a benchmark baseline only |

The libflame and AOCL-libFLAME clones sit beside the OpenBLAS and BLIS checkouts, read-only.
libflame's `master` is 62 commits past `5.2.0` (to `f38dd83e`, 2024-05-20) with no later tag.
AOCL-libFLAME's `5.3.2` is its default branch's head; its dated snapshot tags (`AOCL-20261001`)
are newer but not releases.

### How libflame and AOCL-libFLAME serve a LAPACK name

Both put LAPACK's names over three kinds of code (`src/map/lapack2flamec/README`):

- **FLAME**: `FLA_<name>.c` maps the name onto a FLAME algorithm (`FLA_LU_piv`, `FLA_Chol`,
  `FLA_QR_UT`, ...), under `FLA_ENABLE_LAPACK2FLAME`. `bdsqr`, `gehrd`, `larfg` and `trsyl` map
  only under `FLA_ENABLE_YET_LAPACK2FLAME`, off by default (marked "off" below).
- **ref**: everything else is Reference-LAPACK run through f2c (libflame: 3.5.0,
  `f2c/netlib/`; AOCL-libFLAME: up to 3.12, `NOTICE.txt`). That is the reference's method, not
  a source of its own.
- **AMD**: with `FLA_ENABLE_AMD_OPT`, AOCL-libFLAME sends several real routines to its own
  modified f2c copies instead, with AVX2 and AVX-512 small-size kernels
  (`src/lapack/x86/avx2/`, `avx512/`, dispatched by `src/lapack/x86/front/fla_lapack_x86_common.c`)
  or OpenMP (`FLA_OPENMP_MULTITHREADING`). Only paths seen calling such code are marked AMD.

Two of libflame's names are not what they say: `gesdd` is "redirected to GESVD"
(`FLA_gesdd.c:36`), so it has no divide and conquer, and `gelsd` to `gelss` (`FLA_gelsd.c:51`).

## LAPACK routines

"Needs" says why rstd wants the routine: **has** (in `src/linalg/lapack/` on `review`), an item
(T34 SVD, T110 ReLAPACK, T111 the rest from faer, T70), or **M** for `rstd.linalg`'s MATLAB-like
API (todo.md M1 to M4; M2 names the LINPACK and EISPACK routines each replaces). In the source
columns "own" is the candidate's own implementation, with its file (OpenBLAS's from its root,
ReLAPACK's in `relapack/src/`, faer's in `faer/src/linalg/`, libflame's and AOCL-libFLAME's in
`src/map/lapack2flamec/` and its `f2c/c/`); "-" is none. Evidence for
**has** rows is the accuracy tests design.md's Verification lists for `src/linalg/lapack/`.

### Cholesky and triangular

| Routine | Needs | OpenBLAS | ReLAPACK | faer | libflame | AOCL-libFLAME | Choice | Why | Evidence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `potrf` | has | own, `lapack/potrf/` | own, `relapack/src/dpotrf.c` | own, `cholesky/llt/factor.rs` | FLAME `FLA_Chol` | AMD: `fla_dpotrf_small_avx2`, `_avx512`; `dpotrf_var2`, OpenMP (`FLA_potrf.c:197-220`) | OpenBLAS; Accelerate on macOS | Q10: OpenBLAS's own first; ReLAPACK's a build option (T110) | accuracy tests; not raced |
| `potf2` | has | own, `lapack/potf2/` | - | inside `llt/factor.rs` | FLAME | ref | OpenBLAS | `potrf`'s base case | accuracy tests; not raced |
| `potrs` | has | ref (`exports/gensymbol.pl:498`) | - | own, `llt/solve.rs` | ref | ref | the reference's two `trsm`s over rstd's drivers | the method is two solves; the drivers are rstd's | accuracy tests; not raced |
| `potri` | has | own, `interface/lapack/potri.c` | (`trtri`, `lauum`) | own, `llt/inverse.rs` | FLAME `FLA_Trinv`, `FLA_Ttmm` | AMD: `fla_dpotri_small_avx2` (`FLA_potri.c:177`) | OpenBLAS | Q10 | accuracy tests; not raced |
| `pbtrf` | T110 (ReLAPACK covers it; T110 does not name it) | - | own, `dpbtrf.c` | - | ref | ref | ReLAPACK | only non-reference source | accuracy tests; not raced |
| `trtri`, `trti2` | has | own, `lapack/trtri/`, `trti2/` | own (`trtri`) | own, `triangular_inverse.rs` | FLAME `FLA_Trinv` | ref | OpenBLAS | Q10 | accuracy tests; not raced |
| `lauum`, `lauu2` | has | own, `lapack/lauum/`, `lauu2/` | own (`lauum`) | inside `llt/inverse.rs` (`matmul::triangular`) | FLAME `FLA_Ttmm` | ref | OpenBLAS | Q10 | accuracy tests; not raced |
| `trtrs` | has | own, `lapack/trtrs/` | - | own, `triangular_solve.rs` | ref | ref | OpenBLAS | Q10 | accuracy tests; not raced |

### LU and general solves

| Routine | Needs | OpenBLAS | ReLAPACK | faer | libflame | AOCL-libFLAME | Choice | Why | Evidence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `getrf`, `getf2` | has | own, `lapack/getrf/`, `getf2/` | own (`getrf`), `dgetrf.c` | own, `lu/partial_pivoting/factor.rs` | FLAME `FLA_LU_piv` | AMD: `fla_[sdz]getrf_small_avx2`, `_avx512`; OpenMP `FLA_LU_piv_[sdz]_parallel` (`FLA_getrf.c:203-291`) | OpenBLAS | Q10 | accuracy tests; not raced |
| `getrs` | has | own, `lapack/getrs/` | - | own, `partial_pivoting/solve.rs` | ref | AMD: `fla_dgetrs_small_notrans` (`dgetrs.c`) | OpenBLAS | Q10 | accuracy tests; not raced |
| `gesv` | has | own, `interface/lapack/gesv.c` | - | own (`solvers.rs`) | ref | ref | OpenBLAS | Q10 | accuracy tests; not raced |
| `laswp` | has (`blas/ext/`) | own, `lapack/laswp/` | - | permutations inside `lu/` | ref | ref | OpenBLAS | Q10 | accuracy tests; not raced |
| `getri` | M (`wgedi`, `inv`) | - | - | own, `partial_pivoting/inverse.rs` | ref | AMD: `lapack_getri_small_d`, plain C (`src/lapack/x86/generic/`) | faer | the anthology rule: faer before the reference | not yet raced |
| `gecon` | M (`wgeco`, `rcond`, `cond`) | - | - | - | ref | ref | open | only the reference's estimator (`xLACN2`) exists; needs a decision | none |
| `gbtrf` | T110 | - | own, `dgbtrf.c` | - | ref | ref | the reference's `xGBTF2` | ReLAPACK's recursion gives a wrong factor ([issues](../blas/issues.md#relapack-openblass-relapack)) | accuracy tests; not raced |

### Symmetric indefinite and generalized definite

| Routine | Needs | OpenBLAS | ReLAPACK | faer | libflame | AOCL-libFLAME | Choice | Why | Evidence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `sytrf`, `hetrf` (and `_rook`) | T110; M (`ldl`, symmetric `mldivide`) | - | own, recursive, `dsytrf.c`, `dsytrf_rook.c` | own, Bunch-Kaufman with partial, rook and rook-diagonal pivoting, `cholesky/bunch_kaufman/factor.rs` | ref | ref | the reference's blocked loop over ReLAPACK's panels | ReLAPACK's recursion gives a wrong factor ([issues](../blas/issues.md#relapack-openblass-relapack)); faer's storage is not LAPACK's | accuracy tests; not raced; faer is the other candidate |
| `sytrs`, `sytri` (`hetrs`, `hetri`) | M | - | - | own, `bunch_kaufman/solve.rs`, `inverse.rs` | ref | ref | open | faer's solve reads faer's factor, so it follows the `sytrf` choice | none |
| `sygst`, `hegst` | T110; M (`eig(A, B)`) | - | own, `dsygst.c` | - | FLAME `FLA_Eig_gest` | FLAME | ReLAPACK | T110 | accuracy tests; not raced; FLAME is the other candidate |
| `sygv`, `hegv`, `sygvd` | M (`eig(A, B)`) | - | - | - | ref | ref | `potrf`, `sygst`, then `syevd` as chosen | a driver over the rows above | none |

### Sylvester

| Routine | Needs | OpenBLAS | ReLAPACK | faer | libflame | AOCL-libFLAME | Choice | Why | Evidence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `trsyl` | T110; M (`sylvester`, `lyap`, `sqrtm`, `logm`) | - | own, recursive, `dtrsyl.c` | - (a 2 x 2 `lasy2` inside `evd/schur/real_schur.rs` only) | FLAME `FLA_Sylv`, off | FLAME `FLA_Sylv`, off | ReLAPACK | T110 | accuracy tests; not raced |
| `tgsyl` | T110 (covered; not named) | - | own, `dtgsyl.c` | - | ref | ref | ReLAPACK | only non-reference source | not yet raced |

### QR, LQ and least squares

| Routine | Needs | OpenBLAS | ReLAPACK | faer | libflame | AOCL-libFLAME | Choice | Why | Evidence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `geqrf` | T111 (and inside T34) | - | - (QR is "not covered: extra FLOPs", `relapack/coverage.md`) | own, `qr/no_pivoting/factor.rs` | FLAME `FLA_QR_UT` | AMD: `fla_[sdcz]geqrf_small` (AVX2), OpenMP (`dgeqrf.c:224`) | faer | T111; ported on `wip/faer` for T34 | not yet raced |
| `orgqr`, `ungqr` | T111 | - | - | own, `householder.rs` (the block sequence on `I`) | FLAME `FLA_QR_UT_form_Q` | FLAME | faer | T111 | not yet raced |
| `ormqr`, `unmqr` | T111; M (`wqrsl`) | - | - | own, `householder.rs` | FLAME `FLA_Apply_Q_UT` | AMD: OpenMP (`dormqr.c`) | faer | T111 | not yet raced |
| `geqp3` | M (`wqrdc`, pivoted `qr`, `rank`) | - | - | own, `qr/col_pivoting/factor.rs` | FLAME `FLA_QR_UT_piv` (`FLA_geqpf.c`) | AMD: `dlaqp2`'s norms (AVX2) | faer | the anthology rule | not yet raced |
| `gelqf`, `orglq`, `ormlq` | M (wide minimum-norm solves) | - | - | - | FLAME `FLA_LQ_UT` | AMD: `fla_dgelqf_small` | open | FLAME is the only non-reference source; faer's QR of `Aᴴ` is the alternative | none |
| `gels` | T111; M (rectangular `mldivide`) | - | - | own, `qr/no_pivoting/solve.rs` (`solve_lstsq_in_place`) | ref | ref | faer | T111 | not yet raced |
| `gelsd`, `gelss`, `gelsy` | M (`pinv`, `lsqminnorm`, rank-deficient `mldivide`) | - | - | SVD pseudoinverse, `svd/mod.rs` | `gelsd` is `gelss` (ref) | ref | open | faer's SVD covers `gelss`; `gelsy`'s complete orthogonal factorization has no non-reference source | none |

### SVD

| Routine | Needs | OpenBLAS | ReLAPACK | faer | libflame | AOCL-libFLAME | Choice | Why | Evidence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `gesvd` | T34 (`wsvdc`, `svd`, `rank`, `null`, `orth`) | - (OpenBLAS's `test_kernel_regress.c` needs it) | - | own, `svd/mod.rs`, QR iteration | FLAME `FLA_Svd_ext`: restructured bidiagonal QR (below) | AMD: `fla_dgesvd_small6`, `_xx_small10` and others (AVX2, `dgesvd.c`) | faer (`faer/svd.zig`) | T34 | `A - U Σ Vᴴ` and orthogonality within faer's `test_svd` tolerance, every job and precision, and `test_kernel_regress.c`; four faer fixes (`../blas/issues.md`); not raced |
| `gesdd` | T34 | - | - | own, divide and conquer, `svd/bidiag_svd.rs` | redirected to `gesvd` | ref (`dgesdd.c` over `dbdsdc`) | faer (`faer/svd.zig`) | T34 | as `gesvd`; not raced |
| `gebrd`, `orgbr`, `ormbr` | inside T34 | - | - (not covered) | own, `svd/bidiag.rs` | FLAME `FLA_Bidiag_UT` | FLAME | faer | T34 | with `gesvd` |
| `bdsqr`, `bdsdc` | inside T34 | - | - | own, `svd/bidiag_svd.rs` | FLAME `FLA_Bsvd`, off | AMD: `dlasr` rotations (AVX2), `dbdsqr_small` | faer | T34 | with `gesvd` |

### Eigenvalues

| Routine | Needs | OpenBLAS | ReLAPACK | faer | libflame | AOCL-libFLAME | Choice | Why | Evidence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `sytrd`, `hetrd` | T111 | - | - (not covered) | own, `evd/tridiag.rs` | FLAME `FLA_Tridiag_UT` | AMD: OpenMP (`dsytrd.c:255`) | faer | T111 | not yet raced |
| `orgtr`, `ormtr`, `unmtr` | T111; M (`htribk`) | - | - | own, `householder.rs` | FLAME | FLAME | faer | T111 | not yet raced |
| `syevd`, `heevd` (and `syev`, `steqr`, `stedc`) | T111; M (symmetric `eig`) | - | - | own, `evd/mod.rs`: QR iteration below 128, divide and conquer above (`tridiag_evd.rs`) | ref; `FLA_Hevd` and `FLA_Tevd` exist but no LAPACK name maps to them | ref | faer | T111 | not yet raced |
| `laed3` (with `laed4`) | T70 (not ported) | own, `lapack/laed3/` | - | faer's own merge, `tridiag_evd.rs` (`secular_eq`, `divide_and_conquer`) | ref | ref | faer's, inside `syevd` | OpenBLAS's serves only the reference's `stedc`, which rstd does not port | not yet raced |
| `gehrd`, `orghr`, `ormhr` | T111 (`geev`); M (`corth`, `hess`) | - | - (not covered) | own, `evd/hessenberg.rs` | FLAME `FLA_Hess_UT`, off | ref | faer | T111 | not yet raced |
| `hseqr` | T111; M (`comqr3`) | - | - | own, multishift QR with aggressive early deflation, `evd/schur/` | ref | ref | faer | T111 | not yet raced |
| `trevc` | T111; M (`comqr3`) | - | - | own, `evd/mod.rs` (`evd_from_real_schur_imp`, `_cplx_`) | ref | ref | faer | T111 | not yet raced |
| `geev` | T111; M (`eig`) | - | - | own, `evd/mod.rs` (`evd_real`, `evd_cplx`) | ref | ref | faer | T111 | not yet raced |
| `gees` | M (`schur`) | - | - | own, `evd/schur/` | ref | ref | faer | the anthology rule | not yet raced |
| `trexc`, `trsen` | M (`ordschur`) | - | - | swaps only, `real_schur.rs` (`schur_swap`), crate-internal | ref | ref | open | no reordering driver outside the reference | none |
| `gebal`, `gebak` | M (`eig`'s balancing, `balance`) | - | - | - | ref | ref | open | only the reference; needs a decision | none |
| `ggev`, `gges` (`gghrd`, `hgeqz`, `tgevc`) | M (`eig(A, B)`, `qz`) | - | - | own, `gevd/` (`gen_hessenberg`, `qz_real`, `qz_cplx`) | ref | AMD: `dhgeqz` rotations, `fla_dhrot3` (AVX2, AVX-512) | faer | the anthology rule | not yet raced |

### Auxiliary

| Routine | Needs | OpenBLAS | ReLAPACK | faer | libflame | AOCL-libFLAME | Choice | Why | Evidence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `lange`, `lansy`, `lanhe` | M (`norm`) | - | - | `reductions/` (`norm_max`, `norm_l2` Frobenius, `norm_l1` entrywise; no operator 1- or infinity-norm) | ref | AMD: `lange` kernels (AVX2, AVX-512) | open | rstd's own over its BLAS (`amax`, `asum`, `nrm2`) is the obvious candidate | none |
| `laset`, `lacpy`, `lascl`, `lapy2` | M2 (`wset`, `rset`, `pythag`) | - | - | element loops | ref | ref (AMD `dlacpy`) | rstd's own, the reference's semantics | too small to source | none |

faer also has what no item asks for yet: pivoted Cholesky (`llt_pivoting`, LAPACK's `pstrf`),
`LDLᴴ` without pivoting, and full-pivoting LU.

## BLAS kernel families

Where more than one source serves a family. The per-core detail is in the `coverage.md` files
(one per level-1 routine, [`blas2/coverage.md`](../blas/blas2/coverage.md),
[`blas3/gemm/coverage.md`](../blas/blas3/gemm/coverage.md)); the reasoning in
[BLAS design](../blas/design.md) and [`comparisons.md`](../blas/comparisons.md); the numbers in
`tools/blasbench/results-*.md`.

| Family | Candidates | Choice | Why | Evidence |
| --- | --- | --- | --- | --- |
| Level-3 tile (`gemm` and the kernels over it) | OpenBLAS's per-core assembly (`gemm/coverage.md`), BLIS's microkernels, rstd's `@Vector` tile, Accelerate | the `@Vector` tile on every core; Accelerate on macOS | one source for every core, BLIS's reference-kernel idea with a register model (comparisons.md 9) | `gemm/notes.md` (about 100% of measured peak on Zen 3); flyoff: real `gemm` 0.97 to 1.30 of OpenBLAS on Zen 3, complex 0.79 to 1.01, 0.94 to 1.13 on the M6; Accelerate 6.2 to 17 times (SME) |
| Level-3 loops and threading | OpenBLAS `level3.c`, `level3_thread.c`; BLIS's ways and barriers | `level3.c`'s loops, BLIS's threading | Settled decisions, "Level-3 threading follows BLIS" | flyoff threaded rows, indicative only |
| `beta` on `C` | OpenBLAS's separate pass; BLIS's first `KC` slice | BLIS's | same bits, one pass less (comparisons.md 7) | bits by construction; speed not raced alone |
| `trsm` diagonal | OpenBLAS's inverted diagonal; BLIS's division option | division | `1/d` overflows for denormal `d` (comparisons.md 10) | special-value tests; speed not raced |
| Small and skinny `gemm` | OpenBLAS's small-matrix kernels; BLIS's `sup` path | none yet (T30) | not ported | not yet raced |
| `gemmt` | rstd's driver (OpenBLAS's); ReLAPACK's recursive `xgemmt`; BLIS's `gemmt` | rstd's driver; ReLAPACK's for its LAPACK callers (T110) | T110 | not yet raced |
| Apple blocking | OpenBLAS `VORTEX` `P`, `Q`; BLIS Firestorm `KC` 4096 and 3072 | OpenBLAS's | the extracted tables | not yet raced (comparisons.md 9's sweep) |
| Routed reductions (extremes, `nrm2`, `sdot`; x86_64's peeling reductions) | converted assembly; `@Vector` (`vector.zig`); Accelerate | `@Vector` where the vector unit is IEEE, else portable; Accelerate on macOS only where `apple_accel.zig`'s `sends` says | Settled decisions, "No address dependence"; one contract (comparisons.md X1) | flyoff: `nrm2` up to 11.7 times OpenBLAS on Zen 3, `asum` 0.67 to 2.6; macOS routing table from results-flyoff-macos-aarch64.md |
| `nrm2` algorithm | OpenBLAS's per-arch loops; BLIS's `sumsqv`; Reference-LAPACK 3.10's Blue's algorithm | `f64` sums for `s` and `c`, Blue's for `d` and `z` | accuracy (comparisons.md 2) | flyoff accuracy: BLIS `snrm2`, `scnrm2` and OpenBLAS arm64 `scnrm2` about 1% off at 16M |
| `i?amax` NaN rule | OpenBLAS's; BLIS's first NaN; Reference-LAPACK master | first NaN | comparisons.md 3 | edge tests |
| arm64 `vortex`, `vortexm4` level 1 and 2 | converted assembly; `@Vector` | `@Vector` in every slot but `rot` and `sum` | it matches or beats the assembly on the M6 | results-m6.md; results-m6-libraries.md |
| arm64 slots whose upstream kernel is C | OpenBLAS's C (not ported); `@Vector`; portable loops | `@Vector` | the C kernels are unported; `@Vector` matches them | results-m6-libraries.md: `@Vector` against the real C, 0.82 to 2.05 |
| Other level-1 slots (`axpy`, `dot`, `rot`, `scal`, `copy`, `swap`, `asum`, `sum` per arch) | converted assembly or ported C drivers; `@Vector`; portable | converted where upstream has assembly; `@Vector` only where routed | upstream's own selection per core (design.md, Kernel selection) | flyoff on Zen 3: `copy`, `swap`, complex `rot` 0.94 to 1.05, `axpy` 0.85 to 1.8; `@Vector` against the converted kernels not raced off arm64 |
| `gemv` | converted assembly and C drivers; `@Vector` (every form); portable; Accelerate | converted where present, `@Vector` on arm64 slots with none, portable elsewhere (T22's gaps) | as above | flyoff: Zen 3 `gemv` 0.06 to 0.61 on the portable loop; M6 `sgemv_n` `@Vector` 0.82 to 0.97 of OpenBLAS's C |
| `symv`, `hemv` | x86_64's converted kernels; the portable composition (`symv_k.c`) | converted on x86_64, portable elsewhere | upstream's own selection per core | not yet raced |
| `ger` and rank updates | OpenBLAS's per-core `ger` kernels (Sandy Bridge, A64FX, Neoverse V1); the portable driver | portable | not ported (T22) | not yet raced |
| `laswp_ncopy`, `neg_tcopy` | OpenBLAS only | unused | `getrf` is not fused (design.md, LAPACK) | none |

## Notably better in libflame or AOCL-libFLAME

Candidates for the timing session; none is ported, and each would race the current choice.

1. **AOCL-libFLAME's small-size x86 kernels** (`src/lapack/x86/avx2/`, `avx512/`): `getrf` for
   `s`, `d` and `z` (`d` below 81 on AVX2 and 160 on AVX-512, `FLA_macro_defs.h:294-296`;
   `FLA_getrf.c:203-291`), `potrf` and `potri`
   for `d`, `getrs` untransposed, `getri` (plain C), `geqrf` for every type, `gelqf`, and
   `gesvd`'s small cases (`FLA_DGESVD_SMALL_SIZE_THRESH*`, 16 to 386, `FLA_macro_defs.h:349-351`), with `dlasr`, `larf` and `lange` kernels underneath. `rstd.linalg`'s
   MATLAB-like layer calls small matrices often (M3), and Zen 3 is a host rstd times. They are
   intrinsics; a port would be `@Vector` code, as T20's strategy does for level 1 and 2.
2. **AOCL-libFLAME's OpenMP paths**: `geqrf`, `ormqr`, `sytrd`, `potrf` (`dpotrf_var2`) and
   `getrf` (`FLA_LU_piv_*_parallel`). Each must first be shown to give the same bits for any
   thread count before it can be a candidate.
3. **libflame's restructured tridiagonal and bidiagonal QR** (F. G. Van Zee, R. A. van de Geijn,
   G. Quintana-Ortí, "Restructuring the Tridiagonal and Bidiagonal QR Algorithms for
   Performance", ACM TOMS 40(3), 2014, <https://doi.org/10.1145/2535371>): Givens rotations
   applied in waves, so vectors update at level-3-like speed. It is libflame's `gesvd`
   (`FLA_Svd_ext`) and its `FLA_Hevd` and `FLA_Tevd`, which no LAPACK name reaches. It would race
   faer's QR iteration in `gesvd` and divide and conquer in `syevd`.
4. **libflame's `FLA_Eig_gest` and `FLA_Sylv`**: the only non-reference rivals to ReLAPACK's
   `sygst` and `trsyl` (T110); `FLA_Sylv` is mapped only with `FLA_ENABLE_YET_LAPACK2FLAME`.
5. **libflame's `FLA_LQ_UT`**: the only non-reference `gelqf`.

Not worth a race: libflame's `gesdd` (it is `gesvd`) and `gelsd` (it is `gelss`), and
AOCL-libFLAME's non-LAPACK extensions (`getrfnp`, `getrfnpi`, `spffrt2`, `spffrtx`).
