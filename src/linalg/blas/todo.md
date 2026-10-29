# BLAS to-do

The BLAS backlog in priority order: work goes top to bottom. IDs stay fixed when items move.
**Needs:** names items that must land first; **Blocked:** names an outside cause (hardware, an
upstream fix, a decision). Items are removed once done and recorded (commits, design.md,
issues.md). Upstream problems go to issues.md as they are found; the user triages them after the
BLAS work.

## Feature complete means

The port and its optimization are done, and other rstd work can start, when all of these items hold:

1. **Every promised routine exists**, for every element type it promises: no public routine
   reaches a null or stub slot, and the reference BLAS and OpenBLAS's column-major extensions
   are covered (T30; batch forms optional, T75).
2. **Threading is done**: threading.md's plan, workers on the caller's floating-point state, and
   results bit-identical for any worker count and `Io` on every hosted rig (T01, T02, T03).
3. **Every kernel object a build can link runs the per-core tests** on some rig or machine, or
   its item says **Blocked:** and why (T43, T65, T68).
4. **No undecided kernel gap**: every bold `coverage.md` cell, level 2 included, is ported,
   replaced by a measured `@Vector` kernel, or left portable with the measurement recorded
   (T20, T22, T67).
5. **OpenBLAS's suites pass in full**, skipping only what `_test/from_openblas/__root__.zig`
   lists as untranslatable (T69; `dgesvd` per T34).
6. **The rig matrix is green** in one sweep on `review` (T62), wasm included (T64), with LAPACK
   on the rigs (T08).
7. **It builds where it is used**: the dev toolchain and Zig 0.17.0 (zfls) for zfls's targets,
   plus `rstd-check-macho` and `rstd-check-uefi` (T66).
8. **The docs describe the code**: no "not yet", "planned" or "serial" for landed work (T63), and
   every issues.md workaround has a removal condition here (T52).
9. **Speed is accounted for**: thresholds tuned by benchmark (T05), and the flyoff shows rstd at
   or above OpenBLAS within noise, or the gap has an item (T21, T57, T73).
10. **The public interface is final**: Z1 to Z5 done, or scheduled with M3.

## Running now

- [ ] **T03** Workers adopt the caller's floating-point control state (design.md, Settled
  decisions): capture it on the caller in `threading.team` and `forEachIn`, apply it as each task's
  first step (x87 control word and MXCSR on x86, FPCR on arm64, the others' equivalents). Then a
  `-Dcpu=pentium3` build's threaded results match serial ones (on Windows new threads start at
  53-bit x87 precision), and so do callers who change rounding or flush-to-zero. Running
  (`wip/fpstate`).
- [ ] **T25** Speed follow-ups: a `@Vector` `sum` (`csum`/`zsum` are scalar on the old x86_64 cores;
  it also fills arm64's `sum` cells); measure the scalar fallbacks (extremes and `sdot` on arm,
  pre-VSX POWER, mips64, SPARC) and the `f64` `sdot` (comparisons.md (1)). Running (`wip/vecsum`).
- [ ] **T55** The full i386 rig (the trimmed 32-bit objects, 886324d; the x87 audit's checks) stalls
  in "interface: symm and hemm match the wide reference" and the watchdog kills it, after the
  level-3 threading landed; investigating (`wip/i386stall`).
- [ ] **T33** LAPACK `getrf`, then call `laswp_ncopy` and `neg_tcopy` from it, or drop them as
  upstream's `#if 1` does. Running (`wip/getrf`).
- [ ] **T37** Coverage follow-ups: let `survey.zig` read `axpy`'s hand ports (its `?` cells), give
  the `gemm` table the rstd column, keep its `arches` table in step with each registry's
  `usesVector`. Running (`wip/covfix`).

- [ ] **T57** Accelerate routing on arm64 loses in 37 of 228 flyoff rows (`results-flyoff.md`):
  every `nrm2` (`dznrm2` 0.06 to 0.21 of rstd's own), contiguous `sdot` through `cblas_dsdot` (0.41
  to 0.50), strided `copy` (0.52), strided complex `dot`, `axpy`, `rot`, `swap` (0.77 to 0.96): stop
  routing those to Accelerate. Running (`wip/acceltrim`).
- [ ] **R4** Mac timing rerun under the strict quiet rule, started automatically when the Mac is
  idle; other Mac runs wait while `scratchpad/MAC_TIMING_LOCK` exists (`wip/retime2`).
- [ ] **R5** Fact-check every external claim in the BLAS and LAPACK docs and comments: papers,
  upstream code, issue numbers, hardware facts and untraced numbers each get a real, durable
  citation, a correction from a source read, or a listing as unverified (scratchpad
  `factcheck_brief.md`). Running: `wip/fc-design` (design.md), `wip/fc-issues` (issues.md),
  `wip/fc-docs` (readme, threading, comparisons, notes, matlab README, results prose), `wip/fc-code`
  (comments and constants).

## Queue

- [ ] **T62** Main-merge gate: once T55 lands and R4 has run, one full rig sweep on review; report
  whether dev is clean to merge to main. **Needs:** T55, R4.
- [ ] **T65** The 32-bit x86 cores from SSE2 up (`northwood` to `bobcat`, readme.md, and 886324d's
  fallback objects) never run: the `i386` rig is `pentium3`. Add a freestanding `i386-max` rig
  (`-cpu max`; the TCG `fldt` bug is a paging case, issues.md), or a native `-Dtarget=x86-windows`
  run under WoW64 with `-Dblas-cores=all`, which also tests T03's x87 case. Part of T62's sweep.
- [ ] **T63** Docs still describe the routines as serial after the level-3 threading: design.md
  (`io` "unused until the routines thread", level 3 "Not ported: ... threading", "All routines are
  single-threaded"), readme.md ("Not yet: ... threading", "every routine is serial for now", and
  `GOTO_NUM_THREADS` "has no effect", against threading.md and design.md's fallback ceiling),
  comparisons.md's summary rows (7) and (10) still "Investigate", (9)'s open `params_data.zig`
  question (OpenBLAS `param.h:3574` gives VORTEX the 512/1024 set, which it holds), and
  `results-flyoff.md` ("rstd's routines run one thread"). design.md and issues.md also name
  0.18.0-dev.751 as the toolchain; it is dev.1864. Before T62.
- [ ] **T64** wasm runs levels 1 and 2 only (`build/wasmtest.zig` filters `blas1.`, `blas2.`,
  `simd.`): widen it to every OS-free test (level 3's tile, `ext/`, `_test/`, `from_openblas/`,
  `lapack/`), and add a build without `simd128` for the portable path.
- [ ] **T66** Zig 0.17.0 is a requirement nobody checks: zfls vendors rstd at `minimum_zig_version =
  "0.17.0"`, and rstd keeps 0.17 shims (`threading.zig`, `cpu/_corespec.zig`, the arm-softfp pin).
  State it in design.md's Requirements and add a compile check with the 0.17.0 toolchain (the VS
  Code extension's) for zfls's targets, linux-aarch64 at least, with `-Dblas-cores=all`.
- [ ] **T61** Uring (`Io.Evented`) hangs a plain futex-only barrier of four fibers once it has
  worker threads (`thread_limit = 3`; issues.md); the fiber test runs single-threaded Uring only.
  T02's `Io.Evented` tests start here.
- [ ] **T01** Threading levels 1 and 2 (threading.md, plan); level 3 landed (3b2caec to 05d9c59),
  all but `gemmt`, which stays serial as upstream.
- [ ] **T02** Threading tests beyond phase 2: a test that `Q` never depends on the thread count;
  race stress on aarch64, ppc64le and riscv64; `Io.Evented` past T61's hang, and its macOS
  `Dispatch`, which has no thread limit and is untested (threading.md); a multithreaded freestanding
  build with its own `Io` (the rigs are `single_threaded`); TSan with this toolchain. **Needs:** T61
  for the Uring part.
- [ ] **T56** Builds where `f32` runs on the x87 (no SSE, e.g. `-Dcpu=pentium2`) are untested; the
  half-precision `gemv` model test and the other `f32` bit promises are at risk there.
- [ ] **T05** Thresholds: OpenBLAS's defaults and per-core overrides, tuned by benchmark; then the
  overrides left out of phase 1 (z13 `axpy`, Windows-on-ARM64 `gemv`) if the numbers agree.
  **Needs:** T01, T73.
- [ ] **T06** `tools/blasbench` runs rstd with `Io.failing`; give it a threaded `Io`.
- [ ] **T07** Tie Accelerate into rstd's threading (threading.md, Vendor backends).
- [ ] **T08** LAPACK's own tests on the rigs: add `linalg/lapack` to `src/_test_qemu.zig`, with a
  sparc64 part.
- [ ] **T09** `linux-powerpc64le` (pwr8, `f128` by library calls) is not hit by the FastISel `stfd`
  bug: drop its `.safe`; the other POWER rigs keep it (T45).
- [ ] **T11** A rig for 32-bit POWER on a VSX CPU (QEMU's POWER8 or POWER9 model running a 32-bit
  image), so the power8 to 10 objects a 32-bit build links are tested (design.md, Settled decisions:
  filling upstream gaps).
- [ ] **T68** RVV runs at one VLEN: `riscv64` and `linux-riscv64` are `vlen=512`, so `zvl128b`'s
  VLEN-agnostic kernels (the `axpy` ports, the level-3 tile) never run at 128 or 256, the hardware
  they are named for. Add rig variants at both.
- [ ] **T18** The ISA check cannot see inline assembly (x86's assembler ignores CPU features).
  Detection, as OpenBLAS's, hands some cores to weaker CPUs: `katmai` to any 32-bit x86 without SSE3
  (down to the 486), `vortex` to any Apple CPU (pre-A12 lacks FCMLA), `vortexm4` to non-Apple SME
  without SVE, a forced core's fallback to `neoversen1` on ARMv8.0; decide per case.
- [ ] **T19** `XDouble` as `f128`: `quad_precision` is hard-wired off; add a build option and a rig
  configuration for it. The edge tests and level 3's per-routine driver tests still cover no
  extended type off x86 (the interfaces do).
- [ ] **T30** Level 3: `gemm3m`, `symm3m`, `hemm3m` (stubs kept), a small-matrix path,
  half-precision level 3 (`shgemm`, whose vtable slots are all null), bfloat16 `sbgemm` and
  `sbgemmt` (OpenBLAS's bfloat16 `gemmt`), then `bgemm` (design.md, Implementation plan).
- [ ] **T69** OpenBLAS's comparison tests are not ported: `compare_sgemv_{bgemv,sbgemv,shgemv}.c`
  cover routines rstd has (`gemvHalf`); `compare_sgemm_{bgemm,sbgemm,shgemm}.c` come with T30.
  `blat3.zig` runs `[cz]blat3_3m` through `gemm`: switch it to the real `gemm3m`, with 3m's looser
  imaginary-part bound (comparisons.md (7)). **Needs:** T30 for the level-3 parts.
- [ ] **T73** Benchmarks above level 1 race `gemv` and `gemm` only: add `ger`, `symv`, `trsv`, a
  banded routine, `symm`, `syrk`, `trmm` and `trsm`, which T05 and the speed claims need. Rerun the
  Zen 3 flyoff (Windows and Linux), which predates T43's x86_64 `rot`, `scal`, `gemv` paths and the
  level-3 threading, under a timing rule for that host. **Needs:** T06 for threaded rows.
- [ ] **T20** Kernel strategy: race `@Vector` against the converted assembly per core on a quiet
  machine (only vortex is done, `tools/blasbench/results-m6.md`); keep per-core code only where Zig
  loses.
- [ ] **T21** Zen 3 (flyoff): short real `dot`s (0.30 of OpenBLAS at 1024) and complex `gemm` (0.79
  to 1.01) trail OpenBLAS.
- [ ] **T67** Level 2's `coverage.md` has no rstd column, so T22's level-2 gaps are invisible: make
  `survey.zig l2` print what rstd serves, as levels 1 and 3 do. `survey.zig` also marks the BG/P
  `*_hummer.S` sources as gaps, though LLVM cannot assemble them; list them as unsupported, as
  `tools/blas_asm/arch/power.zig` does `axpy_hummer.S`. **Needs:** T37.
- [ ] **T22** C and intrinsics kernels not yet ported, levels 1 and 2, by T20's strategy, `@Vector`
  first as arm64 did (`vector.zig`); the coverage.md gaps, the largest first. Level 1: riscv64
  outside `axpy`; POWER8 to 10 `asum`, `scal`, `rot`, `swap`, `cdot`; MIPS MSA and zarch; x86_64
  `asum` on Haswell and AVX-512; wasm SIMD128. Level 2: x86_64 `*_4.c` gemv, symv, `sger.c`,
  `dger.c`; riscv64 `gemv_[nt]_rvv.c`, `symv_[UL]_rvv.c`, `zsymv`, `zhemv`; MIPS MSA
  `?gemv_[nt]_msa.c` and `gemv_[nt]_loongson3a.c`; zarch z13/z14 `?gemv_[nt]_4.c`; arm64
  `symv_[UL]_asimd_4x4.c`, SVE `symv_[UL]_sve_v1x4.c` and `ger_sve_v1x3.c`; wasm `gemv_[nt].c`.
  **Needs:** T20, T67.
- [ ] **T23** riscv `zvl256b` has no object (`zvl128b` serves it): read VLEN and port the
  `*_vector.c` kernels (design.md, Open issues). **Needs:** T20.
- [ ] **T24** Skylake-X family `sgemv_n`/`sgemv_t`: those cores add large intrinsics kernels on top
  of Haswell's (`sgemv_t` uses a 1,121-line template); they keep the portable `sgemv` until ported.
  Testable only on AVX-512 hardware (T43).
- [ ] **T26** vortex: a tuned `@Vector` `rot` (the assembly still leads by 2 to 17% in L1); race the
  Zig complex `gemm` against OpenBLAS's vortex assembly on the M6.
- [ ] **T27** gemm tuning: faster Zig packing (whole-vector transposes in `oncopy`, `incopy`),
  `C`-tile prefetch at kernel entry, a larger `KC` on Apple cores (comparisons.md (9)).
- [ ] **T28** Measure whether composed `symv`/`hemv` pay off on cores with assembly gemv (atom,
  opteron, arm64); today they compose only where the core has its own gemv.
- [ ] **T29** The bfloat16 intrinsics (cooperlake, sapphirerapids, Neoverse N2, V1, V2, POWER8 to
  10, riscv zvl256b); `sbdot` could call the core's own `sdot`. **Needs:** T20.
- [ ] **T31** Assembly where `@Vector` cannot reach: arm64 SVE, SME (`vortexm4`), POWER10 MMA, AMX
  and AVX-512 BF16 for `sbgemm`; Skylake-X only if under 85% of its kernel (design.md,
  Implementation plan; `blas3/gemm/notes.md`). **Needs:** T20.
- [ ] **T32** BLIS assembly kernels, decided per kernel where benchmarks favor them (BSD-3, with a
  license section and per-file provenance). **Needs:** T20.
- [ ] **T70** OpenBLAS's own LAPACK past `potrf`, `potrs`, `laswp` and T33's `getrf`: `getf2`,
  `getrs`, `gesv`, `potf2`, `potri`, `lauu2`/`lauum`, `trti2`/`trtri`, `trtrs` (threaded variants
  with T01) and `laed3` with `test_laed3_nan.c` (`interface/lapack/`, `lapack/`). M2 needs `getrs`
  at least. **Needs:** T33.
- [ ] **T34** LAPACK `dgesvd` (for `test_kernel_regress.c`): about 5,000 to 8,000 lines of reference
  LAPACK.
- [ ] **T35** Optionally route `potrf` to Accelerate on macOS; threaded `potrf`
  (`potrf_parallel.c`). **Needs:** T01.
- [ ] **T36** The Accelerate level-2 special-value sweep takes minutes per Mac run (complex, even
  trimmed to every third variant at the largest shape); speed it up.
- [ ] **T38** Check in the C-against-Zig differential harness (stub `common.h`) used for `dot`.
- [ ] **T39** `sparc64`: replace the image split with a real fix, mapping more memory at boot
  through the Open Firmware client interface (claim and map), so one image fits.
- [ ] **T40** Linux rig gaps: `linux-powerpc64le` never runs POWER10 kernels (baseline build),
  `linux-aarch64` leaves out the SVE cores, and under HVF it also skips `vortexm4` (cause unknown).
- [ ] **T71** Linux rigs with no item: loongarch64 (the `AT_HWCAP` LSX/LASX probe), big-endian
  powerpc64 (`AT_PLATFORM`, the `64be` variant), sparc64, 32-bit mips.
- [ ] **T72** UEFI never runs the tests: `rstd-check-uefi` only links. Add a run rig for the OS-free
  suite under OVMF and AAVMF (UEFI's inline `Io`; RandyOS boots through it).
- [ ] **T74** No CI runs a BLAS test: rstd's `.gitlab-ci.yml` is rgui's copy (`EXE_NAME: "rgui"`,
  Zig 0.16.0). A job for `rstd-test`, `rstd-check-macho`, `rstd-check-uefi` and T66's 0.17 check
  would catch drift between rig sweeps; optional if testing stays manual.
- [ ] **T41** `zrot_ppc440.S` (`ppc440fp2`) has a zero-stride patch no rig has run; no rig runs
  `ppc440` (it shares code with `ppcg4`). The 32-bit `ppc440` and `ppc440fp2` kernels rebuilt in
  719a41b are compile-checked only.
- [ ] **T42** Test-only code still auto-vectorizes on the G4 (`blat1` rot, the scal, gemv and
  half-gemv models, `triangular.checkAll`, `blat3`'s complex type). Guard with `util.scalar` the
  loops a test feeds subnormals (AltiVec flushes them); close the rest, which see normal-range data
  only.
- [ ] **T43** Cores no rig runs: the Skylake-X family (AVX-512) and the Bulldozer family (FMA4,
  XOP), including their new x86_64 `rot`, `scal` and `gemv` paths, which are compile-only: run
  `x86_64-whpx` on a host that has them. **Blocked:** hardware.
- [ ] **T44** The AVX-512-off workaround for the self-hosted backend (`vectorBits` in
  `blas3/gemm/kernel.zig`, `zig_backend != .stage2_x86_64`) has no issues.md entry and no recorded
  reason (since 2cf8757); 0.18.0-dev.751 compiles 512-bit `@Vector` multiply-adds there. Run the
  gemm tests with `-fno-llvm` on AVX-512 hardware, then drop it or record it. **Blocked:** hardware.
- [ ] **T45** POWER rigs stay at `.safe` (dcb4628): LLVM's FastISel stores an `f64` held in a
  v-register with `stfd` (issues.md); drop `.safe` once fixed. **Blocked:** upstream.
- [ ] **T46** No `powerpc64` (970) AltiVec scan: LLVM stops with "Invalid global physical register".
  Scan once it builds; rerun the G4 `-femit-asm` scan after LLVM upgrades, since the SLP vectorizer
  is not proven blocked by `util.scalar`. **Blocked:** upstream.
- [ ] **T47** Hardware-only paths never run: the mips `cpucfg` decode, the sparc `sllx` rewrite,
  Linux ppc64 `mfpvr`, the Windows arm64 registry MIDR and feature probes, the arm64 sysfs fallback,
  the FreeBSD arm64 probe, Accelerate's pre-13.3 symbol fallbacks (`$NEWLAPACK`, classic `cblas_*`);
  NEON and AltiVec flushing seen on rigs only. **Blocked:** hardware.
- [ ] **T48** 32-bit SPARC never runs (lld cannot link it; includes the `%icc` index-tie fix).
  **Blocked:** upstream.
- [ ] **T49** Linux rigs: `linux-i386` is a stub (QEMU x87 page-crossing bug); `s390x` and
  `mips64el` wait on `zig fetch` taking single files. **Blocked:** upstream.
- [ ] **T50** The `powerpc-macho` rig is a stub; finish it once the powerpc-macos Zig port lands.
  **Blocked:** upstream.
- [ ] **T51** translate-c: re-enable the package once its bundled Aro builds on current Zig.
  **Blocked:** upstream.
- [ ] **T52** On each Zig update, recheck and drop the workarounds once fixed (T48 and T49 track
  lld's SPARC and `zig fetch`): the powerpc64le two-float return (`zdot` patch), arm soft-float
  struct returns, `f80` under `x86_64_sysv` on Windows, the comptime naked-function call, ppc64
  big-endian ELFv2, the two-lane `f32` vector stored as 16 bytes, sparc64 `_Qp_sqrt`
  (`check.sqrtW`), sparc64 `-Oz`, `SafeAllocator`'s CPU count, the `Thread` override doc, run-step
  timeouts and output, the interrupted-link cache entry, inline-asm errors exiting 0
  (`checkedObject`), the MIPS64r6 select copy and spill (comptime `conj`), MIPS MSA vector `fma`
  (MSA cores differ from the reference's fused bits until fixed), arm CPU models without an FPU (rig
  pins), WHPX passing AMD CPUID (macOS rig on TCG), QEMU `xststdcsp` (portable `nrm2` avoids
  `isInf`). The Zig entries were last checked on dev.751.
- [ ] **T53** macOS VM extras: base-image overlays (`qemu-clone`, `-Dmac-base`), a `Forgejo`
  self-hosted runner, a setup that needs no clicking, less guest load when headless (HID driver,
  `SecurityAgent`), and freeze detection that asks before killing.
- [ ] **T54** Optional: an inline-asm `maddf` for R6 scalar multiply-adds in the main library, if R6
  speed ever matters (LLVM calls `fmaf`).
- [ ] **T75** Optional: `?gemm_batch` and `?gemm_batch_strided` (`interface/gemm_batch.c`,
  `gemm_batch_strided.c`) and the `sbgemm`/`bgemm` batch forms, column-major only; a strided batch
  maps onto M3's fixed-size forms.
- [ ] **T76** Optional, from comparisons.md: `Dot2` (compensated `f32` lanes) if T25 finds the `f64`
  `sdot` costly in cache; BLIS's `JR`/`IR` ways for threads sharing an L2; an A-panel prefetch
  beside T27's `C`-tile one; small-`n` `sgemv_n` tuning (`results-m6.md`).

## After the BLAS is feature complete

- [ ] **M1** Audit `rstd.linalg.matlab` (the matlab88 port of 1980s LINPACK and EISPACK,
  `src/linalg/matlab/`) against the modern BLAS and LAPACK (design.md, Settled decisions: the matlab
  port defers): `blas1.zig`'s level-1 kernels and `@Vector` stride-1 paths, `matmul.zig`, and what
  LU, QR, Cholesky, SVD and the eigen routines call. The default is the modern routine; list only
  where the old one is more accurate or right on an edge case the modern one gets wrong, each with a
  test that shows it. **Needs:** the queue above.
- [ ] **M2** Replace the matlab port's own BLAS and LAPACK-level code with `rstd.linalg.blas` and
  `rstd.linalg.lapack`, keeping only what M1 showed better. Decide first how to bridge its split
  real/imaginary complex layout (`xr`/`xi`) to the BLAS's interleaved `Complex`: convert at the
  boundary or change the port's layout. Its LINPACK and EISPACK routines map to LAPACK's `getrf`,
  `getrs`, `getri`, `gecon` (`wgefa`, `wgesl`, `wgedi`, `wgeco`), `potrf` (`wpofa`), `geqp3`,
  `unmqr`, `trtrs` (`wqrdc`, `wqrsl`), `gesvd` (`wsvdc`, T34), `gehrd` (`corth`), `hseqr` with
  `trevc` (`comqr3`), `hetrd` and `unmtr` (`htridi`, `htribk`); rstd has `potrf`, T33 and T70 add
  OpenBLAS's own, and the rest come from Reference-LAPACK (BSD-3, as OpenBLAS ships it in
  `lapack-netlib/`). Its fills (`wset`, `rset`) are LAPACK `laset`, `pythag` is `lapy2`. Its own
  tests must still pass. **Needs:** M1.

- [ ] **M3** Build `rstd.linalg`'s public layer, the BLAS behind it as an implementation detail
  (design.md, Settled decisions): `linalg.vector` takes fixed-size vectors and small matrices by
  value, with math names, and per size runs inline `@Vector` code or calls the BLAS (crossover from
  M4); under SPIR-V each operation is SPIR-V's own instruction (`OpDot`, `OpMatrixTimesVector`,
  GLSL.std.450). Split `rstd.simd` between them: slice reductions become BLAS routines, element-wise
  helpers (`mix`, `mulAdd`) move to `linalg.vector`; move its callers (the matlab port, the rgl font
  code, shaderlib's `mix` and `sum`, `rgl/shaders/types.zig`'s `zigvec` types).
  `arrayMax`/`arrayMin` use `@reduce`'s NaN rule, not the BLAS's first-NaN one, so moving them
  changes their special-value results; `arrayProduct`, the bitwise reductions and general
  `mapReduce` have no BLAS counterpart. Then grow it toward any-size vectors and matrices, as MATLAB
  takes them. **Needs:** M1, S1.
- [ ] **M4** Time `linalg.vector`'s inline `@Vector` forms against the BLAS calls per size (2, 3, 4,
  then up through sizes such as 28 where either could win) and operation (`dot`, `nrm2`-style
  `length`, `axpy`-style updates, matrix-vector), on the M6 and Zen 3 under their timing rules; set
  M3's crossover from the numbers. **Needs:** M3's first forms.

## SPIR-V: vector math in shaders

Not the BLAS kernels: the fixed-size vector math (`rstd.linalg.vector`, `rstd.simd` until M3) that
`src/rgl/shaders/lib/` (`math.zig`, `sdf.zig`) builds on. Can start any time. SPIR-V's own
instructions are used as they are, not reimplemented or tuned: the aim is one API, not speed.

- [ ] **S1** Audit what the SPIR-V backend emits for shaderlib's math: dot, cross, length, distance,
  normalize, reflect, refract, faceforward, inversesqrt, fma, mix, step, smoothstep, fract, clamp,
  min/max, matrix-vector and matrix-matrix products. For each, whether it lowers to the native
  instruction (`OpDot`, `OpVectorTimesMatrix`, `OpMatrixTimesVector`, or the GLSL.std.450 extended
  set) or to scalar code; check by disassembling the compiled shaders. Known already:
  `blas/fma.zig`'s `hardware` is false on SPIR-V, so `mulAdd` (and `mix`, `linalg.vector.dot`) is an
  unfused `a*b + c`, never `Fma`; `linalg.vector.dot` is a scalar `mulAdd` chain by design, `length`
  is `@sqrt(dot2)`, not `Length`, and `norm` is not `Normalize`. Shaderlib's `Mat2` to `Mat4` are
  arrays of vectors, not `OpTypeMatrix`, so the matrix instructions need a matrix type first; record
  which.
- [ ] **S2** Shaderlib stays the shaders' one API (design.md, Settled decisions): GLSL-friendly
  names that alias or wrap `linalg` calls where `linalg` has the operation (`dot`, `length`,
  matrix-vector products, through M3's SPIR-V forms); `cross` moves to `rstd.linalg.vector` and is
  re-exported; `normalize`, `reflect`, `refract`, `faceforward`, `smoothstep`, `fract`, `step` stay
  in shaderlib, lowered to GLSL.std.450 instructions. Add a build-time test that disassembles the
  shaders and fails if a listed operation comes out as scalar code. **Needs:** S1, M3 (for the
  `linalg`-backed ones).

## Zig-idiomatic interface, last

Decided, not built (design.md, "Planned: a Zig-idiomatic interface").

- [ ] **Z1** Slices for every operand instead of pointer and size arguments; routines check that the
  slices cover what sizes and strides need.
- [ ] **Z2** Exact error sets per routine or family, with plain names (`SliceTooShort`,
  `ZeroStride`, `LeadingDimensionTooSmall`); `Canceled` only if cancelable, `OutOfMemory` only if
  allocating.
- [ ] **Z3** Errors, not panics, for bad input; checks no longer need the reference's order.
- [ ] **Z4** Doc comments say where a routine panics explicitly (only unrecoverable states, above
  all possible memory corruption from assembly).
- [ ] **Z5** User docs without BLAS jargon (`xerbla`, `LDA`, `INCX`).
- [ ] **Z6** The last item: an rstd-wide Zig-idiomatic pass (slices, exact error sets, plain error
  names, documented panics). **Needs:** Z1 to Z5.

## Waiting on the user

Nothing waiting.
