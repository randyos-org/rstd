# BLAS to-do

The BLAS backlog in priority order: work goes top to bottom. IDs stay fixed when items move.
**Needs:** names items that must land first; **Blocked:** names an outside cause (hardware, an
upstream fix, a decision). Items are removed once done and recorded (commits, design.md,
issues.md). Upstream problems go to issues.md as they are found; the user triages them after the
BLAS work.

## Running now

- [ ] **T03** Workers adopt the caller's floating-point control state (design.md, Settled
  decisions): capture it on the caller in `threading.team` and `forEachIn`, apply it as each task's
  first step (x87 control word and MXCSR on x86, FPCR on arm64, the others' equivalents). Then a
  `-Dcpu=pentium3` build's threaded results match serial ones (on Windows new threads start at
  53-bit x87 precision), and so do callers who change rounding or flush-to-zero. Running
  (`wip/fpstate`).
- [ ] **T25** Speed follow-ups: a `@Vector` `sum` (`csum`/`zsum` are scalar on the old x86_64
  cores); measure the scalar fallbacks (extremes and `sdot` on arm, pre-VSX POWER, mips64, SPARC)
  and the `f64` `sdot` (comparisons.md (1)). Running (`wip/vecsum`).
- [ ] **T55** The full i386 rig (the trimmed 32-bit objects, 886324d; the x87 audit's checks) stalls
  in "interface: symm and hemm match the wide reference" and the watchdog kills it, after the
  level-3 threading landed; investigating (`wip/i386stall`).
- [ ] **T58** Run the edge tests on the Mac (`ssh-test-aarch64-macos`): they now put a zero factor
  beside NaN and infinity for every rank update, so they show whether Accelerate's `syr`, `her`,
  `spr`, `hpr` and rank-2 routines skip as the reference does; if not, guard them as `ger` is in
  `apple_accel.zig`. Accelerate's `spr2` (and likely `her2`, `hpr2`) fail it with NaN; guarding
  (`wip/accelrank`).

## Queue

- [ ] **T61** Uring (`Io.Evented`) hangs a plain futex-only barrier of four fibers once it has
  worker threads (`thread_limit = 3`; issues.md); the fiber test runs single-threaded Uring only.
- [ ] **T57** Accelerate routing on arm64 loses in 37 of 228 flyoff rows (`results-flyoff.md`):
  every `nrm2` (`dznrm2` 0.06 to 0.21 of rstd's own), contiguous `sdot` through `cblas_dsdot` (0.41
  to 0.50), strided `copy` (0.52), strided complex `dot`, `axpy`, `rot`, `swap` (0.77 to 0.96): stop
  routing those to Accelerate.
- [ ] **T01** Threading levels 1 and 2 (threading.md, plan); level 3 landed (3b2caec to 05d9c59),
  all but `gemmt`, which stays serial as upstream.
- [ ] **T02** Threading tests beyond phase 2: a test that `Q` never depends on the thread count;
  race stress on aarch64, ppc64le and riscv64; the futex protocol on `Io.Evented`; TSan with this
  toolchain.
- [ ] **T56** Builds where `f32` runs on the x87 (no SSE, e.g. `-Dcpu=pentium2`) are untested; the
  half-precision `gemv` model test and the other `f32` bit promises are at risk there.
- [ ] **T05** Thresholds: OpenBLAS's defaults and per-core overrides, tuned by benchmark; then the
  overrides left out of phase 1 (z13 `axpy`, Windows-on-ARM64 `gemv`) if the numbers agree.
  **Needs:** T01.
- [ ] **T06** `tools/blasbench` runs rstd with `Io.failing`; give it a threaded `Io`.
- [ ] **T07** Tie Accelerate into rstd's threading (threading.md, Vendor backends).
- [ ] **T08** LAPACK's own tests on the rigs: add `linalg/lapack` to `src/_test_qemu.zig`, with a
  sparc64 part.
- [ ] **T09** `linux-powerpc64le` (pwr8, `f128` by library calls) is not hit by the FastISel `stfd`
  bug: drop its `.safe`; the other POWER rigs keep it (T45).
- [ ] **T11** A rig for 32-bit POWER on a VSX CPU (QEMU's POWER8 or POWER9 model running a 32-bit
  image), so the power8 to 10 objects a 32-bit build links are tested (design.md, Settled decisions:
  filling upstream gaps).
- [ ] **T18** The ISA check cannot see inline assembly (x86's assembler ignores CPU features).
  Detection, as OpenBLAS's, hands some cores to weaker CPUs: `katmai` to any 32-bit x86 without SSE3
  (down to the 486), `vortex` to any Apple CPU (pre-A12 lacks FCMLA), `vortexm4` to non-Apple SME
  without SVE, a forced core's fallback to `neoversen1` on ARMv8.0; decide per case.
- [ ] **T19** `XDouble` as `f128`: `quad_precision` is hard-wired off; add a build option and a rig
  configuration for it. The edge tests and level 3's per-routine driver tests still cover no
  extended type off x86 (the interfaces do).
- [ ] **T20** Kernel strategy: race `@Vector` against the converted assembly per core on a quiet
  machine (only vortex is done, `tools/blasbench/results-m6.md`); keep per-core code only where Zig
  loses.
- [ ] **T21** Zen 3 (flyoff): short real `dot`s (0.30 of OpenBLAS at 1024) and complex `gemm` (0.79
  to 1.01) trail OpenBLAS.
- [ ] **T22** C and intrinsics kernels not yet ported, levels 1 and 2, by T20's strategy, `@Vector`
  first as arm64 did (`vector.zig`); the coverage.md gaps, largest first: riscv64 level 1 outside
  `axpy`; POWER8 to 10 `asum`, `scal`, `rot`, `swap`, `cdot`; MIPS MSA and zarch level 1; arm64
  `sum`; x86_64 `asum` on Haswell and AVX-512, `*_4.c` gemv, symv and ger; wasm SIMD128. **Needs:**
  T20.
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
- [ ] **T30** Level 3: `gemm3m`, `symm3m`, `hemm3m` (stubs kept), a small-matrix path, bfloat16 and
  extended level 3, then `bgemm` (design.md, Implementation plan).
- [ ] **T31** Assembly where `@Vector` cannot reach: arm64 SVE, SME (`vortexm4`), POWER10 MMA, AMX
  and AVX-512 BF16 for `sbgemm`; Skylake-X only if under 85% of its kernel (design.md,
  Implementation plan; `blas3/gemm/notes.md`). **Needs:** T20.
- [ ] **T32** BLIS assembly kernels, decided per kernel where benchmarks favor them (BSD-3, with a
  license section and per-file provenance). **Needs:** T20.
- [ ] **T33** LAPACK `getrf`, then call `laswp_ncopy` and `neg_tcopy` from it, or drop them as
  upstream's `#if 1` does.
- [ ] **T34** LAPACK `dgesvd` (for `test_kernel_regress.c`): about 5,000 to 8,000 lines of reference
  LAPACK.
- [ ] **T35** Optionally route `potrf` to Accelerate on macOS; threaded `potrf`
  (`potrf_parallel.c`). **Needs:** T01.
- [ ] **T36** The Accelerate level-2 special-value sweep takes minutes per Mac run (complex, even
  trimmed to every third variant at the largest shape); speed it up.
- [ ] **T37** Coverage follow-ups: let `survey.zig` read `axpy`'s hand ports (its `?` cells), give
  the `gemm` table the rstd column, keep its `arches` table in step with each registry's
  `usesVector`.
- [ ] **T38** Check in the C-against-Zig differential harness (stub `common.h`) used for `dot`.
- [ ] **T39** `sparc64`: replace the image split with a real fix, mapping more memory at boot
  through the Open Firmware client interface (claim and map), so one image fits.
- [ ] **T40** Linux rig gaps: `linux-powerpc64le` never runs POWER10 kernels (baseline build),
  `linux-aarch64` leaves out the SVE cores, and under HVF it also skips `vortexm4` (cause unknown).
- [ ] **T41** `zrot_ppc440.S` (`ppc440fp2`) has a zero-stride patch no rig has run; no rig runs
  `ppc440` (it shares code with `ppcg4`). The 32-bit `ppc440` and `ppc440fp2` kernels rebuilt in
  719a41b are compile-checked only.
- [ ] **T42** Test-only code still auto-vectorizes on the G4 (`blat1` rot, the scal, gemv and
  half-gemv models, `triangular.checkAll`, `blat3`'s complex type); fine for normal-range data.
- [ ] **T43** Cores no rig runs: the Skylake-X family (AVX-512) and the Bulldozer family (FMA4,
  XOP), including their new x86_64 `rot`, `scal` and `gemv` paths, which are compile-only: run
  `x86_64-whpx` on a host that has them. **Blocked:** hardware.
- [ ] **T44** The AVX-512-off workaround for the self-hosted backend (`blas3/gemm/kernel.zig:56`)
  has no issues.md entry and no recorded reason (since 2cf8757); 0.18.0-dev.751 compiles 512-bit
  `@Vector` multiply-adds there. Run the gemm tests with `-fno-llvm` on AVX-512 hardware, then drop
  it or record it. **Blocked:** hardware.
- [ ] **T45** POWER rigs stay at `.safe` (dcb4628): LLVM's FastISel stores an `f64` held in a
  v-register with `stfd` (issues.md); drop `.safe` once fixed. **Blocked:** upstream.
- [ ] **T46** No `powerpc64` (970) AltiVec scan: LLVM stops with "Invalid global physical register".
  Scan once it builds; rerun the G4 `-femit-asm` scan after LLVM upgrades, since the SLP vectorizer
  is not proven blocked by `util.scalar`. **Blocked:** upstream.
- [ ] **T47** Hardware-only paths never run: the mips `cpucfg` decode, the sparc `sllx` rewrite,
  Linux ppc64 `mfpvr`, the Windows arm64 registry MIDR and feature probes, the arm64 sysfs fallback;
  NEON and AltiVec flushing seen on rigs only. **Blocked:** hardware.
- [ ] **T48** 32-bit SPARC never runs (lld cannot link it; includes the `%icc` index-tie fix).
  **Blocked:** upstream.
- [ ] **T49** Linux rigs: `linux-i386` is a stub (QEMU x87 page-crossing bug); `s390x` and
  `mips64el` wait on `zig fetch` taking single files. **Blocked:** upstream.
- [ ] **T50** The `powerpc-macho` rig is a stub; finish it once the powerpc-macos Zig port lands.
  **Blocked:** upstream.
- [ ] **T51** translate-c: re-enable the package once its bundled Aro builds on current Zig.
  **Blocked:** upstream.
- [ ] **T52** On each Zig update, recheck and drop the workarounds once fixed: the powerpc64le
  two-float return (`zdot` patch), arm soft-float struct returns, `f80` under `x86_64_sysv` on
  Windows, the comptime naked-function call, ppc64 big-endian ELFv2, the two-lane `f32` vector
  stored as 16 bytes, sparc64 `_Qp_sqrt` (`check.sqrtW`), sparc64 `-Oz`, lld and 32-bit SPARC,
  `SafeAllocator`'s CPU count, the `Thread` override doc, `zig fetch`, run-step timeouts and output,
  the interrupted-link cache entry.
- [ ] **T53** macOS VM extras: base-image overlays (`qemu-clone`, `-Dmac-base`), a `Forgejo`
  self-hosted runner, a setup that needs no clicking, less guest load when headless (HID driver,
  `SecurityAgent`), and freeze detection that asks before killing.
- [ ] **T54** Optional: an inline-asm `maddf` for R6 scalar multiply-adds in the main library, if R6
  speed ever matters (LLVM calls `fmaf`).

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
