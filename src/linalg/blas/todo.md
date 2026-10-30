# BLAS to-do

The BLAS backlog in priority order: work goes top to bottom. IDs stay fixed when items move.
**Needs:** names items that must land first; **Blocked:** names an outside cause (hardware, an
upstream fix, a decision). Upstream problems go to issues.md as they are found; the user
triages them after the BLAS work.

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

- [ ] **T22** C and intrinsics kernels not yet ported, levels 1 and 2, by T20's strategy, `@Vector`
  first as arm64 did (`vector.zig`); the coverage.md gaps, the largest first. Level 1: riscv64
  outside `axpy`; POWER8 to 10 `asum`, `scal`, `rot`, `swap`, `cdot`; MIPS MSA and zarch; x86_64
  `asum` on Haswell and AVX-512; wasm SIMD128. Level 2 (464 gaps, `blas2/coverage.md`): riscv64 140
  (gemv 82, hemv 32, symv 24, `sbgemv` and `shgemv` on zvl256b); mips64 100 (gemv on I6400, I6500,
  P6600, Loongson 3R3, 3R4); arm64 60 (symv on 11 cores, `sger`/`dger` on A64FX and Neoverse V1,
  `sbgemv`/`bgemv` on N2, V1, V2); POWER 60 (gemv on POWER8 to 10, `sbgemv` on POWER10); x86_64 50
  (symv on 10 cores, `sgemv` on the Skylake-X family (T24), `sger`/`dger` on Sandy Bridge,
  `sbgemv_n` on Cooper Lake and Sapphire Rapids); zarch 30 (gemv on z13, z14); mips 20 (gemv on
  P5600); wasm 4 (gemv). **Needs:** T20. Running (`wip/portrv`: riscv64; `wip/portxa`: x86_64 and
  arm64). riscv64 landed 2026-10-30 (f7b6be16): `@Vector` kernels fill every slot on
  `riscv64_zvl128b`, `riscv64_zvl256b` and `x280` (`vector.Rule.unfilled_on`); level-2 gaps 140 to
  32, all C910V, which is **Blocked** (no `xtheadvector` in LLVM or Zig).
- [ ] **T99** Multithreaded freestanding tests (from T02, see threading.md): start secondary cores
  in the test kernels with their own stacks, give each core thread-local storage (std's `debug.zig`
  makes `panic_stage` thread-local once a build is not `single_threaded`), an `Io` that hands group
  tasks to those cores with futex waits, wakes and cancellation, and start the rigs with more than
  one CPU. Running (`wip/smp`).
- [ ] **T73** Benchmarks above level 1 race `gemv` and `gemm` only: add `ger`, `symv`, `trsv`, a
  banded routine, `symm`, `syrk`, `trmm` and `trsm`, which T05 and the speed claims need. Rerun the
  Zen 3 flyoff (Windows and Linux), which predates T43's x86_64 `rot`, `scal`, `gemv` paths and the
  level-3 threading, under a timing rule for that host. **Needs:** T06 for threaded rows. Reference
  binaries beside OpenBLAS (design.md, the anthology decision): reference LAPACK, and Eigen (MPL
  2.0) built from upstream through Nix as a benchmark-only binary, never ported, for the LAPACK and
  level-3 races. Running (`wip/bench`).
- [ ] **T24** Skylake-X family `sgemv_n`/`sgemv_t`: those cores add large intrinsics kernels on top
  of Haswell's (`sgemv_t` uses a 1,121-line template); they keep the portable `sgemv` until ported.
  Testable only on AVX-512 hardware (T43). Running (`wip/portxa`).
- [ ] **T111** The rest of LAPACK that would otherwise come from reference LAPACK, ported from faer
  (Q10): QR (`geqrf`, `orgqr`/`ungqr`, `ormqr`), least squares (`gels`), the symmetric and Hermitian
  eigensolvers (`syevd`/`heevd`), the general eigensolver (`geev`), `laed3` and OpenBLAS's
  `test_laed3_nan.c` (T70 left them: they need reference LAPACK's `dlaed4`, `dlaed5`, `dlaed6` and
  the test drives `dsyevd`; take faer's divide-and-conquer tridiagonal solver), and whatever else
  rstd.linalg's MATLAB-like API needs. Reference LAPACK comes later only as a benchmark baseline.
  Reuse T34's `lapack/faer/*.zig` (`householder`, `qr`, `bidiag`, `ops`) and add files beside them;
  `qr.zig` already has faer's QR for `geqrf`. Running (`wip/faer2`).
- [ ] **T19** `XDouble` as `f128`: `quad_precision` is hard-wired off; add a build option and a rig
  configuration for it. The edge tests and level 3's per-routine driver tests still cover no
  extended type off x86 (the interfaces do). Running (`wip/quad`).
- [ ] **T116** The generic family follows `parameter.c:451`, an x86-only rule: a static
  `TARGET=GENERIC` OpenBLAS elsewhere gets `p = q = r = 128` (`parameter.c:46-48`), or fixed sets on
  MIPS and LoongArch (`param.h:4386-4415`). Make rstd's generic `p`/`q`/`r` follow the non-x86 rule
  off x86 (from T112). Running (`wip/quad`).
- [ ] **T40** Linux rig gaps: `linux-powerpc64le` never runs POWER10 kernels (baseline build),
  `linux-aarch64` leaves out the SVE cores, and under HVF it also skips `vortexm4` (cause unknown).
  Running (`wip/riggaps`).
- [ ] **T115** Routines only reference LAPACK has (T113, `lapack/sources.md`): `gecon`,
  `gebal`/`gebak`, `trexc`/`trsen`, `gelsy`, plus `lange` (no source chosen) and `gelqf` (only
  libflame's FLAME code besides the reference). The anthology rule picks the best open-source
  implementation, and here that is the reference: port these from reference LAPACK (its license is
  in the accepted BSD class) unless faer turns out to have an equivalent, and record each in the
  table. Running (`wip/refonly`).

## Queue

- [ ] **T62** Main-merge gate: once T55 and T77 land and R4 has run, one full rig sweep on review,
  plus `wasm-test`, `rstd-check-macho`, `rstd-check-uefi` and `rstd-check-zig017` (Zig 0.17.0,
  `-j4`); report whether dev is clean to merge to main. **Needs:** the items running now to land, so
  one sweep covers them (T55, T77, T91 and R4 done).
- [ ] **T08** LAPACK's own tests (`potrf`, `getrf`, `getrs`) on the rigs: add `linalg/lapack` to
  `src/_test_qemu.zig`, with a sparc64 part. Landed (8c574827): `linalg.lapack.` in the rigs'
  default filters and `_test_qemu.zig` (so `rstd-check-uefi` links it), sparc64's LAPACK in the
  `blas3-tri` part, sparc64-debug in four `lapack-*` parts (`getrs`'s every-core test split into
  real and complex to fit the load window; `getrf` and `cholesky` have about 140 and 240 KiB left);
  every listed `qemu-check` passes. Check off once T62's sweep runs them; the 32 MiB rig heap is
  unproven.
- [ ] **T114** Accelerate's `potrf` (T35) is untested on the Mac: which symbol spelling Apple
  exports (`load` tries `dpotrf$NEWLAPACK$ILP64` and `dpotrf_$NEWLAPACK$ILP64`), whether its special
  values match the reference (a quirk entry if not), and in the timing session whether routing beats
  rstd's threaded `potrf`. Sweep filters: `potrf goes to Accelerate`, `potr`.
- [ ] **T117** faer follow-ups (T34): `@Vector` the fused bidiagonalization pass (scalar now, faer's
  four partial sums); keep `gesvd`'s `INFO > 0` count (faer doesn't); run divide and conquer's two
  halves in parallel if rstd's threading gains task trees; race `gesvd`/`gesdd` against OpenBLAS's
  LAPACK, reference LAPACK and Eigen in the timing session (T73's bench). Sweep: watch
  sparc64-debug's `openblas-l1` part (0x3f3800), which now also pulls in `gesvd` (T39).
- [ ] **T118** ReLAPACK follow-ups (T110): `tgsyl`, which needs reference LAPACK's `xTGSY2` (2x2
  blocks, up to 8x8 systems), `xGETC2`, `xGESC2` and `xLATDF` (about 1,500 lines); a blocked `gbtrf`
  from LAPACK's `dgbtrf` (rstd runs the unblocked `xGBTF2`, since ReLAPACK's recursion is wrong);
  race faer's Bunch-Kaufman (different storage) against `sytrf`. Sweep: `relapack` on every rig, and
  with `-Drelapack=true` on one or two.
- [ ] **T119** Level-3 small-matrix follow-ups (T30): `sbgemm`'s small-matrix kernels are not
  ported; the per-core small kernels upstream has (AVX-512, SVE, POWER10, LASX) are not converted;
  the small-path thresholds are upstream's, tuned for its kernels, so time them under T05. Sweep:
  `3m`, `small`, `interface: gemm` on every rig, `aarch64-sve`, `powerpc64le-power10` and
  loongarch64 first (their permits allow the small path), `blat3_3m` on `x86_64` and `i386-max`.
- [ ] **T94** Fact-check items still open after T78: the sparc64 `-Oz` miscompile (needs a reduced
  case on the sparc64 rig); the OpenBIOS 0x400000 to 0x800000 window (our observation; check
  OpenBIOS's source); 32-bit x86 Windows threads starting at x87 0x27f (T03's observation; no
  source, no test records it: add a test); "MKL leaves FTZ to the application" (from memory); Tegra
  2 as VFPv3-D16 (a forum post only); `symmetric.zig`'s "4 elements of slack past each".
- [ ] **T101** rstd's `sdot` from 256K floats runs 51 to 53 GB/s on arm64 where its kernel alone
  runs 98 (R4, results-flyoff.md); bits unchanged. Level-1 threading or `Io.Threaded` (both landed
  since 6d66239) are the suspects; also look at threaded `cgemm` n = 128 in the VM, whose rate
  varies 198 to 249 GFLOP/s run to run.
- [ ] **T102** Strided `saxpy` and `zscal` (inc 2) in the arm64 Linux VM halved between 6d66239 and
  28171adc for both rstd and the converted kernel (`saxpy` 54 to 28 GB/s, `zscal` 139 to 71; R4,
  results-flyoff.md). Find whether it is code or the guest.
- [ ] **T21** Zen 3 (flyoff): short real `dot`s (0.30 of OpenBLAS at 1024) and complex `gemm` (0.79
  to 1.01) trail OpenBLAS.
- [ ] **T20** Kernel strategy: race `@Vector` against the converted assembly per core on a quiet
  machine (only vortex is done, `tools/blasbench/results-m6.md`); keep per-core code only where Zig
  loses. riscv64 races (from T22): `@Vector` against the `asm` `axpy` on zvl128b and x280, and
  fixed-length `@Vector` against VLEN-agnostic RVV on hardware wider than a core's minimum VLEN.
- [ ] **T05** Thresholds: OpenBLAS's defaults and per-core overrides, tuned by benchmark; then the
  overrides left out of phase 1 (z13 `axpy`, Windows-on-ARM64 `gemv`) if the numbers agree.
  **Needs:** T01, T73. T06's smoke run (busy Zen 3, `--threads=4`): threaded level-1 calls were
  slower than serial at the sizes it ran, so the thresholds need this tuning. Also time T07's split
  Accelerate calls against one whole call; and the bench's "Accelerate multithreaded" column is
  single-threaded while rstd's init sets `VECLIB_MAXIMUM_THREADS=1` process-wide, so a fair column
  runs Accelerate in a separate process. From T113's table: race AOCL-libFLAME's small-matrix
  AVX2/AVX-512 kernels (`getrf`, `potrf`, `potri`, `getrs`, `getri`, `geqrf`, `gelqf`, small
  `gesvd`; its OpenMP paths only after a same-bits check), libflame's restructured QR for `gesvd`
  and Hermitian eigenvalues (Van Zee et al., TOMS 2014), `FLA_Hevd`/`FLA_Tevd`, and
  `FLA_Eig_gest`/`FLA_Sylv` against ReLAPACK's `sygst`/`trsyl`.
- [ ] **T88** Time T25's races on an idle machine: portable against `@Vector` for `sum`, `sdot` in
  `f64` and in `f32` lanes, `nrm2` and the extremes, natively and as an x87 build (a scalar-FPU
  stand-in); bench and runner in the scratchpad (`vecsum/_vecsum_bench.zig`, `timeit.ps1`). Also
  race `vortex`'s converted `sum` against `@Vector` on the M6 (kept in `vector.displacing.keep`
  until then).
- [ ] **T26** vortex: a tuned `@Vector` `rot` (the assembly still leads by 2 to 17% in L1); race the
  Zig complex `gemm` against OpenBLAS's vortex assembly on the M6.
- [ ] **T27** gemm tuning: faster Zig packing (whole-vector transposes in `oncopy`, `incopy`),
  `C`-tile prefetch at kernel entry, a larger `KC` on Apple cores (comparisons.md (9)).
- [ ] **T28** Measure whether composed `symv`/`hemv` pay off on cores with assembly gemv (atom,
  opteron, arm64); today they compose only where the core has its own gemv.
- [ ] **T29** The bfloat16 intrinsics (cooperlake, sapphirerapids, Neoverse N2, V1, V2, POWER8 to
  10); `sbdot` could call the core's own `sdot`. **Needs:** T20.
- [ ] **T31** Assembly where `@Vector` cannot reach: arm64 SVE, SME (`vortexm4`), POWER10 MMA, AMX
  and AVX-512 BF16 for `sbgemm`; Skylake-X only if under 85% of its kernel (design.md,
  Implementation plan; `blas3/gemm/notes.md`). **Needs:** T20.
- [ ] **T32** BLIS assembly kernels, decided per kernel where benchmarks favor them (BSD-3, with a
  license section and per-file provenance). **Needs:** T20.
- [ ] **T39** `sparc64`: replace the image split with a real fix, mapping more memory at boot
  through the Open Firmware client interface (claim and map), so one image fits. Until then, ten
  sparc64-debug parts sit inside the window but above 0x3c0000 headroom: `openblas-l1` 0x3f3800,
  `blas2-edge-complex-cores-2` 0x3f2300, `lapack-getrf` 0x3f1cc0, `openblas-cores-l2-complex`
  0x3d99c0, `lapack-cholesky` 0x3d9000, `blas2-edge-complex` 0x3d4880, `blas3-kernel` 0x3c7000,
  `blas2-gemv-cores` 0x3c3000, `openblas-blat2` 0x3c04c0, `blas2-edge-real` 0x3bca80 (T95's agent;
  `scratchpad/sparcsplit/spans.py` reads them); split them before code growth pushes them over.
- [ ] **T71** Linux rigs with no item: loongarch64 (the `AT_HWCAP` LSX/LASX probe), big-endian
  powerpc64 (`AT_PLATFORM`, the `64be` variant), sparc64, 32-bit mips.
- [ ] **T41** `zrot_ppc440.S` (`ppc440fp2`) has a zero-stride patch no rig has run; no rig runs
  `ppc440` (it shares code with `ppcg4`). The 32-bit `ppc440` and `ppc440fp2` kernels rebuilt in
  719a41b are compile-checked only.
- [ ] **T74** No CI runs a BLAS test: rstd's `.gitlab-ci.yml` is rgui's copy (`EXE_NAME: "rgui"`,
  Zig 0.16.0). A job for `rstd-test`, `rstd-check-macho`, `rstd-check-uefi` and `rstd-check-zig017`
  (on Zig 0.17.0) would catch drift between rig sweeps; optional if testing stays manual.
- [ ] **T53** macOS VM extras: base-image overlays (`qemu-clone`, `-Dmac-base`), a `Forgejo`
  self-hosted runner, a setup that needs no clicking, less guest load when headless (HID driver,
  `SecurityAgent`), and freeze detection that asks before killing.
- [ ] **T104** Optional, best effort: the reference's GCC build contracts the MSA kernels' `y0 +=
  da_vec * x0` into `fmadd` on every MSA core, but rstd fuses vectors only on `loongson3r4` (T86).
  The same inline MSA `fmadd` from `fma.zig` could model `p5600`, `i6400`, `i6500` and `p6600` as
  the reference builds them; QEMU rounds neither way, so it rests on reading the assembly.
- [ ] **T76** Optional, from comparisons.md: `Dot2` (compensated `f32` lanes) if T25 finds the `f64`
  `sdot` costly in cache; BLIS's `JR`/`IR` ways for threads sharing an L2; an A-panel prefetch
  beside T27's `C`-tile one; small-`n` `sgemv_n` tuning (`results-m6.md`).
- [ ] **T54** Optional: an inline-asm `maddf` for R6 scalar multiply-adds in the main library, if R6
  speed ever matters (LLVM calls `fmaf`).
- [ ] **T52** On each Zig update, recheck and drop the workarounds once fixed (T48 and T49 track
  lld's SPARC and `zig fetch`): the powerpc64le two-float return (`zdot` patch), arm soft-float
  struct returns, `f80` under `x86_64_sysv` on Windows, the comptime naked-function call, ppc64
  big-endian ELFv2, the two-lane `f32` vector stored as 16 bytes, sparc64 `_Qp_sqrt`
  (`check.sqrtW`), sparc64 `-Oz`, `SafeAllocator`'s CPU count, the `Thread` override doc, run-step
  timeouts and output, the interrupted-link cache entry, inline-asm errors exiting 0
  (`checkedObject`), the MIPS64r6 select copy and spill (comptime `conj`), MIPS MSA vector `fma`
  (MSA cores differ from the reference's fused bits until fixed), arm CPU models without an FPU (rig
  pins), WHPX passing AMD CPUID (macOS rig on TCG), QEMU `xststdcsp` (portable `nrm2` avoids
  `isInf`). Also, wasm bitcode getting its feature flags twice (`asmCheck`'s
  `-disable-llvm-verifier`). The Zig entries were last checked on dev.751. Also: building the
  library for an R6+MSA target CPU (`i6400`) trips LLVM's `Cannot copy registers` assertion
  (`MipsSEInstrInfo.cpp:261`) in `blas2.interface`, `blas3` and `_test`, likely the R6 select bug
  already recorded (fixed upstream after 23.1.x); the rigs' main CPU has no MSA, so it never shows
  there (T89).
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
- [ ] **T105** The power10 objects in a 32-bit POWER build don't link and are untested: the dev
  toolchain's LLVM can't build 32-bit compiler-rt for POWER4 to POWER10 (`64bit` crashes, VSX with
  `fpcvt` fails an `f128` conversion, POWER9 vector and ISA 3.0 features error; Zig 0.17.0 builds
  them all). `powerpc-power9` (T11) runs power8 and power9 only. Any real 32-bit POWER8 to 10 build
  of rstd fails with this toolchain, against design.md's filling-upstream-gaps decision.
  **Blocked:** upstream (LLVM; issues.md, 32-bit compiler-rt).

## Done

- [X] **T70** OpenBLAS's own LAPACK past `potrf`, `potrs`, `laswp`, `getrf` and `getrs` (`getf2` is
  `getrf`'s base case): `gesv`, `potf2`, `potri`, `lauu2`/`lauum`, `trti2`/`trtri`, `trtrs`, and
  `laed3` with `test_laed3_nan.c` (`interface/lapack/`, `lapack/`). Threaded variants
  (`getrf_parallel.c`, the only caller of `laswp_ncopy`; `potrf_parallel.c`, T35) once T01 lands.
  **Needs:** T01 for the threaded ones. Running (`wip/lapack2`). Done 2026-10-30, landed ee46ad6d:
  `gesv`, `potf2`, `potri`, `trtri`/`trti2`, `lauum`/`lauu2`, `trtrs` (`lapack/triangular.zig`),
  `getrs` on workers and `getrf` with `getrf_parallel.c`'s lookahead over `getrf_single.c`'s blocks
  (its thread-dependent block widths would break the bit rule; bit-identical to serial at 2, 3, 8
  workers). Reference wins on `potri`'s zero diagonal and `trtrs`'s argument order (issues.md).
  `laed3` and `test_laed3_nan.c` moved to T111 (faer).
- [X] **T35** Optionally route `potrf` to Accelerate on macOS; threaded `potrf`
  (`potrf_parallel.c`). **Needs:** T01. Running (`wip/lapack2`). Done 2026-10-30, landed e7dd85f0:
  threaded `potrf` (`potrf_[UL]_parallel.c`), the algorithm picked by size (upstream's interface
  thresholds), never worker count; `potrf` routed to Accelerate on macOS (T114 verifies it).
- [X] **T65** The 32-bit x86 cores from SSE2 up (`northwood` to `bobcat`, readme.md, and 886324d's
  fallback objects) never run: the `i386` rig is `pentium3`. Add a freestanding `i386-max` rig
  (`-cpu max`; the TCG `fldt` bug is a paging case, issues.md), or a native `-Dtarget=x86-windows`
  run under WoW64 with `-Dblas-cores=all`, which also tests T03's x87 case. Part of T62's sweep.
  Debug passes on the Mac (477 passed, 0 failed, all 20 cores), landed 6676d41d with `x86.ld`
  keeping the Multiboot header; ReleaseSafe run left. Running (`wip/i386max`). Done 2026-10-30,
  landed 6676d41d and the runner fix: Debug 477 passed, 0 failed; ReleaseSafe 486 passed, 31
  skipped, 0 failed; all 20 32-bit x86 cores in each of 71 per-core tests. The ReleaseSafe hang was
  the rig's runner: `paintStack`, inlined into `qemutest_main`'s 16 KiB release frame, painted over
  main's own loop bound; it is now `noinline` (`test/qemu/runner.zig`).
- [X] **T100** At the smallest `n`, rstd's public level-1 routines (`axpy`, `scal`, `copy`, `swap`,
  `rot`, `iamax`, `dot`, `asum`) now run 7 to 30% behind the kernel they call, on macOS and Linux;
  they matched it at 6d66239 (R4, results-flyoff.md). Find the per-call overhead (likely the level-1
  threading dispatch) and remove it for small `n`. Running (`wip/smallfix`). Done 2026-10-30: the
  per-call threading plan (core-name lookup and rule) cost `dot` 11 to 28 ns; calls at or below
  every core's `serial_max` (computed at compile time) now skip it, with identical results; `ddot` n
  = 4 went from 3.66 to 0.99 of its kernel (native Zen 3, ReleaseFast). The R4-style timing confirms
  it later.
- [X] **T18** The ISA check cannot see inline assembly (x86's assembler ignores CPU features).
  Detection, as OpenBLAS's, hands some cores to weaker CPUs: `katmai` to any 32-bit x86 without SSE3
  (down to the 486), `vortex` to any Apple CPU (pre-A12 lacks FCMLA), `vortexm4` to non-Apple SME
  without SVE, a forced core's fallback to `neoversen1` on ARMv8.0; decide per case. Running
  (`wip/smallfix`). Done 2026-10-30 (design.md, Cores handed to weaker CPUs): `katmai` on pre-SSE3
  32-bit x86 and `vortex` on Apple CPUs kept as OpenBLAS (no target affected); `vortexm4` now needs
  an Apple CPU and non-Apple SME without SVE takes the SVE test; a forced arm64 core falls back only
  within what the detected core has (hit Raspberry Pi 3 and 4).
- [X] **T42** Test-only code still auto-vectorizes on the G4 (`blat1` rot, the scal, gemv and
  half-gemv models, `triangular.checkAll`, `blat3`'s complex type). Guard with `util.scalar` the
  loops a test feeds subnormals (AltiVec flushes them); close the rest, which see normal-range data
  only. Running (`wip/smallfix`). Done 2026-10-30, no change: a ReleaseSafe G4 assembly scan finds
  test loops in AltiVec lanes on normal-range data only; subnormal-fed loops are scalar (issues.md).
- [X] **T38** Check in the C-against-Zig differential harness (stub `common.h`) used for `dot`.
  Running (`wip/smallfix`). Done 2026-10-30: `tools/dotverify` (stub `common.h`, `run.sh`) and
  `_verify_dot.zig`; 8352 cases, 0 mismatches across nehalem, sandybridge, haswell, zen, skylakex
  (AVX2 width) and generic; a swapped-reference control finds 4328.
- [X] **T87** Debug builds round-trip every value through memory, which hides x87 width problems;
  T56's reduction failures showed only in ReleaseSafe. Run the x87 rigs (`i386`, and T65's) and
  native x87 runs in ReleaseSafe as well as the default mode, at least in each full sweep. Running
  (`wip/smallfix`). Done 2026-10-30: rigs `i386-safe` and `i386-max-safe` (ReleaseSafe copies,
  `x87Safe` in `build/qemutest.zig`) join the sweep; `rstd-test-x87` runs the native pentium3 suite
  in Debug and ReleaseSafe.
- [X] **T07** Tie Accelerate into rstd's threading (threading.md, Vendor backends). Found by T06:
  rstd's init makes Accelerate single-threaded for the process (`VECLIB_MAXIMUM_THREADS=1`) and the
  calling thread, so the bench's macOS "Accelerate multithreaded" rows (`--threads=default`) likely
  ran one thread, and `results-m6-libraries.md`'s "gains at most 1.2 times from more threads" may
  rest on that. The bench now prints Accelerate's state after rstd's init; a fair multithreaded
  Accelerate column needs it on its own thread, or this item. Running (`wip/macaccel`). Done
  2026-10-30, landed: Accelerate's own threads stay off (its threaded `gemv` rounds differently from
  one thread, and Apple documents no fix; threading.md, Vendor backends). Above rstd's threshold, a
  routed call splits by size only (even at one worker) into single-threaded Accelerate calls on
  rstd's workers, reductions combined in part order: level 1
  `copy`/`swap`/`scal`/`rot`/`axpy`/`dot`/`dotc`/`dsdot`/`asum`, `gemv`, `ger`, `gemm`, `symm`,
  `hemm`, `trmm`, `trsm`; the rest stay one call. Split calls differ in the last bits from one whole
  call. Mac 41/41.
- [X] **T36** The Accelerate level-2 special-value sweep takes minutes per Mac run (complex, even
  trimmed to every third variant at the largest shape); speed it up. Running (`wip/macaccel`). Done
  2026-10-30, landed be39c986: the Accelerate level-2 sweep runs its routines on a thread per CPU,
  the same cases; 35.5 s to 6.0 s (the cost is the `f128` reference).
- [X] **T109** In `blas2/banded.zig` (about line 1004), "an allocation failure comes back as
  OutOfMemory" saves `const was = blas.useAccelerate(false)`, but `useAccelerate` returns the state
  after the change, so `was` is always false and routing never comes back on: on macOS every later
  test in that process runs without Accelerate. Save `apple_accel.routing()` first, as
  `blas3/threaded.zig`'s `rstdOnly` does (T108), and check the later tests on the Mac. Running
  (`wip/macaccel`). Done 2026-10-30, landed: the test saves `accel.routing()` before turning routing
  off.
- [X] **T112** Since 2cfd70ff/13768d8a (2026-10-29) the portable blocking's buffer is 128 MB
  (`BUFFER_SIZE`), so its `R` is about 262K `f32` columns and anything scaling with `R` grew
  fourfold; T103 found it overflowing the freestanding heap. Check the generic table's 128 MB
  against OpenBLAS's `common_*.h` for the targets that use it, and keep it only if the reference
  does. Running (`wip/sources`). Done 2026-10-30, landed 7c6e8474 (via `wip/sources`): the generator
  had read the generic table's 128 MB from `common_x86_64.h:260`; OpenBLAS's generic target uses its
  own architecture's header (x86 16 MB, arm and arm64 32 MB, POWER 256 MB under `DYNAMIC_ARCH`,
  Alpha and C-SKY 32 MB). A `.generic` install now takes the build target's own size
  (`native_buffer_size`); 128 MB stays only on x86_64 and where OpenBLAS has no header. Sweep:
  `params` and `blas3` on `i386`, `i386-max`, `arm-softfp`, `powerpc-750`. Found T116.
- [X] **T113** The anthology's source table (design.md, Settled decisions, 2026-10-30): for each
  LAPACK routine, and each BLAS kernel family where more than one source exists, list the candidates
  (OpenBLAS's own, ReLAPACK, faer, libflame, AOCL-libFLAME, BLIS, reference LAPACK as baseline), the
  one chosen, why, and its benchmark or accuracy evidence or "not yet raced". Clone libflame and
  AOCL-libFLAME read-only beside OpenBLAS and BLIS, pinned, so agents can compare them. Feeds the
  timing session. Eigen is listed only as a benchmark baseline, never as a source. Running
  (`wip/sources`). Done 2026-10-30, landed 7f9e0fbe: `src/linalg/lapack/sources.md`, about 50 LAPACK
  rows and 18 BLAS kernel families, each candidate checked in source; libflame 5.2.0 (`53bd60c50b`)
  and AOCL-libFLAME 5.3.2 (`3b3251adec`) cloned beside OpenBLAS, both BSD-3 (AOCL also MIT per
  file). libflame's `gesdd` is its `gesvd` and its `gelsd` its `gelss`. ReLAPACK also covers `pbtrf`
  and `tgsyl` (T110). Found T115; timing candidates added to T05.
- [X] **T75** Optional: `?gemm_batch` and `?gemm_batch_strided` (`interface/gemm_batch.c`,
  `gemm_batch_strided.c`) and the `sbgemm`/`bgemm` batch forms, column-major only; a strided batch
  maps onto M3's fixed-size forms. Running (`wip/l3half`). Done 2026-10-30, landed 8dd19974:
  `gemmBatch`, `gemmBatchStrided`, `gemmHalfBatch`, `gemmHalfBatchStrided` (`blas3/batch.zig`),
  column-major; every group checked before any problem runs, as upstream.
- [X] **T34** LAPACK SVD from faer (Q10): `gesvd` and `gesdd` (`test_kernel_regress.c` needs
  `dgesvd`), ported from faer's SVD (MIT, pure Rust) behind LAPACK-style entry points in
  `rstd.linalg.lapack`, its matrix products through rstd's BLAS; results bit-identical for any
  thread count. Running (`wip/faer`). After it lands, add `lapack/sources.md` (T113) to design.md's
  Source tree row for `../lapack/`, which this item rewrites. Done 2026-10-30, landed (through
  73bfefb2): faer `faer-v0.24.4` (`0539947`, Codeberg, MIT) in `src/linalg/lapack/faer/`, one file
  per faer source; `gesvd` (bidiagonal QR iteration) and `gesdd` (divide and conquer, QR below 128)
  in all four precisions, LAPACK's arguments and scaling; products through rstd's level-3 drivers;
  bit-identical at 1, 3, 8 workers. OpenBLAS's `test_kernel_regress.c` now runs and passes. Four
  faer bugs fixed and recorded in issues.md (a root finder stopping early near a pole, reciprocal of
  a denormal norm, a zero subproblem's missing 1, an unwritten block of full `U`). Follow-ups: T117.
- [X] **T110** ReLAPACK (MIT; OpenBLAS's `relapack/`, its `BUILD_RELAPACK` option) for the
  factorizations it covers that OpenBLAS doesn't reimplement (`sytrf`/`hetrf`, `trsyl`,
  `sygst`/`hegst`, `gbtrf`, `gemmt`'s LAPACK uses); for `potrf`, `getrf`, `trtri` and `lauum`,
  OpenBLAS's own stay the default and ReLAPACK's are a build option, as in OpenBLAS (Q10).
  **Needs:** T70 to land (same area). Running (`wip/relapack`). Done 2026-10-30 except `tgsyl`
  (T118), landed through 7b81809e: ReLAPACK's `sygst`/`hegst`, `pbtrf`, `trsyl` and internal `gemmt`
  by default; its `potrf`, `getrf`, `trtri`, `lauum` behind `-Drelapack` (OpenBLAS's stay default;
  match within 64n ulps). ReLAPACK's `sytrf`/`hetrf` (and rook) and `gbtrf` give wrong factors
  upstream (checked against its own C), so rstd runs LAPACK's blocked loop over ReLAPACK's panel
  kernel, and the unblocked `xGBTF2`; six more upstream bugs fixed (issues.md). Base cases (`sytf2`,
  `sygs2`, `pbtf2`, `gbtf2`, `laln2`, `lasy2`) are reference LAPACK ports, faer having none in
  LAPACK's storage (Q11).
- [X] **T23** riscv `zvl256b` has no object (`zvl128b` serves it): read VLEN and port the
  `*_vector.c` kernels (design.md, Open issues). **Needs:** T20. Running (`wip/portrv`). Done
  2026-10-30, landed f7b6be16: `riscv64_zvl256b` has its own object (VLEN from `vlenb`), its
  `@Vector` kernels and gemm tiles at 256 bits.
- [X] **Q11** May rstd port reference LAPACK code where no other open-source candidate has the
  routine? T110 (landed on review) ported ReLAPACK's base cases (`sytf2`, `sygs2`, `pbtf2`, `gbtf2`,
  `laln2`, `lasy2`) from reference LAPACK, since faer has none in LAPACK's storage, and added
  reference LAPACK's BSD notice to `LICENSE`; T115 queues `gecon`, `gebal`/`gebak`, `trexc`/`trsen`,
  `gelsy` the same way, and T118's `tgsyl` would too. Its LICENSE is BSD-3's three conditions plus a
  patent disclaimer, inside the anthology decision's "MIT or BSD-style" rule, but design.md (LAPACK
  section) still says its license is not yet cleared. Recommended: yes, as the anthology's source of
  last resort, and I update that line; if no, T110's base cases come out and those routines wait.
  Answered 2026-10-30: yes, reference LAPACK is always available as a fallback (permissive, no
  copyleft, free for commercial reuse, as MIT); design.md updated. ReLAPACK's base cases are not in
  ReLAPACK: it calls LAPACK's unblocked routines rather than shipping them, hence the reference
  ports. T115 can start.
- [X] **T30** Level 3: `gemm3m`, `symm3m`, `hemm3m` (stubs kept), a small-matrix path,
  half-precision level 3 (`shgemm`, whose vtable slots are all null), bfloat16 `sbgemm` and
  `sbgemmt` (OpenBLAS's bfloat16 `gemmt`), then `bgemm` (design.md, Implementation plan). Running
  (`wip/l3ext`, `wip/l3half`). Half-precision part landed 2026-10-30 (`sbgemm`, `bgemm`, `shgemm` as
  `gemmHalf`, `sbgemmt` as `gemmtHalf`, portable kernels, `blas3/half.zig`); 3m and the small-matrix
  path left (`wip/l3ext`). Done 2026-10-30: 3M (`blas3/gemm3m/`, `gemm3m_level3.c`, `symm3m_k.c`,
  `hemm3m_k.c`; on every target, never routed to Accelerate) and the small-matrix path
  (`blas3/gemm/small_kernel.zig`, upstream's per-core permits; POWER10's single-precision permit
  fixed at 64 cubed, upstream's follows the thread count) landed through 573d70a9. Follow-ups: T119.
- [X] **T69** OpenBLAS's comparison tests are not ported: `compare_sgemv_{bgemv,sbgemv,shgemv}.c`
  cover routines rstd has (`gemvHalf`); `compare_sgemm_{bgemm,sbgemm,shgemm}.c` come with T30.
  `blat3.zig` runs `[cz]blat3_3m` through `gemm`: switch it to the real `gemm3m`, with 3m's looser
  imaginary-part bound (comparisons.md (7)). **Needs:** T30 for the level-3 parts. Running
  (`wip/l3ext`, `wip/l3half`). The `compare_*` ports and `test_bgemm.c` landed 2026-10-30; the
  `blat3` 3m switch left (`wip/l3ext`). Done 2026-10-30: `[cz]blat3_3m` runs the real `gemm3m`, the
  3M-built part gauged at three times the bound (2e5ca907).
- [X] **T72** UEFI never runs the tests: `rstd-check-uefi` only links. Add a run rig for the OS-free
  suite under OVMF and AAVMF (UEFI's inline `Io`; RandyOS boots through it). Running
  (`wip/uefirig`). Done 2026-10-30, landed 318951f1: `x86_64-uefi` (OVMF) and `aarch64-uefi` (AAVMF)
  run the OS-free tests as a UEFI app on rstd's UEFI `Io`, sharing the freestanding runner
  (`test/qemu/uefi.zig`); firmware read-only pflash from Nix's QEMU data, so they run through
  `-Dtest-host` only; smoke runs 105 and 103 passed, 0 failed.

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
  `trevc` (`comqr3`), `hetrd` and `unmtr` (`htridi`, `htribk`); rstd has `potrf`, `getrf` and
  `getrs`; T70 adds OpenBLAS's own, and the rest come from Reference-LAPACK (BSD-3, as OpenBLAS
  ships it in `lapack-netlib/`). Its fills (`wset`, `rset`) are LAPACK `laset`, `pythag` is `lapy2`.
  Its own tests must still pass. **Needs:** M1.

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
