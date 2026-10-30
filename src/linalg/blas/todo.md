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

- [ ] **T65** The 32-bit x86 cores from SSE2 up (`northwood` to `bobcat`, readme.md, and 886324d's
  fallback objects) never run: the `i386` rig is `pentium3`. Add a freestanding `i386-max` rig
  (`-cpu max`; the TCG `fldt` bug is a paging case, issues.md), or a native `-Dtarget=x86-windows`
  run under WoW64 with `-Dblas-cores=all`, which also tests T03's x87 case. Part of T62's sweep.
  Debug passes on the Mac (477 passed, 0 failed, all 20 cores), landed 6676d41d with `x86.ld`
  keeping the Multiboot header; ReleaseSafe run left. Running (`wip/i386max`).
- [ ] **T30** Level 3: `gemm3m`, `symm3m`, `hemm3m` (stubs kept), a small-matrix path,
  half-precision level 3 (`shgemm`, whose vtable slots are all null), bfloat16 `sbgemm` and
  `sbgemmt` (OpenBLAS's bfloat16 `gemmt`), then `bgemm` (design.md, Implementation plan). Running
  (`wip/l3ext`, `wip/l3half`).
- [ ] **T69** OpenBLAS's comparison tests are not ported: `compare_sgemv_{bgemv,sbgemv,shgemv}.c`
  cover routines rstd has (`gemvHalf`); `compare_sgemm_{bgemm,sbgemm,shgemm}.c` come with T30.
  `blat3.zig` runs `[cz]blat3_3m` through `gemm`: switch it to the real `gemm3m`, with 3m's looser
  imaginary-part bound (comparisons.md (7)). **Needs:** T30 for the level-3 parts. Running
  (`wip/l3ext`, `wip/l3half`).
- [ ] **T75** Optional: `?gemm_batch` and `?gemm_batch_strided` (`interface/gemm_batch.c`,
  `gemm_batch_strided.c`) and the `sbgemm`/`bgemm` batch forms, column-major only; a strided batch
  maps onto M3's fixed-size forms. Running (`wip/l3half`).
- [ ] **T70** OpenBLAS's own LAPACK past `potrf`, `potrs`, `laswp`, `getrf` and `getrs` (`getf2` is
  `getrf`'s base case): `gesv`, `potf2`, `potri`, `lauu2`/`lauum`, `trti2`/`trtri`, `trtrs`, and
  `laed3` with `test_laed3_nan.c` (`interface/lapack/`, `lapack/`). Threaded variants
  (`getrf_parallel.c`, the only caller of `laswp_ncopy`; `potrf_parallel.c`, T35) once T01 lands.
  **Needs:** T01 for the threaded ones. Running (`wip/lapack2`).
- [ ] **T35** Optionally route `potrf` to Accelerate on macOS; threaded `potrf`
  (`potrf_parallel.c`). **Needs:** T01. Running (`wip/lapack2`).
- [ ] **T07** Tie Accelerate into rstd's threading (threading.md, Vendor backends). Found by T06:
  rstd's init makes Accelerate single-threaded for the process (`VECLIB_MAXIMUM_THREADS=1`) and the
  calling thread, so the bench's macOS "Accelerate multithreaded" rows (`--threads=default`) likely
  ran one thread, and `results-m6-libraries.md`'s "gains at most 1.2 times from more threads" may
  rest on that. The bench now prints Accelerate's state after rstd's init; a fair multithreaded
  Accelerate column needs it on its own thread, or this item. Running (`wip/macaccel`).
- [ ] **T36** The Accelerate level-2 special-value sweep takes minutes per Mac run (complex, even
  trimmed to every third variant at the largest shape); speed it up. Running (`wip/macaccel`).
- [ ] **T109** In `blas2/banded.zig` (about line 1004), "an allocation failure comes back as
  OutOfMemory" saves `const was = blas.useAccelerate(false)`, but `useAccelerate` returns the state
  after the change, so `was` is always false and routing never comes back on: on macOS every later
  test in that process runs without Accelerate. Save `apple_accel.routing()` first, as
  `blas3/threaded.zig`'s `rstdOnly` does (T108), and check the later tests on the Mac. Running
  (`wip/macaccel`).
- [ ] **T22** C and intrinsics kernels not yet ported, levels 1 and 2, by T20's strategy, `@Vector`
  first as arm64 did (`vector.zig`); the coverage.md gaps, the largest first. Level 1: riscv64
  outside `axpy`; POWER8 to 10 `asum`, `scal`, `rot`, `swap`, `cdot`; MIPS MSA and zarch; x86_64
  `asum` on Haswell and AVX-512; wasm SIMD128. Level 2 (464 gaps, `blas2/coverage.md`): riscv64 140
  (gemv 82, hemv 32, symv 24, `sbgemv` and `shgemv` on zvl256b); mips64 100 (gemv on I6400, I6500,
  P6600, Loongson 3R3, 3R4); arm64 60 (symv on 11 cores, `sger`/`dger` on A64FX and Neoverse V1,
  `sbgemv`/`bgemv` on N2, V1, V2); POWER 60 (gemv on POWER8 to 10, `sbgemv` on POWER10); x86_64 50
  (symv on 10 cores, `sgemv` on the Skylake-X family (T24), `sger`/`dger` on Sandy Bridge,
  `sbgemv_n` on Cooper Lake and Sapphire Rapids); zarch 30 (gemv on z13, z14); mips 20 (gemv on
  P5600); wasm 4 (gemv). **Needs:** T20. Running (`wip/portrv`: riscv64 part only`).
- [ ] **T23** riscv `zvl256b` has no object (`zvl128b` serves it): read VLEN and port the
  `*_vector.c` kernels (design.md, Open issues). **Needs:** T20. Running (`wip/portrv`).
- [ ] **T100** At the smallest `n`, rstd's public level-1 routines (`axpy`, `scal`, `copy`, `swap`,
  `rot`, `iamax`, `dot`, `asum`) now run 7 to 30% behind the kernel they call, on macOS and Linux;
  they matched it at 6d66239 (R4, results-flyoff.md). Find the per-call overhead (likely the level-1
  threading dispatch) and remove it for small `n`. Running (`wip/smallfix`).
- [ ] **T18** The ISA check cannot see inline assembly (x86's assembler ignores CPU features).
  Detection, as OpenBLAS's, hands some cores to weaker CPUs: `katmai` to any 32-bit x86 without SSE3
  (down to the 486), `vortex` to any Apple CPU (pre-A12 lacks FCMLA), `vortexm4` to non-Apple SME
  without SVE, a forced core's fallback to `neoversen1` on ARMv8.0; decide per case. Running
  (`wip/smallfix`).
- [ ] **T42** Test-only code still auto-vectorizes on the G4 (`blat1` rot, the scal, gemv and
  half-gemv models, `triangular.checkAll`, `blat3`'s complex type). Guard with `util.scalar` the
  loops a test feeds subnormals (AltiVec flushes them); close the rest, which see normal-range data
  only. Running (`wip/smallfix`).
- [ ] **T38** Check in the C-against-Zig differential harness (stub `common.h`) used for `dot`.
  Running (`wip/smallfix`).
- [ ] **T87** Debug builds round-trip every value through memory, which hides x87 width problems;
  T56's reduction failures showed only in ReleaseSafe. Run the x87 rigs (`i386`, and T65's) and
  native x87 runs in ReleaseSafe as well as the default mode, at least in each full sweep. Running
  (`wip/smallfix`).
- [ ] **T34** LAPACK SVD from faer (Q10): `gesvd` and `gesdd` (`test_kernel_regress.c` needs
  `dgesvd`), ported from faer's SVD (MIT, pure Rust) behind LAPACK-style entry points in
  `rstd.linalg.lapack`, its matrix products through rstd's BLAS; results bit-identical for any
  thread count. Running (`wip/faer`).

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
- [ ] **T112** Since 2cfd70ff/13768d8a (2026-10-29) the portable blocking's buffer is 128 MB
  (`BUFFER_SIZE`), so its `R` is about 262K `f32` columns and anything scaling with `R` grew
  fourfold; T103 found it overflowing the freestanding heap. Check the generic table's 128 MB
  against OpenBLAS's `common_*.h` for the targets that use it, and keep it only if the reference
  does.
- [ ] **T111** The rest of LAPACK that would otherwise come from reference LAPACK, ported from faer
  (Q10): QR (`geqrf`, `orgqr`/`ungqr`, `ormqr`), least squares (`gels`), the symmetric and Hermitian
  eigensolvers (`syevd`/`heevd`), the general eigensolver (`geev`), and whatever else rstd.linalg's
  MATLAB-like API needs. Reference LAPACK comes later only as a benchmark baseline. **Needs:** T34's
  port layout.
- [ ] **T110** ReLAPACK (MIT; OpenBLAS's `relapack/`, its `BUILD_RELAPACK` option) for the
  factorizations it covers that OpenBLAS doesn't reimplement (`sytrf`/`hetrf`, `trsyl`,
  `sygst`/`hegst`, `gbtrf`, `gemmt`'s LAPACK uses); for `potrf`, `getrf`, `trtri` and `lauum`,
  OpenBLAS's own stay the default and ReLAPACK's are a build option, as in OpenBLAS (Q10).
  **Needs:** T70 to land (same area).
- [ ] **T113** The anthology's source table (design.md, Settled decisions, 2026-10-30): for each
  LAPACK routine, and each BLAS kernel family where more than one source exists, list the candidates
  (OpenBLAS's own, ReLAPACK, faer, libflame, AOCL-libFLAME, BLIS, reference LAPACK as baseline), the
  one chosen, why, and its benchmark or accuracy evidence or "not yet raced". Clone libflame and
  AOCL-libFLAME read-only beside OpenBLAS and BLIS, pinned, so agents can compare them. Feeds the
  timing session. Eigen is listed only as a benchmark baseline, never as a source.
- [ ] **T72** UEFI never runs the tests: `rstd-check-uefi` only links. Add a run rig for the OS-free
  suite under OVMF and AAVMF (UEFI's inline `Io`; RandyOS boots through it).
- [ ] **T99** Multithreaded freestanding tests (from T02, see threading.md): start secondary cores
  in the test kernels with their own stacks, give each core thread-local storage (std's `debug.zig`
  makes `panic_stage` thread-local once a build is not `single_threaded`), an `Io` that hands group
  tasks to those cores with futex waits, wakes and cancellation, and start the rigs with more than
  one CPU.
- [ ] **T19** `XDouble` as `f128`: `quad_precision` is hard-wired off; add a build option and a rig
  configuration for it. The edge tests and level 3's per-routine driver tests still cover no
  extended type off x86 (the interfaces do).
- [ ] **T73** Benchmarks above level 1 race `gemv` and `gemm` only: add `ger`, `symv`, `trsv`, a
  banded routine, `symm`, `syrk`, `trmm` and `trsm`, which T05 and the speed claims need. Rerun the
  Zen 3 flyoff (Windows and Linux), which predates T43's x86_64 `rot`, `scal`, `gemv` paths and the
  level-3 threading, under a timing rule for that host. **Needs:** T06 for threaded rows. Reference
  binaries beside OpenBLAS (design.md, the anthology decision): reference LAPACK, and Eigen (MPL
  2.0) built from upstream through Nix as a benchmark-only binary, never ported, for the LAPACK and
  level-3 races.
- [ ] **T40** Linux rig gaps: `linux-powerpc64le` never runs POWER10 kernels (baseline build),
  `linux-aarch64` leaves out the SVE cores, and under HVF it also skips `vortexm4` (cause unknown).
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
  loses.
- [ ] **T05** Thresholds: OpenBLAS's defaults and per-core overrides, tuned by benchmark; then the
  overrides left out of phase 1 (z13 `axpy`, Windows-on-ARM64 `gemv`) if the numbers agree.
  **Needs:** T01, T73. T06's smoke run (busy Zen 3, `--threads=4`): threaded level-1 calls were
  slower than serial at the sizes it ran, so the thresholds need this tuning.
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
  10, riscv zvl256b); `sbdot` could call the core's own `sdot`. **Needs:** T20.
- [ ] **T31** Assembly where `@Vector` cannot reach: arm64 SVE, SME (`vortexm4`), POWER10 MMA, AMX
  and AVX-512 BF16 for `sbgemm`; Skylake-X only if under 85% of its kernel (design.md,
  Implementation plan; `blas3/gemm/notes.md`). **Needs:** T20.
- [ ] **T32** BLIS assembly kernels, decided per kernel where benchmarks favor them (BSD-3, with a
  license section and per-file provenance). **Needs:** T20.
- [ ] **T24** Skylake-X family `sgemv_n`/`sgemv_t`: those cores add large intrinsics kernels on top
  of Haswell's (`sgemv_t` uses a 1,121-line template); they keep the portable `sgemv` until ported.
  Testable only on AVX-512 hardware (T43).
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

- [X] **Q7** R4 still cannot start, and the load average looks like the wrong test. Since about
  11:50 the Mac has been idle by every other measure (no QEMU, no process over 5%, top processes
  Activity Monitor 4%, `airportd` 2.4%, `sysmond` 1.6%), yet the one-minute load sits at 1.4 to 1.7
  against a limit of about 1.05. From 11:15 to 11:25, with one QEMU core at 99% and Activity Monitor
  open, the load was 1.01 to 1.07: busier, but lower. So load here seems to count short wake-ups
  (Wi-Fi, Activity Monitor's sampling), not CPU use. Proposal: replace "load under 1" with "all
  processes but Activity Monitor and `sysmond` together under 10% of one core over the 10 s sample",
  keeping no QEMU, no process over 5%, two runs agreeing within 3%, and the load before and after
  recorded. Or keep load and raise it to under 2. R4 polls until you choose; T62 waits on R4.
  Details: scratchpad `retime2/poll.log`. R4's agent, after 8 h: on 65 polls only the load blocked
  it, median 1.53 at 99.2 to 99.7% idle with no runnable threads, so the load is a counting
  artifact; it suggests a load limit near 2.0 as the smallest change. Answered 2026-10-29: keep the
  load test with a limit of 2 (Activity Monitor and `sysmond` still excluded, their share added); R4
  restarted, the Mac held for it.
- [X] **Q6** How should T84's level 2 products thread? (a) Per-chunk partial vectors with chunks
  fixed by size, summed in chunk order: keeps the bit rule for any core, as threading.md's plan had
  it, but changes these routines' serial results in every build (large `n` now sums in chunks) and
  costs `n` elements per chunk. (b) Split the output, as for `gemv`: keeps today's serial bits, but
  only if every core's `axpy` gives an element the same bits wherever it sits in a call, which is
  untested across cores; it would need that test on every rig first. Recommended: (a), since the bit
  rule then holds by construction; OpenBLAS's threaded `symv` also gives each thread its own buffer
  and sums them into `y` with `AXPYU_K` (`driver/level2/symv_thread.c:197`, `:259-291`), though in
  an order that follows the thread count. Answered 2026-10-29: (a), fixed chunks; T84 running.
- [X] **Q8** Loongson 3A model (T86)? (A) Match the reference: model `loongson3r3` and `loongson3r4`
  as fused and emit scalar multiply-adds as inline `madd.fmt` (`@mulAdd` is a slow library call
  there); no rig can check the fused bits, since QEMU rounds `madd`. (B) Stay unfused and
  deterministic: add LLVM's `nomadd4` to `loongson3r4` (no `madd` asm in it) and document that
  Loongson users build with it; `loongson3r3` cannot take it (its converted kernels use `madd`), so
  its assembly keeps fusing. Recommended (the agent's): B for `loongson3r4`, accepting
  `loongson3r3`'s fused assembly, since without hardware neither choice can be tested. Answered
  2026-10-29: (A), model `loongson3r3` and `loongson3r4` as fused, as the reference; T86 running.
- [X] **R4** Mac timing rerun under the strict quiet rule, started automatically when the Mac is
  idle; other Mac runs wait while `scratchpad/MAC_TIMING_LOCK` exists (`wip/retime2`). Stopped 16:07
  without timing: in 8 h (97 polls) the Mac never met the rule; on 65 polls only the load blocked
  it, reading 1.1 to 2.3 while 99.2 to 99.7% idle with no runnable threads
  (`retime2/phantom-load.txt`). Ready to run (`wip/retime2` on 28171adc, benches prebuilt) once Q7
  sets the rule. The Mac now runs macOS 27.2 (26B123). Rerunning under the Q7 rule (load under 2),
  the Mac drained for it (`wip/retime2`). Done 2026-10-29 21:52 to 23:56 CDT, landed
  45199c68..0fe0c74f: one-thread numbers match the busy-Mac run (OpenBLAS median ratio 1.000), so
  the load did not skew them; only threaded OpenBLAS `gemm` on macOS gained 7 to 9% at n = 128 and
  256. One noisy row left, Linux threaded `cgemm` n = 128 (rstd 198 to 249 GFLOP/s; T100).
  Regressions found: T100, T101, T102.
- [X] **T84** Thread `symv`, `hemv`, `sbmv`, `hbmv`, `spmv`, `hpmv`, `trmv`, `tbmv`, `tpmv` and
  untransposed `gbmv`, left serial by T01: each output element sums over several kernel calls, so
  splitting the output moves the bits with the split. Also, `sbdot` and the bfloat16 conversions,
  not threaded (T05). Fixed chunks (Q6). Running (`wip/thread2b`). Done 2026-10-30, landed
  380d98ff..eaf33d7d: `symv`, `hemv`, `sbmv`, `hbmv`, `spmv`, `hpmv`, `trmv`, `tbmv`, `tpmv` and
  untransposed `gbmv` fold per-chunk partial vectors in chunk order (`threading.Partials`);
  thresholds from OpenBLAS's interface files, and `sbmv`/`hbmv`/`spmv` (no upstream threshold) take
  their siblings' rules (threading.md, Open issues). Native and x87 111 passed; `linux-aarch64`
  threaded level 2 26/26.
- [X] **T61** Uring (`Io.Evented`) hangs a plain futex-only barrier of four fibers once it has
  worker threads (`thread_limit = 3`; issues.md); the fiber test runs single-threaded Uring only.
  T02's `Io.Evented` tests start here. Running (`wip/uring`). Done 2026-10-30, landed
  efe4a960..3cc6e58c: the hang is std's `Uring.futexWake` passing `maxInt(u32)` to the kernel as -1,
  waking one waiter; `Signal.set` now wakes `maxInt(i32)`, with a test that fails without it.
  `Io.Evented` with worker threads still crashes in std's `Group` (Uring and `Dispatch`) and Uring
  reuses awaited fibers early (issues.md), so the tests keep Evented single-threaded (Q9).
  `-Drstd-tsan` added (T02).
- [X] **T02** Threading tests beyond phase 2: a test that `Q` never depends on the thread count;
  race stress on aarch64, ppc64le and riscv64. TSan is done (T61: `-Drstd-tsan`, no races on the
  `linux-aarch64` rig). `Io.Evented` with worker threads (Uring, and macOS `Dispatch`) is blocked on
  std's `Group` and fiber-reuse bugs (issues.md; Q9); the multithreaded freestanding build moved to
  T99. threading.md's Verification describes the race loops (about 1000 iterations on the
  weak-memory rigs) as existing, and comparisons.md (8) pitfall 4 asks for the `Q` test; neither
  exists yet (docs63): write them, or mark the text planned. Running (`wip/thrtests`). Done
  2026-10-30, landed d050361d..a7b86e11: `Q` is checked thread-count-free per core (`checkCore`) and
  per call (each threaded `gemm`/`trsm`/`trmm` kernel sees the serial run's depths); race stress for
  levels 1 to 3 (1000 rounds, 2 to 8 workers, bit-checked against one) passes 28/28 on
  `linux-aarch64`, `linux-powerpc64le` and `linux-riscv64`; threading.md's Verification now
  describes only existing tests. Freestanding threading is T99; Evented with workers waits on Q9.
- [X] **T86** Model Loongson 3A's fused `madd.fmt` (T82): on `loongson3r3` and `loongson3r4`
  hardware, rstd's Zig code fuses wherever LLVM picked `madd` and rounds elsewhere (`loongson3r4`'s
  vector bodies round, scalar tails fuse), and `loongson3r3`'s converted kernels fuse; no rig can
  see it, since QEMU rounds `madd`. Thread-count independence holds either way. Model both as fused
  (Q8). Running (`wip/loongfuse`). Done 2026-10-30, landed 819b6c12: `loongson3r3` and `loongson3r4`
  model fused, as the reference's `-march=loongson3a`: scalar multiply-adds as inline `madd.fmt`,
  `loongson3r4`'s MSA vectors as inline `fmadd` (addend tied to the destination, avoiding LLVM's
  operand swap); `p5600` assembly is unchanged instruction for instruction. A native test checks
  fused `axpy`; rigs skip it (QEMU rounds `madd`). `mips64r2el` rig 252 passed, 1 failed (T103, not
  from T86).
- [X] **T09** `linux-powerpc64le` (pwr8, `f128` by library calls) is not hit by the FastISel `stfd`
  bug: drop its `.safe`; the other POWER rigs keep it (T45). Running (`wip/rigs2`). Done 2026-10-30,
  landed a4456735: Debug, 783 passed, 34 skipped, 0 failed; issues.md's workaround text updated.
- [X] **T11** A rig for 32-bit POWER on a VSX CPU (QEMU's POWER8 or POWER9 model running a 32-bit
  image), so the power8 to 10 objects a 32-bit build links are tested (design.md, Settled decisions:
  filling upstream gaps). Running (`wip/rigs2`). Done 2026-10-30 for power8 and power9, landed
  002fedcb: `powerpc-power9` is a 32-bit POWER8 build without `64bit` and `fpcvt` on QEMU's POWER9,
  `test/qemu/ppc32.zig` sets MSR[VSX]; 480 passed, 30 skipped, 3 failed (T103). power10 can't link:
  T105.
- [X] **T68** RVV runs at one VLEN: `riscv64` and `linux-riscv64` are `vlen=512`, so `zvl128b`'s
  VLEN-agnostic kernels (the `axpy` ports, the level-3 tile) never run at 128 or 256, the hardware
  they are named for. Add rig variants at both. Running (`wip/rigs2`). Done 2026-10-30, landed
  45a90ff6: `riscv64-vlen128` and `riscv64-vlen256` (same image as `riscv64`), each 482 passed, 30
  skipped, 1 failed (T103).
- [X] **Q9** Should rstd stop using std's `Io.Group` for its worker team on macOS? On `Io.Evented`
  with worker threads, `Group` crashes on both Uring and `Dispatch` (issues.md). A team built on
  plain futures (`io.concurrent` and `await`) passed 20000 stress rounds on `Dispatch`, but on Uring
  futures hit a separate fiber-reuse bug, so it would help macOS only. (A) Wait for upstream:
  Evented with worker threads stays out of the tests, and `Io.Threaded` (the default) is unaffected.
  (B) Add a futures-based team used on `Dispatch` only, a second code path to keep bit-identical.
  Recommended: (A), since the default `Io` works and (B) is a workaround on one OS that upstream's
  fix would remove. Answered 2026-10-30: (B) as a fallback, reverting to (A) once upstream fixes
  `Group`: a futures-based team on `Dispatch` only, chosen by one comptime switch so that the revert
  is one line, with the std-only `Group` repro kept as a check to rerun on each toolchain bump; T106
  running.
- [X] **T98** The iBook G3 (PowerPC 750FX, no AltiVec) is a RandyOS target, but no rig runs it: a
  750 gets the `unknown` core (the portable kernels), and the only 32-bit PowerPC rig is a G4. Add a
  `powerpc-750` rig (QEMU `-cpu 750`, target CPU `750`) and run it in Debug and ReleaseSafe
  (design.md, Settled decisions: target platforms come first). Running (`wip/ppc750`). Done
  2026-10-30, landed 06576620..b32f33d4: `powerpc-750` (RandyOS's `powerpc-freestanding-eabi`, soft
  float, per-rig `blas_cores = .none`) and `powerpc-750-hf` on `mac99`; QEMU `-cpu 750`, as OpenBIOS
  hangs on `750fx`/`750gx` (issues.md). `powerpc-750` Debug and ReleaseSafe 482 passed, 35 skipped,
  0 failed; `-hf` 483/30/0. A soft-float ABI no longer fuses (`fma.zig`: `@mulAdd` there is
  compiler-rt's soft `fma`, 4 times slower). Only the portable kernels run on a 750, as the
  reference. Follow-ups: `arm-softfp` OOM (T103), RandyOS `blas_cores` (T107).
- [X] **T106** Q9's fallback: on `Io.Evented` with worker threads over macOS `Dispatch`, rstd's
  worker team uses plain futures (`io.concurrent` and `await`) instead of std's crashing `Io.Group`
  (issues.md), chosen by one comptime switch; `Group` everywhere else. Re-enable the
  `Dispatch`-with-workers tests; Uring stays out (its futures hit the fiber-reuse bug). Keep the
  std-only `Group` repro as a script to rerun on each toolchain bump, and revert to `Group` once it
  passes. Running (`wip/dispatchteam`). Done 2026-10-30, landed 9d0a7a1d..236d1668:
  `threading.dispatch_futures` runs worker tasks on futures on `Dispatch`, `Group` elsewhere; the
  threaded tests run on `Dispatch` with workers (3 Mac runs; `linux-aarch64` 47 passed). Revert
  check: `zig build ssh-repro-group-aarch64-macos` (fails until std fixes `Group`; then set the
  switch false and delete the futures path). Found: `Dispatch` moves fibers between threads without
  their FP control state (issues.md); level-3 workers now hold the caller's `FpState` only while
  running (`threading.Lent`). Two Mac failures from before: T108.
- [X] **T108** On the Mac (`ssh-test-aarch64-macos`), "the public gemm threads above the threshold"
  and "the public symm, herk, syr2k and trmm thread" fail, already on review 642e1a12 before T106.
  Both use only `Io.Threaded` and pass on Windows and Linux; macOS is a zfls host. Find why
  (Accelerate routing, which is single-threaded, is the first suspect) and fix the code or the tests
  (scratchpad `dispatchteam/`). Running (`wip/macpublic`). Done 2026-10-30, landed bedd6366:
  Accelerate routing, not threading. On macOS the public `gemm`, `symm`, `herk`, `syr2k` and `trmm`
  go to Accelerate, so the tests compared Accelerate's bits with rstd's serial driver and never ran
  rstd's threading there; they now turn routing off (`rstdOnly`) and restore it. Mac 27/27 twice,
  `linux-aarch64` 27/27. Found T109.
- [X] **T107** RandyOS integration (from T98): RandyOS's `build/rstdlib.zig` passes no `blas_cores`,
  so building rstd for `powerpc-freestanding-eabi` (soft float) stops unless it passes `.none`, as
  the `powerpc-750` rig does. Set it in RandyOS for its soft-float targets (iBook G3, 32-bit arm).
  Done 2026-10-30 in rstd instead, landed e4ff4686: `-Dblas-cores` now defaults to `none` where no
  kernel takes the target's calling convention (`blaskern.defaultCores`), so RandyOS's build needs
  no change; asking for kernels there still stops with the error. Arm soft float was never affected
  (its `sf` variant). RandyOS picks it up when its top-level `vendor/rstd` syncs.
- [X] **Q10** Where should LAPACK past OpenBLAS's own routines (T34 `dgesvd`, eigensolvers, QR) come
  from? (A) Reference LAPACK (netlib; OpenBLAS ships it as `lapack-netlib/`): its LICENSE is BSD-3's
  three conditions plus a patent and third-party IP disclaimer, so the same class as OpenBLAS's;
  complete and the de facto standard, but blocked, not tuned. (B) ReLAPACK (MIT, Elmar Peise;
  already in OpenBLAS's tree as `relapack/`, used by its `BUILD_RELAPACK` option): recursive
  versions of the factorizations (Cholesky, LU, `sytrf`, triangular inverse, Sylvester, generalized
  eigenproblem setup), which its paper reports beat reference LAPACK and match or beat OpenBLAS and
  MKL; builds on BLAS plus LAPACK's unblocked kernels, no SVD or eigensolvers. (C) faer (MIT, pure
  Rust): modern, explicitly SIMD and parallel LU, QR, Cholesky, SVD and eigensolvers; a port from
  Rust is a rewrite, not a translation. Recommended: allow (A) as the baseline for anything not
  ported elsewhere (license already in the accepted class), take (B) for the factorizations it
  covers, as OpenBLAS's own option does, and look at (C) for SVD and eigensolvers if (A)'s are too
  slow; each by benchmark, per design.md's matlab-port rule. Answered 2026-10-30: (B) and (C), no
  reference LAPACK up front: ReLAPACK for the factorizations it covers (T110), and faer for
  everything that would otherwise come from reference LAPACK, ported now for performance out of the
  box (T34, T111); reference LAPACK only later, as a benchmark baseline for faer's versions.
- [X] **T103** The per-core `symm`/`hemm` test ("every linked core's symm and hemm match the wide
  reference") runs out of memory on the kernel rigs: `loongson3r4` on `mips64r2el` (already at
  44e1e683, before T86; scratchpad `loongfuse/rig-symm-base.log`), `zvl128b` and `x280` on
  `riscv64`, `riscv64_zvl128b` on `riscv64-vlen128`/`-vlen256`, and on `powerpc-power9` also the
  per-core `gemm` test and the portable `symm` test (scratchpad `rigs2/`). Every kernel rig's heap
  is a fixed 32 MiB whatever the guest memory, so raising guest memory alone won't fix it: size the
  tests to the heap, free per core, or size the heap to the guest. Also, `arm-softfp` (2 tests in
  `blas3.symm`, already at 4b0f4d9b; hard-float `arm` passes; RandyOS's 32-bit arm is soft float;
  scratchpad `ppc750/`). Blocks a clean T62 sweep. Running (`wip/heapfit`): per-rig heap sizes in
  `build/qemutest.zig`, rig runs left. Done 2026-10-30, landed 77409535..3bd5af33: the level-3 case
  past `R` (`_test/blas3/check.zig`) held three `unroll_m + 4`-row matrices past an `R` that
  2cfd70ff/13768d8a's 128 MB buffer made about 1 MB per row, overflowing the 32 MiB heap for
  `unroll_m` of 6 or more (also `x86_64` and `aarch64`); on freestanding rigs it now takes one row
  (about 12.6 MB for any core). Lost there: crossing the real `R` with a full row panel (the
  shrunk-blocking sweep still crosses a shrunk `R`). The runner now reports heap use, the failing
  request and a heap high water. Symm passes on `powerpc-power9`, `mips64r2el`, `riscv64`,
  `riscv64-vlen128`/`-vlen256`, `arm-softfp`, `x86_64`, `aarch64`. Found T112.

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
