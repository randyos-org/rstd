# BLAS to-do

Open BLAS work, kept here so that findings are not lost. Check items off as they land.

## In progress

- [ ] Role strings become `blas.Slot` everywhere: 55 role-string functions and 1138 generated
  `.role` entries (`wip/slot`).
- [ ] Level 3: `trsm` divides by the diagonal instead of pre-inverting it, and `beta` folds into the
  first `KC` slice (`blas3/driver.zig`, `wip/l3`); comparisons.md (7) and (10).
- [X] Land the AltiVec fix: `util.scalar` keeps portable `f32` loops out of AltiVec lanes
  (`wip/altivec`; issues.md, AltiVec entry).
- [ ] Full qemu rig sweep on the Mac on the latest `review`; check every `sparc64-debug` part fits
  OpenBIOS's 4 MB window (the largest was 0x3f15e8 of 0x400000).
- [ ] Small-n level-2 paths from OpenBLAS (`spr.c`, `syr.c`, `syr2.c`, `zsyr.c`, `zher.c` and kin)
  where they change the serial path; port unless measured slower (`wip/l2small`).
- [ ] Flyoff: rstd against prebuilt OpenBLAS and BLIS (and Accelerate), speed and results, on native
  CPUs only: Windows and Linux x86_64 (WHPX VM), macOS and Linux arm64 (HVF VM) (`wip/flyoff`).
- [ ] Move kernels and functions that do not apply to column-major Zig data (CBLAS row-major paths,
  dead slots) and converted kernels no core uses to `blas/unused/` (git-ignored); the generator
  stops emitting them. After `wip/slot` lands.
- [ ] gemm tails: go back to OpenBLAS-style unrolled edge tiles; the `thin` tails were a wash on
  speed and grew the objects (blas3/gemm/notes.md, Tails). After `wip/l3` lands.

## Correctness and verification

- [ ] `zrot_ppc440.S` (`ppc440fp2`) has a zero-stride patch no rig has run; no rig runs `ppc440` or
  `loongson3r4` either (they share code with `ppcg4` and `loongson3r3`).
- [ ] powerpc64le `.safe` workaround for LLVM's `-O0` `fcmpu`: it looks fixed; verify Debug on
  `powerpc64le`, `powerpc64le-power10`, `powerpc64` and `linux-powerpc64le`, then drop it.
- [ ] The `arm` rig may test soft-float: Zig's `cortex_a15` model has no NEON or VFP features, so
  `f32` compiles to library calls; check and add the features (`build/qemutest.zig`). Its build
  also prints `error: <inline asm>: instruction requires: VFP2 (vmsr fpexc)` from
  `test/qemu/arm.zig`, though the step succeeds.
- [ ] No `powerpc64` (970) AltiVec scan yet: LLVM stops with "Invalid global physical register".
  Scan once it builds; also rerun the G4 `-femit-asm` scan after LLVM upgrades, since the SLP
  vectorizer is not proven blocked by `util.scalar`.
- [ ] Test-only code still auto-vectorizes on the G4 (`blat1` rot, the scal, gemv and half-gemv
  models, `triangular.checkAll`, `blat3`'s complex type); fine for normal-range data.
- [ ] Cores no rig runs: the Skylake-X family (AVX-512) and the Bulldozer family (FMA4, XOP). Run
  `x86_64-whpx` on a host that has them (design.md, Open issues).
- [ ] Hardware-only paths never run: the mips `cpucfg` decode, the sparc `sllx` rewrite, Linux ppc64
  `mfpvr`, the Windows arm64 registry MIDR and feature probes, the arm64 sysfs fallback; NEON and
  AltiVec flushing was seen on rigs only.
- [ ] Variants never run: 32-bit SPARC (lld cannot link it; includes the `%icc` index-tie fix) and
  32-bit x86 kernels on Windows and macOS, where complex dot returns a struct (verify).
- [ ] The init race test's `single_threaded` skip has not run on a freestanding target (verify).
- [ ] `f80` portable kernels off x86 (Zig calling convention) are untested.
- [ ] Level 2 has no large-stride cases (`_test/blas2/edge.zig`).
- [ ] `blas2/buffers.md` says `zgemv_t.S` needs a NaN-free buffer; the NaN-poisoned tests never
  showed it. Confirm, or drop the interface's x86_64 clear of that buffer.
- [ ] i386 x87 `xgemv_n` column blocks per core are checked only indirectly; add a test that calls a
  slot built for one core (`blas2/gemv`, `cpu/x86.zig`).
- [ ] Linux rig gaps: `linux-powerpc64le` never runs POWER10 kernels (baseline build),
  `linux-aarch64` leaves out the SVE cores, and under HVF it also skips `vortexm4` (cause unknown).
- [ ] Zig exits 0 on an inline-assembly error and writes an empty object. Confirm the link then
  fails on the missing symbols, or add a symbol check to `build/blaskern.zig` (`c_kernels.md` says
  to check by hand).
- [ ] Rank updates skip zero columns (`ger`, `syr`, `her` and kin); follow reference if it drops the
  skip, and note `ger`'s run splitting may double the cost with strided `x` and many zero `y`
  (comparisons.md (5)).

## Performance

- [ ] `trsm` with a unit diagonal now divides by 1 (exact, but a division per element); skip it
  for `Diag.u` (`blas3/trsm/kernel.zig`).
- [ ] Speed follow-ups: a `@Vector` `sum` (`csum`/`zsum` are scalar on the old x86_64 cores);
  measure the scalar fallbacks (extremes and `sdot` on arm, pre-VSX POWER, mips64, SPARC) and the
  f64 `sdot` (comparisons.md (1)).
- [ ] Kernel strategy: race `@Vector` against the converted assembly per core on a quiet machine
  (only vortex is done, `tools/blasbench/results-m6.md`); keep per-core code only where Zig loses.
- [ ] C and intrinsics kernels not yet ported, levels 1 and 2 (about 1120 and 776 cells: riscv64
  RVV, zarch, POWER8 to 10, mips MSA, x86_64 `*_4.c` gemv, symv and ger, wasm SIMD128): fill them by
  that strategy, `@Vector` first, as arm64 did (`vector.zig`).
- [ ] Same for the bfloat16 intrinsics (cooperlake, sapphirerapids, Neoverse N2, V1, V2, POWER8 to
  10, riscv zvl256b); `sbdot` could call the core's own `sdot`.
- [ ] A tuned `@Vector` `rot` for vortex: the assembly still leads by 2 to 17% in L1.
- [ ] Race the Zig complex gemm against OpenBLAS's vortex assembly on the M6.
- [ ] Assembly where `@Vector` cannot reach: arm64 SVE, SME (`vortexm4`), POWER10 MMA, AMX and
  AVX-512 BF16 for `sbgemm`; Skylake-X only if under 85% of its kernel (design.md, Implementation
  plan 6; `blas3/gemm/notes.md`).
- [ ] BLIS assembly kernels, decided per kernel where benchmarks favor them (BSD-3, with a license
  section and per-file provenance).
- [ ] gemm tuning: faster Zig packing (whole-vector transposes in `oncopy`, `incopy`), `C`-tile
  prefetch at kernel entry, a larger `KC` on Apple cores (comparisons.md (9)).
- [ ] Measure whether composed `symv`/`hemv` pay off on cores with assembly gemv (atom, opteron,
  arm64); today they compose only where the core has its own gemv.

## Features and coverage

- [ ] `gemm3m`, `symm3m`, `hemm3m` (stubs kept), a small-matrix path, bfloat16 and extended level 3
  (design.md, Implementation plan 4).
- [ ] Extension interfaces and portable kernels: `geadd`, `imatcopy`, `omatcopy`, `laswp`, `neg`,
  `gemmt`, `bgemm`; their slots are null and `utest_ext` cases wait on them (plan 5).
- [ ] LAPACK, starting with what the OpenBLAS suites wait on: `potrf`, `potrs`, `dgesvd`.
- [ ] `scabs1`/`dcabs1` (`cabs.S` upstream) are missing.
- [ ] riscv `zvl256b` has no object (`zvl128b` serves it): read VLEN and port the `*_vector.c`
  kernels (design.md, Open issues).
- [X] Complex `max`/`min` stay a compile error, now pointing to `.amax`/`.amin`: reference BLAS and
  BLIS have no signed complex extreme, and OpenBLAS's `cblas_icmax` scans the floats as reals
  (issues.md; d14858d).
- [ ] `sbger_k`, `sbsymv_L`, `sbsymv_U`: dead upstream too (only `common_param.h`'s table names
  them; no interface, driver or kernel). Move their stubs to `unused/` with the rest.
- [ ] `gerv`/`gerd` serve only OpenBLAS's row-major CBLAS `zgerc`: move them to `unused/`; keep
  upstream's row-major branches as comments in our interfaces, never as code.
- [ ] `coverage.md` exists only for `axpy`, `dot`, level 2 and `gemm`; generate the rest with
  `tools/blas_asm/survey.zig`.
- [ ] Level-2 assembly for the 32-bit-only `ppc440fp2` and `ppcg4` still carries 64-bit variant
  keys; prune (verify).
- [ ] arm64 gemv assembly saves and restores `x18`, the Windows and Apple platform register, as
  upstream does; document or guard it.
- [ ] Converted kernels no core uses anymore: `nrm2` and the extremes on every arch, `sdot`, and the
  vortex assembly in non-Apple objects. Move them to `blas/unused/` (git-ignored, kept locally) and
  stop the generator emitting them, once `wip/slot` lands.

## Threading

- [ ] Rest of phase 1: the cached CPU count and ceiling in `blas.config`, thresholds, chunking, the
  claim loop and ordered combine, with unit tests (threading.md, Implementation plan 1).
- [ ] Level-3 threading prototype from BLIS (tag 2.1), close to clean room, against OpenBLAS's
  `level3_thread.c`, before phase 2 (threading.md, Open issues).
- [ ] Phases 2 onward: real `gemm`, level 1, level 2, the rest of level 3 (threading.md).
- [ ] Thresholds: OpenBLAS's defaults and per-core overrides now; tune after benchmarks.
- [ ] Tests: threaded against serial, bit for bit, at 1 to 9 threads; a test that `Q` never depends
  on the thread count; race stress on aarch64, ppc64le, riscv64.
- [ ] Unverified: the futex protocol on `Io.Evented`; TSan with this toolchain.
- [ ] `tools/blasbench` runs rstd with `Io.failing`; revisit once rstd threads.

## Accelerate

- [ ] Part 2: route the public routines to it on macOS: the order enum, 0-based `i?amax`, complex
  arguments by pointer, the `_sub` dots, 32-bit range checks on the LP64 table, rstd's argument
  checks before any call, quirk entries where its results differ (design.md, Apple Accelerate).
- [ ] The off switch stops using the table and never unloads; `singleThreadHere` on every calling
  thread; `load` during init, since `setenv` is not thread-safe; tie it into rstd's threading.

## Tooling and tests

- [ ] The AVX-512-off workaround for the self-hosted backend (`blas3/gemm/kernel.zig:56`) has no
  issues.md entry.
- [ ] `blas1/c_kernels.md` says the AVX-512 and FMA4 cores need `qemu-test-x86_64` with a matching
  `-cpu`; TCG has neither, so point at `x86_64-whpx`.
- [ ] `sparc64`: replace the image split with a real fix, mapping more memory at boot through the
  Open Firmware client interface (claim and map), so one image fits.
- [ ] Linux rigs: `linux-i386` is a stub (QEMU x87 page-crossing bug); `s390x` and `mips64el` wait
  on `zig fetch` taking single files.
- [ ] The `powerpc-macho` rig is a stub; finish it once the powerpc-macos Zig port lands.
- [ ] `-Dfail-fast` has not been seen to stop a real failing rig (verify).
- [ ] Check in the C-against-Zig differential harness (stub `common.h`) used for `dot`.
- [ ] `build/blaskern.zig` has three `catch @panic("OOM")`; use `try`.
- [X] `fmadd` is duplicated in `blas1/axpy/kernel.zig` and `blas1/dot/kernel.zig`; merge
  (`fma.madd`, 3dd9d35).
- [ ] macOS VM extras: base-image overlays (`qemu-clone`, `-Dmac-base`), a `Forgejo` self-hosted
  runner, a setup that needs no clicking, less guest load when headless (HID driver,
  `SecurityAgent`), and freeze detection that asks before killing.

## Upstream

- [ ] Keep recording upstream problems in issues.md as they are found; the user triages them (filing
  included) once the BLAS work is done.
- [ ] Recheck on each Zig update and drop the workaround once fixed: the powerpc64le two-float
  return (`zdot` patch), arm soft-float struct returns, f80 under `x86_64_sysv` on Windows, the
  comptime naked-function call, ppc64 big-endian ELFv2, the two-lane `f32` vector stored as 16 bytes
  (array store in the complex tile), sparc64 `_Qp_sqrt` (`check.sqrtW`; rstd's own f128 square roots
  stay f64-accurate there), sparc64 `-Oz`, lld and 32-bit SPARC.
- [ ] Also recheck: `SafeAllocator`'s CPU count, the `Thread` override doc, `zig fetch`, run-step
  timeouts and output, the interrupted-link cache entry.
- [ ] translate-c: re-enable the package once its bundled Aro builds on current Zig.

## Zig-idiomatic interface

Decided, not built (design.md, "Planned: a Zig-idiomatic interface").

- [ ] Slices for every operand instead of pointer and size arguments; routines check that the slices
  cover what sizes and strides need.
- [ ] Exact error sets per routine or family, with plain names (`SliceTooShort`, `ZeroStride`,
  `LeadingDimensionTooSmall`); `Canceled` only if cancelable, `OutOfMemory` only if allocating.
- [ ] Errors, not panics, for bad input; checks no longer need the reference's order.
- [ ] Doc comments say where a routine panics explicitly (only unrecoverable states, above all
  possible memory corruption from assembly).
- [ ] User docs without BLAS jargon (`xerbla`, `LDA`, `INCX`).

## Zigify the codebase

- [ ] The last item: an rstd-wide Zig-idiomatic pass (slices, exact error sets, plain error names,
  documented panics). The BLAS interface items above are part of it.
