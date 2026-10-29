# BLAS to-do

Open BLAS work, kept here so that findings are not lost. Check items off as they land.

## In progress

- [ ] Optional: an inline-asm `maddf` for R6 scalar multiply-adds in the main library, if R6 speed
  ever matters (LLVM calls `fmaf`).
- [ ] Native `x86-windows-gnu -Dcpu=pentium3`: `thread_proto`'s "prototypes match the serial driver
  bit for bit" fails (`f64 w=2 m=1 n=16 k=13`), already at 6c6def0; likely x87 excess precision
  between the threaded and serial paths.
- [ ] The ISA check cannot see inline assembly (x86's assembler ignores CPU features); and
  detection, as OpenBLAS's, hands some cores to weaker CPUs: `katmai` to any 32-bit x86 without SSE3
  (down to the 486), `vortex` to any Apple CPU (pre-A12 lacks FCMLA), `vortexm4` to non-Apple SME
  without SVE, a forced core's fallback to `neoversen1` on ARMv8.0.
- [ ] Zen 3 (flyoff): short real `dot`s (0.30 of OpenBLAS at 1024) and complex `gemm` (0.79 to 1.01)
  trail OpenBLAS.
- [ ] Skylake-X family `sgemv_n`/`sgemv_t`: those cores add large intrinsics kernels on top of
  Haswell's (`sgemv_t` uses a 1,121-line template); they keep the portable `sgemv` until ported.

## Correctness and verification

- [ ] `zrot_ppc440.S` (`ppc440fp2`) has a zero-stride patch no rig has run; no rig runs `ppc440` or
  `loongson3r4` either (they share code with `ppcg4` and `loongson3r3`).
- [ ] POWER rigs stay at `.safe` (dcb4628): the `-O0` `fcmpu` compare is fixed on 0.18.0-dev.751
  (now `xscmpudp`), but FastISel still stores an f64 held in a v-register with `stfd` (reads f2 for
  v2); seen in `_test/blas1/edge.zig` `contractCase`, 10-line reproducer in issues.md. All four rigs
  pass in Debug today (400/24/0, linux 431/30/0). Drop `.safe` once LLVM fixes the store;
  `linux-powerpc64le` (pwr8, f128 by library calls) looks unaffected and could drop it now.
- [ ] No `powerpc64` (970) AltiVec scan yet: LLVM stops with "Invalid global physical register".
  Scan once it builds; also rerun the G4 `-femit-asm` scan after LLVM upgrades, since the SLP
  vectorizer is not proven blocked by `util.scalar`.
- [ ] Test-only code still auto-vectorizes on the G4 (`blat1` rot, the scal, gemv and half-gemv
  models, `triangular.checkAll`, `blat3`'s complex type); fine for normal-range data.
- [ ] Cores no rig runs: the Skylake-X family (AVX-512) and the Bulldozer family (FMA4, XOP). Run
  `x86_64-whpx` on a host that has them (design.md, Open issues).
- [ ] The new x86_64 `rot`, `scal` and `gemv` paths for the Skylake-X and Bulldozer families are
  compile-only; include them in that `x86_64-whpx` run.
- [ ] Hardware-only paths never run: the mips `cpucfg` decode, the sparc `sllx` rewrite, Linux ppc64
  `mfpvr`, the Windows arm64 registry MIDR and feature probes, the arm64 sysfs fallback; NEON and
  AltiVec flushing was seen on rigs only.
- [ ] Variants never run: 32-bit SPARC (lld cannot link it; includes the `%icc` index-tie fix) and
  32-bit x86 kernels on Windows and macOS, where complex dot returns a struct (verify).
- [ ] The init race test's `single_threaded` skip has not run on a freestanding target (verify).
- [ ] Nothing runs `XDouble` as `f128` (`quad_precision` is hard-wired off): add a build option and
  a rig configuration for it. The edge tests and level 3's per-routine driver tests still cover no
  extended type off x86 (the interfaces do).
- [ ] x87: other multistep arithmetic (the complex `mul`/`sub` updates) can change bits with where
  LLVM spills; one code path stays self-consistent, but two routines doing the same arithmetic may
  not. Audit where rstd promises equal bits across routines on x87 targets.
- [ ] 32-bit builds: `Makefile.system` gives some targets another core's getarch entry (x86 haswell,
  skylakex, cooperlake, sapphirerapids, sandybridge, zen take nehalem's; the bulldozer family
  barcelona's; power8 to 10 power6's; armv8 armv7's); model it together with the matching kernel
  choice. Also, the power64 family assumes little-endian, so a big-endian ppc64 build takes the
  wrong ppc970 unrolls.
- [ ] Linux rig gaps: `linux-powerpc64le` never runs POWER10 kernels (baseline build),
  `linux-aarch64` leaves out the SVE cores, and under HVF it also skips `vortexm4` (cause unknown).
- [ ] Rank updates skip zero columns (`ger`, `syr`, `her` and kin); follow reference if it drops the
  skip, and note `ger`'s run splitting may double the cost with strided `x` and many zero `y`
  (comparisons.md (5)).

## Performance

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
- [ ] Once LAPACK `getrf` lands, call `laswp_ncopy` and `neg_tcopy` from it, or drop them as
  upstream's `#if 1` does.
- [ ] LAPACK, starting with what the OpenBLAS suites wait on: `potrf`, `potrs`, `dgesvd`.
- [ ] riscv `zvl256b` has no object (`zvl128b` serves it): read VLEN and port the `*_vector.c`
  kernels (design.md, Open issues).
- [ ] `coverage.md` exists only for `axpy`, `dot`, level 2 and `gemm`; generate the rest with
  `tools/blas_asm/survey.zig`.
- [ ] Level-2 assembly for the 32-bit-only `ppc440fp2` and `ppcg4` still carries 64-bit variant
  keys; prune (verify).
- [ ] arm64 gemv assembly saves and restores `x18`, the Windows and Apple platform register, as
  upstream does; document or guard it.
- [ ] `tools/blas_asm/arch` patches of sources no longer converted (the extremes', `nrm2`'s, x86_64
  `zdot_atom.S`) never run now, so nothing checks they still match: move them to `unused/` too, or
  keep them as the record `issues.md` points to.

## Threading

- [ ] Threshold overrides left out of phase 1: z13 `axpy` and Windows-on-ARM64 `gemv`
  (`threading.zig` `rule`); add them if benchmarks agree.
- [ ] Unmeasured: B on A's `CNAME` grid, and B multiplying its own chunks as it packs them, before
  the barrier (passes the tests, uncommitted: `scratchpad/agents/l3proto/thread_proto_fused.zig`,
  binaries `l3bench_bgrid.exe`, `l3bench_bfused.exe`); time them when the machine is quiet. The
  prototypes cap at 64 workers (`chunks_max`).
- [ ] Phases 2 onward: real `gemm`, level 1, level 2, the rest of level 3 (threading.md).
- [ ] Thresholds: OpenBLAS's defaults and per-core overrides now; tune after benchmarks.
- [ ] Tests: threaded against serial, bit for bit, at 1 to 9 threads; a test that `Q` never depends
  on the thread count; race stress on aarch64, ppc64le, riscv64.
- [ ] Unverified: the futex protocol on `Io.Evented`; TSan with this toolchain.
- [ ] `tools/blasbench` runs rstd with `Io.failing`; revisit once rstd threads.

## Accelerate

- [ ] The Accelerate level-2 special-value sweep takes minutes per Mac run (complex, even trimmed to
  every third variant at the largest shape); speed it up.
- [ ] Tie Accelerate into rstd's threading (threading.md, Vendor backends).

## Tooling and tests

- [ ] The AVX-512-off workaround for the self-hosted backend (`blas3/gemm/kernel.zig:56`) has no
  issues.md entry. Its reason was never recorded (since 2cf8757); the self-hosted backend compiles
  512-bit `@Vector` multiply-adds without error on 0.18.0-dev.751. Run the gemm tests with
  `-fno-llvm` on AVX-512 hardware (this host is Zen 3; TCG has none), then drop it or record it.
- [ ] Zig 0.18.0-dev.751 gives four x86 models LLVM names LLVM does not know (`skylake_avx512`,
  `icelake_client`, `icelake_server`, `goldmont_plus`: underscores for hyphens), so those cores
  likely build without CPU tuning (issues.md). Check the per-core objects' tuning, and map the names
  in `build/blaskern.zig` if so.
- [ ] `sparc64`: replace the image split with a real fix, mapping more memory at boot through the
  Open Firmware client interface (claim and map), so one image fits.
- [ ] Linux rigs: `linux-i386` is a stub (QEMU x87 page-crossing bug); `s390x` and `mips64el` wait
  on `zig fetch` taking single files.
- [ ] The `powerpc-macho` rig is a stub; finish it once the powerpc-macos Zig port lands.
- [ ] Check in the C-against-Zig differential harness (stub `common.h`) used for `dot`.
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

## Open questions for the user

Items set aside overnight for a decision; each says what is blocked and where its work sits.

- [X] `arm-softfp` rig: on Zig 0.17 its `cortex_a15` had `vfp4`, so it built hardware float under
  the soft-float ABI (true softfp); on 0.18 it is fully soft-float. Pin `vfp4`/`neon`/`d32` so it is
  softfp on both toolchains, or keep it fully soft? Its kernels' `sf` variants run either way
  (553c63a). Decided (2026-10-29): softfp on both toolchains; features pinned as in `arm`
  (design.md, Settled decisions).
- [X] x86_64 `dot`: on the old x86_64 cores the alignment-peeling `ddot`/`cdot`/`zdot` slots get the
  portable kernel, not `@Vector`, because `blas1/dot/x86_64/kernel.zig` never exposed the converted
  list (now `skips`). `readme.md` says `@Vector`. Adding `pub const skips = K.skips;` there matches
  the docs but changes results on those cores. Do it, or fix the docs? Decided: `@Vector`, per the
  2026-10-23 no-address-dependence decision; `dot` missed it by an oversight (design.md, Settled
  decisions).
- [X] Hermitian `Uplo` letters `v`/`m` (`HEMVREV`): upstream reaches them for `hpmv`, `her`, `hpr`,
  `her2`, `hpr2` only from CBLAS row-major (Fortran takes V/M only for `zhemv`/`zhbmv`). Per the
  column-major-only direction, drop them from those five (a public API change), or keep them?
  Decided (2026-10-29): drop everything row-major-only; `v`/`m` go from those five (`wip/colmajor`).
- [X] Testing bar for mechanical refactors: native tests plus `ast-check` let the Slot rewrite land
  with a broken aarch64-sve build (fixed in 275c67b). Add `zig build qemu-check` (compile-only,
  every rig image, about a minute, no QEMU) to that bar? Decided (2026-10-29): targeted only; no
  routine full `qemu-check`; per-rig `qemu-check-<rig>` steps for the touched arches (design.md,
  Settled decisions).
- [X] Mac timing: `audiomxd` (~72% CPU) and `configd` (~42%) have spun for 10 days, holding the M6's
  load near 1.6, so the flyoff relaxed the quiet rule (no `qemu-system`, all else under 5% of a
  core; documented in `results-flyoff-macos-aarch64.md`). OK as the rule going forward, or will you
  reboot the Mac and have its timings rerun? Decided (2026-10-29): the quiet rule stays; the Mac was
  rebooted and its flyoff timings are being rerun (`wip/retime`).
- [X] Level-3 threading protocol for phase 2: (B) BLIS barriers as recommended (half the code, one
  barrier per slice, a one-line deadlock argument), or (A) OpenBLAS's flags (0 to 17% faster on Zen
  3 today; threading.md still describes A)? If B: one barrier with double buffering (twice A's `B`
  memory, as prototyped), or BLIS's two barriers around one buffer? Phase 2 waits on this. Decided
  (2026-10-29): BLIS, with BLIS's two barriers around one buffer (deferring where OpenBLAS and BLIS
  agree); phase 2 started (`wip/l3thread`).
- [X] Allow a short bounded spin before `futexWait` in level-3 threading? The no-spin rule
  (threading.md, Waiting) keeps fiber `Io`s unblocked; a spin might help small sizes. Decided
  (2026-10-29): spin, as OpenBLAS and BLIS both do.
- [X] Accelerate: for `herk`, `her2k`, `her`, `hpr`, `her2`, `hpr2`, Accelerate leaves the
  diagonal's imaginary parts as `beta` times their old value (or NaN with an infinity), where the
  reference zeroes them. Rstd now zeroes them after the call, giving the reference's bits, rather
  than routing those routines to rstd. Keep the fix-up, or route them? Also, x86_64 macOS was
  checked only under Rosetta 2. Decided (2026-10-29): the mathematically correct zero; the fix-up
  stays.
- [X] `armv9sme` spec adds `sme_f64f64` (OpenBLAS's `-march=armv9-a+sve2+sme+sme-f64f64`), an
  optional SME extension `runnable` does not check; no rstd kernel uses SME yet. Drop it from the
  spec, or require it at run time (new probes for Linux, macOS, Windows, freestanding)? It sits in
  the ISA check's `undecided` list meanwhile (bd6dfe8). Decided (2026-10-29): keep it, and require
  it at run time so that no core is chosen that can fault (`wip/smef64`).
