# Linear algebra to-do

The backlog of `rstd.linalg` (the BLAS, LAPACK, the public layer, the matlab port, shaderlib's
vector math) in priority order: work goes top to bottom. IDs stay fixed when items move.
**Needs:** names items that must land first; **Blocked:** names an outside cause (hardware, an
upstream fix, a decision). Upstream problems go to issues.md as they are found; the user
triages them after the BLAS work.

## Feature complete means

**MVP first (user, 2026-10-30):** support tiers as Zig's ([`support.md`](../../support.md)): **tier 1** is `x86_64-windows`, fully working with everything implemented so far, as soon as possible, and SPIR-V (`spirv32-vulkan`), whose work (S1, S2) waits until that CPU BLAS works; **tier 2** is every other RandyOS and zfls platform (RandyOS's x86_64, aarch64, 32-bit arm and PowerPC 750 freestanding targets and its UEFI build chain; zfls's x86_64 and aarch64 Linux, aarch64 Windows and x86_64 and aarch64 macOS 11+, with rstd's Accelerate backend on macOS), whose bugs are also required, after tier 1; **tier 3** is every other platform, best effort; **tier 4** is unsupported or blocked. Tier 3 and 4 bugs, follow-ups, tuning, optional work and Blocked items live in Known issues and enhancements, after the MVP, and the Queue gets no new items.

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

## In progress

One task at a time, worked with the user in the trunk on `review` (user, 2026-10-30). Unfinished
work stays there as unstaged edits, and the item says what is left, so nothing lives only in a chat
or a scratchpad.

### Working on

Nothing yet: next is M2, after a check-in with the user (2026-10-31).

### In the background

Sweeps by `.agent/sweep.sh` (gitignored), results in `.agent/sweeps/<UTC>-<tier>/summary.md`;
started 2026-10-30 22:20 on `review` faeea2aa: tier 1, and tier 2 then tier 3. These two ran from
the trunk, so steps after their first compile what has landed since (T124, M1, kern3, T127); from
tier 3 on, each tier runs from its own worktree, `C:\scratch\gitstd-wt\sweep-<tier>`, pinned
at review's tip when it starts (user, 2026-10-30).

- [ ] **T62** Main-merge gate: once T55 and T77 land and R4 has run, one full rig sweep on review,
  plus `wasm-test`, `rstd-check-macho`, `rstd-check-uefi` and `rstd-check-zig017` (Zig 0.17.0,
  `-j4`); report whether dev is clean to merge to main. **Needs:** the items running now to land, so
  one sweep covers them (T55, T77, T91 and R4 done). MVP order: tier 1 first: native
  `x86_64-windows` (`rstd-test`, with `-Dblas-cores=all`, and the `x86_64-whpx` rig for the x86_64
  cores), fixing every failure. Then tier 2: `x86_64`, `x86_64-uefi`, `x86_64-smp` (T99),
  `linux-x86_64`, `x86_64-macos` (Intel Mac), `aarch64`, `aarch64-uefi`, `linux-aarch64`,
  `ssh-test-aarch64-macos`, `arm`, `arm-softfp`, `linux-arm`, `powerpc-750`, `powerpc-750-hf`, and
  `rstd-check-zig017` (zfls's Zig 0.17.0 targets, `aarch64-windows` included; its tests run only in
  zfls's CI), fixing every failure. Then tier 3 (`i386` and `rstd-test-x87` with their ReleaseSafe
  forms, `wasm-test`, every other rig) once for the record; its failures go to Known issues and
  enhancements. Tier 1 passed 2026-10-31
  (`.agent/sweeps/20261031T032046Z-tier1`): `native` (`-Dblas-cores=all`, faeea2aa) 963 passed, 32
  skipped, 0 failed; `x86_64-whpx` (after T124, M1, `kern3` and T127 landed) 643 passed, 42
  skipped, 0 failed; no inline-asm errors.

#### Build slots

Compiling, one row per slot; a native test compiled on the Mac and run here shows in both tables.

| Slot | Item | Job | Since |
| --- | --- | --- | --- |
| Local 1 | free | | |
| Local 2 | free | | |
| Mac 0.18 | free | | |
| Mac 0.17 | free | | |

#### Running tests

Test executables and QEMU rigs, outside the slots.

| Where | Item | Run | Since |
| --- | --- | --- | --- |
| Mac | T62 tier 2 | step 2 of 26, `x86_64-uefi` (step 1, `x86_64`: 640 passed, 1 failed, T130) | 01:18 |

## Queue

### Testing: sweeps and checks of landed code

Each moves to In progress when the sweep step that runs it starts.

- [ ] **T08** LAPACK's own tests (`potrf`, `getrf`, `getrs`) on the rigs: add `linalg/lapack` to
  `src/_test_qemu.zig`, with a sparc64 part. Landed (8c574827): `linalg.lapack.` in the rigs'
  default filters and `_test_qemu.zig` (so `rstd-check-uefi` links it), sparc64's LAPACK in the
  `blas3-tri` part, sparc64-debug in four `lapack-*` parts (`getrs`'s every-core test split into
  real and complex to fit the load window; `getrf` and `cholesky` have about 140 and 240 KiB left);
  every listed `qemu-check` passes. Check off once T62's sweep runs them; the 32 MiB rig heap is
  unproven. Scheduled: tier 2 sweep steps 18 to 21
  (`T08-<rig>`), started by the sweep itself.
- [ ] **T114** Accelerate's `potrf` (T35) is untested on the Mac: which symbol spelling Apple
  exports (`load` tries `dpotrf$NEWLAPACK$ILP64` and `dpotrf_$NEWLAPACK$ILP64`), whether its special
  values match the reference (a quirk entry if not), and in the timing session whether routing beats
  rstd's threaded `potrf`. Sweep filters: `potrf goes to Accelerate`, `potr`. Tier 2 (Accelerate on
  macOS, user 2026-10-30): the correctness part is required; the timing part waits for the timing
  session. The 13.3 SDK answers the spelling: `dpotrf$NEWLAPACK$ILP64`, no underscore (T47); the
  special values are left. Scheduled: tier 2 sweep
  steps 22 and 23 (`T114-<host>`), started by the sweep itself.

### Code and design

- [ ] **M2** Replace the matlab port's own BLAS and LAPACK-level code with `rstd.linalg.blas` and
  `rstd.linalg.lapack`, keeping only what M1 showed better. Decide first how to bridge its split
  real/imaginary complex layout (`xr`/`xi`) to the BLAS's interleaved `Complex`: convert at the
  boundary or change the port's layout. Its LINPACK and EISPACK routines map to LAPACK's `getrf`,
  `getrs`, `getri`, `gecon` (`wgefa`, `wgesl`, `wgedi`, `wgeco`), `potrf` (`wpofa`), `geqp3`,
  `unmqr`, `trtrs` (`wqrdc`, `wqrsl`), `gesvd` (`wsvdc`, T34), `gehrd` (`corth`), `hseqr` with
  `trevc` (`comqr3`), `hetrd` and `unmtr` (`htridi`, `htribk`); rstd has `potrf`, `getrf` and
  `getrs`; T70 adds OpenBLAS's own, and the rest come from Reference-LAPACK (BSD-3, as OpenBLAS
  ships it in `lapack-netlib/`). Its fills (`wset`, `rset`) are LAPACK `laset`, `pythag` is `lapy2`.
  Its own tests must still pass. **Needs:** M1 (done), T123, T125, T126.
- [ ] **M3** Build `rstd.linalg`'s public layer, the BLAS and LAPACK behind it, as designed in
  design.md's "Planned: the public `rstd.linalg` layer" (settled in M5, 2026-10-30):
  `StaticMatrix(R, C, T)` (comptime shape, a value) and unmanaged `Matrix(T)` (run-time shape),
  real and complex, one operation set over views, inline `@Vector` code up to 8 per dimension and
  the BLAS and LAPACK above by size alone, SPIR-V's own instructions in shaders; the `rstd.simd`
  split and its callers; its migration plan's seven steps, each a commit;
  `src/linalg/readme.md`, the layer's user guide, linking down to `blas/readme.md`. **Needs:** nothing
  (M1, S1, M5 done); M4 confirms the inline limit after step 4.
- [ ] **T21** Zen 3 (flyoff): short real `dot`s (0.30 of OpenBLAS at 1024) and complex `gemm` (0.79
  to 1.01) trail OpenBLAS. Tier 1 (user, 2026-10-30): re-measure after T100, then fix what's left;
  the timing session verifies. Diagnosed 2026-10-30, no code yet: Zen's `sdot` kernel
  (`Real(f32,4).dotSum`) is bound by its conversions (32 `vcvtps2pd` per 64 elements, about 70
  GB/s: the cost of the `f64` rule), and on AVX2 its two `@Vector(32, f64)` accumulators fill all
  16 ymm registers (about 25 moves and 3 spills per iteration): halve them. `cgemm`/`zgemm`: tiles,
  `P`/`Q` and the inner loop match OpenBLAS's, but `cgemm_kernel_8x2_haswell.c` reads two or three
  B panels at once (m8n4/m8n6) with 12 accumulators where ours has 8, exactly FMA latency times
  throughput: port that shape. `blas_bench`'s `kern3` section (59faab6c, any host) times the
  kernel alone against rstd's whole `gemm` on one worker; its first run (busy Zen 3, 2026-10-30)
  has the driver keeping 0.89 to 0.98 of the kernel for `s`, `d` and `z`, and `cgemm`'s whole call
  at 1.04 to 1.08 of its kernel block: the `cgemm` gap is the kernel's. For the timing session:
  `kern3`, `flyoff.sh zen3`'s `win-1t`/`lx-1t` `dot s` rows and the level-3 `c`/`z` rows.
- [ ] **T102** Strided `saxpy` and `zscal` (inc 2) in the arm64 Linux VM halved between 6d66239 and
  28171adc for both rstd and the converted kernel (`saxpy` 54 to 28 GB/s, `zscal` 139 to 71; R4,
  results-flyoff.md). Read 2026-10-30, not measured: nothing on rstd's path changed (one kernel call
  below the threshold, vortex's slot is still the displaced `@Vector` kernel, `vector.zig`'s `axpy`
  and `scal` and the per-core flags unchanged). The bench's `asm` column reads the same function
  rstd calls, so the two moved together, while its own copy of that source (`vec x4`) stayed at 55
  and OpenBLAS's `axpy.S` at 28: the drop follows that one copy in the binary (likely code
  placement) or the guest. Left, in the timing session: `flyoff.sh m6`'s inc-2 `saxpy`/`zscal`
  rows, then align the copy or close this as the guest's.

## Known issues and enhancements (after the MVP)

Not needed for the MVP (see Feature complete means): bugs off the RandyOS targets, best-effort ports and rigs, follow-ups, tuning, optional work and Blocked items. Picked up after the MVP, in this order.

- [ ] **T131** dqds for a bidiagonal's values alone (from T126): Reference-LAPACK's `xLASQ1` to
  `xLASQ6` (`dlasq1.f` to `dlasq6.f`, about 1,070 lines of code), which `xBDSQR` calls when no
  vectors are wanted: faster than the QR iteration and accurate to a few eps relatively. Now
  `gesvd` and `gesdd` with `.n` run `xBDSQR`'s QR sweeps, already accurate (T126's test).
- [ ] **T120** faer follow-ups (split from T117, 2026-10-30): `@Vector` the fused bidiagonalization
  pass (scalar now, faer's four partial sums); run divide-and-conquer's two halves in parallel if
  rstd's threading gains task trees; port faer's multishift QR with aggressive early deflation
  (n >= 75) and blocked Hessenberg reduction (n >= 256), so large `geev` isn't slow; thread the fused
  tridiagonalization and Hessenberg passes with fixed chunks; make `sytrd`, `gehrd`, `hseqr`,
  `trevc` public; faer routines beyond what M2 needs (T123 has those): Cholesky rank-1 update,
  pivoted Cholesky, full-pivoting LU, `ggev`/`gges`; race `gesvd`/`gesdd` against OpenBLAS's LAPACK,
  reference LAPACK and Eigen in the timing session (T73's bench); split sparc64-debug's
  `lapack-faer` part (tier 3; `openblas-l1` at 0x3f3800 also pulls in `gesvd`, T39).
- [ ] **T121** ReLAPACK follow-ups (split from T118, 2026-10-30): a blocked `gbtrf` from LAPACK's
  `dgbtrf` (rstd runs the unblocked `xGBTF2`, since ReLAPACK's recursion is wrong); race faer's
  Bunch-Kaufman (different storage) against `sytrf`. `tgsyl`'s `xLATDF` estimator is a local loop
  in `sylvester.zig` copying the update step of `gecon`'s private `Gecon` (it needs the last
  vector, which `gecon` doesn't return), not `norm.zig`'s `lacn2` over `latrs` as the docs say:
  share it or fix the docs; `tgsyl` not timed yet. Sweep: `relapack` on every rig, and with
  `-Drelapack=true` on one or two.
- [ ] **T122** Small-matrix follow-ups (split from T119, 2026-10-30): the per-core small kernels for
  SVE, POWER10 and LASX (tier 3 or untestable here); the small-path thresholds are upstream, tuned
  for its kernels, so time them under T05; `sbgemm` small kernels off x86_64; from T119's port:
  `nn`'s dot rows are read in place, not through upstream row-major copy (same bits, likely
  slower on AVX-512 hardware), and Cooper Lake's `sbgemm` small path keeps the portable bits, not
  `vdpbf16ps`'s pair rounding (by design, documented). Sweep: `3m`, `small`, `interface: gemm` on
  every rig, `aarch64-sve`, `powerpc64le-power10` and loongarch64 first (their permits allow the
  small path), `blat3_3m` on `x86_64` and `i386-max`.
- [ ] **T47** Hardware-only paths never run: the mips `cpucfg` decode, the sparc `sllx` rewrite,
  Linux ppc64 `mfpvr`, the Windows arm64 registry MIDR and feature probes, the arm64 sysfs fallback,
  the FreeBSD arm64 probe, Accelerate's pre-13.3 symbol fallbacks (`$NEWLAPACK`, classic `cblas_*`);
  NEON and AltiVec flushing seen on rigs only. **Blocked:** hardware. Accelerate's pre-13.3 symbol
  fallbacks are tier 2 (zfls supports macOS 11.0): they need a macOS 11 or 12 machines or VM to run.
  Tier-2 parts done 2026-10-30 (landed b1cad78a's copy): the Windows arm64 registry and feature
  probes, the arm64 sysfs fallback and Accelerate's pre-13.3 fallbacks are tested with recorded
  inputs (Raspberry Pi 3 to 5, Neoverse N1 and N2, Apple M1 under Asahi, Snapdragon X Elite,
  Microsoft SQ1/SQ2, Azure Cobalt 100; macOS 11.3, 12.3 and 13.3 SDK export lists); no bugs found.
  Still open: Linux ppc64 `mfpvr`, the mips `cpucfg` decode, the sparc `sllx` rewrite, the FreeBSD
  arm64 probe, NEON and AltiVec flushing. Possible enhancements: read every Windows
  `CentralProcessor` key, map Cobalt 100 to `neoversen2` (as LLVM), alias Snapdragon X Elite to
  `neoversen1` (as OpenBLAS's build-time detection).
- [ ] **T22** C and intrinsics kernels not yet ported, levels 1 and 2, by T20's strategy, `@Vector`
  first as arm64 did (`vector.zig`); the coverage.md gaps, the largest first. Level 1: riscv64
  outside `axpy`; POWER8 to 10 `asum`, `scal`, `rot`, `swap`, `cdot`; MIPS MSA and zarch; x86_64
  `asum` on Haswell and AVX-512; wasm SIMD128. Level 2 (464 gaps, `blas2/coverage.md`): riscv64 140
  (gemv 82, hemv 32, symv 24, `sbgemv` and `shgemv` on zvl256b); mips64 100 (gemv on I6400, I6500,
  P6600, Loongson 3R3, 3R4); arm64 60 (symv on 11 cores, `sger`/`dger` on A64FX and Neoverse V1,
  `sbgemv`/`bgemv` on N2, V1, V2); POWER 60 (gemv on POWER8 to 10, `sbgemv` on POWER10); x86_64 50
  (symv on 10 cores, `sgemv` on the Skylake-X family (T24), `sger`/`dger` on Sandy Bridge,
  `sbgemv_n` on Cooper Lake and Sapphire Rapids); zarch 30 (gemv on z13, z14); mips 20 (gemv on
  P5600); wasm 4 (gemv). **Needs:** T20. riscv64 landed 2026-10-30 (f7b6be16): `@Vector` kernels
  fill every slot on `riscv64_zvl128b`, `riscv64_zvl256b` and `x280` (`vector.Rule.unfilled_on`);
  level-2 gaps 140 to 32, all C910V, which is **Blocked** (no `xtheadvector` in LLVM or Zig). x86_64
  landed 2026-10-30 (c1f9f988): `symv` on nehalem, sandybridge, haswell/zen, the Skylake-X family
  and Bulldozer family, `asum` on Haswell and AVX-512, `sger`/`dger` on Sandy Bridge, `sbgemv_n` on
  Cooper Lake and Sapphire Rapids (portable bits): x86_64 has no level-2 or `asum` gaps left; the
  AVX-512 and FMA4 kernels' `@Vector` code is tested natively, not on their hardware. arm64 landed
  2026-10-30: the NEON `symv_[LU]_asimd_4x4.c` port on Neoverse N1, Ampere1, `vortex`, `vortexm4`
  (upstream's upper driver drops rows when the start column isn't a multiple of 4; the port keeps
  them, issues.md), and `sbgemv`/`bgemv` on Neoverse N2, V1, V2 through `vector.Half`. Left after
  the MVP (tier 3, or tier-2 cores no rig runs): arm64's SVE `symv` drivers (A64FX, ARMv8 SVE, ARMv9
  SME, Cortex-A510, A710, X2, Neoverse V1) and `ger_sve_v1x3.c` (A64FX, V1); riscv64's C910V
  (**Blocked:** no `xtheadvector`); POWER, MIPS, zarch and wasm.
- [ ] **T73** Benchmarks above level 1 race `gemv` and `gemm` only: add `ger`, `symv`, `trsv`, a
  banded routine, `symm`, `syrk`, `trmm` and `trsm`, which T05 and the speed claims need. Rerun the
  Zen 3 flyoff (Windows and Linux), which predates T43's x86_64 `rot`, `scal`, `gemv` paths and the
  level-3 threading, under a timing rule for that host. **Needs:** T06 for threaded rows. Reference
  binaries beside OpenBLAS (design.md, the anthology decision): reference LAPACK, and Eigen (MPL
  2.0) built from upstream through Nix as a benchmark-only binary, never ported, for the LAPACK and
  level-3 races. Code landed 2026-10-30 (765b30fa..37d4b2ba, `bench-check` 240/240 on every tier 1
  and 2 host): benches for `ger`, `symv`, `trsv`, `gbmv` (kl = ku = 32), `symm`, `syrk`, `trmm`,
  `trsm`, a `lapack` section as table entries with hooks for T34, T110, T111; reference LAPACK
  3.12.1, Eigen 3.4.1 (our shim, built by Nix on the Mac) and AOCL-libFLAME 5.3.2 (x86_64 Linux)
  through Nix; a separate-process Accelerate column; `flyoff.sh zen3` and `flyoff.sh m6` under the
  R4 rule. Left: the first real flyoff runs and Nix builds (untested end to end); confirm the
  Windows idle rule it improvised (processor queue length under 2 for load); `sygst`/`hegst`,
  `trsyl`, `hetrf` entries; rstd calls for `geqrf`/`syevd` once T111 lands; a `-Drelapack` switch.
  Its tiny busy-machine run saw `trmm` at n = 64 at 0.09 to 0.40 of OpenBLAS: look at it in the
  timing session.
- [ ] **T94** Fact-check items still open after T78: the sparc64 `-Oz` miscompile (needs a reduced
  case on the sparc64 rig); the OpenBIOS 0x400000 to 0x800000 window (our observation; check
  OpenBIOS's source); 32-bit x86 Windows threads starting at x87 0x27f (T03's observation; no
  source, no test records it: add a test); "MKL leaves FTZ to the application" (from memory); Tegra
  2 as VFPv3-D16 (a forum post only); `symmetric.zig`'s "4 elements of slack past each".
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
  `FLA_Eig_gest`/`FLA_Sylv` against ReLAPACK's `sygst`/`trsyl`. From T101 (threading.md): OpenBLAS
  threads x86_64 `ddot`/`zdot` and arm64's `dot.c` above 10000 (arm64 with its own neoversev1 and
  a64fx tables) and never `sdot`, `cdot`, `dsdot` or `dot.S`; rstd threads every precision above
  10000. Threaded `cgemm` at n = 128 in the arm64 VM varies 198 to 249 GFLOP/s: each call (about
  75 µs of work) wakes five sleeping `Io.Threaded` pool threads (a task allocation, a mutex and a
  condvar signal each), whose wake time varies under HVF; no spinning server (threading.md), so a
  threshold that keeps such calls serial is the fix to time.
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
  `isInf`), no fiber `Io` tests on SVE targets (std's fiber switch clobbers `ffr`, T40). Also, wasm
  bitcode getting its feature flags twice (`asmCheck`'s
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

- [X] **T99** Multithreaded freestanding tests (from T02, see threading.md): start secondary cores
  in the test kernels with their own stacks, give each core thread-local storage (std's `debug.zig`
  makes `panic_stage` thread-local once a build is not `single_threaded`), an `Io` that hands group
  tasks to those cores with futex waits, wakes and cancellation, and start the rigs with more than
  one CPU. Done 2026-10-30, landed 7cab01a9: rigs `x86_64-smp` (fw_cfg core count, real-mode stub,
  INIT-SIPI-SIPI) and `aarch64-smp` (PSCI `CPU_ON`) start four cores with their own stacks and
  thread-local storage, and `test/qemu/smp.zig`'s `Pool` `Io` hands futures and group tasks to idle
  cores (waits in `mwait`/`wfe`); `threading.TestIo` puts the threaded tests on it. Smoke runs:
  x86_64-smp 30 passed, aarch64-smp 27 passed, 0 failed, 4 of 4 cores up. Left for Known issues: the
  canceled-gemm test is hosted-only, and timed waits return after one check (no clock on the rigs).
- [X] **T24** Skylake-X family `sgemv_n`/`sgemv_t`: those cores add large intrinsics kernels on top
  of Haswell's (`sgemv_t` uses a 1,121-line template); they keep the portable `sgemv` until ported.
  Testable only on AVX-512 hardware (T43). Done 2026-10-30, landed c1f9f988: Skylake-X `sgemv_n`
  (n_128 and n_64 blocks) and `sgemv_t` (the template's 8 shapes) as `@Vector` around Haswell's
  blocks in `gemv/kernel.zig`; their exact code runs natively on Zen 3 (LLVM splits the wide
  vectors), AVX-512 hardware still untested (T43).
- [X] **T116** The generic family follows `parameter.c:451`, an x86-only rule: a static
  `TARGET=GENERIC` OpenBLAS elsewhere gets `p = q = r = 128` (`parameter.c:46-48`), or fixed sets on
  MIPS and LoongArch (`param.h:4386-4415`). Make rstd's generic `p`/`q`/`r` follow the non-x86 rule
  off x86 (from T112). Running (`wip/quad`). Done 2026-10-30, landed d253ea5f: off x86 the generic
  table takes `parameter.c:46-48`'s `p = q = r = 128`, or 32-bit MIPS's and LoongArch's fixed sets
  (`param.h:4388-4415`); x86 keeps its L2 rule (`parameter.c:451`); the generator
  (`tools/blas_asm/params.zig`) emits `generic_off_x86`, every other table unchanged. Native 65/65,
  `qemu-check-powerpc-750` passed. Sweep: `params` and `blas3` on `arm-softfp` and `powerpc-750`.
- [X] **T111** The rest of LAPACK that would otherwise come from reference LAPACK, ported from faer
  (Q10): QR (`geqrf`, `orgqr`/`ungqr`, `ormqr`), least squares (`gels`), the symmetric and Hermitian
  eigensolvers (`syevd`/`heevd`), the general eigensolver (`geev`), `laed3` and OpenBLAS's
  `test_laed3_nan.c` (T70 left them: they need reference LAPACK's `dlaed4`, `dlaed5`, `dlaed6` and
  the test drives `dsyevd`; take faer's divide-and-conquer tridiagonal solver), and whatever else
  rstd.linalg's MATLAB-like API needs. Reference LAPACK comes later only as a benchmark baseline.
  Reuse T34's `lapack/*.zig` (`householder`, `qr`, `bidiag`, `ops`) and add files beside them;
  `qr.zig` already has faer's QR for `geqrf`. Running (`wip/faer2`). Done 2026-10-30, landed
  8ee93919, cd01f0a5: `geqrf`, `orgqr`/`ungqr`, `ormqr`/`unmqr`, `gels`, `syevd`/`heevd`,
  `syev`/`heev` and `geev` from faer v0.24.4 (`qr_solve.zig`, `tridiag*.zig`, `hessenberg.zig`,
  `real_schur.zig`, `complex_schur.zig`, `evd.zig`), `gebal`/`gebak` from reference LAPACK
  (`balance.zig`; faer has none), OpenBLAS's `test_laed3_nan.c` against `syevd`. Five more faer bugs
  fixed (issues.md: a NaN column in staircase QR, complex `rotg`'s `r`, a 2x2 determinant guard that
  broke left eigenvectors at `SMLNUM`, more denormal reciprocals). Native 99 of 100, 1 skipped.
  Follow-ups: T117.
- [X] **T115** Routines only reference LAPACK has (T113, `lapack/sources.md`): `gecon`,
  `gebal`/`gebak`, `trexc`/`trsen`, `gelsy`, plus `lange` (no source chosen) and `gelqf` (only
  libflame's FLAME code besides the reference). The anthology rule picks the best open-source
  implementation, and here that is the reference: port these from reference LAPACK (its license is
  in the accepted BSD class) unless faer turns out to have an equivalent, and record each in the
  table. Running (`wip/refonly`). Done 2026-10-30, landed 82ce5804, 9a9d4f95: `lange`, `gecon` (with
  `lacn2`, `latrs`, `rscl`, Blue's `lassq`), `trexc`, `trsen` (Sylvester solves through ReLAPACK's
  `trsyl`) and `gelsy` from reference LAPACK; `gebal`/`gebak` folded into T111's `balance.zig`
  (every job); `gelqf`/`orglq`/`ormlq` as libflame's UT-transform method on faer's kernels (no
  libflame code); `geqp3` from faer's pivoted QR (`qr_pivot.zig`). Native 96 of 97, 1 skipped;
  `arm-softfp`, `powerpc-750` compile. Known issues: extended-type `nrm2` underflows near the
  threshold (the extreme `gelsy` test skips there); a faer pivoted-QR `1/0` on an all-zero matrix,
  guarded (issues.md). ReLAPACK's `sygst`, `hegst`, `pbtrf`, `trsyl`, `gbtrf` are not exported from
  `lapack/__root__.zig`.
- [X] **T19** `XDouble` as `f128`: `quad_precision` is hard-wired off; add a build option and a rig
  configuration for it. The edge tests and level 3's per-routine driver tests still cover no
  extended type off x86 (the interfaces do). Parked 2026-10-30 for the MVP (tier 3 work off x86):
  wip d2ea5703 on `wip/quad` (worktree removed) has `-Dquad-precision`, a `blas_options` module and
  an `aarch64-quad` rig; only `zig build -h` (0.18) and the native types test (31/31) ran; left:
  0.17 `-h`, wider native tests with quad on and off, `qemu-check-aarch64-quad`. Finishing the
  started work (user, 2026-10-30). Running (`wip/quad`). Done 2026-10-30, landed: `-Dquad-precision`
  makes `XDouble` `f128` through a `blas_options` module, with an `aarch64-quad` rig; native edge
  and level-3 driver tests pass with quad off (72/75) and on (73/75), `qemu-check-aarch64-quad` and
  `arm-softfp` with quad on 64/64, smoke run 23/23 on the Mac. Sweep: `qemu-test-aarch64-quad` in
  full (tier 3).
- [X] **T119** AVX-512 small-matrix kernels (user, 2026-10-30, tier 1): OpenBLAS's Skylake-X
  `sgemm`/`dgemm` small kernels (`nn`/`nt`/`tn`/`tt` and their `b0` forms) and Cooper Lake and
  Sapphire Rapids `sbgemm` small kernels as `@Vector`, tested natively on Zen 3 (no AVX-512 here).
  Done 2026-10-30, landed cb5d5a00, 58f9ae29: `blas3/gemm/x86_64/kernel.zig` ports the Skylake-X
  small kernels (512-bit vectors; the generic kernels' bits except `tn` and `nn`'s last rows, whose
  lanes sum every 16th or 8th product and reduce in upstream's two orders) for Skylake-X, Cooper
  Lake and Sapphire Rapids through a `small_kernel.forCore` hook; Cooper Lake's `sbgemm` small path
  with upstream's permit and portable bits (Sapphire Rapids has none, as upstream); `survey.zig`'s
  small-matrix table in `gemm/coverage.md`. Native 43/43 (small-kernel filters), `-Dblas-cores=all`
  134/136 (2 skipped); `qemu-check-x86_64`, `-arm-softfp`, `-powerpc-750` 71/71 on the M6. The
  AVX-512 code runs only as Zen 3 splits it (T43). Follow-ups: T122. Sweep: `small`, `half gemm`
  and `interface: gemm` on `x86_64` and `linux-x86_64`; `small permits` everywhere.
- [X] **T101** rstd's `sdot` from 256K floats runs 51 to 53 GB/s on arm64 where its kernel alone
  runs 98 (R4, results-flyoff.md). Done 2026-10-30, landed f2b453bc: above n = 10000 a threaded
  `sdot` summed its chunks with `dsdot_k`, which was arm64's `dot.S` as DSDOT (one dependent `fadd`
  chain, about 50 GB/s) and x86_64's `sdot.c` as DSDOT (blocks multiplied and summed in `f32`,
  against the `f64` rule readme.md and threading.md state: a correctness bug too). `dsdot_k` is
  now routed as `sdot_k` (`vector.routed`): on every core exactly `sdot_k`'s `f64` sum before
  rounding; bits identical for any worker count; a per-core test checks `sdot_k` against `dsdot_k`
  rounded. The arm64, loongarch64, riscv64 and x86_64 `dsdot` ports are gone (regenerated). Native
  `-Dblas-cores=all` 129/132 (3 skipped); no rig compiled it yet, and the speed is unmeasured
  (timing session, `flyoff.sh m6` `dot s` rows). Sweep: `qemu-check-aarch64`, `-arm-softfp`,
  `-powerpc-750`; filters `dot` and `level 1` on `aarch64`, `linux-aarch64`, `arm-softfp` and
  `powerpc-750`. The `cgemm` n = 128 variance went to T05.
- [X] **T40** Linux rig gaps: `linux-powerpc64le` never ran POWER10 kernels (baseline build),
  `linux-aarch64` left out the SVE cores, and neither Linux arm64 rig ran `vortexm4` nor `armv9sme`
  (Debian's 6.12 netboot kernel builds without `CONFIG_ARM64_SME`). Done 2026-10-30, landed
  d6aadaf3, 0d257b04, 408246a2, d44510df: rigs `linux-aarch64-sve` (TCG, `neoverse_v1` build, no
  fiber `Io` tests on SVE: std's fiber switch clobbers `ffr`, issues.md) and
  `linux-powerpc64le-power10` (`pwr10`, ReleaseSafe), `qemu-check-<rig>` for the Linux rigs,
  Debian 14's 7.0.13 arm64 netboot kernel (built with SME) for every Linux arm64 rig, and a check
  that each rig selects its core (`the core in use` prints it; `core`/`hvf_core` in
  `qemutest.zig`). Smoke runs on the M6 (`linalg.blas.cpu.`): `power10`, `armv9sme` (TCG) and
  `vortexm4` (HVF) selected, 119 passed, 0 failed each. Sweep: tier 1 compiles the merged
  `threading.zig` `Evented` line natively (built only for SVE so far); `linux-aarch64` (tier 2) in
  full on the new kernel.
- [X] **T117** `geev` and `gesvd` correctness (from T111 and T34): `xTREVC3`-style overflow scaling
  in `geev`'s eigenvector back substitution, what `geev` returns after `NoConvergence`, and
  `gesvd`'s `INFO > 0`. Done 2026-10-30, landed 64f5de92, cf5f88e2: `gesvd` and `geev` take an
  optional last `unconverged: ?*usize` (LAPACK's `INFO > 0`: `gesvd`'s unconverged
  superdiagonals, `min(m, n) - 1` for an infinity or NaN; `geev`'s `i`, with `w[i..n]` the
  converged eigenvalues and no vectors, as `dgeev.f`), reached in tests through `ops.zig`'s
  `test_iter_limit`; `trevc.zig` ports `dtrevc3.f`, `ztrevc3.f` and `dlaln2.f` (overflow test:
  a 20 x 20 triangular matrix whose unscaled vectors pass the largest float); `ztrevc3`'s left
  solve passes `ZLATRS` the wrong columns' norms (issues.md; rstd passes the right ones). Native
  `lapack`, `faer`, `relapack` and `utest: kernel_regress` 158 of 160 (2 skipped);
  `qemu-check-x86_64`, `-arm-softfp`, `-powerpc-750` 74 of 74 on the M6. Follow-ups: T120. Sweep:
  `faer svd`, `faer evd` on every tier-1/2 rig.
- [X] **T118** ReLAPACK's `tgsyl` and the missing exports. Done 2026-10-30, landed 8cd41249,
  b79d8486: `tgsyl` by ReLAPACK's recursion over reference LAPACK's `xTGSY2`, `xGETC2`, `xGESC2`
  and `xLATDF`, every precision, `dif` for `ijob` 0 to 4; two ReLAPACK bugs fixed (issues.md: the
  transposed updates' leading dimension for `F` and `E`; real `dif` for `ijob` 2 and 4 counting
  only the last base case); exported from `lapack/__root__.zig`: `sygst`, `hegst`, `pbtrf`,
  `gbtrf`, `trsyl`, `tgsyl`, `sytrf`, `hetrf`, `sytrfRook`, `hetrfRook`, each with a smoke test.
  Checks as T117's. Follow-ups: T121. Sweep: `relapack tgsyl`, `relapack trsyl`, `lapack:
  ReLAPACK` on every tier-1/2 rig, `x86_64`, `aarch64`, `arm-softfp` and `powerpc-750` first
  (extended types).
- [X] **S1** Audit what the SPIR-V backend emits for shaderlib's math, operation by operation.
  Done 2026-10-30, landed 27981ef4, 3352bb51: `src/rgl/shaders/probe.zig` is one fragment shader
  per operation (`zig build shader-probes -Dvulkan`, one probe with `-Dshader-probe=<op>`; the
  build reads its `Op` enum), disassembled and `spirv-val`-clean in Debug and ReleaseFast; the
  results table is in `src/rgl/shaders/lib/readme.md`. Only the math builtins, `@min`/`@max` and
  `@mulAdd` lower natively; `@Vector` arithmetic stays vector; `dot`, `length`, `distance`,
  `normalize`, `cross`, `reflect`, `refract`, `faceforward`, `inversesqrt`, rstd's `mulAdd` and
  `mix`, `step`, `smoothstep`, `fract` and `clamp` are scalar or hand-built; inline `asm` reaches
  `OpDot` and GLSL.std.450 (`dotAsm`, `crossAsm`); inline asm cannot declare `OpTypeMatrix`
  (confirmed by compiling one) and `comptimePrint` as an `asm` template crashes the backend (both
  issues.md); every inlined call keeps its structured blocks even in ReleaseFast. Fixed on the
  way: shaderlib's `inversesqrt`, `reflect` and `refract` did not compile for vectors. Follow-up:
  S2 (the probes are its fixture).
- [X] **T124** `lapack.geev` on `Complex(f128)` didn't compile (found by M1): faer's complex Schur
  (`eig22`) called `std.math.complex.sqrt`, which has only `f32` and `f64`. Done 2026-10-30, landed
  ff132aa4: `ops.zig`'s `sqrt` keeps std's for `f32`/`f64` (same bits) and runs musl's `csqrt`
  over any other float type; tests: `geev` in `XDouble`, unslotted `f80` and `f128`, real and
  complex, and the square root against `f64`'s; the std gap in issues.md. Native `faer evd` 42/42;
  `qemu-check-arm-softfp`, `-powerpc-750` 24/24 on the M6. Sweep: `faer evd` on every tier-1/2 rig.
- [X] **M1** Audit `rstd.linalg.matlab` (the matlab88 port) against the modern BLAS and LAPACK. Done
  2026-10-30, landed (the audit, now in design.md's "The matlab port
  against the modern BLAS and LAPACK"; `lu.zig`'s test): every routine's modern
  counterpart and verdict, checked by experiments E1 to E9 against an `f128` LAPACK reference. The
  modern routine wins or ties everywhere (balancing helps, E3; `gecon`'s rcond is tighter; `nrm2`
  at underflow; Cholesky's NaN and complex diagonals) except E1 (`wgedi`'s scaled determinant,
  which M2 keeps over `getrf`; the test), E2 (graded Hermitian, T125) and E4 (graded bidiagonal,
  T126), which the modern LAPACK fixes (user, 2026-10-30). M2's layout: convert at the boundary.
  The matlab tests run natively and in the macOS sweeps, not on the rigs (`_test_qemu.zig` takes
  the BLAS and LAPACK only).
- [X] **T127** Flatten `src/linalg/lapack/faer/` into `src/linalg/lapack/` (user, 2026-10-30): the
  folder sorted by source, which no longer held, and provenance lives in each file's doc comment
  and `sources.md`. Done 2026-10-30, landed 9f823ac3: the 16 files moved, imports and paths fixed
  (upstream faer's `faer/src/...` paths and the `faer ...` test names kept). Native LAPACK,
  faer, ReLAPACK and `wgedi` tests 161/163 (2 skipped); `qemu-check-x86_64`, `-arm-softfp`,
  `-powerpc-750` and `bench-check` 319/319 on the M6. Merging related files is Z6's.
- [X] **M5** Review M3's proposal with the user. Done 2026-10-30: all 11 questions settled
  (storage packed column-major; `StaticMatrix(R, C, T)` and unmanaged `Matrix(T)` over one
  operation set on views; inline up to 8 and the BLAS above, by size alone; a root declaration
  for the limit; `norm(v)` with `error.ZeroVector`; `lu` factors with `zero_pivot`; `eig` adapting
  to its input; `linalg.blas`/`lapack` public; complex fixed sizes now; `error.Singular` in
  shaders), folded with M1's audit into design.md ("Planned: the public `rstd.linalg` layer", "The
  matlab port against the modern BLAS and LAPACK"); `design-linalg.md` and `matlab/audit.md`
  dropped (user, 2026-10-31). Unstaged for the user's review.
- [X] **T128** Move `design.md` and `todo.md` from `src/linalg/blas/` to `src/linalg/` (user,
  2026-10-31): they cover LAPACK, the public layer and the matlab port, not only the BLAS. Done
  2026-10-31, unstaged for the user's review: both moved and retitled ("Linear algebra design",
  "Linear algebra to-do"), every relative Markdown link they broke rewritten (in them and in 8
  other files; a check finds none broken), `lu.zig`'s mention, `.agent`'s tooling and the memory
  notes updated. design.md's backticked paths keep their stated `src/linalg/blas/` base.
- [X] **T129** Move `issues.md` from `src/linalg/blas/` to `src/linalg/` (user, 2026-10-31). Done
  2026-10-31, unstaged for the user's review: moved, its intro widened to `rstd.linalg`, its links
  and the links into it rewritten (none broken), the plain mentions in `relapack.zig`,
  `sources.md`, `tools/evented_group_repro.zig`, `build/sshtest.zig`, `.agent` and the memory notes
  updated. `blas/readme.md` stays (the BLAS's user guide); M3 adds `src/linalg/readme.md`.
- [X] **T123** LAPACK routines the matlab port (M2) needs. Done 2026-10-31: pieces 1 and 2 landed
  4b5d9c4a (`getri` from faer's LU inverse; `laset`; `lapy2`; `gehrd` and `orghr`/`unghr` over
  faer's Hessenberg reduction with LAPACK's `tau` and `ilo`/`ihi`); pieces 3 to 5 staged by the
  user for review (`hseqr` over the `lahqr`s; `trevc3`/`trevc`, every `SIDE` and `HOWMNY`;
  `gees` as `dgees.f`/`zgees.f`, sorting by `trsen`), with `sources.md`, design.md and the BLAS
  readme. Native: every LAPACK test 166/168 (2 skipped); `qemu-check-x86_64`, `-arm-softfp`,
  `-powerpc-750` 74/74 on the M6 after each piece. Sweep: `getri`, `laset`, `gehrd`, `hseqr`,
  `trevc3`, `gees` on every tier-1/2 rig.
- [X] **T125** Hermitian and symmetric eigenvalues of graded matrices (M1's E2). Done 2026-10-31,
  unstaged for the user's review (`evd.zig`): `.u` reduces `J A J` (the reversal) as its lower
  triangle, so faer's top-down reduction runs from the bottom, as `xSYTRD`; the vectors' rows are
  reversed after. Tested against `f128` both ways, `syevd` and `syev`, real and complex (without
  the fix `2.8e5` off); every LAPACK and matlab test 236/238 (2 skipped); rig compile checks 74/74.
- [X] **T130** `gesdd`'s singular values out of order on the `x86_64` rig. Done 2026-10-31,
  unstaged for the user's review (`bidiag_svd.zig`, issues.md): faer's `compute_svd_of_m` sorts a
  merge's values with one bubble pass, and heavy deflation left several out of place on
  `haswell`'s rounding (reproduced natively with `-Dcpu=haswell`); a stable insertion sort
  (`sortAscending`) with its own test. `haswell`'s SVD tests 39/39, `zen`'s every LAPACK test
  167/169, `qemu-test-x86_64`'s SVD tests on the rig 53/53, rig compile checks 74/74.
- [X] **T126** Singular values of graded bidiagonals (M1's E4). Done 2026-10-31, unstaged for the
  user's review (`bidiag_svd.zig`, `svd.zig`, issues.md, `sources.md`, design.md):
  `qrAlgorithm` is now `dbdsqr.f` (with `dlasv2.f`, `dlas2.f`), faer's scaling kept around it,
  the rotations applied as made, so `gesvd`, `gesdd`'s small sizes and the divide-and-conquer
  leaves all keep a graded bidiagonal's small values. E4's graded bidiagonal, reversed and random,
  to `256 n eps` against `f128` (faer's QR: 13% off); every LAPACK and matlab test 237/239 (2
  skipped); rig compile checks 74/74. Follow-up: T131 (dqds for values alone).

## After the BLAS is feature complete

- [ ] **M4** Time `linalg.vector`'s inline `@Vector` forms against the BLAS and LAPACK at sizes 4
  to 16 per operation and type (`dot`, `length`, `mul` matrix by vector and by matrix, `solve`), on
  the M6 and Zen 3 under their timing rules: confirm the inline limit of 8 (design.md, the public
  layer's Dispatch) or move it. **Needs:** M3's step 4.

## SPIR-V: vector math in shaders

Not the BLAS kernels: the fixed-size vector math (`rstd.linalg.vector`, `rstd.simd` until M3) that
`src/rgl/shaders/lib/` (`math.zig`, `sdf.zig`) builds on. Can start any time. SPIR-V's own
instructions are used as they are, not reimplemented or tuned: the aim is one API, not speed.

- [ ] **S2** Shaderlib stays the shaders' one API (design.md, Settled decisions): GLSL-friendly
  names that alias or wrap `linalg` calls where `linalg` has the operation (`dot`, `length`,
  matrix-vector products, through M3's SPIR-V forms); `cross` moves to `rstd.linalg.vector` and is
  re-exported; `normalize`, `reflect`, `refract`, `faceforward`, `smoothstep`, `fract`, `step` stay
  in shaderlib, lowered to GLSL.std.450 instructions. Add a build-time test that disassembles the
  shaders and fails if a listed operation comes out as scalar code: `zig build shader-probes`
  (S1's `probe.zig`) with a check over `spirv-dis` output, against `shaders/lib/readme.md`'s table.
  **Needs:** M3 (for the `linalg`-backed ones); S1 done.

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
  names, documented panics), and the file layout: fewer, bigger files where pieces belong together
  (`lapack/`'s `schur.zig`, `real_schur.zig` and `complex_schur.zig`, after T127 flattened it).
  **Needs:** Z1 to Z5.

## Waiting on the user

Nothing waiting.
