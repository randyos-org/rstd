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

T62's tier 2 failures, fixed one at a time on `review` as the sweep finds them (In the background).

### In the background

Sweeps by `.agent/sweep.sh` (gitignored), results in `.agent/sweeps/<UTC>-<tier>/summary.md`;
started 2026-10-30 22:20 on `review` faeea2aa: tier 1, and tier 2 then tier 3. These two ran from
the trunk, so steps after their first compile what has landed since (T124, M1, kern3, T127); from
tier 3 on, each tier runs from its own worktree, `C:\scratch\git\rstd-wt\sweep-<tier>`, pinned
at review's tip when it starts (user, 2026-10-30). M2's trunk edits moved tier 2 there too: it
resumed `--from x86_64-smp` at 1650346b (2026-10-31 07:23), after steps 3 to 7 had compiled the
trunk with M2 half done. With the fixes landed (2d48532d, 4962dbf2, a54e0240, ec15b11c), a fresh full tier 2
started at da2ae06a (2026-10-31 14:05), then tier 3. Split 14:44 (user): the main sweep stops
after `x86_64`; `--tree t2b` (the other x86_64 rigs, `T08-x86_64`), `t2c` (aarch64 rigs,
`T08-aarch64`), `t2d` (arm and PowerPC, their T08 steps) run on the Mac side by side, `t2e`
(`x86_64-macos`, T114, the compile checks) here, each in its own worktree and `-<tree>` folder.
The Accelerate fixes landed as b18d7dae; `t2f` runs `ssh-x86_64-macos` and `ssh-aarch64-macos`
there (14:46).
Stopped 15:40 (user: fail fast, fix, retest only the failure, then resume): `sweep.sh` now passes
`-Dfail-fast` and stops at a failing step (`--keep-going` for the old way), with
`sweeps/STOP-<commit>` stopping the other sweeps of that commit. `aarch64` had passed (645, 51
skipped). `x86_64` failed `trevc3` (selected vectors differ in the last bits from all of them):
on `prescott` and 10 other SSE-era cores, the converted complex `axpy`/`axpyc` assembly (`c`, `z`)
gives other bits when `y` is misaligned (8 mod 16, even at `n` 1), breaking threading.md's "no
address dependence". The address test now covers the element-wise kernels (`edge.zig`) and finds
them. Upstream's own: its aligned-`y` paths add `y` before `alpha_i` times the swapped `x`, the
misaligned ones after `alpha * x`, as do three of `zaxpy_sse.S`'s one-element tails. Fixed by
generator patches (`yLast` and three tail patches, `arch/x86_64.zig`) that add `y` last on every
path, keeping the kernels; the address test (complex `scal` with a zero part too) and `trevc3`
pass natively on every core and on the `x86_64` rig (user: no dropping OpenBLAS kernels without a
very high bar; T133 to T136 evaluate the earlier such choices). Local builds
use the user's experimental Zig 0.18.0-dev.856: a compiler problem there is reported, not worked
around.

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
  skipped, 0 failed; no inline-asm errors. Tier 2 so far
  (`.agent/sweeps/20261031T032047Z-tier2`), to fix: `x86_64` 1 failed (faer SVD special matrices,
  T130, compiled before it landed); `x86_64-uefi` 16 failed (BLAS gemm tests out of the 32 MB heap,
  the linked-core triangular test crashed or hung); `x86_64-smp` 2 failed (threaded level 2 out of
  heap; threaded level 1 bits differ by worker count); `linux-x86_64` a compiler crash (`unreachable`)
  generating `bidiag_svd.lasv2` for `prescott` (T126's code); `ssh-x86_64-macos` a failed
  `expectEqualSlices`, ssh exit 255; `x86_64-macos` names `qemu-test-x86_64-macos`, which no longer
  exists (only setup, install, start, snapshot); `aarch64` passed (644, 51 skipped). Steps 3 to 7
  rerun in the resume; `x86_64` and `x86_64-uefi` need `--only` reruns after. Paused 2026-10-31
  ~07:50 for a reboot, the resume having redone `x86_64-smp` (the same 2 failed) and
  `linux-x86_64` (the same crash). Plan: fix the known failures on top of 7511eef3, then a fresh
  full tier 2 at the fixed tip, then tier 3. Fixes so far (unstaged): `linux-x86_64`'s crash is
  the self-hosted x86_64 backend on `lasv2`'s f80 instance, its `if` expression now a statement
  (issues.md); `x86_64-smp`'s level 1 compared an f80's undefined padding bytes, now value bits
  (`check.appendBits`); its level 2 rank updates kept 16 MB of bytes twice, now digests
  (`check.appendDigest`, `ger` too); `x86_64-uefi` ran the native-only cases (biggest sizes, other
  targets' tiles) because `check.rig()` missed `.uefi`, hence its out-of-memory failures; now a
  rig, so the tile test skips; its AVX cores' gemm tests still need 24 MB plus a 20 MB packing
  buffer, so the x86_64 rigs' heap is 64 MB (T132); its triangular test passed on rerun; `x86_64-macos`:
  `zb.sh` no longer pre-builds a `qemu-check` it lacks, and `sweep.sh` links the trunk's
  `.macos-vm` into the worktree. `x86_64-smp`'s threaded tests then passed (23, 0 failed; level 2 reads elements back through
  `check.fromBits`; heap high water 19.8 of 64 MB). `ssh-x86_64-macos` rerun: 906 of 997 passed, then
  utest's `axpby` at stride 0 sent `incy` 0 to Accelerate's `catlas_saxpby`, which rejects it and
  exits (255); `apple_accel.axpby` now keeps `incy` 0 on rstd (ec15b11c). Its check (utest and
  `apple_accel`, both Macs) found two more: x86_64's `cblas_sswap` rejects a zero `incx` the same
  way, so every zero stride on a written vector now stays on rstd (`copy`, `swap`, `scal`, `rot`,
  `rotm`, `axpy`, `axpby`); arm64's complex `gemv` with `beta` `(2, 0)` and an infinite `y` skips
  the zero part's products (`(inf, 1)` gives no NaN), so such calls now stay on rstd too. The
  recheck: arm64 78 passed; x86_64 failed utest_ext's complex `gemm` with `R` (2.46e-4 over 1e-4),
  its two calls on different backends (Accelerate has no `R`), now both on rstd; then 78 of 78
  on both (b18d7dae). The tier 2 sweep, pinned at da2ae06a, lacks these: its
  `ssh-*-macos` steps need `--only` reruns once they land.

Tier 2 at b165b084 (2026-10-31 evening), so far: `arm` passed (643, 53 skipped). Failures:
`ssh-x86_64-macos` (Rosetta, Accelerate) "trsv and trsm give the same bits": `trsm` gave NaN for
every element but an infinity, on a diagonal system with subnormal entries and NaN in the unused
triangle; then threaded level 2 went 60 minutes without output (the step counted as BROKEN).
Accelerate's `trsm` multiplies by the diagonal's reciprocal, which overflows for a subnormal
(the reference divides; the unused triangle and subnormals in `B` are fine): a diagonal entry too
small to invert now stays on rstd, and the test runs both on rstd. Checked on the M6
(`ZB_WHERE=mac`, new: builds and runs only there). `aarch64-smp` ran out of its 32 MB heap in
threaded gemm (35 MB): aarch64's rig heap is 64 MB, as x86_64's. `x86_64-uefi` (2 h 37 min) ran out of its 64
MB heap in gemm3m's per-core test at `sandybridge` (41.9 MB wanted, 48.2 MB in use; the
`x86_64` rig's peak is 23 MB): `blas3/check.zig` crossed `R` with one row only where
`os.tag == .freestanding`, so UEFI took `unroll_m + 1`; it now takes one too (23 MB peak on the rig);
`t2b`'s other steps (`x86_64-smp`, `linux-x86_64`, `T08-x86_64`) run meanwhile, at 74305c34.
`x86_64` passed (654, 42 skipped, 2 h 58 min); `aarch64-smp` passed; `T114-aarch64-macos` passed
(46, 1 skipped). `ssh-x86_64-macos` (986 passed, 1 failed): T114's own test, Accelerate's
x86_64 `cpotrf_` and `zpotrf_` under Rosetta call a well-conditioned `L Lᴴ` not positive
definite from order 2 (real and arm64 right); complex `potrf` on x86_64 now stays on rstd
(46 of 47, 1 skipped, under Rosetta and on arm64). The native runner's `FAIL (...)` line does not end an SSH step under
`-Dfail-fast` (it ran on to the end). `linux-aarch64` passed (965, 45 skipped), `T08-aarch64` passed (149, 11
skipped). `linux-x86_64` could not unpack its initramfs in 1 GB ("write error", the test killed
by SIGBUS): its test is 711 MB (every x86_64 core's kernels with debug info; 146 MB for an arm
one), so the rig has 2 GB (8ab4a028). Then it booted, but `rgl.image`'s `writeRawFile` test got
NoSpaceLeft writing into `/` (the unpacked initramfs, a tmpfs of half the memory, the test 724 MB
of it): 3 GB (943c7810), rerunning. That test writes into the working directory, not
`TMPDIR` (`/tmp` on the rigs); rgl's to change.
`T114-aarch64-macos`: "potrf: a failed pivot is stored" kept the complex pivot's imaginary part,
because `accel.route(false)` came before the BLAS's init, which routes again, when the filter
made it the first BLAS call; it and `threaded.zig`'s `rstdOnly` now init first (unstaged).

The two Mac suites, compiled here with the experimental Zig 0.18.0-dev.858 (earlier runs were
compiled on the Mac's dev.751-era Zig), fail tests no rstd change touched: `gemm.kernel`'s
"geometry: the most accumulators the registers hold" (expected 1, found 4) and `_corespec`'s "no
kernel object needs a feature" (mips64 `i6400`, `i6500` objects now need `fpu_d fpu_s`), both on
x86_64 and arm64, and on arm64 `faer evd: geev` and `faer svd: gesdd` "one worker's bits at any
worker count". Built by the Mac's Zig the first two pass (x86_64 and arm64) and fail natively
with dev.858: the toolchain's MIPS feature data (`kernel.zig:1020` compares mips r6 with `p5600`;
`_corespec.zig:687`), reported to the user, rstd unchanged. `geev` and `gesdd` failed with the
Mac's Zig too: their tests run one worker unchunked and the others in many chunks, and
Accelerate's bits follow a part's shape, so they now run on rstd's kernels (the real contract,
the same chunks at any worker count, `apple_accel`'s split test covers with routing on).

#### Build slots

Compiling, one row per slot; a native test compiled on the Mac and run here shows in both tables.

| Slot | Item | Job | Since |
| --- | --- | --- | --- |
| Local 1 | T62 | `t2e`: `x86_64-macos`, `check-macho`, `check-uefi`, `check-zig017` (74305c34; Zig 0.18.0-dev.858) | 20:40 |
| Local 2 | T62 | `t2f`: `ssh-x86_64-macos`, `ssh-aarch64-macos`, `T114-x86_64-macos` (23f164c4) | 20:54 |
| Mac 0.18 | T62 | shared by the sweeps' compiles | 19:05 |
| Mac 0.17 | free | | |

#### Running tests

Test executables and QEMU rigs, outside the slots. From 21:10 each queued step has a sweep of
its own (`--tree t2g` to `t2n`, at 23f164c4); `t2b`, `t2d` and `t2f` end after their current step.

| Where | Item | Run | Since |
| --- | --- | --- | --- |
| Mac | T62 | `t2a`: `x86_64-uefi` (23f164c4) | 20:54 |
| Mac | T62 | `t2b`: `x86_64-smp` (74305c34) | 19:05 |
| Mac | T62 | `t2c`: `linux-aarch64` (74305c34; its sweep died on TERM, the step runs on) | 21:02 |
| Mac | T62 | `t2d`: `linux-arm` (b165b084; `arm`, `arm-softfp` passed) | 20:00 |
| Mac | T62 | `t2f`: `ssh-x86_64-macos` (23f164c4) | 20:54 |
| Mac | T62 | `t2g` `linux-x86_64`, `t2h` `T08-x86_64`, `t2i` `T08-aarch64`, `t2j` `powerpc-750`, `t2k` `powerpc-750-hf`, `t2l` `T08-arm-softfp` and `T08-powerpc-750`, `t2m` `ssh-aarch64-macos`, `t2n` `T114-x86_64-macos` | 21:10 |
| Local | T62 | `t2e`: the macOS VM, then the three compile checks (74305c34) | 20:40 |

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

- [ ] **M7** Redesign the fixed-size vectors and matrices around `@Matrix` (user, 2026-10-31),
  ahead of M3, whose types follow from it: `Vec(N, T)` over `@Vector`, `StaticMat(R, C, T)` over
  `@Matrix` where the compiler has it (this toolchain's randy-matrix), and a fallback with the
  same API where it doesn't (Zig 0.17 for zfls, stock 0.18), so code written for either compiles
  on both. Settle: the comptime test for `@Matrix`; storage (`@Matrix` pads columns to the vector
  size, so `@Matrix(3, 2, f32)` is two `@Vector(3, f32)`s, against M5's packed `[C][R]T`) and
  what views and the BLAS see; complex elements; SPIR-V's `OpTypeMatrix` through `@Matrix`
  (`tools/spirv_matrix_probe.zig`); the same bits from both forms; `Mat(T)` beside them. It also
  hosts the vector primitives ported kernels need (pffft's interleave, transpose and half-swap,
  rstd's N10). Revisits M5's decisions 1, 8 and 11 (design.md); M3, M6 and rstd's N1 follow it.
  **Needs:** nothing.
- [ ] **M3** Build `rstd.linalg`'s public layer, the BLAS and LAPACK behind it, as designed in
  design.md's "Planned: the public `rstd.linalg` layer" (settled in M5, 2026-10-30):
  `StaticMat(R, C, T)` (comptime shape, a value) and unmanaged `Mat(T)` (run-time shape),
  real and complex, one operation set over views, inline `@Vector` code up to 8 per dimension and
  the BLAS and LAPACK above by size alone, SPIR-V's own instructions in shaders; the `rstd.simd`
  split and its callers; its migration plan's seven steps, each a commit;
  `src/linalg/readme.md`, the layer's user guide, linking down to `blas/readme.md`. **Needs:** M7
  (M1, S1, M5 done); M4 confirms the inline limit after step 4.
- [ ] **M6** A type-agnostic matlab port, as the BLAS and LAPACK are: `CxMat` and its routines
  generic over the element type (`comptime T`, any float type), not hard-coded `f64`, with a
  caller-chosen result type where MATLAB fixes `double` (an explicit `comptime` type, the nearest
  Zig has to `@floatFromInt`'s inference). Pulled into the Queue by the user, 2026-10-31.
  On `@Matrix` where the compiler has it, as M7 designs. **Needs:** M7, then a decision with M3
  on whether the port keeps its own types or moves onto `Mat(T)`.
- [ ] **T138** Compensated summation (Neumaier) as an option on the BLAS reductions (`sum`,
  `asum`, `dot`), still bit-identical for any thread count; `rstd.stats` uses it through the BLAS.
  Benchmark it against the plain forms on the M6 Mac mini and Zen 3, then choose the defaults.
  Pulled into the Queue by the user, 2026-10-31.
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

- [ ] **T137** `tools/macos_vm.zig` copies the tests in as soon as port 2222 answers; on
  2026-10-31 (T62, `t2e`) the VM came up after 527 s and the first `scp` got "Connection reset by
  peer", so the step ran no test and stopped its sweep; the second try failed alike. `patient`
  retries only ssh's exit 255, `scp` exits 1; the copy now retries on any failure (unstaged,
  checking). Then (2026-11-01) the copy worked, but the run got "permission denied": Windows'
  `scp` gives the copy no execute bit; `chmod +x` before running it (unstaged, checking).
- [ ] **T136** Evaluate the x86_64 `peels` slots (`vector.peeling`: `asum`, `sum`, `ddot`, `cdot*`,
  `zdot*`, `dgemv_n`), whose converted assembly (`asum_sse*.S`, `zasum_*.S`, `zsum_sse*.S`,
  `dot_sse2.S`, `dot_atom.S`, `ddot_bulldozer.S`, `zdot_*.S`, `dgemv_n*.S`) rstd skips for
  `@Vector` kernels on the cores that build it, because its bits follow an operand's address.
  Never benchmarked on those cores. Deviating from OpenBLAS needs a very high bar (user,
  2026-10-31): confirm the dependence upstream, look for a generator patch that keeps the kernel
  (as `zaxpy_sse2.S`'s add order, T62), time both on each core (forced on this host, indicative),
  then decide with the user.
- [ ] **T135** Evaluate the routed extremes (`i?amax`, `i?amin`, `i?max`, `i?min` and values; the
  first NaN wins, `readme.md`, One contract): every architecture runs rstd's kernels, not
  upstream's. On `zen` `isamax` is 0.24 to 0.41 of OpenBLAS below 16M floats (flyoff,
  2026-10-27). Weigh the contract against upstream's NaN behavior and speed per architecture;
  patch upstream where it is wrong, if that keeps it; decide with the user.
- [ ] **T134** Evaluate routed `nrm2` (exact `f64` sum for `s`, `c`; Blue's for `d`, `z`): rstd's is
  faster on `zen` (1.2 to 11.7) and the M6 (4.3 to 17), unmeasured elsewhere. issues.md lists
  upstream kernels wrong on special values: confirm which are wrong rather than only different,
  patch those that can be, time the rest per architecture; decide with the user.
- [ ] **T133** Evaluate routed `sdot` and `dsdot` (one `f64` sum): short `sdot` is 0.30 to 0.65 of
  OpenBLAS on `zen`, 0.37 at 256K floats on the M6. Weigh the accuracy contract against
  upstream's `f32` sums per architecture; decide with the user.
- [ ] **T132** The gemm driver allocates the whole `BUFFER_SIZE` packing buffer (20 MB on the AVX
  x86_64 cores) per call, however small the problem; allocating only what the blocking uses for
  `m`, `n`, `k` would keep the bits and save the memory (found by T62: the x86_64 rigs' heap is 64
  MB for it).
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

Cleared 2026-10-31; landed work is in `review`'s log.

- [X] **M2** The matlab port over `rstd.linalg.blas` and `rstd.linalg.lapack` (7511eef3), reviewed
  by the user 2026-10-31; its rig tests run in T62's tier 2.

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
