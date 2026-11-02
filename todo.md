# rstd to-do

Open work in priority order, top to bottom. **Needs:** items that must land first; **Blocked:** an outside cause.

## In progress

### Working on

- [ ] **T154** "Skylake-X small kernels: upstream's blocks, lanes and reductions, on any x86_64" fails on the Mac under Rosetta with `ROSETTA_ADVERTISE_AVX=1` (`f32`, `nn`, m = n = k = 1; `small.zig`'s `checkSkylakex`) at 5f7ae56f, though it passed there at d3fa4af4 and passes natively on Zen 3. Bisect the commits between (ff7255f9's fusion setting first), then tell rstd's path apart from Rosetta's FMA rounding. Tier 2.

### In the background

- None.

## Queue

- [ ] **M10** The public layer's step 3 (src/linalg/design.md's "Planned: the public `rstd.linalg` layer", Migration plan): `View`, `ConstView` and unmanaged `Mat(T)` (run-time shape), real and complex; the dynamic functions over views, with the run-time switch into the inline kernels up to `L`. **Needs:** nothing (M9 done).

- [ ] **M11** The public layer's step 4: past `L`, the BLAS and LAPACK routines behind the `...Into` forms, chosen by size alone, and the inline limit's `Options` (`rstd_linalg_options.inline_limit`); both paths tested. **Needs:** M10.

- [ ] **M12** The public layer's step 5: its SPIR-V forms, each operation SPIR-V's own instruction where one exists. **Needs:** M11.

- [ ] **S2** Shaderlib stays the shaders' one API (src/linalg/design.md, Settled decisions): GLSL-friendly names that alias or wrap `linalg` calls where `linalg` has the operation (`dot`, `length`, matrix-vector products, through M12's SPIR-V forms); `cross` moves to `rstd.linalg.vector` and is re-exported; `normalize`, `reflect`, `refract`, `faceforward`, `smoothstep`, `fract`, `step` stay in shaderlib, lowered to GLSL.std.450 instructions. Add a build-time test that disassembles the shaders and fails if a listed operation comes out as scalar code: `zig build shader-probes` (S1's `probe.zig`) with a check over `spirv-dis` output, against `src/rgl/shaders/lib/readme.md`'s table. SPIR-V's own instructions are used as they are, not tuned: the aim is one API, not speed. **Needs:** M12 (for the `linalg`-backed ones); S1 done.

- [ ] **M13** The public layer's step 6: the float-slice reductions move to `linalg` (the NaN change). **Needs:** M12.

- [ ] **M16** pffft's vector primitives (interleave, transpose, half-swap) in `Vec` and `StaticMat` on `@Matrix`, as M7 settled (src/linalg/design.md, decisions 13 to 19), for N10. **Needs:** nothing (M9 done).

- [ ] **M15** `src/linalg/readme.md`: the public layer's user guide, linking down to `src/linalg/blas/readme.md`. **Needs:** M11.

- [ ] **Z1** Slices for every operand instead of pointer and size arguments; routines check that the slices cover what sizes and strides need. Z1 to Z6 are decided, not built (src/linalg/design.md, "Planned: a Zig-idiomatic interface"); before M6, so the matlab port moves once.

- [ ] **Z2** Exact error sets per routine or family, with plain names (`SliceTooShort`, `ZeroStride`, `LeadingDimensionTooSmall`); `Canceled` only if cancelable, `OutOfMemory` only if allocating.

- [ ] **Z3** Errors, not panics, for bad input; checks no longer need the reference's order.

- [ ] **Z4** Doc comments say where a routine panics explicitly (only unrecoverable states, above all possible memory corruption from assembly).

- [ ] **Z5** User docs without BLAS jargon (`xerbla`, `LDA`, `INCX`).

- [ ] **M14** The public layer's step 7: the matlab port moves onto the layer (`Mat(T)` and views). **Needs:** M13.

- [ ] **M6** A type-agnostic matlab port, as the BLAS and LAPACK are: `CxMat` and its routines generic over the element type (`comptime T`, any float type), not hard-coded `f64`, with a caller-chosen result type where MATLAB fixes `double` (an explicit `comptime` type, the nearest Zig has to `@floatFromInt`'s inference). Pulled into the Queue by the user, 2026-10-31. On `@Matrix` where the compiler has it, as M7 settled. **Needs:** M14 (the port on `Mat(T)`; whether it keeps its own types is decided there).

- [ ] **N7** `fft.Plan`: complex 1-D, portable, every float type (PocketFFT's mixed radix and Bluestein, from its C++ header at a pinned commit). The correctness baseline every later FFT item is tested against. Add PocketFFT's `LICENSE` section and per-file provenance with the first ported code. **Needs:** nothing (N6 done).

- [ ] **N8** `fft.RealPlan` (PocketFFT's real transforms). **Needs:** N7.

- [ ] **N9** `Layout` batches, run concurrently through `io`, as many at once as the scratch holds. **Needs:** N7.

- [ ] **N10** Audio: port pffft's f32 SIMD kernels faithfully (design.md, Decisions), within one transform, sizes `2^a 3^b 5^c`, and the unordered convolution (`forwardUnordered`, `convolveAccumulate`, `backwardUnordered`): SSE on x86_64, NEON on aarch64 and 32-bit arm, generic code on the PowerPC 750. Where a kernel repeats one rstd has (BLAS, LAPACK), call rstd's. Its vector primitives (interleave, transpose, half-swap) are M16's. Add pffft's `LICENSE` section with its first ported code. **Needs:** N8, M16.

- [ ] **N35** FFT kernels for tier 1 and 2 targets pffft lacks or serves poorly (AVX on x86_64, say), derived from its own or new. **Needs:** N10.

- [ ] **N29** `rstd.signal` window functions: Hann, Hamming, Blackman, Kaiser, sine and KBD, symmetric or periodic. **Needs:** nothing (N6 done).

- [ ] **N30** `signal` STFT and its inverse (weighted overlap-add) over `RealPlan`. **Needs:** N8, N29.

- [ ] **N31** `signal` uniform-partitioned convolution over the unordered forms: fixed latency, no allocation per block. **Needs:** N10.

- [ ] **T138** Compensated summation (Neumaier) as an option on the BLAS reductions (`sum`, `asum`, `dot`), still bit-identical for any thread count; `rstd.stats` uses it through the BLAS. Benchmark it against the plain forms on the M6 Mac mini and Zen 3, then choose the defaults. Pulled into the Queue by the user, 2026-10-31.

- [ ] **N11** `stats.mean`, `variance` and `stdDev` through the BLAS `sum`, with its accuracy options once T138 adds them. **Needs:** nothing (N6 done).

- [ ] **N12** `stats.median` and the quantiles, all nine methods. **Needs:** nothing (N6 done).

- [ ] **N13** `stats.histogram` and `histogramInto`. **Needs:** nothing (N6 done).

- [ ] **N33** Type-agnostic `rstd.stats`, as the BLAS: integer and `Complex` inputs beside floats, the result type chosen by the caller (`mean(io, R, xs)`, the element type from `xs`). **Needs:** N11.

- [ ] **N14** `sort.unique`, `intersect`, `setUnion`, `setDiff` and `isMember`, each NaN distinct. **Needs:** nothing (N6 done).

- [ ] **N32** `Complex` in `rstd.sort`, in MATLAB's order (by magnitude, then angle), used by every ordering operation, `stats.median` included. **Needs:** nothing (N6 done).

- [ ] **N1** N-d array stubs (design.md, "Planned: N-d arrays"). **Needs:** M10, whose `View(T)` and `Mat(T)` they extend.

- [ ] **N15** `stats.cov` and `corrcoef`: `syrk` over the centered columns. **Needs:** M11.

- [ ] **M4** Time `linalg.vector`'s inline `@Vector` forms against the BLAS and LAPACK at sizes 4 to 16 per operation and type (`dot`, `length`, `mul` matrix by vector and by matrix, `solve`), on the M6 and Zen 3 under their timing rules: confirm the inline limit of 8 (src/linalg/design.md, the public layer's Dispatch) or move it. `linalg.matrix`'s LU (`det`, `inv` past 3 by 3, every `solve`; 5f7ae56f) is scalar loops, not the design's inline `@Vector` code: time it too, and vectorize its column updates if it trails. **Needs:** M11.

- [ ] **T21** Zen 3 (flyoff): short real `dot`s (0.30 of OpenBLAS at 1024) and complex `gemm` (0.79 to 1.01) trail OpenBLAS. Tier 1 (user, 2026-10-30): re-measure after T100, then fix what's left; the timing session verifies. Diagnosed 2026-10-30, no code yet: Zen's `sdot` kernel (`Real(f32,4).dotSum`) is bound by its conversions (32 `vcvtps2pd` per 64 elements, about 70 GB/s: the cost of the `f64` rule), and on AVX2 its two `@Vector(32, f64)` accumulators fill all 16 ymm registers (about 25 moves and 3 spills per iteration): halve them. `cgemm`/`zgemm`: tiles, `P`/`Q` and the inner loop match OpenBLAS's, but `cgemm_kernel_8x2_haswell.c` reads two or three B panels at once (m8n4/m8n6) with 12 accumulators where ours has 8, exactly FMA latency times throughput: port that shape. `blas_bench`'s `kern3` section (59faab6c, any host) times the kernel alone against rstd's whole `gemm` on one worker; its first run (busy Zen 3, 2026-10-30) has the driver keeping 0.89 to 0.98 of the kernel for `s`, `d` and `z`, and `cgemm`'s whole call at 1.04 to 1.08 of its kernel block: the `cgemm` gap is the kernel's. For the timing session: `kern3`, `flyoff.sh zen3`'s `win-1t`/`lx-1t` `dot s` rows and the level-3 `c`/`z` rows.

- [ ] **T102** Strided `saxpy` and `zscal` (inc 2) in the arm64 Linux VM halved between 6d66239 and 28171adc for both rstd and the converted kernel (`saxpy` 54 to 28 GB/s, `zscal` 139 to 71; R4, tools/blasbench/results-flyoff.md). Read 2026-10-30, not measured: nothing on rstd's path changed (one kernel call below the threshold, vortex's slot is still the displaced `@Vector` kernel, `vector.zig`'s `axpy` and `scal` and the per-core flags unchanged). The bench's `asm` column reads the same function rstd calls, so the two moved together, while its own copy of that source (`vec x4`) stayed at 55 and OpenBLAS's `axpy.S` at 28: the drop follows that one copy in the binary (likely code placement) or the guest. Left, in the timing session: `flyoff.sh m6`'s inc-2 `saxpy`/`zscal` rows, then align the copy or close this as the guest's.

- [ ] **T114** Accelerate's `potrf` (T35), the timing part: in the timing session, whether routing beats rstd's threaded `potrf` (its correctness is checked: T62's `T114-<host>` steps).

- [ ] **N2** Sparse matrix stubs (design.md, "Planned: sparse matrices"). **Needs:** M10.

- [ ] **T150** Sparse matrices (src/linalg/design.md, "Planned: sparse matrices"; Randy ChatGPT discussions): decide with the user whether the kernels live in `rstd.linalg` or stay in srapnol's `src/sparse.zig` (its DES-012), then build CSC with `A x`, `Aᵀ y`, sparse `axpy` and `dot`, `scal`, submatrix extraction, Gustavson `spgemm` (symbolic and numeric) and sparse triangular solves, bit-identical for any thread count. **Needs:** that decision.

- [ ] **T151** Time CSC `A x` and `Aᵀ y` against dense `gemv` over fill (1% to 60%), size (100 to 10 000) and pattern (random, banded) on the M6 and Zen 3: the crossover for a dense-or-sparse switch (ChatGPT's 20% to 30% cites no source), and whether AVX2's gather beats scalar loads for `Aᵀ y` (Randy ChatGPT discussions). **Needs:** T150.

- [ ] **Z6** The last item: an rstd-wide Zig-idiomatic pass (slices, exact error sets, plain error names, documented panics), and the file layout: fewer, bigger files where pieces belong together (`lapack/`'s `schur.zig`, `real_schur.zig` and `complex_schur.zig`, after T127 flattened it). **Needs:** Z1 to Z5.

## Known issues and enhancements

Picked up after the Queue: tier 3 and 4 bugs, follow-ups, tuning, optional and Blocked work.

- [ ] **T144** `powerpc64le`, `powerpc64le-power10` (ReleaseSafe) and `powerpc64` do not compile: LLVM cannot select a `v1i128` bitcast of an `f128` load in `tridiag_evd.qrAlgorithm` (`tridiag_evd.zig:163`, issues.md). Tier 3; no workaround in rstd. **Blocked:** upstream.

- [ ] **T152** `linux-powerpc64le-power10`: six LAPACK tests fail in `Complex(f32)` at denormal or tiny scales: `getrf`'s denormal pivot (dividing by 1.1e-44 gives `inf`, `nan`), `syevd`, `gees` at scale 1e-30, `gesvd` (NoConvergence), `gels` (Singular), the evd NaN test. The power9 rig (`linux-powerpc64le`) passes them all: likely power10 code (LLVM or QEMU) flushing or overflowing in the complex division. Not traced. Tier 3.

- [ ] **T153** `linalg.matrix`'s `@Matrix` path and its column code (5f7ae56f) can differ in bits: on MIPS with MSA the fork's Sema fuses `M v`, `v M` and `M N` products as `std_float_options` says, while `fma.mulAdd` (the column code) stays unfused there (`fma.compiles`: LLVM compiles MSA's `fmadd` wrong, issues.md). Also, with both operands NaN, the column code adds `term + sum` where Sema adds `sum + term`, so the NaN payload may differ; the bits test ("mul through @Matrix has its column code's bits") uses finite values only. Fix: make the fallback follow Sema there, or keep Sema off `@mulAdd` under MSA in the fork. Tier 3 (fork compiler only).

- [ ] **T145** x87, Debug only (`i386`, `rstd-test-x87`; the ReleaseSafe forms pass): "geev on special matrices" fails `f64` n = 2 kind 3, the rotation `[0 -1; 1 0]`: an eigenvalue's `x[kbig].im` is past `tolerance(F, n, 1)` (`evd.zig:852`). Likely `f80` intermediates in the 2 x 2 standardization; not traced. Tier 3.

- [ ] **T146** x87 (`rstd-test-x87`): matlab's `rat` loses precision (22/7 comes back with numerator 3.14306640625; 3.1416015625 for π, `rat.zig:47`, `:55`), and `base` renders the wrong digits (`base.zig:75`). Not traced; M6 reworks both. Tier 3.

- [ ] **T147** "the registry read works on Windows" (`_detect/arm64.zig:1011`) fails on the 32-bit `rstd-test-x87`: `InstallTime` is missing, likely from WOW64's 32-bit view of `HKLM\SOFTWARE` (`KEY_WOW64_64KEY` reads the native one). Matters for 32-bit arm Windows detection too. Tier 3.

- [ ] **T148** Big-endian MSA (`mips32r5`'s `p5600`): `util.abs` clears the low word's top bit of an `f64` instead of its sign (`|x|` off by 2^-22 relative, `-inf` kept), so `amax`, `amin`, `iamax`, `iamin`, `nrm2`, blat1 and utest fail, and `getrf` pivots wrong (`|L| > 1`). The `f64` to `@Vector(n, u32)` bitcast's lane order on big-endian MSA is the likely cause; a mask in `F`-wide lanes avoids it. `mips64r6` fails the same `getrf` tests (its run ends at T149 before the level-1 ones); the little-endian rigs pass. Tier 3.

- [ ] **T149** MIPS R6 (`mips64r6`, `mips64r6el`): "3m: the three passes' copies and kernel give the complex product" panics with integer overflow (no stack on the rig). Not traced. Tier 3.

- [ ] **T136** Evaluate the x86_64 `peels` slots (`vector.peeling`: `asum`, `sum`, `ddot`, `cdot*`, `zdot*`, `dgemv_n`), whose converted assembly (`asum_sse*.S`, `zasum_*.S`, `zsum_sse*.S`, `dot_sse2.S`, `dot_atom.S`, `ddot_bulldozer.S`, `zdot_*.S`, `dgemv_n*.S`) rstd skips for `@Vector` kernels on the cores that build it, because its bits follow an operand's address. Never benchmarked on those cores. Deviating from OpenBLAS needs a very high bar (user, 2026-10-31): confirm the dependence upstream, look for a generator patch that keeps the kernel (as `zaxpy_sse2.S`'s add order, T62), time both on each core (forced on this host, indicative), then decide with the user.

- [ ] **T135** Evaluate the routed extremes (`i?amax`, `i?amin`, `i?max`, `i?min` and values; the first NaN wins, `src/linalg/blas/readme.md`, One contract): every architecture runs rstd's kernels, not upstream's. On `zen` `isamax` is 0.24 to 0.41 of OpenBLAS below 16M floats (flyoff, 2026-10-27). Weigh the contract against upstream's NaN behavior and speed per architecture; patch upstream where it is wrong, if that keeps it; decide with the user.

- [ ] **T134** Evaluate routed `nrm2` (exact `f64` sum for `s`, `c`; Blue's for `d`, `z`): rstd's is faster on `zen` (1.2 to 11.7) and the M6 (4.3 to 17), unmeasured elsewhere. issues.md lists upstream kernels wrong on special values: confirm which are wrong rather than only different, patch those that can be, time the rest per architecture; decide with the user.

- [ ] **T133** Evaluate routed `sdot` and `dsdot` (one `f64` sum): short `sdot` is 0.30 to 0.65 of OpenBLAS on `zen`, 0.37 at 256K floats on the M6. Weigh the accuracy contract against upstream's `f32` sums per architecture; decide with the user.

- [ ] **T132** The gemm driver allocates the whole `BUFFER_SIZE` packing buffer (20 MB on the AVX x86_64 cores) per call, however small the problem; allocating only what the blocking uses for `m`, `n`, `k` would keep the bits and save the memory (found by T62: the x86_64 rigs' heap is 64 MB for it).

- [ ] **T131** dqds for a bidiagonal's values alone (from T126): Reference-LAPACK's `xLASQ1` to `xLASQ6` (`dlasq1.f` to `dlasq6.f`, about 1,070 lines of code), which `xBDSQR` calls when no vectors are wanted: faster than the QR iteration and accurate to a few eps relatively. Now `gesvd` and `gesdd` with `.n` run `xBDSQR`'s QR sweeps, already accurate (T126's test).

- [ ] **T120** faer follow-ups (split from T117, 2026-10-30): `@Vector` the fused bidiagonalization pass (scalar now, faer's four partial sums); run divide-and-conquer's two halves in parallel if rstd's threading gains task trees; port faer's multishift QR with aggressive early deflation (n >= 75) and blocked Hessenberg reduction (n >= 256), so large `geev` isn't slow; thread the fused tridiagonalization and Hessenberg passes with fixed chunks; make `sytrd`, `gehrd`, `hseqr`, `trevc` public; faer routines beyond what M2 needs (T123 has those): Cholesky rank-1 update, pivoted Cholesky, full-pivoting LU, `ggev`/`gges`; race `gesvd`/`gesdd` against OpenBLAS's LAPACK, reference LAPACK and Eigen in the timing session (T73's bench).

- [ ] **T121** ReLAPACK follow-ups (split from T118, 2026-10-30): a blocked `gbtrf` from LAPACK's `dgbtrf` (rstd runs the unblocked `xGBTF2`, since ReLAPACK's recursion is wrong); race faer's Bunch-Kaufman (different storage) against `sytrf`. `tgsyl`'s `xLATDF` estimator is a local loop in `sylvester.zig` copying the update step of `gecon`'s private `Gecon` (it needs the last vector, which `gecon` doesn't return), not `norm.zig`'s `lacn2` over `latrs` as the docs say: share it or fix the docs; `tgsyl` not timed yet. Sweep: `relapack` on every rig, and with `-Drelapack=true` on one or two.

- [ ] **T122** Small-matrix follow-ups (split from T119, 2026-10-30): the per-core small kernels for SVE, POWER10 and LASX (tier 3 or untestable here); the small-path thresholds are upstream, tuned for its kernels, so time them under T05; `sbgemm` small kernels off x86_64; from T119's port: `nn`'s dot rows are read in place, not through upstream row-major copy (same bits, likely slower on AVX-512 hardware), and Cooper Lake's `sbgemm` small path keeps the portable bits, not `vdpbf16ps`'s pair rounding (by design, documented). Sweep: `3m`, `small`, `interface: gemm` on every rig, `aarch64-sve`, `powerpc64le-power10` and loongarch64 first (their permits allow the small path), `blat3_3m` on `x86_64` and `i386-max`.

- [ ] **T47** Hardware-only paths never run: Linux ppc64 `mfpvr`, the mips `cpucfg` decode, the sparc `sllx` rewrite, the FreeBSD arm64 probe, NEON and AltiVec flushing (rigs only). The tier 2 ones (the Windows arm64 registry and feature probes, the arm64 sysfs fallback, Accelerate's pre-13.3 symbol fallbacks) are tested with recorded inputs (b1cad78a); only a macOS 11 or 12 machine would run the last for real. Possible enhancements: read every Windows `CentralProcessor` key, map Cobalt 100 to `neoversen2` (as LLVM), alias Snapdragon X Elite to `neoversen1` (as OpenBLAS's build-time detection). **Blocked:** hardware.

- [ ] **T22** C and intrinsics kernels not yet ported, levels 1 and 2, by T20's strategy, `@Vector` first (`vector.zig`; the routines' `coverage.md` files have the cells). x86_64 and riscv64's zvl cores have no gaps left; their AVX-512, FMA4 and RVV kernels are tested natively or under TCG, not on that hardware. Left (tier 3, or tier-2 cores no rig runs): arm64's SVE `symv` drivers (A64FX, ARMv8 SVE, ARMv9 SME, Cortex-A510, A710, X2, Neoverse V1) and `ger_sve_v1x3.c` (A64FX, V1); POWER8 to 10 (`asum`, `scal`, `rot`, `swap`, `cdot`, gemv, POWER10's `sbgemv`); MIPS MSA (gemv on I6400, I6500, P6600, Loongson 3R3 and 3R4, P5600); zarch (z13, z14); wasm SIMD128 (gemv); riscv64's C910V (**Blocked:** no `xtheadvector` in LLVM or Zig). **Needs:** T20.

- [ ] **T73** Benchmarks above level 1, landed 2026-10-30 (765b30fa..37d4b2ba): `ger`, `symv`, `trsv`, `gbmv`, `symm`, `syrk`, `trmm`, `trsm` and a `lapack` section, against reference LAPACK, Eigen (benchmark-only, never ported) and AOCL-libFLAME through Nix, an Accelerate column in its own process, `flyoff.sh zen3` and `m6` under the R4 rule. Left: the first real flyoff runs and Nix builds (untested end to end); confirm the Windows idle rule (processor queue length under 2); `sygst`/`hegst`, `trsyl`, `hetrf` entries; rstd's `geqrf`/`syevd` once T111 lands; a `-Drelapack` switch; `trmm` at n = 64 (0.09 to 0.40 of OpenBLAS, busy machine). **Needs:** T06 for threaded rows.

- [ ] **T94** Fact-check items still open after T78: 32-bit x86 Windows threads starting at x87 0x27f (T03's observation; no source, no test records it: add a test); "MKL leaves FTZ to the application" (from memory); Tegra 2 as VFPv3-D16 (a forum post only); `symmetric.zig`'s "4 elements of slack past each".

- [ ] **T20** Kernel strategy: race `@Vector` against the converted assembly per core on a quiet machine (only vortex is done, `tools/blasbench/results-m6.md`); keep per-core code only where Zig loses. riscv64 races (from T22): `@Vector` against the `asm` `axpy` on zvl128b and x280, and fixed-length `@Vector` against VLEN-agnostic RVV on hardware wider than a core's minimum VLEN.

- [ ] **T05** Thresholds: OpenBLAS's defaults and per-core overrides, tuned by benchmark; then the overrides left out of phase 1 (z13 `axpy`, Windows-on-ARM64 `gemv`) if the numbers agree. **Needs:** T01, T73. T06's smoke run (busy Zen 3, `--threads=4`): threaded level-1 calls were slower than serial at the sizes it ran, so the thresholds need this tuning. Also time T07's split Accelerate calls against one whole call; and the bench's "Accelerate multithreaded" column is single-threaded while rstd's init sets `VECLIB_MAXIMUM_THREADS=1` process-wide, so a fair column runs Accelerate in a separate process. From T113's table: race AOCL-libFLAME's small-matrix AVX2/AVX-512 kernels (`getrf`, `potrf`, `potri`, `getrs`, `getri`, `geqrf`, `gelqf`, small `gesvd`; its OpenMP paths only after a same-bits check), libflame's restructured QR for `gesvd` and Hermitian eigenvalues (Van Zee et al., TOMS 2014), `FLA_Hevd`/`FLA_Tevd`, and `FLA_Eig_gest`/`FLA_Sylv` against ReLAPACK's `sygst`/`trsyl`. From T101 (src/linalg/blas/threading.md): OpenBLAS threads x86_64 `ddot`/`zdot` and arm64's `dot.c` above 10000 (arm64 with its own neoversev1 and a64fx tables) and never `sdot`, `cdot`, `dsdot` or `dot.S`; rstd threads every precision above 10000. Threaded `cgemm` at n = 128 in the arm64 VM varies 198 to 249 GFLOP/s: each call (about 75 µs of work) wakes five sleeping `Io.Threaded` pool threads (a task allocation, a mutex and a condvar signal each), whose wake time varies under HVF; no spinning server (src/linalg/blas/threading.md), so a threshold that keeps such calls serial is the fix to time.

- [ ] **T88** Time T25's races on an idle machine: portable against `@Vector` for `sum`, `sdot` in `f64` and in `f32` lanes, `nrm2` and the extremes, natively and as an x87 build (a scalar-FPU stand-in); bench and runner in the scratchpad (`vecsum/_vecsum_bench.zig`, `timeit.ps1`). Also race `vortex`'s converted `sum` against `@Vector` on the M6 (kept in `vector.displacing.keep` until then).

- [ ] **T26** vortex: a tuned `@Vector` `rot` (the assembly still leads by 2 to 17% in L1); race the Zig complex `gemm` against OpenBLAS's vortex assembly on the M6.

- [ ] **T27** gemm tuning: faster Zig packing (whole-vector transposes in `oncopy`, `incopy`), `C`-tile prefetch at kernel entry, a larger `KC` on Apple cores (src/linalg/blas/comparisons.md (9)).

- [ ] **T28** Measure whether composed `symv`/`hemv` pay off on cores with assembly gemv (atom, opteron, arm64); today they compose only where the core has its own gemv.

- [ ] **T29** The bfloat16 intrinsics (cooperlake, sapphirerapids, Neoverse N2, V1, V2, POWER8 to 10); `sbdot` could call the core's own `sdot`. **Needs:** T20.

- [ ] **T31** Assembly where `@Vector` cannot reach: arm64 SVE, SME (`vortexm4`), POWER10 MMA, AMX and AVX-512 BF16 for `sbgemm`; Skylake-X only if under 85% of its kernel (src/linalg/design.md, Implementation plan; `src/linalg/blas/blas3/gemm/notes.md`). **Needs:** T20.

- [ ] **T32** BLIS assembly kernels, decided per kernel where benchmarks favor them (BSD-3, with a license section and per-file provenance). **Needs:** T20.

- [ ] **T39** sparc64 is tier 4 (user, 2026-11-01): its images outgrew OpenBIOS's 0x400000 to 0x800000 window again (tier 3 at 5ddbd422: `blas3-tri` in ReleaseFast, 34 of 46 parts in Debug), so both rigs are stubs (`build/qemutest.zig`) and out of the sweeps. Bringing it back needs a real fix first: map more memory at boot through the Open Firmware client interface (claim and map), so one image fits. Open there too: the `-Oz` miscompile (issues.md, needs a reduced case) and the window itself (our observation; OpenBIOS's source unchecked).

- [ ] **T71** Linux rigs with no item: loongarch64 (the `AT_HWCAP` LSX/LASX probe), big-endian powerpc64 (`AT_PLATFORM`, the `64be` variant), 32-bit mips.

- [ ] **T41** `zrot_ppc440.S` (`ppc440fp2`) has a zero-stride patch no rig has run; no rig runs `ppc440` (it shares code with `ppcg4`). The 32-bit `ppc440` and `ppc440fp2` kernels rebuilt in 719a41b are compile-checked only.

- [ ] **T74** No CI runs a BLAS test: rstd's `.gitlab-ci.yml` is rgui's copy (`EXE_NAME: "rgui"`, Zig 0.16.0). A job for `rstd-test`, `rstd-check-macho`, `rstd-check-uefi` and `rstd-check-zig017` (on Zig 0.17.0) would catch drift between rig sweeps; optional if testing stays manual.

- [ ] **T53** macOS VM extras: base-image overlays (`qemu-clone`, `-Dmac-base`), a `Forgejo` self-hosted runner, a setup that needs no clicking, less guest load when headless (HID driver, `SecurityAgent`), and freeze detection that asks before killing.

- [ ] **T104** Optional, best effort: the reference's GCC build contracts the MSA kernels' `y0 += da_vec * x0` into `fmadd` on every MSA core, but rstd fuses vectors only on `loongson3r4` (T86). The same inline MSA `fmadd` from `fma.zig` could model `p5600`, `i6400`, `i6500` and `p6600` as the reference builds them; QEMU rounds neither way, so it rests on reading the assembly.

- [ ] **T76** Optional, from src/linalg/blas/comparisons.md: `Dot2` (compensated `f32` lanes) if T25 finds the `f64` `sdot` costly in cache; BLIS's `JR`/`IR` ways for threads sharing an L2; an A-panel prefetch beside T27's `C`-tile one; small-`n` `sgemv_n` tuning (`tools/blasbench/results-m6.md`).

- [ ] **T54** Optional: an inline-asm `maddf` for R6 scalar multiply-adds in the main library, if R6 speed ever matters (LLVM calls `fmaf`).

- [ ] **T52** On each Zig update, recheck and drop the workarounds once fixed (T48 and T49 track lld's SPARC and `zig fetch`): the powerpc64le two-float return (`zdot` patch), arm soft-float struct returns, `f80` under `x86_64_sysv` on Windows, the comptime naked-function call, ppc64 big-endian ELFv2, the two-lane `f32` vector stored as 16 bytes, sparc64 `_Qp_sqrt` (`check.sqrtW`), sparc64 `-Oz`, `SafeAllocator`'s CPU count, the `Thread` override doc, run-step timeouts and output, the interrupted-link cache entry, inline-asm errors exiting 0 (`checkedObject`), the MIPS64r6 select copy and spill (comptime `conj`), MIPS MSA vector `fma` (MSA cores differ from the reference's fused bits until fixed), arm CPU models without an FPU (rig pins), WHPX passing AMD CPUID (macOS rig on TCG), QEMU `xststdcsp` (portable `nrm2` avoids `isInf`), no fiber `Io` tests on SVE targets (std's fiber switch clobbers `ffr`, T40). Also, wasm bitcode getting its feature flags twice (`asmCheck`'s `-disable-llvm-verifier`). The Zig entries were last checked on dev.751. Also: building the library for an R6+MSA target CPU (`i6400`) trips LLVM's `Cannot copy registers` assertion (`MipsSEInstrInfo.cpp:261`) in `blas2.interface`, `blas3` and `_test`, likely the R6 select bug already recorded (fixed upstream after 23.1.x); the rigs' main CPU has no MSA, so it never shows there (T89).

- [ ] **T43** Cores no rig runs: the Skylake-X family (AVX-512) and the Bulldozer family (FMA4, XOP), including their new x86_64 `rot`, `scal` and `gemv` paths, which are compile-only: run `x86_64-whpx` on a host that has them. **Blocked:** hardware.

- [ ] **T44** The AVX-512-off workaround for the self-hosted backend (`vectorBits` in `blas3/gemm/kernel.zig`, `zig_backend != .stage2_x86_64`) has no issues.md entry and no recorded reason (since 2cf8757); 0.18.0-dev.751 compiles 512-bit `@Vector` multiply-adds there. Run the gemm tests with `-fno-llvm` on AVX-512 hardware, then drop it or record it. **Blocked:** hardware.

- [ ] **T45** POWER rigs stay at `.safe` (dcb4628): LLVM's FastISel stores an `f64` held in a v-register with `stfd` (issues.md); drop `.safe` once fixed. **Blocked:** upstream.

- [ ] **T46** No `powerpc64` (970) AltiVec scan: LLVM stops with "Invalid global physical register". Scan once it builds; rerun the G4 `-femit-asm` scan after LLVM upgrades, since the SLP vectorizer is not proven blocked by `util.scalar`. **Blocked:** upstream.

- [ ] **T48** 32-bit SPARC never runs (lld cannot link it; includes the `%icc` index-tie fix). **Blocked:** upstream.

- [ ] **T49** Linux rigs: `linux-i386` is a stub (QEMU x87 page-crossing bug); `s390x` and `mips64el` wait on `zig fetch` taking single files. **Blocked:** upstream.

- [ ] **T50** The `powerpc-macho` rig is a stub; finish it once the powerpc-macos Zig port lands. **Blocked:** upstream.

- [ ] **T51** translate-c: re-enable the package once its bundled Aro builds on current Zig. **Blocked:** upstream.

- [ ] **T105** The power10 objects in a 32-bit POWER build don't link and are untested: the dev toolchain's LLVM can't build 32-bit compiler-rt for POWER4 to POWER10 (`64bit` crashes, VSX with `fpcvt` fails an `f128` conversion, POWER9 vector and ISA 3.0 features error; Zig 0.17.0 builds them all). `powerpc-power9` (T11) runs power8 and power9 only. Any real 32-bit POWER8 to 10 build of rstd fails with this toolchain, against src/linalg/design.md's filling-upstream-gaps decision. **Blocked:** upstream (LLVM; issues.md, 32-bit compiler-rt).

- [ ] **T141** The macOS VM (`qemu-test-x86_64-macos`, Skylake-Client under TCG) fails "relapack potrf: A = L Lᴴ from a known L": `NotPositiveDefinite` on every core and the `@Vector` kernels, so Accelerate's AVX2 path (BLAS calls inside ReLAPACK's `potrf`) or TCG; Rosetta passes, AVX advertised or not. Its `cholesky` `potri` test then ran over an hour silent. First on the Intel Mac once it is set up (`potrfCases` now names the type and size); a filtered VM run with routing on and off would tell the two apart.

- [ ] **T140** A pass over all of rstd for std declarations this toolchain marks deprecated (`std.builtin` for `std.lang`, `Arch.isMIPS` for `isMips`, and the rest), using the replacement wherever Zig 0.17 has it too (user, 2026-11-01). New code avoids them already.

- [ ] **N16** FFT: f64 SIMD (pffft's double) and a cache-blocked split for `n` past L2, for numerics. **Needs:** N10.

- [ ] **N17** FFT: SIMD across a batch's transforms, for images. **Needs:** N9.

- [ ] **N18** `fft.TrigPlan` (DCT and DST types 1 to 4) and `Mdct`. **Needs:** N8.

- [ ] **N19** `fft.StaticPlan` and `StaticRealPlan`. **Needs:** N8.

- [ ] **N20** FFT: route vDSP's sizes (`f 2^n`, `f` in 1, 3, 5, 15) to Accelerate on macOS by default, as the BLAS does, converting its real layout and scaling; a way to force native. **Needs:** N8.

- [ ] **N21** FFT: threads within one large transform, by a fixed split. **Needs:** N9.

- [ ] **N22** N-d array bodies. **Needs:** N1.

- [ ] **N23** `fft.transformAxes` and `realForwardAxes`. **Needs:** N22.

- [ ] **N24** `rstd.stats` along a dimension. **Needs:** N22.

- [ ] **N25** `rstd.stats`: an option to omit NaNs.

- [ ] **N26** Sparse bodies: `Csc` first, `mulVec`, `cg` and `gmres`. A likely first user is `linalg/optim`'s sparse Jacobian. **Needs:** N2.

- [ ] **N27** Sparse: BiCGSTAB and MINRES. **Needs:** N26.

- [ ] **N28** Sparse direct factorizations. **Blocked:** a license decision (design.md, Licenses).

- [ ] **N34** FFT: native rewrites of pffft's kernels, kept only where a benchmark shows them faster, or as fast and cleaner. **Needs:** N10.

## Done

Cleared 2026-11-02; landed work is in `review`'s log.

## Waiting on the user

- T150: whether the sparse kernels live in `rstd.linalg` or stay in srapnol's `src/sparse.zig`.
