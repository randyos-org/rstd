# rstd to-do

Open work in priority order, top to bottom. **Needs:** items that must land first; **Blocked:** an outside cause.

## In progress

### Working on

- None.

### In the background

- Housekeeping H1 to H112 (user, 2026-11-07): landed on `housekeeping` with the generator items, `review`'s Intel Mac scripts merged in (now in `tools/macos/`); RandyOS's follow-ups staged on its `main`; waiting for the user's review and a first build.

## Queue

- [ ] **N45** `ssh-test-x86_64-macos` runs on the Intel Mac (`macmini-intel`), not under Rosetta on the M6: an option for the Intel Mac's alias beside `-Dssh-mac` (default `macmini-m6`, help text "Apple Silicon"), the x86_64 step pointed at it; Rosetta stays reachable for comparison. Then T141 can rerun there. **Needs:** N43, N44.

- [ ] **N46** The qemu rigs on the Intel Mac (`-Dtest-host=macmini-intel`): `tools/nix/flake.nix` follows nixpkgs-unstable, which dropped x86_64-darwin (26.11 release notes), so give x86_64-darwin a `nixpkgs-26.05-darwin` input (supported to the end of 2026) and check `tools/test_host.zig`/`run.sh` there (`hvf` for x86_64 guests). **Needs:** N44.

- [ ] **N47** `support.md`: the Intel Mac replaces Rosetta 2 and the macOS VM as x86_64 macOS coverage (lines 35, 64 "no Intel Mac is at hand", 77) and joins the test hosts. **Needs:** N45, N46.

- [ ] **N48** rstd `CLAUDE.md`, This machine: the Intel Mac is set up (alias `macmini-intel`, `tools/macos/configure_mac_host_intel.sh`, `configure_mac_nix.sh`), not "coming"; how its ssh and rig steps are reached; drop "the macOS VM runs nothing" once N49 retires the VM. **Needs:** N45, N46.

- [ ] **N37** msdf `Shape.calcBounds` (`Shape.zig:97-101`) starts `right` and `top` at `floatMin(f64)`, the smallest positive float, not msdfgen's `-LARGE_VALUE`: a glyph wholly below the baseline (`_`) or left of 0 gets top or right 0, so a wrong box. Use `-floatMax`; test with such a shape.

- [ ] **N38** `buildTextAtlas` (`src/rgl/font/render/atlas.zig:159`) looks up each byte of `text` as a code point, so non-ASCII UTF-8 gives wrong glyphs ("é" as two Latin-1 glyphs), while `layoutText` decodes UTF-8 (`layout.zig:112`). Decode with `std.unicode.Utf8View` and return `InvalidUtf8` as `layoutText` does; test a non-ASCII string.

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

- [ ] **T102** Strided `saxpy` and `zscal` (inc 2) in the arm64 Linux VM halved between 6d66239 and 28171adc for both rstd and the converted kernel (`saxpy` 54 to 28 GB/s, `zscal` 139 to 71; R4, tools/blasbench/results/results-flyoff.md). Read 2026-10-30, not measured: nothing on rstd's path changed (one kernel call below the threshold, vortex's slot is still the displaced `@Vector` kernel, `vector.zig`'s `axpy` and `scal` and the per-core flags unchanged). The bench's `asm` column reads the same function rstd calls, so the two moved together, while its own copy of that source (`vec x4`) stayed at 55 and OpenBLAS's `axpy.S` at 28: the drop follows that one copy in the binary (likely code placement) or the guest. Left, in the timing session: `flyoff.sh m6`'s inc-2 `saxpy`/`zscal` rows, then align the copy or close this as the guest's.

- [ ] **T114** Accelerate's `potrf` (T35), the timing part: in the timing session, whether routing beats rstd's threaded `potrf` (its correctness is checked: T62's `T114-<host>` steps).

- [ ] **N2** Sparse matrix stubs (design.md, "Planned: sparse matrices"). **Needs:** M10.

- [ ] **T150** Sparse matrices (src/linalg/design.md, "Planned: sparse matrices"; Randy ChatGPT discussions), in `rstd.linalg`: CSC with `A x`, `Aᵀ y`, sparse `axpy` and `dot`, `scal`, submatrix extraction, Gustavson `spgemm` (symbolic and numeric) and sparse triangular solves, bit-identical for any thread count. **Needs:** nothing.

- [ ] **T151** Time CSC `A x` and `Aᵀ y` against dense `gemv` over fill (1% to 60%), size (100 to 10 000) and pattern (random, banded) on the M6 and Zen 3: the crossover for a dense-or-sparse switch (ChatGPT's 20% to 30% cites no source), and whether AVX2's gather beats scalar loads for `Aᵀ y` (Randy ChatGPT discussions). **Needs:** T150.

- [ ] **Z6** The last item: an rstd-wide Zig-idiomatic pass (slices, exact error sets, plain error names, documented panics), and the file layout: fewer, bigger files where pieces belong together (`lapack/`'s `schur.zig`, `real_schur.zig` and `complex_schur.zig`, after T127 flattened it). **Needs:** Z1 to Z5.

## Known issues and enhancements

Picked up after the Queue: tier 3 and 4 bugs, follow-ups, tuning, optional and Blocked work.

- [ ] **N49** Deprecate the macOS VM, replaced by the Intel Mac: the `x86_64-macos` rig and its `qemu-setup-`, `qemu-install-`, `qemu-start-` and snapshot steps (`build/qemu/macos.zig`, `tools/macos_vm.zig`, `tools/macos/configure_macos_vm.sh`, `install_macos.sh`, `snapshot_macos.sh`), their docs and `.macos-vm`; retire T53 with it. **Needs:** N45, N46.

- [ ] **T158** matlab on an empty real matrix: `lu.zig`'s `Lu.init`, `solve`, `lange1` and `rcond` (so `det`, `inv`, `mldivide`) pass its `lda` of 0 to LAPACK, which answers `InvalidLda`, and `ops.narrow` makes that `unreachable`. Give 0x0 MATLAB's answers (`det([])` 1, `inv([])` empty, `norm` 0) and test them, as `chol` now does.

- [ ] **T159** `build/blaskern.zig`'s `-Dblas-llvm` help text ("unset uses it only where the kernels need it") is stale: unset now always uses LLVM.

- [ ] **N39** GPOS format-1 `ClassDef.classOf` (`tables/gpos/ClassDef.zig`) computes `start + count` and the record offset `2 * (glyph - start)` in `u16`: they overflow (a panic in safe builds) past glyph 65535 or from 32768 entries. Widen them; test such a ClassDef.

- [ ] **N40** `png.decode` leaks the image when un-filtering fails: no `errdefer` frees it. Add one; test a corrupt filter byte under the testing allocator.

- [ ] **N41** `color.zig`: `RGBAf32`'s, `HSLAf32`'s and `YCbCrAf32`'s `equals` bit-cast the 128-bit color to `u32`, a compile error once anything calls them (nothing does yet). Compare the components; test colors differing only in alpha.

- [ ] **N42** `shaders/layout.zig`'s `Globals` is 24 bytes in the shaders (`V2f32`) and 20 in `gpu.zig`'s push constant (`[2]f32`): check the offsets the shaders read against what the host writes.

- [ ] **T144** `powerpc64le`, `powerpc64le-power10` (ReleaseSafe) and `powerpc64` do not compile: LLVM cannot select a `v1i128` bitcast of an `f128` load in `tridiag_evd.qrAlgorithm` (`tridiag_evd.zig:163`, issues.md). Tier 3; no workaround in rstd. **Blocked:** upstream.

- [ ] **T152** `linux-powerpc64le-power10`: six LAPACK tests fail in `Complex(f32)` at denormal or tiny scales: `getrf`'s denormal pivot (dividing by 1.1e-44 gives `inf`, `nan`), `syevd`, `gees` at scale 1e-30, `gesvd` (NoConvergence), `gels` (Singular), the evd NaN test. The power9 rig (`linux-powerpc64le`) passes them all: likely power10 code (LLVM or QEMU) flushing or overflowing in the complex division. Not traced. Tier 3.

- [ ] **T153** `linalg.matrix`'s `@Matrix` path and its column code (5f7ae56f) can differ in bits: on MIPS with MSA the fork's Sema fuses `M v`, `v M` and `M N` products as `std_float_options` says, while `fma.mulAdd` (the column code) stays unfused there (`fma.compiles`: LLVM compiles MSA's `fmadd` wrong, issues.md). Also, with both operands NaN, the column code adds `term + sum` where Sema adds `sum + term`, so the NaN payload may differ; the bits test ("mul through @Matrix has its column code's bits") uses finite values only. Fix: make the fallback follow Sema there, or keep Sema off `@mulAdd` under MSA in the fork. Tier 3 (fork compiler only).

- [ ] **T145** x87, Debug only (`i386`, `rstd-test-x87`; the ReleaseSafe forms pass): "geev on special matrices" fails `f64` n = 2 kind 3, the rotation `[0 -1; 1 0]`: an eigenvalue's `x[kbig].im` is past `tolerance(F, n, 1)` (`evd.zig:852`). Likely `f80` intermediates in the 2 x 2 standardization; not traced. Tier 3.

- [ ] **T146** x87 (`rstd-test-x87`): matlab's `rat` loses precision (22/7 comes back with numerator 3.14306640625; 3.1416015625 for π), and `base` renders the wrong digits (`matlab/rat.zig`, `matlab/base.zig`). Not traced; M6 reworks both. Tier 3.

- [ ] **T147** "the registry read works on Windows" (`cpu/detect/arm64.zig:917`) fails on the 32-bit `rstd-test-x87`: `InstallTime` is missing, likely from WOW64's 32-bit view of `HKLM\SOFTWARE` (`KEY_WOW64_64KEY` reads the native one). Matters for 32-bit arm Windows detection too. Tier 3.

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

- [ ] **T20** Kernel strategy: race `@Vector` against the converted assembly per core on a quiet machine (only vortex is done, `tools/blasbench/results/results-m6.md`); keep per-core code only where Zig loses. riscv64 races (from T22): `@Vector` against the `asm` `axpy` on zvl128b and x280, and fixed-length `@Vector` against VLEN-agnostic RVV on hardware wider than a core's minimum VLEN.

- [ ] **T05** Thresholds: OpenBLAS's defaults and per-core overrides, tuned by benchmark; then the overrides left out of phase 1 (z13 `axpy`, Windows-on-ARM64 `gemv`) if the numbers agree. **Needs:** T01, T73. T06's smoke run (busy Zen 3, `--threads=4`): threaded level-1 calls were slower than serial at the sizes it ran, so the thresholds need this tuning. Also time T07's split Accelerate calls against one whole call; and the bench's "Accelerate multithreaded" column is single-threaded while rstd's init sets `VECLIB_MAXIMUM_THREADS=1` process-wide, so a fair column runs Accelerate in a separate process. From T113's table: race AOCL-libFLAME's small-matrix AVX2/AVX-512 kernels (`getrf`, `potrf`, `potri`, `getrs`, `getri`, `geqrf`, `gelqf`, small `gesvd`; its OpenMP paths only after a same-bits check), libflame's restructured QR for `gesvd` and Hermitian eigenvalues (Van Zee et al., TOMS 2014), `FLA_Hevd`/`FLA_Tevd`, and `FLA_Eig_gest`/`FLA_Sylv` against ReLAPACK's `sygst`/`trsyl`. From T101 (src/linalg/blas/threading.md): OpenBLAS threads x86_64 `ddot`/`zdot` and arm64's `dot.c` above 10000 (arm64 with its own neoversev1 and a64fx tables) and never `sdot`, `cdot`, `dsdot` or `dot.S`; rstd threads every precision above 10000. Threaded `cgemm` at n = 128 in the arm64 VM varies 198 to 249 GFLOP/s: each call (about 75 µs of work) wakes five sleeping `Io.Threaded` pool threads (a task allocation, a mutex and a condvar signal each), whose wake time varies under HVF; no spinning server (src/linalg/blas/threading.md), so a threshold that keeps such calls serial is the fix to time.

- [ ] **T88** Time T25's races on an idle machine: portable against `@Vector` for `sum`, `sdot` in `f64` and in `f32` lanes, `nrm2` and the extremes, natively and as an x87 build (a scalar-FPU stand-in); bench and runner in the scratchpad (`vecsum/_vecsum_bench.zig`, `timeit.ps1`). Also race `vortex`'s converted `sum` against `@Vector` on the M6 (kept in `vector.displacing.keep` until then).

- [ ] **T26** vortex: a tuned `@Vector` `rot` (the assembly still leads by 2 to 17% in L1); race the Zig complex `gemm` against OpenBLAS's vortex assembly on the M6.

- [ ] **T27** gemm tuning: faster Zig packing (whole-vector transposes in `oncopy`, `incopy`), `C`-tile prefetch at kernel entry, a larger `KC` on Apple cores (src/linalg/blas/comparisons.md (9)).

- [ ] **T28** Measure whether composed `symv`/`hemv` pay off on cores with assembly gemv (atom, opteron, arm64); today they compose only where the core has its own gemv.

- [ ] **T29** The bfloat16 intrinsics (cooperlake, sapphirerapids, Neoverse N2, V1, V2, POWER8 to 10); `sbdot` could call the core's own `sdot`. **Needs:** T20.

- [ ] **T31** Assembly where `@Vector` cannot reach: arm64 SVE, SME (`vortexm4`), POWER10 MMA, AMX and AVX-512 BF16 for `sbgemm`; Skylake-X only if under 85% of its kernel (src/linalg/design.md, Implementation plan; `src/linalg/blas/blas3/gemm/notes.md`). **Needs:** T20.

- [ ] **T32** BLIS assembly kernels, decided per kernel where benchmarks favor them (BSD-3, with a license section and per-file provenance). **Needs:** T20.

- [ ] **T39** sparc64 is tier 4 (user, 2026-11-01): its images outgrew OpenBIOS's 0x400000 to 0x800000 window again (tier 3 at 5ddbd422: `blas3-tri` in ReleaseFast, 34 of 46 parts in Debug), so both rigs are stubs (`build/qemu/machines.zig`) and out of the sweeps. Bringing it back needs a real fix first: map more memory at boot through the Open Firmware client interface (claim and map), so one image fits. Open there too: the `-Oz` miscompile (issues.md, needs a reduced case) and the window itself (our observation; OpenBIOS's source unchecked).

- [ ] **T71** Linux rigs with no item: loongarch64 (the `AT_HWCAP` LSX/LASX probe), big-endian powerpc64 (`AT_PLATFORM`, the `64be` variant), 32-bit mips.

- [ ] **T41** `zrot_ppc440.S` (`ppc440fp2`) has a zero-stride patch no rig has run; no rig runs `ppc440` (it shares code with `ppcg4`). The 32-bit `ppc440` and `ppc440fp2` kernels rebuilt in 719a41b are compile-checked only.

- [ ] **T74** No CI runs a BLAS test: rstd's `.gitlab-ci.yml` is rgui's copy (`EXE_NAME: "rgui"`, Zig 0.16.0). A job for `rstd-test`, `rstd-check-macho`, `rstd-check-uefi` and `rstd-check-zig017` (on Zig 0.17.0) would catch drift between rig sweeps; optional if testing stays manual.

- [ ] **T53** macOS VM extras: base-image overlays (`qemu-clone`, `-Dmac-base`), a `Forgejo` self-hosted runner, a setup that needs no clicking, less guest load when headless (HID driver, `SecurityAgent`), and freeze detection that asks before killing.

- [ ] **T104** Optional, best effort: the reference's GCC build contracts the MSA kernels' `y0 += da_vec * x0` into `fmadd` on every MSA core, but rstd fuses vectors only on `loongson3r4` (T86). The same inline MSA `fmadd` from `fma.zig` could model `p5600`, `i6400`, `i6500` and `p6600` as the reference builds them; QEMU rounds neither way, so it rests on reading the assembly.

- [ ] **T76** Optional, from src/linalg/blas/comparisons.md: `Dot2` (compensated `f32` lanes) if T25 finds the `f64` `sdot` costly in cache; BLIS's `JR`/`IR` ways for threads sharing an L2; an A-panel prefetch beside T27's `C`-tile one; small-`n` `sgemv_n` tuning (`tools/blasbench/results/results-m6.md`).

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

## Housekeeping

Readability findings from a read-only audit of every handwritten file (2026-11-07), for the user to promote. Highest payoff first in each group. Moves that rename tests also update `build/qemu/machines.zig`'s filters; generated files change only through `tools/blas_asm`. Of 120k handwritten lines, 15k are comments, in 606 blocks over 4 lines (42 over 20).

### Across rstd

- [ ] **H113** Rest of H2: one `table`/`no_table` for BLAS's and `lapack/interface.zig`'s (not in `util.zig`, which the kernel objects compile: it would import `__root__.zig`).

- [ ] **H115** Slots by enum, part 2: a `portable.Param` enum for the gemm blocking parameter slots (`gemm_p`, `gemm_q`, `gemm_r`, `unroll_*`, `align_k`) in `gemm/params.zig`, `driver.zig`, `driver3m.zig`, `_test/blas3/check.zig`, `lapack/lu.zig` and `lapack/_test/check.zig`.

- [ ] **H116** Slots by enum, part 3: the half-precision and `gemm_small` slots in `half.zig`, `ext/bf16/kernel.zig`, `small.zig`, the half gemv sites in `blas2/interface.zig` and `gemv/__root__.zig`, and symmetric's first and second gemv names.

- [ ] **H117** Slots by enum, part 4: `gemm/kernel.zig`'s installs and its symbol exports by `Slot` name.

- [ ] **H118** Slots by enum, part 5: the ops `blas3/driver.zig` and `blas3/triangular.zig` build from routine and option strings, by an exhaustive mapping.

- [ ] **H119** Slots by enum, part 6: `_test/blas1/edge.zig`'s and `denormal.zig`'s dot, reduction and extreme loops over string lists; `rotSlot` (3 spellings) once.

- [ ] **H114** Rest of H12: one `expectSameBits` (in `_test/bitwise.zig`) for `vector.zig`'s and `matrix.zig`'s.

### BLAS core: `blas/`, `cpu/`, `kernels/`, `ext/`

- [ ] **H147** Rest of H22: `vector.zig`'s `RealSlots` and `ComplexSlots` (150 lines) derive from `portable.slot`; `survey.zig` reads them as text, so it changes too.

- [ ] **H148** Rest of H28: core names as enums, not strings, in `gemm/params.zig`'s `ruleFor`, `gemv/driver.zig`'s `nbmax` and `cform`, the `entryOf(core: []const u8)` functions and `uk_block.find`'s core; `rule()`'s string API (fed by `blas.coreName()`) changes with them.

### BLAS levels 1 and 2

- [ ] **H120** Rest of H36: gemv's half-precision tests to `ext/bf16`.

- [ ] **H38** One set of level-1 wrappers and arithmetic for level 2: `triangular.zig`'s `Ops(T)` and `Num(T)` replace `banded.zig`'s `axpyK` to `scalK` and `mul` to `scaleF`, and `blas2/interface.zig`'s `scaleY`, `gemvCall`, `mvCall`. Set aside 2026-11-07: not provably bit-identical (the two reach the vtable differently and round complex products differently); needs a build and the tests.

- [ ] **H40** One `@Vector` axpy block for `axpy/__root__.zig`'s `SimdBlock`, x86_64's `VectorBlock` and `power/c_kernels.zig`'s `SaxpyVec` (also `wasm/simd.zig`'s strided loop, msa's and wasm's `load`/`store`); check the bits on power. Set aside 2026-11-07: the merged block is not provably the same operations in the same order; needs a build.

- [ ] **H121** Rest of H41: `axpy/{arm,loongarch64,x86,sparc}/kernel.zig` (registry-only wrappers) fold into `kernels/<arch>.zig`.

- [ ] **H122** Rest of H42: `symv/driver.zig`'s one-column tail (twice) and 4-column update (3 times), `gather`/`scatter` against gemv's `small` (floating-point kernel code: check the bits).

- [ ] **H123** Rest of H42: the x86_64 `kernel.zig` feature constants (5 copies) through `util.x86Has`; their unguarded `featureSetHas` gives another value off x86, where the tests also analyze them.

- [ ] **H124** Rest of H43: `blas2/interface.zig`'s flag table by `Trans`'s methods.

- [ ] **H125** Rest of H45: `checkAll` through `check.keep` in every root (changes the tests' flow), one portable-vtable setup in tests (one is static, one on the stack), `blas1/__root__.zig`'s routine list once (by reflection, which changes the installation order).

### BLAS level 3

- [ ] **H49** One driver loop skeleton (js, ls, is) for `driver.zig`, `driver3m.zig` and `half.zig`; verify the bits with the folded-beta and threaded tests. Set aside 2026-11-07: the order of operations needs the tests to confirm.

- [ ] **H53** Kernel helpers once: `gemm/params.zig`'s cache detection (from line 125, with private `cpuid`, `cpucfg`) to `cpu/detect/`; 3M's `Sym` is symm's `Source`; one `sideKernel` for trmm and trsm; one `Source` copy driver; one `TileSource` installer; 3M's kernel loops on the real ones through a store function (lower priority). Set aside 2026-11-07: kernel-object code, needs a build.

- [ ] **H54** One test `Routine` enum (`check.Routine`) and dispatcher for `blas3/driver.zig`'s (504) and `_test/blas3/threaded.zig`'s (223), which swap syrk's `opt`. Set aside 2026-11-07: the two pass different options (syrk's order, trmm and trsm's `conj`, herk's real `alpha`), so one dispatcher changes which cases run.

### LAPACK

- [ ] **H65** `hseqr`, `gees` and `Eigenvalues` move from `evd.zig` into the merged Schur file, ending the schur-evd import cycle. **Needs:** Z6.

- [ ] **H126** Rest of H67: `qr.zig`'s `randomRank` (built by `gemm`: other bits) and the svd and evd special-matrix tests' inline rank-r products onto `rankProduct`.

### linalg public layer, matlab and fft

### rgl

- [ ] **H78** One comptime field reader for sfnt tables replaces 256 hand accessors (`os2.zig`, `head.zig`, `hhea.zig`, `maxp.zig`, the cmap formats). Set aside 2026-11-07: needs a generic big-endian reader over every field type and about 256 call sites checked by a compile; rgui calls `head.unitsPerEm()`, so the accessors stay.

- [ ] **H81** Generic cmap formats: `GroupFormat` for fmt8, 12, 13 and `TrimmedFormat` for fmt0, 6, 10 (9 files to 5); one record binary search for kern, mvar, cmap, CFF and GPOS (6 copies). Set aside 2026-11-07: the formats differ in header layout, accessor names and search, and a shared search needs a context struct per call site, longer than the loops (CFF's FDSelect is an upper-bound search).

- [ ] **H127** Rest of H83: rgui to `font.render`, `font.msdf` and `rasterize`'s glyph bitmaps; then drop `Sfnt`'s `glyphBitmap*` and `font.sfnt.render`/`msdf` aliases, ending the Sfnt-render import cycle.

- [ ] **H128** Rest of H86: `gen_icon.zig` (a build tool's `main`) to `tools/`; rstd's package `.paths` lacks `tools/`, and rgui builds it through the dependency.

- [ ] **H129** Rest of H87: 2-D points on `linalg.vector` (`geometry.zig`'s `PointF32` math, the hand dot and length in sw and vulkan, glyf's `@reduce`); `vector.dot` and `length` may fuse and fall back to hypot, so the bits and `PointF32`'s API change.

- [ ] **H130** Rest of H88: one CFF integer decoder for `Buf` and `Dict.zig`; they treat malformed input differently (`Buf` reads 0 past the end and panics on a bad byte, `Dict` the reverse), so pick one.

- [ ] **H131** Rest of H89: atlas's duplicate placement type, the mirrored `MsdfOptions` and a split of the 190-line `buildTextAtlas`. **Needs:** N38.

- [ ] **H132** Rest of H89: one rectangle type for `geometry.Rect` (`u32` sizes) and `pack.Rect` (`i32`).

- [ ] **H133** Rest of H91: the `rgl.zig` shader becomes `shapes.zig` (`build/__root__.zig`'s shader-probes loop takes one string for the file and the output name).

- [ ] **H134** Rest of H91: one of `create`/`init` (`Generator.create`, `pack.Context.create`, called in `buildTextAtlas`). **Needs:** N38.

### Tests and rigs

- [ ] **H94** `_test/blas1/denormal.zig` folds into `edge.zig` (about 150 duplicated lines), as blas2's does. Set aside 2026-11-07: the fold would run edge's checks, sizes and gap checks on denormal operands, which changes the tests.

- [ ] **H96** One failure collector (`check.keep`) for `edge.zig`'s and `denormal.zig`'s `note` and `blas3/check.zig`'s `Bad`; one rig predicate (`ob.reduced` or `check.rig`); dedupe inside `blas3/check.zig`. Set aside 2026-11-07: `note` and `Bad` print and cap differently; `reduced()` and `rig()` differ on macOS and `.other`.

### Tools, build and the rest of `src/`

- [ ] **H149** Rest of H103: one per-arch table for `gen.zig`, `uk.zig`, `survey.zig` and `params.zig`; one getarch reader (`params.zig` preprocesses `getarch.c`, `openblas.zig` reads its `FORCE_` blocks as text); `coreDefines` as data; `uk.emitAsm` in parts.

- [ ] **H135** Rest of H106: a contestant and row helper, `Call.of`, an `EnumArray` for the library paths, `std.math.Complex` (may change the operation order); check with `bench-check`.

- [ ] **H136** Rest of H107: `tools/rig/` and `tools/macos/macos_vm.zig`, with `tools/remote.zig` as a named module (a file in a subfolder cannot import it).

- [ ] **H137** Rest of H108: an options struct for `addRstdModule`'s two bools.

## Done

Cleared 2026-11-02; landed work is in `review`'s log.

- [X] **N52** The Mac ssh tools keep the Mac awake for their runs (4c9fbba4, 2026-11-08): `ssh_test.zig`'s remote command and `test_host.zig`'s `nix shell … run.sh` run under `caffeinate -i` on a Darwin host, since idle sleep ignores ssh sessions without a tty.

- [X] **N51** Mac hosts save power while staying reachable (eba7ac35, 2026-11-08): both host scripts sleep the display and disks after 10 minutes and the Mac after 60 (Intel also `hibernatemode 0`) with the screen lock off (auto-login covers boot only); an incoming ssh packet dark-wakes either Mac in about a second. Rerun the scripts on the Macs once N52 lands.

- [X] **N50** `configure_mac_host_intel.sh` disables the Spotlight UI agent too (4583eee7, 2026-11-07): with corespotlightd off, Sequoia's `com.apple.Spotlight` retried its IndexAgent endlessly (85% of a core, 40% more in launchd).

- [X] **N43** Intel Mac host scripts (c98a815e, 2026-11-07): `configure_macos.sh` is now `configure_macos_vm.sh`, `configure_mac_host.sh` is `configure_mac_host_arm.sh`, and `configure_mac_host_intel.sh` set up `macmini-intel`.

- [X] **N44** `tools/configure_mac_nix.sh` (df12e3a8, 2026-11-07): Nix and nix-darwin as on the M6; run on the Intel Mac, with nixos-cfg's `macmini-intel` flake on nixpkgs 26.05.

- [X] **H1** One set of BLAS scalar helpers in `util.zig` (`raw`, `isZero`, `isOne`, `zero`, `one`, `first`; one meaning each for `at` and `off`), `apple_accel.zig`'s and `cholesky.zig`'s copies too (`housekeeping`, 2026-11-07).

- [X] **H2** One set of comptime type checks in `util.zig` (`requireFloat`, `requireComplex`, `requireReal`; lapack's alias them); `linalg.vector`'s `isComplex` and `Real` alias `types.zig`'s. Rest: H113 (`housekeeping`, 2026-11-07).

- [X] **H3** One element-arithmetic set in `lapack/ops.zig` for lapack, `matrix.zig`, `dense.zig` and matlab's generators; `triangular.zig`'s `recip` (`ztrti2_U.c`'s ratio) and matlab's `recip` (`scalar.wdiv`) differ and stay (`housekeeping`, 2026-11-07).

- [X] **H4** One `ops.narrow(E, e)` for the error-narrowing copies in lapack, `dense.zig` and matlab (`housekeeping`, 2026-11-07).

- [X] **H5** One stack-or-heap `Scratch` in `util.zig` for levels 2 and 3 (`housekeeping`, 2026-11-07).

- [X] **H6** `threading.forEachPart` for the per-chunk scratch pattern (`housekeeping`, 2026-11-07).

- [X] **H7** `util.scalarArgs` for the two-armed complex kernel calls; two test calls whose arguments are casts stay (`housekeeping`, 2026-11-07).

- [X] **H9** One shared test module `_test/check.zig` beside `_test/guard.zig` (`use`, `restore`, `Poison`), `_test/reference.zig` (the `f128` reference, ulp bounds), `_test/pattern.zig` and `_test/bitwise.zig`; lapack's tests install vtables by `guard.use` and `guard.restore` (`housekeeping`, 2026-11-07).

- [X] **H10** lapack's test helpers and `thresh`, `cplx`, `ratio` in `lapack/_test/check.zig` (`housekeeping`, 2026-11-07).

- [X] **H11** lapack's shared test fixtures in `lapack/_test/fixtures.zig` (`CholProblem`, `LuProblem`, `Tri`, `schur`) and `_test/symmetric.zig` (`Indefinite`, `ldlProduct`, `sandwich`) (`housekeeping`, 2026-11-07).

- [X] **H12** Shared test helpers in `_test/check.zig`; blas3's `Wide` is `Ref`; `dense.zig` takes its `sameBits`. The per-element `sameBits` (`batch.zig`, `half.zig`: they read `f80` differently) and the kernel objects' `fillSeq` stay. Rest: H114 (and `rotSlot`, H119) (`housekeeping`, 2026-11-07).

- [X] **H13** One `WorkerRuns` harness and `eachTestIo` in `_test/io.zig` for qr, lq, lsq, svd, evd and blas3's threaded tests (`housekeeping`, 2026-11-07).

- [X] **H14** Named ulp bounds in `_test/reference.zig` (`housekeeping`, 2026-11-07).

- [X] **H15** Every import at the top of its file; the unused ones dropped (`housekeeping`, 2026-11-07).

- [X] **H16** Semicolons for the comment dashes (minus signs, list bullets and `// ---` separators kept) (`housekeeping`, 2026-11-07).

- [X] **H17** The listed Markdown unwrapped (`housekeeping`, 2026-11-07).

- [X] **H18** One `cpu.Arch` table replaces the 12 `cpu/<arch>.zig` stubs; `select.zig`, `detect/`, `coretypes.zig` and `corespec.zig` lose their underscores (`_auxv.zig` is H27's `osinfo.zig`) (`housekeeping`, 2026-11-07).

- [X] **H19** `registry.Arch` builds each `kernels/<arch>.zig`'s registry (`housekeeping`, 2026-11-07).

- [X] **H20** `coretypes.Family`, `familyOf` and `coretypes.all` for the arch switches and hand lists; `select.archOf` maps the other way and stays (`housekeeping`, 2026-11-07).

- [X] **H21** One `detect` and `canRun` in `select.zig` over each detector's probe; mips's own is `runsHere` (`housekeeping`, 2026-11-07).

- [X] **H22** vector.zig's slot placement (`Rule`, `takes`, `Registry` as `Vector`) is in the registry file; `routed` and `peels` stay, as the coverage pages cite them. Rest: H147 (`housekeeping`, 2026-11-07).

- [X] **H23** `asm_registry.zig` is `registry.zig`, with vector's placement, `Union` for the asum, gemv and symv merges (axpy's power and mips64 ones are not plain unions: the ports win slot by slot), and `uk_block.zig` beside it in `blas/` (`housekeeping`, 2026-11-07).

- [X] **H24** `FpState` in `FpState.zig`; its users import it (`housekeeping`, 2026-11-07).

- [X] **H25** `apple_accel.zig`'s tests in `_test/apple_accel.zig` (hosted only: the "routed calls split" test, which only skipped there, no longer builds on the SMP rigs), its quirk list in design.md's Accelerate section, `portable.letter`; `span` and `mat` compute extents, not `util.fits`'s overflow check, and stay (`housekeeping`, 2026-11-07).

- [X] **H26** Threading's test infrastructure in `_test/io.zig` (`Evented` stays: `usesFutures` needs it); its users import it (`housekeeping`, 2026-11-07).

- [X] **H27** `_auxv.zig` is `cpu/osinfo.zig`, with arm64's registry and sysfs reads (`housekeeping`, 2026-11-07).

- [X] **H28** Core names as enums in threading's per-core thresholds and in the generated uk and skip lists. Rest: H148 (`housekeeping`, 2026-11-07).

- [X] **H29** `coretypes.staticOnly` lists and one `select.testing` for the test doubles; x86's 32-bit-only prongs are not that list (`housekeeping`, 2026-11-07).

- [X] **H30** Half widens through `bf16.widen`; one `vector.splatTo` for the 4 copies (`housekeeping`, 2026-11-07).

- [X] **H31** Dead code dropped: `cpu.unsupported`, x86's `Tier`, `vendor` and `displayFamily`, `auxv.has`, the bf16 aliases. User's call: H138 (`housekeeping`, 2026-11-07).

- [X] **H32** The stale BLAS core comments fixed (`housekeeping`, 2026-11-07).

- [X] **H33** The BLAS core's long comment blocks cut; Vtable's mirrored C `#define` blocks kept, as derived source (`housekeeping`, 2026-11-07).

- [X] **H34** One `blas1/extremum.zig` table for the 8 extremum routines; their folders keep their coverage pages (`housekeeping`, 2026-11-07).

- [X] **H35** `blas1/interface.zig`'s threading tests and special-value tests moved out; `rotg`, `zrotg`, `rotmg` in `blas1/rot/rotg.zig` (`housekeeping`, 2026-11-07).

- [X] **H36** gemv's x87 extended sums and their tests in `gemv/extended.zig`. Rest: H120 (`housekeeping`, 2026-11-07).

- [X] **H37** The rank updates in `blas2/rank.zig`; banded's `Kind` is symmetric's (`housekeeping`, 2026-11-07).

- [X] **H39** `symmetric.zig` takes its types from `types.zig`; `copy_cases` in `_test/blas2/` (`housekeeping`, 2026-11-07).

- [X] **H41** The axpy arch headers are one line each. Rest: H121 (`housekeeping`, 2026-11-07).

- [X] **H42** `util.x86Has` for the gemv and rot drivers' x86 checks. Rest: H122, H123 (`housekeeping`, 2026-11-07).

- [X] **H43** gemv's flag table and `cform` use `Trans`'s methods; `Extreme.op` and `indexOp` name the slots. Rest: H124 (`housekeeping`, 2026-11-07).

- [X] **H44** The routine-level `kernel.zig` drivers are `driver.zig` (`housekeeping`, 2026-11-07).

- [X] **H45** axpy's test names as its siblings', scal's flag named as the portable kernels'. Rest: H125 (`housekeeping`, 2026-11-07).

- [X] **H46** The level-1 and level-2 module docs and long function comments cut (`housekeeping`, 2026-11-07).

- [X] **H47** One gemm and one gemmt argument check, public in `blas3/interface.zig` (`housekeeping`, 2026-11-07).

- [X] **H48** `Blocking.fromSlots` and named `default_p`, `default_q`, `default_r` (`housekeeping`, 2026-11-07).

- [X] **H50** The trmm and trsm drivers and their threaded split in `blas3/triangular.zig` (`housekeeping`, 2026-11-07).

- [X] **H51** `blas3/driver3m.zig`; the kernel registry in `blas3/kernels.zig`; `gemm/small_x86_64.zig`; gemmt's tests in `_test/blas3/`; `blas3/__root__.zig`'s order and doc. The routine folders (hemm, hemm3m, syr2k) stay (`housekeeping`, 2026-11-07).

- [X] **H52** `blas3/threaded.zig`'s tests in `_test/blas3/threaded.zig`. User's call: H139 (`housekeeping`, 2026-11-07).

- [X] **H55** Level-3 small fixes: headers and imports, `batch.zig`'s stale T30, `hv`, 3M's `Desc` takes `Conj`, `inline else` tri dispatch, `small.gemm`, one masked load and store (`housekeeping`, 2026-11-07).

- [X] **H56** Level-3 headers cut to 3 lines, rationale in `gemm/notes.md` (`housekeeping`, 2026-11-07).

- [X] **H57** `cholesky.zig`'s shared argument, vtable and thread helpers in `lapack/interface.zig` (`housekeeping`, 2026-11-07).

- [X] **H58** Real-slice helpers in `ops.zig`; one `lascl` in `auxiliary.zig`; one driver scaling prologue (`housekeeping`, 2026-11-07).

- [X] **H59** BLAS wrappers in `lapack/blas_calls.zig`; `lartg`, `larfg`, `lapy2`, `lapy3` and `laset` in `auxiliary.zig` (`housekeeping`, 2026-11-07).

- [X] **H60** `trsyl` and `tgsyl` with their tests in `trsyl.zig` and `tgsyl.zig`, out of `relapack.zig` (`housekeeping`, 2026-11-07).

- [X] **H61** `sytrf`, `hetrf`, their rook forms and tests in `sytrf.zig`, out of `relapack.zig` (`housekeeping`, 2026-11-07).

- [X] **H62** `gbtrf` and its tests in `gbtrf.zig`; relapack's other tests in `lapack/_test/relapack.zig` (`housekeeping`, 2026-11-07).

- [X] **H63** evd's tests in `lapack/_test/evd.zig`; trevc3's, `sqrt`'s and `rotg`'s beside their code (`housekeeping`, 2026-11-07).

- [X] **H64** One `lasy2` in `sylvester.zig` (as `dlasy2.f`) with `laln2` beside it; `sequenceFactors` in `householder.zig` (`housekeeping`, 2026-11-07).

- [X] **H66** One `eachBlocking` test helper for cholesky and lu, one NaN-column check; cholesky's loops now run type, vtable, triangle, with the same cases (`housekeeping`, 2026-11-07).

- [X] **H67** Shared test generators `rankProduct`, `hermitian` and `expectNear`. Rest: H126 (`housekeeping`, 2026-11-07).

- [X] **H68** latrs's rescaled divide once, `ops.Mat` for `at`, `cabs2`; `norm.Op` is `types.Trans` (`housekeeping`, 2026-11-07).

- [X] **H69** lapack small fixes: root order and aliases, stale docs, dead code, test names; the sparc64 `lapack-faer` filter lists the 12 former `faer_*` modules (it ran none before). User's call: H140 (`housekeeping`, 2026-11-07).

- [X] **H70** lapack's module headers and long docs cut, rationale in design.md (`housekeeping`, 2026-11-07).

- [X] **H71** `dense.zig` public first, one fixed-size dispatcher; `matrix.zig` too (`housekeeping`, 2026-11-07).

- [X] **H72** Not done (user, 2026-11-07): matlab's files stay one per unit, so the merges into `lu.zig`, `eig.zig`, `scalar.zig` and an `operators.zig` were undone (`housekeeping`, 2026-11-07).

- [X] **H73** matlab dedupe: `triu`, one binary loop, `Unary`, one `inline for`. `chol`, `schur` and `hess` of an empty matrix and `hess` and `schur` of a 1x1 return the empty or trivial factors (`chol([])` panicked in `potrf`'s `lda` check before), with tests (`housekeeping`, 2026-11-07).

- [X] **H74** matlab names: `Cx` only, one `Allocator`, `lu.lange1`, `svd.zig`, `vector.mix`'s docs. User's call: H141, H142 (`housekeeping`, 2026-11-07).

- [X] **H75** `matlab/readme.md`'s module table matches the files (`housekeeping`, 2026-11-07).

- [X] **H76** matlab's headers and long docs cut, with `dense.zig`'s, `vector.zig`'s and `fft`'s; HESS and SCHUR cite `MATFN2`, NORM `MATFN3` (`housekeeping`, 2026-11-07).

- [X] **H77** One synthetic test font in `render/_test/font.zig` (`housekeeping`, 2026-11-07).

- [X] **H79** GPOS kerning on the `Lookup`, `PairPos`, `Coverage` and `ClassDef` views; format-1 ClassDef's `u16` offset kept as it was (N39) (`housekeeping`, 2026-11-07).

- [X] **H80** `shaders/layout.zig`: `Globals(V)`, `ShapeKind` and the instance structs, for the shaders and `gpu.zig` (sizes differ: N42) (`housekeeping`, 2026-11-07).

- [X] **H82** glyf's files deduplicated (`appendVertex`, one coordinate decoder, `table_data`); the `tables/glyf/` folder stays (`housekeeping`, 2026-11-07).

- [X] **H83** msdf and render up to `font/`; `Generator.zig`; `glyphBitmap*` in `rasterize.zig`. Rest: H127 (`housekeeping`, 2026-11-07).

- [X] **H84** `stream_types.zig` is `sfnt/types.zig`, with `FixedPoint` for `ShortFrac` and `Fixed`; `aux_types.zig` stays (`housekeeping`, 2026-11-07).

- [X] **H85** msdf coloring's 0- and 1-corner cases shared; HVAR, VVAR and cvar stay their own table files (`housekeeping`, 2026-11-07).

- [X] **H86** render's backends by `inline else`; `sw/` is `sw.zig`; `Image.downsample`; `draw.zig` stays. Rest: H128 (`housekeeping`, 2026-11-07).

- [X] **H87** One `median3`, shaderlib's, for msdf too. Rest: H129 (`housekeeping`, 2026-11-07).

- [X] **H88** The charstring interpreter switches on `Instruction`, with the escape ops apart; `CharstringCtx` private. Rest: H130 (`housekeeping`, 2026-11-07).

- [X] **H89** An options struct for `rasterize`; `pack.zig` up to `rgl/`; `tessellate.zig` and `bounds.zig` stay. Rest: H131, H132 (`housekeeping`, 2026-11-07).

- [X] **H90** rgl's dead code dropped: `textures.zig`, `sdRoundedBoxMultiCorner`, `CharstringCtx.deinit`, `Version16Dot16.format`. User's call: H143; `color.equals`: N41 (`housekeeping`, 2026-11-07).

- [X] **H91** `deinit` everywhere, one `Allocator`, `stepScalar`, PNG's names and split `decode`, `Sfnt.init`'s prints behind `InitOptions.debug`. Rest: H133, H134; user's call: H144 (`housekeeping`, 2026-11-07).

- [X] **H92** The sfnt readme's files and API current; `Sfnt.zig`'s docs trimmed (`housekeeping`, 2026-11-07).

- [X] **H93** rgl's table-view note said once in `tables/__root__.zig`; long comment blocks cut (`housekeeping`, 2026-11-07).

- [X] **H95** Shared `Result`, `bump` and letter parsing for blat2 and blat3; `makeVecZeroMid`. Their `Sweep` and `dat` differ, and stay (`housekeeping`, 2026-11-07).

- [X] **H97** `from_openblas`'s tests named for the suites they run; the level-1 wrappers inlined (`housekeeping`, 2026-11-07).

- [X] **H98** `test/qemu/pc.zig`, `uhi.zig` and `ofw.zig`; ppc32's PAPR constants named (`housekeeping`, 2026-11-07).

- [X] **H99** The `_test/` and `test/` headers cut to 2 or 3 lines, rationale in design.md (`housekeeping`, 2026-11-07).

- [X] **H100** The generated asm files hold only `keys` (core by enum literal, variant, a merged file's parts) and `Body`, with one-line headers; `registry.Entry` answers `serves` and `kernel`; uk.zig's headers on one line too (`housekeeping`, 2026-11-07).

- [X] **H101** One gemv uk index; the generated kernels stay under `<arch>/asm/`, apart from the handwritten `kernel.zig` (`housekeeping`, 2026-11-07).

- [X] **H102** `params_data.zig` without default fields (read from `params.zig`), sizes as shifts, POWER's other byte order as its differences (`ppc970`'s); one unwrapped legend, `blas/coverage.md`, linked from each routine's page (`housekeeping`, 2026-11-07).

- [X] **H103** blas_asm: `openblas.zig` (was `kernel_map.zig`) with `initFromEnv`, one `isName`, `map.format`; one `Patch.apply`; uk's per-target state on `Target`; `gen.convert` in parts. Rest: H149 (`housekeeping`, 2026-11-07).

- [X] **H104** The rigs in `build/qemu/` (`machines.zig`, `macos.zig`, the steps in `__root__.zig`) and the benches in `build/bench.zig`, ending the import cycle; one initramfs packer, `Boot.Linux.commandLine`, `Machine.semihosting()`, `tools.mac_alias` (`housekeeping`, 2026-11-07).

- [X] **H105** `tools.hostTool` for the six tool builders; `tools/remote.zig` for the ssh and spawn code (`housekeeping`, 2026-11-07).

- [X] **H106** `blas_bench.zig` split into it, `kernels`, `libs` and `ext`; one `checkGemm`, `fillRandom`. Rest: H135 (`housekeeping`, 2026-11-07).

- [X] **H107** `tools/probe/`, `tools/macos/`, `tools/blasbench/results/`. Rest: H136; user's call: H145 (`housekeeping`, 2026-11-07).

- [X] **H108** `_test_wasm.zig` gone; `machine.zig` is `bootinfo.zig` (RandyOS's edits staged on its `main`); `buildutils/` keeps its 16 files, with one docs-object helper in `docs.zig`; `_os/` keeps its name: platform code behind the public wrappers (user). Rest: H137 (`housekeeping`, 2026-11-07).

- [X] **H109** The UEFI `Io` vtable built by a loop over its parts; `stubs.zig`; one console and read loop; one `getTimeSeconds`; the stale comments (`housekeeping`, 2026-11-07).

- [X] **H110** Dead code dropped: `needsLlvm`, `patientOptions`' parameter, the unused consts. User's call: H146 (`housekeeping`, 2026-11-07).

- [X] **H111** The tool and build headers cut, usage in the readmes (`housekeeping`, 2026-11-07).

- [X] **H112** The bridge: `portable.Prec`, `slotName`, `slotPtr`, `slotTag`, `gemvOp`; the symv, hemv, ger and gemv installs and tests, level 1's simple edge and denormal sites, ext's matcopy, neg and laswp (one `matcopy.Copy`) and symmetric's copy kernel by op. Rest: H115 to H119 (`housekeeping`, 2026-11-07).

- [X] **N36** msdf `Shape.orientContours` orients contours (2026-11-07): it skips a contour once it has an orientation, as msdfgen's `!orientations[i]` (pin `1c106ed`), and its `y1` searches stop at the first height off `y0`. Tests: a counterclockwise outline and a hole wound like its outline are reversed; a contour closing on a horizontal edge at `y0` is scanned below it.

- [X] **H8** Kernels by slot only (2202fc46, 2026-11-07): the 844 unused per-slot wrappers (`pub fn sasum_k(...) { return blas.kernel(.sasum_k)(...) }`) are gone from the blas1, blas2, blas3 and ext `__root__.zig` files, with their repeated "Argument names are copied" header and the type aliases left unused; the roots that held only wrappers (`blas3/hemm`, `ext/omatcopy`, `ext/imatcopy`) are deleted, and the docs no longer call them vtable wrappers.

- [X] **M10** The public layer's step 3, first part (6648ae73, 2026-11-02): `View`, `ConstView` and unmanaged `Mat(T)` in `matrix.zig`, with `view()`/`constView()`; `dense.zig`'s `dot`, `norm`, `sum`, `max`, `min`, `axpy`, `scale`, `mul`/`mulInto`, `transpose`, `ctranspose` over views, flat in `rstd.linalg`. Up to `L` = 8 they give a `StaticMat`'s bits; past it, BLAS. Complex `sum` adds the parts apart; `max`/`min` real only; `scale` inline at every size (design.md).

- [X] **M17** The public layer's step 3, second part (afde537a, 2026-11-07): `solve`/`solveInto`, `lu`, `chol`, `qr`, `inv`/`invInto`, `det` and `matrixNorm` for views. Up to `L`, `det`, `inv` and square `solve` give a `StaticMat`'s bits; past it, `solve` routes by structure (`trtrs`, `potrf` falling back to `gesv`), rectangular `solve` is `gelsy`, `det` rescales by powers of 2 (design.md).

- [X] **M18** The public layer's step 3, third part (8126600c, 2026-11-07): `svd` (`.values` by `gesvd`, with vectors by `gesdd`), `eig` (a union: Hermitian by `syevd`, else `geev`), `pinv`, `rank` and `cond` for views, on LAPACK at every size; a NaN gives `NoConvergence`. `View.col` no longer slices past an empty view's end.

- [X] **M11** The public layer's step 4 (9eda96a6, 2026-11-07): `rstd.linalg.Options` (`rstd_linalg_options.inline_limit`, default 8, at least 4); a `StaticMat`'s `mul`, `det`, `inv` and `solve` past `L` are compile errors naming `rstd.linalg`'s forms (none under SPIR-V), with the unchecked ones in `linalg.matrix.any_size`; both paths tested at every size 1 to 16 against a wider type.

- [X] **M12** The public layer's step 5 (a1d588fe, 2026-11-07): under SPIR-V, `linalg.vector`'s `dot`, `dot2`, `length`, `distance`, `norm`, `cross`, `mix` and `mulAdd` are `OpDot` and GLSL.std.450's instructions, and `linalg.matrix`'s `det`, `inv` and `outer` are `Determinant`, `MatrixInverse` and `OpOuterProduct` (inline asm on `@Matrix` operands); a shader probe each, all `spirv-val`-clean.

- [X] **T157** Cheaper comptime in `linalg.vector`'s inline forms (2026-11-07): `partDot` takes its real type and complexity from the caller, `dot` uses builtins for its types, `vector.mulAdd` aliases `fma.mulAdd`, whose policy is a table by float width. msdf's `signedDistance` used 985 of its 1000 comptime branches, now 462.

- [X] **T155** `zig fmt`: a test line's indent in `balance.zig` and `arm64.zig`'s alignment (699a8e71, 2026-11-02).

- [X] **T156** `zig fmt`: `gemm3m/driver.zig`'s `@intFromEnum` to `@backingInt` (1e19eb12, 2026-11-07); `zig fmt --check src/linalg` is clean.

- [X] **T154** "Skylake-X small kernels: upstream's blocks, lanes and reductions, on any x86_64" failed on the Mac under Rosetta (f2259e1e, 2026-11-02): not Rosetta. Since ff7255f9 the kernel fuses only with an FMA instruction, which baseline x86_64 (the macOS test target) lacks, while the test's model always fused (`@mulAdd`); it now fuses as `fma.madd` does. Passes on Zen 3 at `-Dcpu=x86_64` and native, and on `ssh-test-x86_64-macos` and `-avx2`.

## Waiting on the user

- [ ] **H138** From H31, keep or delete: `Config.zig`'s unread `thread_timeout`, `goto_debug`, `main_free` and `omp_*` fields; `corespec.undecided`; the commented-out lines at `blas/__root__.zig:24` and `cpu/__root__.zig:18` (which names the gone `_auxv.zig`).

- [ ] **H139** From H52: move the `l3thread bench` (now in `_test/blas3/threaded.zig`) to `tools/blasbench`?

- [ ] **H140** From H69: one trevc normalizer? `normalizeColumn` returns early on a zero column, while `normalize` would scale it by infinity, so merging them changes behavior.

- [ ] **H141** From H74: rename `norm.matNorm1` (left only re-documented, as `matNormInf` mirrors its name)?

- [ ] **H142** From H74: make `vector.mix` `inline` like its siblings? It would raise msdf's comptime branch use (T157).

- [ ] **H143** From H90, keep as a reference note or delete: `shaders/lib/sdf.zig`'s commented-out GLSL, the old commented usage examples in `font/__root__.zig` and `font/sfnt/__root__.zig`, and the commented GLSL to-do list in `shaders/lib/math.zig`.

- [ ] **H144** From H91: should `Sfnt`'s new `InitOptions.debug` default to off? It defaults to on, so `Sfnt.init` prints to stderr as before.

- [ ] **H145** From H107: delete `tools/blasbench/axpy_bench.zig` (it cannot build) and `bench.zig`'s unused half?

- [ ] **H146** From H110, keep or change: `zigtime.zig`'s `@panic("OOM")`, the unused `bytes.zig`, `constants.zig` and `debug.assertMsg` exports, the commented-out code in `renderdoc.zig`, `spirv.zig`, `icon.zig`, `build.zig` and `build.zig.zon`, and `zigstd.zig`'s aliases (127 users).
