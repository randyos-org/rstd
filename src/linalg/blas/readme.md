# BLAS

Rstd's BLAS is a Zig port of OpenBLAS: its interfaces, drivers, per-core kernel selection and
parameters, with OpenBLAS's per-core assembly converted to Zig and its generic kernels rewritten
on `@Vector`. One build carries kernels for every core the target can run and picks one at run
time, as OpenBLAS's `DYNAMIC_ARCH` builds do.

It is a Zig library, not a replacement `libopenblas`: there are no CBLAS (`cblas_*`) or
Fortran-named (`dgemv_`) entry points, no row-major mapping layer and no `xerbla` callback. The
routines are Zig functions over column-major data in the Fortran argument order; bad arguments
return named errors.

How the code is structured, and why, is in [`design.md`](design.md).

## Supported

### Element types

| Type | Letters | Notes |
| --- | --- | --- |
| `f32`, `f64` | `s`, `d` | |
| `types.XDouble` | `q` | x87 `f80` on x86 and x86_64 outside Windows and UEFI, else `f64` |
| `std.math.Complex` of those | `c`, `z`, `x` | |
| any other float type (`f16`, `f80`, `f128`) and its `Complex` | | portable kernels built for the type |
| `types.BF16`, `types.HF16` | `sb`, `b`, `sh` | the routines under [bfloat16 and half precision](#bfloat16-and-half-precision) |

### Routines

- **Level 1:** `copy`, `swap`, `scal`, `dot`, `dotc`, `dsdot`, `sdsdot`, `nrm2`, `asum`, `sum`,
  the index and value extremes (`i?amax`, `i?amin`, `i?max`, `i?min`, `?amax`, `?amin`, `?max`,
  `?min`), `rot`, `rotm`, `rotg`, `zrotg`, `rotmg`, `axpy`, `axpyc`, `axpby`, `cabs1`.
- **Level 2:** `gemv`, `ger`, `gerc`, `symv`, `hemv`, `trmv`, `trsv`, `tbmv`, `tbsv`, `tpmv`,
  `tpsv`, `gbmv`, `sbmv`, `hbmv`, `spmv`, `hpmv`, `spr`, `spr2`, `syr`, `syr2`, `hpr`, `hpr2`,
  `her`, `her2`, including OpenBLAS's complex symmetric forms and `hemv`'s and `hbmv`'s
  transposed (`HEMVREV`) triangles.
- **Level 3:** `gemm`, `gemmt`, `symm`, `hemm`, `syrk`, `herk`, `syr2k`, `her2k`, `trmm`,
  `trsm`, real and complex.
- **bfloat16 and half:** `sbdot`, `bscal`, the conversions, and `bgemv`, `sbgemv`, `shgemv`.
- **Extensions:** `geadd`, `omatcopy` and `imatcopy` (column-major), and LAPACK's `laswp`.
- **LAPACK** (`rstd.linalg.lapack`): `potrs`, and `potrf`, which returns `error.NotPositiveDefinite` for a
  leading minor that is not positive definite, a NaN pivot included, and writes its order to
  its last argument; `getrs`, and `getrf`, which completes the factorization, then returns
  `error.Singular` for an exactly zero `U(i, i)` and writes the first `i` to its last argument.

Not yet: `gemm3m`, `symm3m`, `hemm3m`, the small-matrix kernels, bfloat16 level 3 (`bgemm`,
`sbgemm`), and threading in levels 1 and 2 ([`threading.md`](threading.md)).

### Targets

Every target compiles, with the portable kernels at least. Per-core kernels exist for these
cores; which OpenBLAS source each core uses, and for level 1 and level 3 but `gemm` what rstd
serves in its place, is in the generated `coverage.md` files (every `blas1/<routine>/` and
`blas3/<routine>/`, and `blas2/`), and which routines each architecture covers is in
[`design.md`](design.md#coverage-by-architecture).

| Architecture | Cores |
| --- | --- |
| x86_64 | `prescott`, `atom`, `core2`, `penryn`, `dunnington`, `nehalem`, `opteron`, `opteron_sse3`, `barcelona`, `nano`, `sandybridge`, `bobcat`, `bulldozer`, `piledriver`, `haswell`, `steamroller`, `excavator`, `zen`, `skylakex`, `cooperlake`, `sapphirerapids` |
| x86 | `katmai`, `coppermine`, `northwood`, `prescott`, `banias`, `atom`, `core2`, `penryn`, `dunnington`, `nehalem`, `athlon`, `opteron`, `opteron_sse3`, `barcelona`, `nano`, `bobcat`, `p5`, `p6`, `viac3`, `yonah`; the AVX cores and Zen run `nehalem`'s or `barcelona`'s, as upstream's 32-bit builds do |
| aarch64 | `armv8`, `cortexa53`, `cortexa55`, `cortexa57`, `cortexa72`, `cortexa73`, `falkor`, `thunderx`, `thunderx2t99`, `thunderx3t110`, `tsv110`, `emag8180`, `neoversen1`, `neoversen2`, `neoversev1`, `neoversev2`, `armv8sve`, `a64fx`, `armv9sme`, `vortex`, `vortexm4`, `cortexa76`, `cortexa510`, `cortexa710`, `cortexx1`, `cortexx2`, `ampere1`, `ft2000` |
| arm | `armv6`, `armv7`, `cortexa9`, `cortexa15` |
| powerpc, powerpc64 | `power3`, `power4`, `power5`, `power6`, `power8`, `power9`, `power10`, `ppc970`, `cell`, `ppc440`, `ppc440fp2`, `ppcg4` |
| mips, mips64 | `loongson3r3`, `loongson3r4`, `p5600`, `i6400`, `i6500`, `p6600` |
| loongarch64 | `la64_generic`, `la264`, `la464` |
| riscv64 | `riscv64_zvl128b`, `x280` |
| s390x | `z13`, `z14` |
| sparc, sparc64 | `sparcv7`, `sparc` |
| wasm32, wasm64 | `wasm128_generic` |

The converted assembly needs the target's calling convention to be one it was built for. x32,
MIPS n32 and soft-float ABIs other than arm's have none; building for them stops unless
`-Dblas-cores=none`.

### Not supported by Zig

These OpenBLAS platforms and cores are not converted; the portable kernels serve them where they
are targets at all.

| Platform or core | Why |
| --- | --- |
| alpha | no Zig backend |
| ia64 | no Zig target |
| csky | LLVM's backend emits neither assembly nor objects (Zig's `src/target.zig`, `hasLlvmSupport`); upstream uses the generic C kernels |
| e2k | no Zig target; upstream uses the generic C kernels |
| POWER BG/P (double FPU) | LLVM cannot assemble the "hummer" instructions |
| RISC-V C910V | LLVM has no `xtheadvector` |

## Usage

### Calling a routine

Each routine takes its element type as a comptime parameter. Sizes, strides and leading
dimensions are `isize`; a negative stride names the vector's last element first, as in the
reference. Every routine is `rstd.linalg.blas.<routine>` and takes `(allocator, io, T, ...)`,
or `(io, T, ...)` when it needs no scratch space beyond a small stack buffer. `io` carries the
threading ([`threading.md`](threading.md)): level 3 threads, all but `gemmt`; levels 1 and 2 are
serial for now. `rotg`,
`zrotg`, `rotmg` and `cabs1` work on scalars only and take neither.

| Level | Shape | Errors |
| --- | --- | --- |
| 1 | `axpy(io, T, n, alpha, x, incx, y, incy)` | `Level1Error` |
| 2 | `gemv(allocator, io, T, trans, m, n, alpha, a, lda, x, incx, beta, y, incy)` | `Level2Error` |
| 3 | `gemm(allocator, io, T, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)` | `Level3Error` |
| extensions | `omatcopy(io, T, trans, rows, cols, alpha, a, lda, b, ldb)` | `ExtError` |

Transposes, triangles, diagonals and sides are the letter enums `blas.Trans`, `Uplo`, `Diag` and
`Side`, shared by every level.

```zig
const blas = rstd.linalg.blas;

// C := A * B, all f64, column-major
try blas.gemm(allocator, io, f64, .n, .n, m, n, k, 1.0, a.ptr, lda, b.ptr, ldb, 0.0, c.ptr, ldc);
```

`blas.kernel(.daxpy_k)` reaches a kernel directly, with the reference's kernel arguments.

### Errors

Bad arguments return an error named after the reference's `xerbla` argument (`InvalidTrans`,
`InvalidM`, `InvalidLda`, `InvalidIncx`, ...), checked in the reference's order, so the first bad
argument wins. A vector or matrix too large for the address space returns its stride's or
leading dimension's error. Routines that allocate also return `error.OutOfMemory`, and every
routine that takes `io` may return `error.Canceled` once it threads. Nothing panics on bad
input. `laswp` has no `xerbla` checks upstream; rstd checks that its rows and pivots lie in `A`
before it swaps.

### Initialization

The first kernel call selects the core and fills the kernel table, reading `OPENBLAS_*` from the
process environment where a library can (Windows, and targets that link libc).
`blas.init(&env_map)` does the same with a map the program passes, which is how the variables
reach other targets; it must come before the first kernel call to count. `blas.shutdown()` goes
back to the portable kernels, and the next call selects again; no kernel call may be running.
`blas.coreName()` names the core whose kernels run.

## Numerical behavior

- **Semantics.** Results follow the reference BLAS: a zero `beta` never reads `y` or `C`, a zero
  `alpha` never reads the other operands, and `scal` with a zero `alpha` turns NaN and
  infinities into NaN. The rank updates (`ger`, `syr`, `her`, their packed and rank-2 forms)
  pass over a column whose factor is zero (both, for rank 2), so a NaN or infinity elsewhere in
  `x` or `y` never meets it. `dsdot` accumulates in `f64`; `rotg` computes in `XDouble`; `zrotg`
  follows LAPACK's `zlartg`.
- **Matrix copies.** `omatcopy` and `imatcopy` store zeros for a zero `alpha` without reading
  `A`, and copy its bits for `alpha` 1, complex too, where upstream's complex kernels multiply
  and turn an infinite part's partner into NaN; in place or through a buffer, `imatcopy` gives
  `omatcopy`'s bits. `geadd`'s zero `alpha` and `beta` read nothing, as `axpby`'s.
- **Zero strides.** Level 1 keeps the reference's sequential meaning: a zero input stride repeats
  `x(1)`, a zero output stride applies `n` updates to one element in turn. `axpy` with both
  strides zero is OpenBLAS's closed form, `y(1) += n*alpha*x(1)`. Level 2 rejects
  a zero stride.
- **Across cores.** Different cores run different kernels, with different fusion, accumulation
  order and blocking, so results may differ in the last bits, as in OpenBLAS; beyond rounding
  they agree (below). Pinning `OPENBLAS_CORETYPE` makes machines that can run that core agree.
- **Reproducibility.** For a given build and core, no result depends on where its operands lie in
  memory, nor on the thread count ([`threading.md`](threading.md)).
- **Denormals.** Rstd never changes the floating-point environment, and denormals are honored on
  every core: where a vector unit flushes them (ARMv7 NEON, pre-VSX AltiVec), `f32` stays out of
  its lanes ([`design.md`](design.md#denormals-and-the-floating-point-environment)).
- **Triangular solves.** `trsv`, `tbsv`, `tpsv` and `trsm` divide by the diagonal, complex by
  Smith's algorithm (`util.cdiv`), as Reference-LAPACK master does; OpenBLAS multiplies by a
  reciprocal, infinite for a denormal diagonal element. Level 2 and level 3 give the same
  quotients, except in an x86 build without SSE2 on a core with it: `trsv` divides on the
  x87, that core's `trsm` in SSE2.
- **Extended precision.** Wider than `f64`, the reductions and `gemv` sum in the order of
  OpenBLAS's x87 `q`/`x` kernels and give their bits.

### One contract on every core

Some slots never take a converted kernel (`vector.routed`), so their assembly is not converted:
every core runs the portable or `@Vector` kernels there, so what they do beyond rounding does not
depend on the core ([`comparisons.md`](comparisons.md), X1). Where the `@Vector` kernels go is in
[`design.md`](design.md#vector-level-1-and-level-2-kernels).

- **`nrm2`:** `s` and `c` sum the squares in `f64`, where each is exact and none overflows or
  underflows, and round the root once; `d` and `z` take Blue's algorithm, as Reference-LAPACK
  3.10's `dnrm2.f90` (`util.blue`). It is infinite for an infinity and NaN for a NaN. Upstream
  has three algorithms: x86_64 sums in `f64` (`nrm2_sse.S`), x86 and x86_64 in the x87's `f80`
  (`nrm2.S`), the rest scale a one-pass loop, several of them wrong on special values
  ([`issues.md`](issues.md)).
- **The extremes** (`i?amax`, `i?amin`, `i?max`, `i?min` and their values): the first NaN wins,
  its index and NaN for the value, as in Reference-LAPACK master (`50dc1bf`, `b358788`), BLIS's
  `amaxv` and arXiv:2207.09281; else the first extreme. Complex keys are `|re| + |im|`; where one
  overflows for a finite element, the index compares `|re|/4 + |im|/4`, as `izamax.f90` does,
  so an infinite part still wins. The converted kernels pass NaN over, or pick one by vector
  lane.
- **`sdot`:** products and sum in `f64`, rounded once. Upstream's sum in `f32` on most cores and
  in the x87's `f80` on x86, so the error of the same `sdot` differed by core by far more than
  its last bits.
- **No address dependence.** x86_64's reductions whose converted kernels peel elements until an
  operand is aligned (`vector.peels`: `asum_sse.S`, `asum_sse2.S`, `asum_atom.S`, `zasum_*.S`,
  `zsum_sse*.S`, `dot_sse2.S`, `dot_atom.S`, `ddot_bulldozer.S`, `zdot_*.S` and `dgemv_n*.S`)
  summed the same data otherwise elsewhere in memory. On the cores that use them the `@Vector`
  kernels, which load unaligned and never peel, take their place; the C kernels (`ddot.c` and the
  like) stay.
- **`axpy` with a zero `alpha`** multiplies on every core, so `0 * x` gives NaN for NaN and
  infinities. The interface returns first, as the reference does, but the triangular and banded
  drivers call the kernel with `-x(j)`, which may be zero; kernels that returned at once made
  `trsv` and `tbsv` propagate NaN otherwise on some cores, and are patched or ported without
  that return ([`design.md`](design.md#fixes-to-upstream-sources)).

### bfloat16 and half precision

The kernels are portable Zig, the same on every core.

- **Conversions** (`sbtobf16`, `bf16to`): round to nearest even; denormals flush to a signed zero
  and NaNs come out quiet, in both directions.
- **`bscal`:** `x := alpha*x` through `f32`. Its flag is 0 for a zero `alpha` giving zero, 1 for
  it keeping NaN and infinity as NaN; the test is on the bits of `alpha`, so only +0 counts.
- **`sbdot`:** both operands widened to `f32`, the products and their sum in `f64`, returned as
  `f32`. It differs from the reference's `sdot` on the widened data, and from the AVX-512 BF16
  kernel, by rounding only.
- **Half `gemv`** (`gemvHalf`: `bgemv`, `sbgemv`, `shgemv`): each result sums `a*x` in `f32` in
  index order, unfused, and rounds once on the store (`alpha*acc + beta*y`; `beta == 0` never
  reads `y`). Without a transpose, rows go in blocks of 256 with the columns outermost, so `A`
  is read down its columns. For `bgemv` the scalars and `y` are bfloat16. Arguments are checked
  as `interface/sbgemv.c` does; `alpha == 0` only scales `y` by `beta`. With `m` or `n` zero the
  kernel yields `beta*y`; the interface returns first.
- `sbdot` and the others return at once for `n <= 0`.

## Apple Accelerate

On macOS the public routines over `s`, `d`, `c` and `z` run Apple's Accelerate where it has
the routine, bound at run time with no SDK and run single-threaded, so its results stay
reproducible (`apple_accel.zig`). Rstd's checks come first, so the errors are the same
everywhere. Where its results differ from the reference's on special values, those calls run
rstd's kernels instead, so [one contract](#one-contract-on-every-core) holds there too.
`RSTD_BLAS_ACCELERATE=0` at init, or `blas.useAccelerate(false)`, turns it off. Which routines
and quirks are in [`design.md`](design.md#apple-accelerate).

## Build options

| Option | Values | Effect |
| --- | --- | --- |
| `-Dblas-cores` | `target`, `all`, `none` | which per-core kernel objects are built: only the core the target maps to (the host's own for a native target), every core the target can run, or none (portable kernels only). Default: `target` in debug, `all` otherwise |
| `-Dblas-llvm` | `true`, `false` | forces the kernel objects on or off LLVM; unset uses LLVM where the kernels need it, which today is everywhere |

A consumer using `src/buildutils` `addRstd` passes its optimize mode (`RstdOptions.optimize`) so
the default follows it, and may set `RstdOptions.blas_cores`. A core with no object in the build
runs the closest linked core's kernels, or the portable ones.

## Environment variables

| Variable | Effect |
| --- | --- |
| `OPENBLAS_CORETYPE` | forces a core by name, case-insensitive; one the CPU cannot run is replaced by the closest one it can, and an unknown name leaves detection to decide |
| `OPENBLAS_VERBOSE` | 1 reports fallbacks and unknown names; 2 also the selected core |
| `GET_OPENBLAS_CORETYPE` | POWER: when set, prints the selected core |
| `OPENBLAS_BLOCK_FACTOR` | scales the computed gemm block size, in percent (10 to 200) |
| `OPENBLAS_L2_SIZE` | overrides the detected L2 size used by the blocking rules |
| `OPENBLAS_NUM_THREADS`, `OPENBLAS_DEFAULT_NUM_THREADS` | the thread ceiling, and LA464's blocking set |
| `GOTO_NUM_THREADS` | the thread ceiling when neither of those is set |
| `LD_HWCAP_MASK` | s390x: masks the hardware capabilities before detection |

`OMP_NUM_THREADS`, `OMP_ADAPTIVE`, `OPENBLAS_THREAD_TIMEOUT`, `GOTO_DEBUG` and
`OPENBLAS_MAIN_FREE`/`GOTOBLAS_MAIN_FREE` are read and have no effect. Values longer than 64
bytes are ignored. Details are in [`design.md`](design.md#configuration).

## Tests

| Command | Runs |
| --- | --- |
| `zig build rstd-test` | every test natively; add `-Dblas-cores=all` to test every core the host runs, `-Drstd-test-filter=<name>` to pick tests |
| `zig build qemu-test` | every emulated machine |
| `zig build qemu-test-<machine>` | one machine (`aarch64`, `aarch64-sve`, `arm`, `arm-softfp`, `riscv64`, `loongarch64`, `mips64r2el`, `mips64r6el`, `mips32r2el`, `mips32r5el`, `mips64r2`, `mips64r6`, `mips32r2`, `mips32r5`, `i386`, `x86_64`, `sparc64`, `sparc64-debug`, `powerpc64le`, `powerpc64le-power10`, `powerpc64`, `powerpc`, `powerpc-pseries`, `s390x`); `-Dqemu-test-filter=<name>` to pick tests, `-Dqemu-test-part=<part>` one image of a rig split in parts |
| `zig build qemu-test-x86_64-whpx` | the host CPU under WHPX (x86_64 Windows hosts): the AVX-512, FMA4 and XOP cores |
| `zig build qemu-test-linux-<arch>` | the full suite under Linux (`aarch64`, `x86_64`, `arm`, `riscv64`, `powerpc64le`) |
| `zig build qemu-test-x86_64-macos` | the unit tests in an x86_64 macOS VM, set up once (`build/qemutest.zig`) |
| `zig build ssh-test-aarch64-macos -Dssh-mac=<ssh alias>` | the full suite on an Apple Silicon Mac |
| `zig build ssh-test-x86_64-macos` | the x86_64 macOS suite on that Mac, under Rosetta 2 |
| `zig build ssh-bench-aarch64-macos -Dopenblas=<checkout>` | the kernel benchmarks on that Mac ([`results-m6.md`](../../../tools/blasbench/results-m6.md)) |
| `zig build wasm-test` | level 1 and 2 under Node, with SIMD128 |
| `zig build rstd-check-macho`, `zig build rstd-check-uefi` | build the tests for macOS and UEFI without running them |

Test names follow the file's path (`linalg.blas._test.blas1.edge.test...`), so a filter can pick
a file. The qemu rigs' QEMU runs on another machine, `-Dtest-host=<ssh alias>` (`macmini-m6` by
default when rstd is the root project; `-Dtest-host-off` runs it here). `-Dfail-fast` ends a rig
at its first failing test. Rig logs are in `.zig-cache/qemu-test-<machine>.log`. What each rig
covers, and how the tests are built, is in [`design.md`](design.md#verification).

## Further reading

- [`design.md`](design.md): structure, core selection, kernel provenance and the fixes to
  upstream sources, Accelerate, testing.
- [`threading.md`](threading.md): the threading design, the `io: Io` API and its guardrails.
- [`comparisons.md`](comparisons.md): the numerical design decisions against BLIS, the reference
  and the exception-handling proposal.
- [`issues.md`](issues.md): problems found in OpenBLAS, Zig, QEMU and other tools, with
  workarounds.
- [`blas1/c_kernels.md`](blas1/c_kernels.md): porting OpenBLAS's C level-1 kernels.
- [`blas2/buffers.md`](blas2/buffers.md): the level-2 kernels' scratch buffers.
- [`blas3/gemm/notes.md`](blas3/gemm/notes.md): gemm kernel measurements.
- [`results-m6.md`](../../../tools/blasbench/results-m6.md): Zig against converted assembly on
  Apple M6.
