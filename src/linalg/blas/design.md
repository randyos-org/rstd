# BLAS design

This document specifies how rstd's BLAS is structured: its requirements, the layers from the
public interface down to the per-core kernels, how a core is detected and selected, where every
kernel comes from, and how all of it is tested. It describes the code as it is. Usage is in
[`readme.md`](readme.md), threading in [`threading.md`](threading.md), problems in other
projects in [`issues.md`](issues.md), and the numerical decisions measured against BLIS and the
reference in [`comparisons.md`](comparisons.md).

Upstream references are to the OpenBLAS checkout that `tools/blas_asm` converts, v0.3.34-317
(commit a69062524), unless a file names another; `Vtable.zig` is derived from v0.3.34-127 (commit
9724481b). Paths are relative to `src/linalg/blas/` unless they start with `build/`, `tools/`,
`test/` or `src/`. The toolchain is Zig 0.18.0-dev.751.

## Requirements

- **A Zig-native surgical port.** OpenBLAS's structure is kept where it carries meaning: the
  `gotoblas_t` dispatch table, the `interface/` checks, the `driver/level2` and `driver/level3`
  loops, the per-core kernel selection of `kernel/<arch>/KERNEL.<core>`, and the parameters of
  `param.h`. Everything around it is Zig: no C sources, no build-time code generation, no
  `common.h` macros.
- **Generic, type-dispatched interfaces.** One function per routine takes the element type as a
  comptime parameter (`gemm(allocator, f64, ...)`), not one per precision letter. Any float type
  and its `std.math.Complex` works; the types with vtable slots run the selected core's kernels,
  the others a portable kernel built for the type.
- **Eager per-core install.** Selecting a core fills the whole vtable at once, before the first
  kernel runs. No routine checks or patches the table lazily; a kernel call is one indirect call
  through a slot.
- **Per-core code only where `@Vector` cannot match it.** A kernel of its own for one core
  exists where a generic `@Vector` kernel is measurably slower or cannot express the
  instructions (SVE, MMA, AMX, SME). The level-3 kernels are one `@Vector` source for every core
  ([`blas3/gemm/notes.md`](blas3/gemm/notes.md),
  [`results-m6.md`](../../../tools/blasbench/results-m6.md)).
- **No CBLAS or Fortran ABI.** There are no `cblas_*` or `dgemv_` symbols, no row-major layer
  and no `xerbla` callback. Bad arguments return named errors. This is a decision, so the audit
  gaps for those layers are closed by it, not by code.
- **One numerical contract per routine.** Every core gives a routine's results the same meaning;
  cores differ in rounding only ([`readme.md`](readme.md#one-contract-on-every-core)).
- **Bit-exactness only against converted code.** A Zig kernel standing in for a converted kernel
  of the same core (an x87 `q`/`x` reduction, a C driver around a microkernel, an intrinsics
  block) gives that kernel's bits. Everything else is held to the reference semantics within a
  stated tolerance, not to OpenBLAS's accumulation order.
- **Reproducibility.** For a given build and core, results depend neither on scheduling, the
  number of threads, nor where the operands lie in memory ([`threading.md`](threading.md)).
- **Errors, never panics.** Argument errors and `error.OutOfMemory` return to the caller.
  Nothing allocates behind the caller's back or retries on failure.
- **Allocator conventions.** A routine that needs scratch beyond a small stack buffer takes a
  `std.mem.Allocator`; nothing uses a global allocator.
- **Every target.** The same source serves hosted, freestanding, UEFI and WebAssembly targets.
  A target with no per-core kernels runs the portable ones, with no special cases in callers.
- **Zig 0.17.0 as well as 0.18.** zfls vendors rstd with `minimum_zig_version = "0.17.0"`
  (its `build.zig.zon`) and builds it on Zig 0.17.0 for x86_64 and aarch64 Linux, macOS 11 and
  Windows (its `.github/workflows/ci.yml`). A shim for a 0.17 difference detects the API
  (`@hasDecl`, `@hasField`, `zig_backend`), and checks the version only for a std bug with
  nothing to detect (`threading.Evented`, [`issues.md`](issues.md)). `rstd-check-zig017`, run
  with 0.17.0, builds the unit tests for those targets with every core's kernels.

## Background

OpenBLAS splits a call into four parts:

- `interface/<routine>.c` checks arguments (calling `xerbla` with the first bad argument's
  number), handles quick returns and negative strides, and picks a driver;
- `driver/level2/*.c` and `driver/level3/level3.c` loop over kernels, blocking and packing the
  operands;
- `kernel/<arch>/` holds the kernels, one source per routine and precision for each core,
  chosen by `KERNEL.<core>`, then the architecture's `KERNEL`, then `kernel/Makefile.L1`, `L2`
  and `L3`;
- `gotoblas_t` (`common_param.h`) is the table the drivers call kernels through. A
  `DYNAMIC_ARCH` build compiles the kernel tree once per core with that core's flags, and
  `driver/others/dynamic*.c` picks one table at load time. `kernel/setparam-ref.c` fills each
  table's blocking parameters from `param.h`.

## Design

### Layers

```text
caller
  |
  v
interface      blas1/interface.zig, blas2/interface.zig, blas2/triangular.zig,
               blas2/banded.zig, blas3/interface.zig
               checks, quick returns, negative strides, scratch buffers
  |
  v
drivers        blas2 (triangular, banded, symmetric), blas3/driver.zig
               loops, blocking, packing
  |
  v
dispatch       portable.kernel / portable.fetch, blas.kernel
               element type -> vtable slot, or a portable kernel for types without slots
  |
  v
Vtable         one slot per gotoblas_t kernel, filled when a core is selected
  |
  +--> portable kernels    in the main library (portable.install)
  |
  +--> per-core kernels    in per-core objects (build/blaskern.zig, _kernel.zig):
                           converted assembly, ported C drivers, @Vector kernels
```

Level 1 has no driver: the interface calls the kernel. Every level-3 routine reaches its kernels
through `blas3/driver.zig`.

### Public interface

The routines are Zig functions over column-major data, in the Fortran argument order:

| Level | File | Routines |
| --- | --- | --- |
| 1 | `blas1/interface.zig` | `copy`, `swap`, `scal`, `scalReal`, `dot`, `dotc`, `dsdot`, `sdsdot`, `nrm2`, `asum`, `sum`, `index`, `extreme`, `rot`, `rotm`, `axpy`, `axpyc`, `axpby`, `rotg`, `zrotg`, `rotmg`, `cabs1`, `sbdot`, `bscal`, `sbtobf16`, `bf16to` |
| 2 | `blas2/interface.zig` | `gemv`, `gemvHalf`, `ger`, `gerc`, `symv`, `hemv` |
| 2 | `blas2/triangular.zig` | `trmv`, `trsv`, `tbmv`, `tbsv`, `tpmv`, `tpsv` |
| 2 | `blas2/banded.zig` | `gbmv`, `sbmv`, `hbmv`, `spmv`, `hpmv`, `spr`, `syr`, `hpr`, `her`, `spr2`, `syr2`, `hpr2`, `her2` |
| 3 | `blas3/interface.zig` | `gemm`, `gemmt`, `symm`, `hemm`, `syrk`, `herk`, `syr2k`, `her2k`, `trmm`, `trsm` |
| 3 | `blas3/half.zig` | `gemmHalf`, `gemmtHalf` |
| 3 | `blas3/batch.zig` | `gemmBatch`, `gemmBatchStrided`, `gemmHalfBatch`, `gemmHalfBatchStrided` |
| extensions | `ext/interface.zig` | `geadd`, `omatcopy`, `imatcopy`, `laswp` |

Conventions shared by every interface:

- **Element type.** `T` is a comptime parameter: any float type or its `std.math.Complex`
  (`types.Float(T)` must be a float). `BF16` and the mixed half routines have their own entry
  points (`sbdot`, `bscal`, `gemvHalf`, `gemmHalf`, `gemmtHalf`).
- **Sizes.** Sizes, strides and leading dimensions are `isize`, not the reference's Fortran
  `INTEGER` (`blasint`), since there is no Fortran ABI; `types.BlasLong`, the kernels' index
  type, is `isize` too.
- **Negative strides.** A negative stride names the vector's last element first: the interface
  moves the pointer to it and the kernel walks back from there, as the reference's interfaces do.
- **Zero strides.** Level 1 keeps the reference's sequential meaning: a zero input stride
  repeats `x(1)`, a zero output stride applies `n` updates to one element in turn. `axpy` with
  both strides zero is OpenBLAS's closed form, `y(1) += n*alpha*x(1)`. Level 2
  rejects a zero stride (`InvalidIncx`).
- **Errors.** Each file has an `ArgError` set named after the reference's `xerbla` arguments
  (`InvalidTrans`, `InvalidLda`, `InvalidIncx`, ...), checked in the reference's order, so the
  first bad argument wins. Each level's `Error` adds `Allocator.Error` where it allocates and
  `Io.Cancelable` everywhere; the root exports them as `Level1Error` to `Level3Error`, and
  the extensions' as `ExtError`.
  A vector or matrix that could not fit the address space (`util.fits`, `util.fitsMatrix`)
  returns its stride's or leading dimension's error, so no index computation overflows.
- **Letters.** Transposes, triangles, sides and diagonals are one set of letter enums for every
  level (`types.zig`). Real `r` and `c` transposes mean `n` and `t`. `Trans` has the complex
  forms `o`, `u`, `s` and `d` (conjugated `x`), and `Uplo` the `v` and `m` forms of `hemv`
  and `hbmv` (`HEMVREV`). A letter a routine does not take is an error: `InvalidTrans` for `o`
  to `d` outside complex `gemv` and `gbmv`, `herk`'s `t` and complex `syrk`'s `c`;
  `InvalidUplo` for `v` and `m` outside `hemv` and `hbmv`, the only routines whose Fortran
  interface takes them upstream (`hpmv`, `her`, `hpr`, `her2` and `hpr2` get them from
  row-major CBLAS alone).
- **Leading dimensions.** Level 3 checks them against `max(1, rows)`, as the reference BLAS
  does; OpenBLAS lets `0` through for empty matrices.
- **Quick returns.** Empty problems return before any kernel call; kernels assume `n >= 1`.
  Level 3 does not forward one-column problems to `gemv`.
- **Argument order.** Every routine takes `(allocator, io, T, ...)`, or `(io, T, ...)` without
  an allocator (level 1, `gemvHalf`), as [`threading.md`](threading.md) specifies; the routines
  that thread run their workers on it. `rotg`, `zrotg` and `rotmg` compute from scalars, never
  thread, and take `(T, ...)` and return `void`; `cabs1` takes `(T, z)` and returns
  `|re| + |im|`.

`__root__.zig` re-exports every routine flat (`blas.dot`, `blas.gemv`, `blas.gemm`), with the
argument enums and error sets, beside the setup and dispatch API: `init`, `shutdown`, `vtable`,
`table`, `kernel`, `coreType`, `coreName`, `config`, `Config`, `Vtable`, `Slot`, `types`, `cpu`,
and `vector` for `tools/blasbench`. The interface files stay where they are.

### Planned: a Zig-idiomatic interface

Decided (2026-10-26), not yet implemented; the interface above is the reference's shape, which
users are not assumed to know. Part of the codebase-wide pass to make rstd more Zig-idiomatic.

- **Slices.** Vectors and matrices are slices (`[]const T`, `[]T`), not many-pointers, so a
  short buffer is caught instead of overrun. A size the slice already carries is not a separate
  argument; strides and leading dimensions remain only where they describe a layout inside the
  slice. Each routine checks that its slices cover what its sizes and strides need.
- **Errors, not panics, for bad input.** A wrong length, size or stride, or a letter the routine
  does not take, is an input error the caller can fix and retry: it returns an error.
- **Panics only when unrecoverable.** A broken invariant found mid-calculation, or a state
  where memory may be corrupt (above all from converted assembly kernels), panics. Each
  routine's doc comment says whether it panics explicitly; callees need not be traced.
- **Exact error sets.** A routine's error set holds only what it can return: its own argument
  errors; `error.OutOfMemory` only if it allocates; `error.Canceled` only if it can be
  canceled. Sets are named and per routine or routine family, built from shared pieces,
  rather than one set per level.
- **Plain names.** Errors say what is wrong in words, not as Fortran argument names: for
  example `SliceTooShort`, `NegativeSize`, `ZeroStride`, `LeadingDimensionTooSmall`,
  `UnsupportedTranspose`, instead of `InvalidM`, `InvalidLda`, `InvalidIncx`. The checks no
  longer need to run in the reference's order.
- **Docs.** User-facing docs explain in plain terms and do not lean on BLAS jargon (`xerbla`,
  `LDA`, `INCX`).

### Types

`types.zig` defines the element and index types:

| Type | Definition | Precision letters |
| --- | --- | --- |
| `f32`, `f64` | native | `s`, `d` |
| `XDouble` | `c_longdouble` (x87 `f80`) with `exprecision`, `f128` with `quad_precision`, else `f64` | `q` |
| `Complex32`, `Complex64`, `ComplexXD` | `std.math.Complex` of the above | `c`, `z`, `x` |
| `BF16` | `enum(u16) { _ }`: a storage format, widened to `f32` for arithmetic | `sb`, `b` |
| `HF16` | `f16` | `sh` |
| `BlasLong` | `isize`, the kernels' `BLASLONG` | |

- **`exprecision`** is on for x86 and x86_64 with an x87 unit whose `long double` is 80 bits,
  except on Microsoft-ABI targets (Windows, UEFI), which the reference never builds it for
  (`exprecisionOn`; `Makefile.system:487` sets `NO_EXPRECISION` on `WINNT`). It assumes the x87 control word's extended precision, as Linux and x64
  UEFI leave it (0x037F: Linux v6.12 `arch/x86/kernel/fpu/core.c`, `fpstate_init_fxstate`;
  UEFI Specification 2.10, section 2.3.4,
  <https://uefi.org/specs/UEFI/2.10/02_Overview.html>). Elsewhere, `q`/`x` alias `d`/`z`, as in the reference; `blas_types` is a
  membership list for that reason, since a `switch` cannot hold one type twice.
- **`quad_precision`** is off; turning it on makes `XDouble` Zig's `f128` (compiler-rt
  arithmetic) and is exclusive with `exprecision`.
- **`f80` size.** `f80` keeps Zig's native size (16 bytes on x86_64, 12 on x86), where the
  reference's GCC builds pass `-m128bit-long-double` (`Makefile.system:602`, `:627`). Buffer formulas count bytes and round up to
  whole elements; the gemm blocking keeps the reference's 16-byte figures.
- **Complex results.** Kernels that return a complex value use `CComplex(F)`, an
  `extern struct`, since `std.math.Complex` has automatic layout. At `f80` that struct returns
  through a hidden pointer on SysV rather than in `st(0)`/`st(1)` as the x87 `xdotu`/`xdotc`
  do, so only Zig kernels fill those slots.
- **Helpers.** `isComplex`, `Float` (the part type), `hasSlots` (`f32`, `f64`, `XDouble`),
  `kernelCc` (C convention for types with slots, Zig's own otherwise), `Acc` (`f32` for types
  narrower than 32 bits) and `isExtended` (wider than `f64`).

### The vtable

`Vtable.zig` is OpenBLAS's `gotoblas_t`, mechanically derived from `common_param.h`: 854 kernel
slots, its 87 integer parameters, rstd's `symv_p` and the `init` hook. It leaves out the 47 slots only row-major data
or CBLAS reaches (`?gerv_k`, `?gerd_k`, the `*_direct*` kernels, `?omatcopy_k_r*`,
`?imatcopy_k_r*`) and the `sbger_k`, `sbsymv_L` and `sbsymv_U` that nothing upstream calls
(`kernel/setparam-ref.c:112-113` fills them with `sger_k`, `ssymv_L` and `ssymv_U`); its doc comment lists them (`unused/`, below). The reference guards fields behind `BUILD_*`,
`EXPRECISION`, `ARCH_*` and `SMALL_MATRIX_OPT`; here every field is always present, kernels as
optional pointers that stay `null` where nothing applies. Slots use the C calling convention and
the kernel ABI of the reference (`daxpy_k(n, dummy, dummy, alpha, x, incx, y, incy, dummy,
dummy)`), so a converted kernel installs as is.

`Vtable.zig` imports nothing outside `src/linalg/blas/`, since the per-core objects compile it
on their own; its test lives in `__root__.zig`.

### Initialization

`__root__.zig` keeps one global `Vtable`, the `Config` it was built from, and the selected core.

- **Before selection** the table is `generic_vtable`: the portable kernels and the generic
  target's blocking, built at comptime.
- **`blas.init(env)`** selects the core and fills the table, once; it returns whether this call
  did it. The first kernel call does the same with the process environment where a library can
  read one (`processEnv`: the Windows PEB, libc's `environ`, `_NSGetEnviron()` on Darwin) and
  with none elsewhere, so `blas.init` is how `OPENBLAS_*` reaches freestanding and non-libc
  builds. It must come before the first kernel call to count. Values are copied, so `env` may go
  away after.
- **Handshake.** `phase` moves `uninit` to `running` by compare-exchange; the winner builds the
  table while any other caller spins until it publishes `ready`. A mutex would bring threading
  support into freestanding builds.
- **Building the table.** `cpu.choose` picks the core (below); `cpu.vtableFromCoreType`
  (`cpu/<arch>.zig`) runs `portable.install`, then `blas3/gemm/params.installCore` for the
  core's blocking and `dtb_entries`, then `cpu.kernel.install` for its kernels, then
  `portable.finish`. A non-null `init` slot runs last.
- **Access.** A routine fetches the table once per call (`blas.vtable()`, which selects first if
  needed, or `blas.kernel(.slot)` for one kernel); code below it reads `blas.table()`, which
  asserts the table is ready.
- **`blas.shutdown()`** returns to `generic_vtable` and the unset state; the next kernel call or
  `init` selects again. No kernel call may be running, as for the reference's `gotoblas_quit`.
- **Tests** replace the table with `testUse`/`testRestore`, without selecting a core.

### Type dispatch and the portable kernels

`portable.zig` connects element types to slots:

- `portable.Op` names a kernel less its precision letter (`axpy_k`, `gemv_n`, `trsm_kernel_LT`).
- `slot(T, op)` builds the slot name from `letter(T)` (`saxpy_k`, `isamax_k`, complex `csrot_k`).
  A type with slots that lacks one is a compile error, except where the reference has no
  extended-precision slot (`no_extended_slot`: `axpby`), which runs the portable kernel.
- `kernel(T, vt, op)` returns the vtable's kernel for a type with slots and the portable kernel
  built for `T` otherwise; `fetch(T, op, getter)` does the same with the table from
  `blas.vtable` or `blas.table`.
- `install` fills a vtable with every level's portable kernels (`installGeneric` of `blas1`,
  `blas2`, `blas3`) and the generic target's blocking; `finish` runs after any core's installation task:
  it aliases the bfloat16-named slots to the `f32` kernels (`ext/bf16/kernel.zig`'s `alias`) and
  restores the portable gemm kernels' unrolls where those kernels stayed
  (`blas3/gemm/kernel.zig`'s `settle`).

The portable kernels are Zig algorithms with the reference's documented semantics, not copies of
its C. What they compute is the contract every core keeps
([`readme.md`](readme.md#one-contract-on-every-core)):

- `sdot` and `dsdot` multiply and sum in `f64`; the portable `dot` sums sequentially, not in
  `generic/dot.c`'s unrolled order.
- `nrm2` sums the squares in `f64` for `s` and `c`, and takes Blue's algorithm for `d` and `z`
  (`util.blue`, after Reference-LAPACK 3.10's `dnrm2.f90`); its results on special values are
  LAPACK 3.10's. Sources: J. L. Blue, "A Portable Fortran Program to Find the Euclidean Norm of
  a Vector", ACM TOMS 4(1), 15-23, 1978, <https://doi.org/10.1145/355769.355771>; E. Anderson,
  "Algorithm 978: Safe Scaling in the Level 1 BLAS", ACM TOMS 44(1), 1-28, 2017,
  <https://doi.org/10.1145/3061665>, as `dnrm2.f90`'s header cites them.
- `scal`'s flag argument keeps NaN and infinity for a zero `alpha` when set.
- In the extremes (`i?amax`, `i?amin`, `i?max`, `i?min` and their values) the first NaN wins,
  else the first extreme; complex keys that overflow compare at a quarter, as `izamax.f90` does
  (Reference master, 0234cae25ceb01807a9b84333b732756559be4aa; not in 3.12.1).
- `rotg` computes in `XDouble`, where the reference uses `long double` on x86
  (`interface/rotg.c:26`); `zrotg` follows
  LAPACK's `zlartg`, since `zrotg.c`'s scaled branches are wrong ([`issues.md`](issues.md)).
- Wider than `f64`, the reductions and `gemv` take the order of the x87 `q`/`x` kernels and give
  their bits: `extendedDot` and `extendedZdot` in `blas1/dot/__root__.zig`, `util.extendedSum`
  for `asum` and `nrm2`, and `extendedGemv` and `extendedZgemv` in `blas2/gemv/__root__.zig`
  (`installExtendedCore` sets x86's per-core column block).

Not reproduced, by design: the reference's accumulation order elsewhere (`generic/dot.c`,
`sum.c`), its signed-zero handling for `axpby` with a zero `alpha`, its threading, and
`zrotg.c`'s internals.

### Contraction

`fma.zig` decides whether a multiply-add fuses, the way the reference's compile flags decide it:

- **Hardware.** Fused where the ISA has a scalar fused multiply-add: x86 with FMA, AArch64, ARM
  with VFPv4, PowerPC, RISC-V F/D, LoongArch, s390x, MIPS release 6. Not on SPARC (its FMA extension of SPARC64 VI and `sparc4`, Oracle Solaris Studio 12.4 Fortran guide, `-xarch`, <https://docs.oracle.com/cd/E37069_01/html/E37076/aevia.html>, is left off, as OpenBLAS's gcc builds take `-mcpu=v9` without `-mfmaf`, `Makefile.sparc:7`), WebAssembly
  or earlier MIPS, whose `madd` rounds the product (Loongson 3A's fuses, and LLVM uses it for a
  plain `a*b + c`: `issues.md`). Vector multiply-adds also fuse on MIPS with
  MSA and WebAssembly with `relaxed_simd` (`vector_madd`).
- **Loongson 3A.** The `loongson3r3` and `loongson3r4` objects fuse (`Spec.fuse_fma`), and on MIPS
  before release 6 that means `madd.fmt` fuses (`fused_madd`): their multiply-adds are inline
  `madd.s`/`madd.d`, `@mulAdd` being a libm `fma` call there, and `loongson3r4`'s vectors MSA
  `fmadd` by inline asm, past LLVM's swapped operands. LLVM's own `madd.fmt` for a plain `a*b + c`
  fuses there too. Other MIPS objects and the main library never take the asm.
- **Soft float.** A soft-float ABI (Zig's `eabi`, `gnueabi` and the like, which compute every
  float in software) never fuses: `@mulAdd` there is compiler-rt's `fma`, several soft
  operations, and a C compiler's `-msoft-float` has no fused instruction to contract to. On the
  `powerpc-750` rig the driver's beta-folding test took 754 s fused and 188 s not (build
  included), and the full Debug run outran the watchdog's hour in one test.
- **Per-core policy.** A root may declare `blas_fuse_fma` and `blas_fuse_fma_vector`; the
  per-core objects do (`_kernel.zig`), from the core's `Spec.fuse_fma`. The AMD cores' FMA3/FMA4
  features exist so that the assembler accepts their microkernels; their flags leave contraction off
  (`kernel/Makefile:105` adds no `-mfma` for the Bulldozer family).
- **Width.** Floats wider than `f64` never fuse (`fuseFor`): `@mulAdd(f80)` is a soft `fmal`
  call, and the x87 kernels never fuse.
- **Expression rule.** Contraction follows a C compiler's within one expression (clang's
  `-ffp-contract=on`, "fusion in the same statement", Clang User's Manual,
  <https://clang.llvm.org/docs/UsersManual.html>): `a*b + c*d` is `fma(a, b, c*d)`, a longer chain fuses each product into
  the running sum from the left, and a product in its own statement, or converted before the
  add, stays rounded. `rsimd.mulAdd` delegates here.

The portable kernels follow the build target's policy, so two builds for different CPUs may
round a portable kernel differently; reproducibility is per build and core.

### Denormals and the floating-point environment

Rstd never changes the caller's floating-point environment (MXCSR, FPCR, FPSCR, `VSCR`) and
assumes IEEE default modes. Its worker threads adopt the caller's control state (`util.FpState`:
the x87 control word, MXCSR, FPCR, FPSCR, `frm`, FCSR, FPC, FSR; not `VSCR`), so a threaded call
rounds as its serial path does under any modes ([`threading.md`](threading.md#guardrails)).
Denormals are honored everywhere: where a vector unit flushes `f32` denormals
(ARMv7 NEON, pre-VSX AltiVec, whose `VSCR[NJ]` Linux sets: v6.12
`arch/powerpc/kernel/process.c`, `start_thread`), `util.vectorIeee` is false, and no kernel puts `f32` in its lanes.

- The portable `saxpy` takes its scalar loop there; `vector.zig` has no `s` or `c` roles there.
- The level-3 `f32` tiles are scalar there (`blas3/gemm/kernel.zig`'s `vectorBits`).
- LLVM's loop and SLP vectorizers would still put scalar `f32` code in AltiVec lanes, so the
  AltiVec cores' objects are built without AltiVec (`cpu/_corespec.zig`'s `without`), and the
  portable loops pass each `f32` result they store or accumulate through `util.scalar`, an empty
  asm they cannot widen, compiled out everywhere else ([`issues.md`](issues.md)). ARMv7's
  vectorizers keep off NEON's `f32` lanes themselves (LLVM 21.1.0
  `llvm/lib/Target/ARM/ARMTargetTransformInfo.h:118`, `isFPVectorizationPotentiallyUnsafe`, true
  but on Darwin and MVE; checked by `Transforms/Vectorize/SLPVectorizer.cpp:9214` and `LoopVectorize.cpp:10042`).

[`comparisons.md`](comparisons.md) (10) records how this was found and measured.

### Per-core kernel objects

OpenBLAS compiles its kernel tree once per core with that core's flags. Rstd does the same with
one object per group of cores whose flags are equal:

```text
cpu/_corespec.zig ----------------------------+
  Spec per core: model or features, fuse_fma  |
                                              v
build/blaskern.zig: kernelCores(arch), runsOn(target cpu), group by Spec.eql
  |  one b.addObject per group, root _kernel.zig, .optimize = .fast, LLVM
  |  blaskernel_options: cores, fuse_fma_set, fuse_fma
  |  blaskern_linked: cores linked, fallback core
  v
_kernel.zig: kernels.ForArch(arch).exportCore(core) for each core
  |  exports blaskern_<slot>_<core>, hidden
  v
kernels/<arch>.zig: modules = the routines with kernels on this arch
  |  install(vt, core) -> asm_registry.install -> @extern each symbol
  v
cpu/<arch>.zig: vtableFromCoreType(core)
```

- **Spec.** `cpu/_corespec.zig` is the single description of each core's compile flags: a CPU
  model where the reference passes `-march`/`-mcpu`, else the baseline plus features, from
  `Makefile.<arch>`, the `HAVE_*` tables of `getarch.c` and `kernel/Makefile`; `without` takes
  features out of a model. Only the ISA is described; `-mtune` has no counterpart. The build and
  the runtime both read it, and a test checks every named model and feature exists. A target
  short of some of the baseline's ISA replaces the baseline (`objectQuery`): an x86 one without
  SSE2 (Zig's baseline is the Pentium 4), a macOS one older than Nehalem, a SPARC V9 without
  VIS (the baselines `pentium4`, `nehalem` and `ultrasparc`: Zig `lib/std/Target.zig:2188-2202`). So a Pentium III build's objects compute `f64` on the x87, as its main library does, and
  run on the cores they serve.
- **ISA check.** `machine` in `_corespec.zig` is what each core's CPUs have, from the hardware
  (the LLVM model of the CPU a core is named for, or of the weakest part it covers, plus what
  `runnable` requires), apart from the spec. `excess` is what an object uses beyond that and
  beyond the target; the build stops on any (`checkIsa` in `build/blaskern.zig`), and a test
  tries every architecture with its baselines and each core's own CPU as targets. `undecided`
  lists spec features waiting on a decision; keep it empty. The check sees
  compiler-chosen instructions only: x86 inline assembly is not held to the target's features.
- **Groups.** `build/blaskern.zig` takes `corespec.kernelCores(arch)` (the cores with kernels,
  derived from the registries and checked by a test in `cpu/__root__.zig`), keeps those whose
  kernels run on the target CPU (`corespec.runsOn`), and groups equal specs into one object.
- **Build mode.** The objects always build `-OReleaseFast`, with LLVM (`needsLlvm`: the
  self-hosted x86_64 assembler rejects the converted assembly's register memory operands and
  VEX/EVEX mnemonics, the self-hosted wasm backend emits invalid SIMD, and the other
  architectures have no self-hosted backend that covers them). `-Dblas-llvm` overrides it. They
  link with LLD, except on Mach-O, which Zig's LLD driver cannot write (Zig's `src/link/Lld.zig`
  at 069480a1a takes COFF, ELF and wasm only); the self-hosted linker misplaced
  SPARC symbols when tried (calls landed in zeroed `.text`; our observation, b8f875b, not
  reproduced since and not in `issues.md`).
- **Selection of objects.** `-Dblas-cores` (`BlasCores` in `src/buildutils/rstd.zig`): `target`
  builds only the group the target maps to (the host's detected core for a native target, else
  the group with the most features the target has), `all` every group the target can run,
  `none` no objects. The default follows the consumer's optimize mode: `target` in debug, `all`
  otherwise. `blaskern_linked` tells the runtime which cores are linked and which one the build
  picked as a fallback.
- **ABI variants.** `abi.zig` picks the converted assembly's build variant from the target's C
  calling convention, not its OS: x86_64 `sysv` or `win` (so x86_64 UEFI takes `win`), arm `hf`
  or `sf`, POWER `32`, `32le`, `64` and `64be`, SPARC `32` and `64`. A target with kernels but no
  matching variant (x32, MIPS n32, a soft-float ABI anywhere but arm) stops the build unless
  `-Dblas-cores=none`. `asm_registry.Kernels.exportCore` fails to compile for a core that lacks
  the variant, so the generator and `abi.zig` cannot drift apart.
- **Symbols.** Each kernel is exported as `blaskern_<slot>_<core>` (`util.kernelSymbol`) with
  hidden visibility, so none leaves a shared library rstd is built into. COFF has no visibility;
  a DLL exports nothing without `dllexport`.
- **Install.** `asm_registry.install` points each slot at `@extern` of the symbol. A core with no
  object in the build is served by another (`served`, below). Every kernel is keyed by its
  `Slot`; a rule without a vtable field (`sdsdot_k`) is not converted, nor is one whose slot no
  core would take it for (`tools/blas_asm/roles.zig`'s `skip`, below).
- **Isolation.** The objects compile `_kernel.zig` alone, so everything they reach imports
  nothing outside `src/linalg/blas/`: `Vtable.zig`, `types.zig`, `util.zig`, `fma.zig`,
  `vector.zig`, `cpu/_coretypes.zig`, the registries and the kernels. `kernels/__root__.zig`
  repeats `cpu.GetCpuArch`'s choice of registry for that reason; the two are kept in step.

### Kernel kinds

A slot holds one of six kinds of kernel:

| Kind | Where | Built into |
| --- | --- | --- |
| Portable | each routine's `__root__.zig` (`Portable(T)`), `blas1/extremum.zig`, `blas2/symmetric.zig`, `ext/bf16/kernel.zig`, `blas3/gemm/kernel.zig` | main library |
| Converted assembly | `<routine>/<arch>/asm/*.zig`, generated | per-core object |
| Ported C driver around a converted microkernel | `<routine>/kernel.zig` (driver), `<routine>/<arch>/uk/*.zig` (generated blocks), `blas1/uk_block.zig` | per-core object |
| Ported intrinsics or vector C | `<routine>/<arch>/kernel.zig` and siblings (`axpy/arm64/sve.zig`, `axpy/riscv64/rvv.zig`, `axpy/mips/msa.zig`), shared drivers (`blas2/symv/kernel.zig`), `vector.ported` | per-core object |
| `@Vector` level 1 and 2 | `vector.zig` | per-core object |
| `@Vector` level-3 tile | `blas3/gemm/kernel.zig`, `blas3/symm/kernel.zig`, `blas3/trmm/kernel.zig`, `blas3/trsm/kernel.zig` | main library and every per-core object |

A converted source may call another routine's kernel (upstream, `loongarch64`'s `cnrm2`/`snrm2`
in `lsx` and `lasx` call `camax_k`/`samax_k`; those are skipped now). The generator preprocesses
with `-DDYNAMIC_ARCH`, so the call names one core's kernel (`camax_k_LA264`) rather than the
dispatch trampoline, and `arch.nakedLines` rewrites that to `blaskern_camax_la264`, which the
same per-core object exports.

[`blas1/c_kernels.md`](blas1/c_kernels.md) specifies the C-driver ports: the parts, what the
generator converts, the driver rules (order, contraction, mixed types, block ABI) and how to add
a routine.

### `@Vector` level-1 and level-2 kernels

`vector.zig` holds plain `@Vector` kernels, `U` native vectors a step (four in every object),
for `axpy`, `axpyc`, `dot`, `dotu`, `dotc`, `scal`, `nrm2`, `asum`, `sum`, the extremes, `copy`,
`swap`, `rot` and every `gemv` form, over `s`, `d`, `c` and `z`. Unit strides take the vectors;
other strides a scalar loop. Special values, strides and flags are the portable kernels'.
`nrm2` sums plain squares, and for `d` and `z` takes Blue's algorithm only where they could have
overflowed or underflowed. They need no FMA unit. `vector.Registry` exports them from the
per-core objects; like the tiles, they never make a core need an object of its own.

- **Routed slots** (`vector.routed`): the extremes, `nrm2`, `sdot`, and on x86_64 the reductions
  whose upstream assembly peels elements until an operand is aligned (`vector.peels`). No core
  takes a converted kernel there, so `tools/blas_asm` does not convert one (`roles.zig`'s `skip`,
  below; `asm_registry.Kernels` stops the build if a routed slot is converted), and an
  architecture's registry puts the `@Vector` kernels in where its vector unit is IEEE, else the
  portable ones serve. Where x86_64 routes a peeling slot by core, the generated index lists
  the cores of the assembly it left out (`skipped`).
- **x86_64, x86, loongarch64, mips, sparc, power:** the routed slots (on x86_64 the peeling
  ones only on the cores whose converted kernels peel; their C kernels elsewhere stay). The
  AltiVec cores' objects are built without AltiVec, so theirs are scalar.
- **arm:** the routed slots, but `s` and `c` on the NEON cores (`cortexa9`, `cortexa15`), whose
  lanes flush `f32` (`vector.flushing`).
- **arm64:** the routed slots on every core; on `vortex` and `vortexm4` every slot but `rot`
  and `sum` (unmeasured), where they match or beat the converted assembly on the M6
  ([`results-m6.md`](../../../tools/blasbench/results-m6.md)), so their objects there leave those
  converted kernels out (`vector.displacing`); elsewhere every slot no converted kernel fills
  (below). On Apple's systems, where Accelerate serves most routines, only the routed slots.
- **riscv64:** on the vector cores (`riscv64_zvl128b`, `riscv64_zvl256b`, `x280`) every slot of
  theirs no port fills (`Rule.unfilled_on`); `axpy` keeps its `asm` port. Each core's object
  builds them at its own guaranteed VLEN (`zvl128b`, `zvl256b`, `zvl512b`), where upstream's
  `vsetvl` loops read VLEN at run time; on a wider machine they use part of each register, which
  T20's race weighs. The cores without the ratified vector extension have no object.
- **s390x, wasm:** none.

`vector.ported` holds `@Vector` kernels no `Registry` routes, which an architecture's hand ports
(`.written = .vector`) name where upstream has intrinsics for them: `axpby`, `rotm`, `dsdot`,
`symv`, `hemv`, and the half-precision `gemv`s in the portable kernels' bits (`vector.Half`).
riscv64 uses them; `installAll` puts them in the native tests' table too.

Without a vector unit for the type (SPARC, ARMv7's `f64`, pre-VSX POWER, MIPS without MSA, the
x87) LLVM splits the vectors into scalar operations. The reductions keep their eight independent
sums, against the portable loop's one, and the extremes compare and select where `@max` would
call libm's `fmax` per element (`vector.hasMinMax`, MSA included); [`comparisons.md`](comparisons.md)
(1) has the measurements. On MSA, `@abs` clears sign bits through a mask LLVM cannot fold back
into `fmax_a` ([`issues.md`](issues.md)).

On arm64 these fill the slots whose upstream kernel is C, which the portable scalar loop served
before; their converted kernels stay, unmeasured:

| Cores | Slots with no converted kernel |
| --- | --- |
| all but `tsv110` | `sdot` |
| all | `cdotu`, `cdotc`, `zdotu`, `zdotc` (on `armv8`, `cortexa53`, `cortexa55`, `cortexa57`, `cortexa72`, `cortexa73`, `cortexa76` and `cortexx1` upstream's non-PGI builds take `zdot.S`, but a slot a make conditional leaves two sources for is not converted, `tools/blas_asm/gen.zig:589`) |
| `thunderx` | `ddot` |
| `thunderx2t99`, `thunderx3t110` | `ddot`, `?nrm2`, `?asum`, `?copy`, `i?amax` |
| `neoversen1`, `ampere1`, `vortex`, `vortexm4` | `ddot`, `?asum`, `?copy`, `i?amax`, `sgemv_n`, `ssum`, `dsum`, `zsum` |
| `neoversen2`, `neoversev2` | `ddot`, `?asum`, `?copy`, `i?amax`, `sgemv_n`, `sgemv_t`, `dgemv_n`, `dgemv_t` |
| SVE: `neoversev1`, `armv8sve`, `a64fx`, `armv9sme` | those, and `sswap`, `dswap`, `srot`, `drot`, `zsum` |
| SVE: `cortexa510`, `cortexa710`, `cortexx2` | the same |

### Coverage by architecture

The registries (`kernels/<arch>.zig`) list the routines with per-core kernels; the routed slots
above are not among them:

| Architecture | Level 1 | Level 2 | Level 3 | `@Vector` |
| --- | --- | --- | --- | --- |
| x86_64 | `axpy`, `dot`, `rot`, `scal` (C drivers and assembly); `asum` (`@Vector` ports), `copy`, `sum`, `swap` | `gemv` (C drivers, assembly, Skylake-X's `sgemv_n`), `symv` (assembly, `@Vector` ports), `hemv`, `ger`, `sbgemv` (`@Vector` ports) | tile | routed |
| x86 | `axpy`, `asum`, `copy`, `dot`, `rot`, `scal`, `sum`, `swap` | `gemv` | tile | routed |
| arm | `axpy`, `asum`, `dot`, `rot`, `sum`, `swap` | `gemv` | tile | routed, `d` and `z` on NEON |
| arm64 | `axpy` (assembly, SVE, ThunderX), `asum`, `copy`, `dot`, `rot`, `scal`, `sum`, `swap` | `gemv`; `symv`, `sbgemv`, `bgemv` (`@Vector` ports) | tile | routed and more |
| loongarch64 | `axpy`, `asum`, `axpby`, `copy`, `dot`, `rot`, `scal`, `sum`, `swap` | `gemv`, `symv` | tile | routed |
| mips, mips64 | `axpy` (MSA on 32- and 64-bit), `asum`, `copy`, `dot`, `scal`, `sum`, `swap` | | tile | routed |
| power | `axpy` (assembly, microkernels, C), `asum`, `copy`, `dot`, `rot`, `scal`, `sum`, `swap` | `gemv` | tile | routed |
| sparc | `axpy`, `asum`, `copy`, `dot`, `rot`, `sum`, `swap` | `gemv` | | routed |
| riscv64 | `axpy` (RVV), `axpby`, `dsdot`, `rotm` (`@Vector` ports) | `symv`, `hemv`, `sbgemv`, `shgemv` (`@Vector` ports) | tile | vector cores |
| s390x | `axpy` (microkernels) | | tile | |
| wasm | `axpy` | | tile | |

Everywhere else, and for every routine not listed, the portable kernels serve: `ger`, the
triangular, banded and packed level-2 routines (drivers over level-1 kernels), the bfloat16 and
half-precision kernels (but the half `gemv` ports above), and the `q`/`x` slots,
which no core replaces (a test in `__root__.zig` checks this). SPARC exports no level-3 tiles: no SPARC feature changes them, so an object would
only repeat the portable kernels.

### Level 2

- **Interface.** `blas2/interface.zig` holds the checks, `beta` scaling and pointer set-up of
  `interface/{gemv,zgemv,ger,zger,symv,zsymv,zhemv}.c`. `triangular.zig` and `banded.zig` hold
  those of the triangular, banded, packed and rank-update routines with their
  `driver/level2/*.c` loops, which run over the level-1 and `gemv` kernels (`portable.kernel`).
  The complex symmetric forms (`zsbmv`, `zspmv`, `zspr`, `zspr2`, `zsyr`, `zsyr2`) are OpenBLAS
  extensions and are included.
- **Rank updates.** The rank-update drivers skip a column whose factor is zero, as reference
  BLAS does (OpenBLAS's `zher_k.c`, `zhpr_k.c` and the rank-2 drivers do not). `ger` calls the
  kernel once per run of nonzero `y`, which is one call when `y` has no zero, and copies a
  strided `x` once for all the runs. The skip is decided above the kernels, so every core
  gives it; [`comparisons.md`](comparisons.md) (5) has why it stays. At unit stride and small
  `n` (`syr`, `spr`, `syr2` below 100; `spr2`, `zsyr` below 50: `interface/syr.c:174`,
  `spr.c:172`, `syr2.c:175`, `spr2.c:173`, `zsyr.c:175`) OpenBLAS's interfaces run the driver's
  columns themselves; rstd's driver takes a unit-stride vector as it is, so it is that path,
  keeping the skip (OpenBLAS's `syr2`/`spr2` paths have none).
- **symv and hemv.** `blas2/symmetric.zig` holds the portable kernels, `kernel/generic`'s
  `symv_k.c`, `zsymv_k.c` and `zhemv_k.c`, which compose the running core's copy and `gemv`
  kernels with the same blocking (`SYMV_P`, the vtable's `symv_p`). `blas2/symv/kernel.zig`
  ports the real C drivers of x86_64 (`?symv_[LU].c`) and arm64 (`symv_[LU]_asimd_4x4.c`),
  their 4x4 blocks on `@Vector`s in each core's widths, fusing and reductions (`Shape`).
  arm64's SVE drivers (`symv_[LU]_sve_v1x4.c`) are not ported.
- **Buffers.** The kernels take a scratch `buffer`. The interfaces size it as the reference's
  `STACK_ALLOC` does, on the stack up to `interface.stack_bytes` (2048, 32-aligned, the
  reference's `MAX_STACK_ALLOC`: `Makefile.system:573`, `common_stackalloc.h:58`), else from
  the allocator, 4096-aligned (the kernels round to 4096-byte boundaries in it, `kernel/generic/symv_k.c:51`).
  [`blas2/buffers.md`](blas2/buffers.md) lists what every kernel does with it.
- **Threading.** `gemv`, `ger`, the rank updates and transposed `gbmv` split `y` or `A`'s
  columns into serial kernel calls. The symmetric, triangular, banded and packed products and
  untransposed `gbmv` split their columns into chunks fixed by the size, each into a partial
  vector of its own, folded in chunk order; the triangular ones port the `*_thread.c` kernels
  ([`threading.md`](threading.md#levels-1-and-2-fixed-chunks-dynamic-claiming-ordered-combine)).

### Level 3

`gemm`, `symm`, `hemm`, `syrk`, `herk`, `syr2k`, `her2k`, `trmm` and `trsm` are ported, real and
complex, for every float type, through the same code.

- **Interface.** `blas3/interface.zig` checks the arguments in `xerbla` order, returns for empty
  problems, and calls the drivers with the table for `T` (a type without slots never reads it).
- **Drivers.** `blas3/driver.zig` is `driver/level3/level3.c`: `C` in `R`-wide column blocks,
  the depth in `Q`-deep slices, the rows in `P`-tall blocks; each slice of `op(B)` is packed
  once (`sb`), each block of `op(A)` once per slice (`sa`). One set of loops serves every
  product of two packed operands (`Desc`): `gemm`; `symm` and `hemm`, which pack the symmetric
  operand whole from its triangle (`symm_k.c`, `zhemm_k.c`); and the rank-k updates, which keep
  one triangle of `C` (`level3_syrk.c`, `level3_syr2k.c`). Tiles across the diagonal go through
  a scratch tile and add the triangle's part (`syrk_kernel.c`); the reference's shared
  `sa`/`sb` path, which needs `unroll_m == unroll_n`, is not taken. `trmm` and `trsm` follow
  `trmm_L.c`, `trmm_R.c`, `trsm_L.c` and `trsm_R.c` loop for loop.
- **`beta`.** Each block of `C` is scaled by `beta` just before the first depth slice first
  touches it, while it is in cache, not in a pass of its own (BLIS's first-`KC` `beta`,
  `frame/3/gemm/bli_gemm_blk_var3.c:92`); `trmm`'s
  and `trsm`'s `alpha` goes where each part of `B` is first read. The same multiplies as a pass
  first, so the same bits; a zero `beta` writes zeros without reading `C`. `alpha == 0` and
  `k == 0` keep the separate pass.
- **Chunks.** The product drivers call the kernel on the same `3 * unroll_n` column chunks they
  pack (`Blocking.chunk`), where the reference hands the kernel a slice's whole `sb`. A stated
  unroll that does not match the kernels then costs speed, never correctness. The triangular
  drivers take whole slices, as the reference's do, and rely on the vtable's unrolls being the
  installed kernels' own.
- **Blocking.** `Blocking.of` reads `P`, `Q`, `R` and the unrolls from the vtable, falling back
  to 224, 128 and 4096 where a rule left a value unset or negative (`Q` is the generic table's,
  `params_data.generic`; `P` and `R` are ours, since that table leaves them to `parameter.c`).
  `Blocking.split` is `level3.c`'s split of the remaining rows or depth.
- **Packing space.** `Work(T)` sizes `sa`, `sb` and the scratch tile for the problem, each on a
  64-byte line; up to `driver.stack_bytes` (2048, ours, as `interface.stack_bytes`) it is on the stack, else from the caller's
  allocator.
- **Kernels.** `blas3/gemm/kernel.zig` holds the packing copies, `beta` and the register-tile
  kernel for real `T` and its `Complex` (`incopy`, `itcopy`, `oncopy`, `otcopy`, `beta`,
  `kernel`, complex `kernel_n`, `_l`, `_r`, `_b`). `symm/kernel.zig`, `trmm/kernel.zig` and
  `trsm/kernel.zig` hold the other routines' copies and kernels over the same tile (`sets`),
  each function named as its slot less the precision letter. The `trsm` kernels are
  `kernel/generic`'s `trsm_kernel_*.c` on the gemm layout's panels, dividing by the diagonal
  where upstream packs its reciprocal (`INV`, `kernel/generic/trsm_ltcopy_2.c:43`), which
  overflows for a denormal below about 2^-128 in `f32` or 2^-1024 in `f64` (past `1/max`); `trsv` divides
  alike, complex by Smith's algorithm (`util.cdiv`; R. L. Smith, "Algorithm 116: Complex
  division", Communications of the ACM 5(8), 435, 1962, <https://doi.org/10.1145/368637.368661>), so both give the same quotients
  (`comparisons.md`, 10), but for an x87 build's cores with SSE2, whose objects divide there. The `trmm` kernels sum a
  diagonal block over its triangle only, where the reference's tiles multiply the zeros it
  packs, so an infinity in `B` meets no zero.
- **Packed layout.** `n / u` panels of width `u`, then one of each power of two below `u`
  present in `n % u`; panel `p` at column `c0` with width `w` sits at `m * c0`. For a
  power-of-two `u` this is `kernel/generic`'s layout. A complex element is its two parts, real
  first.
- **Tile.** `geometryFor(T, cpu)` sets the tile from the target's vector width
  (`vectorBits`) and register count: `mv` vectors of `vl` lanes by `nr` columns, the most
  accumulators that fit beside the `mv` loads of A, B (one broadcast, or `nr / vl` vectors where
  LLVM multiplies by lane) and one spare register, with no more loads per step than three
  quarters of the multiply-adds. Complex lanes hold whole elements with two accumulators each,
  after `zgemm_kernel_4x2_haswell.S`; without vectors, or with `f32` lanes that flush, the parts
  are separate scalars. Row tails are unrolled tiles, one per pair of widths, as in OpenBLAS.
- **One source for every core.** The tile kernels are built into the main library as the
  portable kernels (`installGeneric`) and into every per-core object (`Registry`), where the
  object's features set the tile. They never make a core need an object of its own
  (`Registry.has` is false). The objects export their unrolls as data
  (`blaskern_<p>gemm_unroll_<m|n>_<core>`); `cpu/<arch>.zig` sets the core's `param.h` blocking
  before its kernels install, and `portable.finish` restores the portable unrolls where the
  portable kernels stayed, so the vtable's unrolls are always the installed kernels'. The `q`
  and `x` slots keep the portable kernels (`core_precisions`).
- **Parameters.** `blas3/gemm/params.zig` fills the blocking as `kernel/setparam-ref.c` does,
  from `params_data.zig` (generated from `param.h`): `p` of the older x86 cores follows the L2
  size, x86's `r` derives from the packing buffer and its `p` rounds up to the unroll (elsewhere
  `Blocking.of` takes whole panels: loongson3r4's `dgemm` `p` of 44 runs as 40), x86's 3M variants
  take the real ones' `p` and `q` unless the header says otherwise (no other `DYNAMIC_ARCH` build
  sets 3M), and LA464 and the AMD AVX-512 parts pick their sets by cache size and thread count.
  `OPENBLAS_BLOCK_FACTOR` scales a computed `p` as `driver/others/parameter.c:519-527` does;
  OpenBLAS ignores it in dynamic builds, which do not link `parameter.c`
  (`driver/others/Makefile:44`).
  Haswell's and Zen's `dgemm` `q` is 128 under `WINDOWS_ABI`, 256 otherwise (`param.h:1613-1619`,
  `:704-710`); `params_data.zig` takes 256 on every OS. The cut is for the stack, not tuning:
  `dgemm_kernel_4x8_haswell.S:54-70` keeps a `q * 8 * 12`-byte buffer on the stack and touches
  its pages on Windows (`:82-106`), and 47b22763f ("reduced stack usage on windows to 16K")
  halved `q` to hold it near 16 KiB; Zen took Haswell's values in 832a27278. The Zig kernels
  keep no such buffer.
  The packing buffer each core's `r` derives from is its `BUFFER_SIZE` (`common_<arch>.h`),
  preprocessed with `-DDYNAMIC_ARCH` like the kernels (Preprocessing, below), so every POWER
  core takes the dynamic build's 256 MiB (`common_power.h:923-926`), not its static size.
  The generic table (the portable vtable on every target, and the generic family, Zig's
  architectures with no detector) takes the target's own `BUFFER_SIZE`, as OpenBLAS's generic
  target sees its architecture's header: 16 MiB on x86 (`common_x86.h:207`), 32 MiB on arm and
  arm64 (`common_arm.h:141`, `common_arm64.h:183`), 256 MiB on POWER, and so on (T112). Its
  extracted 128 MiB (`common_x86_64.h:260`) was right on x86_64 only. Alpha and C-SKY have
  headers but no family (32 MiB, `common_alpha.h:128`, `common_csky.h:53`); architectures
  OpenBLAS has no port for keep the x86_64 stand-in.
  Only x86 sets the generic `p` and `r` at run time (`parameter.c:451`'s L2 rule, `r` from the
  buffer). A static `TARGET=GENERIC` build elsewhere keeps the header variables at
  `parameter.c:46-48`'s 128 (`p`, `q` and `r`), and MIPS and LoongArch have fixed sets
  (`param.h:4388-4415`). The generic table follows its target's branch (`GenericBranch`,
  `params_data.generic_off_x86`), and so does a core with no table of its own (T116).
  `dtb_entries` per core is `getarch.c`'s `DTB_DEFAULT_ENTRIES`, which `params.zig` also
  extracts into `params_data.zig` (`dtbFor`). A 32-bit family takes `Makefile.system`'s 32-bit
  fallback entries (`Makefile.system:115-158`) (x86's AVX cores Nehalem's, the Bulldozer family and Zen Barcelona's, POWER8
  to 10 POWER6's), as do the 32-bit x86 kernels' `config.h` macros; POWER's tables also follow
  the byte order (`PPC970`'s unrolls).
- **3M.** `gemm3m`, `symm3m` and `hemm3m` (`blas3/gemm3m/`) port `gemm3m_level3.c`,
  `symm3m_k.c` and `hemm3m_k.c`: per depth slice three passes, each packing a real form of
  `op(A)` (both parts summed, the real part, the imaginary part) and of `alpha * op(B)`, then a
  real kernel adding the product to both parts of `C` with a sign each (`ALPHA5` to
  `ALPHA18`). The kernel is the real tile of the parts' type, writing complex `C` (the Haswell
  kernel's lane doubling); the copies are `kernel/generic`'s 3M copies over the gemm layout.
  Both join the gemm tile's `sets`, so every core's object builds them, with their unrolls
  (`<p>gemm3m_unroll_m`, `_n`) exported as the gemm ones are. `P`, `Q`, `R` are the 3M
  parameters where x86's rule sets them, else the real precision's; depth slices split as
  `gemm3m_level3.c` does (halves, not rounded to the unroll). Threads take parts of `C` on whole
  panels, each a serial problem, so the bits are the serial ones. Opt-in, as upstream: one part
  of each element comes from a difference of products and is only normwise stable
  ([`comparisons.md`](comparisons.md), 7); never routed to Accelerate, which has no 3M. Upstream
  exports only `gemm3m`, and only on x86, x86_64, ia64 and MIPS (`interface/Makefile`); here all
  three serve every target.
- **Small matrices.** `gemm` takes OpenBLAS's `SMALL_MATRIX_OPT` path (`blas3/gemm/small.zig`)
  where the core's `gemm_small_matrix_permit` allows: no packing, each element summed over the
  depth in order, then `beta*C + alpha*sum` (`gemm_small_kernel_b0_*` for a zero `beta`, which
  never reads `C`). The kernels are `kernel/generic`'s, the real ones on `@Vector`s of the gemm
  tile, in every core's object; the per-core assembly and intrinsics (Skylake-X, SVE, POWER10,
  LoongArch) are not converted. The permits are upstream's per core (`permitFor`, installed with
  the blocking): Skylake-X and its successors, the `ARMV8SVE` family, POWER10, LA264 and LA464;
  every other core and all complex take the generic permit, which allows nothing. POWER10's
  `sgemm` permit follows the thread count upstream (64^3 with more than one, 100^3 with one);
  here it is 64^3 always, so the path never does. Taken after `driver.gemm`'s quick returns,
  where upstream takes it before them, so a zero `alpha` or `k` still reads neither operand.
  The thresholds are tuned for upstream's kernels, not these (T05).
- **Not ported.** `sbgemm`'s small-matrix kernels, and `level3_thread.c`: level 3 threads by
  BLIS's ways instead ([`threading.md`](threading.md)).

[`blas3/gemm/notes.md`](blas3/gemm/notes.md) has the tile measurements and code generation notes
behind these choices; [`comparisons.md`](comparisons.md) (7, 9, 10) the open performance and
robustness questions (fusing `beta`, `trsm`'s inverted diagonal, prefetch).

### bfloat16 and half precision

The kernels are `ext/bf16/kernel.zig`: the conversions, `bscal`, `sbdot`, the
`bgemv`/`sbgemv`/`shgemv` family and the `bgemm`/`sbgemm`/`shgemm` one. They are portable Zig,
the same on every core in their bits. The cores whose reference has a half `gemv` of its own
install `vector.Half`'s, the portable kernels' bits on vectors: `riscv64_zvl256b` (`sbgemv`,
`shgemv`), Cooper Lake and Sapphire Rapids (`sbgemv_n`), Neoverse N2, V1 and V2 (`sbgemv`,
`bgemv`). Upstream's AVX-512 BF16, BFDOT and BFMLAL sum pairs of products in instructions with
their own rounding, which `@Vector` cannot express and the contract does not want, so only
their shapes carry over; VSX and MMA are not ported. Every bfloat16 conversion
upstream is `kernel/x86_64/tobf16.c` and `bf16to.c` (`kernel/Makefile.L1:311`, `:315`); the AVX-512 BF16 microkernels give the same
bits.

`kernel/setparam-ref.c` fills the bfloat16-named slots that have no kernel of their own with the
`f32` kernels: 17 level-1 ones (`sbamax_k` to `sbswap_k`) and `dsbdot_k` (with `dsdot_k`;
positionally, `setparam-ref.c:107-111` against `common_param.h:86-108`), which
the level-1 interfaces reach, and `sbtrsm_*`, `sbtrmm_*`, `sbsymm_*`, `sbneg_tcopy` and
`sblaswp_ncopy`, which nothing calls yet. `portable.finish` sets them at the end of every
vtable's build.

The half gemms (`blas3/half.zig`) are `interface/gemm.c` built `BFLOAT16`, `BGEMM` and
`HFLOAT16` over `level3.c`'s loops, for one type of `A` and `B` and another of `C`; the
`driver.zig` loops take one type, so `half.zig` has its own copy of them, which uses `driver.zig`'s
`Blocking` and `Work`. The kernels port `kernel/generic/gemmkernel_2x2.c`, `gemm_beta.c` and the
`ncopy`/`tcopy` packing (`Makefile.L3:190-243`) over the `f32` gemm tile; each element's
arithmetic is the same in every lane and tile, so threads split `C` by columns (rows for a tall
one) into problems of their own, each on its own packing space, and get the serial bits. Blocking
is the core's `sb`, `b` and `sh` sizes (`params.zig`), with `Blocking.of`'s fallbacks where its
family sets none; the unrolls are the portable kernels' (`settleGemm`), and `align_k` is unused,
since the portable packing needs no padding. `gemmtHalf` is `interface/sbgemmt.c`: a half `gemv`
per column of the triangle, serial, for all three types (upstream exports `sbgemmt` only).
Upstream compares the bfloat16 scalars by their bits; rstd by value ([`issues.md`](issues.md)).
The `sbgemm` small-matrix slots stay null (T30's small-matrix path), and the AVX-512 BF16, AMX,
NEON, SVE and MMA kernels are T29 and T31.

Their numerical behavior is in [`readme.md`](readme.md#bfloat16-and-half-precision).

### Extensions

OpenBLAS's extensions, column-major only, behind `ext/interface.zig` (`geadd`, `omatcopy`,
`imatcopy`, `laswp`) and `blas3/interface.zig` (`gemmt`). Every core runs the portable kernels:
upstream converts no assembly for these slots, only C and intrinsics (NEON, RVV, LASX, vector).

- **`geadd`** (`kernel/generic/geadd.c`, `zgeadd.c`): `axpby` per column of `op(C)`, or `scal`
  with flag 0 for a zero `alpha`, through the vtable, so it follows the core's `axpby` and
  `scal`. The letters come last, as in the Fortran; `op` never conjugates, so complex `r` and
  `c` are `InvalidTrans` where upstream takes `C` as `T` (`interface/zgeadd.c:65-66`).
- **`omatcopy`, `imatcopy`** (`ext/matcopy.zig`, after `kernel/arm/*omatcopy_c*.c` and
  `kernel/generic/*imatcopy_c*.c`): one scaling for all, a zero `alpha` storing zeros and
  `alpha` 1 copying, complex too (upstream's complex kernels multiply by 1, `issues.md`).
  `imatcopy` runs the in-place kernels where the shape and leading dimension stay, and else goes
  through a buffer from its allocator, as upstream does (`interface/imatcopy.c:128-164`, under
  `NEW_IMATCOPY`, which line 55 always defines); both give the same bits. No `q` or `x` slots, as upstream.
- **`laswp`** (`lapack/laswp/generic/laswp_k_1.c`): the interchanges as `dlaswp.f` defines
  them, for any element type, not a vtable slot upstream either. Rstd adds checks the reference
  has none of: rows and pivots inside `A`, all pivots checked before the first swap.
- **`neg_tcopy`, `laswp_ncopy`**: gemm packing copies (`neg_tcopy_*.c` is `gemm_tcopy_*.c`
  negated; `laswp_ncopy_*.c` swaps and packs as `gemm_oncopy`), so they join the gemm kernels'
  `sets` and take each core's tile. `laswp_ncopy` also swaps rows `k1` to `k2` in `A`, which
  upstream leaves stale, as `getrf_single.c`'s `laswp_plus` and `gemm_oncopy` do in its place.
  Nothing calls either: LAPACK's `getrf` takes `laswp_plus` and `gemm_oncopy`, as
  `getrf_single.c`'s `#if 1` does; `getrf_parallel.c`, its only caller, has the call behind
  `if (0 && ...)` (`:138`, `:277`), and nothing upstream calls `neg_tcopy`.
- **`gemmt`**: upstream takes a `gemv` per column, and for `r` and `c` conjugates the caller's
  `B` in place and back (`interface/gemmt.c:495`, `:697`, the `gemv` at `:569`; `issues.md`). Here it is `gemm`'s operands over `syrk`'s triangle in
  the level-3 driver (`driver.gemmt`), and `B` is only read.
- **Batched `gemm`** (`blas3/batch.zig`; `interface/gemm_batch.c`, `gemm_batch_strided.c`):
  groups of problems sharing sizes and scalars (a `GemmGroup` holds slices of operand pointers;
  unequal lengths are `InvalidCount`), or one problem at fixed strides. Every group is checked
  before any runs, as upstream. Upstream hands whole problems to threads and takes its
  small-matrix kernels up to `100^3`; here they run in order through `gemm` or `gemmHalf`, each
  threading as usual, with the same bits. Upstream has `sbgemm`'s batches; `bgemm`'s and
  `shgemm`'s come with the same code.
- **Accelerate**: none of these is routed; Accelerate's CBLAS has no column-major equivalent of
  them (its `appleblas_?geadd` takes three matrices, `const` `A` and `B` and a written `C`;
  <https://developer.apple.com/documentation/accelerate/appleblas_dgeadd(_:_:_:_:_:_:_:_:_:_:_:_:_:)>).

### LAPACK

`src/linalg/lapack/` (`rstd.linalg.lapack`) holds LAPACK routines over this BLAS, with its
conventions: `(allocator, io, T, ...)`, column-major, the letter enums, `xerbla`-named argument
errors in the reference's order. OpenBLAS's own LAPACK (`lapack/`, `interface/lapack/`) is the
source where it has a routine; past it, faer's (`faer/`), as the anthology decision (Settled
decisions) chooses. Reference-LAPACK (`lapack-netlib/SRC/`) gives the semantics: where a
source differs, the reference wins, and the file says so. LAPACK's per-core kernels
(`laswp_ncopy`, `neg_tcopy`) are vtable slots, so they stay in `blas/` with the other kernels;
the routines live here and call the BLAS's drivers and table directly, which one Zig module
allows.

- **`potrf`** (`cholesky.zig`): `potrf_U_single.c` and `potrf_L_single.c`, the recursive
  blocked factor: blocks of `GEMM_Q` columns, or a quarter of `n` up to `4 * GEMM_Q`
  (`lapack/potrf/potrf_U_single.c:113-114`; `:108` for the `dtb_entries / 2` case), each
  factored recursively, then `trsm` on its panel and `herk` (`syrk`) on the rest through the
  level-3 drivers, which upstream fuses per column block. Up to `dtb_entries / 2` columns it is
  `potf2_[UL].c` and `zpotf2_[UL].c` as plain loops. Two departures, both the reference's
  `xPOTRF2`: a NaN pivot fails (`DISNAN`, `lapack-netlib/SRC/dpotrf2.f:168`), where upstream `ajj <= 0`
  (`lapack/potf2/potf2_U.c:71`) passes it and reports success; and
  the base case divides by the pivot, as `trsm` does, where upstream multiplies by its
  reciprocal (`potf2_U.c:86`). The drivers are called directly, so Accelerate never runs them;
  on macOS the whole `potrf` goes to Accelerate's `?potrf_` while routing is on
  ([Apple Accelerate](#apple-accelerate)). `potf2` is the base case, public.
- **Threads.** Each interface's thread rule is upstream's (`cholesky.zig`'s `threads`):
  `getrf` from `m * n` 40000 (`DOUBLE` 10000; complex 10000), `gesv` from `n * n` alike,
  `potrf` from `n` 128 (`DOUBLE` and complex 64), `lauum`, `getrs` and `trtrs` from 64,
  `trtri` past 150 (complex from 64), `potri` from 180 (complex 200, `DOUBLE` 150). Upstream
  runs its `_parallel.c` algorithm when it has more than one thread; here the size alone
  picks it, so one worker past the threshold runs it too, and the drivers' workers keep the
  bits ([`threading.md`](threading.md)). No result depends on the worker count.
- **`potrf` threaded** (`potrf_U_parallel.c`, `potrf_L_parallel.c`): blocks of half `n` rounded
  up to `unroll_n`, at most `GEMM_Q` (`:92`), each factored recursively until `4 * unroll_n`
  is left for the single form; the panel's `trsm` and the update's `herk` take the workers.
- **`potri`** (`cholesky.zig`, `interface/lapack/potri.c`): `trtri` (non-unit), then `lauum`,
  both single or both parallel by `potri`'s rule. A zero on the factor's diagonal is
  `error.Singular` before `A` is touched, as `xPOTRI` has it through `xTRTRI`; upstream skips
  that check and divides (`issues.md`).
- **`trtri`, `trti2`** (`triangular.zig`): `trtri_U_single.c` and `trtri_L_single.c`, blocks of
  `GEMM_Q` (upper, `trti2` up to `GEMM_Q`; lower, below it), each multiplied by the inverse
  already made (`trmm`), solved by its diagonal block with -1 (`trsm`), then inverted by
  `trti2_[UL].c` and `ztrti2_[UL].c` as loops (`trmv`, then `scal` by `-1 / A(j, j)`; complex
  reciprocals by `ztrti2_U.c`'s ratio). Threaded, `trtri_[UL]_parallel.c`: `trti2` up to
  `dtb_entries` (`:87`), blocks of `GEMM_Q` or a quarter of `n` below `4 * GEMM_Q` (`:93`),
  right-looking: `trsm`, the block recursively, then `gemm` and `trmm` on the rest. A zero on
  the diagonal (non-unit) is `error.Singular` before `A` is touched, found by `==`, where
  upstream asks `amin` and `iamin`. `trti2` has no check, as the reference's.
- **`lauum`, `lauu2`** (`triangular.zig`): `lauum_U_single.c` and `lauum_L_single.c`, blocks of
  `GEMM_Q` or a quarter of `n` up to `4 * GEMM_Q`, `lauu2` up to `dtb_entries`; each block adds
  its off-diagonal part's `herk` (`syrk`) into the part already made, multiplies that part by
  its diagonal block (`trmm`), then recurses, which upstream fuses per column block.
  `lauu2_[UL].c` and `zlauu2_[UL].c` are loops in upstream's order (`scal`, then `dot` and
  `gemv`), not the reference's (`dot` and `gemv` with `beta`). Threaded, `lauum_[UL]_parallel.c`:
  blocks of half `n` rounded up to `unroll_n`, at most `GEMM_Q` (`:91`), down to `2 * unroll_n`.
- **`trtrs`** (`triangular.zig`): `trtrs_single.c` and `ztrtrs_single.c`, `trsm` through the
  driver, `trsv` for one right-hand side; `trtrs_parallel.c` is that `trsm` on threads. It
  takes `trans` `r`, as upstream does, and checks its arguments in the reference's order
  (upstream reports `diag` first, `issues.md`). A zero on the diagonal is `error.Singular`
  before `B` is touched, even with no right-hand side, as the reference.
- **`potrs`** (`cholesky.zig`): `xPOTRS`, which OpenBLAS takes from the reference too
  (`exports/gensymbol.pl:498` lists it with the netlib objects): two
  `trsm`s on the factor, through the drivers.
- **`getrf`** (`lu.zig`): `getrf_single.c`, the recursive blocked factor: blocks of half the
  columns rounded up to `unroll_n`, at most `GEMM_Q`, each factored recursively, then the
  block's interchanges (`laswp`), `trsm` on its rows and `gemm` on the rest through the level-3
  drivers, which upstream fuses per `GEMM_R` columns; the later interchanges reach each block's
  columns last. While that block is at most `2 * unroll_n` wide (`min(m, n)` up to about
  `4 * unroll_n`; `unroll_n` on POWER8: `lapack/getrf/getrf_single.c:76-88`) it is
  `getf2_k.c` and `zgetf2_k.c`, left-looking,
  as loops but for the core's `iamax`, so the pivot is the first NaN, else the first largest
  `|re| + |im|`. Two departures, both the reference's `xGETRF2`: a nonzero pivot below the
  safe minimum (denormal, or NaN) is swapped in and divided by, where upstream `getf2`
  records the swap in `ipiv` but neither swaps nor scales (`issues.md`); and the complex test
  is on `|p|`, where upstream tests each part. A normal pivot's column is scaled by its
  reciprocal, as both do. `laswp_ncopy` and `neg_tcopy` stay unused ([Extensions](#extensions)).
- **`getrf` threaded** (`lu.zig`, `getrf_parallel.c`): the lookahead. The first block is
  factored recursively; then for each block the caller updates the next block's columns
  (rounded up to whole `unroll_n` panels) and factors that block, while the other workers
  update the rest in chunks cut at panel edges (`inner_basic_thread`), each on packing space
  allocated before the tasks start (`threading.team`). Upstream sizes its blocks by its thread
  count (`FORMULA1`, `FORMULA2`, `:502`), which would change the bits; here they are
  `getrf_single.c`'s, and every column gets the serial kernel calls, so the factors are the
  serial ones. The closing interchanges, threaded upstream, are serial.
- **`getrs`** (`lu.zig`): `getrs_single.c` and `zgetrs_single.c`, which are `xGETRS`: the
  interchanges, then `trsm` with `L` and with `U` through the drivers, or for `op(A)`
  transposed the solves first and the interchanges backward; one real right-hand side takes
  `trsv` on the vtable's kernels, as `getrs_single.c:47` does. It takes `trans` `r`, as
  `interface/lapack/zgetrs.c:96` does, and checks that every pivot is a row of `B` (`error.InvalidIpiv`), as
  rstd's `laswp` does; the reference checks none. Threaded, its `trsm`s split the columns, as
  `getrs_parallel.c` gives each thread a column block.
- **`gesv`** (`lu.zig`, `interface/lapack/gesv.c`): `getrf`, then `getrs` unless `U` is
  singular, as `xGESV`, both on `gesv`'s thread rule.
- **faer's routines** (`faer/`): faer v0.24.4 (MIT,
  <https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`), behind Reference-LAPACK's
  interfaces. Each file
  ports one faer source and names it: `ops.zig` (a column-major `Mat` view, element
  arithmetic, and the BLAS calls that replace faer's own `matmul` and triangular kernels: level
  3 through the drivers on the interface's worker rule, level 1 and 2 through the public
  routines), `householder.zig` (`householder.rs` and `jacobi.rs`: reflections, their block
  factors `T = diag(tau) + striu(Vᴴ V)` and block application, real Givens rotations),
  `qr.zig` (`qr/no_pivoting/factor.rs`, faer's staircase QR, with `geqrf`, `orgqr`,
  `ormqr`), `qr_solve.zig` (`qr/no_pivoting/solve.rs`, with `gels`), `bidiag.zig`
  (`svd/bidiag.rs`), `bidiag_svd.zig` (`svd/bidiag_svd.rs`: QR iteration, divide and
  conquer, and the secular root finder the eigensolver shares), `svd.zig` (`svd/mod.rs`, with
  `gesvd` and `gesdd`), `tridiag.zig` (`evd/tridiag.rs`), `tridiag_evd.zig`
  (`evd/tridiag_evd.rs`), `hessenberg.zig` (`evd/hessenberg.rs`), `real_schur.zig` and
  `complex_schur.zig` (`evd/schur/`), `balance.zig` (Reference-LAPACK's `dgebal.f`,
  `dgebak.f`, which faer lacks), and `evd.zig` (`evd/mod.rs`, with `syevd`, `syev`
  and `geev`). Faer's `rayon` parallelism becomes rstd's: the fused bidiagonalization pass
  splits by column chunks (independent columns, so the same bits at any worker count); the
  fused tridiagonalization and Hessenberg passes, whose threads each sum a partial
  product, run serially (faer's `Par::Seq` path), so no sum depends on the worker count;
  the products thread inside the BLAS; and divide and conquer's two halves run in turn, as
  rstd's threading splits flat ranges, not task trees.
- **`gesvd`, `gesdd`** (`faer/svd.zig`): `dgesvd.f`'s and `dgesdd.f`'s arguments, jobs (`Job`:
  `a`, `s`, `o`, `n`), `xerbla` order and scaling of `A` into `[SMLNUM, BIGNUM]` (`dlascl.f`),
  over faer's SVD: QR first past an aspect ratio of `11 / 6`, then bidiagonalization, then the
  bidiagonal's SVD by QR iteration (`gesvd`, as the reference's `xBDSQR`) or divide and conquer
  from 128 up (`gesdd`, as its `xBDSDC`). No `WORK` arguments: the allocator gives workspace.
  Failure is an error: `NoConvergence` for the reference's `INFO > 0` (without its count, which
  faer does not keep) and for an infinite `A`, or a NaN for `gesvd`, since faer refuses a
  bidiagonal that is not finite; `gesdd` returns `InvalidA` for a NaN, the reference's
  `INFO = -4` (`dgesdd.f:595`). Four fixes to faer, each in [`issues.md`](issues.md): its
  secular solver can stop far from a root next to a pole, here bisected; reciprocals of
  denormal norms overflow, here divided by; a zero subproblem's two rows of `U` lack the
  identity's last 1; and its tall path leaves `U`'s top right block as found, here zeroed.
- **ReLAPACK** (`relapack.zig`; Elmar Peise's recursive algorithms, MIT, from OpenBLAS's
  `relapack/src/` at the same pin): each routine halves its problem by `split` (`?REC_SPLIT`, a
  64-byte cache line of elements) down to `CROSSOVER` (24), joining the halves with the level-3
  drivers on the interface's worker rules, so the bits are the same at any worker count.
  `potrf`, `getrf`, `trtri` and `lauum` replace OpenBLAS's only under `-Drelapack`
  (`build/blaskern.zig`'s `lapack_options`), as OpenBLAS's `BUILD_RELAPACK` with
  `RELAPACK_REPLACE`; the default stays OpenBLAS's. ReLAPACK serves `sygst`, `hegst`, `pbtrf`
  and `trsyl` by default; their unblocked kernels are Reference-LAPACK's (`xSYGS2`, `xHEGS2`,
  `xPBTF2`, `xTRSYL` with `xLALN2` and `xLASY2`), since faer has none in LAPACK's storage.
  Where ReLAPACK is wrong ([`issues.md`](issues.md#relapack-openblass-relapack)), the reference
  wins: `sytrf`, `hetrf` and their `_rook` forms run the reference's blocked `xSYTRF` loop over
  ReLAPACK's `rec2` panels (`indefinite.zig`), and `gbtrf` is `xGBTF2`.
- **`geqrf`, `orgqr`, `ormqr`** (`faer/qr.zig`; `ungqr`, `unmqr` the same functions):
  `dgeqrf.f`'s, `dorgqr.f`'s and `dormqr.f`'s arguments and `xerbla` order over faer's
  staircase QR. `tau` is LAPACK's (`H = I - tau v vᴴ`), the reciprocal of faer's, so an
  identity reflection is `tau = 0` with its `v` zeroed. Two departures, both faer's: a
  complex `R`'s diagonal is not made real (`tau` is real instead), and a column within
  `16 (m - j) eps` of the span of those before it takes no reflection (`Q R = A` holds to that
  tolerance). `orgqr` and `ormqr` rebuild the block factors from `V` and `1 / tau`, so the
  reference's complex `tau` works too; `ormqr` takes `trans` `t` and `c` alike for real
  (complex: `c` only, as `zunmqr`). A fix to faer ([`issues.md`](issues.md)): its staircase
  skips a column whose norm is NaN, here reflected, so NaN propagates as in the reference.
- **`gels`** (`faer/qr_solve.zig`): `dgels.f`'s arguments, scaling of `A` and `B`, its four
  cases and where the answer and the residual land, over faer's `solve_lstsq_in_place`
  (`Qᴴ B`, then `trtrs`). The reference factors a wide `A` as `L Q` (`xGELQF`); faer has no
  LQ, so `Aᴴ = Q R`, written back conjugate transposed, which is `xGELQF`'s layout. A zero on
  the factor's diagonal is `error.Singular` with its index in `pivot`, the reference's
  `INFO > 0`.
- **`syevd`, `heevd`, `syev`, `heev`** (`faer/evd.zig`): `dsyevd.f`'s and `dsyev.f`'s
  arguments, `xerbla` order and scaling into `[sqrt(SMLNUM), sqrt(BIGNUM)]`, over faer's
  `self_adjoint_evd`: tridiagonalize, take the real tridiagonal (a complex subdiagonal's
  phases moved into the vectors), decompose it by divide and conquer from 128 up (`syevd`, as
  the reference's `xSTEDC`) or by QR iteration (`syev`, as `xSTEQR`; and for values alone),
  apply the reflections. `NoConvergence` also for a NaN or infinity in `A`, which faer's
  solvers cannot take; with `jobz` `n`, `A` is left as it was.
- **`geev`** (`faer/evd.zig`): `dgeev.f`'s arguments, `xerbla` order, NaN check
  (`InvalidA`, the reference's `INFO = -4`), scaling and normalization (unit 2-norm, largest
  entry real; a real pair turned by `DLARTG`), over faer's `evd_real` and `evd_cplx`:
  Hessenberg reduction, Schur form by `lahqr` (`xLAHQR`'s double-shift and single-shift QR),
  eigenvectors of the triangular factor by back substitution (the left ones from its reversed
  conjugate transpose), times the Schur vectors. Faer has no balancing, so `dgebal.f`'s and
  `dgebak.f`'s (job `B`) come from Reference-LAPACK (`faer/balance.zig`); faer's solvers then
  run over the whole balanced matrix. Departures: back substitution guards small pivots (`eps ‖T‖₁`) but does not scale against
  overflow, as `xTREVC3` does; the eigenvalues are undefined after `NoConvergence`, where the
  reference keeps those it found. Fixes to faer ([`issues.md`](issues.md)): the complex
  `rotg` returns `r = 1` for a zero second entry, here `p`; the tridiagonal QR's and
  `self_adjoint_evd`'s reciprocals of denormals overflow, here divided by; the eigenvector
  solve guards a 2 x 2 block's determinant with a linear bound, here squared. Not ported yet: faer's multishift QR with
  aggressive early deflation (`multishift_qr`, from 75 up) and its blocked Hessenberg
  reduction (from 256 up); the unblocked forms serve every size, as `xLAHQR` and `xGEHD2`
  can.
- **Reference-LAPACK's own** (T115; 3.12.0 as OpenBLAS ships it in `lapack-netlib/SRC/`, the
  fallback where neither faer nor libflame has a method of its own): `lange` with Blue's `xLASSQ`
  and `gecon` with `xLACN2` (its reverse communication a callback) and `xLATRS` (`norm.zig`);
  `gebal` and `gebak`, every job (`faer/balance.zig`, which `geev` calls with job `B`); `trexc` and `trsen` (`schur.zig`: `xLAEXC`
  with `xLASY2`, `xLARFG`, `xLARFX`, `xLANV2`, `xLARTG`; `trsen`'s solves are ReLAPACK's `trsyl`);
  `gelsy` (`lsq.zig`: `xLAIC1`, and `xTZRZF`, `xORMRZ` by their unblocked `xLATRZ`, `xORMR3`) over
  faer's pivoted QR (`faer/qr_pivot.zig`, `geqp3`, which recomputes every norm where the
  reference's `xLAQPS` downdates them). Real and complex share one loop where the `d` and `z`
  sources differ only in a few steps, which the code marks. faer's `schur_swap` and `schur_move`
  are crate-internal translations of `dlaexc.f` and `dtrexc.f`, so `trexc` is the reference
  either way.
- **`gelqf`, `orglq`, `ormlq`** (`lq.zig`): libflame's `FLA_LQ_UT` method, the only LQ of its
  own: Householder blocks with the UT transform, which faer's `householder.zig` already
  implements, run on `Aᴴ` (LQ of `A` is QR of `Aᴴ`), with the reference's storage (`conj(v)` in
  each row, `tau` its own). No libflame code is ported. As faer's `geqrf`, a complex `L`'s
  diagonal stays complex (`tau` is real instead).
- **Not positive definite:** `potrf` returns
  `error.NotPositiveDefinite` and writes the failing minor's order (the reference's `INFO`) to
  an optional `minor: ?*usize`, its last argument, as Zig's diagnostics arguments do; `A` then
  holds the factor's first `k - 1` columns and the failed pivot at `A(k, k)`, as both upstreams
  leave it. The alternative, returning `INFO` as a value, lets a caller use an unfinished factor
  by forgetting a check; an error cannot carry the order itself.
- **Singular:** `getrf` (and `gesv`) completes the factorization, as both upstreams do, then
  returns `error.Singular` with the first zero `U(i, i)`'s `i` in an optional `pivot`. A NaN
  pivot is not zero, as in the reference's `INFO`: its column fills with NaN, and nothing
  fails. `trtri`, `trtrs` and `potri` return it for a zero on the diagonal, the same way.
- **`laed3`:** OpenBLAS's own (`lapack/laed3/`) serves only the reference's `xSTEDC`, which
  rstd does not port: `syevd`'s merge is faer's (`tridiag_evd.zig`), whose `z` from the roots
  takes `sqrt(|prod|)` where `xLAED3` takes `sqrt(-W)` (NaN for a spuriously positive `W`,
  which OpenBLAS's 71f22601 clamps to 0). `test/test_laed3_nan.c` is ported against `syevd`
  (`evd.zig`), at its `n = 8`, through divide and conquer at thresholds 4 and 8, and at
  `n = 200`. The reference routines `dsyevd` reaches that faer has no counterpart for:
  `xSTERF` (root-free QR for values alone; faer runs its QR iteration without vectors),
  `xLAMC3` (an addition forced through memory against extended registers) and `xLAMRG`
  (faer merges inline).
- **Sources.** [`sources.md`](../lapack/sources.md) lists, per routine, the candidates (OpenBLAS,
  ReLAPACK, faer, libflame, AOCL-libFLAME), which implement it, the choice and its evidence.

### Apple Accelerate

On macOS the public routines over `s`, `d`, `c` and `z` go to Accelerate where it has the
routine and matches the reference; `apple_accel.zig` binds it and routes to it.

- **Binding.** `load()` opens the framework with `dlopen` (`dll.openDarwinFramework`), never
  `-framework`, so cross-compiles need no SDK, and fills `Table(I)` with `dll.bind`. The
  declarations are ours, from Apple's documentation; no SDK file is copied. The handle stays
  open for the process, since the table points into it (`src/dll.zig`).
- **Symbols.** macOS 13.3 and later have `cblas_dgemm$NEWLAPACK$ILP64`, with 64-bit integers
  (rstd's `isize`), which `load` takes; else the 32-bit `cblas_dgemm$NEWLAPACK`, then the
  classic `cblas_dgemm` (macOS 13.3 Release Notes, Accelerate: "a new ILP64 interface",
  `ACCELERATE_NEW_LAPACK`, `ACCELERATE_LAPACK_ILP64`,
  <https://developer.apple.com/documentation/macos-release-notes/macos-13_3-release-notes>).
  The spellings are ours, checked by `apple_accel.zig`'s tests on the M6 ("the 32-bit table binds
  the new names first"). macOS 27 binds 150 of the `Table`'s 154 CBLAS fields: there is no
  `cblas_i?amin` (seen on the M6; no test checks it). `axpby` is
  `catlas_?axpby`. LAPACK's `?potrf_` is bound as `dpotrf$NEWLAPACK$ILP64` or
  `dpotrf_$NEWLAPACK$ILP64` (and the 32-bit forms, then the classic `dpotrf_`): which the new
  interface spells is unchecked until a Mac run.
- **Threads.** Threaded, Accelerate's `gemv` with `A` un-transposed rounds differently from one
  thread; `ddot`, `dnrm2`, `dasum`, `A^T` `gemv` and `gemm` gave the same bits at every thread
  count tried (macOS 27, M6). So `load` makes Accelerate single-threaded:
  `VECLIB_MAXIMUM_THREADS=1` (`setMaxThreadsEnv`), read at the process's first BLAS call
  (ACCELERATE(7), macOS 10.14: "may be cached by the library and reused"), and
  `BLASSetThreading` (macOS 15+), per thread, which `singleThreadHere` sets on each thread's
  first routed call (Apple, `BLASSetThreading`: "This setting is per thread",
  <https://developer.apple.com/documentation/accelerate/blassetthreading(_:)>). One thread cost `gemm` 1.15x at 2000 (our timing for 72b2efc; no results file kept).
  Rstd's workers thread it instead: above a routine's threshold a routed call splits into parts
  of its output fixed by the problem size, one single-threaded Accelerate call each
  ([`threading.md`](threading.md#vendor-backends) has which routines and why). The split's
  speed against one Accelerate call is untimed (T05).
- **Loading and the off switch.** `blas.init` loads it (`apple_accel.start`), since `setenv` is
  not thread-safe ("need not be thread-safe", POSIX.1-2024,
  <https://pubs.opengroup.org/onlinepubs/9799919799/functions/setenv.html>), unless `RSTD_BLAS_ACCELERATE` is 0 (`Config.accelerate`).
  `blas.useAccelerate` turns routing off and on later; nothing unloads it. `testUse` turns it off,
  so a test's own table runs its kernels.
- **Routing.** Each public routine makes rstd's checks and quick returns, then calls its function
  in `apple_accel.zig`, which calls the table (column-major, complex scalars by pointer,
  `cblas_?dot?_sub`) or returns false for rstd's kernels. On the LP64 table it does so for any
  size, stride or operand span past `c_int`, rather than truncate.
- **Routed:** `copy`, `swap`, `scal`, `dot`, `dotc`, `dsdot`, `sdsdot`, `nrm2`, `asum`, `rot`,
  `rotm`, `axpy`, `axpby`; `gemv`, `gbmv`, `ger`, `gerc`, the real `symv`, `sbmv`, `spmv`,
  `syr`, `spr`, `syr2`, `spr2`, the Hermitian `hemv`, `hbmv`, `hpmv`, `her`, `hpr`, `her2`,
  `hpr2`, and the six triangular routines; all of level 3, less what is kept for speed (below);
  LAPACK's `potrf` (its `INFO` becomes rstd's error and `minor`; after a failure `A` holds
  Accelerate's partial factor).
  Routed, `sdot` is `cblas_dsdot` rounded, and `snrm2` and `scnrm2` the root of `dsdot(x, x)`,
  so `s` sums in `f64` as rstd's kernels do. `d` and `z` `nrm2`, and `asum`, take `incx > 0` only.
- **Kept for speed (arm64).** Level 1 goes by routine, letter and stride, one table
  (`apple_accel.zig`'s `sends`): Accelerate gets every call, the unit-stride ones only (every
  stride 1), or none. A slot is kept where rstd's kernel ran faster than Accelerate past noise
  (ratio under 0.97) in the flyoff of 2026-10-29
  ([`results-flyoff-macos-aarch64.md`](../../../tools/blasbench/results-flyoff-macos-aarch64.md#accelerate-routing),
  routed against `--rstd-accelerate=0`, strided at 2); a slot within noise or faster stays
  routed. Rerun timings revise the table, not the code. Routed over rstd's own:
  - `nrm2`, none: `dznrm2` 0.06 to 0.21, `dnrm2` 0.25 to 0.48, `snrm2` and `scnrm2` 0.53 to 0.92.
  - `sdot`, none: 0.41 to 0.50 through `cblas_dsdot`, 0.97 strided. `dsdot` and `sdsdot`, with
    a kernel of their own the flyoff did not time, stay routed.
  - `scopy`, `dcopy`, unit only: 0.53 and 0.52 strided; `ccopy` and `zcopy` strided are 2.0 and
    1.8.
  - `swap`, unit only: 0.90, 0.90, 0.94, 0.96 strided.
  - Complex `dot` and `axpy`, unit only: 0.77 and 0.78, 0.82 and 0.81 strided; `dotc` follows
    `dotu`, the only one timed. Real strided `ddot` is 2.4, `axpy` 1.9 and 2.0.
  - `csrot`, unit only: 0.88 strided; `zdrot` strided is 0.97, noise.

  Left routed though slower: contiguous `zcopy` at 4M elements (0.83, where 64K is 4.2), a size
  limit no row places. x86_64 has no such data and keeps every call.
- **Not routed:** what CBLAS lacks: `sum`, the extremes, `axpyc`, `rotg`, `zrotg`, `rotmg`, the
  complex symmetric level 2, the conjugated-`x` letters `o` to `d`, a complex `r`, `HEMVREV`'s
  `v` and `m`; `i?amax`, which does not take the first NaN.
- **Quirks.** Where Accelerate's results differ from the reference's (macOS 27, M6), the call
  runs rstd's kernels. The Accelerate special-value tests (`apple_accel.zig`,
  `_test/blas2/edge.zig`, `blas3/interface.zig`) found each, against rstd's kernels or the wide
  reference.
  - Complex `scal`, `axpy`, `axpby` with a scalar that has a zero part: Accelerate skips that
    part's products, so `(1.5, 0) * (inf, 1)` keeps a finite imaginary part.
  - `axpby` with a zero `alpha`: it reads `x`.
  - The triangular level 2 with NaN or an infinity in `x`, and `trmm` with one in `B`: past
    small orders it multiplies them by the triangle's zeros.
  - `ger`, `gerc` with one in `x`: it does not pass over a zero `y`.
  - `syr`, `spr`, `her`, `hpr` with one in `x` and a zero `x(j)`, and `syr2`, `spr2`, `her2`,
    `hpr2` with one in `x` or `y` and a `j` where both are zero: it does not pass over that
    column, in either triangle.
  - `herk`, `her2k` with one in `C`'s triangle and a nonzero `beta`: past small orders it scales
    `C` by the real `beta` as a complex product.
  - The Hermitian rank updates leave the diagonal's imaginary parts as `beta` times theirs, or
    NaN where an infinity meets them; the reference sets them to zero. They are zeroed after the
    call; `herk` and `her2k` with `beta` 1 and nothing added, which leave `C` as it was, run
    rstd's.
  - x86_64 (under Rosetta 2; no Intel Mac was at hand): `cblas_dsdot` is off by `f32`'s epsilon,
    so `dsdot`, `sdot`, `snrm2` and `scnrm2` run rstd's; `gemm` with NaN or an infinity in `A`,
    which gives NaN at an edge of `m`.

### Core detection and selection

Each architecture family has a light detector, `cpu/_detect/<arch>.zig`, and a wrapper,
`cpu/<arch>.zig`, that adds the kernel registry and `vtableFromCoreType`. Only the target's
detector is analyzed, since detectors hold instructions (`cpuid`, `mrs midr_el1`, `cpucfg`)
other targets cannot assemble.

- **Core types.** `cpu/_coretypes.zig` holds every family's `CoreType` enum, import-free, in the
  order of the reference's `corename[]`. Static-only targets of the reference, which
  `dynamic*.c` never returns, are core types too: x86 `p5`, `p6`, `viac3`, `yonah`; arm64
  `cortexa76`, `cortexa510`, `cortexa710`, `cortexx1`, `cortexx2`, `ampere1`, `ft2000`; riscv
  `c910v`, `u74`, `x280`; mips `mips24k`, `mips1004k`, `p5600`, `i6400`, `i6500`, `p6600`
  (`select.staticOnly`). `OPENBLAS_CORETYPE` can name them; detection and fallbacks never pick
  them.
- **Detector contract.** A thin `probe()` fills a `Caps`; `decide(caps)`, `runnable(core,
  caps)` and the part tables are pure functions the tests feed synthetic values on any host.
  `decide` passes its answer through `select.guard`, so it never returns a core `runnable`
  rejects: the detector's `parent` chain first (arm64: Neoverse V2 to N2 to V1 to `armv8sve` to
  N1), then the closest runnable core. `supportsCore(core)` (`canRun` in `cpu/<arch>.zig`; mips
  has its own `canRun`) asks the running CPU, for tests.

| Family | Probe |
| --- | --- |
| x86, x86_64 | `cpuid` and `xgetbv`, unprivileged, so the same on every OS and freestanding |
| arm64 | Darwin: `sysctl` (`hw.optional.arm.FEAT_SME`, `FEAT_SME_F64F64`: XNU's `ARM_FEATURE_FLAG` list, <https://github.com/apple-oss-distributions/xnu/blob/ac9718fb1af618d5ce8678d0dc6e8a58f252216f/osfmk/arm/arm_features.inc#L118>, xnu-12377.121.6; OpenBLAS reads only `FEAT_SME`, `driver/others/dynamic_arm64.c:612`); Linux and FreeBSD: `MIDR_EL1` through the kernel's emulation (`HWCAP_CPUID`), else per-CPU sysfs `midr_el1` and `cpu_capacity` (Linux v6.12 `Documentation/ABI/testing/sysfs-devices-system-cpu`), plus the hwcaps (`HWCAP2_SME_F64F64` among them); Windows: `MIDR_EL1` from the registry, features from `IsProcessorFeaturePresent` (`PF_ARM_SME_F64F64_INSTRUCTIONS_AVAILABLE`, 85: <https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-isprocessorfeaturepresent>); freestanding and UEFI: the ID registers (`ID_AA64SMFR0_EL1.F64F64`) |
| arm | `AT_HWCAP` VFP and NEON bits and the part in `/proc/cpuinfo`; freestanding and UEFI: CP15 and VFP ID registers |
| power | `AT_PLATFORM`, then the PVR (`mfpvr`, emulated by Linux: v6.12 `arch/powerpc/kernel/traps.c`, `emulate_instruction`), then the ISA bits of `AT_HWCAP`/`AT_HWCAP2`; freestanding: the PVR |
| mips | Loongson `cpucfg` where `AT_HWCAP` advertises it, else `/proc/cpuinfo` (Linux only) |
| loongarch64 | `cpucfg` word 0 and `AT_HWCAP` LSX/LASX; freestanding and UEFI: `cpucfg` word 2 |
| riscv64 | `riscv_hwprobe`, else `AT_HWCAP`'s V bit and an RVV 1.0 check, then `vlenb`; freestanding: `misa`, `mstatus` and `vlenb` |
| s390x | `AT_HWCAP`, masked by `LD_HWCAP_MASK` (`probeMasked`) |
| sparc, wasm | nothing at run time: the target decides |

`cpu.choose` (`cpu/_select.zig`'s `choose`) is `gotoblas_dynamic_init` in two steps:

1. **Which core.** `OPENBLAS_CORETYPE`, matched case-insensitively, if it names a core
   (`unknown` names none) and the CPU can run it; else the detector's answer. A name no core has
   is reported at `OPENBLAS_VERBOSE` 1 and detection decides (the reference falls back to a
   default core: the generic one on arm64, mips64 and s390x, POWER8 on power, Katmai on x86,
   and by AVX level on x86_64, `dynamic.c:1155-1164`). A named core the CPU cannot run is reported and replaced by the closest one it
   can; on arm64, whose `runnable` cannot see the ISA level, only by one the detected core
   covers (`level_unseen`).
2. **Which object.** A core with no kernel object in the build (excluded by `-Dblas-cores`, or
   with no kernels at all) is served by the closest linked core that needs no feature it lacks
   (`select.served`), else by the build's pick for the target, else by the portable kernels.
   The served core's parameters are used too. The result is `Chosen{ core, served }`. On 32-bit
   x86 the AVX cores and Zen have no objects: as in `Makefile.system`'s 32-bit fallbacks, they
   take Nehalem's, the Bulldozer family and Zen Barcelona's (`corespec.objectCore`), whose
   kernels and blocking are theirs there.

"Closest" is decided by instruction set, not by a table: `select.Sets` derives each core's
feature set from its `Spec` (dependencies included, tuning-only features removed), and `closest`
prefers a core whose set is a subset of the wanted one's, with the most features, ties to the
earlier core. The reference's `#define gotoblas_X gotoblas_Y` chains follow no such rule
(Barcelona to Nehalem, `dynamic.c:93`, puts SSSE3 code, `getarch.c:313`'s `HAVE_SSSE3`, on a
CPU without it: `amdfam10` has none, Zig `lib/std/Target/x86.zig:1409`). Tests in `_select.zig` check that
every detectable core has an object that serves it and that no core is served by one that needs
more.

Where the port differs from `dynamic*.c`:

- **x86.** The OS AVX, AVX2 and AVX-512 state checks and the reference's fallbacks and warnings
  are kept. AMD family 10h and 12h are Barcelona (the reference returns NULL, `dynamic.c:924-925`); FMA3 is required
  for the Haswell and Zen kernels; `cpuid` leaves are range-checked; macOS's lazy ZMM state is
  not mistaken for no AVX-512.
- **arm.** The reference has no run-time dispatch; `cpuid_arm.c`'s build-time choice is made at
  run time from the hwcap bits, requiring 32 VFP registers for the ARMv7 kernels. Cortex-A9 and
  A15 are detected by part, as `cpuid_arm.c` does.
- **arm64.** The part table follows `dynamic_arm64.c`. Without the `MIDR_EL1` emulation the
  per-CPU sysfs pass reads `cpu_capacity` (decimal), as `cpuid_arm64.c:267` does;
  `dynamic_arm64.c:395` reads `cpu_capability`, which Linux does not have (Linux v6.12
  `Documentation/ABI/testing/sysfs-devices-system-cpu` lists `cpu_capacity` only). Every core built for
  ARMv8.3 and up needs FCMA, since its `trsm` kernels hold `FCMLA`; the Apple CPUs macOS and
  Linux run on (M1 and later) have it.
  The SME cores (`armv9sme`, `vortexm4`) are built with SME's optional F64F64, as the
  reference's are (`kernel/Makefile:30`, `:41`), so they need it detected: the reference picks
  them on SME alone (`dynamic_arm64.c:602-616`, `support_sme1`). Each probe
  decodes it from its raw values in a pure function, and it counts only alongside SME. A
  Windows without the SME flags answers false, so it never gets an SME core.
- **power.** POWER10 needs ISA 3.1 and MMA. POWER3 to 5, 970, Cell, G4 and 440, static targets
  upstream, are told apart by `AT_PLATFORM`, the PVR or the ISA bits, so none runs a later
  generation's code. A CPU nothing names gets the portable kernels (the reference guesses
  POWER8, `dynamic_power.c:344-346`). AIX is not a Zig target (`std.Target.Os.Tag` has no `aix`).
- **mips.** The reference tests `HWCAP_LOONGSON_CPUCFG` the wrong way round
  (`dynamic_mips64.c:168`) and its `cpucfg` word (`0xc8080118`, `:122`) has both registers zero; both are corrected, unverified on hardware.

### Configuration

`Config.zig` reads every variable the reference reads at init, once, into a value
(`Config.fromEnv`); `blas.config` keeps it. Defaults match the reference's unset behavior. Each
value is read one variable at a time with no allocation, from a map (`blas.init`), a POSIX
`environ` block, or the Windows PEB (`Config.Env`). A value longer than `value_max` (64 bytes, ours: room for any core name)
reads as not parsable and falls back to the default (for a core name, to detection).
Initialization also caches the CPU count in `cpu_count` (1 without threads).

| Variable | Field | Effect here |
| --- | --- | --- |
| `OPENBLAS_CORETYPE` | `core_type_buf` | forces a core (selection, above) |
| `OPENBLAS_VERBOSE` | `verbose` | 1: fallbacks and unknown names; 2: the selected core and the serving object |
| `GET_OPENBLAS_CORETYPE` | `get_openblas_coretype` | POWER only (`dynamic_power.c:352`): set to anything, prints the core to stderr |
| `OPENBLAS_BLOCK_FACTOR` | `block_factor` | scales a computed gemm `p` (clamped to 10 to 200 percent, `parameter.c:521-522`) |
| `OPENBLAS_L2_SIZE` | `l2_size_override` | replaces the detected L2 size (KiB) in the blocking rules |
| `OPENBLAS_NUM_THREADS`, else `OPENBLAS_DEFAULT_NUM_THREADS` | `num_threads` | LA464's blocking set; the thread ceiling of [`threading.md`](threading.md) |
| `GOTO_NUM_THREADS` | `goto_num_threads` | the fallback ceiling of [`threading.md`](threading.md) |
| `OMP_NUM_THREADS`, `OMP_ADAPTIVE`, `OPENBLAS_THREAD_TIMEOUT` | `omp_num_threads`, `omp_adaptive`, `thread_timeout` | parsed, no effect |
| `GOTO_DEBUG` | `goto_debug` | parsed, no effect |
| `OPENBLAS_MAIN_FREE`, `GOTOBLAS_MAIN_FREE` | `main_free` | parsed, no effect |
| `LD_HWCAP_MASK` | `ld_hwcap_mask` | s390x only: ANDed with `AT_HWCAP` before detection (`cpuid_zarch.h:35`) |
| `RSTD_BLAS_ACCELERATE` (rstd's own) | `accelerate` | macOS: 0 keeps the routines off Accelerate ([above](#apple-accelerate)) |

Integers parse with Zig's base-0 rules after trimming whitespace; negative values clamp to 0
except `GOTO_DEBUG`, `OPENBLAS_L2_SIZE` and `LD_HWCAP_MASK`. A bare leading zero is decimal, not
octal as with `strtoul` (`0o` is octal). Diagnostics print as `OpenBLAS: ...` through
`std.debug.print`, never on freestanding targets.

### Kernel provenance

Every per-core kernel is derived from OpenBLAS sources by `tools/blas_asm`, or ported by hand
from the C where noted. Generated files say which tool wrote them and must be regenerated, never
edited.

| Tool | Reads | Writes |
| --- | --- | --- |
| `gen.zig` | the `.S` kernels of level 1 and level 2 | `blas<n>/<routine>/<arch>/asm/*.zig` and their `__root__.zig` index |
| `uk.zig` | the inline-assembly microkernels of the C drivers (`<p><routine>_microk_<core>-2.c`) | `blas1/<routine>/<arch>/uk/*.zig` and their index |
| `params.zig` | `param.h`, through a C program per core | `blas3/gemm/params_data.zig` |
| `survey.zig` | the `KERNEL` files and `Makefile.L1`, `L2`, `L3`; rstd's `asm/` and `uk/` indexes, hand `ports` and registries | the `coverage.md` tables (every `blas1/<routine>/` and `blas3/<routine>/`, and `blas2/`) |
| `kernel_map.zig` | `kernel/<arch>/KERNEL.<core>`, `KERNEL`, `Makefile.L1` | (library) which source each core uses |
| `roles.zig` | the `kernel/Makefile.L*` rules | (library) each role's make variable and flags |
| `arch/<arch>.zig` | | (data) per-architecture defines, variants and patches |

Each tool runs as `OPENBLAS=<checkout> zig run tools/blas_asm/<tool>.zig [-- <names>...]`;
`BLAS_OUT` redirects the output.

#### Kernel selection

Which source a core uses follows `kernel/Makefile`: the core's `KERNEL.<core>` (and what it
includes), then the architecture's `KERNEL`, then `Makefile.L1` and its siblings. A plain
assignment in `KERNEL` overrides the core's choice; only an `ifndef` leaves it.
`kernel_map.zig` evaluates this, deciding `ifdef`/`ifndef` on the variable itself and listing
both sides of any other conditional as alternatives. Consequences of upstream's own selection,
reproduced as they are:

- x86_64 `csum`/`zsum` on the AVX-512 cores are `zsum_sse.S`/`zsum_sse2.S`
  (`kernel/x86_64/KERNEL:504-505`);
- loongarch64 `dsdot` is `dot.S` on LA264 and LA464 (`kernel/loongarch64/KERNEL:238`);
- arm64 `csum` is `csum.S` everywhere (`kernel/arm64/KERNEL:1`);
- MIPS64's MSA complex `axpy`, `rot`, `swap` and `scal` are never selected: `KERNEL` resets
  them to the `arm/` C (`kernel/mips64/KERNEL:1-11`).

#### Conversion of the assembly

`gen.zig` turns each `.S` into one naked function per kernel, whose body is the source's
assembly, unchanged in substance:

- **Preprocessing.** `zig cc -E` runs once per core and variant with the reference's defines
  (`-DDYNAMIC_ARCH`, the core's `getarch.c` macros, `DOUBLE`, `COMPLEX`, `CONJ`), since OpenBLAS
  selects code with `#ifdef <CORE>` and the ABI with `WINDOWS_ABI` and the like.
- **Variants.** One file per build variant (`Arch.variants`): x86_64 `sysv` and `win`, arm `hf`
  and `sf`, POWER `32`, `32le`, `64` (little-endian ELFv2) and `64be` (big-endian; Zig emits
  ELFv2 there too, but on PS3: Zig `src/target.zig:692` at 069480a1a), SPARC `32` and `64`. A variant leaves out the cores it never runs
  (`Variant.without`): POWER's 32-bit-only `ppc440`, `ppc440fp2` and `ppcg4` have no `64` or
  `64be` build.
- **Merged texts.** Cores whose expansions are identical share one file; texts that differ in
  under about a tenth of their lines (`similar`, our threshold: `tools/blas_asm/gen.zig:395`) share one too, with the differing lines under
  `parts`, chosen at comptime.
- **Rewrites.** Directives that declare the function or its section go (Zig places it); tabs
  become spaces; named labels become numeric local labels, so none reach the symbol table; x86
  `.align N`, bytes on ELF and COFF but `2^N` on Mach-O (LLVM 21.1.0: `MCAsmInfo.h`'s
  `AlignmentIsInBytes = true`, `MCAsmInfoDarwin.cpp` sets it false), becomes `.balign` (`align_bytes`);
  `uk.zig` does the same for the microkernels. A variant may rewrite every statement
  (`Variant.rewrite`; the sparc `64` variant's `sll` becomes `sllx`).
- **Output.** Each file holds `keys` (core and variant pairs), `serves(core)`, `kernel(core,
  variant)` and the naked `Body`. File names are the role, the source's suffix, the first core
  where one source splits, and the variant: `daxpy_atom_sysv.zig`,
  `daxpy_sse2_barcelona_win.zig`. The routine's `asm/__root__.zig` lists every file by slot;
  `asm_registry.Kernels` turns it into a module for `kernels/<arch>.zig`.
- **arm64's x18.** One build serves every OS, and x18 is the platform register on Windows (the
  TEB) and Apple (reserved; "Don't use this register"). Sources: Microsoft, "Overview of ARM64 ABI
  conventions", <https://learn.microsoft.com/en-us/cpp/build/arm64-windows-abi-conventions>;
  Apple, "Writing ARM64 code for Apple platforms" ("The platforms reserve register `x18`"),
  <https://developer.apple.com/documentation/xcode/writing-arm64-code-for-apple-platforms>.
  Upstream stopped using it as scratch in 0.3.18, for Apple M1 (`Changelog.txt:1298`); `gemv_n.S`, `gemv_t.S`, `zgemv_n.S` and `zgemv_t.S` still save and
  restore it with the callee-saved registers, writing back the value they found, which is
  harmless. `arch/arm64.zig`'s `check` stops the generator on any other use.

#### Fixes to upstream sources

The per-core tests found bugs in some of OpenBLAS's own assembly, still present on master at
v0.3.34-317. The generator fixes each before converting it: a `Patch` in
`tools/blas_asm/arch/<arch>.zig` is an exact find and replace on the upstream text that must
match `count` times, or the generator stops, with a comment saying what was wrong. The generated
files carry the fixed code. [`issues.md`](issues.md#openblas) holds the issue text for each. No
core takes the extremes' and `nrm2`'s kernels anymore ([`readme.md`](readme.md#one-contract-on-every-core)),
fixed or not, nor x86_64's peeling `zdot_atom.S`, so they are no longer converted (below). Their
patches stay in the same files as `retired`, the record of those fixes: every run checks that
they still match upstream, none applies them, and a full run stops on a live patch that applied
to nothing.

- `x86_64/zdot_atom.S`, Windows: Windows passes the address of a 16-byte result first, shifting
  the arguments along; the kernel read that address as `n`. Fixed the way `zdot_sse2.S` already
  does it.
- `x86_64/zgemv_t_atom.S`, Windows: read the double `alpha_i`, which Windows passes on the
  stack, with a single-precision load; its sibling kernels use the double one.
- `x86_64/scal_sse.S`, `scal_sse2.S`, `scal_atom.S`, Windows: read the flag (whether a zero
  `alpha` keeps NaN) from `incy`'s stack slot, `64(%rsp)`, not the tenth argument's,
  `80(%rsp)`, so a zero `alpha` stored zeros over NaN and infinity where the interface asks for
  NaN. `scal_sse.S` also took a NaN `alpha` for zero (`comiss` sets ZF for NaN; `scal_sse2.S`
  has the `jp`), so with flag 0, as for a level-2 `beta`, NaN became zeros.
- `x86_64/zscal_atom.S`: never read the flag, so a zero `alpha` always stored zeros; it now
  tests it as `scal_atom.S` does. A NaN imaginary part beside a zero real one also counted as
  zero.
- `x86_64/iamax_sse.S`, `iamax_sse2.S`, `izamax_sse.S`, `izamax_sse2.S` (every `i?amax`,
  `i?amin`, `i?max`, `i?min` on x86_64): the index search takes the last element as the match
  unread, but counted it even when it had been read, so for a NaN extreme, which no element
  equals, it returned `n + 1`. It now returns `n` there.
- `x86_64/copy_sse.S`, `copy_sse2.S`, `zcopy_sse.S`, `zcopy_sse2.S` (opteron, opteron_sse3):
  loaded by adding to a zeroed register, so -0 came out +0 and a signaling NaN quiet. They now
  load as on the other cores, so the copy is exact.
- `x86_64/symv_L_sse.S` (`ssymv_L` on atom, barcelona, bobcat, core2, dunnington, nano,
  opteron, penryn, piledriver, prescott): the last two rows of a block, and the last two
  columns, zeroed the upper half of the matrix operand but not of `x`'s, which holds the next
  elements, so an infinity there met a 0 and gave NaN in `y`.
- `power/scal.S`, `zscal.S` (power6, ppc970, cell): on Linux, read the flag argument at ELFv1's
  stack offset. Under ELFv2, and 32-bit, that slot lies past the caller's arguments, so a zero
  `alpha` could multiply where it should store zeros: `cgbmv` with `beta` zero gave NaN. The
  fix reads it through `FRAMESLOT`, which has each ABI's offset; FreeBSD's offsets were already
  right.
- `x86/amax.S`, `iamax.S`, `zamax.S`, `izamax.S`, the minima: took the new element unless it
  compared above the pick, and an unordered compare does not, so a NaN replaced the pick and the
  next element replaced the NaN. They now compare the pick with the new element, as the maxima
  do.
- `arm/amax_vfp.S`, `iamax_vfp.S`, the minima: moved on `lt`, which an unordered compare sets
  too, with the same effect; `mi` is less than only.
- `mips64/iamax.S`, `imax.S`, `izamax.S` and their minima (loongson3r3, loongson3r4), and
  `sparc/iamax.S`, `imax.S`, `izamax.S`: keep a pick per lane, and of equal picks the reduction
  kept the first lane's, which can be the later index. They now take the smaller index of equal
  picks.
- `arm64/axpy.S`, `zaxpy.S`, `daxpy_thunderx2t99.S` (every arm64 core but the SVE ones and
  `vortex`, `vortexm4`) and `loongarch64/axpy_lsx.S`, `axpy_lasx.S`, `caxpy_lsx.S`,
  `caxpy_lasx.S`: returned at once for a zero `alpha`, so NaN and infinities in `x` never
  reached `y`, unlike the other architectures' kernels; the level-2 triangular drivers pass
  `-x(j)`, which may be zero. The comparisons go. The ports of `arm64/axpy_sve.c`,
  `mips/?axpy_msa.c`, `riscv64/axpy_rvv.c`, `zaxpy_rvv.c` and `wasm/axpy.c` drop the same
  return.

The `lsx` and `lasx` kernels under `loongarch64/`, mostly where a stride is not 1:

- `sum_lsx.S`, `sum_lasx.S`, `csum_lsx.S`, `csum_lasx.S`: the double result was returned with a
  single-precision move, leaving its upper half undefined (NaN under qemu). `sum_lasx.S` also
  added doubles as floats.
- `casum_lasx.S`, `csum_lasx.S`: single precision's final sum added the last block of eight
  twice; `casum` also added lane 6 in place of lane 7.
- `axpby_lsx.S`, `caxpby_lsx.S`: with `beta` zero and a stride on `y`, the leftover elements
  were written over the first ones.
- `axpby_lasx.S`: double precision skipped every other block of eight, and with `beta` zero
  never stopped.
- `axpby_lsx.S`, `axpby_lasx.S`, `caxpby_lsx.S`, `caxpby_lasx.S`, with a zero scalar: the
  elements past the last block were `alpha*x + beta*y` whatever the scalars, from `x` (`y`)
  elements the zero-`alpha` (`beta`) blocks left their pointer at, so a NaN or infinity met the
  zero. The tail now reads neither operand of a zero scalar, as `kernel/arm/axpby.c` does. Double
  `axpby` also compared the scalars with a single-precision 0 (`ffint.s.l`, upper half
  undefined), so it never took the zero-scalar blocks; and complex, with one stride 1 and the
  other not, has no such blocks, so a zero scalar now goes to the strided ones, which do.
- `caxpby_lsx.S`, `caxpby_lasx.S`: with `alpha` and `beta` zero and a stride on `y`, zeroed
  only the real parts, and `lasx` single one of them twice.
- `caxpby_lsx.S`: double, with a stride on `x` and none on `y`, stored the upper half of each
  block's second pair over both halves.
- `max_lsx.S`, `min_lsx.S`: with a stride on `x`, compared by magnitude, as `amax` does.
- `camin_lsx.S`, `camin_lasx.S`: with a stride on `x`, the minimum started from `|re(x0)|`
  rather than `|re(x0)| + |im(x0)|`.
- `imax_lasx.S`, `imin_lasx.S`: with a stride on `x`, single precision below 8 compared against
  a first element never loaded, and double counted the indices up twice as fast.
- `iamin_lsx.S`, `iamin_lasx.S`: with a stride on `x`, double kept a minimum by the previous
  comparison, and `lasx` skipped the reduction across lanes.
- `icamin_lsx.S`, `icamin_lasx.S`: `lsx` below 4 elements started the index at 4 rather than 1,
  and with a stride on `x` its single precision started from an absolute value never computed;
  `lasx` single, with a stride on `x`, overwrote its index step, and its reduction over lanes 4
  to 7 mistook its own partial minimum for lane 4's value, returning lane 4's index.
- `snrm2_lsx.S`, `snrm2_lasx.S`, `cnrm2_lsx.S`, `cnrm2_lasx.S`: with a negative stride,
  returned 0. They scale by `?amax_k` of the same vector, which returns 0 for a stride below 1;
  they now walk `x` forwards from its lowest element.

Zero strides, found by the ported `utest/test_swap.c` and the per-core checks, which give every
two-operand kernel one zero stride beside a nonzero one. The reference is sequential, so a zero
stride means `n` operations on one element in turn:

- `sparc/swap.S`, `zswap.S`, `rot.S`, `zrot.S`: the strided loop loads a whole block before
  storing it, so `n` operations on one element acted as one. A zero stride now takes the
  one-at-a-time loop.
- `power/rot.S`, `zrot.S`, `swap.S`, `zswap.S`, `rot_ppc440.S`, `zrot_ppc440.S` (`utest`'s
  `zdrot_inc_0` on cell), loongarch64 `swap.S`, `rot_lsx.S`, `rot_lasx.S`, and mips64 `swap.S`:
  the same bug, and the same fix.
- x86 `swap.S`, `zswap.S` and arm `swap_vfp.S`, `rot_vfp.S`: upstream's own zero-stride cases
  skip the work or stop after one element, right only for some `n`; zero strides there now take
  one element at a time too.
- With only `incy` zero: x86_64 `zaxpy_sse2.S` (whose own zero-stride case fired only with both
  zero), `axpy_atom.S`, `zaxpy_atom.S`; power `axpy.S`, `zaxpy.S`, `axpy_ppc440.S`; sparc
  `axpy.S`, `zaxpy.S` (which now also returns for `n = 0`); loongarch64 `axpy_lsx.S`,
  `axpy_lasx.S`, `caxpy_lsx.S`, `caxpy_lasx.S`; and `cswap_lsx.S`, `cswap_lasx.S` with `incx`
  zero. A zero stride there takes a one-at-a-time loop.
- The C sources ported to Zig had it too, fixed in the port: `riscv64/axpy_rvv.c` started every
  lane at `y[0]` and summed them, adding `y[0]` once a lane, and `zaxpy_rvv.c`'s segments and
  `mips/saxpy_msa.c`, `daxpy_msa.c`'s blocks load `y` before storing it
  (`axpy/riscv64/rvv.zig`, `axpy/mips/msa.zig`).

The `@Vector` ports of C kernels fix two more ([`issues.md`](issues.md#openblas)):

- `x86_64/sasum_microk_haswell-2.c`, `dasum_microk_haswell-2.c`: peeled to a 32-byte boundary,
  so the sums followed the address; the port loads unaligned from the first element, as
  upstream's Skylake-X blocks do (`blas1/asum/x86_64/kernel.zig`).
- `arm64/symv_U_asimd_4x4.c`: with `m - offset` not a multiple of 4 dropped the rows between
  the last multiple of 4 and the first column; the port adds them (`blas2/symv/kernel.zig`).

Patches that work around Zig or QEMU rather than OpenBLAS, each removable once the cause is fixed
([`issues.md`](issues.md)):

- `power/zdot.S`, `zdot_cell.S`: they return the complex result in f1/f2 as the ELFv2 ABI says ("64-Bit ELF V2 ABI Specification", revision 1.5, 2.2.6, homogeneous floating-point aggregates; copy read: <https://kib.kiev.ua/x86docs/POWER/ABI/leabi-20201201.pdf>),
  but Zig reads a two-float struct from r3/r4. The sources already copy it into r3/r4 for the
  Fortran interface; the patch turns that on.
- arm `cdot_vfp.S`, `zdot_vfp.S`, `sf` variant: they return through a pointer in r0, as AAPCS
  says (AAPCS32 2024Q3, Result Return: a Composite Type larger than 4 bytes "is stored in
  memory at an address passed as an extra argument",
  <https://github.com/ARM-software/abi-aa/blob/2024Q3/aapcs32/aapcs32.rst>), but Zig reads the struct from r0 to r3, so the patch drops the pointer and returns there.
- sparc, `64` variant: the sources scale 64-bit strides with `sll`, which under qemu keeps only
  the low 32 bits, so negative strides broke; the variant writes `sllx`. Whether real hardware
  behaves as qemu does is open.

#### Platforms not converted

Alpha, IA-64, C-SKY, E2K, POWER BG/P and RISC-V C910V have no converted kernels; the reasons are
in [`readme.md`](readme.md#not-supported-by-zig). `Arch.unsupported` lists sources the toolchain
cannot assemble.

#### Kernels not converted, and `unused/`

`roles.zig`'s `skip` lists the kernels `gen.zig` and `uk.zig` leave out because no core would
install them: the routed slots (`vector.routed`), everywhere, and on x86_64 the peeling
reductions, whose cores the index keeps (`skipped`). Code that does not apply to rstd's
column-major Zig API is not kept as dead code either: the skipped kernels' converted files, the
vtable slots only row-major data or CBLAS reaches, and the stubs and tests of those slots live in
the git-ignored `unused/`, mirroring their old paths (snippets cut from a file that stays get
one file each, with a header naming the source). Nothing builds from `unused/`. Where an
interface mirrors an upstream one with a CBLAS row-major branch (`ger`'s `gerv_k`, the level-3
`*_direct*` kernels), that branch is a short comment, never code.

### Source tree

| Path | Role |
| --- | --- |
| `__root__.zig` | the module: the global vtable, `init`, `shutdown`, `vtable`, `table`, `kernel`, core queries, test hooks |
| `Vtable.zig` | `gotoblas_t`, derived from `common_param.h` |
| `Config.zig` | the environment variables |
| `types.zig` | element and index types, `XDouble`, `exprecision` |
| `portable.zig` | `Op`, slot names, type dispatch, `install`, `finish` |
| `vector.zig` | the `@Vector` level-1 and level-2 kernels, `routed`, `peels` |
| `asm_registry.zig` | `Kernels` (export and install of converted assembly), `install`, `served`, `isLinked` |
| `abi.zig` | the build variant for a target's calling convention |
| `fma.zig` | the contraction policy |
| `threading.zig` | thresholds, worker count, chunks, the claim loop and ordered combine |
| `util.zig` | negative offsets, unaligned vector loads, `kernelSymbol`, `fits`, `vectorIeee`, `blue`, `extendedSum` |
| `apple_accel.zig` | Accelerate's CBLAS, bound at run time, and the routing to it |
| `_kernel.zig` | root of each per-core object |
| `cpu/__root__.zig` | `choose`, `detect`, `fromName`, the target's `CoreType` |
| `cpu/<arch>.zig` | the family's `vtableFromCoreType`, registry and `canRun` |
| `cpu/_detect/<arch>.zig` | probes and pure decisions |
| `cpu/_coretypes.zig`, `_corespec.zig`, `_select.zig`, `_auxv.zig` | core enums, compile specs, selection, auxiliary vector and OS queries |
| `kernels/__root__.zig`, `kernels/<arch>.zig` | per-architecture kernel registries |
| `blas1/` | level-1 routines, one directory each; `interface.zig`; `extremum.zig` (portable extremes); `uk_block.zig` (C-driver ports); `c_kernels.md` |
| `blas2/` | `gemv`, `ger`, `symv`, `hemv` directories; `interface.zig`, `triangular.zig`, `banded.zig`, `symmetric.zig`; `buffers.md`, `coverage.md` |
| `blas3/` | one directory per routine; `interface.zig`, `driver.zig`; `gemm/` holds the tile kernels, `params.zig`, `params_data.zig`, `notes.md`, `coverage.md` |
| `ext/` | `bf16/` (bfloat16 and half kernels); `interface.zig`; `geadd/`, `matcopy.zig` (the `omatcopy` and `imatcopy` kernels), `laswp/` and `neg/` (their kernels over the gemm tile, in `kernel.zig`) |
| `../lapack/` | LAPACK over these routines: `cholesky.zig` (`potrf`, `potf2`, `potrs`, `potri`, the thread rules), `lu.zig` (`getrf`, `getrs`, `gesv`), `triangular.zig` (`trtri`, `trti2`, `lauum`, `lauu2`, `trtrs`), `relapack.zig` (ReLAPACK's `sytrf`, `hetrf` and `_rook` forms, `sygst`, `hegst`, `pbtrf`, `gbtrf`, `trsyl`, and under `-Drelapack` `potrf`, `getrf`, `trtri`, `lauum`), `indefinite.zig` (their Bunch-Kaufman kernels), `sylvester.zig` (`trsyl`'s); Reference-LAPACK's own: `norm.zig` (`lange`, `gecon`), `schur.zig` (`trexc`, `trsen`), `lsq.zig` (`gelsy`); `lq.zig` (`gelqf`, `orglq`, `ormlq`, libflame's method); `faer/`, faer's ports: `ops.zig`, `householder.zig`, `qr.zig` (`geqrf`, `orgqr`, `ormqr`), `qr_pivot.zig` (`geqp3`), `qr_solve.zig` (`gels`), `bidiag.zig`, `bidiag_svd.zig`, `svd.zig` (`gesvd`, `gesdd`), `tridiag.zig`, `tridiag_evd.zig`, `hessenberg.zig`, `real_schur.zig`, `complex_schur.zig`, `balance.zig` (`gebal`, `gebak`), `evd.zig` (`syevd`, `syev`, `geev`) |
| `blas3/gemmt/` | `gemmt`'s tests; it runs the driver |
| `_test/` | test-only code, wired by `_test/__root__.zig` (below) |
| `unused/` | git-ignored, kept locally: kernels, slots and stubs no column-major Zig caller reaches (above) |
| `readme.md`, `design.md`, `threading.md`, `comparisons.md`, `issues.md` | documentation |

A routine directory holds its `__root__.zig` (the vtable wrappers, its portable kernels and the
tests of every core's kernels), and per architecture `<arch>/asm/` (generated), `<arch>/uk/`
(generated) and `<arch>/kernel.zig` (composition and hand ports) as the routine needs.
Every `blas1/<routine>/` and `blas3/<routine>/`, and `blas2/`, has a generated `coverage.md`:
which OpenBLAS source each core uses, and what rstd serves in its place.

`_test/` holds the code only tests use:

| Path | Role |
| --- | --- |
| `_test/types.zig` | tests of `types.zig`, which the per-core objects compile without `zigstd` |
| `_test/blas1/check.zig` | `eachCore`, the wide reference `Wide`, `Pattern`, `rig`, `sqrtW`, `rootUlps` |
| `_test/blas1/edge.zig`, `_test/blas1/denormal.zig` | level-1 special values and denormals on every core |
| `_test/blas2/edge.zig`, `_test/blas2/guard.zig` | level-2 special values; the buffer guards |
| `_test/blas3/check.zig` | the level-3 harness each routine's tests call |
| `_test/from_openblas/` | OpenBLAS's own suites, ported |

Outside `src/linalg/blas/`:

| Path | Role |
| --- | --- |
| `build/blaskern.zig` | the per-core objects |
| `build/qemutest.zig` | the qemu rigs and the macOS VM |
| `build/sshtest.zig` | the SSH test and benchmark steps on a Mac |
| `build/wasmtest.zig` | the WASI test step |
| `build/__root__.zig` | `-Dblas-cores`, `rstd-test`, `rstd-check-macho`, `rstd-check-uefi`, `rstd-check-zig017` |
| `src/dll.zig` | run-time library loading, for Accelerate |
| `src/_test_qemu.zig` | test root of the OS-free rigs |
| `src/_test_wasm.zig` | the WASI test root: `_test_qemu.zig` |
| `test/qemu/` | boot code, linker scripts and runner of the freestanding and UEFI rigs |
| `test/linux/init.zig` | PID 1 of the Linux rigs |
| `tools/blas_asm/` | the converters |
| `tools/blasbench/` | the benchmarks: kernels, rstd against the libraries and reference binaries, the flyoff ([`readme.md`](../../../tools/blasbench/readme.md)) |
| `tools/test_host.zig`, `tools/nix/` | running the rigs' QEMU on another machine |
| `tools/ssh_test.zig`, `tools/qemu_watch.zig` | the SSH test runner; the rigs' watchdog |

**Cores handed to weaker CPUs** (T18). Where detection, as the reference's, can give a CPU a
core whose object it cannot run, each case is decided by the [settled
decision](#settled-decisions) that target platforms come first: fixed where a RandyOS target or a
zfls host could meet it, else the reference's choice kept and documented.

- **`katmai` on 32-bit x86 without SSE3, down to the 486:** kept. A 32-bit build's last resort
  is Katmai, whose kernels need SSE (`dynamic.c:1156`); 32-bit x86 is neither a RandyOS target
  nor a zfls host.
- **`vortex` on any Apple CPU:** kept for Darwin, which the reference also gives `vortex`
  (`dynamic_arm64.c:352-359`). The A7 to A11 lack FCMA, but run iOS only; every Apple CPU macOS,
  Linux and RandyOS run on is an M1 or later. Elsewhere `vortex` and `vortexm4` now need Apple's
  implementer code, since their objects are built for the M1 and M4.
- **`vortexm4` on a non-Apple CPU with SME and no SVE:** fixed. The reference picks it
  (`dynamic_arm64.c:543-547`), but such a CPU could be a Linux or Windows arm64 zfls host; it
  now takes the SVE test, then generic ARMv8.
- **A forced core's fallback to `neoversen1` on ARMv8.0:** fixed. The reference runs a named
  core unchecked (`force_coretype`); rstd replaces one the CPU cannot run, but arm64's `runnable`
  cannot tell ARMv8.0 from 8.2, so `OPENBLAS_CORETYPE=neoversev1` on a Raspberry Pi 3 or 4 fell
  back to Neoverse N1. A stand-in now must need nothing the detected core does not.

### Threading

[`threading.md`](threading.md) specifies threading, and its plan is implemented: every routine
takes an `io: Io`; levels 1 and 2 split work into chunks fixed by the problem size and combine
partial results in chunk order, so results do not depend on the worker count; level 3 follows
BLIS's `jc x ic` ways over a shared packed `B` on `Io.Group.concurrent` (`blas3/threaded.zig`);
all scratch is allocated before any task starts, and workers adopt the caller's floating-point
control state. `Config` parses the thread variables it reads. `threading.zig` holds the
machinery: the cached CPU count and ceiling, the thresholds, the chunks, the claim loop, the
ordered combine and the partial vectors, and level 3's team and barriers. `nrm2` and `gemmt` stay
serial, for the reasons `threading.md` gives.

## Settled decisions

Questions decided once, with the reason, so they are not reopened without new evidence.

- **No address dependence** (2026-10-23). Results never depend on where operands sit in
  memory. On the x86_64 cores whose converted reductions peel elements until an operand is
  aligned, every such slot takes the `@Vector` kernel, `dot` included: it loads unaligned, and
  its several partial sums bound the error more tightly than one running sum
  ([`readme.md`](readme.md), "No address dependence").
- **Column-major only** (2026-10-29). Rstd has no row-major entry points, and nothing upstream
  reaches only from CBLAS's row-major path survives: the Hermitian `Uplo` letters `v`/`m` go
  from `hpmv`, `her`, `hpr`, `her2` and `hpr2`.
- **`arm-softfp` is softfp** (2026-10-29). The soft-float ABI over the A15's VFPv4 and NEON,
  with the features pinned in `build/qemutest.zig`, so Zig 0.17 (whose model had them) and 0.18
  (whose model has none) build the same rig.
- **Timing on a host macOS keeps busy** (2026-10-29). Quiet the host first: the M6 auto-logs
  in, so `audiomxd` and `configd` rest, and with SIP off its indexing and model daemons are
  disabled (`tools/configure_mac_host.sh`). It runs macOS 27.2 (26B123). A run starts when no
  QEMU runs, the one-minute load is under 2, and no process passes 5% of a core in a
  `top -l 2 -s 5` sample; Activity Monitor and `sysmond`, which the user keeps open, are left out
  of the 5% test and their share is added to the limit. The limit is 2, not 1, because idle this
  macOS reads 1.1 to 2.3 with its CPUs 99% idle and no thread runnable (1.0 to 1.1 with one QEMU
  core busy), so under 1 never came in 8 hours. A watch every 20 s repeats any run another
  process disturbs; each run is made twice and a row kept only where the two agree within 3%; the
  load before and after is recorded, and a ratio within 3% of parity is reported as within noise;
  multithreaded rows are indicative only. Against the earlier busy-Mac runs (load 2.6 to 9.1) the
  one-thread library rates moved by a median of 0.0%.
- **Level-3 threading follows BLIS** (2026-10-29): `jc x ic` ways sharing packed `B` through
  barriers, written from BLIS's published design, not its code (F. G. Van Zee, R. A. van de
  Geijn, "BLIS: A Framework for Rapidly Instantiating BLAS Functionality", ACM TOMS 41(3), 1-33,
  2015, <https://doi.org/10.1145/2764454>; the threading in T. M. Smith, R. van de Geijn,
  M. Smelyanskiy, J. R. Hammond, F. G. Van Zee, "Anatomy of High-Performance Many-Threaded
  Matrix Multiplication", IPDPS 2014, 1049-1059, <https://doi.org/10.1109/IPDPS.2014.110>). Where
  OpenBLAS and BLIS both make a choice, rstd defers to it: workers spin, as both do (inside a
  call without bound: BLIS's `frame/thread/bli_thrcomm.c:264`, OpenBLAS's `driver/level3/level3_thread.c:408`;
  rstd's spin is bounded, then sleeps, [`threading.md`](threading.md)), and BLIS's two barriers guard one `B` buffer rather than the prototype's double buffer.
- **Filling upstream gaps** (2026-10-29). Where OpenBLAS and BLIS agree for a concrete reason
  (numerics, hardware, a measured trade-off), rstd defers to them. Where they leave a combination
  out only for simplicity, because their kernels cannot serve it, and a Zig `@Vector` kernel can,
  rstd fills the gap. So 32-bit POWER builds on POWER8 to 10 keep rstd's `@Vector` kernels: both
  upstreams run older code there only because their POWER8 to 10 kernels are 64-bit (OpenBLAS
  sends 32-bit builds elsewhere, `kernel/power/KERNEL.POWER8:1`: "Big-endian 32bit (AIX) is supported
  through the POWER6 GEMM kernels"; BLIS's POWER9 kernel is inline assembly loading pointers
  with `ld`, `kernels/power9/3/bli_gemm_power9_asm_d12x6.c:66`, its POWER10 ones C on MMA
  builtins, and it has no 32-bit POWER configuration).
- **Hermitian diagonals are real** (2026-10-29). `herk`, `her2k`, `her`, `hpr`, `her2` and
  `hpr2` set the diagonal's imaginary parts to zero, as the reference does: `A = Aᴴ` makes them
  zero, and the stored ones are not part of the matrix. Where Accelerate leaves `beta` times the
  old ones, rstd zeroes them after the call.
- **Optional CPU extensions only when detected** (2026-10-29). A core's object may use an
  extension only if detection (`runnable`) confirms it, so no core is chosen that can fault;
  `armv9sme` takes `sme_f64f64` on that condition.
- **Targeted testing** (2026-10-29). A change is tested where it reaches: filtered tests, and
  only the touched arches' rigs and compile checks (`qemu-check-<rig>`); the all-rig steps run
  periodically over batches, not per change.
- **Target platforms come first; the rest is best effort** (2026-10-29). The platforms that must
  work and be tested are RandyOS's target machines (its `readme.md`: x86_64 PC and Intel Mac,
  aarch64 Raspberry Pi 3 to 5 and Apple silicon, 32-bit arm Raspberry Pi 3, the 32-bit PowerPC
  750FX iBook G3) and the hosts that consume rstd (zfls: Windows, Linux and macOS on x86_64 and
  aarch64). Every other architecture is ported for completeness: where a choice cannot be tested
  there but OpenBLAS and BLIS agree on it, rstd follows them and documents why, rather than
  blocking on a rig; failures or gaps on those architectures are recorded, not release blockers.
  The tiers, platform by platform and core by core, are in [`support.md`](../../../support.md).
- **`rstd.linalg` is the public API; the BLAS is an implementation detail** (2026-10-29). The
  aim is a MATLAB-like library: callers pass vectors and matrices of any size and `linalg` picks
  the implementation; callers never choose a backend. `linalg.vector` takes fixed-size vectors
  and small matrices by value, with math names (`dot`, `length`, `cross`). Per size, it runs
  inline `@Vector` code or calls the BLAS; the crossover is measured, not guessed, and 2-, 3- and
  4-element forms likely stay inline. The choice depends only on size and target, so a call
  gives the same bits every time. The BLAS keeps its slice and stride routines for `linalg` to
  call. `rstd.simd` splits between them: slice reductions become BLAS routines, element-wise
  helpers (`mix`, `mulAdd`) move to `linalg.vector` (todo.md, M3). Under SPIR-V each operation is
  SPIR-V's own instruction where one exists (`OpDot`, GLSL.std.450 `Length`, `Fma`: GLSL.std.450
  1.00 revision 17, <https://registry.khronos.org/SPIR-V/specs/unified1/GLSL.std.450.html>), used as is,
  not rebuilt or tuned: the goal is one API, not GPU speed.
- **The matlab port defers to modern BLAS and LAPACK** (2026-10-29). `rstd.linalg.matlab`
  ports 1980s LINPACK and EISPACK algorithms. Where rstd's BLAS or LAPACK (chosen per routine,
  as the next decision says) has the operation, the modern routine replaces the old one; the old one
  stays only where a test shows it better (more accurate, or right on an edge case the modern
  one gets wrong), and that test is kept with it.
- **`rstd.linalg` is an anthology of the best open-source methods** (2026-10-30). Each routine
  or kernel takes the best open-source implementation of its method, chosen one at a time, as
  OpenBLAS's and BLIS's kernels are mixed today; nothing is all or nothing. For LAPACK the
  candidates are OpenBLAS's own LAPACK, ReLAPACK's recursive factorizations (MIT), faer (MIT,
  <https://codeberg.org/sarah-quinones/faer>), libflame and AOCL-libFLAME (BSD-3). Where only
  reference LAPACK had the routine, faer's port comes first, for speed out of the box. Reference
  LAPACK (modified BSD: permissive, no copyleft) is always available as the fallback where no
  other candidate has a routine or a base case (Q11, user, 2026-10-30), and a benchmark baseline
  otherwise. A choice changes when a benchmark or an accuracy test favors another source.
  Licenses must be MIT or BSD-style, permissive and free for commercial reuse; anything with
  copyleft, or any other license, is asked about first.
  Eigen (MPL 2.0, which would keep any port under MPL) is not a source: it appears only in
  benchmarks, built from upstream as a reference binary, never ported (user, 2026-10-30).
  The candidates and the choice for each routine and kernel family are in
  [`../lapack/sources.md`](../lapack/sources.md) (T113).
- **Bit guarantees are the CPU's** (2026-10-29). The BLAS's promises (the same bits for any
  thread count, one contract on every core, the first-NaN rule) hold on the CPU. Under SPIR-V
  the driver sets its instructions' precision and may fuse `a*b + c`, so Vulkan's precision
  rules apply there.
- **Shaders keep their own API** (2026-10-29). `src/rgl/shaders/lib/` (shaderlib) is the one API
  shader code uses: the standard library for SPIR-V, with GLSL-friendly names that alias or wrap
  `linalg` calls where `linalg` has the operation (`dot`, `length`, the matrix-vector products).
  `cross` moves to `rstd.linalg.vector`, which shaderlib re-exports; `normalize`, `reflect`,
  `refract`, `faceforward`, `smoothstep`, `fract` and `step` stay in shaderlib, as shader
  operations lowered to SPIR-V's own instructions. Shaderlib is the developer-friendly face, not
  a second implementation: everything a shader calls ends in a SPIR-V instruction where one
  exists, through `linalg` or directly.
- **Workers adopt the caller's floating-point control state** (2026-10-29). Threads do not
  inherit it (x86 Windows starts them at 53-bit x87 precision where a MinGW main thread has 64-bit, sources in [`threading.md`](threading.md);
  MXCSR and FPCR differ too if the caller changed rounding or flush to zero), so each task
  starts by applying the state captured on the caller's thread. Rstd still never changes the
  caller's own environment; its workers copy it, and threaded results match serial ones for any
  caller ([`threading.md`](threading.md)).
- **Rank updates keep the reference's zero skip** (2026-10-29). `ger` passes over a column whose
  `y(j)` is zero, `syr`, `her` and their packed forms one whose `x(j)` is, the rank-2 forms one
  where both are, as Reference-LAPACK 3.12.1 and master do, so a NaN or infinity in the other
  vector never meets that zero. OpenBLAS and BLIS disagree with each other and give no reason, so
  the reference decides; revisit when it drops the skip, as its PR #1323 did for the triangular
  routines (<https://github.com/Reference-LAPACK/lapack/pull/1323>, merged 2026-07-13) ([`comparisons.md`](comparisons.md) (5)).
- **Level-2 products fold fixed partial vectors** (2026-10-29). `symv`, `hemv`, `sbmv`, `hbmv`,
  `spmv`, `hpmv`, `trmv`, `tbmv`, `tpmv` and untransposed `gbmv` sum each element of `y` over
  several columns' kernel calls. Above their thresholds the columns split into chunks fixed by
  the size, never the worker count; each chunk adds into a zeroed partial vector of its own and
  the partials are added to `y` in chunk order, as OpenBLAS's `symv_thread.c` sums per-thread
  buffers with `AXPYU_K`. The bit rule then holds by construction on every core, at the cost of
  changing these routines' results from the one-pass serial order for sizes past one chunk,
  which a one-worker call shares; below one chunk they run as before ([`threading.md`](threading.md)).
- **LAPACK failures are errors** (2026-10-29). A matrix that is not positive definite is
  `error.NotPositiveDefinite`, its failing minor written to an optional `minor` argument, so a
  caller cannot go on with half a factor by forgetting a check; a NaN pivot fails, as in the
  reference's `xPOTRF2`, where OpenBLAS's `ajj <= 0` passes it ([Extensions](#extensions)).
  An exactly singular `U` is `error.Singular`, its first zero pivot written to `pivot`; a NaN
  pivot is not zero, as the reference's `xGETRF2` has it.
- **Loongson 3A fuses, as the reference does** (2026-10-29). `loongson3r3` and `loongson3r4`
  are modeled fused: OpenBLAS builds them with `-march=loongson3a` (`Makefile.system:1002`), where
  GCC contracts to `madd.fmt`, which GS464 and GS464E round once (their manuals, `issues.md`).
  QEMU rounds `madd.fmt`, so no rig can check the fused bits (the axpy fused-bit test skips on
  rigs); correctness rests on the manuals and on reading the objects' assembly.
- **Q9: futures on `Dispatch` as a fallback until std's `Group` is fixed, then revert**
  (2026-10-30). On `Io.Evented` with worker threads std's `Group` crashes on `Uring` and
  `Dispatch` ([`issues.md`](issues.md)). On `Dispatch` (macOS) rstd's worker tasks are plain
  futures, chosen by the comptime `threading.dispatch_futures`; every other `Io` keeps `Group`.
  `Uring`'s futures have their own bug, so it stays unsupported with worker threads. Once
  `ssh-repro-group-aarch64-macos` passes on a new toolchain, the switch goes back to `Group` and
  the futures path is deleted ([`threading.md`](threading.md)).

## Alternatives considered

| Alternative | Reason not chosen |
| --- | --- |
| A drop-in `libopenblas` (CBLAS and Fortran symbols, `xerbla`) | A Zig API with named errors is the goal; an ABI layer would add the row-major mapping and callbacks for no Zig caller |
| One function per precision (`dgemm`, `zgemm`) | Duplicates every interface; a comptime type serves every float type, including those without slots |
| A vtable with fields compiled in or out per build | Every field present and `null` where absent keeps one layout for the main library and every object |
| Selecting kernels lazily per routine | Eager installation keeps a kernel call to one indirect call and the table immutable while kernels run |
| A mutex around initialization | Brings threading support into freestanding builds; a three-state handshake suffices |
| The reference's `gotoblas_X` fallback chains | They follow no rule and put SSSE3 code on CPUs without it; the ISA subset rule cannot |
| Fall back to the generic core for an unknown `OPENBLAS_CORETYPE` | Detection gives a faster correct core (reasoned, not measured: generic runs the portable kernels; OpenBLAS's x86_64 picks by ISA there, `driver/others/dynamic.c:1157`); the name is reported |
| Each core's converted kernel in every slot, with quirks listed | Results then mean different things per core (`sdot` precision, NaN in the extremes, `nrm2`); [`comparisons.md`](comparisons.md) (X1) |
| Reproducing the reference's accumulation order in the portable kernels | Not required for correctness; kept only where a Zig kernel replaces converted code of the same core |
| Converting the gemm assembly | The Zig tile runs at the measured peak on AVX2 and within a few percent on Apple cores (`notes.md`, `results-m6.md`); one source serves every core |
| Setting FTZ or Java mode to make vector lanes IEEE | Changes the caller's environment, per thread; scalar code is used where lanes flush |
| `std.math.Complex` returns from kernels | Automatic layout cannot be a C return type; `CComplex` is an `extern struct` |
| Fortran `INTEGER` sizes | No Fortran ABI; `isize` addresses any array that fits in memory |
| Hand-editing converted kernels | Edits would be lost on regeneration; patches on the upstream text fail loudly when they stop matching |

## Verification

### Test layers

| Layer | Files | What it checks |
| --- | --- | --- |
| Per-core kernels | each routine's `__root__.zig` | every linked core's kernels (`check.eachCore`) against a plain loop or a wide reference, over sizes and strides past every unroll, and zero strides, one beside a nonzero one too |
| Edge cases | `_test/blas1/edge.zig`, `_test/blas2/edge.zig` | special values against an `f128` reference, NaN in gaps and past the end, negative and page-sized strides; the reductions give the same bits with their operands shifted by any number of floats (`gemv`'s in `blas2/gemv/__root__.zig`) |
| Denormals | `_test/blas1/denormal.zig` | denormal operands and results through every level-1 kernel, unflushed; `f32` lanes flush where `util.vectorIeee` says they do |
| Level 3 | `_test/blas3/check.zig` (used by each routine's tests) | every side, triangle, transpose and diagonal, shapes around block edges, special values |
| Bit-exact models | `blas1/dot/__root__.zig`, `blas1/asum/__root__.zig`, `blas1/nrm2/__root__.zig`, `blas2/gemv/__root__.zig` | the Zig kernels that replace converted or C code give its bits (the x87 sums, `skylakexModel`, the half-precision `gemv`) |
| Buffers | `_test/blas2/guard.zig` | each level-2 kernel gets exactly the buffer its contract promises, with canaries |
| Interfaces | `blas1/interface.zig`, `blas2/interface.zig`, `blas3/interface.zig` | argument errors and their order, quick returns |
| LAPACK | `../lapack/cholesky.zig` | on every core, shrunk and own blocking: `potrf` of `A = L Lᴴ` from a known `L` (backward within `2 (n + 2)` ulps, forward loosely), the exact failing minor for a negative, NaN or infinite element, nothing outside the triangle read or written; `potrs`'s residual within `8 (n + 2)` ulps of `L`, `Lᴴ` and `X` multiplied in magnitude, a NaN in `B` filling its column alone; `getrf` on square, tall and wide shapes (`P A` within `2 (min(m, n) + 2)` ulps of `L` and `U` multiplied in magnitude, `L` at most 1 in magnitude, `sqrt 2` complex), the first zero pivot with the factorization completed, NaN never singular, pivoting ties, NaN and infinite pivots, denormal pivots divided by; `getrs` for every `op`, its residual within `8 (n + 2)` ulps of `L`, `U` and `X` multiplied in magnitude, a NaN in `B` filling its column alone; `gesv` bit for bit `getrf` then `getrs`, singular before `B` is touched; `potri`'s `A X = I` within `32 (n + 2)` ulps of `A` and `X` multiplied in magnitude; `trtri`'s `A X = I` within `8 (n + 2)` ulps, both triangles and diagonals, a unit diagonal never read; `lauum`'s products within `2 (n + 2)` ulps, the complex diagonal exactly real; `trtrs` for every `op`, within `8 (n + 2)` ulps; each blocked routine's single and parallel algorithms, and the threaded ones (`potrf`, `potri`, `trtri`, `lauum`, `trtrs`, `getrf`'s lookahead, `getrs`, `gesv`) bit for bit one worker's at 2, 3 and 8, inline and on `Io.Threaded`; `potrf` through Accelerate on macOS; ReLAPACK's routines (`relapack.zig`): `sytrf`, `hetrf` and the `_rook` forms, random and zero-diagonal, rebuilt as `L D Lᵀ` from LAPACK's format within `32 (n + 2)` ulps of `A` and the factors multiplied in magnitude, singular at a zero `D` block; `sygst`, `hegst` for every `itype` against `M X Mᴴ` within `32 (n + 2)` ulps; `pbtrf` within `2 (n + 2)` ulps of the factor and its transpose multiplied in magnitude, the minor for a negative or NaN pivot; `gbtrf`'s `P L U` within `2 (min(m, n) + 2)` ulps; `trsyl`'s residual within `32 (m + n + 2)` ulps for every `op` and sign on quasi-triangular input; ReLAPACK's `potrf`, `getrf`, `trtri`, `lauum` by `cholesky.zig`'s and the rest's checks, and with `-Drelapack` against OpenBLAS's within `64 n` ulps; the bits inline and on `Io.Threaded` |
| LAPACK, faer | `../lapack/faer/*.zig` | on the selected core: `gesvd` and `gesdd` in the four precisions, every job, square, tall and wide (past `11 / 6`, and past 128 for divide and conquer), `A - U Σ Vᴴ`, `Uᴴ U - I` and `V Vᴴ - I` within faer's `test_svd` tolerance (`128 eps sqrt(8 max(m, n))`), `Σ` nonincreasing and the same for every job; faer's `test_special` (zero, ones, identity), rank deficient, scaled past `SMLNUM` and `BIGNUM`, a lone denormal; NaN and infinity; the argument errors; one worker's bits at 3 and 8; and the parts: `Q R = A`, `Uᴴ A V = B`, block reflections against single ones, Givens, faer's bidiagonal tests (`test_qr_algorithm`, `test_divide_and_conquer` at thresholds 4 and 8, `test_zink`, the deflation cases). QR: `Q R = A` and `Qᴴ Q - I` within `4 (m + n) eps` of `‖A‖`, square, tall, wide, rank deficient and a zero column; `ormqr` for both sides and `op`s and `orgqr` against one reflection at a time, with geqrf's `tau` and the reference's complex form; `gels` against a known least-squares or minimum-norm answer for `A` and `Aᴴ`, tall and wide, the residual in `B`'s trailing rows, an inconsistent system's `Aᴴ r` (`dqrt17.f`, Reference-LAPACK's threshold 30), rank deficient, zero, scaled, NaN. Eigen: `syevd` and `syev` for both triangles at faer's sizes and past 128, `A V - V Λ` and `Vᴴ V - I` within faer's `128 eps sqrt(8 n)`, values alone the same; faer's `test_special`, clustered, graded and Wilkinson matrices; OpenBLAS's `test_laed3_nan.c`; faer's tridiagonal tests (`test_evd_2_0` to `test_gh_82_mini`) through QR iteration and divide and conquer, and random, clustered, graded and zero tridiagonals; `Qᴴ A Q` tridiagonal and Hessenberg; `lahqr`'s `Z T Zᴴ = H` with standardized blocks; `geev`'s right and left residuals, unit vectors with a real largest entry, pairs positive first, on random, zero, identity, Jordan, rotation, companion and triangular matrices; balancing (`D⁻¹ Pᵀ A P D`, and `D M D⁻¹`'s eigenvalues to `M`'s scale); NaN, infinity, scaling and the argument errors; and one worker's bits at 3 and 8 for `geqrf` with `ormqr`, `syevd` and `geev` |
| LAPACK, reference | `../lapack/norm.zig`, `faer/balance.zig`, `schur.zig`, `lq.zig`, `lsq.zig`, `faer/qr_pivot.zig` | on the selected core, the reference's own test ratios and thresholds (`TESTING/dtest.in`'s 30, `dec.in`'s 20): `gecon` within `dget06.f`'s ratio of `1 / (‖A‖ ‖A⁻¹‖)`, `latrs`'s scaled residual near overflow, `lange`'s norms at the range's ends; `gebal` on `dbal.in`'s cases and `gebak` back; `trexc` every move with `dget36.f`'s residuals, `trsen`'s `s` and `sep` against the Kronecker operator; `gelqf`, `orglq`, `ormlq` by `dlqt01.f`, `dlqt03.f`; `geqp3` by `dqpt01.f`, `dqrt11.f`; `gelsy`'s rank and the SVD's minimum-norm solution; `gelqf` and `gelsy` give the same bits at 1, 3 and 8 workers |
| Ported suites | `_test/from_openblas/` | OpenBLAS's `blat1`, `blat2`, `blat3`, `utest`, `utest/test_extensions` |
| Selection | `cpu/`, `cpu/_select.zig`, `cpu/_detect/*.zig` | decisions on synthetic `Caps`, served objects, names |
| Initialization | `__root__.zig` | the handshake, `shutdown`, environment copies, racing threads |
| Small files | `_test/types.zig`, `fma.zig`, `portable.zig`, `ext/`, `Config.zig`, `abi.zig`, `apple_accel.zig` | types, contraction, dispatch, bfloat16, parsing, variants, the Accelerate binding and routing |

Test names follow the path (`linalg.blas._test.blas1.edge.test...`), which the rigs' filters
use.

### Every core

`_test/blas1/check.zig`'s `eachCore(check)` runs `check` on a vtable holding each linked core's
kernels over the portable ones (`portable.install`, `cpu.kernel.install`, `portable.finish`), for
every core this machine can run (`canRun`). A failing core is named and the rest still run.
Natively, only the target's core is linked in debug unless `-Dblas-cores=all`; the rigs link
every core the target CPU runs. `vector.installAll` puts every `@Vector` kernel in one table, so
they run natively too.

### References and tolerances

- **Wide reference.** Level 1 and 2 compare with an `f128` reference (`check.Wide`); level 3
  builds each routine's dense operator at twice the precision (`f64` for `f32`, else `f128`) and
  accepts `2 (k + 2)` ulps of the sum of magnitudes (this and the LAPACK tests' `2 (n + 2)` and
  `8 (n + 2)` are our margins over the roundings per element, not published bounds), `trsm` its substitution step by step.
  `sqrtW` adds Newton steps where sparc64's `_Qp_sqrt` works in `f64` (Zig
  `lib/compiler_rt/sqrt.zig:493-494` at 069480a1a); `rootUlps` allows for the library's own
  `f128` roots there.
- **Order-free comparisons.** A core may sum in any order, so plain-loop checks allow rounding
  error; complex infinities may meet in another order than the reference's
  (`expectWideAnyOrder`).
- **No flush tolerance.** Denormals must come out as the reference has them on every core.
- **Exactness.** Bits are compared only where a Zig kernel models converted code of the same
  core: the x87 sums, the x86_64 `dot` vector blocks (`skylakexModel`), the half-precision
  `gemv`. The ported C drivers were also compared bit-for-bit against OpenBLAS's C, compiled with
  a stub `common.h` ([`blas1/c_kernels.md`](blas1/c_kernels.md)).

### Poisoning and guards

- **What must not be read holds NaN:** the other triangle, a unit diagonal, a Hermitian
  diagonal's imaginary parts, band corners, padding rows, the gaps of strided vectors, `y` for a
  zero `beta`, the operands for a zero `alpha`.
- **What must not be written holds a sentinel or NaN** that must stay: rows of `C` past `m`,
  what follows its last column, the triangle a rank-k update leaves, the gaps and the elements
  past a vector.
- **Scratch is poisoned.** The test allocator hands out NaN-filled memory (`guard.zig`), and the
  interfaces' stack buffers hold NaN in tests, so a kernel that reads scratch it did not write
  fails.
- **Blocking is shrunk.** Level 3 runs on a vtable whose `P`, `Q` and `R` are shrunk to `2 *
  unroll_m`, 5 and `6 * unroll_n + 1` (ours, `_test/blas3/check.zig:624`), so every block and panel edge comes at small sizes, then
  once at the core's own blocking, past `R` (natively also past twice `P` and `Q`).
- **Foreign tiles.** `blas3/gemm/__root__.zig` installs tiles other targets choose
  (`installTile`: NEON's, 8 x 3 complex) to run every routine over them on any host.

### One contract, no quirks

Every core keeps the portable kernels' contract, so the edge and denormal tests hold every core
to the same expectations. Where a core's converted kernels did otherwise, they are not used
(`vector.routed`) or are patched (above); [`issues.md`](issues.md) records what upstream's do.

### OpenBLAS's suites

`_test/from_openblas/` ports the reference BLAS programs `test/[sdcz]blat[123].f` with their
`.dat` sweeps, the parts of their CBLAS versions (`ctest/`) that differ, and the unit tests
(`utest/`, `utest/test_extensions/`), transcribed mechanically from v0.3.34-317.

- **Tables.** The `.dat` tables run natively; the rigs (`reduced()`) take the example table in
  each program's header (`N` up to 9). One more test per level runs the suites, on that table,
  on every linked core the machine runs.
- **Translation.** Bad arguments map onto the named errors. What cannot be translated is listed
  in `_test/from_openblas/__root__.zig`: `xerbla` and bad letters (only the letters a type
  rejects remain), checks that scalars come back unchanged, row-major CBLAS, `test_fork.c`,
  `test_post_fork*.c`, and the `ARMV8` `sdot`/`ddot` cases that preload a vector register.
- **Pending.** `blat3` runs every level-3 routine through the flat API, and `test_potrs.c`
  and `test_kernel_regress.c`'s `dgesvd` run through `rstd.linalg.lapack`, the latter's
  matrix from glibc's `rand` as on OpenBLAS's Linux CI. `compare.zig` ports `test/compare_sgemv_*.c` and
  `compare_sgemm_*.c`: the half routines against `sgemv` and `sgemm` and a plain loop, to order
  256 natively, 24 on the rigs and per core.

### Rigs

`src/_test_qemu.zig` is the root of the OS-free tests: `simd`, the whole BLAS tree and LAPACK;
`src/_test_wasm.zig` runs it on WASI. The rigs run it, or the full test suite, where no host
here can.

| Rig | How | Covers |
| --- | --- | --- |
| `qemu-test-<machine>` (freestanding) | a test kernel per architecture under `qemu-system-*`: boot code and runner in `test/qemu/`, `single_threaded` but on the `-smp` machines, a 32 MiB heap (ours, `test/qemu/heap.zig:14`), every core the target CPU runs linked, LLD for the linker script | every BLAS and LAPACK test but `abi.zig`'s and `apple_accel.zig`'s, which need an OS (`default_filters`) |
| `qemu-test-<arch>-uefi` | the freestanding rigs' tests and runner as a UEFI application (`test/qemu/uefi.zig`) under QEMU's OVMF and AAVMF, on rstd's UEFI `std.Io`; x86_64 and aarch64 | RandyOS's loader's platforms: PE/COFF, the Microsoft ABI on x86_64, UEFI's CPU state |
| `qemu-test-x86_64-whpx` | the x86_64 kernel on the host CPU under WHPX, `SapphireRapids,+fma4,+xop`; x86_64 Windows hosts only | the AVX-512 cores and the Bulldozer family, which TCG cannot run |
| `qemu-test-linux-<arch>` | a Debian netboot kernel with the hosted test executable as PID 1 (`test/linux/init.zig`), `-smp 4` | auxiliary-vector detection, threads, files; aarch64, x86_64, arm, riscv64, powerpc64le |
| `qemu-test-x86_64-macos` | an installed x86_64 macOS VM under TCG, tests over SSH | Mach-O at run time on an x86 CPU model |
| `ssh-test-aarch64-macos` | the full suite on an Apple Silicon Mac over SSH (`-Dssh-mac`, default `macmini-m6`) | `vortex`, `vortexm4`, SME without SVE, which QEMU's `-cpu max` cannot show; Accelerate |
| `ssh-test-x86_64-macos` | the x86_64 macOS suite on that Mac under Rosetta 2 | Mach-O x86_64, faster than the VM |
| `wasm-test` | WASI under Node, a module with SIMD128 and one without (`wasm-test-simd128`, `wasm-test-portable`) | the freestanding rigs' tests, `abi.zig` and `threading.zig`, on the SIMD core and the portable path |
| `rstd-check-macho` | build only, x86_64 and aarch64 macOS, every core | Mach-O assembly errors |
| `rstd-check-uefi` | build only, x86_64 and aarch64 UEFI, the UEFI rigs' runner, no filter | every OS-free test compiles for UEFI, also those the rigs leave out |
| `rstd-check-zig017` | build and link only, with Zig 0.17.0, zfls's six targets, every core | zfls's toolchain builds rstd (Requirements) |
| `rstd-test-x87` | the native suite as `pentium3` code (`x86-windows-gnu` or `x86-linux-musl`), Debug and ReleaseSafe; x86 hosts | x87 widths with threads and an OS; run in each full sweep, with the `-safe` rigs |

The freestanding and UEFI machines (`build/qemutest.zig`):

| Machine | CPU | Notes |
| --- | --- | --- |
| `aarch64` | `cortex-a72` | |
| `aarch64-sve` | `max` | SVE2, SME2, BF16: every arm64 core but Apple's |
| `aarch64-smp`, `x86_64-smp` | as `aarch64`, `x86_64`; four cores | multithreaded kernels: the threaded tests on their cores (`test/qemu/smp.zig`, [`threading.md`](threading.md) Verification) |
| `arm`, `arm-softfp` | `cortex-a15` | `eabihf`, and `eabi` for the kernels' `sf` variant on the VFP the CPU has |
| `riscv64` | `rv64,v=true,vlen=512` | the VLEN of the widest linked core (x280; `zvl512b` in `Makefile.riscv64:6`) |
| `riscv64-vlen128`, `riscv64-vlen256` | `rv64,v=true,vlen=128`, `256` | `zvl128b`'s and `zvl256b`'s objects on the VLEN they are named for, the VLEN-agnostic `axpy` at both |
| `loongarch64` | `la464` | |
| `mips64r2el` | `Loongson-3A4000` | MSA, which the boot turns on: the `loongson3r4` core's kernels (QEMU v10.1.0 `target/mips/cpu-defs.c.inc` sets `CP0C3_MSAP` for it, `I6400` and `P5600`, not `Loongson-3A1000` or `24Kf`) |
| `mips64r2` | `Loongson-3A1000` | big-endian, no MSA |
| `mips64r6el`, `mips64r6` | `I6400` | the i6400, i6500 and p6600 cores' MSA kernels |
| `mips32r2el`, `mips32r2` | `24Kf` | no 32-bit kernels apply: the portable ones run |
| `mips32r5el`, `mips32r5` | `P5600` | the `p5600` core's MSA kernels |
| `i386` | `pentium3` | |
| `i386-max` | `max` under TCG | built for `nehalem`: the SSE2-and-up cores, each object for its own ISA |
| `i386-safe`, `i386-max-safe` | as `i386`, `i386-max` | ReleaseSafe: Debug stores every value to memory, which hides x87 width bugs (T56's reductions failed only optimized) |
| `x86_64` | `max` under TCG | |
| `sparc64` | `sun4u` | ReleaseFast, in 14 images (OpenBIOS maps 4 MB of each) |
| `sparc64-debug` | `sun4u` | Debug, in 46 images: every test, the parameter tests per family and the every-core OpenBLAS suites per level |
| `powerpc64le`, `powerpc64le-power10` | `power9`, `power10` | ReleaseSafe: LLVM's FastISel bug ([`issues.md`](issues.md)) |
| `powerpc64` | `power10` on `pseries` under VOF, big-endian | every 64-bit POWER core; Zig emits ELFv2 here too, so it boots as the little-endian rigs do; ReleaseSafe |
| `powerpc` | G4 on `mac99`, OpenBIOS | the `ppcg4` core; builds in debug |
| `powerpc-750` | `750` on `mac99`, OpenBIOS | RandyOS's iBook G3 target, `powerpc-freestanding-eabi` with CPU `750`: soft float, so no kernel objects (`-Dblas-cores=none`); the portable kernels. Not `750fx`, which OpenBIOS hangs on ([`issues.md`](issues.md)) |
| `powerpc-750-hf` | `750` on `mac99`, OpenBIOS | the 750 with hard float: every 32-bit core linked, none runnable (the PVR names none), so the portable kernels |
| `powerpc-pseries` | the G4 build on a POWER10 in 32-bit mode | power3 to power6, ppc970 and cell as 32-bit code |
| `powerpc-power9` | a 32-bit POWER8 build, less `64bit` and `fpcvt`, on a POWER9 in 32-bit mode | power8 and power9 as 32-bit code (Settled decisions: filling upstream gaps); power10 cannot link, as LLVM fails on 32-bit ISA 3.0 builds ([`issues.md`](issues.md)) |
| `s390x` | `max` | |
| `x86_64-uefi` | `max` on `q35`, OVMF | baseline `x86_64-uefi` build; the app turns on the AVX and AVX-512 state |
| `aarch64-uefi` | `cortex-a72` on `virt`, AAVMF | baseline `aarch64-uefi` build; HVF on an arm64 Mac, SME there |

Stubs with no step: `sparc` (32-bit; LLD cannot link 32-bit SPARC: LLVM 21.1.0 `lld/ELF/Arch`
has `SPARCV9.cpp` only), `powerpc-macho` and
`linux-i386` (TCG's x87 page-boundary bug). The big-endian rigs are `s390x`, `powerpc64`, and
`mips64r2`, `mips64r6`, `mips32r2` and `mips32r5`. Every family has a rig that builds in debug.
`-Dblas-cores` does not change what the rigs link.

Rig mechanics:

- **Results.** Where the machine has semihosting (aarch64, arm, riscv64, MIPS) the runner exits
  with the number of failed tests; elsewhere the step checks the printed summary for
  "0 failed". Consoles log live to `.zig-cache/qemu-test-<name>.log`. `-Dfail-fast` ends a rig
  at its first failing test.
- **Images.** OpenBIOS on `sun4u` maps 0x400000 to 0x800000 only (seen on the rig, `build/qemutest.zig:333`; OpenBIOS's source not checked), so the sparc64 rigs split
  their tests into images by filter (`parts` in `build/qemutest.zig`), and the build, under
  `qemu-check` too, fails an image whose loaded span is past that window (`load_window`,
  `tools/load_window.zig`). `-Dqemu-test-part=<part>` builds and runs only that image.
- **UEFI.** The app boots from a read-only FAT drive (QEMU's vvfat) with the edk2 firmware
  Nix's QEMU ships, as read-only pflash (`tools/nix/run.sh` finds it, so only on the test
  host). It takes its heap and a 1 MiB stack from boot services, turns off the watchdog,
  turns on the vector state, masks interrupts and then calls no firmware: the handlers save
  only 128-bit SIMD registers, which clobbers SVE and SME state, and a boot service unmasks
  interrupts again. So the console is the UART and exit is `ResetSystem`; the summary carries
  the result. It prints its load address for `addr2line`.
- **Progress.** On the rigs (`check.rig()`: freestanding, or `RSTD_TEST_PROGRESS`, which the
  Linux init and the SSH runs set) `eachCore` names each core as it starts, and the biggest
  cases are left to native runs.
- **Test host.** The rigs' QEMU runs on another machine, `-Dtest-host=<ssh alias>`
  (`macmini-m6` by default when rstd is the root project; `-Dtest-host-off` runs it here); the
  images still build here (`tools/test_host.zig`). Each run copies its files to the host's
  `/tmp`, runs QEMU from the Nix flake in `tools/nix` (nixpkgs pinned by `flake.lock`) and
  streams the console back to the watchdog and live log. At most `-Dtest-host-jobs` (8, our default) rigs run
  there at once, across builds. On an arm64 Mac, `aarch64`, `aarch64-uefi` and `linux-aarch64` run under HVF on
  the host's CPU (Apple, SME2, no SVE: the `vortexm4` kernels run too); `-Dtest-host-tcg` keeps
  them on TCG. `x86_64-whpx` and the macOS VM rig stay local.
- **WHPX.** QEMU's model plus FMA4 and XOP, of which WHPX keeps what the host has; `host` and
  `max` stop in SeaBIOS. TCG implements no AVX-512, FMA4 or XOP (QEMU v10.1.0
  `target/i386/cpu.c`: neither `TCG_7_0_EBX_FEATURES` nor `TCG_EXT3_FEATURES` has them), so the `skylakex`,
  `cooperlake` and `sapphirerapids` kernels, and the Bulldozer family's, run only there, on a
  host that has them.
- **macOS VM.** Set up once (`qemu-setup-`, `qemu-install-x86_64-macos`; the header of
  `build/qemutest.zig` has the steps). It runs under emulation, since WHPX on an AMD host shows
  macOS the AMD CPU, which it does not boot on. If the VM does not answer SSH, the test step
  starts it headless (log in `.macos-vm/x86_64-macos/qemu.log`, VNC on 127.0.0.1:5901), prints
  its serial console (`serial.log`) while it boots, and shuts it down after the tests, pass or
  fail, unless `-Dmac-keep`. `qemu-kill-x86_64-macos` kills a VM a build step started.
- **SSH Mac.** `build/sshtest.zig` cross-compiles and copies the tests to the SSH alias; the
  Mac's address and login stay in the SSH config. `ssh-bench-aarch64-macos
  -Dopenblas=<checkout>` races the kernels there (`tools/blasbench/`); the results are in
  [`results-m6.md`](../../../tools/blasbench/results-m6.md).

### What runs where

| Check | Native (`rstd-test`) | Freestanding rigs | Linux rigs | SSH Mac | wasm |
| --- | --- | --- | --- | --- | --- |
| Cores linked | target core in debug, `-Dblas-cores` | all the CPU runs | all the CPU runs | all | the SIMD core; none without SIMD128 |
| `.dat` sweeps | full | header tables | header tables | full | header tables |
| Detection probes | host | ID registers or instructions | auxv, sysfs, `/proc` | `sysctl` | none |
| Threads | yes | no (`single_threaded`) | yes | yes | no |
| `abi.zig` | yes | no | yes | yes | yes |
| `apple_accel.zig` | yes | no | yes | yes | no |
| Level 3 | yes | yes | yes | yes | yes |
| LAPACK | yes | no | yes | yes | yes |

## Implementation plan

1. **Threading** of levels 1 and 2 ([`threading.md`](threading.md)'s plan); level 3 threads.
2. **Zig-idiomatic interface.** Slices, exact error sets with plain names, and documented
   panics ([above](#planned-a-zig-idiomatic-interface)), with the codebase-wide pass.
3. **Remaining level 3.** Landed: 3M, the small-matrix path, and bfloat16 and half-precision
   level 3 (`blas3/half.zig`); `sbgemm`'s small-matrix kernels are not ported.
4. **Extensions.** Landed ([below](#extensions)), `bgemm` with the half gemms.
5. **Per-core tuning.** Act on [`results-m6.md`](../../../tools/blasbench/results-m6.md)'s and
   [`comparisons.md`](comparisons.md)'s recommendations, and convert assembly where `@Vector`
   cannot reach: arm64 SVE, POWER10 MMA, AMX and AVX-512 BF16 for `sbgemm`, and Skylake-X only
   if measured below 85% of its kernel ([`blas3/gemm/notes.md`](blas3/gemm/notes.md)).

## Open issues

- Upstream, toolchain and emulator problems, with their workarounds, are in
  [`issues.md`](issues.md); open numerical and performance questions in
  [`comparisons.md`](comparisons.md).
- The mips `cpucfg` corrections and the sparc `sllx` rewrite are unverified on hardware; ARMv7
  NEON and AltiVec flushing was seen on the rigs, not on hardware.
- Cores the rigs cannot run (Bulldozer, Piledriver and Steamroller microkernels on a host
  without FMA4 and XOP; the AVX-512 cores on a host without AVX-512) are compiled and
  symbol-checked only.
