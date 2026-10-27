# BLAS design

This document specifies how rstd's BLAS is structured: its requirements, the layers from the
public interface down to the per-core kernels, how a core is detected and selected, where every
kernel comes from, and how all of it is tested. It describes the code as it is. Usage is in
[`readme.md`](readme.md), threading in [`threading.md`](threading.md), problems in other
projects in [`issues.md`](issues.md), and the numerical decisions measured against BLIS and the
reference in [`comparisons.md`](comparisons.md).

Upstream references are to the OpenBLAS checkout that `tools/blas_asm` converts, v0.3.34-317
(commit a69062524), unless a file names another; `Vtable.zig` is derived from v0.3.34 (commit
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
| 1 | `blas1/interface.zig` | `copy`, `swap`, `scal`, `scalReal`, `dot`, `dotc`, `dsdot`, `sdsdot`, `nrm2`, `asum`, `sum`, `index`, `extreme`, `rot`, `rotm`, `axpy`, `axpyc`, `axpby`, `rotg`, `zrotg`, `rotmg`, `sbdot`, `bscal`, `sbtobf16`, `bf16to` |
| 2 | `blas2/interface.zig` | `gemv`, `gemvHalf`, `ger`, `gerc`, `symv`, `hemv` |
| 2 | `blas2/triangular.zig` | `trmv`, `trsv`, `tbmv`, `tbsv`, `tpmv`, `tpsv` |
| 2 | `blas2/banded.zig` | `gbmv`, `sbmv`, `hbmv`, `spmv`, `hpmv`, `spr`, `syr`, `hpr`, `her`, `spr2`, `syr2`, `hpr2`, `her2` |
| 3 | `blas3/interface.zig` | `gemm`, `symm`, `hemm`, `syrk`, `herk`, `syr2k`, `her2k`, `trmm`, `trsm` |

Conventions shared by every interface:

- **Element type.** `T` is a comptime parameter: any float type or its `std.math.Complex`
  (`types.Float(T)` must be a float). `BF16` and the mixed half routines have their own entry
  points (`sbdot`, `bscal`, `gemvHalf`).
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
  `Io.Cancelable` everywhere; the root exports them as `Level1Error` to `Level3Error`.
  A vector or matrix that could not fit the address space (`util.fits`, `util.fitsMatrix`)
  returns its stride's or leading dimension's error, so no index computation overflows.
- **Letters.** Transposes, triangles, sides and diagonals are one set of letter enums for every
  level (`types.zig`). Real `r` and `c` transposes mean `n` and `t`. `Trans` has the complex
  forms `o`, `u`, `s` and `d` (conjugated `x`), and `Uplo` the `v` and `m` forms of `hemv`
  (`HEMVREV`). A letter a routine does not take is an error: `InvalidTrans` for `o` to `d`
  outside complex `gemv` and `gbmv`, `herk`'s `t` and complex `syrk`'s `c`; `InvalidUplo` for
  `v` and `m` outside the Hermitian level-2 routines.
- **Leading dimensions.** Level 3 checks them against `max(1, rows)`, as the reference BLAS
  does; OpenBLAS lets `0` through for empty matrices.
- **Quick returns.** Empty problems return before any kernel call; kernels assume `n >= 1`.
  Level 3 does not forward one-column problems to `gemv`.
- **Argument order.** Every routine takes `(allocator, io, T, ...)`, or `(io, T, ...)` without
  an allocator (level 1, `gemvHalf`), as [`threading.md`](threading.md) specifies. `io` is
  unused until the routines thread. `rotg`, `zrotg` and `rotmg` compute from scalars, never
  thread, and take `(T, ...)` and return `void`.

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
  (`exprecisionOn`). It assumes the x87 control word's extended precision, as Linux and x64
  UEFI leave it. Elsewhere `q`/`x` alias `d`/`z`, as in the reference; `blas_types` is a
  membership list for that reason, since a `switch` cannot hold one type twice.
- **`quad_precision`** is off; turning it on makes `XDouble` Zig's `f128` (compiler-rt
  arithmetic) and is exclusive with `exprecision`.
- **`f80` size.** `f80` keeps Zig's native size (16 bytes on x86_64, 12 on x86), where the
  reference builds with `-m128bit-long-double`. Buffer formulas count bytes and round up to
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
slots, 88 integer parameters and the `init` hook. It leaves out the 50 slots only row-major data
or CBLAS reaches (`?gerv_k`, `?gerd_k`, the `*_direct*` kernels, `?omatcopy_k_r*`,
`?imatcopy_k_r*`) and the `sbger_k`, `sbsymv_L` and `sbsymv_U` that nothing upstream fills or
calls; its doc comment lists them (`unused/`, below). The reference guards fields behind `BUILD_*`,
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
  `blas2`, `blas3`) and the generic target's blocking; `finish` runs after any core's install:
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
  LAPACK 3.10's.
- `scal`'s flag argument keeps NaN and infinity for a zero `alpha` when set.
- In the extremes (`i?amax`, `i?amin`, `i?max`, `i?min` and their values) the first NaN wins,
  else the first extreme; complex keys that overflow compare at a quarter, as `izamax.f90` does.
- `rotg` computes in `XDouble`, where the reference uses `long double` on x86; `zrotg` follows
  LAPACK's `zlartg`, since `zrotg.c`'s scaled branches are wrong.
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
  with VFPv4, PowerPC, RISC-V F/D, LoongArch, s390x, MIPS release 6. Not on SPARC, WebAssembly
  or earlier MIPS, whose `madd` rounds the product. Vector multiply-adds also fuse on MIPS with
  MSA and WebAssembly with `relaxed_simd` (`vector_madd`).
- **Per-core policy.** A root may declare `blas_fuse_fma` and `blas_fuse_fma_vector`; the
  per-core objects do (`_kernel.zig`), from the core's `Spec.fuse_fma`. The AMD cores' FMA3/FMA4
  features exist so the assembler accepts their microkernels; their flags leave contraction off.
- **Width.** Floats wider than `f64` never fuse (`fuseFor`): `@mulAdd(f80)` is a soft `fmal`
  call, and the x87 kernels never fuse.
- **Expression rule.** Contraction follows a C compiler's within one expression (clang's
  `-ffp-contract=on`): `a*b + c*d` is `fma(a, b, c*d)`, a longer chain fuses each product into
  the running sum from the left, and a product in its own statement, or converted before the
  add, stays rounded. `rsimd.mulAdd` delegates here.

The portable kernels follow the build target's policy, so two builds for different CPUs may
round a portable kernel differently; reproducibility is per build and core.

### Denormals and the floating-point environment

rstd never changes the floating-point environment (MXCSR, FPCR, FPSCR, `VSCR`) and assumes IEEE
default modes. Denormals are honored everywhere: where a vector unit flushes `f32` denormals
(ARMv7 NEON, pre-VSX AltiVec), `util.vectorIeee` is false and no kernel puts `f32` in its lanes.

- The portable `saxpy` takes its scalar loop there; `vector.zig` has no `s` or `c` roles there.
- The level-3 `f32` tiles are scalar there (`blas3/gemm/kernel.zig`'s `vectorBits`).
- LLVM's loop and SLP vectorizers would still put scalar `f32` code in AltiVec lanes, so the
  AltiVec cores' objects are built without AltiVec (`cpu/_corespec.zig`'s `without`), and the
  portable loops pass each `f32` result they store or accumulate through `util.scalar`, an empty
  asm they cannot widen, compiled out everywhere else ([`issues.md`](issues.md)). ARMv7's
  vectorizers keep off NEON's `f32` lanes themselves.

[`comparisons.md`](comparisons.md) (10) records how this was found and measured.

### Per-core kernel objects

OpenBLAS compiles its kernel tree once per core with that core's flags. rstd does the same with
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
  the runtime both read it, and a test checks every named model and feature exists.
- **Groups.** `build/blaskern.zig` takes `corespec.kernelCores(arch)` (the cores with kernels,
  derived from the registries and checked by a test in `cpu/__root__.zig`), keeps those whose
  kernels run on the target CPU (`corespec.runsOn`), and groups equal specs into one object.
- **Build mode.** The objects always build `-OReleaseFast`, with LLVM (`needsLlvm`: the
  self-hosted x86_64 assembler rejects the converted assembly's register memory operands and
  VEX/EVEX mnemonics, the self-hosted wasm backend emits invalid SIMD, and the other
  architectures have no self-hosted backend that covers them). `-Dblas-llvm` overrides it. They
  link with LLD, except on Mach-O, which LLD cannot write; the self-hosted linker misplaces
  SPARC symbols.
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
| Ported intrinsics or vector C | `<routine>/<arch>/kernel.zig` and siblings (`axpy/arm64/sve.zig`, `axpy/riscv64/rvv.zig`, `axpy/mips/msa.zig`) | per-core object |
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
for `axpy`, `axpyc`, `dot`, `dotu`, `dotc`, `scal`, `nrm2`, `asum`, the extremes, `copy`,
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
- **x86_64, x86, loongarch64:** the routed slots (on x86_64 the peeling ones only on the cores
  whose converted kernels peel; their C kernels elsewhere stay).
- **power:** the routed slots on POWER8 to 10; AltiVec cores keep the portable kernels.
- **arm64:** the routed slots on every core; on `vortex` and `vortexm4` every slot but `rot`,
  where they match or beat the converted assembly on the M6
  ([`results-m6.md`](../../../tools/blasbench/results-m6.md)), so their objects there leave those
  converted kernels out (`vector.displaces`); elsewhere every slot no converted kernel fills
  (below). On Apple's systems, where Accelerate serves most routines, only the routed slots.
- **arm, mips, sparc, riscv64, s390x, wasm:** none.

On arm64 these fill the slots whose upstream kernel is C, which the portable scalar loop served
before; their converted kernels stay, unmeasured:

| Cores | Slots with no converted kernel |
| --- | --- |
| all but `tsv110` | `sdot` |
| all | `cdotu`, `cdotc`, `zdotu`, `zdotc` |
| `thunderx` | `ddot` |
| `thunderx2t99`, `thunderx3t110` | `ddot`, `?nrm2`, `?asum`, `?copy`, `i?amax` |
| `neoversen1`, `ampere1`, `vortex`, `vortexm4` | `ddot`, `?asum`, `?copy`, `i?amax`, `sgemv_n` |
| `neoversen2`, `neoversev2` | those, and `sgemv_t`, `dgemv_n`, `dgemv_t` |
| SVE: `neoversev1`, `armv8sve`, `a64fx`, `armv9sme` | those, and `sswap`, `dswap`, `srot`, `drot` |
| SVE: `cortexa510`, `cortexa710`, `cortexx2` | the same |

### Coverage by architecture

The registries (`kernels/<arch>.zig`) list the routines with per-core kernels; the routed slots
above are not among them:

| Architecture | Level 1 | Level 2 | Level 3 | `@Vector` |
| --- | --- | --- | --- | --- |
| x86_64 | `axpy`, `dot`, `rot`, `scal` (C drivers and assembly); `asum`, `copy`, `sum`, `swap` | `gemv` (C drivers and assembly), `symv`, `hemv` | tile | routed |
| x86 | `axpy`, `asum`, `copy`, `dot`, `rot`, `scal`, `sum`, `swap` | `gemv` | tile | routed |
| arm | `axpy`, `asum`, `dot`, `rot`, `sum`, `swap` | `gemv` | tile | |
| arm64 | `axpy` (assembly, SVE, ThunderX), `asum`, `copy`, `dot`, `rot`, `scal`, `sum`, `swap` | `gemv` | tile | routed and more |
| loongarch64 | `axpy`, `asum`, `axpby`, `copy`, `dot`, `rot`, `scal`, `sum`, `swap` | `gemv`, `symv` | tile | routed |
| mips, mips64 | `axpy` (MSA on 32- and 64-bit), `asum`, `copy`, `dot`, `scal`, `sum`, `swap` | | tile | |
| power | `axpy` (assembly, microkernels, C), `asum`, `copy`, `dot`, `rot`, `scal`, `sum`, `swap` | `gemv` | tile | routed, POWER8 to 10 |
| sparc | `axpy`, `asum`, `copy`, `dot`, `rot`, `sum`, `swap` | `gemv` | | |
| riscv64 | `axpy` (RVV) | | tile | |
| s390x | `axpy` (microkernels) | | tile | |
| wasm | `axpy` | | tile | |

Everywhere else, and for every routine not listed, the portable kernels serve: `ger`, the
triangular, banded and packed level-2 routines (drivers over level-1 kernels), the bfloat16 and
half-precision kernels, and the `q`/`x` slots, which no core replaces (a test in `__root__.zig`
checks this). SPARC exports no level-3 tiles: no SPARC feature changes them, so an object would
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
  kernel once per run of nonzero `y`, which is one call when `y` has no zero.
  [`comparisons.md`](comparisons.md) (5) weighs dropping the skip. At unit stride and small
  `n` (`syr`, `spr`, `syr2` below 100; `spr2`, `zsyr` below 50) the interface runs the driver's
  columns itself, as OpenBLAS's do, keeping the skip (OpenBLAS's `syr2`/`spr2` paths have none).
- **symv and hemv.** `blas2/symmetric.zig` holds the portable kernels, `kernel/generic`'s
  `symv_k.c`, `zsymv_k.c` and `zhemv_k.c`, which compose the running core's copy and `gemv`
  kernels with the same blocking (`SYMV_P`, the vtable's `symv_p`).
- **Buffers.** The kernels take a scratch `buffer`. The interfaces size it as the reference's
  `STACK_ALLOC` does, on the stack up to `interface.stack_bytes` (2048, 32-aligned, the
  reference's `MAX_STACK_ALLOC`), else from the allocator, 4096-aligned.
  [`blas2/buffers.md`](blas2/buffers.md) lists what every kernel does with it.
- **Threading.** The reference's `*_thread.c` drivers are not ported.

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
  touches it, while it is in cache, not in a pass of its own (BLIS's first-`KC` `beta`); `trmm`'s
  and `trsm`'s `alpha` goes where each part of `B` is first read. The same multiplies as a pass
  first, so the same bits; a zero `beta` writes zeros without reading `C`. `alpha == 0` and
  `k == 0` keep the separate pass.
- **Chunks.** The product drivers call the kernel on the same `3 * unroll_n` column chunks they
  pack (`Blocking.chunk`), where the reference hands the kernel a slice's whole `sb`. A stated
  unroll that does not match the kernels then costs speed, never correctness. The triangular
  drivers take whole slices, as the reference's do, and rely on the vtable's unrolls being the
  installed kernels' own.
- **Blocking.** `Blocking.of` reads `P`, `Q`, `R` and the unrolls from the vtable, falling back
  to the generic target's 224, 128 and 4096 where a rule left a value unset or negative.
  `Blocking.split` is `level3.c`'s split of the remaining rows or depth.
- **Packing space.** `Work(T)` sizes `sa`, `sb` and the scratch tile for the problem, each on a
  64-byte line; up to `driver.stack_bytes` (2048) it is on the stack, else from the caller's
  allocator.
- **Kernels.** `blas3/gemm/kernel.zig` holds the packing copies, `beta` and the register-tile
  kernel for real `T` and its `Complex` (`incopy`, `itcopy`, `oncopy`, `otcopy`, `beta`,
  `kernel`, complex `kernel_n`, `_l`, `_r`, `_b`). `symm/kernel.zig`, `trmm/kernel.zig` and
  `trsm/kernel.zig` hold the other routines' copies and kernels over the same tile (`sets`),
  each function named as its slot less the precision letter. The `trsm` kernels are
  `kernel/generic`'s `trsm_kernel_*.c` on the gemm layout's panels, dividing by the diagonal
  where upstream packs its reciprocal, which overflows for a denormal one; `trsv` divides
  alike, complex by Smith's algorithm (`util.cdiv`), so both give the same quotients
  (`comparisons.md`, 10). The `trmm` kernels sum a
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
  size, x86's `r` derives from the packing buffer, `p` rounds up to the unroll, the 3M variants
  take the real ones' `p` and `q` unless the header says otherwise, and LA464 and the AMD
  AVX-512 parts pick their sets by cache size and thread count. `OPENBLAS_BLOCK_FACTOR` scales a
  computed `p` as `driver/others/parameter.c` does; OpenBLAS ignores it in dynamic builds.
  `dtb_entries` per core is `getarch.c`'s `DTB_DEFAULT_ENTRIES`, which `params.zig` also
  extracts into `params_data.zig` (`dtbFor`).
- **Not ported.** `gemm3m`, `symm3m`, `hemm3m` (their slots and wrappers exist), the
  small-matrix kernels, bfloat16 and half-precision level 3, and threading.

[`blas3/gemm/notes.md`](blas3/gemm/notes.md) has the tile measurements and code generation notes
behind these choices; [`comparisons.md`](comparisons.md) (7, 9, 10) the open performance and
robustness questions (fusing `beta`, `trsm`'s inverted diagonal, prefetch).

### bfloat16 and half precision

The kernels are `ext/bf16/kernel.zig`: the conversions, `bscal`, `sbdot` and the
`bgemv`/`sbgemv`/`shgemv` family. They are portable Zig, the same on every core; no core
installs a kernel of its own for any of them, and the per-core sources the reference selects
(AVX-512 BF16, NEON, SVE, BFDOT, VSX, MMA, RVV) are not ported. Every bfloat16 conversion
upstream is `kernel/x86_64/tobf16.c` and `bf16to.c`; the AVX-512 BF16 microkernels give the same
bits.

`kernel/setparam-ref.c` fills the bfloat16-named slots that have no kernel of their own with the
`f32` kernels: 17 level-1 ones (`sbamax_k` to `sbswap_k`) and `dsbdot_k` (with `dsdot_k`), which
the level-1 interfaces reach, and `sbtrsm_*`, `sbtrmm_*`, `sbsymm_*`, `sbneg_tcopy` and
`sblaswp_ncopy`, which nothing calls yet and which stay null while their `f32`
slots are. `portable.finish` sets them at the end of every vtable's build.

Their numerical behavior is in [`readme.md`](readme.md#bfloat16-and-half-precision).

### Apple Accelerate

On macOS the public routines over `s`, `d`, `c` and `z` go to Accelerate where it has the
routine and matches the reference; `apple_accel.zig` binds it and routes to it.

- **Binding.** `load()` opens the framework with `dlopen` (`dll.openDarwinFramework`), never
  `-framework`, so cross-compiles need no SDK, and fills `Table(I)` with `dll.bind`. The
  declarations are ours, from Apple's documentation; no SDK file is copied. The handle stays
  open for the process, since the table points into it (`src/dll.zig`).
- **Symbols.** macOS 13.3 and later have `cblas_dgemm$NEWLAPACK$ILP64`, with 64-bit integers
  (rstd's `isize`), which `load` takes; else the 32-bit `cblas_dgemm$NEWLAPACK`, then the
  classic `cblas_dgemm`. macOS 27 binds 150 of the 154: there is no `cblas_i?amin`. `axpby` is
  `catlas_?axpby`.
- **Threads.** Threaded, Accelerate's `gemv` with `A` un-transposed rounds differently from one
  thread; `ddot`, `dnrm2`, `dasum`, `A^T` `gemv` and `gemm` gave the same bits at every thread
  count tried (macOS 27, M6). So `load` makes Accelerate single-threaded:
  `VECLIB_MAXIMUM_THREADS=1` (`setMaxThreadsEnv`), read at the process's first BLAS call, and
  `BLASSetThreading` (macOS 15+), per thread, which `singleThreadHere` sets on each thread's
  first routed call. One thread cost `gemm` 1.15x at 2000.
- **Loading and the off switch.** `blas.init` loads it (`apple_accel.start`), since `setenv` is
  not thread-safe, unless `RSTD_BLAS_ACCELERATE` is 0 (`Config.accelerate`).
  `blas.useAccelerate` turns routing off and on later; nothing unloads it. `testUse` turns it off,
  so a test's own table runs its kernels.
- **Routing.** Each public routine makes rstd's checks and quick returns, then calls its function
  in `apple_accel.zig`, which calls the table (column-major, complex scalars by pointer,
  `cblas_?dot?_sub`) or returns false for rstd's kernels. On the LP64 table it does so for any
  size, stride or operand span past `c_int`, rather than truncate.
- **Routed:** `copy`, `swap`, `scal`, `dot`, `dotc`, `dsdot`, `sdsdot`, `nrm2`, `asum`, `rot`,
  `rotm`, `axpy`, `axpby`; `gemv`, `gbmv`, `ger`, `gerc`, the real `symv`, `sbmv`, `spmv`,
  `syr`, `spr`, `syr2`, `spr2`, the Hermitian `hemv`, `hbmv`, `hpmv`, `her`, `hpr`, `her2`,
  `hpr2`, and the six triangular routines; all of level 3. On arm64 `sdot` is `cblas_dsdot`
  rounded, and `snrm2` and `scnrm2` the root of `dsdot(x, x)`, so `s` sums in `f64` as rstd's
  kernels do. `d` and `z` `nrm2`, and `asum`, take `incx > 0` only.
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
  - `herk`, `her2k` with one in `C`'s triangle and a nonzero `beta`: past small orders it scales
    `C` by the real `beta` as a complex product.
  - The Hermitian rank updates leave the diagonal's imaginary parts as `beta` times theirs, or
    NaN where an infinity meets them; the reference sets them to zero. They are zeroed after the
    call; `herk` and `her2k` with `beta` 1 and nothing added, which leave `C` as it was, run
    rstd's.
  - x86_64 (under Rosetta 2; no Intel Mac was at hand): `cblas_dsdot` is off by `f32`'s epsilon,
    so `dsdot`, `sdot`, `snrm2` and `scnrm2` run rstd's; `gemm` with NaN or an infinity in `A`,
    which gives NaN at an edge of `m`.

Left to do: tying its threading into rstd's own ([`threading.md`](threading.md#vendor-backends)).

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
| arm64 | Darwin: `sysctl` (Apple, SME); Linux and FreeBSD: `MIDR_EL1` through the kernel's emulation (`HWCAP_CPUID`), else per-CPU sysfs `midr_el1` and `cpu_capacity`, plus the hwcaps; Windows: `MIDR_EL1` from the registry, features from `IsProcessorFeaturePresent`; freestanding and UEFI: the ID registers |
| arm | `AT_HWCAP` VFP and NEON bits and the part in `/proc/cpuinfo`; freestanding and UEFI: CP15 and VFP ID registers |
| power | `AT_PLATFORM`, then the PVR (`mfpvr`, emulated by Linux), then the ISA bits of `AT_HWCAP`/`AT_HWCAP2`; freestanding: the PVR |
| mips | Loongson `cpucfg` where `AT_HWCAP` advertises it, else `/proc/cpuinfo` (Linux only) |
| loongarch64 | `cpucfg` word 0 and `AT_HWCAP` LSX/LASX; freestanding and UEFI: `cpucfg` word 2 |
| riscv64 | `riscv_hwprobe`, else `AT_HWCAP`'s V bit and an RVV 1.0 check, then `vlenb`; freestanding: `misa`, `mstatus` and `vlenb` |
| s390x | `AT_HWCAP`, masked by `LD_HWCAP_MASK` (`probeMasked`) |
| sparc, wasm | nothing at run time: the target decides |

`cpu.choose` (`cpu/_select.zig`'s `choose`) is `gotoblas_dynamic_init` in two steps:

1. **Which core.** `OPENBLAS_CORETYPE`, matched case-insensitively, if it names a core
   (`unknown` names none) and the CPU can run it; else the detector's answer. A name no core has
   is reported at `OPENBLAS_VERBOSE` 1 and detection decides (the reference falls back to its
   generic core). A named core the CPU cannot run is reported and replaced by the closest one it
   can.
2. **Which object.** A core with no kernel object in the build (excluded by `-Dblas-cores`, or
   with no kernels at all) is served by the closest linked core that needs no feature it lacks
   (`select.served`), else by the build's pick for the target, else by the portable kernels.
   The served core's parameters are used too. The result is `Chosen{ core, served }`.

"Closest" is decided by instruction set, not by a table: `select.Sets` derives each core's
feature set from its `Spec` (dependencies included, tuning-only features removed), and `closest`
prefers a core whose set is a subset of the wanted one's, with the most features, ties to the
earlier core. The reference's `#define gotoblas_X gotoblas_Y` chains follow no such rule
(Barcelona to Nehalem puts SSSE3 code on a CPU without it). Tests in `_select.zig` check that
every detectable core has an object that serves it and that no core is served by one that needs
more.

Where the port differs from `dynamic*.c`:

- **x86.** The OS AVX, AVX2 and AVX-512 state checks and the reference's fallbacks and warnings
  are kept. AMD family 10h and 12h are Barcelona (the reference returns NULL); FMA3 is required
  for the Haswell and Zen kernels; `cpuid` leaves are range-checked; macOS's lazy ZMM state is
  not mistaken for no AVX-512.
- **arm.** The reference has no run-time dispatch; `cpuid_arm.c`'s build-time choice is made at
  run time from the hwcap bits, requiring 32 VFP registers for the ARMv7 kernels. Cortex-A9 and
  A15 are detected by part, as `cpuid_arm.c` does.
- **arm64.** The part table follows `dynamic_arm64.c`. Without the `MIDR_EL1` emulation the
  per-CPU sysfs pass reads `cpu_capacity` (decimal), as `cpuid_arm64.c` does;
  `dynamic_arm64.c` reads `cpu_capability`, which Linux does not have. Every core built for
  ARMv8.3 and up needs FCMA, since its `trsm` kernels hold `FCMLA`; Apple's cores always have it.
- **riscv.** `vlenb` picks `riscv64_zvl256b`, which has no object yet, so `riscv64_zvl128b`
  serves it.
- **power.** POWER10 needs ISA 3.1 and MMA. POWER3 to 5, 970, Cell, G4 and 440, static targets
  upstream, are told apart by `AT_PLATFORM`, the PVR or the ISA bits, so none runs a later
  generation's code. A CPU nothing names gets the portable kernels (the reference guesses
  POWER8). AIX is not a Zig target.
- **mips.** The reference tests `HWCAP_LOONGSON_CPUCFG` the wrong way round and its `cpucfg`
  word has both registers zero; both are corrected, unverified on hardware.

### Configuration

`Config.zig` reads every variable the reference reads at init, once, into a value
(`Config.fromEnv`); `blas.config` keeps it. Defaults match the reference's unset behavior. Each
value is read one variable at a time with no allocation, from a map (`blas.init`), a POSIX
`environ` block, or the Windows PEB (`Config.Env`). A value longer than `value_max` (64 bytes)
reads as unparseable and falls back to the default (for a core name, to detection).
Initialization also caches the CPU count in `cpu_count` (1 without threads).

| Variable | Field | Effect here |
| --- | --- | --- |
| `OPENBLAS_CORETYPE` | `core_type_buf` | forces a core (selection, above) |
| `OPENBLAS_VERBOSE` | `verbose` | 1: fallbacks and unknown names; 2: the selected core and the serving object |
| `GET_OPENBLAS_CORETYPE` | `get_openblas_coretype` | POWER only: set to anything, prints the core to stderr |
| `OPENBLAS_BLOCK_FACTOR` | `block_factor` | scales a computed gemm `p` (clamped to 10 to 200 percent) |
| `OPENBLAS_L2_SIZE` | `l2_size_override` | replaces the detected L2 size (KiB) in the blocking rules |
| `OPENBLAS_NUM_THREADS`, else `OPENBLAS_DEFAULT_NUM_THREADS` | `num_threads` | LA464's blocking set; the thread ceiling of [`threading.md`](threading.md) |
| `GOTO_NUM_THREADS` | `goto_num_threads` | the fallback ceiling of [`threading.md`](threading.md) |
| `OMP_NUM_THREADS`, `OMP_ADAPTIVE`, `OPENBLAS_THREAD_TIMEOUT` | `omp_num_threads`, `omp_adaptive`, `thread_timeout` | parsed, no effect |
| `GOTO_DEBUG` | `goto_debug` | parsed, no effect |
| `OPENBLAS_MAIN_FREE`, `GOTOBLAS_MAIN_FREE` | `main_free` | parsed, no effect |
| `LD_HWCAP_MASK` | `ld_hwcap_mask` | s390x only: ANDed with `AT_HWCAP` before detection |
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
| `survey.zig` | the `KERNEL` files and `Makefile.L1`, `L2`, `L3` | the `coverage.md` tables (`blas1/axpy/`, `blas1/dot/`, `blas2/`, `blas3/gemm/`) |
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

- x86_64 `csum`/`zsum` on the AVX-512 cores are `zsum_sse.S`/`zsum_sse2.S`;
- loongarch64 `dsdot` is `dot.S` on LA264 and LA464;
- arm64 `csum` is `csum.S` everywhere;
- MIPS64's MSA complex `axpy`, `rot`, `swap` and `scal` are never selected: `KERNEL` resets
  them to the `arm/` C.

#### Conversion of the assembly

`gen.zig` turns each `.S` into one naked function per kernel, whose body is the source's
assembly, unchanged in substance:

- **Preprocessing.** `zig cc -E` runs once per core and variant with the reference's defines
  (`-DDYNAMIC_ARCH`, the core's `getarch.c` macros, `DOUBLE`, `COMPLEX`, `CONJ`), since OpenBLAS
  selects code with `#ifdef <CORE>` and the ABI with `WINDOWS_ABI` and the like.
- **Variants.** One file per build variant (`Arch.variants`): x86_64 `sysv` and `win`, arm `hf`
  and `sf`, POWER `32`, `32le`, `64` (little-endian ELFv2) and `64be` (big-endian; Zig emits
  ELFv2 there too), SPARC `32` and `64`.
- **Merged texts.** Cores whose expansions are identical share one file; texts that differ in
  under about a tenth of their lines (`similar`) share one too, with the differing lines under
  `parts`, chosen at comptime.
- **Rewrites.** Directives that declare the function or its section go (Zig places it); tabs
  become spaces; named labels become numeric local labels, so none reach the symbol table; x86
  `.align N`, bytes on ELF and COFF but `2^N` on Mach-O, becomes `.balign` (`align_bytes`);
  `uk.zig` does the same for the microkernels. A variant may rewrite every statement
  (`Variant.rewrite`; the sparc `64` variant's `sll` becomes `sllx`).
- **Output.** Each file holds `keys` (core and variant pairs), `serves(core)`, `kernel(core,
  variant)` and the naked `Body`. File names are the role, the source's suffix, the first core
  where one source splits, and the variant: `daxpy_atom_sysv.zig`,
  `daxpy_sse2_barcelona_win.zig`. The routine's `asm/__root__.zig` lists every file by slot;
  `asm_registry.Kernels` turns it into a module for `kernels/<arch>.zig`.

#### Fixes to upstream sources

The per-core tests found bugs in some of OpenBLAS's own assembly, still present on master at
v0.3.34-317. The generator fixes each before converting it: a `Patch` in
`tools/blas_asm/arch/<arch>.zig` is an exact find and replace on the upstream text that must
match `count` times, or the generator stops, with a comment saying what was wrong. The generated
files carry the fixed code. [`issues.md`](issues.md#openblas) holds the issue text for each. No
core takes the extremes' and `nrm2`'s kernels anymore ([`readme.md`](readme.md#one-contract-on-every-core)),
fixed or not, so they are no longer converted (below) and their patches no longer run.

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
  0. The tail now reads neither operand of a zero scalar, as `kernel/arm/axpby.c` does. Double
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

Patches that work around Zig or QEMU rather than OpenBLAS, each removable once the cause is fixed
([`issues.md`](issues.md)):

- `power/zdot.S`, `zdot_cell.S`: they return the complex result in f1/f2 as the ELFv2 ABI says,
  but Zig reads a two-float struct from r3/r4. The sources already copy it into r3/r4 for the
  Fortran interface; the patch turns that on.
- arm `cdot_vfp.S`, `zdot_vfp.S`, `sf` variant: they return through a pointer in r0, as AAPCS
  says, but Zig reads the struct from r0 to r3, so the patch drops the pointer and returns there.
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
| `ext/` | `bf16/` (bfloat16 and half kernels) and the vtable wrappers of `geadd`, `imatcopy`, `omatcopy`, `laswp`, `neg`, not yet behind an interface |
| `_test/` | test-only code, wired by `_test/__root__.zig` (below) |
| `unused/` | git-ignored, kept locally: kernels, slots and stubs no column-major Zig caller reaches (above) |
| `readme.md`, `design.md`, `threading.md`, `comparisons.md`, `issues.md` | documentation |

A routine directory holds its `__root__.zig` (the vtable wrappers, its portable kernels and the
tests of every core's kernels), and per architecture `<arch>/asm/` (generated), `<arch>/uk/`
(generated) and `<arch>/kernel.zig` (composition and hand ports) as the routine needs.
`blas1/axpy/`, `blas1/dot/`, `blas2/` and `blas3/gemm/` have a generated `coverage.md`: which
OpenBLAS source each core uses.

`_test/` holds the code only tests use:

| Path | Role |
| --- | --- |
| `_test/types.zig` | tests of `types.zig`, which the per-core objects compile without `zigstd` |
| `_test/blas1/check.zig` | `eachCore`, the wide reference `Wide`, `Pattern`, `rig`, `sqrtW` |
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
| `build/__root__.zig` | `-Dblas-cores`, `rstd-test`, `rstd-check-macho`, `rstd-check-uefi` |
| `src/dll.zig` | run-time library loading, for Accelerate |
| `src/_test_qemu.zig` | test root of the OS-free rigs |
| `test/qemu/` | boot code, linker scripts and runner of the freestanding rigs |
| `test/linux/init.zig` | PID 1 of the Linux rigs |
| `test/uefi/runner.zig` | the link-only UEFI test runner |
| `tools/blas_asm/` | the converters |
| `tools/blasbench/` | the kernel benchmarks |
| `tools/test_host.zig`, `tools/nix/` | running the rigs' QEMU on another machine |
| `tools/ssh_test.zig`, `tools/qemu_watch.zig` | the SSH test runner; the rigs' watchdog |

### Threading

All routines are single-threaded. [`threading.md`](threading.md) specifies threading: every
routine takes an `io: Io`; levels 1 and 2 split work into chunks fixed by the problem size and
combine partial results in chunk order, so results do not depend on the worker count; level 3
follows OpenBLAS's shared packed `B` protocol on `Io.Group.concurrent`; all scratch is allocated
before any task starts. `Config` already parses the thread variables it reads. `threading.zig`
holds the machinery of its plan's step 1: the cached CPU count and ceiling, the thresholds, the
chunks, the claim loop and the ordered combine; no routine uses it yet.

## Alternatives considered

| Alternative | Reason not chosen |
| --- | --- |
| A drop-in `libopenblas` (CBLAS and Fortran symbols, `xerbla`) | A Zig API with named errors is the goal; an ABI layer would add the row-major mapping and callbacks for no Zig caller |
| One function per precision (`dgemm`, `zgemm`) | Duplicates every interface; a comptime type serves every float type, including those without slots |
| A vtable with fields compiled in or out per build | Every field present and `null` where absent keeps one layout for the main library and every object |
| Selecting kernels lazily per routine | Eager install keeps a kernel call to one indirect call and the table immutable while kernels run |
| A mutex around initialization | Brings threading support into freestanding builds; a three-state handshake suffices |
| The reference's `gotoblas_X` fallback chains | They follow no rule and put SSSE3 code on CPUs without it; the ISA subset rule cannot |
| Falling back to the generic core for an unknown `OPENBLAS_CORETYPE` | Detection gives a faster correct core; the name is reported |
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
  accepts `2 (k + 2)` ulps of the sum of magnitudes, `trsm` its substitution step by step.
  `sqrtW` adds Newton steps where sparc64's `_Qp_sqrt` works in `f64`.
- **Order-free comparisons.** A core may sum in any order, so plain-loop checks allow rounding
  error; complex infinities may meet in another order than the reference's
  (`expectWideAnyOrder`).
- **No flush tolerance.** Denormals must come out as the reference has them on every core.
- **Exactness.** Bits are compared only where a Zig kernel models converted code of the same
  core: the x87 sums, the x86_64 `dot` vector blocks (`skylakexModel`), the half-precision
  `gemv`. The ported C drivers were also compared bit for bit against OpenBLAS's C, compiled with
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
  unroll_m`, 5 and `6 * unroll_n + 1`, so every block and panel edge comes at small sizes, then
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
- **Pending.** `blat3` runs every level-3 routine through the flat API. Still skipped: LAPACK's `potrf`/`potrs` and the `dgesvd`
  regression (`utest.zig`), and `bgemm`, `geadd`, `imatcopy`, `omatcopy`, `gemmt` and complex
  `gemm` in `utest_ext.zig`, which wait for their interfaces.

### Rigs

`src/_test_qemu.zig` is the root of the OS-free tests: `simd` and the whole BLAS tree. The rigs
run it, or the full test suite, where no host here can.

| Rig | How | Covers |
| --- | --- | --- |
| `qemu-test-<machine>` (freestanding) | a test kernel per architecture under `qemu-system-*`: boot code and runner in `test/qemu/`, `single_threaded`, a 32 MiB heap, every core the target CPU runs linked, LLD for the linker script | every BLAS test but `abi.zig`'s and `apple_accel.zig`'s, which need an OS (`default_filters`) |
| `qemu-test-x86_64-whpx` | the x86_64 kernel on the host CPU under WHPX, `SapphireRapids,+fma4,+xop`; x86_64 Windows hosts only | the AVX-512 cores and the Bulldozer family, which TCG cannot run |
| `qemu-test-linux-<arch>` | a Debian netboot kernel with the hosted test executable as PID 1 (`test/linux/init.zig`), `-smp 4` | auxiliary-vector detection, threads, files; aarch64, x86_64, arm, riscv64, powerpc64le |
| `qemu-test-x86_64-macos` | an installed x86_64 macOS VM under TCG, tests over SSH | Mach-O at run time on an x86 CPU model |
| `ssh-test-aarch64-macos` | the full suite on an Apple silicon Mac over SSH (`-Dssh-mac`, default `macmini-m6`) | `vortex`, `vortexm4`, SME without SVE, which QEMU's `-cpu max` cannot show; Accelerate |
| `ssh-test-x86_64-macos` | the x86_64 macOS suite on that Mac under Rosetta 2 | Mach-O x86_64, faster than the VM |
| `wasm-test` | WASI under Node with SIMD128 | `blas1`, `blas2` and `simd` on the wasm kernels |
| `rstd-check-macho` | build only, x86_64 and aarch64 macOS, every core | Mach-O assembly errors |
| `rstd-check-uefi` | build only, x86_64 and aarch64 UEFI, link-only runner | the OS-free tests compile for UEFI |

The freestanding machines (`build/qemutest.zig`):

| Machine | CPU | Notes |
| --- | --- | --- |
| `aarch64` | `cortex-a72` | |
| `aarch64-sve` | `max` | SVE2, SME2, BF16: every arm64 core but Apple's |
| `arm`, `arm-softfp` | `cortex-a15` | `eabihf`, and `eabi` for the kernels' `sf` variant on the VFP the CPU has |
| `riscv64` | `rv64,v=true,vlen=512` | the VLEN of the widest linked core (x280) |
| `loongarch64` | `la464` | |
| `mips64r2el`, `mips64r2` | `Loongson-3A1000` | little- and big-endian |
| `mips64r6el`, `mips64r6` | `I6400` | |
| `mips32r2el`, `mips32r2` | `24Kf` | no 32-bit kernels apply: the portable ones run |
| `mips32r5el`, `mips32r5` | `P5600` | the `p5600` core's MSA `axpy` |
| `i386` | `pentium3` | |
| `x86_64` | `max` under TCG | |
| `sparc64` | `sun4u` | ReleaseFast, in 14 images (OpenBIOS maps 4 MB of each) |
| `sparc64-debug` | `sun4u` | Debug, in 43 images: every test, the parameter tests per family and the every-core OpenBLAS suites per level |
| `powerpc64le`, `powerpc64le-power10` | `power9`, `power10` | ReleaseSafe: LLVM's FastISel bug ([`issues.md`](issues.md)) |
| `powerpc64` | `power10` on `pseries` under VOF, big-endian | every 64-bit POWER core; Zig emits ELFv2 here too, so it boots as the little-endian rigs do; ReleaseSafe |
| `powerpc` | G4 on `mac99`, OpenBIOS | the `ppcg4` core; builds in debug |
| `powerpc-pseries` | the G4 build on a POWER10 in 32-bit mode | power3 to power6, ppc970 and cell as 32-bit code |
| `s390x` | `max` | |

Stubs with no step: `sparc` (32-bit; LLD cannot link 32-bit SPARC), `powerpc-macho` and
`linux-i386` (TCG's x87 page-boundary bug). The big-endian rigs are `s390x`, `powerpc64`, and
`mips64r2`, `mips64r6`, `mips32r2` and `mips32r5`. Every family has a rig that builds in debug.
`-Dblas-cores` does not change what the rigs link.

Rig mechanics:

- **Results.** Where the machine has semihosting (aarch64, arm, riscv64, MIPS) the runner exits
  with the number of failed tests; elsewhere the step checks the printed summary for
  "0 failed". Consoles log live to `.zig-cache/qemu-test-<name>.log`. `-Dfail-fast` ends a rig
  at its first failing test.
- **Progress.** On the rigs (`check.rig()`: freestanding, or `RSTD_TEST_PROGRESS`, which the
  Linux init and the SSH runs set) `eachCore` names each core as it starts, and the biggest
  cases are left to native runs.
- **Test host.** The rigs' QEMU runs on another machine, `-Dtest-host=<ssh alias>`
  (`macmini-m6` by default when rstd is the root project; `-Dtest-host-off` runs it here); the
  images still build here (`tools/test_host.zig`). Each run copies its files to the host's
  `/tmp`, runs QEMU from the Nix flake in `tools/nix` (nixpkgs pinned by `flake.lock`) and
  streams the console back to the watchdog and live log. At most `-Dtest-host-jobs` (8) rigs run
  there at once, across builds. On an arm64 Mac, `aarch64` and `linux-aarch64` run under HVF on
  the host's CPU (Apple, SME2, no SVE: the `vortexm4` kernels run too); `-Dtest-host-tcg` keeps
  them on TCG. `x86_64-whpx` and the macOS VM rig stay local.
- **WHPX.** QEMU's model plus FMA4 and XOP, of which WHPX keeps what the host has; `host` and
  `max` stop in SeaBIOS. TCG implements no AVX-512, FMA4 or XOP, so the `skylakex`,
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
| Cores linked | target core in debug, `-Dblas-cores` | all the CPU runs | all the CPU runs | all | the SIMD core |
| `.dat` sweeps | full | header tables | header tables | full | not run |
| Detection probes | host | ID registers or instructions | auxv, sysfs, `/proc` | `sysctl` | none |
| Threads | yes | no (`single_threaded`) | yes | yes | no |
| `abi.zig`, `apple_accel.zig` | yes | no | yes | yes | no |
| Level 3 | yes | yes | yes | yes | no |

## Implementation plan

1. **Threading.** [`threading.md`](threading.md)'s plan.
2. **Zig-idiomatic interface.** Slices, exact error sets with plain names, and documented
   panics ([above](#planned-a-zig-idiomatic-interface)), with the codebase-wide pass.
3. **Remaining level 3.** `gemm3m`, `symm3m`, `hemm3m`, the small-matrix kernels, bfloat16 and
   half-precision level 3.
4. **Extensions.** Interfaces for `geadd`, `imatcopy`, `omatcopy`, `laswp`, `neg` and `gemmt`,
   which the pending `utest_ext` cases wait for.
5. **Per-core tuning.** Act on [`results-m6.md`](../../../tools/blasbench/results-m6.md)'s and
   [`comparisons.md`](comparisons.md)'s recommendations, and convert assembly where `@Vector`
   cannot reach: arm64 SVE, POWER10 MMA, AMX and AVX-512 BF16 for `sbgemm`, and Skylake-X only
   if measured below 85% of its kernel ([`blas3/gemm/notes.md`](blas3/gemm/notes.md)).

## Open issues

- Upstream, toolchain and emulator problems, with their workarounds, are in
  [`issues.md`](issues.md); open numerical and performance questions in
  [`comparisons.md`](comparisons.md).
- `riscv64_zvl256b` has no object; `riscv64_zvl128b` serves it.
- The mips `cpucfg` corrections and the sparc `sllx` rewrite are unverified on hardware; ARMv7
  NEON and AltiVec flushing was seen on the rigs, not on hardware.
- Cores the rigs cannot run (Bulldozer, Piledriver and Steamroller microkernels on a host
  without FMA4 and XOP; the AVX-512 cores on a host without AVX-512) are compiled and
  symbol-checked only.
