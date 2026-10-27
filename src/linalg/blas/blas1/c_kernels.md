# Porting OpenBLAS's C level-1 kernels

Some OpenBLAS kernels are C, not assembly: a driver (`kernel/<arch>/ddot.c`)
around a block kernel, which is an inline-assembly microkernel
(`ddot_microk_haswell-2.c`) or C intrinsics, or a driver with no block at all.
`axpy` was the first routine ported this way; `dot` on x86_64 is the second.
Both go through the same machinery, so a new routine adds a table row and its
driver, not a new converter. Where these ports sit among the other kernel kinds
is in [`design.md`](../design.md#kernel-kinds).

## The parts

| Part                  | Where                                   | Made how               |
| --------------------- | --------------------------------------- | ---------------------- |
| Inline-assembly block | `<routine>/<arch>/uk/<p><r>_<core>.zig` | generated, `uk.zig`    |
| Index: core to block  | `<routine>/<arch>/uk/__root__.zig`      | generated, `uk.zig`    |
| Driver                | `<routine>/kernel.zig`                  | hand-ported            |
| Intrinsics block      | `<routine>/<arch>/kernel.zig`           | hand-ported, `@Vector` |
| Composition per arch  | `<routine>/<arch>/kernel.zig`           | `uk.Kernels(...)`      |
| Registry              | `kernels/<arch>.zig` `modules`          | one line per routine   |
| Tests                 | `<routine>/__root__.zig`                | hand-written           |

`blas1/uk_block.zig` holds what is shared: `Entry` (a row of the index),
`find`, the `Block` kinds, and `Kernels`, which exports the core's kernels
(`blaskern_<role>_<core>`) and points the vtable at them. A core whose source is
a whole `.S` kernel still goes through `asm/__root__.zig`; `Kernels` takes both
and prefers the assembly.

## What the generator converts

`tools/blas_asm/uk.zig`, run as

```bash
OPENBLAS=<checkout> zig run tools/blas_asm/uk.zig [-- <routine>...]
```

For each routine in the `Routine` table of `kernel_map.zig` and each target
that lists it, it:

1. converts every `<p><routine>_microk_<core>-2.c` (or, for `zarch`, the
   function in the driver's own file): the `__asm__` statement as written,
   named operands, clobbers at the width of the registers used, `#if CONJ`
   resolved once per form;
2. reads the `#if defined(CORE)` chain at the top of each driver for which
   block each core includes and the driver's block size (`n & -N`);
3. writes the index. A core that uses the driver but no chain link gets
   `.none`; a link whose file has no assembly gets `.vector`.

Which cores use the driver comes from `KERNEL.<core>`, then `KERNEL`, then
`Makefile.L1`, through the roles' make variables (`SDOTKERNEL`, `DSDOTKERNEL`).
It does not edit what it converts, so a mistake upstream is fixed by a patch.

### Patches

A target's `patches` are `arch.Patch` rows, the same type `gen.zig` uses for the
`.S` kernels, matched on the microkernel's file name. `find` must occur
`count` times or the generator stops. List each in
[`design.md`](../design.md#fixes-to-upstream-sources)'s "Fixes to upstream
sources". The x86_64 `dot` blocks the host can run needed none; the
Bulldozer, Piledriver and Steamroller ones are untested at run time.

### Adding a routine

1. A `Routine` and its roles in `kernel_map.zig`: the role's name, precision
   letter, make prefix (`DS` for `DSDOTKERNEL`), and which function of a block
   it calls (the conjugated `axpyc` calls the second).
2. The routine in a target's `routines` in `uk.zig`, then run it.
3. The driver, `<routine>/kernel.zig` (below), with `roles`, `Kernel(role)`,
   `Asm(role, f)` and `entry(role, mask, Block)`.
4. A `<routine>/<arch>/kernel.zig` that is a `uk.Kernels(...)` plus the
   architecture's `Vector` block, and its line in `kernels/<arch>.zig`
   (replacing the routine's plain `registry.Kernels(...)` line).
5. Tests; `survey.zig -- <routine>` for `coverage.md`.

## What is ported by hand: the driver

The driver is the C function around the block, ported statement by statement.
It lives once per routine, not per architecture, since the drivers of one name
are the same C across cores.

- **Shape.** `n <= 0` returns first. The contiguous branch is `inc == 1` on both:
  `n1 = n & -mask`, the block for `n1`, then the scalar tail from `n1`. Anything
  else is the strided loop, with the C's own unroll and its pairing of
  products (`temp1 += m1 + m3`): reordering a sum changes the bits.
- **Block ABI.** The C hands the block a pointer to where the sum goes, zero
  beforehand. `Block.run(n1, x, y)` returns what the block left there: one
  element, or the array of partial sums for complex (`cdot` leaves eight, `zdot`
  four). `Asm` wraps a converted function that way.
- **No block.** The C falls back to a generic loop in the same file, which has
  its own order (`dot += y0*x0 + ... + y7*x7`), not the scalar tail's. That is
  `Fallback`; the index's `.none` means it.
- **Contraction.** The build states whether the core's flags contract
  multiply-adds (`fma.zig`, `fma.fuse`), and clang fuses within one expression:
  `a*b + c` fuses; a product kept in its own statement (`m1 = y*x`) does not; in
  `a*b + c*d` the first product fuses. A `float` product added to a `double`
  is rounded as `float` first, so `sdot`'s tail is unfused where `dsdot`'s is
  not. The driver calls `fmadd(fuse, T, a, b, c)` where, and only where, the C
  statement would contract.
- **Mixed types.** `sdot` sums in `double` and returns `float`; `dsdot` is the
  same file built with `DSDOT` and calls the block in chunks of 32 whatever
  the block size. The role table carries that (`dsdot` is not `sdot`), and
  `uk.Float` gives the element type (`f32` for `dsdot`).
- **Complex.** Results are `types.CComplex(F)` (an `extern struct`, so the C ABI
  returns it as the reference does). `*dotc` differs from `*dotu` only in the
  last two lines.
- **Threads.** The C `CNAME` splits long vectors over threads
  (`#if defined(SMP)`) and adds the partial sums. That is the interface's job
  here; the driver is the single-threaded `dot_compute`, which is also what
  OpenBLAS runs below 10000 elements.
- **Strides.** `util.at` reaches negative offsets; the C indexes with
  `BLASLONG`, so the driver does.

## Kernels written in intrinsics

Where an architecture's block is C intrinsics (`ddot_microk_skylakex-2.c`, AVX-512
cores) there is nothing to convert; the index says `.vector` and the
architecture's `kernel.zig` supplies `Vector(role)`, a block on `@Vector`.

- Keep the C's operation order, not just its value. Reductions are fixed:
  accumulators in a fixed number, folded in a fixed order, halves added, then
  `hadd`. Written the same way, the result is the same bits; written as
  `@reduce`, it is not.
- Mirror the intrinsics: `_mm256_extractf128_pd(a, 1)` is the upper half
  (`@shuffle`), `_mm_hadd_pd(h, h)[0]` is `h[0] + h[1]`, `_mm_hadd_ps` twice is
  `(h0 + h1) + (h2 + h3)`.
- Mirror the preprocessor choice. `#ifdef __AVX512CD__` is the target's
  `avx512cd` feature, read in the per-core object, which is built with the
  core's features. The 512-bit accumulators fold into the 256-bit ones before
  the 256-bit loop.
- Take `fuse` as a parameter of the block (`VectorBlockWith`) so a test can run
  both policies. `acc += a * b` on vector types contracts the same way as on
  scalars.
- Wider vectors than the host has still run, lowered to narrower ones, so the
  512-bit path is testable anywhere against a lane-by-lane model of the C.

Other families the same way: NEON and SVE (`arm_neon.h`, `arm_sve.h`), RVV
(`__riscv_v*`), MSA, z/Architecture vector builtins and POWER `__vector` types
map to `@Vector` and `@shuffle`; where the C uses a fused instruction by name
(`vfmaq`), use `@mulAdd` unconditionally, not `fmadd(fuse, ...)`.

## Tests

- `check.eachCore` runs the routine's check against every linked core the
  machine can run, over `check.sizes` and `check.strides`, so each block size
  and unroll boundary is hit. A failing core is named and the rest still run.
  Build with `-Dblas-cores=all`, or only the host's core is linked.
- Tolerances: a core may sum in any order, so the plain-loop check allows rounding
  error. `dsdot` through an upstream block sums 32 `f32` products in `f32`, so
  it is checked at `f32` precision; the portable one stays strict.
- Where the order is fixed, check the bits: the vector block against a scalar
  model of the C (`skylakexModel` in `dot/__root__.zig`), for each policy.
- A driver can be checked against OpenBLAS's own C by compiling the original
  with a stub `common.h` (`BLASLONG`, `FLOAT`, `CNAME`, the complex struct) and
  `-ffp-contract=on` or `off`, then comparing results bit for bit with the
  ported driver over the same data. Done for `dot` on the cores the host runs
  (6264 cases, no differences).
- Cores the host cannot run (FMA4 and XOP for Bulldozer and Piledriver,
  AVX-512 for the Skylake-X family) are compiled and symbol-checked natively;
  running them needs `qemu-test-x86_64` with a matching `-cpu`.

Zig exits 0 on an inline-assembly error and writes an empty object, so check
that a new object really exports the symbols (`llvm-nm`).
