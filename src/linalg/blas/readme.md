# BLAS

Port of the OpenBLAS kernel tree; the converters are in `tools/blas_asm/`.
Which OpenBLAS source each core uses is in the generated `coverage.md` files
next to each routine.

## Layout

- `Vtable.zig`: the reference's `gotoblas_t`; `blas.kernel(.slot)` calls a slot.
- `blas{1,2,3}/`, `ext/`: one directory per routine. Its `__root__.zig` holds
  the interface (stubs until ported), and `<arch>/asm/` the converted
  assembly; `blas{1,2,3}/interface.zig` are the Zig entry points.
- `kernels/<arch>.zig`: the one list of an architecture's kernel modules, and
  the `install` that points the vtable at a core's kernels.
- `cpu/<arch>.zig`: detection, and `canRun` where the tests need it.
- `blas2/buffers.md`: what each level-2 kernel does with its `buffer`, and
  what the interface gives it.
- `blas1/uk_block.zig`: a routine's C driver around a converted microkernel
  block, for routines whose upstream kernel is C (`axpy`, `dot`);
  `blas1/c_kernels.md` says how such a kernel is ported.

## Initialization

`blas.init(env)` selects the core and fills the vtable with its kernels, once
and thread safe; the first kernel call does it by itself if the program has
not. A routine fetches the table once per call (`blas.vtable()`, or
`blas.kernel` for one kernel) and code under it reads `blas.table()`. On its
own, init reads `OPENBLAS_*` from the process (the Windows PEB, or libc's
`environ`, `_NSGetEnviron()` on Darwin) one variable at a time; Zig gives a
library no environment elsewhere, so there only `blas.init(env)` honors them.
`blas.shutdown()` goes back to the portable kernels; the next call selects
again. `DTB_ENTRIES` per core is in `params_data.zig`
(`tools/blas_asm/params.zig`, from `getarch.c`).

## Not a drop-in

This is a Zig library, not a replacement `libopenblas`: there are no CBLAS
(`cblas_*`) or Fortran-named (`dgemv_`) entry points, no row-major mapping
layer, and no `xerbla` callback. The interfaces are Zig functions
(`blas1/interface.zig`, `blas2/`) that take the Fortran argument order and
column-major data; bad arguments return a named error in the order the
reference's `xerbla` numbers them. Chosen on purpose, so the audit gaps for
those layers are closed by this decision, not by code.

## Level 3

Only real `gemm` is ported so far (`sgemm`, `dgemm`, and `qgemm` for
`XDouble` through the same code); complex, the other routines, the
small-matrix kernels and threading come later.

- **Interface:** `blas3/interface.zig`'s `gemm` checks its arguments in the
  reference's `xerbla` order (leading dimensions against `max(1, rows)`, as
  the reference BLAS does), returns for empty `C`, and calls the driver.
  Sizes and leading dimensions are `isize`; real `r` and `c` transposes are
  `n` and `t`. No `gemv` forwarding for one-column problems.
- **Driver:** `blas3/level3.zig` is `driver/level3/level3.c`'s loops, with
  `P`, `Q`, `R` from the vtable (`gemm/params.zig`). The packing space comes
  from the caller's allocator, sized to the problem (2 KiB or less on the
  stack). `B` is packed in chunks of `3 * unroll_n` columns and the kernel is
  always called per chunk, so a stated unroll that does not match the
  kernels costs speed, never correctness.
- **Kernels:** `blas3/gemm/kernel.zig`, Zig on `@Vector`s, built into the
  main library (the portable kernels) and into every per-core object, so
  each core's features set the vector width and register count, and from
  them the tile (`geometry`). They never make a core need an object of its
  own. The objects export their unrolls as data; `cpu/<arch>.zig` sets the
  core's blocking before its kernels, and `portable.finish` restores the
  portable kernels' unrolls where those stayed, so the vtable's unrolls are
  always the installed kernels' own. `gemm/notes.md` has the measurements.
- **Tests:** `blas3/check.zig` compares with a plain gemm in twice the
  precision, within `2 (k + 2) eps` of the sum of magnitudes, on a vtable
  whose blocking is shrunk so every block and panel edge comes at small
  sizes, and for every linked core this machine runs; then special values in
  each operand, and shapes past the core's own `R` (natively, past twice its
  `P` and `Q` too).

## Converted assembly

`tools/blas_asm/gen.zig` turns each upstream `.S` into a naked function.

- **ABI variants:** a source is preprocessed once per variant: x86_64 `sysv`
  and `win`, arm `hf` and `sf`, POWER `32`, `32le`, `64` (little-endian
  ELFv2) and `64be` (big-endian; Zig emits ELFv2 there too), SPARC `32` and
  `64`.
  `abi.zig` picks the one a target needs from its C calling convention, so
  x86_64 UEFI gets `win`. The build stops for a target whose ABI has none
  (x32, MIPS n32); `-Dblas-cores=none` leaves the portable kernels.
- **Merged texts:** cores whose texts differ in under 10% of lines share one
  file, with the differing lines chosen at comptime.
- **Labels:** named labels become numeric local labels, so none reach the
  symbol table.
- **Alignment:** x86 `.align N` means N bytes on ELF and COFF but 2^N on
  Mach-O, so the x86 sources' `.align` is written as `.balign`, bytes on every
  object format. `uk.zig` does the same for the C drivers' blocks.

## Kernel selection

Which source a core uses follows `kernel/Makefile`: the core's
`KERNEL.<core>` (and what it includes), then the architecture's `KERNEL`, then
`Makefile.L1` and the like. A plain assignment in `KERNEL` overrides the core's
choice; only an `ifndef` leaves it. `kernel_map.zig` evaluates this, deciding
`ifdef`/`ifndef` on the variable itself and listing both sides of any other
conditional as alternatives. Upstream consequences: x86_64 `csum`/`zsum` on the
AVX-512 cores are `zsum_sse.S`/`zsum_sse2.S`; loongarch64 `dsdot` is `dot.S`
on LA264 and LA464; arm64 `csum` is `csum.S` everywhere; MIPS64's MSA complex
`axpy`, `rot`, `swap` and `scal` are never selected, `KERNEL` resets them to
the `arm/` C.

One routine is left portable on purpose: `nrm2` on POWER, Loongson
(`mips64`), `loongarch64` and SPARC. Their kernels are converted but not
wired in (`kernels/<arch>.zig`): they gave NaN for an infinity, 0 or values
far off for denormals (NaN on ppcg4 and ppc440), and the `lsx` and `lasx`
ones 0 for a NaN (`issues.md`). They are the same scaled one-pass loop as
the portable kernel, so nothing is lost.

## Portable kernels

The portable kernels are Zig algorithms with the reference's documented
semantics, not copies of its C: `dsdot` accumulates in `f64` (so does `sdot`),
`nrm2` scales against overflow and underflow and is infinite for infinities
(LAPACK 3.10's; the scaled loop alone gives NaN for two), `scal`'s flag keeps
NaN and Inf, `i?amax` returns the first of equal extremes, and a NaN first
element is never replaced. Not required by the design, so not done: reproducing the
reference's accumulation order (`generic/dot.c`'s unrolled and SIMD partial
sums, `sum.c`'s), its signed-zero handling for `axpby` with a zero `alpha`,
its threading, and `zrotg.c`'s internals (`zrotg` follows LAPACK's `zlartg`;
`zrotg.c`'s scaled branches are wrong). `rotg` computes in `XDouble`, which
is where the reference uses `long double` on x86: the x87 `f80` where
`types.exprecision` is on, else `f64`. Wider than `f64`, the reductions and
`gemv` do take the order of the x87 `q`/`x` kernels, so they give those
kernels' bits (`blas1/x87.zig`, `blas2/gemv/x87.zig`). Tests compare against a
wide-precision reference with tolerances; exact comparison is kept where a Zig
kernel stands in for a converted kernel of the same core. `blas1/edge.zig`,
`blas2/edge.zig` and `blas3/check.zig` put NaN, infinities, signed zeros,
denormals and huge and tiny values through every core's kernels, with NaN
wherever a routine must not read or write; where a core's upstream kernels
do otherwise, `issues.md` says what they do, and the tests check that.

Contraction (`fma.zig`) is decided per target: fused where the ISA has a
scalar fused multiply-add (x86 with FMA, AArch64, ARM with VFPv4, PowerPC,
RISC-V F/D, LoongArch, s390x, MIPS release 6), not on SPARC, WebAssembly or
earlier MIPS, whose `madd` rounds the product. Vector multiply-adds also fuse
on MIPS with MSA and WebAssembly with `relaxed_simd`. A per-core kernel
object overrides the policy with that core's flags.

## Symbols and linking

Kernels are exported as `blaskern_<role>_<core>` with hidden visibility, so
they never leave a shared library rstd is built into. COFF has no visibility;
a DLL exports nothing without `dllexport`.

The qemu rigs link with `lld`, since the self-hosted ELF linker ignores linker
scripts. LLVM-built kernel objects do too, since the self-hosted linker
misplaces their symbols, except on Mach-O, which `lld` cannot write.

## Not supported by Zig

These OpenBLAS platforms and cores are not converted; the portable kernels
serve them where they are targets at all.

| Platform or core          | Why                                                  |
|---------------------------|------------------------------------------------------|
| alpha                     | no Zig backend                                       |
| ia64                      | no Zig target                                        |
| csky, e2k                 | no Zig target; upstream uses the generic C kernels   |
| POWER BG/P (double FPU)   | LLVM cannot assemble the "hummer" instructions       |
| RISC-V C910V              | LLVM has no `xtheadvector`                           |

## bfloat16 and half precision

The kernels are `ext/bf16/kernel.zig` (conversions, `bscal`, `sbdot`, and the
`bgemv`/`sbgemv`/`shgemv` family); the interfaces are `blas1/interface.zig`
(`sbdot`, `bscal`, `sbtobf16`, `bf16to`) and `blas2/interface.zig`
(`gemvHalf`). They are portable Zig, the same on every core: no core installs
a kernel of its own for any of them, and the per-core sources the reference
selects (AVX-512 BF16, NEON, SVE, BFDOT, VSX, MMA, RVV) are not ported.

- **Aliases:** `kernel/setparam-ref.c` fills the bfloat16-named slots that
  have no kernel of their own with the `f32` kernels: 17 level-1 ones
  (`sbamax_k` to `sbswap_k`) and `dsbdot_k` (with `dsdot_k`), which the
  level-1 interfaces reach, and `sbger_k`, `sbsymv_*`, `sbtrsm_*`, `sbtrmm_*`,
  `sbsymm_*`, `sbneg_tcopy` and `sblaswp_ncopy`, which nothing calls yet and
  stay null while their `f32` slots are. `portable.finish` sets them once, at
  the end of every vtable's build (`generic_vtable`, `vtableFromCoreType`).
- **Conversions:** round to nearest even; denormals flush to a signed zero and
  NaNs come out quiet, in both directions (`x86_64/tobf16.c`, `bf16to.c`).
- **`bscal`:** `x := alpha*x` through `f32`. The flag `dummy2` is 0 for a
  zero `alpha` giving zero, 1 for it keeping NaN and infinity as NaN; as in
  `generic/scal.c` that test is on the bits of `alpha`, so only +0 counts.
- **`sbdot`:** both operands widened to `f32`, the products and their sum in
  `f64`, returned as `f32`; no buffer. Differs from the reference's `sdot` on
  the widened data (and from the AVX-512 BF16 kernel) by rounding only.
- **Half `gemv`:** each result sums `a*x` in `f32` in index order, unfused,
  and rounds once on the store (`alpha*acc + beta*y`, `beta == 0` never reads
  `y`). Without a transpose, rows go in blocks of 256 with the columns
  outermost, so `A` is read down its columns. For `bgemv` the scalars and `y`
  are bfloat16. With `m` or `n` zero the kernel yields `beta*y`; the
  interface returns first.
- **Interface:** `gemvHalf` checks arguments as `interface/sbgemv.c` does, and
  for `alpha == 0` only scales `y` by `beta`. The reference's threaded drivers
  are not ported. `sbdot` and the others return at once for `n <= 0`.

## Fixes to upstream sources

The per-core tests in the qemu rigs found bugs in some of OpenBLAS's own
assembly files (still there on OpenBLAS master, v0.3.34-317). The generator
fixes each one before converting it: `patches` in `tools/blas_asm/arch/`, an
exact text edit that must match, with a comment saying what was wrong. The
generated files carry the fixed code.

- `x86_64/zdot_atom.S`, Windows: Windows passes the address of a 16-byte
  result first, shifting the arguments along; the kernel read that address
  as `n`. Fixed the way `zdot_sse2.S` already does it.
- `x86_64/zgemv_t_atom.S`, Windows: read the double `alpha_i`, which Windows
  passes on the stack, with a single-precision load; its sibling kernels use
  the double one.
- `x86_64/scal_sse.S`, `scal_sse2.S`, `scal_atom.S`, Windows: read the flag
  (whether a zero `alpha` keeps NaN) from `incy`'s stack slot, `64(%rsp)`,
  not the tenth argument's, `80(%rsp)`, so a zero `alpha` stored zeros over
  NaN and infinity where the interface asks for NaN. `scal_sse.S` also took
  a NaN `alpha` for zero (`comiss` sets ZF for NaN; `scal_sse2.S` has the
  `jp`), so with flag 0, as for a level-2 `beta`, NaN became zeros.
- `x86_64/zscal_atom.S`: never read the flag at all, so a zero `alpha`
  always stored zeros; it now tests it as `scal_atom.S` does. A NaN
  imaginary part beside a zero real one also counted as zero.
- `x86_64/iamax_sse.S`, `iamax_sse2.S`, `izamax_sse.S`, `izamax_sse2.S`
  (every `i?amax`, `i?amin`, `i?max`, `i?min` on x86_64): the index search
  takes the last element as the match unread, but counted it even when it
  had been read, so for a NaN extreme, which no element equals, returned
  `n + 1`. It now returns `n` there.
- `x86_64/copy_sse.S`, `copy_sse2.S`, `zcopy_sse.S`, `zcopy_sse2.S`
  (opteron, opteron_sse3): loaded by adding to a zeroed register, so -0
  came out +0 and a signaling NaN quiet. They now load as on the other
  cores, so the copy is exact.
- `x86_64/symv_L_sse.S` (`ssymv_L` on atom, barcelona, bobcat, core2,
  dunnington, nano, opteron, penryn, piledriver, prescott): the last two
  rows of a block, and the last two columns, zeroed the upper half of the
  matrix operand but not of `x`'s, which holds the next elements, so an
  infinity there met a 0 and gave NaN in `y`.
- `power/scal.S`, `zscal.S` (power6, ppc970, cell): on Linux, read the flag
  argument (whether a zero `alpha` keeps NaN) at ELFv1's stack offset. Under
  ELFv2, and 32-bit, that slot lies past the caller's arguments, so a zero
  `alpha` could multiply where it should store zeros: `cgbmv` with `beta`
  zero gave NaN. The fix reads it through `FRAMESLOT`, which has each ABI's
  offset; FreeBSD's offsets were already right.
- `x86/amax.S`, `iamax.S`, `zamax.S`, `izamax.S`, the minima: took the new
  element unless it compared above the pick, and an unordered compare does
  not, so a NaN replaced the pick and the next element replaced the NaN.
  They now compare the pick with the new element, as the maxima do.
- `arm/amax_vfp.S`, `iamax_vfp.S`, the minima: moved on `lt`, which an
  unordered compare sets too, with the same effect; `mi` is less than only.
- `mips64/iamax.S`, `imax.S`, `izamax.S` and their minima (loongson3r3,
  loongson3r4), and `sparc/iamax.S`, `imax.S`, `izamax.S`: keep a pick per
  lane, and of equal picks the reduction kept the first lane's, which can
  be the later index. They now take the smaller index of equal picks.

The rest are the `lsx` and `lasx` kernels under `loongarch64/`, mostly where
a stride is not 1:

- `sum_lsx.S`, `sum_lasx.S`, `csum_lsx.S`, `csum_lasx.S`: the double result
  was returned with a single-precision move, leaving its upper half undefined
  (NaN under qemu). `sum_lasx.S` also added doubles as floats.
- `casum_lasx.S`, `csum_lasx.S`: single precision's final sum added the last
  block of eight twice; `casum` also added lane 6 in place of lane 7.
- `axpby_lsx.S`, `caxpby_lsx.S`: with `beta` zero and a stride on `y`,
  the leftover elements were written over the first ones.
- `axpby_lasx.S`: double precision skipped every other block of eight, and
  with `beta` zero never stopped.
- `axpby_lsx.S`, `axpby_lasx.S`, `caxpby_lsx.S`, `caxpby_lasx.S`, with a
  zero scalar: the elements past the last block were `alpha*x + beta*y`
  whatever the scalars, from `x` (`y`) elements the zero-`alpha` (`beta`)
  blocks left their pointer at, so a NaN or infinity met the 0. The tail
  now reads neither operand of a zero scalar, as `kernel/arm/axpby.c` does.
  Double `axpby` also compared the scalars with a single-precision 0
  (`ffint.s.l`, upper half undefined), so it never took the zero-scalar
  blocks; and complex, with one stride 1 and the other not, has no such
  blocks, so a zero scalar now goes to the strided ones, which do.
- `caxpby_lsx.S`, `caxpby_lasx.S`: with `alpha` and `beta` zero and a
  stride on `y`, zeroed only the real parts, and `lasx` single one of them twice.
- `caxpby_lsx.S`: double, with a stride on `x` and none on `y`, stored the
  upper half of each block's second pair over both halves.
- `max_lsx.S`, `min_lsx.S`: with a stride on `x`, compared by magnitude, as
  `amax` does.
- `camin_lsx.S`, `camin_lasx.S`: with a stride on `x`, the minimum started from
  `|re(x0)|` rather than `|re(x0)| + |im(x0)|`.
- `imax_lasx.S`, `imin_lasx.S`: with a stride on `x`, single precision below 8
  compared against a first element never loaded, and double counted the
  indices up twice as fast.
- `iamin_lsx.S`, `iamin_lasx.S`: with a stride on `x`, double kept a minimum by
  the previous comparison, and `lasx` skipped the reduction across lanes.
- `icamin_lsx.S`, `icamin_lasx.S`: `lsx` below 4 elements started the index
  at 4 rather than 1, and with a stride on `x` its single precision started
  from an absolute value never computed; `lasx` single, with a stride on `x`,
  overwrote its index step, and its reduction over lanes 4 to 7 mistook its
  own partial minimum for lane 4's value, returning lane 4's index.
- `snrm2_lsx.S`, `snrm2_lasx.S`, `cnrm2_lsx.S`, `cnrm2_lasx.S`: with a
  negative stride, returned 0. They scale by `?amax_k` of the same vector,
  which returns 0 for a stride below 1; they now walk x forwards from its
  lowest element.

One more patch works around Zig, not OpenBLAS: `power/zdot.S` and
`zdot_cell.S` return the complex result in f1/f2 as the ELFv2 ABI says, but
Zig reads a two-float struct from r3/r4. The sources already copy it into
r3/r4 for the Fortran interface; the patch turns that on. It can go once Zig
follows the ABI. Likewise arm's `cdot_vfp.S` and `zdot_vfp.S` for the
soft-float `sf` variant: they return through a pointer in r0, as AAPCS
says, but Zig reads the struct from r0 to r3, so the patch drops the pointer
and returns there.

And the sparc kernels' `64` variant shifts with `sllx` where the sources say
`sll`: they scale 64-bit strides with `sll`, which under qemu keeps only the
low 32 bits, so negative strides broke. `issues.md` has the details, and the
open question of whether real hardware behaves the same.

The ported `utest/test_swap.c` found one more, in `sparc/swap.S` and
`zswap.S`: with a zero stride, the strided loop loads a whole block before
storing it, so `n` swaps of one element acted as one. A zero stride now
takes the one-at-a-time loop. `power/rot.S`, `zrot.S`, `swap.S`, `zswap.S`,
`rot_ppc440.S` and `zrot_ppc440.S` had the same bug (`utest`'s
`zdrot_inc_0` on cell) and take the same fix, as do loongarch64 `swap.S`,
`rot_lsx.S`, `rot_lasx.S` and mips64 `swap.S`. Upstream's own zero-stride
cases in x86 `swap.S`, `zswap.S` and arm `swap_vfp.S`, `rot_vfp.S` skip the
work or stop after one element, right only for some `n`; zero strides there
now take one element at a time too.

## Kernels that call another routine's kernel

A handful of upstream sources branch into another routine's own kernel rather
than being self-contained (`loongarch64`'s `cnrm2`/`snrm2`, `lsx`/`lasx`: `bl
camax_k`/`bl samax_k`). Preprocessed the way OpenBLAS itself would for a
`DYNAMIC_ARCH` build (`gen.zig` always passes `-DDYNAMIC_ARCH`, since that is
what this project is doing), the call resolves to one core's own kernel (`bl
camax_k_LA264`) rather than the generic dispatch trampoline, which
`nakedLines()` rewrites to this project's own name for that symbol
(`blaskern_camax_la264`) -- already exported by the same isolated per-core
kernel object, so no cross-routine import is needed.

## Building the kernel objects

The per-core kernel objects are built with LLVM. Every architecture needs it
today: the self-hosted x86_64 backend's assembler rejects the register-based
memory operands (`(%rdx)`) and the VEX/EVEX mnemonics of the converted
assembly, the self-hosted wasm backend emits invalid wasm for the SIMD kernel,
and aarch64 and riscv64 have no self-hosted backend that covers them.
`-Dblas-llvm=true|false` overrides the choice; leaving it unset uses LLVM only
where `build/blaskern.zig` (`needsLlvm`) says the kernels need it.

`-Dblas-cores` picks which objects are built: `target`, the debug default,
builds only the core the target maps to (the host's own, detected, for a
native target); `all`, the release default, every core the target can run;
`none`, no objects, leaving the portable kernels. A consumer passes its
optimize mode to rstd so that the default follows it.

## Core selection

`blas.init` calls `cpu.choose` (`cpu/_select.zig`), which is `gotoblas_dynamic_init`
in two steps.

1. **Which core.** `OPENBLAS_CORETYPE` if it names a core (case-insensitive;
   `unknown` is no name) and the CPU can run it, else the detector's answer. A
   name no core has is reported at `OPENBLAS_VERBOSE` 1 and detection decides
   (the reference falls back to its generic core instead). A core the CPU cannot
   run is replaced by the closest one it can.
2. **Which object.** A core with no kernel object in the build (`-Dblas-cores`,
   or a core with no assembly) is served by the closest linked core that needs no
   feature it lacks, by instruction set (`select.served`), else the build's pick,
   else the portable kernels. The reference's `#define gotoblas_X gotoblas_Y`
   chains follow no such rule (Barcelona to Nehalem puts SSSE3 code on a CPU
   without it). The served core's parameters are used too.

Each `cpu/_detect/<arch>.zig` is a thin `probe` (cpuid, auxv, sysfs, sysctl) that
fills a `Caps`, and pure `decide`, `runnable` and table functions the tests feed
synthetic values. `decide` never returns a core `runnable` rejects.

Where the port differs from `dynamic*.c`:

- x86: the OS AVX, AVX2 and AVX-512 state checks and the reference's fallbacks
  and warnings are kept; AMD family 10h and 12h are Barcelona (the reference
  returns NULL), FMA3 is required for the Haswell and Zen kernels, `cpuid` leaves
  are range-checked, and macOS's lazy ZMM state is not mistaken for no AVX-512.
- arm: the reference has no run-time dispatch; `cpuid_arm.c`'s choice is made at
  run time from the hwcap bits, requiring 32 VFP registers for the ARMv7 kernels.
- arm64: the part table follows `dynamic_arm64.c`; Windows reads `MIDR_EL1` from
  the registry and SVE, SVE2, SME and BF16 from `IsProcessorFeaturePresent`,
  FreeBSD uses the auxiliary vector as Linux does. Without the `MIDR_EL1`
  emulation the per-CPU sysfs pass reads `cpu_capacity` (decimal), as
  `cpuid_arm64.c` does; `dynamic_arm64.c` reads `cpu_capability`, which Linux
  does not have.
- riscv: `vlenb` picks ZVL256B (no object yet, so ZVL128B serves it).
- power: POWER10 needs ISA 3.1 and MMA. POWER3 to 5, 970, Cell, G4 and 440
  (static targets upstream) are told by `AT_PLATFORM`, the PVR or the hwcap ISA
  bits, so none runs a later generation's code; a CPU nothing names gets the
  portable kernels (the reference guesses POWER8). AIX is not a Zig target.
- mips: the reference tests `HWCAP_LOONGSON_CPUCFG` the wrong way round and its
  `cpucfg` word has both registers zero; both are corrected here, unverified on
  hardware.
- Static-only `TARGET`s (P5, P6, VIAC3, Yonah; Cortex-A76, A510, A710, X1,
  X2, Ampere1, FT2000) are cores that `OPENBLAS_CORETYPE` can name, with their
  objects; detection, like the reference's, never picks them. Cortex-A9 and A15
  are detected by part, as `cpuid_arm.c` does.

## Qemu rigs

`zig build qemu-test-<machine>` runs the kernel tests under `qemu-system-*`
(see `build/qemutest.zig` for the machines). `zig build rstd-test` runs them
natively and `zig build wasm-test` under node.

`qemu-test-x86_64-whpx` runs them on the host's own CPU under WHPX (an
x86_64 Windows host only), with QEMU's `SapphireRapids` model plus FMA4 and
XOP, of which WHPX keeps what the host has. TCG implements no AVX-512, FMA4
or XOP, so the skylakex, cooperlake and sapphirerapids kernels, and the
bulldozer family's, run only there, on a host that has them.

The big-endian rigs are `s390x`, `powerpc64` (`pseries` under VOF, a
POWER10, so every 64-bit POWER core runs; Zig emits ELFv2 there too, so it
boots as the little-endian rigs do), and `mips64r2`, `mips64r6`, `mips32r2`
and `mips32r5`, the `malta` rigs' big-endian twins. 32-bit MIPS has no
converted kernels, so its rigs test the portable ones.

`powerpc` is 32-bit PowerPC, a G4 under OpenBIOS on `mac99`, for the ppcg4
core. It builds in debug, which the 64-bit POWER rigs cannot (LLVM's
FastISel bug in issues.md needs VSX). `powerpc-pseries` runs the same build
in a POWER10's 32-bit mode under VOF, so power3 to power6, ppc970 and cell
run as 32-bit code.

`arm-softfp` builds for `eabi`, the soft-float convention, so the kernels'
`sf` variant runs, using the VFP the CPU has.

`sparc64-debug` is `sparc64` in debug, in ten images where `sparc64` needs
seven (OpenBIOS maps 4 MB of each), less one test that alone is 6 MB in debug.
`sparc` (32-bit) is a stub: lld cannot link 32-bit SPARC.

The kernel rigs link every core the target CPU can run, whatever the build
mode, and `-Dblas-cores` does not change that (nor for the Linux rigs). Every
family has a rig that builds in debug.

Mach-O needs macOS itself: `qemu-test-x86_64-macos` runs the unit tests in an
x86_64 macOS VM (set up once; see `build/qemutest.zig`), under emulation,
since WHPX on an AMD host shows macOS the AMD CPU, which it does not boot on.
If the VM does not answer SSH, the step starts it headless (log in
`.macos-vm/x86_64-macos/qemu.log`, VNC on 127.0.0.1:5901), prints its serial
console (`serial.log`) while it boots, and shuts it down
after the tests, pass or fail, unless `-Dmac-keep`. `qemu-kill-x86_64-macos`
kills a VM a build step started.
`zig build rstd-check-macho` builds the tests for x86_64 and aarch64 macOS
without running them, which is where Mach-O assembly errors show.

`qemu-test-linux-<arch>` (aarch64, x86_64, arm, riscv64, powerpc64le) boots
a Debian netboot kernel with the hosted test executable as PID 1
(`test/linux/init.zig`), for the tests that need an OS: auxiliary-vector
detection, threads, files.

Apple silicon runs on a real Mac: `zig build remote-test-aarch64-macos`
copies the tests to the SSH alias `-Dremote-mac` and runs them there
(`build/remotetest.zig`); the Mac's address and login stay in the SSH config.
On an M6 this covers `vortex`, `vortexm4` and SME without SVE, which QEMU's
`-cpu max` cannot show. `zig build remote-bench-aarch64-macos
-Dopenblas=<checkout>` races the kernels there (`tools/blasbench/`); the
results are in `tools/blasbench/results-m6.md`.

## OpenBLAS's test suites

`_test_openblas/` ports OpenBLAS's own tests to the Zig interfaces: the
reference BLAS programs `test/[sdcz]blat[123].f` with their `.dat` sweeps,
the parts of their CBLAS versions (`ctest/`) that differ, and the unit tests
(`utest/`, `utest/test_extensions/`), from v0.3.34-317. The `.dat` tables run
natively; the qemu rigs take the example table in each program's header
(`N` up to 9). Every linked core the machine runs gets the suites again, on
that small table. Bad arguments map onto the named errors, and what cannot
be translated (bad letters, row-major CBLAS, forks) is listed in
`_test_openblas/__root__.zig`. Tests of routines not ported yet (level 3
beyond real `gemm`, LAPACK, the matrix copies) skip until they are; the
level-3 ones are written against the signatures they expect.
