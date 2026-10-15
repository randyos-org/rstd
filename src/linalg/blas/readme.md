# BLAS

Port of the OpenBLAS kernel tree; the converters are in `tools/blas_asm/`.
Which OpenBLAS source each core uses is in the generated `coverage.md` files
next to each routine.

## Layout

- `Vtable.zig`: the reference's `gotoblas_t`; `blas.kernel(.slot)` calls a slot.
- `blas{1,2,3}/`, `ext/`: one directory per routine. Its `__root__.zig` holds
  the interface (stubs until ported), and `<arch>/asm/` the converted
  assembly.
- `kernels/<arch>.zig`: the one list of an architecture's kernel modules, and
  the `install` that points the vtable at a core's kernels.
- `cpu/<arch>.zig`: detection, and `canRun` where the tests need it.
- `blas1/uk_block.zig`: a routine's C driver around a converted microkernel
  block, for routines whose upstream kernel is C (`axpy`, `dot`);
  `blas1/c_kernels.md` says how such a kernel is ported.

## Converted assembly

`tools/blas_asm/gen.zig` turns each upstream `.S` into a naked function.

- **ABI variants:** a source is preprocessed once per variant: x86_64 `sysv`
  and `win`, arm `hf` and `sf`, POWER and SPARC `32` and `64`.
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

## Portable kernels

The portable kernels are Zig algorithms with the reference's documented
semantics, not copies of its C: `dsdot` accumulates in `f64` (so does `sdot`),
`nrm2` scales against overflow and underflow, `scal`'s flag keeps NaN and
Inf, `i?amax` returns the first of equal extremes, and a NaN first element is
never replaced. Not required by the design, so not done: reproducing the
reference's accumulation order (`generic/dot.c`'s unrolled and SIMD partial
sums, `sum.c`'s), its signed-zero handling for `axpby` with a zero `alpha`,
its threading, and `zrotg.c`'s internals (`zrotg` follows LAPACK's `zlartg`;
`zrotg.c`'s scaled branches are wrong). `rotg` computes in `XDouble`, which
is where the reference uses `long double` on x86; it is an `f64` alias until
the x87 work lands, after which `rotg` widens with it. Tests compare against a
wide-precision reference with tolerances; exact comparison is kept where a Zig
kernel stands in for a converted kernel of the same core.

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
- `caxpby_lsx.S`, `caxpby_lasx.S`: with `alpha` and `beta` zero and a
  stride on `y`, zeroed only the real parts, and `lasx` single one of them twice.
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

One more patch works around Zig, not OpenBLAS: `power/zdot.S` and
`zdot_cell.S` return the complex result in f1/f2 as the ELFv2 ABI says, but
Zig reads a two-float struct from r3/r4. The sources already copy it into
r3/r4 for the Fortran interface; the patch turns that on. It can go once Zig
follows the ABI.

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

## Qemu rigs

`zig build qemu-test-<machine>` runs the kernel tests under `qemu-system-*`
(see `build/qemutest.zig` for the machines). `zig build rstd-test` runs them
natively and `zig build wasm-test` under node.

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
