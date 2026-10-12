# BLAS

Port of the OpenBLAS kernel tree; the converters are in `tools/blas_asm/`.
Which OpenBLAS source each core uses is in the generated `coverage.md` files
next to each routine.

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

## Qemu rigs

`zig build qemu-test-<machine>` runs the kernel tests under `qemu-system-*`
(see `build/qemutest.zig` for the machines). `zig build rstd-test` runs them
natively and `zig build wasm-test` under node.
