# BLAS flyoff: rstd, OpenBLAS and BLIS on every native CPU

On 2026-09-21 the question was whether rstd could be checked against
prebuilt BLAS libraries, Linux ones included, "to check our work on
crosscompiling". This races rstd's public routines against prebuilt
OpenBLAS and BLIS (and Accelerate on the Mac) on every OS and CPU at hand,
with no emulated instruction pipeline: bare metal, or a VM on the host's
own CPU. It compares results as well as speed: rstd's largest difference
from each library, rstd's own error against an f64 reference, and a hash
of rstd's output, which must match on two machines that run the same
kernels.

Measured 2026-10-27 with `blas_bench.zig`'s `lib1`, `lib2` and `lib3`
sections. The full tables, with each platform's machine, libraries,
threads and times:

- [Windows x86_64](results-flyoff-windows-x86_64.md): bare metal, Ryzen 7
  7730U; official OpenBLAS 0.3.33.
- [Linux x86_64](results-flyoff-linux-x86_64.md): the same CPU in a WHPX
  VM; nixpkgs' OpenBLAS 0.3.33 and BLIS 2.1.
- [macOS arm64](results-flyoff-macos-aarch64.md): bare metal, M6 Mac
  mini; nixpkgs' OpenBLAS 0.3.33 and Accelerate.
- [Linux arm64](results-flyoff-linux-aarch64.md): the M6 in an HVF VM;
  nixpkgs' OpenBLAS 0.3.33.

## Answers

- **Cross-compiled rstd computes the same bits.** Windows and Linux on the
  same Zen 3, with different OS, ABI, libc and loader, select the same
  `zen` kernels and agree on all 228 outputs, bit for bit. On the M6,
  macOS and Linux agree on 162 of 228; the 66 others all lie in the slots
  where rstd deliberately differs by OS (`kernels/arm64.zig`'s
  `usesVector`: on Apple's systems the `@Vector` kernels fill only the
  routed slots, leaving the rest for Accelerate): `dot`, `asum`, complex
  `scal`, `gemv` and two large complex `axpy`s. Across architectures,
  `gemm`, one `@Vector` source for every core, gives the same bits on
  arm64 and x86_64 for every precision at n = 64 and 128, and for all but
  `z` at 256; from 512 the cores' blockings part.
- **rstd agrees with every library to rounding.** No difference passes
  what accumulation order explains (no `!` anywhere), no contestant failed
  a check, and every `i?amax`, `copy` and `swap` matches exactly. Where
  rstd and a library run the same converted kernel the difference is 0
  in nearly every row.
- **Three library results are far less accurate than rstd's**, inside the
  bound but worth knowing (below): BLIS's `snrm2` and `scnrm2`, and
  OpenBLAS's arm64 `scnrm2`, are off by about 1% for 16M floats.
- **Speed, one thread, Zen 3** (the same on Windows and Linux): rstd
  matches OpenBLAS where it runs a converted kernel (`copy`, `swap`,
  complex `rot` 0.94 to 1.05, `axpy` 0.85 to 1.8, real `gemm` 0.97 to 1.30)
  and mostly beats it in `nrm2` (up to 11.7 times, against OpenBLAS's
  scaled loop) and `asum` (0.67 to 2.6). It trails where Zen's upstream
  kernels are C and still unconverted, so rstd runs the portable loop:
  every `gemv` (0.06 to 0.61), real `scal` and `rot` (down to 0.13). Short
  real `dot`s (0.30 at 1024 floats) and complex `gemm` (0.79 to 1.01) are
  slower too. BLIS's `gemm` is within 0.81 to 1.04 of OpenBLAS and its `gemv`
  within 0.33 to 1.19, but its reductions are slow: `nrm2` 0.02 to 0.36,
  `asum` 0.07 to 0.45.
- **Speed, one thread, M6:** on Linux rstd runs vortex's converted and
  `@Vector` kernels and is at or above OpenBLAS almost everywhere: `dot`
  0.58 to 4.4 (short real ones the low end), `gemv` 0.87 to 3.2, `nrm2` 4.4
  to 17, `gemm` 0.94 to 1.13. On macOS, without the `@Vector` kernels,
  `asum`, `copy`, `ddot`, complex `dot` and `sgemv_n` fall back to the
  portable loop (down to 0.06). Accelerate's SME `gemm` is 6.1 to 17 times
  OpenBLAS's.
- **Default threading:** rstd's routines run one thread with any `Io` so
  far ([`threading.md`](../../src/linalg/blas/threading.md) is a design),
  so against OpenBLAS's threads its `gemm` falls to 0.14 to 0.52 from
  n = 256. BLIS stays at one thread unless told otherwise.

## At a glance

GB/s for `dot` and `axpy`, GFLOP/s for `gemv` and `gemm`; one thread
unless marked. `ob` is OpenBLAS, `accel` Accelerate.

| routine | p | n | Win ob | Win rstd | Lx ob | Lx rstd | Lx blis | Mac ob | Mac rstd | Mac accel | La ob | La rstd |
| ------- | - | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: |
| dot | d | 4096 | 134.7 | 134.3 | 124.4 | 125.8 | 117.7 | 170.0 | 34.9 | 508.0 | 169.3 | 216.0 |
| dot | d | 8388608 | 30.0 | 30.1 | 28.6 | 28.6 | 28.6 | 122.9 | 32.7 | 154.7 | 119.2 | 120.0 |
| axpy | d | 4096 | 198.6 | 196.4 | 170.5 | 171.3 | 162.5 | 250.1 | 249.9 | 744.2 | 249.8 | 288.0 |
| axpy | d | 8388608 | 26.8 | 26.8 | 25.1 | 25.6 | 24.9 | 148.2 | 148.2 | 146.3 | 133.8 | 133.8 |
| gemv_n | d | 724 | 22.4 | 6.8 | 20.4 | 6.5 | 19.5 | 18.9 | 19.0 | 198.7 | 18.9 | 28.9 |
| gemv_t | d | 724 | 23.9 | 2.9 | 22.3 | 2.6 | 20.4 | 20.1 | 20.1 | 172.9 | 20.1 | 29.6 |
| gemv_n | d | 2896 | 7.2 | 4.1 | 6.9 | 4.0 | 8.2 | 18.9 | 18.9 | 36.0 | 18.6 | 17.6 |
| gemm | s | 2048 | 114.5 | 113.8 | 107.4 | 105.8 | 105.5 | 127.5 | 121.6 | 1819.1 | 127.1 | 127.0 |
| gemm | d | 256 | 54.0 | 57.2 | 49.6 | 52.6 | 49.1 | 61.0 | 59.6 | 460.6 | 62.1 | 61.9 |
| gemm | d | 2048 | 49.1 | 54.7 | 48.7 | 49.9 | 50.7 | 64.1 | 60.9 | 455.3 | 63.9 | 63.4 |
| gemm | z | 1024 | 56.7 | 46.2 | 53.1 | 42.9 | 51.8 | 56.6 | 62.0 | 470.5 | 56.6 | 63.7 |
| gemm, default threads | d | 2048 | 170.6 | 52.4 | 163.3 | 48.9 | 45.7 | 353.7 | 59.4 | 522.9 | 297.3 | 63.4 |

`Win` is Windows x86_64, `Lx` Linux x86_64, `Mac` macOS arm64 and `La`
Linux arm64. The two x86_64 platforms share a CPU, as do the two arm64
ones. Level 1 runs 5 to 15% slower in the VMs than on bare metal, the
ratios between contestants much less.

## Results that disagree

rstd's differences from the libraries, against its own error from the f64
reference, tell whose result is off. The largest, all inside the bound:

| platform | routine | p | n | rstd-lib (eps) | rstd err (eps) | off |
| -------- | ------- | - | -: | -: | -: | --- |
| Linux x86_64 | `nrm2` vs BLIS | s | 16777216 | 87201.6 | 0.3 | BLIS |
| Linux x86_64 | `nrm2` vs BLIS | c | 8388608 | 87201.6 | 0.3 | BLIS |
| macOS, Linux arm64 | `nrm2` vs OpenBLAS | c | 8388608 | 87201.6 | 0.3 | OpenBLAS |
| macOS, Linux arm64 | `nrm2` vs OpenBLAS | c | 131072 | 220.2 | 0.2 | OpenBLAS |
| Linux x86_64 | `asum` vs BLIS | s | 16777216 | 260.0 | 4.2 | BLIS |
| macOS | `asum` vs both | s | 16777216 | 262.0 | 255.9 | rstd |

- **Single-precision `nrm2` in BLIS, and `scnrm2` in OpenBLAS on arm64,
  lose about 1%** (87,202 eps of `f32` is 1.04%) for 16M floats.
  OpenBLAS's `kernel/arm64/znrm2.S` (`CNRM2KERNEL` of `NEOVERSEN1`, which
  serves the M6 in both OSes) keeps the scaled sum of squares (`SSQ`,
  `SCALE`) in one `f32` register and adds the terms in turn. They average
  1/3, and the sum reaches about 5.6M, where an `f32` ulp is 0.5: each
  addition rounds by a large part of its term. BLIS differs from rstd by
  the very same amount, so it runs the same sum in the same order. rstd
  (0.3 eps), OpenBLAS's x86_64 kernels and Accelerate are exact to
  rounding, as is OpenBLAS's real `snrm2` on arm64 (`nrm2.S`). The bound
  is 2 (n + 8) eps, 33M here, so the bench does not flag it.
- **rstd's portable `sasum` on macOS** (no converted kernel there, and
  `@Vector` left to Accelerate) is off by 256 eps for 16M floats, where
  OpenBLAS and Accelerate are near exact: one sequential `f32`
  accumulator (`blas1/asum`'s generic loop). On Linux arm64 the `@Vector` kernel takes it, and on Zen 3 the converted
  kernel (4.2 eps).
- **`f64` sums** (`dasum`, `dznrm2`, `dgemv_n`) differ by up to 516 eps
  between rstd and each library from accumulation order alone; there the
  reference is itself an `f64` sum, so `rstd err` says only how close
  rstd's order is to a sequential loop.

## Reading the tables

- **Rates and ratios:** each contestant's rate, then each over OpenBLAS's.
  `rstd` is rstd's public routine (`blas.dot`, ...), `rstd runs` whether
  its slot holds a converted kernel (`asm`) or the portable loop. On arm64
  `asm` and `vec x4` are the converted and `@Vector` kernels alone, as in
  [`results-m6-libraries.md`](results-m6-libraries.md).
- **`rstd-<lib>`:** rstd's largest difference from that library's result,
  over every output float, in eps of the type relative to the size of the
  terms behind it (the sum of the products' magnitudes; for `gemm`,
  `|C0| + |A| |B|` per entry, all `n^2` of them). `!` marks one past what
  accumulation order explains: 0 for `copy`, `swap` and `i?amax`, 8 for
  the other element-wise routines, 2 (n + 8) for a sum of n terms.
- **`rstd err`:** rstd's own error from an f64 reference in the same units
  (`gemm`: at 256 sampled entries).
- **`rstd bits`:** the low 32 bits of a Wyhash of rstd's output, for
  comparing machines; the inputs come from the same seeded generator
  everywhere.

## Skipped

- **macOS x86_64:** only through Rosetta 2, which translates the
  instruction stream, or the macOS VM under WHPX, which hangs on this AMD
  host (`EXITBS:START`). Both are emulation for this purpose.
- **Windows arm64:** no hardware.
- **BLIS on Windows:** no official build without an account or license
  click-through. AMD's AOCL for Windows is an installer behind an EULA
  page; `amd/blis` releases carry Linux tarballs only (and none since
  2.2), `flame/blis` none at all. Not built from source, as asked.
- **BLIS on arm64 (macOS and Linux):** nixpkgs' `blis` lists
  `x86_64-linux` as its only platform, so its binary cache has no arm64
  build.

## Rerunning

From the rstd root, with MSYS2's `ssh` and `scp` first on `PATH` (the Mac
is SSH alias `macmini-m6`: `-Dssh-mac`, `-Dtest-host`):

```sh
# Windows x86_64, here: OpenBLAS's official x64 zip, unpacked
zig build bench-native -Dopenblas-lib=<dir>/bin/libopenblas.dll -- lib1 lib2 lib3
# Linux x86_64, a WHPX VM here: nixpkgs' OpenBLAS and BLIS through the Mac's Nix
zig build qemu-bench-linux-x86_64 -Dvm-bench-arg=lib1 -Dvm-bench-arg=lib2 -Dvm-bench-arg=lib3
# macOS arm64, on the Mac: nixpkgs' OpenBLAS, and Accelerate
zig build ssh-bench-aarch64-macos -- lib1 lib2 lib3
# Linux arm64, an HVF VM on the Mac
zig build qemu-bench-linux-aarch64 -Dvm-bench-arg=lib1 -Dvm-bench-arg=lib2 -Dvm-bench-arg=lib3
```

`--threads=default` (in the VMs `-Dvm-bench-arg=--threads=default`) lets
the libraries pick their thread counts. `bench-native` and
`ssh-bench-aarch64-macos` take `-Dblis-lib` as well.

- **The libraries:** `tools/nix/blas.sh`, run on the Mac by
  `tools/ssh_test.zig`, gets nixpkgs' `openblasCompat` (32-bit integers;
  plain `openblas` is ILP64 on x86_64-linux) and, for x86_64-linux,
  `blis` at the nixpkgs of `tools/nix/flake.lock`, substituted from the
  binary cache and never built (`--max-jobs 0`). For Linux it writes a
  `newc` cpio of their runtime closure's shared libraries, with glibc's
  loader where Zig's executables look for it; the build caches it once
  per lock. For macOS it links them under `/tmp/rstd-blas`.
- **The VMs** boot Debian's netboot kernel with an initramfs of
  `test/linux/init.zig`, the bench as `/test` and that cpio, as the Linux
  qemu rigs do; `init` passes the kernel command line's arguments after
  `--` on to the bench. The x86_64 VM runs on WHPX here (`-cpu
  EPYC-Milan`, as `host` stops in SeaBIOS), the arm64 one on HVF on the
  Mac through `tools/test_host.zig` with `-cpu host`; a host without HVF
  fails rather than falls back to TCG.
- **Quiet machines:** each timed run started with no other `zig` or
  `qemu-system` process on its machine and was repeated until none
  appeared while it ran. On the Mac the 1-minute load stayed above 1 all
  night from two spinning system daemons; see
  [macOS arm64](results-flyoff-macos-aarch64.md).
