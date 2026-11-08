# BLAS flyoff: rstd, OpenBLAS and BLIS on every native CPU

On 2026-09-21 the question was whether rstd could be checked against prebuilt BLAS libraries, Linux ones included, "to check our work on crosscompiling". This races rstd's public routines against prebuilt OpenBLAS and BLIS (and Accelerate on the Mac) on every OS and CPU at hand, with no emulated instruction pipeline: bare metal, or a VM on the host's own CPU. It compares results as well as speed: rstd's largest difference from each library, rstd's own error against an f64 reference, and a hash of rstd's output, which must match on two machines that run the same kernels.

Measured with `blas_bench.zig`'s `lib1`, `lib2` and `lib3` sections: x86_64 on 2026-10-27, arm64 again late on 2026-10-29, on a quiet Mac under the [timing rule](#timing-rule) below, at `review` `28171adc`. The full tables, with each platform's machine, libraries, threads and times:

- [Windows x86_64](results-flyoff-windows-x86_64.md): bare metal, Ryzen 7 7730U; official OpenBLAS 0.3.33.
- [Linux x86_64](results-flyoff-linux-x86_64.md): the same CPU in a WHPX VM; nixpkgs' OpenBLAS 0.3.33 and BLIS 2.1.
- [macOS arm64](results-flyoff-macos-aarch64.md): bare metal, M6 Mac mini; nixpkgs' OpenBLAS 0.3.33 and Accelerate; rstd both routed to Accelerate (its default) and on its own kernels.
- [Linux arm64](results-flyoff-linux-aarch64.md): the M6 in an HVF VM; nixpkgs' OpenBLAS 0.3.33.

## Answers

- **Cross-compiled rstd computes the same bits.** Windows and Linux on the same Zen 3, with different OS, ABI, libc and loader, select the same `zen` kernels and agree on all 228 outputs, bit for bit. On the M6, macOS on rstd's own kernels (`--rstd-accelerate=0`) and Linux agree on 166 of 228; the 62 others all lie in the slots where rstd deliberately differs by OS (`kernels/arm64.zig`'s `usesVector`: on Apple's systems the `@Vector` kernels fill only the routed slots, leaving the rest for Accelerate): `dot`, `asum`, complex `scal`, `gemv` and two large complex `axpy`s. Routed to Accelerate, macOS's default, 121 agree: the rest are Accelerate's bits. Across architectures, `gemm`, one `@Vector` source for every core, gives the same bits on arm64 and x86_64 for every precision at n = 64 and 128, and for all but `z` at 256; from 512 the cores' blockings part.
- **rstd agrees with every library to rounding.** No difference passes what accumulation order explains (no `!` anywhere), no contestant failed a check, and every `i?amax`, `copy` and `swap` matches exactly. Where rstd and a library run the same converted kernel the difference is 0 in nearly every row.
- **Three library results are far less accurate than rstd's**, inside the bound but worth knowing (below): BLIS's `snrm2` and `scnrm2`, and OpenBLAS's arm64 `scnrm2`, are off by about 1% for 16M floats.
- **Speed, one thread, Zen 3** (the same on Windows and Linux): rstd matches OpenBLAS where it runs a converted kernel (`copy`, `swap`, complex `rot` 0.94 to 1.05, `axpy` 0.85 to 1.8, real `gemm` 0.97 to 1.30) and mostly beats it in `nrm2` (up to 11.7 times, against OpenBLAS's scaled loop) and `asum` (0.67 to 2.6). It trails where Zen's upstream kernels are C and still unconverted, so rstd runs the portable loop: every `gemv` (0.06 to 0.61), real `scal` and `rot` (down to 0.13). Short real `dot`s (0.30 at 1024 floats) and complex `gemm` (0.79 to 1.01) are slower too. BLIS's `gemm` is within 0.81 to 1.04 of OpenBLAS and its `gemv` within 0.33 to 1.19, but its reductions are slow: `nrm2` 0.02 to 0.36, `asum` 0.07 to 0.45.
- **Speed, one thread, M6:** on Linux rstd runs vortex's converted and `@Vector` kernels and is at or above OpenBLAS almost everywhere: `dot` 0.37 to 4.2 (`sdot` from 256K floats the low end: rstd's routine runs at about half its kernel's rate), `gemv` 0.87 to 3.2, `nrm2` 4.3 to 17, `gemm` 0.94 to 1.13. On macOS rstd runs Accelerate where it routes, so its `gemm` is 6.2 to 17 times OpenBLAS's (Accelerate's SME) and its `gemv` up to 18. `nrm2`, `sdot` and `i?amax` now stay on rstd's kernels; routing loses only in `zcopy` of 4M elements (0.83) and strided `zrot` (0.97); see [macOS arm64](results-flyoff-macos-aarch64.md#accelerate-routing). On its own kernels, without the `@Vector` ones, `asum`, `copy`, `ddot`, complex `dot` and `sgemv_n` fall back to the portable loop (down to 0.05). At the smallest n, rstd's public routines now run 7 to 30% behind the kernels they call, on both OSes.
- **Default threading (indicative):** on x86_64 (2026-10-27) rstd ran one thread, so against OpenBLAS's threads its `gemm` fell to 0.14 to 0.52 from n = 256. In the arm64 VM rstd now threads on its own `Io.Threaded` ([`threading.md`](../../../src/linalg/blas/threading.md)): its `gemm` is 0.21 to 1.15 of OpenBLAS's 6 threads, below 1 at n = 64 and 128 and in real, above from n = 256 in complex. On macOS rstd routes `gemm` to Accelerate, which rstd's init leaves single-threaded, so its 1.28 to 8.7 times OpenBLAS's 12 threads is one Accelerate thread: not a multithreaded comparison. BLIS stays at one thread unless told otherwise.

## At a glance

GB/s for `dot` and `axpy`, GFLOP/s for `gemv` and `gemm`; one thread unless marked. `ob` is OpenBLAS, `accel` Accelerate; `Mac rstd` is rstd routed to Accelerate, its default, `Mac own` rstd on its own kernels (`--rstd-accelerate=0`).

| routine | p | n | Win ob | Win rstd | Lx ob | Lx rstd | Lx blis | Mac ob | Mac rstd | Mac own | Mac accel | La ob | La rstd |
| ------- | - | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: |
| dot | d | 4096 | 134.7 | 134.3 | 124.4 | 125.8 | 117.7 | 169.4 | 504.3 | 33.9 | 507.4 | 169.8 | 210.1 |
| dot | d | 8388608 | 30.0 | 30.1 | 28.6 | 28.6 | 28.6 | 123.4 | 158.9 | 32.4 | 159.9 | 124.3 | 124.5 |
| axpy | d | 4096 | 198.6 | 196.4 | 170.5 | 171.3 | 162.5 | 255.4 | 742.1 | 245.4 | 742.1 | 249.8 | 282.3 |
| axpy | d | 8388608 | 26.8 | 26.8 | 25.1 | 25.6 | 24.9 | 145.6 | 148.4 | 151.3 | 148.5 | 147.7 | 147.0 |
| gemv_n | d | 724 | 22.4 | 6.8 | 20.4 | 6.5 | 19.5 | 19.0 | 199.7 | 19.0 | 199.1 | 19.0 | 28.9 |
| gemv_t | d | 724 | 23.9 | 2.9 | 22.3 | 2.6 | 20.4 | 20.0 | 173.7 | 20.1 | 173.3 | 20.1 | 29.7 |
| gemv_n | d | 2896 | 7.2 | 4.1 | 6.9 | 4.0 | 8.2 | 19.0 | 38.4 | 19.0 | 38.1 | 18.9 | 19.5 |
| gemm | s | 2048 | 114.5 | 113.8 | 107.4 | 105.8 | 105.5 | 127.7 | 1846.0 | 121.8 | 1839.7 | 127.6 | 127.7 |
| gemm | d | 256 | 54.0 | 57.2 | 49.6 | 52.6 | 49.1 | 61.5 | 472.7 | 59.7 | 469.2 | 62.1 | 61.9 |
| gemm | d | 2048 | 49.1 | 54.7 | 48.7 | 49.9 | 50.7 | 64.1 | 437.4 | 60.9 | 437.4 | 64.0 | 63.6 |
| gemm | z | 1024 | 56.7 | 46.2 | 53.1 | 42.9 | 51.8 | 56.6 | 459.2 | 61.0 | 458.6 | 56.7 | 63.9 |
| gemm, default threads (indicative) | d | 2048 | 170.6 | 52.4 | 163.3 | 48.9 | 45.7 | 356.3 | 457.7 | - | 457.2 | 303.1 | 260.8 |

`Win` is Windows x86_64, `Lx` Linux x86_64, `Mac` macOS arm64 and `La` Linux arm64. The two x86_64 platforms share a CPU, as do the two arm64 ones. Level 1 runs 5 to 15% slower in the VMs than on bare metal, the ratios between contestants much less.

rstd over OpenBLAS from the same rows; a ratio from 0.97 to 1.03 is within noise, too close to name a winner.

| routine | p | n | Win | Lx | Mac | Mac own | La |
| ------- | - | -: | --- | --- | --- | --- | --- |
| dot | d | 4096 | within noise | within noise | 2.98 | 0.20 | 1.24 |
| dot | d | 8388608 | within noise | within noise | 1.29 | 0.26 | within noise |
| axpy | d | 4096 | within noise | within noise | 2.91 | 0.96 | 1.13 |
| axpy | d | 8388608 | within noise | within noise | within noise | within noise | within noise |
| gemv_n | d | 724 | 0.30 | 0.32 | 10.51 | within noise | 1.52 |
| gemv_t | d | 724 | 0.12 | 0.12 | 8.68 | within noise | 1.48 |
| gemv_n | d | 2896 | 0.57 | 0.58 | 2.02 | within noise | 1.03 |
| gemm | s | 2048 | within noise | within noise | 14.46 | 0.95 | within noise |
| gemm | d | 256 | 1.06 | 1.06 | 7.69 | within noise | within noise |
| gemm | d | 2048 | 1.11 | within noise | 6.82 | 0.95 | within noise |
| gemm | z | 1024 | 0.81 | 0.81 | 8.11 | 1.08 | 1.13 |
| gemm, default threads (indicative) | d | 2048 | 0.31 | 0.30 | 1.28 (one thread) | - | 0.86 |

## Results that disagree

rstd's differences from the libraries, against its own error from the f64 reference, tell whose result is off. The largest, all inside the bound:

| platform | routine | p | n | rstd-lib (eps) | rstd err (eps) | off |
| -------- | ------- | - | -: | -: | -: | --- |
| Linux x86_64 | `nrm2` vs BLIS | s | 16777216 | 87201.6 | 0.3 | BLIS |
| Linux x86_64 | `nrm2` vs BLIS | c | 8388608 | 87201.6 | 0.3 | BLIS |
| macOS, Linux arm64 | `nrm2` vs OpenBLAS | c | 8388608 | 87201.6 | 0.3 | OpenBLAS |
| macOS, Linux arm64 | `nrm2` vs OpenBLAS | c | 131072 | 220.2 | 0.2 | OpenBLAS |
| Linux x86_64 | `asum` vs BLIS | s | 16777216 | 260.0 | 4.2 | BLIS |

- **Single-precision `nrm2` in BLIS, and `scnrm2` in OpenBLAS on arm64, lose about 1%** (87,202 eps of `f32` is 1.04%) for 16M floats. OpenBLAS's `kernel/arm64/znrm2.S` (`CNRM2KERNEL` of `NEOVERSEN1`, which serves the M6 in both OSes) keeps the scaled sum of squares (`SSQ`, `SCALE`) in one `f32` register and adds the terms in turn. They average 1/3, and the sum reaches about 5.6M, where an `f32` ulp is 0.5: each addition rounds by a large part of its term. BLIS differs from rstd by the very same amount, so it runs the same sum in the same order. rstd (0.3 eps), OpenBLAS's x86_64 kernels and Accelerate are exact to rounding, as is OpenBLAS's real `snrm2` on arm64 (`nrm2.S`). The bound is 2 (n + 8) eps, 33M here, so the bench does not flag it.
- **rstd's portable `sasum` on macOS** (its own kernels: no converted kernel there, and `@Vector` left to Accelerate) was off by 256 eps for 16M floats on 2026-10-27, one running `f32` sum. Since `3a649af` the portable `f32` loop sums in `f64`, and it is 0.2 eps off, at about 6% less speed.
- **`f64` sums** (`dasum`, `dznrm2`, `dgemv_n`) differ by up to 516 eps between rstd and each library from accumulation order alone; there the reference is itself an `f64` sum, so `rstd err` says only how close rstd's order is to a sequential loop.

## Reading the tables

- **Rates and ratios:** each contestant's rate, then each over OpenBLAS's. `rstd` is rstd's public routine (`blas.dot`, ...), `rstd runs` whether its slot holds a converted kernel (`asm`) or the portable loop. On arm64 `asm` and `vec x4` are the converted and `@Vector` kernels alone, as in [`results-m6-libraries.md`](results-m6-libraries.md).
- **Timing:** a rate is the fastest of 15 batches of at least 10 ms (`gemm`: 7 of 20 ms), the contestants' batches interleaved (`blas_bench.zig`'s `l1_timing`, `l3_timing`; as [`results-m6.md`](results-m6.md#method)).
- **`rstd-<lib>`:** rstd's largest difference from that library's result, over every output float, in eps of the type relative to the size of the terms behind it (the sum of the products' magnitudes; for `gemm`, `|C0| + |A| |B|` per entry, all `n^2` of them). `!` marks one past what accumulation order explains: 0 for `copy`, `swap` and `i?amax`, 8 for the other element-wise routines, 2 (n + 8) for a sum of n terms.
- **`rstd err`:** rstd's own error from an f64 reference in the same units (`gemm`: at 256 sampled entries).
- **`rstd bits`:** the low 32 bits of a Wyhash of rstd's output, for comparing machines; the inputs come from the same seeded generator everywhere.

## Skipped

- **macOS x86_64:** only through Rosetta 2, which translates the instruction stream, or the macOS VM under WHPX, which hangs on this AMD host (`EXITBS:START`). Both are emulation for this purpose.
- **Windows arm64:** no hardware.
- **BLIS on Windows:** no official build without an account or license click-through. AMD's AOCL for Windows is an installer behind an EULA page; `amd/blis` releases carry Linux tarballs only (and none since 2.2), `flame/blis` none at all. Not built from source, as asked.
- **BLIS on arm64 (macOS and Linux):** nixpkgs' `blis` lists `x86_64-linux` as its only platform, so its binary cache has no arm64 build.

## Rerunning

Under the timing rule, one command per host: `tools/blasbench/flyoff.sh` ([`readme.md`](readme.md#rerunning-the-flyoff)). One run of each, from the rstd root, with MSYS2's `ssh` and `scp` first on `PATH` (the Mac is SSH alias `macmini-m6`: `-Dssh-mac`, `-Dtest-host`):

```sh
# Windows x86_64, here: OpenBLAS's official x64 zip, unpacked
zig build bench-native -Dopenblas-lib=<dir>/bin/libopenblas.dll -- lib1 lib2 lib3
# Linux x86_64, a WHPX VM here: nixpkgs' OpenBLAS and BLIS through the Mac's Nix
zig build qemu-bench-linux-x86_64 -Dvm-bench-arg=lib1 -Dvm-bench-arg=lib2 -Dvm-bench-arg=lib3
# macOS arm64, on the Mac: nixpkgs' OpenBLAS, and Accelerate
zig build ssh-bench-aarch64-macos -- lib1 lib2 lib3
# the same, rstd on its own kernels rather than Accelerate
zig build ssh-bench-aarch64-macos -- --rstd-accelerate=0 lib1 lib2 lib3
# Linux arm64, an HVF VM on the Mac
zig build qemu-bench-linux-aarch64 -Dvm-bench-arg=lib1 -Dvm-bench-arg=lib2 -Dvm-bench-arg=lib3
```

`--threads=default` (in the VMs `-Dvm-bench-arg=--threads=default`) lets the libraries pick their thread counts, and rstd thread up to the CPUs on the bench's `Io.Threaded`; `--threads=<n>` sets each library's count and rstd's ceiling (`OPENBLAS_NUM_THREADS`) to n, 1 by default. `--rstd-io=failing` runs rstd on `Io.failing`, serial, as the x86_64 runs above did; each table's notes name rstd's `Io`. `bench-native` and `ssh-bench-aarch64-macos` take `-Dblis-lib` as well.

- **The libraries:** `tools/nix/blas.sh`, run on the Mac by `tools/ssh_test.zig`, gets nixpkgs' `openblasCompat` (32-bit integers; plain `openblas` is ILP64 on x86_64-linux) and, for x86_64-linux, `blis` at the nixpkgs of `tools/nix/flake.lock`, substituted from the binary cache and never built (`--max-jobs 0`). For Linux it writes a `newc` cpio of their runtime closure's shared libraries, with glibc's loader where Zig's executables look for it; the build caches it once per lock. For macOS it links them under `/tmp/rstd-blas`.
- **The VMs** boot Debian's netboot kernel with an initramfs of `test/linux/init.zig`, the bench as `/test` and that cpio, as the Linux qemu rigs do; `init` passes the kernel command line's arguments after `--` on to the bench. The x86_64 VM runs on WHPX here (`-cpu EPYC-Milan`, as `host` stops in SeaBIOS), the arm64 one on HVF on the Mac through `tools/test_host.zig` with `-cpu host`; a host without HVF fails rather than falls back to TCG.
- **Quiet machines:** the x86_64 runs started with no other `zig` or `qemu-system` process on this machine and were repeated until none appeared while they ran. The arm64 runs follow the timing rule below.

## Timing rule

The user's rule of 2026-10-29 for the Mac, its indexing daemons (`mediaanalysisd`, `photoanalysisd`, `textunderstandingd`, `spotlightknowledged`, `corespotlightd`, Siri's models) disabled with SIP off (`tools/macos/configure_mac_host_arm.sh`), on macOS 27.2 (26B123):

- A run starts when no `qemu-system` process (by name) runs on the Mac, the 1-minute load is under 2, and no process passes 5% of a core in a `top -l 2 -s 5` sample. Activity Monitor and its `sysmond`, which the user keeps open, are left out of the 5% test and their share is added to the load limit. The limit is 2, not 1: idle, this macOS reads a load of 1.1 to 2.3 with its CPUs 99% idle and no thread runnable, and 1.0 to 1.1 with one QEMU core busy. `uptime` and the top 8 processes are recorded before and after each run. Runs go strictly one after another, the Linux VM never beside a macOS run.
- A watch every 20 s notes any other process over 5%; such a run is set aside and repeated, 3 times at most (12 were, for `XProtect` scans, `mds`, `WindowServer` and other daemons waking).
- Each invocation runs twice. A row counts when two runs agree within 3% on every rate; otherwise its section runs again, up to 3 more times, and a row still apart is marked noisy. A row shows the better rate of its agreeing pair. One arm64 row is noisy: Linux's threaded `cgemm` at n = 128, where rstd's threads gave 198 to 249 GFLOP/s over five runs.
- A ratio from 0.97 to 1.03 is within noise. The `--threads=default` runs are indicative only.
- The benches ran from a script on the Mac itself (2026-10-29, 21:52 to 23:56 CDT), with the binaries `zig build` makes for `ssh-bench-aarch64-macos` and `qemu-bench-linux-aarch64`.

Against the first 2026-10-29 runs (01:35 to 02:35 CDT), taken with the daemons at about three cores and a load of 2.6 to 9.1, the background load had not skewed the one-thread rates: OpenBLAS's median ratio is 1.000 on both OSes (0.97 to 1.06), the converted and `@Vector` kernels' 1.000 to 1.001. Threaded OpenBLAS gained up to 9% on macOS at n = 128 and 256. The rows that moved are code and OS changes, which the platform pages list: rstd's call overhead at the smallest n and its `sdot` (both OSes), Accelerate's `icamax` on macOS 27.2, and the routing trimmed since.
