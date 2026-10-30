# blasbench

Races rstd's BLAS and LAPACK against other libraries, and its arm64 kernels against OpenBLAS's
assembly. `blas_bench.zig` is the bench (its header lists the sections and options), `bench.zig`
the timing loop, `merge.zig` the merge of repeated runs, `flyoff.sh` the timed reruns, `nix/` the
Eigen shim. Results: [`results-flyoff.md`](results-flyoff.md) (every OS and CPU at hand),
[`results-m6.md`](results-m6.md) and [`results-m6-libraries.md`](results-m6-libraries.md) (the M6).

## Sections

- `probe`, `l1`, `l2`, `l3`: arm64 kernels in their vtable slots' ABI (`l3` needs `-Dopenblas`).
- `lib1`, `lib2`, `lib3`: rstd's public routines against the libraries through CBLAS (level 1,
  `gemv`, `gemm`), then the table-driven races through each library's Fortran-ABI routines:
  `ger`, `symv`, `trsv`, `gbmv` in `lib2`; `gemm`, `symm`, `syrk`, `trmm`, `trsm` in `lib3`.
- `lapack`: `potrf`, `getrf`, `getrs`, `gesv`, `trtri`, `potri`, and the routines rstd does not
  have yet, which race the libraries until it does: `gesvd`, `gesdd` (T34), `geqrf`, `syevd` (T111),
  `sytrf`, `gbtrf` (T110).

Adding a routine to `lib2`, `lib3` or `lapack` is adding an entry to `ext_l2`, `ext_l3` or
`lapack_ops`: a struct naming its operands, FLOP count, Fortran-ABI call, rstd's call and a residual
(the comment above `Level` in `blas_bench.zig`). An entry written before rstd has the routine
gates rstd's call on `@hasDecl`, so it joins the race when the routine lands.

Every contestant is checked once on fresh data before timing (its residual, in eps of the type
relative to the size of its terms); rstd's output is hashed (`rstd bits`) and compared with each
library's (`rstd-<lib>`). Rates are the best of several batches, the contestants' batches taken in
turn.

## Contestants

| column | library | where from |
| ------ | ------- | ---------- |
| `openblas` | OpenBLAS, with its LAPACK | nixpkgs' `openblasCompat` (32-bit integers); on Windows the official zip |
| `rstd` | rstd's public routines | this tree |
| `blis` | BLIS (no LAPACK) | nixpkgs, x86_64 Linux only |
| `accelerate` | Apple's Accelerate | macOS' framework |
| `reference` | Reference-LAPACK 3.12.1 over its own reference BLAS | nixpkgs' `lapack-reference` |
| `eigen` | Eigen 3.4.1 through `nix/eigen_shim.cpp` | nixpkgs' `eigen`, built by `nix/eigen.nix` |
| `aocl-libflame` | AOCL-libFLAME 5.3.2 (AMD's libFLAME) | nixpkgs' `amd-libflame`, x86_64 Linux only |

The reference binaries (design.md, the anthology decision) are baselines only:

- **Reference LAPACK** runs one thread; level 3 and LAPACK rows stop at n = 1024 for it.
- **Eigen** (MPL 2.0) is never ported and no Eigen source enters rstd's tree. The shim is rstd's own:
  Eigen's routines behind the Fortran names the bench calls (`gemm`, `trmm`, `trsm`, real `symm` and
  `syrk`, `potrf`, `getrf`, `getrs`, `gesv`, `trtri`, `potri`, `gesdd`, `geqrf`, `syevd`/`heevd`),
  for the arguments the bench passes only. Nix builds it with nixpkgs' Zig for each target, so the
  Mac builds the Linux ones too, with the flyoff hosts' CPUs (`-mcpu=znver3` for x86_64 Linux,
  `apple_m1` for arm64); it runs one thread (no OpenMP).
- **AOCL-libFLAME** is substituted from the binary cache; when the cache lacks it, `blas.sh` says
  so and the column is skipped. It threads by `OMP_NUM_THREADS`, which `--threads` sets.
- **libflame** (FLAME's own) is not in nixpkgs at the pinned revision, so it is not a baseline.

`tools/nix/blas.sh <system> <rev> <bench nix dir>` fetches and builds them on the Mac's Nix, for
macOS under `/tmp/rstd-blas` and for the Linux VMs into their initramfs (`/flyoff`). Windows has
no Nix: its runs race OpenBLAS alone. A library that cannot load is skipped with a warning;
OpenBLAS must load.

**Accelerate's threads.** rstd's init makes Accelerate single-threaded for its process
(`VECLIB_MAXIMUM_THREADS=1`, threading.md), so a run with `--threads` other than 1 then starts the
bench again with `--accelerate-alone=1`: OpenBLAS and Accelerate only, rstd never initialized,
Accelerate on its own threads. Its tables follow the first process's, without rstd's columns.

## Steps

| step | what |
| ---- | ---- |
| `bench-native -Dopenblas-lib=<library>` | here; `-Dblis-lib`, `-Dreference-lib`, `-Deigen-lib`, `-Daocl-libflame-lib` add the others |
| `ssh-bench-aarch64-macos` | on the Mac, nixpkgs' libraries and the reference binaries through its Nix |
| `qemu-bench-linux-<arch>` | in a Linux VM on the host's CPU (`-Dvm-bench-arg` for the bench) |
| `bench-check` | compiles the bench for Windows x86_64, Linux and macOS on x86_64 and arm64, and tests `merge.zig` |
| `flyoff-stage-<host>` | `flyoff.sh`'s files in `zig-out/flyoff-<host>` |

Arguments after `--` go to the bench. A section runs alone when named, as in
`zig build ssh-bench-aarch64-macos -- --threads=1 lapack`.

## Rerunning the flyoff

One command per host, from the rstd root under MSYS2's bash (MSYS2's `ssh` and `scp` first on
`PATH`, `ZIG` the Zig to use, `MAC` the Mac's SSH alias, default `macmini-m6`):

```sh
# Zen 3: Windows here, then the Linux VM under WHPX (its libraries from the Mac's Nix)
tools/blasbench/flyoff.sh zen3 <dir>/bin/libopenblas.dll
# M6: macOS, then the Linux VM under HVF, all on the Mac
tools/blasbench/flyoff.sh m6
```

Each stages its files (`zig build flyoff-stage-<host>`), then runs every configuration under the
timing rule: one thread with `lib1 lib2 lib3 lapack`, and default threads with `lib3 lapack`
(indicative); on macOS also rstd on its own kernels (`--rstd-accelerate=0`). The M6's run lives on
the Mac (`~/rstd-flyoff`, detached under `caffeinate`), so the Windows host may sleep or reboot;
`flyoff.sh m6` waits for it and copies the results back, and `flyoff.sh m6 --collect` collects a
run already finished. Results land in `zig-out/flyoff-<host>/runs`: per configuration
`<cfg>.merged.md` and `<cfg>.report`, and each run's `.log`, `.before`, `.after` and `.watch`.

**The timing rule** (R4; design.md, "Timing on a host macOS keeps busy"; its details in
[`results-flyoff.md`](results-flyoff.md#timing-rule)):

- Quiet the host first: on the Mac `tools/configure_mac_host.sh` (auto-login, indexing daemons off
  with SIP off). A first window starts when the host is idle twice two minutes apart.
- A run starts when no QEMU runs, the one-minute load is under 2, and no other process passes 5%
  of a core in a 5 s sample. Activity Monitor and `sysmond` are left out of the 5% test, their
  share added to the load limit. Windows has no load average: there the processor queue length
  stands in for it, under 2, sampled with `typeperf`.
- A watch every 20 s notes any other process over 5% or another QEMU; such a run is set aside and
  repeated, 3 times at most. `uptime` and the top processes are recorded before and after each run.
- Each configuration runs twice; a row counts when the two agree within 3% on every rate
  (`blas_bench --merge`); sections with rows apart run again, up to 3 more times, and a row still
  apart is marked noisy.
- A ratio from 0.97 to 1.03 is within noise; multithreaded rows are indicative only.
