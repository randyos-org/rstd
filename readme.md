# Randy Standard Library

A Zig std lib supplement

## What's in it

- **`linalg`:** a Zig-native, column-major BLAS and LAPACK, the best open-source method for each routine (OpenBLAS's kernels, ReLAPACK, faer, Reference-LAPACK) ported into Zig, with results bit-identical for any thread count; `linalg.vector` and `linalg.matrix` for small `@Vector` math, and a MATLAB-like layer.
- **`rgl`:** fonts (sfnt, MSDF atlases, text layout and rasterizing) and shaderlib for Zig's SPIR-V shaders.
- **`fft`, `stats`, `sort`:** the FFT, statistics, and sorting and sets over `linalg`; planned or in progress ([`design.md`](design.md)).
- **Smaller pieces:** `crypto` (DER, RSA, X.509), `zip`, `msoffice` (compound files), `dll`, `io`, `log`, `debug`, `mem` and the build helpers in `buildutils`.

It builds with Zig 0.17.0 and with 0.18-dev; zfls and RandyOS check it out as a git submodule.

## Supported platforms

For what each tier promises, how each platform is tested, and the untested cores and caveats, see [`support.md`](support.md).

| Tier | Meaning | Platforms |
| --- | --- | --- |
| 1 | Required first | `x86_64-windows`; `spirv32-vulkan` (rgl's shaders) |
| 2 | Required | RandyOS: x86_64, aarch64, 32-bit arm and PowerPC 750 freestanding, x86_64 and aarch64 UEFI; zfls: x86_64 and aarch64 Linux (musl and glibc) and macOS 11+, aarch64 Windows |
| 3 | Best effort | 32-bit x86, arm Linux, riscv64, powerpc64 and powerpc64le, 32-bit POWER with hard float, MIPS and MIPS64, loongarch64, s390x, `wasm32-wasi` |
| 4 | Unsupported or blocked | sparc64, 32-bit SPARC, x86 Linux, Linux on the other tier 3 ISAs, `powerpc-macos`; everything else compiles with the portable kernels only |

## Documentation

- [`design.md`](design.md): statistics, sorting and sets, the FFT, signal processing, N-d arrays and sparse matrices.
- [`src/linalg/design.md`](src/linalg/design.md): linear algebra's design and decisions; [`src/linalg/blas/readme.md`](src/linalg/blas/readme.md), the BLAS's contract.
- [`support.md`](support.md): platforms and BLAS cores by tier.

## Issues

Work is tracked in the [issues](https://github.com/randyos-org/rstd/issues), and the open `dev` to `main` pull request lists the current milestone's in work order. Labels mark the rest: `backlog` is outside the milestone, `isa: ...` and `tier N` say which platforms an issue concerns (none: all of them), and `upstream: ...` tracks a problem in another project that rstd works around.

GitHub is the canonical repository; [Codeberg](https://codeberg.org/randyos/rstd) mirrors it.

## License

MIT ([`LICENSE`](LICENSE)), with the licenses of the ported code (OpenBLAS's BSD-3 and others) in the same file.
