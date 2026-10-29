# BLAS flyoff: Linux arm64

A Linux VM on the M6's own cores under HVF. Part of the
[BLAS flyoff](results-flyoff.md).

- **Machine:** the M6 Mac mini under QEMU 11.1.2 (nixpkgs, as
  `tools/nix` pins it) with `-accel hvf -M virt -cpu host -smp 6 -m 4G`:
  the host's cores run the guest, never translated. Debian 13's netboot
  kernel and an initramfs, as for Linux x86_64. The clock probe read 4.66
  to 4.78 GHz.
- **Build:** Zig 0.18.0-dev.751+069480a1a, ReleaseFast, `aarch64-linux-gnu`,
  CPU model `apple_m1` as for macOS; rstd detects `vortex`, not
  `vortexm4`: HVF does not pass SME to the guest. rstd at `review`
  `6d66239`.
- **OpenBLAS:** 0.3.33 from nixpkgs at `5adc5083`, `openblasCompat`, the
  same store path as plain `openblas`
  (`/nix/store/hrvb350dvin5id9mjh8ja07g3xqczbqq-openblas-0.3.33`):
  `OpenBLAS 0.3.33 DYNAMIC_ARCH NO_AFFINITY USE_OPENMP neoversen1
  MAX_THREADS=64`.
- **BLIS:** none: nixpkgs builds BLIS for x86_64-linux only.
- **Threads:** one each; default threading: OpenBLAS 6 (the VM's CPUs).
- **When:** 2026-10-29, 02:02 to 02:35 CDT, after the macOS runs and never
  beside them, under the
  [timing rule](results-flyoff-macos-aarch64.md#timing-rule) of the macOS
  runs (the user's, of 2026-10-29: no other `qemu-system` on the Mac,
  1-minute load under 10, two runs agreeing within 3% per row, up to 3
  reruns, else noisy). Level 1 and 2 took 3 runs, default threading 3. No
  row is noisy. Each run's conditions are under
  [run conditions](#run-conditions).

Columns as [`results-flyoff.md`](results-flyoff.md#reading-the-tables) has them.

## Level 1

GB/s; ratios are over openblas.

### Single (`s`)

| routine | p | n | inc | rstd runs | from | openblas | rstd | asm | vec x4 | rstd/ob | asm/ob | x4/ob | rstd-ob | rstd err | rstd bits |
| ------- | - | -: | --: | --------- | ---- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --- |
| axpy | s | 1024 | 1 | asm | vortex | 209.4 | 278.8 | 291.3 | 289.6 | 1.33 | 1.39 | 1.38 | 0.0 | 0.5 | caad32e5 |
| axpy | s | 8192 | 1 | asm | vortex | 211.2 | 288.2 | 289.6 | 289.5 | 1.36 | 1.37 | 1.37 | 0.0 | 0.5 | 43c33837 |
| axpy | s | 262144 | 1 | asm | vortex | 156.8 | 166.2 | 166.1 | 166.1 | 1.06 | 1.06 | 1.06 | 0.0 | 0.5 | f0ad0812 |
| axpy | s | 16777216 | 1 | asm | vortex | 145.6 | 145.6 | 145.0 | 146.4 | 1.00 | 1.00 | 1.01 | 0.0 | 0.5 | 89e4227f |
| axpy | s | 4096 | 2 | asm | vortex | 27.7 | 54.3 | 53.9 | 54.2 | 1.96 | 1.95 | 1.96 | 0.0 | 0.5 | 7d653b35 |
| dot | s | 1024 | 1 | asm | vortex | 158.4 | 95.9 | 96.0 | 95.9 | 0.61 | 0.61 | 0.61 | 0.0 | 0.0 | cb8a2643 |
| dot | s | 8192 | 1 | asm | vortex | 171.7 | 97.7 | 97.7 | 97.7 | 0.57 | 0.57 | 0.57 | 0.0 | 0.0 | a65eae7d |
| dot | s | 262144 | 1 | asm | vortex | 142.8 | 97.7 | 97.8 | 97.8 | 0.68 | 0.68 | 0.68 | 0.1 | 0.0 | a7fc0d7d |
| dot | s | 16777216 | 1 | asm | vortex | 124.0 | 97.6 | 97.5 | 97.6 | 0.79 | 0.79 | 0.79 | 0.0 | 0.0 | 308b8086 |
| dot | s | 4096 | 2 | asm | vortex | 11.8 | 40.2 | 40.3 | 40.5 | 3.41 | 3.42 | 3.43 | 0.1 | 0.0 | f51e3c89 |
| scal | s | 1024 | 1 | asm | vortex | 224.8 | 230.1 | 234.7 | 235.0 | 1.02 | 1.04 | 1.05 | 0.0 | 0.0 | 724ae348 |
| scal | s | 8192 | 1 | asm | vortex | 228.3 | 285.8 | 285.5 | 287.7 | 1.25 | 1.25 | 1.26 | 0.0 | 0.0 | e77a38e6 |
| scal | s | 262144 | 1 | asm | vortex | 160.6 | 158.3 | 158.2 | 158.2 | 0.99 | 0.99 | 0.99 | 0.0 | 0.0 | ea31acf9 |
| scal | s | 16777216 | 1 | asm | vortex | 138.5 | 138.4 | 139.0 | 138.7 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | 23e14c6d |
| scal | s | 4096 | 2 | asm | vortex | 35.1 | 37.2 | 37.2 | 37.3 | 1.06 | 1.06 | 1.06 | 0.0 | 0.0 | 56552b1a |
| nrm2 | s | 1024 | 1 | asm | vortex | 5.1 | 71.0 | 71.1 | 70.9 | 13.92 | 13.94 | 13.90 | 0.0 | 0.2 | 8a414b3b |
| nrm2 | s | 8192 | 1 | asm | vortex | 4.8 | 73.0 | 73.1 | 73.0 | 15.21 | 15.23 | 15.21 | 0.0 | 0.0 | 5fb36603 |
| nrm2 | s | 262144 | 1 | asm | vortex | 4.8 | 73.3 | 73.4 | 73.3 | 15.27 | 15.29 | 15.27 | 0.0 | 0.2 | 488d3ca7 |
| nrm2 | s | 16777216 | 1 | asm | vortex | 4.7 | 73.0 | 73.0 | 73.0 | 15.53 | 15.53 | 15.53 | 0.0 | 0.3 | 7b1207c6 |
| nrm2 | s | 4096 | 2 | asm | vortex | 4.8 | 20.7 | 20.8 | 20.6 | 4.31 | 4.33 | 4.29 | 0.0 | 0.0 | 207c05bc |
| asum | s | 1024 | 1 | asm | vortex | 143.2 | 126.9 | 128.1 | 135.4 | 0.89 | 0.89 | 0.95 | 0.5 | 0.6 | ba9e885c |
| asum | s | 8192 | 1 | asm | vortex | 146.1 | 142.7 | 143.0 | 144.6 | 0.98 | 0.98 | 0.99 | 0.5 | 0.4 | aed43f7e |
| asum | s | 262144 | 1 | asm | vortex | 146.0 | 139.6 | 139.4 | 140.8 | 0.96 | 0.95 | 0.96 | 2.0 | 1.5 | 3a2a771e |
| asum | s | 16777216 | 1 | asm | vortex | 108.4 | 107.5 | 107.5 | 107.4 | 0.99 | 0.99 | 0.99 | 14.0 | 7.8 | 8b557299 |
| asum | s | 4096 | 2 | asm | vortex | 8.5 | 30.0 | 30.0 | 27.6 | 3.53 | 3.53 | 3.25 | 15.6 | 0.1 | 12e3d328 |
| iamax | s | 1024 | 1 | asm | vortex | 81.0 | 76.3 | 78.3 | 83.6 | 0.94 | 0.97 | 1.03 | 0.0 | 0.0 | dc379299 |
| iamax | s | 8192 | 1 | asm | vortex | 126.6 | 91.1 | 91.3 | 92.6 | 0.72 | 0.72 | 0.73 | 0.0 | 0.0 | d2336677 |
| iamax | s | 262144 | 1 | asm | vortex | 141.2 | 95.6 | 95.6 | 96.8 | 0.68 | 0.68 | 0.69 | 0.0 | 0.0 | c22f0fea |
| iamax | s | 16777216 | 1 | asm | vortex | 105.0 | 91.0 | 91.0 | 91.9 | 0.87 | 0.87 | 0.88 | 0.0 | 0.0 | b4a99c3e |
| iamax | s | 4096 | 2 | asm | vortex | 9.3 | 38.0 | 38.1 | 37.6 | 4.09 | 4.10 | 4.04 | 0.0 | 0.0 | 54fe13f7 |
| copy | s | 1024 | 1 | asm | vortex | 132.4 | 270.1 | 276.3 | 271.5 | 2.04 | 2.09 | 2.05 | 0.0 | 0.0 | 238f9f91 |
| copy | s | 8192 | 1 | asm | vortex | 144.7 | 284.8 | 285.6 | 278.9 | 1.97 | 1.97 | 1.93 | 0.0 | 0.0 | 2984510f |
| copy | s | 262144 | 1 | asm | vortex | 140.0 | 147.3 | 149.5 | 148.1 | 1.05 | 1.07 | 1.06 | 0.0 | 0.0 | e162e98c |
| copy | s | 16777216 | 1 | asm | vortex | 120.2 | 120.9 | 120.8 | 120.9 | 1.01 | 1.00 | 1.01 | 0.0 | 0.0 | f6a3c87d |
| copy | s | 4096 | 2 | asm | vortex | 36.4 | 36.0 | 37.0 | 68.7 | 0.99 | 1.02 | 1.89 | 0.0 | 0.0 | 9d9ce9e7 |
| swap | s | 1024 | 1 | asm | vortex | 182.7 | 264.1 | 267.2 | 267.5 | 1.45 | 1.46 | 1.46 | 0.0 | 0.0 | c9422ad8 |
| swap | s | 8192 | 1 | asm | vortex | 126.4 | 223.1 | 222.2 | 222.8 | 1.77 | 1.76 | 1.76 | 0.0 | 0.0 | 43a5fd90 |
| swap | s | 262144 | 1 | asm | vortex | 85.1 | 158.9 | 159.2 | 157.8 | 1.87 | 1.87 | 1.85 | 0.0 | 0.0 | bd1677e4 |
| swap | s | 16777216 | 1 | asm | vortex | 121.9 | 146.4 | 145.2 | 146.3 | 1.20 | 1.19 | 1.20 | 0.0 | 0.0 | d00aab0b |
| swap | s | 4096 | 2 | asm | vortex | 37.8 | 37.7 | 37.8 | 37.6 | 1.00 | 1.00 | 0.99 | 0.0 | 0.0 | 67f6c48a |
| rot | s | 1024 | 1 | asm | vortex | 223.7 | 224.4 | 222.9 | 227.9 | 1.00 | 1.00 | 1.02 | 0.0 | 0.4 | be269a89 |
| rot | s | 8192 | 1 | asm | vortex | 142.6 | 143.3 | 142.6 | 194.4 | 1.00 | 1.00 | 1.36 | 0.0 | 0.5 | 24c0a740 |
| rot | s | 262144 | 1 | asm | vortex | 150.2 | 152.7 | 150.2 | 156.5 | 1.02 | 1.00 | 1.04 | 0.0 | 0.5 | ee260036 |
| rot | s | 16777216 | 1 | asm | vortex | 134.2 | 133.9 | 135.0 | 148.2 | 1.00 | 1.01 | 1.10 | 0.0 | 0.6 | e454dacb |
| rot | s | 4096 | 2 | asm | vortex | 16.6 | 16.6 | 16.6 | 36.3 | 1.00 | 1.00 | 2.19 | 0.0 | 0.5 | 91451f3d |

### Double (`d`)

| routine | p | n | inc | rstd runs | from | openblas | rstd | asm | vec x4 | rstd/ob | asm/ob | x4/ob | rstd-ob | rstd err | rstd bits |
| ------- | - | -: | --: | --------- | ---- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --- |
| axpy | d | 512 | 1 | asm | vortex | 280.1 | 278.7 | 291.2 | 290.1 | 1.00 | 1.04 | 1.04 | 0.0 | 0.0 | 14fe6b1b |
| axpy | d | 4096 | 1 | asm | vortex | 249.8 | 288.2 | 289.7 | 289.5 | 1.15 | 1.16 | 1.16 | 0.0 | 0.0 | b48fd710 |
| axpy | d | 131072 | 1 | asm | vortex | 165.5 | 165.5 | 165.6 | 165.6 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | 781a2973 |
| axpy | d | 8388608 | 1 | asm | vortex | 147.1 | 146.4 | 146.0 | 146.3 | 1.00 | 0.99 | 0.99 | 0.0 | 1.0 | 25c1d122 |
| axpy | d | 2048 | 2 | asm | vortex | 54.9 | 105.4 | 107.3 | 107.3 | 1.92 | 1.95 | 1.95 | 0.0 | 0.0 | 5c4ebaac |
| dot | d | 512 | 1 | asm | vortex | 158.0 | 212.1 | 212.3 | 212.4 | 1.34 | 1.34 | 1.34 | 0.0 | 0.1 | 573cc684 |
| dot | d | 4096 | 1 | asm | vortex | 169.4 | 216.0 | 215.1 | 215.4 | 1.28 | 1.27 | 1.27 | 0.0 | 0.0 | 3de4aa5b |
| dot | d | 131072 | 1 | asm | vortex | 142.9 | 143.9 | 144.2 | 144.1 | 1.01 | 1.01 | 1.01 | 0.0 | 0.2 | 10cd2e18 |
| dot | d | 8388608 | 1 | asm | vortex | 123.6 | 124.2 | 123.9 | 123.7 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | c114232f |
| dot | d | 2048 | 2 | asm | vortex | 24.1 | 85.8 | 86.2 | 85.2 | 3.56 | 3.58 | 3.54 | 0.2 | 0.2 | 31da4e29 |
| scal | d | 512 | 1 | asm | vortex | 246.6 | 230.1 | 235.6 | 235.3 | 0.93 | 0.96 | 0.95 | 0.0 | 0.0 | 71fb20d8 |
| scal | d | 4096 | 1 | asm | vortex | 255.4 | 284.9 | 285.2 | 287.9 | 1.12 | 1.12 | 1.13 | 0.0 | 0.0 | 33d534c5 |
| scal | d | 131072 | 1 | asm | vortex | 160.7 | 157.7 | 157.6 | 157.6 | 0.98 | 0.98 | 0.98 | 0.0 | 0.0 | 6dcd76ed |
| scal | d | 8388608 | 1 | asm | vortex | 139.6 | 139.1 | 139.2 | 139.1 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | 3297b6cc |
| scal | d | 2048 | 2 | asm | vortex | 69.9 | 71.9 | 72.0 | 72.5 | 1.03 | 1.03 | 1.04 | 0.0 | 0.0 | 63c27c6a |
| nrm2 | d | 512 | 1 | asm | vortex | 10.6 | 174.7 | 177.9 | 177.2 | 16.48 | 16.78 | 16.72 | 1.2 | 0.6 | 3e7b7ad9 |
| nrm2 | d | 4096 | 1 | asm | vortex | 9.6 | 165.6 | 166.0 | 163.2 | 17.25 | 17.29 | 17.00 | 4.3 | 0.9 | 21d8e8ce |
| nrm2 | d | 131072 | 1 | asm | vortex | 9.3 | 147.7 | 147.4 | 147.6 | 15.88 | 15.85 | 15.87 | 24.5 | 23.3 | 14746a82 |
| nrm2 | d | 8388608 | 1 | asm | vortex | 9.3 | 105.7 | 105.1 | 107.1 | 11.37 | 11.30 | 11.52 | 304.3 | 80.2 | a6e8b866 |
| nrm2 | d | 2048 | 2 | asm | vortex | 9.5 | 44.5 | 44.7 | 43.8 | 4.68 | 4.71 | 4.61 | 2.4 | 1.2 | b45479ac |
| asum | d | 512 | 1 | asm | vortex | 144.1 | 140.3 | 140.7 | 142.7 | 0.97 | 0.98 | 0.99 | 0.5 | 3.0 | 4b8ce16e |
| asum | d | 4096 | 1 | asm | vortex | 146.3 | 145.2 | 145.4 | 146.1 | 0.99 | 0.99 | 1.00 | 1.5 | 5.5 | 3a292c62 |
| asum | d | 131072 | 1 | asm | vortex | 146.4 | 139.6 | 139.6 | 140.9 | 0.95 | 0.95 | 0.96 | 3.0 | 23.5 | 3fed70f4 |
| asum | d | 8388608 | 1 | asm | vortex | 104.5 | 104.8 | 104.5 | 104.8 | 1.00 | 1.00 | 1.00 | 28.0 | 485.5 | fb441d7c |
| asum | d | 2048 | 2 | asm | vortex | 17.5 | 56.1 | 56.2 | 56.1 | 3.21 | 3.21 | 3.21 | 0.0 | 0.0 | c9a36c37 |
| iamax | d | 512 | 1 | asm | vortex | 111.5 | 92.9 | 92.9 | 94.1 | 0.83 | 0.83 | 0.84 | 0.0 | 0.0 | dc379299 |
| iamax | d | 4096 | 1 | asm | vortex | 131.5 | 94.3 | 94.4 | 95.7 | 0.72 | 0.72 | 0.73 | 0.0 | 0.0 | 6e83751e |
| iamax | d | 131072 | 1 | asm | vortex | 137.4 | 95.7 | 95.7 | 96.9 | 0.70 | 0.70 | 0.71 | 0.0 | 0.0 | c22f0fea |
| iamax | d | 8388608 | 1 | asm | vortex | 103.4 | 90.7 | 90.5 | 91.2 | 0.88 | 0.88 | 0.88 | 0.0 | 0.0 | 027eda4d |
| iamax | d | 2048 | 2 | asm | vortex | 18.8 | 71.3 | 70.9 | 68.3 | 3.79 | 3.77 | 3.63 | 0.0 | 0.0 | 1c02ab35 |
| copy | d | 512 | 1 | asm | vortex | 131.4 | 232.8 | 237.7 | 237.9 | 1.77 | 1.81 | 1.81 | 0.0 | 0.0 | 846f68e3 |
| copy | d | 4096 | 1 | asm | vortex | 143.6 | 285.7 | 287.5 | 287.4 | 1.99 | 2.00 | 2.00 | 0.0 | 0.0 | 03334116 |
| copy | d | 131072 | 1 | asm | vortex | 140.1 | 148.9 | 147.4 | 149.2 | 1.06 | 1.05 | 1.06 | 0.0 | 0.0 | 3eacc605 |
| copy | d | 8388608 | 1 | asm | vortex | 121.3 | 121.8 | 121.9 | 121.8 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | 6758600f |
| copy | d | 2048 | 2 | asm | vortex | 72.4 | 70.3 | 72.3 | 135.6 | 0.97 | 1.00 | 1.87 | 0.0 | 0.0 | 7e8cf68a |
| swap | d | 512 | 1 | asm | vortex | 36.5 | 100.0 | 93.0 | 94.3 | 2.74 | 2.55 | 2.58 | 0.0 | 0.0 | b656d6ca |
| swap | d | 4096 | 1 | asm | vortex | 126.3 | 223.5 | 222.4 | 223.3 | 1.77 | 1.76 | 1.77 | 0.0 | 0.0 | 736fd008 |
| swap | d | 131072 | 1 | asm | vortex | 86.8 | 158.2 | 158.2 | 156.9 | 1.82 | 1.82 | 1.81 | 0.0 | 0.0 | 36fa1c50 |
| swap | d | 8388608 | 1 | asm | vortex | 121.9 | 144.4 | 143.3 | 144.2 | 1.18 | 1.18 | 1.18 | 0.0 | 0.0 | fe45c39e |
| swap | d | 2048 | 2 | asm | vortex | 72.7 | 72.9 | 72.9 | 72.6 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | 08704655 |
| rot | d | 512 | 1 | asm | vortex | 59.0 | 67.9 | 60.6 | 94.3 | 1.15 | 1.03 | 1.60 | 0.0 | 0.7 | 7005d051 |
| rot | d | 4096 | 1 | asm | vortex | 171.2 | 179.0 | 170.9 | 195.7 | 1.05 | 1.00 | 1.14 | 0.0 | 0.8 | 6be51bfb |
| rot | d | 131072 | 1 | asm | vortex | 161.5 | 161.7 | 161.5 | 156.5 | 1.00 | 1.00 | 0.97 | 0.0 | 0.8 | e4fb385e |
| rot | d | 8388608 | 1 | asm | vortex | 140.6 | 140.7 | 141.1 | 144.4 | 1.00 | 1.00 | 1.03 | 0.0 | 0.8 | a426f250 |
| rot | d | 2048 | 2 | asm | vortex | 33.5 | 33.4 | 33.5 | 72.2 | 1.00 | 1.00 | 2.16 | 0.0 | 0.8 | 51abe0da |

### Single complex (`c`)

| routine | p | n | inc | rstd runs | from | openblas | rstd | asm | vec x4 | rstd/ob | asm/ob | x4/ob | rstd-ob | rstd err | rstd bits |
| ------- | - | -: | --: | --------- | ---- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --- |
| axpy | c | 512 | 1 | asm | vortex | 141.5 | 271.8 | 280.1 | 284.4 | 1.92 | 1.98 | 2.01 | 0.0 | 0.9 | f7304c9f |
| axpy | c | 4096 | 1 | asm | vortex | 145.2 | 278.7 | 279.8 | 286.2 | 1.92 | 1.93 | 1.97 | 0.0 | 0.9 | d34874c1 |
| axpy | c | 131072 | 1 | asm | vortex | 146.5 | 165.6 | 165.7 | 165.7 | 1.13 | 1.13 | 1.13 | 0.0 | 1.0 | 609be9d3 |
| axpy | c | 8388608 | 1 | asm | vortex | 140.2 | 143.8 | 144.6 | 144.9 | 1.03 | 1.03 | 1.03 | 1.0 | 1.0 | 2cd831c8 |
| axpy | c | 2048 | 2 | asm | vortex | 108.7 | 105.2 | 104.6 | 105.0 | 0.97 | 0.96 | 0.97 | 0.0 | 0.9 | bedda113 |
| dot | c | 512 | 1 | asm | vortex | 100.1 | 137.5 | 139.8 | 138.3 | 1.37 | 1.40 | 1.38 | 0.0 | 0.0 | 15bdf5a3 |
| dot | c | 4096 | 1 | asm | vortex | 105.3 | 173.0 | 173.8 | 181.7 | 1.64 | 1.65 | 1.73 | 0.0 | 0.0 | 4d0b37c7 |
| dot | c | 131072 | 1 | asm | vortex | 106.6 | 142.6 | 142.5 | 142.3 | 1.34 | 1.34 | 1.33 | 0.0 | 0.0 | d7e49038 |
| dot | c | 8388608 | 1 | asm | vortex | 106.2 | 123.5 | 123.6 | 123.3 | 1.16 | 1.16 | 1.16 | 0.0 | 0.0 | 690e6dfd |
| dot | c | 2048 | 2 | asm | vortex | 9.8 | 41.8 | 41.9 | 40.6 | 4.27 | 4.28 | 4.14 | 0.2 | 0.0 | d1f88c4d |
| scal | c | 512 | 1 | asm | vortex | 94.6 | 189.3 | 190.2 | 189.9 | 2.00 | 2.01 | 2.01 | 0.7 | 0.4 | e857fb3f |
| scal | c | 4096 | 1 | asm | vortex | 96.8 | 191.2 | 191.4 | 190.5 | 1.98 | 1.98 | 1.97 | 0.8 | 0.5 | dedcb883 |
| scal | c | 131072 | 1 | asm | vortex | 97.7 | 157.6 | 157.6 | 157.1 | 1.61 | 1.61 | 1.61 | 0.8 | 0.6 | 4b03607f |
| scal | c | 8388608 | 1 | asm | vortex | 97.3 | 137.1 | 137.2 | 137.9 | 1.41 | 1.41 | 1.42 | 0.8 | 0.6 | a8f7a3ba |
| scal | c | 2048 | 2 | asm | vortex | 41.6 | 71.3 | 71.5 | 70.8 | 1.71 | 1.72 | 1.70 | 0.8 | 0.5 | b0699f13 |
| nrm2 | c | 512 | 1 | asm | vortex | 5.0 | 70.9 | 71.0 | 70.8 | 14.18 | 14.20 | 14.16 | 1.8 | 0.2 | 8a414b3b |
| nrm2 | c | 4096 | 1 | asm | vortex | 4.8 | 73.0 | 73.1 | 73.0 | 15.21 | 15.23 | 15.21 | 12.9 | 0.0 | 5fb36603 |
| nrm2 | c | 131072 | 1 | asm | vortex | 4.7 | 73.3 | 73.3 | 73.3 | 15.60 | 15.60 | 15.60 | 220.2 | 0.2 | 488d3ca7 |
| nrm2 | c | 8388608 | 1 | asm | vortex | 4.6 | 73.0 | 73.0 | 73.0 | 15.87 | 15.87 | 15.87 | 87201.6 | 0.3 | 7b1207c6 |
| nrm2 | c | 2048 | 2 | asm | vortex | 4.8 | 21.6 | 21.7 | 21.9 | 4.50 | 4.52 | 4.56 | 0.9 | 0.4 | 7eb18e21 |
| asum | c | 512 | 1 | asm | vortex | 143.2 | 125.6 | 127.2 | 135.9 | 0.88 | 0.89 | 0.95 | 0.5 | 0.6 | ba9e885c |
| asum | c | 4096 | 1 | asm | vortex | 146.2 | 143.0 | 142.9 | 144.6 | 0.98 | 0.98 | 0.99 | 0.5 | 0.4 | aed43f7e |
| asum | c | 131072 | 1 | asm | vortex | 146.4 | 139.5 | 140.0 | 140.9 | 0.95 | 0.96 | 0.96 | 2.0 | 1.5 | 3a2a771e |
| asum | c | 8388608 | 1 | asm | vortex | 105.2 | 104.0 | 104.4 | 104.1 | 0.99 | 0.99 | 0.99 | 14.0 | 7.8 | 8b557299 |
| asum | c | 2048 | 2 | asm | vortex | 15.3 | 28.8 | 28.8 | 29.2 | 1.88 | 1.88 | 1.91 | 0.5 | 0.3 | 873eff2b |
| iamax | c | 512 | 1 | asm | vortex | 100.1 | 73.6 | 73.6 | 74.5 | 0.74 | 0.74 | 0.74 | 0.0 | 0.0 | 7c523b83 |
| iamax | c | 4096 | 1 | asm | vortex | 134.2 | 111.1 | 111.5 | 113.7 | 0.83 | 0.83 | 0.85 | 0.0 | 0.0 | a22cbbee |
| iamax | c | 131072 | 1 | asm | vortex | 141.6 | 112.6 | 112.7 | 115.4 | 0.80 | 0.80 | 0.81 | 0.0 | 0.0 | ba474698 |
| iamax | c | 8388608 | 1 | asm | vortex | 103.5 | 100.4 | 100.4 | 102.2 | 0.97 | 0.97 | 0.99 | 0.0 | 0.0 | 48ba82fe |
| iamax | c | 2048 | 2 | asm | vortex | 18.8 | 43.9 | 43.9 | 44.1 | 2.34 | 2.34 | 2.35 | 0.0 | 0.0 | 36357704 |
| copy | c | 512 | 1 | asm | vortex | 131.2 | 270.4 | 276.7 | 271.6 | 2.06 | 2.11 | 2.07 | 0.0 | 0.0 | 238f9f91 |
| copy | c | 4096 | 1 | asm | vortex | 143.5 | 285.7 | 285.9 | 286.9 | 1.99 | 1.99 | 2.00 | 0.0 | 0.0 | 2984510f |
| copy | c | 131072 | 1 | asm | vortex | 139.2 | 148.3 | 146.6 | 147.3 | 1.07 | 1.05 | 1.06 | 0.0 | 0.0 | e162e98c |
| copy | c | 8388608 | 1 | asm | vortex | 119.8 | 120.3 | 120.5 | 120.3 | 1.00 | 1.01 | 1.00 | 0.0 | 0.0 | f6a3c87d |
| copy | c | 2048 | 2 | asm | vortex | 72.3 | 70.4 | 72.3 | 136.8 | 0.97 | 1.00 | 1.89 | 0.0 | 0.0 | df37064b |
| swap | c | 512 | 1 | asm | vortex | 36.1 | 94.0 | 94.0 | 95.6 | 2.60 | 2.60 | 2.65 | 0.0 | 0.0 | c9422ad8 |
| swap | c | 4096 | 1 | asm | vortex | 93.9 | 183.9 | 186.9 | 187.1 | 1.96 | 1.99 | 1.99 | 0.0 | 0.0 | 43a5fd90 |
| swap | c | 131072 | 1 | asm | vortex | 83.5 | 157.1 | 157.7 | 156.2 | 1.88 | 1.89 | 1.87 | 0.0 | 0.0 | bd1677e4 |
| swap | c | 8388608 | 1 | asm | vortex | 133.1 | 149.1 | 148.5 | 149.5 | 1.12 | 1.12 | 1.12 | 0.0 | 0.0 | d00aab0b |
| swap | c | 2048 | 2 | asm | vortex | 51.7 | 51.4 | 51.4 | 51.4 | 0.99 | 0.99 | 0.99 | 0.0 | 0.0 | f17b1e4a |
| rot | c | 512 | 1 | asm | vortex | 60.5 | 66.1 | 61.4 | 90.2 | 1.09 | 1.01 | 1.49 | 0.0 | 0.4 | be269a89 |
| rot | c | 4096 | 1 | asm | vortex | 145.5 | 157.9 | 145.5 | 166.9 | 1.09 | 1.00 | 1.15 | 0.0 | 0.5 | 24c0a740 |
| rot | c | 131072 | 1 | asm | vortex | 156.4 | 156.8 | 156.4 | 156.0 | 1.00 | 1.00 | 1.00 | 0.0 | 0.5 | ee260036 |
| rot | c | 8388608 | 1 | asm | vortex | 149.6 | 148.8 | 149.1 | 149.1 | 0.99 | 1.00 | 1.00 | 0.0 | 0.6 | e454dacb |
| rot | c | 2048 | 2 | asm | vortex | 50.7 | 50.6 | 50.8 | 49.3 | 1.00 | 1.00 | 0.97 | 0.0 | 0.5 | 246c280f |

### Double complex (`z`)

| routine | p | n | inc | rstd runs | from | openblas | rstd | asm | vec x4 | rstd/ob | asm/ob | x4/ob | rstd-ob | rstd err | rstd bits |
| ------- | - | -: | --: | --------- | ---- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --- |
| axpy | z | 256 | 1 | asm | vortex | 141.1 | 266.8 | 281.9 | 275.1 | 1.89 | 2.00 | 1.95 | 0.0 | 1.0 | 544d4934 |
| axpy | z | 2048 | 1 | asm | vortex | 145.1 | 272.0 | 282.1 | 267.7 | 1.87 | 1.94 | 1.84 | 0.0 | 1.0 | 4a164f05 |
| axpy | z | 65536 | 1 | asm | vortex | 146.5 | 165.6 | 165.5 | 165.5 | 1.13 | 1.13 | 1.13 | 0.0 | 1.0 | c3b3756c |
| axpy | z | 4194304 | 1 | asm | vortex | 139.3 | 145.2 | 145.6 | 145.9 | 1.04 | 1.05 | 1.05 | 1.0 | 1.0 | 9d706a6e |
| axpy | z | 1024 | 2 | asm | vortex | 215.3 | 207.0 | 199.4 | 204.7 | 0.96 | 0.93 | 0.95 | 0.0 | 1.0 | 8899d47f |
| dot | z | 256 | 1 | asm | vortex | 100.6 | 153.8 | 152.9 | 161.0 | 1.53 | 1.52 | 1.60 | 0.0 | 0.1 | edb8cdc5 |
| dot | z | 2048 | 1 | asm | vortex | 105.3 | 175.4 | 175.9 | 186.3 | 1.67 | 1.67 | 1.77 | 0.1 | 0.1 | 74732dd6 |
| dot | z | 65536 | 1 | asm | vortex | 106.6 | 141.9 | 141.9 | 141.8 | 1.33 | 1.33 | 1.33 | 0.0 | 0.1 | 5d89350f |
| dot | z | 4194304 | 1 | asm | vortex | 105.9 | 122.5 | 122.3 | 122.4 | 1.16 | 1.15 | 1.16 | 0.1 | 0.1 | fc754d0d |
| dot | z | 1024 | 2 | asm | vortex | 19.7 | 84.7 | 85.1 | 82.1 | 4.30 | 4.32 | 4.17 | 0.1 | 0.2 | 9a68f762 |
| scal | z | 256 | 1 | asm | vortex | 94.6 | 184.2 | 190.2 | 187.6 | 1.95 | 2.01 | 1.98 | 0.7 | 0.7 | a86de12e |
| scal | z | 2048 | 1 | asm | vortex | 96.8 | 188.6 | 191.2 | 190.8 | 1.95 | 1.98 | 1.97 | 0.8 | 0.8 | 05e4cc32 |
| scal | z | 65536 | 1 | asm | vortex | 97.7 | 157.5 | 157.5 | 157.1 | 1.61 | 1.61 | 1.61 | 0.8 | 0.8 | b974bd94 |
| scal | z | 4194304 | 1 | asm | vortex | 97.3 | 137.3 | 137.1 | 137.4 | 1.41 | 1.41 | 1.41 | 0.8 | 0.8 | 2a1d74d7 |
| scal | z | 1024 | 2 | asm | vortex | 83.1 | 138.9 | 139.4 | 137.7 | 1.67 | 1.68 | 1.66 | 0.7 | 0.7 | 9fe85bef |
| nrm2 | z | 256 | 1 | asm | vortex | 10.5 | 170.3 | 177.8 | 176.8 | 16.22 | 16.93 | 16.84 | 1.2 | 0.6 | 3e7b7ad9 |
| nrm2 | z | 2048 | 1 | asm | vortex | 9.5 | 161.8 | 165.9 | 163.2 | 17.03 | 17.46 | 17.18 | 4.3 | 1.7 | 21d8e8ce |
| nrm2 | z | 65536 | 1 | asm | vortex | 9.2 | 146.9 | 146.5 | 146.7 | 15.97 | 15.92 | 15.95 | 24.5 | 6.7 | 14746a82 |
| nrm2 | z | 4194304 | 1 | asm | vortex | 9.2 | 106.2 | 106.0 | 106.7 | 11.54 | 11.52 | 11.60 | 304.3 | 28.2 | a6e8b866 |
| nrm2 | z | 1024 | 2 | asm | vortex | 9.7 | 50.0 | 50.2 | 50.9 | 5.15 | 5.18 | 5.25 | 1.8 | 0.0 | 69e44779 |
| asum | z | 256 | 1 | asm | vortex | 144.1 | 140.5 | 141.0 | 142.7 | 0.98 | 0.98 | 0.99 | 0.5 | 1.0 | 4b8ce16e |
| asum | z | 2048 | 1 | asm | vortex | 146.3 | 145.2 | 145.5 | 146.1 | 0.99 | 0.99 | 1.00 | 1.5 | 9.1 | 3a292c62 |
| asum | z | 65536 | 1 | asm | vortex | 145.6 | 138.9 | 138.7 | 140.3 | 0.95 | 0.95 | 0.96 | 3.0 | 1.0 | 3fed70f4 |
| asum | z | 4194304 | 1 | asm | vortex | 103.8 | 103.5 | 103.5 | 103.7 | 1.00 | 1.00 | 1.00 | 28.0 | 266.5 | fb441d7c |
| asum | z | 1024 | 2 | asm | vortex | 34.8 | 60.3 | 60.4 | 61.0 | 1.73 | 1.74 | 1.75 | 1.5 | 1.5 | e5229744 |
| iamax | z | 256 | 1 | asm | vortex | 116.2 | 89.7 | 90.0 | 86.3 | 0.77 | 0.77 | 0.74 | 0.0 | 0.0 | 56bbea09 |
| iamax | z | 2048 | 1 | asm | vortex | 133.0 | 114.6 | 114.6 | 114.9 | 0.86 | 0.86 | 0.86 | 0.0 | 0.0 | b895d41e |
| iamax | z | 65536 | 1 | asm | vortex | 133.6 | 115.3 | 115.4 | 115.7 | 0.86 | 0.86 | 0.87 | 0.0 | 0.0 | ba474698 |
| iamax | z | 4194304 | 1 | asm | vortex | 104.2 | 99.9 | 100.0 | 101.8 | 0.96 | 0.96 | 0.98 | 0.0 | 0.0 | 48ba82fe |
| iamax | z | 1024 | 2 | asm | vortex | 38.6 | 76.2 | 76.2 | 76.9 | 1.97 | 1.97 | 1.99 | 0.0 | 0.0 | 4301daa9 |
| copy | z | 256 | 1 | asm | vortex | 129.9 | 262.5 | 274.0 | 271.6 | 2.02 | 2.11 | 2.09 | 0.0 | 0.0 | 846f68e3 |
| copy | z | 2048 | 1 | asm | vortex | 143.2 | 283.8 | 283.0 | 284.3 | 1.98 | 1.98 | 1.99 | 0.0 | 0.0 | 03334116 |
| copy | z | 65536 | 1 | asm | vortex | 139.3 | 147.1 | 148.4 | 147.7 | 1.06 | 1.07 | 1.06 | 0.0 | 0.0 | 3eacc605 |
| copy | z | 4194304 | 1 | asm | vortex | 120.1 | 120.4 | 120.6 | 120.3 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | 6758600f |
| copy | z | 1024 | 2 | asm | vortex | 144.5 | 137.5 | 141.5 | 228.3 | 0.95 | 0.98 | 1.58 | 0.0 | 0.0 | ba9ea3c3 |
| swap | z | 256 | 1 | asm | vortex | 183.0 | 267.8 | 266.4 | 266.5 | 1.46 | 1.46 | 1.46 | 0.0 | 0.0 | b656d6ca |
| swap | z | 2048 | 1 | asm | vortex | 125.6 | 221.0 | 223.2 | 222.8 | 1.76 | 1.78 | 1.77 | 0.0 | 0.0 | 736fd008 |
| swap | z | 65536 | 1 | asm | vortex | 85.2 | 157.4 | 157.8 | 156.6 | 1.85 | 1.85 | 1.84 | 0.0 | 0.0 | 36fa1c50 |
| swap | z | 4194304 | 1 | asm | vortex | 122.2 | 147.6 | 146.1 | 147.4 | 1.21 | 1.20 | 1.21 | 0.0 | 0.0 | fe45c39e |
| swap | z | 1024 | 2 | asm | vortex | 137.2 | 139.2 | 133.2 | 132.8 | 1.01 | 0.97 | 0.97 | 0.0 | 0.0 | 6c5680c2 |
| rot | z | 256 | 1 | asm | vortex | 250.7 | 252.5 | 253.7 | 227.4 | 1.01 | 1.01 | 0.91 | 0.0 | 0.7 | 7005d051 |
| rot | z | 2048 | 1 | asm | vortex | 189.9 | 197.1 | 190.1 | 195.0 | 1.04 | 1.00 | 1.03 | 0.0 | 0.8 | 6be51bfb |
| rot | z | 65536 | 1 | asm | vortex | 157.4 | 157.6 | 157.5 | 156.4 | 1.00 | 1.00 | 0.99 | 0.0 | 0.8 | e4fb385e |
| rot | z | 4194304 | 1 | asm | vortex | 142.2 | 142.6 | 142.6 | 145.5 | 1.00 | 1.00 | 1.02 | 0.0 | 0.8 | a426f250 |
| rot | z | 1024 | 2 | asm | vortex | 133.0 | 136.1 | 132.6 | 133.2 | 1.02 | 1.00 | 1.00 | 0.0 | 0.8 | e18ade79 |

## Level 2

`gemv`, square, unit strides, `beta` 1: GFLOP/s.

| routine | p | n | rstd runs | from | openblas | rstd | asm | vec x4 | rstd/ob | asm/ob | x4/ob | rstd-ob | rstd err | rstd bits |
| ------- | - | -: | --------- | ---- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --- |
| gemv_n | s | 90 | asm | vortex | 66.3 | 60.2 | 61.1 | 61.4 | 0.91 | 0.92 | 0.93 | 0.6 | 6.9 | ec887ee3 |
| gemv_n | s | 1024 | asm | vortex | 51.7 | 52.1 | 52.2 | 52.1 | 1.01 | 1.01 | 1.01 | 1.7 | 15.8 | b6d4008b |
| gemv_n | s | 4096 | asm | vortex | 33.5 | 36.9 | 36.7 | 36.8 | 1.10 | 1.10 | 1.10 | 1.3 | 19.6 | 384bdb2b |
| gemv_t | s | 90 | asm | vortex | 37.4 | 41.0 | 41.5 | 41.1 | 1.10 | 1.11 | 1.10 | 0.8 | 0.4 | 4cdae9fe |
| gemv_t | s | 1024 | asm | vortex | 40.2 | 53.8 | 53.8 | 53.9 | 1.34 | 1.34 | 1.34 | 0.5 | 0.3 | 2e17dcca |
| gemv_t | s | 4096 | asm | vortex | 38.9 | 38.1 | 37.8 | 38.0 | 0.98 | 0.97 | 0.98 | 0.5 | 0.3 | f28754a7 |
| gemv_n | d | 64 | asm | vortex | 20.0 | 36.3 | 37.0 | 36.8 | 1.81 | 1.85 | 1.84 | 0.0 | 4.1 | ba428633 |
| gemv_n | d | 724 | asm | vortex | 18.9 | 29.0 | 29.0 | 29.0 | 1.53 | 1.53 | 1.53 | 0.0 | 11.9 | d500cfb8 |
| gemv_n | d | 2896 | asm | vortex | 18.9 | 19.4 | 19.3 | 19.4 | 1.03 | 1.02 | 1.03 | 0.0 | 19.6 | ee62fc60 |
| gemv_t | d | 64 | asm | vortex | 19.7 | 35.1 | 35.6 | 35.6 | 1.78 | 1.81 | 1.81 | 0.0 | 0.6 | f0b4d993 |
| gemv_t | d | 724 | asm | vortex | 20.1 | 29.7 | 29.8 | 29.7 | 1.48 | 1.48 | 1.48 | 0.6 | 0.9 | 2f3d70c0 |
| gemv_t | d | 2896 | asm | vortex | 19.6 | 19.5 | 19.4 | 19.4 | 0.99 | 0.99 | 0.99 | 0.4 | 1.3 | 2aa63d1b |
| gemv_n | c | 64 | asm | vortex | 34.6 | 63.1 | 63.6 | 65.5 | 1.82 | 1.84 | 1.89 | 1.3 | 4.6 | d24de882 |
| gemv_n | c | 724 | asm | vortex | 35.8 | 80.8 | 80.7 | 82.8 | 2.26 | 2.25 | 2.31 | 2.5 | 7.7 | dd8cf00f |
| gemv_n | c | 2896 | asm | vortex | 36.3 | 59.8 | 59.9 | 60.5 | 1.65 | 1.65 | 1.67 | 1.2 | 9.0 | ee261e78 |
| gemv_t | c | 64 | asm | vortex | 45.5 | 50.0 | 50.4 | 49.9 | 1.10 | 1.11 | 1.10 | 0.8 | 0.6 | b6e86b45 |
| gemv_t | c | 724 | asm | vortex | 25.7 | 82.5 | 82.9 | 83.9 | 3.21 | 3.23 | 3.26 | 0.5 | 0.3 | b51f3824 |
| gemv_t | c | 2896 | asm | vortex | 22.6 | 66.1 | 66.1 | 65.7 | 2.92 | 2.92 | 2.91 | 0.2 | 0.2 | 33a07a90 |
| gemv_n | z | 45 | asm | vortex | 23.2 | 33.7 | 34.0 | 35.0 | 1.45 | 1.47 | 1.51 | 0.9 | 5.7 | 4721cbb7 |
| gemv_n | z | 512 | asm | vortex | 24.3 | 38.8 | 39.5 | 40.4 | 1.60 | 1.63 | 1.66 | 2.5 | 8.4 | b1d1864c |
| gemv_n | z | 2048 | asm | vortex | 24.2 | 33.5 | 33.6 | 34.2 | 1.38 | 1.39 | 1.41 | 1.9 | 9.3 | 748372da |
| gemv_t | z | 45 | asm | vortex | 27.5 | 25.2 | 25.4 | 25.3 | 0.92 | 0.92 | 0.92 | 0.9 | 0.9 | 96a49e28 |
| gemv_t | z | 512 | asm | vortex | 22.9 | 40.8 | 41.9 | 42.3 | 1.78 | 1.83 | 1.85 | 0.5 | 0.6 | cdda4ebe |
| gemv_t | z | 2048 | asm | vortex | 21.4 | 37.2 | 37.1 | 37.4 | 1.74 | 1.73 | 1.75 | 0.2 | 0.3 | ee2ce6a7 |

## Level 3

`gemm` NN, square, `C += A*B`, one thread: GFLOP/s.

| p | n | openblas | rstd | rstd/ob | rstd-ob | rstd err | rstd bits |
| - | -: | -: | -: | -: | -: | -: | --- |
| s | 64 | 121.0 | 112.9 | 0.93 | 0.0 | 0.6 | 4c2c8b1f |
| s | 128 | 123.3 | 118.8 | 0.96 | 0.0 | 0.6 | ab770b69 |
| s | 256 | 124.6 | 123.5 | 0.99 | 0.0 | 1.0 | 21fdf633 |
| s | 512 | 125.9 | 126.2 | 1.00 | 0.0 | 0.9 | 45469416 |
| s | 1024 | 127.6 | 127.5 | 1.00 | 2.3 | 0.8 | eb97f50f |
| s | 2048 | 127.6 | 127.3 | 1.00 | 1.5 | 0.4 | fc6facb0 |
| d | 64 | 59.9 | 56.7 | 0.95 | 0.0 | 1.1 | 3fe90444 |
| d | 128 | 60.6 | 59.1 | 0.98 | 0.0 | 0.7 | 92dfae0a |
| d | 256 | 62.1 | 62.0 | 1.00 | 0.0 | 0.8 | 0040fb73 |
| d | 512 | 63.5 | 63.0 | 0.99 | 2.1 | 0.6 | a2b5517a |
| d | 1024 | 64.0 | 63.5 | 0.99 | 1.5 | 0.9 | 6898ad11 |
| d | 2048 | 64.0 | 63.5 | 0.99 | 0.8 | 1.3 | 34bb3f10 |
| c | 64 | 109.3 | 115.4 | 1.06 | 0.8 | 0.3 | a7b61b84 |
| c | 128 | 111.6 | 122.3 | 1.10 | 0.9 | 0.3 | 670b5a37 |
| c | 256 | 113.1 | 125.9 | 1.11 | 0.6 | 0.4 | d7a63c6e |
| c | 512 | 114.5 | 128.4 | 1.12 | 1.0 | 0.3 | 92ad3d6b |
| c | 1024 | 114.9 | 128.4 | 1.12 | 0.5 | 0.2 | 7763b208 |
| c | 2048 | 115.0 | 128.7 | 1.12 | 0.3 | 0.1 | 71fbc007 |
| z | 64 | 54.2 | 57.6 | 1.06 | 1.1 | 0.5 | f96643da |
| z | 128 | 54.9 | 61.2 | 1.11 | 0.7 | 0.4 | 250fd4d6 |
| z | 256 | 56.1 | 62.7 | 1.12 | 0.7 | 0.5 | 5e75ac7b |
| z | 512 | 56.7 | 64.1 | 1.13 | 0.7 | 0.6 | d6fe629e |
| z | 1024 | 56.7 | 63.9 | 1.13 | 0.5 | 0.5 | 25a060ab |
| z | 2048 | 56.8 | 64.1 | 1.13 | 0.3 | 0.4 | d4c7d3e9 |

### Default threading

OpenBLAS at its own thread count; rstd one thread. Indicative only: the
most load-sensitive runs.

| p | n | openblas | rstd | rstd/ob | rstd-ob | rstd err | rstd bits |
| - | -: | -: | -: | -: | -: | -: | --- |
| s | 64 | 121.2 | 113.1 | 0.93 | 0.0 | 0.6 | 4c2c8b1f |
| s | 128 | 292.5 | 118.7 | 0.41 | 0.0 | 0.6 | ab770b69 |
| s | 256 | 443.1 | 123.4 | 0.28 | 0.0 | 1.0 | 21fdf633 |
| s | 512 | 497.3 | 126.2 | 0.25 | 0.0 | 0.9 | 45469416 |
| s | 1024 | 561.8 | 127.4 | 0.23 | 2.3 | 0.8 | eb97f50f |
| s | 2048 | 585.4 | 127.0 | 0.22 | 1.5 | 0.4 | fc6facb0 |
| d | 64 | 59.7 | 56.4 | 0.94 | 0.0 | 1.1 | 3fe90444 |
| d | 128 | 194.1 | 59.1 | 0.30 | 0.0 | 0.7 | 92dfae0a |
| d | 256 | 252.3 | 62.0 | 0.25 | 0.0 | 0.8 | 0040fb73 |
| d | 512 | 270.8 | 62.8 | 0.23 | 2.1 | 0.6 | a2b5517a |
| d | 1024 | 289.5 | 63.7 | 0.22 | 1.5 | 0.9 | 6898ad11 |
| d | 2048 | 303.5 | 63.4 | 0.21 | 0.8 | 1.3 | 34bb3f10 |
| c | 64 | 276.3 | 115.1 | 0.42 | 0.8 | 0.3 | a7b61b84 |
| c | 128 | 341.1 | 122.3 | 0.36 | 0.9 | 0.3 | 670b5a37 |
| c | 256 | 427.4 | 126.2 | 0.30 | 0.6 | 0.4 | d7a63c6e |
| c | 512 | 449.3 | 128.3 | 0.29 | 1.0 | 0.3 | 92ad3d6b |
| c | 1024 | 475.4 | 128.4 | 0.27 | 0.5 | 0.2 | 7763b208 |
| c | 2048 | 498.0 | 128.6 | 0.26 | 0.3 | 0.1 | 71fbc007 |
| z | 64 | 171.9 | 57.5 | 0.33 | 1.1 | 0.5 | f96643da |
| z | 128 | 195.5 | 61.1 | 0.31 | 0.7 | 0.4 | 250fd4d6 |
| z | 256 | 225.3 | 62.7 | 0.28 | 0.7 | 0.5 | 5e75ac7b |
| z | 512 | 231.8 | 63.3 | 0.27 | 0.7 | 0.6 | d6fe629e |
| z | 1024 | 241.5 | 63.9 | 0.26 | 0.5 | 0.5 | 25a060ab |
| z | 2048 | 250.1 | 64.1 | 0.26 | 0.3 | 0.4 | d4c7d3e9 |

## Run conditions

The Mac's 1-minute, 5-minute and 15-minute loads, and `top`'s 8 busiest
processes in % of a core, before and after each run (the VM itself is
`qemu-system`, so it never shows before or after). `lx` is the bench at
one thread, `lx-thr` with `--threads=default`; `lx-3` reran levels 1 and
2, `lx-thr-3` level 3.

- **lx-1**, 2026-10-29 02:02:00 CDT to 02:06:57: load 7.39 / 6.38 / 5.61 before, 6.47 / 5.90 / 5.56 after. `taskpolicy -b` on `corespotlightd`, `mediaanalysisd`, `spotlightknowledged.updater`, `textunderstandingd`. Before: `mediaanalysisd` 216.1, `textunderstandingd` 84.7, `modelcatalogd` 15.8, `corespotlightd` 13.8, `modelmanagerd` 6.8, `spotlightknowledged` 3.5, `loginwindow` 2.2, `dasd` 2.1. After: `mediaanalysisd` 214.7, `textunderstandingd` 85.7, `corespotlightd` 18.1, `modelcatalogd` 15.1, `modelmanagerd` 6.5, `loginwindow` 3.8, `spotlightknowledged` 1.7, `hybridsearchd` 1.3.
- **lx-2**, 2026-10-29 02:07:15 CDT to 02:12:02: load 5.85 / 5.79 / 5.53 before, 2.82 / 4.31 / 4.95 after. `taskpolicy -b` on `corespotlightd`, `mediaanalysisd`, `textunderstandingd`. Before: `mediaanalysisd` 215.9, `textunderstandingd` 87.5, `corespotlightd` 19.5, `modelcatalogd` 12.1, `modelmanagerd` 5.3, `spotlightknowledged` 2.6, `SetStoreUpdateService` 2.3, `hybridsearchd` 1.3. After: `textunderstandingd` 91.5, `photolibraryd` 31.9, `corespotlightd` 30.1, `modelcatalogd` 12.1, `SetStoreUpdateService` 5.5, `modelmanagerd` 5.0, `spotlightknowledged` 2.4, `kernel_task` 1.7.
- **lx-3**, 2026-10-29 02:29:56 CDT to 02:33:43: load 9.13 / 6.83 / 5.39 before, 5.31 / 5.82 / 5.27 after. `taskpolicy -b` on `corespotlightd`, `mediaanalysisd`, `textunderstandingd`. Before: `mediaanalysisd` 215.5, `textunderstandingd` 89.5, `corespotlightd` 20.9, `modelcatalogd` 9.7, `loginwindow` 4.7, `modelmanagerd` 4.3, `spotlightknowledged` 1.5, `kernel_task` 1.4. After: `mediaanalysisd` 214.9, `textunderstandingd` 84.0, `corespotlightd` 20.5, `modelcatalogd` 15.1, `modelmanagerd` 6.6, `hybridsearchd` 1.8, `spotlightknowledged` 1.7, `top` 1.3.
- **lx-thr-1**, 2026-10-29 02:12:26 CDT to 02:13:21: load 2.61 / 4.16 / 4.88 before, 3.65 / 4.15 / 4.83 after. `taskpolicy -b` on `corespotlightd`, `textunderstandingd`. Before: `textunderstandingd` 85.4, `corespotlightd` 23.9, `modelcatalogd` 16.5, `photolibraryd` 15.0, `modelmanagerd` 6.8, `SetStoreUpdateService` 2.8, `spotlightknowledged` 2.7, `hybridsearchd` 1.7. After: `ANECompilerService` 135.8, `textunderstandingd` 85.4, `corespotlightd` 19.0, `modelcatalogd` 18.6, `modelmanagerd` 7.2, `mediaanalysisd` 3.6, `hybridsearchd` 2.5, `GenerativeExperiencesSafetyInferenceProvider` 1.2.
- **lx-thr-2**, 2026-10-29 02:13:34 CDT to 02:14:30: load 3.35 / 4.06 / 4.79 before, 4.41 / 4.17 / 4.78 after. `taskpolicy -b` on `corespotlightd`, `textunderstandingd`. Before: `ANECompilerService` 116.7, `textunderstandingd` 89.3, `corespotlightd` 20.3, `modelcatalogd` 16.5, `modelmanagerd` 6.2, `mediaanalysisd` 5.7, `hybridsearchd` 1.9, `kernel_task` 1.4. After: `textunderstandingd` 91.5, `mediaanalysisd` 58.0, `corespotlightd` 22.6, `modelcatalogd` 11.0, `modelmanagerd` 4.6, `photolibraryd` 3.2, `SetStoreUpdateService` 2.7, `spotlightknowledged` 2.6.
- **lx-thr-3**, 2026-10-29 02:34:08 CDT to 02:35:04: load 5.35 / 5.77 / 5.26 before, 5.51 / 5.74 / 5.29 after. `taskpolicy -b` on `corespotlightd`, `mediaanalysisd`, `textunderstandingd`. Before: `mediaanalysisd` 216.1, `textunderstandingd` 88.6, `corespotlightd` 18.0, `modelcatalogd` 12.2, `modelmanagerd` 5.3, `spotlightknowledged` 1.7, `top` 1.3, `kernel_task` 1.0. After: `mediaanalysisd` 216.1, `textunderstandingd` 88.1, `corespotlightd` 25.0, `modelcatalogd` 9.9, `spotlightknowledged` 5.8, `modelmanagerd` 4.5, `SetStoreUpdateService` 1.9, `kernel_task` 1.5.
