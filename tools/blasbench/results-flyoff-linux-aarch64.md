# BLAS flyoff: Linux arm64

A Linux VM on the M6's own cores under HVF. Part of the
[BLAS flyoff](results-flyoff.md).

- **Machine:** the M6 Mac mini (macOS 27.2) under QEMU 11.1.2 (nixpkgs,
  as `tools/nix` pins it) with `-accel hvf -M virt -cpu host -smp 6 -m 4G`:
  the host's cores run the guest, never translated. Debian 13's netboot
  kernel and an initramfs, as for Linux x86_64. The clock probe read 4.77
  to 4.79 GHz, 4.58 at the end of a threaded run.
- **Build:** Zig 0.18.0-dev.751+069480a1a, ReleaseFast, `aarch64-linux-gnu`,
  CPU model `apple_m1` as for macOS; rstd detects `vortex`, not
  `vortexm4`: HVF does not pass SME to the guest. rstd at `review`
  `28171adc`, on the bench's own `Io.Threaded`.
- **OpenBLAS:** 0.3.33 from nixpkgs at `5adc5083`, `openblasCompat`, the
  same store path as plain `openblas`
  (`/nix/store/hrvb350dvin5id9mjh8ja07g3xqczbqq-openblas-0.3.33`):
  `OpenBLAS 0.3.33 DYNAMIC_ARCH NO_AFFINITY USE_OPENMP neoversen1
  MAX_THREADS=64`.
- **BLIS:** none: nixpkgs builds BLIS for x86_64-linux only.
- **Threads:** `--threads=1`: one each, rstd at thread ceiling 1. Default
  threading: OpenBLAS 6 (the VM's CPUs), rstd up to the 6 CPUs.
- **When:** 2026-10-29, 23:14 to 23:56 CDT, after the macOS runs and never
  beside them, under the
  [timing rule](results-flyoff-macos-aarch64.md#timing-rule) of the macOS
  runs (no other `qemu-system` on the Mac, 1-minute load under 2, no other
  process over 5%, two runs agreeing within 3% per row, up to 3 reruns,
  else noisy). The watch set 5 runs aside (`knowledge-agent`, `airportd`,
  `cfprefsd`, `biomesyncd`, `NotificationCenter`, `AddressBookSourceSync`);
  `lx-1` kept its fourth try, where `stickersd` passed 5% once. Levels 1 and
  3 took 3 runs, default threading 5. One row is noisy. Each run's
  conditions are under [run conditions](#run-conditions).

Columns as [`results-flyoff.md`](results-flyoff.md#reading-the-tables) has them.

## Level 1

GB/s; ratios are over openblas.

### Single (`s`)

| routine | p | n | inc | rstd runs | from | openblas | rstd | asm | vec x4 | rstd/ob | asm/ob | x4/ob | rstd-ob | rstd err | rstd bits |
| ------- | - | -: | --: | --------- | ---- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --- |
| axpy | s | 1024 | 1 | asm | vortex | 209.6 | 239.0 | 290.2 | 289.0 | 1.14 | 1.38 | 1.38 | 0.0 | 0.5 | caad32e5 |
| axpy | s | 8192 | 1 | asm | vortex | 211.2 | 282.4 | 289.7 | 289.4 | 1.34 | 1.37 | 1.37 | 0.0 | 0.5 | 43c33837 |
| axpy | s | 262144 | 1 | asm | vortex | 156.4 | 165.7 | 166.0 | 166.0 | 1.06 | 1.06 | 1.06 | 0.0 | 0.5 | f0ad0812 |
| axpy | s | 16777216 | 1 | asm | vortex | 145.1 | 145.8 | 146.1 | 146.3 | 1.00 | 1.01 | 1.01 | 0.0 | 0.5 | 89e4227f |
| axpy | s | 4096 | 2 | asm | vortex | 27.9 | 27.7 | 27.8 | 55.1 | 0.99 | 1.00 | 1.97 | 0.0 | 0.5 | 7d653b35 |
| dot | s | 1024 | 1 | asm | vortex | 159.3 | 94.0 | 96.1 | 96.0 | 0.59 | 0.60 | 0.60 | 0.0 | 0.0 | cb8a2643 |
| dot | s | 8192 | 1 | asm | vortex | 172.0 | 97.2 | 97.8 | 97.7 | 0.57 | 0.57 | 0.57 | 0.0 | 0.0 | a65eae7d |
| dot | s | 262144 | 1 | asm | vortex | 143.5 | 53.0 | 97.8 | 97.7 | 0.37 | 0.68 | 0.68 | 0.1 | 0.0 | a7fc0d7d |
| dot | s | 16777216 | 1 | asm | vortex | 124.4 | 51.6 | 97.5 | 97.6 | 0.41 | 0.78 | 0.78 | 0.0 | 0.0 | 308b8086 |
| dot | s | 4096 | 2 | asm | vortex | 11.8 | 39.8 | 40.4 | 40.5 | 3.37 | 3.42 | 3.43 | 0.1 | 0.0 | f51e3c89 |
| scal | s | 1024 | 1 | asm | vortex | 224.7 | 184.2 | 235.3 | 234.8 | 0.82 | 1.05 | 1.04 | 0.0 | 0.0 | 724ae348 |
| scal | s | 8192 | 1 | asm | vortex | 228.2 | 250.6 | 283.7 | 284.7 | 1.10 | 1.24 | 1.25 | 0.0 | 0.0 | e77a38e6 |
| scal | s | 262144 | 1 | asm | vortex | 160.8 | 158.3 | 158.6 | 158.7 | 0.98 | 0.99 | 0.99 | 0.0 | 0.0 | ea31acf9 |
| scal | s | 16777216 | 1 | asm | vortex | 139.3 | 139.7 | 140.1 | 139.6 | 1.00 | 1.01 | 1.00 | 0.0 | 0.0 | 23e14c6d |
| scal | s | 4096 | 2 | asm | vortex | 35.2 | 36.8 | 37.5 | 37.6 | 1.05 | 1.07 | 1.07 | 0.0 | 0.0 | 56552b1a |
| nrm2 | s | 1024 | 1 | asm | vortex | 5.0 | 71.0 | 71.2 | 70.9 | 14.20 | 14.24 | 14.18 | 0.0 | 0.2 | 8a414b3b |
| nrm2 | s | 8192 | 1 | asm | vortex | 4.8 | 73.0 | 73.2 | 73.1 | 15.21 | 15.25 | 15.23 | 0.0 | 0.0 | 5fb36603 |
| nrm2 | s | 262144 | 1 | asm | vortex | 4.8 | 73.2 | 73.4 | 73.3 | 15.25 | 15.29 | 15.27 | 0.0 | 0.2 | 488d3ca7 |
| nrm2 | s | 16777216 | 1 | asm | vortex | 4.8 | 73.1 | 73.2 | 73.1 | 15.23 | 15.25 | 15.23 | 0.0 | 0.3 | 7b1207c6 |
| nrm2 | s | 4096 | 2 | asm | vortex | 4.8 | 20.8 | 20.8 | 20.7 | 4.33 | 4.33 | 4.31 | 0.0 | 0.0 | 207c05bc |
| asum | s | 1024 | 1 | asm | vortex | 143.2 | 118.5 | 127.6 | 135.2 | 0.83 | 0.89 | 0.94 | 0.5 | 0.6 | ba9e885c |
| asum | s | 8192 | 1 | asm | vortex | 146.0 | 141.7 | 143.0 | 144.6 | 0.97 | 0.98 | 0.99 | 0.5 | 0.4 | aed43f7e |
| asum | s | 262144 | 1 | asm | vortex | 146.0 | 138.5 | 139.5 | 140.8 | 0.95 | 0.96 | 0.96 | 0.0 | 0.5 | d371127c |
| asum | s | 16777216 | 1 | asm | vortex | 110.5 | 110.4 | 110.2 | 110.4 | 1.00 | 1.00 | 1.00 | 6.0 | 0.2 | 40b71501 |
| asum | s | 4096 | 2 | asm | vortex | 8.5 | 29.6 | 30.1 | 27.7 | 3.48 | 3.54 | 3.26 | 15.6 | 0.1 | 12e3d328 |
| iamax | s | 1024 | 1 | asm | vortex | 81.0 | 69.8 | 78.4 | 83.6 | 0.86 | 0.97 | 1.03 | 0.0 | 0.0 | dc379299 |
| iamax | s | 8192 | 1 | asm | vortex | 126.5 | 89.7 | 91.3 | 92.7 | 0.71 | 0.72 | 0.73 | 0.0 | 0.0 | d2336677 |
| iamax | s | 262144 | 1 | asm | vortex | 141.2 | 91.8 | 95.6 | 96.8 | 0.65 | 0.68 | 0.69 | 0.0 | 0.0 | c22f0fea |
| iamax | s | 16777216 | 1 | asm | vortex | 107.6 | 91.7 | 92.0 | 93.0 | 0.85 | 0.86 | 0.86 | 0.0 | 0.0 | b4a99c3e |
| iamax | s | 4096 | 2 | asm | vortex | 9.4 | 37.4 | 38.1 | 37.6 | 3.98 | 4.05 | 4.00 | 0.0 | 0.0 | 54fe13f7 |
| copy | s | 1024 | 1 | asm | vortex | 133.8 | 192.7 | 277.2 | 274.3 | 1.44 | 2.07 | 2.05 | 0.0 | 0.0 | 238f9f91 |
| copy | s | 8192 | 1 | asm | vortex | 145.5 | 243.1 | 287.6 | 284.8 | 1.67 | 1.98 | 1.96 | 0.0 | 0.0 | 2984510f |
| copy | s | 262144 | 1 | asm | vortex | 140.7 | 147.7 | 149.4 | 148.8 | 1.05 | 1.06 | 1.06 | 0.0 | 0.0 | e162e98c |
| copy | s | 16777216 | 1 | asm | vortex | 121.9 | 122.0 | 122.1 | 122.1 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | f6a3c87d |
| copy | s | 4096 | 2 | asm | vortex | 36.5 | 36.5 | 37.3 | 69.0 | 1.00 | 1.02 | 1.89 | 0.0 | 0.0 | 9d9ce9e7 |
| swap | s | 1024 | 1 | asm | vortex | 182.7 | 229.0 | 267.5 | 267.8 | 1.25 | 1.46 | 1.47 | 0.0 | 0.0 | c9422ad8 |
| swap | s | 8192 | 1 | asm | vortex | 126.4 | 217.9 | 222.7 | 223.0 | 1.72 | 1.76 | 1.76 | 0.0 | 0.0 | 43a5fd90 |
| swap | s | 262144 | 1 | asm | vortex | 85.7 | 159.0 | 159.9 | 158.8 | 1.86 | 1.87 | 1.85 | 0.0 | 0.0 | bd1677e4 |
| swap | s | 16777216 | 1 | asm | vortex | 122.5 | 144.0 | 142.9 | 143.9 | 1.18 | 1.17 | 1.17 | 0.0 | 0.0 | d00aab0b |
| swap | s | 4096 | 2 | asm | vortex | 38.0 | 37.8 | 38.1 | 37.7 | 0.99 | 1.00 | 0.99 | 0.0 | 0.0 | 67f6c48a |
| rot | s | 1024 | 1 | asm | vortex | 223.5 | 202.9 | 222.8 | 227.8 | 0.91 | 1.00 | 1.02 | 0.0 | 0.4 | be269a89 |
| rot | s | 8192 | 1 | asm | vortex | 142.9 | 141.9 | 142.9 | 194.5 | 0.99 | 1.00 | 1.36 | 0.0 | 0.5 | 24c0a740 |
| rot | s | 262144 | 1 | asm | vortex | 150.1 | 150.0 | 150.1 | 156.4 | 1.00 | 1.00 | 1.04 | 0.0 | 0.5 | ee260036 |
| rot | s | 16777216 | 1 | asm | vortex | 133.2 | 133.4 | 133.1 | 145.1 | 1.00 | 1.00 | 1.09 | 0.0 | 0.6 | e454dacb |
| rot | s | 4096 | 2 | asm | vortex | 16.7 | 16.8 | 16.8 | 35.4 | 1.01 | 1.01 | 2.12 | 0.0 | 0.5 | 91451f3d |

### Double (`d`)

| routine | p | n | inc | rstd runs | from | openblas | rstd | asm | vec x4 | rstd/ob | asm/ob | x4/ob | rstd-ob | rstd err | rstd bits |
| ------- | - | -: | --: | --------- | ---- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --- |
| axpy | d | 512 | 1 | asm | vortex | 280.1 | 238.6 | 290.0 | 289.4 | 0.85 | 1.04 | 1.03 | 0.0 | 0.0 | 14fe6b1b |
| axpy | d | 4096 | 1 | asm | vortex | 249.8 | 282.3 | 289.7 | 289.8 | 1.13 | 1.16 | 1.16 | 0.0 | 0.0 | b48fd710 |
| axpy | d | 131072 | 1 | asm | vortex | 165.4 | 165.2 | 165.5 | 165.7 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | 781a2973 |
| axpy | d | 8388608 | 1 | asm | vortex | 147.7 | 147.0 | 147.3 | 147.7 | 1.00 | 1.00 | 1.00 | 0.0 | 1.0 | 25c1d122 |
| axpy | d | 2048 | 2 | asm | vortex | 55.0 | 104.5 | 107.8 | 107.5 | 1.90 | 1.96 | 1.95 | 0.0 | 0.0 | 5c4ebaac |
| dot | d | 512 | 1 | asm | vortex | 158.1 | 186.5 | 212.2 | 212.2 | 1.18 | 1.34 | 1.34 | 0.0 | 0.1 | 573cc684 |
| dot | d | 4096 | 1 | asm | vortex | 169.8 | 210.1 | 215.4 | 215.6 | 1.24 | 1.27 | 1.27 | 0.0 | 0.0 | 3de4aa5b |
| dot | d | 131072 | 1 | asm | vortex | 142.8 | 142.5 | 143.4 | 143.4 | 1.00 | 1.00 | 1.00 | 0.1 | 0.3 | 35de1d5e |
| dot | d | 8388608 | 1 | asm | vortex | 124.3 | 124.5 | 124.4 | 124.4 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | f39398e0 |
| dot | d | 2048 | 2 | asm | vortex | 24.2 | 83.8 | 86.3 | 85.2 | 3.46 | 3.57 | 3.52 | 0.2 | 0.2 | 31da4e29 |
| scal | d | 512 | 1 | asm | vortex | 246.4 | 183.9 | 236.4 | 233.7 | 0.75 | 0.96 | 0.95 | 0.0 | 0.0 | 71fb20d8 |
| scal | d | 4096 | 1 | asm | vortex | 255.3 | 248.1 | 285.4 | 287.6 | 0.97 | 1.12 | 1.13 | 0.0 | 0.0 | 33d534c5 |
| scal | d | 131072 | 1 | asm | vortex | 160.6 | 157.5 | 157.7 | 157.8 | 0.98 | 0.98 | 0.98 | 0.0 | 0.0 | 6dcd76ed |
| scal | d | 8388608 | 1 | asm | vortex | 139.9 | 139.8 | 140.0 | 139.9 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | 3297b6cc |
| scal | d | 2048 | 2 | asm | vortex | 70.1 | 70.6 | 72.5 | 72.6 | 1.01 | 1.03 | 1.04 | 0.0 | 0.0 | 63c27c6a |
| nrm2 | d | 512 | 1 | asm | vortex | 10.8 | 174.9 | 177.8 | 176.7 | 16.19 | 16.46 | 16.36 | 1.2 | 0.6 | 3e7b7ad9 |
| nrm2 | d | 4096 | 1 | asm | vortex | 9.7 | 165.6 | 166.2 | 163.2 | 17.07 | 17.13 | 16.82 | 4.3 | 0.9 | 21d8e8ce |
| nrm2 | d | 131072 | 1 | asm | vortex | 9.5 | 147.5 | 147.5 | 147.7 | 15.53 | 15.53 | 15.55 | 24.5 | 23.3 | 14746a82 |
| nrm2 | d | 8388608 | 1 | asm | vortex | 9.4 | 110.8 | 109.8 | 110.6 | 11.79 | 11.68 | 11.77 | 304.3 | 80.2 | a6e8b866 |
| nrm2 | d | 2048 | 2 | asm | vortex | 9.7 | 44.7 | 44.8 | 44.0 | 4.61 | 4.62 | 4.54 | 2.4 | 1.2 | b45479ac |
| asum | d | 512 | 1 | asm | vortex | 144.1 | 127.5 | 140.6 | 142.7 | 0.88 | 0.98 | 0.99 | 0.5 | 3.0 | 4b8ce16e |
| asum | d | 4096 | 1 | asm | vortex | 146.3 | 142.5 | 145.5 | 146.1 | 0.97 | 0.99 | 1.00 | 1.5 | 5.5 | 3a292c62 |
| asum | d | 131072 | 1 | asm | vortex | 146.4 | 138.9 | 139.5 | 140.8 | 0.95 | 0.95 | 0.96 | 0.0 | 26.5 | a59130c7 |
| asum | d | 8388608 | 1 | asm | vortex | 110.2 | 110.0 | 109.9 | 109.9 | 1.00 | 1.00 | 1.00 | 4.5 | 518.0 | 6792d763 |
| asum | d | 2048 | 2 | asm | vortex | 17.7 | 54.8 | 56.3 | 56.2 | 3.10 | 3.18 | 3.18 | 0.0 | 0.0 | c9a36c37 |
| iamax | d | 512 | 1 | asm | vortex | 111.7 | 89.0 | 92.9 | 94.1 | 0.80 | 0.83 | 0.84 | 0.0 | 0.0 | dc379299 |
| iamax | d | 4096 | 1 | asm | vortex | 131.5 | 92.0 | 94.4 | 95.7 | 0.70 | 0.72 | 0.73 | 0.0 | 0.0 | 6e83751e |
| iamax | d | 131072 | 1 | asm | vortex | 137.4 | 92.7 | 95.7 | 96.9 | 0.67 | 0.70 | 0.71 | 0.0 | 0.0 | c22f0fea |
| iamax | d | 8388608 | 1 | asm | vortex | 106.5 | 91.8 | 91.8 | 92.8 | 0.86 | 0.86 | 0.87 | 0.0 | 0.0 | 027eda4d |
| iamax | d | 2048 | 2 | asm | vortex | 18.9 | 69.9 | 71.7 | 68.3 | 3.70 | 3.79 | 3.61 | 0.0 | 0.0 | 1c02ab35 |
| copy | d | 512 | 1 | asm | vortex | 132.5 | 185.6 | 239.2 | 239.8 | 1.40 | 1.81 | 1.81 | 0.0 | 0.0 | 846f68e3 |
| copy | d | 4096 | 1 | asm | vortex | 144.7 | 261.9 | 287.3 | 287.5 | 1.81 | 1.99 | 1.99 | 0.0 | 0.0 | 03334116 |
| copy | d | 131072 | 1 | asm | vortex | 141.2 | 148.1 | 148.4 | 148.5 | 1.05 | 1.05 | 1.05 | 0.0 | 0.0 | 3eacc605 |
| copy | d | 8388608 | 1 | asm | vortex | 122.1 | 122.5 | 122.4 | 122.4 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | 6758600f |
| copy | d | 2048 | 2 | asm | vortex | 72.8 | 70.8 | 73.4 | 135.9 | 0.97 | 1.01 | 1.87 | 0.0 | 0.0 | 7e8cf68a |
| swap | d | 512 | 1 | asm | vortex | 36.5 | 93.0 | 92.9 | 94.3 | 2.55 | 2.55 | 2.58 | 0.0 | 0.0 | b656d6ca |
| swap | d | 4096 | 1 | asm | vortex | 126.5 | 219.8 | 223.3 | 224.1 | 1.74 | 1.77 | 1.77 | 0.0 | 0.0 | 736fd008 |
| swap | d | 131072 | 1 | asm | vortex | 84.8 | 159.2 | 159.1 | 158.2 | 1.88 | 1.88 | 1.87 | 0.0 | 0.0 | 36fa1c50 |
| swap | d | 8388608 | 1 | asm | vortex | 122.4 | 146.0 | 145.0 | 145.9 | 1.19 | 1.18 | 1.19 | 0.0 | 0.0 | fe45c39e |
| swap | d | 2048 | 2 | asm | vortex | 72.7 | 72.2 | 73.2 | 72.6 | 0.99 | 1.01 | 1.00 | 0.0 | 0.0 | 08704655 |
| rot | d | 512 | 1 | asm | vortex | 59.7 | 65.1 | 61.4 | 94.0 | 1.09 | 1.03 | 1.57 | 0.0 | 0.7 | 7005d051 |
| rot | d | 4096 | 1 | asm | vortex | 171.1 | 173.6 | 170.9 | 195.4 | 1.01 | 1.00 | 1.14 | 0.0 | 0.8 | 6be51bfb |
| rot | d | 131072 | 1 | asm | vortex | 160.7 | 160.8 | 160.9 | 156.2 | 1.00 | 1.00 | 0.97 | 0.0 | 0.8 | e4fb385e |
| rot | d | 8388608 | 1 | asm | vortex | 143.4 | 142.9 | 143.1 | 146.8 | 1.00 | 1.00 | 1.02 | 0.0 | 0.8 | a426f250 |
| rot | d | 2048 | 2 | asm | vortex | 33.6 | 33.4 | 33.7 | 72.5 | 0.99 | 1.00 | 2.16 | 0.0 | 0.8 | 51abe0da |

### Single complex (`c`)

| routine | p | n | inc | rstd runs | from | openblas | rstd | asm | vec x4 | rstd/ob | asm/ob | x4/ob | rstd-ob | rstd err | rstd bits |
| ------- | - | -: | --: | --------- | ---- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --- |
| axpy | c | 512 | 1 | asm | vortex | 141.4 | 234.6 | 280.2 | 284.3 | 1.66 | 1.98 | 2.01 | 0.0 | 0.9 | f7304c9f |
| axpy | c | 4096 | 1 | asm | vortex | 145.2 | 276.6 | 279.9 | 285.9 | 1.90 | 1.93 | 1.97 | 0.0 | 0.9 | d34874c1 |
| axpy | c | 131072 | 1 | asm | vortex | 146.4 | 165.2 | 165.6 | 165.6 | 1.13 | 1.13 | 1.13 | 0.0 | 1.0 | 609be9d3 |
| axpy | c | 8388608 | 1 | asm | vortex | 140.0 | 145.6 | 146.1 | 145.7 | 1.04 | 1.04 | 1.04 | 1.0 | 1.0 | 2cd831c8 |
| axpy | c | 2048 | 2 | asm | vortex | 108.7 | 102.7 | 104.6 | 104.9 | 0.94 | 0.96 | 0.97 | 0.0 | 0.9 | bedda113 |
| dot | c | 512 | 1 | asm | vortex | 100.4 | 138.5 | 140.7 | 139.3 | 1.38 | 1.40 | 1.39 | 0.0 | 0.0 | 15bdf5a3 |
| dot | c | 4096 | 1 | asm | vortex | 105.3 | 172.6 | 173.3 | 181.8 | 1.64 | 1.65 | 1.73 | 0.0 | 0.0 | 4d0b37c7 |
| dot | c | 131072 | 1 | asm | vortex | 106.5 | 139.7 | 142.4 | 142.3 | 1.31 | 1.34 | 1.34 | 0.0 | 0.0 | e52efbff |
| dot | c | 8388608 | 1 | asm | vortex | 106.3 | 123.4 | 123.4 | 123.2 | 1.16 | 1.16 | 1.16 | 0.0 | 0.0 | 555f75e6 |
| dot | c | 2048 | 2 | asm | vortex | 9.8 | 41.4 | 42.0 | 40.7 | 4.22 | 4.29 | 4.15 | 0.2 | 0.0 | d1f88c4d |
| scal | c | 512 | 1 | asm | vortex | 94.6 | 159.8 | 190.1 | 190.0 | 1.69 | 2.01 | 2.01 | 0.7 | 0.4 | e857fb3f |
| scal | c | 4096 | 1 | asm | vortex | 96.8 | 187.0 | 191.3 | 190.8 | 1.93 | 1.98 | 1.97 | 0.8 | 0.5 | dedcb883 |
| scal | c | 131072 | 1 | asm | vortex | 97.7 | 157.5 | 157.6 | 157.1 | 1.61 | 1.61 | 1.61 | 0.8 | 0.6 | 4b03607f |
| scal | c | 8388608 | 1 | asm | vortex | 97.3 | 138.7 | 138.7 | 138.8 | 1.43 | 1.43 | 1.43 | 0.8 | 0.6 | a8f7a3ba |
| scal | c | 2048 | 2 | asm | vortex | 41.6 | 69.7 | 71.6 | 70.9 | 1.68 | 1.72 | 1.70 | 0.8 | 0.5 | b0699f13 |
| nrm2 | c | 512 | 1 | asm | vortex | 5.1 | 70.9 | 71.1 | 70.8 | 13.90 | 13.94 | 13.88 | 1.8 | 0.2 | 8a414b3b |
| nrm2 | c | 4096 | 1 | asm | vortex | 4.8 | 73.0 | 73.1 | 73.0 | 15.21 | 15.23 | 15.21 | 12.9 | 0.0 | 5fb36603 |
| nrm2 | c | 131072 | 1 | asm | vortex | 4.8 | 73.3 | 73.3 | 73.3 | 15.27 | 15.27 | 15.27 | 220.2 | 0.2 | 488d3ca7 |
| nrm2 | c | 8388608 | 1 | asm | vortex | 4.7 | 73.1 | 73.2 | 73.2 | 15.55 | 15.57 | 15.57 | 87201.6 | 0.3 | 7b1207c6 |
| nrm2 | c | 2048 | 2 | asm | vortex | 4.8 | 21.8 | 21.9 | 22.2 | 4.54 | 4.56 | 4.62 | 0.9 | 0.4 | 7eb18e21 |
| asum | c | 512 | 1 | asm | vortex | 143.2 | 120.1 | 127.1 | 135.8 | 0.84 | 0.89 | 0.95 | 0.5 | 0.6 | ba9e885c |
| asum | c | 4096 | 1 | asm | vortex | 146.1 | 141.9 | 142.9 | 144.7 | 0.97 | 0.98 | 0.99 | 0.5 | 0.4 | aed43f7e |
| asum | c | 131072 | 1 | asm | vortex | 146.4 | 136.6 | 139.5 | 140.8 | 0.93 | 0.95 | 0.96 | 0.5 | 0.0 | f2bb55a3 |
| asum | c | 8388608 | 1 | asm | vortex | 110.6 | 109.9 | 109.9 | 110.1 | 0.99 | 0.99 | 1.00 | 6.0 | 0.2 | 40b71501 |
| asum | c | 2048 | 2 | asm | vortex | 15.7 | 27.6 | 28.9 | 29.5 | 1.76 | 1.84 | 1.88 | 0.5 | 0.3 | 873eff2b |
| iamax | c | 512 | 1 | asm | vortex | 100.1 | 68.5 | 73.5 | 74.5 | 0.68 | 0.73 | 0.74 | 0.0 | 0.0 | 7c523b83 |
| iamax | c | 4096 | 1 | asm | vortex | 134.2 | 110.4 | 111.4 | 113.7 | 0.82 | 0.83 | 0.85 | 0.0 | 0.0 | a22cbbee |
| iamax | c | 131072 | 1 | asm | vortex | 141.6 | 109.1 | 112.6 | 115.2 | 0.77 | 0.80 | 0.81 | 0.0 | 0.0 | ba474698 |
| iamax | c | 8388608 | 1 | asm | vortex | 106.1 | 102.5 | 102.8 | 105.2 | 0.97 | 0.97 | 0.99 | 0.0 | 0.0 | 48ba82fe |
| iamax | c | 2048 | 2 | asm | vortex | 18.8 | 43.3 | 44.0 | 44.1 | 2.30 | 2.34 | 2.35 | 0.0 | 0.0 | 36357704 |
| copy | c | 512 | 1 | asm | vortex | 132.1 | 193.2 | 276.9 | 274.0 | 1.46 | 2.10 | 2.07 | 0.0 | 0.0 | 238f9f91 |
| copy | c | 4096 | 1 | asm | vortex | 144.6 | 251.6 | 287.5 | 282.9 | 1.74 | 1.99 | 1.96 | 0.0 | 0.0 | 2984510f |
| copy | c | 131072 | 1 | asm | vortex | 140.4 | 150.0 | 151.3 | 149.7 | 1.07 | 1.08 | 1.07 | 0.0 | 0.0 | e162e98c |
| copy | c | 8388608 | 1 | asm | vortex | 121.8 | 122.2 | 122.1 | 122.1 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | f6a3c87d |
| copy | c | 2048 | 2 | asm | vortex | 72.6 | 70.5 | 72.6 | 137.0 | 0.97 | 1.00 | 1.89 | 0.0 | 0.0 | df37064b |
| swap | c | 512 | 1 | asm | vortex | 36.1 | 88.9 | 92.7 | 95.5 | 2.46 | 2.57 | 2.65 | 0.0 | 0.0 | c9422ad8 |
| swap | c | 4096 | 1 | asm | vortex | 93.7 | 181.4 | 186.9 | 187.0 | 1.94 | 1.99 | 2.00 | 0.0 | 0.0 | 43a5fd90 |
| swap | c | 131072 | 1 | asm | vortex | 84.3 | 158.5 | 158.5 | 157.2 | 1.88 | 1.88 | 1.86 | 0.0 | 0.0 | bd1677e4 |
| swap | c | 8388608 | 1 | asm | vortex | 133.9 | 150.9 | 149.7 | 150.6 | 1.13 | 1.12 | 1.12 | 0.0 | 0.0 | d00aab0b |
| swap | c | 2048 | 2 | asm | vortex | 52.6 | 52.0 | 52.5 | 52.4 | 0.99 | 1.00 | 1.00 | 0.0 | 0.0 | f17b1e4a |
| rot | c | 512 | 1 | asm | vortex | 61.7 | 67.6 | 62.4 | 90.9 | 1.10 | 1.01 | 1.47 | 0.0 | 0.4 | be269a89 |
| rot | c | 4096 | 1 | asm | vortex | 145.5 | 156.0 | 145.7 | 166.8 | 1.07 | 1.00 | 1.15 | 0.0 | 0.5 | 24c0a740 |
| rot | c | 131072 | 1 | asm | vortex | 156.3 | 156.2 | 156.4 | 155.8 | 1.00 | 1.00 | 1.00 | 0.0 | 0.5 | ee260036 |
| rot | c | 8388608 | 1 | asm | vortex | 151.2 | 151.0 | 150.9 | 150.9 | 1.00 | 1.00 | 1.00 | 0.0 | 0.6 | e454dacb |
| rot | c | 2048 | 2 | asm | vortex | 51.3 | 51.0 | 51.6 | 49.2 | 0.99 | 1.01 | 0.96 | 0.0 | 0.5 | 246c280f |

### Double complex (`z`)

| routine | p | n | inc | rstd runs | from | openblas | rstd | asm | vec x4 | rstd/ob | asm/ob | x4/ob | rstd-ob | rstd err | rstd bits |
| ------- | - | -: | --: | --------- | ---- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --- |
| axpy | z | 256 | 1 | asm | vortex | 141.4 | 233.3 | 283.3 | 283.3 | 1.65 | 2.00 | 2.00 | 0.0 | 1.0 | 544d4934 |
| axpy | z | 2048 | 1 | asm | vortex | 145.1 | 272.1 | 282.6 | 286.1 | 1.88 | 1.95 | 1.97 | 0.0 | 1.0 | 4a164f05 |
| axpy | z | 65536 | 1 | asm | vortex | 146.5 | 165.4 | 165.5 | 165.5 | 1.13 | 1.13 | 1.13 | 0.0 | 1.0 | c3b3756c |
| axpy | z | 4194304 | 1 | asm | vortex | 140.5 | 146.5 | 147.1 | 146.6 | 1.04 | 1.05 | 1.04 | 1.0 | 1.0 | 9d706a6e |
| axpy | z | 1024 | 2 | asm | vortex | 215.9 | 198.3 | 198.8 | 208.9 | 0.92 | 0.92 | 0.97 | 0.0 | 1.0 | 8899d47f |
| dot | z | 256 | 1 | asm | vortex | 100.6 | 152.3 | 152.7 | 163.1 | 1.51 | 1.52 | 1.62 | 0.0 | 0.1 | edb8cdc5 |
| dot | z | 2048 | 1 | asm | vortex | 105.3 | 172.6 | 175.2 | 186.4 | 1.64 | 1.66 | 1.77 | 0.1 | 0.1 | 74732dd6 |
| dot | z | 65536 | 1 | asm | vortex | 106.5 | 141.0 | 142.5 | 142.1 | 1.32 | 1.34 | 1.33 | 0.0 | 0.1 | fbf27725 |
| dot | z | 4194304 | 1 | asm | vortex | 106.3 | 123.6 | 123.5 | 123.4 | 1.16 | 1.16 | 1.16 | 0.1 | 0.0 | 658fd268 |
| dot | z | 1024 | 2 | asm | vortex | 20.2 | 83.4 | 85.3 | 82.1 | 4.13 | 4.22 | 4.06 | 0.1 | 0.2 | 9a68f762 |
| scal | z | 256 | 1 | asm | vortex | 94.5 | 159.3 | 190.2 | 189.0 | 1.69 | 2.01 | 2.00 | 0.7 | 0.7 | a86de12e |
| scal | z | 2048 | 1 | asm | vortex | 96.8 | 186.5 | 191.4 | 190.4 | 1.93 | 1.98 | 1.97 | 0.8 | 0.8 | 05e4cc32 |
| scal | z | 65536 | 1 | asm | vortex | 97.7 | 157.4 | 157.6 | 157.1 | 1.61 | 1.61 | 1.61 | 0.8 | 0.8 | b974bd94 |
| scal | z | 4194304 | 1 | asm | vortex | 97.4 | 138.6 | 138.8 | 139.1 | 1.42 | 1.43 | 1.43 | 0.8 | 0.8 | 2a1d74d7 |
| scal | z | 1024 | 2 | asm | vortex | 83.0 | 70.7 | 71.9 | 137.7 | 0.85 | 0.87 | 1.66 | 0.7 | 0.7 | 9fe85bef |
| nrm2 | z | 256 | 1 | asm | vortex | 10.6 | 173.5 | 177.7 | 177.4 | 16.37 | 16.76 | 16.74 | 1.2 | 0.6 | 3e7b7ad9 |
| nrm2 | z | 2048 | 1 | asm | vortex | 9.7 | 162.0 | 166.1 | 163.2 | 16.70 | 17.12 | 16.82 | 4.3 | 1.7 | 21d8e8ce |
| nrm2 | z | 65536 | 1 | asm | vortex | 9.3 | 147.5 | 147.6 | 147.7 | 15.86 | 15.87 | 15.88 | 24.5 | 6.7 | 14746a82 |
| nrm2 | z | 4194304 | 1 | asm | vortex | 9.3 | 110.8 | 110.5 | 110.5 | 11.91 | 11.88 | 11.88 | 304.3 | 28.2 | a6e8b866 |
| nrm2 | z | 1024 | 2 | asm | vortex | 9.8 | 51.1 | 51.3 | 52.5 | 5.21 | 5.23 | 5.36 | 1.8 | 0.0 | 69e44779 |
| asum | z | 256 | 1 | asm | vortex | 144.1 | 127.5 | 140.8 | 142.6 | 0.88 | 0.98 | 0.99 | 0.5 | 1.0 | 4b8ce16e |
| asum | z | 2048 | 1 | asm | vortex | 146.3 | 142.6 | 145.5 | 146.1 | 0.97 | 0.99 | 1.00 | 1.5 | 9.1 | 3a292c62 |
| asum | z | 65536 | 1 | asm | vortex | 146.3 | 137.4 | 139.4 | 140.9 | 0.94 | 0.95 | 0.96 | 3.0 | 1.0 | 3fed70f4 |
| asum | z | 4194304 | 1 | asm | vortex | 110.4 | 110.1 | 109.7 | 110.0 | 1.00 | 0.99 | 1.00 | 4.5 | 299.0 | 6792d763 |
| asum | z | 1024 | 2 | asm | vortex | 35.3 | 58.4 | 61.0 | 61.4 | 1.65 | 1.73 | 1.74 | 1.5 | 1.5 | e5229744 |
| iamax | z | 256 | 1 | asm | vortex | 116.2 | 81.3 | 89.7 | 85.9 | 0.70 | 0.77 | 0.74 | 0.0 | 0.0 | 56bbea09 |
| iamax | z | 2048 | 1 | asm | vortex | 133.0 | 112.4 | 114.7 | 114.8 | 0.85 | 0.86 | 0.86 | 0.0 | 0.0 | b895d41e |
| iamax | z | 65536 | 1 | asm | vortex | 133.8 | 112.8 | 115.3 | 115.7 | 0.84 | 0.86 | 0.86 | 0.0 | 0.0 | ba474698 |
| iamax | z | 4194304 | 1 | asm | vortex | 108.0 | 103.9 | 104.0 | 105.2 | 0.96 | 0.96 | 0.97 | 0.0 | 0.0 | 48ba82fe |
| iamax | z | 1024 | 2 | asm | vortex | 38.6 | 75.8 | 76.9 | 76.8 | 1.96 | 1.99 | 1.99 | 0.0 | 0.0 | 4301daa9 |
| copy | z | 256 | 1 | asm | vortex | 130.5 | 193.2 | 277.2 | 274.3 | 1.48 | 2.12 | 2.10 | 0.0 | 0.0 | 846f68e3 |
| copy | z | 2048 | 1 | asm | vortex | 144.1 | 263.6 | 284.2 | 279.3 | 1.83 | 1.97 | 1.94 | 0.0 | 0.0 | 03334116 |
| copy | z | 65536 | 1 | asm | vortex | 140.4 | 147.5 | 149.0 | 147.9 | 1.05 | 1.06 | 1.05 | 0.0 | 0.0 | 3eacc605 |
| copy | z | 4194304 | 1 | asm | vortex | 122.0 | 122.2 | 122.1 | 122.4 | 1.00 | 1.00 | 1.00 | 0.0 | 0.0 | 6758600f |
| copy | z | 1024 | 2 | asm | vortex | 144.8 | 133.9 | 141.4 | 228.2 | 0.92 | 0.98 | 1.58 | 0.0 | 0.0 | ba9ea3c3 |
| swap | z | 256 | 1 | asm | vortex | 183.1 | 228.4 | 266.7 | 266.6 | 1.25 | 1.46 | 1.46 | 0.0 | 0.0 | b656d6ca |
| swap | z | 2048 | 1 | asm | vortex | 125.9 | 217.4 | 223.4 | 223.2 | 1.73 | 1.77 | 1.77 | 0.0 | 0.0 | 736fd008 |
| swap | z | 65536 | 1 | asm | vortex | 85.5 | 158.6 | 158.5 | 157.2 | 1.85 | 1.85 | 1.84 | 0.0 | 0.0 | 36fa1c50 |
| swap | z | 4194304 | 1 | asm | vortex | 122.3 | 146.9 | 145.4 | 146.6 | 1.20 | 1.19 | 1.20 | 0.0 | 0.0 | fe45c39e |
| swap | z | 1024 | 2 | asm | vortex | 138.4 | 135.4 | 133.8 | 133.6 | 0.98 | 0.97 | 0.97 | 0.0 | 0.0 | 6c5680c2 |
| rot | z | 256 | 1 | asm | vortex | 253.6 | 219.7 | 255.1 | 227.2 | 0.87 | 1.01 | 0.90 | 0.0 | 0.7 | 7005d051 |
| rot | z | 2048 | 1 | asm | vortex | 189.9 | 197.0 | 190.1 | 194.9 | 1.04 | 1.00 | 1.03 | 0.0 | 0.8 | 6be51bfb |
| rot | z | 65536 | 1 | asm | vortex | 157.4 | 157.4 | 157.5 | 156.4 | 1.00 | 1.00 | 0.99 | 0.0 | 0.8 | e4fb385e |
| rot | z | 4194304 | 1 | asm | vortex | 144.2 | 144.2 | 144.0 | 147.3 | 1.00 | 1.00 | 1.02 | 0.0 | 0.8 | a426f250 |
| rot | z | 1024 | 2 | asm | vortex | 133.1 | 133.8 | 132.6 | 133.0 | 1.01 | 1.00 | 1.00 | 0.0 | 0.8 | e18ade79 |

## Level 2

`gemv`, square, unit strides, `beta` 1: GFLOP/s.

| routine | p | n | rstd runs | from | openblas | rstd | asm | vec x4 | rstd/ob | asm/ob | x4/ob | rstd-ob | rstd err | rstd bits |
| ------- | - | -: | --------- | ---- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --- |
| gemv_n | s | 90 | asm | vortex | 66.3 | 57.5 | 61.1 | 61.5 | 0.87 | 0.92 | 0.93 | 0.6 | 6.9 | ec887ee3 |
| gemv_n | s | 1024 | asm | vortex | 51.7 | 52.0 | 52.1 | 52.0 | 1.01 | 1.01 | 1.01 | 1.7 | 15.8 | b6d4008b |
| gemv_n | s | 4096 | asm | vortex | 33.0 | 36.3 | 36.3 | 36.3 | 1.10 | 1.10 | 1.10 | 1.3 | 19.6 | 384bdb2b |
| gemv_t | s | 90 | asm | vortex | 37.4 | 39.9 | 41.5 | 40.7 | 1.07 | 1.11 | 1.09 | 0.8 | 0.4 | 4cdae9fe |
| gemv_t | s | 1024 | asm | vortex | 40.3 | 53.8 | 53.9 | 53.9 | 1.33 | 1.34 | 1.34 | 0.5 | 0.3 | 2e17dcca |
| gemv_t | s | 4096 | asm | vortex | 39.1 | 38.4 | 38.3 | 38.4 | 0.98 | 0.98 | 0.98 | 0.5 | 0.3 | f28754a7 |
| gemv_n | d | 64 | asm | vortex | 20.0 | 34.4 | 37.0 | 36.8 | 1.72 | 1.85 | 1.84 | 0.0 | 4.1 | ba428633 |
| gemv_n | d | 724 | asm | vortex | 19.0 | 28.9 | 28.9 | 28.9 | 1.52 | 1.52 | 1.52 | 0.0 | 11.9 | d500cfb8 |
| gemv_n | d | 2896 | asm | vortex | 18.9 | 19.5 | 19.5 | 19.5 | 1.03 | 1.03 | 1.03 | 0.0 | 19.6 | ee62fc60 |
| gemv_t | d | 64 | asm | vortex | 19.7 | 34.5 | 36.0 | 35.7 | 1.75 | 1.83 | 1.81 | 0.0 | 0.6 | f0b4d993 |
| gemv_t | d | 724 | asm | vortex | 20.1 | 29.7 | 29.7 | 29.7 | 1.48 | 1.48 | 1.48 | 0.6 | 0.9 | 2f3d70c0 |
| gemv_t | d | 2896 | asm | vortex | 19.6 | 19.7 | 19.7 | 19.7 | 1.01 | 1.01 | 1.01 | 0.4 | 1.3 | 2aa63d1b |
| gemv_n | c | 64 | asm | vortex | 34.6 | 61.9 | 63.6 | 65.5 | 1.79 | 1.84 | 1.89 | 1.3 | 4.6 | d24de882 |
| gemv_n | c | 724 | asm | vortex | 35.8 | 80.3 | 80.7 | 82.9 | 2.24 | 2.25 | 2.32 | 2.5 | 7.7 | dd8cf00f |
| gemv_n | c | 2896 | asm | vortex | 36.3 | 60.2 | 60.3 | 61.0 | 1.66 | 1.66 | 1.68 | 1.2 | 9.0 | ee261e78 |
| gemv_t | c | 64 | asm | vortex | 45.5 | 49.3 | 50.4 | 49.9 | 1.08 | 1.11 | 1.10 | 0.8 | 0.6 | b6e86b45 |
| gemv_t | c | 724 | asm | vortex | 25.7 | 82.4 | 82.9 | 83.7 | 3.21 | 3.23 | 3.26 | 0.5 | 0.3 | b51f3824 |
| gemv_t | c | 2896 | asm | vortex | 22.6 | 64.5 | 64.5 | 64.6 | 2.85 | 2.85 | 2.86 | 0.2 | 0.2 | 33a07a90 |
| gemv_n | z | 45 | asm | vortex | 23.2 | 33.2 | 34.0 | 35.1 | 1.43 | 1.47 | 1.51 | 0.9 | 5.7 | 4721cbb7 |
| gemv_n | z | 512 | asm | vortex | 24.3 | 39.3 | 39.5 | 40.6 | 1.62 | 1.63 | 1.67 | 2.5 | 8.4 | b1d1864c |
| gemv_n | z | 2048 | asm | vortex | 24.1 | 33.2 | 33.1 | 33.7 | 1.38 | 1.37 | 1.40 | 1.9 | 9.3 | 748372da |
| gemv_t | z | 45 | asm | vortex | 27.5 | 24.9 | 25.4 | 25.3 | 0.91 | 0.92 | 0.92 | 0.9 | 0.9 | 96a49e28 |
| gemv_t | z | 512 | asm | vortex | 22.9 | 41.7 | 41.8 | 42.5 | 1.82 | 1.83 | 1.86 | 0.5 | 0.6 | cdda4ebe |
| gemv_t | z | 2048 | asm | vortex | 21.4 | 37.2 | 37.0 | 37.2 | 1.74 | 1.73 | 1.74 | 0.2 | 0.3 | ee2ce6a7 |

## Level 3

`gemm` NN, square, `C += A*B`, one thread: GFLOP/s.

| p | n | openblas | rstd | rstd/ob | rstd-ob | rstd err | rstd bits |
| - | -: | -: | -: | -: | -: | -: | --- |
| s | 64 | 120.9 | 113.6 | 0.94 | 0.0 | 0.6 | 4c2c8b1f |
| s | 128 | 123.4 | 118.9 | 0.96 | 0.0 | 0.6 | ab770b69 |
| s | 256 | 124.6 | 123.6 | 0.99 | 0.0 | 1.0 | 21fdf633 |
| s | 512 | 126.0 | 126.3 | 1.00 | 0.0 | 0.9 | 45469416 |
| s | 1024 | 127.5 | 128.0 | 1.00 | 2.3 | 0.8 | eb97f50f |
| s | 2048 | 127.6 | 127.7 | 1.00 | 1.5 | 0.4 | fc6facb0 |
| d | 64 | 59.8 | 56.4 | 0.94 | 0.0 | 1.1 | 3fe90444 |
| d | 128 | 60.5 | 59.1 | 0.98 | 0.0 | 0.7 | 92dfae0a |
| d | 256 | 62.1 | 61.9 | 1.00 | 0.0 | 0.8 | 0040fb73 |
| d | 512 | 63.5 | 63.0 | 0.99 | 2.1 | 0.6 | a2b5517a |
| d | 1024 | 64.0 | 63.7 | 1.00 | 1.5 | 0.9 | 6898ad11 |
| d | 2048 | 64.0 | 63.6 | 0.99 | 0.8 | 1.3 | 34bb3f10 |
| c | 64 | 109.2 | 115.5 | 1.06 | 0.8 | 0.3 | a7b61b84 |
| c | 128 | 111.5 | 122.2 | 1.10 | 0.9 | 0.3 | 670b5a37 |
| c | 256 | 113.0 | 126.0 | 1.12 | 0.6 | 0.4 | d7a63c6e |
| c | 512 | 114.5 | 128.3 | 1.12 | 1.0 | 0.3 | 92ad3d6b |
| c | 1024 | 114.9 | 128.4 | 1.12 | 0.5 | 0.2 | 7763b208 |
| c | 2048 | 115.0 | 128.7 | 1.12 | 0.3 | 0.1 | 71fbc007 |
| z | 64 | 54.1 | 57.5 | 1.06 | 1.1 | 0.5 | f96643da |
| z | 128 | 54.8 | 60.9 | 1.11 | 0.7 | 0.4 | 250fd4d6 |
| z | 256 | 56.1 | 62.7 | 1.12 | 0.7 | 0.5 | 5e75ac7b |
| z | 512 | 56.7 | 64.0 | 1.13 | 0.7 | 0.6 | d6fe629e |
| z | 1024 | 56.7 | 63.9 | 1.13 | 0.5 | 0.5 | 25a060ab |
| z | 2048 | 56.8 | 64.2 | 1.13 | 0.3 | 0.4 | d4c7d3e9 |

### Default threading

OpenBLAS at its 6 threads, rstd on `Io.Threaded` up to the 6 CPUs.
Indicative only. rstd's threaded `gemm` reaches 0.86 to 0.95 of OpenBLAS's
in real from n = 512 and passes it in complex from n = 256 (1.05 to 1.15);
at n = 64 and 128 its threads cost more than they bring (0.21 to 0.96).
Its rate at `c` 128 varied from 198 to 249 GFLOP/s over five runs
(OpenBLAS's 340.3 to 340.7), so that row is noisy.

| p | n | openblas | rstd | rstd/ob | rstd-ob | rstd err | rstd bits |
| - | -: | -: | -: | -: | -: | -: | --- |
| s | 64 | 123.1 | 117.9 | 0.96 | 0.0 | 0.6 | 4c2c8b1f |
| s | 128 | 293.5 | 97.0 | 0.33 | 0.0 | 0.6 | ab770b69 |
| s | 256 | 442.7 | 262.4 | 0.59 | 0.0 | 1.0 | 21fdf633 |
| s | 512 | 497.4 | 471.3 | 0.95 | 0.0 | 0.9 | 45469416 |
| s | 1024 | 561.2 | 501.6 | 0.89 | 2.3 | 0.8 | eb97f50f |
| s | 2048 | 587.4 | 527.0 | 0.90 | 1.5 | 0.4 | fc6facb0 |
| d | 64 | 59.8 | 56.3 | 0.94 | 0.0 | 1.1 | 3fe90444 |
| d | 128 | 194.2 | 80.8 | 0.42 | 0.0 | 0.7 | 92dfae0a |
| d | 256 | 252.5 | 184.0 | 0.73 | 0.0 | 0.8 | 0040fb73 |
| d | 512 | 271.0 | 241.2 | 0.89 | 2.1 | 0.6 | a2b5517a |
| d | 1024 | 289.2 | 253.2 | 0.88 | 1.5 | 0.9 | 6898ad11 |
| d | 2048 | 303.1 | 260.8 | 0.86 | 0.8 | 1.3 | 34bb3f10 |
| c | 64 | 270.5 | 57.1 | 0.21 | 0.8 | 0.3 | a7b61b84 |
| c (noisy) | 128 | 340.7 | 249.4 | 0.73 | 0.9 | 0.3 | 670b5a37 |
| c | 256 | 427.7 | 447.2 | 1.05 | 0.6 | 0.4 | d7a63c6e |
| c | 512 | 448.8 | 511.5 | 1.14 | 1.0 | 0.3 | 92ad3d6b |
| c | 1024 | 477.2 | 529.4 | 1.11 | 0.5 | 0.2 | 7763b208 |
| c | 2048 | 494.3 | 570.4 | 1.15 | 0.3 | 0.1 | 71fbc007 |
| z | 64 | 163.2 | 50.9 | 0.31 | 1.1 | 0.5 | f96643da |
| z | 128 | 195.5 | 158.6 | 0.81 | 0.7 | 0.4 | 250fd4d6 |
| z | 256 | 225.2 | 240.8 | 1.07 | 0.7 | 0.5 | 5e75ac7b |
| z | 512 | 231.5 | 249.4 | 1.08 | 0.7 | 0.6 | d6fe629e |
| z | 1024 | 239.5 | 254.0 | 1.06 | 0.5 | 0.5 | 25a060ab |
| z | 2048 | 245.7 | 271.7 | 1.11 | 0.3 | 0.4 | d4c7d3e9 |

## Run conditions

The Mac's 1-minute, 5-minute and 15-minute loads, and `top`'s 8 busiest
processes in % of a core, before and after each kept run (the VM itself is
`qemu-system`, gone by the time the after sample is taken). `lx` is the bench at
one thread, `lx-thr` with `--threads=default`; `lx-3` reran levels 1 and
3, `lx-thr-3` to `-5` level 3.

- **lx-1**, 2026-10-29 23:29:10 CDT to 23:33:49: load 1.13 / 1.26 / 1.49 before, 1.09 / 1.18 / 1.39 after; the watch saw `stickersd` at 9.2% (23:32:36). Before: `Activity Monitor` 3.8, `sysmond` 1.7, `top` 1.4, `remindd` 0.5, `kernel_task` 0.3, `AppleCentauriAlp` 0.1, `SetStoreUpdateService` 0.1, `WindowServer` 0.1. After: `Activity Monitor` 3.8, `sysmond` 1.7, `top` 1.5, `kernel_task` 0.7, `stickersd` 0.3, `WindowServer` 0.1, `AppleCentauriAlp` 0.1, `clocksyncd` 0.1.
- **lx-2**, 2026-10-29 23:38:50 CDT to 23:43:28: load 1.19 / 1.24 / 1.35 before, 1.14 / 1.20 / 1.29 after. Before: `Activity Monitor` 3.8, `sysmond` 1.7, `top` 1.5, `kernel_task` 0.4, `stickersd` 0.2, `AppleCentauriAlp` 0.2, `WindowServer` 0.1, `clocksyncd` 0.1. After: `Activity Monitor` 3.7, `sysmond` 1.8, `top` 1.5, `launchd` 1.1, `kernel_task` 0.9, `Bitwarden Helper` 0.3, `AppleCentauriAlp` 0.3, `clocksyncd` 0.2.
- **lx-3**, 2026-10-29 23:43:38 CDT to 23:47:51: load 1.12 / 1.19 / 1.29 before, 1.14 / 1.17 / 1.24 after. Before: `Activity Monitor` 3.6, `sysmond` 1.7, `top` 1.5, `runningboardd` 0.5, `kernel_task` 0.4, `AppleCentauriAlp` 0.2, `Bitwarden Helper` 0.2, `airportd` 0.1. After: `Activity Monitor` 3.6, `sysmond` 1.8, `top` 1.4, `kernel_task` 0.7, `locationd` 0.3, `AppleCentauriAlp` 0.1, `clocksyncd` 0.1, `ControlCenter` 0.1.
- **lx-thr-1**, 2026-10-29 23:48:02 CDT to 23:48:30: load 1.12 / 1.16 / 1.24 before, 2.12 / 1.37 / 1.31 after. Before: `Activity Monitor` 3.6, `sysmond` 1.7, `top` 1.4, `kernel_task` 0.4, `AppleCentauriAlp` 0.1, `clocksyncd` 0.1, `locationd` 0.1, `ControlCenter` 0.1. After: `Activity Monitor` 4.2, `sysmond` 1.7, `top` 1.5, `launchd` 1.1, `kernel_task` 0.6, `clocksyncd` 0.2, `AppleCentauriAlp` 0.2, `dasd` 0.2.
- **lx-thr-2**, 2026-10-29 23:50:30 CDT to 23:50:58: load 1.90 / 1.68 / 1.45 before, 3.23 / 1.99 / 1.56 after. Before: `Activity Monitor` 3.7, `sysmond` 1.7, `top` 1.5, `launchd` 1.1, `kernel_task` 0.4, `clocksyncd` 0.2, `AppleCentauriAlp` 0.1, `WindowServer` 0.1. After: `Activity Monitor` 4.0, `loginwindow` 2.4, `sysmond` 1.8, `top` 1.5, `launchd` 1.1, `kernel_task` 1.1, `Bitwarden` 0.9, `WindowServer` 0.7.
- **lx-thr-3**, 2026-10-29 23:52:19 CDT to 23:52:47: load 1.74 / 1.84 / 1.54 before, 2.64 / 2.03 / 1.62 after. Before: `Activity Monitor` 3.7, `sysmond` 1.8, `top` 1.4, `kernel_task` 0.4, `AppleCentauriAlp` 0.2, `opendirectoryd` 0.2, `trustd` 0.2, `Bitwarden Helper` 0.1. After: `Activity Monitor` 3.1, `sysmond` 1.7, `top` 1.4, `dasd` 0.7, `kernel_task` 0.6, `Bitwarden Helper` 0.3, `searchpartyd` 0.2, `locationd` 0.2.
- **lx-thr-4**, 2026-10-29 23:54:08 CDT to 23:54:36: load 1.76 / 1.89 / 1.60 before, 3.24 / 2.21 / 1.72 after. Before: `Activity Monitor` 3.9, `sysmond` 1.7, `top` 1.5, `kernel_task` 0.4, `Bitwarden Helper` 0.3, `AppleCentauriAlp` 0.1, `WindowServer` 0.1, `clocksyncd` 0.1. After: `Activity Monitor` 3.8, `sysmond` 1.8, `top` 1.5, `kernel_task` 0.5, `Bitwarden Helper` 0.2, `AppleCentauriAlp` 0.2, `clocksyncd` 0.1, `airportd` 0.1.
- **lx-thr-5**, 2026-10-29 23:55:57 CDT to 23:56:25: load 1.91 / 2.01 / 1.69 before, 2.79 / 2.19 / 1.76 after. Before: `Activity Monitor` 4.0, `loginwindow` 2.5, `Bitwarden Helper` 2.2, `sysmond` 1.7, `top` 1.5, `launchd` 1.1, `WindowServer` 0.7, `kernel_task` 0.7. After: `Activity Monitor` 3.7, `sysmond` 1.8, `top` 1.5, `kernel_task` 0.6, `Bitwarden Helper` 0.3, `AppleCentauriAlp` 0.1, `clocksyncd` 0.1, `WindowServer` 0.1.

## Against the busy Mac

The 2026-10-29 runs (02:02 to 02:35 CDT, `review` `6d66239`) took a host
load of 2.6 to 9.1. Against them, OpenBLAS's median ratio is 1.000 (0.98
to 1.06, the top end large `asum`), the converted and `@Vector` kernels'
1.001: the background load had not skewed the one-thread numbers.
Threaded OpenBLAS is within 5%. What moved is code:

- **rstd's call overhead:** its public routine runs 7 to 30% behind the
  kernel it calls at the smallest n of `axpy`, `scal`, `copy`, `swap`,
  `rot`, `iamax`, `dot` and `asum`, where it matched before; `sdot` from
  256K floats runs 51 to 53 GB/s against the kernel's 98.
- **Strided `saxpy` and `zscal`** (inc 2) halved, rstd's and the converted
  kernel's alike (54 to 28 GB/s, as macOS's; 139 to 71).
- **Default threading:** rstd threads now (it ran one thread before), so its
  `gemm` there rose up to 4.4 times.
- **`rstd bits`** changed for large `dot` and `asum` (13 rows), a new
  summation order.
