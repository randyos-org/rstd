# Support tiers

Which platforms rstd supports, and how far, in the style of Zig's support table
(<https://ziglang.org/download/0.14.0/release-notes.html#Support-Table>). The tiers follow rstd's
consumers, RandyOS and zfls (user, 2026-10-30). Each tier includes the guarantees of the tiers
after it. Rig names are `zig build` steps: `qemu-test-<rig>` unless named in full
([`design.md`](src/linalg/design.md#rigs) has the rigs and their machines).

## Tiers

- **Tier 1: required first.** Everything implemented works and is tested. Bugs block the MVP and every release, and go in the current milestone (the open `dev` to `main` pull request's list of issues). Every kernel object the platform links runs the per-core tests, or is listed below as untested with its item.
- **Tier 2: required, after tier 1.** RandyOS's and zfls's other platforms. Bugs are required fixes, also for the MVP, and go in the current milestone too. Tested on a rig or machine where one exists; untested cores and paths are listed with their items.
- **Tier 3: best effort.** Every other platform with a rig. The rigs run in each sweep for the record, and CI keeps the tier 3 and 4 legs, the tier 3 rigs among them, but skips them for now: work is on tiers 1 and 2 only (user, 2026-11-20), and a run started by hand runs them on request (`.github/workflows/ci.yml`); failures go to issues labeled `backlog`, as do tier 4's, tuning, optional and Blocked work. Choices no rig can check follow OpenBLAS and BLIS where they agree
  ([`design.md`](src/linalg/design.md#settled-decisions), "Target platforms come first").
- **Tier 4: unsupported or blocked.** No rig, blocked upstream, or not worth its upkeep. Every
  target compiles with the portable kernels at least; nothing more is promised. Each entry cites
  its item.

## Platforms

| Target | Tier | Consumer | Tested by |
| --- | --- | --- | --- |
| `x86_64-windows` | 1 | zfls | native `rstd-test` on the Zen 3 host, with `-Dblas-cores=all`; `x86_64-whpx` (the host CPU); `rstd-check-zig017`; zfls's CI |
| `spirv32-vulkan` | 1 | rgl, shaderlib | `-Dvulkan` compiles rgl's shaders (`vulkan_v1_2`); `shader-check` validates them with `spirv-val`; `shader-run` draws with them and runs shaderlib's probes on a Vulkan device against the CPU (`tools/shader_run`); both in CI, on Mesa's lavapipe |
| `x86_64-freestanding-none` | 2 | RandyOS: PC, Intel Mac | `x86_64` and `x86_64-smp` (TCG `max`, HVF on the Intel Mac, or KVM in CI), `x86_64-whpx` |
| `aarch64-freestanding-none` | 2 | RandyOS: Raspberry Pi 3 to 5, Apple Silicon | `aarch64` (HVF on the M6, or `cortex-a72`), `aarch64-smp`; in CI on arm64 runners. `aarch64-sve` (TCG `max`) is tier 3: no RandyOS target has SVE (#619) |
| `arm-freestanding-eabi` | 2 | RandyOS: Raspberry Pi 3, 32-bit | `arm-softfp` (soft-float ABI, A15's VFP), `arm` (`eabihf`); both in CI under TCG |
| `powerpc-freestanding-eabi` | 2 | RandyOS: iBook G3, PowerPC 750FX | `powerpc-750` (soft float: portable kernels only), `powerpc-750-hf`; both in CI under TCG |
| `x86_64-uefi-msvc` | 2 | RandyOS's loader | `x86_64-uefi` (OVMF; HVF on the Intel Mac, or KVM in CI), `rstd-check-uefi` (CI too) |
| `aarch64-uefi-msvc` | 2 | RandyOS's loader (pftf) | `aarch64-uefi` (AAVMF; HVF on the M6, or arm64 runners in CI), `rstd-check-uefi`. `aarch64-uefi-sve` (TCG `max`) is tier 3: no RandyOS target has SVE (#619) |
| `x86_64-linux-musl` | 2 | zfls (release, conda, Alpine) | CI's native `x86_64-linux-musl` legs; locally `linux-x86_64` (TCG `max` or HVF on the Intel Mac; 4 CPUs), `rstd-check-zig017` |
| `aarch64-linux-musl` | 2 | zfls (release, conda, Alpine) | CI's native legs on arm64 runners, which are Cobalt 100 (Neoverse N2, SVE2) and run the `neoversen2` kernels (CI's `blas-core` step, #619); locally `linux-aarch64` (HVF on the M6, or TCG `max`); the SVE and SME cores, for SVE servers such as Graviton 3 and 4, in CI's `aarch64-linux-musl, sve` leg under qemu-user's `max` CPU ([#653](https://github.com/randyos-org/rstd/issues/653)), and locally the `linux-aarch64-sve` rig (TCG `max`); `rstd-check-zig017` |
| `x86_64-linux-gnu`, `aarch64-linux-gnu` | 2 | zfls | CI's native legs (`ubuntu-24.04`, `ubuntu-24.04-arm`); locally the Linux rigs above with `-Drstd-tsan` (glibc from `-Dtest-host`'s Nix); in CI both ThreadSanitizer runs are native legs ([#539](https://github.com/randyos-org/rstd/issues/539), [#652](https://github.com/randyos-org/rstd/issues/652)) |
| `aarch64-windows` | 2 | zfls | CI's native legs (`windows-11-arm`), `rstd-check-zig017` |
| `x86_64-macos` 11.0+ | 2 | zfls | `ssh-test-x86_64-macos` (Intel Mac, i5-8500B, macOS 15.8), CI's native legs (`macos-15-intel`), `ssh-test-x86_64-macos-rosetta` (Rosetta 2 on the M6, for comparison), `rstd-check-macho`, `rstd-check-zig017` |
| `aarch64-macos` 11.0+ | 2 | zfls | `ssh-test-aarch64-macos` (M6, macOS 27.2), CI's native legs (`macos-15`), `rstd-check-macho`, `rstd-check-zig017` |
| `x86-windows-gnu` | 3 | | `rstd-test-x87` (`pentium3` code on the Zen 3 host, Debug and ReleaseSafe) |
| `x86-freestanding-none` | 3 | | `i386`, `i386-max`, `i386-safe`, `i386-max-safe` |
| `arm-linux` | 3 | | `linux-arm`; runs with tier 2 in the MVP sweep, for the arm cores under an OS |
| `riscv64-freestanding-none`, `riscv64-linux` | 3 | | `riscv64`, `riscv64-vlen128`, `riscv64-vlen256`, `linux-riscv64` |
| `powerpc64le-freestanding-none`, `powerpc64le-linux` | 3 | | `powerpc64le`, `powerpc64le-power10`, `linux-powerpc64le` |
| `powerpc64-freestanding-none` | 3 | | `powerpc64` (big-endian, POWER10) |
| `powerpc-freestanding-none` (hard float) | 3 | | `powerpc` (G4), `powerpc-pseries`, `powerpc-power9` |
| `mips64el`, `mips64`, `mipsel`, `mips` (freestanding) | 3 | | `mips64r2el`, `mips64r2`, `mips64r6el`, `mips64r6`, `mips32r2el`, `mips32r2`, `mips32r5el`, `mips32r5` |
| `loongarch64-freestanding-none` | 3 | | `loongarch64` |
| `s390x-freestanding-none` | 3 | | `s390x` |
| `wasm32-wasi` | 3 | | `wasm-test` (`wasm-test-simd128`, `wasm-test-portable`) |
| `sparc64-freestanding-none` | 4 | | `sparc64`, `sparc64-debug`, out of the sweeps: their images outgrow OpenBIOS's 4 MiB window (T39) |
| `sparc-freestanding-none` (32-bit) | 4 | | stub: lld cannot link 32-bit SPARC (T48) |
| `x86-linux` | 4 | | `linux-i386` stub: TCG's x87 page-crossing bug (T49); `rstd-test-x87` runs it on an x86 Linux host, which is not at hand |
| `s390x-linux`, `mips64el-linux` | 4 | | no rig: `zig fetch` of single files (T49) |
| `loongarch64-linux`, `powerpc64-linux`, `sparc64-linux`, 32-bit `mips-linux` | 4 | | no Linux rig (T71) |
| `powerpc-macos` (Mach-O) | 4 | | `powerpc-macho` stub: no powerpc-macos Zig port (T50) |
| alpha, ia64, csky, e2k | 4 | | no Zig target or backend ([`readme.md`](src/math/kernels/blas/readme.md#not-supported-by-zig)) |
| every other target (FreeBSD, riscv32, ...) | 4 | | none; compiles with the portable kernels |

Notes:

- **SPIR-V** (`rstd.linalg` under SPIR-V, rgl's shaders, shaderlib) waited for the CPU BLAS
  MVP; S1 and S2 ([#51](https://github.com/randyos-org/rstd/issues/51)) made shaderlib the shaders' one API over `linalg`, and [#353](https://github.com/randyos-org/rstd/issues/353) validates and runs them. The bit guarantees (same bits
  for any thread count, one contract on every core, the first-NaN rule) are the CPU's only;
  under SPIR-V, Vulkan's precision rules apply (design.md, Settled decisions).
- **macOS 11 and 12** never run: the Intel Mac has macOS 15.8 and the M6 27.2.
- **The Intel Mac** replaces Rosetta 2 and the macOS VM (`x86_64-macos`, TCG; retired by [#98](https://github.com/randyos-org/rstd/issues/98)) as x86_64 macOS coverage; Rosetta 2's steps stay for comparison.
- **Hardware-only paths** (T47): Linux ppc64 `mfpvr`, the mips `cpucfg` decode, the sparc `sllx`
  rewrite. The Windows arm64 registry MIDR and feature probes and the arm64 sysfs fallback are
  tested with recorded inputs (design.md, Core detection and selection), not on hardware.

## Accelerate (tier 2)

The macOS backend (`apple_accel.zig`: the binding, routing, quirks, and LAPACK's `potrf`) is a
tier-2 component of both macOS platforms.

| Path | Tested by |
| --- | --- |
| arm64, macOS 13.3+ ILP64 symbols | `ssh-test-aarch64-macos` |
| x86_64, macOS 13.3+ symbols | `ssh-test-x86_64-macos` (Intel Mac: Accelerate's AVX2 path), `ssh-test-x86_64-macos-rosetta` (its SSE path) |
| both, Mach-O build | `rstd-check-macho` |
| the pre-13.3 fallbacks (`$NEWLAPACK` 32-bit, classic `cblas_*`) | recorded inputs: the binding against the 11.3, 12.3 and 13.3 SDKs' export lists (native tests); no macOS 11 or 12 machine is at hand |
| LAPACK `potrf`'s symbol spelling and special values | none yet (T114) |

## BLAS cores

A core is in its platform's tier. Kernel cores are `kernelCores` in
[`cpu/corespec.zig`](src/math/kernels/blas/cpu/corespec.zig); a core runs where detection
(`runnable`) says the rig's CPU can run it. Unlisted cores have no kernels: the portable ones
serve them.

### Untested in tiers 1 and 2

Compiled and symbol-checked, never run:

| Cores | Platforms | Why | Item |
| --- | --- | --- | --- |
| `skylakex`, `cooperlake`, `sapphirerapids` | every x86_64 one | AVX-512: the Zen 3 host lacks it, so WHPX cannot show it; TCG has none | T43, T44 |
| `bulldozer`, `piledriver`, `steamroller`, `excavator` | every x86_64 one | FMA4 and XOP: the host lacks them; TCG has none | T43 |
| every aarch64 core but the runner's | `aarch64-windows` | CI's `windows-11-arm` legs run the cores their CPU can; detection is tested with recorded inputs | T47 |

### x86_64 (tier 1)

| Core | Runs on |
| --- | --- |
| `prescott`, `atom`, `core2`, `penryn`, `dunnington`, `nehalem`, `opteron`, `opteron_sse3`, `barcelona`, `nano`, `sandybridge`, `bobcat`, `haswell`, `zen` | native Zen 3 with `-Dblas-cores=all`, `x86_64-whpx`, `x86_64`, `x86_64-uefi`, `linux-x86_64`; `ssh-test-x86_64-macos` (the Intel Mac, through `haswell`) |
| `bulldozer`, `piledriver`, `steamroller`, `excavator` | none (untested, above) |
| `skylakex`, `cooperlake`, `sapphirerapids` | none (untested, above) |

### aarch64 (tier 2)

| Core | Runs on |
| --- | --- |
| `armv8`, `cortexa53`, `cortexa55`, `cortexa57`, `cortexa72`, `cortexa73`, `cortexa76`, `cortexx1`, `falkor`, `thunderx`, `thunderx2t99`, `tsv110`, `emag8180`, `neoversen1`, `ft2000` | `aarch64`, `aarch64-sve`, `aarch64-uefi`, `linux-aarch64` |
| `thunderx3t110`, `ampere1` (need FCMA) | `aarch64-sve`, `linux-aarch64`; `aarch64` and `aarch64-uefi` under HVF |
| `a64fx`, `armv8sve`, `neoversev1`, `neoversen2`, `neoversev2`, `cortexa510`, `cortexa710`, `cortexx2` (SVE) | CI's `aarch64-linux-musl, sve` leg (qemu-user), `linux-aarch64-sve`; `aarch64-sve` (tier 3) |
| `armv9sme` (SVE2, SME) | CI's `aarch64-linux-musl, sve` leg (qemu-user), `linux-aarch64-sve`; `aarch64-sve` (tier 3) |
| `vortex` (Apple) | `ssh-test-aarch64-macos`; `aarch64`, `aarch64-uefi` and `linux-aarch64` under HVF |
| `vortexm4` (Apple, SME) | `ssh-test-aarch64-macos`; `aarch64` and `aarch64-uefi` under HVF |

### arm (tier 2)

| Core | Runs on |
| --- | --- |
| `armv6`, `armv7`, `cortexa9`, `cortexa15` | `arm`, `arm-softfp`, `linux-arm` (Cortex-A15, VFPv4 and NEON) |

RandyOS's PowerPC 750 is soft float, so it links no kernel object: the portable kernels serve it,
and every POWER core is tier 3.

### Tier 3 architectures

| Architecture | Core | Runs on |
| --- | --- | --- |
| x86 | `katmai`, `coppermine`, `northwood`, `banias`, `p5`, `p6`, `viac3`, `opteron` | `i386` (`pentium3`), `i386-max` |
| x86 | `prescott`, `yonah`, `atom`, `core2`, `penryn`, `dunnington`, `nehalem`, `opteron_sse3`, `barcelona`, `nano`, `bobcat` | `i386-max` |
| x86 | `athlon` (3DNow!) | unsettled: whether TCG `max` shows 3DNow! is unchecked; the Zen 3 host has none |
| powerpc64, powerpc64le | `power3`, `power4`, `power5`, `power6`, `ppc970`, `cell`, `power8`, `power9` | `powerpc64` (every 64-bit core), `powerpc64le`, `linux-powerpc64le` |
| powerpc64, powerpc64le | `power10` | `powerpc64`, `powerpc64le-power10`, `linux-powerpc64le-power10` |
| powerpc | `power3`, `power4`, `power5`, `power6`, `ppc970`, `cell` | `powerpc-pseries` (32-bit code on POWER10) |
| powerpc | `power8`, `power9` | `powerpc-power9` (32-bit code on POWER9) |
| powerpc | `ppcg4` | `powerpc` (G4) |
| mips64 | `loongson3r3` | `mips64r2el`, `mips64r2` |
| mips64 | `loongson3r4` (MSA) | `mips64r2el`; its fused bits rest on the manuals, as QEMU rounds `madd.fmt` |
| mips64 | `i6400`, `i6500`, `p6600` | `mips64r6el`, `mips64r6` |
| mips | `p5600` | `mips32r5el`, `mips32r5` |
| loongarch64 | `la64_generic`, `la264`, `la464` | `loongarch64` |
| riscv64 | `riscv64_zvl128b` | `riscv64-vlen128` |
| riscv64 | `riscv64_zvl256b` | `riscv64-vlen256` |
| riscv64 | `x280` | `riscv64`, `linux-riscv64` (VLEN 512) |
| s390x | `z13`, `z14` | `s390x` |
| wasm32 | `wasm128_generic` | `wasm-test-simd128` |

### Tier 4 cores

| Core | Why | Item |
| --- | --- | --- |
| `ppc440`, `ppc440fp2` | no rig runs a 440; compile-checked only | T41 |
| `power10` in 32-bit POWER builds | the dev toolchain's LLVM cannot build 32-bit compiler-rt for it, so it does not link | T105 |
| every sparc64 object | sparc64 is tier 4; its rigs are out of the sweeps | T39 |
| every 32-bit SPARC object | lld cannot link 32-bit SPARC | T48 |
| `c910v` (RISC-V) | no `xtheadvector` in LLVM; no object | [`readme.md`](src/math/kernels/blas/readme.md#not-supported-by-zig) |
| POWER BG/P (double FPU) | LLVM cannot assemble its instructions; not converted | [`readme.md`](src/math/kernels/blas/readme.md#not-supported-by-zig) |

## Toolchains

| Toolchain | Used by | Checked by |
| --- | --- | --- |
| Zig 0.18.0-dev.1029+1ebdf2f8f (self-built) | rstd itself, RandyOS | every step above |
| Zig 0.17.1 | zfls (ReleaseSafe, its six targets) | `rstd-check-zig017`: builds and links zfls's targets, every core, without running them; `zig build -h` after build-script edits |
