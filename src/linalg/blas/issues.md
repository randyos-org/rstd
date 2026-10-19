# Known upstream issues

Problems found in other projects while building and testing the BLAS port.
None is filed. For each: what we saw, how we work around it and where,
and the text of an issue we could file later. The toolchain throughout is
Zig 0.18.0-dev (self-built, `C:\scratch\git\zig-build`) unless an entry says
otherwise.

When an issue is fixed upstream, drop the workaround and its entry.

## Zig compiler

### f80 returned through a hidden pointer under `x86_64_sysv` on a Windows target

- **Seen:** a function declared `callconv(.{ .x86_64_sysv = .{} })` in an
  `x86_64-windows` build returns `f80` through a hidden pointer in `rdi`;
  SysV returns it in `st(0)`. Arguments are passed correctly. A call into
  the x87 `qdot` assembly crashed on it.
- **Workaround:** none needed in the library: OpenBLAS never builds its x87
  kernels for Windows, and `types.exprecisionOn` keeps `XDouble` at `f64`
  there. The f80 research harness used a naked thunk.
- **Issue text:**

  > **x86_64: f80 return value goes through memory with `callconv(.x86_64_sysv)` on x86_64-windows**
  >
  > On an `x86_64-windows` target, a function with
  > `callconv(.{ .x86_64_sysv = .{} })` returning `f80` is lowered with a
  > hidden result pointer in `rdi`. The System V ABI returns `long double`
  > in `st(0)` (class X87). The same declaration on `x86_64-linux` returns
  > in `st(0)`. It looks like the Win64 rule for `long double` (by memory)
  > is applied although the function uses the SysV convention. Arguments
  > are passed correctly; only the return differs. Calling a SysV assembly
  > routine that returns in `st(0)` then reads garbage, or crashes when the
  > callee writes through the "pointer". Reproduce: declare
  > `extern fn f() callconv(.{ .x86_64_sysv = .{} }) f80;`, call it, and
  > compare the emitted assembly for `x86_64-windows` and `x86_64-linux`.

### Comptime-known pointer to a naked function, cast to SysV and called, clobbers the caller

- **Seen:** in ReleaseFast and ReleaseSafe, `@ptrCast` of a comptime-known
  pointer to a `callconv(.naked)` function to a SysV function type, then
  called directly, corrupted the caller's callee-saved registers (it
  surfaced as a `@memcpy arguments alias` panic). Debug was fine, and so
  was the same call through a runtime pointer. The call looks emitted with
  the default (Win64) convention, so `rsi`/`rdi` were assumed preserved.
- **Workaround:** the vtable only calls kernels through runtime function
  pointers on their own ABI, which is not affected. Keep it that way: do not
  turn vtable calls into comptime-known calls of the naked kernels.
- **Issue text:**

  > **Calling a comptime-known naked function through `@ptrCast` uses the default calling convention**
  >
  > On x86_64-windows, given `fn k() callconv(.naked) void { asm volatile (...) }`
  > implementing a System V function, the call
  > `@as(*const fn (...) callconv(.{ .x86_64_sysv = .{} }) void, @ptrCast(&k))(args)`
  > in ReleaseFast or ReleaseSafe does not treat `rsi`/`rdi` as clobbered:
  > the caller keeps values there across the call, as the Win64 convention
  > allows, and they come back changed. Loading the same pointer through a
  > `volatile` (so it is not comptime-known) produces a correct call. Debug
  > builds are correct either way. It seems the direct call is lowered
  > against the callee's own convention (`.naked`, so the target default)
  > instead of the pointer type's.

### powerpc64: `.abiversion 2` emitted where the target's convention is ELFv1

- **Seen:** for big-endian `powerpc64`, where `Target.cCallingConvention`
  reports ELFv1, Zig still emits `.abiversion 2` and ELFv2 code.
- **Workaround:** `abi.zig` picks the converted kernels' `64be` variant for
  `powerpc64` regardless of what `cc` says, matching what Zig emits.
- **Issue text:**

  > **powerpc64 (big-endian) objects are ELFv2 although `cCallingConvention` says ELFv1**
  >
  > For `powerpc64-linux-gnu`, `std.Target.cCallingConvention()` returns
  > the ELFv1 convention, but the emitted objects carry `.abiversion 2` and
  > follow ELFv2 (no function descriptors, TOC set up the ELFv2 way). Code
  > that selects assembly by the reported convention then links ELFv1
  > routines into an ELFv2 object. Either the emitted ABI or the reported
  > convention should change so the two agree.

### powerpc64le: a two-float struct is returned in `r3`/`r4`, not `f1`/`f2`

- **Seen:** a function returning an `extern struct { re: f64, im: f64 }` is
  read from `r3`/`r4` by Zig. ELFv2 returns homogeneous floating-point
  aggregates in `f1`/`f2`, which is where OpenBLAS's `power/zdot.S` and
  `zdot_cell.S` put the result.
- **Workaround:** those two kernels are patched to also copy the result to
  `r3`/`r4` (a path the sources already have for the Fortran interface).
  See "Fixes to upstream sources" in `readme.md`.
- **Issue text:**

  > **powerpc64le: homogeneous float aggregates returned in GPRs instead of FPRs**
  >
  > On `powerpc64le-linux`, a C-callconv function returning
  > `extern struct { a: f64, b: f64 }` is expected by Zig in `r3`/`r4`. The
  > ELFv2 ABI (section 2.2.4.2) returns homogeneous floating-point
  > aggregates of up to eight members in `f1` to `f8`. Calling C or assembly
  > that follows the ABI reads the wrong registers. Example: OpenBLAS's
  > `zdot` kernels for POWER return `double _Complex` in `f1`/`f2`.

### arm soft-float: a struct of floats is returned in r0 to r3, not through memory

- **Seen:** on `arm-freestanding-eabi` and `arm-linux-gnueabi`, a
  `callconv(.c)` function returning `extern struct { re: f32, im: f32 }` is
  lowered as returning in r0/r1, and one of two `f64` in r0 to r3. AAPCS
  returns a composite over 4 bytes in memory, through a pointer in r0, as
  clang does for the same C (`zig cc`). `src/codegen/arm/abi.zig`'s
  `classifyType` makes any struct of up to four floats `.byval`, the VFP
  rule, whatever the convention. OpenBLAS's arm `cdot` and `zdot` (`sf`)
  store through r0, so on the `arm-softfp` rig they took `n` as the result's
  address, read `incx` as `x`, and faulted.
- **Workaround:** the generator patches the `sf` `cdot_vfp.S` and
  `zdot_vfp.S` to take no pointer and return in r0 to r3
  (`tools/blas_asm/arch/arm.zig`).
- **Issue text:**

  > **arm: float structs are returned in core registers under the soft-float ABI**
  >
  > With `callconv(.c)` on `arm-linux-gnueabi` (or `-eabi`), a function
  > returning `extern struct { a: f32, b: f32 }` returns in r0/r1, and
  > `extern struct { a: f64, b: f64 }` in r0 to r3. The base AAPCS returns
  > a composite type larger than 4 bytes in memory, at an address passed in
  > r0; clang does so for the equivalent C struct, so calls between Zig and
  > C disagree. `classifyType` in `src/codegen/arm/abi.zig` returns
  > `.byval` for a struct of up to four floats (the VFP variant's
  > homogeneous aggregate rule) without checking for the soft-float
  > convention (`arm_aapcs` rather than `arm_aapcs_vfp`).

### Self-hosted x86_64 backend leaves a value on the x87 stack converting to f80

- **Seen:** with the self-hosted x86_64 backend (the debug default), converting
  `u128` or `f128` to `f80` (a compiler-rt libcall) leaves an entry on the x87
  register stack. Later x87 arithmetic overflows the stack and produces NaN.
- **Workaround:** the tests never convert `u128` or `f128` to `f80`:
  `blas1/check.zig` builds random `f80` values from bits, and compares `f128`
  references at `f128` (`expectCloseRef`).
- **Issue text:**

  > **x86_64 backend: `@floatFromInt(u128)` / `@floatCast(f128)` to `f80` unbalances the x87 stack**
  >
  > With the self-hosted x86_64 backend, `@as(f80, @floatFromInt(x))` for a
  > `u128` `x`, or `@as(f80, @floatCast(y))` for an `f128` `y`, returns the
  > right value but leaves one extra entry on the x87 register stack. After
  > several such conversions, x87 arithmetic in the same function produces
  > NaN (stack overflow, invalid operation). The LLVM backend is correct.
  > Reproduce: convert a few `u128` values to `f80` in a loop, add them, and
  > compare against the LLVM backend; `fstenv` after the call shows the tag
  > word with a register in use.

### LLVM `-Oz` miscompile on sparc64

- **Seen:** since Zig 0.18.0-dev.250, the sparc64 rig built with `.small`
  jumps into unrelated code in the gemm parameter test.
- **Workaround:** the sparc64 rig builds with `.fast` (`build/qemutest.zig`).
- **Issue text:** not ready; needs a reduced case first. It may belong to
  LLVM rather than Zig.

### sparc64: compiler-rt's `_Qp_sqrt` takes the square root in f64

- **Seen:** on the sparc64 rig, `@sqrt` of a runtime `f128` (lowered to
  `_Qp_sqrt`) gives only f64's precision: `@sqrt(2)` came back as
  `0x3fff6a09e667f3bcd000000000000000`. "portable nrm2 kernels at the
  extended type match a plain loop" failed for `Complex(f80)` on it: the
  kernel was right, the `f128` reference was not. `lib/compiler_rt/sqrt.zig`
  has `c.* = sqrt(@floatCast(a.*))`, where `sqrt` is the f64 routine; still
  so on Zig master (added in 2c5f54c898, 2025-04-28).
- **Workaround:** test references take the root with `check.sqrtW`, which on
  sparc64 refines `@sqrt` with two Newton steps in `f128`, scaling by even
  powers of two into f64's range first (a denormal's square is 0 in f64, and
  the steps from 0 give infinity). The library's own
  `f128` square roots (the portable `nrm2` and `rotg` at `f128`) stay at f64
  precision there until Zig is fixed.
- **Issue text:**

  > **compiler-rt: sparc64 `_Qp_sqrt` computes the square root in f64**
  >
  > `lib/compiler_rt/sqrt.zig` defines
  > `fn _Qp_sqrt(c: *f128, a: *f128) callconv(.c) void { c.* = sqrt(@floatCast(a.*)); }`.
  > `sqrt` there is the f64 routine, so `a` is narrowed to f64 and the root
  > widened back: `@sqrt` of a runtime `f128` on sparc64 has 53 bits of
  > precision. For example `@sqrt(@as(f128, 2))` with a runtime 2 returns
  > `0x3fff6a09e667f3bcd000000000000000` instead of
  > `0x3fff6a09e667f3bcc908b2fb1366ea95`. It should call `sqrt_f128`, as
  > `sqrtq` (`_Q_sqrt` on sparc32) does.

### powerpc64le: LLVM's `-O0` compares a VSX register with `fcmpu`

- **Seen:** on the powerpc64le rigs (pwr9 and pwr10, Debug), every gemm test
  failed at f64 with "padding written", though the padding held the sentinel:
  `got != sentinel` came out true for equal values. The f64 was loaded with
  `lxsdx 34, ...` (VSR34, which is v2) and compared with `fcmpu 0, 2, 0`,
  which reads f2 (VSR2). The same compare in a release build is correct. It
  looks like PowerPC FastISel's compare using `FCMPUD` on a register of the
  VSX class `VSFRC` without constraining it to `F8RC`.
- **Workaround:** the powerpc64le rigs build with `.safe`
  (`build/qemutest.zig`), which does not use FastISel.
- **Issue text:**

  > **[PowerPC] FastISel emits `fcmpu` on a VSX-class register at -O0**
  >
  > On powerpc64le with VSX (seen with `-mcpu=pwr9` and `pwr10`) at -O0, an
  > f64 loaded in one block and compared in another was loaded by FastISel
  > with `lxsdx` into a `VSFRC` virtual register, which register allocation
  > put in VSR34 (v2). The compare is `fcmpu cr0, 2, 0`: `FCMPUD` takes
  > `F8RC` operands, so the encoded `2` reads f2 (VSR2), not v2, and the
  > comparison uses a stale value. `PPCFastISel::PPCEmitCmp` appears to pass
  > the source registers to `FCMPUD` without constraining their class (or
  > using `XSCMPUDP` for VSX registers). Seen in a large Zig-generated
  > function; the same IR at -O1 and above is correct.

### lld cannot link 32-bit SPARC

- **Seen:** `zig build rstd-test -Dtarget=sparc-linux-gnu` compiles every
  object, then fails to link: `ld.lld: unsupported e_machine value: 2`
  (EM_SPARC). lld supports only SPARC V9 (EM_SPARCV9). Still so on
  0.18.0-dev.267 (2026-10-17), for a one-function freestanding image too;
  `-fno-lld` instead crashes Zig's own ELF linker (`invalid enum value` in
  `link/Elf/Object.zig`'s `parseCommon`, which reads every object header as
  `Elf64_Ehdr`), and that linker ignores linker scripts anyway.
- **Workaround:** none; 32-bit SPARC is compile-checked only (the kernel
  objects build). Its rig, `sparc` (OpenBIOS on `SS-5`), is a stub in
  `build/qemutest.zig`.
- **Issue text:** an lld feature request rather than a Zig bug, and probably
  known upstream: check LLVM's issue tracker for "lld sparc32" before filing.

## Zig standard library

### `SafeAllocator` asks for the CPU count in single-threaded builds

- **Seen:** `std.heap.SafeAllocator.acquireThread` calls
  `std.Thread.getCpuCount()` unconditionally. On freestanding, which has no
  `Thread` implementation, that is a compile error even with
  `-fsingle-threaded`, so `std.testing.allocator` cannot be used there.
- **Workaround:** the freestanding test runner (`test/qemu/runner.zig`)
  overrides `std_os_options.Thread` with a one-CPU implementation.
- **Issue text:**

  > **`std.heap.SafeAllocator` requires `std.Thread` even when single-threaded**
  >
  > `SafeAllocator.acquireThread` calls `std.Thread.getCpuCount()` to size
  > its thread table, whatever `builtin.single_threaded` says. On a target
  > without a `Thread` implementation (freestanding), this makes
  > `std.testing.allocator` a compile error
  > (`Unsupported operating system freestanding`). In single-threaded
  > builds there is only one thread, so it could use one table entry and
  > skip the CPU count.

### The `Thread` override's documented `getCpuCount` signature does not compile

- **Seen:** `std.Options.OperatingSystem.Thread`'s doc comment gives
  `pub fn getCpuCount() anyerror!usize;`, but `std.Thread.getCpuCount` does
  `return try Impl.getCpuCount();` into `CpuCountError!usize`, so an override
  written from the docs fails with "global error set cannot cast into a
  smaller set".
- **Workaround:** our override returns `std.Thread.CpuCountError!usize`.
- **Issue text:**

  > **Docs: `std_os_options.Thread.getCpuCount` must return `CpuCountError!usize`, not `anyerror!usize`**
  >
  > The doc comment on `std.Options.OperatingSystem.Thread` lists
  > `pub fn getCpuCount() anyerror!usize;`. `std.Thread.getCpuCount` returns
  > `CpuCountError!usize` and `try`s the override, so an override with the
  > documented signature does not compile. Either the doc should say
  > `std.Thread.CpuCountError!usize`, or `std.Thread.getCpuCount` should
  > accept any error set from the override.

## Zig package manager and build system

### `zig fetch` needs an archive and a `Content-Type` header

- **Seen:** `zig fetch` of a single file fails with "unrecognized
  'Content-Type' header: 'text/plain'", and of a `.tar.gz` from
  snapshot.debian.org (which sends no `Content-Type`) with "missing
  'Content-Type' header", although the URL ends in `.tar.gz`.
- **Workaround:** the Linux rigs' kernels come from Debian's netboot
  `.tar.gz` on deb.debian.org (which sends `application/x-gzip`), pinned to a
  dated installer build (`build.zig.zon`). s390x and mips64el have no Linux
  rig: Debian ships their installer kernels only as single files, sent with
  no `Content-Type` (13's s390x `generic/kernel.debian`, 12's mips64el
  `malta/netboot/vmlinuz-*-5kc-malta`). Both boot the rig's initramfs by
  hand: s390x on `s390-ccw-virtio` with `console=ttysclp0`, mips64el on
  `malta` (`-cpu MIPS64R2-generic`, 256 MB; 1 GB hangs) with `ttyS0`.
- **Issue text:**

  > **`zig fetch`: fall back to the URL's extension when `Content-Type` is missing, and allow single files**
  >
  > Fetching `https://snapshot.debian.org/.../netboot.tar.gz` fails with
  > `missing 'Content-Type' header`: the server sends none, though the path
  > says `.tar.gz`. Using the extension when the header is absent (or
  > `application/octet-stream`) would make such mirrors usable. Separately,
  > a package that is a single file (a kernel image, a firmware blob)
  > cannot be fetched at all; storing it as a one-file package would avoid
  > wrapping it in an archive.

### Run steps have no timeout for system commands

- **Seen:** `--test-timeout` applies only to unit tests. A qemu run whose
  firmware halted (OpenBIOS "Stopping execution") never ended, and with it
  the build.
- **Workaround:** `tools/qemu_watch.zig` runs each rig's qemu and ends it on
  firmware halt text or an hour without output.
- **Issue text:**

  > **Build system: a timeout for `Step.Run` system commands**
  >
  > `--test-timeout` covers unit tests run through the test protocol, but a
  > `Step.Run` of an external program has no time limit. An emulator or
  > other long-running tool that hangs blocks the build forever. A per-step
  > `timeout` (and ideally an inactivity limit on its output) would let the
  > build fail instead.

### A run step that checks stdout shows nothing until the program exits

- **Seen:** with `expectStdOutMatch`, a run step captures all of stdout and
  shows it only when the program exits, so a long run gives no progress.
- **Workaround:** the rigs have qemu also log the console to
  `.zig-cache/qemu-test-<name>.log` as it is written.
- **Issue text:**

  > **Build system: stream a run step's output while still checking it**
  >
  > A `Step.Run` with `expectStdOutMatch` (or any `stdio = .check`) buffers
  > the child's whole stdout and prints it only at the end. For long-running
  > steps there is no progress at all. An option to pass the output through
  > as it arrives while still matching it afterward would help.

### Build-script compile errors are hidden by a cached runner

- **Seen:** after a compile error in `build.zig` or a file it imports,
  `zig build` silently runs the previously cached build runner; `zig build -h`
  lacks the new option and prints no error. `--cache-dir` with a fresh
  directory shows the error.
- **Workaround:** check build-script changes with
  `zig build -h --cache-dir .zig-cache-fresh`.
- **Issue text:**

  > **`zig build` runs a stale cached build runner when the build script fails to compile**
  >
  > With the configurer/maker split, introducing a compile error in
  > `build.zig` and running `zig build` (or `zig build -h`) reports nothing
  > and behaves as the old script did. The error appears only with a fresh
  > `--cache-dir`. The build should fail with the compile error.

### A link interrupted midway leaves a cache entry without its executable

- **Seen:** a `zig build` killed while linking left `.zig-cache/o/<digest>/`
  with only `<name>_zcu.obj` and a manifest in `.zig-cache/h/`. Later builds
  treat it as cached, and the run step fails with `FileNotFound` on the exe.
- **Workaround:** move that manifest and `o/` directory out of the cache.
- **Issue text:**

  > **Cache hit does not check that the emitted files exist**
  >
  > Killing `zig build` during the link of an executable leaves a manifest
  > whose `o/` directory holds only the object file. Every later build takes
  > the cache hit and fails to spawn the missing executable
  > (`FileNotFound`). The manifest should be written after the outputs, or
  > a hit should verify them.

## translate-c (codeberg.org/ziglang/translate-c)

### `build.zig` treats `standardOptimizeOption` as non-optional

- **Seen:** at 9b5eb01 (its latest as of 2026-10-12), translate-c's
  `build.zig` uses `b.standardOptimizeOption(...)` as an `Optimize`. On
  current Zig master it returns `?Optimize` (aa8c774e0f), so the package
  does not build.
- **Workaround:** rstd's `translate_c` dependency is commented out and
  `src/buildutils/renderdoc.zig` uses the deprecated built-in `addTranslateC`.
- **Issue text:**

  > **build.zig does not compile on current Zig: `standardOptimizeOption` returns an optional**
  >
  > Since ziglang/zig aa8c774e0f, `std.Build.standardOptimizeOption`
  > returns `?std.builtin.OptimizeMode`. translate-c's `build.zig` passes
  > it where an `OptimizeMode` is required, so depending on the package
  > fails to configure. Using `orelse .Debug` (or the new default) fixes it.

## QEMU

### WHPX passes the host's AMD CPUID through to an Intel-only guest

- **Seen:** on an AMD Ryzen host, the macOS rig under `-accel whpx` stops
  right after `EXITBS:START`, whatever `-cpu ...,vendor=GenuineIntel` says.
  Under TCG it boots.
- **Workaround:** the macOS rig uses TCG.
- **Issue text:** none: QEMU's WHPX CPUID-trapping patches (leaf 0 vendor)
  were already in review in 2026. Check their status before filing.

### x86: TCG has no AVX-512, FMA4 or XOP

- **Seen:** QEMU 11.0.1's TCG warns "TCG doesn't support requested
  feature" for `avx512f`, `avx512dq`, `avx512cd`, `avx512bw` and `avx512vl`
  (`-cpu Skylake-Server-v4`) and for `fma4` and `xop` (`-cpu Opteron_G5`),
  and `-cpu max` has none of them. So the skylakex, cooperlake and
  sapphirerapids kernels, and bulldozer, piledriver, steamroller and
  excavator, never run on the TCG rigs.
- **Workaround:** `x86_64-whpx` runs the tests on the host CPU under WHPX,
  which runs them where the host has the extensions. This host (Ryzen 7
  7730U, Zen 3) has neither, so they are untested here. Under WHPX, `-cpu
  host` and `-cpu max` stop in SeaBIOS ("Unexpected VP exit code 4"); named
  models boot, keeping only the features the host has.
- **Issue text:** none; long-standing missing features, not bugs.

### x86: an x87 80-bit load or store across a page boundary misplaces bytes

- **Seen:** QEMU 11.0.1 under TCG, 32-bit Linux userspace (Debian 12's i386
  kernel, or Debian 13's amd64 kernel running an i386 program). When an
  `fldt` or `fstpt` operand crosses a page, the bytes of the 8-byte
  significand that lie on the second page are read from (or written to)
  the wrong place: a load returns stale bytes, a store leaves the guest's
  bytes unchanged and writes elsewhere. The 2-byte exponent, wholly on the
  second page, is right, and so is a `movq` across the same boundary. The
  same program under WHPX is right. i386's `f80` is 12 bytes, 4-aligned, so
  arrays of it cross pages; the rstd suite failed its `c_longdouble`
  checks and then crashed in `SafeAllocator` (a misplaced store). x86_64's
  16-byte `f80` never crosses.
- **Workaround:** `linux-i386` is a stub (`build/qemutest.zig`).
- **Issue text:**

  > **target/i386: fldt/fstpt across a page boundary access the wrong bytes under TCG**
  >
  > In 32-bit protected mode with paging, an `fldt` whose 10-byte operand
  > starts at page offset 0xff9 to 0xfff loads the significand bytes past the
  > boundary from somewhere other than the guest's next page (they come back
  > stale), while the exponent bytes are right; `fstpt` likewise fails to
  > write them. A `movq` of the same 8 bytes is right. It looks like the
  > 8-byte significand access in the x87 helpers treats the range as lying
  > in the first page (an `X86Access` pointer check on the start rather than
  > the end). Reproduce: in an i386 Linux guest, store 10 known bytes at a
  > stack address ending in 0xffc, `fldt` them and `fstpt` them elsewhere,
  > and compare.

### sparc64: `sll` drops the upper 32 bits of a 64-bit register

- **Seen:** under `qemu-system-sparc64` (sun4u), OpenBLAS's sparc kernels
  fault on any negative stride. They scale the 64-bit stride with
  `sll INCX, BASE_SHIFT, INCX`; the result under qemu is the low 32 bits,
  zero-extended, so `-2` became `0xfffffff8` and the next load landed 4 GB
  away. `sllx` instead fixes it (checked on `scopy`).
- **Workaround:** the generator writes the sparc `64` variant's `sll` as
  `sllx` (`tools/blas_asm/arch/sparc.zig`), which is right whichever reading
  of `sll` is correct.
- **Open question:** the SPARC V9 manual describes `sll` as shifting all 64
  bits, which would make this a qemu bug; if real UltraSPARC hardware
  zero-extends instead, it is an OpenBLAS bug (its kernels would then fail on
  hardware too). Check on hardware, or against the manual's exact text,
  before filing either.
- **Issue text (qemu, if hardware shifts 64 bits):**

  > **target/sparc: SLL on sparc64 zero-extends a 32-bit result instead of shifting all 64 bits**
  >
  > On `qemu-system-sparc64`, `sll %i2, 2, %i2` with `%i2 = -2` leaves
  > `0x00000000fffffff8` in `%i2`. The SPARC V9 manual (A.49) says SLL
  > shifts all 64 bits of `r[rs1]` (only the count is limited to 5 bits),
  > which gives `0xfffffffffffffff8`. Code that scales a negative 64-bit
  > index with `sll` (OpenBLAS's sparc kernels do) then computes an address
  > 4 GB away.
- **Issue text (OpenBLAS, if hardware zero-extends):**

  > **kernel/sparc: strides scaled with `sll` break for negative increments on 64-bit**
  >
  > The sparc level-1 kernels (for example `copy.S`, `nrm2.S`) scale
  > `INCX`/`INCY` with `sll INCX, BASE_SHIFT, INCX`. With 64-bit `BLASLONG`
  > and a negative increment, `sll` produces a zero-extended 32-bit value,
  > so the kernel steps forward by about 4 GB per element. Using `sllx`
  > under `__64BIT__` fixes it.

## OpenBLAS

The kernel bugs we fixed in the converted sources are listed under "Fixes to
upstream sources" in `readme.md`: `x86_64` argument and flag handling
(mostly on Windows), a POWER pair on ELFv2, and the `loongarch64`
`lsx`/`lasx` kernels, mostly with strides other than 1. One
issue covering the `loongarch64` ones could read:

> **loongarch64 LSX/LASX level-1 kernels: wrong results with non-unit strides and some edge cases**
>
> Testing the `kernel/loongarch64/*_lsx.S` and `*_lasx.S` level-1 kernels
> against plain reference loops (under qemu-system-loongarch64) found these
> (each is a separate, small fix):
>
> - `sum_lsx`, `sum_lasx`, `csum_lsx`, `csum_lasx`: the double result is
>   returned with a single-precision move, leaving its upper half undefined;
>   `sum_lasx` also adds doubles as floats.
> - `casum_lasx`, `csum_lasx`: the single-precision final sum adds the last
>   block of eight twice; `casum` also adds lane 6 in place of lane 7.
> - `axpby_lsx`, `caxpby_lsx`: with `beta == 0` and `incy != 1`, the
>   leftover elements are written over the first ones.
> - `axpby_lasx`: double precision skips every other block of eight, and
>   with `beta == 0` never terminates.
> - `caxpby_lsx`, `caxpby_lasx`: with `alpha == beta == 0` and `incy != 1`,
>   only the real parts are zeroed (`lasx` single zeroes one twice).
> - `max_lsx`, `min_lsx`: with `incx != 1`, compare by magnitude, as `amax`
>   does.
> - `camin_lsx`, `camin_lasx`: with `incx != 1`, the minimum starts from
>   `|re(x0)|` rather than `|re(x0)| + |im(x0)|`.
> - `imax_lasx`, `imin_lasx`: with `incx != 1`, single precision below 8
>   elements compares against a first element never loaded, and double
>   advances the index twice as fast.
> - `iamin_lsx`, `iamin_lasx`: with `incx != 1`, double keeps the minimum of
>   the previous comparison, and `lasx` skips the reduction across lanes.
> - `icamin_lsx`, `icamin_lasx`: `lsx` below 4 elements starts the index at 4
>   rather than 1, and with `incx != 1` its single precision starts from a
>   value never computed; `lasx` single with `incx != 1` overwrites its index
>   step, and its reduction over lanes 4 to 7 mistakes its own partial
>   minimum for lane 4's value.
> - `snrm2_lsx`, `snrm2_lasx`, `cnrm2_lsx`, `cnrm2_lasx`: with `incx < 0`
>   (which `interface/nrm2.c` passes on), return 0: the scale comes from
>   `?amax_k` on the same arguments, and it returns 0 for `incx <= 0`.

And for `x86_64`:

> **x86_64 `zdot_atom.S` and `zgemv_t_atom.S`: wrong argument handling on Windows**
>
> - `zdot_atom.S`: Windows passes the address of the 16-byte result first,
>   shifting the arguments; the kernel reads that address as `n`.
>   `zdot_sse2.S` already handles this.
> - `zgemv_t_atom.S`: reads the double `alpha_i`, passed on the stack on
>   Windows, with a single-precision load; its sibling kernels use the
>   double one.

And the `scal` flag on `x86_64`:

> **x86_64 `scal_sse.S`, `scal_sse2.S`, `scal_atom.S` (Windows) and `zscal_atom.S`: the flag argument is lost, and a NaN alpha counts as zero**
>
> `interface/scal.c` and `zscal.c` pass `dummy2 = 1` so that a zero alpha
> multiplies, turning NaN and Inf into NaN.
>
> - `scal_sse.S`, `scal_sse2.S`, `scal_atom.S`: under `WINDOWS_ABI` they
>   load the flag with `movq 64(%rsp)`, which is `inc_y` (the eighth
>   argument); the flag, the tenth, is at `80(%rsp)`. The interface passes
>   `inc_y = 0`, so a zero alpha always stores zeros on Windows.
> - `zscal_atom.S` never reads the flag: a zero alpha always stores zeros.
>   Loading it as `scal_atom.S` does (`24(%rsp)`, `88(%rsp)` on Windows) and
>   branching to the multiply path when it is 1 fixes it.
> - `scal_sse.S` tests `alpha == 0` with `comiss` and `jne` alone; a NaN
>   alpha sets ZF too and takes the zero path (`scal_sse2.S` has
>   `jp .L100`). With flag 0, as the level-2 drivers scale `y` by `beta`, a
>   NaN `beta` gives zeros. `zscal_atom.S` likewise tests `alpha_i == 0`
>   without `jp`, so `(0, NaN)` counts as zero.

And for POWER:

> **power `scal.S` and `zscal.S`: the flag argument is read at ELFv1's offset on Linux ELFv2**
>
> `kernel/power/scal.S` and `zscal.S` load the last argument (`dummy2`,
> the flag that keeps NaN for a zero alpha) with `LDLONG FLAG, 120(SP)` /
> `128(SP)` on Linux, the ELFv1 parameter save area. Under ELFv2 (ppc64le,
> and big-endian ELFv2 such as musl) the slots are 16 bytes lower, as the
> FreeBSD branch already has (104/112), so the kernel reads past its
> arguments and a zero alpha either stores zeros or multiplies depending on
> stack contents (NaN where zeros belong, for example `cgbmv` with
> `beta == 0` on POWER6, PPC970 and CELL). 32-bit SysV passes the flag in
> the first stack word. `LDLONG FLAG, FRAMESLOT(1)(SP)` (`scal.S`) and
> `FRAMESLOT(2)(SP)` (`zscal.S`) under `__64BIT__`, `FRAMESLOT(0)(SP)`
> otherwise, give the right slot for each ABI.

And for the x87 extended kernels (not converted here; the portable `scal`
runs for `XDouble`, keeping NaN as the flag asks):

> **x86_64 x87 `scal.S` and `zscal.S` (`qscal`, `xscal`): a zero alpha ignores the flag**
>
> `interface/scal.c` and `zscal.c` pass the flag (`dummy2 = 1`) that asks a
> zero alpha to multiply, so NaN and Inf in `x` give NaN, as reference BLAS
> does and as the other `scal` kernels do. `kernel/x86_64/scal.S` and
> `zscal.S` never read it: a zero alpha stores zeros, so `qscal` and `xscal`
> turn NaN and Inf into 0.

And in the x87 gemv kernels (the portable gemv takes their order but reads
only the elements it uses):

> **x87 `qgemv_n.S`, `qgemv_t.S`, `xgemv_n.S`, `xgemv_t.S`: a read past the end of `x`**
>
> `kernel/x86_64/` and `kernel/x86/` load the next `x` element before the
> loop tests whether there is one. With `incx == 1`, where `x` is the
> caller's own array, the last block reads `x[m]` (`qgemv_t`, when `n >= 4`)
> or, for an even block length, `x[n]` (`qgemv_n`, when `m >= 4`) or
> `x[m]`/`x[n]` (`xgemv_t`, `xgemv_n`). The value is discarded, but the read
> faults where `x` ends at an unmapped page. Testing the count before the
> load, or loading at the top of each trip, avoids it.

### Special values some converted kernels treat otherwise

The edge-case tests (`blas1/edge.zig` and `blas2/edge.zig`, on every rig)
found these. Reference BLAS's `i?amax` replaces its pick only on a strict
`>`, so NaN is picked only as the first element, and then nothing replaces
it; LAPACK 3.10's `nrm2` is infinite for an infinity; `kernel/arm/axpby.c`
reads no operand whose scalar is zero.

Not fixed, so `blas1/edge.zig`'s `Quirk` names it and the tests check what
the kernels do:

- **Extremes, lane by lane.** `x86_64/amax_sse.S`, `amax_sse2.S`,
  `zamax_sse.S`, `zamax_sse2.S`, `iamax_sse.S`, `iamax_sse2.S`,
  `izamax_sse.S` and `izamax_sse2.S` (every x86_64 core, all eight
  routines), `power/amax.S` and its kin (every POWER core), the
  `loongarch64` extremes and `arm64/amax.S` (`amax`, every arm64 core)
  compare in vector lanes, so what a NaN does depends on where it falls: a
  NaN first may be passed over, a later one picked. `power/amax.S`, with a
  NaN last, returned an element smaller than the maximum;
  `loongarch64/amax_lsx.S` returned 0 for a lone NaN. Left as they are: the
  reference's order needs each SIMD reduction rewritten.

Not used: the `nrm2` kernels below are converted but not wired in, and the
portable `nrm2` runs on those cores (`readme.md`):

- **`nrm2` of an infinity is NaN.** `power/nrm2.S` and its kin (every POWER
  core), `mips64/dnrm2.S` and `znrm2.S` (Loongson), `loongarch64/dnrm2.S`
  and `znrm2.S` (generic), `znrm2_lsx.S` and `znrm2_lasx.S`, and
  `sparc/znrm2.S`: their scaled sums divide the infinity by itself.
- **`nrm2` of a NaN is 0.** `loongarch64/snrm2_lsx.S`, `dnrm2_lsx.S`,
  `cnrm2_lsx.S` and the `lasx` ones, for a lone NaN: the scale is `?amax_k`
  of the vector, which passes NaN over.
- **`nrm2` of denormals.** The same POWER, Loongson, `loongarch64` and
  `sparc` kernels: `dnrm2` gave 0 or a value 15% off, `znrm2` (and
  `loongarch64`'s `snrm2` and `cnrm2`) infinity, where LAPACK's is
  accurate; they square the denormals, or scale by a reciprocal that
  overflows. `power/dnrm2_ppc440.S` and `znrm2_ppc440.S` (ppcg4, ppc440)
  gave NaN: their square root goes through `frsqrte` of a denormal sum.

Fixed in the converted sources (`readme.md`). Issues could read, for the
x86_64 index kernels:

> **x86_64 `iamax_sse.S`, `iamax_sse2.S`, `izamax_sse.S`, `izamax_sse2.S`: `n + 1` for a NaN extreme**
>
> The second pass looks for the element equal to the extreme the first
> found. Its last step (`.L28`, `.L68`, `.L98` in `iamax_sse2.S`, and the
> like) takes the last element as the match unread, and counts it even when
> `M` is even and it was read already. A NaN extreme equals no element, so
> the index returned is `n + 1`, past the end of `x`. Testing `M` odd there
> keeps it in range.

And the x86_64 copies on OPTERON:

> **x86_64 `copy_sse.S`, `copy_sse2.S`, `zcopy_sse.S`, `zcopy_sse2.S`: OPTERON loses the sign of zero**
>
> Under `OPTERON`, `LOAD` is `xorps REG, REG; addps OFFSET(ADDR), REG`
> (`addpd` for double): `+0 + -0` is `+0`, so `-0` is copied as `+0`, and a
> signaling NaN comes out quiet. The plain `movaps`/`movups` the other
> cores use copies the bits.

And x86_64 `symv`:

> **x86_64 `symv_L_sse.S`: an infinity in `x` gives NaN in `y`**
>
> The last two rows of a block (`.L17`) zero the upper half of each column
> operand (`movlhps xtemp2, a1`) but not of `xtemp1`, which still holds the
> next two elements of `x`; the two columns at the end (`.L20`) multiply a
> `movsd` of the matrix by all four lanes of `atemp4`. An infinity in those
> lanes times the zeroed half is NaN, which the horizontal sum adds into
> `y` (`ssymv` lower, `n = 7`, `x[6]` infinite, gives NaN in `y[0]`).
> Zeroing the upper half of `xtemp1` and of `atemp4` fixes both.

And the x86 and arm minima:

> **x86 `amax.S`, `iamax.S`, `zamax.S`, `izamax.S` and arm `amax_vfp.S`, `iamax_vfp.S`: the minima pick NaN**
>
> The x86 minima move the new element in unless `fcomi` says it is above
> the pick (`fcmovnbe`, `cmovb`); unordered is not above, so a NaN replaces
> the pick and the next element replaces the NaN. The arm minima move on
> `lt`, which an unordered `vcmpe` sets too. The maxima keep NaN out. The
> fix, for x86: compare the pick with the new element (`fxch` around
> `fcomi`) and use the maxima's moves; for arm: `mi` in place of `lt`.

And ties, on Loongson and SPARC:

> **mips64 `iamax.S`, `imax.S`, `izamax.S` (and the minima) and sparc `iamax.S`, `imax.S`, `izamax.S`: of equal extremes, a later index**
>
> They keep a pick per lane, then reduce the lanes keeping the first
> lane's pick on a tie, whose index may be later than another lane's.
> Reference BLAS returns the first. Taking the smaller index of equal picks
> at each step of the reduction fixes it.

And loongarch64 `axpby`:

> **loongarch64 `axpby_lsx.S`, `axpby_lasx.S`, `caxpby_lsx.S`, `caxpby_lasx.S`: a zero scalar's operand is read past the blocks**
>
> The tail (`.L998`) computes `alpha*x + beta*y` whatever the scalars, so a
> NaN or infinity meets a zero scalar's 0; `kernel/arm/axpby.c` reads
> neither operand then. The blocked loops for a zero `alpha` (`beta`) do
> not advance `X` (`Y`) either, so the tail reads the wrong elements. A tail
> that tests the scalars, as the blocks do, fixes both. Two more paths read
> the operand of a zero scalar:
>
> - `axpby_lsx.S`, `axpby_lasx.S`, double: the zero they compare the
>   scalars with is made by `ffint.s.l`, whose upper half is undefined (NaN
>   under qemu), so the zero-scalar blocks are never taken; `caxpby_lsx.S`
>   uses `ffint.d.l` for double.
> - `caxpby_lsx.S`, `caxpby_lasx.S`, one stride 1 and the other not (`.L12`,
>   `.L21`): no zero-scalar blocks at all; branching to `.L22`, the strided
>   path, which has them, fixes it.

Ours, not a converted kernel:

- **Rank updates and a zero `x(j)`.** Reference `zher`, `zhpr`, `dsyr2` and
  the rest skip a column whose factor is zero (for rank 2, `x(j)` and `y(j)`
  both). OpenBLAS's `driver/level2/syr_k.c`, `spr_k.c`, `zsyr_k.c` and
  `zspr_k.c` do; `zher_k.c`, `zhpr_k.c` and the rank-2 drivers do not, so
  an infinity elsewhere in `x` meets the zero and gives NaN in that column.
  `banded.zig`'s drivers skip as the reference does. Left: `ger`, whose
  reference skips a zero `y(j)`; the per-core `ger` kernels do not.

And SPARC `swap`, found by the ported `utest/test_swap.c`:

> **sparc `swap.S`, `zswap.S`: a zero stride swaps a block as one element**
>
> The strided loop (`.LL51`) loads a block of elements (8 for `swap.S`, 4
> for `zswap.S`) from `x` and `y` before storing any, so with `incx` or
> `incy` zero the `n` swaps of one element act as one. `utest`'s
> `zswap_inc_0` (`n` = 4, both strides 0) gets `x` and `y` exchanged where
> four swaps leave them as they were; `sswap` and `dswap` go wrong from 8
> elements. Sending a zero stride to the one-at-a-time tail loop (`.LL56`)
> fixes both.
