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

### sparc64: f80 complex nrm2 off by about f64 precision (unconfirmed)

- **Seen:** on the sparc64 rig (`.fast`, software f80), "portable nrm2
  kernels at the extended type match a plain loop" fails for
  `Complex(f80)`: got 0.799018197057005788, want 0.799018197057005830, a
  relative difference of about 5e-17 against an allowance of 2.8e-17, as if
  one side were computed in f64. The same test passes on aarch64 (also
  software f80), freestanding and under Linux, and on x86.
- **Workaround:** none yet; the test stays failing on that rig.
- **Next:** find which side is wrong (kernel or reference), then whether it
  is `.fast` or the sparc backend (a smaller image that fits a debug build
  would tell). No issue text until then.

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
  dated installer build (`build.zig.zon`).
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
upstream sources" in `readme.md`: an `x86_64` pair on Windows, and the
`loongarch64` `lsx`/`lasx` kernels, mostly with strides other than 1. One
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

And for `x86_64`:

> **x86_64 `zdot_atom.S` and `zgemv_t_atom.S`: wrong argument handling on Windows**
>
> - `zdot_atom.S`: Windows passes the address of the 16-byte result first,
>   shifting the arguments; the kernel reads that address as `n`.
>   `zdot_sse2.S` already handles this.
> - `zgemv_t_atom.S`: reads the double `alpha_i`, passed on the stack on
>   Windows, with a single-precision load; its sibling kernels use the
>   double one.
