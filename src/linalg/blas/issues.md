# Known upstream issues

Problems found in other projects while building and testing the BLAS port.
None is filed. For each: what we saw, how we work around it and where,
and the text of an issue we could file later. The toolchain throughout is
Zig 0.18.0-dev (self-built, `C:\scratch\git\zig-build`) unless an entry says
otherwise. The Zig entries were rechecked on 0.18.0-dev.751 (2026-10-26).
That build is upstream master `cb6f0c25f1` (2026-10-26) with the open
<https://codeberg.org/ziglang/zig/pulls/31235> merged, so its `dev` number
is not upstream. Zig source links below are at `cb6f0c25f1`.

When an issue is fixed upstream, drop the workaround and its entry.

## Zig compiler

### f80 returned through a hidden pointer under `x86_64_sysv` on a Windows target

- **Seen:** a function declared `callconv(.{ .x86_64_sysv = .{} })` in an
  `x86_64-windows` build returns `f80` through a hidden pointer in `rdi`;
  SysV returns it in `st(0)`. Arguments are passed correctly. A call into
  the x87 `qdot` assembly crashed on it. Reproduced 2026-10-29 on 0.18.0-dev.751 with the
  issue text's declaration (`lea rdi` before the call on Windows, `fstp` of `st(0)` on Linux).
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
  > is applied, although the function uses the SysV convention. Arguments
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
  was the same call through a runtime pointer (reproduced 2026-10-29 with the issue text's case:
  `g` compiles to `unreachable`). On 0.18.0-dev.751 the LLVM
  IR shows why: `call x86_64_sysvcc void @k(i64 ..)` against
  `define void @k()`. LLVM treats a call whose convention differs from
  the callee's as undefined, and now compiles the whole caller to
  `unreachable`, so it crashes.
- **Workaround:** the vtable only calls kernels through runtime function
  pointers on their own ABI, which is not affected. Keep it that way: do not
  turn vtable calls into comptime-known calls of the naked kernels.
- **Issue text:**

  > **Calling a comptime-known naked function through `@ptrCast` uses the default calling convention**
  >
  > On x86_64-windows, given `fn k() callconv(.naked) void { asm volatile (...) }`
  > implementing a System V function, the call
  > `@as(*const fn (usize, usize) callconv(.{ .x86_64_sysv = .{} }) void, @ptrCast(&k))(a, b)`
  > is emitted as a direct `call x86_64_sysvcc void @k(i64, i64)` to
  > `define void @k()`, which has the default convention. LLVM considers a
  > call with a mismatched convention undefined behavior: in ReleaseFast
  > or ReleaseSafe the calling function compiles to a lone `unreachable`
  > (earlier builds instead kept values in `rsi`/`rdi` across the call).
  > Loading the same pointer through a `volatile`, so it is not
  > comptime-known, produces a correct indirect call; Debug with the
  > self-hosted backend is correct either way. A direct call through a
  > cast pointer should be emitted as an indirect call, or with the
  > callee's own convention matching the pointer type's.

### powerpc64: `.abiversion 2` emitted where the target's convention is ELFv1

- **Seen:** for big-endian `powerpc64`, where `Target.cCallingConvention`
  reports ELFv1, Zig still emits `.abiversion 2` and ELFv2 code.
  `cCallingConvention` gives ELFv1 for PS3 or a GNU ABI
  ([`lib/std/Target.zig:3773`](https://codeberg.org/ziglang/zig/src/commit/cb6f0c25f176dbebf0c6801170c360c2c1198524/lib/std/Target.zig#L3773)),
  but LLVM is told `elfv2` for all but PS3
  ([`src/target.zig:692`](https://codeberg.org/ziglang/zig/src/commit/cb6f0c25f176dbebf0c6801170c360c2c1198524/src/target.zig#L692)).
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
  See [`design.md`](design.md#fixes-to-upstream-sources), "Fixes to
  upstream sources".
- **Issue text:**

  > **powerpc64le: homogeneous float aggregates returned in GPRs instead of FPRs**
  >
  > On `powerpc64le-linux`, a C-callconv function returning
  > `extern struct { a: f64, b: f64 }` is expected by Zig in `r3`/`r4`. The
  > ELFv2 ABI (revision 1.5, section 2.2.6) returns homogeneous floating-point
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

### ThreadSanitizer: faults at start on static musl; macOS needs an SDK

- **Seen:** `-fsanitize-thread` for `x86_64-linux-musl` (static, PIE) links,
  but rstd's test binary dies of `SIGSEGV` before its first test on
  `qemu-test-linux-x86_64`. For `x86_64-linux-gnu` it runs and reports races.
  For `aarch64-macos`, cross-compiled from Windows, libtsan does not build:
  `sys/timeb.h`, `util.h`, `asl.h`, `sys/syscall.h`,
  `CoreFoundation/CFBase.h` and `objc/objc-sync.h` are not found. 0.18.0-dev.751,
  2026-10-29.
- **Workaround:** `-Drstd-tsan` sanitizes `linux-gnu` tests only; the x86_64 and
  aarch64 Linux rigs then build for glibc and take it, with its loader, from Nix
  (`tools/nix/blas.sh`).
- **Issue text:**

  > **-fsanitize-thread: static musl executables crash at startup**
  >
  > A test binary built with `-fsanitize-thread` for `x86_64-linux-musl` links
  > as a static PIE and dies of SIGSEGV before its first test; built for
  > `x86_64-linux-gnu` it runs under TSan. If static TSan is unsupported, the
  > compiler should reject it rather than link it. (Not yet reduced: check a
  > two-thread race program first.)

### Self-hosted x86_64 backend stores a two-lane `f32` vector as 16 bytes

- **Seen:** with the self-hosted x86_64 backend, `@sizeOf(@Vector(2, f32))`
  is 16 (8 with LLVM): `Type.abiAlignment` gives every `stage2_x86_64`
  vector at least 16, and `Type.abiSize` rounds up to it. A store through a `*align(1) @Vector(2, f32)` cast
  from an `[*]f32`, as `util.store` does, writes 16 bytes, so the two
  `f32` after the vector are overwritten. Still so on 0.18.0-dev.751.
- **Workaround:** the complex gemm tile assigns whole arrays
  (`p[0..lanes].* = ...`), not `util.store` (`blas3/gemm/kernel.zig`).
- **Issue text:**

  > **x86_64 backend: `@Vector(2, f32)` is 16 bytes, so a store through a pointer to it overwrites the next 8**
  >
  > With `-fno-llvm` on x86_64, `@sizeOf(@Vector(2, f32))` is 16; the LLVM
  > backend says 8, the vector's natural size. `Type.abiAlignment` gives every
  > `stage2_x86_64` vector a 16-byte minimum, so a store through
  > `*align(4) @Vector(2, f32)` writes 16 bytes. Code that views part of an
  > `[N]f32` as a vector (a common way to load and store lanes of an
  > array) silently corrupts the elements after it. Reproduce:
  >
  > ```zig
  > var buf: [4]f32 = .{ 9, 9, 9, 9 };
  > const p: *align(4) @Vector(2, f32) = @ptrCast(&buf);
  > p.* = .{ 1, 2 };
  > // -fllvm: { 1, 2, 9, 9 }; -fno-llvm: buf[2] and buf[3] overwritten
  > ```
  >
  > Coercing to an array (`buf[0..2].* = v`) stores 8 bytes. Vectors have
  > no guaranteed memory layout, but the size differing between backends,
  > and padding written past a vector that `[2]f32` holds exactly, are
  > surprising; giving the vector its natural size, as LLVM does, would
  > avoid both.

### Self-hosted x86_64 backend: no `(%reg)` operand in inline assembly

- **Seen:** `asm volatile ("stmxcsr (%[p])" :: [p] "r" (&v) : .{ .memory = true })`
  fails with "invalid memory operand: '(%[p])'" under `-fno-llvm`; so do
  `0(%[p])` and `(%rdx)`. Zig 0.17.0 and 0.18.0-dev.751. Found by
  `rstd-check-zig017`: 0.17.0 builds debug x86_64-linux with this backend.
- **Workaround:** `util.zig`'s `FpState` asm takes `m` operands on this
  backend (`fp_mem_ops`), the `(%[p])` form elsewhere (the next entry).
- **Issue text:**

  > **x86_64 backend: inline asm rejects register-indirect memory operands**
  >
  > `asm volatile ("ldmxcsr (%[p])" : : [p] "r" (&v) : .{ .memory = true });`
  > builds with `-fllvm` and fails with `-fno-llvm`: "invalid memory
  > operand: '(%[p])'". `0(%[p])` and a fixed `(%rdx)` fail the same way;
  > an `m` constraint (`ldmxcsr %[v]` with `[v] "m" (v)`) works.

### LLVM backend: an `m` input operand reads its address, not its value

- **Seen:** `asm volatile ("ldmxcsr %[v]" :: [v] "m" (v))` lowers to
  `call void asm "ldmxcsr ${0}", "m"(ptr %1)`, `%1` a spill of `v`. A direct
  `m` takes the pointer as the value, so LLVM spills the pointer and the
  instruction reads the address (`ldmxcsr dword ptr [rbp - 16]`, which holds
  `rbp - 4`); here a `#GP`. `*m` (indirect) is what the spill needs. `=m`
  outputs work. Zig 0.17.0 and 0.18.0-dev.751, x86_64.
- **Workaround:** `util.zig` passes a pointer in a register with LLVM
  (`(%[p])`), and `m` only on the self-hosted backend (the previous entry).
- **Issue text:**

  > **LLVM backend: inline asm `"m"` inputs pass a pointer to the value as the value**
  >
  > ```zig
  > export fn f(v: u32) void {
  >     asm volatile ("ldmxcsr %[v]" : : [v] "m" (v) : .{ .memory = true });
  > }
  > ```
  >
  > With `-fllvm -target x86_64-linux` the IR is
  > `asm sideeffect "ldmxcsr ${0}", "m,..."(ptr %1)`, where `%1` holds `v`.
  > With a direct `m` constraint LLVM treats `%1` as the operand's value and
  > spills it, so `ldmxcsr` loads the pointer, not `v`. Emitting `*m` for a
  > spilled input would load `v`. The self-hosted backend loads `v`.

### LLVM `-Oz` miscompile on sparc64

- **Seen:** since Zig 0.18.0-dev.250, the sparc64 rig built with `.small`
  jumps into unrelated code in the gemm parameter test. Still so on
  0.18.0-dev.751: "every core of the x86_64 family" traps with
  `Unhandled Exception 0x10`; the same image at `.fast` passes.
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
  so on 0.18.0-dev.751 (added in 2c5f54c898, 2025-04-28).
- **Workaround:** test references take the root with `check.sqrtW`, which on
  sparc64 refines `@sqrt` with two Newton steps in `f128`, scaling by even
  powers of two into f64's range first (a denormal's square is 0 in f64, and
  the steps from 0 give infinity). The library's own
  `f128` square roots (the portable `nrm2` and `rotg` at `f128`) stay at f64
  precision there until Zig is fixed; tests allow for them with
  `check.rootUlps`. "level 1 at any float type" failed at `f80` there while
  its reference took a plain `@sqrt` (T90).
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

### powerpc64: FastISel uses FPR-only instructions on VSX v-registers

- **Seen:** on the powerpc64le rigs (pwr9 and pwr10, Debug), every gemm test
  failed at f64 with "padding written", though the padding held the sentinel:
  `got != sentinel` came out true for equal values. The f64 was loaded with
  `lxsdx 34, ...` (VSR34, which is v2) and compared with `fcmpu 0, 2, 0`,
  which reads f2 (VSR2). The same compare in a release build is correct.
- **Status:** still there on 0.18.0-dev.751 (LLVM 23.1.2), rechecked
  2026-10-27. That compare is now right: FastISel compares f64 with
  `xscmpudp` on pwr9 and later, and on pwr8 constrains `fcmpu`'s operands to
  `F8RC`; the 0f0bcde^ check in a standalone function compiles clean. Its
  stores are not: an f64 that another block converts to f128 (`xscvdpqp`,
  native for freestanding pwr9 and later) gets the class `VFRC` (v0 to v31),
  and FastISel stores it with `STFD`, which reads f0 to f31.
  In the three freestanding rigs' Debug images, the MIR after instruction
  selection has FPR instructions on VSX classes only in
  `_test/blas1/edge.zig`'s two `contractCase` instances at f64 (the `say` of
  `got` and `want`), whose asm has `lxsd 2, 184(31)` then `stfd 2, 752(31)`.
  All four rigs pass in Debug, as that store only reaches a failure message;
  one like it could feed a check. Clean at `-O1`, with `-fast-isel=false`,
  and for pwr8 and `powerpc64le-linux` (f128 by library calls).
- **Workaround:** the freestanding pwr9 and pwr10 rigs build with `.safe`
  (`build/qemutest.zig`), which does not use FastISel; `linux-powerpc64le`
  (pwr8) builds in Debug.
- **Issue text:**

  > **[PowerPC] FastISel stores a `VFRC` register with `stfd` at -O0**
  >
  > With `-mcpu=pwr9` or `pwr10` at -O0, an f64 defined in one block and
  > converted to f128 in another gets the register class `VFRC` (from
  > `XSCVDPQP`'s operand). A third block's `store double` of it is selected
  > by FastISel as `STFD`, whose operand class is `F8RC`;
  > `-verify-machineinstrs` reports "Expected a F8RC register, but got a VFRC
  > register". Register allocation puts it in v3 (VSR35), and the encoded
  > `3` stores f3. `PPCFastISel::PPCEmitStore` picks `STFD` for any
  > class but exactly `VSFRC`
  > ([`PPCFastISel.cpp:637,651`](https://github.com/llvm/llvm-project/blob/85ac560262434c9ccfc0c183ec22d4138ed647fb/llvm/lib/Target/PowerPC/PPCFastISel.cpp#L637-L651),
  > 23.1.2); `PPCEmitCmp` once did the same with `FCMPUD`.
  > Zig source, `zig build-obj -O Debug -target powerpc64le-freestanding
  > -mcpu pwr9`:
  >
  > ```zig
  > extern fn sink(p: *const [2]f64) void;
  >
  > export fn f(get: *const fn () callconv(.c) f64, c: bool) void {
  >     const x = get();
  >     const y = get();
  >     if (c and @as(f128, x) != y) {
  >         const args: [2]f64 = .{ x, y };
  >         sink(&args);
  >     }
  > }
  > ```
  >
  > gives, for the `sink` block:
  >
  > ```asm
  > lxsd 2, 144(31)
  > lxsd 3, 128(31)
  > stfd 3, 184(31)
  > stfd 2, 192(31)
  > ```

### powerpc (32-bit): compiler-rt does not build for POWER4 and later

- **Seen:** building a 32-bit POWER10 rig (T11), on 0.18.0-dev.751
  (2026-10-29). Any `powerpc` build that bundles compiler-rt fails in LLVM
  for a server CPU, three ways by feature: with `64bit` (pwr4 to pwr6, g5,
  `pwr8-vsx`; g4 plus `64bit` too) the compiler segfaults in codegen;
  with VSX and `fpcvt` (pwr7, pwr8) `@intFromFloat` of an `f128` into an
  `i32` that goes back to float becomes `PPCISD::FCTIWZ` on the `f128` (an
  "invalid node" error, or "Cannot select" without `stfiwx`; compiler-rt's
  `rem_pio2l.zig:141`); with POWER9's vector features (native `f128`)
  "Cannot select" on `f128` `fcopysign` (`log_f128.zig:79`) or `v1i128`
  bitcasts (`fmod.zig:168`), or a "Potential legalization loop" assertion.
  Without those vector features but with ISA 3.0 (needed for ISA 3.1, which
  the power10 core needs), the tests' own code stops at "Invalid D-Form
  Pseudo-ops on Pre-P9 target" (`PPCInstrInfo.cpp:3252`, an assertion). The
  per-core kernel objects, which do not bundle compiler-rt, build.
  Zig 0.17.0 builds all of these, so it is new with the dev toolchain's LLVM.
- **Workaround:** `powerpc-power9` builds for POWER8 less `64bit` and
  `fpcvt`, on a POWER9: the power8 and power9 objects link and run, power10's
  cannot (T11). A 32-bit build has no workaround that keeps a POWER4 to 10
  model.
- **Issue text:**

  > **[PowerPC] 32-bit compiler-rt fails to compile for pwr4 and later**
  >
  > `zig build-exe -O Debug -target powerpc-freestanding-none -mcpu pwr6
  > e.zig -fcompiler-rt`, with `e.zig` an exported `_start` that loops,
  > segfaults in `ZigLLVMTargetMachineEmitToFile`; so do pwr4, pwr5, g5,
  > and `g4+64bit`; `pwr6-64bit` builds. With `-mcpu pwr8`, `pwr7` or
  > `pwr8-64bit`:
  >
  > ```zig
  > export fn g(x: f128) f64 {
  >     return @floatFromInt(@as(i32, @intFromFloat(x)));
  > }
  > ```
  >
  > fails in ReleaseFast with `LLVM ERROR: invalid node: operand #1 must
  > have type f64, but has type f128` on `PPCISD::STFIWX` of
  > `f128 = PPCISD::FCTIWZ`; `pwr8-vsx` builds it. With `pwr9-64bit` or
  > `pwr10-64bit-fpcvt`, compiler-rt fails with `Cannot select: f128 =
  > fcopysign` (`log_f128.zig:79`). Zig 0.17.0 builds every case.

### powerpc: LLVM vectorizes scalar `f32` loops into AltiVec, which flushes

- **Seen:** with the level-3 tiles scalar on the AltiVec cores, complex
  `trmm` and `trsm` still lost denormals there, in the per-core objects
  (built for speed) and in a ReleaseFast build for the G4; a Debug build
  and every real routine were right. LLVM's loop vectorizer had put the
  split tile's column loop in AltiVec lanes, which flush with `VSCR[NJ]`
  set. ARM's `isFPVectorizationPotentiallyUnsafe` keeps it off NEON's `f32`
  lanes; PowerPC has no such check. LLVM's `pwr3` to `pwr5` models have
  AltiVec, though those chips do not. Still so on 0.18.0-dev.751. That
  build's assembly (`-femit-asm`) also had AltiVec in the portable `saxpy`,
  both `gemm_beta`s and the scalar tile, real `trmm` and `trsm`, the half
  `gemv`, the `cdot` block of `dot/kernel.zig` (SLP, in ReleaseSmall too)
  and the tests' denormal operands, which came out zero; kept denormal,
  `saxpy` failed on the rig.
- **Workaround:** the AltiVec cores' objects are built without it
  (`cpu/_corespec.zig`'s `without`), and those loops pass each `f32` result
  they store or accumulate through `util.scalar`, an empty asm the
  vectorizers cannot widen. A new loop without it can still be vectorized:
  the `-femit-asm` of a `-mcpu=g4` build should have no `vaddfp`, `vmaddfp`
  or `vsubfp` outside `simd`, the other targets' tiles and test models.
  A ReleaseSafe G4 scan (2026-10-30, T42) found test code in AltiVec lanes
  only where the data is normal-range (`check.fill`'s [-2, 2), `blat3`'s
  tables, `threaded.zig`'s 2^±12): the `scal`, half-`gemv`, triangular,
  `geadd`, `matcopy`, level-3 driver and threaded models and
  `std.math.Complex`. Loops fed subnormals (`denormal.zig`'s `scaled`, the
  edge tests) stayed scalar, but for `denormal.zig`'s test that the lanes
  flush; a new one needs `util.scalar`.

### lld cannot link 32-bit SPARC

- **Seen:** `zig build rstd-test -Dtarget=sparc-linux-gnu` compiles every
  object, then fails to link: `ld.lld: unsupported e_machine value: 2`
  (EM_SPARC). lld supports only SPARC V9 (EM_SPARCV9). Still so on
  0.18.0-dev.751 (2026-10-26), for a one-function freestanding image too;
  `-fno-lld` instead crashes Zig's own ELF linker (`invalid enum value` in
  `link/Elf/Object.zig`'s `parseCommon`, which reads every object header as
  `Elf64_Ehdr`), and that linker ignores linker scripts anyway.
- **Workaround:** none; 32-bit SPARC is compile-checked only (the kernel
  objects build). Its rig, `sparc` (OpenBIOS on `SS-5`), is a stub in
  `build/qemutest.zig`.
- **Issue text:** none; an lld feature, not a Zig bug, and pending upstream:
  <https://github.com/llvm/llvm-project/pull/228200> (open) adds 32-bit SPARC
  to lld.

### Inline assembly errors exit 0

- **Seen:** on 0.18.0-dev.751 (LLVM 23.1.2), an instruction LLVM rejects
  in inline assembly (a bad mnemonic; on aarch64, also one the target's
  features lack) prints `error: <inline asm>: ...`, yet `zig build-obj` exits 0. The
  object is written without that line; the symbol is still there.
  `zig build` shows the step's text as "failure", then counts the step a
  success, and caches the object, so a rerun is silent. A per-core BLAS
  object can come out missing instructions. Zig 0.17.0 exits 1 (both rechecked 2026-10-29 with
  the issue text's `bogusinsn` case).
- **Workaround:** `build/blaskern.zig` (`checkedObject`) also emits each
  kernel object's LLVM bitcode, compiles it with `zig cc -c` in a run
  step, and links the object only through a copy that depends on that
  step. Clang fails on the same lines and its run step does not cache a
  failure. Objects built without LLVM, or with an explicit `ofmt`
  (`zig cc` ignores `-ofmt`), are not checked. The kernel rigs' test
  images (boot code and all) get the same check (`asmCheck`) before QEMU
  runs; `zig build qemu-check` builds and checks them all without QEMU.
  `-mcpu` is passed for module-level asm, which LLVM assembles with the
  default CPU, except on x86, whose assembler ignores CPU features.
  0.18.0-dev.751's `std.Target.x86` gives `skylake_avx512`,
  `icelake_client`, `icelake_server` and `goldmont_plus` those LLVM names
  (0.17.0: hyphens). Not a bug: LLVM 23 takes the underscore spellings as
  aliases, tuning included. The same IR compiles to the same code with
  either spelling of each, without a warning; a made-up name warns and
  changes the code. So the `skylakex` object is tuned; nothing to file.
  Emitting bitcode also prints a spurious `error: <inline asm>` for
  module-level asm that needs more than the default CPU's features; the
  object is unaffected.
- **Issue text:**

  > **LLVM inline assembly errors no longer fail the compilation**
  >
  > On 0.18.0-dev.751, `zig build-obj k.zig -target x86_64-linux -OReleaseFast`
  > with `fn f() callconv(.naked) void { asm volatile ("bogusinsn %rax\n ret\n"); }`
  > exported prints `error: <inline asm>:1:2: invalid instruction mnemonic
  > 'bogusinsn'` but exits 0 and writes an object in which `f` is a lone
  > `ret`. The result is cached, so the next build prints nothing. Under
  > `zig build` the step is reported as succeeding. Zig 0.17.0 exits 1, and
  > `zig cc -c` on the same module's bitcode fails as expected. The error
  > should fail the compilation; perhaps check `DiagnosticHandler::HasErrors`
  > after emitting.

### wasm bitcode gets its feature flags twice

- **Seen:** the asm check (above) of the `wasm128_generic` object fails in
  clang's verifier: "module flag identifiers must be unique" for each
  `wasm-feature-*` flag, then "Broken module found". Zig writes the
  bitcode (`-femit-llvm-bc`) after the codegen passes have run on the
  module (`src/zig_llvm.cpp`, `ZigLLVMTargetMachineEmitToFile`), and
  LLVM's WebAssembly codegen has added those flags by then; clang's
  codegen adds them again.
- **Workaround:** `build/blaskern.zig` (`asmCheck`) passes `-Xclang
  -disable-llvm-verifier` for wasm. Inline asm errors still fail it.
- **Issue text:**

  > **`-femit-llvm-bc` with `-femit-bin` writes post-codegen bitcode**
  >
  > `ZigLLVMTargetMachineEmitToFile` runs `codegen_pm` before
  > `WriteBitcodeToFile`, so the bitcode holds codegen's changes to the
  > module. For wasm32 that includes the `wasm-feature-*` module flags;
  > `zig cc -c` on that bitcode adds them again and fails verification
  > ("module flag identifiers must be unique"). Writing the bitcode before
  > `codegen_pm.run` gives the optimized IR that the object came from.

### MIPS64r6: LLVM cannot select between `-x` and `x` in a loop

- **Seen:** since a594fa7 (`trsm` divides by the diagonal), the
  `mips64r6` and `mips64r6el` rigs stopped building: the i6400 kernel
  object (`mips64r6`, no MSA) hit `Assertion failed: Opc && "Cannot copy
  registers"` (`MipsSEInstrInfo.cpp:261`, LLVM 23.1.2). The `ztrsm`
  `solveRight` conjugated with a run-time `conj` (`if (c) -x else x` on an
  `f64`). R6's `sel.d` takes its condition in an FPR; MachineLICM hoisted
  the condition's copy out of the loop as a `COPY` from a 32-bit GPR to a
  64-bit FPR, which `copyPhysReg` has no instruction for. Before a594fa7
  the same select compiled, by luck of the surrounding code. Where the
  copy has to be spilled instead, LLVM asserts elsewhere (next entry).
- **Workaround:** `conj` is comptime in `blas3/trsm/kernel.zig`'s `run`,
  `solveLeft`, `solveRight` and `conjIf`, so no select is left. A sign-bit
  `xor` would also do; `x * (if (c) -1 else 1)` would not: LLVM folds it
  back into the select.
- **Upstream:** reported as <https://github.com/llvm/llvm-project/issues/228690>
  (another select pair); fixed on main by
  <https://github.com/llvm/llvm-project/pull/226838> (2026-10-12, not in 23.1.x),
  which drops `FGR64CC`. Retest after the next LLVM before filing.
- **Issue text:**

  > **[Mips] R6: "Cannot copy registers" for a loop-invariant `select` of
  > `fneg x` and `x`**
  >
  > `llc -mtriple=mipsisa64r6el -mcpu=mips64r6 -O2` (here `zig cc -c -O2`)
  > on
  >
  > ```llvm
  > define void @f(i1 zeroext %c, i64 %n, ptr %p) {
  > entry:
  >   br label %loop
  > loop:
  >   %i = phi i64 [ 0, %entry ], [ %i.next, %loop ]
  >   %a = getelementptr inbounds double, ptr %p, i64 %i
  >   %x = load double, ptr %a, align 8
  >   %nx = fneg double %x
  >   %r = select i1 %c, double %nx, double %x
  >   store double %r, ptr %a, align 8
  >   %i.next = add i64 %i, 1
  >   %done = icmp eq i64 %i.next, %n
  >   br i1 %done, label %exit, label %loop
  > exit:
  >   ret void
  > }
  > ```
  >
  > asserts `Opc && "Cannot copy registers"` in
  > `MipsSEInstrInfo::copyPhysReg`, from `ExpandPostRAPseudos`. Before that
  > pass the entry block has `$d0_64 = COPY killed $at`: the `SEL_D`
  > condition, hoisted out of the loop, copied from a GPR32 into an FGR64.
  > Big-endian `mips64r6` too; `-O0` compiles. A select of two loaded
  > values does not trigger it.

### MIPS64r6: LLVM cannot spill a select's condition

- **Seen:** with MSA (`-mcpu=mips64r6+msa`, once `blaskern.zig`'s
  `groupQuery` stopped dropping a model's `features`), the i6400 kernel
  object asserted `Opc && "Register class not handled!"`
  (`MipsSEInstrInfo.cpp:373`, LLVM 23.1.2); da01d84's did too. The culprit
  was the complex gemm kernel's `run` (`blas3/gemm/kernel.zig`): its
  run-time `conj` picked `acc * sa` or `acc` per accumulator vector. The
  condition, copied from a GPR into R6's `FGR64CC` class as in the entry
  above, lived across the tile's loop, and the greedy allocator spilled it;
  `storeRegToStack` has no store for `FGR64CC`. The same bug as above, met
  on the spill path.
- **Workaround:** `conj` is comptime in the vector path of
  `blas3/gemm/kernel.zig`'s complex kernels (`run`, `rows`, `tile`,
  `vectorTile`); `splitTile`, used only without vectors, keeps it at run
  time, so its scalar tile is still one function.
- **Upstream:** the same <https://github.com/llvm/llvm-project/pull/226838>
  removes `FGR64CC`, so no spill of it is left.
- **Issue text:**

  > **[Mips] R6: "Register class not handled!" spilling a `select`
  > condition (`FGR64CC`)**
  >
  > `llc -mtriple=mipsisa64r6el -mcpu=mips64r6 -O2` (here `zig cc -c -O2`)
  > on
  >
  > ```llvm
  > define void @f(i1 zeroext %c, i64 %n, ptr %p) {
  > entry:
  >   br label %loop
  > loop:
  >   %i = phi i64 [ 0, %entry ], [ %i.next, %loop ]
  >   %a = getelementptr inbounds double, ptr %p, i64 %i
  >   %x = load double, ptr %a, align 8
  >   call void asm sideeffect "", "~{$f0},~{$f1},~{$f2},~{$f3},~{$f4},~{$f5},~{$f6},~{$f7},~{$f8},~{$f9},~{$f10},~{$f11},~{$f12},~{$f13},~{$f14},~{$f15},~{$f16},~{$f17},~{$f18},~{$f19},~{$f20},~{$f21},~{$f22},~{$f23},~{$f24},~{$f25},~{$f26},~{$f27},~{$f28},~{$f29},~{$f30},~{$f31}"()
  >   %nx = fneg double %x
  >   %r = select i1 %c, double %nx, double %x
  >   store double %r, ptr %a, align 8
  >   %i.next = add i64 %i, 1
  >   %done = icmp eq i64 %i.next, %n
  >   br i1 %done, label %exit, label %loop
  > exit:
  >   ret void
  > }
  > ```
  >
  > asserts `Opc && "Register class not handled!"` in
  > `MipsSEInstrInfo::storeRegToStack`, from the greedy allocator's inline
  > spiller. Before greedy the entry block has `%13:fgr64cc = COPY
  > %9:gpr32`, the `SEL_D` condition hoisted out of the loop; it cannot
  > stay in an FPR across the asm, and `FGR64CC` has no spill opcode. Without the asm the copy is kept in a register and
  > `copyPhysReg` asserts "Cannot copy registers" instead. Big-endian
  > `mips64r6`, and `<2 x double>` with `+msa`, too; `-O0` compiles.

### MIPS MSA: LLVM's vector `fma` adds the wrong operand

- **Seen:** once the MSA cores got MSA (`blaskern.zig`'s `groupQuery`),
  the `mips64r6el` rig failed 20 tests on the i6400, i6500 and p6600 cores:
  `axpy`, `gemv`, `ger`, `trmv`, `gbmv` and `gemm` were far off.
  `@mulAdd(V, x, alpha, y)` became `fmadd.w $w2, $w0, $w1` with `x` in
  `$w2`, which computes `x + alpha*y`. Every `@Vector` `@mulAdd` on MSA is
  wrong this way, `f32` and `f64`, release 5 and 6, and so are scalar ones
  the SLP vectorizer packs (`axpy`'s strided loop and the complex gemm
  tile on release 6, whose scalar `@mulAdd`s fuse). The Loongson 3A R4
  object, built with MSA all along, had it too; no rig runs that core.
- **Workaround:** `fma.compiles` is false with MSA, so `fma.madd`,
  `fma.fusible` (and with it `fuseFor` and the gemm tile's `madd`) and
  `vector.zig`'s `mulAdd` multiply and add separately there. The
  reference's GCC build fuses (`fmadd.w`), so MSA cores can differ from it
  in the last bit. The `loongson3r4` object, modeled fused (T86), takes
  `fmadd` as inline asm with the operands right (`fma.fusedMadd`).
- **Issue text:**

  > **[Mips] MSA: `llvm.fma` on vectors selects `fmadd` with the addend as
  > a multiplicand**
  >
  > `llc -mtriple=mipsisa64r6el -mattr=+msa -O2` (here `zig cc -S -O2
  > -mcpu=mips64r6+msa`) on
  >
  > ```llvm
  > define void @f(ptr %a, ptr %b, ptr %c, ptr %r) {
  >   %x = load <4 x float>, ptr %a
  >   %y = load <4 x float>, ptr %b
  >   %z = load <4 x float>, ptr %c
  >   %m = call <4 x float> @llvm.fma.v4f32(<4 x float> %x, <4 x float> %y, <4 x float> %z)
  >   store <4 x float> %m, ptr %r
  >   ret void
  > }
  > declare <4 x float> @llvm.fma.v4f32(<4 x float>, <4 x float>, <4 x float>)
  > ```
  >
  > gives
  >
  > ```text
  > ld.w    $w0, 0($6)
  > ld.w    $w1, 0($5)
  > ld.w    $w2, 0($4)
  > fmadd.w $w2, $w1, $w0
  > st.w    $w2, 0($7)
  > ```
  >
  > `fmadd.w wd, ws, wt` is `wd + ws*wt` (MSA manual MD00866 rev. 1.12,
  > FMADD.df; QEMU agrees), so this
  > stores `x + y*z` where `fma(x, y, z)` is `x*y + z`. The cause:
  > `FMADD_W_DESC` and `FMADD_D_DESC` in `MipsMSAInstrInfo.td` get the
  > pattern `(fma $wd_in, $ws, $wt)` from `MSA_3R_4R_DESC_BASE`. `<2 x double>`
  > (`fmadd.d`) and `-mcpu=p5600 -mattr=+msa` do the same. Separate `fmul`
  > and `fadd` are right.

### MIPS R6: a scalar `fma` is a library call, not `maddf`

- **Seen:** on `mips64r6` and `mips32r6`, freestanding and Linux, with or
  without MSA, a scalar `@mulAdd` (`llvm.fma.f32`/`f64` in the IR) becomes a
  call to `fmaf`/`fma`, though release 6 has the fused `maddf.s`/`maddf.d`.
  Not a feature we leave out: LLVM's MIPS backend has none for it (Zig's
  `mips` feature list has no `fma`; `MipsISelLowering.cpp` sets `ISD::FMA` on
  `f32` and `f64` to `Expand`), and `-mcpu=mips64r6` is all it takes.
  Pending upstream: <https://github.com/llvm/llvm-project/pull/233171> (open)
  selects `MADDF`/`MSUBF` for it.
  Speed only; the call is correctly rounded, so the bits are `maddf`'s.
- **Where:** `fma.hardware` is true on R6, so the main library of an R6
  target (the portable kernels, `rsimd.mulAdd`) calls it per multiply-add.
  The R6 kernel objects have MSA, where `fma.compiles` keeps them unfused.
- **Workaround:** none taken. An inline-asm `maddf` in `fma.madd` and
  `fma.mulAdd` for R6 scalars would do, at the cost of hiding the operation
  from LLVM's scheduling and vectorizers.
- **Issue text** (not needed while that PR is open):

  > **[Mips] R6: scalar `llvm.fma` is expanded to a libcall instead of
  > `maddf.s`/`maddf.d`**
  >
  > `llc -mtriple=mipsisa64r6el -O2` (here `zig build-obj -OReleaseFast
  > -target mips64el-freestanding -mcpu=mips64r6`) on
  >
  > ```llvm
  > define float @f(float %a, float %b, float %c) {
  >   %r = call float @llvm.fma.f32(float %a, float %b, float %c)
  >   ret float %r
  > }
  > declare float @llvm.fma.f32(float, float, float)
  > ```
  >
  > calls `fmaf` (`jalr $25` through `%call16(fmaf)`); `double` calls `fma`.
  > MIPS32r6 and MIPS64r6 have the fused `maddf.fmt` (`fd = fd + fs*ft`, one
  > rounding), so `fma` could select `maddf.s`/`maddf.d` with the addend tied
  > to the destination, as MSA's `fmadd` is used for vector `fma`.

### MIPS before R6: a plain `a*b + c` becomes `madd.fmt`, which Loongson 3A fuses

- **Seen:** on `mips4`, `mips64`, `mips64r2` and `mips64r5` (with MSA too),
  `a*b + c` without fast-math is `madd.d`/`madd.s` and `a*b - c` is
  `msub.d`, even where the product has another use; `c - a*b` stays
  `mul.d`, `sub.d`. The MIPS64 spec rounds the product, so LLVM takes them
  as exact, but Loongson's GS464 and GS464E round once: Loongson 3A1000
  user manual vol. 2 (GS464, V1.4) §7.3.2, and 3A2000/3B2000 vol. 2
  (GS464E, V1.03) and 3A3000/3B3000 vol. 2 (GS464E, V1.2) §2.4.10. GCC
  models it (gcc 15.2 `mips.h` `ISA_HAS_FUSED_MADD4`: `-march=loongson3a`,
  `gs464`, `gs464e`, `gs264e`, `r8000`); LLVM has no Loongson CPU. The 3A4000
  register manual (V1.5) does not say what GS464V does.
- **Where:** every pre-R6 rstd build's main library, modeled unfused
  (`fma.hardware` false): on that silicon it fuses where LLVM picked
  `madd.fmt` and rounds elsewhere. The `loongson3r3` and `loongson3r4`
  objects are modeled fused (T86, design.md's Settled decisions), as
  OpenBLAS's `-march=loongson3a` build is (its `.S` kernels' `madd.fmt`,
  its C under GCC's default `-ffp-contract=fast`), so LLVM's picks there
  fuse as modeled; a `c - a*b` it leaves as `mul.d`, `sub.d` still rounds.
  The rigs run QEMU, which rounds all of it (below).
- **Workaround:** for the Loongson objects, `fma.fused_madd` emits fused
  multiply-adds as inline `madd.fmt`. `+nomadd4` stops the selection, but
  it also rejects inline `madd.fmt` (with exit 0, "Inline assembly errors
  exit 0").
- **Issue text:**

  > **[Mips] No Loongson 3A CPU: strict `fadd(fmul)` selects `madd.fmt`,
  > which GS464/GS464E execute fused**
  >
  > `llc -mtriple=mips64el -mcpu=mips64r2 -O2` (here `zig build-obj
  > -OReleaseFast -target mips64el-linux -mcpu=mips64r2`) on
  >
  > ```llvm
  > define double @f(double %a, double %b, double %c) {
  >   %m = fmul double %a, %b
  >   %r = fadd double %m, %c
  >   ret double %r
  > }
  > ```
  >
  > gives `madd.d $f0, $f14, $f12, $f13`. Right by the MIPS64 spec, where
  > MADD.fmt rounds the product, but Loongson's GS464 and GS464E cores
  > (3A1000 to 3A3000) round once, as their manuals say and GCC models
  > (`ISA_HAS_FUSED_MADD4`). Code built for a generic MIPS64 target thus fuses
  > on them where `madd.fmt` is selected and rounds elsewhere (MSA vectors,
  > `c - a*b`). A Loongson 3A CPU that treats `madd4` as fused (legal for
  > `fma`, not for strict `fadd(fmul)`) would match GCC's `-march=loongson3a`;
  > `-mattr=+nomadd4` is the only control today, and it also rejects
  > `madd.fmt` in inline assembly.

### MIPS MSA: LLVM takes `fmax_a` for `fabs`, which keeps the sign

- **Seen:** once the routed `@Vector` kernels went to mips (T25), the
  `mips64r2el` rig's `samax` on `loongson3r4` (MSA) gave 1.9792 for 1.9802
  at a stride of 3, and its amin test hung: the scalar search never met
  the block's key. `@abs(v)` became
  `fmax_a.w $w0, $w0, $w0`; FMAX_A is the operand of the larger
  magnitude, sign kept (QEMU's `msa_helper.c`, `FMAXMIN_A`; the MSA
  manual's FMAX_A, from memory), so `|x|` came out as `x` and negative keys
  lost. A plain sign mask (`and` with `0x7fffffff`) folds back into `fabs`
  and the same `fmax_a`, as do scalar `@abs`es the SLP vectorizer packs
  (the complex extremes' `|re| + |im|`). Release 5 and 6, `f32` and `f64`.
- **Workaround:** `util.abs` clears the sign bits with a mask whose word
  passes through an empty `asm`, so LLVM cannot see it is one (`and.v`);
  only on MSA. `vector.zig` takes it, as do the complex extremes' keys
  (`extremum.zig`, `interface.zig`'s `extreme` and `cabs1`) and the
  from-OpenBLAS tests' `abs1`. T89 compiled the kernel objects of every
  MSA core and the BLAS tests for an MSA target CPU (`i6400`,
  `mips64r5+msa`, `p5600+msa`): the only other `fmax_a` were
  `extreme`'s `incx == 0` complex key (R6) and `abs1` (packed by SLP).
- **Issue text:**

  > **[Mips] MSA: `fabs` on vectors selects `fmax_a.df`, which returns
  > the signed operand**
  >
  > `zig build-obj -OReleaseFast -target mips64el-freestanding-none
  > -mcpu=mips64r5+msa` (`llc -mattr=+msa`) on
  >
  > ```zig
  > export fn absA(p: *[4]f32) void {
  >     const v: @Vector(4, f32) = p.*;
  >     p.* = @abs(v);
  > }
  > ```
  >
  > gives `fmax_a.w $w0, $w0, $w0`. FMAX_A.df writes whichever operand has
  > the larger absolute value, with its sign, so this stores `v`, not
  > `|v|`; QEMU computes it so. `and` with the sign mask is combined into
  > the same `fabs`. Selecting `bclri.w`/`bclri.d` on bit 31/63 (as
  > `and` with a mask the combiner cannot see already gives `and.v`)
  > would be right.

## Zig standard library

### Zig 0.17.0's `Io.Dispatch` does not compile

- **Seen:** on `aarch64-macos` and `x86_64-macos`, Zig 0.17.0's
  `lib/std/Io/Dispatch.zig:439` sets `.processReplacePath` in its `Io.VTable`
  literal, a field 0.17.0's `Io.VTable` lacks ("no field named
  'processReplacePath'"). Anything that instantiates `Io.Dispatch` fails to
  build, so rstd's `threading.usesFutures` (T106) broke `rstd-check-zig017`'s
  macOS targets, which zfls builds.
- **Workaround:** `threading.Evented` is `void` for `Dispatch` on Zig 0.17, as
  it already was for `Uring`, and `usesFutures` and the futures cap test
  `Evented`, never `Io.Dispatch` itself, so 0.17 never instantiates it. The
  dev toolchain is unaffected. Drop the gate when zfls moves past 0.17.
- **Issue text:** fixed upstream after 0.17.0 (master's `Io.VTable` has the
  field); nothing to file.

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
  smaller set". Both come from the open
  <https://codeberg.org/ziglang/zig/pulls/31235>, not upstream master, so
  this belongs in that PR's review.
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

### arm: most CPU models have no FPU

- **Seen:** `std.Target.arm.cpu.cortex_a15` has no VFP or NEON on
  0.18.0-dev.751 (0.17 gave it `vfp4`, never `neon`), so
  `-mcpu=cortex_a15` with `eabihf` compiles every float operation to
  `__aeabi_*` calls and passes floats in core registers. Of 100 models only
  9 have an FPU feature; `cortex_a5` to `cortex_a17`, `cortex_r5` and
  others lost theirs since 0.17. LLVM takes a core's default FPU from its
  target parser (`cortex-a15`: `neon-vfpv4`), which Zig's models skip.
  Deliberate: Zig's generator now omits features ARM's cores may lack
  (dff54eecb0 and b0ce4f0996, 2026-10-20), so `+vfp4,+neon` must be asked for. A
  `vmsr` in the `arm` rig's boot code then failed to assemble, though the
  build exited 0, so the FPU stayed off and every converted kernel was
  skipped as not runnable.
- **Workaround:** the `arm` rig adds `vfp4`, `neon`, `d32`, `hwdiv` and
  `hwdiv_arm` (`build/qemutest.zig`); the boot code writes its `vmsr` by
  encoding (`test/qemu/arm.zig`).
- **Issue text** (likely by design, per those commits; check before filing):

  > **arm: CPU models lack their default FPU, so `-mcpu=cortex_a15` is soft-float**
  >
  > `zig build-obj -target arm-freestanding-eabihf -mcpu=cortex_a15` of
  > `export fn f(a: f32, b: f32) f32 { return a + b; }` calls
  > `__aeabi_fadd` with the operands in r0 and r1; clang's
  > `--target=armv7a-none-eabihf -mcpu=cortex-a15` emits `vadd.f32`. LLVM
  > keeps each core's default FPU and extensions (`ARMTargetParser.def`)
  > outside its `.td` processor definitions, and the generated
  > `std.Target.arm` models do not include them. On 0.18.0-dev.751 only 9
  > of 100 models have any FPU feature; 0.17 had 20.

### `Io.Uring`: `futexWake` of `maxInt(u32)` wakes one waiter

- **Seen:** `Uring.futexWake` (`lib/std/Io/Uring.zig:2052`) puts `max_waiters`
  in an `IORING_OP_FUTEX_WAKE` as is. The kernel's `futex_wake` takes the count
  as an `int`, so `maxInt(u32)` is -1 and wakes one waiter. `Io.Threaded`
  clamps it to `maxInt(i32)` (`Threaded.zig:1162`). std's own
  `Io.Event.set` passes `maxInt(u32)`. On `qemu-test-linux-x86_64` (Debian 13
  kernel, TCG), with `thread_limit = 3`, `Event.set` woke 1 of 3 waiters, and
  rstd's level-3 barriers hung. With `thread_limit = 0` it passes by luck: the
  fibers' `FUTEX_WAIT`s reach the kernel only after the word changed. 0.18.0-dev.751,
  2026-10-29.
- **Workaround:** `threading.Signal.set` wakes `maxInt(i32)`; its test
  "a signal wakes every sleeper" fails without it.
- **Issue text:**

  > **std.Io.Uring: futexWake(maxInt(u32)) wakes only one waiter, so Event.set hangs**
  >
  > `Uring.futexWake` passes `max_waiters` unchanged as the count of an
  > `IORING_OP_FUTEX_WAKE`. The kernel's `futex_wake` takes it as an `int`, so
  > `std.math.maxInt(u32)`, which `Io.Event.set` passes, becomes -1 and wakes
  > one waiter. `Io.Threaded` clamps to `maxInt(i32)`. With worker threads,
  > waiters on an `Event` past the first stay asleep:
  >
  > ```zig
  > test "Event.set wakes every waiter" {
  >     var ev: std.Io.Evented = undefined;
  >     try ev.init(std.heap.page_allocator, .{ .thread_limit = 3 });
  >     defer ev.deinit();
  >     const io = ev.io();
  >     const W = struct {
  >         fn wait(e: *std.Io.Event, wio: std.Io, n: *std.atomic.Value(u32)) void {
  >             e.waitUncancelable(wio);
  >             _ = n.fetchAdd(1, .acq_rel);
  >         }
  >     };
  >     var event: std.Io.Event = .unset;
  >     var woken: std.atomic.Value(u32) = .init(0);
  >     var group: std.Io.Group = .init;
  >     for (0..3) |_| try group.concurrent(io, W.wait, .{ &event, io, &woken });
  >     try io.sleep(.fromMilliseconds(500), .awake);
  >     event.set(io);
  >     try io.sleep(.fromMilliseconds(500), .awake);
  >     const got = woken.load(.acquire);
  >     io.futexWake(std.Io.Event, &event, std.math.maxInt(i32)); // let the rest go
  >     try group.await(io);
  >     try std.testing.expectEqual(@as(u32, 3), got); // found 1
  > }
  > ```
  >
  > Clamping as `Threaded` does fixes it.

### `Io.Uring` and `Io.Dispatch`: `Group.await` blocks inside a context switch

- **Seen:** with worker threads, `Group.concurrent` then `Group.await`, repeated,
  crashes within 2000 rounds on both `Evented`s. The awaiter's pending
  `.group_await` runs in `SwitchMessage.handle`, on the context switched to,
  and takes the group's futex mutex there. When a task ending on another thread
  holds it, `Uring` futex-waits from the idle context and trips
  `assert` in `currentFiber` (`Uring.zig:124`, from `:1305`, `:1676`, `:1587`
  or `:1603`); `Dispatch` wakes from `Group.unlock` there and faults in
  libdispatch (`Sleeper.init`). `Dispatch.Fiber.destroy` does the same on the
  backing allocator's mutex. Level 3's threaded tests hit it at once on
  `qemu-test-linux-x86_64` and on the M6. 0.18.0-dev.751, 2026-10-30.
- **Workaround:** on `Dispatch`, rstd's worker tasks are plain futures
  (`threading.Futures`, chosen by the comptime `threading.dispatch_futures`),
  which pass the same loop; the threaded tests run there with worker threads.
  `Uring` keeps `Group` and no worker threads in the tests (its futures hit the
  entry below), so threaded BLAS on a multithreaded `Uring` is unsupported
  (threading.md). `threading.initEvented` tells `Evented` its page allocator
  needs no mutex.
- **Revert check:** `zig build ssh-repro-group-aarch64-macos` runs
  `tools/evented_group_repro.zig` (std only) on the Mac, `futures` then
  `group`, 20000 rounds each; rerun it on each toolchain bump. Once `group`
  passes, set `dispatch_futures` false and delete the futures path. The file
  checks `Uring` too on a Linux host; no step runs it there.
- **Issue text:**

  > **std.Io (Uring, Dispatch): Group.await crashes with worker threads**
  >
  > `SwitchMessage.handle` runs `.group_await` and `.group_cancel` after the
  > switch, and `Group.lock`/`unlock` may then futex-wait or wake there. On
  > `Uring` that trips the `currentFiber` assert (no fiber: the idle context);
  > on `Dispatch` it faults in libdispatch. Repeating this crashes within a
  > few thousand rounds on `x86_64-linux` and `aarch64-macos`:
  >
  > ```zig
  > var ev: std.Io.Evented = undefined;
  > try ev.init(std.heap.page_allocator, .{}); // Uring: .thread_limit = 3
  > defer ev.deinit();
  > const io = ev.io();
  > for (0..20000) |_| {
  >     var group: std.Io.Group = .init;
  >     for (0..3) |_| try group.concurrent(io, task, .{&counter});
  >     try group.await(io);
  > }
  > ```

### `Io.Uring`: an awaited fiber is reused before it has switched away

- **Seen:** with worker threads, `io.concurrent` then `Future.await`, repeated,
  panics within 2000 rounds at `Uring.zig:1424` ("switched to dead fiber").
  A finishing task marks itself `finished` (`:1423`) and then switches away on
  its own stack; an awaiter on another thread sees `finished`, `destroy`s the
  fiber to its free list, and a new task takes it while the old one is still
  switching. `Dispatch` passes the same 20000 rounds. 0.18.0-dev.751,
  2026-10-30.
- **Workaround:** as above: no worker threads on `Uring`, so the futures
  fallback is for `Dispatch` only.
- **Issue text:**

  > **std.Io.Uring: Future.await frees a fiber that is still switching away**
  >
  > `AsyncClosure.call` publishes `Fiber.finished` before its final `yield`
  > leaves the fiber's stack. An awaiter on another thread can then `destroy`
  > it and `concurrent` can reuse it, so the old `yield` writes into the new
  > task's fiber. Repeating `io.concurrent` of three trivial tasks and
  > awaiting them, with `thread_limit = 3`, panics "switched to dead fiber"
  > within a few thousand rounds.

### `Io.Dispatch`: a fiber's floating-point control state stays with the thread

- **Seen:** `Dispatch` resumes a fiber that waited, or woke a waiter
  (`futexWake` yields too), on whichever queue thread runs it next, the
  process's main fiber included. The context switch (`lib/std/Io/fiber.zig`)
  lists `fpcr` (x86_64: `mxcsr`) as clobbered rather than saving it, so the
  control state stays with the thread. A worker that adopted the caller's
  rounding toward zero then ran under another thread's, and left its first
  thread rounding toward zero: level 3's rounding test lost its bits on the
  M6. 0.18.0-dev.751, 2026-10-30.
- **Workaround:** `threading.Lent`: a level 3 worker gives its thread's own
  state back around each futex wait and wake and takes the caller's again
  after; the caller re-adopts its own after them (threading.md, Futures on
  `Dispatch`).
- **Issue text:**

  > **std.Io.Dispatch: floating-point control state does not follow a fiber**
  >
  > The fiber context switch (`Io/fiber.zig`) marks the floating-point
  > control register clobbered (`fpcr` on aarch64, `mxcsr` on x86_64)
  > instead of saving it. A fiber that sets rounding toward zero, then
  > waits on a futex, can resume on another `Dispatch` thread under that
  > thread's rounding, while the thread it left keeps toward zero for
  > whatever runs there next. The process's main fiber migrates too. Either
  > carry the control state in the fiber's context, or document that it is
  > per thread and undefined across a yield.

### Zig 0.17.0: `Io.Uring` does not compile

- **Seen:** with the 0.17.0 release, any use of `Io.Evented` on Linux fails:
  `lib/std/Io/Uring.zig:759:14: error: no field named 'processReplacePath'
  in struct 'Io.VTable'`. `Uring.io()` fills a field that 0.17.0's `Io.VTable`
  does not have. 0.18.0-dev.751 builds it (the field is gone from both).
  Found by `rstd-check-zig017` for x86_64 and aarch64 Linux, 2026-10-29.
- **Workaround:** `threading.Evented` is void on Zig 0.17 with `Uring`, so
  `blas3/threaded.zig`'s fiber test skips there. Gated on the version, since
  the field `Uring` sets is private. Drop it once rstd needs 0.18.
- **Issue text** (0.17.0 only; check a 0.17 patch release before filing):

  > **std.Io.Uring does not compile on 0.17.0: `processReplacePath` is not an `Io.VTable` field**
  >
  > `var ev: std.Io.Evented = undefined; ... _ = ev.io();` on
  > `x86_64-linux` with Zig 0.17.0 fails in `lib/std/Io/Uring.zig:759`:
  > "no field named 'processReplacePath' in struct 'Io.VTable'".
  > `Io.VTable` has `processReplace` only.

### aarch64 SVE: `Io/fiber.zig` clobbers `ffr`, which LLVM rejects

- **Seen:** `qemu-check-linux-aarch64-sve` (Debug, CPU `neoverse_v1`) fails compiling the test
  executable: "inline asm clobber list contains reserved registers: FFR". `lib/std/Io/fiber.zig`'s
  aarch64 context switch lists `.ffr = true`; LLVM warns on an SVE target, and Zig makes the
  warning an error. The build still exits 0 (see "Inline assembly errors exit 0"). 0.18.0-dev.751,
  2026-10-30. Baseline aarch64 targets build, so only hosted SVE targets meet it.
- **Where:** the `linux-aarch64-sve` rig (T40), and any hosted SVE build of rstd.
- **Workaround:** `threading.Evented` is void on SVE targets, so the fiber tests are skipped and
  nothing compiles the switch (`src/linalg/blas/threading.zig`). Drop it once std builds there.
- **Issue text:** none yet: first reduce it and check how LLVM expects FFR to be clobbered.

## Zig package manager and build system

### `zig fetch` needs an archive and a `Content-Type` header

- **Seen:** `zig fetch` of a single file fails with "unrecognized
  'Content-Type' header: 'text/plain'", and of a `.tar.gz` from
  snapshot.debian.org (which sends no `Content-Type`) with "missing
  'Content-Type' header", although the URL ends in `.tar.gz`. Single files
  are open upstream as <https://codeberg.org/ziglang/zig/issues/32041>, whose
  first example also gets "missing 'Content-Type' header".
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
  > says `.tar.gz`. Using the extension when the header is absent, as is
  > already done for `application/octet-stream`, would make such mirrors usable. Separately,
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

### A link interrupted midway leaves a cache entry without its executable

- **Seen:** a `zig build` killed while linking left `.zig-cache/o/<digest>/`
  with only `<name>_zcu.obj` and a manifest in `.zig-cache/h/`. Later builds
  treat it as cached, and the run step fails with `FileNotFound` on the EXE.
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

### Its bundled Aro does not build on current Zig

- **Seen:** at first (9b5eb01, 2026-10-12) translate-c's `build.zig` failed
  because `standardOptimizeOption` returned an optional (aa8c774e0f); Zig
  reverted that (21141cb88a). On 0.18.0-dev.751 the package still does not
  build, now in Aro: at 9b5eb01, `build/generate_def.zig` uses the removed
  `std.ArrayListUnmanaged`; at 6867c49 (latest, 2026-10-26),
  `src/aro/Target.zig` and `Compilation.zig` use names Zig renamed or
  dropped (`Os.Tag.amdpal`, s390x `transactional_execution`, spirv
  `v1_0`).
- **Workaround:** rstd's `translate_c` dependency is commented out and
  `src/buildutils/renderdoc.zig` uses the deprecated built-in `addTranslateC`.
- **Issue text:**

  > **translate-c does not build on Zig master: Aro uses removed std and target names**
  >
  > With Zig 0.18.0-dev.751+069480a1a, building the `translate-c` artifact
  > at 6867c49 fails in its Aro dependency: `Target.Os.Tag` has no
  > `amdpal`, `Target.s390x.Feature` no `transactional_execution`, and
  > `Target.spirv.Feature` no `v1_0` (`src/aro/Target.zig`,
  > `src/aro/Compilation.zig`). The pin before it (9b5eb01) fails on
  > `std.ArrayListUnmanaged` in Aro's `build/generate_def.zig`. Updating
  > the Aro pin to one that follows Zig master would fix both.

## QEMU

### WHPX passes the host's AMD CPUID through to an Intel-only guest

- **Seen:** on an AMD Ryzen host, the macOS rig under `-accel whpx` stops
  right after `EXITBS:START`, whatever `-cpu ...,vendor=GenuineIntel` says.
  Under TCG it boots.
- **Workaround:** the macOS rig uses TCG.
- **Issue text:** none: QEMU 11.1.0 intercepts CPUID under WHPX, leaf 0
  included (79b3cb1f01, "whpx: i386: introduce proper cpuid support"); our
  QEMU is 11.0.1. Retest on 11.1 first. Possibly related:
  <https://gitlab.com/qemu-project/qemu/-/issues/4385> (open; a macOS guest on
  an AMD host dies in early boot under WHPX).

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
- **Issue text:** none; long-standing missing features, not bugs. AVX-512 is
  requested in <https://gitlab.com/qemu-project/qemu/-/issues/2878> (open).

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
- **Upstream:** reported as <https://gitlab.com/qemu-project/qemu/-/issues/3582>
  (open; also #4682, closed): `access_ptr` in `target/i386/tcg/access.c` tests
  `offset <= ac->size1 - len`, which wraps when `size1 < len`, so the 8-byte
  significand access goes through the first page's host pointer. Still so on
  master (5624b266cd, 2026-10-28). Nothing to file.
- **Issue text** (superseded by #3582):

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
- **Open question:** the SPARC V9 manual (SAV09R1459912, A.49 "Shift") says
  "SLL and SLLX shift all 64 bits of the value in r[rs1]", so this is a qemu
  bug: `do_shift_r` and `do_shift_i` in `target/sparc/translate.c` (master
  5624b266cd) zero-extend SLL's result from 32 bits. If real UltraSPARC
  hardware zero-extends instead, it is an OpenBLAS bug (its kernels would then
  fail on hardware too). Check on hardware before filing either.
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
  > The sparc level-1 kernels (for example `copy.S`, `snrm2.S`) scale
  > `INCX`/`INCY` with `sll INCX, BASE_SHIFT, INCX`. With 64-bit `BLASLONG`
  > and a negative increment, `sll` produces a zero-extended 32-bit value,
  > so the kernel steps forward by about 4 GB per element. Using `sllx`
  > under `__64BIT__` fixes it.

### powerpc: `xststdcsp` matches a single-precision denormal whatever DCMX asks

- **Seen:** QEMU 11.0.1 TCG, `-cpu power9` and `power10`. LLVM lowers
  `std.math.isInf` on `f32` (ReleaseSafe, `-mcpu=pwr9`) to
  `xststdcsp 0, 1, 48`: DCMX 48 is +Inf and -Inf, and the operand is the
  `f32` as the FPR holds it, in double format, as the ISA wants. LLVM is
  right. Under QEMU, it is true for `f32` denormals (0x00000001, 0x007fffff,
  0x80000001), and right for normals, zeros, infinities, NaN and every `f64`
  (`xststdcdp`). `helper_XSTSTDCSP` in `target/ppc/fpu_helper.c` has
  `match = float64_tstdc(b->VsrD(0), dcmx) || (exp > 0 && exp < 0x381)`: a
  normal double below 2^-126 matches whatever DCMX selects. The 2022 move of
  `XSTSTDC[SDQ]P` to decodetree (da3c53bac3) dropped the DCMX test; the old helper had
  `match = extract32(dcmx, 0 + !sign, 1)` for that range. Unchanged on
  master.
- **Workaround:** the portable `nrm2` does not call `isInf`
  (`blas1/nrm2/__root__.zig`). Other `f32` `isInf` calls are wrong for
  denormals on the powerpc64le rigs.
- **Issue text:**

  > **target/ppc: xststdcsp matches single-precision denormals regardless of DCMX**
  >
  > `helper_XSTSTDCSP` computes
  > `match = float64_tstdc(b->VsrD(0), dcmx) || (exp > 0 && exp < 0x381);`,
  > so a value that is denormal in single precision (a normal double with
  > exponent below 0x381) sets CR[BF].EQ whatever DCMX selects. Before
  > XSTSTDC[SDQ]P moved to decodetree the range set `match` from DCMX bit
  > `0 + !sign`, as for double denormals. Test for infinity
  > (`xststdcsp 0, 1, 48`) with f1 = 2^-149 (`lfs` of 0x00000001) and
  > `beq` is taken, where the ISA matches it only for DCMX bits 5 and 6
  > (+/- denormal). Clang emits exactly this for `__builtin_isinf` on a
  > `float` with `-mcpu=pwr9`. Suggested fix:
  >
  > ```c
  > match = float64_tstdc(b->VsrD(0), dcmx) ||
  >         (exp > 0 && exp < 0x381 && extract32(dcmx, 0 + !sign, 1));
  > ```
  >
  > Reproduce, in a powerpc64le guest built for POWER9 at -O2:
  >
  > ```c
  > __attribute__((noinline)) int f(float x) { return __builtin_isinf(x); }
  > union { unsigned u; float f; } d = { 1 };   /* 2^-149 */
  > /* f(d.f) is 1 under qemu -cpu power9; isinf(2^-149) is 0 */
  > ```

### mips: the Loongson 3A models round `madd.fmt`'s product

- **Seen:** `helper_float_madd_d` and its kin (`target/mips/tcg/fpu_helper.c`,
  v11.0.1 and master) multiply and add with two roundings for every CPU,
  `-cpu Loongson-3A1000` and `Loongson-3A4000` included. The GS464 and
  GS464E manuals say one rounding ("MIPS before R6" above).
- **Where:** the `mips64r2` (3A1000) and `mips64r2el` (3A4000) rigs: the
  converted `loongson3r3` kernels' `madd.fmt` (`saxpy`, `daxpy`, `ddot`,
  the complex `dot`s), the Loongson objects' inline ones and LLVM's round
  there and fuse on the silicon.
- **Workaround:** none needed: the tests that reach them compare against
  a tolerance or with themselves, which hold either way. The rigs cannot
  check fused Loongson bits, so axpy's fused-bit test for these cores
  skips on rigs. `loongson3r4`'s MSA `fmadd` is fused under QEMU too.
- **Issue text:**

  > **target/mips: MADD.fmt and kin round the product on Loongson-3A CPUs,
  > which fuse it**
  >
  > `helper_float_madd_d` computes `float64_mul` then `float64_add` for all
  > pre-R6 CPUs. Loongson's GS464 (3A1000 user manual vol. 2, §7.3.2) and
  > GS464E (3A3000 user manual vol. 2, §2.4.10) do MADD.fmt, MSUB.fmt,
  > NMADD.fmt and NMSUB.fmt with one rounding, as GCC's
  > `-march=loongson3a` assumes when it emits `madd.d` for `fma()`. Under
  > `-cpu Loongson-3A1000`, `madd.d` with `fs = ft = 1 + 2^-52` and
  > `fr = -(1 + 2^-51)` gives 0 where a fused one gives 2^-104. These
  > CPUs could use `float64_muladd`.

### powerpc: OpenBIOS hangs before its banner on `-cpu 750fx`

- **Seen:** QEMU 11.0.1 (OpenBIOS 1.1, Sep 2024) here and 11.1.2
  (OpenBIOS built Jun 2026) on the Mac, `-M mac99` or `g3beige`
  with `-cpu 750fx` or `750gx` prints nothing; `-cpu 750`, `750cxe` and `G4`
  boot the same image. `-d int` shows an ISI at 0xfff08b50 (in the ROM),
  then `rfid` (0x4c000024, not in 32-bit CPUs) at 0x251c, then nothing.
  OpenBIOS knows both PVRs (`arch/ppc/qemu/init.c`, `ppc_defs`, 0x7000 and
  0x7002, e5ac46dd2), so the fault is elsewhere; cause not found.
- **Where:** the `powerpc-750` rigs, for RandyOS's iBook G3 (a 750FX).
- **Workaround:** the rigs run `-cpu 750`. The 750FX adds HID2, four more
  BAT pairs and L2 controls, none of which user code reaches; Zig's `750`
  model is the same for both.
- **Issue text:** none yet: not narrowed to QEMU or OpenBIOS.

## Linux kernel

### arm64: Debian's 6.12 kernel has no SME

- **Seen:** `qemu-test-linux-aarch64` under HVF on the M6 never runs `vortexm4`, and the TCG
  rigs (`-cpu max`) never run `armv9sme`: `AT_HWCAP2` has no `HWCAP2_SME`. QEMU 11.1.2 (our
  Nix pin) does pass SME2 to an HVF guest (QEMU's "hvf/arm: expose FEAT_SME2 to guest if
  available", 2026-03-06). The kernel is Debian's 6.12.107 (netboot 20250803+deb13u7), and 6.12 to 6.15
  have `ARM64_SME` `depends on BROKEN` (`arch/arm64/Kconfig`), so it is built without SME: the
  image has SVE's strings and none of SME's. SME instructions trap there too.
- **Where:** the `linux-aarch64` rigs; SME's Linux detection and kernels go untested there.
- **Workaround:** none; `ssh-test-aarch64-macos` runs `vortexm4` on the M6, and the
  freestanding `aarch64-sve` rig runs `armv9sme`.
- **Issue text:** none: fixed upstream in 6.16 (33c4618d0ac0, "arm64/fpsimd: Allow
  CONFIG_ARM64_SME to be selected"). A netboot kernel of 6.16 or later would run them.

## OpenBLAS

The kernel bugs we fixed in the converted sources are listed under "Fixes to
upstream sources" in [`design.md`](design.md#fixes-to-upstream-sources): `x86_64` argument and flag handling
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

The edge-case tests (`_test/blas1/edge.zig` and `_test/blas2/edge.zig`, on every rig)
found these. In rstd's extremes the first NaN wins ([`readme.md`](readme.md#one-contract-on-every-core)); LAPACK 3.10's `nrm2` is infinite for an
infinity; `kernel/arm/axpby.c` reads no operand whose scalar is zero.

Not used: no core takes a converted extreme kernel. Upstream follows
Reference BLAS 3.12's `i?amax`, which replaces its pick only on a strict
`>`, so NaN is picked only as the first element, or worse:

- **Extremes, lane by lane.** `x86_64/amax_sse.S`, `amax_sse2.S`,
  `zamax_sse.S`, `zamax_sse2.S`, `iamax_sse.S`, `iamax_sse2.S`,
  `izamax_sse.S` and `izamax_sse2.S` (every x86_64 core, all eight
  routines), `power/amax.S` and its kin (every POWER core), the
  `loongarch64` extremes and `arm64/amax.S` (`amax`, every arm64 core)
  compare in vector lanes, so what a NaN does depends on where it falls: a
  NaN first may be passed over, a later one picked. `power/amax.S`, with a
  NaN last, returned an element smaller than the maximum;
  `loongarch64/amax_lsx.S` returned 0 for a lone NaN. Meeting "first NaN"
  needs each SIMD reduction rewritten.
- **The rest.** The x86, arm, arm64 `iamax`, mips64 and SPARC extremes keep
  the strict compare (after the fixes below), so a NaN after the first
  element is passed over.

Not used: no core takes a converted `nrm2` kernel, so `nrm2` is one
algorithm everywhere ([`readme.md`](readme.md#one-contract-on-every-core)). x86_64's
`nrm2_sse.S` (an unscaled `f64` sum for `snrm2` and `cnrm2`) and the x87
`nrm2.S` and `znrm2.S` of x86 and x86_64 (an unscaled `f80` sum) were right
on special values but their own algorithms; the arm and arm64 ones scale a
one-pass loop as those below do. Those below were wrong:

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

Fixed by the generator's patches, the extremes' kept as `retired` now that
they are not converted
([`design.md`](design.md#fixes-to-upstream-sources)). Issues could read, for the
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

Not used: the x86_64 reductions whose bits follow an operand's address.
No core takes them ([`readme.md`](readme.md#one-contract-on-every-core)); an issue
could read:

> **x86_64 `asum_sse.S`, `asum_sse2.S`, `zasum_sse2.S`, `dot_sse2.S`, `zdot_sse2.S`, `zsum_sse2.S`, `dgemv_n.S` and kin: the result depends on the alignment of `x` or `y`**
>
> They peel one to three elements until `x` (or `y`) is 16-byte aligned
> (`testq $SIZE, X` and the like), which changes which lane each element
> is summed in, so the same vector at another address gives other bits:
> `dasum` of 67 elements differs in the last bit with `x` 8 bytes further
> on, `ddot` of 7 with `y` moved. `dgemv_n.S` accumulates in its buffer
> when `y` is misaligned, and adds that to `y` once. Unaligned loads and a
> fixed split would make them reproducible.

And the kernels rstd changed, by patch or port, to keep `axpy` one contract:

> **`arm64/axpy.S`, `zaxpy.S`, `daxpy_thunderx2t99.S`, `axpy_sve.c`, `loongarch64/axpy_lsx.S`, `axpy_lasx.S`, `caxpy_lsx.S`, `caxpy_lasx.S`, `mips/?axpy_msa.c`, `riscv64/axpy_rvv.c`, `zaxpy_rvv.c`, `wasm/axpy.c`: a zero alpha skips `x`**
>
> They return at once for `alpha == 0`, as `kernel/arm/axpy.c` does, while
> the other cores' kernels multiply. `interface/axpy.c` returns first, but
> `driver/level2` calls `AXPYU_K` with `-B[i]`, which may be zero, so on
> these cores an infinity or NaN in the matrix column does not reach the
> solution of `trsv` and `tbsv`, and on the others it does. Dropping the
> early return makes every core multiply; reference master's `dtrsv.f`
> (`6e7b86c`) removed its zero test for the same reason.

A zero `alpha` with a denormal in `y` also showed the ARMv7 NEON and AltiVec
flushing of `comparisons.md` (10): the portable `saxpy`'s `@Vector` block
turned the denormal into 0 on the `arm`, `linux-arm`, `powerpc` and
`powerpc-pseries` rigs. It was the only level-1 or level-2 kernel that did
(`_test/blas1/denormal.zig`, `_test/blas2/edge.zig`); it now takes its scalar loop where
`util.vectorIeee` says the lanes flush, the `@Vector` kernels leave `s` and
`c` there to the portable ones, and the level-3 tiles are scalar there
(`blas3/gemm/kernel.zig`'s `vectorBits`).

Upstream level-1 and level-2 kernels for these cores are scalar: VFP on
arm, the FPU on POWER before POWER8, none of them AltiVec. Its ARMv7
`sgemm_kernel_4x4_vfpv3_neon.S` and the PPC970 and CELL
`gemm_kernel_altivec*.S` and `zgemm_kernel_altivec*.S` compute in NEON and
AltiVec lanes and never set `VSCR`, so they flush `f32` denormals as these
units do.

Ours, not a converted kernel:

- **Rank updates and a zero factor.** Reference `zher`, `zhpr`, `dsyr2` and
  the rest skip a column whose factor is zero (for rank 2, `x(j)` and `y(j)`
  both). OpenBLAS's `driver/level2/syr_k.c`, `spr_k.c`, `zsyr_k.c` and
  `zspr_k.c` do; `zher_k.c`, `zhpr_k.c` and the rank-2 drivers do not, so
  an infinity elsewhere in `x` meets the zero and gives NaN in that column.
  `banded.zig`'s drivers skip as the reference does. So does `ger`, whose
  reference skips a zero `y(j)` and whose per-core kernels do not:
  `interface.zig` calls the kernel once per run of nonzero `y`, which is one
  call when `y` has no zero, and copies a strided `x` once for all of them.
  OpenBLAS's threaded drivers skip where its serial ones do not, and BLIS
  skips only where an `axpyv` kernel returns for a zero scalar
  ([`comparisons.md`](comparisons.md) (5)).

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

And SPARC `rot`, found by the per-core `rot` checks once they took zero
strides:

> **sparc `rot.S`, `zrot.S`: a zero stride acts on one element as a block**
>
> As in `swap.S`: the strided loop (`.LL51`) loads a block (8 for `rot.S`,
> 4 for `zrot.S`) from `x` and `y` before storing any, so with `incx` or
> `incy` zero later elements of the block are read before earlier ones are
> written. `srot` with `n` = 8 and both strides 0 gives `x` = 2.0356 where
> eight rotations give -1.0088. Sending a zero stride to the one-at-a-time
> loop (`.LL56`) fixes both.

And POWER `rot` and `swap`, the same pattern, found by `utest`'s `zdrot_inc_0`
(upstream met it in <https://github.com/OpenMathLib/OpenBLAS/issues/1469>, and
<https://github.com/OpenMathLib/OpenBLAS/pull/1535> moved POWER6 to `arm/zrot.c`;
CELL and the `ppc440` cores kept the assembly):

> **power `rot.S`, `zrot.S`, `swap.S`, `zswap.S`, `rot_ppc440.S`, `zrot_ppc440.S`: a zero stride acts on one element as a block**
>
> The strided loop (`LL(110)`; the `ppc440` files are pipelined from the
> start) loads a block of elements from `x` and `y` before storing any, so
> with `incx` or `incy` zero later elements of the block are read before
> earlier ones are written. `zdrot_inc_0` (`n` = 4, both strides 0) gets
> `x` = 0.3125 for cell's `zrot.S` where four rotations give -0.21484375,
> and so would ppc440fp2's `zrot_ppc440.S`; the others go wrong from 8
> elements (`swap.S` 16). Sending a zero stride to the one-at-a-time loop
> (`LL(160)`) fixes them.

And the same on other architectures, found by the per-core `rot` and `swap`
checks once they took zero strides; OpenBLAS's `interface/rot.c` and
`swap.c` pass a zero stride to the kernel, as rstd's interface does:

> **x86 `swap.S`, `zswap.S`, arm `swap_vfp.S`, `rot_vfp.S`: zero strides skipped or cut short**
>
> Upstream special-cased zero strides to pass `utest`'s `n` = 4 case
> (OpenBLAS #1469, #1575), each in a way right only for some `n`. x86
> `swap.S` returns at once when both are zero, which is right for an even
> `n` only (one swap, `n` = 1, does nothing); with one zero, the strided
> loop loads four elements before storing any. x86 `zswap.S` does the same
> with two. arm `swap_vfp.S` returns at once when either is zero; arm
> `rot_vfp.S` stops after the first element of the tail (`n` = 2 rotates
> once). The x86 fix sends a zero stride to a one-at-a-time loop (`.L35`;
> a new `.L36` for `zswap.S`); the arm one drops the checks, as the
> strided loops already take one element at a time.

And loongarch64 and mips64:

> **loongarch64 `swap.S`, `rot_lsx.S`, `rot_lasx.S` and mips64 `swap.S`: a zero stride acts on one element as a block**
>
> As on POWER: the strided loops load a block of 8 before storing it, so
> from `n` = 8 a zero stride goes wrong (la64_generic, la264, la464,
> loongson3r3). `crot_lsx.S` and `crot_lasx.S` already send zero strides to
> their one-at-a-time loop (`.L996`); the fix does the same for these (the
> mips64 one sets the count in each branch's delay slot).

And with only the written operand's stride zero, found by the per-core checks
once they took one zero stride beside a nonzero one. The reference is
sequential: `axpy` with `incy` = 0 does `y[0] += alpha*x[i]` for each `i` in
turn.

> **x86_64 `zaxpy_sse2.S`, `axpy_atom.S`, `zaxpy_atom.S`: a zero `incy` updates one element as a block**
>
> The strided loops load several `y` elements before storing any, so with
> `incy` = 0 later updates start from the old `y[0]`. `zaxpy_sse2.S` has a
> one-at-a-time loop for zero strides (`.L58`), but its test sends only both
> zero there, despite its comment, and the loop never steps `x`; `zaxpy`
> with `n` = 3, `incy` = 0 gives -3.9665 where the sequential loop gives
> -1.0370 (prescott, core2, nehalem, opteron and the rest). The `atom` files
> have no such loop: `daxpy` goes wrong from `n` = 2. The fix sends a zero
> `incy` to a stepping one-at-a-time loop.

And POWER and SPARC, as in their `rot` and `swap`:

> **power `axpy.S`, `zaxpy.S`, `axpy_ppc440.S`, sparc `axpy.S`, `zaxpy.S`: a zero `incy` updates one element as a block**
>
> The strided loop loads a block of `y` before storing it. `saxpy` with
> `incy` = 0 goes wrong from `n` = 16 on power3 to power6, ppc970, cell and
> ppcg4, and from `n` = 8 on sparc. Zero strides take the one-at-a-time
> loop, as those do; on sparc it now returns for `n` = 0 first, which its
> loop, run once before the count is tested, did not.

And loongarch64:

> **loongarch64 `axpy_lsx.S`, `axpy_lasx.S`, `caxpy_lsx.S`, `caxpy_lasx.S`, `cswap_lsx.S`, `cswap_lasx.S`: a zero stride acts on one element as a block**
>
> The `axpy` loops for a strided `y` (`.L12`, `.L22`) gather a block of `y`
> before storing it; `saxpy` with `incy` = 0 goes wrong from `n` = 8 (la264,
> la464). `cswap`'s `.L22` reads `x` a block ahead of storing it, so
> `incx` = 0 beside a strided `y` goes wrong from `n` = 4. A zero `incy`
> (either zero stride for `cswap`) takes the one-at-a-time loop.

And C sources rstd ports to Zig:

> **riscv64 `axpy_rvv.c`, `zaxpy_rvv.c` and mips `saxpy_msa.c`, `daxpy_msa.c`: a zero `incy` updates one element as a block**
>
> `axpy_rvv.c` has a branch for `incy` = 0 that starts every lane at
> `y[0]` and sums the lanes into it, so it adds `y[0]` once a lane; `saxpy`
> with `n` = 3 gives -1.2692 where the loop gives 1.6555
> (riscv64_zvl128b, x280). `zaxpy_rvv.c`'s strided segments and the MSA
> files' blocks of four and two load `y` before storing it (i6400, i6500,
> p6600 from `n` = 4). Rstd's ports of these (`axpy/riscv64/rvv.zig`,
> `axpy/mips/msa.zig`) take a zero `incy` one element at a time.

### `cblas_icmax` and `cblas_icmin` scan a complex vector as reals

`interface/Makefile` builds the complex `cblas_icmax`, `cblas_icmin`,
`cblas_izmax` and `cblas_izmin` from `imax.c` without `USE_ABS`, and
`imax.c` only picks a complex kernel under `USE_ABS`. So they call the real
`ismax_k`/`idmax_k` (and the `min` ones) on the complex array: `n` floats at
a stride of `incx` floats, real and imaginary parts alike, which covers only
the first `n / 2` elements at unit stride. Reference BLAS and BLIS have no
signed complex extreme (only `i?amax` and `amaxv`, over `|re| + |im|`).
Reported upstream as <https://github.com/OpenMathLib/OpenBLAS/issues/6107>;
<https://github.com/OpenMathLib/OpenBLAS/pull/6129> (2026-10-15, after our
pin) removed the four functions. Nothing to file.

- **Workaround:** none needed; rstd's `extreme` and `index` make `.max` and
  `.min` over a complex type a compile error that points to `.amax` and
  `.amin`.
- **Issue text:** "`cblas_icmax`/`cblas_icmin` (and the `z` forms) run the
  real `?max` kernel over the complex array's floats, so the returned index
  is not an element index and only half the vector is scanned. Either
  define them (by `|re| + |im|`, as `icamax`) or drop them."

### arm64 `znrm2.S`: `scnrm2` keeps an `f32` running sum

`kernel/arm64/znrm2.S` (`CNRM2KERNEL` of `NEOVERSEN1`, which serves Apple
M-series on macOS and Linux) adds the scaled squares (`SSQ`, `SCALE`) in
turn into one `f32` register. Each term averages 1/3 while the sum grows to
about 5.6M, where an `f32` ulp is 0.5, so for 16M floats the result is
about 1% off (87,202 eps); rstd is 0.3 eps off. The real `snrm2` (`nrm2.S`)
and the x86_64 kernels are exact to rounding.

- **Workaround:** none needed; no core takes the converted `nrm2` kernels.
- **Issue text:** "On arm64, `scnrm2` loses about 1% for 8M complex
  elements: `znrm2.S` keeps its single-precision sum of squares in an `f32`
  register. Accumulate in `f64`, as `snrm2` does."

### x86_64: the C kernels of Haswell and Zen round by compiler

`kernel/Makefile` builds the `HASWELL` and `ZEN` kernels with `-mavx2 -mfma`
under clang but only `-mavx2` under GCC, so the C parts of the C kernels
(`gemv`'s row tails and `add_y`, for one) contract multiply-adds under clang
and not under GCC. Rstd follows clang (`cpu/_corespec.zig`). Against a GCC
build (the 0.3.33 Windows release, GCC 9.3: its `libopenblas.dll`'s producer strings read
`GNU C99 9.3.0`), the converted `gemv` blocks
agree bit-for-bit and those C statements differ in the last bit; built
without contraction, rstd's complex `gemv` agreed in all 17,360 cases tried.

- **Workaround:** none needed; the rounding is the core's either way.
- **Issue text:** "On x86_64, `kernel/Makefile` passes `-mfma` with
  `AVX2OPT` for clang but not for GCC, so the C code of the `HASWELL` and
  `ZEN` kernels (for example the tails of `?gemv_n_4.c`) is contracted into
  FMAs by one compiler and not the other, and a GCC and a clang build give
  different last bits. Passing `-mfma` for GCC too, or
  `-ffp-contract=off` for both, would make them agree."

### x86_64 `dgemv_t_4.c`: SSE3 on the `OPTERON` core

`kernel/x86_64/dgemv_t_4.c` defines `dgemv_kernel_4x2` and `dgemv_kernel_4x1`
in inline assembly ending in `haddpd`, an SSE3 instruction. `KERNEL` makes it
the `DGEMVTKERNEL` of every x86_64 core that names no other, `OPTERON`
included, and `dynamic.c` picks `gotoblas_OPTERON` (`DYNAMIC_OLDER=1`)
for an AMD family 0Fh CPU whose CPUID reports no SSE3. So a transposed `dgemv` with `n % 4 != 0`
would fault there. Not seen on hardware.

- **Workaround:** where the core has no SSE3, the ported driver
  (`blas2/gemv/kernel.zig`) does those two blocks' arithmetic in SSE2
  (`Fallback.t4x2`, `t4x1`); a test checks them against the assembly bit for
  bit.
- **Issue text:** "`dgemv_t_4.c`'s `dgemv_kernel_4x2` and `dgemv_kernel_4x1`
  use `haddpd` (SSE3), yet `OPTERON`, which `dynamic.c` picks for family 0Fh
  AMD CPUs without SSE3, builds `dgemv_t_4.c` (`KERNEL`'s default). A transposed
  `dgemv` with `n % 4 != 0` then executes an SSE3 instruction there. Give
  `OPTERON` an SSE2 `DGEMVTKERNEL`, or add the halves with SSE2."

### x86_64 `cscal.c`, `zscal.c`: the tails treat special values otherwise than the blocks

`kernel/x86_64/zscal.c` and `cscal.c` (built without contraction) multiply the
elements their microkernel block covers by the full product,
`(ar*xr - ai*xi, ar*xi + ai*xr)`. The scalar tail (the last `n % 8`
elements, `n % 16` for `cscal`) and the strided loops branch instead:

- in `zscal.c`, with `ai == 0` (and, strided, the interface's flag) they
  compute `(ar*xr, ar*xi)`, so an infinite `xr` leaves a finite imaginary
  part where the product gives NaN (`0*inf`); `cscal.c`'s tail checks for it;
- with `alpha == 0` and the flag, the unit-stride tail tests only its first
  element for NaN and stores that verdict, NaN or 0, in every tail element;
  an infinity gives 0 where the block's product gives NaN;
- `isinf(da_r)` sits inside `if (da_r == 0.0 || isnan(da_r))`, so it never
  holds.

So one element's result depends on where it falls in the vector.

- **Workaround:** rstd's driver (`blas1/scal/kernel.zig`) takes the blocks'
  full product in the tails and strided loops too, which is the portable
  kernels' contract; the bits of finite values are the reference's.
- **Issue text:** "In `kernel/x86_64/zscal.c` and `cscal.c` the tail after
  the `zscal_kernel_8` block, and the strided path, special-case `alpha`:
  in `zscal.c`, with `ai == 0` an infinite real part no longer makes the
  imaginary part NaN as the block's product does; with `alpha == 0` and `dummy2` set, the
  tail checks only its first element for NaN and writes that result (NaN or 0)
  to every remaining element, and an infinity becomes 0 instead of NaN.
  Results then depend on `n % 8` (`n % 16` for `cscal`). The full product,
  as the blocks compute it, would make every element agree."

### x86_64 Haswell `sasum`, `dasum`: the blocks peel to an aligned address

`kernel/x86_64/sasum_microk_haswell-2.c` and `dasum_microk_haswell-2.c` sum
up to seven elements alone until `x` is 32-byte aligned, then group the rest
into their four vector sums, so the same data sums with other roundings
elsewhere in memory. Upstream's Skylake-X blocks dropped the peel for this
reason (their comment: "Deliberately no peeling").

- **Workaround:** rstd's port (`blas1/asum/x86_64/kernel.zig`) loads
  unaligned from the first element, as the Skylake-X blocks do.
- **Issue text:** "`sasum_microk_haswell-2.c` and `dasum_microk_haswell-2.c`
  peel to a 32-byte boundary before their AVX2 loop, so `sasum` and `dasum`
  on Haswell and Zen give results that depend on the buffer's address. The
  Skylake-X kernels removed their peel for this; the same change (unaligned
  loads from `x1[0]`) would make the Haswell ones a function of the data
  alone."

### arm64 `symv_U_asimd_4x4.c`: rows dropped when the first column is not a multiple of 4

The driver calls its block on rows `0..j1` with `j1 = (j/4)*4`, then takes
rows `j..j+3`; where `m - offset` is not a multiple of 4, rows `j1..j-1` of
those four columns are never added. `interface/symv.c` passes
`offset = m` and `symv_thread.c` splits on multiples of 4, so OpenBLAS never
reaches it; a direct kernel call can. `symv_U_sve_v1x4.c` calls its block on
`0..j`.

- **Workaround:** rstd's port (`blas2/symv/kernel.zig`, `ArmU`) adds those
  rows after the block, as the x86_64 driver's `kernel_1x4` does.
- **Issue text:** "In `kernel/arm64/symv_U_asimd_4x4.c`, `j1 = (j / 4) * 4`
  rounds the column start down and `symv_kernel_4x4(0, j1, ...)` is followed
  directly by the diagonal block at `j`, so for `m - offset` not a multiple of
  4 the rows `j1..j-1` miss the contribution of columns `j..j+3`. Unreachable
  through the interface today; a loop over `j1..j`, as `dsymv_U.c` on x86_64
  has, would make the kernel correct for any `offset`."

### x86_64 `ssymv_U_microk_sandy-2.c`, `ssymv_U_microk_haswell-2.c`: floats added as a double

Both blocks end with `vaddsd` (double) where their sums are floats: the low
64 bits of each reduced vector, two copies of the sum, plus `temp2` read with
`vmovss`. `temp2` is 0 whenever the drivers call them, so the low float, the
one stored, is the sum unchanged.

- **Workaround:** none needed; the port stores the sum (`Shape.store`).
- **Issue text:** "`ssymv_U_microk_sandy-2.c` and `ssymv_U_microk_haswell-2.c`
  add `temp2` with `vaddsd`, a double add, to single-precision sums. It is
  harmless because `temp2` is always 0 on entry, but `vaddss` is what the
  other `ssymv` blocks use."

### Complex `omatcopy`, `imatcopy`, `gemmt`: multiplying by 1 makes NaN

The real `omatcopy` and `imatcopy` kernels copy for `alpha` 1 (and store
zeros for 0); the complex ones (`kernel/arm/zomatcopy_c*.c`,
`kernel/generic/zimatcopy_c*.c`) always multiply, so `alpha` `(1, 0)` turns
`(x, Inf)` into `(NaN, Inf)`: `0 * Inf` in the other part. Two interfaces
multiply by `(1, 0)` themselves:

- `interface/zimatcopy.c` with `lda != ldb` (its `NEW_IMATCOPY` in-place path,
  defined by default, needs them equal) copies
  `alpha*op(A)` to a buffer and back with `OMATCOPY_K_CN` at `(1, 0)`, so an
  infinite part of the result makes the other NaN whatever `alpha` was, and
  `zimatcopy` differs from `zomatcopy` there;
- `interface/gemmt.c` conjugates the caller's `B` in place with
  `IMATCOPY_K_CNC` at `(1, 0)` for `transb` `R` or `C`, and back after, so an
  infinite part leaves `B` with NaN beside it; it also writes an input the
  caller may share. <https://github.com/OpenMathLib/OpenBLAS/pull/6060>
  (2026-09-29, after our pin) replaced this `gemmt` with a level-3 driver;
  recheck before filing.

- **Workaround:** rstd's matrix copies copy for `alpha` 1 and store zeros for
  0, real or complex (`ext/matcopy.zig`), and its `gemmt` reads `B` only
  (`blas3/driver.zig`).
- **Issue text:** "The complex `?omatcopy`/`?imatcopy` kernels multiply by
  `alpha` even when it is `(1, 0)`, so an infinite real or imaginary part
  turns the other into NaN (`0 * Inf`). `zimatcopy` with `lda != ldb` copies
  back through `OMATCOPY_K_CN` at `(1, 0)`, so it differs from `zomatcopy`
  for any `alpha`, and `?gemmt` with `transb` `R` or `C` conjugates `B` in
  place the same way, leaving NaN in the caller's `B`. A copy branch for
  `alpha == 1`, as the real kernels have, would fix all three; `gemmt` could
  also pass the conjugation to its `gemv` instead of writing `B`."

### `bgemm`: bfloat16 scalars compared by their bits, and `beta` truncated

Built `BGEMM`, `FLOAT` is `bfloat16`, a `uint16_t` (`common.h:285, 340-341`),
so the float tests on the scalars compare bits: `driver/level3/level3.c`'s
`alpha[0] == ZERO` holds only for +0 and `beta[0] != ONE` always holds, and
`kernel/generic/gemm_beta.c`'s `beta == ZERO` only for +0. A `-0` `alpha`
reads `A` and `B` (a NaN there reaches `C`), a `-0` `beta` multiplies (an
infinite `C` turns NaN), and `beta` 1 multiplies every element. That beta
kernel also narrows by taking the high half of the `f32` bits
(`f32tobfloat16`), a truncation, where the gemm kernel and every other
bfloat16 store round to nearest (`sbstobf16_`).

- **Workaround:** rstd tests the scalars by value, as its `gemvHalf` does, and
  its beta kernel rounds (`blas3/half.zig`, `ext/bf16/kernel.zig`).
- **Issue text:** "`bgemm` compares its bfloat16 `alpha` and `beta` as integers
  (`FLOAT` is `uint16_t`): `alpha == 0` and `beta == 0` miss `-0`, and `beta
  != 1` is always true, so `beta` 1 still rewrites `C`. `gemm_beta.c`'s
  bfloat16 path also truncates `beta*C` to bfloat16 instead of rounding as
  `sbstobf16` and the kernels do. Widening the scalars with `sbf16tos` before
  the tests, and rounding the product, would fix both."

### `potrf` reports success on a NaN pivot

`lapack/potf2/potf2_[UL].c` and `zpotf2_[UL].c`, which
`lapack/potrf/potrf_[UL]_single.c` calls for each diagonal block, test each pivot with
`ajj <= 0`, which is false for NaN, so a NaN on the diagonal (or one that
reaches it) passes, the factor fills with NaN, and `info` is 0. Reference
LAPACK's `xPOTRF2` checks `DISNAN` and returns the minor.

- **Workaround:** rstd's `potrf` fails a NaN pivot as the reference does
  (`src/linalg/lapack/cholesky.zig`).
- **Issue text:** "`?potrf` returns `info = 0` for a matrix whose pivot is
  NaN, because `ajj <= 0` is false for NaN; the factor is all NaN. Reference
  LAPACK's `xPOTRF2` tests `DISNAN(ajj)` and returns the minor; testing
  `!(ajj > 0)` would match it."

### `getf2` records a swap it does not make for a tiny or NaN pivot

`lapack/getf2/getf2_k.c` and `zgetf2_k.c` swap the pivot row in and scale the
column below only when `|p| >= DBL_MIN` (`FLT_MIN`; complex: either part),
yet set `ipiv` first. For a denormal or NaN pivot, columns `0..j` keep their
rows while the later columns and `getrf_single.c`'s interchanges take the
swap, and `L(:, j)` stays unscaled: `P L U` is not `A`, and `info` is 0.
Reference LAPACK's `xGETRF2` swaps for any nonzero pivot and divides by one
below the safe minimum. The threshold came in with
<https://github.com/OpenMathLib/OpenBLAS/pull/3924> (2023, against overflow in
`1/p`; then `1e-305`, now `DBL_MIN`/`FLT_MIN`).

- **Workaround:** rstd's `getrf` swaps and divides as the reference does
  (`src/linalg/lapack/lu.zig`).
- **Issue text:** "`?getrf` returns a wrong factorization when a pivot is
  nonzero but below `DBL_MIN` (`FLT_MIN` for single), or NaN: `getf2_k.c`
  records the interchange in `ipiv` but skips both the row swap and the
  scaling, so `P*L*U != A` with `info = 0`. Reference `xGETRF2` swaps for
  any nonzero pivot and divides element by element when `|p| < sfmin`;
  doing the same in the `else` branch would fix it."

### `potri` divides by a zero on the factor's diagonal

`interface/lapack/potri.c` and `zpotri.c` call `TRTRI_[UL]N_SINGLE` (or
`_PARALLEL`) directly (`potri.c:133`), skipping the zero-diagonal check that
`trtri.c` makes first (`:113`), so a singular factor gives infinities and NaN
with `info = 0`. Reference LAPACK's `xPOTRI` calls `xTRTRI`, which returns `INFO = i` for the
first zero `A(i, i)` before touching `A`. Found by reading the source.

- **Workaround:** rstd's `potri` checks the diagonal first and returns
  `error.Singular` (`src/linalg/lapack/cholesky.zig`).
- **Issue text:** "`?potri` returns `info = 0` and a matrix of infinities and
  NaN when the Cholesky factor has a zero on its diagonal: `potri.c` calls the
  `trtri` kernels without `trtri.c`'s `AMIN_K == 0` check. Reference `xPOTRI`
  returns `INFO = i` for the first zero `A(i, i)`."

### `trtrs` reports a bad `diag` before a bad `uplo` or `trans`

`interface/lapack/trtrs.c` and `ztrtrs.c` set `info = 3` for `diag` last
(`:115`), so a call with bad `uplo` and `diag` reports argument 3 where
Reference LAPACK's `xTRTRS` reports argument 1. Found by reading the source.

- **Workaround:** rstd's `trtrs` checks in the reference's order
  (`src/linalg/lapack/triangular.zig`).
- **Issue text:** "`?trtrs` checks `diag` after `uplo` and `trans`, so its
  `xerbla` names argument 3 where Reference LAPACK's `xTRTRS` names 1 or 2;
  moving the `diag` check above the `trans` one would match it."

### `zrotg.c`: the scaled branches take square roots where `zlartg` squares

`interface/zrotg.c` transcribes Reference LAPACK's `zlartg`, but its scaled
branches set `g2` and `f2` to `sqrt(re^2 + im^2)` (`:137`, `:143`, `:149`)
where `SRC/zlartg.f90` (v3.12.1, `:234`, `:243`, `:251`) has `ABSSQ`, the
squared norm, so `c`, `s` and `r` are wrong whenever an input is tiny or huge
enough to scale. Found by reading the source; not run against OpenBLAS.

- **Workaround:** rstd's `zrotg` ports `zlartg` itself
  (`blas1/interface.zig`), tested against a wide reference over every branch.
- **Issue text:** "`?rotg` for complex inputs (`interface/zrotg.c`) is wrong
  on its scaled paths: `g2` and `f2` are computed as `sqrt(x_r*x_r + x_i*x_i)`
  where LAPACK's `zlartg` uses the squared norm `ABSSQ`, so the rotation is
  wrong when `|a|` or `|b|` is outside `[rtmin, rtmax]`."

## ReLAPACK (OpenBLAS's `relapack/`)

At OpenBLAS `v0.3.34-317`, the same as upstream <https://github.com/HPAC/ReLAPACK> `master`
where the entry says so. The ground truth is a C harness that links ReLAPACK's sources with
naive BLAS loops and Reference-LAPACK's unblocked kernels (`lapack-netlib/SRC/*.c`), and checks
each result against the reference's routine or a residual (2026-10-30). The reference's own
`dsytrf` passes the same check.

### `?sytrf`, `?hetrf`: the factor is wrong past one level of recursion

`RELAPACK_dsytrf_rec2` (`relapack/src/dsytrf_rec2.c`, `xLASYF` less its update) ends by
putting its panel in LAPACK's form, undoing the panel's interchanges in its columns. The caller
(`dsytrf.c`, `RELAPACK_dsytrf_rec`) then updates the rest with those columns, which `xLASYF`
does before the undo, so the update reads rows in the wrong order. Deeper, the recursion defers
each half's update of the columns past its block, yet the pivot search and interchanges of the
next block read and move those columns. `n = 30` (two levels) already fails: `L D Lᵀ` is off by
up to 8 for entries of size 2; at `n = 100`, by 60. Every precision and both triangles, and the
same in upstream ReLAPACK.

- **Workaround:** rstd's `sytrf` and `hetrf` run Reference-LAPACK's blocked `xSYTRF` loop
  (`ILAENV`'s 64 columns) over ReLAPACK's `rec2` panel, its update ahead of the undo
  (`src/linalg/lapack/relapack.zig`, `indefinite.zig`).
- **Issue text:** "`RELAPACK_?sytrf` and `?hetrf` (and the `_rook` forms) return a wrong
  factorization for `n` past the crossover: `?sytrf_rec2` undoes the panel's interchanges
  before the recursion uses the panel to update `A_BR`, where `xLASYF` updates first; and the
  deferred updates of the columns right of each block are read stale by the next block's pivot
  search. `A - L*D*L**T` is O(1) at `n = 30`."

### `?sytrf_rook`, `?hetrf_rook`: the non-rook kernel below the crossover

`dsytrf_rook.c` and the other `_rook` forms call `LAPACK(dsytf2)` (`zhetf2`) for a small whole
matrix, not `dsytf2_rook`, so `ipiv` comes back in the non-rook format, which the rook solvers
read otherwise: at `n = 20`, `A - L D Lᵀ` read as rook is off by 6.

- **Workaround:** rstd's `sytrfRook` and `hetrfRook` call the rook kernel (`indefinite.tf2`).
- **Issue text:** "`RELAPACK_?sytrf_rook` and `?hetrf_rook` call `?sytf2` (`?hetf2`) instead
  of `?sytf2_rook` (`?hetf2_rook`) for `n <= CROSSOVER`, so `ipiv` is not in the rook format."

### `?gbtrf`: the recursion gives a wrong factor

OpenBLAS's `dgbtrf.c` takes `xGBTF2` unless `n <= kl` (`1036299da`, "due to crashes in LAPACK
tests"); where its recursion runs (`m = 50, n = 40, kl = 45, ku = 3`), 8 pivots and the factor
differ from `dgbtf2`'s, by up to 8. Upstream's full recursion differs on six of seven shapes
tried. The recursion also swaps `A_Rr` from its first column (`A_Rrj[i] = A_Rr[ip]`), and the
second half's `dgbtf2` overwrites the first half's `info`.

- **Workaround:** rstd's `gbtrf` is `xGBTF2` throughout (`relapack.zig`).
- **Issue text:** "`RELAPACK_?gbtrf` returns a wrong factorization whenever its recursion runs
  (`n <= kl` in OpenBLAS; most shapes upstream); compare with `?gbtf2` at `m = 50, n = 40,
  kl = 45, ku = 3`."

### `?pbtrf`: the rest factored as dense reads outside the band

`dpbtrf.c` factors the rest with `RELAPACK_dpotrf` whenever `kd > n1`, through the band's
storage taken as a dense matrix of leading dimension `ldAb - 1`; once the rest is wider than
`kd + 1` (`n = 100, kd = 50`), `potrf` reads and writes places of other columns. `zpbtrf.c`
also bounds `n22` by `kd` where `dpbtrf.c` takes `n1`. Found by reading the source.

- **Workaround:** rstd's `pbtrf` recurses on the band there (`relapack.zig`).
- **Issue text:** "`RELAPACK_?pbtrf` calls `?potrf` on `A_BR` as a dense matrix whenever
  `kd > n1`, which is outside the band storage once `n - n1 > kd + 1`."

### Smaller ones

- `?trtri` sets `INFO = i - 1` for a zero `A(i, i)` (`ztrtri.c`), so a zero first diagonal
  element reports success. rstd's `trtri` reports `i`.
- `?getrf` returns as soon as the left half finds a zero pivot (`if (*info) return`), leaving
  the factorization unfinished, where `xGETRF` completes it. rstd's completes it.
- `?sytrf` reports `INFO` as `info1 || info2`, 0 or 1, not the zero block's index. rstd's
  reports the index.
- `?trsyl`, `op(A)` transposed: the second solve's `scale` rescales `C_B` (`dlascl` on `&m1`
  rows of `C_B`), the half solved second, not `C_T`. rstd's rescales `C_T`.
- `?sygst`: `dsygst.c`, `ssygst.c` and `chegst.c` use their workspace only if `lWork > n2 * n1`,
  `zhegst.c` if `>=`; the top level's buffer is exactly `n1 * n2`, so three of four never use
  it there. rstd's takes `>=` for all.

These need no issue of their own beyond a line each in a ReLAPACK report.

## BLIS

### `snrm2` and `scnrm2` keep an `f32` running sum

BLIS 2.1's single-precision `nrm2` (Linux x86_64, `zen3` config) differs
from rstd by the same 87,202 eps as OpenBLAS's arm64 `scnrm2`, so it runs
the same sequential `f32` sum: about 1% off for 16M floats. Its `sasum` is
260 eps off there too, where rstd's converted kernel is 4.2 eps off.

- **Workaround:** none needed; rstd does not call BLIS.
- **Issue text:** "Single-precision `nrm2` (`snrm2`, `scnrm2`) loses about 1% for 16M floats on `zen3`, from a sequential `f32` sum of squares. Accumulate in `f64`, or in several partial sums."

## faer

faer v0.24.4 (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`), as
`src/linalg/lapack/faer/` ports it. Found by rstd's tests of the port, then
reproduced in a standalone copy of it; not run against faer, whose crates
are not fetched here.

### The divide-and-conquer SVD loses a root next to a pole

`secular_eq_root_finder` (`faer/src/linalg/svd/bidiag_svd.rs`) ends its
secant loop once `|f_cur - f_prev| <= eps` and returns the last iterate as
`mu`. For a root within rounding of a pole (a `col0` entry just above the
deflation threshold, so `mu` is about `z² / 2d`), `f` is 1 to rounding
everywhere the bracketing reached, so the loop stops at once, far from the
root; `perturb_col0` then makes `zhat` from that `mu`, and the vectors
reconstruct the bidiagonal with an error of `1e-4` (f64) where
`8 eps` was expected. Random bidiagonals graded over four decades hit it at
`n = 100` (f64) and `n = 64` (f32), threshold 4 or 8.

- **Workaround:** when the secant stops that way, unconverged, the bracket's
  bisection runs, as it does for faer's other exits
  (`bidiag_svd.zig`, `secant`). Graded, clustered and tiny-entry bidiagonals
  up to `n = 257` then reconstruct within 13 eps.
- **Issue text:** "The bidiagonal divide and conquer (`secular_eq_root_finder`)
  can return a root far from the true one when it lies within rounding of a
  pole: the secant loop exits on `|f_cur - f_prev| <= eps` without
  converging and without setting `use_bisection`, so `perturb_col0` builds
  `zhat` from a wrong `mu` and `U S Vᵀ` misses the bidiagonal by `1e-4`
  (f64). Falling back to bisection on that exit fixes it."

### Reciprocals of tiny norms overflow

`bidiag_in_place` (`svd/bidiag.rs`) scales the row `A12` by `norm.recip()`
whenever `norm != 0`, `qr_algorithm` and `divide_and_conquer`
(`svd/bidiag_svd.rs`) scale the bidiagonal by `max.recip()`, and
`compute_bidiag_cplx_svd`'s `normalized` (`svd/mod.rs`) scales an entry by
the reciprocal of its larger part. Below the safe
minimum the reciprocal is infinite, and the row or bidiagonal fills with NaN;
the SVD then reports `NoConvergence`. A rank-1 `f32` matrix (all ones,
`15 x 10`) gets there: its bidiagonal decays by `eps` per step to denormals.

- **Workaround:** below the safe minimum, rstd divides instead
  (`bidiag.zig`, `bidiag_svd.zig`'s `scaleDown`, `svd.zig`'s `normalized`);
  the bits elsewhere are faer's. The eigensolver has the same two:
  `tridiag_evd.rs`'s `qr_algorithm` scales by `max.recip()` and
  `evd/mod.rs`'s `normalized` closure by `x.abs().recip()`; rstd divides
  there too (`tridiag_evd.zig`'s `scaleDown`, `evd.zig`'s `normalized`).
- **Issue text:** "The SVD of a rank-deficient `f32` matrix (e.g. `15 x 10`
  all ones) returns `NoConvergence`: the bidiagonal decays into denormals,
  and `bidiag_in_place` and the bidiagonal solvers multiply by the
  reciprocal of a denormal norm or maximum, which is infinite. Dividing when
  the value is below `min_positive` avoids it."

### `divide_and_conquer`'s zero case leaves the last row's 1 out

For a subproblem whose bidiagonal is all zero, `divide_and_conquer` sets
`MatU::TwoRows` to the identity's first and last rows but stores
`u[(1, n)] = 0` where the identity has 1. When that subproblem is the second
half of a first half, the 0 becomes the parent's coupling `lambda`, and the
values-only path (`TwoRowsStorage`) drifts from the full one (1.5e-13 at
`n = 24`, zero block at 7..12).

- **Workaround:** rstd stores 1 (`bidiag_svd.zig`); both paths then agree.
- **Issue text:** "In `divide_and_conquer`'s `max == 0` branch,
  `MatU::TwoRows` stores `u[(1, n)] = zero()`; the last row of the identity
  has `one()` there. The values-only SVD then differs from the full one when
  a zero block ends a left half."

### The tall path leaves part of `U` as found

`svd` (`svd/mod.rs`) on a matrix taller than `11 / 6` of its width zeroes
`U`'s rows past `n` and sets the trailing diagonal, but never writes rows
`0..n` of columns `n..m` of a full `U` before applying `Q`. Faer's tests
pass zeroed matrices, so it only shows on a caller's own buffer.

- **Workaround:** rstd zeroes that block (`svd.zig`).
- **Issue text:** "`linalg::svd::svd` with a full `U` and `nrows > 11/6 ncols`
  does not initialize `U[0..n, n..m]` before applying the Householder
  sequence, so a caller passing an uninitialized `U` gets garbage in those
  columns."

### The staircase QR skips a NaN column

`qr_in_place_unblocked` (`qr/no_pivoting/factor.rs`) applies a column's
reflection only when `norm > threshold`. A NaN in the column makes that
false, so the column takes no reflection and the row does not advance: a
NaN matrix factors as rank 0 with finite pieces, and a least-squares solve
over it reports a zero pivot instead of returning NaN, as the reference's
`xGEQRF` does.

- **Workaround:** rstd tests `!(norm <= threshold)` (`qr.zig`), so NaN
  takes the reflection and propagates.
- **Issue text:** "`qr_in_place_unblocked` decides whether to reflect a
  column with `norm > threshold`, which is false when the column holds a
  NaN; the column is then skipped as if dependent, and the factorization of
  a NaN matrix comes out rank 0 and finite. `!(norm <= threshold)` keeps the
  NaN visible."

### Complex `rotg` returns `r = 1` for a zero second entry

`JacobiRotation::rotg` (`jacobi.rs`), complex branch, returns
`(c = 1, s = 0, r = 1)` when `q = 0`; the rotation leaves `p` unchanged, so
`r` must be `p` (the reference's `xLARTG` returns `r = f`). The complex
`lahqr` (`evd/schur/complex_schur.rs`) stores that `r` in the bulge
chase's `A(i, i - 1)`, so a bulge that is exactly zero (an underflowed
product) would replace a subdiagonal entry with 1. The real branch returns
`p`.

- **Workaround:** rstd returns `p` (`householder.zig`, `ComplexRotation`).
- **Issue text:** "The complex branch of `JacobiRotation::rotg` returns
  `r = one()` when `q == 0`, but the identity rotation it returns maps `p`
  to `p`; `r` should be `p`, as LAPACK's `zlartg` has it. `lahqr` writes
  this `r` into the Hessenberg matrix during the bulge chase."

### The eigenvector solve guards a determinant with a linear bound

`solve_cplx_shifted_upper_quasi_triangular_system` (`evd/mod.rs`) solves
each 2 x 2 block of the real Schur form shifted by a complex eigenvalue
through its determinant `z = (a - p - i q)² - b c`, quadratic in `A`, and
replaces it by `eps ‖A‖₁` whenever `|z| < eps ‖A‖₁`, a linear bound. At
unit scale that is rare, but a matrix scaled small (as `dgeev.f`'s scaling
to `SMLNUM` leaves one) always has `|z| < eps ‖A‖₁`, so the vector of every
complex pair whose solve meets another 2 x 2 block comes out wrong; scaled
large, the guard never acts. A 12 x 12 random `f64` scaled by `2^-958` gave
a left eigenvector with residual 1.15.

- **Workaround:** rstd guards with `(eps ‖A‖₁)²` (`evd.zig`,
  `solveQuasiComplex`), so the test is scale invariant.
- **Issue text:** "In `solve_cplx_shifted_upper_quasi_triangular_system`'s
  2 x 2 branch, `z` is the shifted block's determinant (quadratic in the
  matrix) but is compared with `eps * norm` (linear), so the guard depends on
  the matrix's scale: a small-scaled matrix always clamps and gets wrong
  eigenvectors. Comparing with `(eps * norm)²` fixes it."

### The pivoted QR scales a zero matrix by an infinite reciprocal

`qr_in_place_unblocked` (`faer/src/linalg/qr/col_pivoting/factor.rs`) scales
`A` by `best.recip()`, `best` its largest column norm, with no check for zero:
a zero matrix becomes NaN (`0 * inf`), and so do its `R` and `tau`s. Found
reading the source while porting `geqp3`.

- **Workaround:** a zero largest norm leaves `A` unscaled (`qr_pivot.zig`);
  the zero-matrix case in its tests covers it.
- **Issue text:** "Column-pivoted `qr_in_place` on an all-zero matrix scales
  it by `1 / 0` and returns NaN in `R` and the Householder factors; skipping
  the scaling when the largest column norm is zero fixes it."
