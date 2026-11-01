//! `qemu-test` step: runs rstd's tests on emulated machines with
//! `qemu-system-*`. Needs `qemu-system-<arch>` on `PATH`. A machine boots one
//! of three ways (`Boot`):
//!
//!   * `.kernel`: a freestanding test kernel per architecture (boot code in
//!     `test/qemu/`, runner in `test/qemu/runner.zig`, test root in
//!     `src/_test_qemu.zig`). Where the machine has semihosting (aarch64, arm,
//!     riscv64, MIPS) the runner exits with the number of failed tests, so a
//!     failing test fails the step. Elsewhere the machine just powers off and the
//!     step checks the printed summary for "0 failed". The console is also
//!     logged as it runs, to `.zig-cache/qemu-test-<name>.log`. A UEFI target
//!     builds the same tests as an application (`test/qemu/uefi.zig`), which
//!     QEMU's own OVMF or AAVMF boots from a FAT drive.
//!   * `.linux`: rstd's full test suite as a static Linux executable, in an
//!     initramfs whose init (`test/linux/init.zig`) mounts procfs, sysfs
//!     and devtmpfs, runs it and reports. The kernel is a Debian netboot one,
//!     a lazy dependency fetched only when rstd is the project being built
//!     (`-Dlinux-vm-off` skips it). No disk, network or
//!     setup; the console is logged live like the kernel rigs'.
//!   * `.macos`: rstd's unit tests in an installed macOS, which the Mach-O
//!     kernels need, copied in and run over SSH. It boots with OSX-KVM's
//!     OpenCore and firmware, a lazy dependency fetched only when rstd is the
//!     project being built (`-Dmacos-vm-off` skips it), so a dependent never
//!     fetches it. The VM lives in `.macos-vm/<name>/` (gitignored).
//!
//! A macOS machine, once:
//!
//!   1. `zig build qemu-setup-<name>`: fetches a recovery image
//!      (`-Dmac-release`, Sonoma by default) and makes the disk.
//!   2. `zig build qemu-install-<name>`, then in the window: "macOS Base
//!      System" in the OpenCore picker; Disk Utility, erase the QEMU HARDDISK
//!      as APFS; reinstall macOS onto it (pick "macOS Installer", then the
//!      disk, at each reboot); create an account; run `tools/configure_macos.sh`
//!      in it (its header says how), which makes the disk the default boot,
//!      turns on Remote Login and turns off what slows an emulated Mac.
//!   3. Authorize an SSH key (`ssh-keygen` makes one if you have none):
//!      `type %USERPROFILE%\.ssh\id_ed25519.pub | ssh -p 2222 <account>@localhost
//!      "mkdir -p ~/.ssh && cat >> ~/.ssh/authorized_keys"`.
//!
//! Then `zig build qemu-start-<name>` boots it in a window, and, while it is
//! up, `zig build qemu-test-<name> -Dmac-user=<account>` runs the tests in it.
//! With the VM shut down, `qemu-snapshot-<name>` saves a restore point (both
//! disks and the NVRAM) and `qemu-restore-<name>` goes back to it
//! (`-Dmac-snapshot=<name>`, default `clean`): worth doing once it is set
//! up, since a hard stop of the VM can leave macOS unable to boot.

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const console = rstdbuild.console;
const blaskern = @import("blaskern.zig");
const buildroot = @import("__root__.zig");
const sshtest = @import("sshtest.zig");

const Machine = struct {
    /// Name for steps and executables.
    name: []const u8,
    query: std.Target.Query,
    /// `qemu-system-<arch>` binary.
    qemu: []const u8,
    /// Acceleration, machine, CPU and devices.
    args: []const []const u8,
    /// Overrides the build's optimize mode; OpenBIOS maps about 4 MB of a sparc64
    /// image, which debug builds of the tests outgrow. sparc64 uses `.fast`, as
    /// `.small` from Zig 0.18.0-dev.250 jumps into unrelated code in the gemm
    /// parameter test (an LLVM `-Oz` miscompile on sparc64, it appears).
    optimize: ?Optimize = null,
    memory: []const u8 = "128M",
    boot: Boot,
    /// Why this rig is not wired up yet; a stub gets no step.
    stub: ?[]const u8 = null,
    /// Runs under WHPX (`args` has the `-accel`), so only on an x86_64
    /// Windows host, and never on `-Dtest-host`.
    whpx: bool = false,
    /// `args` under HVF, on a `-Dtest-host` that is an arm64 Mac, unless
    /// `-Dtest-host-tcg`; null where the rig needs a CPU the host lacks.
    hvf: ?[]const []const u8 = null,
    /// The BLAS core a Linux rig must select under `args`, and under `hvf`
    /// (on the M6); its run fails otherwise.
    core: ?[]const u8 = null,
    hvf_core: ?[]const u8 = null,
    /// Kernel objects linked: `.none` for a soft-float ABI, which no
    /// converted kernel serves (`abi.zig`).
    blas_cores: rstdbuild.BlasCores = .all,
    /// `XDouble` as `f128` (`-Dquad-precision` for this rig alone).
    quad_precision: bool = false,
    /// Test sets each run as a kernel of its own, `<name>-<part>`, for a
    /// machine whose firmware maps too little for all the tests in one image;
    /// `-Dqemu-test-part` picks some, and `-Dqemu-test-filter` builds one
    /// image as usual.
    parts: []const Part = &.{},
    /// CPUs: past 1 the kernel builds multithreaded, starts them all
    /// (`test/qemu/smp.zig`) and runs `smp_filters` unless filtered.
    cpus: u8 = 1,
};

const Part = struct { name: []const u8, filters: []const []const u8 };

/// How the tests reach a machine.
const Boot = union(enum) {
    kernel: Kernel,
    /// An installed macOS, set up once (see the file header).
    macos,
    /// A Linux kernel booted with an initramfs holding rstd's full test suite.
    linux: Linux,

    /// A freestanding test kernel, booted bare with `-kernel`.
    const Kernel = struct {
        linker_script: []const u8,
        /// Use the LLVM backend instead of the self-hosted one for code generation.
        llvm: bool = false,
        code_model: std.lang.CodeModel = .default,
        /// Exit status carries the result; otherwise the printed summary does.
        semihosting: bool = true,
        /// Text qemu always prints to stderr here. The build treats any unchecked
        /// stderr as a diagnostic (a "failed command" echo on a passing step), so
        /// expecting it keeps a passing run quiet.
        stderr_note: ?[]const u8 = null,
        /// The addresses the firmware maps for the image, if it maps too few
        /// for every image: one whose loaded span is past them fails the build.
        load_window: ?Window = null,
        /// A UEFI target's firmware, by its name in QEMU's data directory,
        /// where Nix's QEMU ships it; no linker script then. Test host only:
        /// `tools/nix/run.sh` finds that directory.
        firmware: ?[]const u8 = null,
    };

    const Window = struct { base: u64, size: u64 };

    /// A Debian netboot kernel (a lazy dependency) and an initramfs of
    /// `test/linux/init.zig` as PID 1 and the hosted test executable.
    const Linux = struct {
        /// The `build.zig.zon` dependency, and the kernel's path in it.
        dep: []const u8,
        kernel: []const u8,
        /// The kernel's console, on the machine's first serial port.
        console: []const u8,
    };
};

const machines = base_machines ++ x87Safe(&.{ "i386", "i386-max" }) ++ quad(&.{"aarch64"});

/// `names` again as `<name>-safe`, in ReleaseSafe (T87): Debug stores every
/// value to memory, which hides x87 width bugs (T56's showed in ReleaseSafe only).
fn x87Safe(comptime names: []const []const u8) [names.len]Machine {
    var out: [names.len]Machine = undefined;
    for (names, &out) |name, *o| {
        o.* = for (base_machines) |m| {
            if (std.mem.eql(u8, m.name, name)) break m;
        } else @compileError("no rig " ++ name);
        o.name = name ++ "-safe";
        o.optimize = .safe;
    }
    return out;
}

/// `names` again as `<name>-quad`, with `XDouble` as `f128` (T19): the edge
/// and level-3 driver tests' only extended type off x86.
fn quad(comptime names: []const []const u8) [names.len]Machine {
    var out: [names.len]Machine = undefined;
    for (names, &out) |name, *o| {
        o.* = for (base_machines) |m| {
            if (std.mem.eql(u8, m.name, name)) break m;
        } else @compileError("no rig " ++ name);
        o.name = name ++ "-quad";
        o.quad_precision = true;
    }
    return out;
}

const base_machines = [_]Machine{
    .{
        .name = "aarch64",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .freestanding, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.cortex_a72 } },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a72" },
        // no semihosting: the boot code falls back to the PL011 and PSCI
        .hvf = &.{ "-accel", "hvf", "-M", "virt", "-cpu", "host" },
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/aarch64.ld",
        } },
    },
    .{
        // the threaded tests on four cores, started by PSCI
        .name = "aarch64-smp",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .freestanding, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.cortex_a72 } },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a72" },
        .hvf = &.{ "-accel", "hvf", "-M", "virt", "-cpu", "host" },
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/aarch64.ld",
        } },
        .cpus = 4,
    },
    .{
        // TCG on a Mac too, which has no SVE
        .name = "aarch64-sve",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .freestanding, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.neoverse_v1 } },
        .qemu = "qemu-system-aarch64",
        // max has SVE2, SME2 and BF16, so every core but Apple's runs
        .args = &.{ "-M", "virt", "-cpu", "max,pauth-impdef=on" },
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/aarch64.ld",
        } },
    },
    .{
        .name = "arm",
        // the A15's FPU and divide: Zig's model has neither (LLVM adds them
        // from its target parser), so the bare model is soft-float
        .query = .{
            .cpu_arch = .arm,
            .os_tag = .freestanding,
            .abi = .eabihf,
            .cpu_model = .{ .explicit = &std.Target.arm.cpu.cortex_a15 },
            .cpu_features_add = std.Target.arm.featureSet(&.{ .vfp4, .neon, .d32, .hwdiv, .hwdiv_arm }),
        },
        .qemu = "qemu-system-arm",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a15" },
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/arm.ld",
        } },
    },
    .{
        // the soft-float ABI over the A15's FPU (softfp): floats in integer
        // registers, the kernels' `sf` variant; the features pinned as in
        // `arm`, since Zig 0.18's model has none and 0.17's had them
        .name = "arm-softfp",
        .query = .{
            .cpu_arch = .arm,
            .os_tag = .freestanding,
            .abi = .eabi,
            .cpu_model = .{ .explicit = &std.Target.arm.cpu.cortex_a15 },
            .cpu_features_add = std.Target.arm.featureSet(&.{ .vfp4, .neon, .d32, .hwdiv, .hwdiv_arm }),
        },
        .qemu = "qemu-system-arm",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a15" },
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/arm.ld",
        } },
    },
    .{
        .name = "riscv64",
        .query = .{ .cpu_arch = .riscv64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-riscv64",
        // the vector extension, which the RVV kernels need, at the VLEN of
        // the widest linked core (x280's 512)
        .args = &.{ "-M", "virt", "-cpu", "rv64,v=true,vlen=512", "-smp", "1", "-bios", "none" },
        .boot = .{
            .kernel = .{
                .linker_script = "test/qemu/riscv64.ld",
                // The kernel sits at 0x80000000, out of reach of the default `medlow`.
                .code_model = .medany,
            },
        },
    },
    .{
        // `riscv64` at VLEN 128 and 256: `zvl128b`'s and `zvl256b`'s objects
        // on the hardware they are named for, the VLEN-agnostic `axpy` at both
        .name = "riscv64-vlen128",
        .query = .{ .cpu_arch = .riscv64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-riscv64",
        .args = &.{ "-M", "virt", "-cpu", "rv64,v=true,vlen=128", "-smp", "1", "-bios", "none" },
        .boot = .{ .kernel = .{ .linker_script = "test/qemu/riscv64.ld", .code_model = .medany } },
    },
    .{
        .name = "riscv64-vlen256",
        .query = .{ .cpu_arch = .riscv64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-riscv64",
        .args = &.{ "-M", "virt", "-cpu", "rv64,v=true,vlen=256", "-smp", "1", "-bios", "none" },
        .boot = .{ .kernel = .{ .linker_script = "test/qemu/riscv64.ld", .code_model = .medany } },
    },
    .{
        .name = "loongarch64",
        .query = .{ .cpu_arch = .loongarch64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-loongarch64",
        .args = &.{ "-M", "virt", "-cpu", "la464" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/loongarch64.ld",
            .semihosting = false,
        } },
    },
    .{
        .name = "mips64r2el",
        .query = .{ .cpu_arch = .mips64el, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips64r2 } },
        .qemu = "qemu-system-mips64el",
        // with MSA, so the Loongson 3A R4 core runs; the big-endian rig keeps
        // the 3A1000, a release 2 machine without it
        .args = &.{ "-M", "malta", "-cpu", "Loongson-3A4000" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips64.ld",
        } },
    },
    .{
        .name = "mips64r6el",
        .query = .{ .cpu_arch = .mips64el, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips64r6 } },
        .qemu = "qemu-system-mips64el",
        .args = &.{ "-M", "malta", "-cpu", "I6400" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips64.ld",
        } },
    },
    .{
        .name = "mips32r2el",
        .query = .{ .cpu_arch = .mipsel, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips32r2 } },
        .qemu = "qemu-system-mipsel",
        .args = &.{ "-M", "malta", "-cpu", "24Kf" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips32.ld",
        } },
    },
    .{
        .name = "mips32r5el",
        .query = .{ .cpu_arch = .mipsel, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.p5600 } },
        .qemu = "qemu-system-mipsel",
        .args = &.{ "-M", "malta", "-cpu", "P5600" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips32.ld",
        } },
    },
    .{
        .name = "mips64r2",
        .query = .{ .cpu_arch = .mips64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips64r2 } },
        .qemu = "qemu-system-mips64",
        .args = &.{ "-M", "malta", "-cpu", "Loongson-3A1000" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips64.ld",
        } },
    },
    .{
        .name = "mips64r6",
        .query = .{ .cpu_arch = .mips64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips64r6 } },
        .qemu = "qemu-system-mips64",
        .args = &.{ "-M", "malta", "-cpu", "I6400" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips64.ld",
        } },
    },
    .{
        .name = "mips32r2",
        .query = .{ .cpu_arch = .mips, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips32r2 } },
        .qemu = "qemu-system-mips",
        .args = &.{ "-M", "malta", "-cpu", "24Kf" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips32.ld",
        } },
    },
    .{
        .name = "mips32r5",
        .query = .{ .cpu_arch = .mips, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.p5600 } },
        .qemu = "qemu-system-mips",
        .args = &.{ "-M", "malta", "-cpu", "P5600" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips32.ld",
        } },
    },
    .{
        .name = "i386",
        .query = .{ .cpu_arch = .x86, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.pentium3 } },
        .qemu = "qemu-system-i386",
        .args = &.{ "-M", "pc", "-cpu", "pentium3" },
        .memory = "128M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/x86.ld",
            .semihosting = false,
        } },
    },
    .{
        // the SSE2-and-up cores, which a `pentium3` cannot run; a CPU past the
        // baseline (`pentium4`) also builds each core's object for its own ISA
        .name = "i386-max",
        .query = .{ .cpu_arch = .x86, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.nehalem } },
        .qemu = "qemu-system-i386",
        .args = &.{ "-M", "pc", "-cpu", "max" },
        .memory = "128M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/x86.ld",
            .semihosting = false,
        } },
    },
    .{
        .name = "x86_64",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.haswell } },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "pc", "-cpu", "max" },
        .memory = "128M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/x86_64.ld",
            .llvm = true,
            .semihosting = false,
        } },
    },
    .{
        // the threaded tests on four cores, started by INIT-SIPI-SIPI
        .name = "x86_64-smp",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.haswell } },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "pc", "-cpu", "max" },
        .memory = "128M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/x86_64.ld",
            .llvm = true,
            .semihosting = false,
        } },
        .cpus = 4,
    },
    .{
        // The host CPU under Windows' hypervisor: TCG has no AVX-512, FMA4 or
        // XOP (issues.md), so their cores run only here, on a host that has
        // them. WHPX keeps only the model's features the host has; `host`
        // and `max` stop in SeaBIOS.
        .name = "x86_64-whpx",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.haswell } },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-accel", "whpx", "-M", "pc", "-cpu", "SapphireRapids,+fma4,+xop" },
        .whpx = true,
        .memory = "128M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/x86_64.ld",
            .llvm = true,
            .semihosting = false,
        } },
    },
    .{
        .name = "sparc64",
        .optimize = .fast,
        // OpenBIOS maps 0x400000 to 0x800000 only, and the default set's
        // code alone is about 5 MB; blas2 by itself is 3.3 MB
        .parts = &.{
            .{ .name = "blas1", .filters = &.{ "linalg.blas.blas1.", "linalg.blas._test.blas1." } },
            .{ .name = "blas2-dense", .filters = &.{ "linalg.blas.blas2.__root__.", "linalg.blas.blas2.gemv.", "linalg.blas.blas2.ger.", "linalg.blas.blas2.hemv.", "linalg.blas.blas2.symv." } },
            .{ .name = "blas2-banded", .filters = &.{ "linalg.blas.blas2.banded.", "linalg.blas.blas2.symmetric." } },
            .{ .name = "blas2-triangular", .filters = &.{"linalg.blas.blas2.triangular."} },
            .{ .name = "blas2-interface", .filters = &.{"linalg.blas.blas2.interface."} },
            .{ .name = "blas2-edge", .filters = &.{"linalg.blas._test.blas2.edge."} },
            // one per routine family: together past the window
            .{ .name = "blas3-gemm", .filters = &.{ "linalg.blas.blas3.__root__.", "linalg.blas.blas3.gemm." } },
            .{ .name = "blas3-symm", .filters = &.{ "linalg.blas.blas3.symm.", "linalg.blas.blas3.hemm." } },
            .{ .name = "blas3-rank", .filters = &.{ "linalg.blas.blas3.syrk.", "linalg.blas.blas3.syr2k." } },
            .{ .name = "blas3-tri", .filters = &.{ "linalg.blas.blas3.trmm.", "linalg.blas.blas3.trsm.", "linalg.lapack." } },
            .{ .name = "blas3-interface", .filters = &.{"linalg.blas.blas3.interface."} },
            .{ .name = "openblas", .filters = &.{"linalg.blas._test.from_openblas."} },
            .{ .name = "rest", .filters = &.{ "linalg.blas.__root__.", "linalg.blas.Config.", "linalg.blas.cpu.", "simd." } },
            .{ .name = "small", .filters = &small_filters },
        },
        .query = .{ .cpu_arch = .sparc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.sparc.cpu.v9 } },
        .qemu = "qemu-system-sparc64",
        .args = &.{ "-M", "sun4u" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/sparc64.ld",
            .semihosting = false,
            .load_window = .{ .base = 0x400000, .size = 0x400000 },
        } },
    },
    .{
        // sparc64 in debug, the mode the other rigs build in; its images are
        // larger, so blas3 is split too, and gemm's three test files apart
        .name = "sparc64-debug",
        .optimize = .debug,
        .parts = &.{
            // blas1 in three, as two are past the window
            .{ .name = "blas1", .filters = &.{ "linalg.blas.blas1.__root__.", "linalg.blas.blas1.amax.", "linalg.blas.blas1.amin.", "linalg.blas.blas1.asum.", "linalg.blas.blas1.axpby.", "linalg.blas.blas1.axpy.", "linalg.blas._test.blas1.check.", "linalg.blas.blas1.copy.", "linalg.blas.blas1.dot.", "linalg.blas.blas1.extremum.", "linalg.blas.blas1.iamax.", "linalg.blas.blas1.iamin.", "linalg.blas.blas1.imax.", "linalg.blas.blas1.imin.", "linalg.blas.blas1.max.", "linalg.blas.blas1.min.", "linalg.blas.blas1.nrm2.", "linalg.blas.blas1.rot.", "linalg.blas.blas1.rotm.", "linalg.blas.blas1.scal.", "linalg.blas.blas1.sum.", "linalg.blas.blas1.swap.", "linalg.blas.blas1.uk_block." } },
            .{ .name = "blas1-interface", .filters = &.{"linalg.blas.blas1.interface."} },
            .{ .name = "blas1-edge", .filters = &.{"linalg.blas._test.blas1.edge."} },
            .{ .name = "blas1-denormal", .filters = &.{"linalg.blas._test.blas1.denormal."} },
            // the dense routines, the triangular ones, the interface and the
            // complex edge cases each in two, as one is past the window
            .{ .name = "blas2-dense", .filters = &.{ "linalg.blas.blas2.__root__.", "linalg.blas.blas2.gemv.__root__.test.portable", "linalg.blas.blas2.gemv.__root__.test_" } },
            .{ .name = "blas2-gemv-x87", .filters = &.{ "linalg.blas.blas2.gemv.__root__.test.half-precision gemv:", "linalg.blas.blas2.gemv.__root__.test.extended gemv:" } },
            .{ .name = "blas2-gemv-cores", .filters = &.{"linalg.blas.blas2.gemv.__root__.test.every linked core's"} },
            .{ .name = "blas2-dense-2", .filters = &.{ "linalg.blas.blas2.ger.", "linalg.blas.blas2.hemv.", "linalg.blas.blas2.symv." } },
            .{ .name = "blas2-banded", .filters = &.{
                "linalg.blas.blas2.banded.test.banded",
                "linalg.blas.blas2.banded.test.the small",
                "linalg.blas.blas2.banded.test.strided",
                "linalg.blas.blas2.banded.test.an allocation",
                "linalg.blas.blas2.banded.test.invalid",
                "linalg.blas.blas2.symmetric.",
            } },
            .{ .name = "blas2-triangular", .filters = &.{"linalg.blas.blas2.triangular.test.triangular"} },
            // every core's triangular results in three, as one is past the window
            .{ .name = "blas2-triangular-cores-real", .filters = &.{"linalg.blas.blas2.triangular.test.every linked core's kernels give the same triangular results, real"} },
            .{ .name = "blas2-triangular-cores-complex", .filters = &.{"linalg.blas.blas2.triangular.test.every linked core's kernels give the same triangular results, complex products"} },
            .{ .name = "blas2-triangular-cores-complex-2", .filters = &.{"linalg.blas.blas2.triangular.test.every linked core's kernels give the same triangular results, complex solves"} },
            .{ .name = "blas2-f80", .filters = &.{ "linalg.blas.blas2.banded.test.f80", "linalg.blas.blas2.triangular.test.f80" } },
            .{ .name = "blas2-interface", .filters = &.{
                "linalg.blas.blas2.interface.test.interface: gemv",
                "linalg.blas.blas2.interface.test.interface: ger",
                "linalg.blas.blas2.interface.test.interface: symv",
                "linalg.blas.blas2.interface.test.interface: problems",
            } },
            .{ .name = "blas2-interface-2", .filters = &.{
                "linalg.blas.blas2.interface.test.interface: the half",
                "linalg.blas.blas2.interface.test.interface: half",
                "linalg.blas.blas2.interface.test.every linked core's",
                "linalg.blas.blas2.interface.test_",
            } },
            // f16, f80 and f128: past the window with the rest
            .{ .name = "blas2-interface-any", .filters = &.{"linalg.blas.blas2.interface.test.interface: level 2"} },
            // real, complex, and complex on every core in two, as one is past
            // the window
            .{ .name = "blas2-edge-real", .filters = &.{ "linalg.blas._test.blas2.edge.test.portable level-2 routines, real", "linalg.blas._test.blas2.edge.test.every linked core's level-2 routines, real" } },
            .{ .name = "blas2-edge-complex", .filters = &.{"linalg.blas._test.blas2.edge.test.portable level-2 routines, complex"} },
            .{ .name = "blas2-edge-complex-cores", .filters = &.{"linalg.blas._test.blas2.edge.test.every linked core's level-2 routines, complex products"} },
            .{ .name = "blas2-edge-complex-cores-2", .filters = &.{"linalg.blas._test.blas2.edge.test.every linked core's level-2 routines, complex solves"} },
            .{ .name = "blas3-gemm", .filters = &.{"linalg.blas.blas3.gemm.__root__."} },
            // but "every core of" each family: those are the blas3-params-* parts
            .{ .name = "blas3-params", .filters = &.{
                "linalg.blas.blas3.gemm.params.test.dtb entries",
                "linalg.blas.blas3.gemm.params.test.generic follows",
                "linalg.blas.blas3.gemm.params.test.the block factor",
                "linalg.blas.blas3.gemm.params.test.fixed sets",
                "linalg.blas.blas3.gemm.params.test.old x86",
                "linalg.blas.blas3.gemm.params.test.LA464",
                "linalg.blas.blas3.gemm.params.test.detection gives",
            } },
            .{ .name = "blas3-params-x86_64", .filters = &.{"linalg.blas.blas3.gemm.params.test.every core of the x86_64 family"} },
            .{ .name = "blas3-params-x86", .filters = &.{"linalg.blas.blas3.gemm.params.test.every core of the x86 family"} },
            .{ .name = "blas3-params-arm", .filters = &.{ "linalg.blas.blas3.gemm.params.test.every core of the arm64 family", "linalg.blas.blas3.gemm.params.test.every core of the arm family" } },
            .{ .name = "blas3-params-power", .filters = &.{ "linalg.blas.blas3.gemm.params.test.every core of the power64 family", "linalg.blas.blas3.gemm.params.test.every core of the power32 family" } },
            .{ .name = "blas3-params-rest", .filters = &.{ "linalg.blas.blas3.gemm.params.test.every core of the riscv64 family", "linalg.blas.blas3.gemm.params.test.every core of the s390x family", "linalg.blas.blas3.gemm.params.test.every core of the mips64 family", "linalg.blas.blas3.gemm.params.test.every core of the mips32 family", "linalg.blas.blas3.gemm.params.test.every core of the loongarch64 family", "linalg.blas.blas3.gemm.params.test.every core of the sparc family", "linalg.blas.blas3.gemm.params.test.every core of the wasm family", "linalg.blas.blas3.gemm.params.test.every core of the generic family" } },
            .{ .name = "blas3-kernel", .filters = &.{
                "linalg.blas.blas3.__root__.",
                "linalg.blas.blas3.gemm.kernel.test.packed kernels",
                "linalg.blas.blas3.gemm.kernel.test.geometry",
                "linalg.blas.blas3.gemm.kernel.test.beta",
            } },
            // the complex kernels in three, as the 8-lane tile's alone are
            // past the window
            .{ .name = "blas3-kernel-complex", .filters = &.{"linalg.blas.blas3.gemm.kernel.test.complex packed kernels match"} },
            .{ .name = "blas3-kernel-complex-wide", .filters = &.{"linalg.blas.blas3.gemm.kernel.test.complex packed kernels of an 8-lane tile match plain arithmetic, A plain"} },
            .{ .name = "blas3-kernel-complex-wide-2", .filters = &.{"linalg.blas.blas3.gemm.kernel.test.complex packed kernels of an 8-lane tile match plain arithmetic, A conjugated"} },
            .{ .name = "blas3-symm", .filters = &.{ "linalg.blas.blas3.symm.", "linalg.blas.blas3.hemm." } },
            .{ .name = "blas3-syrk", .filters = &.{"linalg.blas.blas3.syrk."} },
            .{ .name = "blas3-syr2k", .filters = &.{"linalg.blas.blas3.syr2k."} },
            .{ .name = "blas3-trmm", .filters = &.{"linalg.blas.blas3.trmm."} },
            .{ .name = "blas3-trsm", .filters = &.{"linalg.blas.blas3.trsm."} },
            // the interface in five, as two are near the window
            .{ .name = "blas3-interface", .filters = &.{ "linalg.blas.blas3.interface.test.interface: gemm", "linalg.blas.blas3.interface.test.interface: a matrix past" } },
            .{ .name = "blas3-interface-symm", .filters = &.{"linalg.blas.blas3.interface.test.interface: symm"} },
            .{ .name = "blas3-interface-rank", .filters = &.{ "linalg.blas.blas3.interface.test.interface: syrk", "linalg.blas.blas3.interface.test.interface: the rank-k" } },
            .{ .name = "blas3-interface-tri", .filters = &.{"linalg.blas.blas3.interface.test.interface: trmm"} },
            .{ .name = "blas3-interface-tri-unslotted", .filters = &.{"linalg.blas.blas3.interface.test.interface: unslotted trmm"} },
            // the OpenBLAS suites in three
            .{ .name = "openblas-l1", .filters = &.{ "linalg.blas._test.from_openblas.blat1.", "linalg.blas._test.from_openblas.utest.", "linalg.blas._test.from_openblas.utest_ext." } },
            // blat2 in two, past the window with level 3's kernels
            .{ .name = "openblas-blat2", .filters = &.{ "linalg.blas._test.from_openblas.blat2.test_", "linalg.blas._test.from_openblas.blat2.test.blat2: SMVCH", "linalg.blas._test.from_openblas.blat2.test.blat2: xGEMV", "linalg.blas._test.from_openblas.blat2.test.blat2: xGBMV", "linalg.blas._test.from_openblas.blat2.test.blat2: xSYMV", "linalg.blas._test.from_openblas.blat2.test.blat2: xSBMV", "linalg.blas._test.from_openblas.blat2.test.blat2: xSPMV", "linalg.blas._test.from_openblas.blat2.test.blat2: xTRMV", "linalg.blas._test.from_openblas.blat2.test.blat2: xTBMV", "linalg.blas._test.from_openblas.blat2.test.blat2: xTPMV", "linalg.blas._test.from_openblas.blat2.test.blat2: error exits" } },
            .{ .name = "openblas-blat2-b", .filters = &.{ "linalg.blas._test.from_openblas.blat2.test.blat2: xTRSV", "linalg.blas._test.from_openblas.blat2.test.blat2: xTBSV", "linalg.blas._test.from_openblas.blat2.test.blat2: xTPSV", "linalg.blas._test.from_openblas.blat2.test.blat2: xGER,", "linalg.blas._test.from_openblas.blat2.test.blat2: xGERC", "linalg.blas._test.from_openblas.blat2.test.blat2: xSYR,", "linalg.blas._test.from_openblas.blat2.test.blat2: xSPR,", "linalg.blas._test.from_openblas.blat2.test.blat2: xSYR2", "linalg.blas._test.from_openblas.blat2.test.blat2: xSPR2" } },
            .{ .name = "openblas-blat3", .filters = &.{"linalg.blas._test.from_openblas.blat3."} },
            // the suites on every linked core, by level, level 1 by suite
            // (utest_ext's level 3 runs with level 3), level 2 by type
            .{ .name = "openblas-cores-l1-blat1", .filters = &.{"linalg.blas._test.from_openblas.__root__.test.the level-1 suites, blat1"} },
            .{ .name = "openblas-cores-l1-utest", .filters = &.{"linalg.blas._test.from_openblas.__root__.test.the level-1 suites, utest,"} },
            .{ .name = "openblas-cores-l1-utest-2", .filters = &.{"linalg.blas._test.from_openblas.__root__.test.the level-1 suites, utest_ext's"} },
            .{ .name = "openblas-cores-l2-real", .filters = &.{ "linalg.blas._test.from_openblas.__root__.test.the level-2 suite at f32", "linalg.blas._test.from_openblas.__root__.test.the level-2 suite at f64" } },
            .{ .name = "openblas-cores-l2-complex", .filters = &.{"linalg.blas._test.from_openblas.__root__.test.the level-2 suite at complex"} },
            .{ .name = "openblas-cores-l3", .filters = &.{"linalg.blas._test.from_openblas.__root__.test.the level-3 suite"} },
            // LAPACK in five, as one is past the window
            .{ .name = "lapack-cholesky", .filters = &.{"linalg.lapack.cholesky."} },
            .{ .name = "lapack-getrf", .filters = &.{"linalg.lapack.lu.test.getrf"} },
            .{ .name = "lapack-getrs-real", .filters = &.{ "linalg.lapack.lu.test.getrs: op(A) X = B on getrf's factor, every op, NaN in one column of B, real", "linalg.lapack.lu.test.getrs: a known", "linalg.lapack.lu.test.getrs: reports" } },
            .{ .name = "lapack-getrs-complex", .filters = &.{"linalg.lapack.lu.test.getrs: op(A) X = B on getrf's factor, every op, NaN in one column of B, complex"} },
            .{ .name = "lapack-faer", .filters = &.{"linalg.lapack.faer."} },
            .{ .name = "rest", .filters = &.{ "linalg.blas.__root__.", "linalg.blas.Config.", "linalg.blas.cpu.", "simd." } },
            .{ .name = "small", .filters = &small_filters },
        },
        .query = .{ .cpu_arch = .sparc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.sparc.cpu.v9 } },
        .qemu = "qemu-system-sparc64",
        .args = &.{ "-M", "sun4u" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/sparc64.ld",
            .semihosting = false,
            .load_window = .{ .base = 0x400000, .size = 0x400000 },
        } },
    },
    .{
        // OpenBIOS on sun4m, as sparc64's on sun4u
        .name = "sparc",
        .query = .{ .cpu_arch = .sparc, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.sparc.cpu.v8 } },
        .qemu = "qemu-system-sparc",
        .args = &.{ "-M", "SS-5" },
        .memory = "256M",
        .stub = "lld cannot link 32-bit SPARC, and Zig's own ELF linker crashes on it (issues.md)",
        .boot = .{ .kernel = .{
            .linker_script = "",
            .semihosting = false,
        } },
    },
    .{
        .name = "powerpc64le",
        // LLVM PowerPC FastISel (debug) uses FPR instructions on VSX registers (issues.md)
        .optimize = .safe,
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr9 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power9" },
        .memory = "512M",
        .boot = .{
            .kernel = .{
                .linker_script = "test/qemu/ppc64.ld",
                .semihosting = false,
                // the firmware loader's 0x400000 shift, see the linker script
                .stderr_note = "kernel_addr changed",
            },
        },
    },
    .{
        .name = "powerpc64le-power10",
        // LLVM PowerPC FastISel (debug) uses FPR instructions on VSX registers (issues.md)
        .optimize = .safe,
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr10 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power10" },
        .memory = "512M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc64.ld",
            .semihosting = false,
            .stderr_note = "kernel_addr changed",
        } },
    },
    .{
        // big-endian, ELFv2 as Zig emits it; POWER10 runs every 64-bit core
        .name = "powerpc64",
        // LLVM PowerPC FastISel (debug) uses FPR instructions on VSX registers (issues.md)
        .optimize = .safe,
        .query = .{ .cpu_arch = .powerpc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr10 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power10" },
        .memory = "512M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc64.ld",
            .semihosting = false,
            .stderr_note = "kernel_addr changed",
        } },
    },
    .{
        // a G4 under OpenBIOS: the ppcg4 core
        .name = "powerpc",
        .query = .{ .cpu_arch = .powerpc, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.g4 } },
        .qemu = "qemu-system-ppc",
        .args = &.{ "-M", "mac99", "-cpu", "G4" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc32.ld",
            .semihosting = false,
        } },
    },
    .{
        // RandyOS's iBook G3 target (a 750FX, no AltiVec) on `mac99`: soft
        // float, so no kernel objects; the portable kernels run. A 750, as
        // OpenBIOS hangs on QEMU's `750fx` (issues.md); the same ISA
        .name = "powerpc-750",
        .query = .{ .cpu_arch = .powerpc, .os_tag = .freestanding, .abi = .eabi, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.@"750" } },
        .qemu = "qemu-system-ppc",
        .args = &.{ "-M", "mac99", "-cpu", "750" },
        .memory = "256M",
        .blas_cores = .none,
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc32.ld",
            .semihosting = false,
        } },
    },
    .{
        // the 750 with hard float: every 32-bit core linked, none runnable
        .name = "powerpc-750-hf",
        .query = .{ .cpu_arch = .powerpc, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.@"750" } },
        .qemu = "qemu-system-ppc",
        .args = &.{ "-M", "mac99", "-cpu", "750" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc32.ld",
            .semihosting = false,
        } },
    },
    .{
        // the G4 build in a POWER10's 32-bit mode under VOF: the server cores
        // through power6, built 32-bit
        .name = "powerpc-pseries",
        .query = .{ .cpu_arch = .powerpc, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.g4 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power10" },
        .memory = "512M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc64.ld",
            .semihosting = false,
            .stderr_note = "kernel_addr changed",
        } },
    },
    .{
        // a 32-bit VSX build on a POWER9 under VOF: the power8 and power9
        // cores' `@Vector` kernels as 32-bit code. The build drops `64bit`
        // and `fpcvt`, with which LLVM fails on compiler-rt in 32-bit mode;
        // ISA 3.0 fails it too, so power10's kernels cannot link (issues.md).
        .name = "powerpc-power9",
        .query = .{
            .cpu_arch = .powerpc,
            .os_tag = .freestanding,
            .abi = .none,
            .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr8 },
            .cpu_features_sub = std.Target.powerpc.featureSet(&.{ .@"64bit", .fpcvt }),
        },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power9" },
        .memory = "512M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc64.ld",
            .semihosting = false,
            .stderr_note = "kernel_addr changed",
        } },
    },
    .{
        // Mach-O on PowerPC, the one big-endian 32-bit Mach-O target qemu can
        // boot (OpenBIOS on `mac99`); covers the object format the other rigs
        // do not.
        .name = "powerpc-macho",
        .query = .{ .cpu_arch = .powerpc, .os_tag = .freestanding, .abi = .none, .ofmt = .macho, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.g4 } },
        .qemu = "qemu-system-ppc",
        .args = &.{ "-M", "mac99", "-cpu", "G4" },
        .memory = "256M",
        .stub = "needs the powerpc-macos Zig port, a linker script, and a check that OpenBIOS boots the image",
        .boot = .{ .kernel = .{
            .linker_script = "",
            .semihosting = false,
        } },
    },
    .{
        .name = "s390x",
        .query = .{ .cpu_arch = .s390x, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.s390x.cpu.z14 } },
        .qemu = "qemu-system-s390x",
        .args = &.{ "-M", "s390-ccw-virtio", "-cpu", "max" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/s390x.ld",
            .semihosting = false,
        } },
    },
    .{
        // RandyOS's loader's platforms: the tests as a UEFI application under
        // QEMU's OVMF and AAVMF, on rstd's UEFI `std.Io`
        .name = "x86_64-uefi",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .uefi },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "q35", "-cpu", "max" },
        .memory = "1G",
        .boot = .{ .kernel = .{
            .linker_script = "",
            .llvm = true,
            .semihosting = false,
            .firmware = "edk2-x86_64-code.fd",
        } },
    },
    .{
        .name = "aarch64-uefi",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .uefi },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a72" },
        .hvf = &.{ "-accel", "hvf", "-M", "virt", "-cpu", "host" },
        .memory = "1G",
        .boot = .{ .kernel = .{
            .linker_script = "",
            .semihosting = false,
            .firmware = "edk2-aarch64-code.fd",
        } },
    },
    .{
        .name = "linux-aarch64",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .linux },
        .qemu = "qemu-system-aarch64",
        // four CPUs, for the tests that race threads
        .args = &.{ "-M", "virt", "-cpu", "max", "-smp", "4" },
        .hvf = &.{ "-accel", "hvf", "-M", "virt", "-cpu", "host", "-smp", "4" },
        .hvf_core = "vortexm4",
        .memory = "1G",
        .boot = .{ .linux = .{ .dep = "debian_netboot_arm64", .kernel = "debian-installer/arm64/linux", .console = "ttyAMA0" } },
    },
    .{
        // the SVE and SME cores, linked only for an SVE target; TCG, as `aarch64-sve`
        .name = "linux-aarch64-sve",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .linux, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.neoverse_v1 } },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "max,pauth-impdef=on", "-smp", "4" },
        .core = "armv9sme",
        .memory = "1G",
        .boot = .{ .linux = .{ .dep = "debian_netboot_arm64", .kernel = "debian-installer/arm64/linux", .console = "ttyAMA0" } },
    },
    .{
        .name = "linux-x86_64",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .linux },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "pc", "-cpu", "max", "-smp", "4" },
        // the initramfs and its unpacked copy: the test, with every core's
        // kernels and their debug info, is about 700 MB, and `/` (where tests
        // write files) gets half the memory
        .memory = "3G",
        .boot = .{ .linux = .{ .dep = "debian_netboot_amd64", .kernel = "debian-installer/amd64/linux", .console = "ttyS0" } },
    },
    .{
        .name = "linux-arm",
        .query = .{ .cpu_arch = .arm, .os_tag = .linux },
        .qemu = "qemu-system-arm",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a15", "-smp", "4" },
        .memory = "1G",
        .boot = .{ .linux = .{ .dep = "debian_netboot_armhf", .kernel = "debian-installer/armhf/vmlinuz", .console = "ttyAMA0" } },
    },
    .{
        .name = "linux-riscv64",
        .query = .{ .cpu_arch = .riscv64, .os_tag = .linux },
        .qemu = "qemu-system-riscv64",
        // the vector extension at the widest linked core's VLEN, as `riscv64`
        .args = &.{ "-M", "virt", "-cpu", "rv64,v=true,vlen=512", "-smp", "4" },
        .memory = "1G",
        .boot = .{ .linux = .{ .dep = "debian_netboot_riscv64", .kernel = "debian-installer/riscv64/linux", .console = "ttyS0" } },
    },
    .{
        .name = "linux-powerpc64le",
        // pwr8, `f128` by library calls: clear of the FastISel bug the other
        // POWER rigs build `.safe` for (issues.md)
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .linux },
        .qemu = "qemu-system-ppc64",
        // detection sees a POWER10 (the baseline build links up to POWER9's
        // kernels); the caps TCG lacks are off, so qemu warns of none
        .args = &.{ "-M", "pseries,cap-cfpc=broken,cap-sbbc=broken,cap-ibs=broken,cap-ccf-assist=off", "-cpu", "power10", "-smp", "4" },
        .memory = "1G",
        .boot = .{ .linux = .{ .dep = "debian_netboot_ppc64el", .kernel = "debian-installer/ppc64el/vmlinux", .console = "hvc0" } },
    },
    .{
        // a POWER10 build, so power10's kernels link and run
        .name = "linux-powerpc64le-power10",
        // native `f128` from pwr9: FastISel's bug, as the freestanding rigs (issues.md)
        .optimize = .safe,
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .linux, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr10 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,cap-cfpc=broken,cap-sbbc=broken,cap-ibs=broken,cap-ccf-assist=off", "-cpu", "power10", "-smp", "4" },
        .core = "power10",
        .memory = "1G",
        .boot = .{ .linux = .{ .dep = "debian_netboot_ppc64el", .kernel = "debian-installer/ppc64el/vmlinux", .console = "hvc0" } },
    },
    .{
        .name = "linux-i386",
        .query = .{ .cpu_arch = .x86, .os_tag = .linux },
        .qemu = "qemu-system-i386",
        .args = &.{ "-M", "pc", "-cpu", "max", "-smp", "4" },
        .memory = "1G",
        // boots and runs; filtered runs pass
        .stub = "TCG's fldt/fstpt misplace the significand bytes past a page boundary, which i386's 4-aligned f80 arrays cross (issues.md)",
        .boot = .{ .linux = .{ .dep = "debian_netboot_i386", .kernel = "debian-installer/i386/linux", .console = "ttyS0" } },
    },
    .{
        // macOS 13 and up; Zig's default minimum is 15.0
        .name = "x86_64-macos",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .macos, .os_version_min = .{ .semver = .{ .major = 13, .minor = 0, .patch = 0 } } },
        .qemu = "qemu-system-x86_64",
        // `-accel` is added by `macosQemu` (`-Dmac-mttcg`)
        .args = &.{
            "-M",      "q35",
            // OSX-KVM's OpenCore-Boot.sh model, less its KVM-only flags: macOS
            // wants an Intel CPU, from 13 with AVX2. Under WHPX `max` and `host`
            // fault the firmware (randyos's `build/qemu.zig`).
            "-cpu",    "Skylake-Client,-hle,-rtm,vendor=GenuineIntel,+invtsc,+ssse3,+sse4.2,+popcnt,+avx,+aes,+xsave,+xsaveopt",
            // each vCPU is a host thread under `thread=multi`
            "-smp",    "8,cores=4,sockets=1",
            // keyboard and mouse on xHCI: macOS has no USB 1.1 driver
            "-device", "qemu-xhci,id=xhci",
            "-device", "usb-kbd,bus=xhci.0",
            "-device", "usb-tablet,bus=xhci.0",
            "-device", "usb-ehci,id=ehci",
            "-smbios", "type=2",
            "-device", "ich9-ahci,id=sata",
            "-device", "vmware-svga",
        },
        .memory = "6G",
        .boot = .macos,
    },
};

/// Tests that need an OS (files, threads) live in files the test root also
/// reaches, so the runnable ones are picked by name; the runner supplies a heap.
/// Left out: `linalg.blas.abi.`, which resolves targets through `std.testing.io`.
pub const default_filters: []const []const u8 = &([_][]const u8{ "linalg.blas.__root__.", "linalg.blas.Config.", "linalg.blas.cpu.", "linalg.blas.blas1.", "linalg.blas.blas2.", "linalg.blas.blas3.", "linalg.blas._test.", "linalg.lapack.", "simd." } ++ small_filters);

/// The multicore rigs' tests: those that thread.
const smp_filters: []const []const u8 = &.{ "linalg.blas.threading.", "linalg.blas.blas3.threaded.", "linalg.blas._test.blas2.threaded.", "threaded level 1", "at any worker count" };

/// The smaller BLAS test files: types, `fma`, `portable`, `ext` and `util`.
const small_filters = [_][]const u8{ "linalg.blas._test.types.", "linalg.blas.fma.", "linalg.blas.portable.", "linalg.blas.ext.", "linalg.blas.util." };

/// Adds the `qemu-test` step and each machine's `qemu-test-<name>` step.
pub fn addQemuTests(b: *Build, optimize: Optimize) !void {
    const all = b.step("qemu-test", "Run rstd's hardware-independent tests under qemu-system-* (all archs)");
    const check = b.step("qemu-check", "Build every kernel rig's test images and check their inline asm, without running qemu");
    const filters = b.option([]const []const u8, "qemu-test-filter", "Only run tests whose name contains this (repeatable)") orelse &.{};
    const only_parts = b.option([]const []const u8, "qemu-test-part", "Only build and run this part of a rig that splits its tests into images, e.g. blas3-kernel-complex (repeatable)") orelse &.{};
    for (only_parts) |o| if (!hasPart(o)) std.debug.panic("-Dqemu-test-part: no rig has a part named {s}", .{o});
    // On when rstd is the project being built, so a dependent never fetches OSX-KVM.
    const macos_vm = vmOption(b, "macos-vm", "the macOS VM rigs", "OSX-KVM");
    const mac: MacOptions = if (macos_vm) .{
        .user = b.option([]const u8, "mac-user", "macOS account in the VM (default randyos)") orelse "randyos",
        .port = b.option(u16, "mac-port", "Host port forwarded to the VM's ssh") orelse 2222,
        .release = b.option(MacRelease, "mac-release", "macOS release qemu-setup fetches (default sonoma)") orelse .sonoma,
        .gdb = b.option(bool, "mac-gdb", "Let gdb attach to the macOS VM on localhost:1234") orelse false,
        .detach = b.option(bool, "mac-detach", "qemu-install/start return once the macOS VM is up rather than when it stops") orelse false,
        .headless = b.option(bool, "mac-headless", "qemu-install/start show no window, VNC on 127.0.0.1:5901 instead; test runs always start the VM this way") orelse false,
        .keep = b.option(bool, "mac-keep", "Leave the macOS VM a test run started running afterward") orelse false,
        .mttcg = b.option(bool, "mac-mttcg", "Run the macOS VM's vCPUs on host threads of their own: faster in the window, but macOS has hung and misread files under it (default: true in a window, false headless)"),
        .ovmf_vars = b.option([]const u8, "mac-ovmf-vars", "Empty OVMF variable store for qemu-setup (default: QEMU's edk2-i386-vars.fd)"),
        .snapshot = b.option([]const u8, "mac-snapshot", "Name of the restore point qemu-snapshot/restore make or go back to (default clean)") orelse "clean",
    } else .{};
    // a one-time fetch of about 55 MB per kernel
    const linux_vm = vmOption(b, "linux-vm", "the Linux rigs", "their kernels");
    // the rigs run on the Mac when rstd is the project being built; a
    // dependent runs them here unless it names a host
    const host_off = b.option(bool, "test-host-off", "Run the rigs' QEMU here rather than on -Dtest-host") orelse false;
    const host_named = b.option([]const u8, "test-host", "SSH alias of a machine to run the rigs' QEMU on, Nix's from tools/nix (default macmini-m6 when rstd is the root project); images still build here, and x86_64-whpx and the macOS VM rigs stay local");
    const host_alias: ?[]const u8 = if (host_off) null else host_named orelse if (b.isRoot()) "macmini-m6" else null;
    const host_jobs = b.option(u32, "test-host-jobs", "Rigs running at once on -Dtest-host, across builds (default 8)") orelse 8;
    const host_tcg = b.option(bool, "test-host-tcg", "Run every rig on -Dtest-host under TCG, also those HVF could run") orelse false;
    const test_host: ?TestHost = if (host_alias) |alias| .{ .alias = alias, .jobs = host_jobs, .tcg = host_tcg } else null;

    const host = b.graph.host.result;
    const whpx_host = host.os.tag == .windows and host.cpu.arch == .x86_64;
    for (machines) |m| {
        if (m.stub != null or (m.whpx and !whpx_host)) continue;
        const target = b.resolveTargetQuery(m.query);
        // one rig's images and asm check, for changes that touch only its arch
        const rig_check = switch (m.boot) {
            .kernel => blk: {
                const s = b.step(b.fmt("qemu-check-{s}", .{m.name}), b.fmt("Build {s}'s test images and check their inline asm, without running qemu", .{m.name}));
                check.dependOn(s);
                break :blk s;
            },
            // not in `qemu-check`, which builds the kernel rigs
            .linux => b.step(b.fmt("qemu-check-{s}", .{m.name}), b.fmt("Build {s}'s initramfs and check its inline asm, without running qemu", .{m.name})),
            .macos => check,
        };
        if (m.boot == .kernel and m.parts.len != 0 and filters.len == 0) {
            const step = b.step(b.fmt("qemu-test-{s}", .{m.name}), b.fmt("Run rstd's tests under qemu-system on {s}, in {d} images (live logs: .zig-cache/qemu-test-{s}-<part>.log)", .{ m.name, m.parts.len, m.name }));
            for (m.parts) |p| {
                const picked = for (only_parts) |o| {
                    if (std.mem.eql(u8, o, p.name)) break true;
                } else only_parts.len == 0;
                if (!picked) continue;
                var part = m;
                part.name = b.fmt("{s}-{s}", .{ m.name, p.name });
                const run = try addKernelRun(b, part, m.boot.kernel, target, m.optimize orelse optimize, p.filters, test_host, rig_check);
                step.dependOn(&run.step);
                all.dependOn(&run.step);
            }
            continue;
        }
        const run = switch (m.boot) {
            .kernel => |k| try addKernelRun(b, m, k, target, m.optimize orelse optimize, filters, test_host, rig_check),
            .macos => if (macos_vm) try addMacosRun(b, m, target, m.optimize orelse optimize, mac) orelse continue else continue,
            .linux => |l| if (linux_vm) try addLinuxRun(b, m, l, target, m.optimize orelse optimize, filters, test_host, rig_check) orelse continue else continue,
        };
        const about = switch (m.boot) {
            .kernel, .linux => b.fmt("Run rstd's tests under qemu-system on {s} (live log: {s})", .{ m.name, liveLog(b, m) }),
            .macos => b.fmt("Run rstd's tests under qemu-system on {s}", .{m.name}),
        };
        const step = b.step(b.fmt("qemu-test-{s}", .{m.name}), about);
        step.dependOn(&run.step);
        all.dependOn(&run.step);
    }
    if (linux_vm) try addLinuxBenches(b, test_host, whpx_host);
}

fn hasPart(name: []const u8) bool {
    for (machines) |m| for (m.parts) |p| if (std.mem.eql(u8, p.name, name)) return true;
    return false;
}

/// Whether to add a family of VM rigs: on when rstd is the root project, so a
/// dependent never fetches what they need; `-D<name>-off` turns them off, and
/// `-D<name>` on in a dependent. `-off` wins.
fn vmOption(b: *Build, comptime name: []const u8, comptime rigs: []const u8, comptime needs: []const u8) bool {
    const on = b.option(bool, name, "Add " ++ rigs ++ ", fetching " ++ needs ++ " (default: when rstd is the root project)");
    const off = b.option(bool, name ++ "-off", "Leave out " ++ rigs ++ " and do not fetch " ++ needs) orelse false;
    return !off and (on orelse b.isRoot());
}

/// `-Dtest-host`: the machine the rigs' QEMU runs on (`tools/test_host.zig`).
const TestHost = struct { alias: []const u8, jobs: u32, tcg: bool };

/// Where `m` runs: the test host, or here (null).
fn testHostFor(m: Machine, host: ?TestHost) ?TestHost {
    return if (m.whpx) null else host;
}

/// `m`'s qemu and `args`, under `tools/qemu_watch.zig`: a firmware halt or
/// an hour with no output ends it and fails the step, where it would
/// otherwise never end. The slowest test prints nothing for about 20 minutes
/// on the x86 rigs. On a test host, through `tools/test_host.zig`, which
/// copies over the files `fileArg` adds.
fn qemuCommand(b: *Build, m: Machine, host: ?TestHost) *Build.Step.Run {
    const run = watchCommand(b);
    // OpenBIOS's halt, and OVMF's and AAVMF's CPU exceptions
    run.addArgs(&.{ "Stopping execution", "!!!! X64 Exception", "Synchronous Exception at", "--" });
    const h = testHostFor(m, host) orelse {
        run.addArg(m.qemu);
        run.addArgs(m.args);
        return run;
    };
    run.addArtifactArg2(testHostTool(b), .{});
    run.addArgs(&.{ h.alias, liveLog(b, m), b.fmt("{d}", .{h.jobs}) });
    run.addDirectoryArg2(b.path("tools/nix"), .{});
    run.addArgs(&.{ m.qemu, "--tcg" });
    run.addArgs(m.args);
    if (m.hvf) |hvf| if (!h.tcg) {
        run.addArg("--hvf");
        run.addArgs(hvf);
    };
    run.addArg("--");
    return run;
}

/// A file `m`'s qemu reads, copied over first on a test host.
fn fileArg(run: *Build.Step.Run, m: Machine, host: ?TestHost, file: Build.LazyPath) void {
    if (testHostFor(m, host) != null) run.addPrefixedFileArg("local:", file) else run.addFileArg(file);
}

var test_host_tool: ?*Build.Step.Compile = null;

fn testHostTool(b: *Build) *Build.Step.Compile {
    if (test_host_tool) |t| if (t.step.owner == b) return t;
    test_host_tool = b.addExecutable(.{ .name = "test_host", .root_module = b.createModule(.{
        .root_source_file = b.path("tools/test_host.zig"),
        .target = b.graph.host,
    }) });
    return test_host_tool.?;
}

/// `tools/qemu_watch.zig` with its idle limit (an hour) and `-Dfail-fast`;
/// the caller adds halt texts, `--` and the command.
pub fn watchCommand(b: *Build) *Build.Step.Run {
    const run = b.addRunArtifact(watchTool(b));
    run.addArg("60");
    if (failFast(b)) run.addArg("--fail-fast");
    return run;
}

var fail_fast: ?struct { owner: *Build, on: bool } = null;

/// `-Dfail-fast`: end a rig at its first failing test rather than run the rest.
fn failFast(b: *Build) bool {
    if (fail_fast) |f| if (f.owner == b) return f.on;
    const on = b.option(bool, "fail-fast", "End a qemu or SSH rig at its first failing test") orelse false;
    fail_fast = .{ .owner = b, .on = on };
    return on;
}

var watch_tool: ?*Build.Step.Compile = null;

pub fn watchTool(b: *Build) *Build.Step.Compile {
    if (watch_tool) |t| if (t.step.owner == b) return t;
    watch_tool = b.addExecutable(.{ .name = "qemu_watch", .root_module = b.createModule(.{
        .root_source_file = b.path("tools/qemu_watch.zig"),
        .target = b.graph.host,
    }) });
    return watch_tool.?;
}

/// Where `m`'s console is logged live, relative to the build root.
fn liveLog(b: *Build, m: Machine) []const u8 {
    return b.fmt(".zig-cache/qemu-test-{s}.log", .{m.name});
}

/// qemu's console on stdout and, as it is written, in `liveLog`: a run step
/// checking stdout shows it only once qemu exits. Naming the serial port and
/// monitor keeps `-nographic`'s firmware console. On a test host
/// `tools/test_host.zig` writes the log.
fn addConsole(b: *Build, run: *Build.Step.Run, m: Machine, host: ?TestHost) void {
    run.addArgs(&.{ "-m", m.memory, "-nographic", "-no-reboot" });
    const log = if (testHostFor(m, host) != null) "" else b.fmt(",logfile={s}", .{liveLog(b, m)});
    // not `-mon`, deprecated (a warning) from QEMU 11.1
    run.addArgs(&.{ "-chardev", b.fmt("stdio,id=con,mux=on,signal=off{s}", .{log}), "-serial", "chardev:con", "-monitor", "chardev:con" });
}

/// QEMU booting a Debian kernel with an initramfs of `test/linux/init.zig`
/// and rstd's hosted tests, which see real auxv, sysfs, procfs and threads;
/// null until the kernel's package is fetched.
fn addLinuxRun(b: *Build, m: Machine, l: Boot.Linux, target: rstdbuild.ResolvedTarget, optimize: Optimize, filters: []const []const u8, host: ?TestHost, rig_check: *Build.Step) !?*Build.Step.Run {
    const kernel = b.dependencyLazy(l.dep, .{}) catch return null;
    // ThreadSanitizer wants glibc, and its loader, from Nix (`tools/nix/blas.sh`)
    const nix: ?[]const u8 = if (!buildroot.tsan or host == null) null else switch (target.result.cpu.arch) {
        .x86_64 => "x86_64-linux",
        .aarch64 => "aarch64-linux",
        else => null,
    };
    var query = target.query;
    if (nix != null) query.abi = .gnu;
    const t = if (nix != null) b.resolveTargetQuery(query) else target;
    const hvf = if (testHostFor(m, host)) |h| m.hvf != null and !h.tcg else false;
    const core = if (hvf) m.hvf_core else m.core;
    // a filtered run still prints the core
    const core_test = "linalg.blas.__root__.test.the core in use";
    const run_filters = if (core != null and filters.len != 0) try std.mem.concat(b.allocator, []const u8, &.{ filters, &.{core_test} }) else filters;
    const mod = try buildroot.addRstdModule(b, null, t, optimize, .all, false, false);
    const test_exe = try buildroot.addRstdTest(b, mod, t, optimize, run_filters);
    const init_exe = b.addExecutable(.{ .name = b.fmt("init-{s}", .{m.name}), .root_module = b.createModule(.{
        .root_source_file = b.path("test/linux/init.zig"),
        .target = target,
        .optimize = optimize,
    }) });
    const pack = b.addRunArtifact(initramfsTool(b));
    const cpio = pack.addOutputFileArg2("initramfs.cpio", .{});
    pack.addArtifactArg2(init_exe, .{ .prefix = "init=" });
    pack.addArtifactArg2(test_exe, .{ .prefix = "test=" });
    if (nix) |system| pack.addPrefixedFileArg("+", (try sshtest.nixBlas(b, host.?.alias, system, false)).captureStdOut(.{}));
    rig_check.dependOn(&pack.step);

    const run = qemuCommand(b, m, host);
    addConsole(b, run, m, host);
    run.addArg("-kernel");
    fileArg(run, m, host, kernel.path(l.kernel));
    run.addArg("-initrd");
    fileArg(run, m, host, cpio);
    // panic=-1: a dead init restarts at once, which -no-reboot makes an exit
    run.addArgs(&.{ "-append", b.fmt("console={s} rdinit=/init panic=-1 quiet loglevel=3", .{l.console}) });
    run.has_side_effects = true;
    // `passed` in test/linux/init.zig
    run.expectStdOutMatch("linuxtest: /test exited 0");
    // no newline: the serial console sends `\r\n`
    if (core) |c| run.expectStdOutMatch(b.fmt("blas core: {s}", .{c}));
    return run;
}

var load_window_tool: ?*Build.Step.Compile = null;

fn loadWindowTool(b: *Build) *Build.Step.Compile {
    if (load_window_tool) |t| if (t.step.owner == b) return t;
    load_window_tool = b.addExecutable(.{ .name = "load_window", .root_module = b.createModule(.{
        .root_source_file = b.path("tools/load_window.zig"),
        .target = b.graph.host,
    }) });
    return load_window_tool.?;
}

var initramfs_tool: ?*Build.Step.Compile = null;

fn initramfsTool(b: *Build) *Build.Step.Compile {
    if (initramfs_tool) |t| if (t.step.owner == b) return t;
    initramfs_tool = b.addExecutable(.{ .name = "initramfs", .root_module = b.createModule(.{
        .root_source_file = b.path("tools/initramfs.zig"),
        .target = b.graph.host,
    }) });
    return initramfs_tool.?;
}

/// The Linux VMs of `qemu-bench-linux-<arch>`, on the host's own CPU: never
/// emulated, so the HVF arguments stand in for TCG's too, and a host
/// without HVF fails rather than emulates.
const bench_machines = [_]struct { m: Machine, cpu: ?*const std.Target.Cpu.Model, nix: []const u8, blis: bool, aocl: bool }{
    .{
        .m = .{
            .name = "bench-linux-x86_64",
            .query = .{ .cpu_arch = .x86_64, .os_tag = .linux, .abi = .gnu },
            .qemu = "qemu-system-x86_64",
            // `host` stops in SeaBIOS (x86_64-whpx); this host's own family
            // keeps the vendor and family the libraries dispatch on
            .args = &.{ "-accel", "whpx", "-M", "pc", "-cpu", "EPYC-Milan", "-smp", "8" },
            .whpx = true,
            .memory = "4G",
            .boot = .{ .linux = .{ .dep = "debian_netboot_amd64", .kernel = "debian-installer/amd64/linux", .console = "ttyS0" } },
        },
        // built here, for this CPU, as bench-native is
        .cpu = null,
        .nix = "x86_64-linux",
        .blis = true,
        // AOCL-libFLAME: nixpkgs builds it for x86_64-linux only
        .aocl = true,
    },
    .{
        .m = .{
            .name = "bench-linux-aarch64",
            .query = .{ .cpu_arch = .aarch64, .os_tag = .linux, .abi = .gnu },
            .qemu = "qemu-system-aarch64",
            .args = &hvf_bench_args,
            .hvf = &hvf_bench_args,
            .memory = "4G",
            .boot = .{ .linux = .{ .dep = "debian_netboot_arm64", .kernel = "debian-installer/arm64/linux", .console = "ttyAMA0" } },
        },
        // as ssh-bench-aarch64-macos builds for the Mac
        .cpu = &std.Target.aarch64.cpu.apple_m1,
        .nix = "aarch64-linux",
        .blis = false,
        .aocl = false,
    },
};

const hvf_bench_args = [_][]const u8{ "-accel", "hvf", "-M", "virt", "-cpu", "host", "-smp", "6" };

/// `qemu-bench-linux-<arch>`: `tools/blasbench/blas_bench.zig`'s library
/// sections in a Linux VM, against nixpkgs' OpenBLAS and BLIS, which
/// `tools/nix/blas.sh` packs on the Nix machine (`-Dtest-host`). The VM
/// boots like the Linux rigs, the bench as `/test`, with `-Dvm-bench-arg`.
fn addLinuxBenches(b: *Build, host: ?TestHost, whpx_host: bool) !void {
    // on the kernel command line, which the build fixes before `--` is known
    const bench_args = b.option([]const []const u8, "vm-bench-arg", "An argument for the bench in the qemu-bench-linux steps (repeatable), as `--` gives the other bench steps") orelse &.{};
    for (bench_machines) |bm| {
        const m = bm.m;
        if (m.whpx and !whpx_host) continue;
        const l = m.boot.linux;
        const step = b.step(b.fmt("qemu-{s}", .{m.name}), b.fmt("Race rstd's BLAS against nixpkgs' libraries in a Linux VM on the host's CPU (live log: {s})", .{liveLog(b, m)}));
        const kernel = b.dependencyLazy(l.dep, .{}) catch continue;
        var query = m.query;
        query.cpu_model = if (bm.cpu) |c| .{ .explicit = c } else .native;
        const exe = try sshtest.benchExe(b, b.resolveTargetQuery(query), null);
        const init_exe = b.addExecutable(.{ .name = b.fmt("init-{s}", .{m.name}), .root_module = b.createModule(.{
            .root_source_file = b.path("test/linux/init.zig"),
            .target = b.resolveTargetQuery(.{ .cpu_arch = m.query.cpu_arch, .os_tag = .linux }),
            .optimize = .fast,
        }) });
        const pack = b.addRunArtifact(initramfsTool(b));
        const cpio = pack.addOutputFileArg2("initramfs.cpio", .{});
        pack.addArtifactArg2(init_exe, .{ .prefix = "init=" });
        pack.addArtifactArg2(exe, .{ .prefix = "test=" });
        // after `--`, init's arguments, which it passes on
        var append: std.ArrayList(u8) = .empty;
        try append.print(b.allocator, "console={s} rdinit=/init panic=-1 quiet loglevel=3 -- --openblas=/flyoff/libopenblas.so", .{l.console});
        if (bm.blis) try append.appendSlice(b.allocator, " --blis=/flyoff/libblis.so");
        // the reference binaries (`tools/nix/blas.sh`); one it skipped is skipped here
        try append.appendSlice(b.allocator, " --reference=/flyoff/liblapack.so --eigen=/flyoff/libeigen.so");
        if (bm.aocl) try append.appendSlice(b.allocator, " --aocl-libflame=/flyoff/libflame.so");
        const vm_append = try append.clone(b.allocator);

        // the Mac's flyoff packs the libraries there: the initramfs without them
        const host_mac = m.hvf != null;
        if (host_mac) {
            const base = b.addRunArtifact(initramfsTool(b));
            const base_cpio = base.addOutputFileArg2("base.cpio", .{});
            base.addArtifactArg2(init_exe, .{ .prefix = "init=" });
            base.addArtifactArg2(exe, .{ .prefix = "test=" });
            try addFlyoffStage(b, m, kernel.path(l.kernel), base_cpio, "base.cpio", vm_append.items);
        }

        const nix = host orelse {
            step.dependOn(&b.addFail(b.fmt("qemu-{s} gets its libraries from Nix on -Dtest-host", .{m.name})).step);
            if (!host_mac) b.step("flyoff-stage-zen3", "Stage the Zen 3 flyoff's files in zig-out/flyoff-zen3 (tools/blasbench/flyoff.sh)").dependOn(&b.addFail("flyoff-stage-zen3 gets the VM's libraries from Nix on -Dtest-host").step);
            continue;
        };
        // the libraries' closure, fetched once per flake.lock
        const libs = (try sshtest.nixBlas(b, nix.alias, bm.nix, true)).captureStdOut(.{});
        pack.addPrefixedFileArg("+", libs);
        if (!host_mac) try addFlyoffStage(b, m, kernel.path(l.kernel), cpio, "initramfs.cpio", vm_append.items);

        const run = qemuCommand(b, m, host);
        addConsole(b, run, m, host);
        run.addArg("-kernel");
        fileArg(run, m, host, kernel.path(l.kernel));
        run.addArg("-initrd");
        fileArg(run, m, host, cpio);
        for (bench_args) |a| try append.print(b.allocator, " {s}", .{a});
        run.addArgs(&.{ "-append", append.items });
        run.has_side_effects = true;
        run.stdio = .inherit;
        step.dependOn(&run.step);
    }
}

/// `flyoff-stage-<host>`: `tools/blasbench/flyoff.sh`'s files for the host
/// whose CPU `m` runs on, in `zig-out/flyoff-<host>`: the host's own bench
/// (`blas_bench`), the VM's kernel (`linux`) and initramfs (`cpio`, as
/// `cpio_name`), its QEMU arguments (`vm.args`, one a line) and kernel
/// command line up to the bench's own arguments (`vm.append`). The Mac's
/// also gets `tools/nix` (`nix`), `tools/blasbench/nix` (`bench-nix`) and
/// `flyoff.sh`, since it runs there.
fn addFlyoffStage(b: *Build, m: Machine, kernel: Build.LazyPath, cpio: Build.LazyPath, cpio_name: []const u8, vm_append: []const u8) !void {
    const mac = m.hvf != null;
    const name = if (mac) "m6" else "zen3";
    const dir: Build.InstallDir = .{ .custom = b.fmt("flyoff-{s}", .{name}) };
    const step = b.step(b.fmt("flyoff-stage-{s}", .{name}), b.fmt("Stage the {s} flyoff's files in zig-out/flyoff-{s} (tools/blasbench/flyoff.sh)", .{ if (mac) "M6" else "Zen 3", name }));
    const target = if (mac) b.resolveTargetQuery(.{ .cpu_arch = .aarch64, .os_tag = .macos }) else b.graph.host;
    const exe = try sshtest.benchExe(b, target, null);
    step.dependOn(&b.addInstallArtifact(exe, .{ .dest_dir = .{ .override = dir } }).step);
    step.dependOn(&b.addInstallFileWithDir(kernel, dir, "linux").step);
    step.dependOn(&b.addInstallFileWithDir(cpio, dir, cpio_name).step);
    var args: std.ArrayList(u8) = .empty;
    const console_args = [_][]const u8{ "-m", m.memory, "-nographic", "-no-reboot", "-chardev", "stdio,id=con,mux=on,signal=off", "-serial", "chardev:con", "-monitor", "chardev:con" };
    for ([_][]const []const u8{ &.{m.qemu}, m.hvf orelse m.args, &console_args }) |part| {
        for (part) |a| try args.print(b.allocator, "{s}\n", .{a});
    }
    const files = b.addWriteFiles();
    step.dependOn(&b.addInstallFileWithDir(files.add("vm.args", args.items), dir, "vm.args").step);
    step.dependOn(&b.addInstallFileWithDir(files.add("vm.append", vm_append), dir, "vm.append").step);
    step.dependOn(&b.addInstallFileWithDir(b.path("tools/blasbench/flyoff.sh"), dir, "flyoff.sh").step);
    if (mac) {
        step.dependOn(&b.addInstallDirectory(.{ .source_dir = b.path("tools/nix"), .install_dir = dir, .install_subdir = "nix" }).step);
        step.dependOn(&b.addInstallDirectory(.{ .source_dir = b.path("tools/blasbench/nix"), .install_dir = dir, .install_subdir = "bench-nix" }).step);
    }
}

/// QEMU booting `m`'s freestanding test kernel, once its inline asm is checked
/// (`blaskern.asmCheck`), which `check` also does.
fn addKernelRun(b: *Build, m: Machine, k: Boot.Kernel, target: rstdbuild.ResolvedTarget, optimize: Optimize, filters: []const []const u8, host: ?TestHost, check: *Build.Step) !*Build.Step.Run {
    const test_mod = b.createModule(.{
        .root_source_file = b.path("src/_test_qemu.zig"),
        .target = target,
        .optimize = optimize,
        .code_model = k.code_model,
        // threads only where the kernel starts more cores
        .single_threaded = m.cpus == 1,
    });
    if (k.firmware != null) try addUefiImports(b, test_mod);
    // the kernels of the architectures that have any
    try blaskern.addBlasKernels(b, test_mod, target, m.blas_cores, if (m.quad_precision) true else null);
    const test_exe = b.addTest(.{
        .name = b.fmt("rstd-qemu-{s}", .{m.name}),
        .root_module = test_mod,
        .test_runner = .{ .path = b.path("test/qemu/runner.zig"), .mode = .simple },
        .filters = if (filters.len != 0) filters else if (m.cpus > 1) smp_filters else default_filters,
        .use_llvm = if (k.llvm) true else null,
        // every rig needs its linker script, which only lld honors
        .use_lld = true,
    });
    if (k.firmware == null) test_exe.setLinkerScript(b.path(k.linker_script));
    check.dependOn(&test_exe.step);

    const run = qemuCommand(b, m, host);
    if (k.load_window) |w| {
        const fit = b.addRunArtifact(loadWindowTool(b));
        fit.setName(b.fmt("check load window of {s}", .{test_exe.name}));
        fit.addArgs(&.{ m.name, b.fmt("0x{x}", .{w.base}), b.fmt("0x{x}", .{w.size}) });
        fit.addFileArg(test_exe.getEmittedBin());
        _ = fit.addOutputFileArg2("span", .{});
        check.dependOn(&fit.step);
        run.step.dependOn(&fit.step);
    }
    // only x86_64 defaults to the self-hosted backend; `zig cc` ignores `ofmt`
    const llvm = k.llvm or target.result.cpu.arch != .x86_64;
    if (llvm and m.query.ofmt == null) {
        const asm_check = try blaskern.asmCheck(b, test_exe, m.query);
        check.dependOn(asm_check);
        run.step.dependOn(asm_check);
    }
    addConsole(b, run, m, host);
    if (m.cpus > 1) run.addArgs(&.{ "-smp", b.fmt("{d}", .{m.cpus}) });
    if (k.semihosting) run.addArgs(&.{ "-semihosting-config", "enable=on,target=native,chardev=con" });
    if (k.firmware) |fw| {
        // the app at the removable-media boot path of a read-only FAT drive
        // (vvfat), booted first; no NIC, whose network boot would come first
        const esp = b.addWriteFiles();
        const boot = switch (target.result.cpu.arch) {
            .x86_64 => "BOOTX64.EFI",
            .aarch64 => "BOOTAA64.EFI",
            else => unreachable,
        };
        _ = esp.addCopyFile(test_exe.getEmittedBin(), b.fmt("EFI/BOOT/{s}", .{boot}));
        // read-only pflash, without a variable store: OVMF's code is no
        // multiple of 64 KiB, which `-bios` wants
        if (testHostFor(m, host) == null) run.step.dependOn(&b.addFail(b.fmt("{s} runs on -Dtest-host, whose tools/nix/run.sh finds QEMU's firmware", .{m.name})).step);
        run.addArgs(&.{ "-drive", b.fmt("if=pflash,unit=0,format=raw,readonly=on,file=@qemu-data@/{s}", .{fw}), "-nic", "none", "-drive" });
        const drive = "if=none,id=esp,format=raw,readonly=on,file=fat:";
        run.addPrefixedDirectoryArg(if (testHostFor(m, host) != null) drive ++ "local:" else drive, esp.getDirectory());
        run.addArgs(&.{ "-device", "virtio-blk-pci,drive=esp,bootindex=0" });
    } else {
        run.addArg("-kernel");
        fileArg(run, m, host, test_exe.getEmittedBin());
    }
    run.has_side_effects = true;
    // HVF has no semihosting, so no exit status
    const hvf = if (testHostFor(m, host)) |h| m.hvf != null and !h.tcg else false;
    if (!k.semihosting or hvf) run.expectStdOutMatch(", 0 failed");
    if (k.stderr_note) |note| run.expectStdErrMatch(note);
    return run;
}

/// What the UEFI boot code (`test/qemu/uefi.zig`) needs of the test module:
/// the module itself as `qemutest_root`, for rstd's UEFI `std.Io`, and the
/// `build_options` that Io's logger reads, with no scope levels.
pub fn addUefiImports(b: *Build, mod: *Build.Module) !void {
    mod.addImport("qemutest_root", mod);
    const options = b.addOptions();
    options.addOption([]const []const u8, "logger_scope_names", &.{});
    options.addOption([]const @typeInfo(std.log.Level).@"enum".tag_type, "logger_scope_name_levels", &.{});
    mod.addImport("build_options", try rstdbuild.optionsModule(options));
}

/// `m`'s `qemu-setup-`, `qemu-install-` and `qemu-start-<name>` steps, and the
/// test run: rstd's unit tests copied into the running VM and run over SSH.
/// Null until the lazy OSX-KVM dependency is fetched. The VM's files are named
/// relative to the build root, where its commands run, so none is hashed.
fn addMacosRun(b: *Build, m: Machine, target: rstdbuild.ResolvedTarget, optimize: Optimize, mac: MacOptions) !?*Build.Step.Run {
    const kvm = b.dependencyLazy("osx_kvm", .{}) catch return null;
    const vm = b.fmt(".macos-vm/{s}", .{m.name});
    try b.root.createDirPath(b.graph.io, vm);
    // setup skips what is there, so a change there configures again
    b.dependOnDirectoryContents(b.path(vm));
    const paths: MacPaths = .{
        .dmg = b.fmt("{s}/BaseSystem.dmg", .{vm}),
        .img = b.fmt("{s}/BaseSystem.img", .{vm}),
        .disk = b.fmt("{s}/mac_hdd.qcow2", .{vm}),
        .vars = b.fmt("{s}/OVMF_VARS.fd", .{vm}),
        .opencore = b.fmt("{s}/OpenCore.qcow2", .{vm}),
    };

    const setup = b.step(b.fmt("qemu-setup-{s}", .{m.name}), b.fmt("macOS VM {s}: fetch a recovery image, make the disk", .{m.name}));
    var fetch: ?*Build.Step.Run = null;
    if (!try exists(b, paths.dmg)) {
        const run = console(b, &.{"python"});
        run.addFileArg(kvm.path("fetch-macOS-v2.py"));
        run.addArgs(&.{ "-s", @tagName(mac.release) });
        // it saves where it runs; the package is read-only
        run.setCwd(b.path(vm));
        setup.dependOn(&run.step);
        fetch = run;
    }
    if (!try exists(b, paths.img)) {
        const run = console(b, &.{ "qemu-img", "convert", "-f", "dmg", "-O", "raw", paths.dmg, paths.img });
        if (fetch) |f| run.step.dependOn(&f.step);
        setup.dependOn(&run.step);
    }
    if (!try exists(b, paths.disk)) setup.dependOn(&console(b, &.{ "qemu-img", "create", "-f", "qcow2", paths.disk, "96G" }).step);
    // The firmware's settings (NVRAM: the default boot disk, boot-args). It
    // must be the variable store `OVMF_CODE_4M.fd` was built for, 528 KiB;
    // OSX-KVM's 128 KiB ones are for the older 2 MiB firmware, and with one
    // the firmware keeps NVRAM in memory and the choices never persist.
    if (!try exists(b, paths.vars)) {
        if (mac.ovmf_vars orelse ovmfVarsTemplate(b)) |template| {
            const copy = b.addUpdateSourceFiles();
            copy.addCopyFileToSource(.{ .cwd_relative = template }, paths.vars);
            setup.dependOn(&copy.step);
        } else setup.dependOn(&b.addFail("no edk2-i386-vars.fd next to qemu-system-x86_64; pass the 4 MiB OVMF variable store with -Dmac-ovmf-vars=<path>").step);
    }
    // OpenCore's own copy, writable: `configure_macos.sh` edits its config
    // (boot-args) from inside the VM
    const opencore_missing = !try exists(b, paths.opencore);
    if (opencore_missing) {
        const copy = b.addUpdateSourceFiles();
        copy.addCopyFileToSource(kvm.path("OpenCore/OpenCore.qcow2"), paths.opencore);
        setup.dependOn(&copy.step);
    }
    // what keeps the VM from starting
    const not_ready: ?*Build.Step = if (if (rstdbuild.io.fileSize(b.graph.io, try b.root.joinString(b.allocator, paths.vars))) |size| size != ovmf_vars_size else false)
        &b.addFail(b.fmt("{s} is not a 528 KiB OVMF variable store, so NVRAM would not persist; with the VM off, delete it and run qemu-setup-{s}", .{ paths.vars, m.name })).step
    else if (opencore_missing)
        &b.addFail(b.fmt("no {s}; run qemu-setup-{s}", .{ paths.opencore, m.name })).step
    else
        null;
    const windows = b.graph.host.result.os.tag == .windows;

    const tool = b.addExecutable(.{ .name = "macos_vm", .root_module = b.createModule(.{
        .root_source_file = b.path("tools/macos_vm.zig"),
        .target = b.graph.host,
    }) });
    // a second qemu on the same files fails, and deletes the first's pid file
    const guard = b.addRunArtifact(tool);
    guard.setCwd(b.path("."));
    guard.stdio = .inherit;
    guard.has_side_effects = true;
    guard.disable_zig_progress = true;
    guard.addArgs(&.{ "guard", vm });

    inline for (.{ "install", "start" }) |cmd| {
        // Detached on Windows, `start` gives qemu its own console (the
        // monitor's) and returns; elsewhere `-daemonize` (`qemuArgs`).
        const run = console(b, if (mac.detach and windows) &.{ "cmd", "/c", "start", "", m.qemu } else &.{m.qemu});
        try qemuArgs(b, run, m, kvm, paths, mac, .{
            .install = std.mem.eql(u8, cmd, "install"),
            .headless = mac.headless,
            .monitor = if (mac.detach) (if (windows) .stdio else .daemon) else .stdio,
        });
        if (not_ready) |s| run.step.dependOn(s);
        run.step.dependOn(&guard.step);
        b.step(b.fmt("qemu-" ++ cmd ++ "-{s}", .{m.name}), b.fmt("macOS VM {s}: " ++ cmd ++ "{s}", .{ m.name, if (mac.headless) ", headless" else ", in a window" })).dependOn(&run.step);
    }

    // A named restore point: a snapshot of each disk (macOS's and OpenCore's,
    // whose config holds the boot-args) plus a copy of the firmware
    // settings, which live outside them. qemu-img will not touch a disk the
    // running VM holds, so neither step works on a live VM.
    const saved_vars = b.fmt("{s}/OVMF_VARS.{s}.fd", .{ vm, mac.snapshot });
    const keep = b.addUpdateSourceFiles();
    keep.addCopyFileToSource(b.path(paths.vars), saved_vars);
    const put = b.addUpdateSourceFiles();
    put.addCopyFileToSource(b.path(saved_vars), paths.vars);
    for ([_][]const u8{ paths.disk, paths.opencore }) |disk| {
        const snap = console(b, &.{ "qemu-img", "snapshot", "-c", mac.snapshot, disk });
        // qemu-img would add a second snapshot of the same name
        if (hasSnapshot(b, disk, mac.snapshot)) snap.step.dependOn(&b.addFail(b.fmt("{s} already has a snapshot named {s}; pick another with -Dmac-snapshot", .{ disk, mac.snapshot })).step);
        snap.step.dependOn(&guard.step);
        keep.step.dependOn(&snap.step);
        // restore points from before OpenCore had its own disk lack one there
        if (disk.ptr == paths.opencore.ptr and !hasSnapshot(b, disk, mac.snapshot)) continue;
        const back = console(b, &.{ "qemu-img", "snapshot", "-a", mac.snapshot, disk });
        back.step.dependOn(&guard.step);
        put.step.dependOn(&back.step);
    }
    b.step(b.fmt("qemu-snapshot-{s}", .{m.name}), b.fmt("macOS VM {s}: save a restore point (-Dmac-snapshot) of the stopped VM", .{m.name})).dependOn(&keep.step);
    b.step(b.fmt("qemu-restore-{s}", .{m.name}), b.fmt("macOS VM {s}: go back to a restore point (-Dmac-snapshot) of the stopped VM", .{m.name})).dependOn(&put.step);

    const mod = try buildroot.addRstdModule(b, null, target, optimize, .all, false, false);
    // In one program, so the VM it started is stopped when tests fail too.
    const run = b.addRunArtifact(tool);
    run.setCwd(b.path("."));
    run.stdio = .inherit;
    // QEMU, started from it, would hold the progress pipe open
    run.disable_zig_progress = true;
    run.has_side_effects = true;
    run.addArgs(&.{ "test", vm, b.fmt("{d}", .{mac.port}), b.fmt("{s}@localhost", .{mac.user}), if (mac.keep) "true" else "false" });
    const exe = try buildroot.addRstdTest(b, mod, target, optimize, &.{});
    run.addArtifactArg2(exe, .{});
    run.addArg(m.qemu);
    // what that runs, built alone, so a build queue frees its slot before the VM boots
    const check = b.step(b.fmt("qemu-check-{s}", .{m.name}), b.fmt("Build what qemu-test-{s} runs, without starting the VM", .{m.name}));
    check.dependOn(&exe.step);
    check.dependOn(&tool.step);
    try qemuArgs(b, run, m, kvm, paths, mac, .{ .headless = true, .monitor = .socket });
    if (mac.user.len == 0) run.step.dependOn(&b.addFail("the macOS VM test needs -Dmac-user=<account>").step);
    if (not_ready) |s| run.step.dependOn(s);

    const Control = struct { name: []const u8, about: []const u8, login: bool };
    for ([_]Control{
        .{ .name = "up", .about = "start headless unless SSH answers; returns once it does", .login = true },
        .{ .name = "down", .about = "shut down over SSH (passwordless sudo shutdown), else kill", .login = true },
        .{ .name = "kill", .about = "kill its qemu (`qemu.pid`), as pulling the plug", .login = false },
    }) |c| {
        const ctl = b.addRunArtifact(tool);
        ctl.setCwd(b.path("."));
        ctl.stdio = .inherit;
        ctl.disable_zig_progress = true;
        ctl.has_side_effects = true;
        ctl.addArgs(&.{ c.name, vm });
        if (c.login) {
            ctl.addArgs(&.{ b.fmt("{d}", .{mac.port}), b.fmt("{s}@localhost", .{mac.user}) });
            if (mac.user.len == 0) ctl.step.dependOn(&b.addFail(b.fmt("qemu-{s} needs -Dmac-user=<account>", .{c.name})).step);
        }
        if (std.mem.eql(u8, c.name, "up")) {
            ctl.addArg(m.qemu);
            try qemuArgs(b, ctl, m, kvm, paths, mac, .{ .headless = true, .monitor = .socket });
            if (not_ready) |s| ctl.step.dependOn(s);
        }
        b.step(b.fmt("qemu-{s}-{s}", .{ c.name, m.name }), b.fmt("macOS VM {s}: {s}", .{ m.name, c.about })).dependOn(&ctl.step);
    }
    return run;
}

/// The macOS VMs' SSH login (`-Dmac-user`, `-Dmac-port`) and release.
const MacOptions = struct { user: []const u8 = "", port: u16 = 2222, release: MacRelease = .sonoma, gdb: bool = false, detach: bool = false, headless: bool = false, keep: bool = false, mttcg: ?bool = null, ovmf_vars: ?[]const u8 = null, snapshot: []const u8 = "clean" };

/// The releases `fetch-macOS-v2.py` offers that the tests run on: the machine
/// builds them for macOS 13 and up. The tags are the fetcher's `-s` names.
/// Sonoma is the one OSX-KVM recommends.
const MacRelease = enum { ventura, sonoma, sequoia, tahoe };

/// A macOS machine's files in `.macos-vm/<name>/`, relative to the build root.
const MacPaths = struct { dmg: []const u8, img: []const u8, disk: []const u8, vars: []const u8, opencore: []const u8 };

/// How `qemuArgs` sets up the VM.
const QemuMode = struct {
    /// the recovery image too
    install: bool = false,
    /// no window; VNC on 127.0.0.1:5901
    headless: bool = false,
    /// `stdio`: a `(qemu)` prompt in this console; `socket`: none, the
    /// monitor on `monitor.sock` in the VM's directory; `daemon`: that, and
    /// qemu goes to the background
    monitor: enum { stdio, socket, daemon } = .stdio,
};

/// Appends `m`'s command line to `run`, which starts QEMU, plus OSX-KVM's
/// firmware, SMC key and OpenCore, and the disks.
fn qemuArgs(b: *Build, run: *Build.Step.Run, m: Machine, kvm: *Build.Dependency, paths: MacPaths, mac: MacOptions, mode: QemuMode) !void {
    // TCG emulates (slow) but honors the Intel vendor the CPU model sets;
    // under WHPX ("whpx", as randyos's own VM uses) the kernel stopped at
    // EXITBS:START, as if it saw the host's AMD CPUID. A host thread per
    // vCPU (`thread=multi`) keeps the window usable, but under it macOS hung
    // with interrupts off at a reboot and the installer's restart, and
    // decompressing system files failed now and then (a boot loop on a hash
    // mismatch, the login window crashing on a page-in) where the same files
    // read back fine; so one host thread headless, where speed matters less.
    // A 2 GB translation cache keeps macOS's code from being translated again.
    const multi = mac.mttcg orelse !mode.headless;
    run.addArgs(&.{ "-accel", b.fmt("tcg,thread={s},tb-size=2048", .{if (multi) "multi" else "single"}) });
    run.addArgs(m.args);
    run.addArgs(&.{ "-m", m.memory, "-device", b.fmt("isa-applesmc,osk={s}", .{try smcKey(b, kvm)}), "-drive" });
    run.addPrefixedFileArg("if=pflash,format=raw,readonly=on,file=", kvm.path("OVMF_CODE_4M.fd"));
    run.addArgs(&.{
        "-drive",  b.fmt("if=pflash,format=raw,file={s}", .{paths.vars}),
        "-drive",  b.fmt("id=OpenCoreBoot,if=none,format=qcow2,file={s}", .{paths.opencore}),
        "-device", "ide-hd,bus=sata.2,drive=OpenCoreBoot",
        "-drive",  b.fmt("id=MacHDD,if=none,format=qcow2,cache=writeback,file={s}", .{paths.disk}),
        "-device", "ide-hd,bus=sata.4,drive=MacHDD",
        "-netdev", b.fmt("user,id=net0,hostfwd=tcp::{d}-:22", .{mac.port}),
        "-device", "virtio-net-pci,netdev=net0,id=net0,mac=52:54:00:c9:18:27",
    });
    if (mode.install) run.addArgs(&.{
        "-drive",  b.fmt("id=InstallMedia,if=none,format=raw,file={s}", .{paths.img}),
        "-device", "ide-hd,bus=sata.3,drive=InstallMedia",
    });
    if (mode.headless) run.addArgs(&.{ "-display", "none", "-vnc", "127.0.0.1:1" });
    // A `(qemu)` prompt: `info status` says whether it runs, and `RIP` in
    // `info registers` moves while the guest executes. Without a console
    // the monitor is a socket.
    const dir = std.Io.Dir.path.dirname(paths.disk).?;
    // for `qemu-kill` and test runs; QEMU removes it at a clean exit
    run.addArgs(&.{ "-pidfile", b.fmt("{s}/qemu.pid", .{dir}) });
    // COM1: macOS's console with `serial=3` in OpenCore's boot-args
    run.addArgs(&.{ "-serial", b.fmt("file:{s}/serial.log", .{dir}) });
    const sock = b.fmt("{s}/monitor.sock", .{dir});
    switch (mode.monitor) {
        .stdio => run.addArgs(&.{ "-monitor", "stdio" }),
        .socket => run.addArgs(&.{ "-monitor", b.fmt("unix:{s},server,nowait", .{sock}) }),
        .daemon => run.addArgs(&.{ "-daemonize", "-monitor", b.fmt("unix:{s},server,nowait", .{sock}) }),
    }
    // gdb-multiarch, then `target remote localhost:1234`
    if (mac.gdb) run.addArg("-s");
}

/// The Apple SMC key OSX-KVM's boot script passes, read from it.
fn smcKey(b: *Build, kvm: *Build.Dependency) ![]const u8 {
    const path = try kvm.builder.root.joinString(b.allocator, "OpenCore-Boot.sh");
    const script = try std.Io.Dir.cwd().readFileAlloc(b.graph.io, path, b.allocator, .limited(1 << 20));
    const start = (std.mem.find(u8, script, "osk=\"") orelse return error.NoSmcKey) + 5;
    const end = std.mem.findScalarPos(u8, script, start, '"') orelse return error.NoSmcKey;
    return script[start..end];
}

/// Whether the qcow2 image at `path` holds a snapshot named `name`, from its
/// snapshot table (QEMU's `docs/interop/qcow2.rst`); false when unreadable.
fn hasSnapshot(b: *Build, path: []const u8, name: []const u8) bool {
    const io = b.graph.io;
    const file = b.root.openFile(io, path, .{}) catch return false;
    defer file.close(io);
    var header: [72]u8 = undefined;
    if ((file.readPositionalAll(io, &header, 0) catch return false) != header.len) return false;
    if (!std.mem.eql(u8, header[0..4], "QFI\xfb")) return false;
    const count = std.mem.readInt(u32, header[60..64], .big);
    var at = std.mem.readInt(u64, header[64..72], .big);
    for (0..count) |_| {
        // l1 offset and size, id and name lengths, dates, vm state, extra
        var entry: [40]u8 = undefined;
        if ((file.readPositionalAll(io, &entry, at) catch return false) != entry.len) return false;
        const id_len = std.mem.readInt(u16, entry[12..14], .big);
        const name_len = std.mem.readInt(u16, entry[14..16], .big);
        const extra_len = std.mem.readInt(u32, entry[36..40], .big);
        const name_at = at + entry.len + extra_len + id_len;
        if (name_len == name.len) {
            var buf: [256]u8 = undefined;
            if (name_len <= buf.len and (file.readPositionalAll(io, buf[0..name_len], name_at) catch return false) == name_len and
                std.mem.eql(u8, buf[0..name_len], name)) return true;
        }
        // entries are padded to a multiple of 8
        at = std.mem.alignForward(u64, name_at + name_len, 8);
    }
    return false;
}

/// The size of the variable store `OVMF_CODE_4M.fd` uses.
const ovmf_vars_size = 540672;

/// QEMU's own empty variable store for its 4 MiB x86 edk2, which
/// `OVMF_CODE_4M.fd` matches: `../share/qemu/edk2-i386-vars.fd` from
/// `qemu-system-x86_64`, as MSYS2 and the official builds lay it out.
fn ovmfVarsTemplate(b: *Build) ?[]const u8 {
    const qemu = b.findProgram(.{ .names = &.{"qemu-system-x86_64"} }) orelse return null;
    const bin = std.Io.Dir.path.dirname(qemu) orelse return null;
    for ([_][]const u8{ "../share/qemu", "share" }) |rel| {
        const path = b.pathJoin(&.{ bin, rel, "edk2-i386-vars.fd" });
        if (rstdbuild.io.fileExists(b.graph.io, path)) return path;
    }
    return null;
}

/// Whether `path`, relative to the build root, exists.
fn exists(b: *Build, path: []const u8) !bool {
    return rstdbuild.io.fileExists(b.graph.io, try b.root.joinString(b.allocator, path));
}
