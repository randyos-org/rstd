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
//!     logged as it runs, to `.zig-cache/qemu-test-<name>.log`.
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
    /// Windows host.
    whpx: bool = false,
    /// Test sets each run as a kernel of its own, `<name>-<part>`, for a
    /// machine whose firmware maps too little for all the tests in one image;
    /// `-Dqemu-test-filter` builds one image as usual.
    parts: []const Part = &.{},
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
        code_model: std.builtin.CodeModel = .default,
        /// Exit status carries the result; otherwise the printed summary does.
        semihosting: bool = true,
        /// Text qemu always prints to stderr here. The build treats any unchecked
        /// stderr as a diagnostic (a "failed command" echo on a passing step), so
        /// expecting it keeps a passing run quiet.
        stderr_note: ?[]const u8 = null,
    };

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

const machines = [_]Machine{
    .{
        .name = "aarch64",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .freestanding, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.cortex_a72 } },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a72" },
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/aarch64.ld",
        } },
    },
    .{
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
        .query = .{ .cpu_arch = .arm, .os_tag = .freestanding, .abi = .eabihf, .cpu_model = .{ .explicit = &std.Target.arm.cpu.cortex_a15 } },
        .qemu = "qemu-system-arm",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a15" },
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/arm.ld",
        } },
    },
    .{
        // the soft-float ABI: floats in integer registers, the kernels' `sf` variant
        .name = "arm-softfp",
        .query = .{ .cpu_arch = .arm, .os_tag = .freestanding, .abi = .eabi, .cpu_model = .{ .explicit = &std.Target.arm.cpu.cortex_a15 } },
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
        .args = &.{ "-M", "malta", "-cpu", "Loongson-3A1000" },
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
            .{ .name = "blas1", .filters = &.{"linalg.blas.blas1."} },
            .{ .name = "blas2-dense", .filters = &.{ "linalg.blas.blas2.__root__.", "linalg.blas.blas2.gemv.", "linalg.blas.blas2.ger.", "linalg.blas.blas2.hemv.", "linalg.blas.blas2.symv." } },
            .{ .name = "blas2-banded", .filters = &.{ "linalg.blas.blas2.banded.", "linalg.blas.blas2.symmetric." } },
            .{ .name = "blas2-triangular", .filters = &.{"linalg.blas.blas2.triangular."} },
            .{ .name = "blas2-interface", .filters = &.{ "linalg.blas.blas2.interface.", "linalg.blas.blas2.edge." } },
            .{ .name = "blas3", .filters = &.{"linalg.blas.blas3."} },
            .{ .name = "openblas", .filters = &.{"linalg.blas._test_openblas."} },
            .{ .name = "rest", .filters = &.{ "linalg.blas.__root__.", "linalg.blas.Config.", "linalg.blas.cpu.", "simd." } },
        },
        .query = .{ .cpu_arch = .sparc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.sparc.cpu.v9 } },
        .qemu = "qemu-system-sparc64",
        .args = &.{ "-M", "sun4u" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/sparc64.ld",
            .semihosting = false,
        } },
    },
    .{
        // sparc64 in debug, the mode the other rigs build in; its images are
        // larger, so blas3 is split too, and gemm's three test files apart
        .name = "sparc64-debug",
        .optimize = .debug,
        .parts = &.{
            .{ .name = "blas1", .filters = &.{"linalg.blas.blas1."} },
            .{ .name = "blas2-dense", .filters = &.{ "linalg.blas.blas2.__root__.", "linalg.blas.blas2.gemv.", "linalg.blas.blas2.ger.", "linalg.blas.blas2.hemv.", "linalg.blas.blas2.symv." } },
            .{ .name = "blas2-banded", .filters = &.{ "linalg.blas.blas2.banded.", "linalg.blas.blas2.symmetric." } },
            .{ .name = "blas2-triangular", .filters = &.{"linalg.blas.blas2.triangular."} },
            .{ .name = "blas2-interface", .filters = &.{"linalg.blas.blas2.interface."} },
            // in two, as one is past the window
            .{ .name = "blas2-edge-real", .filters = &.{ "linalg.blas.blas2.edge.test.portable level-2 routines, real", "linalg.blas.blas2.edge.test.every linked core's level-2 routines, real" } },
            .{ .name = "blas2-edge-complex", .filters = &.{ "linalg.blas.blas2.edge.test.portable level-2 routines, complex", "linalg.blas.blas2.edge.test.every linked core's level-2 routines, complex" } },
            .{ .name = "blas3-gemm", .filters = &.{"linalg.blas.blas3.gemm.__root__."} },
            // but "every core of every family", 6 MB alone in debug: `sparc64` runs it
            .{ .name = "blas3-params", .filters = &.{
                "linalg.blas.blas3.gemm.params.test.dtb entries",
                "linalg.blas.blas3.gemm.params.test.generic follows",
                "linalg.blas.blas3.gemm.params.test.the block factor",
                "linalg.blas.blas3.gemm.params.test.fixed sets",
                "linalg.blas.blas3.gemm.params.test.old x86",
                "linalg.blas.blas3.gemm.params.test.LA464",
                "linalg.blas.blas3.gemm.params.test.detection gives",
            } },
            .{ .name = "blas3-kernel", .filters = &.{ "linalg.blas.blas3.__root__.", "linalg.blas.blas3.gemm.kernel." } },
            .{ .name = "blas3-interface", .filters = &.{"linalg.blas.blas3.interface."} },
            // the OpenBLAS suites in three; not "every linked core", which
            // holds all of them (4.7 MB in debug): `sparc64` runs it
            .{ .name = "openblas-l1", .filters = &.{ "linalg.blas._test_openblas.blat1.", "linalg.blas._test_openblas.utest.", "linalg.blas._test_openblas.utest_ext." } },
            .{ .name = "openblas-blat2", .filters = &.{"linalg.blas._test_openblas.blat2."} },
            .{ .name = "openblas-blat3", .filters = &.{"linalg.blas._test_openblas.blat3."} },
            .{ .name = "rest", .filters = &.{ "linalg.blas.__root__.", "linalg.blas.Config.", "linalg.blas.cpu.", "simd." } },
        },
        .query = .{ .cpu_arch = .sparc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.sparc.cpu.v9 } },
        .qemu = "qemu-system-sparc64",
        .args = &.{ "-M", "sun4u" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/sparc64.ld",
            .semihosting = false,
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
        // LLVM PowerPC FastISel (debug) compares f64 in the wrong register (issues.md)
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
        // LLVM PowerPC FastISel (debug) compares f64 in the wrong register (issues.md)
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
        // LLVM PowerPC FastISel (debug) compares f64 in the wrong register (issues.md)
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
        .name = "linux-aarch64",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .linux },
        .qemu = "qemu-system-aarch64",
        // four CPUs, for the tests that race threads
        .args = &.{ "-M", "virt", "-cpu", "max", "-smp", "4" },
        .memory = "1G",
        .boot = .{ .linux = .{ .dep = "debian_netboot_arm64", .kernel = "debian-installer/arm64/linux", .console = "ttyAMA0" } },
    },
    .{
        .name = "linux-x86_64",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .linux },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "pc", "-cpu", "max", "-smp", "4" },
        .memory = "1G",
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
        // LLVM PowerPC FastISel (debug) compares f64 in the wrong register (issues.md)
        .optimize = .safe,
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .linux },
        .qemu = "qemu-system-ppc64",
        // detection sees a POWER10 (the baseline build links up to POWER9's
        // kernels); the caps TCG lacks are off, so qemu warns of none
        .args = &.{ "-M", "pseries,cap-cfpc=broken,cap-sbbc=broken,cap-ibs=broken,cap-ccf-assist=off", "-cpu", "power10", "-smp", "4" },
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
const default_filters: []const []const u8 = &.{ "linalg.blas.__root__.", "linalg.blas.Config.", "linalg.blas.cpu.", "linalg.blas.blas1.", "linalg.blas.blas2.", "linalg.blas.blas3.", "simd." };

/// Adds the `qemu-test` step and each machine's `qemu-test-<name>` step.
pub fn addQemuTests(b: *Build, optimize: Optimize) !void {
    const all = b.step("qemu-test", "Run rstd's hardware-independent tests under qemu-system-* (all archs)");
    const filters = b.option([]const []const u8, "qemu-test-filter", "Only run tests whose name contains this (repeatable)") orelse &.{};
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

    const host = b.graph.host.result;
    const whpx_host = host.os.tag == .windows and host.cpu.arch == .x86_64;
    for (machines) |m| {
        if (m.stub != null or (m.whpx and !whpx_host)) continue;
        const target = b.resolveTargetQuery(m.query);
        if (m.boot == .kernel and m.parts.len != 0 and filters.len == 0) {
            const step = b.step(b.fmt("qemu-test-{s}", .{m.name}), b.fmt("Run rstd's tests under qemu-system on {s}, in {d} images (live logs: .zig-cache/qemu-test-{s}-<part>.log)", .{ m.name, m.parts.len, m.name }));
            for (m.parts) |p| {
                var part = m;
                part.name = b.fmt("{s}-{s}", .{ m.name, p.name });
                const run = addKernelRun(b, part, m.boot.kernel, target, m.optimize orelse optimize, p.filters);
                step.dependOn(&run.step);
                all.dependOn(&run.step);
            }
            continue;
        }
        const run = switch (m.boot) {
            .kernel => |k| addKernelRun(b, m, k, target, m.optimize orelse optimize, filters),
            .macos => if (macos_vm) try addMacosRun(b, m, target, m.optimize orelse optimize, mac) orelse continue else continue,
            .linux => |l| if (linux_vm) try addLinuxRun(b, m, l, target, m.optimize orelse optimize, filters) orelse continue else continue,
        };
        const about = switch (m.boot) {
            .kernel, .linux => b.fmt("Run rstd's tests under qemu-system on {s} (live log: {s})", .{ m.name, liveLog(b, m) }),
            .macos => b.fmt("Run rstd's tests under qemu-system on {s}", .{m.name}),
        };
        const step = b.step(b.fmt("qemu-test-{s}", .{m.name}), about);
        step.dependOn(&run.step);
        all.dependOn(&run.step);
    }
}

/// Whether to add a family of VM rigs: on when rstd is the root project, so a
/// dependent never fetches what they need; `-D<name>-off` turns them off, and
/// `-D<name>` on in a dependent. `-off` wins.
fn vmOption(b: *Build, comptime name: []const u8, comptime rigs: []const u8, comptime needs: []const u8) bool {
    const on = b.option(bool, name, "Add " ++ rigs ++ ", fetching " ++ needs ++ " (default: when rstd is the root project)");
    const off = b.option(bool, name ++ "-off", "Leave out " ++ rigs ++ " and do not fetch " ++ needs) orelse false;
    return !off and (on orelse b.isRoot());
}

/// `m`'s qemu, under `tools/qemu_watch.zig`: a firmware halt or an hour with
/// no output ends it and fails the step, where it would otherwise never end.
/// The slowest test prints nothing for about 20 minutes on the x86 rigs.
fn qemuCommand(b: *Build, m: Machine) *Build.Step.Run {
    const run = b.addRunArtifact(watchTool(b));
    run.addArgs(&.{ "60", "Stopping execution", "--", m.qemu });
    return run;
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
/// monitor keeps `-nographic`'s firmware console.
fn addConsole(b: *Build, run: *Build.Step.Run, m: Machine) void {
    run.addArgs(&.{ "-m", m.memory, "-nographic", "-no-reboot" });
    run.addArgs(&.{ "-chardev", b.fmt("stdio,id=con,mux=on,signal=off,logfile={s}", .{liveLog(b, m)}), "-serial", "chardev:con", "-mon", "chardev=con" });
}

/// QEMU booting a Debian kernel with an initramfs of `test/linux/init.zig`
/// and rstd's hosted tests, which see real auxv, sysfs, procfs and threads;
/// null until the kernel's package is fetched.
fn addLinuxRun(b: *Build, m: Machine, l: Boot.Linux, target: rstdbuild.ResolvedTarget, optimize: Optimize, filters: []const []const u8) !?*Build.Step.Run {
    const kernel = b.lazyDependency(l.dep, .{}) orelse return null;
    const mod = try buildroot.addRstdModule(b, null, target, optimize, .all, false, false);
    const test_exe = try buildroot.addRstdTest(b, mod, target, filters);
    const init_exe = b.addExecutable(.{ .name = b.fmt("init-{s}", .{m.name}), .root_module = b.createModule(.{
        .root_source_file = b.path("test/linux/init.zig"),
        .target = target,
        .optimize = optimize,
    }) });
    const tool = b.addExecutable(.{ .name = "initramfs", .root_module = b.createModule(.{
        .root_source_file = b.path("tools/initramfs.zig"),
        .target = b.graph.host,
    }) });
    const pack = b.addRunArtifact(tool);
    const cpio = pack.addOutputFileArg2("initramfs.cpio", .{});
    pack.addPrefixedArtifactArg("init=", init_exe);
    pack.addPrefixedArtifactArg("test=", test_exe);

    const run = qemuCommand(b, m);
    run.addArgs(m.args);
    addConsole(b, run, m);
    run.addArg("-kernel");
    run.addFileArg(kernel.path(l.kernel));
    run.addArg("-initrd");
    run.addFileArg(cpio);
    // panic=-1: a dead init restarts at once, which -no-reboot makes an exit
    run.addArgs(&.{ "-append", b.fmt("console={s} rdinit=/init panic=-1 quiet loglevel=3", .{l.console}) });
    run.has_side_effects = true;
    // `passed` in test/linux/init.zig
    run.expectStdOutMatch("linuxtest: /test exited 0");
    return run;
}

/// QEMU booting `m`'s freestanding test kernel.
fn addKernelRun(b: *Build, m: Machine, k: Boot.Kernel, target: rstdbuild.ResolvedTarget, optimize: Optimize, filters: []const []const u8) *Build.Step.Run {
    const test_mod = b.createModule(.{
        .root_source_file = b.path("src/_test_qemu.zig"),
        .target = target,
        .optimize = optimize,
        .code_model = k.code_model,
        // no threads under a test kernel
        .single_threaded = true,
    });
    // the kernels of the architectures that have any
    blaskern.addBlasKernels(b, test_mod, target, .all);
    const test_exe = b.addTest(.{
        .name = b.fmt("rstd-qemu-{s}", .{m.name}),
        .root_module = test_mod,
        .test_runner = .{ .path = b.path("test/qemu/runner.zig"), .mode = .simple },
        .filters = if (filters.len != 0) filters else default_filters,
        .use_llvm = if (k.llvm) true else null,
        // every rig needs its linker script, which only lld honors
        .use_lld = true,
    });
    test_exe.setLinkerScript(b.path(k.linker_script));

    const run = qemuCommand(b, m);
    run.addArgs(m.args);
    addConsole(b, run, m);
    if (k.semihosting) run.addArgs(&.{ "-semihosting-config", "enable=on,target=native,chardev=con" });
    run.addArg("-kernel");
    run.addArtifactArg(test_exe);
    run.has_side_effects = true;
    if (!k.semihosting) run.expectStdOutMatch(", 0 failed");
    if (k.stderr_note) |note| run.expectStdErrMatch(note);
    return run;
}

/// `m`'s `qemu-setup-`, `qemu-install-` and `qemu-start-<name>` steps, and the
/// test run: rstd's unit tests copied into the running VM and run over SSH.
/// Null until the lazy OSX-KVM dependency is fetched. The VM's files are named
/// relative to the build root, where its commands run, so none is hashed.
fn addMacosRun(b: *Build, m: Machine, target: rstdbuild.ResolvedTarget, optimize: Optimize, mac: MacOptions) !?*Build.Step.Run {
    const kvm = b.lazyDependency("osx_kvm", .{}) orelse return null;
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
    run.addArtifactArg(try buildroot.addRstdTest(b, mod, target, &.{}));
    run.addArg(m.qemu);
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
    const dir = std.fs.path.dirname(paths.disk).?;
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
    const start = (std.mem.indexOf(u8, script, "osk=\"") orelse return error.NoSmcKey) + 5;
    const end = std.mem.indexOfScalarPos(u8, script, start, '"') orelse return error.NoSmcKey;
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
    const bin = std.fs.path.dirname(qemu) orelse return null;
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
