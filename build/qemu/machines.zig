//! The qemu rigs (`__root__.zig`): each machine, how its tests boot, and the test filters they run.

const std = @import("std");
const rstdbuild = @import("../../src/buildutils/__root__.zig");
const Optimize = rstdbuild.Optimize;
const Build = rstdbuild.Build;

pub const Machine = struct {
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
    /// `args` under HVF, on a `-Dtest-host` Mac of the rig's architecture
    /// (`-Dtest-host-arch`), unless `-Dtest-host-tcg`; null where the rig
    /// needs a CPU the host lacks.
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

    /// Whether QEMU's semihosting carries a kernel's exit status, else its printed
    /// summary does: on Arm, RISC-V and MIPS, but not under UEFI firmware.
    pub fn semihosting(m: Machine) bool {
        return m.boot.kernel.firmware == null and switch (m.query.cpu_arch.?) {
            .aarch64, .arm, .riscv64, .mips, .mipsel, .mips64, .mips64el => true,
            else => false,
        };
    }
};

pub const Part = struct { name: []const u8, filters: []const []const u8 };

/// How the tests reach a machine.
pub const Boot = union(enum) {
    kernel: Kernel,
    /// An installed macOS, set up once (`macos.zig`).
    macos,
    /// A Linux kernel booted with an initramfs holding rstd's full test suite.
    linux: Linux,

    /// A freestanding test kernel, booted bare with `-kernel`.
    pub const Kernel = struct {
        linker_script: []const u8,
        /// Use the LLVM backend instead of the self-hosted one for code generation.
        llvm: bool = false,
        code_model: std.lang.CodeModel = .default,
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

    pub const Window = struct { base: u64, size: u64 };

    /// A Debian netboot kernel (a lazy dependency) and an initramfs of
    /// `test/linux/init.zig` as PID 1 and the hosted test executable.
    pub const Linux = struct {
        /// The `build.zig.zon` dependency, and the kernel's path in it.
        dep: []const u8,
        kernel: []const u8,
        /// The kernel's console, on the machine's first serial port.
        console: []const u8,

        /// The kernel command line; `panic=-1`: a dead init restarts at once,
        /// which `-no-reboot` makes an exit.
        pub fn commandLine(l: Linux, b: *Build) []const u8 {
            return b.fmt("console={s} rdinit=/init panic=-1 quiet loglevel=3", .{l.console});
        }
    };
};

/// Tier 4 (user, 2026-11-01): the images outgrow OpenBIOS's window faster
/// than parts split them (rstd issue #131).
const sparc64_stub = "tier 4: its images outgrow OpenBIOS's 4 MiB window (rstd issue #131)";

/// Every rig.
pub const all = base_machines ++ x87Safe(&.{ "i386", "i386-max" }) ++ quad(&.{"aarch64"});

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
        } },
    },
    .{
        .name = "x86_64",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.haswell } },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "pc", "-cpu", "max" },
        .hvf = &.{ "-accel", "hvf", "-M", "pc", "-cpu", "host" },
        .memory = "128M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/x86_64.ld",
            .llvm = true,
        } },
    },
    .{
        // the threaded tests on four cores, started by INIT-SIPI-SIPI
        .name = "x86_64-smp",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.haswell } },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "pc", "-cpu", "max" },
        .hvf = &.{ "-accel", "hvf", "-M", "pc", "-cpu", "host" },
        .memory = "128M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/x86_64.ld",
            .llvm = true,
        } },
        .cpus = 4,
    },
    .{
        // The host CPU under Windows' hypervisor: TCG has no AVX-512, FMA4 or
        // XOP (issue #220), so their cores run only here, on a host that has
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
        } },
    },
    .{
        .name = "sparc64",
        .stub = sparc64_stub,
        .optimize = .fast,
        // OpenBIOS maps 0x400000 to 0x800000 only, and the default set's
        // code alone is about 5 MB; blas2 by itself is 3.3 MB
        .parts = &.{
            .{ .name = "blas1", .filters = &.{ "linalg.blas.blas1.", "linalg.blas._test.blas1.", "linalg.blas._test.check.", "linalg.blas._test.bitwise." } },
            .{ .name = "blas2-dense", .filters = &.{ "linalg.blas.blas2.__root__.", "linalg.blas.blas2.gemv.", "linalg.blas._test.bf16.", "linalg.blas.blas2.ger.", "linalg.blas.blas2.hemv.", "linalg.blas.blas2.symv." } },
            .{ .name = "blas2-banded", .filters = &.{ "linalg.blas.blas2.banded.", "linalg.blas.blas2.symmetric." } },
            .{ .name = "blas2-triangular", .filters = &.{"linalg.blas.blas2.triangular."} },
            .{ .name = "blas2-interface", .filters = &.{"linalg.blas.blas2.interface."} },
            .{ .name = "blas2-edge", .filters = &.{"linalg.blas._test.blas2.edge."} },
            // one per routine family: together past the window
            .{ .name = "blas3-gemm", .filters = &.{ "linalg.blas.blas3.__root__.", "linalg.blas.blas3.gemm." } },
            .{ .name = "blas3-symm", .filters = &.{"linalg.blas.blas3.symm."} },
            .{ .name = "blas3-rank", .filters = &.{ "linalg.blas.blas3.syrk.", "linalg.blas.blas3.syr2k." } },
            .{ .name = "blas3-tri", .filters = &.{ "linalg.blas.blas3.trmm.", "linalg.blas.blas3.trsm.", "linalg.lapack." } },
            .{ .name = "blas3-interface", .filters = &.{"linalg.blas.blas3.interface."} },
            .{ .name = "openblas", .filters = &.{"linalg.blas._test.from_openblas."} },
            .{ .name = "rest", .filters = &.{ "linalg.blas.__root__.", "linalg.blas.Config.", "linalg.blas.cpu." } },
            .{ .name = "small", .filters = &small_filters },
        },
        .query = .{ .cpu_arch = .sparc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.sparc.cpu.v9 } },
        .qemu = "qemu-system-sparc64",
        .args = &.{ "-M", "sun4u" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/sparc64.ld",
            .load_window = .{ .base = 0x400000, .size = 0x400000 },
        } },
    },
    .{
        // sparc64 in debug, the mode the other rigs build in; its images are
        // larger, so blas3 is split too, and gemm's three test files apart
        .name = "sparc64-debug",
        .stub = sparc64_stub,
        .optimize = .debug,
        .parts = &.{
            // blas1 in three, as two are past the window
            .{ .name = "blas1", .filters = &.{ "linalg.blas.blas1.__root__.", "linalg.blas.blas1.asum.", "linalg.blas.blas1.axpby.", "linalg.blas.blas1.axpy.", "linalg.blas._test.check.", "linalg.blas._test.bitwise.", "linalg.blas.blas1.copy.", "linalg.blas.blas1.dot.", "linalg.blas.blas1.extremum.", "linalg.blas._test.blas1.extremum.", "linalg.blas.blas1.nrm2.", "linalg.blas.blas1.rot.", "linalg.blas.blas1.rotm.", "linalg.blas.blas1.scal.", "linalg.blas.blas1.sum.", "linalg.blas.blas1.swap.", "linalg.blas.uk_block." } },
            .{ .name = "blas1-interface", .filters = &.{ "linalg.blas.blas1.interface.", "linalg.blas._test.blas1.threaded." } },
            .{ .name = "blas1-edge", .filters = &.{ "linalg.blas._test.blas1.edge.test.portable level-1 kernels: special", "linalg.blas._test.blas1.edge.test.every linked core's level-1 kernels: special", "linalg.blas._test.blas1.edge.test.portable level-1 reductions", "linalg.blas._test.blas1.edge.test.every linked core's level-1 reductions", "linalg.blas._test.blas1.edge.test.every linked core keeps" } },
            .{ .name = "blas1-denormal", .filters = &.{ "linalg.blas._test.blas1.edge.test.portable level-1 kernels: denormal", "linalg.blas._test.blas1.edge.test.every linked core's level-1 kernels: denormal", "linalg.blas._test.blas1.edge.test.where util.vectorIeee" } },
            // the dense routines, the triangular ones, the interface and the
            // complex edge cases each in two, as one is past the window
            .{ .name = "blas2-dense", .filters = &.{ "linalg.blas.blas2.__root__.", "linalg.blas.blas2.gemv.__root__.test.portable", "linalg.blas.blas2.gemv.__root__.test_" } },
            .{ .name = "blas2-gemv-x87", .filters = &.{ "linalg.blas._test.bf16.gemv.", "linalg.blas.blas2.gemv.extended.test.extended gemv:" } },
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
            .{ .name = "blas3-symm", .filters = &.{"linalg.blas.blas3.symm."} },
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
            .{ .name = "openblas-cores-l1-blat1", .filters = &.{"linalg.blas._test.from_openblas.__root__.test.blat1 on every"} },
            .{ .name = "openblas-cores-l1-utest", .filters = &.{"linalg.blas._test.from_openblas.__root__.test.utest and utest_ext's level-1 cases"} },
            .{ .name = "openblas-cores-l1-utest-2", .filters = &.{"linalg.blas._test.from_openblas.__root__.test.utest_ext's level-2 cases"} },
            .{ .name = "openblas-cores-l2-real", .filters = &.{ "linalg.blas._test.from_openblas.__root__.test.blat2 at f32", "linalg.blas._test.from_openblas.__root__.test.blat2 at f64" } },
            .{ .name = "openblas-cores-l2-complex", .filters = &.{"linalg.blas._test.from_openblas.__root__.test.blat2 at complex"} },
            .{ .name = "openblas-cores-l3", .filters = &.{"linalg.blas._test.from_openblas.__root__.test.blat3, utest's LAPACK cases"} },
            // LAPACK in five, as one is past the window
            .{ .name = "lapack-cholesky", .filters = &.{"linalg.lapack.cholesky."} },
            .{ .name = "lapack-getrf", .filters = &.{"linalg.lapack.lu.test.getrf"} },
            .{ .name = "lapack-getrs-real", .filters = &.{ "linalg.lapack.lu.test.getrs: op(A) X = B on getrf's factor, every op, NaN in one column of B, real", "linalg.lapack.lu.test.getrs: a known", "linalg.lapack.lu.test.getrs: reports" } },
            .{ .name = "lapack-getrs-complex", .filters = &.{"linalg.lapack.lu.test.getrs: op(A) X = B on getrf's factor, every op, NaN in one column of B, complex"} },
            .{ .name = "lapack-faer", .filters = &.{ "linalg.lapack.ops.", "linalg.lapack.interface.", "linalg.lapack.blas_calls.", "linalg.lapack.auxiliary.", "linalg.lapack.householder.", "linalg.lapack.qr.", "linalg.lapack.qr_solve.", "linalg.lapack.tridiag.", "linalg.lapack.tridiag_evd.", "linalg.lapack.hessenberg.", "linalg.lapack.real_schur.", "linalg.lapack.complex_schur.", "linalg.lapack.balance.", "linalg.lapack.bidiag.", "linalg.lapack.bidiag_svd." } },
            .{ .name = "rest", .filters = &.{ "linalg.blas.__root__.", "linalg.blas.Config.", "linalg.blas.cpu." } },
            .{ .name = "small", .filters = &small_filters },
        },
        .query = .{ .cpu_arch = .sparc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.sparc.cpu.v9 } },
        .qemu = "qemu-system-sparc64",
        .args = &.{ "-M", "sun4u" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/sparc64.ld",
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
        .stub = "lld cannot link 32-bit SPARC, and Zig's own ELF linker crashes on it (issue #189)",
        .boot = .{ .kernel = .{
            .linker_script = "",
        } },
    },
    .{
        .name = "powerpc64le",
        // LLVM PowerPC FastISel (debug) uses FPR instructions on VSX registers (issue #186)
        .optimize = .safe,
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr9 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power9" },
        .memory = "512M",
        .boot = .{
            .kernel = .{
                .linker_script = "test/qemu/ppc64.ld",
                // the firmware loader's 0x400000 shift, see the linker script
                .stderr_note = "kernel_addr changed",
            },
        },
    },
    .{
        .name = "powerpc64le-power10",
        // LLVM PowerPC FastISel (debug) uses FPR instructions on VSX registers (issue #186)
        .optimize = .safe,
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr10 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power10" },
        .memory = "512M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc64.ld",
            .stderr_note = "kernel_addr changed",
        } },
    },
    .{
        // big-endian, ELFv2 as Zig emits it; POWER10 runs every 64-bit core
        .name = "powerpc64",
        // LLVM PowerPC FastISel (debug) uses FPR instructions on VSX registers (issue #186)
        .optimize = .safe,
        .query = .{ .cpu_arch = .powerpc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr10 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power10" },
        .memory = "512M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc64.ld",
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
        } },
    },
    .{
        // RandyOS's iBook G3 target (a 750FX, no AltiVec) on `mac99`: soft
        // float, so no kernel objects; the portable kernels run. A 750, as
        // OpenBIOS hangs on QEMU's `750fx` (issue #225); the same ISA
        .name = "powerpc-750",
        .query = .{ .cpu_arch = .powerpc, .os_tag = .freestanding, .abi = .eabi, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.@"750" } },
        .qemu = "qemu-system-ppc",
        .args = &.{ "-M", "mac99", "-cpu", "750" },
        .memory = "256M",
        .blas_cores = .none,
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc32.ld",
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
            .stderr_note = "kernel_addr changed",
        } },
    },
    .{
        // a 32-bit VSX build on a POWER9 under VOF: the power8 and power9
        // cores' `@Vector` kernels as 32-bit code. The build drops `64bit`
        // and `fpcvt`, with which LLVM fails on compiler-rt in 32-bit mode;
        // ISA 3.0 fails it too, so power10's kernels cannot link (issue #187).
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
        } },
    },
    .{
        // RandyOS's loader's platforms: the tests as a UEFI application under
        // QEMU's OVMF and AAVMF, on rstd's UEFI `std.Io`
        .name = "x86_64-uefi",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .uefi },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "q35", "-cpu", "max" },
        .hvf = &.{ "-accel", "hvf", "-M", "q35", "-cpu", "host" },
        .memory = "1G",
        .boot = .{ .kernel = .{
            .linker_script = "",
            .llvm = true,
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
        .hvf = &.{ "-accel", "hvf", "-M", "pc", "-cpu", "host", "-smp", "4" },
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
        // POWER rigs build `.safe` for (issue #186)
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
        // native `f128` from pwr9: FastISel's bug, as the freestanding rigs (issue #186)
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
        .stub = "TCG's fldt/fstpt misplace the significand bytes past a page boundary, which i386's 4-aligned f80 arrays cross (issue #221)",
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
pub const default_filters: []const []const u8 = &([_][]const u8{ "linalg.blas.__root__.", "linalg.blas.Config.", "linalg.blas.cpu.", "linalg.blas.blas1.", "linalg.blas.blas2.", "linalg.blas.blas3.", "linalg.blas._test.", "linalg.lapack." } ++ small_filters);

/// The multicore rigs' tests: those that thread.
pub const smp_filters: []const []const u8 = &.{ "linalg.blas.threading.", "linalg.blas.blas3.threaded.", "linalg.blas._test.blas3.threaded.", "linalg.blas._test.blas2.threaded.", "threaded level 1", "at any worker count" };

/// The smaller BLAS test files: types, `fma`, `portable`, `ext`, `util` and `FpState`.
const small_filters = [_][]const u8{ "linalg.blas._test.types.", "linalg.blas.fma.", "linalg.blas.portable.", "linalg.blas.ext.", "linalg.blas.util.", "linalg.blas.FpState." };
