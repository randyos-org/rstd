//! The qemu rigs (`__root__.zig`): each machine, how its tests boot, and the test filters they run.

pub const std = @import("std");
const rstdbuild = @import("../../src/buildutils/__root__.zig");
const arm_rigs = @import("machines/arm.zig");
const riscv_mips_rigs = @import("machines/riscv_mips.zig");
const x86_rigs = @import("machines/x86.zig");
const sparc_rigs = @import("machines/sparc.zig");
const power_rigs = @import("machines/power.zig");
const uefi_rigs = @import("machines/uefi.zig");
const linux_rigs = @import("machines/linux.zig");
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
    /// parameter test (an LLVM `-Oz` miscompile on sparc64, it appears; issue #183).
    optimize: ?Optimize = null,
    memory: []const u8 = "256M",
    /// `memory` under `-Drstd-tsan`, whose shadow memory and bigger test
    /// outgrow it (ENOMEM in faer's evd tests at 3G, issue #408)
    tsan_memory: ?[]const u8 = null,
    /// A Linux rig's `parts` under `-Drstd-tsan`, each a VM of its own, and a
    /// `rest` part of every test they leave out.
    tsan_parts: []const Part = &.{},
    /// As `tsan_parts`, in every run, kernel rigs too.
    split: []const Part = &.{},
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
    /// `args` under KVM, with `-Dqemu-kvm` on a Linux host of the rig's
    /// architecture (CI's runners); null where the rig needs a CPU the host
    /// lacks, or has not been tried under KVM.
    kvm: ?[]const []const u8 = null,
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

    /// nixpkgs' kernel (`tools/nix/kernel.sh`) and an initramfs of
    /// `test/linux/init.zig` as PID 1 and the hosted test executable.
    pub const Linux = struct {
        /// The kernel's Nix system.
        nix: []const u8,
        /// The kernel's console, on the machine's first serial port.
        console: []const u8,

        /// The kernel command line; `panic=-1`: a dead init restarts at once,
        /// which `-no-reboot` makes an exit.
        pub fn commandLine(l: Linux, b: *Build) []const u8 {
            return b.fmt("console={s} rdinit=/init panic=-1 quiet loglevel=3", .{l.console});
        }
    };
};

/// The Linux rigs' ThreadSanitizer runs, each part a VM: in one, TCG on CI's
/// arm64 runners passes the job's time limit, and x86_64 runs out of memory
/// after hours (issue #412); split as the kernel rigs' default set.
pub const tsan_parts = [_]Part{
    .{ .name = "blas1-2", .filters = &.{ "math.kernels.blas.blas1.", "math.kernels.blas.blas2.", "math.kernels.blas._test.blas1.", "math.kernels.blas._test.blas2." } },
    .{ .name = "blas3", .filters = &.{ "math.kernels.blas.blas3.", "math.kernels.blas._test.blas3." } },
    .{ .name = "lapack", .filters = &.{"math.kernels.lapack."} },
};

/// The SVE rigs' runs, each part a VM or image: in one, TCG on CI's arm64
/// runners takes 3.3 to 4.6 hours (issue #569). Balanced on `aarch64-sve`'s
/// test times in run 44702146209, 22 to 33 minutes a part.
pub const sve_parts = [_]Part{
    .{ .name = "blas1", .filters = &.{ "math.kernels.blas.blas1.", "math.kernels.blas._test.blas1.", "math.kernels.blas._test.from_openblas." } },
    .{ .name = "blas2-edge", .filters = &.{"math.kernels.blas._test.blas2.edge."} },
    .{ .name = "blas2", .filters = &.{
        "math.kernels.blas.blas2.",
        "math.kernels.blas._test.blas2.banded.",
        "math.kernels.blas._test.blas2.interface.",
        "math.kernels.blas._test.blas2.interface_any.",
        "math.kernels.blas._test.blas2.symcopy.",
        "math.kernels.blas._test.blas2.threaded.",
        "math.kernels.blas._test.blas2.triangular.",
    } },
    .{ .name = "blas3-tri", .filters = &.{ "math.kernels.blas.blas3.trmm.", "math.kernels.blas.blas3.trsm.", "math.kernels.blas.blas3.syrk.", "math.kernels.blas.blas3.syr2k.", "math.kernels.blas._test.blas3." } },
    .{ .name = "blas3", .filters = &.{ "math.kernels.blas.blas3.gemm.", "math.kernels.blas.blas3.gemm3m.", "math.kernels.blas.blas3.gemmt.", "math.kernels.blas.blas3.hemm.", "math.kernels.blas.blas3.hemm3m.", "math.kernels.blas.blas3.symm.", "math.kernels.blas.blas3.symm3m." } },
    .{ .name = "lapack-lu", .filters = &.{ "math.kernels.lapack._test.lu.", "math.kernels.lapack._test.lu_solve.", "math.kernels.lapack._test.cholesky." } },
};

/// Tier 4 (user, 2026-11-01): the images outgrow OpenBIOS's window faster
/// than parts split them (rstd issue #131).
pub const sparc64_stub = "tier 4: its images outgrow OpenBIOS's 4 MiB window (rstd issue #131)";

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

const base_machines = arm_rigs.list ++ riscv_mips_rigs.list ++ x86_rigs.list ++ sparc_rigs.list ++ power_rigs.list ++ uefi_rigs.list ++ linux_rigs.list;

/// Tests that need an OS (files, threads) live in files the test root also
/// reaches, so the runnable ones are picked by name; the runner supplies a heap.
/// Left out: `math.kernels.blas.abi.`, which resolves targets through `std.testing.io`.
pub const default_filters: []const []const u8 = &([_][]const u8{ "math.kernels.blas.__root__.", "math.kernels.blas.Config.", "math.kernels.blas.cpu.", "math.kernels.blas.blas1.", "math.kernels.blas.blas2.", "math.kernels.blas.blas3.", "math.kernels.blas._test.", "math.kernels.lapack." } ++ small_filters);

/// The multicore rigs' tests: those that thread.
pub const smp_filters: []const []const u8 = &.{ "math.kernels.blas.threading.", "math.kernels.blas.blas3.threaded.", "math.kernels.blas._test.blas3.threaded.", "math.kernels.blas._test.blas2.threaded.", "threaded level 1", "at any worker count" };

/// The smaller BLAS test files: types, `fma`, `portable`, `ext`, `util` and `FpState`.
pub const small_filters = [_][]const u8{ "math.kernels.blas._test.types.", "math.kernels.blas.fma.", "math.kernels.blas.portable.", "math.kernels.blas.ext.", "math.kernels.blas.util.", "math.kernels.blas.FpState." };
