//! `qemu-test` step: runs rstd's hardware-independent tests on emulated
//! machines with `qemu-system-*`, one freestanding test kernel per
//! architecture (boot code in `src/_qemutest/`, runner in `src/_test_qemu_runner.zig`, test root in
//! `src/_test_qemu.zig`). Needs `qemu-system-<arch>` on `PATH`.
//!
//! Where the machine has semihosting (aarch64, arm, riscv64, mips64el) the
//! runner exits with the number of failed tests, so a failing test fails the
//! step. Elsewhere (loongarch64, s390x, powerpc64le) the machine just powers
//! off and the step checks the printed summary for "0 failed".

const std = @import("std");
const rstd = @import("../src/__root__.zig");
const rstdbuild = rstd.buildutils;
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const blaskern = @import("blaskern.zig");

const Machine = struct {
    /// Name for step and executable.
    name: []const u8,
    query: std.Target.Query,
    /// `qemu-system-<arch>` binary.
    qemu: []const u8,
    /// Machine and CPU selection plus whatever boots a bare `-kernel` ELF.
    args: []const []const u8,
    linker_script: []const u8,
    /// Use the LLVM backend instead of the self-hosted one for code generation.
    llvm: bool = false,
    /// Overrides the build's optimize mode; OpenBIOS maps about 4 MB of a sparc64
    /// image, which debug builds of the tests outgrow. sparc64 uses `.fast`, as
    /// `.small` from Zig 0.18.0-dev.250 jumps into unrelated code in the gemm
    /// parameter test (an LLVM `-Oz` miscompile on sparc64, it appears).
    optimize: ?Optimize = null,
    code_model: std.builtin.CodeModel = .default,
    memory: []const u8 = "128M",
    /// Exit status carries the result; otherwise the printed summary does.
    semihosting: bool = true,
    /// Text qemu always prints to stderr here. The build treats any unchecked
    /// stderr as a diagnostic (a "failed command" echo on a passing step), so
    /// expecting it keeps a passing run quiet.
    stderr_note: ?[]const u8 = null,
};

const machines = [_]Machine{
    .{
        .name = "aarch64",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .freestanding, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.cortex_a72 } },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a72" },
        .linker_script = "src/_qemutest/aarch64.ld",
    },
    .{
        .name = "aarch64-sve",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .freestanding, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.neoverse_v1 } },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "neoverse-v1" },
        .linker_script = "src/_qemutest/aarch64.ld",
    },
    .{
        .name = "arm",
        .query = .{ .cpu_arch = .arm, .os_tag = .freestanding, .abi = .eabihf, .cpu_model = .{ .explicit = &std.Target.arm.cpu.cortex_a15 } },
        .qemu = "qemu-system-arm",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a15" },
        .linker_script = "src/_qemutest/arm.ld",
    },
    .{
        .name = "riscv64",
        .query = .{ .cpu_arch = .riscv64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-riscv64",
        // the vector extension, which the RVV kernels need
        .args = &.{ "-M", "virt", "-cpu", "rv64,v=true,vlen=256", "-smp", "1", "-bios", "none" },
        .linker_script = "src/_qemutest/riscv64.ld",
        // The kernel sits at 0x80000000, out of reach of the default `medlow`.
        .code_model = .medany,
    },
    .{
        .name = "loongarch64",
        .query = .{ .cpu_arch = .loongarch64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-loongarch64",
        .args = &.{ "-M", "virt", "-cpu", "la464" },
        .linker_script = "src/_qemutest/loongarch64.ld",
        .memory = "256M",
        .semihosting = false,
    },
    .{
        .name = "mips64r2el",
        .query = .{ .cpu_arch = .mips64el, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips64r2 } },
        .qemu = "qemu-system-mips64el",
        .args = &.{ "-M", "malta", "-cpu", "Loongson-3A1000" },
        .linker_script = "src/_qemutest/mips64.ld",
        .memory = "256M",
    },
    .{
        .name = "mips64r6el",
        .query = .{ .cpu_arch = .mips64el, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips64r6 } },
        .qemu = "qemu-system-mips64el",
        .args = &.{ "-M", "malta", "-cpu", "I6400" },
        .linker_script = "src/_qemutest/mips64.ld",
        .memory = "256M",
    },
    .{
        .name = "mips32r2el",
        .query = .{ .cpu_arch = .mipsel, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips32r2 } },
        .qemu = "qemu-system-mipsel",
        .args = &.{ "-M", "malta", "-cpu", "24Kf" },
        .linker_script = "src/_qemutest/mips32.ld",
        .memory = "256M",
    },
    .{
        .name = "mips32r5el",
        .query = .{ .cpu_arch = .mipsel, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.p5600 } },
        .qemu = "qemu-system-mipsel",
        .args = &.{ "-M", "malta", "-cpu", "P5600" },
        .linker_script = "src/_qemutest/mips32.ld",
        .memory = "256M",
    },
    .{
        .name = "i386",
        .query = .{ .cpu_arch = .x86, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.pentium3 } },
        .qemu = "qemu-system-i386",
        .args = &.{ "-M", "pc", "-cpu", "pentium3" },
        .linker_script = "src/_qemutest/x86.ld",
        .memory = "128M",
        .semihosting = false,
    },
    .{
        .name = "x86_64",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.haswell } },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "pc", "-cpu", "max" },
        .linker_script = "src/_qemutest/x86_64.ld",
        .llvm = true,
        .memory = "128M",
        .semihosting = false,
    },
    .{
        .name = "sparc64",
        .optimize = .fast,
        .query = .{ .cpu_arch = .sparc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.sparc.cpu.v9 } },
        .qemu = "qemu-system-sparc64",
        .args = &.{ "-M", "sun4u" },
        .linker_script = "src/_qemutest/sparc64.ld",
        .memory = "256M",
        .semihosting = false,
    },
    .{
        .name = "powerpc64le",
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr9 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power9" },
        .linker_script = "src/_qemutest/ppc64le.ld",
        .memory = "512M",
        .semihosting = false,
        // the firmware loader's 0x400000 shift, see the linker script
        .stderr_note = "kernel_addr changed",
    },
    .{
        .name = "powerpc64le-power10",
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr10 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power10" },
        .linker_script = "src/_qemutest/ppc64le.ld",
        .memory = "512M",
        .semihosting = false,
        .stderr_note = "kernel_addr changed",
    },
    .{
        .name = "s390x",
        .query = .{ .cpu_arch = .s390x, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.s390x.cpu.z14 } },
        .qemu = "qemu-system-s390x",
        .args = &.{ "-M", "s390-ccw-virtio", "-cpu", "max" },
        .linker_script = "src/_qemutest/s390x.ld",
        .memory = "256M",
        .semihosting = false,
    },
};

/// Tests that need an OS (allocator, files, threads) live in files the test
/// root also reaches, so the runnable ones are picked by name instead.
const default_filters: []const []const u8 = &.{ "linalg.blas.blas1.axpy", "linalg.blas.blas3.gemm.params", "simd." };

/// Adds the `qemu-test` step and its per-arch `qemu-test-<arch>` steps.
pub fn addQemuTests(b: *Build, optimize: Optimize) void {
    const all = b.step("qemu-test", "Run rstd's hardware-independent tests under qemu-system-* (all archs)");
    const filters = b.option([]const []const u8, "qemu-test-filter", "Only run tests whose name contains this (repeatable)") orelse &.{};

    for (machines) |m| {
        const target = b.resolveTargetQuery(m.query);
        const test_mod = b.createModule(.{
            .root_source_file = b.path("src/_test_qemu.zig"),
            .target = target,
            .optimize = m.optimize orelse optimize,
            .code_model = m.code_model,
        });
        // the kernels of the architectures that have any
        blaskern.addBlasKernels(b, test_mod, target, .all);
        const test_exe = b.addTest(.{
            .name = b.fmt("rstd-qemu-{s}", .{m.name}),
            .root_module = test_mod,
            .test_runner = .{ .path = b.path("src/_test_qemu_runner.zig"), .mode = .simple },
            .filters = if (filters.len != 0) filters else default_filters,
            .use_llvm = if (m.llvm) true else null,
            // every rig needs its linker script, which only lld honors
            .use_lld = true,
        });
        test_exe.setLinkerScript(b.path(m.linker_script));

        const run = b.addSystemCommand(&.{m.qemu});
        run.addArgs(m.args);
        run.addArgs(&.{ "-m", m.memory, "-nographic", "-no-reboot" });
        if (m.semihosting) run.addArgs(&.{ "-semihosting-config", "enable=on,target=native" });
        run.addArg("-kernel");
        run.addArtifactArg(test_exe);
        run.has_side_effects = true;
        if (!m.semihosting) run.expectStdOutMatch(", 0 failed");
        if (m.stderr_note) |note| run.expectStdErrMatch(note);

        const step = b.step(b.fmt("qemu-test-{s}", .{m.name}), b.fmt("Run rstd's tests under qemu-system on {s}", .{m.name}));
        step.dependOn(&run.step);
        all.dependOn(&run.step);
    }
}
