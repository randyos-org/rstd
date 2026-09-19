//! `qemu-test` step: runs rstd's hardware-independent tests on emulated
//! machines with `qemu-system-*`, one freestanding test kernel per
//! architecture (boot code in `src/qemutest/`, runner in `src/qemutest_runner.zig`, test root in
//! `src/test_qemu.zig`). Needs `qemu-system-<arch>` on `PATH`.
//!
//! Where the machine has semihosting (aarch64, arm, riscv64, mips64el) the
//! runner exits with the number of failed tests, so a failing test fails the
//! step. Elsewhere (loongarch64, s390x, powerpc64le) the machine just powers
//! off and the step checks the printed summary for "0 failed".

const std = @import("std");
const rstd = @import("../src/__root__.zig");
const rstdbuild = rstd.buildutils;
const Build = rstdbuild.Build;
const OptimizeMode = rstdbuild.OptimizeMode;

const Machine = struct {
    /// Name for step and executable.
    name: []const u8,
    query: std.Target.Query,
    /// `qemu-system-<arch>` binary.
    qemu: []const u8,
    /// Machine and CPU selection plus whatever boots a bare `-kernel` ELF.
    args: []const []const u8,
    linker_script: []const u8,
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
        .linker_script = "src/qemutest/aarch64.ld",
    },
    .{
        .name = "arm",
        .query = .{ .cpu_arch = .arm, .os_tag = .freestanding, .abi = .eabihf, .cpu_model = .{ .explicit = &std.Target.arm.cpu.cortex_a15 } },
        .qemu = "qemu-system-arm",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a15" },
        .linker_script = "src/qemutest/arm.ld",
    },
    .{
        .name = "riscv64",
        .query = .{ .cpu_arch = .riscv64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-riscv64",
        .args = &.{ "-M", "virt", "-cpu", "rv64", "-smp", "1", "-bios", "none" },
        .linker_script = "src/qemutest/riscv64.ld",
        // The kernel sits at 0x80000000, out of reach of the default `medlow`.
        .code_model = .medany,
    },
    .{
        .name = "loongarch64",
        .query = .{ .cpu_arch = .loongarch64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-loongarch64",
        .args = &.{ "-M", "virt", "-cpu", "la464" },
        .linker_script = "src/qemutest/loongarch64.ld",
        .memory = "256M",
        .semihosting = false,
    },
    .{
        .name = "mips64el",
        .query = .{ .cpu_arch = .mips64el, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips64r6 } },
        .qemu = "qemu-system-mips64el",
        .args = &.{ "-M", "malta", "-cpu", "I6400" },
        .linker_script = "src/qemutest/mips64.ld",
        .memory = "256M",
    },
    .{
        .name = "powerpc64le",
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr9 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power9" },
        .linker_script = "src/qemutest/ppc64le.ld",
        .memory = "512M",
        .semihosting = false,
        // the firmware loader's 0x400000 shift, see the linker script
        .stderr_note = "kernel_addr changed",
    },
    .{
        .name = "s390x",
        .query = .{ .cpu_arch = .s390x, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.s390x.cpu.z14 } },
        .qemu = "qemu-system-s390x",
        .args = &.{ "-M", "s390-ccw-virtio", "-cpu", "max" },
        .linker_script = "src/qemutest/s390x.ld",
        .memory = "256M",
        .semihosting = false,
    },
};

/// Tests that need an OS (allocator, files, threads) live in files the test
/// root also reaches, so the runnable ones are picked by name instead.
const default_filters: []const []const u8 = &.{ "linalg.blas.blas1.axpy", "simd." };

/// Adds the `qemu-test` step and its per-arch `qemu-test-<arch>` steps.
pub fn addQemuTests(b: *Build, optimize: OptimizeMode) void {
    const all = b.step("qemu-test", "Run rstd's hardware-independent tests under qemu-system-* (all archs)");
    const filters = b.option([]const []const u8, "qemu-test-filter", "Only run tests whose name contains this (repeatable)") orelse &.{};

    for (machines) |m| {
        const target = b.resolveTargetQuery(m.query);
        const test_exe = b.addTest(.{
            .name = b.fmt("rstd-qemu-{s}", .{m.name}),
            .root_module = b.createModule(.{
                .root_source_file = b.path("src/test_qemu.zig"),
                .target = target,
                .optimize = optimize,
                .code_model = m.code_model,
            }),
            .test_runner = .{ .path = b.path("src/qemutest_runner.zig"), .mode = .simple },
            .filters = if (filters.len != 0) filters else default_filters,
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
