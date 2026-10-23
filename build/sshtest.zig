//! `ssh-test-<target>` steps: rstd's full test suite, cross-compiled here
//! and run on a real machine over SSH (`tools/ssh_test.zig`), for targets
//! no emulator covers. The machine is an SSH alias, so `ssh <alias>` must
//! work without a password from the shell the build runs in; its address
//! and login live in that SSH config, not here.

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const buildroot = @import("__root__.zig");
const qemutest = @import("qemutest.zig");

/// Adds the SSH steps; `filters` is `-Drstd-test-filter`.
pub fn addSshTests(b: *Build, optimize: Optimize, filters: []const []const u8) !void {
    // Apple silicon: the cores only real hardware runs (vortex, vortexm4, SME)
    const mac = b.option([]const u8, "ssh-mac", "SSH alias of the Apple silicon Mac for ssh-test-aarch64-macos (default macmini-m6)") orelse "macmini-m6";
    const target = b.resolveTargetQuery(.{ .cpu_arch = .aarch64, .os_tag = .macos });
    const mod = try buildroot.addRstdModule(b, null, target, optimize, .all, false, false);
    const test_exe = try buildroot.addRstdTest(b, mod, target, filters);

    const tool = b.addExecutable(.{ .name = "ssh_test", .root_module = b.createModule(.{
        .root_source_file = b.path("tools/ssh_test.zig"),
        .target = b.graph.host,
    }) });
    // under the qemu rigs' watchdog: an hour with no output ends it
    const run = b.addRunArtifact(qemutest.watchTool(b));
    run.addArgs(&.{ "60", "--" });
    run.addArtifactArg(tool);
    run.addArg(mac);
    run.addArtifactArg(test_exe);
    run.stdio = .inherit;
    run.has_side_effects = true;
    b.step("ssh-test-aarch64-macos", b.fmt("Run rstd's tests on the Mac at SSH alias {s} (-Dssh-mac)", .{mac})).dependOn(&run.step);

    try addSshBench(b, mac, target, tool);
}

/// `ssh-bench-aarch64-macos`: `tools/blasbench/blas_bench.zig`, always
/// ReleaseFast, with OpenBLAS's vortex gemm `.S` assembled from the
/// `-Dopenblas` checkout; arguments after `--` go to the bench.
fn addSshBench(b: *Build, mac: []const u8, target: rstdbuild.ResolvedTarget, tool: *rstdbuild.CompileStep) !void {
    const step = b.step("ssh-bench-aarch64-macos", b.fmt("Race rstd's BLAS kernels on the Mac at SSH alias {s} (needs -Dopenblas)", .{mac}));
    const openblas = b.option([]const u8, "openblas", "OpenBLAS checkout whose gemm assembly ssh-bench-aarch64-macos races") orelse {
        step.dependOn(&b.addFail("ssh-bench-aarch64-macos needs -Dopenblas=<OpenBLAS checkout>").step);
        return;
    };
    const optimize: Optimize = .fast;
    const rstd = try buildroot.addRstdModule(b, null, target, optimize, .all, false, false);
    const mod = b.createModule(.{
        .root_source_file = b.path("tools/blasbench/blas_bench.zig"),
        .target = target,
        .optimize = optimize,
        .link_libc = true,
    });
    mod.addImport("rstd", rstd);

    // `config.h` as `tools/blas_asm` writes it
    const inc = b.addWriteFiles();
    _ = inc.add("config.h", "#define ARCH_ARM64 1\n#define BINARY64 1\n");
    const kernel_dir = b.pathJoin(&.{ openblas, "kernel", "arm64" });
    mod.addIncludePath(inc.getDirectory());
    mod.addIncludePath(.{ .cwd_relative = openblas });
    mod.addIncludePath(.{ .cwd_relative = kernel_dir });
    const sources = [_][3][]const u8{
        .{ "sgemm_kernel_16x4.S", "sgemm_kernel", "" },
        .{ "sgemm_tcopy_16.S", "sgemm_itcopy", "" },
        .{ "sgemm_ncopy_4.S", "sgemm_oncopy", "" },
        .{ "dgemm_kernel_8x4.S", "dgemm_kernel", "-DDOUBLE" },
        .{ "dgemm_tcopy_8.S", "dgemm_itcopy", "-DDOUBLE" },
        .{ "dgemm_ncopy_4.S", "dgemm_oncopy", "-DDOUBLE" },
    };
    for (sources) |s| {
        const name = b.fmt("_ob_{s}", .{s[1]});
        var flags: std.ArrayList([]const u8) = .empty;
        try flags.appendSlice(b.allocator, &vortex_defines);
        try flags.appendSlice(b.allocator, &.{
            b.fmt("-DASMNAME={s}", .{name}),
            b.fmt("-DASMFNAME={s}_", .{name}),
            b.fmt("-DNAME={s}_", .{name}),
            b.fmt("-DCNAME={s}", .{name}),
        });
        if (s[2].len > 0) try flags.append(b.allocator, s[2]);
        mod.addCSourceFile(.{
            .file = .{ .cwd_relative = b.pathJoin(&.{ kernel_dir, s[0] }) },
            .flags = flags.items,
            .language = .assembly_with_preprocessor,
        });
    }
    const exe = b.addExecutable(.{ .name = "blas_bench", .root_module = mod });

    const run = b.addRunArtifact(qemutest.watchTool(b));
    run.addArgs(&.{ "60", "--" });
    run.addArtifactArg(tool);
    run.addArg(mac);
    run.addArtifactArg(exe);
    run.addPassthruArgs();
    run.stdio = .inherit;
    run.has_side_effects = true;
    step.dependOn(&run.step);
}

/// The preprocessor flags of OpenBLAS's build for VORTEX: `getarch.c`'s
/// `ARCHCONFIG`, then what its build passes every arm64 `.S`.
const vortex_defines = [_][]const u8{
    "-w",
    "-DASSEMBLER",
    "-DARCH_ARM64",
    "-DOS_DARWIN",
    "-DDYNAMIC_ARCH",
    "-DVORTEX",
    "-DL1_DATA_SIZE=32768",
    "-DL1_DATA_LINESIZE=64",
    "-DL2_SIZE=262144",
    "-DL2_LINESIZE=64",
    "-DDTB_DEFAULT_ENTRIES=64",
    "-DDTB_SIZE=4096",
    "-DL2_ASSOCIATIVE=32",
    "-DHAVE_VFPV4",
    "-DHAVE_VFPV3",
    "-DHAVE_VFP",
    "-DHAVE_NEON",
    "-DARMV8",
};
