//! `ssh-test-<target>` steps: rstd's full test suite, cross-compiled here
//! and run on a real machine over SSH (`tools/ssh_test.zig`), for targets
//! no emulator covers. The machine is an SSH alias, so `ssh <alias>` must
//! work without a password from the shell the build runs in; its address
//! and login live in that SSH config, not here.
//!
//! Also the BLAS bench (`tools/blasbench/blas_bench.zig`) on real CPUs:
//! `ssh-bench-aarch64-macos` on the Mac, `bench-native` here; the Linux VMs'
//! `qemu-bench-linux-<arch>` are in `qemutest.zig`.

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const buildroot = @import("__root__.zig");
const qemutest = @import("qemutest.zig");

/// Adds the SSH steps and `bench-native`; `filters` is `-Drstd-test-filter`.
pub fn addSshTests(b: *Build, optimize: Optimize, filters: []const []const u8) !void {
    const mac = b.option([]const u8, "ssh-mac", "SSH alias of the Apple silicon Mac for the ssh-test steps (default macmini-m6)") orelse "macmini-m6";
    const tool = sshTool(b);
    // Apple silicon: the cores only real hardware runs (vortex, vortexm4, SME)
    const arm = b.resolveTargetQuery(.{ .cpu_arch = .aarch64, .os_tag = .macos });
    try addSshTest(b, optimize, filters, mac, tool, arm, "on the Mac");
    // x86_64 macOS under Rosetta 2, which the Mac must have installed
    const x86 = b.resolveTargetQuery(.{ .cpu_arch = .x86_64, .os_tag = .macos });
    try addSshTest(b, optimize, filters, mac, tool, x86, "on the Mac under Rosetta 2");

    const libs: BenchLibs = .{
        .openblas = b.option([]const u8, "openblas-lib", "OpenBLAS library for bench-native and ssh-bench-aarch64-macos, as a path where the bench runs (the Mac's default: nixpkgs' through Nix)"),
        .blis = b.option([]const u8, "blis-lib", "BLIS library for bench-native and ssh-bench-aarch64-macos, as a path where the bench runs"),
    };
    try addSshBench(b, mac, arm, tool, libs);
    try addNativeBench(b, libs);
}

var ssh_tool: ?*rstdbuild.CompileStep = null;

/// `tools/ssh_test.zig`, built once per build.
pub fn sshTool(b: *Build) *rstdbuild.CompileStep {
    if (ssh_tool) |t| if (t.step.owner == b) return t;
    ssh_tool = b.addExecutable(.{ .name = "ssh_test", .root_module = b.createModule(.{
        .root_source_file = b.path("tools/ssh_test.zig"),
        .target = b.graph.host,
    }) });
    return ssh_tool.?;
}

/// `ssh-test-<arch>-macos`: the full suite for `target`, run on `mac`.
fn addSshTest(b: *Build, optimize: Optimize, filters: []const []const u8, mac: []const u8, tool: *rstdbuild.CompileStep, target: rstdbuild.ResolvedTarget, where: []const u8) !void {
    const mod = try buildroot.addRstdModule(b, null, target, optimize, .all, false, false);
    const test_exe = try buildroot.addRstdTest(b, mod, target, optimize, filters);
    // under the qemu rigs' watchdog: an hour with no output ends it
    const run = qemutest.watchCommand(b);
    run.addArg("--");
    run.addArtifactArg2(tool, .{});
    run.addArg(mac);
    run.addArtifactArg2(test_exe, .{});
    run.stdio = .inherit;
    run.has_side_effects = true;
    const name = b.fmt("ssh-test-{t}-macos", .{target.result.cpu.arch});
    b.step(name, b.fmt("Run rstd's {t} tests {s} at SSH alias {s} (-Dssh-mac)", .{ target.result.cpu.arch, where, mac })).dependOn(&run.step);
}

/// The libraries the bench's `lib*` sections race, by their paths where it
/// runs.
const BenchLibs = struct { openblas: ?[]const u8, blis: ?[]const u8 };

/// `tools/blasbench/blas_bench.zig` for `target`, always ReleaseFast, with
/// OpenBLAS's vortex gemm `.S` assembled from the `openblas` checkout when
/// given (arm64 only).
pub fn benchExe(b: *Build, target: rstdbuild.ResolvedTarget, openblas: ?[]const u8) !*rstdbuild.CompileStep {
    const optimize: Optimize = .fast;
    const rstd = try buildroot.addRstdModule(b, null, target, optimize, .all, false, false);
    const mod = b.createModule(.{
        .root_source_file = b.path("tools/blasbench/blas_bench.zig"),
        .target = target,
        .optimize = optimize,
        .link_libc = true,
    });
    mod.addImport("rstd", rstd);
    const options = b.addOptions();
    options.addOption(bool, "ob_asm", openblas != null);
    mod.addOptions("bench_options", options);
    if (openblas) |checkout| try addVortexAsm(b, mod, checkout);
    return b.addExecutable(.{ .name = "blas_bench", .root_module = mod });
}

/// `tools/nix/blas.sh <system>` on `host` (`ssh_test`), at the nixpkgs of
/// `tools/nix/flake.lock`: nixpkgs' OpenBLAS and BLIS for `system`.
pub fn nixBlas(b: *Build, host: []const u8, system: []const u8) !*Build.Step.Run {
    const lock = try std.Io.Dir.cwd().readFileAlloc(b.graph.io, try b.root.joinString(b.allocator, "tools/nix/flake.lock"), b.allocator, .limited(1 << 16));
    const Lock = struct { nodes: struct { nixpkgs: struct { locked: struct { rev: []const u8 } } } };
    const parsed = try std.json.parseFromSliceLeaky(Lock, b.allocator, lock, .{ .ignore_unknown_fields = true });
    const run = b.addRunArtifact(sshTool(b));
    run.addArg(host);
    run.addFileArg(b.path("tools/nix/blas.sh"));
    run.addArgs(&.{ system, parsed.nodes.nixpkgs.locked.rev });
    return run;
}

/// `ssh-bench-aarch64-macos`: `tools/blasbench/blas_bench.zig` on the Mac,
/// with OpenBLAS's vortex gemm `.S` from the `-Dopenblas` checkout when
/// given. The libraries are `BenchLibs`, OpenBLAS by default nixpkgs'
/// (`tools/nix/blas.sh`), and Accelerate (`-Daccelerate`). Arguments after
/// `--` go to the bench.
fn addSshBench(b: *Build, mac: []const u8, target: rstdbuild.ResolvedTarget, tool: *rstdbuild.CompileStep, libs: BenchLibs) !void {
    const step = b.step("ssh-bench-aarch64-macos", b.fmt("Race rstd's BLAS against the libraries on the Mac at SSH alias {s}, and its kernels with -Dopenblas", .{mac}));
    const openblas = b.option([]const u8, "openblas", "OpenBLAS checkout whose vortex gemm assembly ssh-bench-aarch64-macos races");
    const exe = try benchExe(b, target, openblas);

    const run = qemutest.watchCommand(b);
    run.addArg("--");
    run.addArtifactArg2(tool, .{});
    run.addArg(mac);
    run.addArtifactArg2(exe, .{});
    if (libs.openblas) |p| run.addArg(b.fmt("--openblas={s}", .{p})) else {
        const fetch = try nixBlas(b, mac, "aarch64-darwin");
        fetch.has_side_effects = true;
        fetch.stdio = .inherit;
        run.step.dependOn(&fetch.step);
        run.addArg("--openblas=/tmp/rstd-blas/aarch64-darwin/openblasCompat/lib/libopenblas.dylib");
    }
    if (libs.blis) |p| run.addArg(b.fmt("--blis={s}", .{p}));
    if (b.option([]const u8, "accelerate", "Accelerate on the Mac for ssh-bench-aarch64-macos (default the system framework)")) |p|
        run.addArg(b.fmt("--accelerate={s}", .{p}));
    run.addPassthruArgs();
    run.stdio = .inherit;
    run.has_side_effects = true;
    step.dependOn(&run.step);
}

/// `bench-native`: `tools/blasbench/blas_bench.zig` here, on this machine's
/// own CPU, against `-Dopenblas-lib` and `-Dblis-lib`. Arguments after `--`
/// go to the bench.
fn addNativeBench(b: *Build, libs: BenchLibs) !void {
    const step = b.step("bench-native", "Race rstd's BLAS against -Dopenblas-lib and -Dblis-lib on this machine");
    const openblas = libs.openblas orelse {
        step.dependOn(&b.addFail("bench-native needs -Dopenblas-lib=<OpenBLAS library>").step);
        return;
    };
    const run = b.addRunArtifact(try benchExe(b, b.graph.host, null));
    run.addArg(b.fmt("--openblas={s}", .{openblas}));
    if (libs.blis) |p| run.addArg(b.fmt("--blis={s}", .{p}));
    run.addPassthruArgs();
    run.stdio = .inherit;
    run.has_side_effects = true;
    step.dependOn(&run.step);
}

/// OpenBLAS's vortex gemm kernels and copies from `checkout`, as
/// `tools/blas_asm` has them, named `ob_<routine>`.
fn addVortexAsm(b: *Build, mod: *rstdbuild.Module, checkout: []const u8) !void {
    // `config.h` as `tools/blas_asm` writes it
    const inc = b.addWriteFiles();
    _ = inc.add("config.h", "#define ARCH_ARM64 1\n#define BINARY64 1\n");
    const kernel_dir = b.pathJoin(&.{ checkout, "kernel", "arm64" });
    mod.addIncludePath(inc.getDirectory());
    mod.addIncludePath(.{ .cwd_relative = checkout });
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
