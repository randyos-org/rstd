//! The bench steps: `tools/blasbench/blas_bench.zig` on the Mac over SSH (`ssh-bench-aarch64-macos`), here (`bench-native`),
//! in a Linux VM on the host's CPU (`qemu-bench-linux-<arch>`), staged for the flyoffs, and compiled for every host (`bench-check`).

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const buildroot = @import("__root__.zig");
const qemutest = @import("qemu/__root__.zig");
const machines = @import("qemu/machines.zig");
const tools = @import("tools.zig");

const Machine = machines.Machine;
const TestHost = qemutest.TestHost;

/// `ssh-bench-aarch64-macos` on `mac` (`target`, through `tool`), `bench-native` and `bench-check`.
pub fn addBenches(b: *Build, mac: []const u8, target: rstdbuild.ResolvedTarget, tool: *rstdbuild.CompileStep) !void {
    const libs: BenchLibs = .{
        .openblas = b.option([]const u8, "openblas-lib", "OpenBLAS library for bench-native and ssh-bench-aarch64-macos, as a path where the bench runs (the Mac's default: nixpkgs' through Nix)"),
        .blis = b.option([]const u8, "blis-lib", "BLIS library for bench-native and ssh-bench-aarch64-macos, as a path where the bench runs"),
        .reference = b.option([]const u8, "reference-lib", "Reference-LAPACK's liblapack for bench-native, as a path where the bench runs"),
        .eigen = b.option([]const u8, "eigen-lib", "The bench's Eigen shim (tools/blasbench/nix) for bench-native, as a path where the bench runs"),
        .aocl_libflame = b.option([]const u8, "aocl-libflame-lib", "AOCL-libFLAME for bench-native, as a path where the bench runs"),
    };
    try addSshBench(b, mac, target, tool, libs);
    try addNativeBench(b, libs);
    addBenchCheck(b);
}

/// The libraries the bench's `lib*` and `lapack` sections race, by their
/// paths where it runs.
const BenchLibs = struct {
    openblas: ?[]const u8,
    blis: ?[]const u8,
    reference: ?[]const u8 = null,
    eigen: ?[]const u8 = null,
    aocl_libflame: ?[]const u8 = null,

    /// `--<option>=<path>` for each library given but OpenBLAS.
    fn addArgs(libs: BenchLibs, b: *Build, run: *Build.Step.Run) void {
        inline for (.{ .{ "blis", "blis" }, .{ "reference", "reference" }, .{ "eigen", "eigen" }, .{ "aocl_libflame", "aocl-libflame" } }) |f|
            if (@field(libs, f[0])) |p| run.addArg(b.fmt("--" ++ f[1] ++ "={s}", .{p}));
    }
};

/// `bench-check`: the bench compiled for each host it races on, run nowhere.
fn addBenchCheck(b: *Build) void {
    const step = b.step("bench-check", "Compile tools/blasbench/blas_bench.zig for Windows x86_64, Linux and macOS on x86_64 and arm64");
    const queries = [_]std.Target.Query{
        .{ .cpu_arch = .x86_64, .os_tag = .windows, .abi = .gnu },
        .{ .cpu_arch = .x86_64, .os_tag = .linux, .abi = .gnu },
        .{ .cpu_arch = .aarch64, .os_tag = .linux, .abi = .gnu, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.apple_m1 } },
        .{ .cpu_arch = .aarch64, .os_tag = .macos },
        .{ .cpu_arch = .x86_64, .os_tag = .macos },
    };
    for (queries) |q| {
        const exe = benchExe(b, b.resolveTargetQuery(q), null) catch |e| {
            step.dependOn(&b.addFail(b.fmt("bench-check: {t}", .{e})).step);
            return;
        };
        step.dependOn(&exe.step);
    }
    const merge_test = b.addTest(.{ .root_module = b.createModule(.{
        .root_source_file = b.path("tools/blasbench/merge.zig"),
        .target = b.graph.host,
    }) });
    step.dependOn(&b.addRunArtifact(merge_test).step);
}

/// `tools/blasbench/blas_bench.zig` for `target`, always ReleaseFast, with
/// OpenBLAS's vortex gemm `.S` assembled from the `openblas` checkout when
/// given (arm64 only).
fn benchExe(b: *Build, target: rstdbuild.ResolvedTarget, openblas: ?[]const u8) !*rstdbuild.CompileStep {
    const optimize: Optimize = .fast;
    const rstd = try buildroot.addRstdModule(b, null, target, optimize, .all, .{});
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

/// Where `tools/nix/blas.sh` links the Mac's libraries.
const mac_libs = "/tmp/rstd-blas/aarch64-darwin/";

/// `ssh-bench-aarch64-macos`: `tools/blasbench/blas_bench.zig` on the Mac,
/// with OpenBLAS's vortex gemm `.S` from the `-Dopenblas` checkout when
/// given. The libraries are `BenchLibs`, OpenBLAS by default nixpkgs'
/// (`tools/nix/blas.sh`), and Accelerate (`-Daccelerate`). Arguments after
/// `--` go to the bench.
fn addSshBench(b: *Build, mac: []const u8, target: rstdbuild.ResolvedTarget, tool: *rstdbuild.CompileStep, libs: BenchLibs) !void {
    const step = b.step("ssh-bench-aarch64-macos", b.fmt("Race rstd's BLAS against the libraries on the Mac at SSH alias {s}, and its kernels with -Dopenblas", .{mac}));
    const openblas = b.option([]const u8, "openblas", "OpenBLAS checkout whose vortex gemm assembly ssh-bench-aarch64-macos races");
    const exe = try benchExe(b, target, openblas);
    const check = b.step("ssh-bench-aarch64-macos-build", "Build what ssh-bench-aarch64-macos runs, without running it");
    check.dependOn(&exe.step);
    check.dependOn(&tool.step);

    const run = qemutest.watchCommand(b);
    run.addArg("--");
    run.addArtifactArg2(tool, .{});
    run.addArg(mac);
    run.addArtifactArg2(exe, .{});
    if (libs.openblas) |p| run.addArg(b.fmt("--openblas={s}", .{p})) else {
        const fetch = try tools.nixBlas(b, mac, "aarch64-darwin", true);
        fetch.has_side_effects = true;
        fetch.stdio = .inherit;
        run.step.dependOn(&fetch.step);
        run.addArg("--openblas=" ++ mac_libs ++ "openblasCompat/lib/libopenblas.dylib");
        run.addArg("--reference=" ++ mac_libs ++ "lapack-reference/lib/liblapack.dylib");
        run.addArg("--eigen=" ++ mac_libs ++ "eigen/lib/libeigen_shim.dylib");
    }
    libs.addArgs(b, run);
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
    const exe = try benchExe(b, b.graph.host, null);
    b.step("bench-native-build", "Build what bench-native runs, without running it").dependOn(&exe.step);
    const openblas = libs.openblas orelse {
        step.dependOn(&b.addFail("bench-native needs -Dopenblas-lib=<OpenBLAS library>").step);
        return;
    };
    const run = b.addRunArtifact(exe);
    run.addArg(b.fmt("--openblas={s}", .{openblas}));
    libs.addArgs(b, run);
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
pub fn addLinuxBenches(b: *Build, host: ?TestHost, whpx_host: bool) !void {
    // on the kernel command line, which the build fixes before `--` is known
    const bench_args = b.option([]const []const u8, "vm-bench-arg", "An argument for the bench in the qemu-bench-linux steps (repeatable), as `--` gives the other bench steps") orelse &.{};
    for (bench_machines) |bm| {
        const m = bm.m;
        if (m.whpx and !whpx_host) continue;
        const l = m.boot.linux;
        const step = b.step(b.fmt("qemu-{s}", .{m.name}), b.fmt("Race rstd's BLAS against nixpkgs' libraries in a Linux VM on the host's CPU (live log: {s})", .{qemutest.liveLog(b, m)}));
        const kernel = b.dependencyLazy(l.dep, .{}) catch continue;
        var query = m.query;
        query.cpu_model = if (bm.cpu) |c| .{ .explicit = c } else .native;
        const exe = try benchExe(b, b.resolveTargetQuery(query), null);
        const init_exe = b.addExecutable(.{ .name = b.fmt("init-{s}", .{m.name}), .root_module = b.createModule(.{
            .root_source_file = b.path("test/linux/init.zig"),
            .target = b.resolveTargetQuery(.{ .cpu_arch = m.query.cpu_arch, .os_tag = .linux }),
            .optimize = .fast,
        }) });
        const pack: qemutest.Initramfs = .init(b, "initramfs.cpio", init_exe, exe);
        // after `--`, init's arguments, which it passes on
        var append: std.ArrayList(u8) = .empty;
        try append.print(b.allocator, "{s} -- --openblas=/flyoff/libopenblas.so", .{l.commandLine(b)});
        if (bm.blis) try append.appendSlice(b.allocator, " --blis=/flyoff/libblis.so");
        // the reference binaries (`tools/nix/blas.sh`); one it skipped is skipped here
        try append.appendSlice(b.allocator, " --reference=/flyoff/liblapack.so --eigen=/flyoff/libeigen.so");
        if (bm.aocl) try append.appendSlice(b.allocator, " --aocl-libflame=/flyoff/libflame.so");
        const vm_append = try append.clone(b.allocator);

        // the Mac's flyoff packs the libraries there: the initramfs without them
        const host_mac = m.hvf != null;
        if (host_mac) {
            const base: qemutest.Initramfs = .init(b, "base.cpio", init_exe, exe);
            try addFlyoffStage(b, m, kernel.path(l.kernel), base.cpio, "base.cpio", vm_append.items);
        }

        const nix = host orelse {
            step.dependOn(&b.addFail(b.fmt("qemu-{s} gets its libraries from Nix on -Dtest-host", .{m.name})).step);
            if (!host_mac) b.step("flyoff-stage-zen3", "Stage the Zen 3 flyoff's files in zig-out/flyoff-zen3 (tools/blasbench/flyoff.sh)").dependOn(&b.addFail("flyoff-stage-zen3 gets the VM's libraries from Nix on -Dtest-host").step);
            continue;
        };
        // the libraries' closure, fetched once per flake.lock
        const libs = (try tools.nixBlas(b, nix.alias, bm.nix, true)).captureStdOut(.{});
        pack.run.addPrefixedFileArg("+", libs);
        if (!host_mac) try addFlyoffStage(b, m, kernel.path(l.kernel), pack.cpio, "initramfs.cpio", vm_append.items);

        const run = qemutest.qemuCommand(b, m, host);
        qemutest.addConsole(b, run, m, host);
        run.addArg("-kernel");
        qemutest.fileArg(run, m, host, kernel.path(l.kernel));
        run.addArg("-initrd");
        qemutest.fileArg(run, m, host, pack.cpio);
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
    const exe = try benchExe(b, target, null);
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
