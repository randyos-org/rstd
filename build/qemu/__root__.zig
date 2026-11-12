//! `qemu-test` steps: rstd's tests on `qemu-system-*` machines (`machines.zig`): a freestanding kernel (`test/qemu/`, root
//! `src/_test_qemu.zig`; on UEFI an application), Linux with an initramfs (`test/linux/init.zig`), or an installed macOS over SSH
//! (`macos.zig`). Consoles log live to `.zig-cache/qemu-test-<name>.log`; without semihosting the summary carries the result.

const std = @import("std");
const rstdbuild = @import("../../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const blaskern = @import("../blaskern.zig");
const buildroot = @import("../__root__.zig");
const tools = @import("../tools.zig");
const machines = @import("machines.zig");
const macos = @import("macos.zig");

const Machine = machines.Machine;
const Boot = machines.Boot;

/// What the VM benches (`bench.zig`) take from the rigs' options.
pub const Rigs = struct { test_host: ?TestHost, whpx_host: bool, linux_vm: bool };

/// Adds the `qemu-test` step and each machine's `qemu-test-<name>` step.
pub fn addQemuTests(b: *Build, optimize: Optimize) !Rigs {
    const all = b.step("qemu-test", "Run rstd's hardware-independent tests under qemu-system-* (all archs)");
    const check = b.step("qemu-check", "Build every kernel rig's test images and check their inline asm, without running qemu");
    const filters = b.option([]const []const u8, "qemu-test-filter", "Only run tests whose name contains this (repeatable)") orelse &.{};
    const only_parts = b.option([]const []const u8, "qemu-test-part", "Only build and run this part of a rig that splits its tests into images, e.g. blas3-kernel-complex (repeatable)") orelse &.{};
    for (only_parts) |o| if (!hasPart(o)) std.debug.panic("-Dqemu-test-part: no rig has a part named {s}", .{o});
    // On when rstd is the project being built, so a dependent never fetches OSX-KVM.
    const macos_vm = vmOption(b, "macos-vm", "the macOS VM rigs", "OSX-KVM");
    const mac: macos.Options = if (macos_vm) .{
        .user = b.option([]const u8, "mac-user", "macOS account in the VM (default randyos)") orelse "randyos",
        .port = b.option(u16, "mac-port", "Host port forwarded to the VM's ssh") orelse 2222,
        .release = b.option(macos.Release, "mac-release", "macOS release qemu-setup fetches (default sonoma)") orelse .sonoma,
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
    const host_named = b.option([]const u8, "test-host", "SSH alias of a machine to run the rigs' QEMU on, Nix's from tools/nix (default " ++ tools.mac_alias ++ " when rstd is the root project); images still build here, and x86_64-whpx and the macOS VM rigs stay local");
    const host_alias: ?[]const u8 = if (host_off) null else host_named orelse if (b.isRoot()) tools.mac_alias else null;
    const host_jobs = b.option(u32, "test-host-jobs", "Rigs running at once on -Dtest-host, across builds (default 8)") orelse 8;
    const host_tcg = b.option(bool, "test-host-tcg", "Run every rig on -Dtest-host under TCG, also those HVF could run") orelse false;
    const test_host: ?TestHost = if (host_alias) |alias| .{ .alias = alias, .jobs = host_jobs, .tcg = host_tcg } else null;

    const host = b.graph.host.result;
    const whpx_host = host.os.tag == .windows and host.cpu.arch == .x86_64;
    for (machines.all) |m| {
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
            .macos => if (macos_vm) try macos.addRun(b, m, target, m.optimize orelse optimize, filters, mac) orelse continue else continue,
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
    return .{ .test_host = test_host, .whpx_host = whpx_host, .linux_vm = linux_vm };
}

fn hasPart(name: []const u8) bool {
    for (machines.all) |m| for (m.parts) |p| if (std.mem.eql(u8, p.name, name)) return true;
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
pub const TestHost = struct { alias: []const u8, jobs: u32, tcg: bool };

/// Where `m` runs: the test host, or here (null).
fn testHostFor(m: Machine, host: ?TestHost) ?TestHost {
    return if (m.whpx) null else host;
}

/// `m`'s qemu and `args`, under `tools/qemu_watch.zig`: a firmware halt or
/// an hour with no output ends it and fails the step, where it would
/// otherwise never end. The slowest test prints nothing for about 20 minutes
/// on the x86 rigs. On a test host, through `tools/test_host.zig`, which
/// copies over the files `fileArg` adds.
pub fn qemuCommand(b: *Build, m: Machine, host: ?TestHost) *Build.Step.Run {
    const run = watchCommand(b);
    // OpenBIOS's halt, and OVMF's and AAVMF's CPU exceptions
    run.addArgs(&.{ "Stopping execution", "!!!! X64 Exception", "Synchronous Exception at", "--" });
    const h = testHostFor(m, host) orelse {
        run.addArg(m.qemu);
        run.addArgs(m.args);
        return run;
    };
    run.addArtifactArg2(tools.hostTool(b, .test_host), .{});
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
pub fn fileArg(run: *Build.Step.Run, m: Machine, host: ?TestHost, file: Build.LazyPath) void {
    if (testHostFor(m, host) != null) run.addPrefixedFileArg("local:", file) else run.addFileArg(file);
}

/// `tools/qemu_watch.zig` with its idle limit (an hour) and `-Dfail-fast`;
/// the caller adds halt texts, `--` and the command.
pub fn watchCommand(b: *Build) *Build.Step.Run {
    const run = b.addRunArtifact(tools.hostTool(b, .qemu_watch));
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

/// Where `m`'s console is logged live, relative to the build root.
pub fn liveLog(b: *Build, m: Machine) []const u8 {
    return b.fmt(".zig-cache/qemu-test-{s}.log", .{m.name});
}

/// qemu's console on stdout and, as it is written, in `liveLog`: a run step
/// checking stdout shows it only once qemu exits. Naming the serial port and
/// monitor keeps `-nographic`'s firmware console. On a test host
/// `tools/test_host.zig` writes the log.
pub fn addConsole(b: *Build, run: *Build.Step.Run, m: Machine, host: ?TestHost) void {
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
    const mod = try buildroot.addRstdModule(b, null, t, optimize, .all, .{});
    const test_exe = try buildroot.addRstdTest(b, mod, t, optimize, run_filters);
    const init_exe = b.addExecutable(.{ .name = b.fmt("init-{s}", .{m.name}), .root_module = b.createModule(.{
        .root_source_file = b.path("test/linux/init.zig"),
        .target = target,
        .optimize = optimize,
    }) });
    const pack: Initramfs = .init(b, "initramfs.cpio", init_exe, test_exe);
    if (nix) |system| pack.run.addPrefixedFileArg("+", (try tools.nixBlas(b, host.?.alias, system, false)).captureStdOut(.{}));
    rig_check.dependOn(&pack.run.step);

    const run = qemuCommand(b, m, host);
    addConsole(b, run, m, host);
    run.addArg("-kernel");
    fileArg(run, m, host, kernel.path(l.kernel));
    run.addArg("-initrd");
    fileArg(run, m, host, pack.cpio);
    run.addArgs(&.{ "-append", l.commandLine(b) });
    run.has_side_effects = true;
    // `passed` in test/linux/init.zig
    run.expectStdOutMatch("linuxtest: /test exited 0");
    // no newline: the serial console sends `\r\n`
    if (core) |c| run.expectStdOutMatch(b.fmt("blas core: {s}", .{c}));
    return run;
}

/// An initramfs (`tools/initramfs.zig`) of `init` as `/init` and `test` as
/// `/test`; `run` takes more files.
pub const Initramfs = struct {
    run: *Build.Step.Run,
    cpio: Build.LazyPath,

    pub fn init(b: *Build, name: []const u8, init_exe: *Build.Step.Compile, test_exe: *Build.Step.Compile) Initramfs {
        const run = b.addRunArtifact(tools.hostTool(b, .initramfs));
        const cpio = run.addOutputFileArg2(name, .{});
        run.addArtifactArg2(init_exe, .{ .prefix = "init=" });
        run.addArtifactArg2(test_exe, .{ .prefix = "test=" });
        return .{ .run = run, .cpio = cpio };
    }
};

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
        .filters = if (filters.len != 0) filters else if (m.cpus > 1) machines.smp_filters else machines.default_filters,
        .use_llvm = if (k.llvm) true else null,
        // every rig needs its linker script, which only lld honors
        .use_lld = true,
    });
    if (k.firmware == null) test_exe.setLinkerScript(b.path(k.linker_script));
    check.dependOn(&test_exe.step);

    const run = qemuCommand(b, m, host);
    if (k.load_window) |w| {
        const fit = b.addRunArtifact(tools.hostTool(b, .load_window));
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
    const semihosting = m.semihosting();
    if (semihosting) run.addArgs(&.{ "-semihosting-config", "enable=on,target=native,chardev=con" });
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
    if (!semihosting or hvf) run.expectStdOutMatch(", 0 failed");
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
