//! `qemu-test` step: runs rstd's tests on emulated machines with
//! `qemu-system-*`. Needs `qemu-system-<arch>` on `PATH`. A machine boots one
//! of two ways (`Boot`):
//!
//!   * `.kernel`: a freestanding test kernel per architecture (boot code in
//!     `src/_qemutest/`, runner in `src/_test_qemu_runner.zig`, test root in
//!     `src/_test_qemu.zig`). Where the machine has semihosting (aarch64, arm,
//!     riscv64, mips64el) the runner exits with the number of failed tests, so a
//!     failing test fails the step. Elsewhere the machine just powers off and the
//!     step checks the printed summary for "0 failed".
//!   * `.macos`: rstd's unit tests in an installed macOS, which the Mach-O
//!     kernels need, copied in and run over SSH. It boots with OSX-KVM's
//!     OpenCore and firmware, a lazy dependency fetched only when rstd is the
//!     project being built (`-Dmacos-vm=false` skips it), so a dependent never
//!     fetches it. The VM lives in `.macos-vm/<name>/` (gitignored).
//!
//! A macOS machine, once:
//!
//!   1. `zig build qemu-setup-<name>`: fetches a recovery image
//!      (`-Dmac-release`, Sonoma by default) and makes the disk.
//!   2. `zig build qemu-install-<name>`, then in the window: "macOS Base
//!      System" in the OpenCore picker; Disk Utility, erase the QEMU HARDDISK
//!      as APFS; reinstall macOS onto it (pick "macOS Installer", then the
//!      disk, at each reboot); create an account; System Settings > General >
//!      Sharing, turn on Remote Login; at the next picker, select the disk and
//!      press Ctrl+Enter to make it the default.
//!   3. Authorize an SSH key (`ssh-keygen` makes one if you have none):
//!      `type %USERPROFILE%\.ssh\id_ed25519.pub | ssh -p 2222 <account>@localhost
//!      "mkdir -p ~/.ssh && cat >> ~/.ssh/authorized_keys"`.
//!
//! Then `zig build qemu-start-<name>` boots it in a window, and, while it is
//! up, `zig build qemu-test-<name> -Dmac-user=<account>` runs the tests in it.

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
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
};

/// How the tests reach a machine.
const Boot = union(enum) {
    kernel: Kernel,
    /// An installed macOS, set up once (see the file header).
    macos,

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
};

const machines = [_]Machine{
    .{
        .name = "aarch64",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .freestanding, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.cortex_a72 } },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a72" },
        .boot = .{ .kernel = .{
            .linker_script = "src/_qemutest/aarch64.ld",
        } },
    },
    .{
        .name = "aarch64-sve",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .freestanding, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.neoverse_v1 } },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "neoverse-v1" },
        .boot = .{ .kernel = .{
            .linker_script = "src/_qemutest/aarch64.ld",
        } },
    },
    .{
        .name = "arm",
        .query = .{ .cpu_arch = .arm, .os_tag = .freestanding, .abi = .eabihf, .cpu_model = .{ .explicit = &std.Target.arm.cpu.cortex_a15 } },
        .qemu = "qemu-system-arm",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a15" },
        .boot = .{ .kernel = .{
            .linker_script = "src/_qemutest/arm.ld",
        } },
    },
    .{
        .name = "riscv64",
        .query = .{ .cpu_arch = .riscv64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-riscv64",
        // the vector extension, which the RVV kernels need
        .args = &.{ "-M", "virt", "-cpu", "rv64,v=true,vlen=256", "-smp", "1", "-bios", "none" },
        .boot = .{
            .kernel = .{
                .linker_script = "src/_qemutest/riscv64.ld",
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
            .linker_script = "src/_qemutest/loongarch64.ld",
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
            .linker_script = "src/_qemutest/mips64.ld",
        } },
    },
    .{
        .name = "mips64r6el",
        .query = .{ .cpu_arch = .mips64el, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips64r6 } },
        .qemu = "qemu-system-mips64el",
        .args = &.{ "-M", "malta", "-cpu", "I6400" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "src/_qemutest/mips64.ld",
        } },
    },
    .{
        .name = "mips32r2el",
        .query = .{ .cpu_arch = .mipsel, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips32r2 } },
        .qemu = "qemu-system-mipsel",
        .args = &.{ "-M", "malta", "-cpu", "24Kf" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "src/_qemutest/mips32.ld",
        } },
    },
    .{
        .name = "mips32r5el",
        .query = .{ .cpu_arch = .mipsel, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.p5600 } },
        .qemu = "qemu-system-mipsel",
        .args = &.{ "-M", "malta", "-cpu", "P5600" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "src/_qemutest/mips32.ld",
        } },
    },
    .{
        .name = "i386",
        .query = .{ .cpu_arch = .x86, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.pentium3 } },
        .qemu = "qemu-system-i386",
        .args = &.{ "-M", "pc", "-cpu", "pentium3" },
        .memory = "128M",
        .boot = .{ .kernel = .{
            .linker_script = "src/_qemutest/x86.ld",
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
            .linker_script = "src/_qemutest/x86_64.ld",
            .llvm = true,
            .semihosting = false,
        } },
    },
    .{
        .name = "sparc64",
        .optimize = .fast,
        .query = .{ .cpu_arch = .sparc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.sparc.cpu.v9 } },
        .qemu = "qemu-system-sparc64",
        .args = &.{ "-M", "sun4u" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "src/_qemutest/sparc64.ld",
            .semihosting = false,
        } },
    },
    .{
        .name = "powerpc64le",
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr9 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power9" },
        .memory = "512M",
        .boot = .{
            .kernel = .{
                .linker_script = "src/_qemutest/ppc64le.ld",
                .semihosting = false,
                // the firmware loader's 0x400000 shift, see the linker script
                .stderr_note = "kernel_addr changed",
            },
        },
    },
    .{
        .name = "powerpc64le-power10",
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr10 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power10" },
        .memory = "512M",
        .boot = .{ .kernel = .{
            .linker_script = "src/_qemutest/ppc64le.ld",
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
            .linker_script = "src/_qemutest/s390x.ld",
            .semihosting = false,
        } },
    },
    .{
        // macOS 13 and up; Zig's default minimum is 15.0
        .name = "x86_64-macos",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .macos, .os_version_min = .{ .semver = .{ .major = 13, .minor = 0, .patch = 0 } } },
        .qemu = "qemu-system-x86_64",
        .args = &.{
            // TCG emulates (slow) but honors the Intel vendor below; under WHPX
            // ("whpx", as randyos's own VM uses) the kernel stopped at
            // EXITBS:START, as if it saw the host's AMD CPUID. Once the
            // installer hung here with every CPU spinning in the kernel,
            // interrupts off, possibly after the host slept; if it recurs
            // with the host awake, use `thread=single` (one host thread).
            // A 2 GB translation cache keeps macOS's code from being flushed
            // and translated again.
            "-accel",  "tcg,thread=multi,tb-size=2048",
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

/// Tests that need an OS (allocator, files, threads) live in files the test
/// root also reaches, so the runnable ones are picked by name instead.
const default_filters: []const []const u8 = &.{ "linalg.blas.blas1.", "linalg.blas.blas3.gemm.params", "simd." };

/// Adds the `qemu-test` step and each machine's `qemu-test-<name>` step.
pub fn addQemuTests(b: *Build, optimize: Optimize) !void {
    const all = b.step("qemu-test", "Run rstd's hardware-independent tests under qemu-system-* (all archs)");
    const filters = b.option([]const []const u8, "qemu-test-filter", "Only run tests whose name contains this (repeatable)") orelse &.{};
    // On when rstd is the project being built, so a dependent never fetches OSX-KVM.
    const macos_vm = b.option(bool, "macos-vm", "Add the macOS VM rigs, fetching OSX-KVM (default: when rstd is the root project)") orelse b.isRoot();
    const mac: MacOptions = if (macos_vm) .{
        .user = b.option([]const u8, "mac-user", "macOS account in the VM") orelse "",
        .port = b.option(u16, "mac-port", "Host port forwarded to the VM's ssh") orelse 2222,
        .release = b.option(MacRelease, "mac-release", "macOS release qemu-setup fetches (default sonoma)") orelse .sonoma,
        .gdb = b.option(bool, "mac-gdb", "Let gdb attach to the macOS VM on localhost:1234") orelse false,
    } else .{};

    for (machines) |m| {
        if (m.stub != null) continue;
        const target = b.resolveTargetQuery(m.query);
        const run = switch (m.boot) {
            .kernel => |k| addKernelRun(b, m, k, target, m.optimize orelse optimize, filters),
            .macos => if (macos_vm) try addMacosRun(b, m, target, m.optimize orelse optimize, mac) orelse continue else continue,
        };
        const step = b.step(b.fmt("qemu-test-{s}", .{m.name}), b.fmt("Run rstd's tests under qemu-system on {s}", .{m.name}));
        step.dependOn(&run.step);
        all.dependOn(&run.step);
    }
}

/// QEMU booting `m`'s freestanding test kernel.
fn addKernelRun(b: *Build, m: Machine, k: Boot.Kernel, target: rstdbuild.ResolvedTarget, optimize: Optimize, filters: []const []const u8) *Build.Step.Run {
    const test_mod = b.createModule(.{
        .root_source_file = b.path("src/_test_qemu.zig"),
        .target = target,
        .optimize = optimize,
        .code_model = k.code_model,
    });
    // the kernels of the architectures that have any
    blaskern.addBlasKernels(b, test_mod, target, .all);
    const test_exe = b.addTest(.{
        .name = b.fmt("rstd-qemu-{s}", .{m.name}),
        .root_module = test_mod,
        .test_runner = .{ .path = b.path("src/_test_qemu_runner.zig"), .mode = .simple },
        .filters = if (filters.len != 0) filters else default_filters,
        .use_llvm = if (k.llvm) true else null,
        // every rig needs its linker script, which only lld honors
        .use_lld = true,
    });
    test_exe.setLinkerScript(b.path(k.linker_script));

    const run = b.addSystemCommand(&.{m.qemu});
    run.addArgs(m.args);
    run.addArgs(&.{ "-m", m.memory, "-nographic", "-no-reboot" });
    if (k.semihosting) run.addArgs(&.{ "-semihosting-config", "enable=on,target=native" });
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
    };

    const setup = b.step(b.fmt("qemu-setup-{s}", .{m.name}), b.fmt("macOS VM {s}: fetch a recovery image, make the disk", .{m.name}));
    var fetch: ?*Build.Step.Run = null;
    if (!exists(b, paths.dmg)) {
        const run = console(b, &.{"python"});
        run.addFileArg(kvm.path("fetch-macOS-v2.py"));
        run.addArgs(&.{ "-s", @tagName(mac.release) });
        // it saves where it runs; the package is read-only
        run.setCwd(b.path(vm));
        setup.dependOn(&run.step);
        fetch = run;
    }
    if (!exists(b, paths.img)) {
        const run = console(b, &.{ "qemu-img", "convert", "-f", "dmg", "-O", "raw", paths.dmg, paths.img });
        if (fetch) |f| run.step.dependOn(&f.step);
        setup.dependOn(&run.step);
    }
    if (!exists(b, paths.disk)) setup.dependOn(&console(b, &.{ "qemu-img", "create", "-f", "qcow2", paths.disk, "96G" }).step);
    if (!exists(b, paths.vars)) {
        const copy = b.addUpdateSourceFiles();
        copy.addCopyFileToSource(kvm.path("OVMF_VARS-1920x1080.fd"), paths.vars);
        setup.dependOn(&copy.step);
    }

    inline for (.{ "install", "start" }) |cmd| {
        const run = try macosQemu(b, m, kvm, paths, mac, std.mem.eql(u8, cmd, "install"));
        b.step(b.fmt("qemu-" ++ cmd ++ "-{s}", .{m.name}), b.fmt("macOS VM {s}: " ++ cmd ++ ", in a window", .{m.name})).dependOn(&run.step);
    }

    const mod = try buildroot.addRstdModule(b, null, target, optimize, .all, false, false);
    const login = b.fmt("{s}@localhost", .{mac.user});
    const port = b.fmt("{d}", .{mac.port});
    const ssh_options: []const []const u8 = &.{ "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "-o", "StrictHostKeyChecking=accept-new" };
    const copy = b.addSystemCommand(&.{ "scp", "-q", "-P", port });
    copy.addArgs(ssh_options);
    copy.addArtifactArg(try buildroot.addRstdTest(b, mod, target));
    copy.addArg(b.fmt("{s}:/tmp/rstd-test", .{login}));
    copy.has_side_effects = true;
    if (mac.user.len == 0) copy.step.dependOn(&b.addFail("the macOS VM test needs -Dmac-user=<account>").step);
    const run = console(b, &.{ "ssh", "-p", port });
    run.addArgs(ssh_options);
    run.addArgs(&.{ login, "/tmp/rstd-test" });
    run.step.dependOn(&copy.step);
    return run;
}

/// The macOS VMs' SSH login (`-Dmac-user`, `-Dmac-port`) and release.
const MacOptions = struct { user: []const u8 = "", port: u16 = 2222, release: MacRelease = .sonoma, gdb: bool = false };

/// The releases `fetch-macOS-v2.py` offers that the tests run on: the machine
/// builds them for macOS 13 and up. The tags are the fetcher's `-s` names.
/// Sonoma is the one OSX-KVM recommends.
const MacRelease = enum { ventura, sonoma, sequoia, tahoe };

/// A macOS machine's files in `.macos-vm/<name>/`, relative to the build root.
const MacPaths = struct { dmg: []const u8, img: []const u8, disk: []const u8, vars: []const u8 };

/// QEMU on `m`'s command line plus OSX-KVM's firmware, SMC key and OpenCore,
/// and the disks; `install` adds the recovery image.
fn macosQemu(b: *Build, m: Machine, kvm: *Build.Dependency, paths: MacPaths, mac: MacOptions, install: bool) !*Build.Step.Run {
    const run = console(b, &.{m.qemu});
    run.addArgs(m.args);
    run.addArgs(&.{ "-m", m.memory, "-device", b.fmt("isa-applesmc,osk={s}", .{try smcKey(b, kvm)}), "-drive" });
    run.addPrefixedFileArg("if=pflash,format=raw,readonly=on,file=", kvm.path("OVMF_CODE_4M.fd"));
    run.addArgs(&.{ "-drive", b.fmt("if=pflash,format=raw,file={s}", .{paths.vars}), "-drive" });
    run.addPrefixedFileArg("id=OpenCoreBoot,if=none,snapshot=on,format=qcow2,file=", kvm.path("OpenCore/OpenCore.qcow2"));
    run.addArgs(&.{
        "-device", "ide-hd,bus=sata.2,drive=OpenCoreBoot",
        "-drive",  b.fmt("id=MacHDD,if=none,format=qcow2,cache=writeback,file={s}", .{paths.disk}),
        "-device", "ide-hd,bus=sata.4,drive=MacHDD",
        "-netdev", b.fmt("user,id=net0,hostfwd=tcp::{d}-:22", .{mac.port}),
        "-device", "virtio-net-pci,netdev=net0,id=net0,mac=52:54:00:c9:18:27",
    });
    if (install) run.addArgs(&.{
        "-drive",  b.fmt("id=InstallMedia,if=none,format=raw,file={s}", .{paths.img}),
        "-device", "ide-hd,bus=sata.3,drive=InstallMedia",
    });
    // A `(qemu)` prompt in this console: `info status` says whether it runs,
    // and `RIP` in `info registers` moves while the guest executes.
    run.addArgs(&.{ "-monitor", "stdio" });
    // gdb-multiarch, then `target remote localhost:1234`
    if (mac.gdb) run.addArg("-s");
    return run;
}

/// The Apple SMC key OSX-KVM's boot script passes, read from it.
fn smcKey(b: *Build, kvm: *Build.Dependency) ![]const u8 {
    const path = try kvm.builder.root.joinString(b.allocator, "OpenCore-Boot.sh");
    const script = try std.Io.Dir.cwd().readFileAlloc(b.graph.io, path, b.allocator, .limited(1 << 20));
    const start = (std.mem.indexOf(u8, script, "osk=\"") orelse return error.NoSmcKey) + 5;
    const end = std.mem.indexOfScalarPos(u8, script, start, '"') orelse return error.NoSmcKey;
    return script[start..end];
}

fn exists(b: *Build, path: []const u8) bool {
    b.root.access(b.graph.io, path, .{}) catch return false;
    return true;
}

/// A command on the console (it may ask questions or open a window), from the
/// build root.
fn console(b: *Build, argv: []const []const u8) *Build.Step.Run {
    const run = b.addSystemCommand(argv);
    run.setCwd(b.path("."));
    run.stdio = .inherit;
    return run;
}
