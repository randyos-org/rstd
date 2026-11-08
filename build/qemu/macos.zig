//! The macOS VM rig (`x86_64-macos` in `machines.zig`): its setup, control and test steps; set up once by `tools/macos/readme.md`.

const std = @import("std");
const rstdbuild = @import("../../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const console = rstdbuild.console;
const buildroot = @import("../__root__.zig");
const tools = @import("../tools.zig");
const machines = @import("machines.zig");

const Machine = machines.Machine;

/// `m`'s `qemu-setup-`, `qemu-install-` and `qemu-start-<name>` steps, and the
/// test run: rstd's unit tests (`-Dqemu-test-filter` narrows them) copied
/// into the running VM and run over SSH.
/// Null until the lazy OSX-KVM dependency is fetched. The VM's files are named
/// relative to the build root, where its commands run, so none is hashed.
pub fn addRun(b: *Build, m: Machine, target: rstdbuild.ResolvedTarget, optimize: Optimize, filters: []const []const u8, mac: Options) !?*Build.Step.Run {
    const kvm = b.dependencyLazy("osx_kvm", .{}) catch return null;
    const vm = b.fmt(".macos-vm/{s}", .{m.name});
    try b.root.createDirPath(b.graph.io, vm);
    // setup skips what is there, so a change there configures again
    b.dependOnDirectoryContents(b.path(vm));
    const paths: Paths = .{
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
    // OpenCore's own copy, writable: `configure_macos_vm.sh` edits its config
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

    const tool = tools.hostTool(b, .macos_vm);
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
    const exe = try buildroot.addRstdTest(b, mod, target, optimize, filters);
    run.addArtifactArg2(exe, .{});
    run.addArg(m.qemu);
    // what that runs, built alone, so a build queue frees its slot before the VM boots
    const check = b.step(b.fmt("qemu-check-{s}", .{m.name}), b.fmt("Build what qemu-test-{s} runs, without starting the VM", .{m.name}));
    check.dependOn(&exe.step);
    check.dependOn(&tool.step);
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
pub const Options = struct { user: []const u8 = "", port: u16 = 2222, release: Release = .sonoma, gdb: bool = false, detach: bool = false, headless: bool = false, keep: bool = false, mttcg: ?bool = null, ovmf_vars: ?[]const u8 = null, snapshot: []const u8 = "clean" };

/// The releases `fetch-macOS-v2.py` offers that the tests run on: the machine
/// builds them for macOS 13 and up. The tags are the fetcher's `-s` names.
/// Sonoma is the one OSX-KVM recommends.
pub const Release = enum { ventura, sonoma, sequoia, tahoe };

/// A macOS machine's files in `.macos-vm/<name>/`, relative to the build root.
const Paths = struct { dmg: []const u8, img: []const u8, disk: []const u8, vars: []const u8, opencore: []const u8 };

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
fn qemuArgs(b: *Build, run: *Build.Step.Run, m: Machine, kvm: *Build.Dependency, paths: Paths, mac: Options, mode: QemuMode) !void {
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
    const dir = std.Io.Dir.path.dirname(paths.disk).?;
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
    const start = (std.mem.find(u8, script, "osk=\"") orelse return error.NoSmcKey) + 5;
    const end = std.mem.findScalarPos(u8, script, start, '"') orelse return error.NoSmcKey;
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
    const bin = std.Io.Dir.path.dirname(qemu) orelse return null;
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
