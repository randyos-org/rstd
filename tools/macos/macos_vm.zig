//! The macOS VM for the tests: `test` starts QEMU headless (detached, log in `<dir>/qemu.log`) when SSH does not answer, runs the tests
//! and shuts down a VM it started unless `keep`; `up`, `down`, `kill` (`<dir>/qemu.pid`) and `guard` (fails while it runs); usage in `run`.

const std = @import("std");
const sysinfo = @import("builtin");
const Io = std.Io;
const remote = @import("remote");

/// How often a shutdown is checked for a frozen VM, and how long nothing
/// may change before it counts as one.
const freeze_poll_s = 15;
const freeze_s = 60;
const poll_s = 10;
/// How often a boot that prints nothing says it is still going.
const beat_s = 30;
const windows = sysinfo.target.os.tag == .windows;

const Vm = struct {
    io: Io,
    arena: std.mem.Allocator,
    dir: []const u8,
    port: []const u8 = "",
    login: []const u8 = "",

    fn path(vm: Vm, name: []const u8) ![]const u8 {
        return vm.arena.print("{s}/{s}", .{ vm.dir, name });
    }
};

pub fn main(init: std.process.Init) void {
    // Exits rather than returns: returning tears down `init.io`, which
    // waits for the task watching a qemu that keeps running.
    std.process.exit(run(init) catch |err| {
        std.debug.print("macOS VM: {t}\n", .{err});
        std.process.exit(1);
    });
}

fn run(init: std.process.Init) !u8 {
    const arena = init.arena.allocator();
    const args = try init.minimal.args.toSlice(arena);
    const usage =
        \\usage: macos_vm test <dir> <port> <login> <keep> <test exe> <qemu argv...>
        \\       macos_vm up <dir> <port> <login> <qemu argv...>
        \\       macos_vm down <dir> <port> <login>
        \\       macos_vm kill <dir>
        \\       macos_vm guard <dir>
        \\
    ;
    const Command = enum { @"test", up, down, kill, guard };
    const command = if (args.len >= 3) std.meta.stringToEnum(Command, args[1]) else null;
    const min: usize = switch (command orelse .kill) {
        .@"test" => 8,
        .up => 6,
        .down => 5,
        .kill, .guard => 3,
    };
    if (command == null or args.len < min) {
        std.debug.print(usage, .{});
        return 2;
    }
    var vm: Vm = .{ .io = init.io, .arena = arena, .dir = args[2] };
    if (min >= 5) {
        vm.port = args[3];
        vm.login = args[4];
    }
    switch (command.?) {
        .kill => try kill(vm),
        .guard => if (try qemuPid(vm)) |pid| if (try running(vm, pid)) {
            std.debug.print("macOS VM: already running (qemu pid {s}; if headless, VNC on 127.0.0.1:5901); qemu-down or qemu-kill it first\n", .{pid});
            return 1;
        },
        .down => try shutdown(vm),
        .up => try up(vm, args[5..], false),
        .@"test" => {
            try up(vm, args[7..], true);
            const ok = try runTests(vm, args[6]);
            // ours to stop: the marker outlives an interrupted run
            if (!std.mem.eql(u8, args[5], "true") and try fileExists(vm, "autostarted")) shutdown(vm) catch {};
            return if (ok) 0 else 1;
        },
    }
    return 0;
}

/// Returns once SSH answers, starting QEMU first when it does not;
/// `mark` has a test run stop it afterward. A QEMU already running is
/// waited for, never started again: the second would fail on the disk lock
/// and, exiting, delete the first one's pid file.
fn up(vm: Vm, qemu_argv: []const []const u8, mark: bool) !void {
    const io = vm.io;
    if (try ssh(vm, &.{"true"}, .quiet)) return;
    if (try qemuPid(vm)) |pid| if (try running(vm, pid)) {
        std.debug.print("macOS VM: qemu (pid {s}) is running but SSH does not answer yet; waiting\n", .{pid});
        return waitUp(vm, null);
    };

    const log_path = try vm.path("qemu.log");
    std.debug.print("macOS VM: not answering on port {s}; starting it headless (log: {s})\n", .{ vm.port, log_path });
    const log = try Io.Dir.cwd().createFile(io, log_path, .{});
    var qemu = try std.process.spawn(io, .{
        .argv = qemu_argv,
        .stdin = .ignore,
        .stdout = .{ .file = log },
        .stderr = .{ .file = log },
        // out of reach of the terminal's Ctrl+C, which `abort` handles
        .create_no_window = true,
        .pgid = if (windows) null else 0,
    });
    log.close(io);
    if (mark) try Io.Dir.cwd().writeFile(io, .{ .sub_path = try vm.path("autostarted"), .data = "" });

    starting = .{ .vm = vm, .id = qemu.id.? };
    defer starting = null;
    onInterrupt();
    var exited: std.atomic.Value(bool) = .init(false);
    // qemu keeps running; `watch` is never awaited unless it exits
    var watch = try io.concurrent(waitQemu, .{ &qemu, io, &exited });
    waitUp(vm, &exited) catch |err| {
        if (exited.load(.acquire)) _ = watch.await(io);
        return err;
    };
}

/// Waits for SSH while the VM boots, printing its serial console. A VM
/// that looks frozen is killed only if the user says so.
fn waitUp(vm: Vm, exited: ?*std.atomic.Value(bool)) !void {
    const io = vm.io;
    const serial_path = try vm.path("serial.log");
    var serial_at: u64 = 0;
    var last: ?u64 = null;
    var idle_since: u64 = 0;
    // wall clock: an ssh try or monitor query takes a while of its own
    const started: Io.Timestamp = .now(io, .awake);
    var next_ssh: u64 = poll_s;
    var next_check: u64 = freeze_poll_s;
    var next_beat: u64 = beat_s;
    while (true) {
        try io.sleep(.fromSeconds(1), .awake);
        echoSerial(io, serial_path, &serial_at);
        const elapsed: u64 = @intCast(started.untilNow(io, .awake).toSeconds());
        if (exited) |e| if (e.load(.acquire)) {
            Io.Dir.cwd().deleteFile(io, try vm.path("autostarted")) catch {};
            std.debug.print("macOS VM: qemu exited before SSH answered; see {s}\n", .{try vm.path("qemu.log")});
            return error.VmFailed;
        };
        if (elapsed >= next_ssh) {
            if (try ssh(vm, &.{"true"}, .quiet)) {
                std.debug.print("macOS VM: up after {d} s\n", .{started.untilNow(io, .awake).toSeconds()});
                return;
            }
            next_ssh = @intCast(started.untilNow(io, .awake).toSeconds() + poll_s);
        }
        if (elapsed < next_check) continue;
        const now = activity(vm) catch null;
        const at: u64 = @intCast(started.untilNow(io, .awake).toSeconds());
        if (now == null or now != last) idle_since = at;
        last = now;
        next_check = at + freeze_poll_s;
        if (at >= next_beat) {
            std.debug.print("macOS VM: booting, {d} s; {s}\n", .{ at, if (at - idle_since < freeze_poll_s) "vCPUs and disks active" else "no change in vCPUs or disks for a while" });
            next_beat = at + beat_s;
        }
        if (at - idle_since < freeze_s) continue;
        std.debug.print("macOS VM: looks frozen: no vCPU has moved and no disk I/O for {d} s (see the screen on VNC 127.0.0.1:5901)\n", .{at - idle_since});
        if (try ask(io, "kill qemu?")) {
            try kill(vm);
            return error.VmFailed;
        }
        idle_since = at;
    }
}

/// The QEMU `up` started, until SSH answers: Ctrl+C before then aborts the
/// start, so it is killed rather than left booting with no one waiting.
var starting: ?struct { vm: Vm, id: std.process.Child.Id } = null;

extern "kernel32" fn SetConsoleCtrlHandler(handler: ?*const fn (u32) callconv(.winapi) i32, add: i32) callconv(.winapi) i32;
extern "kernel32" fn TerminateProcess(process: std.os.windows.HANDLE, code: u32) callconv(.winapi) i32;

fn onInterrupt() void {
    if (windows) {
        _ = SetConsoleCtrlHandler(&abortWindows, 1);
    } else {
        const act: std.posix.Sigaction = .{
            .handler = .{ .handler = abortPosix },
            .mask = std.posix.sigemptyset(),
            .flags = 0,
        };
        std.posix.sigaction(.INT, &act, null);
    }
}

fn abortWindows(_: u32) callconv(.winapi) i32 {
    const s = starting orelse return 0;
    _ = TerminateProcess(s.id, 1);
    std.debug.print("\nmacOS VM: start aborted; qemu killed\n", .{});
    // a handler thread, so the files can go too
    for ([_][]const u8{ "qemu.pid", "autostarted" }) |name| {
        const path = s.vm.path(name) catch continue;
        Io.Dir.cwd().deleteFile(s.vm.io, path) catch {};
    }
    std.process.exit(130);
}

fn abortPosix(_: std.posix.SIG) callconv(.c) void {
    // a signal handler: only the kill; `kill` clears the stale pid file
    if (starting) |s| std.posix.kill(s.id, .TERM) catch {};
    std.process.exit(130);
}

/// Prints what QEMU wrote to `path` (`-serial file:`) past `at`.
fn echoSerial(io: Io, path: []const u8, at: *u64) void {
    const file = Io.Dir.cwd().openFile(io, path, .{}) catch return;
    defer file.close(io);
    var buf: [4096]u8 = undefined;
    while (true) {
        const n = file.readPositionalAll(io, &buf, at.*) catch return;
        if (n == 0) return;
        printPlain(buf[0..n]);
        at.* += n;
    }
}

/// Prints `text` without its escape sequences: the firmware's clear
/// screen, cursor moves and window resizes would act on this terminal.
fn printPlain(text: []const u8) void {
    var i: usize = 0;
    var start: usize = 0;
    while (i < text.len) {
        if (text[i] != 0x1b) {
            i += 1;
            continue;
        }
        std.debug.print("{s}", .{text[start..i]});
        i += 1;
        // CSI: `ESC [`, parameters, then a final byte in 0x40..0x7e
        if (i < text.len and text[i] == '[') {
            i += 1;
            while (i < text.len and (text[i] < 0x40 or text[i] > 0x7e)) i += 1;
        }
        i = @min(i + 1, text.len);
        start = i;
    }
    std.debug.print("{s}", .{text[start..]});
}

/// The size of the serial log so far; 0 when there is none.
fn serialLength(io: Io, path: []const u8) u64 {
    const file = Io.Dir.cwd().openFile(io, path, .{}) catch return 0;
    defer file.close(io);
    return file.length(io) catch 0;
}

fn waitQemu(qemu: *std.process.Child, io: Io, exited: *std.atomic.Value(bool)) void {
    _ = qemu.wait(io) catch {};
    exited.store(true, .release);
}

/// Copies `exe` into the VM and runs it there, its output in this console.
fn runTests(vm: Vm, exe: []const u8) !bool {
    const dest = try vm.arena.print("{s}:/tmp/rstd-test", .{vm.login});
    var argv: std.ArrayList([]const u8) = .empty;
    try argv.appendSlice(vm.arena, &.{ "scp", "-q", "-P", vm.port });
    try argv.appendSlice(vm.arena, &patient_options);
    try argv.appendSlice(vm.arena, &.{ exe, dest });
    // `scp` exits 1, not 255, when the connection drops, which `patient`
    // takes as final; the first copy after boot has been reset (T137)
    for (0..patient_tries) |attempt| {
        if (try patient(vm, argv.items, .inherit) == 0) break;
        if (attempt + 1 == patient_tries) {
            std.debug.print("macOS VM: copying the tests in failed\n", .{});
            return false;
        }
        std.debug.print("macOS VM: copying the tests in failed; trying again in {d} s\n", .{patient_pause_s});
        try vm.io.sleep(.fromSeconds(patient_pause_s), .awake);
    }
    // Windows' scp gives the copy no execute bit
    return try sshPatient(vm, &.{ "chmod", "+x", "/tmp/rstd-test", "&&", "/tmp/rstd-test" }, .inherit) == 0;
}

/// Shuts the VM down from inside, which needs passwordless `sudo shutdown`
/// (`tools/macos/configure_macos_vm.sh`), then kills QEMU if it is still there.
fn shutdown(vm: Vm) !void {
    std.debug.print("macOS VM: shutting down\n", .{});
    // `-l` lists the rule without running anything
    if (try sshPatient(vm, &.{ "sudo", "-n", "-l", "/sbin/shutdown" }, .quiet) != 0) {
        std.debug.print("macOS VM: no passwordless sudo shutdown (run tools/macos/configure_macos_vm.sh in the VM); left running\n", .{});
        return error.NoShutdown;
    }
    // the serial console from here on, not the boot again
    const serial_path = try vm.path("serial.log");
    var serial_at: u64 = serialLength(vm.io, serial_path);
    // the connection drops as it goes down, so the result says little
    _ = try execCode(vm.io, try sshArgv(vm, &patient_options, &.{ "sudo", "-n", "/sbin/shutdown", "-h", "now" }), .quiet);
    const pid = try qemuPid(vm) orelse return;
    // Never killed for taking long: it may be writing the disk out. Only a
    // VM that looks frozen is, and only if the user says so.
    var waited: u32 = 0;
    var last: ?u64 = null;
    var same: u32 = 0;
    while (try running(vm, pid)) : (waited += 1) {
        try vm.io.sleep(.fromSeconds(1), .awake);
        echoSerial(vm.io, serial_path, &serial_at);
        if (waited % freeze_poll_s != 0) continue;
        const now = activity(vm) catch null;
        same = if (now != null and now == last) same + 1 else 0;
        last = now;
        if (same * freeze_poll_s < freeze_s) continue;
        std.debug.print("macOS VM: looks frozen: no vCPU has moved and no disk I/O for {d} s\n", .{freeze_s});
        if (try ask(vm.io, "kill qemu? the disk may need repair or a restore")) return kill(vm);
        same = 0;
    }
    echoSerial(vm.io, serial_path, &serial_at);
    std.debug.print("macOS VM: down after {d} s\n", .{waited});
    Io.Dir.cwd().deleteFile(vm.io, try vm.path("qemu.pid")) catch {};
    Io.Dir.cwd().deleteFile(vm.io, try vm.path("autostarted")) catch return;
}

/// A digest of what the VM is doing: every vCPU's instruction pointer and
/// the disks' operation counts, from the QEMU monitor. Equal digests a
/// while apart mean nothing ran: a hang, or a halt with interrupts off.
fn activity(vm: Vm) !?u64 {
    var hash: std.hash.Wyhash = .init(0);
    const regs = try monitor(vm, "info registers -a") orelse return null;
    var lines = std.mem.splitScalar(u8, regs, '\n');
    while (lines.next()) |line| if (std.mem.find(u8, line, "RIP=")) |at| hash.update(line[at..@min(line.len, at + 20)]);
    const disks = try monitor(vm, "info blockstats") orelse return null;
    lines = std.mem.splitScalar(u8, disks, '\n');
    while (lines.next()) |line| {
        var fields = std.mem.tokenizeScalar(u8, line, ' ');
        while (fields.next()) |f| if (std.mem.startsWith(u8, f, "rd_operations=") or std.mem.startsWith(u8, f, "wr_operations=")) hash.update(f);
    }
    return hash.final();
}

/// The output of HMP `command` on the VM's monitor socket; null when
/// there is no monitor to ask.
fn monitor(vm: Vm, command: []const u8) !?[]const u8 {
    const address = Io.net.UnixAddress.init(try vm.path("monitor.sock")) catch return null;
    const stream = address.connect(vm.io) catch return null;
    defer stream.close(vm.io);
    var got: std.ArrayList(u8) = .empty;
    // the greeting, then the command's echo and output, each up to a prompt
    if (!try untilPrompt(vm, stream, &got)) return null;
    got.clearRetainingCapacity();
    var buf: [256]u8 = undefined;
    var w = stream.writer(vm.io, &buf);
    w.interface.print("{s}\n", .{command}) catch return null;
    w.interface.flush() catch return null;
    if (!try untilPrompt(vm, stream, &got)) return null;
    return got.items;
}

const prompt = "(qemu) ";

fn untilPrompt(vm: Vm, stream: Io.net.Stream, got: *std.ArrayList(u8)) !bool {
    var buf: [4096]u8 = undefined;
    while (std.mem.find(u8, got.items, prompt) == null) {
        var bufs = [_][]u8{&buf};
        const n = stream.read(vm.io, &bufs) catch return false;
        if (n == 0) return false;
        try got.appendSlice(vm.arena, buf[0..n]);
    }
    return true;
}

/// Asks a yes/no question on the console; no unless answered yes.
fn ask(io: Io, question: []const u8) !bool {
    std.debug.print("{s} [y/N] ", .{question});
    var buf: [256]u8 = undefined;
    var r = Io.File.stdin().readerStreaming(io, &buf);
    const line = (r.interface.takeDelimiter('\n') catch return false) orelse return false;
    const answer = std.mem.trim(u8, line, " \t\r");
    return answer.len != 0 and (answer[0] == 'y' or answer[0] == 'Y');
}

/// Ends the QEMU in `qemu.pid`, if it is still running.
fn kill(vm: Vm) !void {
    const io = vm.io;
    if (try qemuPid(vm)) |pid| {
        if (try running(vm, pid)) {
            std.debug.print("macOS VM: killing qemu (pid {s})\n", .{pid});
            const argv: []const []const u8 = if (windows) &.{ "taskkill", "/F", "/PID", pid } else &.{ "kill", pid };
            _ = try exec(io, argv, .quiet);
            while (try running(vm, pid)) try io.sleep(.fromMilliseconds(200), .awake);
        } else std.debug.print("macOS VM: not running\n", .{});
    } else std.debug.print("macOS VM: no {s}; not started by a build step\n", .{try vm.path("qemu.pid")});
    // QEMU removes its pid file at a clean exit, not when killed
    Io.Dir.cwd().deleteFile(io, try vm.path("qemu.pid")) catch {};
    Io.Dir.cwd().deleteFile(io, try vm.path("autostarted")) catch return;
}

/// The pid in `qemu.pid`, as text.
fn qemuPid(vm: Vm) !?[]const u8 {
    const text = Io.Dir.cwd().readFileAlloc(vm.io, try vm.path("qemu.pid"), vm.arena, .limited(64)) catch |err| switch (err) {
        error.FileNotFound => return null,
        else => |e| return e,
    };
    const pid = std.mem.trim(u8, text, " \r\n");
    _ = std.fmt.parseInt(u32, pid, 10) catch return null;
    return pid;
}

/// Whether `pid` is a running QEMU; a stale pid file's pid may since have
/// gone to something else.
fn running(vm: Vm, pid: []const u8) !bool {
    const argv: []const []const u8 = if (windows)
        &.{ "tasklist", "/NH", "/FO", "CSV", "/FI", try vm.arena.print("PID eq {s}", .{pid}) }
    else
        &.{ "ps", "-p", pid, "-o", "comm=" };
    const result = std.process.run(vm.arena, vm.io, .{ .argv = argv }) catch return false;
    return std.mem.find(u8, result.stdout, "qemu-system") != null;
}

/// Polling a VM that is still booting: a short connect timeout.
const poll_options = [_][]const u8{ "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "-o", "StrictHostKeyChecking=accept-new" };

/// A VM that answered but is saturated (the login window and USB HID drivers
/// can hold its load average near 70 under emulation) takes minutes to
/// answer: a long connect timeout, and a try again while ssh reports the
/// connection failed (exit 255, not the command's own status).
const patient_timeout_s = 180;
const patient_tries = 10;
const patient_pause_s = 15;

const patient_options = [_][]const u8{ "-o", "BatchMode=yes", "-o", std.fmt.comptimePrint("ConnectTimeout={d}", .{patient_timeout_s}), "-o", "ServerAliveInterval=30", "-o", "StrictHostKeyChecking=accept-new" };

/// `ssh <options> login command`.
fn sshArgv(vm: Vm, options: []const []const u8, command: []const []const u8) ![]const []const u8 {
    return remote.sshArgv(vm.arena, vm.port, options, vm.login, command);
}

/// Runs `command` in the VM over SSH; whether it succeeded. QEMU's port
/// forward accepts before sshd is up, so a TCP connect alone says nothing.
fn ssh(vm: Vm, command: []const []const u8, output: Output) !bool {
    return try execCode(vm.io, try sshArgv(vm, &poll_options, command), output) == 0;
}

/// `command` over SSH, patiently; its exit status.
fn sshPatient(vm: Vm, command: []const []const u8, output: Output) !u8 {
    return patient(vm, try sshArgv(vm, &patient_options, command), output);
}

/// Runs an ssh or scp `argv`, again while it could not connect (255).
fn patient(vm: Vm, argv: []const []const u8, output: Output) !u8 {
    var code: u8 = 255;
    for (0..patient_tries) |attempt| {
        code = try execCode(vm.io, argv, output);
        if (code != 255) return code;
        if (attempt + 1 < patient_tries) {
            std.debug.print("macOS VM: no answer yet (too busy?); trying again in {d} s\n", .{patient_pause_s});
            try vm.io.sleep(.fromSeconds(patient_pause_s), .awake);
        }
    }
    return code;
}

const Output = enum { quiet, inherit };

/// Runs `argv` to completion; whether it exited 0.
fn exec(io: Io, argv: []const []const u8, output: Output) !bool {
    return try execCode(io, argv, output) == 0;
}

/// Runs `argv` to completion; its exit status.
fn execCode(io: Io, argv: []const []const u8, output: Output) !u8 {
    const stdio: std.process.SpawnOptions.StdIo = if (output == .quiet) .ignore else .inherit;
    return remote.exec(io, "", argv, stdio, stdio);
}

fn fileExists(vm: Vm, name: []const u8) !bool {
    Io.Dir.cwd().access(vm.io, try vm.path(name), .{}) catch return false;
    return true;
}
