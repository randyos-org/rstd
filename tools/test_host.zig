//! Runs a qemu rig's QEMU on another machine over SSH (`-Dtest-host`):
//! copies the rig's files and `tools/nix` (the flake that pins QEMU, and
//! `run.sh`) into a new directory in the host's /tmp, runs QEMU there
//! through the flake, and passes its console to this stdout and the live
//! log. `run.sh` waits for one of `jobs` slots, ends QEMU when this ssh's
//! stdin closes (this process died) and removes the directory.
//!
//!   test_host <host> <live log> <jobs> <nix dir> <qemu> --tcg <args>...
//!     [--hvf <args>...] -- <args>...
//!
//! An argument `local:<path>` is a file copied over, passed as its copy's
//! path. A macOS arm64 host takes the `--hvf` arguments when given, any
//! other the `--tcg` ones. `ssh` and `scp` come from PATH with their own
//! config; nothing about the machine is stored here.

const std = @import("std");
const Io = std.Io;

pub fn main(init: std.process.Init) void {
    std.process.exit(run(init) catch |err| {
        std.debug.print("test_host: {t}\n", .{err});
        std.process.exit(1);
    });
}

const options = [_][]const u8{ "-o", "BatchMode=yes", "-o", "ConnectTimeout=30", "-o", "ServerAliveInterval=30" };

const local_prefix = "local:";

fn run(init: std.process.Init) !u8 {
    const io = init.io;
    const arena = init.arena.allocator();
    const args = try init.minimal.args.toSlice(arena);
    if (args.len < 8 or !std.mem.eql(u8, args[6], "--tcg")) return usage();
    const host = args[1];
    const log_path = args[2];
    const jobs = args[3];
    const nix_dir = args[4];
    const qemu = args[5];
    const dashes = indexOf(args, 7, "--") orelse return usage();
    const hvf_at = indexOf(args[0..dashes], 7, "--hvf");
    const tcg = args[7 .. hvf_at orelse dashes];
    const hvf: ?[]const []const u8 = if (hvf_at) |h| args[h + 1 .. dashes] else null;
    const rest = args[dashes + 1 ..];

    // the host's kind, and a directory of this run's own
    const probe = try capture(io, arena, try sshArgv(arena, host, "uname -sm && mktemp -d /tmp/rstd-qemu.XXXXXX"));
    var lines = std.mem.tokenizeAny(u8, probe, "\r\n");
    const kind = lines.next() orelse return error.HostProbe;
    const dir = lines.next() orelse return error.HostProbe;
    const use_hvf = hvf != null and std.mem.eql(u8, kind, "Darwin arm64");

    // the files to copy, by name: one directory holds them all
    var scp: std.ArrayList([]const u8) = .empty;
    try scp.append(arena, "scp");
    try scp.appendSlice(arena, &options);
    try scp.appendSlice(arena, &.{ "-q", "-r", nix_dir });
    var names: std.ArrayList([]const u8) = .empty;
    try names.append(arena, std.Io.Dir.path.basename(nix_dir));
    var qemu_argv: std.ArrayList([]const u8) = .empty;
    try qemu_argv.append(arena, qemu);
    try qemu_argv.appendSlice(arena, if (use_hvf) hvf.? else tcg);
    for (rest) |a| {
        if (!std.mem.startsWith(u8, a, local_prefix)) {
            try qemu_argv.append(arena, a);
            continue;
        }
        const path = a[local_prefix.len..];
        const name = std.Io.Dir.path.basename(path);
        for (names.items) |n| if (std.mem.eql(u8, n, name)) return error.DuplicateFileName;
        try names.append(arena, name);
        try scp.append(arena, path);
        try qemu_argv.append(arena, try arena.print("{s}/{s}", .{ dir, name }));
    }
    try scp.append(arena, try arena.print("{s}:{s}/", .{ host, dir }));
    if (try exec(io, scp.items) != 0) {
        std.debug.print("test_host: copying the rig's files to {s} failed\n", .{host});
        _ = exec(io, try sshArgv(arena, host, try arena.print("rm -rf '{s}'", .{dir}))) catch {};
        return 1;
    }

    // each argument single-quoted for the remote shell
    var command: std.ArrayList(u8) = .empty;
    const nix = names.items[0];
    try command.print(arena, "cd '{0s}' && exec nix shell ./{1s}#qemu ./{1s}#flock -c bash {1s}/run.sh '{2s}'", .{ dir, nix, jobs });
    for (qemu_argv.items) |a| {
        if (std.mem.findScalar(u8, a, '\'') != null) return error.QuoteInArgument;
        try command.print(arena, " '{s}'", .{a});
    }
    const log = try Io.Dir.cwd().createFile(io, log_path, .{});
    defer log.close(io);
    var out: Io.File.Writer = .init(log, io, &.{});
    const stdout = Io.File.stdout();
    const banner = try arena.print("test-host: {s} ({s}, {s})\n", .{ host, kind, if (use_hvf) "hvf" else "tcg" });
    stdout.writeStreamingAll(io, banner) catch {};
    try out.interface.writeAll(banner);

    // stdin stays open while this process lives: `run.sh` ends QEMU when it closes
    var ssh = try spawn(io, try sshArgv(arena, host, command.items), .pipe, .pipe);
    var buf: [4096]u8 = undefined;
    while (true) {
        const n = ssh.stdout.?.readStreaming(io, &.{&buf}) catch break;
        if (n == 0) continue;
        stdout.writeStreamingAll(io, buf[0..n]) catch {};
        out.interface.writeAll(buf[0..n]) catch {};
    }
    return switch (try ssh.wait(io)) {
        .exited => |code| code,
        else => 255,
    };
}

fn usage() u8 {
    std.debug.print("usage: test_host <host> <live log> <jobs> <nix dir> <qemu> --tcg <args>... [--hvf <args>...] -- <args>...\n", .{});
    return 2;
}

fn indexOf(args: []const []const u8, from: usize, arg: []const u8) ?usize {
    for (args[from..], from..) |a, i| if (std.mem.eql(u8, a, arg)) return i;
    return null;
}

fn sshArgv(arena: std.mem.Allocator, host: []const u8, command: []const u8) ![]const []const u8 {
    var argv: std.ArrayList([]const u8) = .empty;
    try argv.append(arena, "ssh");
    try argv.appendSlice(arena, &options);
    try argv.appendSlice(arena, &.{ host, command });
    return argv.items;
}

fn spawn(io: Io, argv: []const []const u8, stdin: std.process.SpawnOptions.StdIo, stdout: std.process.SpawnOptions.StdIo) !std.process.Child {
    return std.process.spawn(io, .{ .argv = argv, .stdin = stdin, .stdout = stdout, .stderr = .inherit }) catch |err| switch (err) {
        error.FileNotFound => {
            std.debug.print("test_host: {s} is not on PATH\n", .{argv[0]});
            return error.FileNotFound;
        },
        else => |e| return e,
    };
}

/// Runs `argv` with this process's stdout; its exit status (255 when it did
/// not exit).
fn exec(io: Io, argv: []const []const u8) !u8 {
    var child = try spawn(io, argv, .ignore, .inherit);
    return switch (try child.wait(io)) {
        .exited => |code| code,
        else => 255,
    };
}

/// `argv`'s stdout, failing unless it exits 0.
fn capture(io: Io, arena: std.mem.Allocator, argv: []const []const u8) ![]const u8 {
    var child = try spawn(io, argv, .ignore, .pipe);
    var reader = child.stdout.?.readerStreaming(io, &.{});
    const text = try reader.interface.allocRemaining(arena, .limited(1 << 16));
    return switch (try child.wait(io)) {
        .exited => |code| if (code == 0) text else error.HostProbe,
        else => error.HostProbe,
    };
}
