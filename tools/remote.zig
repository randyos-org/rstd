//! The ssh and spawn code of `ssh_test`, `test_host` and `macos_vm`.

const std = @import("std");
const Io = std.Io;
const StdIo = std.process.SpawnOptions.StdIo;

/// `ssh_test`'s and `test_host`'s ssh and scp options.
pub const options = [_][]const u8{ "-o", "BatchMode=yes", "-o", "ConnectTimeout=30", "-o", "ServerAliveInterval=30" };

/// `ssh [-p port] <options> <host> <command>`.
pub fn sshArgv(arena: std.mem.Allocator, port: ?[]const u8, opts: []const []const u8, host: []const u8, command: []const []const u8) ![]const []const u8 {
    var argv: std.ArrayList([]const u8) = .empty;
    try argv.append(arena, "ssh");
    if (port) |p| try argv.appendSlice(arena, &.{ "-p", p });
    try argv.appendSlice(arena, opts);
    try argv.append(arena, host);
    try argv.appendSlice(arena, command);
    return argv.items;
}

/// Spawns `argv`; a missing program prints `<who><program> is not on PATH`.
pub fn spawn(io: Io, who: []const u8, argv: []const []const u8, stdin: StdIo, stdout: StdIo, stderr: StdIo) !std.process.Child {
    return std.process.spawn(io, .{ .argv = argv, .stdin = stdin, .stdout = stdout, .stderr = stderr }) catch |err| switch (err) {
        error.FileNotFound => {
            std.debug.print("{s}{s} is not on PATH\n", .{ who, argv[0] });
            return error.FileNotFound;
        },
        else => |e| return e,
    };
}

/// Runs `argv` to completion, stdin ignored; its exit status.
pub fn exec(io: Io, who: []const u8, argv: []const []const u8, stdout: StdIo, stderr: StdIo) !u8 {
    var child = try spawn(io, who, argv, .ignore, stdout, stderr);
    return status(try child.wait(io));
}

/// An exit status; 255 when the child did not exit.
pub fn status(term: std.process.Child.Term) u8 {
    return switch (term) {
        .exited => |code| code,
        else => 255,
    };
}
