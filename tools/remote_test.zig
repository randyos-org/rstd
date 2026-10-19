//! Runs a test executable on another machine over SSH: copies it into the
//! host's /tmp, runs it there with `RSTD_TEST_PROGRESS` set (so tests that
//! run every core name each one) and any further arguments, removes it, and
//! exits with its status. `ssh` and `scp` come from PATH with their own
//! config, so `host` is an SSH alias or name; nothing about the machine is
//! stored here.
//!
//!   remote_test <host> <test exe> [args...]

const std = @import("std");
const Io = std.Io;

pub fn main(init: std.process.Init) void {
    std.process.exit(run(init) catch |err| {
        std.debug.print("remote_test: {t}\n", .{err});
        std.process.exit(1);
    });
}

fn run(init: std.process.Init) !u8 {
    const io = init.io;
    const arena = init.arena.allocator();
    const args = try init.minimal.args.toSlice(arena);
    if (args.len < 3) {
        std.debug.print("usage: remote_test <host> <test exe> [args...]\n", .{});
        return 2;
    }
    const host = args[1];
    // unique, so runs from two checkouts do not collide
    const remote = try std.fmt.allocPrint(arena, "/tmp/rstd-test-{d}", .{Io.Clock.awake.now(io).toNanoseconds()});
    const options = [_][]const u8{ "-o", "BatchMode=yes", "-o", "ConnectTimeout=30", "-o", "ServerAliveInterval=30" };

    var scp: std.ArrayList([]const u8) = .empty;
    try scp.append(arena, "scp");
    try scp.appendSlice(arena, &options);
    try scp.appendSlice(arena, &.{ "-q", args[2], try std.fmt.allocPrint(arena, "{s}:{s}", .{ host, remote }) });
    if (try exec(io, scp.items) != 0) {
        std.debug.print("remote_test: copying the tests to {s} failed\n", .{host});
        return 1;
    }

    // each argument single-quoted for the remote shell
    var extra: std.ArrayList(u8) = .empty;
    for (args[3..]) |a| {
        if (std.mem.indexOfScalar(u8, a, '\'') != null) return error.QuoteInArgument;
        try extra.print(arena, " '{s}'", .{a});
    }
    const command = try std.fmt.allocPrint(arena, "chmod +x {0s} && RSTD_TEST_PROGRESS=1 {0s}{1s}; s=$?; rm -f {0s}; exit $s", .{ remote, extra.items });
    var ssh: std.ArrayList([]const u8) = .empty;
    try ssh.append(arena, "ssh");
    try ssh.appendSlice(arena, &options);
    try ssh.appendSlice(arena, &.{ host, command });
    return exec(io, ssh.items);
}

/// Runs `argv` with this process's stdio; its exit status (255 when it did
/// not exit).
fn exec(io: Io, argv: []const []const u8) !u8 {
    var child = std.process.spawn(io, .{ .argv = argv, .stdin = .ignore, .stdout = .inherit, .stderr = .inherit }) catch |err| switch (err) {
        error.FileNotFound => {
            std.debug.print("remote_test: {s} is not on PATH\n", .{argv[0]});
            return error.FileNotFound;
        },
        else => |e| return e,
    };
    return switch (try child.wait(io)) {
        .exited => |code| code,
        else => 255,
    };
}
