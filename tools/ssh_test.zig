//! Runs a test executable on another machine over SSH: copies it into the
//! host's /tmp, runs it there with `RSTD_TEST_PROGRESS` set (so tests that
//! run every core name each one) and any further arguments, removes it, and
//! exits with its status. `ssh` and `scp` come from PATH with their own
//! config, so `host` is an SSH alias or name; nothing about the machine is
//! stored here.
//!
//!   ssh_test <host> <test exe> [args...]
//!
//! An argument `local:<path>` is a file or directory copied over too, next
//! to the executable, passed as its copy's path and removed after; one
//! `env:<name>=<value>` sets that variable for the run instead.

const std = @import("std");
const Io = std.Io;

const local_prefix = "local:";
const env_prefix = "env:";

pub fn main(init: std.process.Init) void {
    std.process.exit(run(init) catch |err| {
        std.debug.print("ssh_test: {t}\n", .{err});
        std.process.exit(1);
    });
}

fn run(init: std.process.Init) !u8 {
    const io = init.io;
    const arena = init.arena.allocator();
    const args = try init.minimal.args.toSlice(arena);
    if (args.len < 3) {
        std.debug.print("usage: ssh_test <host> <test exe> [args...]\n", .{});
        return 2;
    }
    const host = args[1];
    // unique, so runs from two checkouts do not collide
    const remote = try arena.print("/tmp/rstd-test-{d}", .{Io.Clock.awake.now(io).toNanoseconds()});
    const options = [_][]const u8{ "-o", "BatchMode=yes", "-o", "ConnectTimeout=30", "-o", "ServerAliveInterval=30" };

    var scp: std.ArrayList([]const u8) = .empty;
    try scp.append(arena, "scp");
    try scp.appendSlice(arena, &options);
    try scp.appendSlice(arena, &.{ "-q", args[2], try arena.print("{s}:{s}", .{ host, remote }) });
    if (try exec(io, scp.items) != 0) {
        std.debug.print("ssh_test: copying the tests to {s} failed\n", .{host});
        return 1;
    }

    // each argument single-quoted for the remote shell; `local:` ones copied
    var extra: std.ArrayList(u8) = .empty;
    var copies: std.ArrayList(u8) = .empty;
    var env: std.ArrayList(u8) = .empty;
    for (args[3..], 0..) |a, i| {
        if (std.mem.startsWith(u8, a, env_prefix)) {
            const set = a[env_prefix.len..];
            const eq = std.mem.findScalar(u8, set, '=') orelse return error.EnvWithoutValue;
            if (std.mem.findScalar(u8, set, '\'') != null) return error.QuoteInArgument;
            try env.print(arena, "{s}='{s}' ", .{ set[0..eq], set[eq + 1 ..] });
            continue;
        }
        const arg = if (std.mem.startsWith(u8, a, local_prefix)) blk: {
            const path = a[local_prefix.len..];
            const copy = try arena.print("{s}.{d}", .{ remote, i });
            var cp: std.ArrayList([]const u8) = .empty;
            try cp.append(arena, "scp");
            try cp.appendSlice(arena, &options);
            try cp.appendSlice(arena, &.{ "-q", "-r", path, try arena.print("{s}:{s}", .{ host, copy }) });
            if (try exec(io, cp.items) != 0) {
                std.debug.print("ssh_test: copying {s} to {s} failed\n", .{ path, host });
                return 1;
            }
            try copies.print(arena, " {s}", .{copy});
            break :blk copy;
        } else a;
        if (std.mem.findScalar(u8, arg, '\'') != null) return error.QuoteInArgument;
        try extra.print(arena, " '{s}'", .{arg});
    }
    const command = try arena.print("chmod +x {0s} && {3s}RSTD_TEST_PROGRESS=1 {0s}{1s}; s=$?; rm -rf {0s}{2s}; exit $s", .{ remote, extra.items, copies.items, env.items });
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
            std.debug.print("ssh_test: {s} is not on PATH\n", .{argv[0]});
            return error.FileNotFound;
        },
        else => |e| return e,
    };
    return switch (try child.wait(io)) {
        .exited => |code| code,
        else => 255,
    };
}
