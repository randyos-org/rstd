//! `ssh_test <host> <test exe> [args...]`: runs the tests in the host's /tmp with `RSTD_TEST_PROGRESS` set and exits with their status.
//! `local:<path>` arguments are copied over too, `env:<name>=<value>` ones set for the run; ssh and scp come from PATH with their config.

const std = @import("std");
const Io = std.Io;
const remote = @import("remote.zig");

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
    const remote_exe = try arena.print("/tmp/rstd-test-{d}", .{Io.Clock.awake.now(io).toNanoseconds()});

    var scp: std.ArrayList([]const u8) = .empty;
    try scp.append(arena, "scp");
    try scp.appendSlice(arena, &remote.options);
    try scp.appendSlice(arena, &.{ "-q", args[2], try arena.print("{s}:{s}", .{ host, remote_exe }) });
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
            const copy = try arena.print("{s}.{d}", .{ remote_exe, i });
            var cp: std.ArrayList([]const u8) = .empty;
            try cp.append(arena, "scp");
            try cp.appendSlice(arena, &remote.options);
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
    // a Mac sleeps through an ssh session without a tty: `caffeinate -i` there
    const command = try arena.print("chmod +x {0s} && {3s}RSTD_TEST_PROGRESS=1 $([ \"$(uname)\" = Darwin ] && echo caffeinate -i) {0s}{1s}; s=$?; rm -rf {0s}{2s}; exit $s", .{ remote_exe, extra.items, copies.items, env.items });
    return exec(io, try remote.sshArgv(arena, null, &remote.options, host, &.{command}));
}

/// Runs `argv` with this process's stdout and stderr; its exit status.
fn exec(io: Io, argv: []const []const u8) !u8 {
    return remote.exec(io, "ssh_test: ", argv, .inherit, .inherit);
}
