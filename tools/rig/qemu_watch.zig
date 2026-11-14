//! `qemu_watch <idle minutes> [--fail-fast] <halt text>... -- <qemu argv>`: runs a rig's qemu, stdout and stderr passed through, and ends it on a
//! firmware halt text, on silence that long, or with `--fail-fast` at the first failing test. Exits with qemu's status, or 1 after ending it.

const std = @import("std");
const Io = std.Io;

const State = struct {
    /// When qemu last wrote, in seconds of the awake clock.
    last: std.atomic.Value(i64),
    halted: std.atomic.Value(bool) = .init(false),
    /// When the first test failed, in seconds of the awake clock; 0 if none has.
    failed: std.atomic.Value(i64) = .init(0),
    /// Streams that have ended, of two.
    done: std.atomic.Value(u8) = .init(0),
};

/// Both test runners' failure line: `FAIL <error>` and `...FAIL (<error>)`.
const fail_text = "FAIL ";

/// Seconds a failure's trace gets to print before the run is ended.
const fail_grace = 3;

pub fn main(init: std.process.Init) void {
    // exits rather than returns: returning would wait for a reader still blocked on qemu
    std.process.exit(run(init) catch |err| {
        std.debug.print("qemu_watch: {t}\n", .{err});
        std.process.exit(1);
    });
}

fn run(init: std.process.Init) !u8 {
    const io = init.io;
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    const dashes = for (args, 0..) |a, i| {
        if (std.mem.eql(u8, a, "--")) break i;
    } else return error.NoQemuCommand;
    if (dashes < 2 or dashes + 1 == args.len) return error.NoQemuCommand;
    const idle = try std.fmt.parseInt(i64, args[1], 10) * 60;
    const fail_fast = dashes > 2 and std.mem.eql(u8, args[2], "--fail-fast");
    const halts = args[@as(usize, if (fail_fast) 3 else 2)..dashes];

    var qemu = try std.process.spawn(io, .{ .argv = args[dashes + 1 ..], .stdin = .ignore, .stdout = .pipe, .stderr = .pipe });
    var state: State = .{ .last = .init(now(io)) };
    // both: test runners report on stderr
    var out = try io.concurrent(passThrough, .{ io, qemu.stdout.?, Io.File.stdout(), halts, fail_fast, &state });
    var err = try io.concurrent(passThrough, .{ io, qemu.stderr.?, Io.File.stderr(), halts, fail_fast, &state });

    while (state.done.load(.acquire) < 2) {
        try io.sleep(.fromSeconds(2), .awake);
        const failed = state.failed.load(.acquire);
        const why: ?[]const u8 = if (state.halted.load(.acquire))
            "the firmware halted"
        else if (failed != 0 and now(io) - failed >= fail_grace)
            "a test failed (-Dfail-fast)"
        else if (now(io) - state.last.load(.acquire) > idle)
            "no output for too long"
        else
            null;
        if (why) |w| {
            qemu.kill(io);
            std.debug.print("\nqemu_watch: ended qemu: {s}\n", .{w});
            return 1;
        }
    }
    out.await(io);
    err.await(io);
    return switch (try qemu.wait(io)) {
        .exited => |code| code,
        else => 1,
    };
}

fn now(io: Io) i64 {
    return Io.Clock.awake.now(io).toSeconds();
}

/// Copies one of qemu's streams to `out`, noting when it last wrote, any halt
/// text and, for `fail_fast`, the first failure; a chunk boundary may split
/// them, hence the kept tail.
fn passThrough(io: Io, from: Io.File, out: Io.File, halts: []const []const u8, fail_fast: bool, state: *State) void {
    defer _ = state.done.fetchAdd(1, .release);
    var buf: [4096 + 256]u8 = undefined;
    var kept: usize = 0;
    while (true) {
        const n = from.readStreaming(io, &.{buf[kept..]}) catch return;
        if (n == 0) continue;
        out.writeStreamingAll(io, buf[kept..][0..n]) catch {};
        state.last.store(now(io), .release);
        const seen = buf[0 .. kept + n];
        for (halts) |h| if (std.mem.find(u8, seen, h) != null) state.halted.store(true, .release);
        if (fail_fast and state.failed.load(.acquire) == 0 and std.mem.find(u8, seen, fail_text) != null)
            state.failed.store(now(io), .release);
        kept = @min(seen.len, 256);
        std.mem.copyForwards(u8, buf[0..kept], seen[seen.len - kept ..]);
    }
}
