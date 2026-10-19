//! Runs a rig's qemu with its output passed through, and ends it when the
//! guest can no longer finish: firmware that halts rather than exits (OpenBIOS
//! "Stopping execution"), or no output at all for too long.
//!
//!   qemu_watch <idle minutes> <halt text>... -- <qemu argv>
//!
//! Exits with qemu's status, or 1 after ending it.

const std = @import("std");
const Io = std.Io;

const State = struct {
    /// When qemu last wrote, in seconds of the awake clock.
    last: std.atomic.Value(i64),
    halted: std.atomic.Value(bool) = .init(false),
    done: std.atomic.Value(bool) = .init(false),
};

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
    const halts = args[2..dashes];

    var qemu = try std.process.spawn(io, .{ .argv = args[dashes + 1 ..], .stdin = .ignore, .stdout = .pipe, .stderr = .inherit });
    var state: State = .{ .last = .init(now(io)) };
    var pump = try io.concurrent(passThrough, .{ io, qemu.stdout.?, halts, &state });

    while (!state.done.load(.acquire)) {
        try io.sleep(.fromSeconds(2), .awake);
        const why: ?[]const u8 = if (state.halted.load(.acquire))
            "the firmware halted"
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
    pump.await(io);
    return switch (try qemu.wait(io)) {
        .exited => |code| code,
        else => 1,
    };
}

fn now(io: Io) i64 {
    return Io.Clock.awake.now(io).toSeconds();
}

/// Copies qemu's output to this stdout, noting when it last wrote and any halt
/// text, which a chunk boundary may split, hence the kept tail.
fn passThrough(io: Io, from: Io.File, halts: []const []const u8, state: *State) void {
    defer state.done.store(true, .release);
    const out = Io.File.stdout();
    var buf: [4096 + 256]u8 = undefined;
    var kept: usize = 0;
    while (true) {
        const n = from.readStreaming(io, &.{buf[kept..]}) catch return;
        if (n == 0) continue;
        out.writeStreamingAll(io, buf[kept..][0..n]) catch {};
        state.last.store(now(io), .release);
        const seen = buf[0 .. kept + n];
        for (halts) |h| if (std.mem.indexOf(u8, seen, h) != null) state.halted.store(true, .release);
        kept = @min(seen.len, 256);
        std.mem.copyForwards(u8, buf[0..kept], seen[seen.len - kept ..]);
    }
}
