//! std's `Io.Group` on `Io.Evented` with worker threads, std only
//! (`src/linalg/issues.md`, "`Group.await` blocks inside a context
//! switch"): three trivial tasks started and awaited, round after round. On
//! 0.18.0-dev.751 `group` crashes within a few thousand rounds on `Dispatch`
//! and `Uring`; `futures`, the same loop on plain futures, passes on
//! `Dispatch`. Rerun on each toolchain bump (`zig build
//! ssh-repro-group-aarch64-macos`); once `group` passes, set
//! `threading.dispatch_futures` false and delete its futures path.
//!
//!     evented_group_repro group|futures [rounds]

const std = @import("std");
const Io = std.Io;

pub fn main(init: std.process.Init) !void {
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len < 2) usage();
    const futures = if (std.mem.eql(u8, args[1], "futures")) true else if (std.mem.eql(u8, args[1], "group")) false else usage();
    const rounds = if (args.len > 2) try std.fmt.parseInt(usize, args[2], 10) else 20000;

    var ev: Io.Evented = undefined;
    var opts: Io.Evented.InitOptions = .{};
    // `Uring` starts no worker threads by default; `Dispatch` sizes its own pool
    if (@hasField(Io.Evented.InitOptions, "thread_limit")) opts.thread_limit = 3;
    try ev.init(std.heap.page_allocator, opts);
    defer ev.deinit();
    const io = ev.io();

    var counter: std.atomic.Value(u32) = .init(0);
    for (0..rounds) |round| {
        if (round % 2000 == 0) std.debug.print("{s} {d}\n", .{ args[1], round });
        if (futures) {
            var list: [3]Io.Future(void) = undefined;
            for (&list) |*f| f.* = try io.concurrent(task, .{&counter});
            for (&list) |*f| f.await(io);
        } else {
            var group: Io.Group = .init;
            for (0..3) |_| try group.concurrent(io, task, .{&counter});
            try group.await(io);
        }
    }
    std.debug.print("{s}: {d} rounds passed\n", .{ args[1], rounds });
}

fn task(counter: *std.atomic.Value(u32)) void {
    _ = counter.fetchAdd(1, .monotonic);
}

fn usage() noreturn {
    std.debug.print("usage: evented_group_repro group|futures [rounds]\n", .{});
    std.process.exit(2);
}
