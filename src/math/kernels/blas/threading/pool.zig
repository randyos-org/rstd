//! Worker threads that spin between calls, as OpenBLAS's
//! `driver/others/blas_server.c` (#710): started as calls first need them,
//! each waits on its own slot, spinning `2^OPENBLAS_THREAD_TIMEOUT` ticks of
//! `rpcc` before it sleeps. A call's `Crew` hands its tasks to idle workers
//! and waits for them. Used on `Io.Threaded`, for x86_64 and arm64 with OS
//! threads; elsewhere teams start through the `Io` (`tasks.zig`).

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../__root__.zig");

const Io = std.Io;
const chunks_max = @import("chunks.zig").chunks_max;

/// Whether this target has the pool: OS threads and an `rpcc`.
pub const built = !sysinfo.single_threaded and switch (sysinfo.target.cpu.arch) {
    .x86_64, .aarch64 => true,
    else => false,
} and switch (sysinfo.target.os.tag) {
    .windows, .linux, .macos => true,
    else => false,
};

/// The most pool workers, over every call at once.
const workers_max = chunks_max;

/// The `Io.Threaded` behind `io`, if it is one.
fn threaded(io: Io) ?*Io.Threaded {
    if (io.vtable != Io.Threaded.global_single_threaded.io().vtable) return null;
    return @ptrCast(@alignCast(io.userdata));
}

/// The pool workers a call on `io` may take: its `concurrent_limit` on
/// `Io.Threaded`, which never changes after `init`; 0 on any other `Io`.
pub fn allowance(io: Io) u32 {
    if (comptime !built) return 0;
    const t = threaded(io) orelse return 0;
    const n = t.concurrent_limit.toInt() orelse return workers_max;
    return @intCast(@min(n, workers_max));
}

/// OpenBLAS's `rpcc`: `rdtsc` (`common_x86_64.h`); on arm64 the virtual
/// count shifted by `cntfrq`'s leading zeros, about 3 GHz (`common_arm64.h`).
fn rpcc() u64 {
    switch (sysinfo.target.cpu.arch) {
        .x86_64 => {
            var lo: u32 = undefined;
            var hi: u32 = undefined;
            asm volatile ("rdtsc"
                : [lo] "={eax}" (lo),
                  [hi] "={edx}" (hi),
            );
            return (@as(u64, hi) << 32) | lo;
        },
        .aarch64 => {
            const count = asm volatile ("isb\nmrs %[r], cntvct_el0"
                : [r] "=r" (-> u64),
            );
            const freq = asm volatile ("mrs %[r], cntfrq_el0"
                : [r] "=r" (-> u64),
            );
            return count << @intCast(@clz(@as(u32, @truncate(freq))));
        },
        else => unreachable,
    }
}

/// `common.h`'s `YIELDING` on x86_64 and ARMV8: eight `nop`s.
fn yielding() void {
    asm volatile ("nop\nnop\nnop\nnop\nnop\nnop\nnop\nnop");
}

/// Ticks a worker or a waiting caller spins before it sleeps:
/// `blas_server.c`'s `thread_timeout`.
fn timeout() u64 {
    return @as(u64, 1) << blas.config.thread_timeout;
}

/// Spins on `done` until it is true or the timeout passes; whether it is.
fn spinFor(ctx: anytype, comptime done: fn (@TypeOf(ctx)) bool) bool {
    const start = rpcc();
    const limit = timeout();
    while (!done(ctx)) {
        yielding();
        if (rpcc() -% start > limit) return done(ctx);
    }
    return true;
}

/// The pool's own futex calls: `Io.Threaded`'s are the OS's, from any thread.
fn sleepIo() Io {
    return Io.Threaded.global_single_threaded.io();
}

/// One pool worker's mailbox (`blas_server.c`'s `thread_status`).
const Slot = struct {
    /// The task handed over; null while idle.
    task: std.atomic.Value(?*Entry) align(std.atomic.cache_line) = .init(null),
    /// Bumped with each task: the word a sleeping worker waits on.
    seq: std.atomic.Value(u32) = .init(0),
    sleeping: std.atomic.Value(bool) = .init(false),
    /// Held by a call from `claim` until its task ends.
    claimed: std.atomic.Value(bool) = .init(false),

    fn give(s: *Slot, e: *Entry) void {
        // with `serve`'s order, either it sees the task or this sees it asleep
        s.task.store(e, .seq_cst);
        _ = s.seq.fetchAdd(1, .seq_cst);
        if (s.sleeping.load(.seq_cst)) sleepIo().futexWake(u32, &s.seq.raw, 1);
    }

    fn hasTask(s: *Slot) bool {
        return s.task.load(.acquire) != null;
    }

    /// The next task, spinning first, then asleep until one comes.
    fn next(s: *Slot) *Entry {
        while (!spinFor(s, hasTask)) {
            const seq = s.seq.load(.seq_cst);
            s.sleeping.store(true, .seq_cst);
            if (s.task.load(.seq_cst) == null) sleepIo().futexWaitUncancelable(u32, &s.seq.raw, seq);
            s.sleeping.store(false, .monotonic);
        }
        return s.task.load(.acquire).?;
    }
};

var slots: [workers_max]Slot = @splat(.{});
var spawned: std.atomic.Value(u32) = .init(0);
var spawning: std.atomic.Value(bool) = .init(false);

/// An idle worker, started if none is; null when the pool is full or a
/// thread does not start.
fn claim() ?*Slot {
    if (claimSpawned()) |s| return s;
    while (spawning.swap(true, .acquire)) std.atomic.spinLoopHint();
    defer spawning.store(false, .release);
    // one may have started meanwhile, or finished its task
    if (claimSpawned()) |s| return s;
    const n = spawned.load(.monotonic);
    if (n == workers_max) return null;
    const s = &slots[n];
    s.claimed.store(true, .monotonic);
    const thread = std.Thread.spawn(.{}, serve, .{s}) catch {
        s.claimed.store(false, .monotonic);
        return null;
    };
    thread.detach();
    spawned.store(n + 1, .release);
    return s;
}

fn claimSpawned() ?*Slot {
    for (slots[0..spawned.load(.acquire)]) |*s| {
        if (s.claimed.cmpxchgStrong(false, true, .acquire, .monotonic) == null) return s;
    }
    return null;
}

/// The call whose task this pool worker runs; null off the pool.
threadlocal var current: ?*const Crew = null;

/// A pool worker's loop (`blas_thread_server`).
fn serve(s: *Slot) void {
    while (true) {
        const e = s.next();
        const crew = e.crew;
        current = crew;
        e.run(e);
        current = null;
        s.task.store(null, .monotonic);
        s.claimed.store(false, .release);
        // the last touch: the caller may return once all have finished
        crew.finish();
    }
}

/// Whether this thread is a pool worker running a task.
pub fn onPool() bool {
    return built and current != null;
}

/// Whether this pool worker's call has stopped waiting for it (`Crew.cancel`).
pub fn aborted() bool {
    if (comptime !built) return false;
    const c = current orelse return false;
    return c.abort.load(.monotonic);
}

/// How long a pool worker sleeps between looks at `aborted`.
pub const abort_poll: Io.Timeout = .{ .duration = .{ .raw = .fromMilliseconds(1), .clock = .awake } };

/// `io.recancel()`, except on a pool worker, which the `Io` never canceled.
pub fn recancel(io: Io) void {
    if (onPool()) return;
    io.recancel();
}

/// The most bytes of a task's arguments, stored in its entry.
const args_max = 64;

/// A task handed to a worker; it lives in its `Crew`.
const Entry = struct {
    run: *const fn (*Entry) void,
    crew: *Crew,
    args: [args_max]u8 align(16),
};

/// One call's pool tasks, awaited and canceled as `Io.Group`'s are: `await`
/// waits cancelably for them all, then `cancel` (the caller's `defer`) stops
/// the rest and waits for them.
pub const Crew = struct {
    entries: [cap]Entry = undefined,
    len: u32 = 0,
    /// `allowance`.
    limit: u32 = 0,
    /// Tasks finished, with `waiting` set while the caller sleeps on it.
    done: std.atomic.Value(u32) = .init(0),
    abort: std.atomic.Value(bool) = .init(false),

    const cap = if (built) workers_max else 0;
    const waiting: u32 = 1 << 31;

    fn Task(comptime function: anytype) type {
        return struct {
            const Args = std.meta.ArgsTuple(@TypeOf(function));

            fn run(e: *Entry) void {
                const args: *const Args = @ptrCast(@alignCast(&e.args));
                @call(.auto, function, args.*);
            }
        };
    }

    /// `function(args)` on an idle pool worker.
    pub fn concurrent(c: *Crew, function: anytype, args: std.meta.ArgsTuple(@TypeOf(function))) Io.ConcurrentError!void {
        const A = @TypeOf(args);
        comptime std.debug.assert(@sizeOf(A) <= args_max and @alignOf(A) <= 16);
        if (c.len == @min(c.limit, cap)) return error.ConcurrencyUnavailable;
        const s = claim() orelse return error.ConcurrencyUnavailable;
        const e = &c.entries[c.len];
        e.run = Task(function).run;
        e.crew = c;
        @as(*A, @ptrCast(&e.args)).* = args;
        c.len += 1;
        s.give(e);
    }

    fn finish(c: *Crew) void {
        // a wake after the caller left reaches only an address, harmlessly
        if (c.done.fetchAdd(1, .release) & waiting != 0) sleepIo().futexWake(u32, &c.done.raw, 1);
    }

    fn allDone(c: *Crew) bool {
        return c.done.load(.acquire) & ~waiting == c.len;
    }

    /// Until every task has finished: spinning first, as `exec_blas`, then asleep.
    fn wait(c: *Crew, io: Io, comptime cancelable: bool) Io.Cancelable!void {
        if (spinFor(c, allDone)) return;
        while (true) {
            const d = c.done.load(.acquire);
            if (d & ~waiting == c.len) return;
            if (d & waiting == 0 and c.done.cmpxchgWeak(d, d | waiting, .acquire, .acquire) != null) continue;
            if (cancelable) try io.futexWait(u32, &c.done.raw, d | waiting) else io.futexWaitUncancelable(u32, &c.done.raw, d | waiting);
        }
    }

    pub fn await(c: *Crew, io: Io) Io.Cancelable!void {
        return c.wait(io, true);
    }

    pub fn cancel(c: *Crew, io: Io) void {
        if (c.len == 0) return;
        c.abort.store(true, .monotonic);
        c.wait(io, false) catch unreachable;
        c.len = 0;
        c.done.store(0, .monotonic);
        c.abort.store(false, .monotonic);
    }
};
