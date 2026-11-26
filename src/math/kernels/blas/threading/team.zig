//! Level 3's workers, which wait on each other (`threading.zig`): the team,
//! its hold on the caller's control state, and the signal and barrier.

const std = @import("std");
const sysinfo = @import("builtin");
const FpState = @import("../FpState.zig");
const Qos = @import("Qos.zig");
const Tasks = @import("tasks.zig").Tasks;
const pool = @import("pool.zig");

const Io = std.Io;

/// Runs `job` on up to `want` workers that wait on each other, the caller
/// as worker 0: `job.plan(n)` sizes the work for the `n` that started and
/// returns how many take part, `job.work(io, i, lent)` runs worker `i`.
/// Workers hold the caller's control state (`Lent`); returns when all end.
pub fn team(io: Io, job: anytype, want: u32) Io.Cancelable!void {
    if (comptime sysinfo.single_threaded) {
        _ = job.plan(1);
        return job.work(io, 0, null);
    }
    const Job = @TypeOf(job.*);
    const S = struct {
        fn worker(j: *Job, go: *Signal, i: u32, fp: FpState, qos: Qos, wio: Io) void {
            qos.raise();
            var lent: Lent = .{ .want = fp };
            lent.take();
            defer lent.give();
            const active = go.waitWhile(wio, 0, &lent) catch return pool.recancel(wio);
            if (i < active) j.work(wio, i, &lent) catch pool.recancel(wio);
        }
    };
    const fp: FpState = .capture();
    const qos: Qos = .capture();
    // the caller's own: a fiber `Io` may resume it on another thread too
    var mine: Lent = .{ .want = fp };
    var go: Signal = .{};
    var tasks: Tasks = .init(io);
    defer {
        // first, canceling workers still waiting on `go`
        tasks.cancel(io);
        // the caller keeps its state where it resumed
        mine.take();
    }
    var started: u32 = 1;
    while (started < want) : (started += 1) tasks.concurrent(io, S.worker, .{ job, &go, started, fp, qos, io }) catch break;
    go.set(io, job.plan(started), &mine);
    try job.work(io, 0, &mine);
    try tasks.await(io);
}

/// A worker's hold on the caller's control state (`FpState`), in force on
/// its thread only while it runs: a fiber `Io` (`Dispatch`) may resume a
/// waiter on another thread, which must not keep it, nor the first its own.
/// The caller's starts given (its thread's own is `want`).
pub const Lent = struct {
    want: FpState,
    /// The thread's own state while `want` replaces it.
    prior: ?FpState = null,

    /// `want` in force on this thread.
    pub fn take(l: *Lent) void {
        l.prior = l.want.adopt();
    }

    /// The thread's own state back, before a wait that may move the worker.
    pub fn give(l: *Lent) void {
        if (l.prior) |p| p.apply();
        l.prior = null;
    }
};

/// Spins of a wait before it sleeps; its time per spin (`pause` on Zen 3)
/// is unmeasured, so the wait's length is too. BLIS and OpenBLAS spin without
/// bound; the bound lets a fiber `Io` run the peer a waiter waits on, and a
/// canceled waiter reach its futex.
pub const spin_limit = 1 << 12;

/// A word workers wait on to change. They spin first, as BLIS's
/// `bli_thrcomm_barrier_atomic` and OpenBLAS's `level3_thread.c` do, then
/// sleep on a futex; setters wake only when someone sleeps. A worker sleeps
/// and wakes others without its `Lent` state, and takes it again after.
pub const Signal = struct {
    value: std.atomic.Value(u32) = .init(0),
    sleepers: std.atomic.Value(u32) = .init(0),

    /// Waits until the value is other than `old`; returns it.
    pub fn waitWhile(s: *Signal, io: Io, old: u32, lent: ?*Lent) Io.Cancelable!u32 {
        for (0..spin_limit) |_| {
            const x = s.value.load(.acquire);
            if (x != old) return x;
            std.atomic.spinLoopHint();
        }
        while (true) {
            // with `set`'s order, either it sees this sleeper or this sees its value
            _ = s.sleepers.fetchAdd(1, .seq_cst);
            defer _ = s.sleepers.fetchSub(1, .monotonic);
            const x = s.value.load(.seq_cst);
            if (x != old) return x;
            if (lent) |l| l.give();
            defer if (lent) |l| l.take();
            // a pool worker looks for its call's cancel now and then
            if (pool.onPool()) {
                if (pool.aborted()) return error.Canceled;
                try io.futexWaitTimeout(u32, &s.value.raw, old, pool.abort_poll);
            } else try io.futexWait(u32, &s.value.raw, old);
        }
    }

    /// Sets the value; a wake may move the setter too (`Dispatch`).
    pub fn set(s: *Signal, io: Io, v: u32, lent: ?*Lent) void {
        s.value.store(v, .seq_cst);
        if (s.sleepers.load(.seq_cst) == 0) return;
        if (lent) |l| l.give();
        defer if (lent) |l| l.take();
        // all: `Uring` passes the count on as a C `int`, where `maxInt(u32)`
        // wakes one (issue #207)
        io.futexWake(u32, &s.value.raw, std.math.maxInt(i32));
    }
};

/// A reusable sense-reversing barrier of `parties` workers
/// (`bli_thrcomm_barrier_atomic`), waiting as `Signal` does.
pub const Barrier = struct {
    parties: u32,
    arrived: std.atomic.Value(u32) align(std.atomic.cache_line) = .init(0),
    phase: Signal align(std.atomic.cache_line) = .{},

    pub fn wait(b: *Barrier, io: Io, lent: ?*Lent) Io.Cancelable!void {
        if (b.parties <= 1) return;
        const p = b.phase.value.load(.monotonic);
        if (b.arrived.fetchAdd(1, .acq_rel) + 1 == b.parties) {
            b.arrived.store(0, .monotonic);
            b.phase.set(io, p +% 1, lent);
            return;
        }
        _ = try b.phase.waitWhile(io, p, lent);
    }
};
