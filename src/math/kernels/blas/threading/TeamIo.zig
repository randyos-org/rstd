//! A caller-owned worker team (issue #561): `init` starts its workers once,
//! on the caller's `Io`, and they wait between jobs; `deinit` ends them. Its
//! `io()` is that `Io` again, so it goes wherever an `io` does; rstd's
//! threading (`region.zig`) finds the team behind it and runs a region on
//! the waiting workers, where any other `Io` starts tasks as before. Bits are
//! as ever the plan's: chunks and folds do not know who ran them.

const std = @import("std");
const sysinfo = @import("builtin");
const Signal = @import("team.zig").Signal;
const Tasks = @import("tasks.zig").Tasks;

const Io = std.Io;
const Team = @This();

/// The team's `Io` vtable: the caller's, with `crashHandler` (called only on
/// `std.Options.debug_io`) as the marker `of` looks for.
vtable: Io.VTable,
inner: Io,
tasks: Tasks,
/// Workers started, the caller not among them.
started: u32 = 0,
/// A word of its own for each worker, bumped to wake it for a job or to stop:
/// the caller wakes the participants only, and each worker sees each of its
/// own wakes, so none runs a job meant for the next generation.
wake: [workers_max]Wake = @splat(.{}),
/// Participants that have finished the job.
fin: Signal = .{},
stop: std.atomic.Value(bool) = .init(false),
busy: std.atomic.Value(bool) = .init(false),
/// The job: its function and context.
run_fn: ?*const fn (*anyopaque, Io, u32) void = null,
run_ctx: ?*anyopaque = null,

/// The most workers a team starts, the caller not among them.
pub const workers_max = 64;

const Wake = struct { signal: Signal align(std.atomic.cache_line) = .{} };

fn marker(_: ?*anyopaque) void {}

/// Starts up to `workers - 1` workers on `parent`, as many as it admits. `t`
/// must stay where it is until `deinit`.
pub fn init(t: *Team, parent: Io, workers: u32) void {
    t.* = .{ .vtable = parent.vtable.*, .inner = parent, .tasks = .init(parent) };
    t.vtable.crashHandler = &marker;
    if (comptime sysinfo.single_threaded) return;
    while (t.started < workers_max and t.started + 1 < workers) : (t.started += 1) {
        t.tasks.concurrent(parent, worker, .{ t, t.started + 1 }) catch break;
    }
}

/// Ends the workers; no job may be running.
pub fn deinit(t: *Team) void {
    t.stop.store(true, .release);
    for (t.wake[0..t.started]) |*w| w.signal.set(t.inner, w.signal.value.load(.monotonic) +% 1, null);
    t.tasks.await(t.inner) catch t.tasks.cancel(t.inner);
    t.* = undefined;
}

/// The `Io` of this team: the caller's own functions on the caller's data.
pub fn io(t: *Team) Io {
    return .{ .userdata = t.inner.userdata, .vtable = &t.vtable };
}

/// The team behind `i`, if it is a team's `Io`.
pub fn of(i: Io) ?*Team {
    if (i.vtable.crashHandler != &marker) return null;
    return @alignCast(@fieldParentPtr("vtable", @constCast(i.vtable)));
}

/// Workers and caller: the most a job can use.
pub fn size(t: *const Team) u32 {
    return t.started + 1;
}

fn worker(t: *Team, i: u32) void {
    const wio = t.inner;
    var seen: u32 = 0;
    while (true) {
        seen = t.wake[i - 1].signal.waitWhile(wio, seen, null) catch return wio.recancel();
        if (t.stop.load(.acquire)) return;
        // woken only for a job it takes part in, and the caller waits for it
        t.run_fn.?(t.run_ctx.?, wio, i);
        _ = t.fin.value.fetchAdd(1, .seq_cst);
        if (t.fin.sleepers.load(.seq_cst) != 0) wio.futexWake(u32, &t.fin.value.raw, std.math.maxInt(i32));
    }
}

/// Takes the team for one job of up to `want` participants, the caller one;
/// how many may take part, or null where the team is busy (a job inside a
/// job) or has too few workers to matter. `end` gives it back.
pub fn begin(t: *Team, want: u32) ?u32 {
    const take = @min(want, t.size());
    if (take < 2) return null;
    if (t.busy.swap(true, .acquire)) return null;
    return take;
}

pub fn end(t: *Team) void {
    t.busy.store(false, .release);
}

/// Runs `f(ctx, io, i)` on workers 1 to `take - 1` (`take` from `begin`, or
/// fewer) and `caller(ctx)` here, returning once every participant is done,
/// also after a cancel, whose chunks the workers then finish.
pub fn go(
    t: *Team,
    take: u32,
    ctx: *anyopaque,
    comptime f: fn (*anyopaque, Io, u32) void,
    comptime caller: fn (*anyopaque) Io.Cancelable!void,
) Io.Cancelable!void {
    if (take < 2) return caller(ctx);
    const need = take - 1;
    t.run_fn = f;
    t.run_ctx = ctx;
    t.fin.value.store(0, .release);
    for (t.wake[0..need]) |*w| w.signal.set(t.inner, w.signal.value.load(.monotonic) +% 1, null);
    const first = caller(ctx);
    // the caller's own wait may be canceled; the workers' chunks still end
    var canceled = false;
    var v = t.fin.value.load(.acquire);
    while (v != need) {
        v = if (canceled) blk: {
            std.atomic.spinLoopHint();
            break :blk t.fin.value.load(.acquire);
        } else t.fin.waitWhile(t.inner, v, null) catch blk: {
            canceled = true;
            break :blk t.fin.value.load(.acquire);
        };
    }
    if (canceled) return error.Canceled;
    return first;
}

/// `begin`, `go` and `end`; null where `begin` is.
pub fn run(
    t: *Team,
    want: u32,
    ctx: *anyopaque,
    comptime f: fn (*anyopaque, Io, u32) void,
    comptime caller: fn (*anyopaque) Io.Cancelable!void,
) ?Io.Cancelable!void {
    const take = t.begin(want) orelse return null;
    defer t.end();
    return t.go(take, ctx, f, caller);
}
