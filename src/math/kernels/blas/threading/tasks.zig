//! A region's worker tasks (`threading.zig`): an `Io.Group`, or plain futures
//! on `Io.Dispatch`.

const std = @import("std");
const sysinfo = @import("builtin");
const Signal = @import("team.zig").Signal;

const Io = std.Io;
const chunks_max = @import("chunks.zig").chunks_max;

/// Plain futures instead of `Io.Group` on `Io.Dispatch`, whose `Group` crashes
/// with worker threads (issue #208). Once `ssh-repro-group-aarch64-macos`
/// passes, set this false and delete `Futures`.
pub const dispatch_futures = true;

/// `Io.Evented` for the fiber tests; void where it does not build: Zig 0.17's
/// `Uring`, and its `Dispatch`, whose vtable names a field 0.17's `Io.VTable`
/// lacks (version-gated, std exposes nothing to detect); SVE targets, whose
/// LLVM rejects the fiber switch's `ffr` clobber (issue #212).
pub const Evented = if ((Io.Evented == Io.Uring or Io.Evented == Io.Dispatch) and sysinfo.zig_version.major == 0 and sysinfo.zig_version.minor == 17) void else if (sysinfo.cpu.arch.isAarch64() and std.Target.aarch64.featureSetHas(sysinfo.cpu.features, .sve)) void else Io.Evented;

/// Whether `io`'s worker tasks run on `Futures`.
pub fn usesFutures(io: Io) bool {
    // `Evented`, not `Io.Evented`: Zig 0.17.0 cannot instantiate `Dispatch`
    if (comptime !dispatch_futures or Evented != Io.Dispatch) return false;
    const S = struct {
        var ev: Io.Dispatch = undefined;
    };
    // `Dispatch.io` returns one constant vtable
    return io.vtable == S.ev.io().vtable;
}

/// A region's worker tasks: an `Io.Group`, or `Futures` (`usesFutures`).
pub const Tasks = struct {
    group: Io.Group = .init,
    futures: Futures = .{},
    on_futures: bool,

    pub fn init(io: Io) Tasks {
        return .{ .on_futures = usesFutures(io) };
    }

    pub fn async(t: *Tasks, io: Io, function: anytype, args: std.meta.ArgsTuple(@TypeOf(function))) void {
        if (t.on_futures) return t.futures.async(io, function, args);
        t.group.async(io, function, args);
    }

    pub fn concurrent(t: *Tasks, io: Io, function: anytype, args: std.meta.ArgsTuple(@TypeOf(function))) Io.ConcurrentError!void {
        if (t.on_futures) return t.futures.concurrent(io, function, args);
        return t.group.concurrent(io, function, args);
    }

    pub fn await(t: *Tasks, io: Io) Io.Cancelable!void {
        if (t.on_futures) return t.futures.await(io);
        return t.group.await(io);
    }

    pub fn cancel(t: *Tasks, io: Io) void {
        if (t.on_futures) return t.futures.cancel(io);
        t.group.cancel(io);
    }
};

/// Tasks as plain futures, awaited and canceled as `Io.Group`'s are: `await`
/// waits cancelably for them all, then `cancel` (the caller's `defer`)
/// cancels and awaits each.
const Futures = struct {
    list: [cap]Io.Future(void) = undefined,
    len: usize = 0,
    /// Tasks not yet finished, plus one the caller holds until `await`.
    left: std.atomic.Value(u32) = .init(1),
    done: Signal = .{},

    const cap = if (dispatch_futures and Evented == Io.Dispatch) chunks_max else 0;

    fn Task(comptime function: anytype) type {
        return struct {
            fn run(fs: *Futures, io: Io, args: std.meta.ArgsTuple(@TypeOf(function))) void {
                @call(.auto, function, args);
                if (fs.left.fetchSub(1, .acq_rel) == 1) fs.done.set(io, 1, null);
            }
        };
    }

    /// Past `cap`, inline, as `Group.async` may run it.
    fn async(fs: *Futures, io: Io, function: anytype, args: std.meta.ArgsTuple(@TypeOf(function))) void {
        if (cap == 0 or fs.len == cap) return @call(.auto, function, args);
        _ = fs.left.fetchAdd(1, .monotonic);
        fs.list[fs.len] = io.async(Task(function).run, .{ fs, io, args });
        fs.len += 1;
    }

    fn concurrent(fs: *Futures, io: Io, function: anytype, args: std.meta.ArgsTuple(@TypeOf(function))) Io.ConcurrentError!void {
        if (cap == 0 or fs.len == cap) return error.ConcurrencyUnavailable;
        _ = fs.left.fetchAdd(1, .monotonic);
        fs.list[fs.len] = io.concurrent(Task(function).run, .{ fs, io, args }) catch |err| {
            _ = fs.left.fetchSub(1, .monotonic);
            return err;
        };
        fs.len += 1;
    }

    fn await(fs: *Futures, io: Io) Io.Cancelable!void {
        if (fs.left.fetchSub(1, .acq_rel) != 1) _ = try fs.done.waitWhile(io, 0, null);
        for (fs.list[0..fs.len]) |*f| f.await(io);
        fs.len = 0;
    }

    fn cancel(fs: *Futures, io: Io) void {
        for (fs.list[0..fs.len]) |*f| f.cancel(io);
        fs.len = 0;
    }
};
