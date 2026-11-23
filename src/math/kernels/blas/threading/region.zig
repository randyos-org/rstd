//! The claim loop and the ordered combine (`threading.zig`).

const std = @import("std");
const sysinfo = @import("builtin");
const FpState = @import("../FpState.zig");
const chunking = @import("chunks.zig");
const Plan = @import("../threading.zig").Plan;
const Tasks = @import("tasks.zig").Tasks;

const Allocator = std.mem.Allocator;
const Io = std.Io;
const chunks_max = chunking.chunks_max;
const Chunk = chunking.Chunk;
const Chunks = chunking.Chunks;

/// `f(ctx, chunk)` for every chunk of `p`, on up to `p.workers` workers that
/// claim chunks in turn; the caller is worker 0, and the rest start through
/// `Group.async` (`Tasks`), so an `Io` that admits fewer runs them inline,
/// and adopt the caller's floating-point control state (`FpState.zig`).
/// Returns once every task has finished. `error.Canceled` leaves chunks undone.
pub fn forEach(io: Io, p: *const Plan, ctx: anytype, comptime f: fn (@TypeOf(ctx), Chunk) void) Io.Cancelable!void {
    return forEachIn(io, p, ctx, f, null);
}

/// `forEach` for a reduction: `partial(ctx, chunk)` per chunk, folded in chunk
/// order on the caller, so the result does not depend on who computed what.
pub fn reduce(
    io: Io,
    p: *const Plan,
    comptime R: type,
    ctx: anytype,
    comptime partial: fn (@TypeOf(ctx), Chunk) R,
    comptime fold: fn (R, R) R,
) Io.Cancelable!R {
    return reduceIn(io, p, R, ctx, partial, fold, null);
}

/// `forEach` with scratch: one allocation, cut into `len(ctx, chunk)`
/// elements per chunk in chunk order; `f` gets its chunk's part.
pub fn forEachPart(
    allocator: Allocator,
    io: Io,
    p: *const Plan,
    comptime E: type,
    comptime alignment: usize,
    ctx: anytype,
    comptime len: fn (@TypeOf(ctx), Chunk) usize,
    comptime f: fn (@TypeOf(ctx), Chunk, []E) void,
) (Allocator.Error || Io.Cancelable)!void {
    var at: [chunks_max + 1]usize = undefined;
    at[0] = 0;
    for (0..p.chunks.count) |i| at[i + 1] = at[i] + len(ctx, p.chunks.get(i));
    const buf = try allocator.alignedAlloc(E, .fromByteUnits(alignment), at[p.chunks.count]);
    defer allocator.free(buf);
    const Parts = struct {
        ctx: @TypeOf(ctx),
        buf: [*]E,
        at: *const [chunks_max + 1]usize,

        fn one(s: @This(), c: Chunk) void {
            f(s.ctx, c, s.buf[s.at[c.index]..s.at[c.index + 1]]);
        }
    };
    return forEach(io, p, Parts{ .ctx = ctx, .buf = buf.ptr, .at = &at }, Parts.one);
}

/// `forEach`, claiming chunks in `order` if given (tests).
pub fn forEachIn(io: Io, p: *const Plan, ctx: anytype, comptime f: fn (@TypeOf(ctx), Chunk) void, order: ?[]const u8) Io.Cancelable!void {
    var region: Region(@TypeOf(ctx), f) = .{ .ctx = ctx, .chunks = &p.chunks, .order = order };
    if (comptime !sysinfo.single_threaded) {
        const w = @min(p.workers, p.chunks.count);
        if (w > 1) return region.parallel(io, w);
    }
    return region.claim(io);
}

pub fn reduceIn(
    io: Io,
    p: *const Plan,
    comptime R: type,
    ctx: anytype,
    comptime partial: fn (@TypeOf(ctx), Chunk) R,
    comptime fold: fn (R, R) R,
    order: ?[]const u8,
) Io.Cancelable!R {
    const Into = struct {
        ctx: @TypeOf(ctx),
        partials: *[chunks_max]R,

        fn one(into: @This(), c: Chunk) void {
            into.partials[c.index] = partial(into.ctx, c);
        }
    };
    var partials: [chunks_max]R = undefined;
    try forEachIn(io, p, Into{ .ctx = ctx, .partials = &partials }, Into.one, order);
    var acc = partials[0];
    for (partials[1..p.chunks.count]) |x| acc = fold(acc, x);
    return acc;
}

/// One parallel region's shared state; tasks get a pointer to it.
fn Region(comptime Ctx: type, comptime f: fn (Ctx, Chunk) void) type {
    return struct {
        const Self = @This();

        ctx: Ctx,
        chunks: *const Chunks,
        order: ?[]const u8,
        next: std.atomic.Value(usize) = .init(0),
        /// The caller's, which `work` adopts.
        fp: FpState = undefined,

        fn claim(r: *Self, io: Io) Io.Cancelable!void {
            while (true) {
                try checkCancel(io);
                const k = r.next.fetchAdd(1, .monotonic);
                if (k >= r.chunks.count) return;
                f(r.ctx, r.chunks.get(if (r.order) |o| o[k] else k));
            }
        }

        /// A worker task. It re-arms a cancellation it took: a task run inline
        /// shares the caller's, and `Group` drops the error.
        fn work(r: *Self, io: Io) void {
            const own = r.fp.adopt();
            defer if (own) |s| s.apply();
            r.claim(io) catch |err| switch (err) {
                error.Canceled => io.recancel(),
            };
        }

        fn parallel(r: *Self, io: Io, w: usize) Io.Cancelable!void {
            r.fp = .capture();
            var tasks: Tasks = .init(io);
            // runs before the caller's frame goes; every task has then run
            defer {
                tasks.cancel(io);
                // the caller's state where a fiber `Io` resumed it (`Lent`)
                _ = r.fp.adopt();
            }
            // `concurrent`, not `async`: a helper run inline would claim every
            // chunk before the next started (issue #685)
            for (1..w) |_| tasks.concurrent(io, Self.work, .{ r, io }) catch break;
            try r.claim(io);
            try tasks.await(io);
        }
    };
}

/// `io.checkCancel()`, or nothing on an `Io` without cancellation, whose
/// `checkCancel` is unreachable (`Io.failing`, rstd's UEFI `Io`).
pub fn checkCancel(io: Io) Io.Cancelable!void {
    if (io.vtable.checkCancel == &Io.unreachableCheckCancel) return;
    return io.checkCancel();
}
