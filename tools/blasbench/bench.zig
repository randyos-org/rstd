//! Times calls against each other: `race` grows each one's batch past the clock's cost, then takes the
//! fastest of several batches, the contestants' batches in turn (`blas_bench`'s sections).

const std = @import("std");
const Io = std.Io;

/// How a race is timed: each batch runs at least `min_ns`, and the fastest
/// of `batches` counts.
pub const Timing = struct {
    min_ns: i96,
    batches: usize,
};

/// One contestant: `run(ctx)` makes one call.
pub const Call = struct {
    ctx: *anyopaque,
    run: *const fn (*anyopaque) void,

    /// `ctx`'s call: its type's `run`.
    pub fn of(ctx: anytype) Call {
        return .{ .ctx = ctx, .run = @typeInfo(@TypeOf(ctx)).pointer.child.run };
    }
};

/// The most contestants one race takes.
pub const max_calls = 16;

/// Seconds per call of each of `calls` into `out`, from each one's fastest
/// batch; the batches go round the contestants in turn.
pub fn race(io: Io, t: Timing, calls: []const Call, out: []f64) error{TooManyCalls}!void {
    if (calls.len > max_calls or out.len < calls.len) return error.TooManyCalls;
    var reps: [max_calls]usize = undefined;
    for (calls, 0..) |c, i| {
        c.run(c.ctx); // warm up
        reps[i] = 1;
        while (batch(io, c, reps[i]) < t.min_ns) reps[i] *= 2;
        out[i] = std.math.inf(f64);
    }
    for (0..t.batches) |_| {
        for (calls, 0..) |c, i| {
            const s = @as(f64, @floatFromInt(batch(io, c, reps[i]))) / 1e9 / @as(f64, @floatFromInt(reps[i]));
            out[i] = @min(out[i], s);
        }
    }
}

/// Nanoseconds for `reps` calls of `c`.
fn batch(io: Io, c: Call, reps: usize) i96 {
    const t0 = nanos(io);
    for (0..reps) |_| c.run(c.ctx);
    return nanos(io) - t0;
}

fn nanos(io: Io) i96 {
    return Io.Clock.awake.now(io).toNanoseconds();
}
