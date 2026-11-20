//! A call's fixed chunks (`threading.zig`): a range split by its length, the
//! count asked for, shape and quantum only.

const std = @import("std");

/// The most chunks per call; a reduction keeps one partial each on the stack.
pub const chunks_max = 64;

/// How cost spreads over a range's units.
pub const Shape = enum {
    /// Every unit costs the same.
    even,
    /// Unit `i` costs about `i`, as a triangle's growing rows.
    rising,
    /// Unit `i` costs about `len - i`.
    falling,
};

/// Units `start..end` of a range, chunk `index` of its split.
pub const Chunk = struct { index: usize, start: usize, end: usize };

/// A range split into balanced, non-empty chunks by its length, the chunk
/// count asked for, shape and quantum only.
pub const Chunks = struct {
    count: usize,
    bounds: [chunks_max + 1]usize,

    /// `want` chunks (clamped to 1 and `chunks_max`) of `len` units, inner
    /// bounds rounded down to multiples of `quantum`; empty ones dropped. An
    /// empty range is one empty chunk.
    pub fn init(len: usize, want: usize, shape: Shape, quantum: usize) Chunks {
        const q = @max(quantum, 1);
        const c = std.math.clamp(want, 1, chunks_max);
        var r: Chunks = .{ .count = 0, .bounds = undefined };
        r.bounds[0] = 0;
        for (1..c + 1) |k| {
            const b = if (k == c) len else blk: {
                const x = cut(len, c, k, shape);
                break :blk x - x % q;
            };
            if (b > r.bounds[r.count]) {
                r.count += 1;
                r.bounds[r.count] = b;
            }
        }
        if (r.count == 0) {
            r.count = 1;
            r.bounds[1] = 0;
        }
        return r;
    }

    pub fn get(r: *const Chunks, i: usize) Chunk {
        return .{ .index = i, .start = r.bounds[i], .end = r.bounds[i + 1] };
    }

    /// Bound `k` of `c` before rounding: `k/c` of the cost.
    fn cut(len: usize, c: usize, k: usize, shape: Shape) usize {
        return switch (shape) {
            .even => @intCast(fraction(len, c, k)),
            .rising => @intCast(std.math.sqrt(fraction(@as(u128, len) * len, c, k))),
            .falling => len - @as(usize, @intCast(std.math.sqrt(fraction(@as(u128, len) * len, c, c - k)))),
        };
    }

    /// `floor(x * k / c)` for `k <= c`, without overflow.
    fn fraction(x: u128, c: usize, k: usize) u128 {
        return (x / c) * k + (x % c) * k / c;
    }
};
