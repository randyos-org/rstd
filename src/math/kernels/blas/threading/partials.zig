//! Per-chunk partial vectors (`threading.zig`), for level 2's folded products.

const std = @import("std");
const types = @import("../types.zig");
const util = @import("../util.zig");
const threading = @import("../threading.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const chunks_max = threading.chunks_max;
const Chunk = threading.Chunk;
const Plan = threading.Plan;
const plan = threading.plan;
const forEach = threading.forEach;

/// Rows `lo..hi` of a vector.
pub const Rows = struct { lo: usize, hi: usize };

/// Per-chunk partial vectors, for a level-2 product whose elements of `y` sum
/// over several chunks' columns: chunk `i` adds its part into a zeroed vector
/// of its own over the rows it reaches, and `fold` adds those to `y` in chunk
/// order. Chunks follow the size alone, so the bits do too, at any worker
/// count (OpenBLAS's `symv_thread.c` keeps one vector per thread instead).
pub fn Partials(comptime T: type) type {
    return struct {
        const Self = @This();
        /// Each vector starts on a cache line.
        const line = 64 / std.math.gcd(64, @sizeOf(T));

        buf: []align(64) T,
        rows: [chunks_max]Rows,
        at: [chunks_max]usize,
        count: usize,
        /// The rows `fold` reaches: the highest `hi`.
        len: usize,

        /// Room for `p`'s chunks, chunk `c` over `rows(ctx, c)`.
        pub fn init(allocator: Allocator, p: *const Plan, ctx: anytype, comptime rows: fn (@TypeOf(ctx), Chunk) Rows) Allocator.Error!Self {
            var s: Self = .{ .buf = &.{}, .rows = undefined, .at = undefined, .count = p.chunks.count, .len = 0 };
            var total: usize = 0;
            for (0..s.count) |i| {
                const r = rows(ctx, p.chunks.get(i));
                s.rows[i] = r;
                s.at[i] = total;
                s.len = @max(s.len, r.hi);
                total += std.mem.alignForward(usize, r.hi - r.lo, line);
            }
            s.buf = try allocator.alignedAlloc(T, .@"64", total);
            return s;
        }

        pub fn deinit(s: *const Self, allocator: Allocator) void {
            allocator.free(s.buf);
        }

        /// Chunk `i`'s vector, zeroed; element `r - rows[i].lo` is row `r`.
        pub fn take(s: *const Self, i: usize) []T {
            const v = s.buf[s.at[i]..][0 .. s.rows[i].hi - s.rows[i].lo];
            @memset(v, std.mem.zeroes(T));
            return v;
        }

        /// `y[r * incy] += ` each chunk's row `r`, chunks in order, each sum
        /// rounded to its type; rows split over workers, which cannot change
        /// a bit.
        pub fn fold(s: *const Self, io: Io, y: [*]T, incy: isize) Io.Cancelable!void {
            var work: u64 = 0;
            for (s.rows[0..s.count]) |r| work += r.hi - r.lo;
            // split as `axpy` is
            const p = plan(.axpy, T, work, s.len, .even, line);
            return forEach(io, &p, Fold{ .s = s, .y = y, .incy = incy }, Fold.rows);
        }

        const Fold = struct {
            s: *const Self,
            y: [*]T,
            incy: isize,

            fn rows(f: Fold, c: Chunk) void {
                for (f.s.rows[0..f.s.count], f.s.at[0..f.s.count]) |r, at| {
                    const lo = @max(r.lo, c.start);
                    const hi = @min(r.hi, c.end);
                    if (lo >= hi) continue;
                    for (lo..hi) |row| {
                        const d = util.at(T, f.y, @as(i64, @intCast(row)) * f.incy);
                        d.* = add(d.*, f.s.buf[at + row - r.lo]);
                    }
                }
            }
        };

        fn add(a: T, b: T) T {
            if (comptime types.isComplex(T)) {
                const F = types.Float(T);
                return .init(util.rounded(F, a.re + b.re), util.rounded(F, a.im + b.im));
            }
            return util.rounded(T, a + b);
        }
    };
}
