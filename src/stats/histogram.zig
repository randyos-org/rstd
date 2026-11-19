//! Histograms: counts over even bins or given edges.

const std = @import("std");
const util = @import("util.zig");

const Allocator = std.mem.Allocator;
const Elem = util.Elem;
const to = util.to;
const checkReal = util.checkReal;
const checkFloat = util.checkFloat;

/// `count` bins spread evenly over `xs`'s range, or the given sorted edges.
pub fn Bins(comptime F: type) type {
    return union(enum) { count: usize, edges: []const F };
}

/// `counts[k]` holds `edges[k] <= x < edges[k + 1]`; the last bin includes
/// its right edge, as MATLAB's `histcounts` and NumPy's `histogram`.
pub fn Histogram(comptime F: type) type {
    return struct {
        counts: []u64,
        edges: []F,

        pub fn deinit(self: *@This(), allocator: Allocator) void {
            allocator.free(self.counts);
            allocator.free(self.edges);
            self.* = undefined;
        }
    };
}

pub const HistogramError = error{ EmptyInput, OutOfMemory, ZeroBins, EdgesNotSorted };
pub const HistogramIntoError = error{ ZeroBins, EdgesNotSorted, SliceTooShort };

/// `.count` bins need a value in `xs` that is not NaN, for their range: its
/// least to its greatest, widened by 0.5 each way when they are equal, as
/// NumPy's.
pub fn histogram(allocator: Allocator, comptime F: type, xs: anytype, bins: Bins(F)) HistogramError!Histogram(F) {
    comptime checkFloat(F);
    comptime checkReal(Elem(@TypeOf(xs)));
    const edges: []F = switch (bins) {
        .count => |n| blk: {
            if (n == 0) return error.ZeroBins;
            var lo = std.math.inf(F);
            var hi = -std.math.inf(F);
            for (xs) |x0| {
                const x = to(F, x0);
                if (std.math.isNan(x)) continue;
                lo = @min(lo, x);
                hi = @max(hi, x);
            }
            if (lo > hi) return error.EmptyInput;
            if (lo == hi) {
                lo -= 0.5;
                hi += 0.5;
            }
            const e = try allocator.alloc(F, n + 1);
            // NumPy's `linspace`: `lo + i step`, the last edge `hi` exactly
            const step = (hi - lo) / @as(F, @floatFromInt(n));
            for (e[0..n], 0..) |*v, i| v.* = lo + @as(F, @floatFromInt(i)) * step;
            e[n] = hi;
            break :blk e;
        },
        .edges => |given| blk: {
            if (given.len < 2) return error.ZeroBins;
            try sortedEdges(F, given);
            break :blk try allocator.dupe(F, given);
        },
    };
    errdefer allocator.free(edges);
    const counts = try allocator.alloc(u64, edges.len - 1);
    @memset(counts, 0);
    histogramInto(F, counts, xs, edges) catch |e| switch (e) {
        // checked above
        error.ZeroBins, error.EdgesNotSorted, error.SliceTooShort => unreachable,
    };
    return .{ .counts = counts, .edges = edges };
}

/// Adds to `counts[0 .. edges.len - 1]`; values outside the edges, and NaNs,
/// are not counted.
pub fn histogramInto(comptime F: type, counts: []u64, xs: anytype, edges: []const F) HistogramIntoError!void {
    comptime checkFloat(F);
    comptime checkReal(Elem(@TypeOf(xs)));
    if (edges.len < 2) return error.ZeroBins;
    if (counts.len < edges.len - 1) return error.SliceTooShort;
    try sortedEdges(F, edges);
    const last = edges.len - 1;
    for (xs) |x0| {
        const x = to(F, x0);
        if (!(x >= edges[0] and x <= edges[last])) continue;
        // the last `k` with `edges[k] <= x`, the right edge in the last bin
        var lo: usize = 0;
        var hi: usize = last;
        while (hi - lo > 1) {
            const mid = lo + (hi - lo) / 2;
            if (edges[mid] <= x) lo = mid else hi = mid;
        }
        counts[lo] += 1;
    }
}

/// Edges ascend (equal neighbors allowed, as NumPy's), and none is NaN.
fn sortedEdges(comptime F: type, edges: []const F) error{EdgesNotSorted}!void {
    for (edges[0 .. edges.len - 1], edges[1..]) |a, b| if (!(a <= b)) return error.EdgesNotSorted;
}
