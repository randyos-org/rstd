//! Histogram counts over sorted edges, NumPy's half-open bins.

const std = @import("std");
const strided = @import("../blas/strided.zig");

const Allocator = std.mem.Allocator;
const ConstStrided = strided.ConstStrided;

/// `count` bins spread evenly over the data's range, or the given sorted edges.
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

pub const HistogramIntoError = error{ ZeroBins, EdgesNotSorted, SliceTooShort };

/// Adds `x`'s counts to `counts[0 .. edges.len - 1]`; values outside the
/// edges, and NaNs, are not counted.
pub fn histogramInto(comptime T: type, edges: []const T, x: ConstStrided(T), counts: []u64) HistogramIntoError!void {
    try checkBins(T, edges, counts.len);
    for (0..x.n) |i| if (bin(T, edges, x.get(i))) |k| {
        counts[k] += 1;
    };
}

/// `edges` bound at least one bin, ascend (equal neighbors allowed, as
/// NumPy's) with no NaN, and `counts` holds a count per bin.
pub fn checkBins(comptime T: type, edges: []const T, counts: usize) HistogramIntoError!void {
    if (@typeInfo(T) != .float) @compileError("rstd.math.kernels.stats: " ++ @typeName(T) ++ " is not a float type");
    if (edges.len < 2) return error.ZeroBins;
    if (counts < edges.len - 1) return error.SliceTooShort;
    for (edges[0 .. edges.len - 1], edges[1..]) |a, b| if (!(a <= b)) return error.EdgesNotSorted;
}

/// The bin of `x` in `edges` (checked by `checkBins`): the last `k` with
/// `edges[k] <= x`, the right edge in the last bin; null outside, or NaN.
pub fn bin(comptime T: type, edges: []const T, x: T) ?usize {
    const last = edges.len - 1;
    if (!(x >= edges[0] and x <= edges[last])) return null;
    var lo: usize = 0;
    var hi: usize = last;
    while (hi - lo > 1) {
        const mid = lo + (hi - lo) / 2;
        if (edges[mid] <= x) lo = mid else hi = mid;
    }
    return lo;
}
