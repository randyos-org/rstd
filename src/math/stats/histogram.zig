//! Histograms: counts over even bins or given edges, over
//! `kernels.stats`' bins.

const std = @import("std");
const kernels_stats = @import("../kernels/stats/__root__.zig");
const util = @import("util.zig");

const Allocator = std.mem.Allocator;
const Elem = util.Elem;
const to = util.to;
const checkReal = util.checkReal;
const checkFloat = util.checkFloat;

pub const Bins = kernels_stats.Bins;
pub const Histogram = kernels_stats.Histogram;
pub const HistogramError = error{ EmptyInput, OutOfMemory, ZeroBins, EdgesNotSorted };
pub const HistogramIntoError = kernels_stats.HistogramIntoError;

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
            kernels_stats.checkBins(F, given, given.len -| 1) catch |e| switch (e) {
                error.ZeroBins => return error.ZeroBins,
                error.EdgesNotSorted => return error.EdgesNotSorted,
                error.SliceTooShort => unreachable,
            };
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
/// are not counted. Floats of type `F` count through
/// `kernels.stats.histogramInto`, other elements here, by the same bins.
pub fn histogramInto(comptime F: type, counts: []u64, xs: anytype, edges: []const F) HistogramIntoError!void {
    const E = Elem(@TypeOf(xs));
    comptime checkFloat(F);
    comptime checkReal(E);
    if (comptime E == F) {
        const sl: []const F = xs;
        return kernels_stats.histogramInto(F, edges, .of(sl), counts);
    }
    try kernels_stats.checkBins(F, edges, counts.len);
    for (xs) |x| if (kernels_stats.bin(F, edges, to(F, x))) |k| {
        counts[k] += 1;
    };
}
