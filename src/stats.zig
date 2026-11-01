//! Descriptive statistics over slices (rstd's design.md, "Statistics"). Stubs
//! for now: bodies return `error.NotImplemented`, which leaves each error set
//! when its body lands, and their tests skip until then. `Running` is real.
//!
//! NaN propagates, as in MATLAB. `cov` and `corrcoef` come with the public
//! layer's views, the forms along a dimension with N-d arrays.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Io = std.Io;

/// The variance's divisor: `n - 1` (MATLAB's default) or `n`. A sample of
/// one has variance 0 either way, as MATLAB's `var`.
pub const Normalize = enum { sample, population };

pub const MomentError = error{ EmptyInput, Canceled, NotImplemented };

/// The mean, through the BLAS `sum`.
pub fn mean(io: Io, comptime F: type, xs: []const F) MomentError!F {
    comptime checkFloat(F);
    _ = io;
    if (xs.len == 0) return error.EmptyInput;
    return error.NotImplemented;
}

/// Two passes: the mean, then the squared deviations.
pub fn variance(io: Io, comptime F: type, xs: []const F, normalize: Normalize) MomentError!F {
    comptime checkFloat(F);
    _ = io;
    _ = normalize;
    if (xs.len == 0) return error.EmptyInput;
    return error.NotImplemented;
}

pub fn stdDev(io: Io, comptime F: type, xs: []const F, normalize: Normalize) MomentError!F {
    return @sqrt(try variance(io, F, xs, normalize));
}

/// Streaming mean and variance (Welford's update). Its bits differ from
/// `mean`'s and `variance`'s.
pub fn Running(comptime F: type) type {
    comptime checkFloat(F);
    return struct {
        count: u64 = 0,
        avg: F = 0,
        /// Sum of squared deviations from `avg`.
        m2: F = 0,

        const Self = @This();

        pub fn push(self: *Self, x: F) void {
            self.count += 1;
            const d = x - self.avg;
            self.avg += d / @as(F, @floatFromInt(self.count));
            self.m2 += d * (x - self.avg);
        }

        pub fn mean(self: Self) error{EmptyInput}!F {
            if (self.count == 0) return error.EmptyInput;
            return self.avg;
        }

        pub fn variance(self: Self, normalize: Normalize) error{EmptyInput}!F {
            if (self.count == 0) return error.EmptyInput;
            const d = switch (normalize) {
                .sample => if (self.count == 1) return 0 else self.count - 1,
                .population => self.count,
            };
            return self.m2 / @as(F, @floatFromInt(d));
        }
    };
}

pub const MedianError = error{ EmptyInput, OutOfMemory, NotImplemented };

/// Selects on a copy; `xs` is untouched.
pub fn median(allocator: Allocator, comptime F: type, xs: []const F) MedianError!F {
    comptime checkFloat(F);
    _ = allocator;
    if (xs.len == 0) return error.EmptyInput;
    return error.NotImplemented;
}

/// Hyndman and Fan's nine (The American Statistician 50(4), 1996), in their
/// order, under NumPy's names. `linear` is type 7 (NumPy, R), `hazen` type 5
/// (MATLAB), `weibull` type 6 (SPSS).
pub const QuantileMethod = enum {
    inverted_cdf,
    averaged_inverted_cdf,
    closest_observation,
    interpolated_inverted_cdf,
    hazen,
    weibull,
    linear,
    median_unbiased,
    normal_unbiased,
};

pub const QuantileOptions = struct { method: QuantileMethod = .linear };

pub const QuantileError = error{ EmptyInput, OutOfMemory, ProbabilityOutOfRange, NotImplemented };
pub const QuantileIntoError = QuantileError || error{SliceTooShort};

/// One quantile of `xs` per `p`, each in [0, 1]; owned by the caller.
pub fn quantile(allocator: Allocator, comptime F: type, xs: []const F, p: []const F, options: QuantileOptions) QuantileError![]F {
    const out = try allocator.alloc(F, p.len);
    errdefer allocator.free(out);
    try quantiles(allocator, F, out, xs, p, options);
    return out;
}

/// `quantile` into `out[0..p.len]`; the allocator holds the sorted copy.
pub fn quantileInto(allocator: Allocator, comptime F: type, out: []F, xs: []const F, p: []const F, options: QuantileOptions) QuantileIntoError!void {
    if (out.len < p.len) return error.SliceTooShort;
    return quantiles(allocator, F, out, xs, p, options);
}

fn quantiles(allocator: Allocator, comptime F: type, out: []F, xs: []const F, p: []const F, options: QuantileOptions) QuantileError!void {
    comptime checkFloat(F);
    _ = allocator;
    _ = out;
    _ = options;
    if (xs.len == 0) return error.EmptyInput;
    for (p) |q| if (!(q >= 0 and q <= 1)) return error.ProbabilityOutOfRange;
    return error.NotImplemented;
}

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

pub const HistogramError = error{ EmptyInput, OutOfMemory, ZeroBins, EdgesNotSorted, NotImplemented };
pub const HistogramIntoError = error{ ZeroBins, EdgesNotSorted, SliceTooShort, NotImplemented };

/// `.count` bins need a nonempty `xs` for their range.
pub fn histogram(allocator: Allocator, comptime F: type, xs: []const F, bins: Bins(F)) HistogramError!Histogram(F) {
    comptime checkFloat(F);
    _ = allocator;
    switch (bins) {
        .count => |n| {
            if (n == 0) return error.ZeroBins;
            if (xs.len == 0) return error.EmptyInput;
        },
        .edges => |e| if (e.len < 2) return error.ZeroBins,
    }
    return error.NotImplemented;
}

/// Adds to `counts[0 .. edges.len - 1]`; values outside the edges, and NaNs,
/// are not counted.
pub fn histogramInto(comptime F: type, counts: []u64, xs: []const F, edges: []const F) HistogramIntoError!void {
    comptime checkFloat(F);
    _ = xs;
    if (edges.len < 2) return error.ZeroBins;
    if (counts.len < edges.len - 1) return error.SliceTooShort;
    return error.NotImplemented;
}

fn checkFloat(comptime F: type) void {
    if (@typeInfo(F) != .float) @compileError("rstd.stats: " ++ @typeName(F) ++ " is not a float type");
}

// ---- tests ----

const testing = std.testing;

/// A stub's tests skip until its body lands.
fn skip(err: anyerror) anyerror {
    return if (err == error.NotImplemented) error.SkipZigTest else err;
}

const sample = [_]f64{ 2, 4, 4, 4, 5, 5, 7, 9 };

test "mean, variance, stdDev" {
    const m = mean(testing.io, f64, &sample) catch |e| return skip(e);
    try testing.expectEqual(5, m);
    try testing.expectEqual(4, try variance(testing.io, f64, &sample, .population));
    try testing.expectApproxEqRel(32.0 / 7.0, try variance(testing.io, f64, &sample, .sample), 1e-15);
    try testing.expectEqual(2, try stdDev(testing.io, f64, &sample, .population));
    try testing.expectEqual(0, try variance(testing.io, f64, sample[0..1], .sample));
    try testing.expectError(error.EmptyInput, mean(testing.io, f64, &.{}));
}

test "Running" {
    var r: Running(f64) = .{};
    try testing.expectError(error.EmptyInput, r.mean());
    for (sample) |x| r.push(x);
    try testing.expectApproxEqRel(5, try r.mean(), 1e-15);
    try testing.expectApproxEqRel(4, try r.variance(.population), 1e-14);
    try testing.expectApproxEqRel(32.0 / 7.0, try r.variance(.sample), 1e-14);

    var one: Running(f32) = .{};
    one.push(3);
    try testing.expectEqual(0, try one.variance(.sample));
}

test "median" {
    const odd = median(testing.allocator, f64, &.{ 3, 1, 2 }) catch |e| return skip(e);
    try testing.expectEqual(2, odd);
    try testing.expectEqual(2.5, try median(testing.allocator, f64, &.{ 4, 1, 3, 2 }));
    try testing.expect(std.math.isNan(try median(testing.allocator, f64, &.{ 1, std.math.nan(f64), 2 })));
}

test "quantile methods" {
    // design.md's worked example: [1, 2, 3, 4] at 0.25
    const xs = [_]f64{ 4, 2, 1, 3 };
    const p = [_]f64{0.25};
    var out: [1]f64 = undefined;
    quantileInto(testing.allocator, f64, &out, &xs, &p, .{}) catch |e| return skip(e);
    try testing.expectApproxEqAbs(1.75, out[0], 1e-15);
    try quantileInto(testing.allocator, f64, &out, &xs, &p, .{ .method = .hazen });
    try testing.expectApproxEqAbs(1.5, out[0], 1e-15);
    try quantileInto(testing.allocator, f64, &out, &xs, &p, .{ .method = .weibull });
    try testing.expectApproxEqAbs(1.25, out[0], 1e-15);
    try testing.expectError(error.ProbabilityOutOfRange, quantileInto(testing.allocator, f64, &out, &xs, &.{1.5}, .{}));

    const owned = try quantile(testing.allocator, f64, &xs, &.{ 0, 0.5, 1 }, .{});
    defer testing.allocator.free(owned);
    try testing.expectEqualSlices(f64, &.{ 1, 2.5, 4 }, owned);
}

test "histogram" {
    const xs = [_]f64{ 0, 0.5, 1, 1.5, 2 };
    var counts: [2]u64 = @splat(0);
    histogramInto(f64, &counts, &xs, &.{ 0, 1, 2 }) catch |e| return skip(e);
    try testing.expectEqualSlices(u64, &.{ 2, 3 }, &counts); // 2 falls in the last bin

    var h = try histogram(testing.allocator, f64, &xs, .{ .count = 4 });
    defer h.deinit(testing.allocator);
    try testing.expectEqualSlices(u64, &.{ 1, 1, 1, 2 }, h.counts);
    try testing.expectEqualSlices(f64, &.{ 0, 0.5, 1, 1.5, 2 }, h.edges);
    try testing.expectError(error.EdgesNotSorted, histogramInto(f64, &counts, &xs, &.{ 1, 0, 2 }));
}
