//! Descriptive statistics over slices (rstd's design.md, "Statistics"). The
//! moments and `Running` are done; the rest are stubs whose bodies return
//! `error.NotImplemented`, which leaves each error set when its body lands,
//! and their tests skip until then.
//!
//! NaN propagates, as in MATLAB. `cov` and `corrcoef` come with the public
//! layer's views, the forms along a dimension with N-d arrays.

const std = @import("std");
const blas = @import("linalg/blas/__root__.zig");
const blas_types = @import("linalg/blas/types.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;

/// The variance's divisor: `n - 1` (MATLAB's default) or `n`. A sample of
/// one has variance 0 either way, as MATLAB's `var`.
pub const Normalize = enum { sample, population };

pub const MomentError = error{ EmptyInput, Canceled };

/// How the sums are taken: as the BLAS's kernels do, or compensated
/// (Neumaier's, `blas.Summation`).
pub const Summation = blas.Summation;

/// The mean, through the BLAS `sum`.
pub fn mean(io: Io, comptime F: type, xs: []const F) MomentError!F {
    return meanWith(io, F, xs, .plain);
}

/// `mean`, summed as `s` says.
pub fn meanWith(io: Io, comptime F: type, xs: []const F, s: Summation) MomentError!F {
    comptime checkFloat(F);
    if (xs.len == 0) return error.EmptyInput;
    return try total(io, F, xs, s) / @as(F, @floatFromInt(xs.len));
}

/// Two passes: the mean, then the squared deviations.
pub fn variance(io: Io, comptime F: type, xs: []const F, normalize: Normalize) MomentError!F {
    return varianceWith(io, F, xs, normalize, .plain);
}

/// `variance`, both passes summed as `s` says; the deviations' squares
/// summed in order, in `f64` or wider.
pub fn varianceWith(io: Io, comptime F: type, xs: []const F, normalize: Normalize, s: Summation) MomentError!F {
    const m = try meanWith(io, F, xs, s);
    const W = Wide(F);
    var acc: Acc(W) = .{};
    for (xs) |x| {
        const d: W = @as(W, x) - @as(W, m);
        acc.add(d * d, s);
    }
    const n = xs.len;
    // a sample of one: 0, or NaN from a NaN or an infinity
    const d = if (normalize == .sample and n > 1) n - 1 else n;
    return @floatCast(acc.value() / @as(W, @floatFromInt(d)));
}

pub fn stdDev(io: Io, comptime F: type, xs: []const F, normalize: Normalize) MomentError!F {
    return @sqrt(try variance(io, F, xs, normalize));
}

/// `stdDev`, summed as `s` says.
pub fn stdDevWith(io: Io, comptime F: type, xs: []const F, normalize: Normalize, s: Summation) MomentError!F {
    return @sqrt(try varianceWith(io, F, xs, normalize, s));
}

/// The sum of `xs`: the BLAS's for its types, else in order here.
fn total(io: Io, comptime F: type, xs: []const F, s: Summation) error{Canceled}!F {
    if (comptime F == f32 or F == f64 or F == blas_types.XDouble) return blas.sumWith(io, F, .of(xs), s);
    var acc: Acc(Wide(F)) = .{};
    for (xs) |x| acc.add(x, s);
    return @floatCast(acc.value());
}

fn Wide(comptime F: type) type {
    return if (@bitSizeOf(F) < 64) f64 else F;
}

/// A running sum, plain or with Neumaier's compensation.
fn Acc(comptime W: type) type {
    return struct {
        sum: W = 0,
        c: W = 0,

        fn add(a: *@This(), x: W, s: Summation) void {
            const t = a.sum + x;
            if (s == .compensated) {
                if (@abs(a.sum) >= @abs(x)) a.c += (a.sum - t) + x else a.c += (x - t) + a.sum;
            }
            a.sum = t;
        }

        fn value(a: @This()) W {
            return a.sum + a.c;
        }
    };
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

pub const MedianError = error{ EmptyInput, OutOfMemory };

/// Sorts a copy; `xs` is untouched. NaN if `xs` holds one, as MATLAB's.
pub fn median(allocator: Allocator, comptime F: type, xs: []const F) MedianError!F {
    comptime checkFloat(F);
    if (xs.len == 0) return error.EmptyInput;
    const s = try sortedCopy(allocator, F, xs);
    defer allocator.free(s);
    const n = s.len;
    if (std.math.isNan(s[n - 1])) return s[n - 1];
    if (n % 2 == 1) return s[n / 2];
    return lerp(F, s[n / 2 - 1], s[n / 2], @as(Wide(F), 0.5));
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

pub const QuantileError = error{ EmptyInput, OutOfMemory, ProbabilityOutOfRange };
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
    if (xs.len == 0) return error.EmptyInput;
    for (p) |q| if (!(q >= 0 and q <= 1)) return error.ProbabilityOutOfRange;
    const s = try sortedCopy(allocator, F, xs);
    defer allocator.free(s);
    // a NaN, sorted last, makes every quantile NaN, as NumPy's
    if (std.math.isNan(s[s.len - 1])) {
        @memset(out[0..p.len], s[s.len - 1]);
        return;
    }
    for (out[0..p.len], p) |*o, q| o.* = quantileOf(F, s, q, options.method);
}

/// NumPy's `_quantile` on the sorted `s` (numpy/lib/_function_base_impl.py,
/// v2.1): the method's virtual index, then the two neighbors and `gamma`.
fn quantileOf(comptime F: type, s: []const F, p: F, method: QuantileMethod) F {
    const W = if (@bitSizeOf(F) < 64) f64 else F;
    const n: W = @floatFromInt(s.len);
    const q: W = p;
    const last = s.len - 1;
    switch (method) {
        .inverted_cdf, .closest_observation => {
            // the discontinuous ones pick an order statistic
            const vi = if (method == .inverted_cdf) n * q - 1 else n * q - 1 - 0.5;
            const prev = @floor(vi);
            const exact = vi == prev;
            const take_prev = if (method == .inverted_cdf) exact else exact and @mod(prev, 2) == 0;
            const k = if (take_prev) prev else prev + 1;
            return s[index(W, k, last)];
        },
        .averaged_inverted_cdf => {
            const vi = n * q - 1;
            const prev = @floor(vi);
            const gamma: W = if (vi == prev) 0.5 else 1;
            return lerp(F, s[index(W, prev, last)], s[index(W, prev + 1, last)], gamma);
        },
        else => {
            const ab: [2]W = switch (method) {
                .interpolated_inverted_cdf => .{ 0, 1 },
                .hazen => .{ 0.5, 0.5 },
                .weibull => .{ 0, 0 },
                .linear => .{ 1, 1 },
                .median_unbiased => .{ 1.0 / 3.0, 1.0 / 3.0 },
                .normal_unbiased => .{ 3.0 / 8.0, 3.0 / 8.0 },
                .inverted_cdf, .averaged_inverted_cdf, .closest_observation => unreachable,
            };
            const vi = n * q + (ab[0] + q * (1 - ab[0] - ab[1])) - 1;
            // `_get_indexes`: past either end, that end
            if (vi < 0) return s[0];
            if (vi >= @as(W, @floatFromInt(last))) return s[last];
            const prev = @floor(vi);
            return lerp(F, s[index(W, prev, last)], s[index(W, prev + 1, last)], vi - prev);
        },
    }
}

/// A virtual index clamped to `0..last`.
fn index(comptime W: type, k: W, last: usize) usize {
    if (!(k > 0)) return 0;
    const u: usize = @intFromFloat(@min(k, @as(W, @floatFromInt(last))));
    return @min(u, last);
}

/// NumPy's `_lerp`: from `a` for `t < 0.5`, from `b` past it, so the ends
/// are exact.
fn lerp(comptime F: type, a: F, b: F, t: anytype) F {
    const W = @TypeOf(t);
    const aw: W = a;
    const bw: W = b;
    const d = bw - aw;
    return @floatCast(if (t >= 0.5) bw - d * (1 - t) else aw + d * t);
}

/// `xs` sorted ascending, NaN last, owned by the caller.
fn sortedCopy(allocator: Allocator, comptime F: type, xs: []const F) Allocator.Error![]F {
    const s = try allocator.dupe(F, xs);
    std.mem.sort(F, s, {}, struct {
        fn lt(_: void, a: F, b: F) bool {
            if (std.math.isNan(a)) return false;
            if (std.math.isNan(b)) return true;
            return a < b;
        }
    }.lt);
    return s;
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

test "moments: compensated sums, the types the BLAS lacks, NaN" {
    const io = testing.io;
    const hard = [_]f64{ 1e16, 1, -1e16, 1 };
    try testing.expectEqual(0.5, try meanWith(io, f64, &hard, .compensated));
    try testing.expectEqual(try mean(io, f64, &hard), try meanWith(io, f64, &hard, .plain));
    const h = [_]f16{ 1, 2, 3, 4 };
    try testing.expectEqual(2.5, try mean(io, f16, &h));
    try testing.expectEqual(1.25, try variance(io, f16, &h, .population));
    const f = [_]f32{ 2, 4, 4, 4, 5, 5, 7, 9 };
    try testing.expectEqual(2, try stdDevWith(io, f32, &f, .population, .compensated));
    const nan = [_]f64{ 1, std.math.nan(f64), 3 };
    try testing.expect(std.math.isNan(try mean(io, f64, &nan)));
    try testing.expect(std.math.isNan(try variance(io, f64, &nan, .sample)));
    try testing.expect(std.math.isNan(try variance(io, f64, nan[1..2], .sample)));
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

test "quantile: all nine methods on [1, 2, 3, 4], NumPy's values" {
    const al = testing.allocator;
    const xs = [_]f64{ 3, 1, 4, 2 };
    const want = [_][2]f64{
        .{ 1, 2 }, // inverted_cdf at 0.25, 0.5
        .{ 1.5, 2.5 }, // averaged_inverted_cdf
        .{ 1, 2 }, // closest_observation
        .{ 1, 2 }, // interpolated_inverted_cdf
        .{ 1.5, 2.5 }, // hazen
        .{ 1.25, 2.5 }, // weibull
        .{ 1.75, 2.5 }, // linear
        .{ 1.0 + 5.0 / 12.0, 2.5 }, // median_unbiased
        .{ 1.4375, 2.5 }, // normal_unbiased
    };
    for (std.enums.values(QuantileMethod), want) |m, w| {
        var out: [2]f64 = undefined;
        try quantileInto(al, f64, &out, &xs, &.{ 0.25, 0.5 }, .{ .method = m });
        try testing.expectApproxEqAbs(w[0], out[0], 1e-15);
        try testing.expectApproxEqAbs(w[1], out[1], 1e-15);
    }
    var out: [1]f64 = undefined;
    try quantileInto(al, f64, &out, &.{ 1, std.math.nan(f64) }, &.{0.5}, .{});
    try testing.expect(std.math.isNan(out[0]));
    try testing.expectEqual(5, try median(al, f16, &.{ 9, 5, 1 }));
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
