//! Tests of `rstd.stats` on slices: the moments, `Running`, the median and quantiles, and histograms.

const std = @import("std");
const stats = @import("../__root__.zig");
const mean = stats.mean;
const meanWith = stats.meanWith;
const variance = stats.variance;
const stdDev = stats.stdDev;
const stdDevWith = stats.stdDevWith;
const Running = stats.Running;
const median = stats.median;
const QuantileMethod = stats.QuantileMethod;
const quantile = stats.quantile;
const quantileInto = stats.quantileInto;
const histogram = stats.histogram;
const histogramInto = stats.histogramInto;

const testing = std.testing;

const sample = [_]f64{ 2, 4, 4, 4, 5, 5, 7, 9 };

test "mean, variance, stdDev" {
    const m = try mean(testing.io, f64, &sample);
    try testing.expectEqual(5, m);
    try testing.expectEqual(4, try variance(testing.io, f64, &sample, .population));
    try testing.expectApproxEqRel(32.0 / 7.0, try variance(testing.io, f64, &sample, .sample), 1e-15);
    try testing.expectEqual(2, try stdDev(testing.io, f64, &sample, .population));
    try testing.expectEqual(0, try variance(testing.io, f64, sample[0..1], .sample));
    try testing.expectError(error.EmptyInput, mean(testing.io, f64, &[_]f64{}));
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
    const odd = try median(testing.allocator, f64, &[_]f64{ 3, 1, 2 });
    try testing.expectEqual(2, odd);
    try testing.expectEqual(2.5, try median(testing.allocator, f64, &[_]f64{ 4, 1, 3, 2 }));
    try testing.expect(std.math.isNan(try median(testing.allocator, f64, &[_]f64{ 1, std.math.nan(f64), 2 })));
}

test "quantile methods" {
    // design.md's worked example: [1, 2, 3, 4] at 0.25
    const xs = [_]f64{ 4, 2, 1, 3 };
    const p = [_]f64{0.25};
    var out: [1]f64 = undefined;
    try quantileInto(testing.allocator, f64, &out, &xs, &p, .{});
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
    try quantileInto(al, f64, &out, &[_]f64{ 1, std.math.nan(f64) }, &.{0.5}, .{});
    try testing.expect(std.math.isNan(out[0]));
    try testing.expectEqual(5, try median(al, f16, &[_]f16{ 9, 5, 1 }));
}

test "histogram" {
    const xs = [_]f64{ 0, 0.5, 1, 1.5, 2 };
    var counts: [2]u64 = @splat(0);
    try histogramInto(f64, &counts, &xs, &.{ 0, 1, 2 });
    try testing.expectEqualSlices(u64, &.{ 2, 3 }, &counts); // 2 falls in the last bin

    var h = try histogram(testing.allocator, f64, &xs, .{ .count = 4 });
    defer h.deinit(testing.allocator);
    try testing.expectEqualSlices(u64, &.{ 1, 1, 1, 2 }, h.counts);
    try testing.expectEqualSlices(f64, &.{ 0, 0.5, 1, 1.5, 2 }, h.edges);
    try testing.expectError(error.EdgesNotSorted, histogramInto(f64, &counts, &xs, &.{ 1, 0, 2 }));
}

test "histogram: equal values, NaN and out-of-range values, edges given" {
    const al = testing.allocator;
    var one = try histogram(al, f64, &[_]f64{ 3, 3, std.math.nan(f64) }, .{ .count = 2 });
    defer one.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 2.5, 3, 3.5 }, one.edges);
    try testing.expectEqualSlices(u64, &.{ 0, 2 }, one.counts);
    var given = try histogram(al, f32, &[_]f32{ -1, 0, 0.5, 1, 9 }, .{ .edges = &.{ 0, 1 } });
    defer given.deinit(al);
    try testing.expectEqualSlices(u64, &.{3}, given.counts);
    try testing.expectError(error.EmptyInput, histogram(al, f64, &[_]f64{std.math.nan(f64)}, .{ .count = 3 }));
    try testing.expectError(error.ZeroBins, histogram(al, f64, &[_]f64{1}, .{ .count = 0 }));
    try testing.expectError(error.EdgesNotSorted, histogram(al, f64, &[_]f64{1}, .{ .edges = &.{ 0, std.math.nan(f64) } }));
}

test "type-agnostic: integers and Complex in, the caller's result type out" {
    const io = testing.io;
    const al = testing.allocator;
    const ints = [_]i32{ 2, 4, 4, 4, 5, 5, 7, 9 };
    try testing.expectEqual(5, try mean(io, f64, &ints));
    try testing.expectEqual(2, try stdDev(io, f32, &ints, .population));
    try testing.expectEqual(4.5, try median(al, f64, &ints));
    var q: [1]f64 = undefined;
    try quantileInto(al, f64, &q, &[_]u8{ 1, 2, 3, 4 }, &.{0.25}, .{});
    try testing.expectEqual(1.75, q[0]);
    var h = try histogram(al, f64, &[_]i64{ 0, 1, 1, 2 }, .{ .count = 2 });
    defer h.deinit(al);
    try testing.expectEqualSlices(u64, &.{ 1, 3 }, h.counts);
    const C = std.math.Complex(f64);
    const cs = [_]C{ .init(1, 2), .init(3, -2) };
    const m = try mean(io, C, &cs);
    try testing.expectEqual(2, m.re);
    try testing.expectEqual(0, m.im);
    // |x - mean|²: (1 + 4) and (1 + 4)
    try testing.expectEqual(5, try variance(io, f64, &cs, .population));
}

test "median of Complex: MATLAB's order" {
    const C = std.math.Complex(f64);
    const xs = [_]C{ .init(3, 0), .init(0, 1), .init(-2, 0) };
    const m = try median(testing.allocator, C, &xs);
    try testing.expectEqual(-2, m.re);
    try testing.expectEqual(0, m.im);
}
