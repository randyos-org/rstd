//! Descriptive statistics over slices (rstd's design.md, "Statistics"): the
//! moments through the BLAS, `Running`, the median and quantiles, and
//! histograms.
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

/// The mean of `xs` (a slice or array of integers, floats or `Complex`), as
/// `R`: a float for real elements, a `Complex` for complex ones. Floats of
/// `R`'s own type sum through the BLAS `sum`.
pub fn mean(io: Io, comptime R: type, xs: anytype) MomentError!R {
    return meanWith(io, R, xs, .plain);
}

/// `mean`, summed as `s` says.
pub fn meanWith(io: Io, comptime R: type, xs: anytype, s: Summation) MomentError!R {
    const E = Elem(@TypeOf(xs));
    comptime checkPair(R, E);
    if (xs.len == 0) return error.EmptyInput;
    const n: Real(R) = @floatFromInt(xs.len);
    if (comptime isComplex(E)) {
        const W = Wide(Real(R));
        var re: Acc(W) = .{};
        var im: Acc(W) = .{};
        for (xs) |x| {
            re.add(to(W, x.re), s);
            im.add(to(W, x.im), s);
        }
        return .{ .re = @as(Real(R), @floatCast(re.value())) / n, .im = @as(Real(R), @floatCast(im.value())) / n };
    }
    return try total(io, R, xs, s) / n;
}

/// Two passes: the mean, then the squared deviations (`|x - mean|²` for
/// complex `xs`). `R` is a float.
pub fn variance(io: Io, comptime R: type, xs: anytype, normalize: Normalize) MomentError!R {
    return varianceWith(io, R, xs, normalize, .plain);
}

/// `variance`, both passes summed as `s` says; the deviations' squares
/// summed in order, in `f64` or wider.
pub fn varianceWith(io: Io, comptime R: type, xs: anytype, normalize: Normalize, s: Summation) MomentError!R {
    const E = Elem(@TypeOf(xs));
    comptime checkFloat(R);
    const W = Wide(R);
    var acc: Acc(W) = .{};
    if (comptime isComplex(E)) {
        const m = try meanWith(io, std.math.Complex(R), xs, s);
        for (xs) |x| {
            const dr = to(W, x.re) - @as(W, m.re);
            const di = to(W, x.im) - @as(W, m.im);
            acc.add(dr * dr + di * di, s);
        }
    } else {
        const m = try meanWith(io, R, xs, s);
        for (xs) |x| {
            const d = to(W, x) - @as(W, m);
            acc.add(d * d, s);
        }
    }
    const n = xs.len;
    // a sample of one: 0, or NaN from a NaN or an infinity
    const d = if (normalize == .sample and n > 1) n - 1 else n;
    return @floatCast(acc.value() / @as(W, @floatFromInt(d)));
}

pub fn stdDev(io: Io, comptime R: type, xs: anytype, normalize: Normalize) MomentError!R {
    return @sqrt(try variance(io, R, xs, normalize));
}

/// `stdDev`, summed as `s` says.
pub fn stdDevWith(io: Io, comptime R: type, xs: anytype, normalize: Normalize, s: Summation) MomentError!R {
    return @sqrt(try varianceWith(io, R, xs, normalize, s));
}

/// The sum of real `xs` as `R`: the BLAS's where the elements are `R` and
/// it has a `sum` for them, else in order here.
fn total(io: Io, comptime R: type, xs: anytype, s: Summation) error{Canceled}!R {
    const E = Elem(@TypeOf(xs));
    if (comptime E == R and (R == f32 or R == f64 or R == blas_types.XDouble)) {
        const sl: []const R = xs;
        return blas.sumWith(io, R, .of(sl), s);
    }
    const W = Wide(R);
    var acc: Acc(W) = .{};
    for (xs) |x| acc.add(to(W, x), s);
    return @floatCast(acc.value());
}

fn Wide(comptime F: type) type {
    return if (@bitSizeOf(F) < 64) f64 else F;
}

/// The element type of a slice, array or pointer to one.
fn Elem(comptime S: type) type {
    return std.meta.Elem(S);
}

/// `std.math.Complex` and other `re`, `im` structs.
fn isComplex(comptime T: type) bool {
    return @typeInfo(T) == .@"struct" and @hasField(T, "re") and @hasField(T, "im");
}

fn Real(comptime R: type) type {
    return if (isComplex(R)) @FieldType(R, "re") else R;
}

/// A real element as the float `W`.
inline fn to(comptime W: type, x: anytype) W {
    return switch (@typeInfo(@TypeOf(x))) {
        .int, .comptime_int => @floatFromInt(x),
        .float, .comptime_float => @floatCast(x),
        else => @compileError("rstd.stats: " ++ @typeName(@TypeOf(x)) ++ " is not a real number type"),
    };
}

/// `R` a float for real `E` (an integer or a float), a `Complex` of a float
/// for complex `E`.
fn checkPair(comptime R: type, comptime E: type) void {
    if (isComplex(E)) {
        if (!isComplex(R)) @compileError("rstd.stats: complex " ++ @typeName(E) ++ " needs a Complex result, not " ++ @typeName(R));
        checkFloat(Real(R));
    } else {
        checkReal(E);
        checkFloat(R);
    }
}

fn checkReal(comptime E: type) void {
    switch (@typeInfo(E)) {
        .int, .float => {},
        else => @compileError("rstd.stats: " ++ @typeName(E) ++ " is not an integer or float type"),
    }
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

/// The median of real `xs` as the float `R`. Sorts a copy; `xs` is
/// untouched. NaN if `xs` holds one, as MATLAB's.
pub fn median(allocator: Allocator, comptime R: type, xs: anytype) MedianError!R {
    const E = Elem(@TypeOf(xs));
    comptime checkOrdered(E);
    comptime checkFloat(R);
    if (xs.len == 0) return error.EmptyInput;
    const s = try sortedCopy(allocator, E, xs);
    defer allocator.free(s);
    const n = s.len;
    if (comptime @typeInfo(E) == .float) if (std.math.isNan(s[n - 1])) return to(R, s[n - 1]);
    if (n % 2 == 1) return to(R, s[n / 2]);
    return lerp(R, to(R, s[n / 2 - 1]), to(R, s[n / 2]), @as(Wide(R), 0.5));
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
pub fn quantile(allocator: Allocator, comptime R: type, xs: anytype, p: []const R, options: QuantileOptions) QuantileError![]R {
    const out = try allocator.alloc(R, p.len);
    errdefer allocator.free(out);
    try quantiles(allocator, R, out, xs, p, options);
    return out;
}

/// `quantile` into `out[0..p.len]`; the allocator holds the sorted copy.
pub fn quantileInto(allocator: Allocator, comptime R: type, out: []R, xs: anytype, p: []const R, options: QuantileOptions) QuantileIntoError!void {
    if (out.len < p.len) return error.SliceTooShort;
    return quantiles(allocator, R, out, xs, p, options);
}

fn quantiles(allocator: Allocator, comptime R: type, out: []R, xs: anytype, p: []const R, options: QuantileOptions) QuantileError!void {
    const E = Elem(@TypeOf(xs));
    comptime checkOrdered(E);
    comptime checkFloat(R);
    if (xs.len == 0) return error.EmptyInput;
    for (p) |q| if (!(q >= 0 and q <= 1)) return error.ProbabilityOutOfRange;
    const s = try sortedCopy(allocator, E, xs);
    defer allocator.free(s);
    // a NaN, sorted last, makes every quantile NaN, as NumPy's
    if (comptime @typeInfo(E) == .float) if (std.math.isNan(s[s.len - 1])) {
        @memset(out[0..p.len], to(R, s[s.len - 1]));
        return;
    };
    for (out[0..p.len], p) |*o, q| o.* = quantileOf(R, E, s, q, options.method);
}

/// NumPy's `_quantile` on the sorted `s` (numpy/lib/_function_base_impl.py,
/// v2.1): the method's virtual index, then the two neighbors and `gamma`.
fn quantileOf(comptime R: type, comptime E: type, s: []const E, p: R, method: QuantileMethod) R {
    const W = Wide(R);
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
            return to(R, s[index(W, k, last)]);
        },
        .averaged_inverted_cdf => {
            const vi = n * q - 1;
            const prev = @floor(vi);
            const gamma: W = if (vi == prev) 0.5 else 1;
            return lerp(R, to(R, s[index(W, prev, last)]), to(R, s[index(W, prev + 1, last)]), gamma);
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
            if (vi < 0) return to(R, s[0]);
            if (vi >= @as(W, @floatFromInt(last))) return to(R, s[last]);
            const prev = @floor(vi);
            return lerp(R, to(R, s[index(W, prev, last)]), to(R, s[index(W, prev + 1, last)]), vi - prev);
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
fn sortedCopy(allocator: Allocator, comptime E: type, xs: anytype) Allocator.Error![]E {
    const s = try allocator.alloc(E, xs.len);
    for (s, xs) |*d, x| d.* = x;
    std.mem.sort(E, s, {}, struct {
        fn lt(_: void, a: E, b: E) bool {
            if (comptime @typeInfo(E) == .float) {
                if (std.math.isNan(a)) return false;
                if (std.math.isNan(b)) return true;
            }
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

/// Sortable: integers and floats; `Complex` waits for its order (N32, #76).
fn checkOrdered(comptime E: type) void {
    if (isComplex(E)) @compileError("rstd.stats: ordering " ++ @typeName(E) ++ " needs MATLAB's complex order (N32, #76)");
    checkReal(E);
}

fn checkFloat(comptime F: type) void {
    if (@typeInfo(F) != .float) @compileError("rstd.stats: " ++ @typeName(F) ++ " is not a float type");
}

// ---- tests ----

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
