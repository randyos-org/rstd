//! The median and quantiles (Hyndman and Fan's nine) of a slice, over a
//! sorted copy.

const std = @import("std");
const math = @import("../math.zig");
const sort = @import("../sort.zig");
const util = @import("util.zig");

const Allocator = std.mem.Allocator;
const Wide = util.Wide;
const Elem = util.Elem;
const isComplex = util.isComplex;
const Real = util.Real;
const to = util.to;
const checkPair = util.checkPair;
const checkReal = util.checkReal;
const checkFloat = util.checkFloat;

pub const MedianError = error{ EmptyInput, OutOfMemory };

/// The median of `xs` as `R` (a `Complex` for complex `xs`, in MATLAB's
/// order, `rstd.sort`'s). Sorts a copy; `xs` is untouched. NaN if `xs`
/// holds one, as MATLAB's.
pub fn median(allocator: Allocator, comptime R: type, xs: anytype) MedianError!R {
    const E = Elem(@TypeOf(xs));
    comptime checkPair(R, E);
    if (xs.len == 0) return error.EmptyInput;
    const s = try sortedCopy(allocator, E, xs);
    defer allocator.free(s);
    const n = s.len;
    if (comptime isComplex(E)) {
        const F = Real(R);
        // a NaN, sorted last, is the median
        const last = s[n - 1];
        if (std.math.isNan(last.re) or std.math.isNan(last.im)) return .{ .re = to(F, last.re), .im = to(F, last.im) };
        const hi = s[n / 2];
        if (n % 2 == 1) return .{ .re = to(F, hi.re), .im = to(F, hi.im) };
        const lo = s[n / 2 - 1];
        const t = @as(Wide(F), 0.5);
        return .{ .re = lerp(F, to(F, lo.re), to(F, hi.re), t), .im = lerp(F, to(F, lo.im), to(F, hi.im), t) };
    }
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
/// `math.floor`, as `@floor` keeps fractions on x87 (issue #481).
fn quantileOf(comptime R: type, comptime E: type, s: []const E, p: R, method: QuantileMethod) R {
    const W = Wide(R);
    const n: W = @floatFromInt(s.len);
    const q: W = p;
    const last = s.len - 1;
    switch (method) {
        .inverted_cdf, .closest_observation => {
            // the discontinuous ones pick an order statistic
            const vi = if (method == .inverted_cdf) n * q - 1 else n * q - 1 - 0.5;
            const prev = math.floor(vi);
            const exact = vi == prev;
            const take_prev = if (method == .inverted_cdf) exact else exact and @mod(prev, 2) == 0;
            const k = if (take_prev) prev else prev + 1;
            return to(R, s[index(W, k, last)]);
        },
        .averaged_inverted_cdf => {
            const vi = n * q - 1;
            const prev = math.floor(vi);
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
            const prev = math.floor(vi);
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

/// `xs` sorted ascending by `rstd.sort` (NaN last, complex in MATLAB's
/// order), owned by the caller.
fn sortedCopy(allocator: Allocator, comptime E: type, xs: anytype) Allocator.Error![]E {
    const s = try allocator.alloc(E, xs.len);
    for (s, xs) |*d, x| d.* = x;
    sort.sort(E, s, .ascend);
    return s;
}

/// The quantiles' elements: integers and floats, as MATLAB's `quantile`.
fn checkOrdered(comptime E: type) void {
    if (isComplex(E)) @compileError("rstd.stats: quantiles of complex " ++ @typeName(E) ++ " are not defined");
    checkReal(E);
}
