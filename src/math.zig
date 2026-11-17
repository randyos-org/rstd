//! Scalar math: exact `floor`, `ceil` and `round` where the builtins are not.

const std = @import("std");
const sysinfo = @import("builtin");

/// `f64` math on x87: compiler-rt's `f64` `@floor`, `@ceil` and `@round`
/// keep fraction bits there (issue #264); `@trunc` masks bits, so it is exact.
pub const x87_f64 = sysinfo.target.cpu.arch == .x86 and !std.Target.x86.featureSetHas(sysinfo.cpu.features, .sse2);

/// `@floor`, exact on x87 too (issue #264).
pub fn floor(x: anytype) @TypeOf(x) {
    if (comptime !x87_f64 or @TypeOf(x) != f64) return @floor(x);
    const t = @trunc(x);
    return if (x < t) t - 1 else t;
}

/// `@ceil`, exact on x87 too (issue #264).
pub fn ceil(x: anytype) @TypeOf(x) {
    if (comptime !x87_f64 or @TypeOf(x) != f64) return @ceil(x);
    const t = @trunc(x);
    return if (x > t) t + 1 else t;
}

/// `@round` (halves away from zero), exact on x87 too (issue #264).
pub fn round(x: anytype) @TypeOf(x) {
    if (comptime !x87_f64 or @TypeOf(x) != f64) return @round(x);
    const t = @trunc(x);
    if (@abs(x - t) < 0.5) return t;
    return if (x < 0) t - 1 else t + 1;
}

pub fn Fraction(comptime F: type) type {
    return struct { numerator: F, denominator: F };
}

/// MATLAB's `rat` (`rat.m`): `x` as a fraction from its continued fraction,
/// within `max(tol, eps(x))`; a null `tol` is MATLAB's `1e-6 * |x|`. The
/// denominator is positive; ±inf gives ±1/0 and NaN 0/0. Ported from
/// matlab88's `RAT` (issue #331).
pub fn rationalApprox(x: anytype, tol: ?@TypeOf(x)) Fraction(@TypeOf(x)) {
    const F = @TypeOf(x);
    if (@typeInfo(F) != .float) @compileError("rationalApprox: " ++ @typeName(F) ++ " is not a float type");
    if (std.math.isNan(x)) return .{ .numerator = 0, .denominator = 0 };
    if (std.math.isInf(x)) return .{ .numerator = std.math.sign(x), .denominator = 0 };
    const ax = @abs(x);
    const eps = std.math.nextAfter(F, ax, std.math.inf(F)) - ax;
    const limit = @max(tol orelse (1e-6 * ax), eps);
    // the last two convergents: n / d and the one before
    var n: F = 1;
    var d: F = 0;
    var n_prev: F = 0;
    var d_prev: F = 1;
    var y = x;
    while (true) {
        const term = round(y);
        y -= term;
        const n_next = n * term + n_prev;
        const d_next = d * term + d_prev;
        n_prev = n;
        d_prev = d;
        n = n_next;
        d = d_next;
        // an exact remainder of 0 ends it too, as 1/0 would not
        if (x == 0 or y == 0 or @abs(n / d - x) < limit) break;
        y = 1 / y;
    }
    return if (d < 0) .{ .numerator = -n, .denominator = -d } else .{ .numerator = n, .denominator = d };
}

test "rationalApprox: MATLAB's convergents and special values" {
    const pi = rationalApprox(@as(f64, std.math.pi), null);
    try std.testing.expectEqual(Fraction(f64){ .numerator = 355, .denominator = 113 }, pi);
    try std.testing.expectEqual(Fraction(f64){ .numerator = 22, .denominator = 7 }, rationalApprox(@as(f64, std.math.pi), 0.01));
    try std.testing.expectEqual(Fraction(f64){ .numerator = -3, .denominator = 4 }, rationalApprox(@as(f64, -0.75), null));
    try std.testing.expectEqual(Fraction(f64){ .numerator = 5, .denominator = 1 }, rationalApprox(@as(f64, 5), null));
    try std.testing.expectEqual(Fraction(f64){ .numerator = 0, .denominator = 1 }, rationalApprox(@as(f64, 0), null));
    try std.testing.expectEqual(Fraction(f32){ .numerator = 1, .denominator = 10 }, rationalApprox(@as(f32, 0.1), null));
    try std.testing.expectEqual(Fraction(f64){ .numerator = -1, .denominator = 0 }, rationalApprox(-std.math.inf(f64), null));
    try std.testing.expectEqual(Fraction(f64){ .numerator = 0, .denominator = 0 }, rationalApprox(std.math.nan(f64), null));
    // within eps when tol is 0: the fraction is x itself
    const third = rationalApprox(@as(f64, 1.0 / 3.0), 0);
    try std.testing.expectEqual(@as(f64, 1.0 / 3.0), third.numerator / third.denominator);
}

test "floor, ceil and round are exact, on x87 too" {
    var x: f64 = 22.0 / 7.0;
    _ = &x;
    try std.testing.expectEqual(@as(f64, 3), floor(x));
    try std.testing.expectEqual(@as(f64, 4), ceil(x));
    try std.testing.expectEqual(@as(f64, 3), round(x));
    try std.testing.expectEqual(@as(f64, -4), floor(-x));
    try std.testing.expectEqual(@as(f64, -3), ceil(-x));
    try std.testing.expectEqual(@as(f64, -3), round(-x));
    for ([_]f64{ 0.5, -0.5, 2.5, -2.5, 0.49999999999999994, 1e300, -0.0 }, [_]f64{ 1, -1, 3, -3, 0, 1e300, -0.0 }) |v, want| {
        try std.testing.expectEqual(want, round(v));
    }
}
