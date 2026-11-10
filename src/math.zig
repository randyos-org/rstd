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
