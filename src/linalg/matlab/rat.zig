//! RAT, from `matlab88.f90`'s `RAT`.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");

pub const Fraction = struct { numerator: f64, denominator: f64 };

/// RAT: `x` as a continued fraction of up to `d.len` rounded terms, kept in
/// `d`, stopping once the remainder is within `1/max_term`.
pub fn rat(x: f64, max_term: f64, d: []f64) Fraction {
    var z = x;
    var k: usize = 0;
    for (0..d.len) |i| {
        k = i;
        d[k] = @round(z);
        z -= d[k];
        if (@abs(z) * max_term <= 1.0) break;
        z = 1.0 / z;
    }

    var t = d[k];
    var s: f64 = 1.0;
    if (k >= 1) {
        var ib: usize = 1;
        while (ib <= k) : (ib += 1) {
            const i = k - ib;
            const prev_t = t;
            t = d[i] * t + s;
            s = prev_t;
        }
    }
    if (s < 0.0) {
        t = -t;
        s = -s;
    }
    return .{ .numerator = t, .denominator = s };
}

// ---- tests ------------------------------------------------------------------

test "rat recovers a simple exact fraction" {
    var d: [10]f64 = @splat(0.0);
    const f = rat(22.0 / 7.0, 100.0, &d);
    try zigstd.testExpectApproxEqAbs(22.0, f.numerator, 1e-9);
    try zigstd.testExpectApproxEqAbs(7.0, f.denominator, 1e-9);
}

test "rat reconstructs the original value" {
    var d: [10]f64 = @splat(0.0);
    const x = std.math.pi;
    const f = rat(x, 1000.0, &d);
    try zigstd.testExpectApproxEqAbs(x, f.numerator / f.denominator, 1e-6);
}
