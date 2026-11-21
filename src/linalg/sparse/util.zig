//! What the sparse routines share.

const std = @import("std");

/// The first position in ascending `xs` whose value is at least `v`.
pub fn lowerBound(xs: []const u32, v: u64) usize {
    var lo: usize = 0;
    var hi = xs.len;
    while (lo < hi) {
        const mid = lo + (hi - lo) / 2;
        if (xs[mid] < v) lo = mid + 1 else hi = mid;
    }
    return lo;
}

/// `n` as a `u32`, or `error.TooLarge`.
pub fn narrow(n: usize) error{TooLarge}!u32 {
    return std.math.cast(u32, n) orelse error.TooLarge;
}
