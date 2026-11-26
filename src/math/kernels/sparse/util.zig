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

/// Whether a column-major `rows x cols` block with leading dimension `ld`
/// fits in `len` elements.
pub fn denseFits(len: usize, ld: usize, rows: usize, cols: usize) bool {
    if (ld < @max(rows, 1)) return false;
    if (rows == 0 or cols == 0) return true;
    return (cols - 1) * ld + rows <= len;
}

/// Column `j` of such a block, empty when `rows` is 0.
pub fn denseCol(a: anytype, ld: usize, j: usize, rows: usize) @TypeOf(a) {
    return if (rows == 0) a[0..0] else a[j * ld ..][0..rows];
}
