//! The reversal `J` (`J x` is `x` backward) that `rq.zig` and `ql.zig` build
//! their factorizations on: `J A J` turns an RQ into an LQ, a QL into a QR.

const std = @import("std");
const types = @import("../blas/types.zig");
const ops = @import("ops.zig");

const Mat = ops.Mat;
const Side = types.Side;

/// `dst := J src J`, or `(J src J)ᴴ` with `adj`.
pub fn flip(comptime T: type, comptime adj: bool, src: Mat(T), dst: Mat(T)) void {
    const r = src.rows;
    const c = src.cols;
    for (0..c) |j| for (0..r) |i| {
        const v = src.get(i, j);
        if (adj) dst.at(c - 1 - j, r - 1 - i).* = ops.conj(T, v) else dst.at(r - 1 - i, c - 1 - j).* = v;
    };
}

/// `dst := J src`.
pub fn reverse(comptime T: type, src: []const T, dst: []T) void {
    for (dst, 0..) |*d, i| d.* = src[src.len - 1 - i];
}

/// `C`'s rows (`l`) or columns (`r`) reversed in place.
pub fn reverseAlong(comptime T: type, side: Side, c: Mat(T)) void {
    switch (side) {
        .l => for (0..c.cols) |j| std.mem.reverse(T, c.col(j)[0..c.rows]),
        .r => for (0..c.cols / 2) |j| for (0..c.rows) |i| {
            std.mem.swap(T, c.at(i, j), c.at(i, c.cols - 1 - j));
        },
    }
}
