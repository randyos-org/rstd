//! The extremum kernels (`blas1/extremum.zig`) against a plain loop: the
//! portable ones and every linked core's.

const std = @import("std");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../check.zig");
const extremum = @import("../../blas1/extremum.zig");

const BlasLong = types.BlasLong;

/// Slot `i` of routine `r` (`extremum.routines`) of `vt`.
fn checkSlot(vt: *const Vtable, comptime r: anytype, comptime i: usize) !void {
    const p = extremum.precisions[i];
    const T = if (p[1] == 2) std.math.Complex(p[0]) else p[0];
    const name = comptime @tagName(r[2][i]);
    try check.checkExtreme(T, r[0], name[0 .. name.len - 2], @field(vt, name).?);
}

/// The single and double slots of every routine.
fn checkAll(vt: *const Vtable) !void {
    var failed: ?anyerror = null;
    inline for (extremum.routines) |r| inline for (.{ 0, 1, 3, 4 }) |i| {
        if (i < r[2].len) check.keep(&failed, checkSlot(vt, r, i));
    };
    if (failed) |e| return e;
}

test "portable extremum kernels match a plain loop" {
    const vt = check.portableVtable();
    inline for (extremum.routines) |r| {
        inline for (.{ 0, 1, 3, 4, 2, 5 }) |i| {
            if (i < r[2].len) try checkSlot(vt, r, i);
        }
    }
}

test "every linked core's extremum kernels match a plain loop, where this machine runs them" {
    try check.eachCore(checkAll);
}

test "portable iamax: the first of equal maxima, 1-based; the first NaN wins" {
    const idamax = extremum.Extremum(f64, 1, extremum.amax).index;
    const izamax = extremum.Extremum(f64, 2, extremum.amax).index;
    var ties = [_]f64{ 1, -3, 3, 2 };
    try std.testing.expectEqual(@as(BlasLong, 2), idamax(4, &ties, 1));
    var nans = [_]f64{ 5, std.math.inf(f64), std.math.nan(f64), 7, std.math.nan(f64) };
    try std.testing.expectEqual(@as(BlasLong, 3), idamax(5, &nans, 1));
    // complex: an infinite part beats a finite element whose |re| + |im| overflows
    const big = std.math.floatMax(f64);
    var z = [_]f64{ big, big, 1, std.math.inf(f64), big, 0 };
    try std.testing.expectEqual(@as(BlasLong, 2), izamax(3, &z, 1));
    z[3] = 1;
    try std.testing.expectEqual(@as(BlasLong, 1), izamax(3, &z, 1));
}
