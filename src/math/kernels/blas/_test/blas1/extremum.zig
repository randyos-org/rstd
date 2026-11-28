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

/// Index slot `i` of routine `r` on `n` keys whose extreme is at both `a`
/// and `b` (`a < b`), the rest a field of lesser keys: the first, `a`, wins
/// wherever the two fall in the vector kernels' blocks and tail (issue #846).
fn checkTie(vt: *const Vtable, comptime r: anytype, comptime i: usize, n: usize, a: usize, b: usize) !void {
    const F = extremum.precisions[i][0];
    const w = extremum.precisions[i][1];
    const find: extremum.Find = r[0];
    const kern = @field(vt, @tagName(r[2][i])).?;
    var x: [2 * 257]F = undefined;
    for (0..n) |k| {
        const t: F = @sin(@as(F, @floatFromInt(k)));
        // keys in [-1, 1] (or [1, 2] in size for `amin`)
        const v: F = if (find.abs and !find.max) 1 + @abs(t) else t;
        if (w == 2) {
            x[2 * k] = v / 2;
            x[2 * k + 1] = -v / 2;
        } else x[k] = v;
    }
    // the extreme twice, by another sign or order of parts at `b`
    const e: F = if (find.abs) (if (find.max) 7 else 0.25) else (if (find.max) 7 else -7);
    for ([_]usize{ a, b }, [_]F{ 1, -1 }) |at, sign| {
        if (w == 2) {
            x[2 * at] = if (sign > 0) e * 3 / 7 else -e * 4 / 7;
            x[2 * at + 1] = if (sign > 0) -e * 4 / 7 else e * 3 / 7;
        } else x[at] = if (find.abs) sign * e else e;
    }
    const got = kern(@intCast(n), &x, 1);
    if (got != @as(BlasLong, @intCast(a + 1))) {
        check.say("{s} n {d}: ties at {d} and {d}, index {d}\n", .{ @tagName(r[2][i]), n, a, b, got - 1 });
        return error.TestUnexpectedResult;
    }
}

fn checkTies(vt: *const Vtable) !void {
    var failed: ?anyerror = null;
    inline for (extremum.routines) |r| if (r[1]) inline for (.{ 0, 1, 3, 4 }) |i| {
        if (i < r[2].len) for ([_]usize{ 64, 130, 257 }) |n| for ([_][2]usize{ .{ 0, 1 }, .{ 3, 100 }, .{ 70, 129 }, .{ 1, 256 }, .{ 200, 256 } }) |ab| {
            if (ab[1] < n) check.keep(&failed, checkTie(vt, r, i, n, ab[0], ab[1]));
        };
    };
    if (failed) |e| return e;
}

test "every linked core's index extremes take the first of equal extremes, in any block or the tail" {
    try check.eachCore(checkTies);
}
