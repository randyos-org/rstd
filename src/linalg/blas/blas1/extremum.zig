//! The largest or smallest element, or its 1-based index: the portable
//! kernels of `amax`, `amin`, `iamax`, `iamin`, `imax`, `imin`, `max` and
//! `min` (`kernel/arm/{,i}{,z}a{max,min}.c`, `{i,}{max,min}.c`), and their
//! test.
//!
//! The scan starts at element 0 and replaces it only on a strict `>` (`<`),
//! so ties keep the first. With `abs` an element counts by `|x|`, by
//! `|re| + |im|` for complex (`CABS1`, not the modulus). `n <= 0` or
//! `inc_x <= 0` gives 0.

const std = @import("std");
const types = @import("../types.zig");
const util = @import("../util.zig");
const check = @import("check.zig");

const BlasLong = types.BlasLong;

pub const Find = struct {
    /// Compare `|x|` rather than `x`.
    abs: bool,
    /// The largest rather than the smallest.
    max: bool,
};

/// The kernels over elements of `w` floats (2 for complex).
pub fn Extremum(comptime F: type, comptime w: comptime_int, comptime find: Find) type {
    return struct {
        fn key(x: [*]F, at: BlasLong) F {
            if (w == 2) return @abs(util.at(F, x, at).*) + @abs(util.at(F, x, at + 1).*);
            const v = util.at(F, x, at).*;
            return if (find.abs) @abs(v) else v;
        }

        /// The index of the extreme element, from 0; `n >= 1`, `inc_x >= 1`.
        fn scan(n: BlasLong, x: [*]F, inc_x: BlasLong) struct { BlasLong, F } {
            var best: BlasLong = 0;
            var best_key = key(x, 0);
            var i: BlasLong = 1;
            while (i < n) : (i += 1) {
                const k = key(x, i * inc_x * w);
                if (if (find.max) k > best_key else k < best_key) {
                    best = i;
                    best_key = k;
                }
            }
            return .{ best, best_key };
        }

        pub fn value(n: BlasLong, x: [*]F, inc_x: BlasLong) callconv(.c) F {
            if (n <= 0 or inc_x <= 0) return 0;
            return scan(n, x, inc_x)[1];
        }

        pub fn index(n: BlasLong, x: [*]F, inc_x: BlasLong) callconv(.c) BlasLong {
            if (n <= 0 or inc_x <= 0) return 0;
            return scan(n, x, inc_x)[0] + 1;
        }
    };
}

/// `kern` (a `value` or `index` kernel, by its return type) against a plain
/// loop; `T` is the element type, real or complex.
pub fn checkKernel(comptime T: type, comptime find: Find, routine: []const u8, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    for (check.sizes) |n| {
        for ([_]usize{ 1, 2, 3 }) |inc| {
            var x: [check.capacity]T = undefined;
            check.fill(T, &x, n + 1);
            var best: usize = 0;
            var best_key: f64 = undefined;
            for (0..n) |i| {
                const v = x[i * inc];
                const k: f64 = if (comptime complex) @abs(@as(f64, v.re)) + @abs(@as(f64, v.im)) else if (find.abs) @abs(@as(f64, v)) else v;
                if (i == 0 or (if (find.max) k > best_key else k < best_key)) {
                    best = i;
                    best_key = k;
                }
            }
            const got = kern(@intCast(n), @as([*]F, @ptrCast(&x)), @intCast(inc));
            const ok = if (@TypeOf(got) == BlasLong)
                check.expectClose(f64, @floatFromInt(best + 1), @floatFromInt(got), 0)
            else
                check.expectClose(F, @floatCast(best_key), got, check.tolerance(F, 1));
            ok catch |e| return check.case(e, routine, T, n, .{ inc, 0 });
        }
    }
}
