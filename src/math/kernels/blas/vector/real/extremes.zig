//! The real `@Vector` extremes: the largest and smallest values and magnitudes, and their indices.

const std = @import("std");
const types = @import("../../types.zig");
const extremum = @import("../../blas1/extremum.zig");
const shared = @import("../util.zig");
const real = @import("../real.zig");

const BlasLong = types.BlasLong;
const Find = extremum.Find;
const abs = shared.abs;
const hasMinMax = shared.hasMinMax;

pub fn Extremes(comptime F: type, comptime U: comptime_int) type {
    const Base = real.Real(F, U);
    return struct {
        const V = Base.V;
        const W = Base.W;
        const ld = Base.ld;

        /// The extremes by `find`, value and 1-based index, as the portable
        /// kernels (`extremum.zig`), to which they fall back for a NaN.
        fn Extreme(comptime find: Find) type {
            return struct {
                const P = extremum.Extremum(F, 1, find);
                const Key = struct {
                    pub inline fn key(p: [*]const F, s: usize, i: usize) F {
                        return if (find.abs) abs(p[i * s]) else p[i * s];
                    }
                    pub inline fn block(p: [*]const F, s: usize, i: usize, chk: *V) F {
                        var m: [4]V = undefined;
                        inline for (0..4) |q| {
                            if (s == 1) {
                                m[q] = ld(p + i + q * W);
                            } else {
                                var r: [W]F = undefined;
                                inline for (0..W) |k| r[k] = p[(i + q * W + k) * s];
                                m[q] = r;
                            }
                            if (find.abs) m[q] = abs(m[q]);
                            chk.* += if (find.abs) m[q] else abs(m[q]);
                        }
                        return pick(find, m);
                    }
                };

                pub fn value(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
                    if (n <= 0 or incx <= 0) return 0;
                    const r = search(@intCast(n), x, @intCast(incx), find, Key) orelse return P.value(n, x, incx);
                    return r[1];
                }

                pub fn index(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) BlasLong {
                    if (n <= 0 or incx <= 0) return 0;
                    const r = search(@intCast(n), x, @intCast(incx), find, Key) orelse return P.index(n, x, incx);
                    return @intCast(r[0] + 1);
                }
            };
        }

        pub const amax = Extreme(.{ .abs = true, .max = true }).value;
        pub const iamax = Extreme(.{ .abs = true, .max = true }).index;
        pub const amin = Extreme(.{ .abs = true, .max = false }).value;
        pub const iamin = Extreme(.{ .abs = true, .max = false }).index;
        pub const max = Extreme(.{ .abs = false, .max = true }).value;
        pub const imax = Extreme(.{ .abs = false, .max = true }).index;
        pub const min = Extreme(.{ .abs = false, .max = false }).value;
        pub const imin = Extreme(.{ .abs = false, .max = false }).index;

        /// The extreme of four vectors of keys: by `@max` or `@min` where
        /// they are instructions (`hasMinMax`), else by compare and select,
        /// which agrees on any keys but NaN, for which `search` gives up.
        pub inline fn pick(comptime find: Find, m: [4]V) F {
            if (comptime !hasMinMax(F)) {
                const a: [W]F = sel(find, sel(find, m[0], m[1]), sel(find, m[2], m[3]));
                var r = a[0];
                inline for (a[1..]) |k| r = sel(find, k, r);
                return r;
            }
            if (find.max) return @reduce(.Max, @max(@max(m[0], m[1]), @max(m[2], m[3])));
            return @reduce(.Min, @min(@min(m[0], m[1]), @min(m[2], m[3])));
        }

        /// `a` where it is `find`'s better than `b`, else `b`.
        inline fn sel(comptime find: Find, a: anytype, b: @TypeOf(a)) @TypeOf(a) {
            const T = @TypeOf(a);
            const take = if (find.max) a > b else a < b;
            return if (@typeInfo(T) == .vector) @select(F, take, a, b) else if (take) a else b;
        }

        inline fn better(comptime find: Find, k: F, best: F) bool {
            return if (find.max) k > best else k < best;
        }

        /// `find`'s first extreme of `len` keys at a stride of `s` floats,
        /// `Key.key` giving one and `Key.block` the extreme of the `4 * W`
        /// from `i`, adding their magnitudes to `chk`, which only a NaN key
        /// makes NaN. One pass keeps the first block with the extreme, no
        /// index kept per element; only that block is searched for the
        /// index. Its index from 0 and its key, or null for a NaN key.
        pub fn search(len: usize, x: [*]const F, s: usize, comptime find: Find, comptime Key: type) ?struct { usize, F } {
            const B = 4 * W;
            var best = Key.key(x, s, 0);
            if (best != best) return null;
            var pos: usize = 0;
            var i: usize = 0;
            if (len >= B) {
                var chk: V = @splat(0);
                var blk: usize = 0;
                while (i + B <= len) : (i += B) {
                    const m = Key.block(x, s, i, &chk);
                    if (better(find, m, best)) {
                        best = m;
                        blk = i;
                    }
                }
                if (std.math.isNan(@reduce(.Add, chk))) return null;
                pos = blk;
                while (Key.key(x, s, pos) != best) pos += 1;
            }
            while (i < len) : (i += 1) {
                const k = Key.key(x, s, i);
                if (k != k) return null;
                if (better(find, k, best)) {
                    best = k;
                    pos = i;
                }
            }
            return .{ pos, Key.key(x, s, pos) };
        }
    };
}
