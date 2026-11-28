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
                    pub inline fn keys(p: [*]const F, s: usize, i: usize, chk: *V) [4]V {
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
                        return m;
                    }
                    pub inline fn extreme(m: [4]V) F {
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

        /// The keys of `a` or `b` that are `find`'s better, lane by lane.
        inline fn keep(comptime find: Find, a: V, b: V) V {
            if (comptime !hasMinMax(F)) return sel(find, a, b);
            return if (find.max) @max(a, b) else @min(a, b);
        }

        /// `find`'s first extreme of `len` keys at a stride of `s` floats:
        /// `Key.key` gives one, `Key.keys` the `4 * W` from `i` as vectors,
        /// adding their magnitudes to `chk`, which only a NaN key makes NaN,
        /// and `Key.extreme` the extreme of such vectors. The blocks keep the
        /// extremes lane by lane, reduced once at the end, as
        /// `kernel/x86_64/iamax_sse.S` keeps them in registers; a second pass
        /// finds the first block holding the extreme, and the key in it
        /// (issue #846). Its index from 0 and its key, or null for a NaN key.
        pub fn search(len: usize, x: [*]const F, s: usize, comptime find: Find, comptime Key: type) ?struct { usize, F } {
            // unit strides apart (2 floats for complex), so their loop holds no gather
            return switch (s) {
                1 => searchAt(len, x, 1, find, Key),
                2 => searchAt(len, x, 2, find, Key),
                else => searchAt(len, x, s, find, Key),
            };
        }

        inline fn searchAt(len: usize, x: [*]const F, s: usize, comptime find: Find, comptime Key: type) ?struct { usize, F } {
            const B = 4 * W;
            var best = Key.key(x, s, 0);
            if (best != best) return null;
            const blocks = len / B * B;
            if (blocks > 0) {
                var chk: V = @splat(0);
                var acc = Key.keys(x, s, 0, &chk);
                var i: usize = B;
                while (i < blocks) : (i += B) {
                    const m = Key.keys(x, s, i, &chk);
                    inline for (0..4) |q| acc[q] = keep(find, m[q], acc[q]);
                }
                if (std.math.isNan(@reduce(.Add, chk))) return null;
                best = Key.extreme(acc);
            }
            var pos: ?usize = if (blocks == 0) 0 else null;
            for (blocks..len) |i| {
                const k = Key.key(x, s, i);
                if (k != k) return null;
                if (better(find, k, best)) {
                    best = k;
                    pos = i;
                }
            }
            if (pos) |p| return .{ p, best };
            // the first block with the extreme, then its first key
            const want: V = @splat(best);
            var chk: V = @splat(0);
            var i: usize = 0;
            while (i < blocks) : (i += B) {
                const m = Key.keys(x, s, i, &chk);
                // the lanes' extremes hold `best` where any key does
                const e = keep(find, keep(find, m[0], m[1]), keep(find, m[2], m[3]));
                if (@reduce(.Or, e == want)) break;
            }
            // from 0 where the vectors' keys round otherwise than `Key.key`'s
            var p: usize = if (i < blocks) i else 0;
            while (p < len and Key.key(x, s, p) != best) p += 1;
            if (p == len) return null;
            return .{ p, best };
        }
    };
}
