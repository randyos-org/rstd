//! The complex `@Vector` extremes: the largest and smallest magnitudes and their indices.

const std = @import("std");
const types = @import("../../types.zig");
const util = @import("../../util.zig");
const extremum = @import("../../blas1/extremum.zig");
const shared = @import("../util.zig");
const Real = @import("../real.zig").Real;
const complex = @import("../complex.zig");

const BlasLong = types.BlasLong;
const Find = extremum.Find;
const lanes = shared.lanes;
const abs = shared.abs;

pub fn Extremes(comptime F: type, comptime U: comptime_int) type {
    const Base = complex.Complex(F, U);
    return struct {
        const V = Base.V;
        const W = Base.W;
        const ld = Base.ld;

        /// The even (`odd` false) or odd floats of `a`, then of `b`.
        inline fn parts(a: V, b: V, comptime odd: bool) V {
            const mask = comptime blk: {
                var m: [W]i32 = undefined;
                for (0..W) |k| {
                    const s = 2 * k + @intFromBool(odd);
                    m[k] = if (s < W) @intCast(s) else ~@as(i32, @intCast(s - W));
                }
                break :blk m;
            };
            return @shuffle(F, a, b, mask);
        }

        /// The extremes of `|re| + |im|`, as `Real`'s. The index falls back
        /// to the portable kernel for an infinite key too, which may be a
        /// finite element's that overflowed.
        fn Extreme(comptime find: Find) type {
            return struct {
                const P = extremum.Extremum(F, 2, find);
                const S = Real(F, W / lanes(F));
                const Key = struct {
                    // rounded: `search` finds the block's stored key again
                    pub inline fn key(p: [*]const F, s: usize, i: usize) F {
                        return util.rounded(F, abs(p[i * s]) + abs(p[i * s + 1]));
                    }
                    pub inline fn block(p: [*]const F, s: usize, i: usize, chk: *V) F {
                        // `4 * W` inlined keys: 256 at x280's 512 bits
                        @setEvalBranchQuota(100_000);
                        var m: [4]V = undefined;
                        inline for (0..4) |q| {
                            if (s == 2) {
                                const a = abs(ld(p + 2 * (i + q * W)));
                                const b = abs(ld(p + 2 * (i + q * W) + W));
                                m[q] = parts(a, b, false) + parts(a, b, true);
                            } else {
                                var r: [W]F = undefined;
                                inline for (0..W) |k| r[k] = key(p, s, i + q * W + k);
                                m[q] = r;
                            }
                            chk.* += m[q];
                        }
                        return util.rounded(F, S.pick(find, m));
                    }
                };

                pub fn search(n: BlasLong, x: [*]F, incx: BlasLong) ?struct { usize, F } {
                    return S.search(@intCast(n), x, 2 * @as(usize, @intCast(incx)), find, Key);
                }

                pub fn value(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
                    if (n <= 0 or incx <= 0) return 0;
                    const r = search(n, x, incx) orelse return P.value(n, x, incx);
                    return r[1];
                }

                pub fn index(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) BlasLong {
                    if (n <= 0 or incx <= 0) return 0;
                    const r = search(n, x, incx) orelse return P.index(n, x, incx);
                    if (std.math.isInf(r[1])) return P.index(n, x, incx);
                    return @intCast(r[0] + 1);
                }
            };
        }

        pub const amax = Extreme(.{ .abs = true, .max = true }).value;
        pub const iamax = Extreme(.{ .abs = true, .max = true }).index;
        pub const amin = Extreme(.{ .abs = true, .max = false }).value;
        pub const iamin = Extreme(.{ .abs = true, .max = false }).index;
    };
}
