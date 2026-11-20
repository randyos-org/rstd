//! The real `@Vector` kernels (`vector.zig`).

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../types.zig");
const util = @import("../util.zig");
const extremum = @import("../blas1/extremum.zig");
const shared = @import("util.zig");

const BlasLong = types.BlasLong;
const Find = extremum.Find;
const lanes = shared.lanes;
const splatTo = shared.splatTo;
const at = shared.at;
const abs = shared.abs;
const mulAdd = shared.mulAdd;
const Pick = shared.Pick;
const hasMinMax = shared.hasMinMax;
const root = shared.root;

/// The real kernels over `F`; reductions keep two vectors of sums.
pub fn Real(comptime F: type, comptime U: comptime_int) type {
    return struct {
        pub const W = U * lanes(F);
        const V = @Vector(W, F);
        /// Four strided elements, for four sums.
        const G = @Vector(4, F);

        inline fn ld(p: [*]const F) V {
            return p[0..W].*;
        }
        inline fn st(p: [*]F, v: V) void {
            p[0..W].* = v;
        }
        inline fn gather(p: [*]F, inc: BlasLong, i: BlasLong) G {
            return .{ at(F, p, i * inc).*, at(F, p, (i + 1) * inc).*, at(F, p, (i + 2) * inc).*, at(F, p, (i + 3) * inc).* };
        }

        pub fn axpy(n: BlasLong, _: BlasLong, _: BlasLong, a: F, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F, _: BlasLong) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) {
                const len: usize = @intCast(n);
                const av: V = @splat(a);
                var i: usize = 0;
                while (i + W <= len) : (i += W) st(y + i, mulAdd(V, av, ld(x + i), ld(y + i)));
                while (i < len) : (i += 1) y[i] = mulAdd(F, a, x[i], y[i]);
                return 0;
            }
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const yp = at(F, y, i * incy);
                yp.* = mulAdd(F, a, at(F, x, i * incx).*, yp.*);
            }
            return 0;
        }

        /// `dot`'s, `nrm2`'s and `sum`'s sums: `f64` for `f32`, as the
        /// portable kernels (products of `f32` are exact there), else `F`.
        pub const D = if (F == f32) f64 else F;
        const DG = @Vector(4, D);
        /// `dotSum`'s and `squares`' accumulator width: `W`, or half of it where
        /// their two `D` accumulators would fill more than half the vector
        /// registers (T21: Zen's `sdot` held all 16 ymm in two `@Vector(32, f64)`).
        pub const AW = blk: {
            const bits = lanes(F) * @bitSizeOf(F);
            const acc_regs = 2 * W * @bitSizeOf(D) / bits;
            break :blk if (D != F and W >= 2 and 2 * acc_regs > util.registers(bits, sysinfo.target.cpu)) W / 2 else W;
        };
        const AV = @Vector(AW, D);

        inline fn ldA(p: [*]const F) @Vector(AW, F) {
            return p[0..AW].*;
        }

        pub fn dot(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) F {
            if (n <= 0) return 0;
            return @floatCast(dotSum(n, x, incx, y, incy));
        }

        /// `dsdot`: `dot`'s sum, unrounded (`kernel/riscv64/dot_rvv.c` with
        /// `DSDOT`); for `f32` only.
        pub fn dsdot(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) D {
            if (n <= 0) return 0;
            return dotSum(n, x, incx, y, incy);
        }

        fn dotSum(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) D {
            if (incx == 1 and incy == 1) {
                const len: usize = @intCast(n);
                var acc: [2]AV = @splat(@splat(0));
                var i: usize = 0;
                while (i + 2 * AW <= len) : (i += 2 * AW) {
                    inline for (0..2) |q| acc[q] = mulAdd(AV, ldA(x + i + q * AW), ldA(y + i + q * AW), acc[q]);
                }
                if (i + AW <= len) {
                    acc[0] = mulAdd(AV, ldA(x + i), ldA(y + i), acc[0]);
                    i += AW;
                }
                var s = @reduce(.Add, acc[0] + acc[1]);
                while (i < len) : (i += 1) s = mulAdd(D, x[i], y[i], s);
                return s;
            }
            var acc: DG = @splat(0);
            var i: BlasLong = 0;
            while (i + 4 <= n) : (i += 4) acc = mulAdd(DG, gather(x, incx, i), gather(y, incy, i), acc);
            var s = @reduce(.Add, acc);
            while (i < n) : (i += 1) s = mulAdd(D, at(F, x, i * incx).*, at(F, y, i * incy).*, s);
            return s;
        }

        const Scale = enum { mul, zero, keep };

        /// `a*v`; for a zero `a`, +0, or NaN for NaN and infinity (`keep`).
        inline fn scaled(comptime T: type, comptime how: Scale, a: F, v: T) T {
            if (comptime @typeInfo(T) == .vector) return switch (how) {
                .mul => @as(T, @splat(a)) * v,
                .zero => @splat(0),
                .keep => @select(F, abs(v) <= @as(T, @splat(std.math.floatMax(F))), @as(T, @splat(0)), @as(T, @splat(std.math.nan(F)))),
            };
            return switch (how) {
                .mul => a * v,
                .zero => 0,
                .keep => if (std.math.isFinite(v)) 0 else std.math.nan(F),
            };
        }

        fn scalBy(comptime how: Scale, a: F, len: usize, x: [*]F, sx: usize) void {
            var i: usize = 0;
            if (sx == 1) {
                while (i + W <= len) : (i += W) st(x + i, scaled(V, how, a, ld(x + i)));
            }
            while (i < len) : (i += 1) x[i * sx] = scaled(F, how, a, x[i * sx]);
        }

        /// A zero `alpha` stores +0, or with `flag` 1 NaN for NaN and infinity.
        pub fn scal(n: BlasLong, _: BlasLong, _: BlasLong, a: F, x: [*]F, incx: BlasLong, _: [*]F, _: BlasLong, _: [*]F, flag: BlasLong) callconv(.c) i32 {
            if (n <= 0 or incx <= 0) return 0;
            const how: Scale = if (a != 0) .mul else if (flag == 0) .zero else .keep;
            switch (how) {
                inline else => |h| scalBy(h, a, @intCast(n), x, @intCast(incx)),
            }
            return 0;
        }

        /// The plain sum of squares of `x[0..len]`, in `D`.
        pub fn squares(len: usize, x: [*]F) D {
            var acc: [2]AV = @splat(@splat(0));
            var i: usize = 0;
            while (i + 2 * AW <= len) : (i += 2 * AW) {
                inline for (0..2) |q| {
                    const v: AV = ldA(x + i + q * AW);
                    acc[q] = mulAdd(AV, v, v, acc[q]);
                }
            }
            var s = @reduce(.Add, acc[0] + acc[1]);
            while (i < len) : (i += 1) s = mulAdd(D, x[i], x[i], s);
            return s;
        }

        pub fn nrm2(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            if (n <= 0 or incx == 0) return 0;
            if (n == 1) return abs(x[0]);
            if (incx == 1) return root(F, 1, squares(@intCast(n), x), n, x, incx);
            var acc: DG = @splat(0);
            var i: BlasLong = 0;
            while (i + 4 <= n) : (i += 4) {
                const v: DG = gather(x, incx, i);
                acc = mulAdd(DG, v, v, acc);
            }
            var s = @reduce(.Add, acc);
            while (i < n) : (i += 1) {
                const v: D = at(F, x, i * incx).*;
                s = mulAdd(D, v, v, s);
            }
            return root(F, 1, s, n, x, incx);
        }

        /// The sum of `|x[0..len]|` in `D`, in `plainSum`'s order (`f64` for
        /// `f32`, as `sum` and `dot`: issue #572).
        pub fn absSum(len: usize, x: [*]F) D {
            var acc: [2]AV = @splat(@splat(0));
            var i: usize = 0;
            while (i + 2 * AW <= len) : (i += 2 * AW) {
                inline for (0..2) |q| acc[q] += abs(@as(AV, ldA(x + i + q * AW)));
            }
            var s = @reduce(.Add, acc[0] + acc[1]);
            while (i < len) : (i += 1) s += abs(@as(D, x[i]));
            return s;
        }

        pub fn asum(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            if (n <= 0 or incx <= 0) return 0;
            if (incx == 1) return @floatCast(absSum(@intCast(n), x));
            var acc: DG = @splat(0);
            var i: BlasLong = 0;
            while (i + 4 <= n) : (i += 4) acc += abs(@as(DG, gather(x, incx, i)));
            var s = @reduce(.Add, acc);
            while (i < n) : (i += 1) s += abs(@as(D, at(F, x, i * incx).*));
            return @floatCast(s);
        }

        /// The sum of `x[0..len]` in `D`, as `dot` sums: float `i` of each
        /// `2 * AW` into lane `i % AW` of sum `i / AW % 2`; then the two lane by
        /// lane, the lanes in order, and the floats left in order.
        pub fn plainSum(len: usize, x: [*]F) D {
            var acc: [2]AV = @splat(@splat(0));
            var i: usize = 0;
            while (i + 2 * AW <= len) : (i += 2 * AW) {
                inline for (0..2) |q| acc[q] += @as(AV, ldA(x + i + q * AW));
            }
            var s = @reduce(.Add, acc[0] + acc[1]);
            while (i < len) : (i += 1) s += x[i];
            return s;
        }

        /// Summed in `D`, rounded once; strided, element `i` into lane
        /// `i % 4` of four sums, then the lanes and the rest in order.
        pub fn sum(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            if (n <= 0 or incx <= 0) return 0;
            if (incx == 1) return @floatCast(plainSum(@intCast(n), x));
            var acc: DG = @splat(0);
            var i: BlasLong = 0;
            while (i + 4 <= n) : (i += 4) acc += @as(DG, gather(x, incx, i));
            var s = @reduce(.Add, acc);
            while (i < n) : (i += 1) s += at(F, x, i * incx).*;
            return @floatCast(s);
        }

        /// The extremes by `find`, value and 1-based index, as the portable
        /// kernels (`extremum.zig`), to which they fall back for a NaN.
        fn Extreme(comptime find: Find) type {
            return struct {
                const P = extremum.Extremum(F, 1, find);
                const Key = struct {
                    inline fn key(p: [*]const F, s: usize, i: usize) F {
                        return if (find.abs) abs(p[i * s]) else p[i * s];
                    }
                    inline fn block(p: [*]const F, s: usize, i: usize, chk: *V) F {
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

        pub fn copy(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) {
                const len: usize = @intCast(n);
                var i: usize = 0;
                while (i + W <= len) : (i += W) st(y + i, ld(x + i));
                while (i < len) : (i += 1) y[i] = x[i];
                return 0;
            }
            var i: BlasLong = 0;
            while (i < n) : (i += 1) at(F, y, i * incy).* = at(F, x, i * incx).*;
            return 0;
        }

        pub fn swap(n: BlasLong, _: BlasLong, _: BlasLong, _: F, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F, _: BlasLong) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) {
                const len: usize = @intCast(n);
                var i: usize = 0;
                while (i + W <= len) : (i += W) {
                    const t = ld(x + i);
                    st(x + i, ld(y + i));
                    st(y + i, t);
                }
                while (i < len) : (i += 1) util.swap(F, &x[i], &y[i]);
                return 0;
            }
            var i: BlasLong = 0;
            while (i < n) : (i += 1) util.swap(F, at(F, x, i * incx), at(F, y, i * incy));
            return 0;
        }

        /// `x := c*x + s*y`, `y := c*y - s*x`.
        pub fn rot(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, c: F, s: F) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) {
                const len: usize = @intCast(n);
                const cv: V = @splat(c);
                const sv: V = @splat(s);
                var i: usize = 0;
                while (i + W <= len) : (i += W) {
                    const xv = ld(x + i);
                    const yv = ld(y + i);
                    st(x + i, mulAdd(V, cv, xv, sv * yv));
                    st(y + i, mulAdd(V, cv, yv, -(sv * xv)));
                }
                while (i < len) : (i += 1) rot1(&x[i], &y[i], c, s);
                return 0;
            }
            var i: BlasLong = 0;
            while (i < n) : (i += 1) rot1(at(F, x, i * incx), at(F, y, i * incy), c, s);
            return 0;
        }

        inline fn rot1(xp: *F, yp: *F, c: F, s: F) void {
            const xv = xp.*;
            const yv = yp.*;
            xp.* = mulAdd(F, c, xv, s * yv);
            yp.* = mulAdd(F, c, yv, -(s * xv));
        }

        /// `y += alpha*A*x`, four columns a pass over a unit-stride `y`.
        pub fn gemv_n(m: BlasLong, n: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F) callconv(.c) i32 {
            if (m <= 0 or n <= 0) return 0;
            const rows: usize = @intCast(m);
            const cols: usize = @intCast(n);
            const ld_a: usize = @intCast(lda);
            var j: usize = 0;
            if (incy == 1) {
                while (j + 4 <= cols) : (j += 4) {
                    var t: [4]F = undefined;
                    var c: [4][*]F = undefined;
                    inline for (0..4) |q| {
                        t[q] = alpha * at(F, x, @as(BlasLong, @intCast(j + q)) * incx).*;
                        c[q] = a + (j + q) * ld_a;
                    }
                    var i: usize = 0;
                    while (i + W <= rows) : (i += W) {
                        var yv = ld(y + i);
                        inline for (0..4) |q| yv = mulAdd(V, ld(c[q] + i), @splat(t[q]), yv);
                        st(y + i, yv);
                    }
                    while (i < rows) : (i += 1) {
                        inline for (0..4) |q| y[i] = mulAdd(F, c[q][i], t[q], y[i]);
                    }
                }
            }
            while (j < cols) : (j += 1) {
                const t = alpha * at(F, x, @as(BlasLong, @intCast(j)) * incx).*;
                const c = a + j * ld_a;
                var i: usize = 0;
                if (incy == 1) {
                    while (i + W <= rows) : (i += W) st(y + i, mulAdd(V, ld(c + i), @splat(t), ld(y + i)));
                }
                while (i < rows) : (i += 1) {
                    const yp = at(F, y, @as(BlasLong, @intCast(i)) * incy);
                    yp.* = mulAdd(F, c[i], t, yp.*);
                }
            }
            return 0;
        }

        /// `y += alpha*A'*x`, four columns' sums a pass over a unit-stride
        /// `x`. With no rows, `y += alpha*0`, as the portable kernel.
        pub fn gemv_t(m: BlasLong, n: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F) callconv(.c) i32 {
            if (n <= 0) return 0;
            const rows: usize = @intCast(@max(m, 0));
            const cols: usize = @intCast(n);
            const ld_a: usize = @intCast(lda);
            var j: usize = 0;
            if (incx == 1) {
                while (j + 4 <= cols) : (j += 4) {
                    var c: [4][*]F = undefined;
                    var acc: [4]V = @splat(@splat(0));
                    inline for (0..4) |q| c[q] = a + (j + q) * ld_a;
                    var i: usize = 0;
                    while (i + W <= rows) : (i += W) {
                        const xv = ld(x + i);
                        inline for (0..4) |q| acc[q] = mulAdd(V, ld(c[q] + i), xv, acc[q]);
                    }
                    inline for (0..4) |q| {
                        var s = @reduce(.Add, acc[q]);
                        var r = i;
                        while (r < rows) : (r += 1) s = mulAdd(F, c[q][r], x[r], s);
                        const yp = at(F, y, @as(BlasLong, @intCast(j + q)) * incy);
                        yp.* = mulAdd(F, alpha, s, yp.*);
                    }
                }
            }
            while (j < cols) : (j += 1) {
                const c = a + j * ld_a;
                var s: F = 0;
                var i: usize = 0;
                if (incx == 1) {
                    var acc: V = @splat(0);
                    while (i + W <= rows) : (i += W) acc = mulAdd(V, ld(c + i), ld(x + i), acc);
                    s = @reduce(.Add, acc);
                }
                while (i < rows) : (i += 1) s = mulAdd(F, c[i], at(F, x, @as(BlasLong, @intCast(i)) * incx).*, s);
                const yp = at(F, y, @as(BlasLong, @intCast(j)) * incy);
                yp.* = mulAdd(F, alpha, s, yp.*);
            }
            return 0;
        }

        /// `y := alpha*x + beta*y`, `kernel/riscv64/axpby_rvv.c` and
        /// `axpby_vector_v2.c`: `beta*y`, then `alpha*x` fused in. A zero
        /// `alpha` does not read `x`, a zero `beta` not `y`.
        pub fn axpby(n: BlasLong, alpha: F, x: [*]F, incx: BlasLong, beta: F, y: [*]F, incy: BlasLong) callconv(.c) i32 {
            if (n <= 0) return 0;
            switch (Pick.of(alpha != 0, beta != 0)) {
                inline else => |p| axpbyBy(p, @intCast(n), alpha, x, incx, beta, y, incy),
            }
            return 0;
        }

        fn axpbyBy(comptime p: Pick, len: usize, a: F, x: [*]F, incx: BlasLong, b: F, y: [*]F, incy: BlasLong) void {
            var i: usize = 0;
            if (incx == 1 and incy == 1) {
                while (i + W <= len) : (i += W) st(y + i, Pick.apply(p, V, a, b, x + i, y + i));
            }
            while (i < len) : (i += 1) {
                const k: BlasLong = @intCast(i);
                const yp = at(F, y, k * incy);
                yp.* = Pick.apply(p, F, a, b, at(F, x, k * incx), yp);
            }
        }

        /// `x := h11*x + h12*y`, `y := h21*x + h22*y`, `kernel/riscv64/rotm_rvv.c`:
        /// `param[0]` -2 the identity, -1 the full matrix (`param[1..5]` =
        /// h11, h21, h12, h22), 0 a unit diagonal, 1 `[h11 1; -1 h22]`, fused
        /// as its `vfmacc` and `vfmsac`. A negative stride starts from the
        /// far end.
        pub fn rotm(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, param: [*]F) callconv(.c) i32 {
            const flag = param[0];
            if (n <= 0 or flag == -2) return 0;
            const h = [4]F{ param[1], param[3], param[2], param[4] };
            const form: Rotm = if (flag < 0) .full else if (flag == 0) .unit else .off;
            switch (form) {
                inline else => |f| rotmBy(f, n, x, incx, y, incy, h),
            }
            return 0;
        }

        const Rotm = enum { full, unit, off };

        inline fn rotm1(comptime f: Rotm, comptime T: type, h: [4]F, w: T, z: T) [2]T {
            return switch (f) {
                .full => .{ mulAdd(T, splatTo(T, h[1]), z, splatTo(T, h[0]) * w), mulAdd(T, splatTo(T, h[3]), z, splatTo(T, h[2]) * w) },
                .unit => .{ mulAdd(T, splatTo(T, h[1]), z, w), mulAdd(T, splatTo(T, h[2]), w, z) },
                .off => .{ mulAdd(T, splatTo(T, h[0]), w, z), mulAdd(T, splatTo(T, h[3]), z, -w) },
            };
        }

        fn rotmBy(comptime f: Rotm, n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, h: [4]F) void {
            var i: BlasLong = 0;
            if (incx == 1 and incy == 1) {
                const len: usize = @intCast(n);
                var e: usize = 0;
                while (e + W <= len) : (e += W) {
                    const r = rotm1(f, V, h, ld(x + e), ld(y + e));
                    st(x + e, r[0]);
                    st(y + e, r[1]);
                }
                i = @intCast(e);
            }
            const ix: BlasLong = if (incx < 0) (1 - n) * incx else 0;
            const iy: BlasLong = if (incy < 0) (1 - n) * incy else 0;
            while (i < n) : (i += 1) {
                const px = at(F, x, ix + i * incx);
                const py = at(F, y, iy + i * incy);
                const r = rotm1(f, F, h, px.*, py.*);
                px.* = r[0];
                py.* = r[1];
            }
        }

        /// `y += alpha*S*x`, `S` symmetric from one stored triangle, over the
        /// columns `offset` selects: the last `offset` for the upper
        /// triangle, the first for the lower. `kernel/riscv64/symv_L_rvv.c`,
        /// `symv_U_rvv.c`, `symv_L_vector.c`, `symv_U_vector.c`: per column
        /// `j`, `alpha*x[j]` times the column into `y`, its dot with `x`
        /// times `alpha` into `y[j]`; the diagonal first.
        fn Symv(comptime lower: bool) type {
            return struct {
                fn k(m: BlasLong, offset: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F) callconv(.c) i32 {
                    if (m <= 0 or offset <= 0) return 0;
                    const ld_a: usize = @intCast(lda);
                    const unit = incx == 1 and incy == 1;
                    var j: BlasLong = if (lower) 0 else m - offset;
                    const end: BlasLong = if (lower) offset else m;
                    while (j < end) : (j += 1) {
                        const col = a + @as(usize, @intCast(j)) * ld_a;
                        const t = alpha * at(F, x, j * incx).*;
                        const yj = at(F, y, j * incy);
                        yj.* = mulAdd(F, t, col[@intCast(j)], yj.*);
                        var i: BlasLong = if (lower) j + 1 else 0;
                        const hi: BlasLong = if (lower) m else j;
                        var s: F = 0;
                        if (unit) {
                            var acc: V = @splat(0);
                            var e: usize = @intCast(i);
                            const stop: usize = @intCast(hi);
                            while (e + W <= stop) : (e += W) {
                                const av = ld(col + e);
                                st(y + e, mulAdd(V, @splat(t), av, ld(y + e)));
                                acc = mulAdd(V, ld(x + e), av, acc);
                            }
                            s = @reduce(.Add, acc);
                            i = @intCast(e);
                        }
                        while (i < hi) : (i += 1) {
                            const av = col[@intCast(i)];
                            const yp = at(F, y, i * incy);
                            yp.* = mulAdd(F, t, av, yp.*);
                            s = mulAdd(F, at(F, x, i * incx).*, av, s);
                        }
                        yj.* = mulAdd(F, alpha, s, yj.*);
                    }
                    return 0;
                }
            };
        }

        pub const symv_U = Symv(false).k;
        pub const symv_L = Symv(true).k;
    };
}
