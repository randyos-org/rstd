//! `@Vector` kernels for level-1 and level-2 slots, in the kernel ABI, `U`
//! native vectors a step; unit strides vectorize, others take a scalar loop,
//! and the special values are the portable kernels' (design.md). At `U = 4`
//! they mostly match or beat OpenBLAS's arm64 kernels on Apple's cores
//! (`tools/blasbench/results/results-m6.md`). `registry.Vector` places them.
//! Imports stay inside `linalg/blas/`, for the per-core objects.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("types.zig");
const util = @import("util.zig");
const slots = @import("slots.zig");
const Vtable = @import("Vtable.zig");
const fma = @import("fma.zig");
const extremum = @import("blas1/extremum.zig");
const bf16 = @import("ext/bf16/kernel.zig");

const BlasLong = types.BlasLong;
const Find = extremum.Find;
const Slot = util.Slot;

/// Lanes of one native vector of `F`: a NEON register on AArch64, where
/// `std.simd` counts SVE as 256 bits, which LLVM splits into NEON pairs.
pub fn lanes(comptime F: type) comptime_int {
    if (comptime sysinfo.target.cpu.arch.isAarch64()) return 16 / @sizeOf(F);
    return std.simd.suggestVectorLength(F) orelse 1;
}

/// `e` as `T`: itself where `T` is a float, else in every lane of `T`.
pub inline fn splatTo(comptime T: type, e: anytype) T {
    return if (@typeInfo(T) == .vector) @splat(e) else e;
}

inline fn at(comptime F: type, p: [*]F, i: BlasLong) *F {
    return util.at(F, p, i);
}

const abs = util.abs;

/// `a*b + c`, fused wherever an instruction does it (`fma.fusible`), whatever
/// the setting; elsewhere rounded twice, as `@mulAdd` there is a software
/// `fma` call per element.
inline fn mulAdd(comptime T: type, a: T, b: T, c: T) T {
    return if (comptime fma.fusible(T)) fma.fusedMadd(T, a, b, c) else a * b + c;
}

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
        const DV = @Vector(W, D);
        const DG = @Vector(4, D);

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
                var acc: [2]DV = @splat(@splat(0));
                var i: usize = 0;
                while (i + 2 * W <= len) : (i += 2 * W) {
                    inline for (0..2) |q| acc[q] = mulAdd(DV, ld(x + i + q * W), ld(y + i + q * W), acc[q]);
                }
                if (i + W <= len) {
                    acc[0] = mulAdd(DV, ld(x + i), ld(y + i), acc[0]);
                    i += W;
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
            var acc: [2]DV = @splat(@splat(0));
            var i: usize = 0;
            while (i + 2 * W <= len) : (i += 2 * W) {
                inline for (0..2) |q| {
                    const v: DV = ld(x + i + q * W);
                    acc[q] = mulAdd(DV, v, v, acc[q]);
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

        /// The sum of `|x[0..len]|`.
        pub fn absSum(len: usize, x: [*]F) F {
            var acc: [2]V = @splat(@splat(0));
            var i: usize = 0;
            while (i + 2 * W <= len) : (i += 2 * W) {
                inline for (0..2) |q| acc[q] += abs(ld(x + i + q * W));
            }
            var s = @reduce(.Add, acc[0] + acc[1]);
            while (i < len) : (i += 1) s += abs(x[i]);
            return s;
        }

        pub fn asum(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            if (n <= 0 or incx <= 0) return 0;
            if (incx == 1) return absSum(@intCast(n), x);
            var acc: G = @splat(0);
            var i: BlasLong = 0;
            while (i + 4 <= n) : (i += 4) acc += abs(gather(x, incx, i));
            var s = @reduce(.Add, acc);
            while (i < n) : (i += 1) s += abs(at(F, x, i * incx).*);
            return s;
        }

        /// The sum of `x[0..len]` in `D`, as `dot` sums: float `i` of each
        /// `2 * W` into lane `i % W` of sum `i / W % 2`; then the two lane by
        /// lane, the lanes in order, and the floats left in order.
        pub fn plainSum(len: usize, x: [*]F) D {
            var acc: [2]DV = @splat(@splat(0));
            var i: usize = 0;
            while (i + 2 * W <= len) : (i += 2 * W) {
                inline for (0..2) |q| acc[q] += @as(DV, ld(x + i + q * W));
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
                while (i < len) : (i += 1) std.mem.swap(F, &x[i], &y[i]);
                return 0;
            }
            var i: BlasLong = 0;
            while (i < n) : (i += 1) std.mem.swap(F, at(F, x, i * incx), at(F, y, i * incy));
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

/// Which operands `axpby` reads: a zero scalar's is not.
const Pick = enum {
    zero,
    ax,
    by,
    both,

    fn of(a: bool, b: bool) Pick {
        return if (a) (if (b) .both else .ax) else if (b) .by else .zero;
    }

    /// `alpha*x + beta*y` from `x` and `y` at `xp` and `yp` as `T`s, a
    /// vector or one float; `beta*y` first, `alpha*x` fused in.
    inline fn apply(comptime p: Pick, comptime T: type, a: anytype, b: @TypeOf(a), xp: anytype, yp: anytype) T {
        const F = @TypeOf(a);
        const s = struct {
            inline fn get(q: anytype) T {
                return if (T == F) q.* else q[0..@typeInfo(T).vector.len].*;
            }
        };
        return switch (p) {
            .zero => splatTo(T, @as(F, 0)),
            .ax => splatTo(T, a) * s.get(xp),
            .by => splatTo(T, b) * s.get(yp),
            .both => mulAdd(T, splatTo(T, a), s.get(xp), splatTo(T, b) * s.get(yp)),
        };
    }
};

/// Whether `@max` and `@min` of `F` are instructions here: elsewhere each
/// is a call to libm's `fmax` or `fmin`, lane by lane (ARMv7 VFP, SPARC,
/// pre-VSX POWER, MIPS before release 6, MSA too, the x87).
fn hasMinMax(comptime F: type) bool {
    const cpu = sysinfo.target.cpu;
    return switch (cpu.arch) {
        .x86, .x86_64 => !util.onX87(F),
        .arm, .armeb, .thumb, .thumbeb => cpu.has(.arm, .fp_armv8),
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => cpu.has(.powerpc, .vsx),
        .mips, .mipsel, .mips64, .mips64el => cpu.hasAny(.mips, &.{ .mips32r6, .mips64r6 }),
        .sparc, .sparc64 => false,
        else => true,
    };
}

/// `nrm2` from `s`, the plain sum of squares of `n` elements of `w` floats.
/// Summed in `f64` for `f32`, it is exact enough as it is. Else `@sqrt(s)`
/// where it neither overflowed nor could have lost more than `eps^2` of
/// itself to underflow, and Blue's algorithm (`util.blue`) where it could have.
fn root(comptime F: type, comptime w: comptime_int, s: anytype, n: BlasLong, x: [*]F, incx: BlasLong) F {
    if (comptime @TypeOf(s) != F) return @floatCast(@sqrt(s));
    const lo = std.math.floatMin(F) / std.math.floatEps(F) * @as(F, @floatFromInt(w * n));
    if (s >= lo and s <= std.math.floatMax(F)) return @sqrt(s);
    return util.blue(F, w, n, x, incx);
}

/// The complex kernels over interleaved `F` pairs.
pub fn Complex(comptime F: type, comptime U: comptime_int) type {
    return struct {
        const R = Real(F, U);
        /// Floats a vector; whole elements.
        pub const W = if (R.W % 2 == 0) R.W else 2 * R.W;
        const V = @Vector(W, F);
        /// One element.
        const V2 = @Vector(2, F);
        const Result = types.CComplex(F);

        inline fn ld(p: [*]const F) V {
            return p[0..W].*;
        }
        inline fn st(p: [*]F, v: V) void {
            p[0..W].* = v;
        }
        inline fn ld2(p: [*]F, i: BlasLong) V2 {
            return @as(*[2]F, @ptrCast(at(F, p, 2 * i))).*;
        }
        inline fn st2(p: [*]F, i: BlasLong, v: V2) void {
            @as(*[2]F, @ptrCast(at(F, p, 2 * i))).* = v;
        }

        /// `(re, im, re, im, ...)` as a `T`.
        inline fn pair(comptime T: type, re: F, im: F) T {
            var r: [@typeInfo(T).vector.len]F = undefined;
            inline for (0..r.len) |k| r[k] = if (k % 2 == 0) re else im;
            return r;
        }

        /// Each element's parts exchanged.
        inline fn swapped(comptime T: type, v: T) T {
            const len = @typeInfo(T).vector.len;
            const mask = comptime blk: {
                var m: [len]i32 = undefined;
                for (0..len) |k| m[k] = @intCast(k ^ 1);
                break :blk m;
            };
            return @shuffle(F, v, undefined, mask);
        }

        /// Multiplying by `(tr, ti)`, the operand conjugated for `conj`: lane
        /// factors on it and on it with its parts exchanged.
        fn By(comptime T: type) type {
            return struct {
                r: T,
                i: T,

                inline fn init(tr: F, ti: F, comptime conj: bool) @This() {
                    return .{ .r = pair(T, tr, if (conj) -tr else tr), .i = pair(T, if (conj) ti else -ti, ti) };
                }
                /// `y + t*v`.
                inline fn madd(b: @This(), v: T, y: T) T {
                    return mulAdd(T, b.r, v, mulAdd(T, b.i, swapped(T, v), y));
                }
                /// `t*v`.
                inline fn mul(b: @This(), v: T) T {
                    return mulAdd(T, b.i, swapped(T, v), b.r * v);
                }
            };
        }

        /// The sums of `a*x` lane by lane and of `a*x` with `x`'s parts
        /// exchanged, as `Dots.finish` makes `sum(op(a)*op(x))` of them.
        fn Dots(comptime ca: bool, comptime cx: bool) type {
            return struct {
                rr: V = @splat(0),
                ri: V = @splat(0),
                rr2: V2 = @splat(0),
                ri2: V2 = @splat(0),

                inline fn add(d: *@This(), av: V, xv: V) void {
                    d.rr = mulAdd(V, av, xv, d.rr);
                    d.ri = mulAdd(V, av, swapped(V, xv), d.ri);
                }
                inline fn add2(d: *@This(), av: V2, xv: V2) void {
                    d.rr2 = mulAdd(V2, av, xv, d.rr2);
                    d.ri2 = mulAdd(V2, av, swapped(V2, xv), d.ri2);
                }
                /// `re = ar*xr - sa*sx*ai*xi`, `im = sx*ar*xi + sa*ai*xr`,
                /// `sa` and `sx` the signs conjugation gives the imaginary parts.
                fn finish(d: @This()) V2 {
                    const sa: F = if (ca) -1 else 1;
                    const sx: F = if (cx) -1 else 1;
                    const re = @reduce(.Add, d.rr * pair(V, 1, -sa * sx)) + @reduce(.Add, d.rr2 * pair(V2, 1, -sa * sx));
                    const im = @reduce(.Add, d.ri * pair(V, sx, sa)) + @reduce(.Add, d.ri2 * pair(V2, sx, sa));
                    return .{ re, im };
                }
            };
        }

        fn Axpy(comptime conj: bool) type {
            return struct {
                fn k(n: BlasLong, _: BlasLong, _: BlasLong, ar: F, ai: F, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F, _: BlasLong) callconv(.c) i32 {
                    if (n <= 0) return 0;
                    var i: BlasLong = 0;
                    if (incx == 1 and incy == 1) {
                        const by = By(V).init(ar, ai, conj);
                        const len: usize = @intCast(2 * n);
                        var f: usize = 0;
                        while (f + W <= len) : (f += W) st(y + f, by.madd(ld(x + f), ld(y + f)));
                        i = @intCast(f / 2);
                    }
                    const by = By(V2).init(ar, ai, conj);
                    while (i < n) : (i += 1) st2(y, i * incy, by.madd(ld2(x, i * incx), ld2(y, i * incy)));
                    return 0;
                }
            };
        }

        /// `y += alpha*x`.
        pub const axpy = Axpy(false).k;
        /// `y += alpha*conj(x)`.
        pub const axpyc = Axpy(true).k;

        /// `sum(op(x)*y)`, `x` conjugated for `conj`.
        fn dot(comptime conj: bool, n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) Result {
            var d: [2]Dots(conj, false) = .{ .{}, .{} };
            var i: BlasLong = 0;
            if (n > 0 and incx == 1 and incy == 1) {
                const len: usize = @intCast(2 * n);
                var f: usize = 0;
                while (f + 2 * W <= len) : (f += 2 * W) {
                    inline for (0..2) |q| d[q].add(ld(x + f + q * W), ld(y + f + q * W));
                }
                i = @intCast(f / 2);
            }
            while (i + 2 <= n) : (i += 2) {
                inline for (0..2) |q| d[q].add2(ld2(x, (i + q) * incx), ld2(y, (i + q) * incy));
            }
            if (i < n) d[0].add2(ld2(x, i * incx), ld2(y, i * incy));
            const s = d[0].finish() + d[1].finish();
            return .{ .re = s[0], .im = s[1] };
        }

        pub fn dotu(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) Result {
            return dot(false, n, x, incx, y, incy);
        }
        pub fn dotc(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) Result {
            return dot(true, n, x, incx, y, incy);
        }

        /// A zero `alpha` with `flag` 0 stores zeros; any other multiplies.
        pub fn scal(n: BlasLong, _: BlasLong, _: BlasLong, ar: F, ai: F, x: [*]F, incx: BlasLong, _: [*]F, _: BlasLong, _: [*]F, flag: BlasLong) callconv(.c) i32 {
            if (n <= 0 or incx <= 0) return 0;
            if (flag == 0 and ar == 0 and ai == 0) {
                var i: BlasLong = 0;
                while (i < n) : (i += 1) st2(x, i * incx, @splat(0));
                return 0;
            }
            var i: BlasLong = 0;
            if (incx == 1) {
                const by = By(V).init(ar, ai, false);
                const len: usize = @intCast(2 * n);
                var f: usize = 0;
                while (f + W <= len) : (f += W) st(x + f, by.mul(ld(x + f)));
                i = @intCast(f / 2);
            }
            const by = By(V2).init(ar, ai, false);
            while (i < n) : (i += 1) st2(x, i * incx, by.mul(ld2(x, i * incx)));
            return 0;
        }

        pub fn nrm2(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            if (n <= 0 or incx == 0) return 0;
            if (incx == 1) return root(F, 2, R.squares(@intCast(2 * n), x), n, x, incx);
            const D2 = @Vector(2, R.D);
            var acc: [2]D2 = @splat(@splat(0));
            var i: BlasLong = 0;
            while (i + 2 <= n) : (i += 2) {
                inline for (0..2) |q| {
                    const v: D2 = ld2(x, (i + q) * incx);
                    acc[q] = mulAdd(D2, v, v, acc[q]);
                }
            }
            if (i < n) {
                const v: D2 = ld2(x, i * incx);
                acc[0] = mulAdd(D2, v, v, acc[0]);
            }
            return root(F, 2, @reduce(.Add, acc[0] + acc[1]), n, x, incx);
        }

        /// The sum of `|re| + |im|`.
        pub fn asum(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            if (n <= 0 or incx <= 0) return 0;
            if (incx == 1) return R.absSum(@intCast(2 * n), x);
            var acc: [2]V2 = @splat(@splat(0));
            var i: BlasLong = 0;
            while (i + 2 <= n) : (i += 2) {
                inline for (0..2) |q| acc[q] += abs(ld2(x, (i + q) * incx));
            }
            if (i < n) acc[0] += abs(ld2(x, i * incx));
            return @reduce(.Add, acc[0] + acc[1]);
        }

        /// The sum of `re + im`: `Real.plainSum`'s over the `2n` floats;
        /// strided, element `i` into sum `i % 2`, by parts, then those lane by
        /// lane and the parts in order.
        pub fn sum(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            if (n <= 0 or incx <= 0) return 0;
            if (incx == 1) return @floatCast(R.plainSum(@intCast(2 * n), x));
            const D2 = @Vector(2, R.D);
            var acc: [2]D2 = @splat(@splat(0));
            var i: BlasLong = 0;
            while (i + 2 <= n) : (i += 2) {
                inline for (0..2) |q| acc[q] += @as(D2, ld2(x, (i + q) * incx));
            }
            if (i < n) acc[0] += @as(D2, ld2(x, i * incx));
            return @floatCast(@reduce(.Add, acc[0] + acc[1]));
        }

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
                    inline fn key(p: [*]const F, s: usize, i: usize) F {
                        return util.rounded(F, abs(p[i * s]) + abs(p[i * s + 1]));
                    }
                    inline fn block(p: [*]const F, s: usize, i: usize, chk: *V) F {
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

                fn search(n: BlasLong, x: [*]F, incx: BlasLong) ?struct { usize, F } {
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

        pub fn copy(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) return R.copy(2 * n, x, 1, y, 1);
            var i: BlasLong = 0;
            while (i < n) : (i += 1) st2(y, i * incy, ld2(x, i * incx));
            return 0;
        }

        pub fn swap(n: BlasLong, _: BlasLong, _: BlasLong, _: F, _: F, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F, _: BlasLong) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) return R.swap(2 * n, 0, 0, 0, x, 1, y, 1, x, 0);
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const t = ld2(x, i * incx);
                st2(x, i * incx, ld2(y, i * incy));
                st2(y, i * incy, t);
            }
            return 0;
        }

        /// `csrot`, `zdrot`: both parts rotated by the real `c` and `s`.
        pub fn rot(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, c: F, s: F) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) return R.rot(2 * n, x, 1, y, 1, c, s);
            const cv: V2 = @splat(c);
            const sv: V2 = @splat(s);
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const xv = ld2(x, i * incx);
                const yv = ld2(y, i * incy);
                st2(x, i * incx, mulAdd(V2, cv, xv, sv * yv));
                st2(y, i * incy, mulAdd(V2, cv, yv, -(sv * xv)));
            }
            return 0;
        }

        /// `y += alpha*op(A)*x`: transposed for `trans`, `A` conjugated for
        /// `ca`, `x` for `cx`, as `gemv/`'s portable `Zgemv`. Without a
        /// transpose, four columns a pass over a unit-stride `y`; with one,
        /// four columns' sums a pass over a unit-stride `x`.
        fn Gemv(comptime trans: bool, comptime ca: bool, comptime cx: bool) type {
            return struct {
                fn k(m: BlasLong, n: BlasLong, _: BlasLong, ar: F, ai: F, a: [*]F, lda: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F) callconv(.c) i32 {
                    if (n <= 0 or (!trans and m <= 0)) return 0;
                    const rows: BlasLong = @max(m, 0);
                    const ld_a: usize = 2 * @as(usize, @intCast(lda));
                    const alpha = By(V2).init(ar, ai, false);
                    // the vectors' stride: `y`'s, or with a transpose `x`'s
                    const unit = (if (trans) incx else incy) == 1;
                    var j: BlasLong = 0;
                    while (j < n) {
                        const g: BlasLong = if (unit and j + 4 <= n) 4 else 1;
                        var c: [4][*]F = undefined;
                        for (0..@intCast(g)) |q| c[q] = a + (@as(usize, @intCast(j)) + q) * ld_a;
                        if (trans) {
                            var d: [4]Dots(ca, cx) = @splat(.{});
                            var i: BlasLong = 0;
                            if (g == 4) {
                                const len: usize = @intCast(2 * rows);
                                var f: usize = 0;
                                while (f + W <= len) : (f += W) {
                                    const xv = ld(x + f);
                                    inline for (0..4) |q| d[q].add(ld(c[q] + f), xv);
                                }
                                i = @intCast(f / 2);
                            }
                            while (i < rows) : (i += 1) {
                                const xv = ld2(x, i * incx);
                                for (0..@intCast(g)) |q| d[q].add2(ld2(c[q], i), xv);
                            }
                            for (0..@intCast(g)) |q| {
                                const at_y = (j + @as(BlasLong, @intCast(q))) * incy;
                                st2(y, at_y, alpha.madd(d[q].finish(), ld2(y, at_y)));
                            }
                        } else {
                            var t: [4]By(V) = undefined;
                            var t2: [4]By(V2) = undefined;
                            for (0..@intCast(g)) |q| {
                                const tv = By(V2).init(ar, ai, cx).mul(ld2(x, (j + @as(BlasLong, @intCast(q))) * incx));
                                t[q] = .init(tv[0], tv[1], ca);
                                t2[q] = .init(tv[0], tv[1], ca);
                            }
                            var i: BlasLong = 0;
                            if (g == 4) {
                                const len: usize = @intCast(2 * rows);
                                var f: usize = 0;
                                while (f + W <= len) : (f += W) {
                                    var yv = ld(y + f);
                                    inline for (0..4) |q| yv = t[q].madd(ld(c[q] + f), yv);
                                    st(y + f, yv);
                                }
                                i = @intCast(f / 2);
                            }
                            while (i < rows) : (i += 1) {
                                var yv = ld2(y, i * incy);
                                for (0..@intCast(g)) |q| yv = t2[q].madd(ld2(c[q], i), yv);
                                st2(y, i * incy, yv);
                            }
                        }
                        j += g;
                    }
                    return 0;
                }
            };
        }

        pub const gemv_n = Gemv(false, false, false).k;
        pub const gemv_t = Gemv(true, false, false).k;
        pub const gemv_r = Gemv(false, true, false).k;
        pub const gemv_c = Gemv(true, true, false).k;
        pub const gemv_o = Gemv(false, false, true).k;
        pub const gemv_u = Gemv(true, false, true).k;
        pub const gemv_s = Gemv(false, true, true).k;
        pub const gemv_d = Gemv(true, true, true).k;

        /// `y := alpha*x + beta*y`, `kernel/riscv64/zaxpby_rvv.c`: `beta*y`,
        /// then `alpha*x` fused in. A zero `alpha` does not read `x`, a zero
        /// `beta` not `y`.
        pub fn axpby(n: BlasLong, ar: F, ai: F, x: [*]F, incx: BlasLong, br: F, bi: F, y: [*]F, incy: BlasLong) callconv(.c) i32 {
            if (n <= 0) return 0;
            switch (Pick.of(ar != 0 or ai != 0, br != 0 or bi != 0)) {
                inline else => |p| axpbyBy(p, n, ar, ai, x, incx, br, bi, y, incy),
            }
            return 0;
        }

        inline fn axpby1(comptime p: Pick, comptime T: type, a: By(T), b: By(T), xv: anytype, yv: anytype) T {
            return switch (p) {
                .zero => @splat(0),
                .ax => a.mul(xv.*),
                .by => b.mul(yv.*),
                .both => a.madd(xv.*, b.mul(yv.*)),
            };
        }

        fn axpbyBy(comptime p: Pick, n: BlasLong, ar: F, ai: F, x: [*]F, incx: BlasLong, br: F, bi: F, y: [*]F, incy: BlasLong) void {
            var i: BlasLong = 0;
            if (incx == 1 and incy == 1) {
                const a = By(V).init(ar, ai, false);
                const b = By(V).init(br, bi, false);
                const len: usize = @intCast(2 * n);
                var f: usize = 0;
                while (f + W <= len) : (f += W) st(y + f, axpby1(p, V, a, b, x[f..][0..W], y[f..][0..W]));
                i = @intCast(f / 2);
            }
            const a = By(V2).init(ar, ai, false);
            const b = By(V2).init(br, bi, false);
            while (i < n) : (i += 1) {
                const yp: *[2]F = @ptrCast(at(F, y, 2 * i * incy));
                yp.* = axpby1(p, V2, a, b, @as(*[2]F, @ptrCast(at(F, x, 2 * i * incx))), yp);
            }
        }

        /// `y += alpha*S*x` over the columns `offset` selects, `S` from one
        /// stored triangle by `kind`; per column, the diagonal first (its
        /// imaginary part unread for Hermitian). `kernel/riscv64/zsymv_L_rvv.c`,
        /// `zsymv_U_rvv.c`, `zhemv_LM_rvv.c`, `zhemv_UV_rvv.c`, `zhemv_*_vector.c`.
        fn Symv(comptime lower: bool, comptime kind: SymKind) type {
            // `y` takes the column conjugated for `HEMVREV` (`V`, `M`), the
            // dot for the plain Hermitian (`U`, `L`)
            const cy = kind == .hermitian_rev;
            const cd = kind == .hermitian;
            return struct {
                fn k(m: BlasLong, offset: BlasLong, ar: F, ai: F, a: [*]F, lda: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F) callconv(.c) i32 {
                    if (m <= 0 or offset <= 0) return 0;
                    const alpha = By(V2).init(ar, ai, false);
                    const ld_a: usize = 2 * @as(usize, @intCast(lda));
                    const unit = incx == 1 and incy == 1;
                    var j: BlasLong = if (lower) 0 else m - offset;
                    const end: BlasLong = if (lower) offset else m;
                    while (j < end) : (j += 1) {
                        const col = a + @as(usize, @intCast(j)) * ld_a;
                        const t = alpha.mul(ld2(x, j * incx));
                        const dy: V2 = if (kind == .symmetric)
                            By(V2).init(t[0], t[1], false).mul(ld2(col, j))
                        else
                            t * @as(V2, @splat(col[2 * @as(usize, @intCast(j))]));
                        st2(y, j * incy, ld2(y, j * incy) + dy);
                        var i: BlasLong = if (lower) j + 1 else 0;
                        const hi: BlasLong = if (lower) m else j;
                        var d: Dots(cd, false) = .{};
                        if (unit) {
                            const by = By(V).init(t[0], t[1], cy);
                            var f: usize = 2 * @as(usize, @intCast(i));
                            const stop: usize = 2 * @as(usize, @intCast(hi));
                            while (f + W <= stop) : (f += W) {
                                const av = ld(col + f);
                                st(y + f, by.madd(av, ld(y + f)));
                                d.add(av, ld(x + f));
                            }
                            i = @intCast(f / 2);
                        }
                        const by = By(V2).init(t[0], t[1], cy);
                        while (i < hi) : (i += 1) {
                            const av = ld2(col, i);
                            st2(y, i * incy, by.madd(av, ld2(y, i * incy)));
                            d.add2(av, ld2(x, i * incx));
                        }
                        st2(y, j * incy, alpha.madd(d.finish(), ld2(y, j * incy)));
                    }
                    return 0;
                }
            };
        }

        pub const symv_U = Symv(false, .symmetric).k;
        pub const symv_L = Symv(true, .symmetric).k;
        pub const hemv_U = Symv(false, .hermitian).k;
        pub const hemv_L = Symv(true, .hermitian).k;
        pub const hemv_V = Symv(false, .hermitian_rev).k;
        pub const hemv_M = Symv(true, .hermitian_rev).k;
    };
}

/// How `S` comes from its stored triangle (`blas2/symmetric.zig`'s `Kind`):
/// symmetric, Hermitian, or the Hermitian one's transpose (`HEMVREV`).
const SymKind = enum { symmetric, hermitian, hermitian_rev };

/// The half-precision gemvs (`sbgemv`, `bgemv`, `shgemv`) in
/// `ext/bf16/kernel.zig`'s bits: `f32` sums in index order, unfused, one
/// rounding on the store. Only the shape of upstream's kernels is kept
/// (`kernel/riscv64/sbgemv_?_vector.c`, `kernel/x86_64/sbgemv_n_microk_cooperlake.c`,
/// `kernel/arm64/sbgemv_n_neon.c`, `sbgemv_t_bfdot.c`, `bgemv_n_sve_v3x4.c`),
/// whose fused or paired products round otherwise.
pub fn Half(comptime E: type, comptime Y: type) type {
    return struct {
        const N = lanes(f32);
        const V = @Vector(N, f32);

        inline fn wide(e: E) f32 {
            return bf16.wide(E, e);
        }

        /// `N` elements from `p`, widened as `bf16.widen` does: a bfloat16
        /// denormal flushes to a signed zero and a NaN comes out quiet.
        inline fn widen(p: [*]const E) V {
            const w = bf16.widen(E, N, p);
            return if (N == 1) @splat(w) else w;
        }

        /// `alpha*acc + beta*y`, rounded once each; `beta == 0` does not read `y`.
        inline fn result(comptime T: type, alpha: f32, beta: f32, acc: T, yv: anytype) T {
            if (beta == 0) return splatTo(T, alpha) * acc;
            return splatTo(T, alpha) * acc + splatTo(T, beta) * @as(T, yv.*);
        }

        /// `y[at] = alpha*acc + beta*y[at]` in `f32`, rounded to `Y`.
        inline fn put(yp: *Y, alpha: f32, beta: f32, acc: f32) void {
            const wy = if (beta == 0) 0 else bf16.wide(Y, yp.*);
            yp.* = bf16.narrow(Y, result(f32, alpha, beta, acc, &wy));
        }

        /// Rows `r` to `r + B` from column 0 on, `B / N` vectors.
        fn rowBlock(comptime B: usize, r: usize, cols: usize, alpha: f32, a: [*]E, ld_a: usize, x: [*]E, incx: BlasLong, beta: f32, y: [*]Y, incy: BlasLong) void {
            var acc: [B / N]V = @splat(@splat(0));
            for (0..cols) |q| {
                const xv: V = @splat(wide(util.at(E, x, @as(BlasLong, @intCast(q)) * incx).*));
                inline for (&acc, 0..) |*s, v| s.* = s.* + widen(a + r + v * N + q * ld_a) * xv;
            }
            inline for (acc, 0..) |s, v| {
                const r0 = r + v * N;
                if (Y == f32 and incy == 1) {
                    const yp: *[N]f32 = (y + r0)[0..N];
                    yp.* = result(V, alpha, beta, s, yp);
                } else {
                    const sa: [N]f32 = s;
                    inline for (sa, 0..) |e, l| {
                        put(util.at(Y, y, @as(BlasLong, @intCast(r0 + l)) * incy), alpha, beta, e);
                    }
                }
            }
        }

        /// `y = alpha*A*x + beta*y`, a lane per row of `A`.
        pub fn gemv_n(m: BlasLong, n: BlasLong, alpha_y: Y, a: [*]E, lda: BlasLong, x: [*]E, incx: BlasLong, beta_y: Y, y: [*]Y, incy: BlasLong) callconv(.c) i32 {
            const alpha = bf16.wide(Y, alpha_y);
            const beta = bf16.wide(Y, beta_y);
            const rows: usize = @intCast(@max(m, 0));
            const cols: usize = @intCast(@max(n, 0));
            const ld_a: usize = @intCast(lda);
            var r: usize = 0;
            while (r + 4 * N <= rows) : (r += 4 * N) rowBlock(4 * N, r, cols, alpha, a, ld_a, x, incx, beta, y, incy);
            while (r + N <= rows) : (r += N) rowBlock(N, r, cols, alpha, a, ld_a, x, incx, beta, y, incy);
            while (r < rows) : (r += 1) {
                var acc: f32 = 0;
                for (0..cols) |q| acc = acc + wide(a[r + q * ld_a]) * wide(util.at(E, x, @as(BlasLong, @intCast(q)) * incx).*);
                put(util.at(Y, y, @as(BlasLong, @intCast(r)) * incy), alpha, beta, acc);
            }
            return 0;
        }

        /// `y = alpha*A'*x + beta*y`, a lane per column of `A`: `N` columns'
        /// `N` rows at a time, transposed, so each column sums in row order.
        pub fn gemv_t(m: BlasLong, n: BlasLong, alpha_y: Y, a: [*]E, lda: BlasLong, x: [*]E, incx: BlasLong, beta_y: Y, y: [*]Y, incy: BlasLong) callconv(.c) i32 {
            const alpha = bf16.wide(Y, alpha_y);
            const beta = bf16.wide(Y, beta_y);
            const rows: usize = @intCast(@max(m, 0));
            const cols: usize = @intCast(@max(n, 0));
            const ld_a: usize = @intCast(lda);
            var c: usize = 0;
            while (c + N <= cols) : (c += N) {
                var acc: V = @splat(0);
                var q: usize = 0;
                while (q + N <= rows) : (q += N) {
                    var t: [N]V = undefined;
                    inline for (&t, 0..) |*v, l| v.* = widen(a + q + (c + l) * ld_a);
                    t = transpose(t);
                    inline for (t, 0..) |v, l| acc = acc + v * @as(V, @splat(wide(util.at(E, x, @as(BlasLong, @intCast(q + l)) * incx).*)));
                }
                while (q < rows) : (q += 1) {
                    var row: [N]f32 = undefined;
                    inline for (&row, 0..) |*e, l| e.* = wide(a[q + (c + l) * ld_a]);
                    acc = acc + @as(V, row) * @as(V, @splat(wide(util.at(E, x, @as(BlasLong, @intCast(q)) * incx).*)));
                }
                const sa: [N]f32 = acc;
                inline for (sa, 0..) |e, l| {
                    put(util.at(Y, y, @as(BlasLong, @intCast(c + l)) * incy), alpha, beta, e);
                }
            }
            while (c < cols) : (c += 1) {
                var acc: f32 = 0;
                for (0..rows) |q| acc = acc + wide(a[q + c * ld_a]) * wide(util.at(E, x, @as(BlasLong, @intCast(q)) * incx).*);
                put(util.at(Y, y, @as(BlasLong, @intCast(c)) * incy), alpha, beta, acc);
            }
            return 0;
        }

        /// `t` transposed: halves of rows trade their off-diagonal halves,
        /// then quarters, down to single lanes.
        inline fn transpose(t: [N]V) [N]V {
            var v = t;
            comptime var h = N / 2;
            inline while (h >= 1) : (h /= 2) {
                const lo = comptime mask(h, false);
                const hi = comptime mask(h, true);
                inline for (0..N) |r| {
                    if (r & h != 0) continue;
                    const p = v[r];
                    const q = v[r + h];
                    v[r] = @shuffle(f32, p, q, lo);
                    v[r + h] = @shuffle(f32, p, q, hi);
                }
            }
            return v;
        }

        /// The upper (`second`) or lower row of a pair after trading the
        /// lanes with bit `h` set in the first for those clear in the second.
        fn mask(comptime h: usize, comptime second: bool) [N]i32 {
            var m: [N]i32 = undefined;
            for (0..N) |l| {
                const own = l & h == 0;
                m[l] = if (!second)
                    (if (own) @intCast(l) else ~@as(i32, @intCast(l - h)))
                else
                    (if (own) @intCast(l + h) else ~@as(i32, @intCast(l)));
            }
            return m;
        }
    };
}

// ---- the slots they fill -----------------------------------------------------

/// Whether `slot`'s converted kernels give way to the portable or these on
/// every core, for one contract (`readme.md`, "One contract"): `nrm2`, the
/// extremes, `sdot` and `dsdot`, and the `peels` slots. `tools/blas_asm`
/// leaves their assembly out (`roles.zig`'s `skip`, kept in step).
pub fn routed(comptime slot: Slot) bool {
    return switch (slot) {
        // the extremes, `nrm2`, `sdot` and `dsdot`
        .samax_k,
        .damax_k,
        .camax_k,
        .zamax_k,
        .samin_k,
        .damin_k,
        .camin_k,
        .zamin_k,
        .smax_k,
        .dmax_k,
        .smin_k,
        .dmin_k,
        .isamax_k,
        .idamax_k,
        .icamax_k,
        .izamax_k,
        .isamin_k,
        .idamin_k,
        .icamin_k,
        .izamin_k,
        .ismax_k,
        .idmax_k,
        .ismin_k,
        .idmin_k,
        .snrm2_k,
        .dnrm2_k,
        .cnrm2_k,
        .znrm2_k,
        .sdot_k,
        .dsdot_k,
        => true,
        else => peels(slot),
    };
}

/// The x86_64 reductions whose converted kernels peel elements until an
/// operand is aligned, so their bits follow its address: `asum_sse.S`,
/// `asum_sse2.S`, `zasum_sse*.S`, `zsum_sse*.S`, `dot_sse2.S`,
/// `ddot_bulldozer.S`, `zdot_sse*.S`, `dgemv_n.S` and their `atom` and
/// `bulldozer` kin. Their C kernels on other cores stay.
pub fn peels(comptime slot: Slot) bool {
    if (sysinfo.target.cpu.arch != .x86_64) return false;
    return std.mem.indexOfScalar(Slot, &peeling, slot) != null;
}

/// The slots `peels` names.
pub const peeling = [_]Slot{ .sasum_k, .dasum_k, .casum_k, .zasum_k, .csum_k, .zsum_k, .ddot_k, .cdotu_k, .cdotc_k, .zdotu_k, .zdotc_k, .dgemv_n };

/// The functions of `Real`, each filling the slot of its op (`slotOf`).
const real_fns = [_][:0]const u8{ "axpy", "dot", "dsdot", "scal", "nrm2", "asum", "sum", "copy", "swap", "rot", "gemv_n", "gemv_t", "amax", "iamax", "amin", "iamin", "max", "imax", "min", "imin" };

/// The functions of `Complex`, likewise.
const complex_fns = [_][:0]const u8{ "axpy", "axpyc", "dotu", "dotc", "scal", "nrm2", "asum", "sum", "copy", "swap", "amax", "iamax", "amin", "iamin", "rot", "gemv_n", "gemv_t", "gemv_r", "gemv_c", "gemv_o", "gemv_u", "gemv_s", "gemv_d" };

/// The slot function `name` over `T` fills: its `slots.Op` at `T`'s letter, or
/// null for `dsdot` past `f32`.
fn slotOf(comptime T: type, comptime name: []const u8) ?Slot {
    if (std.mem.eql(u8, name, "dsdot")) return if (T == f32) .dsdot_k else null;
    const op: slots.Op = if (@hasField(slots.Op, name ++ "_k")) @field(slots.Op, name ++ "_k") else @field(slots.Op, name);
    return slots.slotTag(.from(slots.letter(T).?), op);
}

/// A kernel: the slot it fills, and the function of `Real` or `Complex` over
/// `F` that is it.
const Kernel = struct { slot: Slot, F: type, complex: bool, f: [:0]const u8 };

/// Every function here as a kernel, whatever its type's lanes do with
/// denormals.
pub const every: []const Kernel = blk: {
    var r: []const Kernel = &.{};
    for ([_]type{ f32, f64 }) |F| for (real_fns) |name| {
        if (slotOf(F, name)) |s| r = r ++ &[_]Kernel{.{ .slot = s, .F = F, .complex = false, .f = name }};
    };
    for ([_]type{ f32, f64 }) |F| for (complex_fns) |name| {
        r = r ++ &[_]Kernel{.{ .slot = slotOf(std.math.Complex(F), name).?, .F = F, .complex = true, .f = name }};
    };
    break :blk r;
};

/// Every kernel here: none in a type whose lanes flush denormals
/// (`util.vectorIeee`), whose slots keep the portable kernels.
pub const kernels: []const Kernel = blk: {
    var r: []const Kernel = &.{};
    for (every) |k| {
        if (util.vectorIeee(k.F)) r = r ++ &[_]Kernel{k};
    }
    break :blk r;
};

/// The kernel of `k`, `U` vectors a step.
pub fn kernel(comptime k: Kernel, comptime U: comptime_int) util.SlotFn(k.slot) {
    return &@field(if (k.complex) Complex(k.F, U) else Real(k.F, U), k.f);
}

/// The symbol of `slot`'s kernel here for `core`, apart from a converted one's.
pub fn symbol(comptime slot: Slot, comptime core: anytype) []const u8 {
    return util.kernelSymbol(slot, core) ++ "_vec";
}

/// Points every slot these kernels fill at this build's own, four vectors a
/// step, `ported`'s too, for the tests.
pub fn installAll(vt: *Vtable) void {
    inline for (kernels) |k| @field(vt, @tagName(k.slot)) = kernel(k, 4);
    inline for (ported.slots) |s| {
        if (comptime ported.builds(s)) @field(vt, @tagName(s)) = ported.entry(s);
    }
}

/// The kernels here no `Registry` routes, outside `real_slots` and
/// `complex_slots`: an architecture's hand ports name them (`entry`, as a
/// port's module) where upstream has intrinsics for them.
pub const ported = struct {
    pub const slots = [_]Slot{
        .saxpby_k, .daxpby_k, .caxpby_k, .zaxpby_k,
        .srotm_k,  .drotm_k,  .ssymv_U,  .ssymv_L,
        .dsymv_U,  .dsymv_L,  .csymv_U,  .csymv_L,
        .zsymv_U,  .zsymv_L,  .chemv_U,  .chemv_L,
        .chemv_V,  .chemv_M,  .zhemv_U,  .zhemv_L,
        .zhemv_V,  .zhemv_M,  .sbgemv_n, .sbgemv_t,
        .shgemv_n, .shgemv_t, .bgemv_n,  .bgemv_t,
    };

    /// `slot`'s kernel, four vectors a step.
    pub fn entry(comptime slot: Slot) util.SlotFn(slot) {
        const R = Real(f32, 4);
        const D = Real(f64, 4);
        const C = Complex(f32, 4);
        const Z = Complex(f64, 4);
        return switch (slot) {
            .saxpby_k => &R.axpby,
            .daxpby_k => &D.axpby,
            .caxpby_k => &C.axpby,
            .zaxpby_k => &Z.axpby,
            .srotm_k => &R.rotm,
            .drotm_k => &D.rotm,
            .ssymv_U => &R.symv_U,
            .ssymv_L => &R.symv_L,
            .dsymv_U => &D.symv_U,
            .dsymv_L => &D.symv_L,
            .csymv_U => &C.symv_U,
            .csymv_L => &C.symv_L,
            .zsymv_U => &Z.symv_U,
            .zsymv_L => &Z.symv_L,
            .chemv_U => &C.hemv_U,
            .chemv_L => &C.hemv_L,
            .chemv_V => &C.hemv_V,
            .chemv_M => &C.hemv_M,
            .zhemv_U => &Z.hemv_U,
            .zhemv_L => &Z.hemv_L,
            .zhemv_V => &Z.hemv_V,
            .zhemv_M => &Z.hemv_M,
            .sbgemv_n => &Half(types.BF16, f32).gemv_n,
            .sbgemv_t => &Half(types.BF16, f32).gemv_t,
            .shgemv_n => &Half(types.HF16, f32).gemv_n,
            .shgemv_t => &Half(types.HF16, f32).gemv_t,
            .bgemv_n => &Half(types.BF16, types.BF16).gemv_n,
            .bgemv_t => &Half(types.BF16, types.BF16).gemv_t,
            else => @compileError(@tagName(slot) ++ " has no kernel in vector.ported"),
        };
    }

    /// Whether this build makes `slot`'s kernel for `installAll`: none over
    /// `f32` whose lanes flush (`util.vectorIeee`), and the half-precision
    /// ones only where a port uses them or the native tests run, since
    /// `f16` vectors stop LLVM on some targets.
    pub fn builds(comptime slot: Slot) bool {
        const name = @tagName(slot);
        const half = for ([_][]const u8{ "sbgemv", "shgemv", "bgemv" }) |f| {
            if (std.mem.startsWith(u8, name, f)) break true;
        } else false;
        if (half) return switch (sysinfo.target.cpu.arch) {
            .x86_64, .aarch64, .riscv64 => true,
            else => false,
        };
        const single = name[0] == 's' or name[0] == 'c';
        return !single or util.vectorIeee(f32);
    }
};
