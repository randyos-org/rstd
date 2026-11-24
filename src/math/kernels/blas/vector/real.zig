//! The real `@Vector` kernels (`vector.zig`).

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../types.zig");
const util = @import("../util.zig");
const extremum = @import("../blas1/extremum.zig");
const shared = @import("util.zig");
const real_reductions = @import("real/reductions.zig");
const real_extremes = @import("real/extremes.zig");
const real_level2 = @import("real/level2.zig");

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
        const reductions = real_reductions.Reductions(F, U);
        pub const absSum = reductions.absSum;
        pub const asum = reductions.asum;
        pub const asum_wide = reductions.asum_wide;
        pub const dot = reductions.dot;
        pub const dsdot = reductions.dsdot;
        pub const nrm2 = reductions.nrm2;
        pub const plainSum = reductions.plainSum;
        pub const squares = reductions.squares;
        pub const sum = reductions.sum;
        pub const sum_wide = reductions.sum_wide;
        pub const unitSum = reductions.unitSum;

        const extremes = real_extremes.Extremes(F, U);
        pub const amax = extremes.amax;
        pub const amin = extremes.amin;
        pub const iamax = extremes.iamax;
        pub const iamin = extremes.iamin;
        pub const imax = extremes.imax;
        pub const imin = extremes.imin;
        pub const max = extremes.max;
        pub const min = extremes.min;
        pub const pick = extremes.pick;
        pub const search = extremes.search;

        const level2 = real_level2.Level2(F, U);
        pub const gemv_n = level2.gemv_n;
        pub const gemv_t = level2.gemv_t;
        pub const symv_L = level2.symv_L;
        pub const symv_U = level2.symv_U;

        pub const W = U * lanes(F);
        pub const V = @Vector(W, F);
        /// Four strided elements, for four sums.
        const G = @Vector(4, F);

        pub inline fn ld(p: [*]const F) V {
            return p[0..W].*;
        }
        pub inline fn st(p: [*]F, v: V) void {
            p[0..W].* = v;
        }
        pub inline fn gather(p: [*]F, inc: BlasLong, i: BlasLong) G {
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
        pub const DG = @Vector(4, D);
        /// `dotSum`'s and `squares`' accumulator width: `W`, or half of it where
        /// their two `D` accumulators would fill more than half the vector
        /// registers (T21: Zen's `sdot` held all 16 ymm in two `@Vector(32, f64)`).
        pub const AW = blk: {
            const bits = lanes(F) * @bitSizeOf(F);
            const acc_regs = 2 * W * @bitSizeOf(D) / bits;
            break :blk if (D != F and W >= 2 and 2 * acc_regs > util.registers(bits, sysinfo.target.cpu)) W / 2 else W;
        };
        pub const AV = @Vector(AW, D);

        pub inline fn ldA(p: [*]const F) @Vector(AW, F) {
            return p[0..AW].*;
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

        /// A native vector of `F`.
        pub const L = lanes(F);
        pub const N = @Vector(L, F);

        /// `v`'s upper half added to its lower until `n` lanes are left: a tree
        /// `log2` adds deep, where `@reduce`'s ordered sum chains one add a
        /// lane, which at n = 100 was most of a `sum` (issue #594).
        pub fn halved(comptime n: comptime_int, v: anytype) @Vector(n, @typeInfo(@TypeOf(v)).vector.child) {
            const E = @typeInfo(@TypeOf(v)).vector.child;
            const m = @typeInfo(@TypeOf(v)).vector.len;
            if (m == n) return v;
            const h = m / 2;
            const lo = @shuffle(E, v, undefined, std.simd.iota(i32, h));
            const hi = @shuffle(E, v, undefined, std.simd.iota(i32, h) + @as(@Vector(h, i32), @splat(h)));
            return halved(n, lo + hi);
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
    };
}
