//! The complex `@Vector` kernels (`vector.zig`).

const std = @import("std");
const types = @import("../types.zig");
const util = @import("../util.zig");
const extremum = @import("../blas1/extremum.zig");
const shared = @import("util.zig");
const Real = @import("real.zig").Real;
const complex_reductions = @import("complex/reductions.zig");
const complex_extremes = @import("complex/extremes.zig");
const complex_level2 = @import("complex/level2.zig");

const BlasLong = types.BlasLong;
const Find = extremum.Find;
const lanes = shared.lanes;
const at = shared.at;
const abs = shared.abs;
const mulAdd = shared.mulAdd;
const Pick = shared.Pick;
const root = shared.root;

/// The complex kernels over interleaved `F` pairs.
pub fn Complex(comptime F: type, comptime U: comptime_int) type {
    return struct {
        const reductions = complex_reductions.Reductions(F, U);
        pub const Dots = reductions.Dots;
        pub const asum = reductions.asum;
        pub const asum_wide = reductions.asum_wide;
        pub const dotc = reductions.dotc;
        pub const dotc_wide = reductions.dotc_wide;
        pub const dotu = reductions.dotu;
        pub const dotu_wide = reductions.dotu_wide;
        pub const nrm2 = reductions.nrm2;
        pub const sum = reductions.sum;
        pub const sum_wide = reductions.sum_wide;

        const extremes = complex_extremes.Extremes(F, U);
        pub const amax = extremes.amax;
        pub const amin = extremes.amin;
        pub const iamax = extremes.iamax;
        pub const iamin = extremes.iamin;

        const level2 = complex_level2.Level2(F, U);
        pub const gemv_c = level2.gemv_c;
        pub const gemv_d = level2.gemv_d;
        pub const gemv_n = level2.gemv_n;
        pub const gemv_o = level2.gemv_o;
        pub const gemv_r = level2.gemv_r;
        pub const gemv_s = level2.gemv_s;
        pub const gemv_t = level2.gemv_t;
        pub const gemv_u = level2.gemv_u;
        pub const hemv_L = level2.hemv_L;
        pub const hemv_M = level2.hemv_M;
        pub const hemv_U = level2.hemv_U;
        pub const hemv_V = level2.hemv_V;
        pub const symv_L = level2.symv_L;
        pub const symv_U = level2.symv_U;

        pub const R = Real(F, U);
        /// Floats a vector; whole elements.
        pub const W = if (R.W % 2 == 0) R.W else 2 * R.W;
        pub const V = @Vector(W, F);
        /// One element.
        pub const V2 = @Vector(2, F);
        pub const Result = types.CComplex(F);

        pub inline fn ld(p: [*]const F) V {
            return p[0..W].*;
        }
        pub inline fn st(p: [*]F, v: V) void {
            p[0..W].* = v;
        }
        pub inline fn ld2(p: [*]F, i: BlasLong) V2 {
            return @as(*[2]F, @ptrCast(at(F, p, 2 * i))).*;
        }
        pub inline fn st2(p: [*]F, i: BlasLong, v: V2) void {
            @as(*[2]F, @ptrCast(at(F, p, 2 * i))).* = v;
        }

        /// `(re, im, re, im, ...)` as a `T`, a vector of `F` or of `R.D`.
        pub inline fn pair(comptime T: type, re: F, im: F) T {
            var r: [@typeInfo(T).vector.len]@typeInfo(T).vector.child = undefined;
            inline for (0..r.len) |k| r[k] = if (k % 2 == 0) re else im;
            return r;
        }

        /// Each element's parts exchanged.
        pub inline fn swapped(comptime T: type, v: T) T {
            const len = @typeInfo(T).vector.len;
            const mask = comptime blk: {
                var m: [len]i32 = undefined;
                for (0..len) |k| m[k] = @intCast(k ^ 1);
                break :blk m;
            };
            return @shuffle(@typeInfo(T).vector.child, v, undefined, mask);
        }

        /// Multiplying by `(tr, ti)`, the operand conjugated for `conj`: lane
        /// factors on it and on it with its parts exchanged.
        pub fn By(comptime T: type) type {
            return struct {
                r: T,
                i: T,

                pub inline fn init(tr: F, ti: F, comptime conj: bool) @This() {
                    return .{ .r = pair(T, tr, if (conj) -tr else tr), .i = pair(T, if (conj) ti else -ti, ti) };
                }
                /// `y + t*v`.
                pub inline fn madd(b: @This(), v: T, y: T) T {
                    return mulAdd(T, b.r, v, mulAdd(T, b.i, swapped(T, v), y));
                }
                /// `t*v`.
                pub inline fn mul(b: @This(), v: T) T {
                    return mulAdd(T, b.i, swapped(T, v), b.r * v);
                }
            };
        }

        fn Axpy(comptime conj: bool) type {
            return struct {
                pub fn k(n: BlasLong, _: BlasLong, _: BlasLong, ar: F, ai: F, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F, _: BlasLong) callconv(.c) i32 {
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
    };
}
