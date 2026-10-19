//! `@Vector` candidates for level-1 and level-2 kernels, in the reference's
//! kernel ABI so they race in a vtable slot's place (`blas_bench.zig`). One
//! step is `U` native vectors: `U = 1` is the plain kernel, `U = 4` the wider
//! unroll. Reductions keep two such accumulators. Unit strides take the
//! vector path, other positive strides a scalar loop; negative strides are
//! not handled (bench only).

const std = @import("std");

/// The real kernels over `F`.
pub fn Real(comptime F: type, comptime U: comptime_int) type {
    return struct {
        const L = std.simd.suggestVectorLength(F) orelse 1;
        pub const W = U * L;
        const V = @Vector(W, F);

        inline fn ld(p: [*]const F) V {
            return p[0..W].*;
        }
        inline fn st(p: [*]F, v: V) void {
            p[0..W].* = v;
        }
        inline fn ustride(inc: isize) usize {
            return @intCast(inc);
        }

        /// Four strided elements from `i` on, for four accumulators.
        const G = @Vector(4, F);
        inline fn gather(p: [*]const F, s: usize, i: usize) G {
            return .{ p[i * s], p[(i + 1) * s], p[(i + 2) * s], p[(i + 3) * s] };
        }

        pub fn axpy(n: isize, _: isize, _: isize, a: F, x: [*]F, incx: isize, y: [*]F, incy: isize, _: [*]F, _: isize) callconv(.c) i32 {
            if (n <= 0 or incx <= 0 or incy <= 0) return 0;
            const len: usize = @intCast(n);
            var i: usize = 0;
            if (incx == 1 and incy == 1) {
                const av: V = @splat(a);
                while (i + W <= len) : (i += W) st(y + i, @mulAdd(V, av, ld(x + i), ld(y + i)));
                while (i < len) : (i += 1) y[i] = @mulAdd(F, a, x[i], y[i]);
            } else {
                const sx = ustride(incx);
                const sy = ustride(incy);
                while (i < len) : (i += 1) y[i * sy] = @mulAdd(F, a, x[i * sx], y[i * sy]);
            }
            return 0;
        }

        pub fn dot(n: isize, x: [*]F, incx: isize, y: [*]F, incy: isize) callconv(.c) F {
            if (n <= 0 or incx <= 0 or incy <= 0) return 0;
            const len: usize = @intCast(n);
            var i: usize = 0;
            var s: F = 0;
            if (incx == 1 and incy == 1) {
                var a0: V = @splat(0);
                var a1: V = @splat(0);
                while (i + 2 * W <= len) : (i += 2 * W) {
                    a0 = @mulAdd(V, ld(x + i), ld(y + i), a0);
                    a1 = @mulAdd(V, ld(x + i + W), ld(y + i + W), a1);
                }
                if (i + W <= len) {
                    a0 = @mulAdd(V, ld(x + i), ld(y + i), a0);
                    i += W;
                }
                s = @reduce(.Add, a0 + a1);
                while (i < len) : (i += 1) s = @mulAdd(F, x[i], y[i], s);
            } else {
                const sx = ustride(incx);
                const sy = ustride(incy);
                var acc: G = @splat(0);
                while (i + 4 <= len) : (i += 4) acc = @mulAdd(G, gather(x, sx, i), gather(y, sy, i), acc);
                s = @reduce(.Add, acc);
                while (i < len) : (i += 1) s = @mulAdd(F, x[i * sx], y[i * sy], s);
            }
            return s;
        }

        pub fn scal(n: isize, _: isize, _: isize, a: F, x: [*]F, incx: isize, _: [*]F, _: isize, _: [*]F, flag: isize) callconv(.c) i32 {
            if (n <= 0 or incx <= 0) return 0;
            const len: usize = @intCast(n);
            // a zero alpha zeroes, unless the flag asks NaN and Inf to stay NaN
            const zero = a == 0 and flag == 0;
            var i: usize = 0;
            if (incx == 1) {
                const av: V = @splat(a);
                while (i + W <= len) : (i += W) st(x + i, if (zero) @splat(0) else av * ld(x + i));
                while (i < len) : (i += 1) x[i] = if (zero) 0 else a * x[i];
            } else {
                const sx = ustride(incx);
                while (i < len) : (i += 1) x[i * sx] = if (zero) 0 else a * x[i * sx];
            }
            return 0;
        }

        /// Unscaled sum of squares; the scaled loop only where it could have
        /// overflowed or lost to underflow.
        pub fn nrm2(n: isize, x: [*]F, incx: isize) callconv(.c) F {
            if (n <= 0 or incx <= 0) return 0;
            const len: usize = @intCast(n);
            var i: usize = 0;
            var s: F = 0;
            if (incx == 1) {
                var a0: V = @splat(0);
                var a1: V = @splat(0);
                while (i + 2 * W <= len) : (i += 2 * W) {
                    const v0 = ld(x + i);
                    const v1 = ld(x + i + W);
                    a0 = @mulAdd(V, v0, v0, a0);
                    a1 = @mulAdd(V, v1, v1, a1);
                }
                s = @reduce(.Add, a0 + a1);
                while (i < len) : (i += 1) s = @mulAdd(F, x[i], x[i], s);
            } else {
                const sx = ustride(incx);
                var acc: G = @splat(0);
                while (i + 4 <= len) : (i += 4) {
                    const v = gather(x, sx, i);
                    acc = @mulAdd(G, v, v, acc);
                }
                s = @reduce(.Add, acc);
                while (i < len) : (i += 1) s = @mulAdd(F, x[i * sx], x[i * sx], s);
            }
            const lo = std.math.floatMin(F) / std.math.floatEps(F) * @as(F, @floatFromInt(len));
            if (s >= lo and s <= std.math.floatMax(F)) return @sqrt(s);
            return scaledNrm2(len, x, ustride(incx));
        }

        fn scaledNrm2(len: usize, x: [*]F, sx: usize) F {
            var scale: F = 0;
            var ssq: F = 1;
            for (0..len) |i| {
                const v = @abs(x[i * sx]);
                if (v == 0) continue;
                if (scale < v) {
                    const r = scale / v;
                    ssq = 1 + ssq * r * r;
                    scale = v;
                } else {
                    const r = v / scale;
                    ssq += r * r;
                }
            }
            return scale * @sqrt(ssq);
        }

        pub fn asum(n: isize, x: [*]F, incx: isize) callconv(.c) F {
            if (n <= 0 or incx <= 0) return 0;
            const len: usize = @intCast(n);
            var i: usize = 0;
            var s: F = 0;
            if (incx == 1) {
                var a0: V = @splat(0);
                var a1: V = @splat(0);
                while (i + 2 * W <= len) : (i += 2 * W) {
                    a0 += @abs(ld(x + i));
                    a1 += @abs(ld(x + i + W));
                }
                s = @reduce(.Add, a0 + a1);
                while (i < len) : (i += 1) s += @abs(x[i]);
            } else {
                const sx = ustride(incx);
                var acc: G = @splat(0);
                while (i + 4 <= len) : (i += 4) acc += @abs(gather(x, sx, i));
                s = @reduce(.Add, acc);
                while (i < len) : (i += 1) s += @abs(x[i * sx]);
            }
            return s;
        }

        /// The first largest `|x|`, from one (`firstMax`).
        pub fn iamax(n: isize, x: [*]F, incx: isize) callconv(.c) isize {
            if (n <= 0 or incx <= 0) return 0;
            return firstMax(@intCast(n), x, ustride(incx), struct {
                inline fn key(p: [*]const F, s: usize, i: usize) F {
                    return @abs(p[i * s]);
                }
                inline fn block(p: [*]const F, s: usize, i: usize) F {
                    var m: [4]V = undefined;
                    if (s == 1) {
                        inline for (0..4) |q| m[q] = @abs(ld(p + i + q * W));
                    } else {
                        inline for (0..4) |q| {
                            var r: [W]F = undefined;
                            inline for (0..W) |k| r[k] = key(p, s, i + q * W + k);
                            m[q] = r;
                        }
                    }
                    return @reduce(.Max, @max(@max(m[0], m[1]), @max(m[2], m[3])));
                }
            });
        }

        /// The first largest `Key.key(x, s, i)`, from one, `Key.block` giving
        /// the largest of the `4 * W` keys from `i`. One pass keeps the first
        /// block whose maximum is largest, with no index kept per element;
        /// only that block is searched for the index. NaN keys are passed
        /// over, as by the reference's comparisons, but for a NaN first one.
        pub fn firstMax(len: usize, x: [*]const F, s: usize, comptime Key: type) isize {
            const B = 4 * W;
            const k0 = Key.key(x, s, 0);
            if (std.math.isNan(k0)) return 1;
            var best = k0;
            var at: usize = 0;
            var i: usize = 0;
            if (len >= B) {
                var blk: usize = 0;
                while (i + B <= len) : (i += B) {
                    const m = Key.block(x, s, i);
                    if (m > best) {
                        best = m;
                        blk = i;
                    }
                }
                at = blk;
                while (Key.key(x, s, at) != best) at += 1;
            }
            while (i < len) : (i += 1) {
                const k = Key.key(x, s, i);
                if (k > best) {
                    best = k;
                    at = i;
                }
            }
            return @intCast(at + 1);
        }

        pub fn copy(n: isize, x: [*]F, incx: isize, y: [*]F, incy: isize) callconv(.c) i32 {
            if (n <= 0 or incx <= 0 or incy <= 0) return 0;
            const len: usize = @intCast(n);
            var i: usize = 0;
            if (incx == 1 and incy == 1) {
                while (i + W <= len) : (i += W) st(y + i, ld(x + i));
                while (i < len) : (i += 1) y[i] = x[i];
            } else {
                const sx = ustride(incx);
                const sy = ustride(incy);
                while (i < len) : (i += 1) y[i * sy] = x[i * sx];
            }
            return 0;
        }

        pub fn swap(n: isize, _: isize, _: isize, _: F, x: [*]F, incx: isize, y: [*]F, incy: isize, _: [*]F, _: isize) callconv(.c) i32 {
            if (n <= 0 or incx <= 0 or incy <= 0) return 0;
            const len: usize = @intCast(n);
            var i: usize = 0;
            if (incx == 1 and incy == 1) {
                while (i + W <= len) : (i += W) {
                    const t = ld(x + i);
                    st(x + i, ld(y + i));
                    st(y + i, t);
                }
                while (i < len) : (i += 1) std.mem.swap(F, &x[i], &y[i]);
            } else {
                const sx = ustride(incx);
                const sy = ustride(incy);
                while (i < len) : (i += 1) std.mem.swap(F, &x[i * sx], &y[i * sy]);
            }
            return 0;
        }

        pub fn rot(n: isize, x: [*]F, incx: isize, y: [*]F, incy: isize, c: F, s: F) callconv(.c) i32 {
            if (n <= 0 or incx <= 0 or incy <= 0) return 0;
            const len: usize = @intCast(n);
            var i: usize = 0;
            if (incx == 1 and incy == 1) {
                const cv: V = @splat(c);
                const sv: V = @splat(s);
                while (i + W <= len) : (i += W) {
                    const xv = ld(x + i);
                    const yv = ld(y + i);
                    st(x + i, @mulAdd(V, cv, xv, sv * yv));
                    st(y + i, @mulAdd(V, cv, yv, -(sv * xv)));
                }
                while (i < len) : (i += 1) rot1(&x[i], &y[i], c, s);
            } else {
                const sx = ustride(incx);
                const sy = ustride(incy);
                while (i < len) : (i += 1) rot1(&x[i * sx], &y[i * sy], c, s);
            }
            return 0;
        }

        inline fn rot1(xp: *F, yp: *F, c: F, s: F) void {
            const xv = xp.*;
            const yv = yp.*;
            xp.* = @mulAdd(F, c, xv, s * yv);
            yp.* = @mulAdd(F, c, yv, -(s * xv));
        }

        /// `y += alpha*A*x`, four columns per pass over `y`.
        pub fn gemv_n(m: isize, n: isize, _: isize, alpha: F, a: [*]F, lda: isize, x: [*]F, incx: isize, y: [*]F, incy: isize, _: [*]F) callconv(.c) i32 {
            if (m <= 0 or n <= 0 or incx <= 0 or incy <= 0) return 0;
            const rows: usize = @intCast(m);
            const cols: usize = @intCast(n);
            const ld_a: usize = @intCast(lda);
            const sx = ustride(incx);
            const sy = ustride(incy);
            var j: usize = 0;
            if (sy == 1) {
                while (j + 4 <= cols) : (j += 4) {
                    var t: [4]F = undefined;
                    var c: [4][*]F = undefined;
                    inline for (0..4) |q| {
                        t[q] = alpha * x[(j + q) * sx];
                        c[q] = a + (j + q) * ld_a;
                    }
                    var i: usize = 0;
                    while (i + W <= rows) : (i += W) {
                        var yv = ld(y + i);
                        inline for (0..4) |q| yv = @mulAdd(V, ld(c[q] + i), @splat(t[q]), yv);
                        st(y + i, yv);
                    }
                    while (i < rows) : (i += 1) {
                        inline for (0..4) |q| y[i] = @mulAdd(F, c[q][i], t[q], y[i]);
                    }
                }
            }
            while (j < cols) : (j += 1) {
                const t = alpha * x[j * sx];
                const c = a + j * ld_a;
                var i: usize = 0;
                if (sy == 1) {
                    while (i + W <= rows) : (i += W) st(y + i, @mulAdd(V, ld(c + i), @splat(t), ld(y + i)));
                }
                while (i < rows) : (i += 1) y[i * sy] = @mulAdd(F, c[i], t, y[i * sy]);
            }
            return 0;
        }

        /// `y += alpha*A'*x`, four columns' dot products per pass over `x`.
        pub fn gemv_t(m: isize, n: isize, _: isize, alpha: F, a: [*]F, lda: isize, x: [*]F, incx: isize, y: [*]F, incy: isize, _: [*]F) callconv(.c) i32 {
            if (m <= 0 or n <= 0 or incx <= 0 or incy <= 0) return 0;
            const rows: usize = @intCast(m);
            const cols: usize = @intCast(n);
            const ld_a: usize = @intCast(lda);
            const sx = ustride(incx);
            const sy = ustride(incy);
            var j: usize = 0;
            if (sx == 1) {
                while (j + 4 <= cols) : (j += 4) {
                    var c: [4][*]F = undefined;
                    var acc: [4]V = @splat(@splat(0));
                    inline for (0..4) |q| c[q] = a + (j + q) * ld_a;
                    var i: usize = 0;
                    while (i + W <= rows) : (i += W) {
                        const xv = ld(x + i);
                        inline for (0..4) |q| acc[q] = @mulAdd(V, ld(c[q] + i), xv, acc[q]);
                    }
                    inline for (0..4) |q| {
                        var s = @reduce(.Add, acc[q]);
                        var r = i;
                        while (r < rows) : (r += 1) s = @mulAdd(F, c[q][r], x[r], s);
                        y[(j + q) * sy] = @mulAdd(F, alpha, s, y[(j + q) * sy]);
                    }
                }
            }
            while (j < cols) : (j += 1) {
                const c = a + j * ld_a;
                var s: F = 0;
                var i: usize = 0;
                if (sx == 1) {
                    var acc: V = @splat(0);
                    while (i + W <= rows) : (i += W) acc = @mulAdd(V, ld(c + i), ld(x + i), acc);
                    s = @reduce(.Add, acc);
                }
                while (i < rows) : (i += 1) s = @mulAdd(F, c[i], x[i * sx], s);
                y[j * sy] = @mulAdd(F, alpha, s, y[j * sy]);
            }
            return 0;
        }
    };
}

/// The complex kernels over interleaved `F` pairs; `C` is the `extern`
/// struct a complex dot returns.
pub fn Complex(comptime F: type, comptime U: comptime_int, comptime C: type) type {
    return struct {
        const R = Real(F, U);
        const W = R.W;
        const V = @Vector(W, F);

        inline fn ld(p: [*]const F) V {
            return p[0..W].*;
        }
        inline fn st(p: [*]F, v: V) void {
            p[0..W].* = v;
        }

        /// Each pair's parts exchanged.
        inline fn swapped(v: V) V {
            const mask = comptime blk: {
                var m: [W]i32 = undefined;
                for (0..W) |k| m[k] = @intCast(k ^ 1);
                break :blk m;
            };
            return @shuffle(F, v, undefined, mask);
        }

        /// `(-v, v, -v, v, ...)`.
        inline fn signed(v: F) V {
            var r: [W]F = undefined;
            inline for (0..W) |k| r[k] = if (k % 2 == 0) -v else v;
            return r;
        }

        /// The even (`odd = false`) or odd elements of `a` then `b`.
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

        inline fn ustride(inc: isize) usize {
            return @intCast(inc);
        }

        /// One element as a vector, for the strided loops.
        const V2 = @Vector(2, F);
        inline fn ld2(p: [*]const F) V2 {
            return p[0..2].*;
        }
        inline fn st2(p: [*]F, v: V2) void {
            p[0..2].* = v;
        }
        inline fn swap2(v: V2) V2 {
            return @shuffle(F, v, undefined, [2]i32{ 1, 0 });
        }

        /// One element times `(ar, ai)`, added to `(yr, yi)`.
        inline fn cmla(yr: *F, yi: *F, ar: F, ai: F, xr: F, xi: F) void {
            const r = yr.* + (ar * xr - ai * xi);
            yi.* = yi.* + (ar * xi + ai * xr);
            yr.* = r;
        }

        pub fn axpy(n: isize, _: isize, _: isize, ar: F, ai: F, x: [*]F, incx: isize, y: [*]F, incy: isize, _: [*]F, _: isize) callconv(.c) i32 {
            if (n <= 0 or incx <= 0 or incy <= 0) return 0;
            const len: usize = @intCast(n);
            var i: usize = 0;
            if (incx == 1 and incy == 1) {
                const vr: V = @splat(ar);
                const vi = signed(ai);
                while (2 * i + W <= 2 * len) : (i += W / 2) {
                    const xv = ld(x + 2 * i);
                    st(y + 2 * i, @mulAdd(V, vi, swapped(xv), @mulAdd(V, vr, xv, ld(y + 2 * i))));
                }
                while (i < len) : (i += 1) cmla(&y[2 * i], &y[2 * i + 1], ar, ai, x[2 * i], x[2 * i + 1]);
            } else {
                const sx = 2 * ustride(incx);
                const sy = 2 * ustride(incy);
                const vr: V2 = @splat(ar);
                const vi: V2 = .{ -ai, ai };
                while (i < len) : (i += 1) {
                    const xv = ld2(x + i * sx);
                    st2(y + i * sy, @mulAdd(V2, vi, swap2(xv), @mulAdd(V2, vr, xv, ld2(y + i * sy))));
                }
            }
            return 0;
        }

        pub fn dotu(n: isize, x: [*]F, incx: isize, y: [*]F, incy: isize) callconv(.c) C {
            var re: F = 0;
            var im: F = 0;
            if (n <= 0 or incx <= 0 or incy <= 0) return .{ .re = 0, .im = 0 };
            const len: usize = @intCast(n);
            var i: usize = 0;
            if (incx == 1 and incy == 1) {
                // lane products: (xr*yr, xi*yi) and (xr*yi, xi*yr)
                var rr: [2]V = @splat(@splat(0));
                var ri: [2]V = @splat(@splat(0));
                while (2 * i + 2 * W <= 2 * len) : (i += W) {
                    inline for (0..2) |q| {
                        const xv = ld(x + 2 * i + q * W);
                        const yv = ld(y + 2 * i + q * W);
                        rr[q] = @mulAdd(V, xv, yv, rr[q]);
                        ri[q] = @mulAdd(V, xv, swapped(yv), ri[q]);
                    }
                }
                re = -@reduce(.Add, (rr[0] + rr[1]) * signed(1));
                im = @reduce(.Add, ri[0] + ri[1]);
                while (i < len) : (i += 1) cmla(&re, &im, x[2 * i], x[2 * i + 1], y[2 * i], y[2 * i + 1]);
            } else {
                const sx = 2 * ustride(incx);
                const sy = 2 * ustride(incy);
                var rr: [2]V2 = @splat(@splat(0));
                var ri: [2]V2 = @splat(@splat(0));
                while (i + 2 <= len) : (i += 2) {
                    inline for (0..2) |q| {
                        const xv = ld2(x + (i + q) * sx);
                        const yv = ld2(y + (i + q) * sy);
                        rr[q] = @mulAdd(V2, xv, yv, rr[q]);
                        ri[q] = @mulAdd(V2, xv, swap2(yv), ri[q]);
                    }
                }
                const r = rr[0] + rr[1];
                re = r[0] - r[1];
                im = @reduce(.Add, ri[0] + ri[1]);
                while (i < len) : (i += 1) cmla(&re, &im, x[i * sx], x[i * sx + 1], y[i * sy], y[i * sy + 1]);
            }
            return .{ .re = re, .im = im };
        }

        pub fn scal(n: isize, _: isize, _: isize, ar: F, ai: F, x: [*]F, incx: isize, _: [*]F, _: isize, _: [*]F, flag: isize) callconv(.c) i32 {
            if (n <= 0 or incx <= 0) return 0;
            const len: usize = @intCast(n);
            const zero = ar == 0 and ai == 0 and flag == 0;
            var i: usize = 0;
            const sx: usize = if (incx == 1) 2 else 2 * ustride(incx);
            if (incx == 1 and !zero) {
                const vr: V = @splat(ar);
                const vi = signed(ai);
                while (2 * i + W <= 2 * len) : (i += W / 2) {
                    const xv = ld(x + 2 * i);
                    st(x + 2 * i, @mulAdd(V, vi, swapped(xv), vr * xv));
                }
            } else if (!zero) {
                const vr: V2 = @splat(ar);
                const vi: V2 = .{ -ai, ai };
                while (i < len) : (i += 1) {
                    const xv = ld2(x + i * sx);
                    st2(x + i * sx, @mulAdd(V2, vi, swap2(xv), vr * xv));
                }
            }
            while (i < len) : (i += 1) {
                const xr = x[i * sx];
                const xi = x[i * sx + 1];
                x[i * sx] = if (zero) 0 else ar * xr - ai * xi;
                x[i * sx + 1] = if (zero) 0 else ar * xi + ai * xr;
            }
            return 0;
        }

        pub fn nrm2(n: isize, x: [*]F, incx: isize) callconv(.c) F {
            if (n <= 0 or incx <= 0) return 0;
            if (incx == 1) return R.nrm2(2 * n, x, 1);
            const len: usize = @intCast(n);
            const sx = 2 * ustride(incx);
            var acc: [2]V2 = @splat(@splat(0));
            var i: usize = 0;
            while (i + 2 <= len) : (i += 2) {
                inline for (0..2) |q| {
                    const v = ld2(x + (i + q) * sx);
                    acc[q] = @mulAdd(V2, v, v, acc[q]);
                }
            }
            var s = @reduce(.Add, acc[0] + acc[1]);
            while (i < len) : (i += 1) s += x[i * sx] * x[i * sx] + x[i * sx + 1] * x[i * sx + 1];
            const lo = std.math.floatMin(F) / std.math.floatEps(F) * @as(F, @floatFromInt(2 * len));
            if (s >= lo and s <= std.math.floatMax(F)) return @sqrt(s);
            // the parts as two strided real vectors, scaled
            return std.math.hypot(R.nrm2(n, x, 2 * incx), R.nrm2(n, x + 1, 2 * incx));
        }

        pub fn asum(n: isize, x: [*]F, incx: isize) callconv(.c) F {
            if (n <= 0 or incx <= 0) return 0;
            if (incx == 1) return R.asum(2 * n, x, 1);
            const len: usize = @intCast(n);
            const sx = 2 * ustride(incx);
            var acc: [2]V2 = @splat(@splat(0));
            var i: usize = 0;
            while (i + 2 <= len) : (i += 2) {
                inline for (0..2) |q| acc[q] += @abs(ld2(x + (i + q) * sx));
            }
            var s = @reduce(.Add, acc[0] + acc[1]);
            while (i < len) : (i += 1) s += @abs(x[i * sx]) + @abs(x[i * sx + 1]);
            return s;
        }

        /// The first largest `|re| + |im|`, as `Real.iamax`.
        pub fn iamax(n: isize, x: [*]F, incx: isize) callconv(.c) isize {
            if (n <= 0 or incx <= 0) return 0;
            return R.firstMax(@intCast(n), x, 2 * ustride(incx), struct {
                inline fn key(p: [*]const F, s: usize, i: usize) F {
                    return @abs(p[i * s]) + @abs(p[i * s + 1]);
                }
                inline fn block(p: [*]const F, s: usize, i: usize) F {
                    var m: [4]V = undefined;
                    if (s == 2) {
                        inline for (0..4) |q| {
                            const a = @abs(ld(p + 2 * (i + q * W)));
                            const b = @abs(ld(p + 2 * (i + q * W) + W));
                            m[q] = parts(a, b, false) + parts(a, b, true);
                        }
                    } else {
                        inline for (0..4) |q| {
                            var r: [W]F = undefined;
                            inline for (0..W) |k| r[k] = key(p, s, i + q * W + k);
                            m[q] = r;
                        }
                    }
                    return @reduce(.Max, @max(@max(m[0], m[1]), @max(m[2], m[3])));
                }
            });
        }

        pub fn copy(n: isize, x: [*]F, incx: isize, y: [*]F, incy: isize) callconv(.c) i32 {
            if (n <= 0 or incx <= 0 or incy <= 0) return 0;
            if (incx == 1 and incy == 1) return R.copy(2 * n, x, 1, y, 1);
            const sx = 2 * ustride(incx);
            const sy = 2 * ustride(incy);
            for (0..@intCast(n)) |i| st2(y + i * sy, ld2(x + i * sx));
            return 0;
        }

        pub fn swap(n: isize, _: isize, _: isize, _: F, _: F, x: [*]F, incx: isize, y: [*]F, incy: isize, _: [*]F, _: isize) callconv(.c) i32 {
            if (n <= 0 or incx <= 0 or incy <= 0) return 0;
            if (incx == 1 and incy == 1) return R.swap(2 * n, 0, 0, 0, x, 1, y, 1, x, 0);
            const sx = 2 * ustride(incx);
            const sy = 2 * ustride(incy);
            for (0..@intCast(n)) |i| {
                const t = ld2(x + i * sx);
                st2(x + i * sx, ld2(y + i * sy));
                st2(y + i * sy, t);
            }
            return 0;
        }

        pub fn rot(n: isize, x: [*]F, incx: isize, y: [*]F, incy: isize, c: F, s: F) callconv(.c) i32 {
            if (n <= 0 or incx <= 0 or incy <= 0) return 0;
            if (incx == 1 and incy == 1) return R.rot(2 * n, x, 1, y, 1, c, s);
            const sx = 2 * ustride(incx);
            const sy = 2 * ustride(incy);
            const cv: V2 = @splat(c);
            const sv: V2 = @splat(s);
            for (0..@intCast(n)) |i| {
                const xv = ld2(x + i * sx);
                const yv = ld2(y + i * sy);
                st2(x + i * sx, @mulAdd(V2, cv, xv, sv * yv));
                st2(y + i * sy, @mulAdd(V2, cv, yv, -(sv * xv)));
            }
            return 0;
        }

        /// `y += alpha*A*x`, four columns per pass over `y`.
        pub fn gemv_n(m: isize, n: isize, _: isize, ar: F, ai: F, a: [*]F, lda: isize, x: [*]F, incx: isize, y: [*]F, incy: isize, _: [*]F) callconv(.c) i32 {
            if (m <= 0 or n <= 0 or incx <= 0 or incy <= 0) return 0;
            const rows: usize = @intCast(m);
            const cols: usize = @intCast(n);
            const ld_a: usize = 2 * @as(usize, @intCast(lda));
            const sx = 2 * ustride(incx);
            const sy = 2 * ustride(incy);
            var j: usize = 0;
            while (j < cols) {
                const g: usize = if (sy == 2 and j + 4 <= cols) 4 else 1;
                var tr: [4]F = undefined;
                var ti: [4]F = undefined;
                var c: [4][*]F = undefined;
                for (0..g) |q| {
                    const xr = x[(j + q) * sx];
                    const xi = x[(j + q) * sx + 1];
                    tr[q] = ar * xr - ai * xi;
                    ti[q] = ar * xi + ai * xr;
                    c[q] = a + (j + q) * ld_a;
                }
                var i: usize = 0;
                if (g == 4) {
                    while (2 * i + W <= 2 * rows) : (i += W / 2) {
                        var yv = ld(y + 2 * i);
                        inline for (0..4) |q| {
                            const av = ld(c[q] + 2 * i);
                            yv = @mulAdd(V, signed(ti[q]), swapped(av), @mulAdd(V, @splat(tr[q]), av, yv));
                        }
                        st(y + 2 * i, yv);
                    }
                }
                while (i < rows) : (i += 1) {
                    for (0..g) |q| cmla(&y[i * sy], &y[i * sy + 1], tr[q], ti[q], c[q][2 * i], c[q][2 * i + 1]);
                }
                j += g;
            }
            return 0;
        }

        /// `y += alpha*A'*x` (no conjugate), four columns per pass over `x`.
        pub fn gemv_t(m: isize, n: isize, _: isize, ar: F, ai: F, a: [*]F, lda: isize, x: [*]F, incx: isize, y: [*]F, incy: isize, _: [*]F) callconv(.c) i32 {
            if (m <= 0 or n <= 0 or incx <= 0 or incy <= 0) return 0;
            const rows: usize = @intCast(m);
            const cols: usize = @intCast(n);
            const ld_a: usize = 2 * @as(usize, @intCast(lda));
            const sx = 2 * ustride(incx);
            const sy = 2 * ustride(incy);
            var j: usize = 0;
            while (j < cols) {
                const g: usize = if (sx == 2 and j + 4 <= cols) 4 else 1;
                var c: [4][*]F = undefined;
                for (0..g) |q| c[q] = a + (j + q) * ld_a;
                var re: [4]F = @splat(0);
                var im: [4]F = @splat(0);
                var i: usize = 0;
                if (g == 4) {
                    var rr: [4]V = @splat(@splat(0));
                    var ri: [4]V = @splat(@splat(0));
                    while (2 * i + W <= 2 * rows) : (i += W / 2) {
                        const xv = ld(x + 2 * i);
                        const xs = swapped(xv);
                        inline for (0..4) |q| {
                            const av = ld(c[q] + 2 * i);
                            rr[q] = @mulAdd(V, av, xv, rr[q]);
                            ri[q] = @mulAdd(V, av, xs, ri[q]);
                        }
                    }
                    inline for (0..4) |q| {
                        re[q] = -@reduce(.Add, rr[q] * signed(1));
                        im[q] = @reduce(.Add, ri[q]);
                    }
                }
                while (i < rows) : (i += 1) {
                    for (0..g) |q| cmla(&re[q], &im[q], c[q][2 * i], c[q][2 * i + 1], x[i * sx], x[i * sx + 1]);
                }
                for (0..g) |q| cmla(&y[(j + q) * sy], &y[(j + q) * sy + 1], ar, ai, re[q], im[q]);
                j += g;
            }
            return 0;
        }
    };
}
