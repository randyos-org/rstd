//! The real `@Vector` level-2 kernels: `gemv` both ways and `symv` from either triangle.

const types = @import("../../types.zig");
const shared = @import("../util.zig");
const real = @import("../real.zig");

const BlasLong = types.BlasLong;
const at = shared.at;
const mulAdd = shared.mulAdd;

pub fn Level2(comptime F: type, comptime U: comptime_int) type {
    const Base = real.Real(F, U);
    return struct {
        const V = Base.V;
        const W = Base.W;
        const ld = Base.ld;
        const st = Base.st;

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

        /// `y += alpha*S*x`, `S` symmetric from one stored triangle, over the
        /// columns `offset` selects: the last `offset` for the upper
        /// triangle, the first for the lower. `kernel/riscv64/symv_L_rvv.c`,
        /// `symv_U_rvv.c`, `symv_L_vector.c`, `symv_U_vector.c`: per column
        /// `j`, `alpha*x[j]` times the column into `y`, its dot with `x`
        /// times `alpha` into `y[j]`; the diagonal first.
        fn Symv(comptime lower: bool) type {
            return struct {
                pub fn k(m: BlasLong, offset: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F) callconv(.c) i32 {
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
