//! The complex `@Vector` level-2 kernels: `gemv` under every transpose and conjugation, and `symv`, `hemv` from either triangle.

const types = @import("../../types.zig");
const complex = @import("../complex.zig");

const BlasLong = types.BlasLong;

pub fn Level2(comptime F: type, comptime U: comptime_int) type {
    const Base = complex.Complex(F, U);
    return struct {
        const By = Base.By;
        const Dots = Base.Dots;
        const V = Base.V;
        const V2 = Base.V2;
        const W = Base.W;
        const ld = Base.ld;
        const ld2 = Base.ld2;
        const st = Base.st;
        const st2 = Base.st2;

        /// `y += alpha*op(A)*x`: transposed for `trans`, `A` conjugated for
        /// `ca`, `x` for `cx`, as `gemv/`'s portable `Zgemv`. Without a
        /// transpose, four columns a pass over a unit-stride `y`; with one,
        /// four columns' sums a pass over a unit-stride `x`.
        fn Gemv(comptime trans: bool, comptime ca: bool, comptime cx: bool) type {
            return struct {
                pub fn k(m: BlasLong, n: BlasLong, _: BlasLong, ar: F, ai: F, a: [*]F, lda: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F) callconv(.c) i32 {
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
                pub fn k(m: BlasLong, offset: BlasLong, ar: F, ai: F, a: [*]F, lda: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F) callconv(.c) i32 {
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
