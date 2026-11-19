//! The complex `@Vector` kernels (`vector.zig`).

const std = @import("std");
const types = @import("../types.zig");
const util = @import("../util.zig");
const extremum = @import("../blas1/extremum.zig");
const shared = @import("util.zig");
const Real = @import("real.zig").Real;

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
