//! The complex `@Vector` reductions: the dot products, plain and conjugated, in `F` and wide, `nrm2`, and the absolute and plain sums.

const types = @import("../../types.zig");
const shared = @import("../util.zig");
const complex = @import("../complex.zig");

const BlasLong = types.BlasLong;
const abs = shared.abs;
const mulAdd = shared.mulAdd;
const root = shared.root;

pub fn Reductions(comptime F: type, comptime U: comptime_int) type {
    const Base = complex.Complex(F, U);
    return struct {
        const R = Base.R;
        const Result = Base.Result;
        const V = Base.V;
        const V2 = Base.V2;
        const W = Base.W;
        const ld = Base.ld;
        const ld2 = Base.ld2;
        const pair = Base.pair;
        const swapped = Base.swapped;

        /// The sums of `a*x` lane by lane and of `a*x` with `x`'s parts
        /// exchanged, as `Dots.finish` makes `sum(op(a)*op(x))` of them.
        pub fn Dots(comptime ca: bool, comptime cx: bool) type {
            return struct {
                rr: V = @splat(0),
                ri: V = @splat(0),
                rr2: V2 = @splat(0),
                ri2: V2 = @splat(0),

                pub inline fn add(d: *@This(), av: V, xv: V) void {
                    d.rr = mulAdd(V, av, xv, d.rr);
                    d.ri = mulAdd(V, av, swapped(V, xv), d.ri);
                }
                pub inline fn add2(d: *@This(), av: V2, xv: V2) void {
                    d.rr2 = mulAdd(V2, av, xv, d.rr2);
                    d.ri2 = mulAdd(V2, av, swapped(V2, xv), d.ri2);
                }
                /// `re = ar*xr - sa*sx*ai*xi`, `im = sx*ar*xi + sa*ai*xr`,
                /// `sa` and `sx` the signs conjugation gives the imaginary parts.
                pub fn finish(d: @This()) V2 {
                    const sa: F = if (ca) -1 else 1;
                    const sx: F = if (cx) -1 else 1;
                    const re = @reduce(.Add, d.rr * pair(V, 1, -sa * sx)) + @reduce(.Add, d.rr2 * pair(V2, 1, -sa * sx));
                    const im = @reduce(.Add, d.ri * pair(V, sx, sa)) + @reduce(.Add, d.ri2 * pair(V2, sx, sa));
                    return .{ re, im };
                }
            };
        }

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

        /// `.fast`'s, in `F`'s lanes.
        pub fn dotu(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) Result {
            return dot(false, n, x, incx, y, incy);
        }
        pub fn dotc(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) Result {
            return dot(true, n, x, incx, y, incy);
        }

        /// Floats a vector of `wideDot`'s sums: `R.AW`, whole elements.
        const DW = if (R.AW % 2 == 0) R.AW else 2 * R.AW;

        /// `Dots(conj, false)`'s sums in `R.D` (`f64` for `f32`, whose
        /// products are exact there), `DW` floats a vector.
        fn WideDots(comptime conj: bool) type {
            const A = R.D;
            const AV = @Vector(DW, A);
            const A2 = @Vector(2, A);
            return struct {
                rr: AV = @splat(0),
                ri: AV = @splat(0),
                rr2: A2 = @splat(0),
                ri2: A2 = @splat(0),

                pub inline fn add(d: *@This(), av: @Vector(DW, F), xv: @Vector(DW, F)) void {
                    const a: AV = av;
                    const x: AV = xv;
                    d.rr = mulAdd(AV, a, x, d.rr);
                    d.ri = mulAdd(AV, a, swapped(AV, x), d.ri);
                }
                pub inline fn add2(d: *@This(), av: V2, xv: V2) void {
                    const a: A2 = av;
                    const x: A2 = xv;
                    d.rr2 = mulAdd(A2, a, x, d.rr2);
                    d.ri2 = mulAdd(A2, a, swapped(A2, x), d.ri2);
                }
                /// As `Dots.finish`, `x` conjugated for `conj`.
                pub fn finish(d: @This()) A2 {
                    const sa: F = if (conj) -1 else 1;
                    const re = @reduce(.Add, d.rr * pair(AV, 1, -sa)) + @reduce(.Add, d.rr2 * pair(A2, 1, -sa));
                    const im = @reduce(.Add, d.ri * pair(AV, 1, sa)) + @reduce(.Add, d.ri2 * pair(A2, 1, sa));
                    return .{ re, im };
                }
            };
        }

        /// `dot`'s order with `WideDots`, rounded once (issue #584).
        fn wideDot(comptime conj: bool, n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) Result {
            var d: [2]WideDots(conj) = .{ .{}, .{} };
            var i: BlasLong = 0;
            if (n > 0 and incx == 1 and incy == 1) {
                const len: usize = @intCast(2 * n);
                var f: usize = 0;
                while (f + 2 * DW <= len) : (f += 2 * DW) {
                    inline for (0..2) |q| d[q].add((x + f + q * DW)[0..DW].*, (y + f + q * DW)[0..DW].*);
                }
                i = @intCast(f / 2);
            }
            while (i + 2 <= n) : (i += 2) {
                inline for (0..2) |q| d[q].add2(ld2(x, (i + q) * incx), ld2(y, (i + q) * incy));
            }
            if (i < n) d[0].add2(ld2(x, i * incx), ld2(y, i * incy));
            const s = d[0].finish() + d[1].finish();
            return .{ .re = @floatCast(s[0]), .im = @floatCast(s[1]) };
        }

        /// `.plain`'s, in `R.D`.
        pub fn dotu_wide(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) Result {
            return wideDot(false, n, x, incx, y, incy);
        }
        pub fn dotc_wide(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) Result {
            return wideDot(true, n, x, incx, y, incy);
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

        /// The sum of `|re| + |im|` (`absolute`) or of `re + im` in `A`:
        /// `Real.unitSum`'s over the `2n` floats; strided, element `i` into
        /// sum `i % 2`, by parts, then those lane by lane and the parts in order.
        fn sums(comptime A: type, comptime absolute: bool, n: BlasLong, x: [*]F, incx: BlasLong) F {
            if (n <= 0 or incx <= 0) return 0;
            if (incx == 1) return @floatCast(R.unitSum(A, absolute, @intCast(2 * n), x));
            const A2 = @Vector(2, A);
            var acc: [2]A2 = @splat(@splat(0));
            var i: BlasLong = 0;
            while (i + 2 <= n) : (i += 2) {
                inline for (0..2) |q| {
                    const v: A2 = ld2(x, (i + q) * incx);
                    acc[q] += if (absolute) abs(v) else v;
                }
            }
            if (i < n) {
                const v: A2 = ld2(x, i * incx);
                acc[0] += if (absolute) abs(v) else v;
            }
            return @floatCast(@reduce(.Add, acc[0] + acc[1]));
        }

        /// `.fast`'s, in `F`'s lanes.
        pub fn asum(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            return sums(F, true, n, x, incx);
        }

        /// `.plain`'s, in `R.D` (`f64` for `f32`: issue #572).
        pub fn asum_wide(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            return sums(R.D, true, n, x, incx);
        }

        /// `.fast`'s, in `F`'s lanes.
        pub fn sum(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            return sums(F, false, n, x, incx);
        }

        /// `.plain`'s, in `R.D`.
        pub fn sum_wide(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            return sums(R.D, false, n, x, incx);
        }
    };
}
