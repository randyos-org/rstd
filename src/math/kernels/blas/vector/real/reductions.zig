//! The real `@Vector` reductions: the dot products, the sums of squares and `nrm2`, and the plain, absolute and `.fast` sums.

const types = @import("../../types.zig");
const shared = @import("../util.zig");
const real = @import("../real.zig");

const BlasLong = types.BlasLong;
const at = shared.at;
const abs = shared.abs;
const mulAdd = shared.mulAdd;
const root = shared.root;

pub fn Reductions(comptime F: type, comptime U: comptime_int) type {
    const Base = real.Real(F, U);
    return struct {
        const AV = Base.AV;
        const AW = Base.AW;
        const D = Base.D;
        const DG = Base.DG;
        const L = Base.L;
        const N = Base.N;
        const V = Base.V;
        const W = Base.W;
        const gather = Base.gather;
        const halved = Base.halved;
        const ld = Base.ld;
        const ldA = Base.ldA;

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
            // no fold of untouched sums, which are 0 (issue #594: n = 16)
            var s: D = if (i == 0) 0 else halved(1, acc[0] + acc[1])[0];
            while (i < len) : (i += 1) s += abs(@as(D, x[i]));
            return s;
        }

        /// `.fast`'s `asum` and `sum` (issue #584), in `F`'s own lanes, as
        /// `blas1/fast.zig` sums: float `i` of each `2 * W` into lane `i % W`
        /// of sum `i / W % 2`, a last `W` into the first, then native vectors
        /// into a third; the first two added and halved to `L` lanes, the
        /// third added, halved to one, and the floats left added in order.
        fn laneSum(comptime absolute: bool, len: usize, x: [*]F) F {
            var acc: [2]V = @splat(@splat(0));
            var i: usize = 0;
            while (i + 2 * W <= len) : (i += 2 * W) {
                inline for (0..2) |q| {
                    const v = ld(x + i + q * W);
                    acc[q] += if (absolute) abs(v) else v;
                }
            }
            if (i + W <= len) {
                const v = ld(x + i);
                acc[0] += if (absolute) abs(v) else v;
                i += W;
            }
            var t: N = @splat(0);
            while (i + L <= len) : (i += L) {
                const v: N = x[i..][0..L].*;
                t += if (absolute) abs(v) else v;
            }
            var s = halved(1, halved(L, acc[0] + acc[1]) + t)[0];
            while (i < len) : (i += 1) s += if (absolute) abs(x[i]) else x[i];
            return s;
        }

        /// The unit-stride sum of `|x|` (`absolute`) or `x` over `len` floats
        /// in `A`: `D`, `.plain`'s, or `F`, `.fast`'s.
        pub fn unitSum(comptime A: type, comptime absolute: bool, len: usize, x: [*]F) A {
            if (A != D) return laneSum(absolute, len, x);
            return if (absolute) absSum(len, x) else plainSum(len, x);
        }

        /// `asum` (`absolute`) or `sum` in `A`, rounded once; strided,
        /// element `i` into lane `i % 4` of four sums, then the lanes and the
        /// rest in order.
        fn sums(comptime A: type, comptime absolute: bool, n: BlasLong, x: [*]F, incx: BlasLong) F {
            if (n <= 0 or incx <= 0) return 0;
            if (incx == 1) return @floatCast(unitSum(A, absolute, @intCast(n), x));
            const AG = @Vector(4, A);
            var acc: AG = @splat(0);
            var i: BlasLong = 0;
            while (i + 4 <= n) : (i += 4) {
                const g: AG = gather(x, incx, i);
                acc += if (absolute) abs(g) else g;
            }
            var s = @reduce(.Add, acc);
            while (i < n) : (i += 1) {
                const v: A = at(F, x, i * incx).*;
                s += if (absolute) abs(v) else v;
            }
            return @floatCast(s);
        }

        /// `.fast`'s, in `F`'s lanes.
        pub fn asum(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            return sums(F, true, n, x, incx);
        }

        /// `.plain`'s, in `D` (`f64` for `f32`: issue #572).
        pub fn asum_wide(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            return sums(D, true, n, x, incx);
        }

        /// The sum of `x[0..len]` in `D`: float `i` of each
        /// `2 * AW` into lane `i % AW` of sum `i / AW % 2`; then the two lane by
        /// lane, halved to one (`halved`), and the floats left in order.
        pub fn plainSum(len: usize, x: [*]F) D {
            var acc: [2]AV = @splat(@splat(0));
            var i: usize = 0;
            while (i + 2 * AW <= len) : (i += 2 * AW) {
                inline for (0..2) |q| acc[q] += @as(AV, ldA(x + i + q * AW));
            }
            var s: D = if (i == 0) 0 else halved(1, acc[0] + acc[1])[0];
            while (i < len) : (i += 1) s += x[i];
            return s;
        }

        /// `.fast`'s, in `F`'s lanes.
        pub fn sum(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            return sums(F, false, n, x, incx);
        }

        /// `.plain`'s, in `D`.
        pub fn sum_wide(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            return sums(D, false, n, x, incx);
        }
    };
}
