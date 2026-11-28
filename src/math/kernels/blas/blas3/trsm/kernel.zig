//! The `trsm` copies (`kernel/generic/trsm_{u,l}{n,t}copy_*.c`, `ztrsm_*`) and
//! kernels (`trsm_kernel_{LT,LN,RN,RT}.c`) over the gemm tile; `solve` divides
//! by the diagonal, where upstream packs its reciprocal (`comparisons.md`, 10).

const std = @import("std");
const types = @import("../../types.zig");
const util = @import("../../util.zig");
const gk = @import("../gemm/kernel.zig");

const BlasLong = types.BlasLong;

pub fn Kernels(comptime T: type, comptime geo: gk.Geometry) type {
    const complex = types.isComplex(T);
    const F = types.Float(T);
    const cc = types.kernelCc(F);
    const E = if (complex) [2]F else F;
    const G = gk.Kernels(T, geo);
    const mr = geo.mr();
    const nr = geo.nr;
    return struct {
        // ---- arithmetic on elements -----------------------------------------

        inline fn mul(x: E, y: E) E {
            if (comptime !complex) return util.scalar(E, x * y);
            return util.scalar(E, .{ x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0] });
        }
        // `c` comptime: a run-time one selects between `-x` and `x`, which
        // LLVM cannot compile in a loop on MIPS64r6 (issue #192)
        fn conjIf(x: E, comptime c: bool) E {
            if (comptime !complex) return x;
            return if (c) .{ x[0], -x[1] } else x;
        }
        inline fn sub(x: E, y: E) E {
            if (comptime !complex) return util.scalar(E, x - y);
            return util.scalar(E, .{ x[0] - y[0], x[1] - y[1] });
        }
        /// `x / y`; complex by `util.cdiv`, as `trsv` divides.
        inline fn div(x: E, y: E) E {
            if (comptime !complex) return x / y;
            return util.cdiv(F, x[0], x[1], y[0], y[1]);
        }

        // ---- copies --------------------------------------------------------

        fn Copy(comptime u: usize, comptime upper: bool, comptime across: bool, comptime unit: bool) type {
            return struct {
                fn f(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, offset: BlasLong, b: [*]F) callconv(cc) i32 {
                    gk.packTriangle(E, upper, across, unit, u, @intCast(m), @intCast(n), @ptrCast(b), @ptrCast(a), @intCast(lda), offset);
                    return 0;
                }
            };
        }

        pub const trsm_iunucopy = Copy(mr, true, false, true).f;
        pub const trsm_iunncopy = Copy(mr, true, false, false).f;
        pub const trsm_iutucopy = Copy(mr, true, true, true).f;
        pub const trsm_iutncopy = Copy(mr, true, true, false).f;
        pub const trsm_ilnucopy = Copy(mr, false, false, true).f;
        pub const trsm_ilnncopy = Copy(mr, false, false, false).f;
        pub const trsm_iltucopy = Copy(mr, false, true, true).f;
        pub const trsm_iltncopy = Copy(mr, false, true, false).f;
        pub const trsm_ounucopy = Copy(nr, true, false, true).f;
        pub const trsm_ounncopy = Copy(nr, true, false, false).f;
        pub const trsm_outucopy = Copy(nr, true, true, true).f;
        pub const trsm_outncopy = Copy(nr, true, true, false).f;
        pub const trsm_olnucopy = Copy(nr, false, false, true).f;
        pub const trsm_olnncopy = Copy(nr, false, false, false).f;
        pub const trsm_oltucopy = Copy(nr, false, true, true).f;
        pub const trsm_oltncopy = Copy(nr, false, true, false).f;

        // ---- kernels -------------------------------------------------------

        fn Kernel(comptime left: bool, comptime forward: bool, comptime conj: bool) type {
            return struct {
                fn real(m: BlasLong, n: BlasLong, k: BlasLong, tag: F, a: [*]F, b: [*]F, c: [*]F, ldc: BlasLong, offset: BlasLong) callconv(cc) i32 {
                    run(left, forward, conj, tag > 0, m, n, k, @ptrCast(a), @ptrCast(b), @ptrCast(c), ldc, offset);
                    return 0;
                }
                fn cplx(m: BlasLong, n: BlasLong, k: BlasLong, tag: F, _: F, a: [*]F, b: [*]F, c: [*]F, ldc: BlasLong, offset: BlasLong) callconv(cc) i32 {
                    run(left, forward, conj, tag > 0, m, n, k, @ptrCast(a), @ptrCast(b), @ptrCast(c), ldc, offset);
                    return 0;
                }
                const f = if (complex) cplx else real;
            };
        }

        pub const trsm_kernel_LT = Kernel(true, true, false).f;
        pub const trsm_kernel_LN = Kernel(true, false, false).f;
        pub const trsm_kernel_RN = Kernel(false, true, false).f;
        pub const trsm_kernel_RT = Kernel(false, false, false).f;
        pub const trsm_kernel_LC = if (complex) Kernel(true, true, true).f else {};
        pub const trsm_kernel_LR = if (complex) Kernel(true, false, true).f else {};
        pub const trsm_kernel_RR = if (complex) Kernel(false, true, true).f else {};
        pub const trsm_kernel_RC = if (complex) Kernel(false, false, true).f else {};

        /// `C -= A * B` on one panel pair, `depth` deep, conjugating the
        /// triangle's side.
        fn update(left: bool, conj: bool, m: usize, n: usize, depth: usize, a: [*]E, b: [*]E, c: [*]E, ldc: usize) void {
            if (depth == 0) return;
            const minus_one: [2]F = .{ -1, 0 };
            gk.sideKernel(G, F, complex, left, conj, m, n, depth, &minus_one, a, b, c, ldc);
        }

        /// The start of panel `p` from the end of `n` lines in panels of `u`:
        /// the tails smallest first, then the whole panels.
        fn backPanel(n: usize, u: usize, from_end: usize) [2]usize {
            const rem = n % u;
            const full = n - rem;
            var seen: usize = 0;
            var w: usize = 1;
            while (w < u) : (w *= 2) {
                if (rem & w != 0) {
                    if (seen == from_end) return .{ full + (rem & ~(2 * w - 1)), w };
                    seen += 1;
                }
            }
            const p = from_end - seen;
            return .{ full - (p + 1) * u, u };
        }

        fn panelCount(n: usize, u: usize) usize {
            return n / u + @popCount(n % u);
        }

        fn run(left: bool, forward: bool, comptime conj: bool, unit: bool, bm: BlasLong, bn: BlasLong, bk: BlasLong, a: [*]E, b: [*]E, c: [*]E, ldc: BlasLong, offset: BlasLong) void {
            const m: usize = @intCast(bm);
            const n: usize = @intCast(bn);
            const k: usize = @intCast(bk);
            const ld: usize = @intCast(ldc);
            if (left) {
                // every column panel solves the rows from their diagonal blocks
                var c0: usize = 0;
                while (c0 < n) {
                    const wn = gk.panelWidth(n - c0, nr);
                    const pb = b + c0 * k;
                    if (forward) {
                        // `trsm_kernel_LT.c`: from the top
                        var kk: isize = offset;
                        var r0: usize = 0;
                        while (r0 < m) {
                            const w = gk.panelWidth(m - r0, mr);
                            const pa = a + r0 * k;
                            const pc = c + r0 + c0 * ld;
                            const kku: usize = @intCast(kk);
                            update(true, conj, w, wn, kku, pa, pb, pc, ld);
                            solveLeft(true, conj, unit, w, wn, pa + kku * w, pb + kku * wn, pc, ld);
                            kk += @intCast(w);
                            r0 += w;
                        }
                    } else {
                        // `trsm_kernel_LN.c`: from the bottom
                        var kk: isize = @as(isize, @intCast(m)) + offset;
                        for (0..panelCount(m, mr)) |p| {
                            const rw = backPanel(m, mr, p);
                            const w = rw[1];
                            const pa = a + rw[0] * k;
                            const pc = c + rw[0] + c0 * ld;
                            const kku: usize = @intCast(kk);
                            update(true, conj, w, wn, k - kku, pa + kku * w, pb + kku * wn, pc, ld);
                            solveLeft(false, conj, unit, w, wn, pa + (kku - w) * w, pb + (kku - w) * wn, pc, ld);
                            kk -= @intCast(w);
                        }
                    }
                    c0 += wn;
                }
            } else {
                // every column panel, in the order of the solve, solves every row panel
                var kk: isize = if (forward) -offset else @as(isize, @intCast(n)) - offset;
                for (0..panelCount(n, nr)) |p| {
                    const cw = if (forward) forwardPanel(n, nr, p) else backPanel(n, nr, p);
                    const wn = cw[1];
                    const pb = b + cw[0] * k;
                    const kku: usize = @intCast(kk);
                    var r0: usize = 0;
                    while (r0 < m) {
                        const w = gk.panelWidth(m - r0, mr);
                        const pa = a + r0 * k;
                        const pc = c + r0 + cw[0] * ld;
                        if (forward) {
                            // `trsm_kernel_RN.c`
                            update(false, conj, w, wn, kku, pa, pb, pc, ld);
                            solveRight(true, conj, unit, w, wn, pa + kku * w, pb + kku * wn, pc, ld);
                        } else {
                            // `trsm_kernel_RT.c`
                            update(false, conj, w, wn, k - kku, pa + kku * w, pb + kku * wn, pc, ld);
                            solveRight(false, conj, unit, w, wn, pa + (kku - wn) * w, pb + (kku - wn) * wn, pc, ld);
                        }
                        r0 += w;
                    }
                    if (forward) kk += @intCast(wn) else kk -= @intCast(wn);
                }
            }
        }

        fn forwardPanel(n: usize, u: usize, p: usize) [2]usize {
            var c0: usize = 0;
            for (0..p) |_| c0 += gk.panelWidth(n - c0, u);
            return .{ c0, gk.panelWidth(n - c0, u) };
        }

        /// `solve` of `trsm_kernel_LT.c` (`forward`) and `_LN.c`: the `m x n`
        /// block of `C` by the diagonal block `a` (`m x m`), into `C` and the
        /// packed `b`; a `unit` diagonal is not divided by.
        pub fn solveLeft(forward: bool, comptime conj: bool, unit: bool, m: usize, n: usize, a: [*]E, b: [*]E, c: [*]E, ldc: usize) void {
            if (comptime vector_solve) if (m == mr) return if (forward) solveLeftVector(true, unit, n, a, b, c, ldc) else solveLeftVector(false, unit, n, a, b, c, ldc);
            solveLeftScalar(forward, conj, unit, m, n, a, b, c, ldc);
        }

        /// `solveLeft` an element at a time.
        pub fn solveLeftScalar(forward: bool, comptime conj: bool, unit: bool, m: usize, n: usize, a: [*]E, b: [*]E, c: [*]E, ldc: usize) void {
            for (0..m) |s| {
                const i = if (forward) s else m - 1 - s;
                const d = conjIf(a[i * m + i], conj);
                for (0..n) |j| {
                    const x = if (unit) c[i + j * ldc] else div(c[i + j * ldc], d);
                    b[i * n + j] = x;
                    c[i + j * ldc] = x;
                    const lo = if (forward) i + 1 else 0;
                    const hi = if (forward) m else i;
                    for (lo..hi) |r| c[r + j * ldc] = sub(c[r + j * ldc], mul(x, conjIf(a[i * m + r], conj)));
                }
            }
        }

        /// `solve` of `trsm_kernel_RN.c` (`forward`) and `_RT.c`: the `m x n`
        /// block of `C` by the diagonal block `b` (`n x n`), into `C` and the
        /// packed `a`; a `unit` diagonal is not divided by.
        pub fn solveRight(forward: bool, comptime conj: bool, unit: bool, m: usize, n: usize, a: [*]E, b: [*]E, c: [*]E, ldc: usize) void {
            if (comptime vector_solve) if (m == mr) return solveRightVector(forward, unit, n, a, b, c, ldc);
            solveRightScalar(forward, conj, unit, m, n, a, b, c, ldc);
        }

        /// `solveRight` an element at a time.
        pub fn solveRightScalar(forward: bool, comptime conj: bool, unit: bool, m: usize, n: usize, a: [*]E, b: [*]E, c: [*]E, ldc: usize) void {
            for (0..n) |s| {
                const i = if (forward) s else n - 1 - s;
                const d = conjIf(b[i * n + i], conj);
                for (0..m) |j| {
                    const x = if (unit) c[j + i * ldc] else div(c[j + i * ldc], d);
                    a[i * m + j] = x;
                    c[j + i * ldc] = x;
                    const lo = if (forward) i + 1 else 0;
                    const hi = if (forward) n else i;
                    for (lo..hi) |q| c[j + q * ldc] = sub(c[j + q * ldc], mul(x, conjIf(b[i * n + q], conj)));
                }
            }
        }

        /// Whether a whole tile solves `mr` lanes at a time: real, on vectors
        /// that keep denormals (`geo.vl`). Each element takes the scalar
        /// solve's operations in its order, so the bits are the same (#720).
        pub const vector_solve = !complex and geo.vl > 1;
        pub const tile_rows = mr;
        const V = @Vector(mr, F);
        const Lane = std.math.IntFittingRange(0, mr);

        /// `solveLeft` on a whole tile (`mr` rows): `C`'s columns in vector
        /// registers, the solved row's element taken by lane and put back by
        /// a select, the rows past it updated under a mask (a scalar store
        /// into a column reloaded as a vector stalled: 0.34 to 0.11 of
        /// OpenBLAS at n = 64).
        pub fn solveLeftVector(comptime forward: bool, unit: bool, n: usize, a: [*]const F, b: [*]F, c: [*]F, ldc: usize) void {
            @setEvalBranchQuota(100_000);
            var cols: [nr]V = undefined;
            inline for (0..nr) |j| if (j < n) {
                cols[j] = c[j * ldc ..][0..mr].*;
            };
            const lanes = std.simd.iota(Lane, mr);
            inline for (0..mr) |s| {
                const i = if (forward) s else mr - 1 - s;
                const acol: V = a[i * mr ..][0..mr].*;
                const at: @Vector(mr, Lane) = @splat(@intCast(i));
                const rest = if (forward) lanes > at else lanes < at;
                const here = lanes == at;
                const d = a[i * mr + i];
                inline for (0..nr) |j| if (j < n) {
                    const col = cols[j];
                    const x = if (unit) col[i] else col[i] / d;
                    b[i * n + j] = x;
                    const xs: V = @splat(x);
                    cols[j] = @select(F, here, xs, @select(F, rest, col - xs * acol, col));
                };
            }
            inline for (0..nr) |j| if (j < n) {
                c[j * ldc ..][0..mr].* = cols[j];
            };
        }

        /// `solveRight` on a whole tile (`mr` rows): each of `C`'s columns
        /// solved as one vector, then taken from the later ones.
        pub fn solveRightVector(forward: bool, unit: bool, n: usize, a: [*]F, b: [*]const F, c: [*]F, ldc: usize) void {
            var cols: [nr]V = undefined;
            inline for (0..nr) |j| if (j < n) {
                cols[j] = c[j * ldc ..][0..mr].*;
            };
            if (forward) {
                inline for (0..nr) |i| if (i < n) solveColumn(i, true, unit, n, a, b, &cols);
            } else inline for (0..nr) |t| {
                const i = nr - 1 - t;
                if (i < n) solveColumn(i, false, unit, n, a, b, &cols);
            }
            inline for (0..nr) |j| if (j < n) {
                c[j * ldc ..][0..mr].* = cols[j];
            };
        }

        /// Column `i` of `solveRightVector`, then its share of the others.
        inline fn solveColumn(comptime i: usize, comptime forward: bool, unit: bool, n: usize, a: [*]F, b: [*]const F, cols: *[nr]V) void {
            const x = if (unit) cols[i] else cols[i] / @as(V, @splat(b[i * n + i]));
            a[i * mr ..][0..mr].* = x;
            cols[i] = x;
            inline for (0..nr) |q| if ((if (forward) q > i else q < i) and q < n) {
                cols[q] = cols[q] - x * @as(V, @splat(b[i * n + q]));
            };
        }
    };
}
