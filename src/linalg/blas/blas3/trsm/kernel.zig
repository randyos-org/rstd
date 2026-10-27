//! The `trsm` kernels and packing copies over the gemm tile.
//!
//! The copies (`kernel/generic/trsm_{u,l}{n,t}copy_*.c`, `ztrsm_*`) pack a
//! block of the triangle from its corner: `m` deep, `n` lines, line `j` at
//! `offset + j` of the depth; element `(j, i)` is `A(j, i)` (`t`) or `A(i,
//! j)` (`n`) where `A` stores it, one on a unit diagonal, and left unwritten
//! elsewhere. Upstream packs the diagonal's reciprocal (`compinv` for
//! complex); here `solve` divides by it, as `trsv` does and
//! Reference-LAPACK master's `dtrsm.f` (`comparisons.md`, 10).
//!
//! The kernels are `kernel/generic/trsm_kernel_{LT,LN,RN,RT}.c` on the gemm
//! layout's panels (whole panels of `unroll_m` or `unroll_n`, then the
//! powers of two): each panel takes the gemm kernel's update by what is
//! solved, then `solve`s its diagonal block, writing the solution to `C` and
//! to the packed right side the later updates read. `LT` and `RN` go
//! forward, `LN` and `RT` backward; `LR`, `LC`, `RR` and `RC` are `LN`, `LT`,
//! `RN` and `RT` conjugating the triangle. Each slot is a thin wrapper: one
//! copy and one kernel body per precision, the options run-time values.

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
    const one: E = if (complex) .{ 1, 0 } else 1;
    return struct {
        // ---- arithmetic on elements -----------------------------------------

        inline fn mul(x: E, y: E) E {
            if (comptime !complex) return util.scalar(E, x * y);
            return util.scalar(E, .{ x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0] });
        }
        fn conjIf(x: E, c: bool) E {
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

        const Source = struct {
            a: [*]const E,
            ld: usize,
            offset: isize,
            across: bool,
            upper: bool,
            unit: bool,

            pub fn at(s: Source, j: usize, i: usize) ?E {
                const jj = s.offset + @as(isize, @intCast(j));
                const ii: isize = @intCast(i);
                const v = s.a[if (s.across) j + i * s.ld else i + j * s.ld];
                if (ii == jj) return if (s.unit) one else v;
                // the element's row and column in the triangle
                const row = if (s.across) jj else ii;
                const col = if (s.across) ii else jj;
                if (if (s.upper) row > col else row < col) return null;
                return v;
            }
        };

        fn copy(u: usize, upper: bool, across: bool, unit: bool, m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, offset: BlasLong, b: [*]F) i32 {
            const src: Source = .{ .a = @ptrCast(a), .ld = @intCast(lda), .offset = offset, .across = across, .upper = upper, .unit = unit };
            gk.packWith(E, u, @intCast(m), @intCast(n), @ptrCast(b), src);
            return 0;
        }

        fn Copy(comptime u: usize, comptime upper: bool, comptime across: bool, comptime unit: bool) type {
            return struct {
                fn f(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, offset: BlasLong, b: [*]F) callconv(cc) i32 {
                    return copy(u, upper, across, unit, m, n, a, lda, offset, b);
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
                fn real(m: BlasLong, n: BlasLong, k: BlasLong, _: F, a: [*]F, b: [*]F, c: [*]F, ldc: BlasLong, offset: BlasLong) callconv(cc) i32 {
                    run(left, forward, conj, m, n, k, @ptrCast(a), @ptrCast(b), @ptrCast(c), ldc, offset);
                    return 0;
                }
                fn cplx(m: BlasLong, n: BlasLong, k: BlasLong, _: F, _: F, a: [*]F, b: [*]F, c: [*]F, ldc: BlasLong, offset: BlasLong) callconv(cc) i32 {
                    run(left, forward, conj, m, n, k, @ptrCast(a), @ptrCast(b), @ptrCast(c), ldc, offset);
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
            const args = .{ @as(BlasLong, @intCast(m)), @as(BlasLong, @intCast(n)), @as(BlasLong, @intCast(depth)) };
            if (comptime complex) {
                const f = if (!conj) &G.gemm_kernel_n else if (left) &G.gemm_kernel_l else &G.gemm_kernel_r;
                _ = f(args[0], args[1], args[2], -1, 0, @ptrCast(a), @ptrCast(b), @ptrCast(c), @intCast(ldc));
            } else {
                _ = G.gemm_kernel(args[0], args[1], args[2], -1, a, b, c, @intCast(ldc));
            }
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

        fn run(left: bool, forward: bool, conj: bool, bm: BlasLong, bn: BlasLong, bk: BlasLong, a: [*]E, b: [*]E, c: [*]E, ldc: BlasLong, offset: BlasLong) void {
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
                            solveLeft(true, conj, w, wn, pa + kku * w, pb + kku * wn, pc, ld);
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
                            solveLeft(false, conj, w, wn, pa + (kku - w) * w, pb + (kku - w) * wn, pc, ld);
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
                            solveRight(true, conj, w, wn, pa + kku * w, pb + kku * wn, pc, ld);
                        } else {
                            // `trsm_kernel_RT.c`
                            update(false, conj, w, wn, k - kku, pa + kku * w, pb + kku * wn, pc, ld);
                            solveRight(false, conj, w, wn, pa + (kku - wn) * w, pb + (kku - wn) * wn, pc, ld);
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
        /// packed `b`.
        fn solveLeft(forward: bool, conj: bool, m: usize, n: usize, a: [*]E, b: [*]E, c: [*]E, ldc: usize) void {
            for (0..m) |s| {
                const i = if (forward) s else m - 1 - s;
                const d = conjIf(a[i * m + i], conj);
                for (0..n) |j| {
                    const x = div(c[i + j * ldc], d);
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
        /// packed `a`.
        fn solveRight(forward: bool, conj: bool, m: usize, n: usize, a: [*]E, b: [*]E, c: [*]E, ldc: usize) void {
            for (0..n) |s| {
                const i = if (forward) s else n - 1 - s;
                const d = conjIf(b[i * n + i], conj);
                for (0..m) |j| {
                    const x = div(c[j + i * ldc], d);
                    a[i * m + j] = x;
                    c[j + i * ldc] = x;
                    const lo = if (forward) i + 1 else 0;
                    const hi = if (forward) n else i;
                    for (lo..hi) |q| c[j + q * ldc] = sub(c[j + q * ldc], mul(x, conjIf(b[i * n + q], conj)));
                }
            }
        }
    };
}
