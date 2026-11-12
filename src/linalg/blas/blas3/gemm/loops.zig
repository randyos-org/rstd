//! The real gemm kernel's loops, `C += alpha * A * B` on real `A` packed in
//! panels of `mr` and `B` in panels of `nr`, shared by the real kernel and
//! 3M's (`../gemm3m/kernel.zig`), which differ only in how a sum meets `C`.

const types = @import("../../types.zig");
const gk = @import("kernel.zig");

const Geometry = gk.Geometry;
const Lanes = gk.Lanes;
const madd = gk.madd;
const splat = gk.splat;
const tails = gk.tails;
const panelWidth = gk.panelWidth;

/// The loops over real `T` with tile `geo`, ending in `Out`'s store:
/// - `C`, the type of `C`'s elements, `scale` of them per tile element;
/// - `Alpha`, and `prep(lanes, alpha)`, `alpha` as one tile's stores take it;
/// - `put(p, sum, alpha)`, one scalar sum into `C` at `p`;
/// - `putLanes(lanes, p, sums, prepped)`, `lanes` sums into `C` at `p`.
pub fn Loops(comptime T: type, comptime geo: Geometry, comptime Out: type) type {
    const vl = geo.vl;
    // `f16` sums in `f32`; vectors are only ever `f32` or `f64`, summed as is
    const A = types.Acc(T);
    if (A != T and vl > 1) @compileError("a widened sum needs scalar lanes");
    const mr = geo.mr();
    const nr = geo.nr;
    const C = Out.C;
    const s: usize = Out.scale;
    return struct {
        /// Columns outermost, so a B panel stays in L1 while the A panels
        /// stream from L2.
        pub inline fn run(m: usize, n: usize, k: usize, alpha: Out.Alpha, ba: [*]const T, bb: [*]const T, c: [*]C, ld: usize) void {
            if (vl == 1) return scalarPanels(m, n, k, alpha, ba, bb, c, ld);
            var q0: usize = 0;
            while (q0 + nr <= n) : (q0 += nr) rows(nr, m, k, alpha, ba, bb + q0 * k, c + s * q0 * ld, ld);
            inline for (comptime tails(nr)) |w| {
                if ((n - q0) & w != 0) {
                    rows(w, m, k, alpha, ba, bb + q0 * k, c + s * q0 * ld, ld);
                    q0 += w;
                }
            }
        }

        /// Without vectors: the whole tiles unrolled, every partial one through
        /// one tile shaped by values, so debug builds of wide scalar tiles
        /// (SPARC's `2 x 14`) stay small.
        fn scalarPanels(m: usize, n: usize, k: usize, alpha: Out.Alpha, ba: [*]const T, bb: [*]const T, c: [*]C, ld: usize) void {
            var q0: usize = 0;
            while (q0 < n) {
                const wn = panelWidth(n - q0, nr);
                var r0: usize = 0;
                while (r0 < m) {
                    const wm = panelWidth(m - r0, mr);
                    const p = c + s * (r0 + q0 * ld);
                    if (wm == mr and wn == nr) {
                        tile(mr, nr, k, alpha, ba + r0 * k, bb + q0 * k, p, ld);
                    } else {
                        partTile(wm, wn, k, alpha, ba + r0 * k, bb + q0 * k, p, ld);
                    }
                    r0 += wm;
                }
                q0 += wn;
            }
        }

        /// A partial scalar tile, `wm x wn` of the tile's most.
        fn partTile(wm: usize, wn: usize, k: usize, alpha: Out.Alpha, pa: [*]const T, pb: [*]const T, c: [*]C, ld: usize) void {
            var acc: [nr][mr]A = @splat(@splat(0));
            var a = pa;
            var b = pb;
            for (0..k) |_| {
                for (0..wn) |j| {
                    const y: A = b[j];
                    for (0..wm) |i| acc[j][i] = madd(A, 1, @as(A, a[i]), y, acc[j][i]);
                }
                a += wm;
                b += wn;
            }
            for (0..wn) |j| for (0..wm) |i| Out.put(c + s * (i + j * ld), acc[j][i], alpha);
        }

        /// Every row panel against one B panel of width `wn`: the full tile
        /// inline, each row tail an unrolled tile of its own, as OpenBLAS's
        /// edge tiles are.
        fn rows(comptime wn: usize, m: usize, k: usize, alpha: Out.Alpha, ba: [*]const T, pb: [*]const T, c: [*]C, ld: usize) void {
            var r0: usize = 0;
            while (r0 + mr <= m) : (r0 += mr) @call(.always_inline, tile, .{ mr, wn, k, alpha, ba + r0 * k, pb, c + s * r0, ld });
            inline for (comptime tails(mr)) |w| {
                if ((m - r0) & w != 0) {
                    @call(.never_inline, tile, .{ w, wn, k, alpha, ba + r0 * k, pb, c + s * r0, ld });
                    r0 += w;
                }
            }
        }

        /// One `wm x wn` tile: `wm / lanes` vectors per column in registers
        /// over all of `k`, then each into `C` through `Out`.
        fn tile(comptime wm: usize, comptime wn: usize, k: usize, alpha: Out.Alpha, pa: [*]const T, pb: [*]const T, c: [*]C, ld: usize) void {
            @setEvalBranchQuota(100_000);
            const lanes = @min(vl, wm);
            const nv = wm / lanes;
            const V = Lanes(A, lanes);
            var acc: [wn][nv]V = undefined;
            inline for (0..wn) |j| inline for (0..nv) |v| {
                acc[j][v] = splat(A, lanes, 0);
            };
            var a = pa;
            var b = pb;
            for (0..k) |_| {
                var av: [nv]V = undefined;
                inline for (0..nv) |v| av[v] = wide(lanes, a + v * lanes);
                inline for (0..wn) |j| {
                    const bv = splat(A, lanes, b[j]);
                    inline for (0..nv) |v| acc[j][v] = madd(A, lanes, av[v], bv, acc[j][v]);
                }
                a += wm;
                b += wn;
            }
            const al = Out.prep(lanes, alpha);
            inline for (0..wn) |j| inline for (0..nv) |v| {
                Out.putLanes(lanes, c + s * (j * ld + v * lanes), acc[j][v], al);
            };
        }

        /// `lanes` elements at `p` as `A`.
        pub inline fn wide(comptime lanes: usize, p: [*]const T) Lanes(A, lanes) {
            return if (A == T) gk.load(T, lanes, p) else p[0];
        }
    };
}
