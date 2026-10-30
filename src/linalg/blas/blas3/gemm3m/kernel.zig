//! The 3M kernels for complex `T` (Higham, "Stability of a Method for
//! Multiplying Complex Matrices with Three Real Matrix Multiplications", SIAM
//! J. Matrix Anal. Appl. 13(3), 681-687, 1992, <https://doi.org/10.1137/0613043>):
//! the copies that pack one part of an operand, or the sum of both, as real
//! panels (`kernel/generic/zgemm3m_ncopy_*.c`, `zgemm3m_tcopy_*.c`,
//! `zsymm3m_ucopy_*.c`, `zsymm3m_lcopy_*.c`, `zhemm3m_ucopy_*.c`,
//! `zhemm3m_lcopy_*.c`), the `o` copies scaling by `alpha` first, and the
//! kernel that adds a real product to both parts of `C`
//! (`kernel/x86_64/zgemm3m_kernel_4x4_haswell.c`). One of the gemm tile's
//! `sets` (`gemm/kernel.zig`), so every core's object builds them.
//!
//! The tile is the real gemm tile of `T`'s parts (`realTile`); the copies
//! pack in its layout, `m` deep and `n` lines.

const std = @import("std");
const types = @import("../../types.zig");
const gk = @import("../gemm/kernel.zig");

const BlasLong = types.BlasLong;
const isComplex = types.isComplex;

/// What a 3M copy packs of each element: both parts summed, or one.
pub const Part = enum { b, r, i };

/// The real tile under complex tile `geo`: the target's real tile for its
/// own complex one (every build but the tests' `installTile`), else `geo`'s
/// lanes as real ones.
pub fn realTile(comptime T: type, comptime geo: gk.Geometry) gk.Geometry {
    if (comptime std.meta.eql(geo, gk.geometry(T))) return gk.geometry(types.Float(T));
    return .{ .vl = if (geo.split) 1 else 2 * geo.vl, .mv = geo.mv, .nr = geo.nr };
}

/// The 3M kernel set of complex `T` under its complex tile `geo`; none for
/// a real `T`. Functions are named as their slots less the precision letter.
pub fn Kernels(comptime T: type, comptime geo: gk.Geometry) type {
    if (comptime !isComplex(T)) return struct {};
    return ThreeM(T, realTile(T, geo));
}

fn ThreeM(comptime T: type, comptime rg: gk.Geometry) type {
    const F = types.Float(T);
    const cc = types.kernelCc(F);
    const A = types.Acc(F);
    const vl = rg.vl;
    if (A != F and vl > 1) @compileError("a widened sum needs scalar lanes");
    const mr = rg.mr();
    const nr = rg.nr;
    const E = [2]F;
    return struct {
        pub const unroll_m: i32 = mr;
        pub const unroll_n: i32 = nr;

        fn Vec(comptime n: usize) type {
            return if (n == 1) F else @Vector(n, F);
        }

        /// `C += alpha * A * B` for real `A` packed `bm x bk` and `B` packed
        /// `bk x bn`, complex `C`: each sum times `alpha_r` into the real
        /// part, times `alpha_i` into the imaginary one. The real kernel's loops
        /// (`gemm/kernel.zig`).
        pub fn gemm3m_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha_r: F, alpha_i: F, ba: [*]F, bb: [*]F, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
            const m: usize = @intCast(bm);
            const n: usize = @intCast(bn);
            const k: usize = @intCast(bk);
            const ld: usize = @intCast(ldc);
            if (vl == 1) {
                scalarPanels(m, n, k, alpha_r, alpha_i, ba, bb, c, ld);
                return 0;
            }
            var q0: usize = 0;
            while (q0 + nr <= n) : (q0 += nr) rows(nr, m, k, alpha_r, alpha_i, ba, bb + q0 * k, c + 2 * q0 * ld, ld);
            inline for (comptime gk.tails(nr)) |w| {
                if ((n - q0) & w != 0) {
                    rows(w, m, k, alpha_r, alpha_i, ba, bb + q0 * k, c + 2 * q0 * ld, ld);
                    q0 += w;
                }
            }
            return 0;
        }

        /// Without vectors, the whole tiles unrolled and the partial ones
        /// shaped by values, as the real kernel does.
        fn scalarPanels(m: usize, n: usize, k: usize, ar: F, ai: F, ba: [*]const F, bb: [*]const F, c: [*]F, ld: usize) void {
            var q0: usize = 0;
            while (q0 < n) {
                const wn = gk.panelWidth(n - q0, nr);
                var r0: usize = 0;
                while (r0 < m) {
                    const wm = gk.panelWidth(m - r0, mr);
                    const p = c + 2 * (r0 + q0 * ld);
                    // a widened sum (`f16`) takes the partial tile throughout
                    if (A == F and wm == mr and wn == nr) tile(mr, nr, k, ar, ai, ba + r0 * k, bb + q0 * k, p, ld) else partTile(wm, wn, k, ar, ai, ba + r0 * k, bb + q0 * k, p, ld);
                    r0 += wm;
                }
                q0 += wn;
            }
        }

        /// `C(i, j) += alpha * s`, one part from each part of `alpha`.
        inline fn put(p: [*]F, s: A, ar: F, ai: F) void {
            const re = gk.madd(A, 1, s, ar, p[0]);
            const im = gk.madd(A, 1, s, ai, p[1]);
            p[0] = if (A == F) re else @floatCast(re);
            p[1] = if (A == F) im else @floatCast(im);
        }

        fn partTile(wm: usize, wn: usize, k: usize, ar: F, ai: F, pa: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            var acc: [nr][mr]A = @splat(@splat(0));
            var a = pa;
            var b = pb;
            for (0..k) |_| {
                for (0..wn) |j| {
                    const y: A = b[j];
                    for (0..wm) |i| acc[j][i] = gk.madd(A, 1, @as(A, a[i]), y, acc[j][i]);
                }
                a += wm;
                b += wn;
            }
            for (0..wn) |j| for (0..wm) |i| put(c + 2 * (i + j * ld), acc[j][i], ar, ai);
        }

        fn rows(comptime wn: usize, m: usize, k: usize, ar: F, ai: F, ba: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            var r0: usize = 0;
            while (r0 + mr <= m) : (r0 += mr) @call(.always_inline, tile, .{ mr, wn, k, ar, ai, ba + r0 * k, pb, c + 2 * r0, ld });
            inline for (comptime gk.tails(mr)) |w| {
                if ((m - r0) & w != 0) {
                    @call(.never_inline, tile, .{ w, wn, k, ar, ai, ba + r0 * k, pb, c + 2 * r0, ld });
                    r0 += w;
                }
            }
        }

        /// `[0, 0, 1, 1, ...]`: each lane twice, for both parts of `C`.
        fn twice(comptime len: usize) @Vector(2 * len, i32) {
            var s: [2 * len]i32 = undefined;
            for (0..2 * len) |i| s[i] = @intCast(i / 2);
            return s;
        }

        /// One `wm x wn` tile: the real kernel's sums, then each lane
        /// doubled against `C`'s parts, as the Haswell kernel's `vunpcklpd`.
        fn tile(comptime wm: usize, comptime wn: usize, k: usize, ar: F, ai: F, pa: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            @setEvalBranchQuota(100_000);
            if (A != F) unreachable;
            const lanes: usize = @min(vl, wm);
            const nv = wm / lanes;
            const V = Vec(lanes);
            var acc: [wn][nv]V = undefined;
            inline for (0..wn) |j| inline for (0..nv) |v| {
                acc[j][v] = if (lanes == 1) 0 else @splat(0);
            };
            var a = pa;
            var b = pb;
            for (0..k) |_| {
                var av: [nv]V = undefined;
                inline for (0..nv) |v| av[v] = if (lanes == 1) a[v] else a[v * lanes ..][0..lanes].*;
                inline for (0..wn) |j| {
                    const bv: V = if (lanes == 1) b[j] else @splat(b[j]);
                    inline for (0..nv) |v| acc[j][v] = gk.madd(F, lanes, av[v], bv, acc[j][v]);
                }
                a += wm;
                b += wn;
            }
            const W = @Vector(2 * lanes, F);
            var al: W = undefined;
            inline for (0..lanes) |i| {
                al[2 * i] = ar;
                al[2 * i + 1] = ai;
            }
            inline for (0..wn) |j| inline for (0..nv) |v| {
                const p = c + 2 * (j * ld + v * lanes);
                if (lanes == 1) {
                    put(p, acc[j][v], ar, ai);
                } else {
                    const d = @shuffle(F, acc[j][v], undefined, twice(lanes));
                    // whole arrays, not `util.store`: the self-hosted x86_64
                    // backend stores a two-lane `f32` vector as sixteen bytes
                    const cv: W = p[0 .. 2 * lanes].*;
                    p[0 .. 2 * lanes].* = gk.madd(F, 2 * lanes, d, al, cv);
                }
            };
        }

        // ---- the copies ----------------------------------------------------

        /// A general operand: line `j`, depth `l` is `a[l * ld + j]`
        /// `across` (`tcopy`), else `a[j * ld + l]` (`ncopy`).
        const General = struct {
            a: [*]const E,
            ld: usize,
            across: bool,

            fn at(s: General, j: usize, l: usize) E {
                return if (s.across) s.a[l * s.ld + j] else s.a[j * s.ld + l];
            }
        };

        /// Element `(x + j, y + l)` of the matrix whose `upper` (or lower)
        /// triangle `a` stores, conjugated for `herm` where it comes from the
        /// other side, with a real diagonal (`symm/kernel.zig`'s source).
        const Sym = struct {
            a: [*]const E,
            ld: usize,
            x: usize,
            y: usize,
            upper: bool,
            herm: bool,

            fn at(s: Sym, j: usize, l: usize) E {
                const x = s.x + j;
                const y = s.y + l;
                const mirrored = if (s.upper) x > y else x < y;
                const v = if (mirrored) s.a[y + x * s.ld] else s.a[x + y * s.ld];
                if (!s.herm) return v;
                if (x == y) return .{ v[0], 0 };
                return if (mirrored) .{ v[0], -v[1] } else v;
            }
        };

        /// `src`'s elements as `cmult` packs them.
        fn Parted(comptime S: type) type {
            return struct {
                s: S,
                part: Part,
                alpha: ?E,
                conj: bool,

                pub fn at(p: @This(), j: usize, l: usize) F {
                    return cmult(p.part, p.alpha, p.conj, p.s.at(j, l));
                }
            };
        }

        /// `CMULT` of the reference's copies: `v`, or `alpha * v` (`conj`:
        /// `alpha * conj(v)`, the `hemm3m` `o` copies' order), then `part`
        /// of it. Products contract as `fma.zig` says a C compiler does.
        fn cmult(part: Part, alpha: ?E, conj: bool, v: E) F {
            const x = v[0];
            const y = v[1];
            const w: E = if (alpha) |al| (if (conj)
                .{ gk.madd(F, 1, al[0], x, al[1] * y), gk.madd(F, 1, al[1], x, -(al[0] * y)) }
            else
                .{ gk.madd(F, 1, al[0], x, -(al[1] * y)), gk.madd(F, 1, al[1], x, al[0] * y) }) else v;
            return switch (part) {
                .b => w[0] + w[1],
                .r => w[0],
                .i => w[1],
            };
        }

        /// One function per kind of copy: the options are values.
        fn general(u: usize, across: bool, part: Part, alpha: ?E, m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) i32 {
            const src: Parted(General) = .{ .s = .{ .a = @ptrCast(a), .ld = @intCast(lda), .across = across }, .part = part, .alpha = alpha, .conj = false };
            gk.packWith(F, u, @intCast(m), @intCast(n), b, src);
            return 0;
        }

        fn sym(u: usize, upper: bool, herm: bool, part: Part, alpha: ?E, m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) i32 {
            const s: Sym = .{ .a = @ptrCast(a), .ld = @intCast(lda), .x = @intCast(posX), .y = @intCast(posY), .upper = upper, .herm = herm };
            const src: Parted(Sym) = .{ .s = s, .part = part, .alpha = alpha, .conj = herm and alpha != null };
            gk.packWith(F, u, @intCast(m), @intCast(n), b, src);
            return 0;
        }

        // `gemm3m`: rows of `op(A)` (`i`, panels of `mr`), columns of `op(B)`
        // times `alpha` (`o`, panels of `nr`); `n` reads down columns, `t` across

        pub fn gemm3m_incopyb(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            return general(mr, false, .b, null, m, n, a, lda, b);
        }
        pub fn gemm3m_incopyr(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            return general(mr, false, .r, null, m, n, a, lda, b);
        }
        pub fn gemm3m_incopyi(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            return general(mr, false, .i, null, m, n, a, lda, b);
        }
        pub fn gemm3m_itcopyb(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            return general(mr, true, .b, null, m, n, a, lda, b);
        }
        pub fn gemm3m_itcopyr(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            return general(mr, true, .r, null, m, n, a, lda, b);
        }
        pub fn gemm3m_itcopyi(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            return general(mr, true, .i, null, m, n, a, lda, b);
        }
        pub fn gemm3m_oncopyb(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return general(nr, false, .b, .{ alpha_r, alpha_i }, m, n, a, lda, b);
        }
        pub fn gemm3m_oncopyr(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return general(nr, false, .r, .{ alpha_r, alpha_i }, m, n, a, lda, b);
        }
        pub fn gemm3m_oncopyi(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return general(nr, false, .i, .{ alpha_r, alpha_i }, m, n, a, lda, b);
        }
        pub fn gemm3m_otcopyb(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return general(nr, true, .b, .{ alpha_r, alpha_i }, m, n, a, lda, b);
        }
        pub fn gemm3m_otcopyr(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return general(nr, true, .r, .{ alpha_r, alpha_i }, m, n, a, lda, b);
        }
        pub fn gemm3m_otcopyi(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return general(nr, true, .i, .{ alpha_r, alpha_i }, m, n, a, lda, b);
        }

        // `symm3m` and `hemm3m`: line `j`, depth `l` is element `(posX + j,
        // posY + l)` of the whole matrix, from its upper (`u`) or lower (`l`)
        // triangle

        pub fn symm3m_iucopyb(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return sym(mr, true, false, .b, null, m, n, a, lda, posX, posY, b);
        }
        pub fn symm3m_ilcopyb(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return sym(mr, false, false, .b, null, m, n, a, lda, posX, posY, b);
        }
        pub fn symm3m_iucopyr(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return sym(mr, true, false, .r, null, m, n, a, lda, posX, posY, b);
        }
        pub fn symm3m_ilcopyr(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return sym(mr, false, false, .r, null, m, n, a, lda, posX, posY, b);
        }
        pub fn symm3m_iucopyi(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return sym(mr, true, false, .i, null, m, n, a, lda, posX, posY, b);
        }
        pub fn symm3m_ilcopyi(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return sym(mr, false, false, .i, null, m, n, a, lda, posX, posY, b);
        }
        pub fn symm3m_oucopyb(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return sym(nr, true, false, .b, .{ alpha_r, alpha_i }, m, n, a, lda, posX, posY, b);
        }
        pub fn symm3m_olcopyb(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return sym(nr, false, false, .b, .{ alpha_r, alpha_i }, m, n, a, lda, posX, posY, b);
        }
        pub fn symm3m_oucopyr(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return sym(nr, true, false, .r, .{ alpha_r, alpha_i }, m, n, a, lda, posX, posY, b);
        }
        pub fn symm3m_olcopyr(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return sym(nr, false, false, .r, .{ alpha_r, alpha_i }, m, n, a, lda, posX, posY, b);
        }
        pub fn symm3m_oucopyi(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return sym(nr, true, false, .i, .{ alpha_r, alpha_i }, m, n, a, lda, posX, posY, b);
        }
        pub fn symm3m_olcopyi(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return sym(nr, false, false, .i, .{ alpha_r, alpha_i }, m, n, a, lda, posX, posY, b);
        }

        pub fn hemm3m_iucopyb(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return sym(mr, true, true, .b, null, m, n, a, lda, posX, posY, b);
        }
        pub fn hemm3m_ilcopyb(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return sym(mr, false, true, .b, null, m, n, a, lda, posX, posY, b);
        }
        pub fn hemm3m_iucopyr(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return sym(mr, true, true, .r, null, m, n, a, lda, posX, posY, b);
        }
        pub fn hemm3m_ilcopyr(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return sym(mr, false, true, .r, null, m, n, a, lda, posX, posY, b);
        }
        pub fn hemm3m_iucopyi(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return sym(mr, true, true, .i, null, m, n, a, lda, posX, posY, b);
        }
        pub fn hemm3m_ilcopyi(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return sym(mr, false, true, .i, null, m, n, a, lda, posX, posY, b);
        }
        pub fn hemm3m_oucopyb(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return sym(nr, true, true, .b, .{ alpha_r, alpha_i }, m, n, a, lda, posX, posY, b);
        }
        pub fn hemm3m_olcopyb(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return sym(nr, false, true, .b, .{ alpha_r, alpha_i }, m, n, a, lda, posX, posY, b);
        }
        pub fn hemm3m_oucopyr(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return sym(nr, true, true, .r, .{ alpha_r, alpha_i }, m, n, a, lda, posX, posY, b);
        }
        pub fn hemm3m_olcopyr(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return sym(nr, false, true, .r, .{ alpha_r, alpha_i }, m, n, a, lda, posX, posY, b);
        }
        pub fn hemm3m_oucopyi(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return sym(nr, true, true, .i, .{ alpha_r, alpha_i }, m, n, a, lda, posX, posY, b);
        }
        pub fn hemm3m_olcopyi(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: F, alpha_i: F, b: [*]F) callconv(cc) i32 {
            return sym(nr, false, true, .i, .{ alpha_r, alpha_i }, m, n, a, lda, posX, posY, b);
        }
    };
}

// ---- tests -----------------------------------------------------------------

fn fillSeq(comptime F: type, buf: []F, seed: usize) void {
    for (buf, 0..) |*v, i| v.* = @floatFromInt(@as(i64, @intCast((i * 7 + seed * 13) % 23)) - 11);
}

/// `C += alpha * A * B` through the three passes of `gemm3m_level3.c` (`NN`)
/// against plain complex arithmetic, for every panel remainder, `A` and `B`
/// packed both ways; exact on small integers.
fn checkProduct(comptime T: type, comptime geo: gk.Geometry) !void {
    const K = Kernels(T, geo);
    const F = types.Float(T);
    const mr: usize = K.unroll_m;
    const nr: usize = K.unroll_n;
    const al = std.testing.allocator;
    const sizes = [_]usize{ 1, 2, 3, 5, mr - 1, mr, mr + 1, 2 * mr + 3, nr + 1, 2 * nr + 1 };
    const ld = std.mem.max(usize, &sizes);
    const kmax = 13;
    const buf = try al.alloc(F, 2 * (4 * ld * kmax + 3 * ld * ld) + 2 * ld * kmax);
    defer al.free(buf);
    const a = buf[0 .. 2 * ld * kmax];
    const at = buf[2 * ld * kmax ..][0 .. 2 * ld * kmax];
    const b = buf[4 * ld * kmax ..][0 .. 2 * ld * kmax];
    const bt = buf[6 * ld * kmax ..][0 .. 2 * ld * kmax];
    const c = buf[8 * ld * kmax ..][0 .. 2 * ld * ld];
    const want = buf[8 * ld * kmax + 2 * ld * ld ..][0 .. 2 * ld * ld];
    const got = buf[8 * ld * kmax + 4 * ld * ld ..][0 .. 2 * ld * ld];
    const pa = buf[8 * ld * kmax + 6 * ld * ld ..][0 .. ld * kmax];
    const pb = buf[9 * ld * kmax + 6 * ld * ld ..][0 .. ld * kmax];
    fillSeq(F, a, 1);
    fillSeq(F, b, 2);
    fillSeq(F, c, 3);
    for (0..ld) |i| for (0..kmax) |l| for (0..2) |p| {
        at[2 * (l + i * kmax) + p] = a[2 * (i + l * ld) + p];
        bt[2 * (i + l * ld) + p] = b[2 * (l + i * kmax) + p];
    };
    const ar: F = 2;
    const ai: F = -1;
    for (sizes) |m| for (sizes) |n| for ([_]usize{ 1, 4, kmax }) |k| {
        @memcpy(want, c);
        for (0..n) |j| for (0..m) |i| {
            var sr: F = 0;
            var si: F = 0;
            for (0..k) |l| {
                const x = a[2 * (i + l * ld)];
                const y = a[2 * (i + l * ld) + 1];
                const u = b[2 * (l + j * kmax)];
                const v = b[2 * (l + j * kmax) + 1];
                sr += x * u - y * v;
                si += x * v + y * u;
            }
            want[2 * (i + j * ld)] += ar * sr - ai * si;
            want[2 * (i + j * ld) + 1] += ar * si + ai * sr;
        };
        inline for (.{ false, true }) |trans| {
            @memcpy(got, c);
            const lm: BlasLong = @intCast(ld);
            const passes = .{ .{ "b", 0, 1 }, .{ "r", 1, -1 }, .{ "i", -1, -1 } };
            inline for (passes) |p| {
                if (trans) {
                    _ = @field(K, "gemm3m_incopy" ++ p[0])(@intCast(k), @intCast(m), at.ptr, kmax, pa.ptr);
                    _ = @field(K, "gemm3m_otcopy" ++ p[0])(@intCast(k), @intCast(n), bt.ptr, lm, ar, ai, pb.ptr);
                } else {
                    _ = @field(K, "gemm3m_itcopy" ++ p[0])(@intCast(k), @intCast(m), a.ptr, lm, pa.ptr);
                    _ = @field(K, "gemm3m_oncopy" ++ p[0])(@intCast(k), @intCast(n), b.ptr, kmax, ar, ai, pb.ptr);
                }
                _ = K.gemm3m_kernel(@intCast(m), @intCast(n), @intCast(k), p[1], p[2], pa.ptr, pb.ptr, got.ptr, lm);
            }
            if (!std.mem.eql(F, want, got)) return error.TestUnexpectedResult;
        }
    };
}

test "3m: the three passes' copies and kernel give the complex product, every panel remainder" {
    try checkProduct(types.Complex64, gk.geometry(types.Complex64));
    try checkProduct(types.Complex32, gk.geometry(types.Complex32));
    try checkProduct(types.Complex64, .{ .vl = 1, .mv = 3, .nr = 5, .split = true });
    try checkProduct(types.Complex32, .{ .vl = 2, .mv = 2, .nr = 3 });
}

test "3m: the real tile is the target's, or the given complex tile's lanes" {
    const F = f64;
    const T = types.Complex64;
    try std.testing.expectEqual(gk.geometry(F), realTile(T, gk.geometry(T)));
    try std.testing.expectEqual(gk.Geometry{ .vl = 4, .mv = 2, .nr = 3 }, realTile(T, .{ .vl = 2, .mv = 2, .nr = 3 }));
    try std.testing.expectEqual(gk.Geometry{ .vl = 1, .mv = 3, .nr = 5 }, realTile(T, .{ .vl = 1, .mv = 3, .nr = 5, .split = true }));
}

test "3m: symmetric and Hermitian copies pack the whole matrix's element, scaled for o" {
    const T = types.Complex64;
    const K = Kernels(T, gk.geometry(T));
    const n = 5;
    var a: [2 * n * n]f64 = undefined;
    fillSeq(f64, &a, 4);
    var out: [n * n]f64 = undefined;
    const ar: f64 = 2;
    const ai: f64 = -3;
    inline for (.{ false, true }) |herm| inline for (.{ true, false }) |upper| inline for (.{ "b", "r", "i" }) |part| inline for (.{ "i", "o" }) |io| {
        const name = (if (herm) "hemm3m_" else "symm3m_") ++ io ++ (if (upper) "u" else "l") ++ "copy" ++ part;
        const f = @field(K, name);
        // the block from line 1, depth 2: 3 lines of depth 3
        if (comptime io[0] == 'i') _ = f(3, 3, &a, n, 1, 2, &out) else _ = f(3, 3, &a, n, 1, 2, ar, ai, &out);
        // packed in panels of the unroll, which `packWith` lays out as the
        // gemm copies do: check each element by recomputing its offset
        const u: usize = if (io[0] == 'i') K.unroll_m else K.unroll_n;
        var c0: usize = 0;
        while (c0 < 3) {
            const w = gk.panelWidth(3 - c0, u);
            for (0..3) |l| for (0..w) |cj| {
                const x = 1 + c0 + cj;
                const y = 2 + l;
                const mirrored = if (upper) x > y else x < y;
                const at = if (mirrored) y + x * n else x + y * n;
                var re = a[2 * at];
                var im = a[2 * at + 1];
                if (herm) {
                    if (x == y) im = 0 else if (mirrored) im = -im;
                }
                if (io[0] == 'o') {
                    // `o` packs `alpha` times the element of `op(B)`: for
                    // `hemm` the conjugate of the row copy's
                    if (herm) im = -im;
                    const r2 = ar * re - ai * im;
                    im = ai * re + ar * im;
                    re = r2;
                }
                const want = if (part[0] == 'b') re + im else if (part[0] == 'r') re else im;
                try std.testing.expectEqual(want, out[3 * c0 + l * w + cj]);
            };
            c0 += w;
        }
    };
}
