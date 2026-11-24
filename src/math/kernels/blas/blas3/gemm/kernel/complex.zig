//! The complex gemm kernels: `beta` and the register tile over packed panels, under each conjugation, split or contracted.

const types = @import("../../../types.zig");
const util = @import("../../../util.zig");
const kroot = @import("../kernel.zig");
const kpack = @import("pack.zig");

const BlasLong = types.BlasLong;
const Conj = kroot.Conj;
const Geometry = kroot.Geometry;
const madd = kroot.madd;
const panelWidth = kroot.panelWidth;
const tails = kroot.tails;
const pack = kpack.pack;

pub fn ComplexKernels(comptime T: type, comptime geo: Geometry) type {
    const F = types.Float(T);
    const cc = types.kernelCc(F);
    const A = types.Acc(F);
    if (A != F and !geo.split) @compileError("a widened sum needs split parts");
    const vl = geo.vl;
    const mr = geo.mr();
    const nr = geo.nr;
    // one element's parts, as the copies move them
    const P = [2]F;
    return struct {
        pub const unroll_m: i32 = mr;
        pub const unroll_n: i32 = nr;

        /// `C := beta*C` on `m x n`; a zero `beta` writes zeros without
        /// reading `C` (`kernel/generic/zgemm_beta.c`).
        pub fn gemm_beta(m: BlasLong, n: BlasLong, _: BlasLong, br: F, bi: F, _: [*]F, _: BlasLong, _: [*]F, _: BlasLong, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
            const height: usize = @intCast(m);
            const ld: usize = @intCast(ldc);
            for (0..@intCast(n)) |j| {
                const col = c[2 * j * ld ..][0 .. 2 * height];
                if (br == 0 and bi == 0) {
                    @memset(col, 0);
                } else {
                    for (0..height) |i| {
                        const x = col[2 * i];
                        const y = col[2 * i + 1];
                        col[2 * i] = util.scalar(F, br * x - bi * y);
                        col[2 * i + 1] = util.scalar(F, br * y + bi * x);
                    }
                }
            }
            return 0;
        }

        pub fn gemm_incopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            pack(P, mr, false, m, n, @ptrCast(a), lda, @ptrCast(b));
            return 0;
        }
        pub fn gemm_itcopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            pack(P, mr, true, m, n, @ptrCast(a), lda, @ptrCast(b));
            return 0;
        }
        pub fn gemm_oncopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            pack(P, nr, false, m, n, @ptrCast(a), lda, @ptrCast(b));
            return 0;
        }
        pub fn gemm_otcopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            pack(P, nr, true, m, n, @ptrCast(a), lda, @ptrCast(b));
            return 0;
        }

        pub fn gemm_kernel_n(bm: BlasLong, bn: BlasLong, bk: BlasLong, ar: F, ai: F, ba: [*]F, bb: [*]F, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
            return run(.n, bm, bn, bk, ar, ai, ba, bb, c, ldc);
        }
        pub fn gemm_kernel_l(bm: BlasLong, bn: BlasLong, bk: BlasLong, ar: F, ai: F, ba: [*]F, bb: [*]F, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
            return run(.l, bm, bn, bk, ar, ai, ba, bb, c, ldc);
        }
        pub fn gemm_kernel_r(bm: BlasLong, bn: BlasLong, bk: BlasLong, ar: F, ai: F, ba: [*]F, bb: [*]F, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
            return run(.r, bm, bn, bk, ar, ai, ba, bb, c, ldc);
        }
        pub fn gemm_kernel_b(bm: BlasLong, bn: BlasLong, bk: BlasLong, ar: F, ai: F, ba: [*]F, bb: [*]F, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
            return run(.b, bm, bn, bk, ar, ai, ba, bb, c, ldc);
        }

        /// `C += alpha * op(A) * op(B)`, `op` conjugating as `conj` says;
        /// the loops of the real kernel. `conj` comptime: a run-time one
        /// selects between vectors, which LLVM cannot compile on MIPS64r6
        /// (issue #192).
        fn run(comptime conj: Conj, bm: BlasLong, bn: BlasLong, bk: BlasLong, ar: F, ai: F, ba: [*]F, bb: [*]F, c: [*]F, ldc: BlasLong) i32 {
            const m: usize = @intCast(bm);
            const n: usize = @intCast(bn);
            const k: usize = @intCast(bk);
            const ld: usize = @intCast(ldc);
            if (geo.split) {
                // one tile function, shaped by values: the scalar fallback's
                // tiles would otherwise be one per width pair
                var p0: usize = 0;
                while (p0 < n) {
                    const wn = panelWidth(n - p0, nr);
                    var r0: usize = 0;
                    while (r0 < m) {
                        const wm = panelWidth(m - r0, mr);
                        splitTile(conj, wm, wn, k, ar, ai, ba + 2 * r0 * k, bb + 2 * p0 * k, c + 2 * (r0 + p0 * ld), ld);
                        r0 += wm;
                    }
                    p0 += wn;
                }
                return 0;
            }
            var q0: usize = 0;
            while (q0 + nr <= n) : (q0 += nr) rows(conj, nr, m, k, ar, ai, ba, bb + 2 * q0 * k, c + 2 * q0 * ld, ld);
            inline for (comptime tails(nr)) |w| {
                if ((n - q0) & w != 0) {
                    rows(conj, w, m, k, ar, ai, ba, bb + 2 * q0 * k, c + 2 * q0 * ld, ld);
                    q0 += w;
                }
            }
            return 0;
        }

        fn rows(comptime conj: Conj, comptime wn: usize, m: usize, k: usize, ar: F, ai: F, ba: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            var r0: usize = 0;
            while (r0 + mr <= m) : (r0 += mr) @call(.always_inline, tile, .{ conj, mr, wn, k, ar, ai, ba + 2 * r0 * k, pb, c + 2 * r0, ld });
            inline for (comptime tails(mr)) |w| {
                if ((m - r0) & w != 0) {
                    @call(.never_inline, tile, .{ conj, w, wn, k, ar, ai, ba + 2 * r0 * k, pb, c + 2 * r0, ld });
                    r0 += w;
                }
            }
        }

        fn tile(comptime conj: Conj, comptime wm: usize, comptime wn: usize, k: usize, ar: F, ai: F, pa: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            if (geo.split) unreachable;
            if (geo.contracted) contractedTile(conj, wm, wn, k, ar, ai, pa, pb, c, ld) else vectorTile(conj, wm, wn, k, ar, ai, pa, pb, c, ld);
        }

        /// One `wm x wn` tile with one accumulator per vector of A and column,
        /// as `cgemm_kernel_8x2_haswell.c` and `zgemm_kernel_4x2_haswell.c`
        /// (OpenBLAS a69062524, BSD-3, LICENSE): A times B's real parts, then A
        /// with its parts swapped times B's imaginary ones, each by `fmaddsub`
        /// (`fmsubadd` where one operand conjugates), which leaves the parts'
        /// sums with the signs `vectorTile`'s `sb` sets right.
        fn contractedTile(comptime conj: Conj, comptime wm: usize, comptime wn: usize, k: usize, ar: F, ai: F, pa: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            @setEvalBranchQuota(100_000);
            const lanes: usize = 2 * @as(usize, @min(vl, wm));
            const nv = 2 * wm / lanes;
            const V = @Vector(lanes, F);
            const same = comptime conj.ofA() == conj.ofB();
            const swap = comptime swapMask(lanes);
            var acc: [wn][nv]V = undefined;
            inline for (0..wn) |j| inline for (0..nv) |v| {
                acc[j][v] = @splat(0);
            };
            var a = pa;
            var b = pb;
            for (0..k) |_| {
                var av: [nv]V = undefined;
                inline for (0..nv) |v| av[v] = a[v * lanes ..][0..lanes].*;
                inline for (0..wn) |j| {
                    const br: V = @splat(b[2 * j]);
                    inline for (0..nv) |v| acc[j][v] = addSub(same, lanes, av[v], br, acc[j][v]);
                }
                inline for (0..nv) |v| av[v] = @shuffle(F, av[v], undefined, swap);
                inline for (0..wn) |j| {
                    const bi: V = @splat(b[2 * j + 1]);
                    inline for (0..nv) |v| acc[j][v] = addSub(same, lanes, av[v], bi, acc[j][v]);
                }
                a += 2 * wm;
                b += 2 * wn;
            }
            const sb = signs(conj, lanes);
            inline for (0..wn) |j| inline for (0..nv) |v| {
                put(lanes, acc[j][v] * sb, ar, ai, c + 2 * j * ld + v * lanes);
            };
        }

        /// `fmaddsub` where `same` (`a b - c` in the real parts, `a b + c` in
        /// the imaginary ones), else `fmsubadd`: LLVM's x86 backend fuses this
        /// shuffle of two fused multiply-adds into one.
        inline fn addSub(comptime same: bool, comptime lanes: usize, a: @Vector(lanes, F), b: @Vector(lanes, F), c: @Vector(lanes, F)) @Vector(lanes, F) {
            const plus = @mulAdd(@Vector(lanes, F), a, b, c);
            const minus = @mulAdd(@Vector(lanes, F), a, b, -c);
            const mask = comptime blk: {
                var m: [lanes]i32 = undefined;
                for (0..lanes) |i| m[i] = if (i % 2 == 0) @intCast(i) else ~@as(i32, @intCast(i));
                break :blk m;
            };
            return if (same) @shuffle(F, minus, plus, mask) else @shuffle(F, plus, minus, mask);
        }

        /// The signs that make a tile's sums `op(A) op(B)`: with x, y from A
        /// and u, v from B, `vectorTile`'s swapped imaginary sums are (yv, xv),
        /// `contractedTile`'s real parts hold `-(xu - yv)` or `xu + yv`.
        fn signs(comptime conj: Conj, comptime lanes: usize) @Vector(lanes, F) {
            return switch (conj) {
                .n => alternate(lanes, -1, 1),
                .l => alternate(lanes, 1, 1),
                .r => alternate(lanes, 1, -1),
                .b => alternate(lanes, -1, -1),
            };
        }

        /// `C += alpha res` on one vector of elements at `p`.
        inline fn put(comptime lanes: usize, res: @Vector(lanes, F), ar: F, ai: F, p: [*]F) void {
            const V = @Vector(lanes, F);
            const alr: V = @splat(ar);
            const ali = alternate(lanes, -ai, ai);
            const rs = @shuffle(F, res, undefined, swapMask(lanes));
            // whole arrays, not `util.store`: the self-hosted x86_64 backend
            // stores a two-lane `f32` vector as sixteen bytes
            const cv: V = p[0..lanes].*;
            p[0..lanes].* = madd(F, lanes, res, alr, madd(F, lanes, rs, ali, cv));
        }

        /// `[1, 0, 3, 2, ...]`: each element's parts swapped.
        fn swapMask(comptime len: usize) @Vector(len, i32) {
            var s: [len]i32 = undefined;
            for (0..len) |i| s[i] = @intCast(i ^ 1);
            return s;
        }

        /// `x` in the real parts, `y` in the imaginary ones.
        fn alternate(comptime len: usize, x: F, y: F) @Vector(len, F) {
            var s: [len]F = undefined;
            inline for (0..len) |i| s[i] = if (i % 2 == 0) x else y;
            return s;
        }

        /// One `wm x wn` tile on vectors of whole elements: `acc1` sums A
        /// times the real parts of B, `acc2` A times the imaginary ones; the
        /// products' cross terms come from swapping `acc2`'s parts at the end.
        fn vectorTile(comptime conj: Conj, comptime wm: usize, comptime wn: usize, k: usize, ar: F, ai: F, pa: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            @setEvalBranchQuota(100_000);
            const lanes: usize = 2 * @as(usize, @min(vl, wm));
            const nv = 2 * wm / lanes;
            const V = @Vector(lanes, F);
            var acc1: [wn][nv]V = undefined;
            var acc2: [wn][nv]V = undefined;
            inline for (0..wn) |j| inline for (0..nv) |v| {
                acc1[j][v] = @splat(0);
                acc2[j][v] = @splat(0);
            };
            var a = pa;
            var b = pb;
            for (0..k) |_| {
                var av: [nv]V = undefined;
                inline for (0..nv) |v| av[v] = a[v * lanes ..][0..lanes].*;
                inline for (0..wn) |j| {
                    const br: V = @splat(b[2 * j]);
                    const bi: V = @splat(b[2 * j + 1]);
                    inline for (0..nv) |v| {
                        acc1[j][v] = madd(F, lanes, av[v], br, acc1[j][v]);
                        acc2[j][v] = madd(F, lanes, av[v], bi, acc2[j][v]);
                    }
                }
                a += 2 * wm;
                b += 2 * wn;
            }
            const swap = comptime swapMask(lanes);
            // `(x + iy)(u + iv)` with x, y from A and u, v from B: acc1 holds
            // (xu, yu), the swapped acc2 (yv, xv)
            const sb = signs(conj, lanes);
            const sa = comptime alternate(lanes, 1, -1);
            inline for (0..wn) |j| inline for (0..nv) |v| {
                const s = @shuffle(F, acc2[j][v], undefined, swap);
                const first = if (conj.ofA()) acc1[j][v] * sa else acc1[j][v];
                put(lanes, madd(F, lanes, s, sb, first), ar, ai, c + 2 * j * ld + v * lanes);
            };
        }

        /// The same with each part a scalar, summed in `A`.
        fn splitTile(conj: Conj, wm: usize, wn: usize, k: usize, ar: F, ai: F, pa: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            // the tile's most, of which `wm x wn` are used
            var re: [nr][mr]A = @splat(@splat(0));
            var im: [nr][mr]A = @splat(@splat(0));
            var a = pa;
            var b = pb;
            // a conjugated operand's imaginary parts change sign
            const ya: A = if (conj.ofA()) -1 else 1;
            const yb: A = if (conj.ofB()) -1 else 1;
            // loops over the tile, not unrolled: a debug build of the unrolled
            // tile overran its frame on MIPS release 6
            for (0..k) |_| {
                var x: [mr]A = undefined;
                var y: [mr]A = undefined;
                for (0..wm) |i| {
                    x[i] = a[2 * i];
                    y[i] = @as(A, a[2 * i + 1]) * ya;
                }
                for (0..wn) |j| {
                    const u: A = b[2 * j];
                    const v: A = @as(A, b[2 * j + 1]) * yb;
                    for (0..wm) |i| {
                        re[j][i] = madd(A, 1, -y[i], v, madd(A, 1, x[i], u, re[j][i]));
                        im[j][i] = madd(A, 1, x[i], v, madd(A, 1, y[i], u, im[j][i]));
                    }
                }
                a += 2 * wm;
                b += 2 * wn;
            }
            const war: A = ar;
            const wai: A = ai;
            for (0..wn) |j| for (0..wm) |i| {
                const p = c + 2 * (i + j * ld);
                const nre = madd(A, 1, war, re[j][i], madd(A, 1, -wai, im[j][i], p[0]));
                const nim = madd(A, 1, war, im[j][i], madd(A, 1, wai, re[j][i], p[1]));
                p[0] = @floatCast(nre);
                p[1] = @floatCast(nim);
            };
        }
    };
}
