//! The `trmm` copies (`kernel/generic/trmm_{u,l}{n,t}copy_*.c`, `ztrmm_*`) and
//! kernels (`trmmkernel_*.c`) over the gemm tile; the diagonal block sums its
//! triangle alone, so an infinity meets no packed zero, unlike upstream.

const std = @import("std");
const types = @import("../../types.zig");
const util = @import("../../util.zig");
const gk = @import("../gemm/kernel.zig");

const BlasLong = types.BlasLong;

pub fn Kernels(comptime T: type, comptime geo: gk.Geometry) type {
    const complex = types.isComplex(T);
    const F = types.Float(T);
    const cc = types.kernelCc(F);
    const A = types.Acc(F);
    // one element as the kernels move it
    const E = if (complex) [2]F else F;
    const G = gk.Kernels(T, geo);
    const mr = geo.mr();
    const nr = geo.nr;
    return struct {
        // ---- copies --------------------------------------------------------

        /// The copy `trmm_<io><ul><nt><un>copy`.
        fn Copy(comptime u: usize, comptime upper: bool, comptime across: bool, comptime unit: bool) type {
            return struct {
                fn f(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
                    const ld: usize = @intCast(lda);
                    const x: usize = @intCast(posX);
                    const y: usize = @intCast(posY);
                    const base = @as([*]const E, @ptrCast(a)) + if (across) y + x * ld else x + y * ld;
                    gk.packTriangle(E, upper, across, unit, u, @intCast(m), @intCast(n), @ptrCast(b), base, ld, @as(isize, @intCast(y)) - @as(isize, @intCast(x)));
                    return 0;
                }
            };
        }

        pub const trmm_iunucopy = Copy(mr, true, false, true).f;
        pub const trmm_iunncopy = Copy(mr, true, false, false).f;
        pub const trmm_iutucopy = Copy(mr, true, true, true).f;
        pub const trmm_iutncopy = Copy(mr, true, true, false).f;
        pub const trmm_ilnucopy = Copy(mr, false, false, true).f;
        pub const trmm_ilnncopy = Copy(mr, false, false, false).f;
        pub const trmm_iltucopy = Copy(mr, false, true, true).f;
        pub const trmm_iltncopy = Copy(mr, false, true, false).f;
        pub const trmm_ounucopy = Copy(nr, true, false, true).f;
        pub const trmm_ounncopy = Copy(nr, true, false, false).f;
        pub const trmm_outucopy = Copy(nr, true, true, true).f;
        pub const trmm_outncopy = Copy(nr, true, true, false).f;
        pub const trmm_olnucopy = Copy(nr, false, false, true).f;
        pub const trmm_olnncopy = Copy(nr, false, false, false).f;
        pub const trmm_oltucopy = Copy(nr, false, true, true).f;
        pub const trmm_oltncopy = Copy(nr, false, true, false).f;

        // ---- kernels -------------------------------------------------------

        /// The kernel of a side; `from`: depths from the row (or column) on
        /// (`LN`, `RT`), else up to it; `conj` conjugates the triangle.
        fn Kernel(comptime left: bool, comptime from: bool, comptime conj: bool) type {
            return struct {
                fn real(m: BlasLong, n: BlasLong, k: BlasLong, alpha: F, a: [*]F, b: [*]F, c: [*]F, ldc: BlasLong, offset: BlasLong) callconv(cc) i32 {
                    run(left, from, conj, m, n, k, @ptrCast(&alpha), a, b, c, ldc, offset);
                    return 0;
                }
                fn cplx(m: BlasLong, n: BlasLong, k: BlasLong, ar: F, ai: F, a: [*]F, b: [*]F, c: [*]F, ldc: BlasLong, offset: BlasLong) callconv(cc) i32 {
                    const alpha: [2]F = .{ ar, ai };
                    run(left, from, conj, m, n, k, &alpha, a, b, c, ldc, offset);
                    return 0;
                }
                const f = if (complex) cplx else real;
            };
        }

        pub const trmm_kernel_LN = Kernel(true, true, false).f;
        pub const trmm_kernel_LT = Kernel(true, false, false).f;
        pub const trmm_kernel_RN = Kernel(false, false, false).f;
        pub const trmm_kernel_RT = Kernel(false, true, false).f;
        pub const trmm_kernel_LR = if (complex) Kernel(true, true, true).f else {};
        pub const trmm_kernel_LC = if (complex) Kernel(true, false, true).f else {};
        pub const trmm_kernel_RR = if (complex) Kernel(false, false, true).f else {};
        pub const trmm_kernel_RC = if (complex) Kernel(false, true, true).f else {};

        fn run(left: bool, from: bool, conj: bool, bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: [*]const F, ba: [*]F, bb: [*]F, bc: [*]F, ldc: BlasLong, offset: BlasLong) void {
            const m: usize = @intCast(bm);
            const n: usize = @intCast(bn);
            const k: isize = bk;
            const ku: usize = @intCast(k);
            const ld: usize = @intCast(ldc);
            const a: [*]E = @ptrCast(ba);
            const b: [*]E = @ptrCast(bb);
            const c: [*]E = @ptrCast(bc);
            const zero: E = if (complex) .{ 0, 0 } else 0;
            var c0: usize = 0;
            while (c0 < n) {
                const wn = gk.panelWidth(n - c0, nr);
                var r0: usize = 0;
                while (r0 < m) {
                    const w = gk.panelWidth(m - r0, mr);
                    const pa = a + r0 * ku;
                    const pb = b + c0 * ku;
                    const pc = c + r0 + c0 * ld;
                    // the diagonal block's depths, and the rest the triangle takes
                    const off: isize = if (left) offset + @as(isize, @intCast(r0)) else @as(isize, @intCast(c0)) - offset;
                    const width: isize = @intCast(if (left) w else wn);
                    const d0 = clip(off, k);
                    const d1 = clip(off + width, k);
                    const rest: [2]usize = if (from) .{ d1, ku } else .{ 0, d0 };
                    for (0..wn) |j| @memset(pc[j * ld ..][0..w], zero);
                    if (rest[1] > rest[0]) gk.sideKernel(G, F, complex, left, conj, w, wn, rest[1] - rest[0], alpha, pa + rest[0] * w, pb + rest[0] * wn, pc, ld);
                    if (vector_diag and w == mr) {
                        diagonal(left, from, wn, d0, d1, off, @ptrCast(pa), @ptrCast(pb), @ptrCast(pc), ld, alpha[0]);
                    } else diagonalScalar(left, from, conj, w, wn, d0, d1, off, pa, pb, pc, ld, alpha);
                    r0 += w;
                }
                c0 += wn;
            }
        }

        /// The diagonal block of a `w x wn` tile at depths `d0` to `d1`, its
        /// triangle alone, into `pc`: `alpha` times the sum of each element's
        /// products, in depth order.
        pub fn diagonalScalar(left: bool, from: bool, conj: bool, w: usize, wn: usize, d0: usize, d1: usize, off: isize, pa: [*]const E, pb: [*]const E, pc: [*]E, ld: usize, alpha: [*]const F) void {
            for (0..wn) |j| for (0..w) |i| {
                // the row (left) or column of `op(A)` this element is
                const at: isize = off + @as(isize, @intCast(if (left) i else j));
                var sr: A = 0;
                var si: A = 0;
                for (d0..d1) |d| {
                    const di: isize = @intCast(d);
                    // left: `D >= L` from, `D <= L` up to; right likewise
                    if (if (from) di < at else di > at) continue;
                    const x = pa[d * w + i];
                    const y = pb[d * wn + j];
                    if (comptime complex) {
                        const xi: A = if (conj and left) -@as(A, x[1]) else x[1];
                        const yi: A = if (conj and !left) -@as(A, y[1]) else y[1];
                        sr += @as(A, x[0]) * y[0] - xi * yi;
                        si += @as(A, x[0]) * yi + xi * y[0];
                    } else sr = gk.madd(A, 1, x, y, sr);
                }
                const p = &pc[i + j * ld];
                if (comptime complex) {
                    const ar: A = alpha[0];
                    const ai: A = alpha[1];
                    p.* = .{ util.scalar(F, @floatCast(p[0] + (ar * sr - ai * si))), util.scalar(F, @floatCast(p[1] + (ar * si + ai * sr))) };
                } else p.* = util.scalar(E, @floatCast(p.* + @as(A, alpha[0]) * sr));
            };
        }

        /// Whether a whole tile's diagonal block runs `mr` lanes at a time:
        /// real, on vectors that keep denormals (`geo.vl`), summed in `F`.
        pub const vector_diag = !complex and geo.vl > 1 and A == F;
        pub const tile_rows = mr;

        /// The diagonal block of a whole tile, `mr` rows by `wn` columns, as
        /// `run`'s scalar loop sums it: a lane outside the triangle keeps its
        /// sum, so an infinity meets no zero and the bits are its own, a NaN's
        /// sign and payload aside (LLVM keeps neither). Each element sums its
        /// depths in order. Left, the block's depths go in segments of `vl`:
        /// in segment `q`, half `q` of the tile's rows is cut by the triangle,
        /// the halves on its kept side are whole and the others are skipped
        /// (#719: the masked product of every lane was 62% of `f32` trmm's
        /// time at n = 64 on Zen 3).
        pub fn diagonal(left: bool, from: bool, wn: usize, d0: usize, d1: usize, off: isize, pa: [*]const F, pb: [*]const F, pc: [*]F, ld: usize, alpha: F) void {
            var s: [nr][mv]H = @splat(@splat(@splat(0)));
            if (left) {
                if (from) leftBlock(true, wn, d0, d1, off, pa, pb, &s) else leftBlock(false, wn, d0, d1, off, pa, pb, &s);
            } else if (from) rightBlock(true, wn, d0, d1, off, pa, pb, &s) else rightBlock(false, wn, d0, d1, off, pa, pb, &s);
            inline for (0..nr) |j| if (j < wn) inline for (0..mv) |h| {
                const c = pc + j * ld + h * vl;
                c[0..vl].* = @as(H, c[0..vl].*) + @as(H, @splat(alpha)) * s[j][h];
            };
        }

        const vl = geo.vl;
        const mv = geo.mv;
        // row numbers as lanes as wide as `F`'s: one compare per vector
        const Ix = if (@bitSizeOf(F) == 64) i64 else i32;
        const H = @Vector(vl, F);
        const IH = @Vector(vl, Ix);

        /// Left: row `i` of the tile is row `off + i` of `op(A)`, kept at depth
        /// `d` when `d >= off + i` (`from`), else when `d <= off + i`.
        fn leftBlock(comptime from: bool, wn: usize, d0: usize, d1: usize, off: isize, pa: [*]const F, pb: [*]const F, s: *[nr][mv]H) void {
            inline for (0..mv) |q| {
                const base = off + @as(isize, @intCast(q * vl));
                const lo = clamp(base, d0, d1);
                const hi = clamp(base + vl, d0, d1);
                const rows = @as(IH, @splat(@as(Ix, @intCast(base)))) + std.simd.iota(Ix, vl);
                for (lo..hi) |d| {
                    const at: IH = @splat(@as(Ix, @intCast(d)));
                    const keep = if (from) at >= rows else at <= rows;
                    inline for (0..mv) |h| if (h == q or (if (from) h < q else h > q)) {
                        const x: H = pa[d * mr + h * vl ..][0..vl].*;
                        inline for (0..nr) |j| if (j < wn) {
                            const f = madd(x, @splat(pb[d * wn + j]), s[j][h]);
                            s[j][h] = if (h == q) @select(F, keep, f, s[j][h]) else f;
                        };
                    };
                }
            }
        }

        /// Right: column `j` of the tile is column `off + j`, kept at depth `d`
        /// when `d >= off + j` (`from`), else when `d <= off + j`.
        fn rightBlock(comptime from: bool, wn: usize, d0: usize, d1: usize, off: isize, pa: [*]const F, pb: [*]const F, s: *[nr][mv]H) void {
            for (d0..d1) |d| {
                const di: isize = @intCast(d);
                inline for (0..mv) |h| {
                    const x: H = pa[d * mr + h * vl ..][0..vl].*;
                    inline for (0..nr) |j| if (j < wn and (if (from) di >= off + j else di <= off + j)) {
                        s[j][h] = madd(x, @splat(pb[d * wn + j]), s[j][h]);
                    };
                }
            }
        }

        /// `gk.madd` once per type: inline, its branch count charges every
        /// unrolled call.
        fn madd(x: H, y: H, c: H) H {
            return gk.madd(F, vl, x, y, c);
        }

        fn clamp(v: isize, lo: usize, hi: usize) usize {
            return @intCast(std.math.clamp(v, @as(isize, @intCast(lo)), @as(isize, @intCast(hi))));
        }

        fn clip(v: isize, k: isize) usize {
            return @intCast(@max(0, @min(v, k)));
        }
    };
}
