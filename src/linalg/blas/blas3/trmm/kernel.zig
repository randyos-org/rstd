//! The `trmm` kernels and packing copies over the gemm tile.
//!
//! The copies (`kernel/generic/trmm_{u,l}{n,t}copy_*.c`, `ztrmm_*`) pack a
//! block of `op(A)`: `m` deep from `posX`, `n` lines from `posY`; line `L`,
//! depth `D` is `A(L, D)` (`t`) or `A(D, L)` (`n`) where `A` stores it, one
//! on a unit diagonal, and left unwritten elsewhere.
//!
//! The kernels follow `trmmkernel_*.c` (`TRMMKERNEL`): `C := alpha * A * B`
//! on packed panels, only over the depths the triangle has. On the left, the
//! row panel at `r` is row `offset + r` of `op(A)`, and takes depths from it
//! (`LN`) or up to it (`LT`); on the right the column panel at `c` is column
//! `c - offset`, and takes depths up to it (`RN`) or from it (`RT`). `R` and
//! `C` conjugate the triangle. Unlike the reference, whose tiles multiply the
//! zeros a diagonal block packs, the diagonal block sums its triangle alone,
//! so an infinity meets no zero.
//!
//! Each slot is a thin wrapper: one copy and one kernel body per precision,
//! the options run-time values, keep the code small.

const types = @import("../../types.zig");
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

        const Source = struct {
            a: [*]const E,
            ld: usize,
            x: usize,
            y: usize,
            across: bool,
            upper: bool,
            unit: bool,

            pub fn at(s: Source, j: usize, l: usize) ?E {
                const line = s.y + j;
                const depth = s.x + l;
                const row = if (s.across) line else depth;
                const col = if (s.across) depth else line;
                if (if (s.upper) row > col else row < col) return null;
                if (s.unit and row == col) return if (complex) .{ 1, 0 } else 1;
                return s.a[row + col * s.ld];
            }
        };

        fn copy(u: usize, upper: bool, across: bool, unit: bool, m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) i32 {
            const src: Source = .{ .a = @ptrCast(a), .ld = @intCast(lda), .x = @intCast(posX), .y = @intCast(posY), .across = across, .upper = upper, .unit = unit };
            gk.packWith(E, u, @intCast(m), @intCast(n), @ptrCast(b), src);
            return 0;
        }

        /// The copy `trmm_<io><ul><nt><un>copy`.
        fn Copy(comptime u: usize, comptime upper: bool, comptime across: bool, comptime unit: bool) type {
            return struct {
                fn f(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
                    return copy(u, upper, across, unit, m, n, a, lda, posX, posY, b);
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

        /// The gemm kernel conjugating the triangle's side.
        fn gemm(left: bool, conj: bool, m: usize, n: usize, k: usize, alpha: [*]const F, a: [*]E, b: [*]E, c: [*]E, ldc: usize) void {
            const bm: BlasLong = @intCast(m);
            const bn: BlasLong = @intCast(n);
            const bk: BlasLong = @intCast(k);
            if (comptime complex) {
                const f = if (!conj) &G.gemm_kernel_n else if (left) &G.gemm_kernel_l else &G.gemm_kernel_r;
                _ = f(bm, bn, bk, alpha[0], alpha[1], @ptrCast(a), @ptrCast(b), @ptrCast(c), @intCast(ldc));
            } else {
                _ = G.gemm_kernel(bm, bn, bk, alpha[0], a, b, c, @intCast(ldc));
            }
        }

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
                    if (rest[1] > rest[0]) gemm(left, conj, w, wn, rest[1] - rest[0], alpha, pa + rest[0] * w, pb + rest[0] * wn, pc, ld);
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
                            } else sr += @as(A, x) * y;
                        }
                        const p = &pc[i + j * ld];
                        if (comptime complex) {
                            const ar: A = alpha[0];
                            const ai: A = alpha[1];
                            p.* = .{ @floatCast(p[0] + (ar * sr - ai * si)), @floatCast(p[1] + (ar * si + ai * sr)) };
                        } else p.* = @floatCast(p.* + @as(A, alpha[0]) * sr);
                    };
                    r0 += w;
                }
                c0 += wn;
            }
        }

        fn clip(v: isize, k: isize) usize {
            return @intCast(@max(0, @min(v, k)));
        }
    };
}
