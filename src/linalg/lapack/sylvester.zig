//! The unblocked kernel of `trsyl` (`relapack.zig`): `op(A) X + isgn X op(B)
//! = scale C` for upper triangular `A` and `B`, real ones quasi-triangular
//! (Schur form, 2 x 2 diagonal blocks).
//!
//! `trsy2` is ReLAPACK's `RELAPACK_?trsyl_rec2` (`relapack/src/dtrsyl_rec2.c`,
//! `ztrsyl_rec2.c` and the rest; OpenBLAS v0.3.34-317), Reference-LAPACK's
//! `xTRSYL` (`lapack-netlib/SRC/dtrsyl.f`, `ztrsyl.f` in that checkout)
//! through f2c. Its 2 x 2 blocks call `xLALN2` (`laln2`, the real `NW = 1`
//! paths of `dlaln2.f`) and `xLASY2` (`lasy2`, the 2 x 2 by 2 x 2 path of
//! `dlasy2.f`, the only one `xTRSYL` takes). The four transposition cases
//! are one loop: the order, the sums' ranges and which side transposes
//! follow `op`, as the reference's four copies do. The sums are `?DOT` in
//! order; complex division is Smith's (`util.cdiv`) where the reference
//! takes `ZLADIV`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const chol = @import("cholesky.zig");

const Float = types.Float;
const isComplex = types.isComplex;
const scalar = chol.scalar;
const add = chol.add;
const sub = chol.sub;
const mul = chol.mul;

/// What `trsy2` returns: `scale` and whether a near-singular system was
/// perturbed (the reference's `INFO = 1`).
pub fn Solved(comptime F: type) type {
    return struct { scale: F, perturbed: bool };
}

/// `|re| + |im|`, `|x|` for real.
fn abs1(comptime T: type, v: T) Float(T) {
    return if (comptime isComplex(T)) @abs(v.re) + @abs(v.im) else @abs(v);
}

/// `xLANGE('M')`: the largest `|a(i, j)|` of the `n x n` `a`, the modulus
/// for complex.
fn maxAbs(comptime T: type, n: usize, a: [*]const T, lda: usize) Float(T) {
    var m: Float(T) = 0;
    for (0..n) |j| for (0..n) |i| {
        const v = a[i + j * lda];
        const x = if (comptime isComplex(T)) std.math.hypot(v.re, v.im) else @abs(v);
        if (x > m or std.math.isNan(x)) m = x;
    };
    return m;
}

/// The diagonal blocks of the `n x n` upper quasi-triangular `a`: their
/// first indices, 1 x 1 unless a subdiagonal is nonzero (real only).
fn blockOf(comptime T: type, n: usize, a: [*]const T, lda: usize, i: usize) usize {
    if (comptime isComplex(T)) return 1;
    return if (i + 1 < n and a[i + 1 + i * lda] != 0) 2 else 1;
}

/// `xTRSYL`: `ta`, `tb` transpose (complex: conjugate-transpose) `A`, `B`.
pub fn trsy2(comptime T: type, ta: bool, tb: bool, isgn: isize, m: usize, n: usize, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, c: [*]T, ldc: usize) Solved(Float(T)) {
    const F = Float(T);
    var out: Solved(F) = .{ .scale = 1, .perturbed = false };
    if (m == 0 or n == 0) return out;
    const eps = std.math.floatEps(F);
    const smlnum = std.math.floatMin(F) * @as(F, @floatFromInt(m * n)) / eps;
    const bignum = 1 / smlnum;
    const smin = @max(smlnum, eps * maxAbs(T, m, a, lda), eps * maxAbs(T, n, b, ldb));
    const sgn: F = @floatFromInt(isgn);
    // the blocks of `B` (columns of `C`) and `A` (rows), in order
    var l: usize = if (tb) n else 0;
    while (if (tb) l > 0 else l < n) {
        // backward from the last block for `tb`, forward otherwise
        const l1, const l2 = blk: {
            if (!tb) {
                const s = blockOf(T, n, b, ldb, l);
                break :blk .{ l, l + s - 1 };
            }
            const s: usize = if (comptime !isComplex(T)) (if (l >= 2 and b[l - 1 + (l - 2) * ldb] != 0) 2 else 1) else 1;
            break :blk .{ l - s, l - 1 };
        };
        var k: usize = if (ta) 0 else m;
        while (if (ta) k < m else k > 0) {
            const k1, const k2 = blk: {
                if (ta) {
                    const s = blockOf(T, m, a, lda, k);
                    break :blk .{ k, k + s - 1 };
                }
                const s: usize = if (comptime !isComplex(T)) (if (k >= 2 and a[k - 1 + (k - 2) * lda] != 0) 2 else 1) else 1;
                break :blk .{ k - s, k - 1 };
            };
            const ctx: Ctx(T) = .{ .ta = ta, .tb = tb, .m = m, .n = n, .a = a, .lda = lda, .b = b, .ldb = ldb, .c = c, .ldc = ldc, .k1 = k1, .k2 = k2, .l1 = l1, .l2 = l2, .sgn = sgn };
            var x: [4]T = undefined;
            var scaloc: F = 1;
            if (l1 == l2 and k1 == k2) {
                const vec = ctx.rhs(k1, l1);
                var a11 = ctx.diag(k1, l1);
                var da11 = abs1(T, a11);
                if (da11 <= smin) {
                    a11 = scalar(T, smin);
                    da11 = smin;
                    out.perturbed = true;
                }
                const db = abs1(T, vec);
                if (da11 < 1 and db > 1 and db > bignum * da11) scaloc = 1 / db;
                x[0] = divide(T, chol.mulReal(T, vec, scaloc), a11);
            } else if (comptime !isComplex(T)) {
                const r = if (l1 == l2) blk: {
                    // a 2 x 1 block of `A`
                    const v = [2]F{ ctx.rhs(k1, l1), ctx.rhs(k2, l1) };
                    break :blk laln2(F, ta, smin, a + k1 + k1 * lda, lda, v, -sgn * b[l1 + l1 * ldb]);
                } else if (k1 == k2) blk: {
                    // a 1 x 2 block of `B`, from the right
                    const v = [2]F{ sgn * ctx.rhs(k1, l1), sgn * ctx.rhs(k1, l2) };
                    break :blk laln2(F, !tb, smin, b + l1 + l1 * ldb, ldb, v, -sgn * a[k1 + k1 * lda]);
                } else blk: {
                    const v = [4]F{ ctx.rhs(k1, l1), ctx.rhs(k2, l1), ctx.rhs(k1, l2), ctx.rhs(k2, l2) };
                    break :blk lasy2(F, ta, tb, sgn, a + k1 + k1 * lda, lda, b + l1 + l1 * ldb, ldb, v);
                };
                if (r.perturbed) out.perturbed = true;
                scaloc = r.scale;
                x = r.x;
            }
            if (scaloc != 1) {
                for (0..n) |j| for (0..m) |i| {
                    c[i + j * ldc] = chol.mulReal(T, c[i + j * ldc], scaloc);
                };
                out.scale *= scaloc;
            }
            if (l1 == l2 and k1 == k2) {
                c[k1 + l1 * ldc] = x[0];
            } else if (l1 == l2) {
                c[k1 + l1 * ldc] = x[0];
                c[k2 + l1 * ldc] = x[1];
            } else if (k1 == k2) {
                c[k1 + l1 * ldc] = x[0];
                c[k1 + l2 * ldc] = x[1];
            } else {
                c[k1 + l1 * ldc] = x[0];
                c[k2 + l1 * ldc] = x[1];
                c[k1 + l2 * ldc] = x[2];
                c[k2 + l2 * ldc] = x[3];
            }
            k = if (ta) k2 + 1 else k1;
        }
        l = if (tb) l1 else l2 + 1;
    }
    return out;
}

/// One block's view of the problem.
fn Ctx(comptime T: type) type {
    return struct {
        ta: bool,
        tb: bool,
        m: usize,
        n: usize,
        a: [*]const T,
        lda: usize,
        b: [*]const T,
        ldb: usize,
        c: [*]T,
        ldc: usize,
        k1: usize,
        k2: usize,
        l1: usize,
        l2: usize,
        sgn: Float(T),

        const Self = @This();

        fn cj(v: T) T {
            return if (comptime isComplex(T)) v.conjugate() else v;
        }

        /// `C(i, j) - (suml + sgn sumr)`, the solved rows (`op(A)`) and
        /// columns (`op(B)`) summed as the reference's `?DOT`s.
        fn rhs(s: Self, i: usize, j: usize) T {
            var suml = scalar(T, 0);
            if (s.ta) {
                // `A(p, i)` conjugated, rows before the block
                for (0..s.k1) |p| suml = add(T, suml, mul(T, cj(s.a[p + i * s.lda]), s.c[p + j * s.ldc]));
            } else {
                for (s.k2 + 1..s.m) |p| suml = add(T, suml, mul(T, s.a[i + p * s.lda], s.c[p + j * s.ldc]));
            }
            var sumr = scalar(T, 0);
            if (s.tb) {
                // `ZDOTC(C(i, l2+1:), B(j, l2+1:))` conjugated back
                for (s.l2 + 1..s.n) |q| sumr = add(T, sumr, mul(T, cj(s.c[i + q * s.ldc]), s.b[j + q * s.ldb]));
                sumr = cj(sumr);
            } else {
                for (0..s.l1) |q| sumr = add(T, sumr, mul(T, s.c[i + q * s.ldc], s.b[q + j * s.ldb]));
            }
            return sub(T, s.c[i + j * s.ldc], add(T, suml, chol.mulReal(T, sumr, s.sgn)));
        }

        /// `op(A)(i, i) + sgn op(B)(j, j)`, as each case of `ZTRSYL` forms it.
        fn diag(s: Self, i: usize, j: usize) T {
            const aii = s.a[i + i * s.lda];
            const bjj = s.b[j + j * s.ldb];
            if (s.ta and s.tb) return cj(add(T, aii, chol.mulReal(T, bjj, s.sgn)));
            return add(T, if (s.ta) cj(aii) else aii, chol.mulReal(T, if (s.tb) cj(bjj) else bjj, s.sgn));
        }
    };
}

fn divide(comptime T: type, x: T, d: T) T {
    if (comptime !isComplex(T)) return x / d;
    const z = util.cdiv(Float(T), x.re, x.im, d.re, d.im);
    return .init(z[0], z[1]);
}

/// A small solve's result: `x` (column-major), its scale, and whether a
/// pivot was perturbed.
pub fn Small(comptime F: type) type {
    return struct { x: [4]F, scale: F, perturbed: bool };
}

/// `DLALN2` with `NA = 2`, `NW = 1`, `CA`, `D1`, `D2` 1: `(A' - wr I) x = s v`
/// for the 2 x 2 `A` at `a` (`A'` its transpose with `ltrans`), by complete
/// pivoting; `s <= 1` keeps `x` from overflowing.
fn laln2(comptime F: type, ltrans: bool, smin: F, a: [*]const F, lda: usize, v: [2]F, wr: F) Small(F) {
    // `IPIVOT`, `RSWAP`, `ZSWAP` (`dlaln2.f`)
    const ipivot = [4][4]usize{ .{ 0, 1, 2, 3 }, .{ 1, 0, 3, 2 }, .{ 2, 3, 0, 1 }, .{ 3, 2, 1, 0 } };
    const rswap = [4]bool{ false, true, false, true };
    const zswap = [4]bool{ false, false, true, true };
    const smlnum = 2 * std.math.floatMin(F);
    const bignum = 1 / smlnum;
    const smini = @max(smin, smlnum);
    var out: Small(F) = .{ .x = @splat(0), .scale = 1, .perturbed = false };
    // `CR`, column-major
    var cr: [4]F = undefined;
    cr[0] = a[0] - wr;
    cr[3] = a[1 + lda] - wr;
    if (ltrans) {
        cr[2] = a[1];
        cr[1] = a[lda];
    } else {
        cr[1] = a[1];
        cr[2] = a[lda];
    }
    var cmax: F = 0;
    var icmax: usize = 0;
    for (cr, 0..) |x, j| if (@abs(x) > cmax) {
        cmax = @abs(x);
        icmax = j;
    };
    if (cmax < smini) {
        const bnorm = @max(@abs(v[0]), @abs(v[1]));
        if (smini < 1 and bnorm > 1 and bnorm > bignum * smini) out.scale = 1 / bnorm;
        const temp = out.scale / smini;
        out.x[0] = temp * v[0];
        out.x[1] = temp * v[1];
        out.perturbed = true;
        return out;
    }
    const ur11 = cr[icmax];
    const cr21 = cr[ipivot[icmax][1]];
    const ur12 = cr[ipivot[icmax][2]];
    const cr22 = cr[ipivot[icmax][3]];
    const ur11r = 1 / ur11;
    const lr21 = ur11r * cr21;
    var ur22 = cr22 - ur12 * lr21;
    if (@abs(ur22) < smini) {
        ur22 = smini;
        out.perturbed = true;
    }
    const br1 = if (rswap[icmax]) v[1] else v[0];
    var br2 = if (rswap[icmax]) v[0] else v[1];
    br2 = br2 - lr21 * br1;
    const bbnd = @max(@abs(br1 * (ur22 * ur11r)), @abs(br2));
    if (bbnd > 1 and @abs(ur22) < 1 and bbnd >= bignum * @abs(ur22)) out.scale = 1 / bbnd;
    const xr2 = (br2 * out.scale) / ur22;
    const xr1 = (out.scale * br1) * ur11r - xr2 * (ur11r * ur12);
    out.x[0] = if (zswap[icmax]) xr2 else xr1;
    out.x[1] = if (zswap[icmax]) xr1 else xr2;
    const xnorm = @max(@abs(xr1), @abs(xr2));
    if (xnorm > 1 and cmax > 1 and xnorm > bignum / cmax) {
        const temp = cmax / bignum;
        out.x[0] *= temp;
        out.x[1] *= temp;
        out.scale *= temp;
    }
    return out;
}

/// `DLASY2` with `N1 = N2 = 2`: `op(TL) X + sgn X op(TR) = s V` for 2 x 2
/// blocks, as a 4 x 4 system by complete pivoting.
pub fn lasy2(comptime F: type, ltranl: bool, ltranr: bool, sgn: F, tl: [*]const F, ldtl: usize, tr: [*]const F, ldtr: usize, v: [4]F) Small(F) {
    const eps = std.math.floatEps(F);
    const smlnum = std.math.floatMin(F) / eps;
    var out: Small(F) = .{ .x = @splat(0), .scale = 1, .perturbed = false };
    const l11 = tl[0];
    const l21 = tl[1];
    const l12 = tl[ldtl];
    const l22 = tl[1 + ldtl];
    const r11 = tr[0];
    const r21 = tr[1];
    const r12 = tr[ldtr];
    const r22 = tr[1 + ldtr];
    var smin = @max(@max(@abs(r11), @abs(r12)), @max(@abs(r21), @abs(r22)));
    smin = @max(smin, @max(@max(@abs(l11), @abs(l12)), @max(@abs(l21), @abs(l22))));
    smin = @max(eps * smin, smlnum);
    // `T16`, row-major here: `t[row][col]`
    var t: [4][4]F = @splat(@splat(0));
    t[0][0] = l11 + sgn * r11;
    t[1][1] = l22 + sgn * r11;
    t[2][2] = l11 + sgn * r22;
    t[3][3] = l22 + sgn * r22;
    if (ltranl) {
        t[0][1] = l21;
        t[1][0] = l12;
        t[2][3] = l21;
        t[3][2] = l12;
    } else {
        t[0][1] = l12;
        t[1][0] = l21;
        t[2][3] = l12;
        t[3][2] = l21;
    }
    if (ltranr) {
        t[0][2] = sgn * r12;
        t[1][3] = sgn * r12;
        t[2][0] = sgn * r21;
        t[3][1] = sgn * r21;
    } else {
        t[0][2] = sgn * r21;
        t[1][3] = sgn * r21;
        t[2][0] = sgn * r12;
        t[3][1] = sgn * r12;
    }
    var bt = v;
    var jpiv: [3]usize = undefined;
    for (0..3) |i| {
        var xmax: F = 0;
        var ipsv: usize = i;
        var jpsv: usize = i;
        for (i..4) |ip| for (i..4) |jp| if (@abs(t[ip][jp]) >= xmax) {
            xmax = @abs(t[ip][jp]);
            ipsv = ip;
            jpsv = jp;
        };
        if (ipsv != i) {
            std.mem.swap([4]F, &t[ipsv], &t[i]);
            std.mem.swap(F, &bt[i], &bt[ipsv]);
        }
        if (jpsv != i) for (0..4) |r| std.mem.swap(F, &t[r][jpsv], &t[r][i]);
        jpiv[i] = jpsv;
        if (@abs(t[i][i]) < smin) {
            out.perturbed = true;
            t[i][i] = smin;
        }
        for (i + 1..4) |j| {
            t[j][i] = t[j][i] / t[i][i];
            bt[j] = bt[j] - t[j][i] * bt[i];
            for (i + 1..4) |k| t[j][k] = t[j][k] - t[j][i] * t[i][k];
        }
    }
    if (@abs(t[3][3]) < smin) {
        out.perturbed = true;
        t[3][3] = smin;
    }
    const big = @max(@max(@abs(bt[0]), @abs(bt[1])), @max(@abs(bt[2]), @abs(bt[3])));
    for (0..4) |i| if ((8 * smlnum) * @abs(bt[i]) > @abs(t[i][i])) {
        out.scale = (1.0 / 8.0) / big;
        for (&bt) |*x| x.* *= out.scale;
        break;
    };
    var tmp: [4]F = undefined;
    for (0..4) |ii| {
        const k = 3 - ii;
        const temp = 1 / t[k][k];
        tmp[k] = bt[k] * temp;
        for (k + 1..4) |j| tmp[k] = tmp[k] - (temp * t[k][j]) * tmp[j];
    }
    for (0..3) |ii| {
        const i = 2 - ii;
        if (jpiv[i] != i) std.mem.swap(F, &tmp[i], &tmp[jpiv[i]]);
    }
    out.x = tmp;
    return out;
}
