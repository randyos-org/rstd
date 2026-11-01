//! The unblocked kernel of `trsyl` (`relapack.zig`): `op(A) X + isgn X op(B)
//! = scale C` for upper triangular `A` and `B`, real ones quasi-triangular
//! (Schur form, 2 x 2 diagonal blocks).
//!
//! `trsy2` is ReLAPACK's `RELAPACK_?trsyl_rec2` (`relapack/src/dtrsyl_rec2.c`,
//! `ztrsyl_rec2.c` and the rest; OpenBLAS v0.3.34-317), Reference-LAPACK's
//! `xTRSYL` (`lapack-netlib/SRC/dtrsyl.f`, `ztrsyl.f` in that checkout)
//! through f2c. Its 2 x 2 blocks call `xLALN2` (`trevc.zig`'s, the
//! real `NW = 1` paths of `dlaln2.f`) and `xLASY2` (`lasy2`, the 2 x 2 by
//! 2 x 2 path of `dlasy2.f`, the only one `xTRSYL` takes). The four
//! transposition cases are one loop: the order, the sums' ranges and which
//! side transposes follow `op`, as the reference's four copies do. The sums
//! are `?DOT` in order; complex division is Smith's (`util.cdiv`) where the
//! reference takes `ZLADIV`.
//!
//! `tgsy2` is `tgsyl`'s base case, Reference-LAPACK's `xTGSY2` with
//! `xGETC2`, `xGESC2` and `xLATDF` (`dtgsy2.f`, `ztgsy2.f` and the rest),
//! `xLATDF`'s `xGECON` from `norm.zig`'s `lacn2` and `latrs`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const chol = @import("cholesky.zig");
const trevc = @import("trevc.zig");
const norm = @import("norm.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
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

/// `DLALN2` with `NA = 2`, `NW = 1` (`trevc.laln2`): `(A' - wr I) x = s v`
/// for the 2 x 2 `A` at `a` (`A'` its transpose with `ltrans`).
fn laln2(comptime F: type, ltrans: bool, smin: F, a: [*]const F, lda: usize, v: [2]F, wr: F) Small(F) {
    const r = trevc.laln2(F, ltrans, 2, 1, smin, a, lda, .{ .{ v[0], 0 }, .{ v[1], 0 } }, wr, 0);
    return .{ .x = .{ r.x[0][0], r.x[1][0], 0, 0 }, .scale = r.scale, .perturbed = r.perturbed };
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

// ---- tgsyl's kernel ----------------------------------------------------------------

/// The `dif` estimate's running sum of squares (`RDSUM`, `RDSCAL`, as
/// `xLASSQ` keeps them) and the count of diagonal block pairs solved
/// (`PQ`).
pub fn Dif(comptime F: type) type {
    return struct { sum: F = 1, scale: F = 0, pq: usize = 0 };
}

/// `|x|`, the modulus for complex (Fortran's `ABS`).
fn absT(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) std.math.hypot(x.re, x.im) else @abs(x);
}

fn conjT(comptime T: type, x: T) T {
    return if (comptime isComplex(T)) x.conjugate() else x;
}

/// `xTGSY2` (`lapack-netlib/SRC/dtgsy2.f`, `ztgsy2.f`), ReLAPACK's base case
/// for `tgsyl`: `A R - L B = scale C`, `D R - L E = scale F` (with `trans`,
/// `Aᴴ R + Dᴴ L = scale C`, `R Bᴴ + L Eᴴ = -scale F`) for `R` in `C` and `L`
/// in `F` (`m x n`), `A`, `B` upper quasi-triangular (complex: triangular),
/// `D`, `E` upper triangular. Each pair of diagonal blocks is a system of 2,
/// 4 or 8 unknowns, `Z x = b` (`Zᴴ` with `trans`) with `Z`'s Kronecker
/// form (`dtgsyl.f`'s (2)), by `getc2` and `gesc2`; with `ijob` 1 or 2 (not
/// `trans`), `latdf` instead adds the block's part of the `dif` estimate to
/// `dif`, unscaled. The reference's four real block shapes are one loop.
pub fn tgsy2(al: Allocator, io: Io, comptime T: type, trans: bool, ijob: u2, m: usize, n: usize, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, c: [*]T, ldc: usize, d: [*]const T, ldd: usize, e: [*]const T, lde: usize, f: [*]T, ldf: usize, dif: *Dif(Float(T))) norm.Error!Solved(Float(T)) {
    const F = Float(T);
    var out: Solved(F) = .{ .scale = 1, .perturbed = false };
    var p: usize = 0;
    var i: usize = 0;
    while (i < m) : (i += blockOf(T, m, a, lda, i)) p += 1;
    var q: usize = 0;
    i = 0;
    while (i < n) : (i += blockOf(T, n, b, ldb, i)) q += 1;
    dif.pq += p * q;
    const zero = scalar(T, 0);
    // rows by `A`'s blocks backward (forward with `trans`), columns by
    // `B`'s forward (backward)
    var jl: usize = if (trans) n else 0;
    while (if (trans) jl > 0 else jl < n) {
        const js: usize, const nb: usize = if (!trans) .{ jl, blockOf(T, n, b, ldb, jl) } else lastBlock(T, b, ldb, jl);
        var il: usize = if (trans) 0 else m;
        while (if (trans) il < m else il > 0) {
            const is: usize, const mb: usize = if (trans) .{ il, blockOf(T, m, a, lda, il) } else lastBlock(T, a, lda, il);
            const k = mb * nb;
            const zdim = 2 * k;
            // `Z = [ kron(I, A_ii)  -kron(B_jjᵀ, I) ; kron(I, D_ii)  -kron(E_jjᵀ, I) ]`
            var z: [64]T = @splat(zero);
            for (0..nb) |jj| for (0..mb) |ii| {
                const r = ii + mb * jj;
                for (0..mb) |kk| {
                    z[r + 8 * (kk + mb * jj)] = a[is + ii + (is + kk) * lda];
                    if (kk >= ii) z[k + r + 8 * (kk + mb * jj)] = d[is + ii + (is + kk) * ldd];
                }
                for (0..nb) |kk| {
                    z[r + 8 * (k + ii + mb * kk)] = chol.neg(T, b[js + kk + (js + jj) * ldb]);
                    if (kk <= jj) z[k + r + 8 * (k + ii + mb * kk)] = chol.neg(T, e[js + kk + (js + jj) * lde]);
                }
            };
            if (trans) for (0..zdim) |col| for (col..zdim) |row| {
                const t = z[row + 8 * col];
                z[row + 8 * col] = conjT(T, z[col + 8 * row]);
                z[col + 8 * row] = conjT(T, t);
            };
            var rhs: [8]T = undefined;
            for (0..nb) |jj| for (0..mb) |ii| {
                rhs[ii + mb * jj] = c[is + ii + (js + jj) * ldc];
                rhs[k + ii + mb * jj] = f[is + ii + (js + jj) * ldf];
            };
            var ipiv: [8]usize = undefined;
            var jpiv: [8]usize = undefined;
            if (getc2(T, zdim, &z, &ipiv, &jpiv) != 0) out.perturbed = true;
            if (!trans and ijob != 0) {
                try latdf(al, io, T, ijob, zdim, &z, &rhs, dif, &ipiv, &jpiv);
            } else {
                const scaloc = gesc2(T, zdim, &z, &rhs, &ipiv, &jpiv);
                if (scaloc != 1) {
                    for (0..n) |col| for (0..m) |row| {
                        c[row + col * ldc] = chol.mulReal(T, c[row + col * ldc], scaloc);
                        f[row + col * ldf] = chol.mulReal(T, f[row + col * ldf], scaloc);
                    };
                    out.scale *= scaloc;
                }
            }
            for (0..nb) |jj| for (0..mb) |ii| {
                c[is + ii + (js + jj) * ldc] = rhs[ii + mb * jj];
                f[is + ii + (js + jj) * ldf] = rhs[k + ii + mb * jj];
            };
            // R(i, j) and L(i, j) into the rest
            const ie = is + mb;
            const je = js + nb;
            if (!trans) {
                for (0..nb) |jj| for (0..is) |row| {
                    var sc = zero;
                    var sf = zero;
                    for (0..mb) |kk| {
                        sc = add(T, sc, mul(T, a[row + (is + kk) * lda], rhs[kk + mb * jj]));
                        sf = add(T, sf, mul(T, d[row + (is + kk) * ldd], rhs[kk + mb * jj]));
                    }
                    c[row + (js + jj) * ldc] = sub(T, c[row + (js + jj) * ldc], sc);
                    f[row + (js + jj) * ldf] = sub(T, f[row + (js + jj) * ldf], sf);
                };
                for (je..n) |col| for (0..mb) |ii| {
                    var sc = zero;
                    var sf = zero;
                    for (0..nb) |kk| {
                        const l = rhs[k + ii + mb * kk];
                        sc = add(T, sc, mul(T, l, b[js + kk + col * ldb]));
                        sf = add(T, sf, mul(T, l, e[js + kk + col * lde]));
                    }
                    c[is + ii + col * ldc] = add(T, c[is + ii + col * ldc], sc);
                    f[is + ii + col * ldf] = add(T, f[is + ii + col * ldf], sf);
                };
            } else {
                for (0..js) |col| for (0..mb) |ii| {
                    var s = zero;
                    for (0..nb) |kk| {
                        s = add(T, s, mul(T, rhs[ii + mb * kk], conjT(T, b[col + (js + kk) * ldb])));
                        s = add(T, s, mul(T, rhs[k + ii + mb * kk], conjT(T, e[col + (js + kk) * lde])));
                    }
                    f[is + ii + col * ldf] = add(T, f[is + ii + col * ldf], s);
                };
                for (0..nb) |jj| for (ie..m) |row| {
                    var s = zero;
                    for (0..mb) |kk| {
                        s = add(T, s, mul(T, conjT(T, a[is + kk + row * lda]), rhs[kk + mb * jj]));
                        s = add(T, s, mul(T, conjT(T, d[is + kk + row * ldd]), rhs[k + kk + mb * jj]));
                    }
                    c[row + (js + jj) * ldc] = sub(T, c[row + (js + jj) * ldc], s);
                };
            }
            il = if (trans) ie else is;
        }
        jl = if (trans) js else js + nb;
    }
    return out;
}

/// The diagonal block of the quasi-triangular `a` that ends before `end`:
/// its first index and size.
fn lastBlock(comptime T: type, a: [*]const T, lda: usize, end: usize) [2]usize {
    const s: usize = if (comptime !isComplex(T)) (if (end >= 2 and a[end - 1 + (end - 2) * lda] != 0) 2 else 1) else 1;
    return .{ end - s, s };
}

/// `xGETC2` (`dgetc2.f`, `zgetc2.f`): `z` (`n x n`, leading dimension 8) to
/// `P L U Q` by complete pivoting, pivots below `SMIN` raised to it; the
/// last such pivot's index (from 1), or 0.
fn getc2(comptime T: type, n: usize, z: *[64]T, ipiv: *[8]usize, jpiv: *[8]usize) usize {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const smlnum = std.math.floatMin(F) / eps;
    var info: usize = 0;
    if (n == 1) {
        ipiv[0] = 0;
        jpiv[0] = 0;
        if (absT(T, z[0]) < smlnum) {
            info = 1;
            z[0] = scalar(T, smlnum);
        }
        return info;
    }
    var smin: F = 0;
    for (0..n - 1) |i| {
        var xmax: F = 0;
        var ipv = i;
        var jpv = i;
        for (i..n) |jp| for (i..n) |ip| if (absT(T, z[ip + 8 * jp]) >= xmax) {
            xmax = absT(T, z[ip + 8 * jp]);
            ipv = ip;
            jpv = jp;
        };
        if (i == 0) smin = @max(eps * xmax, smlnum);
        if (ipv != i) for (0..n) |col| std.mem.swap(T, &z[ipv + 8 * col], &z[i + 8 * col]);
        ipiv[i] = ipv;
        if (jpv != i) for (0..n) |row| std.mem.swap(T, &z[row + 8 * jpv], &z[row + 8 * i]);
        jpiv[i] = jpv;
        if (absT(T, z[i + 8 * i]) < smin) {
            info = i + 1;
            z[i + 8 * i] = scalar(T, smin);
        }
        for (i + 1..n) |j| z[j + 8 * i] = divide(T, z[j + 8 * i], z[i + 8 * i]);
        for (i + 1..n) |col| for (i + 1..n) |row| {
            z[row + 8 * col] = sub(T, z[row + 8 * col], mul(T, z[row + 8 * i], z[i + 8 * col]));
        };
    }
    if (absT(T, z[n - 1 + 8 * (n - 1)]) < smin) {
        info = n;
        z[n - 1 + 8 * (n - 1)] = scalar(T, smin);
    }
    ipiv[n - 1] = n - 1;
    jpiv[n - 1] = n - 1;
    return info;
}

/// `xGESC2` (`dgesc2.f`, `zgesc2.f`): `rhs` to `getc2`'s `z`'s solution
/// times the returned scale, `<= 1`, chosen against overflow.
fn gesc2(comptime T: type, n: usize, z: *const [64]T, rhs: *[8]T, ipiv: *const [8]usize, jpiv: *const [8]usize) Float(T) {
    const F = Float(T);
    const smlnum = std.math.floatMin(F) / std.math.floatEps(F);
    for (0..n - 1) |i| std.mem.swap(T, &rhs[i], &rhs[ipiv[i]]);
    for (0..n - 1) |i| for (i + 1..n) |j| {
        rhs[j] = sub(T, rhs[j], mul(T, z[j + 8 * i], rhs[i]));
    };
    var scale: F = 1;
    // `I?AMAX`: `|re| + |im|` for complex
    var big: usize = 0;
    for (0..n) |i| if (abs1(T, rhs[i]) > abs1(T, rhs[big])) {
        big = i;
    };
    if (2 * smlnum * absT(T, rhs[big]) > absT(T, z[n - 1 + 8 * (n - 1)])) {
        const temp = 0.5 / absT(T, rhs[big]);
        for (rhs[0..n]) |*x| x.* = chol.mulReal(T, x.*, temp);
        scale *= temp;
    }
    backSolve(T, n, z, rhs);
    unswap(T, n, rhs, jpiv);
    return scale;
}

/// `xLASWP` with `INCX = -1`: the swaps undone, last first.
fn unswap(comptime T: type, n: usize, x: *[8]T, piv: *const [8]usize) void {
    var i = n - 1;
    while (i > 0) {
        i -= 1;
        std.mem.swap(T, &x[i], &x[piv[i]]);
    }
}

/// `U x = rhs` in place, as `xGESC2` divides: `rhs(i) *= 1 / U(i, i)`,
/// then `rhs(i) -= rhs(j) (U(i, j) / U(i, i))`.
fn backSolve(comptime T: type, n: usize, z: *const [64]T, rhs: *[8]T) void {
    var i = n;
    while (i > 0) {
        i -= 1;
        const temp = divide(T, scalar(T, 1), z[i + 8 * i]);
        rhs[i] = mul(T, rhs[i], temp);
        for (i + 1..n) |j| rhs[i] = sub(T, rhs[i], mul(T, rhs[j], mul(T, z[i + 8 * j], temp)));
    }
}

/// The real part of `conj(x)ᵀ y` over `n` entries, `?DOT`'s or `ZDOTC`'s.
fn redot(comptime T: type, n: usize, x: [*]const T, y: [*]const T) Float(T) {
    var s: Float(T) = 0;
    for (0..n) |i| s += if (comptime isComplex(T)) x[i].re * y[i].re + x[i].im * y[i].im else x[i] * y[i];
    return s;
}

/// `xLATDF` (`dlatdf.f`, `zlatdf.f`): `Z x = b` for `getc2`'s `z` with `b`
/// chosen, from `rhs`, to make `x` large, then `x`'s sum of squares into
/// `dif`. `ijob` 1: `b`'s entries `±1` by look-ahead; 2: `rhs` plus or
/// minus an approximate null vector of `Z` (`xGECON`'s last `lacn2` vector).
fn latdf(al: Allocator, io: Io, comptime T: type, ijob: u2, n: usize, z: *const [64]T, rhs: *[8]T, dif: *Dif(Float(T)), ipiv: *const [8]usize, jpiv: *const [8]usize) norm.Error!void {
    const F = Float(T);
    const one = scalar(T, 1);
    if (ijob == 2) {
        var xm: [8]T = @splat(scalar(T, 0));
        try nullVector(al, io, T, n, z, &xm);
        unswap(T, n, &xm, ipiv);
        const temp = 1 / @sqrt(redot(T, n, &xm, &xm));
        var xp: [8]T = undefined;
        for (0..n) |k| {
            xm[k] = chol.mulReal(T, xm[k], temp);
            xp[k] = add(T, xm[k], rhs[k]);
            rhs[k] = sub(T, rhs[k], xm[k]);
        }
        _ = gesc2(T, n, z, rhs, ipiv, jpiv);
        _ = gesc2(T, n, z, &xp, ipiv, jpiv);
        var sp: F = 0;
        var sr: F = 0;
        for (0..n) |k| {
            sp += abs1(T, xp[k]);
            sr += abs1(T, rhs[k]);
        }
        if (sp > sr) rhs.* = xp;
        norm.lassq(T, n, rhs, 1, &dif.scale, &dif.sum);
        return;
    }
    for (0..n - 1) |i| std.mem.swap(T, &rhs[i], &rhs[ipiv[i]]);
    // the L part, each `b(j)` +1 or -1 by look-ahead
    var pmone: F = -1;
    for (0..n - 1) |j| {
        const col = z[8 * j ..].ptr;
        const bp = add(T, rhs[j], one);
        const bm = sub(T, rhs[j], one);
        var splus: F = 1 + redot(T, n - j - 1, col + j + 1, col + j + 1);
        const sminu = redot(T, n - j - 1, col + j + 1, rhs[j + 1 ..].ptr);
        splus *= if (comptime isComplex(T)) rhs[j].re else rhs[j];
        if (splus > sminu) {
            rhs[j] = bp;
        } else if (sminu > splus) {
            rhs[j] = bm;
        } else {
            // equal sums: -1 the first time, +1 after
            rhs[j] = add(T, rhs[j], scalar(T, pmone));
            pmone = 1;
        }
        const temp = chol.neg(T, rhs[j]);
        for (j + 1..n) |k| rhs[k] = add(T, rhs[k], mul(T, temp, col[k]));
    }
    // the U part, with look-ahead for b(n)
    var work = rhs.*;
    work[n - 1] = add(T, rhs[n - 1], one);
    rhs[n - 1] = sub(T, rhs[n - 1], one);
    var splus: F = 0;
    var sminu: F = 0;
    var i = n;
    while (i > 0) {
        i -= 1;
        const temp = divide(T, one, z[i + 8 * i]);
        work[i] = mul(T, work[i], temp);
        rhs[i] = mul(T, rhs[i], temp);
        for (i + 1..n) |k| {
            work[i] = sub(T, work[i], mul(T, work[k], mul(T, z[i + 8 * k], temp)));
            rhs[i] = sub(T, rhs[i], mul(T, rhs[k], mul(T, z[i + 8 * k], temp)));
        }
        splus += absT(T, work[i]);
        sminu += absT(T, rhs[i]);
    }
    if (splus > sminu) rhs.* = work;
    unswap(T, n, rhs, jpiv);
    norm.lassq(T, n, rhs, 1, &dif.scale, &dif.sum);
}

/// `xGECON('I', N, Z, LDZ, ONE, ...)` as `xLATDF` uses it: the vector
/// `lacn2` leaves in `xGECON`'s `WORK(N+1:2N)`, `(L U)⁻¹ w` for the `w` it
/// found, from `getc2`'s factors (its swaps left out, as there). An
/// overflowing solve ends it early, as `xGECON`'s `GO TO 20`.
fn nullVector(al: Allocator, io: Io, comptime T: type, n: usize, z: *const [64]T, v: *[8]T) norm.Error!void {
    const F = Float(T);
    const Inv = struct {
        al: Allocator,
        io: Io,
        n: usize,
        z: [*]const T,
        cnorm: [16]F = undefined,
        normin: bool = false,

        const Self = @This();
        const Stop = norm.Error || error{Overflow};

        // the infinity norm's `KASE1 = 2`: `KASE = 2` is `U⁻¹ L⁻¹`
        fn apply(s: *Self, adjoint: bool, x: [*]T) Stop!void {
            const op: norm.Op = if (comptime isComplex(T)) .c else .t;
            const lower: [*]F = &s.cnorm;
            const upper: [*]F = s.cnorm[8..].ptr;
            var sl: F = undefined;
            var su: F = undefined;
            if (adjoint) {
                sl = try norm.latrs(s.al, s.io, T, false, .n, true, s.normin, s.n, s.z, 8, x, lower);
                su = try norm.latrs(s.al, s.io, T, true, .n, false, s.normin, s.n, s.z, 8, x, upper);
            } else {
                su = try norm.latrs(s.al, s.io, T, true, op, false, s.normin, s.n, s.z, 8, x, upper);
                sl = try norm.latrs(s.al, s.io, T, false, op, true, s.normin, s.n, s.z, 8, x, lower);
            }
            const scale = sl * su;
            s.normin = true;
            if (scale != 1) {
                const ix = norm.iamax(T, s.n, x, 1);
                if (scale < norm.abs1(T, x[ix]) * std.math.floatMin(F) or scale == 0) return error.Overflow;
                try norm.rscl(s.io, T, s.n, scale, x, 1);
            }
        }
    };
    var ctx: Inv = .{ .al = al, .io = io, .n = n, .z = z };
    var x: [8]T = undefined;
    _ = norm.lacn2(al, T, Inv.Stop, n, v, &x, &ctx, Inv.apply) catch |e| switch (e) {
        error.Overflow => {},
        error.OutOfMemory, error.Canceled => |y| return y,
    };
}
