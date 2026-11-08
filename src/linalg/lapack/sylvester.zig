//! The kernels of `trsyl` (`trsyl.zig`) and `tgsyl` (`tgsyl.zig`):
//! Reference-LAPACK's `xTRSYL` (`trsy2`), `xLALN2`, `xLASY2`, and `xTGSY2`
//! with `xGETC2`, `xGESC2` and `xLATDF` (`lapack-netlib/SRC/dtrsyl.f`,
//! `dlaln2.f`, `dlasy2.f`, `dtgsy2.f` and their forms, OpenBLAS
//! v0.3.34-317). Complex division is Smith's (`util.cdiv`) where the
//! reference takes `ZLADIV` or `DLADIV`.

const std = @import("std");
const types = @import("../blas/types.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const norm = @import("norm.zig");
const util = @import("../blas/util.zig");
const driver = @import("../blas/blas3/driver.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const scalar = ops.scalar;
const add = ops.add;
const sub = ops.sub;
const mul = ops.mul;

const Error = driver.Error;

/// What `trsy2` returns: `scale` and whether a near-singular system was
/// perturbed (the reference's `INFO = 1`).
pub fn Solved(comptime F: type) type {
    return struct { scale: F, perturbed: bool };
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
                var da11 = ops.abs1(T, a11);
                if (da11 <= smin) {
                    a11 = scalar(T, smin);
                    da11 = smin;
                    out.perturbed = true;
                }
                const db = ops.abs1(T, vec);
                if (da11 < 1 and db > 1 and db > bignum * da11) scaloc = 1 / db;
                x[0] = ops.div(T, ops.mulReal(T, vec, scaloc), a11);
            } else if (comptime !isComplex(T)) {
                const r = if (l1 == l2) blk: {
                    // a 2 x 1 block of `A`
                    const v = [2]F{ ctx.rhs(k1, l1), ctx.rhs(k2, l1) };
                    break :blk laln2Small(F, ta, smin, a + k1 + k1 * lda, lda, v, -sgn * b[l1 + l1 * ldb]);
                } else if (k1 == k2) blk: {
                    // a 1 x 2 block of `B`, from the right
                    const v = [2]F{ sgn * ctx.rhs(k1, l1), sgn * ctx.rhs(k1, l2) };
                    break :blk laln2Small(F, !tb, smin, b + l1 + l1 * ldb, ldb, v, -sgn * a[k1 + k1 * lda]);
                } else blk: {
                    const v = [4]F{ ctx.rhs(k1, l1), ctx.rhs(k2, l1), ctx.rhs(k1, l2), ctx.rhs(k2, l2) };
                    break :blk lasy2Square(F, ta, tb, sgn, a + k1 + k1 * lda, lda, b + l1 + l1 * ldb, ldb, v);
                };
                if (r.perturbed) out.perturbed = true;
                scaloc = r.scale;
                x = r.x;
            }
            if (scaloc != 1) {
                for (0..n) |j| for (0..m) |i| {
                    c[i + j * ldc] = ops.mulReal(T, c[i + j * ldc], scaloc);
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
            return sub(T, s.c[i + j * s.ldc], add(T, suml, ops.mulReal(T, sumr, s.sgn)));
        }

        /// `op(A)(i, i) + sgn op(B)(j, j)`, as each case of `ZTRSYL` forms it.
        fn diag(s: Self, i: usize, j: usize) T {
            const aii = s.a[i + i * s.lda];
            const bjj = s.b[j + j * s.ldb];
            if (s.ta and s.tb) return cj(add(T, aii, ops.mulReal(T, bjj, s.sgn)));
            return add(T, if (s.ta) cj(aii) else aii, ops.mulReal(T, if (s.tb) cj(bjj) else bjj, s.sgn));
        }
    };
}

/// A small solve's result: `x` (column-major), its scale, and whether a
/// pivot was perturbed.
pub fn Small(comptime F: type) type {
    return struct { x: [4]F, scale: F, perturbed: bool };
}

/// `laln2` with `na = 2`, `nw = 1`: `(A' - wr I) x = s v`
/// for the 2 x 2 `A` at `a` (`A'` its transpose with `ltrans`).
fn laln2Small(comptime F: type, ltrans: bool, smin: F, a: [*]const F, lda: usize, v: [2]F, wr: F) Small(F) {
    const r = laln2(F, ltrans, 2, 1, smin, a, lda, .{ .{ v[0], 0 }, .{ v[1], 0 } }, wr, 0);
    return .{ .x = .{ r.x[0][0], r.x[1][0], 0, 0 }, .scale = r.scale, .perturbed = r.perturbed };
}

/// What `laln2` returns: `x[i][k]` the solution's row `i`, real part
/// (`k = 0`) or imaginary (`k = 1`); `scale <= 1`; `xnorm` its infinity
/// norm; `perturbed` when a pivot was raised to `smin` (`INFO = 1`).
pub fn Laln2(comptime F: type) type {
    return struct { x: [2][2]F, scale: F, xnorm: F, perturbed: bool };
}

/// `DLALN2`: `(A - (wr + i wi)) X = scale B` (`Aᵀ` with `ltrans`) for the
/// `na x na` `A` at `a` (`na` 1 or 2), `X` and `B` real (`nw = 1`) or
/// complex (`nw = 2`, `b[i][1]` the imaginary parts), by complete pivoting;
/// pivots below `smin` become `smin`.
pub fn laln2(comptime F: type, ltrans: bool, na: usize, nw: usize, smin: F, a: [*]const F, lda: usize, b: [2][2]F, wr: F, wi: F) Laln2(F) {
    const ipivot = [4][4]usize{ .{ 0, 1, 2, 3 }, .{ 1, 0, 3, 2 }, .{ 2, 3, 0, 1 }, .{ 3, 2, 1, 0 } };
    const rswap = [4]bool{ false, true, false, true };
    const zswap = [4]bool{ false, false, true, true };
    const smlnum = 2 * std.math.floatMin(F);
    const bignum = 1 / smlnum;
    const smini = @max(smin, smlnum);
    var out: Laln2(F) = .{ .x = @splat(@splat(0)), .scale = 1, .xnorm = 0, .perturbed = false };
    if (na == 1) {
        var csr = a[0] - wr;
        var csi: F = if (nw == 1) 0 else -wi;
        var cnorm = @abs(csr) + @abs(csi);
        if (cnorm < smini) {
            csr = smini;
            csi = 0;
            cnorm = smini;
            out.perturbed = true;
        }
        const bnorm = @abs(b[0][0]) + (if (nw == 1) 0 else @abs(b[0][1]));
        if (cnorm < 1 and bnorm > 1 and bnorm > bignum * cnorm) out.scale = 1 / bnorm;
        if (nw == 1) {
            out.x[0][0] = (b[0][0] * out.scale) / csr;
            out.xnorm = @abs(out.x[0][0]);
        } else {
            const z = util.cdiv(F, out.scale * b[0][0], out.scale * b[0][1], csr, csi);
            out.x[0] = z;
            out.xnorm = @abs(z[0]) + @abs(z[1]);
        }
        return out;
    }
    // `CR`, `CI` column-major as `CRV`, `CIV`
    var crv: [4]F = .{ a[0] - wr, a[1], a[lda], a[1 + lda] - wr };
    if (ltrans) std.mem.swap(F, &crv[1], &crv[2]);
    const civ: [4]F = if (nw == 1) @splat(0) else .{ -wi, 0, 0, -wi };
    var cmax: F = 0;
    var icmax: usize = 0;
    for (0..4) |j| if (@abs(crv[j]) + @abs(civ[j]) > cmax) {
        cmax = @abs(crv[j]) + @abs(civ[j]);
        icmax = j;
    };
    if (cmax < smini) {
        const bnorm = @max(@abs(b[0][0]) + @abs(b[0][1]), @abs(b[1][0]) + @abs(b[1][1]));
        if (smini < 1 and bnorm > 1 and bnorm > bignum * smini) out.scale = 1 / bnorm;
        const temp = out.scale / smini;
        for (0..2) |i| for (0..nw) |k| {
            out.x[i][k] = temp * b[i][k];
        };
        out.xnorm = temp * bnorm;
        out.perturbed = true;
        return out;
    }
    const piv = ipivot[icmax];
    if (nw == 1) {
        const ur11 = crv[icmax];
        const cr21 = crv[piv[1]];
        const ur12 = crv[piv[2]];
        const cr22 = crv[piv[3]];
        const ur11r = 1 / ur11;
        const lr21 = ur11r * cr21;
        var ur22 = cr22 - ur12 * lr21;
        if (@abs(ur22) < smini) {
            ur22 = smini;
            out.perturbed = true;
        }
        const br1 = if (rswap[icmax]) b[1][0] else b[0][0];
        var br2 = if (rswap[icmax]) b[0][0] else b[1][0];
        br2 = br2 - lr21 * br1;
        const bbnd = @max(@abs(br1 * (ur22 * ur11r)), @abs(br2));
        if (bbnd > 1 and @abs(ur22) < 1 and bbnd >= bignum * @abs(ur22)) out.scale = 1 / bbnd;
        const xr2 = (br2 * out.scale) / ur22;
        const xr1 = (out.scale * br1) * ur11r - xr2 * (ur11r * ur12);
        out.x[0][0] = if (zswap[icmax]) xr2 else xr1;
        out.x[1][0] = if (zswap[icmax]) xr1 else xr2;
        out.xnorm = @max(@abs(xr1), @abs(xr2));
    } else {
        const ur11 = crv[icmax];
        const ui11 = civ[icmax];
        const cr21 = crv[piv[1]];
        const ci21 = civ[piv[1]];
        const ur12 = crv[piv[2]];
        const ui12 = civ[piv[2]];
        const cr22 = crv[piv[3]];
        const ci22 = civ[piv[3]];
        var ur11r: F = undefined;
        var ui11r: F = undefined;
        var lr21: F = undefined;
        var li21: F = undefined;
        var ur12s: F = undefined;
        var ui12s: F = undefined;
        var ur22: F = undefined;
        var ui22: F = undefined;
        if (icmax == 0 or icmax == 3) {
            // the pivoted off-diagonals are real
            if (@abs(ur11) > @abs(ui11)) {
                const temp = ui11 / ur11;
                ur11r = 1 / (ur11 * (1 + temp * temp));
                ui11r = -temp * ur11r;
            } else {
                const temp = ur11 / ui11;
                ui11r = -1 / (ui11 * (1 + temp * temp));
                ur11r = -temp * ui11r;
            }
            lr21 = cr21 * ur11r;
            li21 = cr21 * ui11r;
            ur12s = ur12 * ur11r;
            ui12s = ur12 * ui11r;
            ur22 = cr22 - ur12 * lr21;
            ui22 = ci22 - ur12 * li21;
        } else {
            // the pivoted diagonals are real
            ur11r = 1 / ur11;
            ui11r = 0;
            lr21 = cr21 * ur11r;
            li21 = ci21 * ur11r;
            ur12s = ur12 * ur11r;
            ui12s = ui12 * ur11r;
            ur22 = cr22 - ur12 * lr21 + ui12 * li21;
            ui22 = -ur12 * li21 - ui12 * lr21;
        }
        const u22abs = @abs(ur22) + @abs(ui22);
        if (u22abs < smini) {
            ur22 = smini;
            ui22 = 0;
            out.perturbed = true;
        }
        const r1: usize, const r2: usize = if (rswap[icmax]) .{ 1, 0 } else .{ 0, 1 };
        var br1 = b[r1][0];
        var bi1 = b[r1][1];
        var br2 = b[r2][0];
        var bi2 = b[r2][1];
        br2 = br2 - lr21 * br1 + li21 * bi1;
        bi2 = bi2 - li21 * br1 - lr21 * bi1;
        const bbnd = @max((@abs(br1) + @abs(bi1)) * (u22abs * (@abs(ur11r) + @abs(ui11r))), @abs(br2) + @abs(bi2));
        if (bbnd > 1 and u22abs < 1 and bbnd >= bignum * u22abs) {
            out.scale = 1 / bbnd;
            br1 *= out.scale;
            bi1 *= out.scale;
            br2 *= out.scale;
            bi2 *= out.scale;
        }
        const x2 = util.cdiv(F, br2, bi2, ur22, ui22);
        const xr1 = ur11r * br1 - ui11r * bi1 - ur12s * x2[0] + ui12s * x2[1];
        const xi1 = ui11r * br1 + ur11r * bi1 - ui12s * x2[0] - ur12s * x2[1];
        const x1 = [2]F{ xr1, xi1 };
        out.x[0] = if (zswap[icmax]) x2 else x1;
        out.x[1] = if (zswap[icmax]) x1 else x2;
        out.xnorm = @max(@abs(xr1) + @abs(xi1), @abs(x2[0]) + @abs(x2[1]));
    }
    // further scaling if `‖A‖ ‖X‖` would overflow
    if (out.xnorm > 1 and cmax > 1 and out.xnorm > bignum / cmax) {
        const temp = cmax / bignum;
        for (&out.x) |*row| for (row) |*v| {
            v.* *= temp;
        };
        out.xnorm *= temp;
        out.scale *= temp;
    }
    return out;
}

/// `LOCU12`, `LOCL21`, `LOCU22`, `XSWPIV`, `BSWPIV` (`dlasy2.f`), 0-based.
const locu12 = [4]usize{ 2, 3, 0, 1 };
const locl21 = [4]usize{ 1, 0, 3, 2 };
const locu22 = [4]usize{ 3, 2, 1, 0 };
const xswpiv = [4]bool{ false, false, true, true };
const bswpiv = [4]bool{ false, true, false, true };

/// `DLASY2` (`dlasy2.f`): `op(TL) X + sgn X op(TR) = scale B` for `TL`
/// `n1 x n1` and `TR` `n2 x n2` (each 1 or 2). `x` is column-major 2 x 2
/// (`X(1,1)`, `X(2,1)`, `X(1,2)`, `X(2,2)`).
pub fn lasy2(comptime F: type, ltranl: bool, ltranr: bool, sgn: F, n1: usize, n2: usize, tl: [*]const F, ldtl: usize, tr: [*]const F, ldtr: usize, b: [*]const F, ldb: usize) Small(F) {
    const eps = std.math.floatEps(F);
    const smlnum = std.math.floatMin(F) / eps;
    var out: Small(F) = .{ .x = @splat(0), .scale = 1, .perturbed = false };
    if (n1 == 2 and n2 == 2) return lasy2Square(F, ltranl, ltranr, sgn, tl, ldtl, tr, ldtr, .{ b[0], b[1], b[ldb], b[1 + ldb] });
    if (n1 == 1 and n2 == 1) {
        var tau1 = tl[0] + sgn * tr[0];
        var bet = @abs(tau1);
        if (bet <= smlnum) {
            tau1 = smlnum;
            bet = smlnum;
            out.perturbed = true;
        }
        const gam = @abs(b[0]);
        if (smlnum * gam > bet) out.scale = 1 / gam;
        out.x[0] = (b[0] * out.scale) / tau1;
        return out;
    }
    var tmp: [4]F = undefined;
    var btmp: [2]F = undefined;
    var smin: F = undefined;
    if (n1 == 1) {
        smin = @max(eps * @max(@abs(tl[0]), @abs(tr[0]), @abs(tr[ldtr]), @abs(tr[1]), @abs(tr[1 + ldtr])), smlnum);
        tmp[0] = tl[0] + sgn * tr[0];
        tmp[3] = tl[0] + sgn * tr[1 + ldtr];
        tmp[1] = sgn * (if (ltranr) tr[1] else tr[ldtr]);
        tmp[2] = sgn * (if (ltranr) tr[ldtr] else tr[1]);
        btmp = .{ b[0], b[ldb] };
    } else {
        smin = @max(eps * @max(@abs(tr[0]), @abs(tl[0]), @abs(tl[ldtl]), @abs(tl[1]), @abs(tl[1 + ldtl])), smlnum);
        tmp[0] = tl[0] + sgn * tr[0];
        tmp[3] = tl[1 + ldtl] + sgn * tr[0];
        tmp[1] = if (ltranl) tl[ldtl] else tl[1];
        tmp[2] = if (ltranl) tl[1] else tl[ldtl];
        btmp = .{ b[0], b[1] };
    }
    // `IDAMAX` of `TMP`
    var ipiv: usize = 0;
    for (1..4) |i| if (@abs(tmp[i]) > @abs(tmp[ipiv])) {
        ipiv = i;
    };
    var p11 = tmp[ipiv];
    if (@abs(p11) <= smin) {
        out.perturbed = true;
        p11 = smin;
    }
    const p12 = tmp[locu12[ipiv]];
    const l21 = tmp[locl21[ipiv]] / p11;
    var p22 = tmp[locu22[ipiv]] - p12 * l21;
    if (@abs(p22) <= smin) {
        out.perturbed = true;
        p22 = smin;
    }
    if (bswpiv[ipiv]) {
        const t = btmp[1];
        btmp[1] = btmp[0] - l21 * t;
        btmp[0] = t;
    } else {
        btmp[1] = btmp[1] - l21 * btmp[0];
    }
    if ((2 * smlnum) * @abs(btmp[1]) > @abs(p22) or (2 * smlnum) * @abs(btmp[0]) > @abs(p11)) {
        out.scale = 0.5 / @max(@abs(btmp[0]), @abs(btmp[1]));
        btmp[0] *= out.scale;
        btmp[1] *= out.scale;
    }
    var x2: [2]F = undefined;
    x2[1] = btmp[1] / p22;
    x2[0] = btmp[0] / p11 - (p12 / p11) * x2[1];
    if (xswpiv[ipiv]) std.mem.swap(F, &x2[0], &x2[1]);
    out.x[0] = x2[0];
    if (n1 == 1) out.x[2] = x2[1] else out.x[1] = x2[1];
    return out;
}
/// `lasy2` with `n1 = n2 = 2`: `op(TL) X + sgn X op(TR) = s V` for 2 x 2
/// blocks, as a 4 x 4 system by complete pivoting.
fn lasy2Square(comptime F: type, ltranl: bool, ltranr: bool, sgn: F, tl: [*]const F, ldtl: usize, tr: [*]const F, ldtr: usize, v: [4]F) Small(F) {
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
pub fn tgsy2(al: Allocator, io: Io, comptime T: type, trans: bool, ijob: u2, m: usize, n: usize, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, c: [*]T, ldc: usize, d: [*]const T, ldd: usize, e: [*]const T, lde: usize, f: [*]T, ldf: usize, dif: *Dif(Float(T))) ops.Error!Solved(Float(T)) {
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
                    z[r + 8 * (k + ii + mb * kk)] = ops.neg(T, b[js + kk + (js + jj) * ldb]);
                    if (kk <= jj) z[k + r + 8 * (k + ii + mb * kk)] = ops.neg(T, e[js + kk + (js + jj) * lde]);
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
                        c[row + col * ldc] = ops.mulReal(T, c[row + col * ldc], scaloc);
                        f[row + col * ldf] = ops.mulReal(T, f[row + col * ldf], scaloc);
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
        for (i + 1..n) |j| z[j + 8 * i] = ops.div(T, z[j + 8 * i], z[i + 8 * i]);
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
    for (0..n) |i| if (ops.abs1(T, rhs[i]) > ops.abs1(T, rhs[big])) {
        big = i;
    };
    if (2 * smlnum * absT(T, rhs[big]) > absT(T, z[n - 1 + 8 * (n - 1)])) {
        const temp = 0.5 / absT(T, rhs[big]);
        for (rhs[0..n]) |*x| x.* = ops.mulReal(T, x.*, temp);
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
        const temp = ops.div(T, scalar(T, 1), z[i + 8 * i]);
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
fn latdf(al: Allocator, io: Io, comptime T: type, ijob: u2, n: usize, z: *const [64]T, rhs: *[8]T, dif: *Dif(Float(T)), ipiv: *const [8]usize, jpiv: *const [8]usize) ops.Error!void {
    const F = Float(T);
    const one = scalar(T, 1);
    if (ijob == 2) {
        var xm: [8]T = @splat(scalar(T, 0));
        try nullVector(al, io, T, n, z, &xm);
        unswap(T, n, &xm, ipiv);
        const temp = 1 / @sqrt(redot(T, n, &xm, &xm));
        var xp: [8]T = undefined;
        for (0..n) |k| {
            xm[k] = ops.mulReal(T, xm[k], temp);
            xp[k] = add(T, xm[k], rhs[k]);
            rhs[k] = sub(T, rhs[k], xm[k]);
        }
        _ = gesc2(T, n, z, rhs, ipiv, jpiv);
        _ = gesc2(T, n, z, &xp, ipiv, jpiv);
        var sp: F = 0;
        var sr: F = 0;
        for (0..n) |k| {
            sp += ops.abs1(T, xp[k]);
            sr += ops.abs1(T, rhs[k]);
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
        const temp = ops.neg(T, rhs[j]);
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
        const temp = ops.div(T, one, z[i + 8 * i]);
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
fn nullVector(al: Allocator, io: Io, comptime T: type, n: usize, z: *const [64]T, v: *[8]T) ops.Error!void {
    const F = Float(T);
    const Inv = struct {
        al: Allocator,
        io: Io,
        n: usize,
        z: [*]const T,
        cnorm: [16]F = undefined,
        normin: bool = false,

        const Self = @This();
        const Stop = ops.Error || error{Overflow};

        // the infinity norm's `KASE1 = 2`: `KASE = 2` is `U⁻¹ L⁻¹`
        fn apply(s: *Self, adjoint: bool, x: [*]T) Stop!void {
            const op: types.Trans = if (comptime isComplex(T)) .c else .t;
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
                const ix = blas_calls.iamax(T, s.n, x, 1);
                if (scale < ops.abs1(T, x[ix]) * std.math.floatMin(F) or scale == 0) return error.Overflow;
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

/// `?LASCL('G', 1, s)`: `C *= s`, part by part; `trsyl`'s and `tgsyl`'s.
pub fn scaleBlock(comptime T: type, m: usize, n: usize, s: Float(T), c: [*]T, ldc: usize) void {
    for (0..n) |j| for (0..m) |i| {
        c[i + j * ldc] = ops.mulReal(T, c[i + j * ldc], s);
    };
}
