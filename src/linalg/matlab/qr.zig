//! QR decomposition (with optional column pivoting) of a general complex
//! matrix, and application of the result to compute `Q*y`, `ctrans(Q)*y`,
//! least-squares solutions, residuals and projections - ported from
//! `matlab88.f90`'s WQRDC/WQRSL (LINPACK's ZQRDC/ZQRSL). See `blas1.zig`
//! for the shared matrix/vector conventions this port uses.

const std = @import("std");
const blas1 = @import("blas1.zig");
const Cx = blas1.Cx;
const ix = blas1.ix;

/// WQRDC: Householder-reduce the complex `n`x`p` matrix `xr`/`xi`
/// (column-major, `ldx >= n`) to `Q*R`, with optional column pivoting.
///
/// `jpvt` (length `p`) controls pivoting when `job != 0` (ignored when
/// `job == 0`): on entry, `jpvt[k] > 0` freezes column `k` at the front,
/// `jpvt[k] < 0` freezes it at the back, `jpvt[k] == 0` leaves it free to
/// be reordered by decreasing norm. On return, `jpvt[k]` is the original
/// (1-based) index of the column now in position `k+1`.
///
/// On return, the upper triangle of `xr`/`xi` holds `R`; below the
/// diagonal it holds information (together with `qrauxr`/`qrauxi`, length
/// `p`) needed to recover the unitary `Q` - see `wqrsl`. `workr`/`worki`
/// (length `p`) are scratch, unused when `job == 0`.
pub fn wqrdc(
    xr: []f64,
    xi: []f64,
    ldx: usize,
    n: usize,
    p: usize,
    qrauxr: []f64,
    qrauxi: []f64,
    jpvt: []isize,
    workr: []f64,
    worki: []f64,
    job: usize,
) void {
    var pl: usize = 1;
    var pu: usize = 0;

    if (job != 0) {
        // Rearrange the columns according to jpvt.
        var j: usize = 1;
        while (j <= p) : (j += 1) {
            const swapj = jpvt[j - 1] > 0;
            const negj = jpvt[j - 1] < 0;
            jpvt[j - 1] = @intCast(j);
            if (negj) jpvt[j - 1] = -@as(isize, @intCast(j));
            if (swapj) {
                if (j != pl) {
                    blas1.wswap(n, xr[ix(ldx, 1, pl)..], xi[ix(ldx, 1, pl)..], 1, xr[ix(ldx, 1, j)..], xi[ix(ldx, 1, j)..], 1);
                }
                jpvt[j - 1] = jpvt[pl - 1];
                jpvt[pl - 1] = @intCast(j);
                pl += 1;
            }
        }
        pu = p;
        var jj: usize = 1;
        while (jj <= p) : (jj += 1) {
            const j2 = p - jj + 1;
            if (jpvt[j2 - 1] < 0) {
                jpvt[j2 - 1] = -jpvt[j2 - 1];
                if (j2 != pu) {
                    blas1.wswap(n, xr[ix(ldx, 1, pu)..], xi[ix(ldx, 1, pu)..], 1, xr[ix(ldx, 1, j2)..], xi[ix(ldx, 1, j2)..], 1);
                    const jp = jpvt[pu - 1];
                    jpvt[pu - 1] = jpvt[j2 - 1];
                    jpvt[j2 - 1] = jp;
                }
                pu -= 1;
            }
        }
    }

    // Compute the norms of the free columns.
    if (pu >= pl) {
        var j: usize = pl;
        while (j <= pu) : (j += 1) {
            qrauxr[j - 1] = blas1.wnrm2(n, xr[ix(ldx, 1, j)..], xi[ix(ldx, 1, j)..], 1);
            qrauxi[j - 1] = 0.0;
            workr[j - 1] = qrauxr[j - 1];
            worki[j - 1] = qrauxi[j - 1];
        }
    }

    // Perform the Householder reduction of x.
    const lup = @min(n, p);
    var l: usize = 1;
    while (l <= lup) : (l += 1) {
        if (l >= pl and l < pu) {
            // Locate the column of largest norm and pivot it into place.
            var maxnrm: f64 = 0.0;
            var maxj = l;
            var j: usize = l;
            while (j <= pu) : (j += 1) {
                if (qrauxr[j - 1] > maxnrm) {
                    maxnrm = qrauxr[j - 1];
                    maxj = j;
                }
            }
            if (maxj != l) {
                blas1.wswap(n, xr[ix(ldx, 1, l)..], xi[ix(ldx, 1, l)..], 1, xr[ix(ldx, 1, maxj)..], xi[ix(ldx, 1, maxj)..], 1);
                qrauxr[maxj - 1] = qrauxr[l - 1];
                qrauxi[maxj - 1] = qrauxi[l - 1];
                workr[maxj - 1] = workr[l - 1];
                worki[maxj - 1] = worki[l - 1];
                const jp = jpvt[maxj - 1];
                jpvt[maxj - 1] = jpvt[l - 1];
                jpvt[l - 1] = jp;
            }
        }
        qrauxr[l - 1] = 0.0;
        qrauxi[l - 1] = 0.0;
        if (l == n) continue;

        // Compute the Householder transformation for column l.
        var nrmxlr = blas1.wnrm2(n - l + 1, xr[ix(ldx, l, l)..], xi[ix(ldx, l, l)..], 1);
        var nrmxli: f64 = 0.0;
        if (blas1.cabs1(nrmxlr, nrmxli) == 0.0) continue;

        if (blas1.cabs1(xr[ix(ldx, l, l)], xi[ix(ldx, l, l)]) != 0.0) {
            const s = blas1.wsign(Cx.init(nrmxlr, nrmxli), Cx.init(xr[ix(ldx, l, l)], xi[ix(ldx, l, l)]));
            nrmxlr = s.re;
            nrmxli = s.im;
        }
        const t0 = blas1.wdiv(Cx.init(1.0, 0.0), Cx.init(nrmxlr, nrmxli));
        blas1.wscal(n - l + 1, t0, xr[ix(ldx, l, l)..], xi[ix(ldx, l, l)..], 1);
        xr[ix(ldx, l, l)] = 1.0 + xr[ix(ldx, l, l)];

        // Apply the transformation to the remaining columns, updating norms.
        const lp1 = l + 1;
        if (p >= lp1) {
            var j: usize = lp1;
            while (j <= p) : (j += 1) {
                var t = Cx.init(
                    -blas1.wdotcr(n - l + 1, xr[ix(ldx, l, l)..], xi[ix(ldx, l, l)..], 1, xr[ix(ldx, l, j)..], xi[ix(ldx, l, j)..], 1),
                    -blas1.wdotci(n - l + 1, xr[ix(ldx, l, l)..], xi[ix(ldx, l, l)..], 1, xr[ix(ldx, l, j)..], xi[ix(ldx, l, j)..], 1),
                );
                t = blas1.wdiv(t, Cx.init(xr[ix(ldx, l, l)], xi[ix(ldx, l, l)]));
                blas1.waxpy(n - l + 1, t, xr[ix(ldx, l, l)..], xi[ix(ldx, l, l)..], 1, xr[ix(ldx, l, j)..], xi[ix(ldx, l, j)..], 1);

                if (j < pl or j > pu or blas1.cabs1(qrauxr[j - 1], qrauxi[j - 1]) == 0.0) continue;

                const nrm = blas1.pythag(xr[ix(ldx, l, j)], xi[ix(ldx, l, j)]);
                var tt = 1.0 - (nrm / qrauxr[j - 1]) * (nrm / qrauxr[j - 1]);
                tt = @max(tt, 0.0);
                const scale2 = tt;
                const stability = 1.0 + 0.05 * tt * (qrauxr[j - 1] / workr[j - 1]) * (qrauxr[j - 1] / workr[j - 1]);
                if (stability != 1.0) {
                    qrauxr[j - 1] *= @sqrt(scale2);
                    qrauxi[j - 1] *= @sqrt(scale2);
                } else {
                    qrauxr[j - 1] = blas1.wnrm2(n - l, xr[ix(ldx, l + 1, j)..], xi[ix(ldx, l + 1, j)..], 1);
                    qrauxi[j - 1] = 0.0;
                    workr[j - 1] = qrauxr[j - 1];
                    worki[j - 1] = qrauxi[j - 1];
                }
            }
        }

        // Save the transformation.
        qrauxr[l - 1] = xr[ix(ldx, l, l)];
        qrauxi[l - 1] = xi[ix(ldx, l, l)];
        xr[ix(ldx, l, l)] = -nrmxlr;
        xi[ix(ldx, l, l)] = -nrmxli;
    }
}

/// WQRSL: apply `wqrdc`'s output to compute coordinate transformations,
/// projections, and least-squares solutions. For `k <= min(n, p)`, let
/// `xk` be the matrix formed from columns `jpvt[0]..jpvt[k-1]` of the
/// original matrix passed to `wqrdc` (or its first `k` columns, if no
/// pivoting was done); `wqrdc` produced a unitary `Q` and upper triangular
/// `R` with `xk = Q * [R; 0]`, encoded in `xr`/`xi`/`qrauxr`/`qrauxi`.
///
/// `job` selects what to compute, as a decimal `abcde`: `a != 0` computes
/// `Q*y` into `qyr`/`qyi`; any of `b,c,d,e != 0` computes `ctrans(Q)*y`
/// into `qtyr`/`qtyi` (needed internally by `c`/`d`/`e` even if not
/// returned); `c != 0` computes the least-squares solution `b` (minimizing
/// `norm(y - xk*b)`) into `br`/`bi`; `d != 0` computes the residual
/// `y - xk*b` into `rsdr`/`rsdi`; `e != 0` computes the projection `xk*b`
/// into `xbr`/`xbi`. Passing the same array for `y`/`qty`/one of
/// `b`,`rsd`,`xb` (and dummies for the rest) is fine, matching the
/// original's storage-sharing conventions - only the last-computed value
/// for a shared array survives.
///
/// Returns `info`: 0, unless `b` was requested and `R` is exactly
/// singular, in which case it's the 1-based index of the first zero
/// diagonal element of `R` (and `b` is left unaltered).
pub fn wqrsl(
    // Not const: the algorithm temporarily overwrites each diagonal entry
    // `x(j,j)` with `qraux(j)` (the Householder vector's leading element,
    // which wqrdc replaced `x(j,j)` with `-nrmxl` in storage) to apply the
    // reflector, then restores the original value before moving on - a net
    // no-op from the caller's perspective, but a real write while it lasts.
    xr: []f64,
    xi: []f64,
    ldx: usize,
    n: usize,
    k: usize,
    qrauxr: []const f64,
    qrauxi: []const f64,
    yr: []const f64,
    yi: []const f64,
    qyr: []f64,
    qyi: []f64,
    qtyr: []f64,
    qtyi: []f64,
    br: []f64,
    bi: []f64,
    rsdr: []f64,
    rsdi: []f64,
    xbr: []f64,
    xbi: []f64,
    job: usize,
) usize {
    var info: usize = 0;

    const cqy = job / 10000 != 0;
    const cqty = job % 10000 != 0;
    const cb = (job % 1000) / 100 != 0;
    const cr = (job % 100) / 10 != 0;
    const cxb = job % 10 != 0;
    const ju = @min(k, n -| 1);

    if (ju == 0) {
        // Special action when n == 1.
        if (cqy) {
            qyr[0] = yr[0];
            qyi[0] = yi[0];
        }
        if (cqty) {
            qtyr[0] = yr[0];
            qtyi[0] = yi[0];
        }
        if (cxb) {
            xbr[0] = yr[0];
            xbi[0] = yi[0];
        }
        if (cb) {
            if (blas1.cabs1(xr[ix(ldx, 1, 1)], xi[ix(ldx, 1, 1)]) != 0.0) {
                const d = blas1.wdiv(Cx.init(yr[0], yi[0]), Cx.init(xr[ix(ldx, 1, 1)], xi[ix(ldx, 1, 1)]));
                br[0] = d.re;
                bi[0] = d.im;
            } else {
                info = 1;
            }
        }
        if (cr) {
            rsdr[0] = 0.0;
            rsdi[0] = 0.0;
        }
        return info;
    }

    // Set up to compute qy or qty.
    if (cqy) blas1.wcopy(n, yr, yi, 1, qyr, qyi, 1);
    if (cqty) blas1.wcopy(n, yr, yi, 1, qtyr, qtyi, 1);

    if (cqy) {
        // Compute q*y.
        var jj: usize = 1;
        while (jj <= ju) : (jj += 1) {
            const j = ju - jj + 1;
            if (blas1.cabs1(qrauxr[j - 1], qrauxi[j - 1]) == 0.0) continue;
            const temp = Cx.init(xr[ix(ldx, j, j)], xi[ix(ldx, j, j)]);
            xr[ix(ldx, j, j)] = qrauxr[j - 1];
            xi[ix(ldx, j, j)] = qrauxi[j - 1];
            var t = Cx.init(
                -blas1.wdotcr(n - j + 1, xr[ix(ldx, j, j)..], xi[ix(ldx, j, j)..], 1, qyr[j - 1 ..], qyi[j - 1 ..], 1),
                -blas1.wdotci(n - j + 1, xr[ix(ldx, j, j)..], xi[ix(ldx, j, j)..], 1, qyr[j - 1 ..], qyi[j - 1 ..], 1),
            );
            t = blas1.wdiv(t, Cx.init(xr[ix(ldx, j, j)], xi[ix(ldx, j, j)]));
            blas1.waxpy(n - j + 1, t, xr[ix(ldx, j, j)..], xi[ix(ldx, j, j)..], 1, qyr[j - 1 ..], qyi[j - 1 ..], 1);
            xr[ix(ldx, j, j)] = temp.re;
            xi[ix(ldx, j, j)] = temp.im;
        }
    }

    if (cqty) {
        // Compute ctrans(q)*y.
        var j: usize = 1;
        while (j <= ju) : (j += 1) {
            if (blas1.cabs1(qrauxr[j - 1], qrauxi[j - 1]) == 0.0) continue;
            const temp = Cx.init(xr[ix(ldx, j, j)], xi[ix(ldx, j, j)]);
            xr[ix(ldx, j, j)] = qrauxr[j - 1];
            xi[ix(ldx, j, j)] = qrauxi[j - 1];
            var t = Cx.init(
                -blas1.wdotcr(n - j + 1, xr[ix(ldx, j, j)..], xi[ix(ldx, j, j)..], 1, qtyr[j - 1 ..], qtyi[j - 1 ..], 1),
                -blas1.wdotci(n - j + 1, xr[ix(ldx, j, j)..], xi[ix(ldx, j, j)..], 1, qtyr[j - 1 ..], qtyi[j - 1 ..], 1),
            );
            t = blas1.wdiv(t, Cx.init(xr[ix(ldx, j, j)], xi[ix(ldx, j, j)]));
            blas1.waxpy(n - j + 1, t, xr[ix(ldx, j, j)..], xi[ix(ldx, j, j)..], 1, qtyr[j - 1 ..], qtyi[j - 1 ..], 1);
            xr[ix(ldx, j, j)] = temp.re;
            xi[ix(ldx, j, j)] = temp.im;
        }
    }

    // Set up to compute b, rsd, or xb.
    if (cb) blas1.wcopy(k, qtyr, qtyi, 1, br, bi, 1);
    const kp1 = k + 1;
    if (cxb) blas1.wcopy(k, qtyr, qtyi, 1, xbr, xbi, 1);
    if (cr and k < n) blas1.wcopy(n - k, qtyr[kp1 - 1 ..], qtyi[kp1 - 1 ..], 1, rsdr[kp1 - 1 ..], rsdi[kp1 - 1 ..], 1);
    if (cxb and kp1 <= n) {
        var i: usize = kp1;
        while (i <= n) : (i += 1) {
            xbr[i - 1] = 0.0;
            xbi[i - 1] = 0.0;
        }
    }
    if (cr) {
        var i: usize = 1;
        while (i <= k) : (i += 1) {
            rsdr[i - 1] = 0.0;
            rsdi[i - 1] = 0.0;
        }
    }

    if (cb) {
        // Compute b.
        var jj: usize = 1;
        b_loop: while (jj <= k) : (jj += 1) {
            const j = k - jj + 1;
            if (blas1.cabs1(xr[ix(ldx, j, j)], xi[ix(ldx, j, j)]) == 0.0) {
                info = j;
                break :b_loop;
            }
            const d = blas1.wdiv(Cx.init(br[j - 1], bi[j - 1]), Cx.init(xr[ix(ldx, j, j)], xi[ix(ldx, j, j)]));
            br[j - 1] = d.re;
            bi[j - 1] = d.im;
            if (j != 1) {
                blas1.waxpy(j - 1, Cx.init(-d.re, -d.im), xr[ix(ldx, 1, j)..], xi[ix(ldx, 1, j)..], 1, br, bi, 1);
            }
        }
    }

    if (cr or cxb) {
        // Compute rsd and/or xb as required.
        var jj: usize = 1;
        while (jj <= ju) : (jj += 1) {
            const j = ju - jj + 1;
            if (blas1.cabs1(qrauxr[j - 1], qrauxi[j - 1]) == 0.0) continue;
            const temp = Cx.init(xr[ix(ldx, j, j)], xi[ix(ldx, j, j)]);
            xr[ix(ldx, j, j)] = qrauxr[j - 1];
            xi[ix(ldx, j, j)] = qrauxi[j - 1];
            if (cr) {
                var t = Cx.init(
                    -blas1.wdotcr(n - j + 1, xr[ix(ldx, j, j)..], xi[ix(ldx, j, j)..], 1, rsdr[j - 1 ..], rsdi[j - 1 ..], 1),
                    -blas1.wdotci(n - j + 1, xr[ix(ldx, j, j)..], xi[ix(ldx, j, j)..], 1, rsdr[j - 1 ..], rsdi[j - 1 ..], 1),
                );
                t = blas1.wdiv(t, Cx.init(xr[ix(ldx, j, j)], xi[ix(ldx, j, j)]));
                blas1.waxpy(n - j + 1, t, xr[ix(ldx, j, j)..], xi[ix(ldx, j, j)..], 1, rsdr[j - 1 ..], rsdi[j - 1 ..], 1);
            }
            if (cxb) {
                var t = Cx.init(
                    -blas1.wdotcr(n - j + 1, xr[ix(ldx, j, j)..], xi[ix(ldx, j, j)..], 1, xbr[j - 1 ..], xbi[j - 1 ..], 1),
                    -blas1.wdotci(n - j + 1, xr[ix(ldx, j, j)..], xi[ix(ldx, j, j)..], 1, xbr[j - 1 ..], xbi[j - 1 ..], 1),
                );
                t = blas1.wdiv(t, Cx.init(xr[ix(ldx, j, j)], xi[ix(ldx, j, j)]));
                blas1.waxpy(n - j + 1, t, xr[ix(ldx, j, j)..], xi[ix(ldx, j, j)..], 1, xbr[j - 1 ..], xbi[j - 1 ..], 1);
            }
            xr[ix(ldx, j, j)] = temp.re;
            xi[ix(ldx, j, j)] = temp.im;
        }
    }

    return info;
}

const testing = std.testing;

test "wqrdc + wqrsl reconstruct a real 3x2 matrix as Q*R with orthonormal Q" {
    // X = [[1,1],[0,1],[1,0]] (column-major, no pivoting).
    const x_orig = [_]f64{ 1.0, 0.0, 1.0, 1.0, 1.0, 0.0 };
    var xr = x_orig;
    var xi: [6]f64 = @splat(0.0);
    var qrauxr: [2]f64 = @splat(0.0);
    var qrauxi: [2]f64 = @splat(0.0);
    var jpvt: [2]isize = @splat(0);
    var workr: [2]f64 = @splat(0.0);
    var worki: [2]f64 = @splat(0.0);

    wqrdc(&xr, &xi, 3, 3, 2, &qrauxr, &qrauxi, &jpvt, &workr, &worki, 0);

    // Extract R (2x2 upper triangular) before x's storage is reused below.
    var r = [_]f64{ 0.0, 0.0, 0.0, 0.0 }; // column-major 2x2
    r[ix(2, 1, 1)] = xr[ix(3, 1, 1)];
    r[ix(2, 1, 2)] = xr[ix(3, 1, 2)];
    r[ix(2, 2, 2)] = xr[ix(3, 2, 2)];

    // Recover Q (3x3) by applying wqrsl's Q*y to each standard basis vector.
    var q: [9]f64 = @splat(0.0); // column-major 3x3
    for (0..3) |col| {
        var e = [_]f64{ 0.0, 0.0, 0.0 };
        e[col] = 1.0;
        var ei_: [3]f64 = @splat(0.0);
        var qy: [3]f64 = @splat(0.0);
        var qyi: [3]f64 = @splat(0.0);
        var dummy: [3]f64 = @splat(0.0);
        var dummyi: [3]f64 = @splat(0.0);
        const info = wqrsl(
            &xr,
            &xi,
            3,
            3,
            2,
            &qrauxr,
            &qrauxi,
            &e,
            &ei_,
            &qy,
            &qyi,
            &dummy,
            &dummyi,
            &dummy,
            &dummyi,
            &dummy,
            &dummyi,
            &dummy,
            &dummyi,
            10000,
        );
        try testing.expectEqual(@as(usize, 0), info);
        for (0..3) |row| q[ix(3, row + 1, col + 1)] = qy[row];
    }

    // Q columns should be orthonormal.
    for (0..3) |a| {
        for (0..3) |b| {
            var dot: f64 = 0.0;
            for (0..3) |row| dot += q[ix(3, row + 1, a + 1)] * q[ix(3, row + 1, b + 1)];
            const expected: f64 = if (a == b) 1.0 else 0.0;
            try testing.expectApproxEqAbs(expected, dot, 1e-9);
        }
    }

    // Q * [R; 0] should reconstruct the original X.
    for (0..3) |row| {
        for (0..2) |col| {
            var sum: f64 = 0.0;
            for (0..2) |k| sum += q[ix(3, row + 1, k + 1)] * r[ix(2, k + 1, col + 1)];
            try testing.expectApproxEqAbs(x_orig[ix(3, row + 1, col + 1)], sum, 1e-9);
        }
    }
}

test "wqrsl solves an exact least-squares fit y = 1 + 2x" {
    // Design matrix [1, x] for x = 0,1,2; y = 1,3,5 lies exactly on the
    // line, so the least-squares solution should be b = [1, 2] with a
    // (near) zero residual.
    var xr = [_]f64{ 1.0, 1.0, 1.0, 0.0, 1.0, 2.0 }; // column-major 3x2
    var xi: [6]f64 = @splat(0.0);
    var qrauxr: [2]f64 = @splat(0.0);
    var qrauxi: [2]f64 = @splat(0.0);
    var jpvt: [2]isize = @splat(0);
    var workr: [2]f64 = @splat(0.0);
    var worki: [2]f64 = @splat(0.0);
    wqrdc(&xr, &xi, 3, 3, 2, &qrauxr, &qrauxi, &jpvt, &workr, &worki, 0);

    var yr = [_]f64{ 1.0, 3.0, 5.0 };
    var yi: [3]f64 = @splat(0.0);
    var qty: [3]f64 = @splat(0.0);
    var qtyi: [3]f64 = @splat(0.0);
    var b: [2]f64 = @splat(0.0);
    var bi: [2]f64 = @splat(0.0);
    var rsd: [3]f64 = @splat(0.0);
    var rsdi: [3]f64 = @splat(0.0);
    var dummy: [3]f64 = @splat(0.0);
    var dummyi: [3]f64 = @splat(0.0);

    const info = wqrsl(
        &xr,
        &xi,
        3,
        3,
        2,
        &qrauxr,
        &qrauxi,
        &yr,
        &yi,
        &dummy,
        &dummyi,
        &qty,
        &qtyi,
        &b,
        &bi,
        &rsd,
        &rsdi,
        &dummy,
        &dummyi,
        110,
    );
    try testing.expectEqual(@as(usize, 0), info);
    try testing.expectApproxEqAbs(1.0, b[0], 1e-9);
    try testing.expectApproxEqAbs(2.0, b[1], 1e-9);
    for (rsd) |r_| try testing.expectApproxEqAbs(0.0, r_, 1e-8);
}
