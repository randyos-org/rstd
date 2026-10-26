//! Singular value decomposition of a general complex matrix, ported from
//! `matlab88.f90`'s WSVDC (LINPACK's ZSVDC). See `blas1.zig` for the
//! shared matrix/vector conventions this port uses.

const std = @import("std");
const blas1 = @import("blas1.zig");
const CxMat = @import("CxMat.zig");
const Cx = blas1.Cx;
const ix = blas1.ix;

/// WSVDC: reduce the complex `n`x`p` matrix `x` (`n = x.m`, `p = x.n`) by
/// unitary transformations `U` and `V` to diagonal form. The diagonal
/// entries of `S` are the singular values of `X`; the columns of `U` are
/// the corresponding left singular vectors, and the columns of `V` the
/// right singular vectors. `x` is destroyed.
///
/// `sr`/`si` (length `min(n+1, p)`) receive the singular values in
/// descending order (only the first `min(n, p)` are meaningful); `er`/`ei`
/// (length `p`) are scratch, ordinarily left zero; `workr`/`worki` (length
/// `n`) are scratch.
///
/// `job` selects what to compute, as a two-digit decimal `AB`:
/// - `A == 0`: don't compute left singular vectors.
/// - `A == 1`: return all `n` left singular vectors in `u`
///   (`u.lda >= n`, `min(n, p)` or more columns).
/// - `A >= 2`: return only the first `min(n, p)` left singular vectors.
/// - `B == 0`: don't compute right singular vectors.
/// - `B == 1`: return the right singular vectors in `v` (`v.lda >= p`).
///
/// Returns `info`: the singular values (and vectors)
/// `s[info]..s[min(n,p)-1]` (0-based) are correct; `info == 0` means all of
/// them are. In any case, `ctrans(U)*X*V` is the bidiagonal matrix with
/// `S` on its diagonal and `E` on its superdiagonal, so `X` and that
/// bidiagonal matrix always share singular values.
pub fn wsvdc(
    x: CxMat,
    sr: []f64,
    si: []f64,
    er: []f64,
    ei: []f64,
    u: CxMat,
    v: CxMat,
    workr: []f64,
    worki: []f64,
    job: usize,
) usize {
    const xr = x.r;
    const xi = x.i;
    const ldx = x.lda;
    const n = x.m;
    const p = x.n;
    const ur = u.r;
    const ui = u.i;
    const ldu = u.lda;
    const vr = v.r;
    const vi = v.i;
    const ldv = v.lda;

    const maxit = 75;
    const small = 1.0 / std.math.pow(f64, 2.0, 48.0);

    const jobu = (job % 100) / 10;
    var ncu = n;
    if (jobu > 1) ncu = @min(n, p);
    const wantu = jobu != 0;
    const wantv = job % 10 != 0;

    var info: usize = 0;
    const nct = @min(n -| 1, p);
    const nrt: usize = if (p >= 2) @min(p - 2, n) else 0;
    const lu = @max(nct, nrt);

    // Reduce x to bidiagonal form, storing the diagonal in s and the
    // super-diagonal in e.
    if (lu >= 1) {
        var l: usize = 1;
        while (l <= lu) : (l += 1) {
            const lp1 = l + 1;
            if (l <= nct) {
                sr[l - 1] = blas1.wnrm2(n - l + 1, xr[ix(ldx, l, l)..], xi[ix(ldx, l, l)..], 1);
                si[l - 1] = 0.0;
                if (blas1.cabs1(sr[l - 1], si[l - 1]) != 0.0) {
                    if (blas1.cabs1(xr[ix(ldx, l, l)], xi[ix(ldx, l, l)]) != 0.0) {
                        const s = blas1.wsign(Cx.init(sr[l - 1], si[l - 1]), Cx.init(xr[ix(ldx, l, l)], xi[ix(ldx, l, l)]));
                        sr[l - 1] = s.re;
                        si[l - 1] = s.im;
                    }
                    const t = blas1.wdiv(Cx.init(1.0, 0.0), Cx.init(sr[l - 1], si[l - 1]));
                    blas1.wscal(n - l + 1, t, xr[ix(ldx, l, l)..], xi[ix(ldx, l, l)..], 1);
                    xr[ix(ldx, l, l)] = 1.0 + xr[ix(ldx, l, l)];
                }
                sr[l - 1] = -sr[l - 1];
                si[l - 1] = -si[l - 1];
            }

            if (p >= lp1) {
                var j: usize = lp1;
                while (j <= p) : (j += 1) {
                    if (l <= nct and blas1.cabs1(sr[l - 1], si[l - 1]) != 0.0) {
                        var t = Cx.init(
                            -blas1.wdotcr(n - l + 1, xr[ix(ldx, l, l)..], xi[ix(ldx, l, l)..], 1, xr[ix(ldx, l, j)..], xi[ix(ldx, l, j)..], 1),
                            -blas1.wdotci(n - l + 1, xr[ix(ldx, l, l)..], xi[ix(ldx, l, l)..], 1, xr[ix(ldx, l, j)..], xi[ix(ldx, l, j)..], 1),
                        );
                        t = blas1.wdiv(t, Cx.init(xr[ix(ldx, l, l)], xi[ix(ldx, l, l)]));
                        blas1.waxpy(n - l + 1, t, xr[ix(ldx, l, l)..], xi[ix(ldx, l, l)..], 1, xr[ix(ldx, l, j)..], xi[ix(ldx, l, j)..], 1);
                    }
                    er[j - 1] = xr[ix(ldx, l, j)];
                    ei[j - 1] = -xi[ix(ldx, l, j)];
                }
            }

            if (wantu and l <= nct) {
                var i: usize = l;
                while (i <= n) : (i += 1) {
                    ur[ix(ldu, i, l)] = xr[ix(ldx, i, l)];
                    ui[ix(ldu, i, l)] = xi[ix(ldx, i, l)];
                }
            }

            if (l <= nrt) {
                er[l - 1] = blas1.wnrm2(p - l, er[lp1 - 1 ..], ei[lp1 - 1 ..], 1);
                ei[l - 1] = 0.0;
                if (blas1.cabs1(er[l - 1], ei[l - 1]) != 0.0) {
                    if (blas1.cabs1(er[lp1 - 1], ei[lp1 - 1]) != 0.0) {
                        const s = blas1.wsign(Cx.init(er[l - 1], ei[l - 1]), Cx.init(er[lp1 - 1], ei[lp1 - 1]));
                        er[l - 1] = s.re;
                        ei[l - 1] = s.im;
                    }
                    const t = blas1.wdiv(Cx.init(1.0, 0.0), Cx.init(er[l - 1], ei[l - 1]));
                    blas1.wscal(p - l, t, er[lp1 - 1 ..], ei[lp1 - 1 ..], 1);
                    er[lp1 - 1] = 1.0 + er[lp1 - 1];
                }
                er[l - 1] = -er[l - 1];

                if (!(lp1 > n or blas1.cabs1(er[l - 1], ei[l - 1]) == 0.0)) {
                    var i: usize = lp1;
                    while (i <= n) : (i += 1) {
                        workr[i - 1] = 0.0;
                        worki[i - 1] = 0.0;
                    }
                    var j: usize = lp1;
                    while (j <= p) : (j += 1) {
                        blas1.waxpy(n - l, Cx.init(er[j - 1], ei[j - 1]), xr[ix(ldx, lp1, j)..], xi[ix(ldx, lp1, j)..], 1, workr[lp1 - 1 ..], worki[lp1 - 1 ..], 1);
                    }
                    j = lp1;
                    while (j <= p) : (j += 1) {
                        const t = blas1.wdiv(Cx.init(-er[j - 1], -ei[j - 1]), Cx.init(er[lp1 - 1], ei[lp1 - 1]));
                        blas1.waxpy(n - l, Cx.init(t.re, -t.im), workr[lp1 - 1 ..], worki[lp1 - 1 ..], 1, xr[ix(ldx, lp1, j)..], xi[ix(ldx, lp1, j)..], 1);
                    }
                }
                if (wantv) {
                    var i: usize = lp1;
                    while (i <= p) : (i += 1) {
                        vr[ix(ldv, i, l)] = er[i - 1];
                        vi[ix(ldv, i, l)] = ei[i - 1];
                    }
                }
            }
        }
    }

    // Set up the final bidiagonal matrix of order m.
    const m = @min(p, n + 1);
    const nctp1 = nct + 1;
    const nrtp1 = nrt + 1;
    if (nct < p) {
        sr[nctp1 - 1] = xr[ix(ldx, nctp1, nctp1)];
        si[nctp1 - 1] = xi[ix(ldx, nctp1, nctp1)];
    }
    if (n < m) {
        sr[m - 1] = 0.0;
        si[m - 1] = 0.0;
    }
    if (nrtp1 < m) {
        er[nrtp1 - 1] = xr[ix(ldx, nrtp1, m)];
        ei[nrtp1 - 1] = xi[ix(ldx, nrtp1, m)];
    }
    er[m - 1] = 0.0;
    ei[m - 1] = 0.0;

    // If required, generate u.
    if (wantu) {
        if (ncu >= nctp1) {
            var j: usize = nctp1;
            while (j <= ncu) : (j += 1) {
                var i: usize = 1;
                while (i <= n) : (i += 1) {
                    ur[ix(ldu, i, j)] = 0.0;
                    ui[ix(ldu, i, j)] = 0.0;
                }
                ur[ix(ldu, j, j)] = 1.0;
                ui[ix(ldu, j, j)] = 0.0;
            }
        }
        if (nct >= 1) {
            var ll: usize = 1;
            while (ll <= nct) : (ll += 1) {
                const l = nct - ll + 1;
                if (blas1.cabs1(sr[l - 1], si[l - 1]) != 0.0) {
                    const lp1 = l + 1;
                    if (ncu >= lp1) {
                        var j: usize = lp1;
                        while (j <= ncu) : (j += 1) {
                            var t = Cx.init(
                                -blas1.wdotcr(n - l + 1, ur[ix(ldu, l, l)..], ui[ix(ldu, l, l)..], 1, ur[ix(ldu, l, j)..], ui[ix(ldu, l, j)..], 1),
                                -blas1.wdotci(n - l + 1, ur[ix(ldu, l, l)..], ui[ix(ldu, l, l)..], 1, ur[ix(ldu, l, j)..], ui[ix(ldu, l, j)..], 1),
                            );
                            t = blas1.wdiv(t, Cx.init(ur[ix(ldu, l, l)], ui[ix(ldu, l, l)]));
                            blas1.waxpy(n - l + 1, t, ur[ix(ldu, l, l)..], ui[ix(ldu, l, l)..], 1, ur[ix(ldu, l, j)..], ui[ix(ldu, l, j)..], 1);
                        }
                    }
                    blas1.wrscal(n - l + 1, -1.0, ur[ix(ldu, l, l)..], ui[ix(ldu, l, l)..], 1);
                    ur[ix(ldu, l, l)] = 1.0 + ur[ix(ldu, l, l)];
                    if (l >= 2) {
                        var i: usize = 1;
                        while (i <= l - 1) : (i += 1) {
                            ur[ix(ldu, i, l)] = 0.0;
                            ui[ix(ldu, i, l)] = 0.0;
                        }
                    }
                } else {
                    var i: usize = 1;
                    while (i <= n) : (i += 1) {
                        ur[ix(ldu, i, l)] = 0.0;
                        ui[ix(ldu, i, l)] = 0.0;
                    }
                    ur[ix(ldu, l, l)] = 1.0;
                    ui[ix(ldu, l, l)] = 0.0;
                }
            }
        }
    }

    // If required, generate v.
    if (wantv) {
        var ll: usize = 1;
        while (ll <= p) : (ll += 1) {
            const l = p - ll + 1;
            const lp1 = l + 1;
            if (l <= nrt and blas1.cabs1(er[l - 1], ei[l - 1]) != 0.0) {
                var j: usize = lp1;
                while (j <= p) : (j += 1) {
                    var t = Cx.init(
                        -blas1.wdotcr(p - l, vr[ix(ldv, lp1, l)..], vi[ix(ldv, lp1, l)..], 1, vr[ix(ldv, lp1, j)..], vi[ix(ldv, lp1, j)..], 1),
                        -blas1.wdotci(p - l, vr[ix(ldv, lp1, l)..], vi[ix(ldv, lp1, l)..], 1, vr[ix(ldv, lp1, j)..], vi[ix(ldv, lp1, j)..], 1),
                    );
                    t = blas1.wdiv(t, Cx.init(vr[ix(ldv, lp1, l)], vi[ix(ldv, lp1, l)]));
                    blas1.waxpy(p - l, t, vr[ix(ldv, lp1, l)..], vi[ix(ldv, lp1, l)..], 1, vr[ix(ldv, lp1, j)..], vi[ix(ldv, lp1, j)..], 1);
                }
            }
            var i: usize = 1;
            while (i <= p) : (i += 1) {
                vr[ix(ldv, i, l)] = 0.0;
                vi[ix(ldv, i, l)] = 0.0;
            }
            vr[ix(ldv, l, l)] = 1.0;
            vi[ix(ldv, l, l)] = 0.0;
        }
    }

    // Transform s and e so that they are real.
    {
        var i: usize = 1;
        transform_loop: while (i <= m) : (i += 1) {
            var tr = blas1.pythag(sr[i - 1], si[i - 1]);
            if (tr != 0.0) {
                const rr = sr[i - 1] / tr;
                const ri = si[i - 1] / tr;
                sr[i - 1] = tr;
                si[i - 1] = 0.0;
                if (i < m) {
                    const d = blas1.wdiv(Cx.init(er[i - 1], ei[i - 1]), Cx.init(rr, ri));
                    er[i - 1] = d.re;
                    ei[i - 1] = d.im;
                }
                if (wantu) blas1.wscal(n, Cx.init(rr, ri), ur[ix(ldu, 1, i)..], ui[ix(ldu, 1, i)..], 1);
            }
            if (i == m) break :transform_loop;
            tr = blas1.pythag(er[i - 1], ei[i - 1]);
            if (tr != 0.0) {
                const d = blas1.wdiv(Cx.init(tr, 0.0), Cx.init(er[i - 1], ei[i - 1]));
                er[i - 1] = tr;
                ei[i - 1] = 0.0;
                const sm = blas1.wmul(Cx.init(sr[i], si[i]), d);
                sr[i] = sm.re;
                si[i] = sm.im;
                if (wantv) blas1.wscal(p, d, vr[ix(ldv, 1, i + 1)..], vi[ix(ldv, 1, i + 1)..], 1);
            }
        }
    }

    // Main iteration loop for the singular values.
    const mm = m;
    var m_var = m;
    var iter: usize = 0;
    outer_iter: while (true) {
        if (m_var == 0) break :outer_iter;
        if (iter >= maxit) {
            info = m_var;
            break :outer_iter;
        }

        // Inspect for negligible elements in s and e.
        var l: usize = 0;
        find_l: {
            var ll: usize = 1;
            while (ll <= m_var) : (ll += 1) {
                l = m_var - ll;
                if (l == 0) break :find_l;
                const negl = @abs(sr[l - 1]) + @abs(sr[l]);
                const ztest = negl + @abs(er[l - 1]) / 2.0;
                if (small * ztest == small * negl) {
                    er[l - 1] = 0.0;
                    break :find_l;
                }
            }
        }

        var kase: u8 = undefined;
        if (l == m_var - 1) {
            kase = 4;
        } else {
            var ls: usize = l;
            find_ls: {
                const lp1 = l + 1;
                const mp1 = m_var + 1;
                var llls: usize = lp1;
                while (llls <= mp1) : (llls += 1) {
                    ls = (lp1 + m_var) - llls; // reordered vs. the Fortran m_var-llls+lp1 to avoid unsigned underflow
                    if (ls == l) break :find_ls;
                    var negl: f64 = 0.0;
                    if (ls != m_var) negl += @abs(er[ls - 1]);
                    if (ls != l + 1) negl += @abs(er[ls - 2]);
                    const ztest = negl + @abs(sr[ls - 1]) / 2.0;
                    if (small * ztest == small * negl) {
                        sr[ls - 1] = 0.0;
                        break :find_ls;
                    }
                }
            }
            if (ls == l) {
                kase = 3;
            } else if (ls == m_var) {
                kase = 1;
            } else {
                kase = 2;
                l = ls;
            }
        }
        l += 1;

        switch (kase) {
            1 => {
                // Deflate negligible sr[m_var - 1].
                const mm1 = m_var - 1;
                var f = er[mm1 - 1];
                er[mm1 - 1] = 0.0;
                var kk: usize = l;
                while (kk <= mm1) : (kk += 1) {
                    const k = mm1 - kk + l;
                    var t1 = sr[k - 1];
                    const g = blas1.rrotg(t1, f);
                    t1 = g.r;
                    sr[k - 1] = t1;
                    if (k != l) {
                        f = -(g.s * er[k - 2]);
                        er[k - 2] = g.c * er[k - 2];
                    }
                    if (wantv) {
                        blas1.rrot(p, vr[ix(ldv, 1, k)..], 1, vr[ix(ldv, 1, m_var)..], 1, g.c, g.s);
                        blas1.rrot(p, vi[ix(ldv, 1, k)..], 1, vi[ix(ldv, 1, m_var)..], 1, g.c, g.s);
                    }
                }
            },
            2 => {
                // Split at negligible sr[l - 1].
                var f = er[l - 2];
                er[l - 2] = 0.0;
                var k: usize = l;
                while (k <= m_var) : (k += 1) {
                    var t1 = sr[k - 1];
                    const g = blas1.rrotg(t1, f);
                    t1 = g.r;
                    sr[k - 1] = t1;
                    f = -(g.s * er[k - 1]);
                    er[k - 1] = g.c * er[k - 1];
                    if (wantu) {
                        blas1.rrot(n, ur[ix(ldu, 1, k)..], 1, ur[ix(ldu, 1, l - 1)..], 1, g.c, g.s);
                        blas1.rrot(n, ui[ix(ldu, 1, k)..], 1, ui[ix(ldu, 1, l - 1)..], 1, g.c, g.s);
                    }
                }
            },
            3 => {
                // One QR step.
                const scale = @max(@max(@max(@max(@abs(sr[m_var - 1]), @abs(sr[m_var - 2])), @abs(er[m_var - 2])), @abs(sr[l - 1])), @abs(er[l - 1]));
                const sm = sr[m_var - 1] / scale;
                const smm1 = sr[m_var - 2] / scale;
                const emm1 = er[m_var - 2] / scale;
                const sl = sr[l - 1] / scale;
                const el = er[l - 1] / scale;
                const b = ((smm1 + sm) * (smm1 - sm) + emm1 * emm1) / 2.0;
                const c = (sm * emm1) * (sm * emm1);
                var shift: f64 = 0.0;
                if (b != 0.0 or c != 0.0) {
                    shift = @sqrt(b * b + c);
                    if (b < 0.0) shift = -shift;
                    shift = c / (b + shift);
                }
                var f = (sl + sm) * (sl - sm) - shift;
                var g = sl * el;
                const mm1 = m_var - 1;
                var k: usize = l;
                while (k <= mm1) : (k += 1) {
                    const g1 = blas1.rrotg(f, g);
                    if (k != l) er[k - 2] = g1.r;
                    f = g1.c * sr[k - 1] + g1.s * er[k - 1];
                    er[k - 1] = g1.c * er[k - 1] - g1.s * sr[k - 1];
                    g = g1.s * sr[k];
                    sr[k] = g1.c * sr[k];
                    if (wantv) {
                        blas1.rrot(p, vr[ix(ldv, 1, k)..], 1, vr[ix(ldv, 1, k + 1)..], 1, g1.c, g1.s);
                        blas1.rrot(p, vi[ix(ldv, 1, k)..], 1, vi[ix(ldv, 1, k + 1)..], 1, g1.c, g1.s);
                    }
                    const g2 = blas1.rrotg(f, g);
                    sr[k - 1] = g2.r;
                    f = g2.c * er[k - 1] + g2.s * sr[k];
                    sr[k] = -(g2.s * er[k - 1]) + g2.c * sr[k];
                    g = g2.s * er[k];
                    er[k] = g2.c * er[k];
                    if (wantu and k < n) {
                        blas1.rrot(n, ur[ix(ldu, 1, k)..], 1, ur[ix(ldu, 1, k + 1)..], 1, g2.c, g2.s);
                        blas1.rrot(n, ui[ix(ldu, 1, k)..], 1, ui[ix(ldu, 1, k + 1)..], 1, g2.c, g2.s);
                    }
                }
                er[m_var - 2] = f;
                iter += 1;
            },
            4 => {
                // Convergence.
                if (sr[l - 1] < 0.0) {
                    sr[l - 1] = -sr[l - 1];
                    if (wantv) blas1.wrscal(p, -1.0, vr[ix(ldv, 1, l)..], vi[ix(ldv, 1, l)..], 1);
                }
                var lcur = l;
                while (lcur != mm and sr[lcur - 1] < sr[lcur]) {
                    std.mem.swap(f64, &sr[lcur - 1], &sr[lcur]);
                    if (wantv and lcur < p) blas1.wswap(p, vr[ix(ldv, 1, lcur)..], vi[ix(ldv, 1, lcur)..], 1, vr[ix(ldv, 1, lcur + 1)..], vi[ix(ldv, 1, lcur + 1)..], 1);
                    if (wantu and lcur < n) blas1.wswap(n, ur[ix(ldu, 1, lcur)..], ui[ix(ldu, 1, lcur)..], 1, ur[ix(ldu, 1, lcur + 1)..], ui[ix(ldu, 1, lcur + 1)..], 1);
                    lcur += 1;
                }
                iter = 0;
                m_var -= 1;
            },
            else => unreachable,
        }
    }

    return info;
}

const zigstd = @import("../../zigstd.zig");

test "wsvdc reconstructs a known real 2x2 matrix via U*S*ctrans(V)" {
    // X = [[1,1],[0,1]] (the classic shear matrix; singular values
    // (1+sqrt5)/2 and (sqrt5-1)/2).
    const x = [_]f64{ 1.0, 0.0, 1.0, 1.0 }; // column-major, preserved for the check
    var xr = x;
    var xi: [4]f64 = @splat(0.0);
    var sr: [2]f64 = @splat(0.0);
    var si: [2]f64 = @splat(0.0);
    var er: [2]f64 = @splat(0.0);
    var ei: [2]f64 = @splat(0.0);
    var ur: [4]f64 = @splat(0.0);
    var ui: [4]f64 = @splat(0.0);
    var vr: [4]f64 = @splat(0.0);
    var vi: [4]f64 = @splat(0.0);
    var workr: [2]f64 = @splat(0.0);
    var worki: [2]f64 = @splat(0.0);

    const info = wsvdc(CxMat.initPacked(&xr, &xi, 2, 2), &sr, &si, &er, &ei, CxMat.initPacked(&ur, &ui, 2, 2), CxMat.initPacked(&vr, &vi, 2, 2), &workr, &worki, 11);
    try zigstd.testExpectEqual(@as(usize, 0), info);

    const golden = (1.0 + @sqrt(5.0)) / 2.0;
    try zigstd.testExpectTrue(sr[0] >= sr[1]);
    try zigstd.testExpectApproxEqAbs(golden, sr[0], 1e-9);
    try zigstd.testExpectApproxEqAbs(1.0 / golden, sr[1], 1e-9);
    try zigstd.testExpectApproxEqAbs(0.0, si[0], 1e-9);
    try zigstd.testExpectApproxEqAbs(0.0, si[1], 1e-9);

    for (0..2) |row| {
        for (0..2) |col| {
            var re: f64 = 0.0;
            var im: f64 = 0.0;
            for (0..2) |k| {
                // U(row,k) * S(k) * conj(V(col,k))
                const ur_ = ur[ix(2, row + 1, k + 1)];
                const ui_ = ui[ix(2, row + 1, k + 1)];
                const vr_ = vr[ix(2, col + 1, k + 1)];
                const vi_ = vi[ix(2, col + 1, k + 1)];
                const u_times_s = Cx.init(ur_ * sr[k], ui_ * sr[k]);
                const term = blas1.wmul(u_times_s, Cx.init(vr_, -vi_));
                re += term.re;
                im += term.im;
            }
            try zigstd.testExpectApproxEqAbs(x[ix(2, row + 1, col + 1)], re, 1e-8);
            try zigstd.testExpectApproxEqAbs(0.0, im, 1e-8);
        }
    }
}

test "wsvdc handles a genuinely complex diagonal matrix" {
    // X = diag(i, 2): singular values are |i|=1 and |2|=2.
    var xr = [_]f64{ 0.0, 0.0, 0.0, 2.0 };
    var xi = [_]f64{ 1.0, 0.0, 0.0, 0.0 };
    var sr: [2]f64 = @splat(0.0);
    var si: [2]f64 = @splat(0.0);
    var er: [2]f64 = @splat(0.0);
    var ei: [2]f64 = @splat(0.0);
    var ur: [4]f64 = @splat(0.0);
    var ui: [4]f64 = @splat(0.0);
    var vr: [4]f64 = @splat(0.0);
    var vi: [4]f64 = @splat(0.0);
    var workr: [2]f64 = @splat(0.0);
    var worki: [2]f64 = @splat(0.0);

    const info = wsvdc(CxMat.initPacked(&xr, &xi, 2, 2), &sr, &si, &er, &ei, CxMat.initPacked(&ur, &ui, 2, 2), CxMat.initPacked(&vr, &vi, 2, 2), &workr, &worki, 11);
    try zigstd.testExpectEqual(@as(usize, 0), info);
    try zigstd.testExpectApproxEqAbs(2.0, sr[0], 1e-9);
    try zigstd.testExpectApproxEqAbs(1.0, sr[1], 1e-9);
}
