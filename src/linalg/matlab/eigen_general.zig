//! Eigenvalues/eigenvectors of a general complex matrix by the (shifted)
//! QR method, ported from `matlab88.f90`'s COMQR3 - a modification of
//! EISPACK's COMQR2 (itself the unitary analogue of COMLR2/Francis's QR
//! algorithm) that adds a `job` parameter to select how much work to do.
//! See `blas1.zig` for the shared matrix/vector conventions this port
//! uses, and `hessenberg.zig` (`corth`) for the routine that produces the
//! Hessenberg form + transformation info this one consumes.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1.zig");
const CxMat = @import("CxMat.zig");
const Cx = blas1.Cx;
const ix = blas1.ix;

/// COMQR3: find the eigenvalues (`wr`/`wi`) of the complex upper
/// Hessenberg matrix `hr`/`hi` (column-major, leading dimension `nm`,
/// order `n`, rows/columns `low..igh` significant - pass `low=1, igh=n` if
/// `corth` was not preceded by a balancing step) and, depending on `job`,
/// its Schur form and/or eigenvectors:
///
/// - `job == 0`: `hr`/`hi` becomes the Schur (triangular) form; `zr`/`zi`
///   unused.
/// - `job == 1`: Schur form in `hr`/`hi`, plus the unitary similarity in
///   `zr`/`zi`.
/// - `job == 2`: same as EISPACK's COMQR2 - eigenvalues in `wr`/`wi` and
///   (unnormalized) eigenvectors of the *original* general matrix in
///   `zr`/`zi`, recovered via `ortr`/`orti` (`corth`'s output).
/// - `job == 3`: only accumulate the Hessenberg-reduction similarity into
///   `zr`/`zi` (from `ortr`/`orti`) and return - no QR iteration at all.
///
/// `ortr`/`orti` (indexed `low..igh`) must be `corth`'s output; pass zeros
/// there if `hr`/`hi` is a Hessenberg matrix in its own right (no prior
/// reduction to back out).
///
/// Returns 0 on success, or the 1-based index of an eigenvalue that failed
/// to converge after `30*n` iterations total (eigenvalues past that index
/// are still correct; none of the eigenvectors are valid in that case).
pub fn comqr3(
    self: CxMat,
    low: usize,
    igh: usize,
    ortr: []f64,
    orti: []f64,
    wr: []f64,
    wi: []f64,
    z: CxMat,
    job: usize,
) usize {
    const nm = self.lda;
    const n = self.n;
    const hr = self.r;
    const hi = self.i;
    const zr = z.r;
    const zi = z.i;

    var skip_realification = false;

    if (job != 0) {
        var row: usize = 1;
        while (row <= n) : (row += 1) {
            var col: usize = 1;
            while (col <= n) : (col += 1) {
                zr[ix(nm, row, col)] = if (row == col) 1.0 else 0.0;
                zi[ix(nm, row, col)] = 0.0;
            }
        }

        if (igh > low + 1) {
            // Form the accumulated transformation from corth's output.
            var ii: usize = 1;
            while (ii <= igh - low - 1) : (ii += 1) {
                const i = igh - ii;
                if ((ortr[i - 1] != 0.0 or orti[i - 1] != 0.0) and
                    (hr[ix(nm, i, i - 1)] != 0.0 or hi[ix(nm, i, i - 1)] != 0.0))
                {
                    // norm below is negative of h formed in corth.
                    const norm = hr[ix(nm, i, i - 1)] * ortr[i - 1] + hi[ix(nm, i, i - 1)] * orti[i - 1];
                    var k: usize = i + 1;
                    while (k <= igh) : (k += 1) {
                        ortr[k - 1] = hr[ix(nm, k, i - 1)];
                        orti[k - 1] = hi[ix(nm, k, i - 1)];
                    }
                    var j: usize = i;
                    while (j <= igh) : (j += 1) {
                        var sr: f64 = 0.0;
                        var si: f64 = 0.0;
                        var k2: usize = i;
                        while (k2 <= igh) : (k2 += 1) {
                            sr += ortr[k2 - 1] * zr[ix(nm, k2, j)] + orti[k2 - 1] * zi[ix(nm, k2, j)];
                            si += ortr[k2 - 1] * zi[ix(nm, k2, j)] - orti[k2 - 1] * zr[ix(nm, k2, j)];
                        }
                        sr /= norm;
                        si /= norm;
                        var k3: usize = i;
                        while (k3 <= igh) : (k3 += 1) {
                            zr[ix(nm, k3, j)] += sr * ortr[k3 - 1] - si * orti[k3 - 1];
                            zi[ix(nm, k3, j)] += sr * orti[k3 - 1] + si * ortr[k3 - 1];
                        }
                    }
                }
            }
            if (job == 3) return 0;
        } else if (igh < low + 1) {
            skip_realification = true;
        }
    }

    if (!skip_realification) {
        // Create real subdiagonal elements.
        var i: usize = low + 1;
        while (i <= igh) : (i += 1) {
            const ll = @min(i + 1, igh);
            if (hi[ix(nm, i, i - 1)] == 0.0) continue;
            const norm = blas1.pythag(hr[ix(nm, i, i - 1)], hi[ix(nm, i, i - 1)]);
            const yr = hr[ix(nm, i, i - 1)] / norm;
            const yi = hi[ix(nm, i, i - 1)] / norm;
            hr[ix(nm, i, i - 1)] = norm;
            hi[ix(nm, i, i - 1)] = 0.0;

            var j: usize = i;
            while (j <= n) : (j += 1) {
                const si = yr * hi[ix(nm, i, j)] - yi * hr[ix(nm, i, j)];
                hr[ix(nm, i, j)] = yr * hr[ix(nm, i, j)] + yi * hi[ix(nm, i, j)];
                hi[ix(nm, i, j)] = si;
            }
            var j2: usize = 1;
            while (j2 <= ll) : (j2 += 1) {
                const si = yr * hi[ix(nm, j2, i)] + yi * hr[ix(nm, j2, i)];
                hr[ix(nm, j2, i)] = yr * hr[ix(nm, j2, i)] - yi * hi[ix(nm, j2, i)];
                hi[ix(nm, j2, i)] = si;
            }
            if (job != 0) {
                var j3: usize = low;
                while (j3 <= igh) : (j3 += 1) {
                    const si = yr * zi[ix(nm, j3, i)] + yi * zr[ix(nm, j3, i)];
                    zr[ix(nm, j3, i)] = yr * zr[ix(nm, j3, i)] - yi * zi[ix(nm, j3, i)];
                    zi[ix(nm, j3, i)] = si;
                }
            }
        }
    }

    // Store roots isolated outside [low, igh].
    var idx: usize = 1;
    while (idx <= n) : (idx += 1) {
        if (idx >= low and idx <= igh) continue;
        wr[idx - 1] = hr[ix(nm, idx, idx)];
        wi[idx - 1] = hi[ix(nm, idx, idx)];
    }

    var en = igh;
    var tr: f64 = 0.0;
    var ti: f64 = 0.0;
    var itn: usize = 30 * n;

    outer_en: while (en >= low) {
        var its: usize = 0;
        const enm1 = en - 1;

        qr_step: while (true) {
            // Look for a single small subdiagonal element.
            var l: usize = low;
            search: {
                var ll: usize = low;
                while (ll <= en) : (ll += 1) {
                    l = en + low - ll;
                    if (l == low) break :search;
                    const xr = @abs(hr[ix(nm, l - 1, l - 1)]) + @abs(hi[ix(nm, l - 1, l - 1)]) +
                        @abs(hr[ix(nm, l, l)]) + @abs(hi[ix(nm, l, l)]);
                    const yr = xr + @abs(hr[ix(nm, l, l - 1)]);
                    if (xr == yr) break :search;
                }
            }

            if (l == en) {
                // A root is found.
                hr[ix(nm, en, en)] += tr;
                wr[en - 1] = hr[ix(nm, en, en)];
                hi[ix(nm, en, en)] += ti;
                wi[en - 1] = hi[ix(nm, en, en)];
                if (en == low) break :outer_en;
                en -= 1;
                continue :outer_en;
            }
            if (itn == 0) return en;

            var sr: f64 = undefined;
            var si: f64 = undefined;
            if (its == 10 or its == 20) {
                // Exceptional shift.
                sr = @abs(hr[ix(nm, en, enm1)]) + @abs(hr[ix(nm, enm1, en - 1)]);
                si = 0.0;
            } else {
                sr = hr[ix(nm, en, en)];
                si = hi[ix(nm, en, en)];
                const xr = hr[ix(nm, enm1, en)] * hr[ix(nm, en, enm1)];
                const xi = hi[ix(nm, enm1, en)] * hr[ix(nm, en, enm1)];
                if (xr != 0.0 or xi != 0.0) {
                    const yr = (hr[ix(nm, enm1, enm1)] - sr) / 2.0;
                    const yi = (hi[ix(nm, enm1, enm1)] - si) / 2.0;
                    var zz = std.math.complex.sqrt(Cx.init(yr * yr - yi * yi + xr, 2.0 * yr * yi + xi));
                    if (yr * zz.re + yi * zz.im < 0.0) zz = Cx.init(-zz.re, -zz.im);
                    const d = blas1.wdiv(Cx.init(xr, xi), Cx.init(yr + zz.re, yi + zz.im));
                    sr -= d.re;
                    si -= d.im;
                }
            }

            var i: usize = low;
            while (i <= en) : (i += 1) {
                hr[ix(nm, i, i)] -= sr;
                hi[ix(nm, i, i)] -= si;
            }
            tr += sr;
            ti += si;
            its += 1;
            itn -= 1;

            // Reduce to triangle (rows).
            const lp1 = l + 1;
            var row: usize = lp1;
            while (row <= en) : (row += 1) {
                var s = hr[ix(nm, row, row - 1)];
                hr[ix(nm, row, row - 1)] = 0.0;
                var norm = @abs(hr[ix(nm, row - 1, row - 1)]) + @abs(hi[ix(nm, row - 1, row - 1)]) + @abs(s);
                norm *= @sqrt((hr[ix(nm, row - 1, row - 1)] / norm) * (hr[ix(nm, row - 1, row - 1)] / norm) +
                    (hi[ix(nm, row - 1, row - 1)] / norm) * (hi[ix(nm, row - 1, row - 1)] / norm) +
                    (s / norm) * (s / norm));
                const xr = hr[ix(nm, row - 1, row - 1)] / norm;
                wr[row - 2] = xr;
                const xi = hi[ix(nm, row - 1, row - 1)] / norm;
                wi[row - 2] = xi;
                hr[ix(nm, row - 1, row - 1)] = norm;
                hi[ix(nm, row - 1, row - 1)] = 0.0;
                hi[ix(nm, row, row - 1)] = s / norm;
                s = hi[ix(nm, row, row - 1)];

                var j: usize = row;
                while (j <= n) : (j += 1) {
                    const yr = hr[ix(nm, row - 1, j)];
                    const yi = hi[ix(nm, row - 1, j)];
                    const zzr = hr[ix(nm, row, j)];
                    const zzi = hi[ix(nm, row, j)];
                    hr[ix(nm, row - 1, j)] = xr * yr + xi * yi + s * zzr;
                    hi[ix(nm, row - 1, j)] = xr * yi - xi * yr + s * zzi;
                    hr[ix(nm, row, j)] = xr * zzr - xi * zzi - s * yr;
                    hi[ix(nm, row, j)] = xr * zzi + xi * zzr - s * yi;
                }
            }

            var diag_si = hi[ix(nm, en, en)];
            var diag_sr: f64 = undefined;
            if (diag_si != 0.0) {
                const norm = blas1.pythag(hr[ix(nm, en, en)], diag_si);
                diag_sr = hr[ix(nm, en, en)] / norm;
                diag_si /= norm;
                hr[ix(nm, en, en)] = norm;
                hi[ix(nm, en, en)] = 0.0;
                if (en != n) {
                    var j: usize = en + 1;
                    while (j <= n) : (j += 1) {
                        const yr = hr[ix(nm, en, j)];
                        const yi = hi[ix(nm, en, j)];
                        hr[ix(nm, en, j)] = diag_sr * yr + diag_si * yi;
                        hi[ix(nm, en, j)] = diag_sr * yi - diag_si * yr;
                    }
                }
            } else {
                diag_sr = undefined;
            }

            // Inverse operation (columns).
            var col: usize = lp1;
            while (col <= en) : (col += 1) {
                const xr = wr[col - 2];
                const xi = wi[col - 2];
                var rh: usize = 1;
                while (rh <= col) : (rh += 1) {
                    const yr = hr[ix(nm, rh, col - 1)];
                    var yi: f64 = 0.0;
                    const zzr = hr[ix(nm, rh, col)];
                    const zzi = hi[ix(nm, rh, col)];
                    if (rh != col) {
                        yi = hi[ix(nm, rh, col - 1)];
                        hi[ix(nm, rh, col - 1)] = xr * yi + xi * yr + hi[ix(nm, col, col - 1)] * zzi;
                    }
                    hr[ix(nm, rh, col - 1)] = xr * yr - xi * yi + hi[ix(nm, col, col - 1)] * zzr;
                    hr[ix(nm, rh, col)] = xr * zzr + xi * zzi - hi[ix(nm, col, col - 1)] * yr;
                    hi[ix(nm, rh, col)] = xr * zzi - xi * zzr - hi[ix(nm, col, col - 1)] * yi;
                }
                if (job != 0) {
                    var rz: usize = low;
                    while (rz <= igh) : (rz += 1) {
                        const yr = zr[ix(nm, rz, col - 1)];
                        const yi = zi[ix(nm, rz, col - 1)];
                        const zzr = zr[ix(nm, rz, col)];
                        const zzi = zi[ix(nm, rz, col)];
                        zr[ix(nm, rz, col - 1)] = xr * yr - xi * yi + hi[ix(nm, col, col - 1)] * zzr;
                        zi[ix(nm, rz, col - 1)] = xr * yi + xi * yr + hi[ix(nm, col, col - 1)] * zzi;
                        zr[ix(nm, rz, col)] = xr * zzr + xi * zzi - hi[ix(nm, col, col - 1)] * yr;
                        zi[ix(nm, rz, col)] = xr * zzi - xi * zzr - hi[ix(nm, col, col - 1)] * yi;
                    }
                }
            }

            if (diag_si == 0.0) continue :qr_step;

            var rh2: usize = 1;
            while (rh2 <= en) : (rh2 += 1) {
                const yr = hr[ix(nm, rh2, en)];
                const yi = hi[ix(nm, rh2, en)];
                hr[ix(nm, rh2, en)] = diag_sr * yr - diag_si * yi;
                hi[ix(nm, rh2, en)] = diag_sr * yi + diag_si * yr;
            }
            if (job != 0) {
                var rz2: usize = low;
                while (rz2 <= igh) : (rz2 += 1) {
                    const yr = zr[ix(nm, rz2, en)];
                    const yi = zi[ix(nm, rz2, en)];
                    zr[ix(nm, rz2, en)] = diag_sr * yr - diag_si * yi;
                    zi[ix(nm, rz2, en)] = diag_sr * yi + diag_si * yr;
                }
            }
            continue :qr_step;
        }
    }

    // All roots found: backsubstitute for eigenvectors of the original
    // matrix (job == 2 only - EISPACK's COMQR2 behavior).
    if (job != 2) return 0;

    var norm: f64 = 0.0;
    {
        var i: usize = 1;
        while (i <= n) : (i += 1) {
            var j: usize = i;
            while (j <= n) : (j += 1) {
                const t = @abs(hr[ix(nm, i, j)]) + @abs(hi[ix(nm, i, j)]);
                if (t > norm) norm = t;
            }
        }
    }
    if (n == 1 or norm == 0.0) return 0;

    var nn: usize = 2;
    while (nn <= n) : (nn += 1) {
        const en2 = n + 2 - nn;
        const xr = wr[en2 - 1];
        const xi = wi[en2 - 1];
        hr[ix(nm, en2, en2)] = 1.0;
        hi[ix(nm, en2, en2)] = 0.0;
        const enm1_2 = en2 - 1;

        var ii: usize = 1;
        while (ii <= enm1_2) : (ii += 1) {
            const i = en2 - ii;
            var zzr: f64 = 0.0;
            var zzi: f64 = 0.0;
            var j: usize = i + 1;
            while (j <= en2) : (j += 1) {
                zzr += hr[ix(nm, i, j)] * hr[ix(nm, j, en2)] - hi[ix(nm, i, j)] * hi[ix(nm, j, en2)];
                zzi += hr[ix(nm, i, j)] * hi[ix(nm, j, en2)] + hi[ix(nm, i, j)] * hr[ix(nm, j, en2)];
            }
            var yr = xr - wr[i - 1];
            var yi = xi - wi[i - 1];
            if (yr == 0.0 and yi == 0.0) {
                yr = norm;
                while (true) {
                    yr /= 100.0;
                    const shifted = norm + yr;
                    if (shifted != norm) continue;
                    break;
                }
                yi = 0.0;
            }
            const d = blas1.wdiv(Cx.init(zzr, zzi), Cx.init(yr, yi));
            hr[ix(nm, i, en2)] = d.re;
            hi[ix(nm, i, en2)] = d.im;
            const t = @abs(d.re) + @abs(d.im);
            if (t == 0.0) continue;
            if (t + 1.0 / t > t) continue;
            var j2: usize = i;
            while (j2 <= en2) : (j2 += 1) {
                hr[ix(nm, j2, en2)] /= t;
                hi[ix(nm, j2, en2)] /= t;
            }
        }
    }

    // Vectors of isolated roots.
    const nm1 = n - 1;
    {
        var i: usize = 1;
        while (i <= nm1) : (i += 1) {
            if (i >= low and i <= igh) continue;
            var j: usize = i + 1;
            while (j <= n) : (j += 1) {
                zr[ix(nm, i, j)] = hr[ix(nm, i, j)];
                zi[ix(nm, i, j)] = hi[ix(nm, i, j)];
            }
        }
    }

    // Multiply by the transformation matrix to give eigenvectors of the
    // original full matrix.
    var jj: usize = low;
    while (jj <= nm1) : (jj += 1) {
        const j = n + low - jj;
        const m = @min(j, igh);
        var i: usize = low;
        while (i <= igh) : (i += 1) {
            var zzr: f64 = 0.0;
            var zzi: f64 = 0.0;
            var k: usize = low;
            while (k <= m) : (k += 1) {
                zzr += zr[ix(nm, i, k)] * hr[ix(nm, k, j)] - zi[ix(nm, i, k)] * hi[ix(nm, k, j)];
                zzi += zr[ix(nm, i, k)] * hi[ix(nm, k, j)] + zi[ix(nm, i, k)] * hr[ix(nm, k, j)];
            }
            zr[ix(nm, i, j)] = zzr;
            zi[ix(nm, i, j)] = zzi;
        }
    }

    return 0;
}

const hessenberg = @import("hessenberg.zig");

test "comqr3 finds eigenvalues of a trivial 2x2 upper-triangular Hessenberg matrix" {
    // Already Hessenberg (in fact triangular): eigenvalues are the
    // diagonal, 2 and 3.
    var hr = [_]f64{ 2.0, 0.0, 1.0, 3.0 }; // column-major
    var hi: [4]f64 = @splat(0.0);
    var ortr: [2]f64 = @splat(0.0);
    var orti: [2]f64 = @splat(0.0);
    var wr: [2]f64 = @splat(0.0);
    var wi: [2]f64 = @splat(0.0);
    var zr: [4]f64 = @splat(0.0);
    var zi: [4]f64 = @splat(0.0);

    const ierr = comqr3(CxMat.initPacked(&hr, &hi, 2, 2), 1, 2, &ortr, &orti, &wr, &wi, CxMat.initPacked(&zr, &zi, 2, 2), 0);
    try zigstd.testExpectEqual(@as(usize, 0), ierr);

    const lo = @min(wr[0], wr[1]);
    const hi_val = @max(wr[0], wr[1]);
    try zigstd.testExpectApproxEqAbs(2.0, lo, 1e-10);
    try zigstd.testExpectApproxEqAbs(3.0, hi_val, 1e-10);
    try zigstd.testExpectApproxEqAbs(0.0, wi[0], 1e-10);
    try zigstd.testExpectApproxEqAbs(0.0, wi[1], 1e-10);
}

test "comqr3 recovers eigenvalues and eigenvectors of a symmetric 2x2 matrix" {
    // A = [[2,1],[1,2]]: eigenvalues 1 and 3.
    var hr = [_]f64{ 2.0, 1.0, 1.0, 2.0 }; // column-major
    var hi: [4]f64 = @splat(0.0);
    var ortr: [2]f64 = @splat(0.0);
    var orti: [2]f64 = @splat(0.0);
    var wr: [2]f64 = @splat(0.0);
    var wi: [2]f64 = @splat(0.0);
    var zr: [4]f64 = @splat(0.0);
    var zi: [4]f64 = @splat(0.0);

    const a = [_]f64{ 2.0, 1.0, 1.0, 2.0 }; // original A, untouched, for verification
    const ierr = comqr3(CxMat.initPacked(&hr, &hi, 2, 2), 1, 2, &ortr, &orti, &wr, &wi, CxMat.initPacked(&zr, &zi, 2, 2), 2);
    try zigstd.testExpectEqual(@as(usize, 0), ierr);

    for (0..2) |col| {
        const lambda = wr[col];
        try zigstd.testExpectApproxEqAbs(0.0, wi[col], 1e-9);
        // A*v - lambda*v should vanish (v = column `col` of z, real since A and
        // lambda are real here).
        for (0..2) |row| {
            var av: f64 = 0.0;
            for (0..2) |k| av += a[ix(2, row + 1, k + 1)] * zr[ix(2, k + 1, col + 1)];
            const residual = av - lambda * zr[ix(2, row + 1, col + 1)];
            try zigstd.testExpectApproxEqAbs(0.0, residual, 1e-8);
        }
    }

    const lo = @min(wr[0], wr[1]);
    const hi_val = @max(wr[0], wr[1]);
    try zigstd.testExpectApproxEqAbs(1.0, lo, 1e-9);
    try zigstd.testExpectApproxEqAbs(3.0, hi_val, 1e-9);
}

test "corth + comqr3 recover eigenpairs of a genuinely non-symmetric 3x3 matrix" {
    const a = [_]f64{ 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 }; // column-major

    var hr = a;
    var hi: [9]f64 = @splat(0.0);
    var ortr: [3]f64 = @splat(0.0);
    var orti: [3]f64 = @splat(0.0);
    hessenberg.corth(CxMat.initPacked(&hr, &hi, 3, 3), 1, 3, &ortr, &orti);

    var wr: [3]f64 = @splat(0.0);
    var wi: [3]f64 = @splat(0.0);
    var zr: [9]f64 = @splat(0.0);
    var zi: [9]f64 = @splat(0.0);
    const ierr = comqr3(CxMat.initPacked(&hr, &hi, 3, 3), 1, 3, &ortr, &orti, &wr, &wi, CxMat.initPacked(&zr, &zi, 3, 3), 2);
    try zigstd.testExpectEqual(@as(usize, 0), ierr);

    // Trace is preserved by similarity.
    try zigstd.testExpectApproxEqAbs(16.0, wr[0] + wr[1] + wr[2], 1e-8);

    for (0..3) |col| {
        // A*v - lambda*v should vanish; lambda/v may be complex in general,
        // so carry both real and imaginary parts through.
        const lr = wr[col];
        const li = wi[col];
        for (0..3) |row| {
            var av_r: f64 = 0.0;
            var av_i: f64 = 0.0;
            for (0..3) |k| {
                av_r += a[ix(3, row + 1, k + 1)] * zr[ix(3, k + 1, col + 1)];
                av_i += a[ix(3, row + 1, k + 1)] * zi[ix(3, k + 1, col + 1)];
            }
            const vr = zr[ix(3, row + 1, col + 1)];
            const vi = zi[ix(3, row + 1, col + 1)];
            const lv_r = lr * vr - li * vi;
            const lv_i = lr * vi + li * vr;
            try zigstd.testExpectApproxEqAbs(0.0, av_r - lv_r, 1e-7);
            try zigstd.testExpectApproxEqAbs(0.0, av_i - lv_i, 1e-7);
        }
    }
}
