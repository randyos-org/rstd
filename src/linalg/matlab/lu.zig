//! Complex LU factorization/solve/determinant/inverse, ported from
//! `matlab88.f90`'s WGECO/WGEFA/WGESL/WGEDI (themselves LINPACK's
//! ZGECO/ZGEFA/ZGESL/ZGEDI, split into real/imaginary array pairs). See
//! `blas1.zig` for the shared matrix/vector conventions this port uses,
//! and `CxMat.zig` for the `self`-shaped calling convention matrix-level
//! functions like these now use - `ipvt`/`br`/`bi`/`workr`/`worki` stay
//! loose parameters since they're vector-shaped, not matrix-shaped.

const std = @import("std");
const blas1 = @import("blas1.zig");
const CxMat = @import("CxMat.zig");
const Cx = blas1.Cx;
const ix = blas1.ix;

/// WGEFA: factor a complex N×N matrix by Gaussian elimination with partial
/// pivoting. `self` (column-major, leading dimension `self.lda`) is
/// overwritten with the packed L\U factors; `ipvt` (length N) receives the
/// pivot vector. Returns 0 on success, or the 1-based column index of the
/// first exactly-zero pivot (`wgesl`/`wgedi` will divide by zero if called
/// after that - use `wgeco`'s `rcond` for a reliable singularity check).
pub fn wgefa(self: CxMat, ipvt: []usize) usize {
    const ar = self.r;
    const ai = self.i;
    const lda = self.lda;
    const n = self.n;

    var info: usize = 0;
    if (n >= 2) {
        var k: usize = 1;
        while (k <= n - 1) : (k += 1) {
            const l = blas1.iwamax(n - k + 1, ar[ix(lda, k, k)..], ai[ix(lda, k, k)..], 1) + k - 1;
            ipvt[k - 1] = l;
            if (blas1.cabs1(ar[ix(lda, l, k)], ai[ix(lda, l, k)]) == 0.0) {
                info = k;
            } else {
                if (l != k) {
                    std.mem.swap(f64, &ar[ix(lda, l, k)], &ar[ix(lda, k, k)]);
                    std.mem.swap(f64, &ai[ix(lda, l, k)], &ai[ix(lda, k, k)]);
                }
                const t = blas1.wdiv(Cx.init(-1.0, 0.0), Cx.init(ar[ix(lda, k, k)], ai[ix(lda, k, k)]));
                blas1.wscal(n - k, t, ar[ix(lda, k + 1, k)..], ai[ix(lda, k + 1, k)..], 1);

                var j: usize = k + 1;
                while (j <= n) : (j += 1) {
                    const tr = ar[ix(lda, l, j)];
                    const ti = ai[ix(lda, l, j)];
                    if (l != k) {
                        std.mem.swap(f64, &ar[ix(lda, l, j)], &ar[ix(lda, k, j)]);
                        std.mem.swap(f64, &ai[ix(lda, l, j)], &ai[ix(lda, k, j)]);
                    }
                    blas1.waxpy(
                        n - k,
                        Cx.init(tr, ti),
                        ar[ix(lda, k + 1, k)..],
                        ai[ix(lda, k + 1, k)..],
                        1,
                        ar[ix(lda, k + 1, j)..],
                        ai[ix(lda, k + 1, j)..],
                        1,
                    );
                }
            }
        }
    }
    ipvt[n - 1] = n;
    if (blas1.cabs1(ar[ix(lda, n, n)], ai[ix(lda, n, n)]) == 0.0) info = n;
    return info;
}

/// WGESL: solve `A*x = b` (`job == 0`) or `ctrans(A)*x = b` (`job != 0`)
/// using the factors from `wgefa`/`wgeco` (`self`). `b` is overwritten
/// with `x`.
pub fn wgesl(self: CxMat, ipvt: []const usize, br: []f64, bi: []f64, job: usize) void {
    const ar = self.r;
    const ai = self.i;
    const lda = self.lda;
    const n = self.n;

    if (job == 0) {
        if (n >= 2) {
            var k: usize = 1;
            while (k <= n - 1) : (k += 1) {
                const l = ipvt[k - 1];
                const tr = br[l - 1];
                const ti = bi[l - 1];
                if (l != k) {
                    br[l - 1] = br[k - 1];
                    bi[l - 1] = bi[k - 1];
                    br[k - 1] = tr;
                    bi[k - 1] = ti;
                }
                blas1.waxpy(n - k, Cx.init(tr, ti), ar[ix(lda, k + 1, k)..], ai[ix(lda, k + 1, k)..], 1, br[k..], bi[k..], 1);
            }
        }
        var kb: usize = 1;
        while (kb <= n) : (kb += 1) {
            const k = n + 1 - kb;
            const d = blas1.wdiv(Cx.init(br[k - 1], bi[k - 1]), Cx.init(ar[ix(lda, k, k)], ai[ix(lda, k, k)]));
            br[k - 1] = d.re;
            bi[k - 1] = d.im;
            blas1.waxpy(k - 1, Cx.init(-d.re, -d.im), ar[ix(lda, 1, k)..], ai[ix(lda, 1, k)..], 1, br, bi, 1);
        }
        return;
    }

    var k: usize = 1;
    while (k <= n) : (k += 1) {
        const tr = br[k - 1] - blas1.wdotcr(k - 1, ar[ix(lda, 1, k)..], ai[ix(lda, 1, k)..], 1, br, bi, 1);
        const ti = bi[k - 1] - blas1.wdotci(k - 1, ar[ix(lda, 1, k)..], ai[ix(lda, 1, k)..], 1, br, bi, 1);
        const d = blas1.wdiv(Cx.init(tr, ti), Cx.init(ar[ix(lda, k, k)], -ai[ix(lda, k, k)]));
        br[k - 1] = d.re;
        bi[k - 1] = d.im;
    }
    if (n >= 2) {
        var kb: usize = 1;
        while (kb <= n - 1) : (kb += 1) {
            const k2 = n - kb;
            br[k2 - 1] += blas1.wdotcr(n - k2, ar[ix(lda, k2 + 1, k2)..], ai[ix(lda, k2 + 1, k2)..], 1, br[k2..], bi[k2..], 1);
            bi[k2 - 1] += blas1.wdotci(n - k2, ar[ix(lda, k2 + 1, k2)..], ai[ix(lda, k2 + 1, k2)..], 1, br[k2..], bi[k2..], 1);
            const l = ipvt[k2 - 1];
            if (l != k2) {
                std.mem.swap(f64, &br[l - 1], &br[k2 - 1]);
                std.mem.swap(f64, &bi[l - 1], &bi[k2 - 1]);
            }
        }
    }
}

/// WGECO: factor a complex N×N matrix (`self`, via `wgefa`) and estimate
/// its reciprocal condition number. `zr`/`zi` (length N, scratch) receive
/// an approximate null vector when `A` is nearly singular. If `rcond` is
/// not needed, `wgefa` alone is faster.
pub fn wgeco(self: CxMat, ipvt: []usize, zr: []f64, zi: []f64) f64 {
    const ar = self.r;
    const ai = self.i;
    const lda = self.lda;
    const n = self.n;

    var anorm: f64 = 0.0;
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        anorm = @max(anorm, blas1.wasum(n, ar[ix(lda, 1, j)..], ai[ix(lda, 1, j)..], 1));
    }

    _ = wgefa(self, ipvt);

    // Solve ctrans(U)*w = e, choosing e to grow w as much as possible.
    var ek = Cx.init(1.0, 0.0);
    blas1.wset(n, Cx.init(0.0, 0.0), zr, zi, 1);
    var k: usize = 1;
    while (k <= n) : (k += 1) {
        const zk = Cx.init(zr[k - 1], zi[k - 1]);
        ek = blas1.wsign(ek, Cx.init(-zk.re, -zk.im));
        const akk = Cx.init(ar[ix(lda, k, k)], ai[ix(lda, k, k)]);
        if (blas1.cabs1(ek.re - zk.re, ek.im - zk.im) > blas1.cabs1(akk.re, akk.im)) {
            const s = blas1.cabs1(akk.re, akk.im) / blas1.cabs1(ek.re - zk.re, ek.im - zk.im);
            blas1.wrscal(n, s, zr, zi, 1);
            ek = Cx.init(s * ek.re, s * ek.im);
        }
        var wk = Cx.init(ek.re - zr[k - 1], ek.im - zi[k - 1]);
        var wkm = Cx.init(-ek.re - zr[k - 1], -ek.im - zi[k - 1]);
        var s = blas1.cabs1(wk.re, wk.im);
        var sm = blas1.cabs1(wkm.re, wkm.im);
        if (blas1.cabs1(akk.re, akk.im) == 0.0) {
            wk = Cx.init(1.0, 0.0);
            wkm = Cx.init(1.0, 0.0);
        } else {
            wk = blas1.wdiv(wk, Cx.init(akk.re, -akk.im));
            wkm = blas1.wdiv(wkm, Cx.init(akk.re, -akk.im));
        }
        if (k + 1 <= n) {
            var jc: usize = k + 1;
            while (jc <= n) : (jc += 1) {
                const akj_conj = Cx.init(ar[ix(lda, k, jc)], -ai[ix(lda, k, jc)]);
                const t = blas1.wmul(wkm, akj_conj);
                sm += blas1.cabs1(zr[jc - 1] + t.re, zi[jc - 1] + t.im);
                const u = blas1.wmul(wk, akj_conj);
                zr[jc - 1] += u.re;
                zi[jc - 1] += u.im;
                s += blas1.cabs1(zr[jc - 1], zi[jc - 1]);
            }
            if (s < sm) {
                const d = Cx.init(wkm.re - wk.re, wkm.im - wk.im);
                wk = wkm;
                var jc2: usize = k + 1;
                while (jc2 <= n) : (jc2 += 1) {
                    const akj_conj = Cx.init(ar[ix(lda, k, jc2)], -ai[ix(lda, k, jc2)]);
                    const u = blas1.wmul(d, akj_conj);
                    zr[jc2 - 1] += u.re;
                    zi[jc2 - 1] += u.im;
                }
            }
        }
        zr[k - 1] = wk.re;
        zi[k - 1] = wk.im;
    }
    var s_norm = 1.0 / blas1.wasum(n, zr, zi, 1);
    blas1.wrscal(n, s_norm, zr, zi, 1);

    // Solve ctrans(L)*y = w.
    var kb: usize = 1;
    while (kb <= n) : (kb += 1) {
        const k2 = n + 1 - kb;
        if (k2 < n) {
            zr[k2 - 1] += blas1.wdotcr(n - k2, ar[ix(lda, k2 + 1, k2)..], ai[ix(lda, k2 + 1, k2)..], 1, zr[k2..], zi[k2..], 1);
            zi[k2 - 1] += blas1.wdotci(n - k2, ar[ix(lda, k2 + 1, k2)..], ai[ix(lda, k2 + 1, k2)..], 1, zr[k2..], zi[k2..], 1);
        }
        if (blas1.cabs1(zr[k2 - 1], zi[k2 - 1]) > 1.0) {
            const s = 1.0 / blas1.cabs1(zr[k2 - 1], zi[k2 - 1]);
            blas1.wrscal(n, s, zr, zi, 1);
        }
        const l2 = ipvt[k2 - 1];
        std.mem.swap(f64, &zr[l2 - 1], &zr[k2 - 1]);
        std.mem.swap(f64, &zi[l2 - 1], &zi[k2 - 1]);
    }
    s_norm = 1.0 / blas1.wasum(n, zr, zi, 1);
    blas1.wrscal(n, s_norm, zr, zi, 1);

    // Solve L*v = y.
    var ynorm: f64 = 1.0;
    var k3: usize = 1;
    while (k3 <= n) : (k3 += 1) {
        const l3 = ipvt[k3 - 1];
        const t = Cx.init(zr[l3 - 1], zi[l3 - 1]);
        std.mem.swap(f64, &zr[l3 - 1], &zr[k3 - 1]);
        std.mem.swap(f64, &zi[l3 - 1], &zi[k3 - 1]);
        if (k3 < n) {
            blas1.waxpy(n - k3, t, ar[ix(lda, k3 + 1, k3)..], ai[ix(lda, k3 + 1, k3)..], 1, zr[k3..], zi[k3..], 1);
        }
        if (blas1.cabs1(zr[k3 - 1], zi[k3 - 1]) > 1.0) {
            const s = 1.0 / blas1.cabs1(zr[k3 - 1], zi[k3 - 1]);
            blas1.wrscal(n, s, zr, zi, 1);
            ynorm *= s;
        }
    }
    s_norm = 1.0 / blas1.wasum(n, zr, zi, 1);
    blas1.wrscal(n, s_norm, zr, zi, 1);
    ynorm *= s_norm;

    // Solve U*z = v.
    var kb2: usize = 1;
    while (kb2 <= n) : (kb2 += 1) {
        const k4 = n + 1 - kb2;
        const akk = Cx.init(ar[ix(lda, k4, k4)], ai[ix(lda, k4, k4)]);
        if (blas1.cabs1(zr[k4 - 1], zi[k4 - 1]) > blas1.cabs1(akk.re, akk.im)) {
            const s = blas1.cabs1(akk.re, akk.im) / blas1.cabs1(zr[k4 - 1], zi[k4 - 1]);
            blas1.wrscal(n, s, zr, zi, 1);
            ynorm *= s;
        }
        if (blas1.cabs1(akk.re, akk.im) != 0.0) {
            const d = blas1.wdiv(Cx.init(zr[k4 - 1], zi[k4 - 1]), akk);
            zr[k4 - 1] = d.re;
            zi[k4 - 1] = d.im;
        } else {
            zr[k4 - 1] = 1.0;
            zi[k4 - 1] = 0.0;
        }
        const t = Cx.init(-zr[k4 - 1], -zi[k4 - 1]);
        blas1.waxpy(k4 - 1, t, ar[ix(lda, 1, k4)..], ai[ix(lda, 1, k4)..], 1, zr, zi, 1);
    }
    s_norm = 1.0 / blas1.wasum(n, zr, zi, 1);
    blas1.wrscal(n, s_norm, zr, zi, 1);
    ynorm *= s_norm;

    return if (anorm != 0.0) ynorm / anorm else 0.0;
}

/// Determinant as `mantissa * 10^exponent10`, with `1 <= |mantissa| < 10`
/// (or `mantissa == 0`) - kept scaled like the original rather than
/// multiplied out, since `exponent10` can be far larger than `f64` can hold
/// in the mantissa alone for large ill-scaled matrices.
pub const Determinant = struct { mantissa: Cx, exponent10: f64 };

/// WGEDI: compute the determinant and/or inverse of a matrix (`self`) from
/// the factors produced by `wgefa`/`wgeco`. `job`: 11 = both, 01 = inverse
/// only, 10 = determinant only (`det_out` is unused and may be `null` when
/// the determinant isn't requested). `workr`/`worki` are length-N scratch.
/// Dividing by zero here means the input had an exactly-zero pivot; check
/// `wgefa`'s `info` or `wgeco`'s `rcond` first if that's a possibility.
pub fn wgedi(
    self: CxMat,
    ipvt: []const usize,
    workr: []f64,
    worki: []f64,
    job: usize,
    det_out: ?*Determinant,
) void {
    const ar = self.r;
    const ai = self.i;
    const lda = self.lda;
    const n = self.n;

    if (job / 10 != 0) {
        var det = Cx.init(1.0, 0.0);
        var exp10: f64 = 0.0;
        var i: usize = 1;
        while (i <= n) : (i += 1) {
            if (ipvt[i - 1] != i) det = Cx.init(-det.re, -det.im);
            det = blas1.wmul(Cx.init(ar[ix(lda, i, i)], ai[ix(lda, i, i)]), det);
            if (blas1.cabs1(det.re, det.im) == 0.0) break;
            while (blas1.cabs1(det.re, det.im) < 1.0) {
                det = Cx.init(det.re * 10.0, det.im * 10.0);
                exp10 -= 1.0;
            }
            while (blas1.cabs1(det.re, det.im) >= 10.0) {
                det = Cx.init(det.re / 10.0, det.im / 10.0);
                exp10 += 1.0;
            }
        }
        if (det_out) |out| out.* = .{ .mantissa = det, .exponent10 = exp10 };
    }

    if (job % 10 == 0) return;

    // Invert U in place.
    var k: usize = 1;
    while (k <= n) : (k += 1) {
        const akk = Cx.init(ar[ix(lda, k, k)], ai[ix(lda, k, k)]);
        const inv = blas1.wdiv(Cx.init(1.0, 0.0), akk);
        ar[ix(lda, k, k)] = inv.re;
        ai[ix(lda, k, k)] = inv.im;
        blas1.wscal(k - 1, Cx.init(-inv.re, -inv.im), ar[ix(lda, 1, k)..], ai[ix(lda, 1, k)..], 1);
        if (k + 1 <= n) {
            var j: usize = k + 1;
            while (j <= n) : (j += 1) {
                const t = Cx.init(ar[ix(lda, k, j)], ai[ix(lda, k, j)]);
                ar[ix(lda, k, j)] = 0.0;
                ai[ix(lda, k, j)] = 0.0;
                blas1.waxpy(k, t, ar[ix(lda, 1, k)..], ai[ix(lda, 1, k)..], 1, ar[ix(lda, 1, j)..], ai[ix(lda, 1, j)..], 1);
            }
        }
    }

    // Form inverse(U) * inverse(L).
    if (n < 2) return;
    var kb: usize = 1;
    while (kb <= n - 1) : (kb += 1) {
        const k2 = n - kb;
        const kp1 = k2 + 1;
        var row: usize = kp1;
        while (row <= n) : (row += 1) {
            workr[row - 1] = ar[ix(lda, row, k2)];
            worki[row - 1] = ai[ix(lda, row, k2)];
            ar[ix(lda, row, k2)] = 0.0;
            ai[ix(lda, row, k2)] = 0.0;
        }
        var j2: usize = kp1;
        while (j2 <= n) : (j2 += 1) {
            const t = Cx.init(workr[j2 - 1], worki[j2 - 1]);
            blas1.waxpy(n, t, ar[ix(lda, 1, j2)..], ai[ix(lda, 1, j2)..], 1, ar[ix(lda, 1, k2)..], ai[ix(lda, 1, k2)..], 1);
        }
        const l = ipvt[k2 - 1];
        if (l != k2) {
            blas1.wswap(n, ar[ix(lda, 1, k2)..], ai[ix(lda, 1, k2)..], 1, ar[ix(lda, 1, l)..], ai[ix(lda, 1, l)..], 1);
        }
    }
}

const zigstd = @import("../../zigstd.zig");

test "wgesl solves a known complex 2x2 system" {
    // A = [[1+i, 2], [3, 4-i]], b chosen so x = [1, 1+i].
    var ar = [_]f64{ 1.0, 3.0, 2.0, 4.0 }; // column-major
    var ai = [_]f64{ 1.0, 0.0, 0.0, -1.0 };
    var ipvt = [_]usize{ 0, 0 };
    var zr = [_]f64{ 0.0, 0.0 };
    var zi = [_]f64{ 0.0, 0.0 };
    const mat = CxMat.initPacked(&ar, &ai, 2, 2);
    const rcond = wgeco(mat, &ipvt, &zr, &zi);
    try zigstd.testExpectTrue(rcond > 0.0);

    // b = A * [1, 1+i]
    var br = [_]f64{ 1.0 * 1.0 - 1.0 * 0.0 + 2.0 * 1.0, 3.0 * 1.0 + 4.0 * 1.0 - (-1.0) * 1.0 };
    var bi = [_]f64{ 1.0 * 0.0 + 1.0 * 1.0 + 2.0 * 1.0, 3.0 * 0.0 + 4.0 * 1.0 + (-1.0) * 1.0 };
    wgesl(mat, &ipvt, &br, &bi, 0);
    try zigstd.testExpectApproxEqAbs(1.0, br[0], 1e-10);
    try zigstd.testExpectApproxEqAbs(0.0, bi[0], 1e-10);
    try zigstd.testExpectApproxEqAbs(1.0, br[1], 1e-10);
    try zigstd.testExpectApproxEqAbs(1.0, bi[1], 1e-10);
}

test "wgedi inverts and finds the determinant of a diagonal matrix" {
    // A = diag(2, 3, -4) (real).
    var ar = [_]f64{ 2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, -4.0 };
    var ai: [9]f64 = @splat(0.0);
    var ipvt = [_]usize{ 0, 0, 0 };
    const mat = CxMat.initPacked(&ar, &ai, 3, 3);
    const info = wgefa(mat, &ipvt);
    try zigstd.testExpectEqual(@as(usize, 0), info);

    var det: Determinant = undefined;
    var workr: [3]f64 = @splat(0.0);
    var worki: [3]f64 = @splat(0.0);
    wgedi(mat, &ipvt, &workr, &worki, 11, &det);

    const value = det.mantissa.re * std.math.pow(f64, 10.0, det.exponent10);
    try zigstd.testExpectApproxEqAbs(-24.0, value, 1e-9);
    try zigstd.testExpectApproxEqAbs(0.0, det.mantissa.im, 1e-9);

    // Inverse of diag(2,3,-4) is diag(0.5, 1/3, -0.25).
    try zigstd.testExpectApproxEqAbs(0.5, ar[ix(3, 1, 1)], 1e-10);
    try zigstd.testExpectApproxEqAbs(1.0 / 3.0, ar[ix(3, 2, 2)], 1e-10);
    try zigstd.testExpectApproxEqAbs(-0.25, ar[ix(3, 3, 3)], 1e-10);
}
