//! The real QZ's 2 x 2 kernels, from faer v0.24.4 (`faer/src/linalg/gevd/qz_real/mod.rs`,
//! MIT): `svd_triu_2x2` (Reference-LAPACK's `xLASV2`) and `generalized_eigval_2x2` (`xLAG2`).

const std = @import("std");

/// faer's `copy_sign`: `a` or `-a` by `b`'s sign; every caller's `a` is
/// nonnegative, so it is Fortran's `SIGN(a, b)`.
pub fn copySign(comptime T: type, a: T, b: T) T {
    return if (b >= 0) a else -a;
}

/// `xLASV2` (faer's `svd_triu_2x2`): the SVD of `[f g; 0 h]`.
pub fn Svd2(comptime T: type) type {
    return struct { ssmin: T, ssmax: T, snr: T, csr: T, snl: T, csl: T };
}

pub fn svdTriu2x2(comptime T: type, f: T, g: T, h: T) Svd2(T) {
    const eps = std.math.floatEps(T);
    var ft = f;
    var fa = @abs(ft);
    var ht = h;
    var ha = @abs(h);
    var pmax: u2 = 1;
    const swap = ha > fa;
    if (swap) {
        pmax = 3;
        std.mem.swap(T, &ft, &ht);
        std.mem.swap(T, &fa, &ha);
    }
    const gt = g;
    const ga = @abs(gt);
    var ssmin: T = 0;
    var ssmax: T = 0;
    var clt: T = 0;
    var crt: T = 0;
    var slt: T = 0;
    var srt: T = 0;
    if (ga == 0) {
        ssmin = ha;
        ssmax = fa;
        clt = 1;
        crt = 1;
        slt = 0;
        srt = 0;
    } else {
        var gasmal = true;
        if (ga > fa) {
            pmax = 2;
            if (fa / ga < eps) {
                gasmal = false;
                ssmax = ga;
                ssmin = if (ha > 1) fa / (ga / ha) else (fa / ga) * ha;
                clt = 1;
                slt = ht / gt;
                srt = 1;
                crt = ft / gt;
            }
        }
        if (gasmal) {
            const d = fa - ha;
            const l: T = if (d == fa) 1 else d / fa;
            const m = gt / ft;
            var t = 2 - l;
            const mm = m * m;
            const tt = t * t;
            const ss = @sqrt(tt + mm);
            const r: T = if (l == 0) @abs(m) else @sqrt(l * l + mm);
            const a = 0.5 * (ss + r);
            ssmin = ha / a;
            ssmax = fa * a;
            if (mm == 0) {
                t = if (l == 0) copySign(T, 2, ft) * copySign(T, 1, gt) else gt / copySign(T, d, ft) + m / t;
            } else {
                t = (m / (ss + t) + m / (r + l)) * (1 + a);
            }
            const ll = @sqrt(t * t + 4);
            crt = 2 / ll;
            srt = t / ll;
            clt = (crt + srt * m) / a;
            slt = (ht / ft) * srt / a;
        }
    }
    const csl = if (swap) srt else clt;
    const snl = if (swap) crt else slt;
    const csr = if (swap) slt else crt;
    const snr = if (swap) clt else srt;
    const tsign: T = switch (pmax) {
        1 => copySign(T, 1, csr) * copySign(T, 1, csl) * copySign(T, 1, f),
        2 => copySign(T, 1, snr) * copySign(T, 1, csl) * copySign(T, 1, g),
        else => copySign(T, 1, snr) * copySign(T, 1, snl) * copySign(T, 1, h),
    };
    ssmax = copySign(T, ssmax, tsign);
    ssmin = copySign(T, ssmin, tsign * copySign(T, 1, f) * copySign(T, 1, h));
    return .{ .ssmin = ssmin, .ssmax = ssmax, .snr = snr, .csr = csr, .snl = snl, .csl = csl };
}

/// `xLAG2` (faer's `generalized_eigval_2x2`): the eigenvalues of the 2 x 2
/// pencil `(A, B)` (row-major `[a11, a12, a21, a22]`, `B` upper
/// triangular), `wr / scale` with an imaginary part `wi / scale`.
pub fn Eig2(comptime T: type) type {
    return struct { scale1: T, scale2: T, wr1: T, wr2: T, wi: T };
}

pub fn eig2x2(comptime T: type, av: [4]T, bv: [4]T) Eig2(T) {
    const safmin = std.math.floatMin(T) * 128;
    const safmax = 1 / safmin;
    const rtmin = @sqrt(safmin);
    const rtmax = 1 / rtmin;
    const half: T = 0.5;
    const anorm = @max(@max(safmin, @abs(av[0]) + @abs(av[2])), @abs(av[1]) + @abs(av[3]));
    const ascale = 1 / anorm;
    const a11 = av[0] * ascale;
    const a12 = av[1] * ascale;
    const a21 = av[2] * ascale;
    const a22 = av[3] * ascale;
    var b11 = bv[0];
    const b12 = bv[1];
    var b22 = bv[3];
    const bmin = rtmin * @max(@max(@max(@abs(b11), @abs(b12)), @abs(b22)), rtmin);
    if (@abs(b11) < bmin) b11 = copySign(T, bmin, b11);
    if (@abs(b22) < bmin) b22 = copySign(T, bmin, b22);
    const bnorm = @max(@max(safmin, @abs(b11)), @abs(b12) + @abs(b22));
    const bsize = @max(@abs(b11), @abs(b22));
    const bscale = 1 / bsize;
    const sb11 = b11 * bscale;
    const sb12 = b12 * bscale;
    const sb22 = b22 * bscale;
    const binv11 = 1 / sb11;
    const binv22 = 1 / sb22;
    const s1 = a11 * binv11;
    const s2 = a22 * binv22;
    var as12: T = undefined;
    var ss: T = undefined;
    var abi22: T = undefined;
    var pp: T = undefined;
    var shift: T = undefined;
    if (@abs(s1) <= @abs(s2)) {
        as12 = a12 - s1 * sb12;
        const as22 = a22 - s1 * sb22;
        ss = a21 * (binv11 * binv22);
        abi22 = as22 * binv22 - ss * sb12;
        pp = half * abi22;
        shift = s1;
    } else {
        as12 = a12 - s2 * sb12;
        const as11 = a11 - s2 * sb11;
        ss = a21 * (binv11 * binv22);
        abi22 = -ss * sb12;
        pp = half * (as11 * binv11 + abi22);
        shift = s2;
    }
    const qq = ss * as12;
    var discr: T = undefined;
    var r: T = undefined;
    if (@abs(pp * rtmin) >= 1) {
        discr = (pp * rtmin) * (pp * rtmin) + qq * safmin;
        r = @sqrt(@abs(discr)) * rtmax;
    } else if (pp * pp + @abs(qq) <= safmin) {
        discr = (rtmax * pp) * (rtmax * pp) + qq * safmax;
        r = @sqrt(@abs(discr)) * rtmin;
    } else {
        discr = pp * pp + qq;
        r = @sqrt(@abs(discr));
    }
    var wr1: T = undefined;
    var wr2: T = undefined;
    var wi: T = undefined;
    if (discr >= 0 or r == 0) {
        const sum = pp + copySign(T, r, pp);
        const diff = pp - copySign(T, r, pp);
        const wbig = shift + sum;
        var wsmall = shift + diff;
        if (@abs(wbig) * half > @max(safmin, @abs(wsmall))) {
            const wdet = (a11 * a22 - a12 * a21) * (binv11 * binv22);
            wsmall = wdet / wbig;
        }
        if (pp > abi22) {
            wr1 = @min(wbig, wsmall);
            wr2 = @max(wbig, wsmall);
        } else {
            wr1 = @max(wbig, wsmall);
            wr2 = @min(wbig, wsmall);
        }
        wi = 0;
    } else {
        wr1 = shift + pp;
        wr2 = wr1;
        wi = r;
    }
    const c1 = bsize * (safmin * @max(ascale, 1));
    const c2 = safmin * @max(bnorm, 1);
    const c3 = bsize * safmin;
    const c4: T = if (ascale <= 1 and bsize <= 1) @min(1, (ascale / safmin) * bsize) else 1;
    const c5: T = if (ascale <= 1 or bsize <= 1) @min(1, ascale * bsize) else 1;
    const fuzzy1: T = 1.00001;
    const wabs = @abs(wr1) + @abs(wi);
    const wsize = @max(safmin, @max(c1, @max(fuzzy1 * (wabs * c2 + c3), @min(c4, half * @max(c5, wabs)))));
    var scale1: T = undefined;
    var scale2: T = 1;
    if (wsize != 1) {
        const wscale = 1 / wsize;
        scale1 = if (wsize > 1) (@max(ascale, bsize) * wscale) * @min(ascale, bsize) else (@min(ascale, bsize) * wscale) * @max(ascale, bsize);
        wr1 *= wscale;
        if (wi != 0) {
            wi *= wscale;
            wr2 = wr1;
            scale2 = scale1;
        }
    } else {
        scale1 = ascale * bsize;
        scale2 = scale1;
    }
    if (wi == 0) {
        const wsize2 = @max(@max(safmin, c1), @max(fuzzy1 * (@abs(wr2) * c2 + c3), @min(c4, half * @max(@abs(wr2), c5))));
        if (wsize2 != 1) {
            const wscale = 1 / wsize2;
            scale2 = if (wsize2 > 1) (@max(ascale, bsize) * wscale) * @min(ascale, bsize) else (@min(ascale, bsize) * wscale) * @max(ascale, bsize);
            wr2 *= wscale;
        } else scale2 = ascale * bsize;
    }
    return .{ .scale1 = scale1, .scale2 = scale2, .wr1 = wr1, .wr2 = wr2, .wi = wi };
}
