//! The secular equation's roots, Reference-LAPACK's `xLAED4`, `xLAED5` and
//! `xLAED6` (`lapack-netlib/SRC/dlaed4.f`, `dlaed5.f`, `dlaed6.f`, OpenBLAS
//! v0.3.34-794-gdbf2c24c9; modified BSD), for `stebr.zig`. The `s` and `d`
//! routines are the same; indices are 0-based.

const std = @import("std");

/// `MAXIT` (`dlaed4.f`).
const maxit = 30;
/// `MAXIT` (`dlaed6.f`).
const maxit6 = 40;

/// `DLAED4`: the `i`-th root `λ` of `1 / rho + Σ z(j)² / (d(j) - λ)`, `d`
/// increasing, `rho > 0`, `‖z‖ = 1`; `delta(j) = d(j) - λ`, or for `n` 1
/// and 2 the root's unit eigenvector. `NoConvergence` after `MAXIT`
/// iterations (`INFO = 1`).
pub fn laed4(comptime F: type, i: usize, d: []const F, z: []const F, delta: []F, rho: F) error{NoConvergence}!F {
    const n = d.len;
    if (n == 1) {
        delta[0] = 1;
        return d[0] + rho * z[0] * z[0];
    }
    if (n == 2) return laed5(F, i, d[0..2].*, z[0..2].*, delta[0..2], rho);
    if (i == n - 1) return lastRoot(F, d, z, delta, rho);
    return innerRoot(F, i, d, z, delta, rho);
}

/// `ψ` over `j < ii`, `φ` over `j > ii` (descending), their derivatives,
/// and the error bound's partial sum (`dlaed4.f` loops 150 and 160).
fn Sums(comptime F: type) type {
    return struct {
        psi: F = 0,
        dpsi: F = 0,
        phi: F = 0,
        dphi: F = 0,
        erretm: F = 0,

        fn of(z: []const F, delta: []const F, ii: usize) @This() {
            var s: @This() = .{};
            for (0..ii) |j| {
                const temp = z[j] / delta[j];
                s.psi += z[j] * temp;
                s.dpsi += temp * temp;
                s.erretm += s.psi;
            }
            s.erretm = @abs(s.erretm);
            var j = delta.len;
            while (j > ii + 1) {
                j -= 1;
                const temp = z[j] / delta[j];
                s.phi += z[j] * temp;
                s.dphi += temp * temp;
                s.erretm += s.phi;
            }
            return s;
        }
    };
}

/// `I = N` (`dlaed4.f`, the first branch): origin `d(n - 1)`.
fn lastRoot(comptime F: type, d: []const F, z: []const F, delta: []F, rho: F) error{NoConvergence}!F {
    const n = d.len;
    const eps = std.math.floatEps(F) / 2; // `DLAMCH('Epsilon')`
    const rhoinv = 1 / rho;
    const ii = n - 2;
    const dn = d[n - 1];
    const midpt = rho / 2;
    for (delta, d) |*x, dj| x.* = (dj - dn) - midpt;
    var psi: F = 0;
    for (0..n - 2) |j| psi += z[j] * z[j] / delta[j];
    var c = rhoinv + psi;
    var w = c + z[ii] * z[ii] / delta[ii] + z[n - 1] * z[n - 1] / delta[n - 1];
    var tau: F = undefined;
    var dltlb: F = undefined;
    var dltub: F = undefined;
    if (w <= 0) {
        const temp = z[n - 2] * z[n - 2] / (d[n - 1] - d[n - 2] + rho) + z[n - 1] * z[n - 1] / rho;
        tau = if (c <= temp) rho else lastGuess(F, d, z, c);
        // D(N) + RHO/2 <= LAMBDA(N) < D(N) + TAU <= D(N) + RHO
        dltlb = midpt;
        dltub = rho;
    } else {
        tau = lastGuess(F, d, z, c);
        // D(N) < D(N) + TAU < LAMBDA(N) < D(N) + RHO/2
        dltlb = 0;
        dltub = midpt;
    }
    for (delta, d) |*x, dj| x.* = (dj - dn) - tau;
    var s = lastSums(F, z, delta, rhoinv, tau);
    w = s.w;
    if (@abs(w) <= eps * s.erretm) return dn + tau;
    if (w <= 0) dltlb = @max(dltlb, tau) else dltub = @min(dltub, tau);
    // the first step
    c = w - delta[n - 2] * s.dpsi - delta[n - 1] * s.dphi;
    var a = (delta[n - 2] + delta[n - 1]) * w - delta[n - 2] * delta[n - 1] * (s.dpsi + s.dphi);
    var b = delta[n - 2] * delta[n - 1] * w;
    if (c < 0) c = @abs(c);
    var eta: F = if (c == 0)
        // Li Ren-Cang's update
        -w / (s.dpsi + s.dphi)
    else if (a >= 0)
        (a + @sqrt(@abs(a * a - 4 * b * c))) / (2 * c)
    else
        2 * b / (a - @sqrt(@abs(a * a - 4 * b * c)));
    // roundoff can give eta w > 0: one Newton step instead
    if (w * eta > 0) eta = -w / (s.dpsi + s.dphi);
    if (tau + eta > dltub or tau + eta < dltlb) {
        const eta1 = -w / (s.dpsi + s.dphi);
        const temp = tau + eta1;
        const eta2 = if (w < 0) (dltub - tau) / 2 else (dltlb - tau) / 2;
        // a Newton step within the bounds: the geometric mean, its sign
        eta = if (dltlb <= temp and temp <= dltub) std.math.copysign(@as(F, 1), eta1) * @sqrt(@abs(eta1)) * @sqrt(@abs(eta2)) else eta2;
    }
    for (delta) |*x| x.* -= eta;
    tau += eta;
    s = lastSums(F, z, delta, rhoinv, tau);
    w = s.w;
    for (3..maxit + 1) |_| {
        if (@abs(w) <= eps * s.erretm) return dn + tau;
        if (w <= 0) dltlb = @max(dltlb, tau) else dltub = @min(dltub, tau);
        c = w - delta[n - 2] * s.dpsi - delta[n - 1] * s.dphi;
        a = (delta[n - 2] + delta[n - 1]) * w - delta[n - 2] * delta[n - 1] * (s.dpsi + s.dphi);
        b = delta[n - 2] * delta[n - 1] * w;
        eta = if (a >= 0)
            (a + @sqrt(@abs(a * a - 4 * b * c))) / (2 * c)
        else
            2 * b / (a - @sqrt(@abs(a * a - 4 * b * c)));
        if (w * eta > 0) eta = -w / (s.dpsi + s.dphi);
        if (tau + eta > dltub or tau + eta < dltlb) {
            eta = if (w < 0) (dltub - tau) / 2 else (dltlb - tau) / 2;
        }
        for (delta) |*x| x.* -= eta;
        tau += eta;
        s = lastSums(F, z, delta, rhoinv, tau);
        w = s.w;
    }
    return error.NoConvergence;
}

/// `dlaed4.f`'s initial `TAU` for `I = N` from the last two poles.
fn lastGuess(comptime F: type, d: []const F, z: []const F, c: F) F {
    const n = d.len;
    const del = d[n - 1] - d[n - 2];
    const a = -c * del + z[n - 2] * z[n - 2] + z[n - 1] * z[n - 1];
    const b = z[n - 1] * z[n - 1] * del;
    return if (a < 0) 2 * b / (@sqrt(a * a + 4 * b * c) - a) else (a + @sqrt(a * a + 4 * b * c)) / (2 * c);
}

/// `dlaed4.f` loops 40, 60 and 80: `ψ` over the first `n - 1`, `φ` the
/// last, the function `w` and its error bound.
fn lastSums(comptime F: type, z: []const F, delta: []const F, rhoinv: F, tau: F) struct { dpsi: F, dphi: F, w: F, erretm: F } {
    const n = delta.len;
    const s = Sums(F).of(z, delta, n - 1);
    const temp = z[n - 1] / delta[n - 1];
    const phi = z[n - 1] * temp;
    const dphi = temp * temp;
    return .{
        .dpsi = s.dpsi,
        .dphi = dphi,
        .w = rhoinv + phi + s.psi,
        .erretm = 8 * (-phi - s.psi) + s.erretm - phi + rhoinv + @abs(tau) * (s.dpsi + dphi),
    };
}

/// `I < N` (`dlaed4.f`, the second branch): origin `d(i)` or `d(i + 1)`,
/// whichever is nearer.
fn innerRoot(comptime F: type, i: usize, d: []const F, z: []const F, delta: []F, rho: F) error{NoConvergence}!F {
    const n = d.len;
    const eps = std.math.floatEps(F) / 2; // `DLAMCH('Epsilon')`
    const rhoinv = 1 / rho;
    const ip1 = i + 1;
    const del = d[ip1] - d[i];
    const midpt = del / 2;
    for (delta, d) |*x, dj| x.* = (dj - d[i]) - midpt;
    var psi: F = 0;
    for (0..i) |j| psi += z[j] * z[j] / delta[j];
    var phi: F = 0;
    var jj = n;
    while (jj > i + 2) {
        jj -= 1;
        phi += z[jj] * z[jj] / delta[jj];
    }
    var c = rhoinv + psi + phi;
    var w = c + z[i] * z[i] / delta[i] + z[ip1] * z[ip1] / delta[ip1];
    var tau: F = undefined;
    var dltlb: F = undefined;
    var dltub: F = undefined;
    // w > 0: d(i) < the root < (d(i) + d(i + 1)) / 2, origin d(i)
    const orgati = w > 0;
    if (orgati) {
        const a = c * del + z[i] * z[i] + z[ip1] * z[ip1];
        const b = z[i] * z[i] * del;
        tau = if (a > 0) 2 * b / (a + @sqrt(@abs(a * a - 4 * b * c))) else (a - @sqrt(@abs(a * a - 4 * b * c))) / (2 * c);
        dltlb = 0;
        dltub = midpt;
    } else {
        const a = c * del - z[i] * z[i] - z[ip1] * z[ip1];
        const b = z[ip1] * z[ip1] * del;
        tau = if (a < 0) 2 * b / (a - @sqrt(@abs(a * a + 4 * b * c))) else -(a + @sqrt(@abs(a * a + 4 * b * c))) / (2 * c);
        dltlb = -midpt;
        dltub = 0;
    }
    const ii = if (orgati) i else ip1;
    const org = d[ii];
    for (delta, d) |*x, dj| x.* = (dj - org) - tau;
    var s = Sums(F).of(z, delta, ii);
    w = rhoinv + s.phi + s.psi;
    // w is the function without its ii-th term
    const swtch3 = ii != 0 and ii != n - 1 and (if (orgati) w < 0 else w > 0);
    var temp = z[ii] / delta[ii];
    var dw = s.dpsi + s.dphi + temp * temp;
    temp = z[ii] * temp;
    w += temp;
    var erretm = 8 * (s.phi - s.psi) + s.erretm + 2 * rhoinv + 3 * @abs(temp) + @abs(tau) * dw;
    if (@abs(w) <= eps * erretm) return org + tau;
    if (w <= 0) dltlb = @max(dltlb, tau) else dltub = @min(dltub, tau);
    // the first step
    var eta: F = undefined;
    if (!swtch3) {
        c = if (orgati)
            w - delta[ip1] * dw - (d[i] - d[ip1]) * sq(F, z[i] / delta[i])
        else
            w - delta[i] * dw - (d[ip1] - d[i]) * sq(F, z[ip1] / delta[ip1]);
        var a = (delta[i] + delta[ip1]) * w - delta[i] * delta[ip1] * dw;
        const b = delta[i] * delta[ip1] * w;
        if (c == 0) {
            if (a == 0) a = if (orgati)
                z[i] * z[i] + delta[ip1] * delta[ip1] * (s.dpsi + s.dphi)
            else
                z[ip1] * z[ip1] + delta[i] * delta[i] * (s.dpsi + s.dphi);
            eta = b / a;
        } else eta = quadStep(F, a, b, c);
    } else {
        // interpolation by the three nearest poles
        eta = try laed6(F, true, orgati, threePoleRho(F, rhoinv, s, d, z, delta, ii, orgati, false), delta[ii - 1 ..][0..3].*, threePoleZ(F, s, z, delta, ii, orgati, false), w);
    }
    if (w * eta >= 0) eta = -w / dw;
    if (tau + eta > dltub or tau + eta < dltlb) {
        eta = if (w < 0) (dltub - tau) / 2 else (dltlb - tau) / 2;
    }
    var prew = w;
    for (delta) |*x| x.* -= eta;
    s = Sums(F).of(z, delta, ii);
    temp = z[ii] / delta[ii];
    dw = s.dpsi + s.dphi + temp * temp;
    temp = z[ii] * temp;
    w = rhoinv + s.phi + s.psi + temp;
    erretm = 8 * (s.phi - s.psi) + s.erretm + 2 * rhoinv + 3 * @abs(temp) + @abs(tau + eta) * dw;
    var swtch = if (orgati) -w > @abs(prew) / 10 else w > @abs(prew) / 10;
    tau += eta;
    for (3..maxit + 1) |_| {
        if (@abs(w) <= eps * erretm) return org + tau;
        if (w <= 0) dltlb = @max(dltlb, tau) else dltub = @min(dltub, tau);
        if (!swtch3) {
            if (!swtch) {
                c = if (orgati)
                    w - delta[ip1] * dw - (d[i] - d[ip1]) * sq(F, z[i] / delta[i])
                else
                    w - delta[i] * dw - (d[ip1] - d[i]) * sq(F, z[ip1] / delta[ip1]);
            } else {
                const t = z[ii] / delta[ii];
                if (orgati) s.dpsi += t * t else s.dphi += t * t;
                c = w - delta[i] * s.dpsi - delta[ip1] * s.dphi;
            }
            var a = (delta[i] + delta[ip1]) * w - delta[i] * delta[ip1] * dw;
            const b = delta[i] * delta[ip1] * w;
            if (c == 0) {
                if (a == 0) a = if (swtch)
                    delta[i] * delta[i] * s.dpsi + delta[ip1] * delta[ip1] * s.dphi
                else if (orgati)
                    z[i] * z[i] + delta[ip1] * delta[ip1] * (s.dpsi + s.dphi)
                else
                    z[ip1] * z[ip1] + delta[i] * delta[i] * (s.dpsi + s.dphi);
                eta = b / a;
            } else eta = quadStep(F, a, b, c);
        } else {
            eta = try laed6(F, false, orgati, threePoleRho(F, rhoinv, s, d, z, delta, ii, orgati, swtch), delta[ii - 1 ..][0..3].*, threePoleZ(F, s, z, delta, ii, orgati, swtch), w);
        }
        if (w * eta >= 0) eta = -w / dw;
        if (tau + eta > dltub or tau + eta < dltlb) {
            const eta1 = -w / dw;
            const t = tau + eta1;
            const eta2 = if (w < 0) (dltub - tau) / 2 else (dltlb - tau) / 2;
            // a Newton step within the bounds: the geometric mean, its sign
            eta = if (dltlb <= t and t <= dltub) std.math.copysign(@as(F, 1), eta1) * @sqrt(@abs(eta1)) * @sqrt(@abs(eta2)) else eta2;
        }
        for (delta) |*x| x.* -= eta;
        tau += eta;
        prew = w;
        s = Sums(F).of(z, delta, ii);
        temp = z[ii] / delta[ii];
        dw = s.dpsi + s.dphi + temp * temp;
        temp = z[ii] * temp;
        w = rhoinv + s.phi + s.psi + temp;
        erretm = 8 * (s.phi - s.psi) + s.erretm + 2 * rhoinv + 3 * @abs(temp) + @abs(tau) * dw;
        if (w * prew > 0 and @abs(w) > @abs(prew) / 10) swtch = !swtch;
    }
    return error.NoConvergence;
}

fn sq(comptime F: type, x: F) F {
    return x * x;
}

/// `dlaed4.f`'s step from `a`, `b`, `c` for `I < N`, `c` nonzero.
fn quadStep(comptime F: type, a: F, b: F, c: F) F {
    return if (a <= 0) (a - @sqrt(@abs(a * a - 4 * b * c))) / (2 * c) else 2 * b / (a + @sqrt(@abs(a * a - 4 * b * c)));
}

/// The three-pole interpolation's `C` (`DLAED6`'s `RHO`), `dlaed4.f`.
fn threePoleRho(comptime F: type, rhoinv: F, s: Sums(F), d: []const F, z: []const F, delta: []const F, ii: usize, orgati: bool, swtch: bool) F {
    const temp = rhoinv + s.psi + s.phi;
    if (swtch) return temp - delta[ii - 1] * s.dpsi - delta[ii + 1] * s.dphi;
    if (orgati) {
        const t1 = sq(F, z[ii - 1] / delta[ii - 1]);
        return temp - delta[ii + 1] * (s.dpsi + s.dphi) - (d[ii - 1] - d[ii + 1]) * t1;
    }
    const t1 = sq(F, z[ii + 1] / delta[ii + 1]);
    return temp - delta[ii - 1] * (s.dpsi + s.dphi) - (d[ii + 1] - d[ii - 1]) * t1;
}

/// The three-pole interpolation's `ZZ`, `dlaed4.f`.
fn threePoleZ(comptime F: type, s: Sums(F), z: []const F, delta: []const F, ii: usize, orgati: bool, swtch: bool) [3]F {
    const mid = z[ii] * z[ii];
    if (swtch) return .{ delta[ii - 1] * delta[ii - 1] * s.dpsi, mid, delta[ii + 1] * delta[ii + 1] * s.dphi };
    if (orgati) {
        const t1 = sq(F, z[ii - 1] / delta[ii - 1]);
        return .{ z[ii - 1] * z[ii - 1], mid, delta[ii + 1] * delta[ii + 1] * ((s.dpsi - t1) + s.dphi) };
    }
    const t1 = sq(F, z[ii + 1] / delta[ii + 1]);
    return .{ delta[ii - 1] * delta[ii - 1] * (s.dpsi + (s.dphi - t1)), mid, z[ii + 1] * z[ii + 1] };
}

/// `DLAED5`: the `i`-th (0 or 1) root of the 2 x 2 secular equation,
/// `delta` its unit eigenvector.
pub fn laed5(comptime F: type, i: usize, d: [2]F, z: [2]F, delta: *[2]F, rho: F) F {
    const del = d[1] - d[0];
    var dlam: F = undefined;
    if (i == 0) {
        const w = 1 + 2 * rho * (z[1] * z[1] - z[0] * z[0]) / del;
        if (w > 0) {
            const b = del + rho * (z[0] * z[0] + z[1] * z[1]);
            const c = rho * z[0] * z[0] * del;
            // b > 0, always
            const tau = 2 * c / (b + @sqrt(@abs(b * b - 4 * c)));
            dlam = d[0] + tau;
            delta[0] = -z[0] / tau;
            delta[1] = z[1] / (del - tau);
        } else {
            const b = -del + rho * (z[0] * z[0] + z[1] * z[1]);
            const c = rho * z[1] * z[1] * del;
            const tau = if (b > 0) -2 * c / (b + @sqrt(b * b + 4 * c)) else (b - @sqrt(b * b + 4 * c)) / 2;
            dlam = d[1] + tau;
            delta[0] = -z[0] / (del + tau);
            delta[1] = -z[1] / tau;
        }
    } else {
        const b = -del + rho * (z[0] * z[0] + z[1] * z[1]);
        const c = rho * z[1] * z[1] * del;
        const tau = if (b > 0) (b + @sqrt(b * b + 4 * c)) / 2 else 2 * c / (-b + @sqrt(b * b + 4 * c));
        dlam = d[1] + tau;
        delta[0] = -z[0] / (del + tau);
        delta[1] = -z[1] / tau;
    }
    const temp = @sqrt(delta[0] * delta[0] + delta[1] * delta[1]);
    delta[0] /= temp;
    delta[1] /= temp;
    return dlam;
}

/// `DLAED6`: the root `tau` nearest the origin of `rho + Σ z(i) / (d(i) -
/// tau)`, three poles, between `d(1)` and `d(2)` (`orgati`) or `d(0)` and
/// `d(1)`; `finit` the function at 0, `first` for `dlaed4.f`'s second
/// iteration (`KNITER = 2`). Gragg, Thornton and Warner's cubically
/// convergent scheme. `NoConvergence` after `MAXIT` (`INFO = 1`).
pub fn laed6(comptime F: type, first: bool, orgati: bool, rho: F, d: [3]F, z: [3]F, finit: F) error{NoConvergence}!F {
    var lbd = if (orgati) d[1] else d[0];
    var ubd = if (orgati) d[2] else d[1];
    if (finit < 0) lbd = 0 else ubd = 0;
    var tau: F = 0;
    if (first) {
        var a: F = undefined;
        var b: F = undefined;
        var c: F = undefined;
        if (orgati) {
            const temp = (d[2] - d[1]) / 2;
            c = rho + z[0] / ((d[0] - d[1]) - temp);
            a = c * (d[1] + d[2]) + z[1] + z[2];
            b = c * d[1] * d[2] + z[1] * d[2] + z[2] * d[1];
        } else {
            const temp = (d[0] - d[1]) / 2;
            c = rho + z[2] / ((d[2] - d[1]) - temp);
            a = c * (d[0] + d[1]) + z[0] + z[1];
            b = c * d[0] * d[1] + z[0] * d[1] + z[1] * d[0];
        }
        const temp = @max(@max(@abs(a), @abs(b)), @abs(c));
        a /= temp;
        b /= temp;
        c /= temp;
        tau = if (c == 0) b / a else quadStep(F, a, b, c);
        if (tau < lbd or tau > ubd) tau = (lbd + ubd) / 2;
        if (d[0] == tau or d[1] == tau or d[2] == tau) {
            tau = 0;
        } else {
            const t = finit + tau * z[0] / (d[0] * (d[0] - tau)) + tau * z[1] / (d[1] * (d[1] - tau)) + tau * z[2] / (d[2] * (d[2] - tau));
            if (t <= 0) lbd = tau else ubd = tau;
            if (@abs(finit) <= @abs(t)) tau = 0;
        }
    }
    const eps = std.math.floatEps(F) / 2; // `DLAMCH('Epsilon')`
    // `BASE ** INT(LOG(SAFMIN) / LOG(BASE) / 3)`, recomputed per call there
    const small1 = comptime std.math.ldexp(@as(F, 1), @divTrunc(std.math.floatExponentMin(F), 3));
    const sminv1 = 1 / small1;
    const small2 = small1 * small1;
    const sminv2 = sminv1 * sminv1;
    // scale against overflow in 1 / temp³
    const near = if (orgati) @min(@abs(d[1] - tau), @abs(d[2] - tau)) else @min(@abs(d[0] - tau), @abs(d[1] - tau));
    var dscale = d;
    var zscale = z;
    var sclinv: ?F = null;
    if (near <= small1) {
        // by the radix power nearest 1 / SAFMIN^(2/3), or ^(1/3); safe, as
        // d, z and tau are O(1)
        const sclfac = if (near <= small2) sminv2 else sminv1;
        sclinv = if (near <= small2) small2 else small1;
        for (&dscale, &zscale) |*x, *y| {
            x.* *= sclfac;
            y.* *= sclfac;
        }
        tau *= sclfac;
        lbd *= sclfac;
        ubd *= sclfac;
    }
    var fc: F = 0;
    var df: F = 0;
    var ddf: F = 0;
    for (dscale, zscale) |x, y| {
        const temp = 1 / (x - tau);
        const temp1 = y * temp;
        const temp2 = temp1 * temp;
        const temp3 = temp2 * temp;
        fc += temp1 / x;
        df += temp2;
        ddf += temp3;
    }
    var f = finit + tau * fc;
    done: {
        if (@abs(f) <= 0) break :done;
        if (f <= 0) lbd = tau else ubd = tau;
        // monotone: up for finit < 0, down for finit > 0
        for (2..maxit6 + 1) |_| {
            const temp1 = if (orgati) dscale[1] - tau else dscale[0] - tau;
            const temp2 = if (orgati) dscale[2] - tau else dscale[1] - tau;
            var a = (temp1 + temp2) * f - temp1 * temp2 * df;
            var b = temp1 * temp2 * f;
            var c = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
            const temp = @max(@max(@abs(a), @abs(b)), @abs(c));
            a /= temp;
            b /= temp;
            c /= temp;
            var eta = if (c == 0) b / a else quadStep(F, a, b, c);
            if (f * eta >= 0) eta = -f / df;
            tau += eta;
            if (tau < lbd or tau > ubd) tau = (lbd + ubd) / 2;
            fc = 0;
            var erretm: F = 0;
            df = 0;
            ddf = 0;
            for (dscale, zscale) |x, y| {
                if (x - tau == 0) break :done;
                const t = 1 / (x - tau);
                const t1 = y * t;
                const t2 = t1 * t;
                const t3 = t2 * t;
                const t4 = t1 / x;
                fc += t4;
                erretm += @abs(t4);
                df += t2;
                ddf += t3;
            }
            f = finit + tau * fc;
            erretm = 8 * (@abs(finit) + @abs(tau) * erretm) + @abs(tau) * df;
            if (@abs(f) <= 4 * eps * erretm or (ubd - lbd) <= 4 * eps * @abs(tau)) break :done;
            if (f <= 0) lbd = tau else ubd = tau;
        }
        return error.NoConvergence;
    }
    if (sclinv) |s| tau *= s;
    return tau;
}
