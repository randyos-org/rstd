//! dqds's step, `xLASQ3` (deflation, then a transform with a shift, failed
//! shifts retried), and its shift, `xLASQ4`, Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dlasq3.f`, `dlasq4.f`, OpenBLAS v0.3.34-317;
//! modified BSD).

const std = @import("std");
const dqds = @import("lasq_dqds.zig");

const Qd = dqds.Qd;
const fmin = dqds.fmin;

/// `CBIAS` (`dlasq2.f`, `dlasq3.f`): how lopsided a block may be before it
/// is flipped.
pub const cbias = 1.5;

/// `DLASQ4` (`dlasq4.f`)'s constants, as written there (`THIRD` is 0.333).
const cnst1 = 0.563;
const cnst2 = 1.01;
const cnst3 = 1.05;
const qurtr = 0.25;
const third = 0.333;
const half = 0.5;
const hundrd = 100;

/// `DLASQ4` (`dlasq4.f`): the shift `tau` for the next transform and its
/// type `ttype`, from the minimal `d`s; `n0in` the block's end before
/// deflation. An early return leaves `tau` as it was, as the reference's.
pub fn lasq4(comptime F: type, s: *Qd(F), n0in: isize) void {
    const z = struct {
        fn f(q: *Qd(F), i: isize) F {
            return q.at(i).*;
        }
    }.f;
    const lo = s.lo;
    const n0 = s.n0;
    const pp = s.pp;
    if (s.dmin <= 0) {
        s.tau = -s.dmin;
        s.ttype = -1;
        return;
    }
    const nn = 4 * n0 + pp;
    var sh: F = 0;
    if (n0in == n0) {
        if (s.dmin == s.dn or s.dmin == s.dn1) {
            var b1 = @sqrt(z(s, nn - 3)) * @sqrt(z(s, nn - 5));
            var b2 = @sqrt(z(s, nn - 7)) * @sqrt(z(s, nn - 9));
            var a2 = z(s, nn - 7) + z(s, nn - 5);
            if (s.dmin == s.dn and s.dmin1 == s.dn1) {
                // cases 2 and 3
                const gap2 = s.dmin2 - a2 - s.dmin2 * qurtr;
                const gap1 = if (gap2 > 0 and gap2 > b2) a2 - s.dn - (b2 / gap2) * b2 else a2 - s.dn - (b1 + b2);
                if (gap1 > 0 and gap1 > b1) {
                    sh = @max(s.dn - (b1 / gap1) * b1, half * s.dmin);
                    s.ttype = -2;
                } else {
                    sh = 0;
                    if (s.dn > b1) sh = s.dn - b1;
                    if (a2 > b1 + b2) sh = @min(sh, a2 - (b1 + b2));
                    sh = @max(sh, third * s.dmin);
                    s.ttype = -3;
                }
            } else {
                // case 4
                s.ttype = -4;
                sh = qurtr * s.dmin;
                var gam: F = undefined;
                var np: isize = undefined;
                if (s.dmin == s.dn) {
                    gam = s.dn;
                    a2 = 0;
                    if (z(s, nn - 5) > z(s, nn - 7)) return;
                    b2 = z(s, nn - 5) / z(s, nn - 7);
                    np = nn - 9;
                } else {
                    np = nn - 2 * pp;
                    gam = s.dn1;
                    if (z(s, np - 4) > z(s, np - 2)) return;
                    a2 = z(s, np - 4) / z(s, np - 2);
                    if (z(s, nn - 9) > z(s, nn - 11)) return;
                    b2 = z(s, nn - 9) / z(s, nn - 11);
                    np = nn - 13;
                }
                // the approximation of the nearby eigenvalue
                a2 += b2;
                var k4 = np;
                while (k4 >= 4 * lo - 1 + pp) : (k4 -= 4) {
                    if (b2 == 0) break;
                    b1 = b2;
                    if (z(s, k4) > z(s, k4 - 2)) return;
                    b2 *= z(s, k4) / z(s, k4 - 2);
                    a2 += b2;
                    if (hundrd * @max(b2, b1) < a2 or cnst1 < a2) break;
                }
                a2 *= cnst3;
                if (a2 < cnst1) sh = gam * (1 - @sqrt(a2)) / (1 + a2);
            }
        } else if (s.dmin == s.dn2) {
            // case 5
            s.ttype = -5;
            sh = qurtr * s.dmin;
            const np = nn - 2 * pp;
            var b1 = z(s, np - 2);
            var b2 = z(s, np - 6);
            const gam = s.dn2;
            if (z(s, np - 8) > b2 or z(s, np - 4) > b1) return;
            var a2 = (z(s, np - 8) / b2) * (1 + z(s, np - 4) / b1);
            if (n0 - lo > 2) {
                b2 = z(s, nn - 13) / z(s, nn - 15);
                a2 += b2;
                var k4 = nn - 17;
                while (k4 >= 4 * lo - 1 + pp) : (k4 -= 4) {
                    if (b2 == 0) break;
                    b1 = b2;
                    if (z(s, k4) > z(s, k4 - 2)) return;
                    b2 *= z(s, k4) / z(s, k4 - 2);
                    a2 += b2;
                    if (hundrd * @max(b2, b1) < a2 or cnst1 < a2) break;
                }
                a2 *= cnst3;
            }
            if (a2 < cnst1) sh = gam * (1 - @sqrt(a2)) / (1 + a2);
        } else {
            // case 6, no information to guide us
            if (s.ttype == -6) {
                s.g += third * (1 - s.g);
            } else if (s.ttype == -18) {
                s.g = qurtr * third;
            } else s.g = qurtr;
            sh = s.g * s.dmin;
            s.ttype = -6;
        }
    } else if (n0in == n0 + 1) {
        // one eigenvalue just deflated
        if (s.dmin1 == s.dn1 and s.dmin2 == s.dn2) {
            // case 7
            s.ttype = -7;
            sh = third * s.dmin1;
            if (z(s, nn - 5) > z(s, nn - 7)) return;
            var b1 = z(s, nn - 5) / z(s, nn - 7);
            var b2 = b1;
            if (b2 != 0) {
                var k4 = 4 * n0 - 9 + pp;
                while (k4 >= 4 * lo - 1 + pp) : (k4 -= 4) {
                    const a2 = b1;
                    if (z(s, k4) > z(s, k4 - 2)) return;
                    b1 *= z(s, k4) / z(s, k4 - 2);
                    b2 += b1;
                    if (hundrd * @max(b1, a2) < b2) break;
                }
            }
            b2 = @sqrt(cnst3 * b2);
            const a2 = s.dmin1 / (1 + b2 * b2);
            const gap2 = half * s.dmin2 - a2;
            if (gap2 > 0 and gap2 > b2 * a2) {
                sh = @max(sh, a2 * (1 - cnst2 * a2 * (b2 / gap2) * b2));
            } else {
                sh = @max(sh, a2 * (1 - cnst2 * b2));
                s.ttype = -8;
            }
        } else {
            // case 9
            sh = qurtr * s.dmin1;
            if (s.dmin1 == s.dn1) sh = half * s.dmin1;
            s.ttype = -9;
        }
    } else if (n0in == n0 + 2) {
        // two eigenvalues deflated: use dmin2, dn2 for dmin and dn
        if (s.dmin2 == s.dn2 and 2 * z(s, nn - 5) < z(s, nn - 7)) {
            // case 10
            s.ttype = -10;
            sh = third * s.dmin2;
            if (z(s, nn - 5) > z(s, nn - 7)) return;
            var b1 = z(s, nn - 5) / z(s, nn - 7);
            var b2 = b1;
            if (b2 != 0) {
                var k4 = 4 * n0 - 9 + pp;
                while (k4 >= 4 * lo - 1 + pp) : (k4 -= 4) {
                    if (z(s, k4) > z(s, k4 - 2)) return;
                    b1 *= z(s, k4) / z(s, k4 - 2);
                    b2 += b1;
                    if (hundrd * b1 < b2) break;
                }
            }
            b2 = @sqrt(cnst3 * b2);
            const a2 = s.dmin2 / (1 + b2 * b2);
            const gap2 = z(s, nn - 7) + z(s, nn - 9) - @sqrt(z(s, nn - 11)) * @sqrt(z(s, nn - 9)) - a2;
            if (gap2 > 0 and gap2 > b2 * a2) {
                sh = @max(sh, a2 * (1 - cnst2 * a2 * (b2 / gap2) * b2));
            } else sh = @max(sh, a2 * (1 - cnst2 * b2));
        } else {
            sh = qurtr * s.dmin2;
            s.ttype = -11;
        }
    } else if (n0in > n0 + 2) {
        // case 12, more than two eigenvalues deflated: no information
        sh = 0;
        s.ttype = -12;
    }
    s.tau = sh;
}

/// `DLASQ3` (`dlasq3.f`): deflates the block's converged end (one or two
/// eigenvalues at a time), flips a lopsided block, then runs one dqds
/// transform with `lasq4`'s shift, retrying a failed shift smaller, or
/// `lasq6` (no shift) as the last resort; `sigma` gains the shift.
pub fn lasq3(comptime F: type, s: *Qd(F)) void {
    const z = struct {
        fn p(q: *Qd(F), i: isize) *F {
            return q.at(i);
        }
    }.p;
    const n0in = s.n0;
    const eps = std.math.floatEps(F);
    const tol = eps * hundrd;
    const tol2 = tol * tol;
    // deflation, `10` to `40`
    while (true) {
        if (s.n0 < s.lo) return;
        var two = false;
        if (s.n0 != s.lo) {
            const nn = 4 * s.n0 + s.pp;
            if (s.n0 == s.lo + 1) {
                two = true;
            } else if (z(s, nn - 5).* > tol2 * (s.sigma + z(s, nn - 3).*) and z(s, nn - 2 * s.pp - 4).* > tol2 * z(s, nn - 7).*) {
                if (z(s, nn - 9).* > tol2 * s.sigma and z(s, nn - 2 * s.pp - 8).* > tol2 * z(s, nn - 11).*) break;
                two = true;
            }
        }
        if (!two) {
            // one eigenvalue
            z(s, 4 * s.n0 - 3).* = z(s, 4 * s.n0 + s.pp - 3).* + s.sigma;
            s.n0 -= 1;
            continue;
        }
        // two eigenvalues
        const nn = 4 * s.n0 + s.pp;
        if (z(s, nn - 3).* > z(s, nn - 7).*) std.mem.swap(F, z(s, nn - 3), z(s, nn - 7));
        var t = half * ((z(s, nn - 7).* - z(s, nn - 3).*) + z(s, nn - 5).*);
        if (z(s, nn - 5).* > z(s, nn - 3).* * tol2 and t != 0) {
            var sv = z(s, nn - 3).* * (z(s, nn - 5).* / t);
            if (sv <= t) {
                sv = z(s, nn - 3).* * (z(s, nn - 5).* / (t * (1 + @sqrt(1 + sv / t))));
            } else sv = z(s, nn - 3).* * (z(s, nn - 5).* / (t + @sqrt(t) * @sqrt(t + sv)));
            t = z(s, nn - 7).* + (sv + z(s, nn - 5).*);
            z(s, nn - 3).* *= z(s, nn - 7).* / t;
            z(s, nn - 7).* = t;
        }
        z(s, 4 * s.n0 - 7).* = z(s, nn - 7).* + s.sigma;
        z(s, 4 * s.n0 - 3).* = z(s, nn - 3).* + s.sigma;
        s.n0 -= 2;
    }
    // `50`: reverse the qd array if warranted
    if (s.pp == 2) s.pp = 0;
    const lo = s.lo;
    const n0 = s.n0;
    const pp = s.pp;
    if (s.dmin <= 0 or n0 < n0in) {
        if (cbias * z(s, 4 * lo + pp - 3).* < z(s, 4 * n0 + pp - 3).*) {
            const ipn4 = 4 * (lo + n0);
            var j4 = 4 * lo;
            while (j4 <= 2 * (lo + n0 - 1)) : (j4 += 4) {
                std.mem.swap(F, z(s, j4 - 3), z(s, ipn4 - j4 - 3));
                std.mem.swap(F, z(s, j4 - 2), z(s, ipn4 - j4 - 2));
                std.mem.swap(F, z(s, j4 - 1), z(s, ipn4 - j4 - 5));
                std.mem.swap(F, z(s, j4), z(s, ipn4 - j4 - 4));
            }
            if (n0 - lo <= 4) {
                z(s, 4 * n0 + pp - 1).* = z(s, 4 * lo + pp - 1).*;
                z(s, 4 * n0 - pp).* = z(s, 4 * lo - pp).*;
            }
            s.dmin2 = @min(s.dmin2, z(s, 4 * n0 + pp - 1).*);
            z(s, 4 * n0 + pp - 1).* = @min(z(s, 4 * n0 + pp - 1).*, z(s, 4 * lo + pp - 1).*, z(s, 4 * lo + pp + 3).*);
            z(s, 4 * n0 - pp).* = @min(z(s, 4 * n0 - pp).*, z(s, 4 * lo - pp).*, z(s, 4 * lo - pp + 4).*);
            s.qmax = @max(s.qmax, z(s, 4 * lo + pp - 3).*, z(s, 4 * lo + pp + 1).*);
            s.dmin = -0.0;
        }
    }
    // the shift, then the transform, `70` to `90`
    lasq4(F, s, n0in);
    const unshifted = while (true) {
        dqds.lasq5(F, s, eps);
        s.ndiv += n0 - lo + 2;
        s.iter += 1;
        if (s.dmin >= 0 and s.dmin1 >= 0) break false;
        if (s.dmin < 0 and s.dmin1 > 0 and z(s, 4 * (n0 - 1) - pp).* < tol * (s.sigma + s.dn1) and @abs(s.dn) < tol * s.sigma) {
            // convergence hidden by a negative dn
            z(s, 4 * (n0 - 1) - pp + 2).* = 0;
            s.dmin = 0;
            break false;
        }
        if (s.dmin < 0) {
            // tau too big: select a new shift and try again
            s.nfail += 1;
            if (s.ttype < -22) {
                s.tau = 0;
            } else if (s.dmin1 > 0) {
                s.tau = (s.tau + s.dmin) * (1 - 2 * eps);
                s.ttype -= 11;
            } else {
                s.tau *= qurtr;
                s.ttype -= 12;
            }
            continue;
        }
        if (std.math.isNan(s.dmin)) {
            if (s.tau == 0) break true;
            s.tau = 0;
            continue;
        }
        break true;
    };
    if (unshifted) {
        // `80`: a risk of underflow
        dqds.lasq6(F, s);
        s.ndiv += n0 - lo + 2;
        s.iter += 1;
        s.tau = 0;
    }
    // `90`
    if (s.tau < s.sigma) {
        s.desig += s.tau;
        const t = s.sigma + s.desig;
        s.desig -= t - s.sigma;
        s.sigma = t;
    } else {
        const t = s.sigma + s.tau;
        s.desig = s.sigma - (t - s.tau) + s.desig;
        s.sigma = t;
    }
}
