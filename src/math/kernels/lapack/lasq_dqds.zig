//! dqds's state and transforms, `xLASQ5` (with a shift) and `xLASQ6` (none,
//! guarding underflow), Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dlasq5.f`, `dlasq6.f`, OpenBLAS v0.3.34-317;
//! modified BSD). Indices are the reference's, 1-based and signed.

const std = @import("std");

/// `dlasq2.f`'s working variables, shared with `xLASQ3` to `xLASQ6` (the
/// reference passes each by reference): `z` the qd array (`4 n`), the
/// block `i0 ..= n0`, ping-pong `pp`, the shifts and the minimal `d`s.
pub fn Qd(comptime F: type) type {
    return struct {
        z: []F,
        lo: isize = 1,
        n0: isize = 0,
        pp: isize = 0,
        dmin: F = 0,
        dmin1: F = 0,
        dmin2: F = 0,
        dn: F = 0,
        dn1: F = 0,
        dn2: F = 0,
        sigma: F = 0,
        desig: F = 0,
        qmax: F = 0,
        g: F = 0,
        tau: F = 0,
        ttype: isize = 0,
        nfail: isize = 0,
        iter: isize = 0,
        ndiv: isize = 0,

        /// `Z(i)`.
        pub fn at(s: *@This(), i: isize) *F {
            return &s.z[@intCast(i - 1)];
        }
    };
}

/// Fortran's `MIN` as the IEEE path needs it: a NaN propagates, so
/// `xLASQ3`'s `DISNAN(DMIN)` sees it (Zig's `@min` would drop it).
pub fn fmin(comptime F: type, a: F, b: F) F {
    if (std.math.isNan(a) or std.math.isNan(b)) return std.math.nan(F);
    return @min(a, b);
}

/// `DLASQ5` (`dlasq5.f`), its IEEE branches (`IEEE` is
/// `ILAENV(10, ...) = 1` on IEEE hardware, which every target is): one
/// dqds transform with shift `tau` (0 when below `eps (sigma + tau) / 2`,
/// then each small `d` flushed to 0), from ping-pong `pp`.
pub fn lasq5(comptime F: type, s: *Qd(F), eps: F) void {
    const lo = s.lo;
    const n0 = s.n0;
    const pp = s.pp;
    if (n0 - lo - 1 <= 0) return;
    const dthresh = eps * (s.sigma + s.tau);
    if (s.tau < dthresh * 0.5) s.tau = 0;
    const tau = s.tau;
    const flush = tau == 0;
    var j4 = 4 * lo + pp - 3;
    var emin = s.at(j4 + 4).*;
    var d = s.at(j4).* - tau;
    s.dmin = d;
    s.dmin1 = -s.at(j4).*;
    j4 = 4 * lo;
    while (j4 <= 4 * (n0 - 3)) : (j4 += 4) {
        if (pp == 0) {
            s.at(j4 - 2).* = d + s.at(j4 - 1).*;
            const temp = s.at(j4 + 1).* / s.at(j4 - 2).*;
            d = d * temp - tau;
            if (flush and d < dthresh) d = 0;
            s.dmin = fmin(F, s.dmin, d);
            s.at(j4).* = s.at(j4 - 1).* * temp;
            emin = fmin(F, s.at(j4).*, emin);
        } else {
            s.at(j4 - 3).* = d + s.at(j4).*;
            const temp = s.at(j4 + 2).* / s.at(j4 - 3).*;
            d = d * temp - tau;
            if (flush and d < dthresh) d = 0;
            s.dmin = fmin(F, s.dmin, d);
            s.at(j4 - 1).* = s.at(j4).* * temp;
            emin = fmin(F, s.at(j4 - 1).*, emin);
        }
    }
    // the last two steps unrolled
    s.dn2 = d;
    s.dmin2 = s.dmin;
    j4 = 4 * (n0 - 2) - pp;
    var j4p2 = j4 + 2 * pp - 1;
    s.at(j4 - 2).* = s.dn2 + s.at(j4p2).*;
    s.at(j4).* = s.at(j4p2 + 2).* * (s.at(j4p2).* / s.at(j4 - 2).*);
    s.dn1 = s.at(j4p2 + 2).* * (s.dn2 / s.at(j4 - 2).*) - tau;
    s.dmin = fmin(F, s.dmin, s.dn1);
    s.dmin1 = s.dmin;
    j4 += 4;
    j4p2 = j4 + 2 * pp - 1;
    s.at(j4 - 2).* = s.dn1 + s.at(j4p2).*;
    s.at(j4).* = s.at(j4p2 + 2).* * (s.at(j4p2).* / s.at(j4 - 2).*);
    s.dn = s.at(j4p2 + 2).* * (s.dn1 / s.at(j4 - 2).*) - tau;
    s.dmin = fmin(F, s.dmin, s.dn);
    s.at(j4 + 2).* = s.dn;
    s.at(4 * n0 - pp).* = emin;
}

/// One `DLASQ6` step: `q` the new `q`, `zn` and `zd` the pair it divides,
/// `d` the running `d`; a zero `q` restarts `d` and `emin`.
fn step6(comptime F: type, s: *Qd(F), qi: isize, ei: isize, next_q: isize, e_old: isize, d: *F, emin: *F, safmin: F) void {
    s.at(qi).* = d.* + s.at(e_old).*;
    if (s.at(qi).* == 0) {
        s.at(ei).* = 0;
        d.* = s.at(next_q).*;
        s.dmin = d.*;
        emin.* = 0;
    } else if (safmin * s.at(next_q).* < s.at(qi).* and safmin * s.at(qi).* < s.at(next_q).*) {
        const temp = s.at(next_q).* / s.at(qi).*;
        s.at(ei).* = s.at(e_old).* * temp;
        d.* *= temp;
    } else {
        s.at(ei).* = s.at(next_q).* * (s.at(e_old).* / s.at(qi).*);
        d.* = s.at(next_q).* * (d.* / s.at(qi).*);
    }
}

/// `DLASQ6` (`dlasq6.f`): one dqd transform (no shift) from ping-pong `pp`,
/// guarding against underflow.
pub fn lasq6(comptime F: type, s: *Qd(F)) void {
    const lo = s.lo;
    const n0 = s.n0;
    const pp = s.pp;
    if (n0 - lo - 1 <= 0) return;
    const safmin = std.math.floatMin(F);
    var j4 = 4 * lo + pp - 3;
    var emin = s.at(j4 + 4).*;
    var d = s.at(j4).*;
    s.dmin = d;
    j4 = 4 * lo;
    while (j4 <= 4 * (n0 - 3)) : (j4 += 4) {
        if (pp == 0) {
            step6(F, s, j4 - 2, j4, j4 + 1, j4 - 1, &d, &emin, safmin);
            s.dmin = fmin(F, s.dmin, d);
            emin = fmin(F, emin, s.at(j4).*);
        } else {
            step6(F, s, j4 - 3, j4 - 1, j4 + 2, j4, &d, &emin, safmin);
            s.dmin = fmin(F, s.dmin, d);
            emin = fmin(F, emin, s.at(j4 - 1).*);
        }
    }
    // the last two steps unrolled
    s.dn2 = d;
    s.dmin2 = s.dmin;
    j4 = 4 * (n0 - 2) - pp;
    var j4p2 = j4 + 2 * pp - 1;
    var dn1 = s.dn2;
    step6(F, s, j4 - 2, j4, j4p2 + 2, j4p2, &dn1, &emin, safmin);
    s.dn1 = dn1;
    s.dmin = fmin(F, s.dmin, s.dn1);
    s.dmin1 = s.dmin;
    j4 += 4;
    j4p2 = j4 + 2 * pp - 1;
    var dn = s.dn1;
    step6(F, s, j4 - 2, j4, j4p2 + 2, j4p2, &dn, &emin, safmin);
    s.dn = dn;
    s.dmin = fmin(F, s.dmin, s.dn);
    s.at(j4 + 2).* = s.dn;
    s.at(4 * n0 - pp).* = emin;
}
