//! `xLASY2`: the Sylvester equation of 1 x 1 and 2 x 2 real blocks, Reference-LAPACK's
//! (`lapack-netlib/SRC/dlasy2.f`, OpenBLAS v0.3.34-317; modified BSD).

const std = @import("std");
const mem = @import("../../mem.zig");
const laln2_solve = @import("laln2.zig");

const Small = laln2_solve.Small;

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
    if (xswpiv[ipiv]) mem.swap(F, &x2[0], &x2[1]);
    out.x[0] = x2[0];
    if (n1 == 1) out.x[2] = x2[1] else out.x[1] = x2[1];
    return out;
}
/// `lasy2` with `n1 = n2 = 2`: `op(TL) X + sgn X op(TR) = s V` for 2 x 2
/// blocks, as a 4 x 4 system by complete pivoting.
pub fn lasy2Square(comptime F: type, ltranl: bool, ltranr: bool, sgn: F, tl: [*]const F, ldtl: usize, tr: [*]const F, ldtr: usize, v: [4]F) Small(F) {
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
            mem.swap([4]F, &t[ipsv], &t[i]);
            mem.swap(F, &bt[i], &bt[ipsv]);
        }
        if (jpsv != i) for (0..4) |r| mem.swap(F, &t[r][jpsv], &t[r][i]);
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
        if (jpiv[i] != i) mem.swap(F, &tmp[i], &tmp[jpiv[i]]);
    }
    out.x = tmp;
    return out;
}
