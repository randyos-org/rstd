//! `stebr.zig`'s deflation, Reference-LAPACK's `xLAED8_BR`
//! (`lapack-netlib/SRC/dlaed8_br.f`, Reference-LAPACK PR 1424, OpenBLAS
//! v0.3.34-794-gdbf2c24c9; modified BSD): `xLAED8` carrying two rows of the
//! eigenvectors (`stebr.zig`'s boundary rows) instead of the whole blocks.
//! Indices are 0-based.

const std = @import("std");
const auxiliary = @import("auxiliary.zig");

/// The merge's scratch, `n` long each (`q2` `2 n`).
pub fn Scratch(comptime F: type) type {
    return struct {
        dlambda: []F,
        w: []F,
        q2: []F,
        indxp: []usize,
        indx: []usize,
        perm: []usize,
    };
}

/// `DLAED8_BR`: sorts the two halves' eigenvalues `d` (`cut` and `n -
/// cut`, each in its `indxq` order) into one increasing order and deflates
/// the rank-one update `rho z zᵀ`: small `z` entries, and close pairs by a
/// Givens rotation applied to the boundary rows `q` (`2 x n`, columns by
/// `d`'s physical index; null at the root, the reference's `NBR = 0`).
/// Returns `k`, the undeflated count: `w.dlambda[0..k]` and `w.w[0..k]` the
/// secular equation's poles and weights, `q`'s first `k` columns (in
/// `w.q2`) their boundary rows; `d[k..]` and `q[k..]` the deflated, in
/// decreasing order. `rho` becomes `|2 rho|`, `z` is overwritten.
///
/// The reference also records the rotations (`GIVCOL`, `GIVNUM`), which
/// nothing reads in this variant; they are left out.
pub fn laed8Br(comptime F: type, d: []F, q: ?[]F, indxq: []usize, rho: *F, cut: usize, z: []F, s: Scratch(F)) usize {
    const n = d.len;
    if (n == 0) return 0;
    const n1 = cut;
    const n2 = n - n1;
    if (rho.* < 0) for (z[n1..n]) |*x| {
        x.* *= -1;
    };
    const t = 1 / @sqrt(@as(F, 2));
    for (s.indx[0..n], 0..) |*x, j| x.* = j;
    for (z[0..n]) |*x| x.* *= t;
    rho.* = @abs(2 * rho.*);
    // the eigenvalues into increasing order
    for (indxq[cut..n]) |*x| x.* += cut;
    for (0..n) |i| {
        s.dlambda[i] = d[indxq[i]];
        s.w[i] = z[indxq[i]];
    }
    auxiliary.lamrg(F, n1, n2, s.dlambda, true, true, s.indx);
    for (0..n) |i| {
        d[i] = s.dlambda[s.indx[i]];
        z[i] = s.w[s.indx[i]];
    }
    const imax = iamax(F, z[0..n]);
    const jmax = iamax(F, d);
    const eps = std.math.floatEps(F) / 2; // `DLAMCH('Epsilon')`
    const tol = 8 * eps * @abs(d[jmax]);
    // everything deflates: the halves' order is the merged one
    if (rho.* * @abs(z[imax]) <= tol) {
        for (0..n) |j| {
            s.perm[j] = indxq[s.indx[j]];
            if (q) |qq| copyCol(F, s.q2, j, qq, s.perm[j]);
        }
        if (q) |qq| @memcpy(qq[0 .. 2 * n], s.q2[0 .. 2 * n]);
        return 0;
    }
    // small z(j) deflate to the back (k2 down from n); a pair of close
    // eigenvalues deflates the one a rotation zeroes
    var k: usize = 0;
    var k2: usize = n;
    var jlam: usize = 0;
    var j: usize = 0;
    while (true) : (j += 1) {
        if (rho.* * @abs(z[j]) <= tol) {
            k2 -= 1;
            s.indxp[k2] = j;
            // unreachable past z(imax), kept as the reference
            if (j == n - 1) return finish(F, d, q, indxq, s, k);
        } else {
            jlam = j;
            break;
        }
    }
    while (true) {
        j += 1;
        if (j >= n) break;
        if (rho.* * @abs(z[j]) <= tol) {
            k2 -= 1;
            s.indxp[k2] = j;
            continue;
        }
        var sn = z[jlam];
        var c = z[j];
        const tau = auxiliary.lapy2(F, c, sn);
        const td = d[j] - d[jlam];
        c /= tau;
        sn = -sn / tau;
        if (@abs(td * c * sn) <= tol) {
            // deflate by the rotation
            z[j] = tau;
            z[jlam] = 0;
            if (q) |qq| rot(F, qq, indxq[s.indx[jlam]], indxq[s.indx[j]], c, sn);
            const tt = d[jlam] * c * c + d[j] * sn * sn;
            d[j] = d[jlam] * sn * sn + d[j] * c * c;
            d[jlam] = tt;
            k2 -= 1;
            var i: usize = 1;
            while (true) {
                if (k2 + i < n and d[jlam] < d[s.indxp[k2 + i]]) {
                    s.indxp[k2 + i - 1] = s.indxp[k2 + i];
                    s.indxp[k2 + i] = jlam;
                    i += 1;
                } else {
                    s.indxp[k2 + i - 1] = jlam;
                    break;
                }
            }
            jlam = j;
        } else {
            s.w[k] = z[jlam];
            s.dlambda[k] = d[jlam];
            s.indxp[k] = jlam;
            k += 1;
            jlam = j;
        }
    }
    s.w[k] = z[jlam];
    s.dlambda[k] = d[jlam];
    s.indxp[k] = jlam;
    k += 1;
    return finish(F, d, q, indxq, s, k);
}

/// `dlaed8_br.f` from label 170: `dlambda` and the boundary rows into the
/// deflation's order, the deflated back into `d` and `q`.
fn finish(comptime F: type, d: []F, q: ?[]F, indxq: []const usize, s: Scratch(F), k: usize) usize {
    const n = d.len;
    for (0..n) |j| {
        const jp = s.indxp[j];
        s.dlambda[j] = d[jp];
        s.perm[j] = indxq[s.indx[jp]];
        if (q) |qq| copyCol(F, s.q2, j, qq, s.perm[j]);
    }
    if (k < n) {
        @memcpy(d[k..n], s.dlambda[k..n]);
        if (q) |qq| @memcpy(qq[2 * k .. 2 * n], s.q2[2 * k .. 2 * n]);
    }
    return k;
}

fn copyCol(comptime F: type, dst: []F, j: usize, src: []const F, p: usize) void {
    dst[2 * j] = src[2 * p];
    dst[2 * j + 1] = src[2 * p + 1];
}

/// `DROT` on columns `a` and `b` of the `2 x n` rows (`drot.f`).
fn rot(comptime F: type, q: []F, a: usize, b: usize, c: F, s: F) void {
    for (0..2) |r| {
        const x = q[2 * a + r];
        const y = q[2 * b + r];
        q[2 * a + r] = c * x + s * y;
        q[2 * b + r] = c * y - s * x;
    }
}

/// `IDAMAX`: the first largest `|x(i)|`.
fn iamax(comptime F: type, x: []const F) usize {
    var best: usize = 0;
    var m: F = -1;
    for (x, 0..) |v, i| if (@abs(v) > m) {
        m = @abs(v);
        best = i;
    };
    return best;
}
