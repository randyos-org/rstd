//! `stebr.zig`'s secular solve, Reference-LAPACK's `xLAED9_BR` and
//! `xLAED4_BR` (`lapack-netlib/SRC/dlaed9_br.f`, `dlaed4_br.f`,
//! Reference-LAPACK PR 1424, OpenBLAS v0.3.34-794-gdbf2c24c9; modified
//! BSD): the merged eigenvalues, and the parent's two boundary rows of their
//! eigenvectors streamed one secular vector at a time, never a `k x k`
//! block. Indices are 0-based.

const std = @import("std");
const laed4_solver = @import("laed4.zig");

const laed4 = laed4_solver.laed4;

/// The secular solve's scratch, `k` long each (`near` `2 k`).
pub fn Scratch(comptime F: type) type {
    return struct {
        delta: []F,
        wprod: []F,
        wsgn: []F,
        tau: []F,
        near: []F,
    };
}

/// `DLAED9_BR`: the `k` roots of the secular equation with poles
/// `dlambda`, weights `w` and `rho` into `d[0..k]`; with `qout` (the
/// reference's `NBR = 2`), the boundary rows `qin` (`2 x k`) times the
/// roots' unit eigenvectors into `qout`, and `w` replaced by the weights
/// Gu and Eisenstat's formula recovers from the roots (`xLAED3`'s).
pub fn laed9Br(comptime F: type, d: []F, qin: []const F, qout: ?[]F, rho: F, dlambda: []const F, w: []F, s: Scratch(F)) error{NoConvergence}!void {
    const k = d.len;
    if (k == 0) return;
    const dl = dlambda[0..k];
    const delta = s.delta[0..k];
    // the root merge: no parent reads the rows
    const q = qout orelse {
        for (0..k) |j| d[j] = try laed4(F, j, dl, w[0..k], delta, rho);
        return;
    };
    // `DLAED4` (`DLAED5`) returns the unit eigenvectors in `delta`, as
    // `DLAED9`'s copy to `S`
    if (k <= 2) {
        for (0..k) |j| {
            d[j] = try laed4(F, j, dl, w[0..k], delta, rho);
            var acc1: F = 0;
            var acc2: F = 0;
            for (0..k) |i| {
                acc1 += qin[2 * i] * delta[i];
                acc2 += qin[2 * i + 1] * delta[i];
            }
            q[2 * j] = acc1;
            q[2 * j + 1] = acc2;
        }
        return;
    }
    // every root, the product for the updated w accumulated from each
    // `delta` as it comes, while `DLAED4_BR`'s accurate values are at hand
    const wsgn = s.wsgn[0..k];
    const wprod = s.wprod[0..k];
    const tau = s.tau[0..k];
    const near = s.near[0 .. 2 * k];
    @memcpy(wsgn, w[0..k]);
    @memset(wprod, 1);
    for (0..k) |j| {
        const r = try laed4Br(F, j, dl, wsgn, delta, rho);
        d[j] = r.dlam;
        tau[j] = r.tau;
        near[2 * j] = r.near1;
        near[2 * j + 1] = r.near2;
        wprod[j] *= delta[j];
        for (0..j) |i| wprod[i] *= delta[i] / (dl[i] - dl[j]);
        for (j + 1..k) |i| wprod[i] *= delta[i] / (dl[i] - dl[j]);
    }
    for (w[0..k], wprod, wsgn) |*x, p, sg| x.* = std.math.copysign(@sqrt(-p), sg);
    // each secular vector rebuilt from its origin and tau, the cached
    // pole-adjacent values where cancellation is worst, through the rows
    for (0..k) |j| {
        // the origin, and the two poles beside the root
        const iorg = if (j == k - 1) k - 1 else if (tau[j] >= 0) j else j + 1;
        const n1 = if (j == k - 1) k - 2 else j;
        const n2 = n1 + 1;
        var temp: F = 0;
        var acc1: F = 0;
        var acc2: F = 0;
        for (0..k) |i| {
            const den = if (i == n1) near[2 * j] else if (i == n2) near[2 * j + 1] else (dl[i] - dl[iorg]) - tau[j];
            const dlt = w[i] / den;
            temp += dlt * dlt;
            acc1 += qin[2 * i] * dlt;
            acc2 += qin[2 * i + 1] * dlt;
        }
        temp = @sqrt(temp);
        q[2 * j] = acc1 / temp;
        q[2 * j + 1] = acc2 / temp;
    }
}

/// `DLAED4_BR`: `DLAED4`'s root, with its origin's offset `tau` (`λ -
/// d(org)`, `org` the nearer of the root's two poles, `n - 1` for the last)
/// and the two pole-adjacent `delta` values. For `n >= 3` (`laed9Br`'s
/// only use; the reference's `n` 1 and 2 branches are left out).
fn laed4Br(comptime F: type, i: usize, d: []const F, z: []const F, delta: []F, rho: F) error{NoConvergence}!struct { dlam: F, tau: F, near1: F, near2: F } {
    const n = d.len;
    const dlam = try laed4(F, i, d, z, delta, rho);
    const org = if (i == n - 1) n - 1 else if (dlam < 0.5 * (d[i] + d[i + 1])) i else i + 1;
    // `-delta(org)`, not `dlam - d(org)`: `DLAED4`'s own rounded
    // denominator, without a rounding through dlam
    const lo = if (i == n - 1) n - 2 else i;
    return .{ .dlam = dlam, .tau = -delta[org], .near1 = delta[lo], .near2 = delta[lo + 1] };
}
