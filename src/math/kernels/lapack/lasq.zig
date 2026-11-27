//! dqds for a bidiagonal's singular values alone, `xLASQ1` and `xLASQ2`,
//! Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dlasq1.f`, `dlasq2.f`,
//! OpenBLAS v0.3.34-317; modified BSD), over `lasq_steps.zig` and
//! `lasq_dqds.zig`. Indices are the reference's, 1-based and signed.

const std = @import("std");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");
const bdsqr = @import("bdsqr.zig");
const dqds = @import("lasq_dqds.zig");
const steps = @import("lasq_steps.zig");

const Qd = dqds.Qd;

/// How `lasq1` ended: the values found (`d` sorted, nonincreasing); dqds
/// stopped short (`INFO = 2`) with `d` and `e` the bidiagonal it reached,
/// for the QR iteration to finish; or nothing usable (`INFO` 1 or 3, a NaN
/// or a bad split), `d` and `e` as given but for `d`'s signs.
pub const Outcome = enum { converged, partial, failed };

/// `DLASQ1` (`dlasq1.f`): the singular values of the `n x n` bidiagonal
/// with diagonal `d` and off-diagonal `e` (`n - 1`) into `d`,
/// nonincreasing, by dqds (Fernando and Parlett; Parlett and Marques), to
/// high relative accuracy. `work` holds `4 n`.
pub fn lasq1(comptime F: type, d: []F, e: []F, work: []F) Outcome {
    const n = d.len;
    if (n == 0) return .converged;
    if (n == 1) {
        d[0] = @abs(d[0]);
        return .converged;
    }
    if (n == 2) {
        const v = bdsqr.las2(F, d[0], e[0], d[1]);
        d[0] = v.ssmax;
        d[1] = v.ssmin;
        return .converged;
    }
    var sigmx: F = 0;
    for (0..n - 1) |i| {
        d[i] = @abs(d[i]);
        sigmx = @max(sigmx, @abs(e[i]));
    }
    d[n - 1] = @abs(d[n - 1]);
    if (sigmx == 0) {
        std.mem.sort(F, d, {}, std.sort.desc(F));
        return .converged;
    }
    // Fortran's `MAX` keeping a NaN, as `DISNAN(SIGMX)` then tests
    for (d) |x| sigmx = if (std.math.isNan(x) or std.math.isNan(sigmx)) std.math.nan(F) else @max(sigmx, x);
    if (std.math.isNan(sigmx)) return .failed;
    // the squares of `d` and `e`, scaled to avoid underflow
    const eps = std.math.floatEps(F);
    const safmin = std.math.floatMin(F);
    const scale = @sqrt(eps / safmin);
    for (0..n) |i| work[2 * i] = d[i];
    for (0..n - 1) |i| work[2 * i + 1] = e[i];
    const w: ops.Mat(F) = .{ .ptr = work.ptr, .rows = 2 * n - 1, .cols = 1, .ld = 2 * n - 1 };
    auxiliary.lascl(F, false, sigmx, scale, w);
    for (work[0 .. 2 * n - 1]) |*x| x.* *= x.*;
    work[2 * n - 1] = 0;
    switch (lasq2(F, n, work)) {
        .converged => {
            for (0..n) |i| d[i] = @sqrt(work[i]);
            auxiliary.lascl(F, false, scale, sigmx, .{ .ptr = d.ptr, .rows = n, .cols = 1, .ld = n });
            return .converged;
        },
        .partial => {
            for (0..n) |i| {
                d[i] = @sqrt(work[2 * i]);
                if (i < n - 1) e[i] = @sqrt(work[2 * i + 1]);
            }
            auxiliary.lascl(F, false, scale, sigmx, .{ .ptr = d.ptr, .rows = n, .cols = 1, .ld = n });
            auxiliary.lascl(F, false, scale, sigmx, .{ .ptr = e.ptr, .rows = n - 1, .cols = 1, .ld = n - 1 });
            return .partial;
        },
        .failed => return .failed,
    }
}

/// `TOL` (`dlasq2.f`): `eps` times this is the convergence tolerance.
const hundrd = 100;

/// `DLASQ2` (`dlasq2.f`): the eigenvalues of the qd array `z` (`4 n`; `q`
/// and `e` interleaved in its first `2 n`) into `z`'s first `n`,
/// nonincreasing; `.partial` leaves the reached `q` and `e` interleaved
/// for `lasq1`. A negative entry, which the squares `lasq1` passes cannot
/// have, is `.failed` (the reference's `INFO = -(200 + k)`).
fn lasq2(comptime F: type, n: usize, zs: []F) Outcome {
    const eps = std.math.floatEps(F);
    const safmin = std.math.floatMin(F);
    const tol = eps * hundrd;
    const tol2 = tol * tol;
    var s: Qd(F) = .{ .z = zs };
    const z = struct {
        fn p(q: *Qd(F), i: isize) *F {
            return q.at(i);
        }
    }.p;
    const nn: isize = @intCast(n);
    // `n` is at least 3 here; `lasq1` takes 1 and 2 itself
    z(&s, 2 * nn).* = 0;
    var emin = z(&s, 2).*;
    var qmax: F = 0;
    var d: F = 0;
    var e: F = 0;
    var k: isize = 1;
    while (k <= 2 * (nn - 1)) : (k += 2) {
        if (z(&s, k).* < 0 or z(&s, k + 1).* < 0) return .failed;
        d += z(&s, k).*;
        e += z(&s, k + 1).*;
        qmax = @max(qmax, z(&s, k).*);
        emin = @min(emin, z(&s, k + 1).*);
    }
    if (z(&s, 2 * nn - 1).* < 0) return .failed;
    d += z(&s, 2 * nn - 1).*;
    qmax = @max(qmax, z(&s, 2 * nn - 1).*);
    // diagonal already
    if (e == 0) {
        k = 2;
        while (k <= nn) : (k += 1) z(&s, k).* = z(&s, 2 * k - 1).*;
        std.mem.sort(F, zs[0..n], {}, std.sort.desc(F));
        z(&s, 2 * nn - 1).* = d;
        return .converged;
    }
    const trace = d + e;
    if (trace == 0) {
        z(&s, 2 * nn - 1).* = 0;
        return .converged;
    }
    // rearrange for locality: `z = (q1, qq1, e1, ee1, q2, ...)`
    k = 2 * nn;
    while (k >= 2) : (k -= 2) {
        z(&s, 2 * k).* = 0;
        z(&s, 2 * k - 1).* = z(&s, k).*;
        z(&s, 2 * k - 2).* = 0;
        z(&s, 2 * k - 3).* = z(&s, k - 1).*;
    }
    s.lo = 1;
    s.n0 = nn;
    // reverse the qd array if warranted
    if (steps.cbias * z(&s, 4 * s.lo - 3).* < z(&s, 4 * s.n0 - 3).*) {
        const ipn4 = 4 * (s.lo + s.n0);
        var k4 = 4 * s.lo;
        while (k4 <= 2 * (s.lo + s.n0 - 1)) : (k4 += 4) {
            std.mem.swap(F, z(&s, k4 - 3), z(&s, ipn4 - k4 - 3));
            std.mem.swap(F, z(&s, k4 - 1), z(&s, ipn4 - k4 - 5));
        }
    }
    // initial split checking via dqd and Li's test
    var pp: isize = 0;
    for (0..2) |_| {
        d = z(&s, 4 * s.n0 + pp - 3).*;
        var k4 = 4 * (s.n0 - 1) + pp;
        while (k4 >= 4 * s.lo + pp) : (k4 -= 4) {
            if (z(&s, k4 - 1).* <= tol2 * d) {
                z(&s, k4 - 1).* = -0.0;
                d = z(&s, k4 - 3).*;
            } else d = z(&s, k4 - 3).* * (d / (d + z(&s, k4 - 1).*));
        }
        // dqd maps `z` to `zz` and back, Li's test
        emin = z(&s, 4 * s.lo + pp + 1).*;
        d = z(&s, 4 * s.lo + pp - 3).*;
        k4 = 4 * s.lo + pp;
        while (k4 <= 4 * (s.n0 - 1) + pp) : (k4 += 4) {
            z(&s, k4 - 2 * pp - 2).* = d + z(&s, k4 - 1).*;
            if (z(&s, k4 - 1).* <= tol2 * d) {
                z(&s, k4 - 1).* = -0.0;
                z(&s, k4 - 2 * pp - 2).* = d;
                z(&s, k4 - 2 * pp).* = 0;
                d = z(&s, k4 + 1).*;
            } else if (safmin * z(&s, k4 + 1).* < z(&s, k4 - 2 * pp - 2).* and safmin * z(&s, k4 - 2 * pp - 2).* < z(&s, k4 + 1).*) {
                const temp = z(&s, k4 + 1).* / z(&s, k4 - 2 * pp - 2).*;
                z(&s, k4 - 2 * pp).* = z(&s, k4 - 1).* * temp;
                d *= temp;
            } else {
                z(&s, k4 - 2 * pp).* = z(&s, k4 + 1).* * (z(&s, k4 - 1).* / z(&s, k4 - 2 * pp - 2).*);
                d = z(&s, k4 + 1).* * (d / z(&s, k4 - 2 * pp - 2).*);
            }
            emin = @min(emin, z(&s, k4 - 2 * pp).*);
        }
        z(&s, 4 * s.n0 - pp - 2).* = d;
        // the new qmax
        qmax = z(&s, 4 * s.lo - pp - 2).*;
        k4 = 4 * s.lo - pp + 2;
        while (k4 <= 4 * s.n0 - pp - 2) : (k4 += 4) qmax = @max(qmax, z(&s, k4).*);
        pp = 1 - pp;
    }
    s.iter = 2;
    s.ndiv = 2 * (s.n0 - s.lo);
    var iwhila: isize = 1;
    while (iwhila <= nn + 1) : (iwhila += 1) {
        if (s.n0 < 1) return finish(F, &s, n, trace);
        // split into unreduced blocks, the last one's start `i0`
        s.desig = 0;
        s.sigma = if (s.n0 == nn) 0 else -z(&s, 4 * s.n0 - 1).*;
        if (s.sigma < 0) return .failed;
        var emax: F = 0;
        emin = if (s.n0 > s.lo) @abs(z(&s, 4 * s.n0 - 5).*) else 0;
        var qmin = z(&s, 4 * s.n0 - 3).*;
        qmax = qmin;
        var k4 = 4 * s.n0;
        while (k4 >= 8) : (k4 -= 4) {
            if (z(&s, k4 - 5).* <= 0) break;
            if (qmin >= 4 * emax) {
                qmin = @min(qmin, z(&s, k4 - 3).*);
                emax = @max(emax, z(&s, k4 - 5).*);
            }
            qmax = @max(qmax, z(&s, k4 - 7).* + z(&s, k4 - 5).*);
            emin = @min(emin, z(&s, k4 - 5).*);
        } else k4 = 4;
        s.lo = @divTrunc(k4, 4);
        s.pp = 0;
        if (s.n0 - s.lo > 1) {
            var dee = z(&s, 4 * s.lo - 3).*;
            var deemin = dee;
            var kmin = s.lo;
            k4 = 4 * s.lo + 1;
            while (k4 <= 4 * s.n0 - 3) : (k4 += 4) {
                dee = z(&s, k4).* * (dee / (dee + z(&s, k4 - 2).*));
                if (dee <= deemin) {
                    deemin = dee;
                    kmin = @divTrunc(k4 + 3, 4);
                }
            }
            if ((kmin - s.lo) * 2 < s.n0 - kmin and deemin <= 0.5 * z(&s, 4 * s.n0 - 3).*) {
                const ipn4 = 4 * (s.lo + s.n0);
                s.pp = 2;
                k4 = 4 * s.lo;
                while (k4 <= 2 * (s.lo + s.n0 - 1)) : (k4 += 4) {
                    std.mem.swap(F, z(&s, k4 - 3), z(&s, ipn4 - k4 - 3));
                    std.mem.swap(F, z(&s, k4 - 2), z(&s, ipn4 - k4 - 2));
                    std.mem.swap(F, z(&s, k4 - 1), z(&s, ipn4 - k4 - 5));
                    std.mem.swap(F, z(&s, k4), z(&s, ipn4 - k4 - 4));
                }
            }
        }
        // the initial shift's lower bound
        s.dmin = -@max(0, qmin - 2 * @sqrt(qmin) * @sqrt(emax));
        s.qmax = qmax;
        const nbig = 100 * (s.n0 - s.lo + 1);
        var iwhilb: isize = 1;
        const done = while (iwhilb <= nbig) : (iwhilb += 1) {
            if (s.lo > s.n0) break true;
            steps.lasq3(F, &s);
            s.pp = 1 - s.pp;
            // when the end is near, check for splits
            if (s.pp == 0 and s.n0 - s.lo >= 3) {
                if (z(&s, 4 * s.n0).* <= tol2 * s.qmax or z(&s, 4 * s.n0 - 1).* <= tol2 * s.sigma) {
                    var splt = s.lo - 1;
                    s.qmax = z(&s, 4 * s.lo - 3).*;
                    emin = z(&s, 4 * s.lo - 1).*;
                    var oldemn = z(&s, 4 * s.lo).*;
                    k4 = 4 * s.lo;
                    while (k4 <= 4 * (s.n0 - 3)) : (k4 += 4) {
                        if (z(&s, k4).* <= tol2 * z(&s, k4 - 3).* or z(&s, k4 - 1).* <= tol2 * s.sigma) {
                            z(&s, k4 - 1).* = -s.sigma;
                            splt = @divTrunc(k4, 4);
                            s.qmax = 0;
                            emin = z(&s, k4 + 3).*;
                            oldemn = z(&s, k4 + 4).*;
                        } else {
                            s.qmax = @max(s.qmax, z(&s, k4 + 1).*);
                            emin = @min(emin, z(&s, k4 - 1).*);
                            oldemn = @min(oldemn, z(&s, k4).*);
                        }
                    }
                    z(&s, 4 * s.n0 - 1).* = emin;
                    z(&s, 4 * s.n0).* = oldemn;
                    s.lo = splt + 1;
                }
            }
        } else false;
        if (!done) {
            restore(F, &s, n);
            return .partial;
        }
    }
    return .failed;
}

/// `dlasq2.f`'s `INFO = 2` exit: too many iterations; restore the shift
/// `sigma` into each unfinished block's `d`s and `e`s and lay them out
/// interleaved. The reference adds the shifts to block `i0 ..= n0` each
/// time and never moves `I1` past a split, so with more than one block it
/// loops forever; this walks the blocks `i1 ..= n1` back to the first, as
/// its comments describe (issue #774).
fn restore(comptime F: type, s: *Qd(F), n: usize) void {
    const z = struct {
        fn p(q: *Qd(F), i: isize) *F {
            return q.at(i);
        }
    }.p;
    var first = s.lo;
    var n1 = s.n0;
    var sigma = s.sigma;
    while (true) {
        var tempq = z(s, 4 * first - 3).*;
        z(s, 4 * first - 3).* += sigma;
        var k = first + 1;
        while (k <= n1) : (k += 1) {
            const tempe = z(s, 4 * k - 5).*;
            z(s, 4 * k - 5).* *= tempq / z(s, 4 * k - 7).*;
            tempq = z(s, 4 * k - 3).*;
            z(s, 4 * k - 3).* += sigma + tempe - z(s, 4 * k - 5).*;
        }
        if (first <= 1) break;
        n1 = first - 1;
        first = n1;
        while (first >= 2 and z(s, 4 * first - 5).* >= 0) first -= 1;
        sigma = -z(s, 4 * n1 - 1).*;
    }
    const nn: isize = @intCast(n);
    var k: isize = 1;
    while (k <= nn) : (k += 1) {
        z(s, 2 * k - 1).* = z(s, 4 * k - 3).*;
        // only the block `1 ..= n0` is unfinished; the other `e`s are
        // essentially zero, though other data may sit there
        z(s, 2 * k).* = if (k < s.n0) z(s, 4 * k - 1).* else 0;
    }
}

/// `dlasq2.f`'s `170`: every eigenvalue found; move them to `z`'s front,
/// sorted, and the statistics after.
fn finish(comptime F: type, s: *Qd(F), n: usize, trace: F) Outcome {
    const nn: isize = @intCast(n);
    var k: isize = 2;
    while (k <= nn) : (k += 1) s.at(k).* = s.at(4 * k - 3).*;
    std.mem.sort(F, s.z[0..n], {}, std.sort.desc(F));
    var e: F = 0;
    k = nn;
    while (k >= 1) : (k -= 1) e += s.at(k).*;
    s.at(2 * nn + 1).* = trace;
    s.at(2 * nn + 2).* = e;
    s.at(2 * nn + 3).* = @floatFromInt(s.iter);
    s.at(2 * nn + 4).* = @as(F, @floatFromInt(s.ndiv)) / @as(F, @floatFromInt(nn * nn));
    s.at(2 * nn + 5).* = hundrd * @as(F, @floatFromInt(s.nfail)) / @as(F, @floatFromInt(s.iter));
    return .converged;
}
