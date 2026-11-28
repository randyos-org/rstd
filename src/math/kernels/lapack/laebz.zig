//! `xLAEBZ`: bisection on Sturm counts of a symmetric tridiagonal matrix,
//! Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dlaebz.f`, OpenBLAS
//! v0.3.34-794-gdbf2c24c9; modified BSD), for `stebz.zig`. Only the reference's scalar
//! loop: `ILAENV(1, 'DSTEBZ')` is 1, so `dstebz.f` passes `NBMIN = 0` and
//! the vectorized branch never runs. Indices are 0-based.

const std = @import("std");

/// `IJOB`: count the eigenvalues in each interval (1), bisect each into
/// intervals holding one eigenvalue cluster (2), or close in on the points
/// where the count reaches `nval` (3).
pub const Job = enum { count, bisect, find };

/// `dlaebz.f`'s arrays: interval `j` is `(ab[j][0], ab[j][1]]`, `nab[j]`
/// the Sturm counts at its ends, `c[j]` its next point; `find` targets
/// `nval[j]`. At most `ab.len` (`MMAX`) intervals.
pub fn Intervals(comptime F: type) type {
    return struct { ab: [][2]F, nab: [][2]isize, c: []F, nval: ?[]isize = null };
}

/// `count`'s eigenvalue total, `MOUT`; or for `bisect` and `find` the
/// intervals now held (`MOUT`) and how many did not converge in `nitmax`
/// steps (`INFO`), the last of them.
pub const Outcome = struct { mout: usize, unconverged: usize = 0 };

/// `DLAEBZ` on the tridiagonal of diagonal `d` and squared off-diagonal
/// `e2`, over the first `minp` intervals of `iv`. Panics if bisection would
/// hold more than `iv.ab.len` intervals, which cannot happen when the
/// intervals hold at most that many eigenvalues, as `stebz`'s do.
pub fn laebz(comptime F: type, job: Job, nitmax: usize, minp: usize, abstol: F, reltol: F, pivmin: F, d: []const F, e2: []const F, iv: Intervals(F)) Outcome {
    if (job == .count) {
        var mout: isize = 0;
        for (0..minp) |ji| {
            for (0..2) |jp| iv.nab[ji][jp] = @intCast(countAt(F, d, e2, pivmin, iv.ab[ji][jp]));
            mout += iv.nab[ji][1] - iv.nab[ji][0];
        }
        return .{ .mout = @intCast(mout) };
    }
    var kf: usize = 0;
    var kl: usize = minp;
    if (job == .bisect) for (0..minp) |ji| {
        iv.c[ji] = 0.5 * (iv.ab[ji][0] + iv.ab[ji][1]);
    };
    for (0..nitmax) |_| {
        var klnew = kl;
        for (kf..kl) |ji| {
            const x = iv.c[ji];
            var k: isize = @intCast(sturm(F, d, e2, pivmin, x));
            if (job == .bisect) {
                k = @min(iv.nab[ji][1], @max(iv.nab[ji][0], k));
                if (k == iv.nab[ji][1]) {
                    iv.ab[ji][1] = x;
                } else if (k == iv.nab[ji][0]) {
                    iv.ab[ji][0] = x;
                } else {
                    // `INFO = MMAX + 1`, out of room
                    std.debug.assert(klnew < iv.ab.len);
                    iv.ab[klnew] = .{ x, iv.ab[ji][1] };
                    iv.nab[klnew] = .{ k, iv.nab[ji][1] };
                    klnew += 1;
                    iv.ab[ji][1] = x;
                    iv.nab[ji][1] = k;
                }
            } else {
                if (k <= iv.nval.?[ji]) {
                    iv.ab[ji][0] = x;
                    iv.nab[ji][0] = k;
                }
                if (k >= iv.nval.?[ji]) {
                    iv.ab[ji][1] = x;
                    iv.nab[ji][1] = k;
                }
            }
        }
        kl = klnew;
        // the converged intervals to the front
        var kfnew = kf;
        for (kf..kl) |ji| {
            const width = @abs(iv.ab[ji][1] - iv.ab[ji][0]);
            const size = @max(@abs(iv.ab[ji][1]), @abs(iv.ab[ji][0]));
            if (width < @max(abstol, pivmin, reltol * size) or iv.nab[ji][0] >= iv.nab[ji][1]) {
                if (ji > kfnew) {
                    std.mem.swap([2]F, &iv.ab[ji], &iv.ab[kfnew]);
                    std.mem.swap([2]isize, &iv.nab[ji], &iv.nab[kfnew]);
                    if (job == .find) std.mem.swap(isize, &iv.nval.?[ji], &iv.nval.?[kfnew]);
                }
                kfnew += 1;
            }
        }
        kf = kfnew;
        for (kf..kl) |ji| iv.c[ji] = 0.5 * (iv.ab[ji][0] + iv.ab[ji][1]);
        if (kf >= kl) break;
    }
    return .{ .mout = kl, .unconverged = kl - kf };
}

/// The eigenvalues at or below `x`, `count`'s form: a pivot under `pivmin`
/// in magnitude becomes `-pivmin`.
fn countAt(comptime F: type, d: []const F, e2: []const F, pivmin: F, x: F) usize {
    var t = d[0] - x;
    if (@abs(t) < pivmin) t = -pivmin;
    var k: usize = @intFromBool(t <= 0);
    for (1..d.len) |j| {
        t = d[j] - e2[j - 1] / t - x;
        if (@abs(t) < pivmin) t = -pivmin;
        k += @intFromBool(t <= 0);
    }
    return k;
}

/// The same, the bisection loop's form: a pivot at or below `pivmin`
/// counts and becomes at most `-pivmin`.
fn sturm(comptime F: type, d: []const F, e2: []const F, pivmin: F, x: F) usize {
    var t = d[0] - x;
    var k: usize = 0;
    if (t <= pivmin) {
        k = 1;
        t = @min(t, -pivmin);
    }
    for (1..d.len) |j| {
        t = d[j] - e2[j - 1] / t - x;
        if (t <= pivmin) {
            k += 1;
            t = @min(t, -pivmin);
        }
    }
    return k;
}
