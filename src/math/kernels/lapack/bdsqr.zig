//! `xBDSQR`: the SVD of a real bidiagonal matrix by implicit zero-shift QR, Reference-LAPACK's
//! (`lapack-netlib/SRC/dbdsqr.f`, `dlas2.f`, `dlasv2.f`, OpenBLAS v0.3.34-317; modified BSD), as
//! `bidiag_svd.zig`'s `qrAlgorithm` (T126).

const std = @import("std");
const mem = @import("../../../mem.zig");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");

const lartg = auxiliary.lartg;
const Mat = ops.Mat;

/// `xBDSQR` (`dbdsqr.f`; Reference-LAPACK 3.12.0 as OpenBLAS v0.3.34-317
/// ships it in `lapack-netlib/SRC/`, BSD-3), in place of faer's
/// `qr_algorithm` (T126): Demmel and Kahan's QR iteration, with zero-shift
/// sweeps where a shift would cost relative accuracy and their relative
/// convergence tests, so a graded bidiagonal keeps its small singular
/// values (faer's absolute `|d| <= eps` test lost them, M1's E4). `diag` and
/// `subdiag` (`n` each; `subdiag[n-1]` ignored) hold faer's lower bidiagonal
/// `B`; its transpose, upper, is the reference's, so the reference's `U` is
/// `v` here and its `VT`'s rows are `u`'s columns. The rotations go on the
/// vectors as they are made, in the reference's order, where it stores them
/// for `DLASR`. The values come out nonnegative and nonincreasing.
/// `ops.test_iter_limit`, when set, caps the sweeps.
pub fn qrAlgorithm(comptime F: type, diag: []F, subdiag: []F, u: ?Mat(F), v: ?Mat(F)) error{NoConvergence}!void {
    const n = diag.len;
    if (n == 0) return;
    // 1-based, as the reference: `d(i)` is `D(I)`, `e(i)` is `E(I)`
    const B = struct {
        dg: []F,
        sd: []F,
        fn d(b: @This(), i: usize) *F {
            return &b.dg[i - 1];
        }
        fn e(b: @This(), i: usize) *F {
            return &b.sd[i - 1];
        }
    };
    const b: B = .{ .dg = diag, .sd = subdiag };
    // faer's scaling to the largest entry, as `xGESVD` scales before
    // `xBDSQR`, which assumes it (its `MAXITR n² UNFL` floor)
    const max = @max(ops.absMax(F, diag), ops.absMax(F, subdiag[0 .. n - 1]));
    if (max == 0) return;
    ops.scaleDown(F, diag, max);
    ops.scaleDown(F, subdiag[0 .. n - 1], max);
    if (n > 1) bdsqrSweeps(F, B, b, n, u, v) catch |err| return unscaled(F, diag, subdiag[0 .. n - 1], max, err);
    // `dbdsqr.f:707-748`: nonnegative, then decreasing by selection
    for (1..n + 1) |i| {
        if (b.d(i).* == 0) b.d(i).* = 0;
        if (b.d(i).* < 0) {
            b.d(i).* = -b.d(i).*;
            if (u) |um| for (um.col(i - 1)[0..um.rows]) |*x| {
                x.* = -x.*;
            };
        }
    }
    for (1..n) |i| {
        var isub: usize = 1;
        var smin = b.d(1).*;
        for (2..n + 2 - i) |j| if (b.d(j).* <= smin) {
            isub = j;
            smin = b.d(j).*;
        };
        const last = n + 1 - i;
        if (isub != last) {
            b.d(isub).* = b.d(last).*;
            b.d(last).* = smin;
            if (u) |um| ops.swapCols(F, um, isub - 1, last - 1);
            if (v) |vm| ops.swapCols(F, vm, isub - 1, last - 1);
        }
    }
    for (diag) |*x| x.* *= max;
    for (subdiag[0 .. n - 1]) |*x| x.* *= max;
}

/// `dbdsqr.f:327-705`: the QR iteration on `n > 1`.
fn bdsqrSweeps(comptime F: type, comptime B: type, b: B, n: usize, u: ?Mat(F), v: ?Mat(F)) error{NoConvergence}!void {
    // `DLAMCH('E')`, the unit roundoff; `DLAMCH('S')`
    const eps = std.math.floatEps(F) / 2;
    const unfl = std.math.floatMin(F);
    const maxitr = 6;
    const nf: F = @floatFromInt(n);
    // `TOLMUL`, in `f64`: `std.math.pow` has no `f128` (every `eps` fits)
    const tolmul: F = @floatCast(@max(10, @min(100, std.math.pow(f64, @floatCast(eps), -0.125))));
    const tol = tolmul * eps;
    var smax: F = 0;
    for (1..n + 1) |i| smax = @max(smax, @abs(b.d(i).*));
    for (1..n) |i| smax = @max(smax, @abs(b.e(i).*));
    // `TOL >= 0`: relative accuracy, and its threshold from a lower bound
    // on the smallest singular value
    var sminoa = @abs(b.d(1).*);
    if (sminoa != 0) {
        var mu = sminoa;
        for (2..n + 1) |i| {
            mu = @abs(b.d(i).*) * (mu / (mu + @abs(b.e(i - 1).*)));
            sminoa = @min(sminoa, mu);
            if (sminoa == 0) break;
        }
    }
    sminoa /= @sqrt(nf);
    const thresh = @max(tol * sminoa, maxitr * (nf * (nf * unfl)));
    const maxitdivn = maxitr * n;
    const sweep_cap = ops.iterLimit(std.math.maxInt(usize));
    var sweeps: usize = 0;
    var iterdivn: usize = 0;
    var iter: usize = 0;
    var oldll: usize = 0;
    var oldm: usize = 0;
    var idir: u2 = 0;
    var m = n;
    // label 60
    outer: while (m > 1) {
        if (iter >= n) {
            iter -= n;
            iterdivn += 1;
            if (iterdivn >= maxitdivn) return error.NoConvergence;
        }
        // the bottom block, `ll..m`, unreduced
        smax = @abs(b.d(m).*);
        var ll: usize = 0;
        for (1..m) |lll| {
            const l = m - lll;
            const abss = @abs(b.d(l).*);
            const abse = @abs(b.e(l).*);
            if (abse <= thresh) {
                ll = l;
                break;
            }
            smax = @max(smax, @max(abss, abse));
        }
        if (ll > 0) {
            b.e(ll).* = 0;
            if (ll == m - 1) {
                m -= 1;
                continue :outer;
            }
        }
        ll += 1;
        if (ll == m - 1) {
            const sv = lasv2(F, b.d(m - 1).*, b.e(m - 1).*, b.d(m).*);
            b.d(m - 1).* = sv.ssmax;
            b.e(m - 1).* = 0;
            b.d(m).* = sv.ssmin;
            if (u) |um| rotCols(F, um, m - 2, sv.csr, sv.snr);
            if (v) |vm| rotCols(F, vm, m - 2, sv.csl, sv.snl);
            m -= 2;
            continue :outer;
        }
        // a new block: chase toward its smaller end
        if (ll > oldm or m < oldll) idir = if (@abs(b.d(ll).*) >= @abs(b.d(m).*)) 1 else 2;
        var smin: F = 0;
        if (idir == 1) {
            if (@abs(b.e(m - 1).*) <= tol * @abs(b.d(m).*)) {
                b.e(m - 1).* = 0;
                continue :outer;
            }
            var mu = @abs(b.d(ll).*);
            smin = mu;
            for (ll..m) |l| {
                if (@abs(b.e(l).*) <= tol * mu) {
                    b.e(l).* = 0;
                    continue :outer;
                }
                mu = @abs(b.d(l + 1).*) * (mu / (mu + @abs(b.e(l).*)));
                smin = @min(smin, mu);
            }
        } else {
            if (@abs(b.e(ll).*) <= tol * @abs(b.d(ll).*)) {
                b.e(ll).* = 0;
                continue :outer;
            }
            var mu = @abs(b.d(m).*);
            smin = mu;
            var l = m - 1;
            while (true) : (l -= 1) {
                if (@abs(b.e(l).*) <= tol * mu) {
                    b.e(l).* = 0;
                    continue :outer;
                }
                mu = @abs(b.d(l).*) * (mu / (mu + @abs(b.e(l).*)));
                smin = @min(smin, mu);
                if (l == ll) break;
            }
        }
        oldll = ll;
        oldm = m;
        // zero shift where a shift would cost relative accuracy
        var shift: F = 0;
        if (!(nf * tol * (smin / smax) <= @max(eps, 0.01 * tol))) {
            const sll = if (idir == 1) @abs(b.d(ll).*) else @abs(b.d(m).*);
            shift = if (idir == 1) las2(F, b.d(m - 1).*, b.e(m - 1).*, b.d(m).*) else las2(F, b.d(ll).*, b.e(ll).*, b.d(ll + 1).*);
            if (sll > 0 and (shift / sll) * (shift / sll) < eps) shift = 0;
        }
        iter += m - ll;
        if (sweeps >= sweep_cap) return error.NoConvergence;
        sweeps += 1;
        if (shift == 0) {
            if (idir == 1) {
                var cs: F = 1;
                var oldcs: F = 1;
                var oldsn: F = 0;
                for (ll..m) |i| {
                    const r1 = lartg(F, b.d(i).* * cs, b.e(i).*);
                    cs = r1.c;
                    if (i > ll) b.e(i - 1).* = oldsn * r1.r;
                    const r2 = lartg(F, oldcs * r1.r, b.d(i + 1).* * r1.s);
                    oldcs = r2.c;
                    oldsn = r2.s;
                    b.d(i).* = r2.r;
                    if (u) |um| rotCols(F, um, i - 1, cs, r1.s);
                    if (v) |vm| rotCols(F, vm, i - 1, oldcs, oldsn);
                }
                const h = b.d(m).* * cs;
                b.d(m).* = h * oldcs;
                b.e(m - 1).* = h * oldsn;
                if (@abs(b.e(m - 1).*) <= thresh) b.e(m - 1).* = 0;
            } else {
                var cs: F = 1;
                var oldcs: F = 1;
                var oldsn: F = 0;
                var i = m;
                while (i > ll) : (i -= 1) {
                    const r1 = lartg(F, b.d(i).* * cs, b.e(i - 1).*);
                    cs = r1.c;
                    if (i < m) b.e(i).* = oldsn * r1.r;
                    const r2 = lartg(F, oldcs * r1.r, b.d(i - 1).* * r1.s);
                    oldcs = r2.c;
                    oldsn = r2.s;
                    b.d(i).* = r2.r;
                    if (u) |um| rotCols(F, um, i - 2, oldcs, -oldsn);
                    if (v) |vm| rotCols(F, vm, i - 2, cs, -r1.s);
                }
                const h = b.d(ll).* * cs;
                b.d(ll).* = h * oldcs;
                b.e(ll).* = h * oldsn;
                if (@abs(b.e(ll).*) <= thresh) b.e(ll).* = 0;
            }
        } else if (idir == 1) {
            var f = (@abs(b.d(ll).*) - shift) * (sign(F, 1, b.d(ll).*) + shift / b.d(ll).*);
            var g = b.e(ll).*;
            for (ll..m) |i| {
                const r1 = lartg(F, f, g);
                if (i > ll) b.e(i - 1).* = r1.r;
                f = r1.c * b.d(i).* + r1.s * b.e(i).*;
                b.e(i).* = r1.c * b.e(i).* - r1.s * b.d(i).*;
                g = r1.s * b.d(i + 1).*;
                b.d(i + 1).* = r1.c * b.d(i + 1).*;
                const r2 = lartg(F, f, g);
                b.d(i).* = r2.r;
                f = r2.c * b.e(i).* + r2.s * b.d(i + 1).*;
                b.d(i + 1).* = r2.c * b.d(i + 1).* - r2.s * b.e(i).*;
                if (i < m - 1) {
                    g = r2.s * b.e(i + 1).*;
                    b.e(i + 1).* = r2.c * b.e(i + 1).*;
                }
                if (u) |um| rotCols(F, um, i - 1, r1.c, r1.s);
                if (v) |vm| rotCols(F, vm, i - 1, r2.c, r2.s);
            }
            b.e(m - 1).* = f;
            if (@abs(b.e(m - 1).*) <= thresh) b.e(m - 1).* = 0;
        } else {
            var f = (@abs(b.d(m).*) - shift) * (sign(F, 1, b.d(m).*) + shift / b.d(m).*);
            var g = b.e(m - 1).*;
            var i = m;
            while (i > ll) : (i -= 1) {
                const r1 = lartg(F, f, g);
                if (i < m) b.e(i).* = r1.r;
                f = r1.c * b.d(i).* + r1.s * b.e(i - 1).*;
                b.e(i - 1).* = r1.c * b.e(i - 1).* - r1.s * b.d(i).*;
                g = r1.s * b.d(i - 1).*;
                b.d(i - 1).* = r1.c * b.d(i - 1).*;
                const r2 = lartg(F, f, g);
                b.d(i).* = r2.r;
                f = r2.c * b.e(i - 1).* + r2.s * b.d(i - 1).*;
                b.d(i - 1).* = r2.c * b.d(i - 1).* - r2.s * b.e(i - 1).*;
                if (i > ll + 1) {
                    g = r2.s * b.e(i - 2).*;
                    b.e(i - 2).* = r2.c * b.e(i - 2).*;
                }
                if (u) |um| rotCols(F, um, i - 2, r2.c, -r2.s);
                if (v) |vm| rotCols(F, vm, i - 2, r1.c, -r1.s);
            }
            b.e(ll).* = f;
            if (@abs(b.e(ll).*) <= thresh) b.e(ll).* = 0;
        }
    }
}

/// One `DLASR`/`DROT` plane on columns `j`, `j + 1` (0-based): the
/// reference's rows of `VT` or columns of `U` alike.
fn rotCols(comptime F: type, m: Mat(F), j: usize, c: F, s: F) void {
    const x = m.col(j);
    const y = m.col(j + 1);
    for (0..m.rows) |r| {
        const t = y[r];
        y[r] = c * t - s * x[r];
        x[r] = s * t + c * x[r];
    }
}

/// Fortran's `SIGN(a, b)`: `|a|` with `b`'s sign.
fn sign(comptime F: type, a: F, b: F) F {
    return if (b >= 0) @abs(a) else -@abs(a);
}

/// `DLAS2` (`dlas2.f`): the smaller singular value of `[f g; 0 h]`.
fn las2(comptime F: type, f: F, g: F, h: F) F {
    const fa = @abs(f);
    const ga = @abs(g);
    const ha = @abs(h);
    const fhmn = @min(fa, ha);
    const fhmx = @max(fa, ha);
    if (fhmn == 0) return 0;
    if (ga < fhmx) {
        const as = 1 + fhmn / fhmx;
        const at = (fhmx - fhmn) / fhmx;
        const au = (ga / fhmx) * (ga / fhmx);
        const c = 2 / (@sqrt(as * as + au) + @sqrt(at * at + au));
        return fhmn * c;
    }
    const au = fhmx / ga;
    if (au == 0) return (fhmn * fhmx) / ga;
    const as = 1 + fhmn / fhmx;
    const at = (fhmx - fhmn) / fhmx;
    const c = 1 / (@sqrt(1 + (as * au) * (as * au)) + @sqrt(1 + (at * au) * (at * au)));
    const ssmin = (fhmn * c) * au;
    return ssmin + ssmin;
}

/// `DLASV2` (`dlasv2.f`): the SVD of `[f g; 0 h]`, `ssmax` and `ssmin` signed
/// so the rotations `(csl, snl)`, `(csr, snr)` diagonalize it.
fn lasv2(comptime F: type, f: F, g: F, h: F) struct { ssmin: F, ssmax: F, snr: F, csr: F, snl: F, csl: F } {
    const eps = std.math.floatEps(F) / 2;
    var ft = f;
    var fa = @abs(ft);
    var ht = h;
    var ha = @abs(h);
    var pmax: u2 = 1;
    const swap = ha > fa;
    if (swap) {
        pmax = 3;
        mem.swap(F, &ft, &ht);
        mem.swap(F, &fa, &ha);
    }
    const gt = g;
    const ga = @abs(gt);
    var ssmin: F = undefined;
    var ssmax: F = undefined;
    var clt: F = undefined;
    var crt: F = undefined;
    var slt: F = undefined;
    var srt: F = undefined;
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
            const dd = fa - ha;
            var l = if (dd == fa) 1 else dd / fa;
            const mm_ = gt / ft;
            var t = 2 - l;
            const mm = mm_ * mm_;
            const tt = t * t;
            const s = @sqrt(tt + mm);
            const r = if (l == 0) @abs(mm_) else @sqrt(l * l + mm);
            const a = 0.5 * (s + r);
            ssmin = ha / a;
            ssmax = fa * a;
            if (mm == 0) {
                t = if (l == 0) sign(F, 2, ft) * sign(F, 1, gt) else gt / sign(F, dd, ft) + mm_ / t;
            } else {
                t = (mm_ / (s + t) + mm_ / (r + l)) * (1 + a);
            }
            l = @sqrt(t * t + 4);
            crt = 2 / l;
            srt = t / l;
            clt = (crt + srt * mm_) / a;
            slt = (ht / ft) * srt / a;
        }
    }
    const csl = if (swap) srt else clt;
    const snl = if (swap) crt else slt;
    const csr = if (swap) slt else crt;
    const snr = if (swap) clt else srt;
    const tsign = switch (pmax) {
        1 => sign(F, 1, csr) * sign(F, 1, csl) * sign(F, 1, f),
        2 => sign(F, 1, snr) * sign(F, 1, csl) * sign(F, 1, g),
        else => sign(F, 1, snr) * sign(F, 1, snl) * sign(F, 1, h),
    };
    return .{
        .ssmax = sign(F, ssmax, tsign),
        .ssmin = sign(F, ssmin, tsign * sign(F, 1, f) * sign(F, 1, h)),
        .snr = snr,
        .csr = csr,
        .snl = snl,
        .csl = csl,
    };
}

/// The scaling undone, and `err`.
fn unscaled(comptime F: type, diag: []F, subdiag: []F, max: F, err: error{NoConvergence}) error{NoConvergence} {
    for (diag) |*x| x.* *= max;
    for (subdiag) |*x| x.* *= max;
    return err;
}
