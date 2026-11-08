//! The SVD of a real lower bidiagonal matrix, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/svd/bidiag_svd.rs`'s divide and conquer and secular
//! solver, with Reference-LAPACK's `xBDSQR` as `qrAlgorithm` (T126).

const std = @import("std");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const auxiliary = @import("auxiliary.zig");
const hh = @import("householder.zig");
const lartg = auxiliary.lartg;

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;

/// What the solvers return: `NoConvergence` once `qrAlgorithm` has run its
/// iterations out.
pub const Error = ops.Error || error{NoConvergence};

/// `nbits`: the significand's bits.
fn nbits(comptime F: type) usize {
    return std.math.floatMantissaBits(F) + 1;
}

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
        std.mem.swap(F, &ft, &ht);
        std.mem.swap(F, &fa, &ha);
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
                // a statement: as an `if` expression, f80 crashes the self-hosted x86_64 backend (issues.md)
                if (l == 0) {
                    t = sign(F, 2, ft) * sign(F, 1, gt);
                } else {
                    t = gt / sign(F, dd, ft) + mm_ / t;
                }
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

/// `U` as `divide_and_conquer` takes it: whole, or only its first and last
/// rows (`2 x (n + 1)`), or none but scratch for those two rows.
pub fn MatU(comptime F: type) type {
    return union(enum) {
        full: Mat(F),
        two_rows: Mat(F),
        two_rows_storage: Mat(F),
    };
}

/// `divide_and_conquer`: the SVD of the `(n + 1) x n` lower bidiagonal
/// matrix, singular values into `diag`, nonincreasing; `subdiag` is
/// scratch. `v`, when given, is `n x n`.
pub fn divideAndConquer(al: Allocator, io: Io, comptime F: type, diag: []F, subdiag: []F, u: MatU(F), v: ?Mat(F), threshold: usize) Error!void {
    const R = hh.Rotation(F);
    const thr = @max(threshold, 4);
    const n = diag.len;
    if (n < thr) {
        const ua = try ops.alloc(F, al, n + 1, n + 1);
        defer ops.free(F, al, ua);
        ua.identity();
        if (v) |vm| vm.identity();
        var val = subdiag[n - 1];
        subdiag[n - 1] = 0;
        var i = n;
        while (i > 0) {
            i -= 1;
            const rot = R.make(diag[i], val);
            diag[i] = rot.c * diag[i] - rot.s * val;
            if (i > 0) {
                val = rot.s * subdiag[i - 1];
                subdiag[i - 1] = rot.c * subdiag[i - 1];
            }
            rot.adjoint().right(n + 1, ua.col(i), ua.col(n));
        }
        try qrAlgorithm(F, diag, subdiag, ua, v);
        switch (u) {
            .full => |um| um.sub(0, 0, n + 1, n + 1).copyFrom(ua),
            .two_rows => |um| for (0..n + 1) |j| {
                um.at(0, j).* = ua.get(0, j);
                um.at(1, j).* = ua.get(n, j);
            },
            .two_rows_storage => {},
        }
        return;
    }
    const max = @max(ops.absMax(F, diag), ops.absMax(F, subdiag));
    if (max == 0) {
        switch (u) {
            .full => |um| um.identity(),
            .two_rows => |um| {
                um.fill(0);
                um.at(0, 0).* = 1;
                // faer stores 0 here: the last row of the identity is eₙ
                um.at(1, n).* = 1;
            },
            .two_rows_storage => {},
        }
        if (v) |vm| vm.identity();
        return;
    }
    ops.scaleDown(F, diag, max);
    ops.scaleDown(F, subdiag, max);
    const k = n / 2;
    const rem = n - k - 1;
    const alpha = diag[k];
    const beta = subdiag[k];
    {
        const u_top: MatU(F), const u_bot: MatU(F) = switch (u) {
            .full => |um| .{ .{ .full = um.sub(0, 1, k + 1, k + 1) }, .{ .full = um.sub(k + 1, k + 1, rem + 1, rem + 1) } },
            .two_rows, .two_rows_storage => |um| .{ .{ .two_rows = um.sub(0, 0, 2, k + 1) }, .{ .two_rows = um.sub(0, k + 1, 2, rem + 1) } },
        };
        const v_top: ?Mat(F) = if (v) |vm| vm.sub(0, 1, k, k) else null;
        const v_bot: ?Mat(F) = if (v) |vm| vm.sub(k + 1, k + 1, rem, rem) else null;
        try divideAndConquer(al, io, F, diag[0..k], subdiag[0..k], u_top, v_top, threshold);
        try divideAndConquer(al, io, F, diag[k + 1 ..], subdiag[k + 1 ..], u_bot, v_bot, threshold);
        switch (u_top) {
            .two_rows => |um| {
                var i = k;
                while (i > 0) {
                    i -= 1;
                    ops.swapCols(F, um, i, i + 1);
                }
            },
            .full, .two_rows_storage => {},
        }
    }
    if (v) |vm| vm.at(k, 0).* = 1;
    {
        var i = k;
        while (i > 0) {
            i -= 1;
            diag[i + 1] = diag[i];
        }
    }
    const lambda = switch (u) {
        .full => |um| um.get(k, k + 1),
        .two_rows, .two_rows_storage => |um| um.get(1, 0),
    };
    const phi = switch (u) {
        .full => |um| um.get(k + 1, n),
        .two_rows, .two_rows_storage => |um| um.get(0, n),
    };
    const al_ = alpha * lambda;
    const bp = beta * phi;
    const r0 = @sqrt(al_ * al_ + bp * bp);
    const c0: F, const s0: F = if (r0 == 0) .{ 1, 0 } else .{ al_ / r0, bp / r0 };
    const col0 = subdiag;
    diag[0] = r0;
    col0[0] = r0;
    switch (u) {
        .full => |um| {
            for (1..k + 1) |j| col0[j] = alpha * um.get(k, j);
            for (k + 1..n) |j| col0[j] = beta * um.get(k + 1, j);
            for (0..k + 1) |i| {
                const a = um.get(i, k + 1);
                um.at(i, 0).* = c0 * a;
                um.at(i, n).* = -s0 * a;
            }
            for (k + 1..n + 1) |i| {
                const xn = um.get(i, n);
                um.at(i, 0).* = s0 * xn;
                um.at(i, n).* = c0 * xn;
            }
        },
        .two_rows => |um| {
            for (1..k + 1) |j| {
                col0[j] = alpha * um.get(1, j);
                um.at(1, j).* = 0;
            }
            for (k + 1..n) |j| {
                col0[j] = beta * um.get(0, j);
                um.at(0, j).* = 0;
            }
            const q10 = um.get(0, 0);
            const q21 = um.get(1, n);
            um.at(0, 0).* = c0 * q10;
            um.at(0, n).* = -s0 * q10;
            um.at(1, 0).* = s0 * q21;
            um.at(1, n).* = c0 * q21;
        },
        .two_rows_storage => |um| {
            for (1..k + 1) |j| col0[j] = alpha * um.get(1, j);
            for (k + 1..n) |j| col0[j] = beta * um.get(0, j);
        },
    }
    const perm = try al.alloc(usize, n);
    defer al.free(perm);
    for (perm, 0..) |*p, i| p.* = i;
    const jacobi_coeff = try al.alloc(R, n);
    defer al.free(jacobi_coeff);
    @memset(jacobi_coeff, .{ .c = 0, .s = 0 });
    const jacobi_idx = try al.alloc(usize, n);
    defer al.free(jacobi_idx);
    for (jacobi_idx, 0..) |*p, i| p.* = i;
    const jacobi_0i: usize, const jacobi_ij: usize = blk: {
        const transpositions = try al.alloc(usize, n);
        defer al.free(transpositions);
        for (transpositions, 0..) |*p, i| p.* = i;
        break :blk try deflate(al, F, diag, col0, jacobi_coeff, jacobi_idx, transpositions, perm, k);
    };
    const um = try ops.alloc(F, al, n + 1, n + 1);
    defer ops.free(F, al, um);
    um.fill(0);
    const vm: ?Mat(F) = if (v != null) try ops.alloc(F, al, n, n) else null;
    defer if (vm) |x| ops.free(F, al, x);
    if (vm) |x| x.fill(0);
    try svdOfM(al, io, F, um, vm, diag, col0, perm);
    @memset(col0, 0);
    {
        var t = jacobi_ij;
        while (t > 0) {
            t -= 1;
            const rot = jacobi_coeff[jacobi_0i + t];
            const ii = jacobi_idx[jacobi_0i + t];
            const ai = perm[ii - 1];
            const aj = perm[ii];
            rot.left(n, um.ptr + aj, um.ld, um.ptr + ai, um.ld);
            if (vm) |x| rot.left(n, x.ptr + aj, x.ld, x.ptr + ai, x.ld);
        }
        t = jacobi_0i;
        while (t > 0) {
            t -= 1;
            const rot = jacobi_coeff[t];
            const ii = jacobi_idx[t];
            rot.left(n, um.ptr + ii, um.ld, um.ptr, um.ld);
        }
    }
    if (v) |vl| try updateV(al, io, F, vl, vm.?, k, rem);
    switch (u) {
        .two_rows_storage => {},
        .two_rows => |ul| {
            const cu = try ops.alloc(F, al, 2, n + 1);
            defer ops.free(F, al, cu);
            try blas_calls.gemm(al, io, F, .n, .n, 1, ul, um, 0, cu);
            ul.copyFrom(cu);
        },
        .full => |ul| try updateU(al, io, F, ul, um, k, rem),
    }
    for (diag) |*x| x.* *= max;
}

/// `update_v`: `V := V Vm`, by its blocks.
fn updateV(al: Allocator, io: Io, comptime F: type, v: Mat(F), vm: Mat(F), k: usize, rem: usize) ops.Error!void {
    const n = v.rows;
    const cv = try ops.alloc(F, al, n, n);
    defer ops.free(F, al, cv);
    try blas_calls.gemm(al, io, F, .n, .n, 1, v.sub(0, 1, k, k), vm.sub(1, 0, k, n), 0, cv.sub(0, 0, k, n));
    try blas_calls.gemm(al, io, F, .n, .n, 1, v.sub(k + 1, k + 1, rem, rem), vm.sub(k + 1, 0, rem, n), 0, cv.sub(k + 1, 0, rem, n));
    try blas_calls.gemm(al, io, F, .n, .n, 1, v.sub(k, 0, 1, 1), vm.sub(0, 0, 1, n), 0, cv.sub(k, 0, 1, n));
    v.copyFrom(cv);
}

/// `update_u`: `U := U Um`, by its blocks.
fn updateU(al: Allocator, io: Io, comptime F: type, u: Mat(F), um: Mat(F), k: usize, rem: usize) ops.Error!void {
    const n = u.rows - 1;
    const cu = try ops.alloc(F, al, n + 1, n + 1);
    defer ops.free(F, al, cu);
    const c1 = cu.sub(0, 0, k + 1, n + 1);
    const c2 = cu.sub(k + 1, 0, rem + 1, n + 1);
    try blas_calls.gemm(al, io, F, .n, .n, 1, u.sub(0, 0, k + 1, k + 1), um.sub(0, 0, k + 1, n + 1), 0, c1);
    try blas_calls.gemm(al, io, F, .n, .n, 1, u.sub(0, n, k + 1, 1), um.sub(n, 0, 1, n + 1), 1, c1);
    try blas_calls.gemm(al, io, F, .n, .n, 1, u.sub(k + 1, k + 1, rem + 1, rem + 1), um.sub(k + 1, 0, rem + 1, n + 1), 0, c2);
    try blas_calls.gemm(al, io, F, .n, .n, 1, u.sub(k + 1, 0, rem + 1, 1), um.sub(0, 0, 1, n + 1), 1, c2);
    u.copyFrom(cu);
}

/// `compute_svd_of_m`: the SVD of the deflated `M` (`col0` in its first
/// column, `diag` on its diagonal), vectors into `um` and `vm` (zeroed).
fn svdOfM(al: Allocator, io: Io, comptime F: type, um: Mat(F), vm: ?Mat(F), diag: []F, col0: []const F, outer_perm: []const usize) ops.Error!void {
    const n = diag.len;
    diag[0] = 0;
    var actual_n = n;
    while (actual_n > 1 and diag[actual_n - 1] == 0) actual_n -= 1;
    const perm_buf = try al.alloc(usize, actual_n);
    defer al.free(perm_buf);
    var o: usize = 0;
    for (0..actual_n) |i| if (@abs(col0[i]) != 0) {
        perm_buf[o] = i;
        o += 1;
    };
    const perm = perm_buf[0..o];
    const scratch = try al.alloc(F, 2 * o + 4 * n);
    defer al.free(scratch);
    const col0_perm = scratch[0..o];
    const diag_perm = scratch[o .. 2 * o];
    const shifts = scratch[2 * o ..][0..n];
    const mus = scratch[2 * o + n ..][0..n];
    const s = scratch[2 * o + 2 * n ..][0..n];
    const zhat = scratch[2 * o + 3 * n ..][0..n];
    for (perm, 0..) |p, kk| {
        col0_perm[kk] = col0[p];
        diag_perm[kk] = diag[p];
    }
    const sec: Secular(F) = .{ .col0_perm = col0_perm, .diag_perm = diag_perm };
    singularValues(io, F, shifts, mus, s, diag, col0, sec);
    perturbCol0(F, zhat, col0, diag, perm, s, shifts, mus);
    const col_perm = try al.alloc(usize, 2 * actual_n);
    defer al.free(col_perm);
    const col_perm_inv = col_perm[actual_n..];
    for (col_perm[0..actual_n], 0..) |*p, i| p.* = i;
    sortAscending(F, s[0..actual_n], col_perm[0..actual_n]);
    for (col_perm[0..actual_n], 0..) |p, i| col_perm_inv[p] = i;
    singularVectors(io, F, um, vm, zhat, diag, perm, outer_perm, col_perm_inv, actual_n, shifts, mus);
    for (0..actual_n) |idx| diag[idx] = s[actual_n - idx - 1];
    for (actual_n..n) |idx| diag[idx] = s[idx];
}

/// `s` ascending, `perm` moved with it. Faer makes one bubble pass,
/// assuming deflation leaves at most the last value out of place; heavy
/// deflation can leave more (T130, issues.md), so a stable insertion sort,
/// the same order wherever one pass was enough.
fn sortAscending(comptime F: type, s: []F, perm: []usize) void {
    for (1..s.len) |i| {
        var j = i;
        while (j > 0 and s[j - 1] > s[j]) : (j -= 1) {
            std.mem.swap(F, &s[j - 1], &s[j]);
            std.mem.swap(usize, &perm[j - 1], &perm[j]);
        }
    }
}

/// `compute_singular_vectors`.
fn singularVectors(io: Io, comptime F: type, um: Mat(F), vm: ?Mat(F), zhat: []const F, diag: []const F, perm: []const usize, outer_perm: []const usize, col_perm_inv: []const usize, actual_n: usize, shifts: []const F, mus: []const F) void {
    const n = diag.len;
    for (0..n) |k| {
        const actual_k = if (k >= actual_n) k else actual_n - col_perm_inv[k] - 1;
        um.at(n, actual_k).* = 0;
        const u = um.col(actual_k);
        const v: ?[*]F = if (vm) |x| x.col(actual_k) else null;
        if (zhat[k] == 0) {
            u[outer_perm[k]] = 1;
            if (v) |vc| vc[outer_perm[k]] = 1;
            continue;
        }
        const mu = mus[k];
        const shift = shifts[k];
        for (perm) |i| u[outer_perm[i]] = (zhat[i] / ((diag[i] - shift) - mu)) / (diag[i] + (shift + mu));
        const nu = 1 / blas_calls.nrm2(io, F, n, u, 1);
        for (u[0..n]) |*x| x.* *= nu;
        if (v) |vc| {
            vc[outer_perm[0]] = -1;
            for (perm[1..]) |i| vc[outer_perm[i]] = ((diag[i] * zhat[i]) / ((diag[i] - shift) - mu)) / (diag[i] + (shift + mu));
            const nv = 1 / blas_calls.nrm2(io, F, n, vc, 1);
            for (vc[0..n]) |*x| x.* *= nv;
        }
    }
    um.at(n, n).* = 1;
}

/// `perturb_col0`: `col0` recomputed from the computed singular values,
/// so the vectors come out orthogonal (Gu and Eisenstat).
fn perturbCol0(comptime F: type, zhat: []F, col0: []const F, diag: []const F, perm: []const usize, s: []const F, shifts: []const F, mus: []const F) void {
    const n = diag.len;
    const m = perm.len;
    if (m == 0) {
        @memset(zhat, 0);
        return;
    }
    const last_idx = perm[m - 1];
    for (0..n) |k| {
        if (col0[k] == 0) {
            zhat[k] = 0;
            continue;
        }
        const dk = diag[k];
        var prod = (s[last_idx] + dk) * (mus[last_idx] + (shifts[last_idx] - dk));
        for (perm, 0..) |i, l| {
            if (i == k) continue;
            if (i >= k and l == 0) {
                prod = 0;
                break;
            }
            const j = if (i < k) i else if (l > 0) perm[l - 1] else i;
            const term = ((s[j] + dk) / (diag[i] + dk)) * ((mus[j] + (shifts[j] - dk)) / (diag[i] - dk));
            prod = prod * term;
        }
        const tmp = @sqrt(prod);
        zhat[k] = if (col0[k] > 0) tmp else -tmp;
    }
}

/// The secular function's data: the deflated `col0` and `diag` entries.
fn Secular(comptime F: type) type {
    return struct {
        col0_perm: []const F,
        diag_perm: []const F,

        const Self = @This();

        /// `secular_eq`: `1 + Σ c² / ((d - shift - mu) (d + shift + mu))`.
        fn eval(e: Self, shift: F, mu: F) F {
            var res: F = 1;
            for (e.col0_perm, e.diag_perm) |c, d| res += (c / ((d - shift) - mu)) * (c / ((d + shift) + mu));
            return res;
        }

        /// `batch_secular_eq`: four at once, its own way.
        fn batch(e: Self, shift: [4]F, mu: [4]F) [4]F {
            var res: [4]F = @splat(1);
            for (e.col0_perm, e.diag_perm) |c, d| {
                for (&res, mu, shift) |*r, m, s| r.* += (c * c) / (((d - s) - m) * ((d + s) + m));
            }
            return res;
        }
    };
}

/// `compute_singular_values`: each root as `shift + mu`, `shift` the nearer
/// pole, for accuracy.
fn singularValues(io: Io, comptime F: type, shifts: []F, mus: []F, s: []F, diag: []const F, col0: []const F, sec: Secular(F)) void {
    const n = diag.len;
    var actual_n = n;
    while (actual_n > 1 and col0[actual_n - 1] == 0) actual_n -= 1;
    for (0..n) |k| {
        s[k] = 0;
        shifts[k] = 0;
        mus[k] = 0;
        if (col0[k] == 0 or actual_n == 1) {
            s[k] = if (k == 0) col0[0] else diag[k];
            shifts[k] = s[k];
            mus[k] = 0;
            continue;
        }
        const last_k = k == actual_n - 1;
        const left = diag[k];
        const right = if (last_k) diag[actual_n - 1] + blas_calls.nrm2(io, F, n, @constCast(col0.ptr), 1) else blk: {
            var l = k + 1;
            while (col0[l] == 0) l += 1;
            break :blk diag[l];
        };
        const r = rootFinder(F, sec, left, right, last_k);
        s[k] = r[0] + r[1];
        shifts[k] = r[0];
        mus[k] = r[1];
    }
}

/// `secular_eq_root_finder`: the root in `(left, right)` of the increasing
/// secular function `sec` (its `eval(shift, mu)` and four at once,
/// `batch`), as `.{ shift, mu }`: bracketing and bisection at geometric
/// midpoints, then the secant method, then bisection if it leaves the
/// bracket. The tridiagonal eigensolver shares it (`tridiag_evd.zig`).
pub fn rootFinder(comptime F: type, sec: anytype, left: F, right: F, last: bool) [2]F {
    const two: F = 2;
    const half: F = 0.5;
    const epsilon = std.math.floatEps(F);
    const inf = std.math.inf(F);
    const mid = left + (right - left) * half;
    const f0 = sec.batch(.{ 0, left, left, right }, .{ mid, if (last) right - left else (right - left) * half, (right - left) * half, (left - right) * half });
    var f_mid = f0[0];
    const f_max = f0[1];
    const f_mid_left_shift = f0[2];
    const f_mid_right_shift = f0[3];
    var shift: F, const mu0: F = if (last or f_mid > 0) .{ left, (right - left) * half } else .{ right, (left - right) * half };
    if (f_mid_left_shift <= 0 and f_mid_right_shift > 0) return .{ shift, mu0 };
    if (!last) {
        if (shift == left) {
            if (f_mid_left_shift < 0) {
                shift = right;
                f_mid = f_mid_right_shift;
            }
        } else if (f_mid_right_shift > 0) {
            shift = left;
            f_mid = f_mid_left_shift;
        }
    }
    var left_shifted: F = undefined;
    var f_left: F = undefined;
    var right_shifted: F = undefined;
    var f_right: F = undefined;
    if (shift == left) {
        left_shifted = 0;
        f_left = -inf;
        right_shifted = if (last) right - left else (right - left) * half;
        f_right = if (last) f_max else f_mid;
    } else {
        left_shifted = (left - right) * half;
        f_left = f_mid;
        right_shifted = 0;
        f_right = inf;
    }
    var iteration_count: usize = 0;
    var f_prev = f_mid;
    const half0 = half;
    const half1 = half0 * half0;
    const half2 = half1 * half1;
    const half3 = half2 * half2;
    const base = if (shift == left) right_shifted else left_shifted;
    const mu_values: [4]F = .{ base * half3, base * half2, base * half1, base * half0 };
    const f_values = sec.batch(@splat(shift), mu_values);
    if (shift == left) {
        var i: usize = 0;
        for (mu_values, f_values, 0..) |m, f, idx| if (f < 0) {
            left_shifted = m;
            f_left = f;
            i = idx + 1;
        };
        if (i < 4) {
            right_shifted = mu_values[i];
            f_right = f_values[i];
        }
    } else {
        var i: usize = 0;
        for (mu_values, f_values, 0..) |m, f, idx| if (f > 0) {
            right_shifted = m;
            f_right = f;
            i = idx + 1;
        };
        if (i < 4) {
            left_shifted = mu_values[i];
            f_left = f_values[i];
        }
    }
    while (right_shifted - left_shifted > two * epsilon * @max(@abs(left_shifted), @abs(right_shifted))) {
        const mid_arith = (left_shifted + right_shifted) * half;
        var mid_geo = @sqrt(@abs(left_shifted)) * @sqrt(@abs(right_shifted));
        if (left_shifted < 0) mid_geo = -mid_geo;
        const mid_shifted = if (mid_geo == 0) mid_arith else mid_geo;
        const fm = sec.eval(shift, mid_shifted);
        if (fm == 0) {
            return .{ shift, mid_shifted };
        } else if (fm > 0) {
            right_shifted = mid_shifted;
            f_prev = f_right;
            f_right = fm;
        } else {
            left_shifted = mid_shifted;
            f_prev = f_left;
            f_left = fm;
        }
        if (iteration_count == 4) break;
        iteration_count += 1;
    }
    const args: [4]F = if (left_shifted == 0)
        .{ right_shifted * two, right_shifted, f_prev, f_right }
    else if (right_shifted == 0)
        .{ left_shifted * two, left_shifted, f_prev, f_left }
    else
        .{ left_shifted, right_shifted, f_left, f_right };
    const sc = secant(F, sec, shift, left, right, args[0], args[1], args[2], args[3]);
    var mu_cur = sc.mu_cur;
    if (sc.left) |lc| if (sc.right) |rc| if (lc < rc) {
        if (lc > left_shifted) left_shifted = lc;
        if (rc < right_shifted) right_shifted = rc;
    };
    if (sc.use_bisection) {
        while (right_shifted - left_shifted > two * epsilon * @max(@abs(left_shifted), @abs(right_shifted))) {
            const m = (left_shifted + right_shifted) * half;
            const fm = sec.eval(shift, m);
            if (fm == 0) {
                break;
            } else if (fm > 0) {
                right_shifted = m;
            } else {
                left_shifted = m;
            }
        }
        mu_cur = (left_shifted + right_shifted) * half;
    }
    return .{ shift, mu_cur };
}

fn SecantResult(comptime F: type) type {
    return struct { use_bisection: bool, mu_cur: F, left: ?F, right: ?F };
}

/// `rootFinder`'s `secant` closure.
fn secant(comptime F: type, sec: anytype, shift: F, left: F, right: F, mu_cur_in: F, mu_prev_in: F, f_cur_in: F, f_prev_in: F) SecantResult(F) {
    const epsilon = std.math.floatEps(F);
    var mu_cur = mu_cur_in;
    var mu_prev = mu_prev_in;
    var f_cur = f_cur_in;
    var f_prev = f_prev_in;
    if (@abs(f_prev) < @abs(f_cur)) {
        std.mem.swap(F, &f_prev, &f_cur);
        std.mem.swap(F, &mu_prev, &mu_cur);
    }
    var left_candidate: ?F = null;
    var right_candidate: ?F = null;
    var use_bisection = false;
    const same_sign = (f_prev > 0) == (f_cur > 0);
    if (!same_sign) {
        if (mu_cur < mu_prev) {
            left_candidate = mu_cur;
            right_candidate = mu_prev;
        } else {
            left_candidate = mu_prev;
            right_candidate = mu_cur;
        }
    }
    while (f_cur != 0 and @abs(mu_cur - mu_prev) > 8 * epsilon * @max(@abs(mu_cur), @abs(mu_prev)) and @abs(f_cur - f_prev) > epsilon and !use_bisection) {
        const a = (f_cur - f_prev) * (mu_prev * mu_cur) / (mu_prev - mu_cur);
        const b = f_cur - a / mu_cur;
        const mu_zero = -a / b;
        const f_zero = sec.eval(shift, mu_zero);
        if (f_zero < 0) {
            left_candidate = mu_zero;
        } else {
            right_candidate = mu_zero;
        }
        mu_prev = mu_cur;
        f_prev = f_cur;
        mu_cur = mu_zero;
        f_cur = f_zero;
        if (shift == left and (mu_cur < 0 or mu_cur > right - left)) use_bisection = true;
        if (shift == right and (mu_cur > 0 or mu_cur < left - right)) use_bisection = true;
        if (@abs(f_cur) > @abs(f_prev)) {
            var k: F = 1;
            for (0..4) |_| {
                const mu_opposite = -a / (k * f_zero + b);
                const f_opposite = sec.eval(shift, mu_opposite);
                if (f_zero < 0 and f_opposite >= 0) {
                    right_candidate = mu_opposite;
                    break;
                }
                if (f_zero > 0 and f_opposite <= 0) {
                    left_candidate = mu_opposite;
                    break;
                }
                k *= 2;
            }
            use_bisection = true;
        }
    }
    // Faer stops here with whatever the secant reached when `f` stops changing
    // by more than `eps`, which happens far from a root that lies within
    // rounding of a pole (a tiny `col0` entry): then `mu` is wrong by orders
    // of magnitude, and `perturbCol0` with it. Bisect instead (issues.md).
    if (!use_bisection and f_cur != 0 and !(@abs(f_cur - f_prev) > epsilon) and
        @abs(mu_cur - mu_prev) > 8 * epsilon * @max(@abs(mu_cur), @abs(mu_prev))) use_bisection = true;
    return .{ .use_bisection = use_bisection, .mu_cur = mu_cur, .left = left_candidate, .right = right_candidate };
}

/// `deflate`: Gu and Eisenstat's deflations; returns how many rotations
/// each kind recorded (`jacobi_0i` of `deflation_43`, then `jacobi_ij` of
/// `deflation_44`), and `perm` the order `svdOfM` reads.
fn deflate(al: Allocator, comptime F: type, diag: []F, col0: []F, jacobi_coeff: []hh.Rotation(F), jacobi_idx: []usize, transpositions: []usize, perm: []usize, k: usize) Allocator.Error![2]usize {
    const n = diag.len;
    var jacobi_0i: usize = 0;
    var jacobi_ij: usize = 0;
    const max_diag = ops.absMax(F, diag[1..]);
    const max_col0 = ops.absMax(F, col0);
    const mx = @max(max_diag, max_col0);
    const eps = std.math.floatEps(F);
    const sml = std.math.floatMin(F);
    const eps_strict = @max(eps * max_diag, sml);
    const eps_coarse = 8 * eps * mx;
    if (diag[0] < eps_coarse) {
        diag[0] = eps_coarse;
        col0[0] = eps_coarse;
    }
    for (1..n) |i| if (@abs(col0[i]) < eps_strict) {
        col0[i] = 0;
    };
    for (1..n) |i| if (diag[i] < eps_coarse) {
        if (deflation43(F, diag, col0, i)) |rot| {
            jacobi_coeff[jacobi_0i] = rot;
            jacobi_idx[jacobi_0i] = i;
            jacobi_0i += 1;
        }
    };
    var total_deflation = true;
    for (1..n) |i| if (!(@abs(col0[i]) < sml)) {
        total_deflation = false;
        break;
    };
    perm[0] = 0;
    var p: usize = 1;
    for (1..n) |i| if (@abs(diag[i]) < sml) {
        perm[p] = i;
        p += 1;
    };
    {
        var i: usize = 1;
        var j: usize = k + 1;
        for (perm[p..]) |*slot| {
            if (i >= k + 1) {
                slot.* = j;
                j += 1;
            } else if (j >= n) {
                slot.* = i;
                i += 1;
            } else if (diag[i] < diag[j]) {
                slot.* = j;
                j += 1;
            } else {
                slot.* = i;
                i += 1;
            }
        }
    }
    if (total_deflation) {
        for (1..n) |i| {
            const pi = perm[i];
            if (@abs(diag[pi]) < sml or diag[pi] > diag[0]) {
                perm[i - 1] = perm[i];
            } else {
                perm[i - 1] = 0;
                break;
            }
        }
    }
    const real = try al.alloc(usize, 2 * n);
    defer al.free(real);
    const real_ind = real[0..n];
    const real_col = real[n..];
    for (real_ind, real_col, 0..) |*a, *b, i| {
        a.* = i;
        b.* = i;
    }
    for ((if (total_deflation) @as(usize, 0) else 1)..n) |i| {
        const pi = perm[n - if (total_deflation) i + 1 else i];
        const j = real_col[pi];
        std.mem.swap(F, &diag[i], &diag[j]);
        if (i != 0 and j != 0) std.mem.swap(F, &col0[i], &col0[j]);
        transpositions[i] = j;
        const real_i = real_ind[i];
        real_col[real_i] = j;
        real_col[pi] = i;
        real_ind[j] = real_i;
        real_ind[i] = pi;
    }
    col0[0] = diag[0];
    for (perm, 0..) |*x, i| x.* = i;
    for (transpositions, 0..) |j, i| std.mem.swap(usize, &perm[i], &perm[j]);
    var i = n - 1;
    while (i > 0 and (@abs(diag[i]) < sml or @abs(col0[i]) < sml)) i -= 1;
    while (i > 1) {
        const ib = i - 1;
        if (diag[i] - diag[ib] < eps_strict) {
            if (deflation44(F, diag, col0, ib, i)) |rot| {
                jacobi_coeff[jacobi_0i + jacobi_ij] = rot;
                jacobi_idx[jacobi_0i + jacobi_ij] = i;
                jacobi_ij += 1;
            }
        }
        i = ib;
    }
    return .{ jacobi_0i, jacobi_ij };
}

/// `deflation_43`: a tiny `diag[i]` folds `col0[i]` into `col0[0]`.
fn deflation43(comptime F: type, diag: []F, col0: []F, i: usize) ?hh.Rotation(F) {
    const p = col0[0];
    const q = col0[i];
    if (p == 0 and q == 0) return null;
    const rot = hh.Rotation(F).make(p, q);
    const r = rot.c * p - rot.s * q;
    col0[0] = r;
    diag[0] = r;
    col0[i] = 0;
    return rot;
}

/// `deflation_44`: equal `diag[i]` and `diag[j]` fold `col0[j]` into
/// `col0[i]`.
fn deflation44(comptime F: type, diag: []F, col0: []F, i: usize, j: usize) ?hh.Rotation(F) {
    const p = col0[i];
    const q = col0[j];
    if (p == 0 and q == 0) {
        diag[i] = diag[j];
        return null;
    }
    const rot = hh.Rotation(F).make(p, q);
    const r = rot.c * p - rot.s * q;
    col0[i] = r;
    col0[j] = 0;
    diag[i] = diag[j];
    return rot;
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");

/// `U [S; 0] Vᵀ` against the `(n + 1) x n` (or with `square`, `n x n`)
/// lower bidiagonal matrix, within faer's `test_divide_and_conquer` (and
/// `test_qr_algorithm`) tolerance: `128 eps` times `max(|d|, |e|) sqrt n`,
/// times 10 for divide and conquer; and `Uᵀ U = I`, `Vᵀ V = I` within it.
fn expectBidiagSvd(comptime F: type, diag: []const F, sub: []const F, d: []const F, u: Mat(F), v: Mat(F), square: bool, factor: F) !void {
    const n = diag.len;
    const rows = if (square) n else n + 1;
    const scale = @max(ops.absMax(F, diag), ops.absMax(F, sub));
    const tol = 128 * std.math.floatEps(F) * @max(scale, 1) * @sqrt(@as(F, @floatFromInt(n))) * factor;
    for (0..n) |j| for (0..rows) |i| {
        var x: F = 0;
        for (0..n) |kk| x += u.get(i, kk) * d[kk] * v.get(j, kk);
        const want: F = if (i == j) diag[j] else if (i == j + 1) sub[j] else 0;
        if (!(@abs(x - want) <= tol)) {
            check.say("bidiag svd n={d}: ({d},{d}) = {e}, want {e}\n", .{ n, i, j, x, want });
            return error.TestUnexpectedResult;
        }
    };
    for (0..n) |j| for (0..n) |i| {
        var x: F = 0;
        for (0..n) |kk| x += v.get(kk, i) * v.get(kk, j);
        try testing.expect(@abs(x - @as(F, if (i == j) 1 else 0)) <= tol);
    };
    for (0..rows) |j| for (0..rows) |i| {
        var x: F = 0;
        for (0..rows) |kk| x += u.get(kk, i) * u.get(kk, j);
        try testing.expect(@abs(x - @as(F, if (i == j) 1 else 0)) <= tol);
    };
    for (1..n) |i| try testing.expect(d[i - 1] >= d[i]);
    for (d) |x| try testing.expect(x >= 0);
}

/// A bidiagonal with entries spread over several magnitudes, and clusters.
fn testBidiag(comptime F: type, diag: []F, sub: []F, seed: u64) void {
    var prng = std.Random.DefaultPrng.init(seed);
    const r = prng.random();
    for (diag, sub, 0..) |*d, *e, i| {
        d.* = (r.float(F) * 2 - 1) * std.math.pow(F, 10, @floatFromInt(@as(i32, @intCast(i % 5)) - 2));
        e.* = (r.float(F) * 2 - 1);
        if (i % 7 == 3) d.* = diag[i - 1];
    }
}

test "faer bidiag svd: qr_algorithm and divide_and_conquer reconstruct B, both thresholds" {
    const al = testing.allocator;
    inline for (.{ f64, f32 }) |F| {
        for ([_]usize{ 1, 2, 3, 7, 20, 64, 129 }) |n| {
            const diag = try al.alloc(F, n);
            defer al.free(diag);
            const sub = try al.alloc(F, n);
            defer al.free(sub);
            testBidiag(F, diag, sub, n);
            const d = try al.alloc(F, n);
            defer al.free(d);
            const e = try al.alloc(F, n);
            defer al.free(e);
            // qr_algorithm on the square part
            {
                @memcpy(d, diag);
                @memcpy(e, sub);
                e[n - 1] = 0;
                const sq = try al.dupe(F, e);
                defer al.free(sq);
                const u = try ops.alloc(F, al, n, n);
                defer ops.free(F, al, u);
                const v = try ops.alloc(F, al, n, n);
                defer ops.free(F, al, v);
                u.identity();
                v.identity();
                try qrAlgorithm(F, d, e, u, v);
                for (e) |x| try testing.expectEqual(@as(F, 0), x);
                try expectBidiagSvd(F, diag, sq, d, u, v, true, 1);
            }
            // divide and conquer, faer's threshold of 4 and 8, and 128
            for ([_]usize{ 4, 8, 128 }) |thr| {
                @memcpy(d, diag);
                @memcpy(e, sub);
                const u = try ops.alloc(F, al, n + 1, n + 1);
                defer ops.free(F, al, u);
                const v = try ops.alloc(F, al, n, n);
                defer ops.free(F, al, v);
                try divideAndConquer(al, check.io, F, d, e, .{ .full = u }, v, thr);
                expectBidiagSvd(F, diag, sub, d, u, v, false, 10) catch |err| {
                    check.say("{s} n={d} threshold {d}\n", .{ @typeName(F), n, thr });
                    return err;
                };
                // the values alone, and without U, agree
                const d2 = try al.dupe(F, diag);
                defer al.free(d2);
                const e2 = try al.dupe(F, sub);
                defer al.free(e2);
                const u_bot = try ops.alloc(F, al, 2, n + 1);
                defer ops.free(F, al, u_bot);
                try divideAndConquer(al, check.io, F, d2, e2, .{ .two_rows_storage = u_bot }, null, thr);
                for (d, d2) |x, y| try testing.expect(@abs(x - y) <= 1024 * std.math.floatEps(F) * @max(d[0], 1));
            }
        }
    }
}

test "faer bidiag svd: a zero block deflates, with or without U; zink's bidiagonal keeps its smallest value" {
    const al = testing.allocator;
    // the zero block is the second half of the first half, so its `TwoRows`
    // last row (faer: 0, the identity: eₙ) becomes the coupling `lambda`
    const n = 24;
    var diag: [n]f64 = undefined;
    var sub: [n]f64 = undefined;
    for (0..n) |i| {
        const zero = i >= 7 and i < 12;
        diag[i] = if (zero) 0 else 1 + @as(f64, @floatFromInt(i));
        sub[i] = if (zero) 0 else 0.5;
    }
    var d = diag;
    var e = sub;
    const u = try ops.alloc(f64, al, n + 1, n + 1);
    defer ops.free(f64, al, u);
    const v = try ops.alloc(f64, al, n, n);
    defer ops.free(f64, al, v);
    try divideAndConquer(al, check.io, f64, &d, &e, .{ .full = u }, v, 4);
    try expectBidiagSvd(f64, &diag, &sub, &d, u, v, false, 10);
    var d2 = diag;
    var e2 = sub;
    const u_bot = try ops.alloc(f64, al, 2, n + 1);
    defer ops.free(f64, al, u_bot);
    try divideAndConquer(al, check.io, f64, &d2, &e2, .{ .two_rows_storage = u_bot }, null, 4);
    // within a few ulps; faer's 0 is off by 1.5e-13 at d[5]
    for (d, d2) |x, y| try testing.expectApproxEqAbs(x, y, 8 * std.math.floatEps(f64) * d[0]);
    // faer's `test_zink` (svd/mod.rs)
    var zd = [_]f64{ -9.931833701529301, -10.920807536026027, -52.33647796311243, 2.3685025127736967, 2.421701994236093, -0.5051763005624579, -0.04808263896606017, -0.003875251886338955, -0.0006413264967716465, -0.003381944152463707, 2.981152313236375e-5, 5.4290648208388795e-6, -6.329275972084404e-7, -6.879142344209158e-7, -5.265228263479126e-9, -2.941999902335516e-9, -1.3060984997930294e-10, 7.07516117218088e-12, 1.8657003929029376e-12, -6.216080089659131e-14 };
    var zs = [_]f64{ -57.8029649868477, 17.67263066467847, 8.884153814270894, -9.01998231080713, -1.028638150814966, 0.22247719217200435, 0.016389886745811315, -0.004090989452162578, 0.00036818904090536926, -0.0031394146217732367, -7.571300829706796e-6, 3.0045718718618155e-6, 2.1329796886727743e-6, 9.259701025627789e-8, 2.2291214755992877e-9, -2.3017207713252894e-9, 6.807967994979358e-11, 2.1677299575405587e-12, -3.07282771050034e-13, 0.0 };
    const zu = try ops.alloc(f64, al, 2, zd.len + 1);
    defer ops.free(f64, al, zu);
    try divideAndConquer(al, check.io, f64, &zd, &zs, .{ .two_rows_storage = zu }, null, 8);
    try testing.expect(zd[zd.len - 1] != 0);
}

test "faer bidiag svd: deflation 4.3 and 4.4, faer's cases" {
    // `test_deflation43`, `test_deflation44`, `test_both_deflation`
    const al = testing.allocator;
    const cases = [_]struct { diag: [7]f32, want: [2]usize }{
        .{ .diag = .{ 1, 5, 3, 1e-7, 4, 2, 2e-7 }, .want = .{ 2, 0 } },
        .{ .diag = .{ 1, 5, 3, 1, 4, 2, 1 }, .want = .{ 0, 1 } },
        .{ .diag = .{ 1, 5, 3, 2, 4, 2, 0 }, .want = .{ 1, 1 } },
    };
    for (cases) |c| {
        var diag = c.diag;
        var col: [7]f32 = @splat(1);
        var coeff: [7]hh.Rotation(f32) = undefined;
        var idx: [7]usize = undefined;
        var perm: [7]usize = undefined;
        var tr: [7]usize = .{ 0, 1, 2, 3, 4, 5, 6 };
        const got = try deflate(al, f32, &diag, &col, &coeff, &idx, &tr, &perm, 3);
        try testing.expectEqual(c.want, got);
        // M' in the permuted places, then the rotations undone, is M: the
        // diagonal and a first column of ones
        var m: [7][7]f32 = @splat(@splat(0));
        for (0..7) |i| {
            m[perm[i]][perm[i]] = diag[i];
            m[perm[i]][perm[0]] = col[i];
        }
        var t = got[1];
        while (t > 0) {
            t -= 1;
            const rot = coeff[got[0] + t];
            const pi = perm[idx[got[0] + t] - 1];
            const pj = perm[idx[got[0] + t]];
            rot.left(7, &m[pj], 1, &m[pi], 1);
            rot.left(7, @as([*]f32, @ptrCast(&m)) + pj, 7, @as([*]f32, @ptrCast(&m)) + pi, 7);
        }
        t = got[0];
        while (t > 0) {
            t -= 1;
            coeff[t].left(7, &m[idx[t]], 1, &m[0], 1);
        }
        for (0..7) |i| for (0..7) |j| {
            const want: f32 = if (j == 0) 1 else if (i == j) c.diag[i] else 0;
            try testing.expect(@abs(m[i][j] - want) <= 128 * std.math.floatEps(f32) * @max(1, @abs(want)));
        };
    }
}

test "faer bidiag svd: the merged values sort fully, several out of place, ties kept in order (T130)" {
    // one bubble pass, faer's, would leave 1 ahead of 0.5 here
    var s = [_]f64{ 3, 1, 0.5, 2, 2, 0 };
    var perm = [_]usize{ 0, 1, 2, 3, 4, 5 };
    sortAscending(f64, &s, &perm);
    try testing.expectEqualSlices(f64, &.{ 0, 0.5, 1, 2, 2, 3 }, &s);
    try testing.expectEqualSlices(usize, &.{ 5, 2, 1, 3, 4, 0 }, &perm);
    // already sorted, and faer's case (one value a step ahead): as its pass leaves them
    var t = [_]f64{ 0, 1, 2 };
    var q = [_]usize{ 0, 1, 2 };
    sortAscending(f64, &t, &q);
    try testing.expectEqualSlices(usize, &.{ 0, 1, 2 }, &q);
    var u = [_]f64{ 0.5, 1, 3, 2 };
    var r = [_]usize{ 0, 1, 2, 3 };
    sortAscending(f64, &u, &r);
    try testing.expectEqualSlices(usize, &.{ 0, 1, 3, 2 }, &r);
}
