//! The real Schur form of a Hessenberg matrix, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/evd/schur/real_schur.rs`'s `lahqr` (`xLAHQR`'s
//! double-shift QR); `multishift_qr` is not ported yet.

const std = @import("std");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");
const hh = @import("householder.zig");

const Io = std.Io;
const Mat = ops.Mat;

/// `lahqr_eig22`: the eigenvalues of `[a00 a01; a10 a11]` as
/// `.{ .{ re, im }, .{ re, im } }`.
pub fn eig22(comptime F: type, a00_: F, a01_: F, a10_: F, a11_: F) [2][2]F {
    const s = @abs(a00_) + @abs(a01_) + @abs(a10_) + @abs(a11_);
    if (s == 0) return .{ .{ 0, 0 }, .{ 0, 0 } };
    const a00 = a00_ / s;
    const a01 = a01_ / s;
    const a10 = a10_ / s;
    const a11 = a11_ / s;
    const tr = (a00 + a11) * 0.5;
    const det = (a00 - tr) * (a00 - tr) + a01 * a10;
    if (det >= 0) {
        const rtdisc = @sqrt(det);
        return .{ .{ s * (tr + rtdisc), 0 }, .{ s * (tr - rtdisc), 0 } };
    }
    const rtdisc = @sqrt(-det);
    return .{ .{ s * tr, s * rtdisc }, .{ s * tr, -(s * rtdisc) } };
}

/// What `schur22` returns: the standardized block, its eigenvalues and the
/// rotation that made it.
pub fn Schur22(comptime F: type) type {
    return struct { a: F, b: F, c: F, d: F, s1: [2]F, s2: [2]F, cs: F, sn: F };
}

/// `lahqr_schur22` (`xLANV2`): the Schur factorization of the real 2 x 2
/// `[a b; c d]`.
pub fn schur22(comptime F: type, a_: F, b_: F, c_: F, d_: F) Schur22(F) {
    var a = a_;
    var b = b_;
    var c = c_;
    var d = d_;
    const multpl: F = 4;
    const eps = std.math.floatEps(F);
    const safmin = std.math.floatMin(F);
    const safmn2 = @sqrt(safmin / eps);
    const safmx2 = 1 / safmn2;
    var cs: F = undefined;
    var sn: F = undefined;
    if (c == 0) {
        cs = 1;
        sn = 0;
    } else if (b == 0) {
        cs = 0;
        sn = 1;
        std.mem.swap(F, &d, &a);
        b = -c;
        c = 0;
    } else if (a == d and (b > 0) != (c > 0)) {
        cs = 1;
        sn = 0;
    } else {
        var temp = a - d;
        var p = temp * 0.5;
        const bcmax = @max(@abs(b), @abs(c));
        var bcmin = @min(@abs(b), @abs(c));
        if ((b > 0) != (c > 0)) bcmin = -bcmin;
        var scale = @max(@abs(p), bcmax);
        var z = (p / scale) * p + (bcmax / scale) * bcmin;
        if (z >= multpl * eps) {
            var t = @sqrt(scale) * @sqrt(z);
            if (p < 0) t = -t;
            z = p + t;
            a = d + z;
            d -= (bcmax / z) * bcmin;
            const tau = @sqrt(c * c + z * z);
            cs = z / tau;
            sn = c / tau;
            b = b - c;
            c = 0;
        } else {
            var sigma = b + c;
            for (0..20) |_| {
                scale = @max(@abs(temp), @abs(sigma));
                if (scale >= safmx2) {
                    sigma *= safmn2;
                    temp *= safmn2;
                    continue;
                }
                if (scale <= safmn2) {
                    sigma *= safmx2;
                    temp *= safmx2;
                    continue;
                }
                break;
            }
            p = temp * 0.5;
            var tau = @sqrt(sigma * sigma + temp * temp);
            cs = @sqrt((1 + @abs(sigma) / tau) * 0.5);
            sn = -p / (tau * cs);
            if (sigma < 0) sn = -sn;
            const aa = a * cs + b * sn;
            const bb = -a * sn + b * cs;
            const cc = c * cs + d * sn;
            const dd = -c * sn + d * cs;
            a = aa * cs + cc * sn;
            b = bb * cs + dd * sn;
            c = -aa * sn + cc * cs;
            d = -bb * sn + dd * cs;
            temp = (a + d) * 0.5;
            a = temp;
            d = temp;
            if (c != 0 and b != 0 and (b > 0) == (c > 0)) {
                const sab = @sqrt(@abs(b));
                const sac = @sqrt(@abs(c));
                p = if (c > 0) sab * sac else -sab * sac;
                tau = 1 / @sqrt(@abs(b + c));
                a = temp + p;
                d = temp - p;
                b -= c;
                c = 0;
                const cs1 = sab * tau;
                const sn1 = sac * tau;
                temp = cs * cs1 - sn * sn1;
                sn = cs * sn1 + sn * cs1;
                cs = temp;
            }
        }
    }
    if (c != 0) {
        const t = @sqrt(@abs(b)) * @sqrt(@abs(c));
        return .{ .a = a, .b = b, .c = c, .d = d, .s1 = .{ a, t }, .s2 = .{ d, -t }, .cs = cs, .sn = sn };
    }
    return .{ .a = a, .b = b, .c = c, .d = d, .s1 = .{ a, 0 }, .s2 = .{ d, 0 }, .cs = cs, .sn = sn };
}

/// `lahqr_shiftcolumn`: the first column of `(H - s1)(H - s2)` for the
/// leading 2 x 2 or 3 x 3 `h`, scaled.
fn shiftColumn(comptime F: type, h: Mat(F), nr: usize, v: *[3]F, s1: [2]F, s2: [2]F) void {
    if (nr == 2) {
        const s = @abs(h.get(0, 0) - s2[0]) + @abs(s2[1]) + @abs(h.get(1, 0));
        if (s == 0) {
            v[0] = 0;
            v[1] = 0;
            return;
        }
        const h10s = h.get(1, 0) / s;
        v[0] = h10s * h.get(0, 1) + (h.get(0, 0) - s1[0]) * ((h.get(0, 0) - s2[0]) / s) - s1[1] * (s2[1] / s);
        v[1] = h10s * (h.get(0, 0) + h.get(1, 1) - s1[0] - s2[0]);
        return;
    }
    const s = @abs(h.get(0, 0) - s2[0]) + @abs(s2[1]) + @abs(h.get(1, 0)) + @abs(h.get(2, 0));
    if (s == 0) {
        v.* = .{ 0, 0, 0 };
        return;
    }
    const h10s = h.get(1, 0) / s;
    const h20s = h.get(2, 0) / s;
    v[0] = (h.get(0, 0) - s1[0]) * ((h.get(0, 0) - s2[0]) / s) - s1[1] * (s2[1] / s) + h.get(0, 1) * h10s + h.get(0, 2) * h20s;
    v[1] = h10s * (h.get(0, 0) + h.get(1, 1) - s1[0] - s2[0]) + h.get(1, 2) * h20s;
    v[2] = h20s * (h.get(0, 0) + h.get(2, 2) - s1[0] - s2[0]) + h.get(2, 1) * h10s;
}

/// The reflection taking `v[0..nr]` to `beta e0`: `v[0]` becomes `beta`,
/// the rest its essential part; returns LAPACK's `tau` (faer's `1 / tau`).
fn reflect(io: Io, comptime F: type, v: *[3]F, nr: usize) F {
    const info = hh.make(io, F, &v[0], nr - 1, v[1..].ptr, 1, null, 0);
    return 1 / info.tau;
}

/// `lahqr`: the real Schur form of `a`'s rows and columns `ilo..ihi`
/// (Hessenberg), the eigenvalues into `w_re`, `w_im` (`n`); with `want_t`
/// the whole of `T`, else only the active block; `z` (`n x n`), when given,
/// times the rotations. Returns 0, or `ihi`'s current value (the reference's
/// `INFO`) when it runs its iterations out.
pub fn lahqr(io: Io, comptime F: type, want_t: bool, a: Mat(F), z: ?Mat(F), w_re: []F, w_im: []F, ilo: usize, ihi: usize) usize {
    const n = a.rows;
    const nh = ihi - ilo;
    const eps = std.math.floatEps(F);
    const small_num = std.math.floatMin(F) / eps;
    const non_convergence_limit = 10;
    const dat1: F = 0.75;
    const dat2: F = -0.4375;
    if (nh == 0) return 0;
    if (nh == 1) {
        w_re[ilo] = a.get(ilo, ilo);
        w_im[ilo] = 0;
    }
    const itmax = ops.iterLimit(30 * @max(10, nh));
    var k_defl: usize = 0;
    var istop = ihi;
    var istart = ilo;
    var v: [3]F = .{ 0, 0, 0 };
    var iter: usize = 0;
    while (iter <= itmax) : (iter += 1) {
        if (iter == itmax) return istop;
        if (istart + 1 >= istop) {
            if (istart + 1 == istop) {
                w_re[istart] = a.get(istart, istart);
                w_im[istart] = 0;
            }
            break;
        }
        const istart_m: usize = if (want_t) 0 else istart;
        const istop_m: usize = if (want_t) n else istop;
        // a negligible subdiagonal
        var i = istop - 1;
        while (i > istart) : (i -= 1) {
            if (@abs(a.get(i, i - 1)) < small_num) {
                a.at(i, i - 1).* = 0;
                istart = i;
                break;
            }
            var tst = @abs(a.get(i - 1, i - 1)) + @abs(a.get(i, i));
            if (tst == 0) {
                if (i >= ilo + 2) tst += @abs(a.get(i - 1, i - 2));
                if (i + 1 < ihi) tst += @abs(a.get(i + 1, i));
            }
            if (@abs(a.get(i, i - 1)) <= eps * tst) {
                const ab = @max(@abs(a.get(i, i - 1)), @abs(a.get(i - 1, i)));
                const ba = @min(@abs(a.get(i, i - 1)), @abs(a.get(i - 1, i)));
                const aa = @max(@abs(a.get(i, i)), @abs(a.get(i, i) - a.get(i - 1, i - 1)));
                const bb = @min(@abs(a.get(i, i)), @abs(a.get(i, i) - a.get(i - 1, i - 1)));
                const s = aa + ab;
                if (ba * (ab / s) <= @max(eps * (bb * (aa / s)), small_num)) {
                    a.at(i, i - 1).* = 0;
                    istart = i;
                    break;
                }
            }
        }
        if (istart + 2 >= istop) {
            if (istart + 1 == istop) {
                k_defl = 0;
                w_re[istart] = a.get(istart, istart);
                w_im[istart] = 0;
                istop = istart;
                istart = ilo;
                continue;
            }
            if (istart + 2 == istop) {
                const r = schur22(F, a.get(istart, istart), a.get(istart, istart + 1), a.get(istart + 1, istart), a.get(istart + 1, istart + 1));
                const rot: hh.Rotation(F) = .{ .c = r.cs, .s = r.sn };
                a.at(istart, istart).* = r.a;
                a.at(istart, istart + 1).* = r.b;
                a.at(istart + 1, istart).* = r.c;
                a.at(istart + 1, istart + 1).* = r.d;
                w_re[istart] = r.s1[0];
                w_im[istart] = r.s1[1];
                w_re[istart + 1] = r.s2[0];
                w_im[istart + 1] = r.s2[1];
                if (want_t) {
                    // rows istart, istart + 1 right of the block; columns above it
                    if (istart + 2 < istop_m) rot.adjoint().left(istop_m - istart - 2, a.ptr + istart + (istart + 2) * a.ld, a.ld, a.ptr + istart + 1 + (istart + 2) * a.ld, a.ld);
                    rot.right(istart - istart_m, a.col(istart) + istart_m, a.col(istart + 1) + istart_m);
                }
                if (z) |zm| rot.right(zm.rows, zm.col(istart), zm.col(istart + 1));
                k_defl = 0;
                istop = istart;
                istart = ilo;
                continue;
            }
        }
        var a00: F = undefined;
        var a01: F = undefined;
        var a10: F = undefined;
        var a11: F = undefined;
        k_defl += 1;
        if (k_defl % non_convergence_limit == 0) {
            var s = @abs(a.get(istop - 1, istop - 2));
            if (istop > ilo + 2) s += @abs(a.get(istop - 2, istop - 3));
            a00 = dat1 * s + a.get(istop - 1, istop - 1);
            a01 = dat2 * s;
            a10 = s;
            a11 = a00;
        } else {
            a00 = a.get(istop - 2, istop - 2);
            a10 = a.get(istop - 1, istop - 2);
            a01 = a.get(istop - 2, istop - 1);
            a11 = a.get(istop - 1, istop - 1);
        }
        var sh = eig22(F, a00, a01, a10, a11);
        if (sh[0][1] == 0 and sh[1][1] == 0) {
            const last = a.get(istop - 1, istop - 1);
            if (@abs(sh[0][0] - last) <= @abs(sh[1][0] - last)) sh[1] = sh[0] else sh[0] = sh[1];
        }
        // a second small subdiagonal to start the sweep at
        var istart2 = istart;
        if (istart + 3 < istop) {
            var ii = istop - 2;
            while (ii > istart + 1) {
                ii -= 1;
                shiftColumn(F, a.sub(ii, ii, 3, 3), 3, &v, sh[0], sh[1]);
                const tau = reflect(io, F, &v, 3);
                const refsum = tau * a.get(ii, ii - 1) + v[1] * a.get(ii + 1, ii - 1);
                if (@abs(a.get(ii + 1, ii - 1) - refsum * v[1]) + @abs(refsum * v[2]) <= eps * ((@abs(a.get(ii, ii - 1)) + @abs(a.get(ii, ii + 1))) + @abs(a.get(ii + 1, ii + 2)))) {
                    istart2 = ii;
                    break;
                }
            }
        }
        for (istart2..istop - 1) |k| {
            const nr = @min(3, istop - k);
            var t1: F = undefined;
            if (k == istart2) {
                shiftColumn(F, a.sub(k, k, nr, nr), nr, &v, sh[0], sh[1]);
                t1 = reflect(io, F, &v, nr);
                if (k > istart) a.at(k, k - 1).* *= 1 - t1;
            } else {
                v[0] = a.get(k, k - 1);
                v[1] = a.get(k + 1, k - 1);
                if (nr == 3) v[2] = a.get(k + 2, k - 1);
                t1 = reflect(io, F, &v, nr);
                a.at(k, k - 1).* = v[0];
                a.at(k + 1, k - 1).* = 0;
                if (nr == 3) a.at(k + 2, k - 1).* = 0;
            }
            const v2 = v[1];
            const t2 = t1 * v2;
            if (nr == 3) {
                const v3 = v[2];
                const t3 = t1 * v3;
                for (k..istop_m) |j| {
                    const sum = a.get(k, j) + v2 * a.get(k + 1, j) + v3 * a.get(k + 2, j);
                    a.at(k, j).* -= sum * t1;
                    a.at(k + 1, j).* -= sum * t2;
                    a.at(k + 2, j).* -= sum * t3;
                }
                for (istart_m..@min(k + 4, istop)) |j| {
                    const sum = a.get(j, k) + v2 * a.get(j, k + 1) + v3 * a.get(j, k + 2);
                    a.at(j, k).* -= sum * t1;
                    a.at(j, k + 1).* -= sum * t2;
                    a.at(j, k + 2).* -= sum * t3;
                }
                if (z) |zm| for (0..zm.rows) |j| {
                    const sum = zm.get(j, k) + v2 * zm.get(j, k + 1) + v3 * zm.get(j, k + 2);
                    zm.at(j, k).* -= sum * t1;
                    zm.at(j, k + 1).* -= sum * t2;
                    zm.at(j, k + 2).* -= sum * t3;
                };
            } else {
                for (k..istop_m) |j| {
                    const sum = a.get(k, j) + v2 * a.get(k + 1, j);
                    a.at(k, j).* -= sum * t1;
                    a.at(k + 1, j).* -= sum * t2;
                }
                for (istart_m..@min(k + 3, istop)) |j| {
                    const sum = a.get(j, k) + v2 * a.get(j, k + 1);
                    a.at(j, k).* -= sum * t1;
                    a.at(j, k + 1).* -= sum * t2;
                }
                if (z) |zm| for (0..zm.rows) |j| {
                    const sum = zm.get(j, k) + v2 * zm.get(j, k + 1);
                    zm.at(j, k).* -= sum * t1;
                    zm.at(j, k + 1).* -= sum * t2;
                };
            }
        }
    }
    return 0;
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");

test "faer real schur: lahqr, Z T Zᵀ = H with T standardized quasi-triangular" {
    const al = testing.allocator;
    inline for (.{ f64, f32 }) |F| {
        for ([_]usize{ 1, 2, 3, 4, 7, 20, 64, 100 }) |n| for (0..3) |kind| {
            const h0 = try ops.alloc(F, al, n, n);
            defer ops.free(F, al, h0);
            check.fill(F, h0.ptr[0 .. n * n], n + 17 * kind);
            switch (kind) {
                0 => {},
                // a random orthogonal-like rotation pattern: all complex pairs
                1 => for (0..n) |j| for (0..n) |i| {
                    h0.at(i, j).* = if (i == j + 1) 1 else if (j == i + 1) -1 else 0;
                },
                // zero
                else => h0.fill(0),
            }
            for (0..n) |j| for (@min(j + 2, n)..n) |i| {
                h0.at(i, j).* = 0;
            };
            const t = try ops.alloc(F, al, n, n);
            defer ops.free(F, al, t);
            t.copyFrom(h0);
            const zm = try ops.alloc(F, al, n, n);
            defer ops.free(F, al, zm);
            zm.identity();
            const wr = try al.alloc(F, n);
            defer al.free(wr);
            const wi = try al.alloc(F, n);
            defer al.free(wi);
            try testing.expectEqual(@as(usize, 0), lahqr(check.io, F, true, t, zm, wr, wi, 0, n));
            // Z T Zᵀ against H, faer's test_real_schur tolerance (`128 eps sqrt(8 n)`)
            const tol = 128 * std.math.floatEps(F) * @sqrt(8 * @as(F, @floatFromInt(n))) * @max(auxiliary.normMax(F, h0), 1);
            for (0..n) |j| for (0..n) |i| {
                var s: F = 0;
                for (0..n) |p| for (0..n) |q| {
                    if (p <= q + 1) s += zm.get(i, p) * t.get(p, q) * zm.get(j, q);
                };
                try testing.expect(@abs(s - h0.get(i, j)) <= tol);
            };
            // quasi-triangular, standardized, eigenvalues on the blocks
            var k: usize = 0;
            while (k < n) {
                if (k + 1 < n and t.get(k + 1, k) != 0) {
                    try testing.expectEqual(t.get(k, k), t.get(k + 1, k + 1));
                    try testing.expect((t.get(k, k + 1) > 0) != (t.get(k + 1, k) > 0));
                    try testing.expect(wi[k] > 0 and wi[k + 1] == -wi[k]);
                    if (k + 2 < n) try testing.expectEqual(@as(F, 0), t.get(k + 2, k + 1));
                    k += 2;
                } else {
                    try testing.expectEqual(@as(F, 0), wi[k]);
                    try testing.expectEqual(t.get(k, k), wr[k]);
                    k += 1;
                }
            }
        };
    }
}
