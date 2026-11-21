//! The complex Schur form of a Hessenberg matrix, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/evd/schur/complex_schur.rs`'s `lahqr` (`xLAHQR`'s
//! single-shift QR, with faer's exceptional shift); `multishift_qr` is not
//! ported yet.

const std = @import("std");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");
const hh = @import("householder.zig");

const Mat = ops.Mat;
const Float = ops.Float;

/// `lahqr_eig22`: the eigenvalues of `[a00 a01; a10 a11]`.
pub fn eig22(comptime T: type, a00_: T, a01_: T, a10_: T, a11_: T) [2]T {
    const s = ops.abs1(T, a00_) + ops.abs1(T, a01_) + ops.abs1(T, a10_) + ops.abs1(T, a11_);
    if (s == 0) return .{ ops.zero(T), ops.zero(T) };
    const s_inv = 1 / s;
    const a00 = ops.mulReal(T, a00_, s_inv);
    const a01 = ops.mulReal(T, a01_, s_inv);
    const a10 = ops.mulReal(T, a10_, s_inv);
    const a11 = ops.mulReal(T, a11_, s_inv);
    const tr = ops.mulReal(T, ops.add(T, a00, a11), 0.5);
    const d = ops.sub(T, a00, tr);
    const det = ops.add(T, ops.mul(T, d, d), ops.mul(T, a01, a10));
    const rtdisc = ops.sqrt(T, det);
    return .{ ops.mulReal(T, ops.add(T, tr, rtdisc), s), ops.mulReal(T, ops.sub(T, tr, rtdisc), s) };
}

/// `lahqr`: the Schur form of `a`'s rows and columns `ilo..ihi`
/// (Hessenberg), the eigenvalues into `w` (`n`); with `want_t` the whole of
/// `T`, else only the active block; `z` (`n x n`), when given, times the
/// rotations. Returns 0, or `ihi`'s current value (the reference's `INFO`)
/// when it runs its iterations out.
pub fn lahqr(comptime T: type, want_t: bool, a: Mat(T), z: ?Mat(T), w: []T, ilo: usize, ihi: usize) usize {
    const F = Float(T);
    const R = hh.ComplexRotation(T);
    const n = a.rows;
    const nh = ihi - ilo;
    const eps = std.math.floatEps(F);
    const smlnum = std.math.floatMin(F) / eps;
    const non_convergence_limit = 10;
    const dat1: F = 0.75;
    const dat2: F = -0.4375;
    if (nh == 0) return 0;
    if (nh == 1) w[ilo] = a.get(ilo, ilo);
    const itmax = ops.iterLimit(@max(30, (std.math.floatMantissaBits(F) + 1) / 2) *| @max(10, nh) *| nh);
    var k_defl: usize = 0;
    var istop = ihi;
    var istart = ilo;
    var iter: usize = 0;
    while (iter < itmax) : (iter += 1) {
        if (iter + 1 == itmax) return istop;
        if (istart + 1 >= istop) {
            if (istart + 1 == istop) w[istart] = a.get(istart, istart);
            break;
        }
        const istart_m: usize = if (want_t) 0 else istart;
        const istop_m: usize = if (want_t) n else istop;
        var i = istop - 1;
        while (i > istart) : (i -= 1) {
            if (ops.abs1(T, a.get(i, i - 1)) < smlnum) {
                a.at(i, i - 1).* = ops.zero(T);
                istart = i;
                break;
            }
            var tst = ops.abs1(T, a.get(i - 1, i - 1)) + ops.abs1(T, a.get(i, i));
            if (tst == 0) {
                if (i >= ilo + 2) tst += ops.abs1(T, a.get(i - 1, i - 2));
                if (i + 1 < ihi) tst += ops.abs1(T, a.get(i + 1, i));
            }
            if (ops.abs1(T, a.get(i, i - 1)) <= eps * tst) {
                const ab = @max(ops.abs1(T, a.get(i, i - 1)), ops.abs1(T, a.get(i - 1, i)));
                const ba = @min(ops.abs1(T, a.get(i, i - 1)), ops.abs1(T, a.get(i - 1, i)));
                const diff = ops.abs1(T, ops.sub(T, a.get(i, i), a.get(i - 1, i - 1)));
                const aa = @max(ops.abs1(T, a.get(i, i)), diff);
                const bb = @min(ops.abs1(T, a.get(i, i)), diff);
                const s = aa + ab;
                if (ba * (ab / s) <= @max(eps * (bb * (aa / s)), smlnum)) {
                    a.at(i, i - 1).* = ops.zero(T);
                    istart = i;
                    break;
                }
            }
        }
        if (istart + 1 >= istop) {
            k_defl = 0;
            w[istart] = a.get(istart, istart);
            istop = istart;
            istart = ilo;
            continue;
        }
        var a00: T = undefined;
        var a01: T = undefined;
        var a10: T = undefined;
        var a11: T = undefined;
        k_defl += 1;
        if (k_defl % non_convergence_limit == 0) {
            var s = ops.abs(T, a.get(istop - 1, istop - 2));
            if (istop > ilo + 2) s += ops.abs(T, a.get(istop - 2, istop - 3));
            a00 = ops.add(T, ops.scalar(T, dat1 * s), a.get(istop - 1, istop - 1));
            a10 = ops.scalar(T, dat2 * s);
            a01 = ops.scalar(T, s);
            a11 = a00;
        } else {
            a00 = a.get(istop - 2, istop - 2);
            a10 = a.get(istop - 1, istop - 2);
            a01 = a.get(istop - 2, istop - 1);
            a11 = a.get(istop - 1, istop - 1);
        }
        const sh = eig22(T, a00, a01, a10, a11);
        const last = a.get(istop - 1, istop - 1);
        const s1 = if (ops.abs1(T, ops.sub(T, sh[0], last)) > ops.abs1(T, ops.sub(T, sh[1], last))) sh[1] else sh[0];
        // a second small subdiagonal to start the sweep at
        var istart2 = istart;
        if (istart + 2 < istop) {
            var ii = istop - 1;
            while (ii > istart + 1) {
                ii -= 1;
                const rot = R.rotg(ops.sub(T, a.get(ii, ii), s1), a.get(ii + 1, ii))[0];
                if (ops.abs1(T, ops.mul(T, ops.conj(T, rot.s), a.get(ii, ii - 1))) <= eps * (ops.abs1(T, a.get(ii, ii - 1)) + ops.abs1(T, a.get(ii, ii + 1)))) {
                    istart2 = ii;
                    break;
                }
            }
        }
        for (istart2..istop - 1) |k| {
            var rot: R = undefined;
            if (k == istart2) {
                rot = R.rotg(ops.sub(T, a.get(k, k), s1), a.get(k + 1, k))[0];
                if (k > istart) a.at(k, k - 1).* = ops.mulReal(T, a.get(k, k - 1), rot.c);
            } else {
                const g = R.rotg(a.get(k, k - 1), a.get(k + 1, k - 1));
                rot = g[0];
                a.at(k, k - 1).* = g[1];
                a.at(k + 1, k - 1).* = ops.zero(T);
            }
            rot.left(istop_m - k, a.ptr + k + k * a.ld, a.ld, a.ptr + k + 1 + k * a.ld, a.ld);
            const top = @min(k + 3, istop);
            rot.adjoint().right(top - istart_m, a.col(k) + istart_m, a.col(k + 1) + istart_m);
            if (z) |zm| rot.adjoint().right(zm.rows, zm.col(k), zm.col(k + 1));
        }
    }
    return 0;
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const types = @import("../blas/types.zig");

test "faer complex schur: lahqr, Z T Zᴴ = H with T upper triangular" {
    const al = testing.allocator;
    inline for (.{ types.Complex64, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 3, 4, 7, 20, 64, 100 }) |n| for (0..2) |kind| {
            const h0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, h0);
            check.fill(T, h0.ptr[0 .. n * n], n + 13 * kind);
            if (kind == 1) h0.fill(ops.zero(T));
            for (0..n) |j| for (@min(j + 2, n)..n) |i| {
                h0.at(i, j).* = ops.zero(T);
            };
            const t = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, t);
            t.copyFrom(h0);
            const zm = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, zm);
            zm.identity();
            const w = try al.alloc(T, n);
            defer al.free(w);
            try testing.expectEqual(@as(usize, 0), lahqr(T, true, t, zm, w, 0, n));
            const tol = 128 * std.math.floatEps(F) * @sqrt(8 * @as(F, @floatFromInt(n))) * @max(auxiliary.normMax(T, h0), 1);
            for (0..n) |j| for (0..n) |i| {
                var s = ops.zero(T);
                for (0..n) |p| for (p..n) |q| {
                    s = ops.add(T, s, ops.mul(T, ops.mul(T, zm.get(i, p), t.get(p, q)), ops.conj(T, zm.get(j, q))));
                };
                try testing.expect(ops.abs(T, ops.sub(T, s, h0.get(i, j))) <= tol);
            };
            for (0..n) |k| {
                if (k + 1 < n) try testing.expect(ops.isZero(T, t.get(k + 1, k)));
                try testing.expectEqual(t.get(k, k), w[k]);
            }
        };
    }
}
