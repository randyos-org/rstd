//! The real Schur form's 2 x 2 kernels, `xLANV2` and `xLAEXC`, Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dlanv2.f`, `dlaexc.f`, `dlarfx.f`, OpenBLAS v0.3.34-317; modified BSD).

const std = @import("std");
const mem = @import("../../../mem.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const auxiliary = @import("auxiliary.zig");
const sylvester = @import("sylvester.zig");

const Io = std.Io;
const Error = ops.Error;

/// `DLARFX` (`dlarfx.f`), order-3 `v`: `C := H C` (`left`, `C` 3 x `n`)
/// or `C := C H` (`C` `m x 3`), `H = I - tau v vᵀ`, each sum in `v`'s order.
fn larfx(comptime F: type, left: bool, m: usize, n: usize, v: [3]F, tau: F, c: [*]F, ldc: usize) void {
    if (tau == 0) return;
    const t = [3]F{ tau * v[0], tau * v[1], tau * v[2] };
    if (left) {
        for (0..n) |j| {
            const col = c + j * ldc;
            const sum = v[0] * col[0] + v[1] * col[1] + v[2] * col[2];
            for (0..3) |i| col[i] -= sum * t[i];
        }
    } else {
        for (0..m) |j| {
            const sum = v[0] * c[j] + v[1] * c[j + ldc] + v[2] * c[j + 2 * ldc];
            for (0..3) |k| c[j + k * ldc] -= sum * t[k];
        }
    }
}

/// `MULTPL` (`dlanv2.f`).
const multpl = 4;

/// `DLANV2` (`dlanv2.f`): the standardized Schur form of the 2 x 2
/// `[a b; c d]`, in place (`c = 0`, or `a = d` and `b c < 0`), its rotation
/// `[cs sn; -sn cs]` and eigenvalues `(rt1r, rt1i)`, `(rt2r, rt2i)`.
pub fn lanv2(comptime F: type, a: *F, b: *F, c: *F, d: *F) struct { rt1r: F, rt1i: F, rt2r: F, rt2i: F, cs: F, sn: F } {
    const eps = std.math.floatEps(F);
    // `DLAMCH('B') ** INT(LOG(SAFMIN / EPS) / LOG(2) / 2)`
    const safmn2 = comptime std.math.ldexp(@as(F, 1), @divTrunc(std.math.floatExponentMin(F) + std.math.floatFractionalBits(F), 2));
    const safmx2 = 1 / safmn2;
    var cs: F = 1;
    var sn: F = 0;
    if (c.* == 0) {} else if (b.* == 0) {
        cs = 0;
        sn = 1;
        mem.swap(F, a, d);
        b.* = -c.*;
        c.* = 0;
    } else if (a.* - d.* == 0 and std.math.signbit(b.*) != std.math.signbit(c.*)) {} else {
        var temp = a.* - d.*;
        var p = 0.5 * temp;
        const bcmax = @max(@abs(b.*), @abs(c.*));
        const bcmis = @min(@abs(b.*), @abs(c.*)) * sign(F, b.*) * sign(F, c.*);
        var scale = @max(@abs(p), bcmax);
        var z = (p / scale) * p + (bcmax / scale) * bcmis;
        if (z >= multpl * eps) {
            // real eigenvalues
            z = p + std.math.copysign(@sqrt(scale) * @sqrt(z), p);
            a.* = d.* + z;
            d.* = d.* - (bcmax / z) * bcmis;
            const tau = auxiliary.lapy2(F, c.*, z);
            cs = z / tau;
            sn = c.* / tau;
            b.* = b.* - c.*;
            c.* = 0;
        } else {
            // complex or nearly equal: equalize the diagonal
            var count: usize = 0;
            var sigma = b.* + c.*;
            while (true) {
                count += 1;
                scale = @max(@abs(temp), @abs(sigma));
                if (scale >= safmx2) {
                    sigma *= safmn2;
                    temp *= safmn2;
                    if (count <= 20) continue;
                }
                if (scale <= safmn2) {
                    sigma *= safmx2;
                    temp *= safmx2;
                    if (count <= 20) continue;
                }
                break;
            }
            p = 0.5 * temp;
            var tau = auxiliary.lapy2(F, sigma, temp);
            cs = @sqrt(0.5 * (1 + @abs(sigma) / tau));
            sn = -(p / (tau * cs)) * sign(F, sigma);
            const aa = a.* * cs + b.* * sn;
            const bb = -a.* * sn + b.* * cs;
            const cc = c.* * cs + d.* * sn;
            const dd = -c.* * sn + d.* * cs;
            a.* = aa * cs + cc * sn;
            b.* = (bb * cs) + (dd * sn);
            c.* = -(aa * sn) + (cc * cs);
            d.* = -bb * sn + dd * cs;
            temp = 0.5 * (a.* + d.*);
            a.* = temp;
            d.* = temp;
            if (c.* != 0) {
                if (b.* != 0) {
                    if (sign(F, b.*) == sign(F, c.*)) {
                        // real eigenvalues after all
                        const sab = @sqrt(@abs(b.*));
                        const sac = @sqrt(@abs(c.*));
                        p = std.math.copysign(sab * sac, c.*);
                        tau = 1 / @sqrt(@abs(b.* + c.*));
                        a.* = temp + p;
                        d.* = temp - p;
                        b.* = b.* - c.*;
                        c.* = 0;
                        const cs1 = sab * tau;
                        const sn1 = sac * tau;
                        const t = cs * cs1 - sn * sn1;
                        sn = cs * sn1 + sn * cs1;
                        cs = t;
                    }
                } else {
                    b.* = -c.*;
                    c.* = 0;
                    const t = cs;
                    cs = -sn;
                    sn = t;
                }
            }
        }
    }
    const im: F = if (c.* == 0) 0 else @sqrt(@abs(b.*)) * @sqrt(@abs(c.*));
    return .{ .rt1r = a.*, .rt1i = im, .rt2r = d.*, .rt2i = -im, .cs = cs, .sn = sn };
}

/// Fortran's `SIGN(ONE, x)`.
fn sign(comptime F: type, x: F) F {
    return std.math.copysign(@as(F, 1), x);
}

/// `DLAEXC` (`dlaexc.f`): swaps the adjacent diagonal blocks of orders `n1`
/// and `n2` (1 or 2) at 0-based `j1` of the real Schur form `t`, updating
/// `q` (`n x n`) when given. Returns false where the swap would leave `t`
/// too far from Schur form (the reference's `INFO = 1`); `t` is then
/// untouched.
pub fn laexc(io: Io, comptime F: type, n: usize, t: [*]F, ldt: usize, q: ?[*]F, ldq: usize, j1: usize, n1: usize, n2: usize) Error!bool {
    if (n == 0 or n1 == 0 or n2 == 0 or j1 + n1 >= n) return true;
    const j2 = j1 + 1;
    const j3 = j1 + 2;
    const j4 = j1 + 3;
    const El = struct {
        fn at(p: [*]F, ld: usize, i: usize, j: usize) *F {
            return &p[i + j * ld];
        }
    };
    if (n1 == 1 and n2 == 1) {
        const t11 = El.at(t, ldt, j1, j1).*;
        const t22 = El.at(t, ldt, j2, j2).*;
        const g = auxiliary.lartg(F, El.at(t, ldt, j1, j2).*, t22 - t11);
        if (j3 < n) try blas_calls.rot(io, F, n - j1 - 2, t + j1 + j3 * ldt, ldt, t + j2 + j3 * ldt, ldt, g.c, g.s);
        try blas_calls.rot(io, F, j1, t + j1 * ldt, 1, t + j2 * ldt, 1, g.c, g.s);
        El.at(t, ldt, j1, j1).* = t22;
        El.at(t, ldt, j2, j2).* = t11;
        if (q) |qq| try blas_calls.rot(io, F, n, qq + j1 * ldq, 1, qq + j2 * ldq, 1, g.c, g.s);
        return true;
    }
    const ldd = 4;
    const nd = n1 + n2;
    var d: [16]F = @splat(0);
    var dnorm: F = 0;
    for (0..nd) |j| for (0..nd) |i| {
        d[i + j * ldd] = El.at(t, ldt, j1 + i, j1 + j).*;
        const v = @abs(d[i + j * ldd]);
        if (dnorm < v or std.math.isNan(v)) dnorm = v;
    };
    const eps = std.math.floatEps(F);
    const smlnum = std.math.floatMin(F) / eps;
    const limit = @max(10 * eps * dnorm, smlnum);
    const sy = sylvester.lasy2(F, false, false, -1, n1, n2, &d, ldd, d[n1 + n1 * ldd ..].ptr, ldd, d[n1 * ldd ..].ptr, ldd);
    const x11 = sy.x[0];
    const x21 = sy.x[1];
    const x12 = sy.x[2];
    const x22 = sy.x[3];
    const scale = sy.scale;
    if (n1 == 1) {
        // n1 = 1, n2 = 2
        var u = [3]F{ scale, x11, x12 };
        const tau = try auxiliary.larfg(io, F, &u[2], .of(u[0..2]));
        u[2] = 1;
        const t11 = El.at(t, ldt, j1, j1).*;
        larfx(F, true, 3, 3, u, tau, &d, ldd);
        larfx(F, false, 3, 3, u, tau, &d, ldd);
        if (@max(@abs(d[2]), @abs(d[2 + ldd]), @abs(d[2 + 2 * ldd] - t11)) > limit) return false;
        larfx(F, true, 3, n - j1, u, tau, t + j1 + j1 * ldt, ldt);
        larfx(F, false, j2 + 1, 3, u, tau, t + j1 * ldt, ldt);
        El.at(t, ldt, j3, j1).* = 0;
        El.at(t, ldt, j3, j2).* = 0;
        El.at(t, ldt, j3, j3).* = t11;
        if (q) |qq| larfx(F, false, n, 3, u, tau, qq + j1 * ldq, ldq);
    } else if (n2 == 1) {
        // n1 = 2, n2 = 1
        var u = [3]F{ -x11, -x21, scale };
        const tau = try auxiliary.larfg(io, F, &u[0], .of(u[1..]));
        u[0] = 1;
        const t33 = El.at(t, ldt, j3, j3).*;
        larfx(F, true, 3, 3, u, tau, &d, ldd);
        larfx(F, false, 3, 3, u, tau, &d, ldd);
        if (@max(@abs(d[1]), @abs(d[2]), @abs(d[0] - t33)) > limit) return false;
        larfx(F, false, j3 + 1, 3, u, tau, t + j1 * ldt, ldt);
        larfx(F, true, 3, n - j1 - 1, u, tau, t + j1 + j2 * ldt, ldt);
        El.at(t, ldt, j1, j1).* = t33;
        El.at(t, ldt, j2, j1).* = 0;
        El.at(t, ldt, j3, j1).* = 0;
        if (q) |qq| larfx(F, false, n, 3, u, tau, qq + j1 * ldq, ldq);
    } else {
        // n1 = n2 = 2
        var r1 = [3]F{ -x11, -x21, scale };
        const tau1 = try auxiliary.larfg(io, F, &r1[0], .of(r1[1..]));
        r1[0] = 1;
        const temp = -tau1 * (x12 + r1[1] * x22);
        var r2 = [3]F{ -temp * r1[1] - x22, -temp * r1[2], scale };
        const tau2 = try auxiliary.larfg(io, F, &r2[0], .of(r2[1..]));
        r2[0] = 1;
        larfx(F, true, 3, 4, r1, tau1, &d, ldd);
        larfx(F, false, 4, 3, r1, tau1, &d, ldd);
        larfx(F, true, 3, 4, r2, tau2, d[1..].ptr, ldd);
        larfx(F, false, 4, 3, r2, tau2, d[ldd..].ptr, ldd);
        if (@max(@abs(d[2]), @abs(d[2 + ldd]), @abs(d[3]), @abs(d[3 + ldd])) > limit) return false;
        larfx(F, true, 3, n - j1, r1, tau1, t + j1 + j1 * ldt, ldt);
        larfx(F, false, j4 + 1, 3, r1, tau1, t + j1 * ldt, ldt);
        larfx(F, true, 3, n - j1, r2, tau2, t + j2 + j1 * ldt, ldt);
        larfx(F, false, j4 + 1, 3, r2, tau2, t + j2 * ldt, ldt);
        El.at(t, ldt, j3, j1).* = 0;
        El.at(t, ldt, j3, j2).* = 0;
        El.at(t, ldt, j4, j1).* = 0;
        El.at(t, ldt, j4, j2).* = 0;
        if (q) |qq| {
            larfx(F, false, n, 3, r1, tau1, qq + j1 * ldq, ldq);
            larfx(F, false, n, 3, r2, tau2, qq + j2 * ldq, ldq);
        }
    }
    // standardize the 2 x 2 blocks where they now sit
    if (n2 == 2) try standardize(io, F, n, t, ldt, q, ldq, j1);
    if (n1 == 2) try standardize(io, F, n, t, ldt, q, ldq, j1 + n2);
    return true;
}

/// `dlaexc.f`'s closing `DLANV2` on the 2 x 2 block at `j`, its rotation
/// applied to the rest of `t` and to `q`.
fn standardize(io: Io, comptime F: type, n: usize, t: [*]F, ldt: usize, q: ?[*]F, ldq: usize, j: usize) Error!void {
    const k = j + 1;
    const r = lanv2(F, &t[j + j * ldt], &t[j + k * ldt], &t[k + j * ldt], &t[k + k * ldt]);
    if (j + 2 < n) try blas_calls.rot(io, F, n - j - 2, t + j + (j + 2) * ldt, ldt, t + k + (j + 2) * ldt, ldt, r.cs, r.sn);
    try blas_calls.rot(io, F, j, t + j * ldt, 1, t + k * ldt, 1, r.cs, r.sn);
    if (q) |qq| try blas_calls.rot(io, F, n, qq + j * ldq, 1, qq + k * ldq, 1, r.cs, r.sn);
}
