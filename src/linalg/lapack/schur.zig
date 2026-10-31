//! Reordering a Schur form, from Reference-LAPACK 3.12.0 as OpenBLAS
//! v0.3.34-317 ships it (`lapack-netlib/SRC/`, modified BSD): `xTREXC`
//! moves one eigenvalue (a real 2 x 2 block's pair) along the diagonal,
//! `xTRSEN` moves a selected cluster to the top with its condition numbers.
//! The real swaps are `dlaexc.f` (with `dlasy2.f`, `dlarfg.f`, `dlarfx.f`,
//! `dlanv2.f`, `dlartg.f90`), the complex ones `ztrexc.f`'s rotations
//! (`zlartg` is rstd's `zrotg`).
//!
//! faer's `schur_swap` and `schur_move` (`evd/schur/real_schur.rs`,
//! `complex_schur.rs`) are crate-internal line-for-line translations of
//! `dlaexc.f` and `dtrexc.f`, and faer has no `trsen`; libflame's are the
//! reference through f2c. So the reference is the source.
//!
//! `trsen`'s Sylvester solves are rstd's `trsyl` (ReLAPACK's recursion,
//! `relapack.zig`), its estimates `norm.lacn2`. Rotations go through rstd's
//! BLAS where it has them (`rot`, real), the rest are loops.

const std = @import("std");
const blas = @import("../blas/__root__.zig");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("faer/ops.zig");
const norm = @import("norm.zig");
const sylvester = @import("sylvester.zig");
const relapack = @import("relapack.zig");
const chol = @import("cholesky.zig");
const evd = @import("faer/evd.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Error = norm.Error;

/// `COMPQ`: update `Q` (`v`) or not (`n`).
pub const Compq = enum { n, v };

/// `trsen`'s eigenvalues, as `geev`'s: real `wr` and `wi`, complex `w`.
pub const Eigenvalues = evd.Eigenvalues;

// ---- small real kernels ------------------------------------------------------------

/// `DLARTG` (`dlartg.f90`): `[c s; -s c] [f; g] = [r; 0]`, `r` with `f`'s
/// sign.
pub fn lartg(comptime F: type, f: F, g: F) struct { c: F, s: F, r: F } {
    const safmin = std.math.floatMin(F);
    const safmax = 1 / safmin;
    const rtmin = @sqrt(safmin);
    const rtmax = @sqrt(safmax / 2);
    const f1 = @abs(f);
    const g1 = @abs(g);
    if (g == 0) return .{ .c = 1, .s = 0, .r = f };
    if (f == 0) return .{ .c = 0, .s = std.math.copysign(@as(F, 1), g), .r = g1 };
    if (f1 > rtmin and f1 < rtmax and g1 > rtmin and g1 < rtmax) {
        const d = @sqrt(f * f + g * g);
        const r = std.math.copysign(d, f);
        return .{ .c = f1 / d, .s = g / r, .r = r };
    }
    const u = @min(safmax, @max(safmin, f1, g1));
    const fs = f / u;
    const gs = g / u;
    const d = @sqrt(fs * fs + gs * gs);
    const r = std.math.copysign(d, f);
    return .{ .c = @abs(fs) / d, .s = gs / r, .r = r * u };
}

/// `DROT` through rstd's BLAS.
fn rot(io: Io, comptime F: type, n: usize, x: [*]F, incx: usize, y: [*]F, incy: usize, c: F, s: F) Error!void {
    if (n == 0) return;
    blas.rot(io, F, @intCast(n), x, @intCast(incx), y, @intCast(incy), c, s) catch |e| return norm.narrow(e);
}

/// `ZROT` (`zrot.f`): `x := c x + s y`, `y := c y - conj(s) x`, `c` real.
fn zrot(comptime T: type, n: usize, x: [*]T, incx: usize, y: [*]T, incy: usize, c: Float(T), s: T) void {
    for (0..n) |i| {
        const xi = x[i * incx];
        const yi = y[i * incy];
        x[i * incx] = ops.add(T, ops.mulReal(T, xi, c), ops.mul(T, s, yi));
        y[i * incy] = ops.sub(T, ops.mulReal(T, yi, c), ops.mul(T, ops.conj(T, s), xi));
    }
}

/// `?LARFG` (`dlarfg.f`, `zlarfg.f`): the reflector `H = I - tau v vᴴ`,
/// `v(0) = 1`, with `Hᴴ [alpha; x] = [beta; 0]`, `beta` real; `beta`
/// replaces `alpha`, `v`'s tail `x` (`n - 1`); returns `tau`. Complex
/// division is Smith's (`norm.div`) where `zlarfg.f` takes `ZLADIV`.
pub fn larfg(io: Io, comptime T: type, n: usize, alpha: *T, x: [*]T, inc: usize) Error!T {
    const F = Float(T);
    const cx = comptime isComplex(T);
    if (n == 0 or (!cx and n == 1)) return ops.zero(T);
    var xnorm = ops.nrm2(io, T, n - 1, x, inc);
    var alphr = ops.re(T, alpha.*);
    var alphi: F = if (cx) alpha.*.im else 0;
    if (xnorm == 0 and alphi == 0) return ops.zero(T);
    var beta = -std.math.copysign(if (cx) lapy3(F, alphr, alphi, xnorm) else lapy2(F, alphr, xnorm), alphr);
    // `DLAMCH('S') / DLAMCH('E')`, `E` being half of `eps`
    const safmin = std.math.floatMin(F) / (std.math.floatEps(F) / 2);
    var knt: usize = 0;
    if (@abs(beta) < safmin) {
        const rsafmn = 1 / safmin;
        while (true) {
            knt += 1;
            try norm.scal(io, T, n - 1, rsafmn, x, inc);
            beta *= rsafmn;
            alphi *= rsafmn;
            alphr *= rsafmn;
            if (!(@abs(beta) < safmin and knt < 20)) break;
        }
        xnorm = ops.nrm2(io, T, n - 1, x, inc);
        beta = -std.math.copysign(if (cx) lapy3(F, alphr, alphi, xnorm) else lapy2(F, alphr, xnorm), alphr);
    }
    var tau: T = undefined;
    if (cx) {
        tau = .init((beta - alphr) / beta, -alphi / beta);
        const s = norm.div(T, ops.one(T), .init(alphr - beta, alphi));
        try norm.scalComplex(io, T, n - 1, s, x, inc);
    } else {
        tau = (beta - alphr) / beta;
        try norm.scal(io, T, n - 1, 1 / (alphr - beta), x, inc);
    }
    for (0..knt) |_| beta *= safmin;
    alpha.* = ops.scalar(T, beta);
    return tau;
}

/// `DLAPY2` (`dlapy2.f`): `sqrt(x² + y²)` without overflow, NaN kept.
pub fn lapy2(comptime F: type, x: F, y: F) F {
    if (std.math.isNan(x)) return x;
    if (std.math.isNan(y)) return y;
    const w = @max(@abs(x), @abs(y));
    const z = @min(@abs(x), @abs(y));
    if (z == 0 or w > std.math.floatMax(F)) return w;
    return w * @sqrt(1 + (z / w) * (z / w));
}

/// `DLAPY3` (`dlapy3.f`): `sqrt(x² + y² + z²)` without overflow.
fn lapy3(comptime F: type, x: F, y: F, z: F) F {
    const xa = @abs(x);
    const ya = @abs(y);
    const za = @abs(z);
    const w = @max(xa, ya, za);
    if (w == 0 or w > std.math.floatMax(F)) return xa + ya + za;
    return w * @sqrt((xa / w) * (xa / w) + (ya / w) * (ya / w) + (za / w) * (za / w));
}

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
        std.mem.swap(F, a, d);
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
            const tau = lapy2(F, c.*, z);
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
            var tau = lapy2(F, sigma, temp);
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

/// `LOCU12`, `LOCL21`, `LOCU22`, `XSWPIV`, `BSWPIV` (`dlasy2.f`), 0-based.
const locu12 = [4]usize{ 2, 3, 0, 1 };
const locl21 = [4]usize{ 1, 0, 3, 2 };
const locu22 = [4]usize{ 3, 2, 1, 0 };
const xswpiv = [4]bool{ false, false, true, true };
const bswpiv = [4]bool{ false, true, false, true };

/// `DLASY2` (`dlasy2.f`) with `LTRANL`, `LTRANR` false and `ISGN = -1`, as
/// `dlaexc.f` calls it: `TL X - X TR = scale B` for `TL` `n1 x n1`, `TR`
/// `n2 x n2` (each 1 or 2). `x` is column-major 2 x 2 (`x[0]`, `x[1]`,
/// `x[2]`, `x[3]` are `X(1,1)`, `X(2,1)`, `X(1,2)`, `X(2,2)`); the 2 x 2 by
/// 2 x 2 case is `sylvester.lasy2`.
fn lasy2(comptime F: type, n1: usize, n2: usize, tl: [*]const F, ldtl: usize, tr: [*]const F, ldtr: usize, b: [*]const F, ldb: usize) sylvester.Small(F) {
    const sgn: F = -1;
    const eps = std.math.floatEps(F);
    const smlnum = std.math.floatMin(F) / eps;
    var out: sylvester.Small(F) = .{ .x = @splat(0), .scale = 1, .perturbed = false };
    if (n1 == 2 and n2 == 2) return sylvester.lasy2(F, false, false, sgn, tl, ldtl, tr, ldtr, .{ b[0], b[1], b[ldb], b[1 + ldb] });
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
        tmp[1] = sgn * tr[ldtr];
        tmp[2] = sgn * tr[1];
        btmp = .{ b[0], b[ldb] };
    } else {
        smin = @max(eps * @max(@abs(tr[0]), @abs(tl[0]), @abs(tl[ldtl]), @abs(tl[1]), @abs(tl[1 + ldtl])), smlnum);
        tmp[0] = tl[0] + sgn * tr[0];
        tmp[3] = tl[1 + ldtl] + sgn * tr[0];
        tmp[1] = tl[1];
        tmp[2] = tl[ldtl];
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
    if (xswpiv[ipiv]) std.mem.swap(F, &x2[0], &x2[1]);
    out.x[0] = x2[0];
    if (n1 == 1) out.x[2] = x2[1] else out.x[1] = x2[1];
    return out;
}

/// `DLAEXC` (`dlaexc.f`): swaps the adjacent diagonal blocks of orders `n1`
/// and `n2` (1 or 2) at 0-based `j1` of the real Schur form `t`, updating
/// `q` (`n x n`) when given. Returns false where the swap would leave `t`
/// too far from Schur form (the reference's `INFO = 1`); `t` is then
/// untouched.
fn laexc(io: Io, comptime F: type, n: usize, t: [*]F, ldt: usize, q: ?[*]F, ldq: usize, j1: usize, n1: usize, n2: usize) Error!bool {
    if (n == 0 or n1 == 0 or n2 == 0 or j1 + n1 >= n) return true;
    const j2 = j1 + 1;
    const j3 = j1 + 2;
    const j4 = j1 + 3;
    const T = struct {
        fn at(p: [*]F, ld: usize, i: usize, j: usize) *F {
            return &p[i + j * ld];
        }
    };
    if (n1 == 1 and n2 == 1) {
        const t11 = T.at(t, ldt, j1, j1).*;
        const t22 = T.at(t, ldt, j2, j2).*;
        const g = lartg(F, T.at(t, ldt, j1, j2).*, t22 - t11);
        if (j3 < n) try rot(io, F, n - j1 - 2, t + j1 + j3 * ldt, ldt, t + j2 + j3 * ldt, ldt, g.c, g.s);
        try rot(io, F, j1, t + j1 * ldt, 1, t + j2 * ldt, 1, g.c, g.s);
        T.at(t, ldt, j1, j1).* = t22;
        T.at(t, ldt, j2, j2).* = t11;
        if (q) |qq| try rot(io, F, n, qq + j1 * ldq, 1, qq + j2 * ldq, 1, g.c, g.s);
        return true;
    }
    const ldd = 4;
    const nd = n1 + n2;
    var d: [16]F = @splat(0);
    var dnorm: F = 0;
    for (0..nd) |j| for (0..nd) |i| {
        d[i + j * ldd] = T.at(t, ldt, j1 + i, j1 + j).*;
        const v = @abs(d[i + j * ldd]);
        if (dnorm < v or std.math.isNan(v)) dnorm = v;
    };
    const eps = std.math.floatEps(F);
    const smlnum = std.math.floatMin(F) / eps;
    const limit = @max(10 * eps * dnorm, smlnum);
    const sy = lasy2(F, n1, n2, &d, ldd, d[n1 + n1 * ldd ..].ptr, ldd, d[n1 * ldd ..].ptr, ldd);
    const x11 = sy.x[0];
    const x21 = sy.x[1];
    const x12 = sy.x[2];
    const x22 = sy.x[3];
    const scale = sy.scale;
    if (n1 == 1) {
        // n1 = 1, n2 = 2
        var u = [3]F{ scale, x11, x12 };
        const tau = try larfg(io, F, 3, &u[2], &u, 1);
        u[2] = 1;
        const t11 = T.at(t, ldt, j1, j1).*;
        larfx(F, true, 3, 3, u, tau, &d, ldd);
        larfx(F, false, 3, 3, u, tau, &d, ldd);
        if (@max(@abs(d[2]), @abs(d[2 + ldd]), @abs(d[2 + 2 * ldd] - t11)) > limit) return false;
        larfx(F, true, 3, n - j1, u, tau, t + j1 + j1 * ldt, ldt);
        larfx(F, false, j2 + 1, 3, u, tau, t + j1 * ldt, ldt);
        T.at(t, ldt, j3, j1).* = 0;
        T.at(t, ldt, j3, j2).* = 0;
        T.at(t, ldt, j3, j3).* = t11;
        if (q) |qq| larfx(F, false, n, 3, u, tau, qq + j1 * ldq, ldq);
    } else if (n2 == 1) {
        // n1 = 2, n2 = 1
        var u = [3]F{ -x11, -x21, scale };
        const tau = try larfg(io, F, 3, &u[0], u[1..].ptr, 1);
        u[0] = 1;
        const t33 = T.at(t, ldt, j3, j3).*;
        larfx(F, true, 3, 3, u, tau, &d, ldd);
        larfx(F, false, 3, 3, u, tau, &d, ldd);
        if (@max(@abs(d[1]), @abs(d[2]), @abs(d[0] - t33)) > limit) return false;
        larfx(F, false, j3 + 1, 3, u, tau, t + j1 * ldt, ldt);
        larfx(F, true, 3, n - j1 - 1, u, tau, t + j1 + j2 * ldt, ldt);
        T.at(t, ldt, j1, j1).* = t33;
        T.at(t, ldt, j2, j1).* = 0;
        T.at(t, ldt, j3, j1).* = 0;
        if (q) |qq| larfx(F, false, n, 3, u, tau, qq + j1 * ldq, ldq);
    } else {
        // n1 = n2 = 2
        var r1 = [3]F{ -x11, -x21, scale };
        const tau1 = try larfg(io, F, 3, &r1[0], r1[1..].ptr, 1);
        r1[0] = 1;
        const temp = -tau1 * (x12 + r1[1] * x22);
        var r2 = [3]F{ -temp * r1[1] - x22, -temp * r1[2], scale };
        const tau2 = try larfg(io, F, 3, &r2[0], r2[1..].ptr, 1);
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
        T.at(t, ldt, j3, j1).* = 0;
        T.at(t, ldt, j3, j2).* = 0;
        T.at(t, ldt, j4, j1).* = 0;
        T.at(t, ldt, j4, j2).* = 0;
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
    if (j + 2 < n) try rot(io, F, n - j - 2, t + j + (j + 2) * ldt, ldt, t + k + (j + 2) * ldt, ldt, r.cs, r.sn);
    try rot(io, F, j, t + j * ldt, 1, t + k * ldt, 1, r.cs, r.sn);
    if (q) |qq| try rot(io, F, n, qq + j * ldq, 1, qq + k * ldq, 1, r.cs, r.sn);
}

// ---- xTREXC ------------------------------------------------------------------------

/// What `trexc` returns: `SwapRejected` for the reference's `INFO = 1`
/// (real only): two adjacent blocks were too close to swap stably, `T`
/// part reordered and `ilst` the moving block's current row.
pub const TrexcError = error{ InvalidN, InvalidLdt, InvalidLdq, InvalidIfst, InvalidIlst, SwapRejected } || Io.Cancelable;

/// `?TREXC` (`dtrexc.f`, `ztrexc.f`): moves the diagonal block of the
/// `n x n` Schur form `T` at row `ifst` (1-based) to row `ilst` by
/// orthogonal (unitary) similarity, `T := Zᴴ T Z`, and with `compq` `v`
/// `Q := Q Z`. Real: `T` is quasi-triangular with standardized 2 x 2
/// blocks; `ifst` moves to its block's first row and `ilst` to the moved
/// block's final first row, as the reference's. Complex: both unchanged.
pub fn trexc(allocator: Allocator, io: Io, comptime T: type, compq: Compq, n: isize, t: [*]T, ldt: isize, q: [*]T, ldq: isize, ifst: *isize, ilst: *isize) TrexcError!void {
    comptime chol.floatType(T, "trexc");
    _ = allocator;
    const wantq = compq == .v;
    if (n < 0) return error.InvalidN;
    if (ldt < @max(1, n)) return error.InvalidLdt;
    if (ldq < 1 or (wantq and ldq < @max(1, n))) return error.InvalidLdq;
    if ((ifst.* < 1 or ifst.* > n) and n > 0) return error.InvalidIfst;
    if ((ilst.* < 1 or ilst.* > n) and n > 0) return error.InvalidIlst;
    if (n <= 1) return;
    if (!util.fitsMatrix(n, n, ldt, @sizeOf(T))) return error.InvalidLdt;
    if (wantq and !util.fitsMatrix(n, n, ldq, @sizeOf(T))) return error.InvalidLdq;
    const un: usize = @intCast(n);
    const qq: ?[*]T = if (wantq) q else null;
    if (comptime isComplex(T)) {
        moveComplex(T, un, t, @intCast(ldt), qq, @intCast(ldq), @intCast(ifst.*), @intCast(ilst.*));
        return;
    }
    var f: usize = @intCast(ifst.*);
    var l: usize = @intCast(ilst.*);
    defer {
        ifst.* = @intCast(f);
        ilst.* = @intCast(l);
    }
    const done = moveReal(io, T, un, t, @intCast(ldt), qq, @intCast(ldq), &f, &l) catch |e| switch (e) {
        error.OutOfMemory => unreachable,
        error.Canceled => return error.Canceled,
    };
    if (!done) return error.SwapRejected;
}

/// `ztrexc.f` after its checks (1-based `ifst`, `ilst`).
fn moveComplex(comptime T: type, n: usize, t: [*]T, ldt: usize, q: ?[*]T, ldq: usize, ifst: usize, ilst: usize) void {
    if (ifst == ilst) return;
    // 0-based swaps of k and k + 1, in the reference's order
    var k: usize = if (ifst < ilst) ifst - 1 else ifst - 2;
    while (true) {
        const t11 = t[k + k * ldt];
        const t22 = t[k + 1 + (k + 1) * ldt];
        var r = t[k + (k + 1) * ldt];
        var cs: Float(T) = undefined;
        var sn: T = undefined;
        blas.zrotg(T, &r, ops.sub(T, t22, t11), &cs, &sn);
        if (k + 2 < n) zrot(T, n - k - 2, t + k + (k + 2) * ldt, ldt, t + k + 1 + (k + 2) * ldt, ldt, cs, sn);
        zrot(T, k, t + k * ldt, 1, t + (k + 1) * ldt, 1, cs, ops.conj(T, sn));
        t[k + k * ldt] = t22;
        t[k + 1 + (k + 1) * ldt] = t11;
        if (q) |qq| zrot(T, n, qq + k * ldq, 1, qq + (k + 1) * ldq, 1, cs, ops.conj(T, sn));
        if (ifst < ilst) {
            if (k + 2 == ilst) break;
            k += 1;
        } else {
            if (k + 1 == ilst) break;
            k -= 1;
        }
    }
}

/// `dtrexc.f` after its checks: `ifst` and `ilst` (1-based) as the
/// reference leaves them; false for `INFO = 1`.
fn moveReal(io: Io, comptime F: type, n: usize, t: [*]F, ldt: usize, q: ?[*]F, ldq: usize, ifst: *usize, ilst: *usize) Error!bool {
    // `T(i, j)`, 1-based
    const at = struct {
        fn f(p: [*]F, ld: usize, i: usize, j: usize) F {
            return p[(i - 1) + (j - 1) * ld];
        }
    }.f;
    if (ifst.* > 1 and at(t, ldt, ifst.*, ifst.* - 1) != 0) ifst.* -= 1;
    var nbf: usize = 1;
    if (ifst.* < n and at(t, ldt, ifst.* + 1, ifst.*) != 0) nbf = 2;
    if (ilst.* > 1 and at(t, ldt, ilst.*, ilst.* - 1) != 0) ilst.* -= 1;
    var nbl: usize = 1;
    if (ilst.* < n and at(t, ldt, ilst.* + 1, ilst.*) != 0) nbl = 2;
    if (ifst.* == ilst.*) return true;
    const ex = struct {
        fn f(i: Io, nn: usize, tt: [*]F, ld: usize, qq: ?[*]F, lq: usize, j1: usize, n1: usize, n2: usize) Error!bool {
            return laexc(i, F, nn, tt, ld, qq, lq, j1 - 1, n1, n2);
        }
    }.f;
    var here = ifst.*;
    if (ifst.* < ilst.*) {
        // forward
        if (nbf == 2 and nbl == 1) ilst.* -= 1;
        if (nbf == 1 and nbl == 2) ilst.* += 1;
        while (true) {
            if (nbf == 1 or nbf == 2) {
                // the block is whole
                var nbnext: usize = 1;
                if (here + nbf + 1 <= n and at(t, ldt, here + nbf + 1, here + nbf) != 0) nbnext = 2;
                if (!try ex(io, n, t, ldt, q, ldq, here, nbf, nbnext)) {
                    ilst.* = here;
                    return false;
                }
                here += nbnext;
                if (nbf == 2 and at(t, ldt, here + 1, here) == 0) nbf = 3;
            } else {
                // a 2 x 2 block split in two 1 x 1s: move each
                var nbnext: usize = 1;
                if (here + 3 <= n and at(t, ldt, here + 3, here + 2) != 0) nbnext = 2;
                if (!try ex(io, n, t, ldt, q, ldq, here + 1, 1, nbnext)) {
                    ilst.* = here;
                    return false;
                }
                if (nbnext == 1) {
                    _ = try ex(io, n, t, ldt, q, ldq, here, 1, nbnext);
                    here += 1;
                } else {
                    if (at(t, ldt, here + 2, here + 1) == 0) nbnext = 1;
                    if (nbnext == 2) {
                        if (!try ex(io, n, t, ldt, q, ldq, here, 1, nbnext)) {
                            ilst.* = here;
                            return false;
                        }
                        here += 2;
                    } else {
                        _ = try ex(io, n, t, ldt, q, ldq, here, 1, 1);
                        _ = try ex(io, n, t, ldt, q, ldq, here + 1, 1, 1);
                        here += 2;
                    }
                }
            }
            if (!(here < ilst.*)) break;
        }
    } else {
        // backward
        while (true) {
            if (nbf == 1 or nbf == 2) {
                var nbnext: usize = 1;
                if (here >= 3 and at(t, ldt, here - 1, here - 2) != 0) nbnext = 2;
                if (!try ex(io, n, t, ldt, q, ldq, here - nbnext, nbnext, nbf)) {
                    ilst.* = here;
                    return false;
                }
                here -= nbnext;
                if (nbf == 2 and at(t, ldt, here + 1, here) == 0) nbf = 3;
            } else {
                var nbnext: usize = 1;
                if (here >= 3 and at(t, ldt, here - 1, here - 2) != 0) nbnext = 2;
                if (!try ex(io, n, t, ldt, q, ldq, here - nbnext, nbnext, 1)) {
                    ilst.* = here;
                    return false;
                }
                if (nbnext == 1) {
                    _ = try ex(io, n, t, ldt, q, ldq, here, nbnext, 1);
                    here -= 1;
                } else {
                    if (at(t, ldt, here, here - 1) == 0) nbnext = 1;
                    if (nbnext == 2) {
                        if (!try ex(io, n, t, ldt, q, ldq, here - 1, 2, 1)) {
                            ilst.* = here;
                            return false;
                        }
                        here -= 2;
                    } else {
                        _ = try ex(io, n, t, ldt, q, ldq, here, 1, 1);
                        _ = try ex(io, n, t, ldt, q, ldq, here - 1, 1, 1);
                        here -= 2;
                    }
                }
            }
            if (!(here > ilst.*)) break;
        }
    }
    ilst.* = here;
    return true;
}

// ---- xTRSEN ------------------------------------------------------------------------

/// `JOB`: no condition numbers (`n`), the cluster's (`e`, `s`), the
/// invariant subspace's (`v`, `sep`), or both (`b`).
pub const TrsenJob = enum { n, e, v, b };

/// What `trsen` returns: `SwapRejected` for the reference's `INFO = 1`
/// (real only): a swap was rejected, `T` part reordered, `s` and `sep`
/// zero, the eigenvalues written.
pub const TrsenError = error{ InvalidN, InvalidLdt, InvalidLdq, SwapRejected } || Allocator.Error || Io.Cancelable;

/// `?TRSEN` (`dtrsen.f`, `ztrsen.f`): reorders the `n x n` Schur form `T`
/// so the eigenvalues `select` picks (`n`; real: either of a pair picks
/// both) lead it, `T := Zᴴ T Z`, `Q := Q Z` with `compq` `v`. Writes the
/// reordered eigenvalues to `w`, their count to `m`, and as `job` asks the
/// reciprocal condition numbers of the cluster (`s`) and of the invariant
/// subspace (`sep`, an estimate).
pub fn trsen(allocator: Allocator, io: Io, comptime T: type, job: TrsenJob, compq: Compq, select: [*]const bool, n: isize, t: [*]T, ldt: isize, q: [*]T, ldq: isize, w: Eigenvalues(T), m: *isize, s: *Float(T), sep: *Float(T)) TrsenError!void {
    comptime chol.floatType(T, "trsen");
    const wants = job == .e or job == .b;
    const wantsp = job == .v or job == .b;
    const wantq = compq == .v;
    if (n < 0) return error.InvalidN;
    if (ldt < @max(1, n)) return error.InvalidLdt;
    if (ldq < 1 or (wantq and ldq < n)) return error.InvalidLdq;
    if (!util.fitsMatrix(n, n, ldt, @sizeOf(T))) return error.InvalidLdt;
    if (wantq and !util.fitsMatrix(n, n, ldq, @sizeOf(T))) return error.InvalidLdq;
    const un: usize = @intCast(n);
    const lt: usize = @intCast(ldt);
    const qq: ?[*]T = if (wantq) q else null;
    // the cluster's size, a real pair counted whole
    var cnt: usize = 0;
    {
        var k: usize = 0;
        while (k < un) : (k += 1) {
            if (!isComplex(T) and k + 1 < un and !ops.isZero(T, t[k + 1 + k * lt])) {
                if (select[k] or select[k + 1]) cnt += 2;
                k += 1;
            } else if (select[k]) cnt += 1;
        }
    }
    m.* = @intCast(cnt);
    defer eigenvalues(T, un, t, lt, w);
    if (cnt == un or cnt == 0) {
        if (wants) s.* = 1;
        if (wantsp) sep.* = try norm.matNorm(allocator, T, .o, un, un, t, lt);
        return;
    }
    // move the selected blocks up, in order
    var ks: usize = 0;
    var k: usize = 0;
    while (k < un) : (k += 1) {
        const pair = !isComplex(T) and k + 1 < un and !ops.isZero(T, t[k + 1 + k * lt]);
        const swap = select[k] or (pair and select[k + 1]);
        if (swap) {
            ks += 1;
            if (k + 1 != ks) {
                if (comptime isComplex(T)) {
                    moveComplex(T, un, t, lt, qq, @intCast(ldq), k + 1, ks);
                } else {
                    var kk = k + 1;
                    if (!try moveReal(io, T, un, t, lt, qq, @intCast(ldq), &kk, &ks)) {
                        if (wants) s.* = 0;
                        if (wantsp) sep.* = 0;
                        return error.SwapRejected;
                    }
                }
            }
            if (pair) ks += 1;
        }
        if (pair) k += 1;
    }
    const n1 = cnt;
    const n2 = un - cnt;
    const nn = n1 * n2;
    if (!wants and !wantsp) return;
    const work = try allocator.alloc(T, 2 * nn);
    defer allocator.free(work);
    const ctx: Sep(T) = .{ .al = allocator, .io = io, .n1 = n1, .n2 = n2, .t = t, .ldt = lt, .x = work.ptr };
    if (wants) {
        // s from the solution of T11 X - X T22 = scale T12
        for (0..n2) |j| @memcpy(work[j * n1 ..][0..n1], t[(n1 + j) * lt ..][0..n1]);
        const scale = try ctx.solve(false, work.ptr);
        const rnorm = try norm.matNorm(allocator, T, .f, n1, n2, work.ptr, n1);
        s.* = if (rnorm == 0) 1 else scale / (@sqrt(scale * scale / rnorm + rnorm) * @sqrt(rnorm));
    }
    if (wantsp) {
        var last = Sep(T).Last{ .ctx = ctx };
        const est = try norm.lacn2(allocator, T, Error, nn, work.ptr + nn, work.ptr, &last, Sep(T).Last.apply);
        sep.* = last.scale / est;
    }
}

/// `trsen`'s Sylvester operator `X ↦ T11 X - X T22` on `n1 x n2` `X`.
fn Sep(comptime T: type) type {
    return struct {
        al: Allocator,
        io: Io,
        n1: usize,
        n2: usize,
        t: [*]T,
        ldt: usize,
        x: [*]T,

        const Self = @This();

        /// Solves for `X` in `x` (the operator's inverse, or its adjoint's),
        /// returning `scale`; a perturbed solve counts, as the reference
        /// ignores its `IERR`.
        fn solve(s: Self, adjoint: bool, x: [*]T) Error!Float(T) {
            var scale: Float(T) = 1;
            const tr: types.Trans = if (adjoint) (if (comptime isComplex(T)) .c else .t) else .n;
            const n1: isize = @intCast(s.n1);
            const n2: isize = @intCast(s.n2);
            const ld: isize = @intCast(s.ldt);
            relapack.trsyl(s.al, s.io, T, tr, tr, -1, n1, n2, s.t, ld, s.t + s.n1 + s.n1 * s.ldt, ld, x, n1, &scale) catch |e| switch (e) {
                error.Perturbed => {},
                error.OutOfMemory => return error.OutOfMemory,
                error.Canceled => return error.Canceled,
                else => unreachable,
            };
            return scale;
        }

        /// `lacn2`'s callback, keeping the last solve's `scale`.
        const Last = struct {
            ctx: Self,
            scale: Float(T) = 1,

            fn apply(l: *Last, adjoint: bool, x: [*]T) Error!void {
                l.scale = try l.ctx.solve(adjoint, x);
            }
        };
    };
}

/// `T`'s eigenvalues into `w`: the diagonal, and real pairs `a ± i b` with
/// `b = sqrt(|T(k, k+1)|) sqrt(|T(k+1, k)|)`.
fn eigenvalues(comptime T: type, n: usize, t: [*]const T, ldt: usize, w: Eigenvalues(T)) void {
    if (comptime isComplex(T)) {
        for (0..n) |k| w[k] = t[k + k * ldt];
        return;
    }
    for (0..n) |k| {
        w.re[k] = t[k + k * ldt];
        w.im[k] = 0;
    }
    for (0..n -| 1) |k| if (t[k + 1 + k * ldt] != 0) {
        w.im[k] = @sqrt(@abs(t[k + (k + 1) * ldt])) * @sqrt(@abs(t[k + 1 + k * ldt]));
        w.im[k + 1] = -w.im[k];
    };
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/blas1/check.zig");
const lu = @import("lu.zig");
const svd = @import("faer/svd.zig");

const test_types = chol.test_types;

/// `THRESH` of the reference's eigenvalue-condition tests
/// (`TESTING/dec.in`, `dchkec.f`'s `dget36` and `dget38`).
const thresh = 20;

/// An `n x n` Schur form: upper triangular with diagonal `1, 2, ...`
/// (complex: plus `i` times as much), real ones with a standardized 2 x 2
/// block `[k b; c k]`, `b c < 0`, wherever `pairs` has a bit set and room
/// remains; random above.
fn schurForm(comptime T: type, al: Allocator, n: usize, pairs: u32, seed: u64) ![]T {
    const F = Float(T);
    const t = try al.alloc(T, n * n);
    check.fill(T, t, seed);
    for (0..n) |j| for (j..n) |i| {
        t[i + j * n] = ops.zero(T);
    };
    var k: usize = 0;
    while (k < n) {
        const v: F = @floatFromInt(k + 1);
        if (!isComplex(T) and k + 1 < n and (pairs >> @intCast(k % 32)) & 1 != 0) {
            t[k + k * n] = ops.scalar(T, v);
            t[k + 1 + (k + 1) * n] = ops.scalar(T, v);
            t[k + (k + 1) * n] = ops.scalar(T, 1.5);
            t[k + 1 + k * n] = ops.scalar(T, -0.75);
            k += 2;
        } else {
            t[k + k * n] = if (comptime isComplex(T)) .init(v, v / 2) else v;
            k += 1;
        }
    }
    return t;
}

/// `Q`'s orthonormality and `Q T Qᴴ = T0` within the threshold of `n eps`
/// (`dget36.f`'s two residuals), and `T` still in Schur form.
fn expectSimilar(comptime T: type, n: usize, t0: []const T, t: []const T, q: []const T) !void {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const nf: F = @floatFromInt(@max(n, 1));
    var tnorm: F = 0;
    for (t0) |v| tnorm = @max(tnorm, ops.abs(T, v));
    for (0..n) |i| for (0..n) |j| {
        var g = ops.zero(T);
        var r = ops.zero(T);
        for (0..n) |k| {
            g = ops.add(T, g, ops.mul(T, ops.conj(T, q[k + i * n]), q[k + j * n]));
            var s = ops.zero(T);
            for (0..n) |l| s = ops.add(T, s, ops.mul(T, t[k + l * n], ops.conj(T, q[j + l * n])));
            r = ops.add(T, r, ops.mul(T, q[i + k * n], s));
        }
        if (i == j) g = ops.sub(T, g, ops.one(T));
        try testing.expect(ops.abs(T, g) <= thresh * nf * eps);
        try testing.expect(ops.abs(T, ops.sub(T, r, t0[i + j * n])) <= thresh * nf * eps * tnorm);
    };
    for (0..n) |j| for (j + 1..n) |i| {
        const v = t[i + j * n];
        if (comptime isComplex(T)) {
            try testing.expect(ops.isZero(T, v));
        } else if (i > j + 1) {
            try testing.expect(ops.isZero(T, v));
        } else if (!ops.isZero(T, v)) {
            // a standardized block, alone on the subdiagonal
            try testing.expect(t[j + j * n] == t[i + i * n]);
            try testing.expect(std.math.signbit(ops.re(T, v)) != std.math.signbit(ops.re(T, t[j + i * n])));
            if (j > 0) try testing.expect(ops.isZero(T, t[j + (j - 1) * n]));
        }
    };
}

/// The eigenvalue (real: of the block) starting at 0-based `k`.
fn eigAt(comptime T: type, n: usize, t: []const T, k: usize) [2]f64 {
    const v = t[k + k * n];
    if (comptime isComplex(T)) return .{ @floatCast(v.re), @floatCast(v.im) };
    if (k + 1 < n and t[k + 1 + k * n] != 0) return .{ @floatCast(v), @floatCast(@sqrt(@abs(t[k + (k + 1) * n])) * @sqrt(@abs(t[k + 1 + k * n]))) };
    return .{ @floatCast(v), 0 };
}

test "lapack trexc: moves every block to every place, Q T Q' = T0, Schur form kept" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 3, 5, 8 }) |n| for ([_]u32{ 0, 0b1, 0b10, 0b1001, 0b100101 }) |pairs| {
            const t0 = try schurForm(T, al, n, pairs, n * 3 + pairs);
            defer al.free(t0);
            const t = try al.alloc(T, n * n);
            defer al.free(t);
            const q = try al.alloc(T, n * n);
            defer al.free(q);
            for (1..n + 1) |f| for (1..n + 1) |l| {
                @memcpy(t, t0);
                @memset(q, ops.zero(T));
                for (0..n) |i| q[i * (n + 1)] = ops.one(T);
                var ifst: isize = @intCast(f);
                var ilst: isize = @intCast(l);
                const ni: isize = @intCast(n);
                try trexc(al, check.io, T, .v, ni, t.ptr, ni, q.ptr, ni, &ifst, &ilst);
                try expectSimilar(T, n, t0, t, q);
                // the moved eigenvalue sits at ilst
                const want = eigAt(T, n, t0, @intCast(ifst - 1));
                const got = eigAt(T, n, t, @intCast(ilst - 1));
                const tol = 100 * @as(f64, @floatCast(std.math.floatEps(F))) * @as(f64, @floatFromInt(n));
                try testing.expectApproxEqAbs(want[0], got[0], tol);
                try testing.expectApproxEqAbs(@abs(want[1]), @abs(got[1]), tol);
            };
        };
    }
}

test "lapack trexc: n = 0 and 1, compq n, bad arguments" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        var t = [_]T{ops.scalar(T, 3)};
        var q = [_]T{ops.scalar(T, 7)};
        var f: isize = 1;
        var l: isize = 1;
        try trexc(al, check.io, T, .v, 1, &t, 1, &q, 1, &f, &l);
        try trexc(al, check.io, T, .n, 0, &t, 1, &q, 1, &f, &l);
        try testing.expect(ops.re(T, t[0]) == 3 and ops.re(T, q[0]) == 7);
    }
    var t: [4]f64 = .{ 1, 0, 2, 3 };
    var q: [4]f64 = undefined;
    var f: isize = 1;
    var l: isize = 2;
    try trexc(al, check.io, f64, .n, 2, &t, 2, &q, 1, &f, &l);
    try testing.expect(t[0] == 3 and t[3] == 1 and @abs(t[2]) == 2 and t[1] == 0);
    try testing.expectError(error.InvalidN, trexc(al, check.io, f64, .n, -1, &t, 1, &q, 1, &f, &l));
    try testing.expectError(error.InvalidLdt, trexc(al, check.io, f64, .n, 2, &t, 1, &q, 1, &f, &l));
    try testing.expectError(error.InvalidLdq, trexc(al, check.io, f64, .v, 2, &t, 2, &q, 1, &f, &l));
    f = 3;
    try testing.expectError(error.InvalidIfst, trexc(al, check.io, f64, .n, 2, &t, 2, &q, 1, &f, &l));
    f = 1;
    l = 0;
    try testing.expectError(error.InvalidIlst, trexc(al, check.io, f64, .n, 2, &t, 2, &q, 1, &f, &l));
}

/// `I ⊗ T11 - T22ᵀ ⊗ I`, the Sylvester operator `trsen` measures, on
/// column-major `vec(X)`.
fn kron(comptime T: type, al: Allocator, n: usize, t: []const T, n1: usize) ![]T {
    const n2 = n - n1;
    const nn = n1 * n2;
    const k = try al.alloc(T, nn * nn);
    @memset(k, ops.zero(T));
    for (0..n2) |j| for (0..n1) |i| {
        const row = i + j * n1;
        for (0..n1) |p| k[row + (p + j * n1) * nn] = ops.add(T, k[row + (p + j * n1) * nn], t[i + p * n]);
        for (0..n2) |r| k[row + (i + r * n1) * nn] = ops.sub(T, k[row + (i + r * n1) * nn], t[n1 + r + (n1 + j) * n]);
    };
    return k;
}

test "lapack trsen: the selected cluster leads, s against the Kronecker solve, sep against its smallest singular value" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const eps = std.math.floatEps(F);
        for ([_]usize{ 2, 4, 7 }) |n| for ([_]u32{ 0, 0b101 }) |pairs| for ([_]u32{ 0b1, 0b10, 0b110, 0b1010010 }) |sel| {
            const t0 = try schurForm(T, al, n, pairs, n * 5 + sel);
            defer al.free(t0);
            const t = try al.alloc(T, n * n);
            defer al.free(t);
            @memcpy(t, t0);
            const q = try al.alloc(T, n * n);
            defer al.free(q);
            @memset(q, ops.zero(T));
            for (0..n) |i| q[i * (n + 1)] = ops.one(T);
            var select: [8]bool = undefined;
            for (0..n) |i| select[i] = (sel >> @intCast(i)) & 1 != 0;
            var wbuf: [16]T = undefined;
            const w: Eigenvalues(T) = if (comptime isComplex(T)) &wbuf else .{ .re = &wbuf, .im = wbuf[8..].ptr };
            var m: isize = -1;
            var s: F = -1;
            var sep: F = -1;
            const ni: isize = @intCast(n);
            try trsen(al, check.io, T, .b, .v, &select, ni, t.ptr, ni, q.ptr, ni, w, &m, &s, &sep);
            try expectSimilar(T, n, t0, t, q);
            const n1: usize = @intCast(m);
            // the selected eigenvalues, in their old order, lead
            var k: usize = 0;
            var got: usize = 0;
            while (k < n) {
                const e = eigAt(T, n, t0, k);
                const size: usize = if (!isComplex(T) and e[1] != 0) 2 else 1;
                if (select[k] or (size == 2 and select[k + 1])) {
                    const g = eigAt(T, n, t, got);
                    try testing.expectApproxEqAbs(e[0], g[0], 1e3 * @as(f64, @floatCast(eps)));
                    got += size;
                }
                k += size;
            }
            try testing.expectEqual(got, n1);
            if (n1 == 0 or n1 == n) {
                try testing.expect(s == 1);
                continue;
            }
            const n2 = n - n1;
            const nn = n1 * n2;
            const kk = try kron(T, al, n, t, n1);
            defer al.free(kk);
            // s = 1 / sqrt(1 + |X|_F^2) for K vec(X) = vec(T12)
            const x = try al.alloc(T, nn);
            defer al.free(x);
            for (0..n2) |j| for (0..n1) |i| {
                x[i + j * n1] = t[i + (n1 + j) * n];
            };
            const kf = try al.alloc(T, nn * nn);
            defer al.free(kf);
            @memcpy(kf, kk);
            const ipiv = try al.alloc(isize, nn);
            defer al.free(ipiv);
            const nni: isize = @intCast(nn);
            try lu.gesv(al, check.io, T, nni, 1, kf.ptr, nni, ipiv.ptr, x.ptr, nni, null);
            var xf: F = 0;
            for (x) |v| xf += ops.abs2(T, v);
            const s_true = 1 / @sqrt(1 + xf);
            try testing.expectApproxEqRel(s_true, s, 1e4 * eps);
            // sep within sqrt(nn) and the threshold of sigma_min(K)
            @memcpy(kf, kk);
            const sv = try al.alloc(F, nn);
            defer al.free(sv);
            var dummy: [1]T = undefined;
            try svd.gesvd(al, check.io, T, .n, .n, nni, nni, kf.ptr, nni, sv.ptr, &dummy, 1, &dummy, 1);
            const smin = sv[nn - 1];
            const root = @sqrt(@as(F, @floatFromInt(nn)));
            if (!(sep >= smin / root * (1 - 1e3 * eps) and sep <= thresh * root * smin)) {
                check.say("trsen {s} n {d} select {b}: sep {e}, sigma_min {e}\n", .{ @typeName(T), n, sel, @as(f64, @floatCast(sep)), @as(f64, @floatCast(smin)) });
                return error.TestUnexpectedResult;
            }
        };
    }
}

test "lapack trsen: none or all selected, a block-diagonal T, n = 0, bad arguments" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        // diag(1, 4): s = 1 (T12 = 0), sep = |1 - 4| exactly
        var t = [_]T{ ops.scalar(T, 1), ops.zero(T), ops.zero(T), ops.scalar(T, 4) };
        var q: [4]T = undefined;
        var wbuf: [4]T = undefined;
        const w: Eigenvalues(T) = if (comptime isComplex(T)) &wbuf else .{ .re = &wbuf, .im = wbuf[2..].ptr };
        var m: isize = -1;
        var s: F = -1;
        var sep: F = -1;
        try trsen(al, check.io, T, .b, .n, &[_]bool{ false, true }, 2, &t, 2, &q, 1, w, &m, &s, &sep);
        try testing.expect(m == 1 and s == 1 and sep == 3);
        try testing.expect(ops.re(T, t[0]) == 4 and ops.re(T, t[3]) == 1);
        try trsen(al, check.io, T, .b, .n, &[_]bool{ false, false }, 2, &t, 2, &q, 1, w, &m, &s, &sep);
        try testing.expect(m == 0 and s == 1 and sep == 4);
        try trsen(al, check.io, T, .v, .n, &[_]bool{ true, true }, 2, &t, 2, &q, 1, w, &m, &s, &sep);
        try testing.expect(m == 2 and sep == 4);
        try trsen(al, check.io, T, .b, .n, &[_]bool{}, 0, &t, 1, &q, 1, w, &m, &s, &sep);
        try testing.expect(m == 0 and s == 1);
    }
    var t: [4]f64 = undefined;
    var q: [4]f64 = undefined;
    var wr: [2]f64 = undefined;
    var wi: [2]f64 = undefined;
    var m: isize = 0;
    var s: f64 = 0;
    var sep: f64 = 0;
    const sel = [_]bool{ true, false };
    try testing.expectError(error.InvalidN, trsen(al, check.io, f64, .n, .n, &sel, -1, &t, 1, &q, 1, .{ .re = &wr, .im = &wi }, &m, &s, &sep));
    try testing.expectError(error.InvalidLdt, trsen(al, check.io, f64, .n, .n, &sel, 2, &t, 1, &q, 1, .{ .re = &wr, .im = &wi }, &m, &s, &sep));
    try testing.expectError(error.InvalidLdq, trsen(al, check.io, f64, .n, .v, &sel, 2, &t, 2, &q, 1, .{ .re = &wr, .im = &wi }, &m, &s, &sep));
}
