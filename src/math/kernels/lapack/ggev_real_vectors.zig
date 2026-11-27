//! The eigenvectors of a real generalized Schur form, faer v0.24.4's
//! `qz_to_gevd_real`, `solve_complex_shifted_1x1` and
//! `solve_complex_shifted_2x2` (`faer/src/linalg/gevd/mod.rs`, MIT). Two
//! places follow Reference-LAPACK (`dtgevc.f`, `dlaln2.f`; OpenBLAS
//! v0.3.34-317) instead (issue #766): a divisor below `smin` is replaced by
//! `smin` alone (`xLALN2`'s `SMINI` identity; faer divides by both), with
//! `xLALN2`'s floor on the 2 x 2 pivot, and a pair's overflow guard takes `B`'s
//! column norm (`WORK(N+J+1)`), where faer repeats `A`'s.

const std = @import("std");
const qz_2x2 = @import("qz_2x2.zig");

const Allocator = std.mem.Allocator;
const C = std.math.Complex;

/// Left (`u_left`, through `Q`) and right (`u_right`, through `Z`)
/// eigenvectors of the real Schur pencil `(S, P)`: one column per real
/// eigenvalue, a pair's real and imaginary parts in two columns, each vector
/// of unit 2-norm.
pub fn eigenvectors(allocator: Allocator, comptime T: type, n: usize, a: []const T, lda: usize, b: []const T, ldb: usize, u_left: ?[]T, ldl: usize, u_right: ?[]T, ldr: usize) Allocator.Error!void {
    const work = try allocator.alloc(T, 2 * n + 4 * n);
    defer allocator.free(work);
    const s: State(T) = .{
        .n = n,
        .a = a,
        .lda = lda,
        .b = b,
        .ldb = ldb,
        .acolnorm = work[0..n],
        .bcolnorm = work[n .. 2 * n],
        .rhs = work[2 * n .. 4 * n],
        .tmp = work[4 * n ..],
    };
    s.norms();
    if (u_left) |u| s.left(u, ldl);
    if (u_right) |u| s.right(u, ldr);
}

fn State(comptime T: type) type {
    return struct {
        n: usize,
        a: []const T,
        lda: usize,
        b: []const T,
        ldb: usize,
        acolnorm: []T,
        bcolnorm: []T,
        /// `n x 2`: a vector's real and imaginary parts
        rhs: []T,
        tmp: []T,

        const Self = @This();
        const ulp = std.math.floatEps(T);
        const safmin = std.math.floatMin(T);

        fn av(s: Self, i: usize, j: usize) T {
            return s.a[i + j * s.lda];
        }

        fn bv(s: Self, i: usize, j: usize) T {
            return s.b[i + j * s.ldb];
        }

        fn r(s: Self, i: usize, jw: usize) *T {
            return &s.rhs[i + jw * s.n];
        }

        fn norms(s: Self) void {
            var j: usize = 0;
            while (j < s.n) {
                const pair = j + 1 < s.n and s.av(j + 1, j) != 0;
                for (j..j + @as(usize, if (pair) 2 else 1)) |jc| {
                    var sa: T = 0;
                    for (0..j) |i| sa += @abs(s.av(i, jc));
                    s.acolnorm[jc] = sa;
                }
                j += if (pair) 2 else 1;
            }
            for (0..s.n) |jj| {
                var sb: T = 0;
                for (0..jj) |i| sb += @abs(s.bv(i, jj));
                s.bcolnorm[jj] = sb;
            }
        }

        fn anorm(s: Self) T {
            var r_: T = 0;
            var j: usize = 0;
            while (j < s.n) {
                const pair = j + 1 < s.n and s.av(j + 1, j) != 0;
                if (!pair) {
                    r_ = @max(r_, s.acolnorm[j] + @abs(s.av(j, j)));
                    j += 1;
                } else {
                    for (j..j + 2) |jc| r_ = @max(r_, s.acolnorm[jc] + @abs(s.av(j, jc)) + @abs(s.av(j + 1, jc)));
                    j += 2;
                }
            }
            return r_;
        }

        fn bnorm(s: Self) T {
            var r_: T = 0;
            for (0..s.n) |j| r_ = @max(r_, s.bcolnorm[j] + @abs(s.bv(j, j)));
            return r_;
        }

        const Coefs = struct { acoef: T, bcoefr: T, bcoefi: T };

        /// A real eigenvalue's scaled `(beta, alpha)`, kept from overflow.
        fn realCoefs(s: Self, je: usize, as: T, bs: T, an: T, bn: T) Coefs {
            const small = safmin * @as(T, @floatFromInt(s.n)) / ulp;
            const big = 1 / small;
            const temp = @max(safmin, @max(@abs(s.av(je, je)) * as, @abs(s.bv(je, je)) * bs));
            const salfar = s.av(je, je) * temp * as;
            const sbeta = s.bv(je, je) * temp * bs;
            var acoef = sbeta * as;
            var bcoefr = salfar * bs;
            var scale: T = 1;
            const lsa = @abs(sbeta) >= safmin and @abs(acoef) < small;
            const lsb = @abs(salfar) >= safmin and @abs(bcoefr) < small;
            if (lsa) scale = (small / @abs(sbeta)) * @min(an, big);
            if (lsb) scale = @max(scale, (small / @abs(salfar)) * @min(bn, big));
            if (lsa or lsb) {
                scale = @min(scale, 1 / (safmin * @max(1, @max(@abs(acoef), @abs(bcoefr)))));
                acoef = if (lsa) as * (scale * sbeta) else scale * acoef;
                bcoefr = if (lsb) bs * (scale * salfar) else scale * bcoefr;
            }
            return .{ .acoef = acoef, .bcoefr = bcoefr, .bcoefi = 0 };
        }

        /// A complex pair's `(scale, wr, wi)` from its 2 x 2 block, scaled.
        fn pairCoefs(s: Self, je: usize, as: T, bs: T) Coefs {
            const e = qz_2x2.eig2x2(T, .{ s.av(je, je), s.av(je, je + 1), s.av(je + 1, je), s.av(je + 1, je + 1) }, .{ s.bv(je, je), s.bv(je, je + 1), s.bv(je + 1, je), s.bv(je + 1, je + 1) });
            var acoef = e.scale1;
            var bcoefr = e.wr1;
            var bcoefi = e.wi;
            const acoefa = @abs(acoef);
            const bcoefa = @abs(bcoefr) + @abs(bcoefi);
            var scale: T = 1;
            if (acoefa * ulp < safmin and acoefa >= safmin) scale = (safmin / ulp) / acoefa;
            if (bcoefa * ulp < safmin and bcoefa >= safmin) scale = @max(scale, (safmin / ulp) / bcoefa);
            if (safmin * acoefa > as) scale = as / (safmin * acoefa);
            if (safmin * bcoefa > bs) scale = @min(scale, bs / (safmin * bcoefa));
            if (scale != 1) {
                acoef *= scale;
                bcoefr *= scale;
                bcoefi *= scale;
            }
            return .{ .acoef = acoef, .bcoefr = bcoefr, .bcoefi = bcoefi };
        }

        fn unit(u: []T, ld: usize, n: usize, je: usize) void {
            for (0..n) |i| u[i + je * ld] = if (i == je) 1 else 0;
        }

        /// `u(:, je..je + nw) = u(:, cols) rhs(cols, 0..nw)`, unit 2-norm.
        fn backTransform(s: Self, u: []T, ld: usize, je: usize, nw: usize, from: usize, to: usize) void {
            const n = s.n;
            var ss: T = 0;
            for (0..nw) |jw| for (0..n) |i| {
                var acc: T = 0;
                for (from..to) |k| acc += u[i + k * ld] * s.r(k, jw).*;
                s.tmp[i + jw * n] = acc;
                ss = std.math.hypot(ss, acc);
            };
            const inv = 1 / ss;
            for (0..nw) |jw| for (0..n) |i| {
                u[i + (je + jw) * ld] = s.tmp[i + jw * n] * inv;
            };
        }

        fn left(s: Self, u: []T, ld: usize) void {
            const n = s.n;
            const an = s.anorm();
            const bn = s.bnorm();
            const as = 1 / @max(an, safmin);
            const bs = 1 / @max(bn, safmin);
            const bignum = 1 / (safmin * @as(T, @floatFromInt(n)));
            var je: usize = 0;
            while (je < n) {
                const pair = je + 1 < n and s.av(je + 1, je) != 0;
                const nw: usize = if (pair) 2 else 1;
                if (!pair and @max(@abs(s.av(je, je)), @abs(s.bv(je, je))) < safmin) {
                    unit(u, ld, n, je);
                    je += 1;
                    continue;
                }
                @memset(s.rhs, 0);
                var c: Coefs = undefined;
                var xmax: T = 1;
                if (!pair) {
                    c = s.realCoefs(je, as, bs, an, bn);
                    s.r(je, 0).* = 1;
                } else {
                    c = s.pairCoefs(je, as, bs);
                    const temp = c.acoef * s.av(je + 1, je);
                    const temp2r = c.acoef * s.av(je, je) - c.bcoefr * s.bv(je, je);
                    const temp2i = -c.bcoefi * s.bv(je, je);
                    if (@abs(temp) > @abs(temp2r) + @abs(temp2i)) {
                        s.r(je, 0).* = 1;
                        s.r(je + 1, 0).* = -temp2r / temp;
                        s.r(je + 1, 1).* = -temp2i / temp;
                    } else {
                        s.r(je + 1, 0).* = 1;
                        const t2 = c.acoef * s.av(je, je + 1);
                        s.r(je, 0).* = (c.bcoefr * s.bv(je + 1, je + 1) - c.acoef * s.av(je + 1, je + 1)) / t2;
                        s.r(je, 1).* = c.bcoefi * s.bv(je + 1, je + 1) / t2;
                    }
                    xmax = @max(@abs(s.r(je, 0).*) + @abs(s.r(je, 1).*), @abs(s.r(je + 1, 0).*) + @abs(s.r(je + 1, 1).*));
                }
                const acoefa = @abs(c.acoef);
                const bcoefa = @abs(c.bcoefr) + @abs(c.bcoefi);
                const dmin = @max(safmin, @max(ulp * acoefa * an, ulp * bcoefa * bn));
                var j = je + nw;
                while (j < n) {
                    const cplx = j + 1 < n and s.av(j + 1, j) != 0;
                    const na: usize = if (cplx) 2 else 1;
                    const xscale = 1 / xmax;
                    var temp = @max(@max(s.acolnorm[j], s.bcolnorm[j]), acoefa * s.acolnorm[j] + bcoefa * s.bcolnorm[j]);
                    // `dtgevc.f`'s `WORK(N+J+1)`, `B`'s norm (faer repeats `A`'s)
                    if (cplx) temp = @max(@max(temp, acoefa * s.acolnorm[j + 1] + bcoefa * s.bcolnorm[j + 1]), @max(s.acolnorm[j + 1], s.bcolnorm[j + 1]));
                    if (temp > bignum * xscale) {
                        for (0..nw) |jw| for (je..j) |jr| {
                            s.r(jr, jw).* *= xscale;
                        };
                        xmax *= xscale;
                    }
                    var sums: [2][2]T = .{ .{ 0, 0 }, .{ 0, 0 } };
                    var sump: [2][2]T = .{ .{ 0, 0 }, .{ 0, 0 } };
                    for (0..nw) |jw| for (0..na) |ja| {
                        for (je..j) |jr| {
                            sums[ja][jw] += s.av(jr, j + ja) * s.r(jr, jw).*;
                            sump[ja][jw] += s.bv(jr, j + ja) * s.r(jr, jw).*;
                        }
                    };
                    for (0..na) |ja| {
                        if (pair) {
                            s.r(j + ja, 0).* = -c.acoef * sums[ja][0] + c.bcoefr * sump[ja][0] - c.bcoefi * sump[ja][1];
                            s.r(j + ja, 1).* = -c.acoef * sums[ja][1] + c.bcoefr * sump[ja][1] + c.bcoefi * sump[ja][0];
                        } else s.r(j + ja, 0).* = -c.acoef * sums[ja][0] + c.bcoefr * sump[ja][0];
                    }
                    // `Sᵀ`'s diagonal block: the left vectors solve the transpose
                    if (cplx) {
                        solve2x2(T, dmin, c.acoef, .{ s.av(j, j), s.av(j + 1, j), s.av(j, j + 1), s.av(j + 1, j + 1) }, s.bv(j, j), s.bv(j + 1, j + 1), s, j, nw, c.bcoefr, c.bcoefi);
                    } else solve1x1(T, dmin, c.acoef, s.av(j, j), s.bv(j, j), s, j, nw, c.bcoefr, c.bcoefi);
                    j += na;
                }
                s.backTransform(u, ld, je, nw, je, n);
                // LAPACK's `VL`: `u = v(:, je) + i v(:, je + 1)` with `uᴴ A = λ uᴴ B`;
                // faer's pair is its conjugate
                if (pair) for (0..n) |i| {
                    u[i + (je + 1) * ld] = -u[i + (je + 1) * ld];
                };
                je += nw;
            }
        }

        fn right(s: Self, u: []T, ld: usize) void {
            const n = s.n;
            const an = s.anorm();
            const bn = s.bnorm();
            const as = 1 / @max(an, safmin);
            const bs = 1 / @max(bn, safmin);
            var je = n;
            while (je > 0) {
                je -= 1;
                const pair = je > 0 and s.av(je, je - 1) != 0;
                const nw: usize = if (pair) 2 else 1;
                je -= nw - 1;
                if (!pair and @max(@abs(s.av(je, je)), @abs(s.bv(je, je))) < safmin) {
                    unit(u, ld, n, je);
                    continue;
                }
                @memset(s.rhs, 0);
                var c: Coefs = undefined;
                if (!pair) {
                    c = s.realCoefs(je, as, bs, an, bn);
                    s.r(je, 0).* = 1;
                    for (0..je) |jr| s.r(jr, 0).* = c.bcoefr * s.bv(jr, je) - c.acoef * s.av(jr, je);
                } else {
                    c = s.pairCoefs(je, as, bs);
                    const temp = c.acoef * s.av(je + 1, je);
                    const temp2r = c.acoef * s.av(je + 1, je + 1) - c.bcoefr * s.bv(je + 1, je + 1);
                    const temp2i = -c.bcoefi * s.bv(je + 1, je + 1);
                    if (@abs(temp) > @abs(temp2r) + @abs(temp2i)) {
                        s.r(je + 1, 0).* = 1;
                        s.r(je, 0).* = -temp2r / temp;
                        s.r(je, 1).* = -temp2i / temp;
                    } else {
                        s.r(je, 0).* = 1;
                        const t2 = c.acoef * s.av(je, je + 1);
                        s.r(je + 1, 0).* = (c.bcoefr * s.bv(je, je) - c.acoef * s.av(je, je)) / t2;
                        s.r(je + 1, 1).* = c.bcoefi * s.bv(je, je) / t2;
                    }
                    const creala = c.acoef * s.r(je, 0).*;
                    const cimaga = c.acoef * s.r(je, 1).*;
                    const crealb = c.bcoefr * s.r(je, 0).* - c.bcoefi * s.r(je, 1).*;
                    const cimagb = c.bcoefi * s.r(je, 0).* + c.bcoefr * s.r(je, 1).*;
                    const cre2a = c.acoef * s.r(je + 1, 0).*;
                    const cim2a = c.acoef * s.r(je + 1, 1).*;
                    const cre2b = c.bcoefr * s.r(je + 1, 0).* - c.bcoefi * s.r(je + 1, 1).*;
                    const cim2b = c.bcoefi * s.r(je + 1, 0).* + c.bcoefr * s.r(je + 1, 1).*;
                    for (0..je) |jr| {
                        s.r(jr, 0).* = -creala * s.av(jr, je) + crealb * s.bv(jr, je) - cre2a * s.av(jr, je + 1) + cre2b * s.bv(jr, je + 1);
                        s.r(jr, 1).* = -cimaga * s.av(jr, je) + cimagb * s.bv(jr, je) - cim2a * s.av(jr, je + 1) + cim2b * s.bv(jr, je + 1);
                    }
                }
                const acoefa = @abs(c.acoef);
                const bcoefa = @abs(c.bcoefr) + @abs(c.bcoefi);
                const dmin = @max(safmin, @max(ulp * acoefa * an, ulp * bcoefa * bn));
                var j = je;
                while (j > 0) {
                    j -= 1;
                    const cplx = j > 0 and s.av(j, j - 1) != 0;
                    const na: usize = if (cplx) 2 else 1;
                    j -= na - 1;
                    if (cplx) {
                        solve2x2(T, dmin, c.acoef, .{ s.av(j, j), s.av(j, j + 1), s.av(j + 1, j), s.av(j + 1, j + 1) }, s.bv(j, j), s.bv(j + 1, j + 1), s, j, nw, c.bcoefr, c.bcoefi);
                    } else solve1x1(T, dmin, c.acoef, s.av(j, j), s.bv(j, j), s, j, nw, c.bcoefr, c.bcoefi);
                    for (0..na) |ja| {
                        if (pair) {
                            const creala = c.acoef * s.r(j + ja, 0).*;
                            const cimaga = c.acoef * s.r(j + ja, 1).*;
                            const crealb = c.bcoefr * s.r(j + ja, 0).* - c.bcoefi * s.r(j + ja, 1).*;
                            const cimagb = c.bcoefi * s.r(j + ja, 0).* + c.bcoefr * s.r(j + ja, 1).*;
                            for (0..j) |jr| {
                                s.r(jr, 0).* -= creala * s.av(jr, j + ja) - crealb * s.bv(jr, j + ja);
                                s.r(jr, 1).* -= cimaga * s.av(jr, j + ja) - cimagb * s.bv(jr, j + ja);
                            }
                        } else {
                            const creala = c.acoef * s.r(j + ja, 0).*;
                            const crealb = c.bcoefr * s.r(j + ja, 0).*;
                            for (0..j) |jr| s.r(jr, 0).* -= creala * s.av(jr, j + ja) - crealb * s.bv(jr, j + ja);
                        }
                    }
                }
                s.backTransform(u, ld, je, nw, 0, je + nw);
            }
        }
    };
}

/// `x := x / (ca a - (wr + i wi) d0)` for row `j` of `rhs` (`nw` columns,
/// the second the imaginary part): faer's `solve_complex_shifted_1x1`, below
/// `smin` the divisor `smin` alone (issue #766).
fn solve1x1(comptime T: type, smin_in: T, ca: T, a00: T, d0: T, s: anytype, j: usize, nw: usize, wr: T, wi: T) void {
    const smin = @max(smin_in, std.math.floatMin(T) * 2);
    const cr = ca * a00 - wr * d0;
    const ci = -wi * d0;
    if (@max(@abs(cr), @abs(ci)) < smin) {
        for (0..nw) |jw| s.r(j, jw).* /= smin;
        return;
    }
    if (nw == 1) {
        s.r(j, 0).* *= 1 / cr;
        return;
    }
    const inv = C(T).init(cr, ci).reciprocal();
    const re = s.r(j, 0).*;
    const im = s.r(j, 1).*;
    s.r(j, 0).* = re * inv.re - im * inv.im;
    s.r(j, 1).* = re * inv.im + im * inv.re;
}

/// `x := (ca A - (wr + i wi) diag(d0, d1))⁻¹ x` for rows `j`, `j + 1`
/// (`a` row-major), by complete pivoting: faer's
/// `solve_complex_shifted_2x2`, below `smin` the divisor `smin I` alone and
/// the second pivot floored at `smin` (`xLALN2`; issue #766).
fn solve2x2(comptime T: type, smin_in: T, ca: T, a: [4]T, d0: T, d1: T, s: anytype, j: usize, nw: usize, wr: T, wi: T) void {
    const Z = C(T);
    const smin = @max(smin_in, std.math.floatMin(T) * 2);
    const wim: T = if (nw == 1) 0 else wi;
    var m: [2][2]Z = .{
        .{ Z.init(ca * a[0] - wr * d0, -wim * d0), Z.init(ca * a[1], 0) },
        .{ Z.init(ca * a[2], 0), Z.init(ca * a[3] - wr * d1, -wim * d1) },
    };
    var cmax: T = 0;
    var pi: usize = 0;
    var pj: usize = 0;
    for (0..2) |i| for (0..2) |k| {
        const v = std.math.hypot(m[i][k].re, m[i][k].im);
        if (v > cmax) {
            cmax = v;
            pi = i;
            pj = k;
        }
    };
    if (cmax < smin) {
        for (0..nw) |jw| for (0..2) |i| {
            s.r(j + i, jw).* /= smin;
        };
        return;
    }
    var x: [2]Z = .{ Z.init(s.r(j, 0).*, if (nw == 2) s.r(j, 1).* else 0), Z.init(s.r(j + 1, 0).*, if (nw == 2) s.r(j + 1, 1).* else 0) };
    // move the pivot to (0, 0): rows by `pi`, columns by `pj`
    if (pi == 1) {
        std.mem.swap([2]Z, &m[0], &m[1]);
        std.mem.swap(Z, &x[0], &x[1]);
    }
    if (pj == 1) {
        std.mem.swap(Z, &m[0][0], &m[0][1]);
        std.mem.swap(Z, &m[1][0], &m[1][1]);
    }
    const l = m[1][0].div(m[0][0]);
    var piv2 = m[1][1].sub(l.mul(m[0][1]));
    if (std.math.hypot(piv2.re, piv2.im) < smin) piv2 = Z.init(smin, 0);
    const y1 = x[1].sub(l.mul(x[0]));
    const x1 = y1.div(piv2);
    const x0 = x[0].sub(m[0][1].mul(x1)).div(m[0][0]);
    // undo the column swap
    const sol: [2]Z = if (pj == 1) .{ x1, x0 } else .{ x0, x1 };
    for (0..2) |i| {
        s.r(j + i, 0).* = sol[i].re;
        if (nw == 2) s.r(j + i, 1).* = sol[i].im;
    }
}
