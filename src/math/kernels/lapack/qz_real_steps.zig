//! The real QZ's steps on a pencil (`qz_real.zig`'s `Pencil`): a 2 x 2 block's
//! standardization and the single- and double-shift sweeps, faer v0.24.4's
//! (`faer/src/linalg/gevd/qz_real/mod.rs`, MIT), as `qz_real.zig` says.

const std = @import("std");
const householder = @import("householder.zig");
const gghrd_reduce = @import("gghrd.zig");
const qz_real = @import("qz_real.zig");
const qz_2x2 = @import("qz_2x2.zig");

const givens = gghrd_reduce.givens;
const rotRows = gghrd_reduce.rotRows;
const rotCols = gghrd_reduce.rotCols;
const svdTriu2x2 = qz_2x2.svdTriu2x2;
const eig2x2 = qz_2x2.eig2x2;

pub fn Steps(comptime T: type) type {
    return struct {
        const Self = qz_real.Pencil(T);

        pub fn block(_: Self, m: []const T, ld: usize, k: usize) [4]T {
            return .{ m[k + k * ld], m[k + (k + 1) * ld], m[k + 1 + k * ld], m[k + 1 + (k + 1) * ld] };
        }

        /// The 2 x 2 block at `ifirst = ilast - 1`: `T`'s block made diagonal
        /// (`xLASV2`) and nonnegative; a pair of complex eigenvalues taken
        /// (true), or, the block's eigenvalues now real, false.
        pub fn standardize(s: Self, ilast: usize, ifirst: usize, ifrstm: usize, ilastm: usize) bool {
            const n = s.n;
            const safmin = std.math.floatMin(T);
            const sv = svdTriu2x2(T, s.tp(ilast - 1, ilast - 1).*, s.tp(ilast - 1, ilast).*, s.tp(ilast, ilast).*);
            var b22 = sv.ssmin;
            var b11 = sv.ssmax;
            var sr = sv.snr;
            var cr = sv.csr;
            const sl = sv.snl;
            const cl = sv.csl;
            if (b11 < 0) {
                cr = -cr;
                sr = -sr;
                b11 = -b11;
                b22 = -b22;
            }
            rotRows(T, cl, sl, s.h, s.ldh, ilast - 1, ilast, ilast - 1, ilastm + 1);
            rotCols(T, cr, sr, s.h, s.ldh, ilast - 1, ilast, ifrstm, ilast + 1);
            if (ilast < ilastm) rotRows(T, cl, sl, s.t, s.ldt, ilast - 1, ilast, ilast + 1, ilastm + 1);
            if (ifrstm + 1 < ilast) rotCols(T, cr, sr, s.t, s.ldt, ilast - 1, ilast, ifrstm, ifirst);
            if (s.q) |q| rotCols(T, cl, sl, q, s.ldq, ilast - 1, ilast, 0, n);
            if (s.z) |z| rotCols(T, cr, sr, z, s.ldz, ilast - 1, ilast, 0, n);
            s.tp(ilast - 1, ilast - 1).* = b11;
            s.tp(ilast - 1, ilast).* = 0;
            s.tp(ilast, ilast - 1).* = 0;
            s.tp(ilast, ilast).* = b22;
            if (b22 < 0) {
                for (ifrstm..ilast + 1) |j| {
                    s.hp(j, ilast).* = -s.hp(j, ilast).*;
                    s.tp(j, ilast).* = -s.tp(j, ilast).*;
                }
                if (s.z) |z| for (0..n) |j| {
                    z[j + ilast * s.ldz] = -z[j + ilast * s.ldz];
                };
                b22 = -b22;
            }
            const e = eig2x2(T, block(s, s.h, s.ldh, ilast - 1), block(s, s.t, s.ldt, ilast - 1));
            const s1 = e.scale1;
            const wr = e.wr1;
            const wi = e.wi;
            if (wi == 0) return false;
            const s1inv = 1 / s1;
            const a11 = s.hp(ilast - 1, ilast - 1).*;
            const a21 = s.hp(ilast, ilast - 1).*;
            const a12 = s.hp(ilast - 1, ilast).*;
            const a22 = s.hp(ilast, ilast).*;
            const c11r = s1 * a11 - wr * b11;
            const c11i = -wi * b11;
            const c12 = s1 * a12;
            const c21 = s1 * a21;
            const c22r = s1 * a22 - wr * b22;
            const c22i = -wi * b22;
            var cz: T = undefined;
            var szr: T = undefined;
            var szi: T = undefined;
            if (@abs(c11r) + @abs(c11i) + @abs(c12) > @abs(c21) + @abs(c22r) + @abs(c22i)) {
                const t1 = std.math.hypot(c12, std.math.hypot(c11r, c11i));
                cz = c12 / t1;
                szr = -c11r / t1;
                szi = -c11i / t1;
            } else {
                cz = std.math.hypot(c22r, c22i);
                if (cz <= safmin) {
                    cz = 0;
                    szr = 1;
                    szi = 0;
                } else {
                    const tempr = c22r / cz;
                    const tempi = c22i / cz;
                    const t1 = std.math.hypot(cz, c21);
                    cz /= t1;
                    szr = -c21 * tempr / t1;
                    szi = c21 * tempi / t1;
                }
            }
            const an = @abs(a11) + @abs(a12) + @abs(a21) + @abs(a22);
            const bn = @abs(b11) + @abs(b22);
            const wabs = @abs(wr) + @abs(wi);
            var cq: T = undefined;
            var sqr: T = undefined;
            var sqi: T = undefined;
            if (s1 * an > wabs * bn) {
                cq = cz * b11;
                sqr = szr * b22;
                sqi = -szi * b22;
            } else {
                const a1r = cz * a11 + szr * a12;
                const a1i = szi * a12;
                const a2r = cz * a21 + szr * a22;
                const a2i = szi * a22;
                cq = std.math.hypot(a1r, a1i);
                if (cq <= safmin) {
                    cq = 0;
                    sqr = 1;
                    sqi = 0;
                } else {
                    const tempr = a1r / cq;
                    const tempi = a1i / cq;
                    sqr = tempr * a2r + tempi * a2i;
                    sqi = tempi * a2r - tempr * a2i;
                }
            }
            const t1 = std.math.hypot(cq, std.math.hypot(sqr, sqi));
            cq /= t1;
            sqr /= t1;
            sqi /= t1;
            const tempr = sqr * szr - sqi * szi;
            const tempi = sqr * szi + sqi * szr;
            const b1r = cq * cz * b11 + tempr * b22;
            const b1i = tempi * b22;
            const b1a = std.math.hypot(b1r, b1i);
            const b2r = cq * cz * b22 + tempr * b11;
            const b2i = -tempi * b11;
            const b2a = std.math.hypot(b2r, b2i);
            s.beta[ilast - 1] = b1a;
            // `dhgeqz.f`'s `B2A` (faer repeats `b1a`; issue #764)
            s.beta[ilast] = b2a;
            s.alphar[ilast - 1] = (wr * b1a) * s1inv;
            s.alphai[ilast - 1] = (wi * b1a) * s1inv;
            s.alphar[ilast] = (wr * b2a) * s1inv;
            s.alphai[ilast] = -(wi * b2a) * s1inv;
            return true;
        }

        pub fn singleShiftSweep(s: Self, s1: T, istart: usize, wr: T, ilast: usize, ilastm: usize, ifrstm: usize) void {
            const n = s.n;
            const g0 = givens(T, s1 * s.hp(istart, istart).* - wr * s.tp(istart, istart).*, s1 * s.hp(istart + 1, istart).*);
            var c = g0.c;
            var sn = g0.s;
            for (istart..ilast) |j| {
                if (j > istart) {
                    const g = givens(T, s.hp(j, j - 1).*, s.hp(j + 1, j - 1).*);
                    c = g.c;
                    sn = g.s;
                    s.hp(j, j - 1).* = g.r;
                    s.hp(j + 1, j - 1).* = 0;
                }
                rotRows(T, c, sn, s.h, s.ldh, j, j + 1, j, ilastm + 1);
                rotRows(T, c, sn, s.t, s.ldt, j, j + 1, j, ilastm + 1);
                if (s.q) |q| rotCols(T, c, sn, q, s.ldq, j, j + 1, 0, n);
                const g = givens(T, s.tp(j + 1, j + 1).*, s.tp(j + 1, j).*);
                c = g.c;
                sn = g.s;
                s.tp(j + 1, j + 1).* = g.r;
                s.tp(j + 1, j).* = 0;
                rotCols(T, c, sn, s.h, s.ldh, j + 1, j, ifrstm, @min(j + 2, ilast) + 1);
                rotCols(T, c, sn, s.t, s.ldt, j + 1, j, ifrstm, j + 1);
                if (s.z) |z| rotCols(T, c, sn, z, s.ldz, j + 1, j, 0, n);
            }
        }

        /// `I - tau v vᵀ` on rows `j..j + 3` of `m`, columns `from..to`.
        pub fn reflectRows(m: []T, ld: usize, j: usize, tau: T, v2: T, v3: T, from: usize, to: usize) void {
            const t2 = tau * v2;
            const t3 = tau * v3;
            for (from..to) |jc| {
                const temp = m[j + jc * ld] + v2 * m[j + 1 + jc * ld] + v3 * m[j + 2 + jc * ld];
                m[j + jc * ld] -= temp * tau;
                m[j + 1 + jc * ld] -= temp * t2;
                m[j + 2 + jc * ld] -= temp * t3;
            }
        }

        /// The same on columns `j..j + 3`, rows `from..to`.
        pub fn reflectCols(m: []T, ld: usize, j: usize, tau: T, v2: T, v3: T, from: usize, to: usize) void {
            const t2 = tau * v2;
            const t3 = tau * v3;
            for (from..to) |jr| {
                const temp = m[jr + j * ld] + v2 * m[jr + (j + 1) * ld] + v3 * m[jr + (j + 2) * ld];
                m[jr + j * ld] -= temp * tau;
                m[jr + (j + 1) * ld] -= temp * t2;
                m[jr + (j + 2) * ld] -= temp * t3;
            }
        }

        /// faer's `make_householder_in_place` on `v`, as LAPACK's `tau`.
        pub fn reflector(s: Self, v: *[3]T) T {
            const info = householder.make(s.io, T, &v[0], 2, v[1..].ptr, 1, null, 0);
            return 1 / info.tau;
        }

        pub fn doubleShiftSweep(s: Self, ascale: T, ilast: usize, bscale: T, ifirst: usize, ilastm: usize, ifrstm: usize) void {
            const n = s.n;
            const safmin = std.math.floatMin(T);
            const hv = struct {
                fn f(p: Self, i: usize, j: usize) T {
                    return p.hp(i, j).*;
                }
            }.f;
            const tv = struct {
                fn f(p: Self, i: usize, j: usize) T {
                    return p.tp(i, j).*;
                }
            }.f;
            const ad11 = (ascale * hv(s, ilast - 1, ilast - 1)) / (bscale * tv(s, ilast - 1, ilast - 1));
            const ad21 = (ascale * hv(s, ilast, ilast - 1)) / (bscale * tv(s, ilast - 1, ilast - 1));
            const ad12 = (ascale * hv(s, ilast - 1, ilast)) / (bscale * tv(s, ilast, ilast));
            const ad22 = (ascale * hv(s, ilast, ilast)) / (bscale * tv(s, ilast, ilast));
            const upper12 = tv(s, ilast - 1, ilast) / tv(s, ilast, ilast);
            const ad11l = (ascale * hv(s, ifirst, ifirst)) / (bscale * tv(s, ifirst, ifirst));
            const ad21l = (ascale * hv(s, ifirst + 1, ifirst)) / (bscale * tv(s, ifirst, ifirst));
            const ad12l = (ascale * hv(s, ifirst, ifirst + 1)) / (bscale * tv(s, ifirst + 1, ifirst + 1));
            const ad22l = (ascale * hv(s, ifirst + 1, ifirst + 1)) / (bscale * tv(s, ifirst + 1, ifirst + 1));
            const ad32l = (ascale * hv(s, ifirst + 2, ifirst + 1)) / (bscale * tv(s, ifirst + 1, ifirst + 1));
            const upper12l = tv(s, ifirst, ifirst + 1) / tv(s, ifirst + 1, ifirst + 1);
            var v: [3]T = .{
                (ad11 - ad11l) * (ad22 - ad11l) - ad12 * ad21 + ad21 * upper12 * ad11l + (ad12l - ad11l * upper12l) * ad21l,
                ((ad22l - ad11l) - ad21l * upper12l - (ad11 - ad11l) - (ad22 - ad11l) + ad21 * upper12) * ad21l,
                ad32l * ad21l,
            };
            var tau = reflector(s, &v);
            v[0] = 1;
            const istart = ifirst;
            for (istart..ilast - 1) |j| {
                if (j > istart) {
                    v = .{ hv(s, j, j - 1), hv(s, j + 1, j - 1), hv(s, j + 2, j - 1) };
                    tau = reflector(s, &v);
                    s.hp(j, j - 1).* = v[0];
                    v[0] = 1;
                    s.hp(j + 1, j - 1).* = 0;
                    s.hp(j + 2, j - 1).* = 0;
                }
                reflectRows(s.h, s.ldh, j, tau, v[1], v[2], j, ilastm + 1);
                reflectRows(s.t, s.ldt, j, tau, v[1], v[2], j, ilastm + 1);
                if (s.q) |q| reflectCols(q, s.ldq, j, tau, v[1], v[2], 0, n);
                // zero `T(j + 1, j)` and `T(j + 2, j)` with a reflection from the right
                var ilpivt = false;
                const temp = @max(@abs(tv(s, j + 1, j + 1)), @abs(tv(s, j + 1, j + 2)));
                const temp2 = @max(@abs(tv(s, j + 2, j + 1)), @abs(tv(s, j + 2, j + 2)));
                var scale: T = undefined;
                var e1: T = undefined;
                var e2: T = undefined;
                solve: {
                    var w11: T = undefined;
                    var w12: T = undefined;
                    var w21: T = undefined;
                    var w22: T = undefined;
                    if (@max(temp, temp2) < safmin) {
                        scale = 0;
                        e1 = 1;
                        e2 = 0;
                        break :solve;
                    } else if (temp >= temp2) {
                        w11 = tv(s, j + 1, j + 1);
                        w21 = tv(s, j + 2, j + 1);
                        w12 = tv(s, j + 1, j + 2);
                        w22 = tv(s, j + 2, j + 2);
                        e1 = tv(s, j + 1, j);
                        e2 = tv(s, j + 2, j);
                    } else {
                        w21 = tv(s, j + 1, j + 1);
                        w11 = tv(s, j + 2, j + 1);
                        w22 = tv(s, j + 1, j + 2);
                        w12 = tv(s, j + 2, j + 2);
                        e2 = tv(s, j + 1, j);
                        e1 = tv(s, j + 2, j);
                    }
                    if (@abs(w12) > @abs(w11)) {
                        ilpivt = true;
                        std.mem.swap(T, &w11, &w12);
                        std.mem.swap(T, &w21, &w22);
                    }
                    const ratio = w21 / w11;
                    e2 -= ratio * e1;
                    w22 -= ratio * w12;
                    scale = 1;
                    if (@abs(w22) < safmin) {
                        scale = 0;
                        e2 = 1;
                        e1 = -w12 / w11;
                        break :solve;
                    }
                    if (@abs(w22) < @abs(e2)) scale = @abs(w22 / e2);
                    if (@abs(w11) < @abs(e1)) scale = @min(scale, @abs(w11 / e1));
                    e2 = (scale * e2) / w22;
                    e1 = (scale * e1 - w12 * e2) / w11;
                }
                if (ilpivt) std.mem.swap(T, &e1, &e2);
                const t1 = @sqrt(scale * scale + e1 * e1 + e2 * e2);
                tau = 1 + scale / t1;
                const vs = -1 / (scale + t1);
                v = .{ 1, vs * e1, vs * e2 };
                reflectCols(s.h, s.ldh, j, tau, v[1], v[2], ifrstm, @min(j + 3, ilast) + 1);
                reflectCols(s.t, s.ldt, j, tau, v[1], v[2], ifrstm, j + 3);
                if (s.z) |z| reflectCols(z, s.ldz, j, tau, v[1], v[2], 0, n);
                s.tp(j + 1, j).* = 0;
                s.tp(j + 2, j).* = 0;
            }
            // the last step: rotations, the bulge's 2 x 2 tail
            const j = ilast - 1;
            const g = givens(T, hv(s, j, j - 1), hv(s, j + 1, j - 1));
            s.hp(j, j - 1).* = g.r;
            s.hp(j + 1, j - 1).* = 0;
            rotRows(T, g.c, g.s, s.h, s.ldh, j, j + 1, j, ilastm + 1);
            rotRows(T, g.c, g.s, s.t, s.ldt, j, j + 1, j, ilastm + 1);
            if (s.q) |q| rotCols(T, g.c, g.s, q, s.ldq, j, j + 1, 0, n);
            const g2 = givens(T, tv(s, j + 1, j + 1), tv(s, j + 1, j));
            s.tp(j + 1, j + 1).* = g2.r;
            s.tp(j + 1, j).* = 0;
            rotCols(T, g2.c, g2.s, s.h, s.ldh, j + 1, j, ifrstm, ilast + 1);
            rotCols(T, g2.c, g2.s, s.t, s.ldt, j + 1, j, ifrstm, ilast);
            if (s.z) |z| rotCols(T, g2.c, g2.s, z, s.ldz, j + 1, j, 0, n);
        }
    };
}
