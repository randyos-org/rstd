//! The real QZ algorithm, `?HGEQZ` for real pencils: faer v0.24.4's
//! `hessenberg_to_qz_unblocked` with its `double_shift_sweep`,
//! `single_shift_sweep`, `svd_triu_2x2` (`xLASV2`) and
//! `generalized_eigval_2x2` (`xLAG2`) (`faer/src/linalg/gevd/qz_real/mod.rs`,
//! MIT), at every size; its multishift sweeps and aggressive early deflation
//! are issue #765. The block standardization and sweeps are in
//! `qz_real_steps.zig`, the 2 x 2 kernels in `qz_2x2.zig`. One place follows
//! Reference-LAPACK's `dhgeqz.f` (OpenBLAS v0.3.34-317) instead (issue
//! #764): a standardized complex pair's second `beta` is `B2A` (`dhgeqz.f`
//! line 1043), where faer repeats `B1A`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const qz_2x2 = @import("qz_2x2.zig");
const qz_real_steps = @import("qz_real_steps.zig");
const gghrd_reduce = @import("gghrd.zig");
const qz_complex = @import("qz_complex.zig");

const Io = std.Io;
const givens = gghrd_reduce.givens;
const rotRows = gghrd_reduce.rotRows;
const rotCols = gghrd_reduce.rotCols;

pub const SchurVectors = qz_complex.SchurVectors;
pub const SchurJob = qz_complex.SchurJob;
pub const HgeqzError = qz_complex.HgeqzError;

/// `?HGEQZ`, real: the real generalized Schur form `(S, P) = (Qᵀ H Z, Qᵀ T
/// Z)` of the Hessenberg-triangular `(H, T)` (`gghrd`'s): `P` upper
/// triangular with a nonnegative diagonal, `S` quasi-triangular with a 2 x 2
/// block per complex pair (`P`'s matching block diagonal). The eigenvalues
/// are `(alphar[j] + i alphai[j]) / beta[j]`, a pair together, the positive
/// imaginary part first. Otherwise as `hgeqzComplex`.
pub fn hgeqzReal(io: Io, comptime T: type, job: SchurJob, compq: SchurVectors, compz: SchurVectors, n: usize, lo: usize, hi: usize, h: []T, ldh: usize, t: []T, ldt: usize, alphar: []T, alphai: []T, beta: []T, q: []T, ldq: usize, z: []T, ldz: usize) HgeqzError!void {
    comptime ops.requireReal(T, "hgeqzReal");
    try interface.blockIn(lo, hi, n);
    try util.holds(h.len, n, n, ldh);
    try util.holds(t.len, n, n, ldt);
    if (alphar.len < n or alphai.len < n or beta.len < n) return error.SliceTooShort;
    if (compq != .none) try util.holds(q.len, n, n, ldq);
    if (compz != .none) try util.holds(z.len, n, n, ldz);
    if (compq == .from_identity) identity(T, n, q, ldq);
    if (compz == .from_identity) identity(T, n, z, ldz);
    if (n == 0) return;
    const s: Pencil(T) = .{
        .io = io,
        .n = n,
        .h = h,
        .ldh = ldh,
        .t = t,
        .ldt = ldt,
        .q = if (compq != .none) q else null,
        .ldq = ldq,
        .z = if (compz != .none) z else null,
        .ldz = ldz,
        .alphar = alphar,
        .alphai = alphai,
        .beta = beta,
        .eigvals_only = job == .eigenvalues,
    };
    for (hi..n) |j| s.outside(j);
    if (hi > lo) try s.iterate(lo, hi - 1);
    for (0..lo) |j| s.outside(j);
}

fn identity(comptime T: type, n: usize, m: []T, ld: usize) void {
    for (0..n) |j| for (0..n) |i| {
        m[i + j * ld] = if (i == j) 1 else 0;
    };
}

pub fn Pencil(comptime T: type) type {
    return struct {
        io: Io,
        n: usize,
        h: []T,
        ldh: usize,
        t: []T,
        ldt: usize,
        q: ?[]T,
        ldq: usize,
        z: ?[]T,
        ldz: usize,
        alphar: []T,
        alphai: []T,
        beta: []T,
        eigvals_only: bool,

        const Self = @This();
        const steps = qz_real_steps.Steps(T);

        pub fn hp(s: Self, i: usize, j: usize) *T {
            return &s.h[i + j * s.ldh];
        }

        pub fn tp(s: Self, i: usize, j: usize) *T {
            return &s.t[i + j * s.ldt];
        }

        /// Negates column `j` of `H` and `T`, rows `from..=j` (only the
        /// diagonal with `eigvals_only`), and of `Z`.
        pub fn negate(s: Self, j: usize, from: usize) void {
            if (!s.eigvals_only) {
                for (from..j + 1) |i| {
                    s.hp(i, j).* = -s.hp(i, j).*;
                    s.tp(i, j).* = -s.tp(i, j).*;
                }
            } else {
                s.hp(j, j).* = -s.hp(j, j).*;
                s.tp(j, j).* = -s.tp(j, j).*;
            }
            if (s.z) |z| for (0..s.n) |i| {
                z[i + j * s.ldz] = -z[i + j * s.ldz];
            };
        }

        /// A real eigenvalue at `j`, `T(j, j)` made nonnegative.
        fn real(s: Self, j: usize, from: usize) void {
            if (s.tp(j, j).* < 0) s.negate(j, from);
            s.alphar[j] = s.hp(j, j).*;
            s.alphai[j] = 0;
            s.beta[j] = s.tp(j, j).*;
        }

        fn outside(s: Self, j: usize) void {
            s.real(j, 0);
        }

        fn hasNan(m: []const T, ld: usize, lo: usize, hi: usize) bool {
            for (lo..hi + 1) |j| for (lo..hi + 1) |i| {
                if (std.math.isNan(m[i + j * ld])) return true;
            };
            return false;
        }

        fn maxAbs(m: []const T, ld: usize, lo: usize, hi: usize) T {
            var r: T = 0;
            for (lo..hi + 1) |j| for (lo..hi + 1) |i| {
                r = @max(r, @abs(m[i + j * ld]));
            };
            return r;
        }

        fn iterate(s: Self, ilo: usize, ihi: usize) HgeqzError!void {
            const n = s.n;
            if (hasNan(s.h, s.ldh, ilo, ihi) or hasNan(s.t, s.ldt, ilo, ihi)) return error.NoConvergence;
            const half: T = 0.5;
            const ulp = std.math.floatEps(T);
            const safmin = std.math.floatMin(T);
            const safmax = 1 / safmin;
            var ilast = ihi;
            var ifrstm: usize = if (s.eigvals_only) ilo else 0;
            var ilastm: usize = if (s.eigvals_only) ihi else n - 1;
            var iiter: usize = 0;
            const maxit = 30 * (ihi + 1 - ilo);
            var eshift: T = 0;
            const anorm = maxAbs(s.h, s.ldh, ilo, ihi);
            const bnorm = maxAbs(s.t, s.ldt, ilo, ihi);
            const atol = @max(safmin, ulp * anorm);
            const btol = @max(safmin, ulp * bnorm);
            const ascale = 1 / @max(safmin, anorm);
            const bscale = 1 / @max(safmin, bnorm);
            var ifirst: usize = ilo;

            for (0..maxit) |_| {
                const Step = enum { deflate, zero_t, sweep };
                const step: Step = find: {
                    if (ilast == ilo) break :find .deflate;
                    if (@abs(s.hp(ilast, ilast - 1).*) <= @max(safmin, ulp * (@abs(s.hp(ilast, ilast).*) + @abs(s.hp(ilast - 1, ilast - 1).*)))) {
                        s.hp(ilast, ilast - 1).* = 0;
                        break :find .deflate;
                    }
                    if (@abs(s.tp(ilast, ilast).*) <= btol) {
                        s.tp(ilast, ilast).* = 0;
                        break :find .zero_t;
                    }
                    var j = ilast;
                    while (j > ilo) {
                        j -= 1;
                        var ilazro: bool = undefined;
                        if (j == ilo) {
                            ilazro = true;
                        } else if (@abs(s.hp(j, j - 1).*) <= @max(safmin, ulp * (@abs(s.hp(j, j).*) + @abs(s.hp(j - 1, j - 1).*)))) {
                            s.hp(j, j - 1).* = 0;
                            ilazro = true;
                        } else ilazro = false;
                        if (@abs(s.tp(j, j).*) < btol) {
                            s.tp(j, j).* = 0;
                            var ilazr2 = false;
                            if (!ilazro) {
                                var temp = @abs(s.hp(j, j - 1).*);
                                var temp2 = @abs(s.hp(j, j).*);
                                const tempr = @max(temp, temp2);
                                if (tempr < 1 and tempr != 0) {
                                    temp /= tempr;
                                    temp2 /= tempr;
                                }
                                if (temp * (ascale * @abs(s.hp(j + 1, j).*)) <= temp2 * (ascale * atol)) ilazr2 = true;
                            }
                            if (ilazro or ilazr2) {
                                for (j..ilast) |jch| {
                                    const g = givens(T, s.hp(jch, jch).*, s.hp(jch + 1, jch).*);
                                    s.hp(jch, jch).* = g.r;
                                    s.hp(jch + 1, jch).* = 0;
                                    rotRows(T, g.c, g.s, s.h, s.ldh, jch, jch + 1, jch + 1, ilastm + 1);
                                    rotRows(T, g.c, g.s, s.t, s.ldt, jch, jch + 1, jch + 1, ilastm + 1);
                                    if (s.q) |q| rotCols(T, g.c, g.s, q, s.ldq, jch, jch + 1, 0, n);
                                    if (ilazr2) s.hp(jch, jch - 1).* *= g.c;
                                    ilazr2 = false;
                                    if (@abs(s.tp(jch + 1, jch + 1).*) >= btol) {
                                        if (jch + 1 >= ilast) break :find .deflate;
                                        ifirst = jch + 1;
                                        break :find .sweep;
                                    }
                                    s.tp(jch + 1, jch + 1).* = 0;
                                }
                                break :find .zero_t;
                            }
                            // chase the zero down to `T(ilast, ilast)`
                            for (j..ilast) |jch| {
                                const g = givens(T, s.tp(jch, jch + 1).*, s.tp(jch + 1, jch + 1).*);
                                s.tp(jch, jch + 1).* = g.r;
                                s.tp(jch + 1, jch + 1).* = 0;
                                if (jch + 1 < ilastm) rotRows(T, g.c, g.s, s.t, s.ldt, jch, jch + 1, jch + 2, ilastm + 1);
                                rotRows(T, g.c, g.s, s.h, s.ldh, jch, jch + 1, jch - 1, ilastm + 1);
                                if (s.q) |q| rotCols(T, g.c, g.s, q, s.ldq, jch, jch + 1, 0, n);
                                const g2 = givens(T, s.hp(jch + 1, jch).*, s.hp(jch + 1, jch - 1).*);
                                s.hp(jch + 1, jch).* = g2.r;
                                s.hp(jch + 1, jch - 1).* = 0;
                                rotCols(T, g2.c, g2.s, s.h, s.ldh, jch, jch - 1, ifrstm, jch + 1);
                                rotCols(T, g2.c, g2.s, s.t, s.ldt, jch, jch - 1, ifrstm, jch);
                                if (s.z) |z| rotCols(T, g2.c, g2.s, z, s.ldz, jch, jch - 1, 0, n);
                            }
                            break :find .zero_t;
                        } else if (ilazro) {
                            ifirst = j;
                            break :find .sweep;
                        }
                    }
                    // the reference's "impossible" drop-through (`INFO = N + 1`)
                    return error.NoConvergence;
                };

                if (step == .zero_t) {
                    const g = givens(T, s.hp(ilast, ilast).*, s.hp(ilast, ilast - 1).*);
                    s.hp(ilast, ilast).* = g.r;
                    s.hp(ilast, ilast - 1).* = 0;
                    rotCols(T, g.c, g.s, s.h, s.ldh, ilast, ilast - 1, ifrstm, ilast);
                    rotCols(T, g.c, g.s, s.t, s.ldt, ilast, ilast - 1, ifrstm, ilast);
                    if (s.z) |z| rotCols(T, g.c, g.s, z, s.ldz, ilast, ilast - 1, 0, n);
                }
                if (step != .sweep) {
                    s.real(ilast, ifrstm);
                    if (ilast == ilo) return;
                    ilast -= 1;
                    iiter = 0;
                    eshift = 0;
                    if (s.eigvals_only) {
                        ilastm = ilast;
                        if (ifrstm > ilast) ifrstm = ilo;
                    }
                    continue;
                }

                iiter += 1;
                if (s.eigvals_only) ifrstm = ifirst;
                var s1: T = undefined;
                var wr: T = undefined;
                const double = pick: {
                    if (iiter % 10 == 0) {
                        // an exceptional shift
                        const mx: T = @floatFromInt(maxit);
                        if (mx * safmin * @abs(s.hp(ilast, ilast - 1).*) < @abs(s.tp(ilast - 1, ilast - 1).*)) {
                            eshift = s.hp(ilast, ilast - 1).* / s.tp(ilast - 1, ilast - 1).*;
                        } else eshift += 1 / (safmin * mx);
                        s1 = 1;
                        wr = eshift;
                        break :pick false;
                    }
                    const e = qz_2x2.eig2x2(T, steps.block(s, s.h, s.ldh, ilast - 1), steps.block(s, s.t, s.ldt, ilast - 1));
                    s1 = e.scale1;
                    var s2 = e.scale2;
                    wr = e.wr1;
                    var wr2 = e.wr2;
                    const tl = s.tp(ilast, ilast).*;
                    const hl = s.hp(ilast, ilast).*;
                    if (@abs((wr / s1) * tl - hl) > @abs((wr2 / s2) * tl - hl)) {
                        std.mem.swap(T, &wr, &wr2);
                        std.mem.swap(T, &s1, &s2);
                    }
                    break :pick e.wi != 0;
                };
                if (!double) {
                    var temp = @min(ascale, 1) * (half * safmax);
                    var scale: T = if (s1 > temp) temp / s1 else 1;
                    temp = @min(bscale, 1) * (half * safmax);
                    if (@abs(wr) > temp) scale = @min(scale, temp / @abs(wr));
                    s1 *= scale;
                    wr *= scale;
                    var istart = ifirst;
                    var j = ilast - 1;
                    while (j > ifirst) : (j -= 1) {
                        var tmp = @abs(s1 * s.hp(j, j - 1).*);
                        var tmp2 = @abs(s1 * s.hp(j, j).* - wr * s.tp(j, j).*);
                        const tempr = @max(tmp, tmp2);
                        if (tempr < 1 and tempr != 0) {
                            tmp /= tempr;
                            tmp2 /= tempr;
                        }
                        if (@abs((ascale * s.hp(j + 1, j).*) * tmp) <= (ascale * atol) * tmp2) {
                            istart = j;
                            break;
                        }
                    }
                    steps.singleShiftSweep(s, s1, istart, wr, ilast, ilastm, ifrstm);
                    continue;
                }
                if (ifirst + 1 == ilast) {
                    // a 2 x 2 block with complex eigenvalues: standardize it
                    if (steps.standardize(s, ilast, ifirst, ifrstm, ilastm)) {
                        if (ifirst <= ilo) return;
                        ilast = ifirst - 1;
                        iiter = 0;
                        eshift = 0;
                        if (s.eigvals_only) {
                            ilastm = ilast;
                            if (ifrstm > ilast) ifrstm = ilo;
                        }
                    }
                    continue;
                }
                steps.doubleShiftSweep(s, ascale, ilast, bscale, ifirst, ilastm, ifrstm);
            }
            return error.NoConvergence;
        }
    };
}
