//! The complex QZ algorithm, `?HGEQZ` for complex pencils: faer v0.24.4's
//! `hessenberg_to_qz_unblocked` (`faer/src/linalg/gevd/qz_cplx/mod.rs`,
//! MIT), single-shift, at every size; its multishift sweeps and aggressive
//! early deflation are issue #762. Three places follow Reference-LAPACK's
//! `zhgeqz.f` (OpenBLAS v0.3.34-317) instead, where faer departs from it
//! (issue #763): the chase of a zero `T(j, j)` rotates `Q`'s columns, not
//! its rows (`zhgeqz.f` line 582), and the columns outside `lo..hi` scale
//! `J - 1` entries of `T` and `J` of `H` (lines 455, 456), one more each.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const gghrd_reduce = @import("gghrd.zig");

const Float = types.Float;
const scalar = ops.scalar;
const add = ops.add;
const sub = ops.sub;
const mul = ops.mul;
const div = ops.div;
const conj = ops.conj;
const mulReal = ops.mulReal;
const abs1 = ops.abs1;
const givens = gghrd_reduce.givens;
const rotRows = gghrd_reduce.rotRows;
const rotCols = gghrd_reduce.rotCols;

const schur = @import("schur.zig");
pub const SchurVectors = schur.SchurVectors;
/// The eigenvalues only (`H` and `T` then unspecified), or the generalized
/// Schur form too: `hseqr`'s `JOB`.
pub const SchurJob = schur.SchurJob;

/// What the QZ returns: `BlockOutOfRange` unless `lo <= hi <= n`, a slice
/// too short or a leading dimension too small; `NoConvergence` when the
/// iteration limit (`30` per row of the block) runs out, a NaN included.
pub const HgeqzError = interface.BlockError || interface.ShapeError || error{NoConvergence};

/// `?HGEQZ`, complex: the generalized Schur form `(S, P) = (Qᴴ H Z, Qᴴ T
/// Z)` of the Hessenberg-triangular `(H, T)` (`gghrd`'s), both upper
/// triangular with `P`'s diagonal real and nonnegative, the eigenvalues
/// `alpha[j] / beta[j]` (`n` each); rows and columns outside `lo..hi` are
/// already triangular (`ggbal`'s). `q` and `z` take or are multiplied by
/// the transformations, as `gghrd`'s.
pub fn hgeqzComplex(comptime T: type, job: SchurJob, compq: SchurVectors, compz: SchurVectors, n: usize, lo: usize, hi: usize, h: []T, ldh: usize, t: []T, ldt: usize, alpha: []T, beta: []T, q: []T, ldq: usize, z: []T, ldz: usize) HgeqzError!void {
    comptime ops.requireComplex(T, "hgeqzComplex");
    try interface.blockIn(lo, hi, n);
    try util.holds(h.len, n, n, ldh);
    try util.holds(t.len, n, n, ldt);
    if (alpha.len < n or beta.len < n) return error.SliceTooShort;
    if (compq != .none) try util.holds(q.len, n, n, ldq);
    if (compz != .none) try util.holds(z.len, n, n, ldz);
    if (compq == .from_identity) identity(T, n, q, ldq);
    if (compz == .from_identity) identity(T, n, z, ldz);
    if (n == 0) return;
    const s: Pencil(T) = .{
        .n = n,
        .h = h,
        .ldh = ldh,
        .t = t,
        .ldt = ldt,
        .q = if (compq != .none) q else null,
        .ldq = ldq,
        .z = if (compz != .none) z else null,
        .ldz = ldz,
        .alpha = alpha,
        .beta = beta,
        .eigvals_only = job == .eigenvalues,
    };
    for (hi..n) |j| s.fixColumn(j);
    if (hi > lo) try s.iterate(lo, hi - 1);
    for (0..lo) |j| s.fixColumn(j);
}

fn identity(comptime T: type, n: usize, m: []T, ld: usize) void {
    for (0..n) |j| for (0..n) |i| {
        m[i + j * ld] = scalar(T, if (i == j) 1 else 0);
    };
}

fn Pencil(comptime T: type) type {
    return struct {
        n: usize,
        h: []T,
        ldh: usize,
        t: []T,
        ldt: usize,
        q: ?[]T,
        ldq: usize,
        z: ?[]T,
        ldz: usize,
        alpha: []T,
        beta: []T,
        eigvals_only: bool,

        const Self = @This();
        const F = Float(T);

        fn hp(s: Self, i: usize, j: usize) *T {
            return &s.h[i + j * s.ldh];
        }

        fn tp(s: Self, i: usize, j: usize) *T {
            return &s.t[i + j * s.ldt];
        }

        /// Scales column `j` so `T(j, j)` is real and nonnegative, rows
        /// `from..` of `T` above the diagonal and of `H` to it, then takes
        /// the eigenvalue.
        fn normalize(s: Self, j: usize, from: usize) void {
            const absb = ops.abs(T, s.tp(j, j).*);
            if (absb > std.math.floatMin(F)) {
                const signbc = mulReal(T, conj(T, s.tp(j, j).*), 1 / absb);
                s.tp(j, j).* = scalar(T, absb);
                if (!s.eigvals_only) {
                    for (from..j) |i| s.tp(i, j).* = mul(T, s.tp(i, j).*, signbc);
                    for (from..j + 1) |i| s.hp(i, j).* = mul(T, s.hp(i, j).*, signbc);
                } else s.hp(j, j).* = mul(T, s.hp(j, j).*, signbc);
                if (s.z) |z| for (0..s.n) |i| {
                    z[i + j * s.ldz] = mul(T, z[i + j * s.ldz], signbc);
                };
            } else s.tp(j, j).* = scalar(T, 0);
            s.alpha[j] = s.hp(j, j).*;
            s.beta[j] = s.tp(j, j).*;
        }

        /// A column outside the block: the reference's `ZSCAL(J - 1, ...)`
        /// and `ZSCAL(J, ...)`, from row 0.
        fn fixColumn(s: Self, j: usize) void {
            s.normalize(j, 0);
        }

        fn maxAbs(m: []const T, ld: usize, lo: usize, hi: usize) F {
            var r: F = 0;
            for (lo..hi + 1) |j| for (lo..hi + 1) |i| {
                r = @max(r, ops.abs(T, m[i + j * ld]));
            };
            return r;
        }

        fn iterate(s: Self, ilo: usize, ihi: usize) HgeqzError!void {
            const n = s.n;
            const ulp = std.math.floatEps(F);
            const safmin = std.math.floatMin(F);
            var ilast = ihi;
            var ifrstm: usize = if (s.eigvals_only) ilo else 0;
            var ilastm: usize = if (s.eigvals_only) ihi else n - 1;
            var iiter: usize = 0;
            const maxit = 30 * (ihi + 1 - ilo);
            var eshift = scalar(T, 0);
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
                    if (abs1(T, s.hp(ilast, ilast - 1).*) <= atol) {
                        s.hp(ilast, ilast - 1).* = scalar(T, 0);
                        break :find .deflate;
                    }
                    if (abs1(T, s.tp(ilast, ilast).*) <= btol) {
                        s.tp(ilast, ilast).* = scalar(T, 0);
                        break :find .zero_t;
                    }
                    var j = ilast;
                    while (j > ilo) {
                        j -= 1;
                        var ilazro: bool = undefined;
                        if (j == ilo) {
                            ilazro = true;
                        } else if (abs1(T, s.hp(j, j - 1).*) < atol) {
                            s.hp(j, j - 1).* = scalar(T, 0);
                            ilazro = true;
                        } else ilazro = false;
                        if (abs1(T, s.tp(j, j).*) < btol) {
                            s.tp(j, j).* = scalar(T, 0);
                            var ilazr2 = false;
                            if (!ilazro) {
                                if (abs1(T, s.hp(j, j - 1).*) * (ascale * abs1(T, s.hp(j + 1, j).*)) < abs1(T, s.hp(j, j).*) * (ascale * atol)) ilazr2 = true;
                            }
                            if (ilazro or ilazr2) {
                                for (j..ilast) |jch| {
                                    const g = givens(T, s.hp(jch, jch).*, s.hp(jch + 1, jch).*);
                                    s.hp(jch, jch).* = g.r;
                                    s.hp(jch + 1, jch).* = scalar(T, 0);
                                    rotRows(T, g.c, g.s, s.h, s.ldh, jch, jch + 1, jch + 1, ilastm + 1);
                                    rotRows(T, g.c, g.s, s.t, s.ldt, jch, jch + 1, jch + 1, ilastm + 1);
                                    // `zhgeqz.f`'s columns (faer rotates rows here; issue #763)
                                    if (s.q) |q| rotCols(T, g.c, conj(T, g.s), q, s.ldq, jch, jch + 1, 0, n);
                                    if (ilazr2) s.hp(jch, jch - 1).* = mulReal(T, s.hp(jch, jch - 1).*, g.c);
                                    ilazr2 = false;
                                    if (abs1(T, s.tp(jch + 1, jch + 1).*) >= btol) {
                                        if (jch + 1 >= ilast) break :find .deflate;
                                        ifirst = jch + 1;
                                        break :find .sweep;
                                    }
                                    s.tp(jch + 1, jch + 1).* = scalar(T, 0);
                                }
                                break :find .zero_t;
                            }
                            // chase the zero down to `T(ilast, ilast)`
                            for (j..ilast) |jch| {
                                const g = givens(T, s.tp(jch, jch + 1).*, s.tp(jch + 1, jch + 1).*);
                                s.tp(jch, jch + 1).* = g.r;
                                s.tp(jch + 1, jch + 1).* = scalar(T, 0);
                                if (jch + 1 < ilastm) rotRows(T, g.c, g.s, s.t, s.ldt, jch, jch + 1, jch + 2, ilastm + 1);
                                rotRows(T, g.c, g.s, s.h, s.ldh, jch, jch + 1, jch - 1, ilastm + 1);
                                if (s.q) |q| rotCols(T, g.c, conj(T, g.s), q, s.ldq, jch, jch + 1, 0, n);
                                const g2 = givens(T, s.hp(jch + 1, jch).*, s.hp(jch + 1, jch - 1).*);
                                s.hp(jch + 1, jch).* = g2.r;
                                s.hp(jch + 1, jch - 1).* = scalar(T, 0);
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
                    // `T(ilast, ilast)` is zero: a rotation zeroes `H(ilast, ilast - 1)`
                    const g = givens(T, s.hp(ilast, ilast).*, s.hp(ilast, ilast - 1).*);
                    s.hp(ilast, ilast).* = g.r;
                    s.hp(ilast, ilast - 1).* = scalar(T, 0);
                    rotCols(T, g.c, g.s, s.h, s.ldh, ilast, ilast - 1, ifrstm, ilast);
                    rotCols(T, g.c, g.s, s.t, s.ldt, ilast, ilast - 1, ifrstm, ilast);
                    if (s.z) |z| rotCols(T, g.c, g.s, z, s.ldz, ilast, ilast - 1, 0, n);
                }
                if (step != .sweep) {
                    // `H(ilast, ilast - 1)` is negligible: `ilast` deflates
                    s.normalize(ilast, ifrstm);
                    if (ilast == ilo) return;
                    ilast -= 1;
                    iiter = 0;
                    eshift = scalar(T, 0);
                    if (s.eigvals_only) {
                        ilastm = ilast;
                        if (ifrstm > ilast) ifrstm = ilo;
                    }
                    continue;
                }

                iiter += 1;
                if (s.eigvals_only) ifrstm = ifirst;
                const shift = if (iiter % 10 != 0) wilkinson(s, ilast, ascale, bscale) else blk: {
                    // an exceptional shift
                    eshift = add(T, eshift, mul(T, mulReal(T, s.hp(ilast, ilast - 1).*, ascale), recip(mulReal(T, s.tp(ilast - 1, ilast - 1).*, bscale))));
                    break :blk eshift;
                };
                // the start of the sweep: the lowest small subdiagonal of `H - shift T`
                var istart = ifirst;
                var ctemp = sub(T, mulReal(T, s.hp(ifirst, ifirst).*, ascale), mul(T, shift, mulReal(T, s.tp(ifirst, ifirst).*, bscale)));
                var j = ilast - 1;
                while (j > ifirst) : (j -= 1) {
                    const c = sub(T, mulReal(T, s.hp(j, j).*, ascale), mul(T, shift, mulReal(T, s.tp(j, j).*, bscale)));
                    var temp = abs1(T, c);
                    var temp2 = ascale * abs1(T, s.hp(j + 1, j).*);
                    const tempr = @max(temp, temp2);
                    if (tempr < 1 and tempr != 0) {
                        temp /= tempr;
                        temp2 /= tempr;
                    }
                    if (abs1(T, s.hp(j, j - 1).*) * temp2 <= temp * atol) {
                        istart = j;
                        ctemp = c;
                        break;
                    }
                }
                const g0 = givens(T, ctemp, mulReal(T, s.hp(istart + 1, istart).*, ascale));
                var c = g0.c;
                var sn = g0.s;
                for (istart..ilast) |jj| {
                    if (jj > istart) {
                        const g = givens(T, s.hp(jj, jj - 1).*, s.hp(jj + 1, jj - 1).*);
                        c = g.c;
                        sn = g.s;
                        s.hp(jj, jj - 1).* = g.r;
                        s.hp(jj + 1, jj - 1).* = scalar(T, 0);
                    }
                    rotRows(T, c, sn, s.h, s.ldh, jj, jj + 1, jj, ilastm + 1);
                    rotRows(T, c, sn, s.t, s.ldt, jj, jj + 1, jj, ilastm + 1);
                    if (s.q) |q| rotCols(T, c, conj(T, sn), q, s.ldq, jj, jj + 1, 0, n);
                    const g = givens(T, s.tp(jj + 1, jj + 1).*, s.tp(jj + 1, jj).*);
                    c = g.c;
                    sn = g.s;
                    s.tp(jj + 1, jj + 1).* = g.r;
                    s.tp(jj + 1, jj).* = scalar(T, 0);
                    rotCols(T, c, sn, s.h, s.ldh, jj + 1, jj, ifrstm, @min(jj + 3, ilast + 1));
                    rotCols(T, c, sn, s.t, s.ldt, jj + 1, jj, ifrstm, jj + 1);
                    if (s.z) |z| rotCols(T, c, sn, z, s.ldz, jj + 1, jj, 0, n);
                }
            }
            return error.NoConvergence;
        }

        fn recip(x: T) T {
            return div(T, scalar(T, 1), x);
        }

        /// The eigenvalue of the trailing 2 x 2 of `(H, T)` (scaled) nearer
        /// its corner, as faer and `zhgeqz.f` choose it.
        fn wilkinson(s: Self, ilast: usize, ascale: F, bscale: F) T {
            const tl = mulReal(T, s.tp(ilast, ilast).*, bscale);
            const tm = mulReal(T, s.tp(ilast - 1, ilast - 1).*, bscale);
            const upper12 = mul(T, mulReal(T, s.tp(ilast - 1, ilast).*, bscale), recip(tl));
            const ad11 = mul(T, mulReal(T, s.hp(ilast - 1, ilast - 1).*, ascale), recip(tm));
            const ad21 = mul(T, mulReal(T, s.hp(ilast, ilast - 1).*, ascale), recip(tm));
            const ad12 = mul(T, mulReal(T, s.hp(ilast - 1, ilast).*, ascale), recip(tl));
            const ad22 = mul(T, mulReal(T, s.hp(ilast, ilast).*, ascale), recip(tl));
            const abi22 = sub(T, ad22, mul(T, upper12, ad21));
            const t1 = mulReal(T, add(T, ad11, abi22), 0.5);
            const rtdisc = ops.sqrt(T, sub(T, add(T, mul(T, t1, t1), mul(T, ad12, ad21)), mul(T, ad11, ad22)));
            const diff = sub(T, t1, abi22);
            const temp = diff.re * rtdisc.re + diff.im * rtdisc.im;
            return if (temp <= 0) add(T, t1, rtdisc) else sub(T, t1, rtdisc);
        }
    };
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const pencil = @import("_test/pencil.zig");

/// `gghrd` then the QZ on a random pencil (`B` upper triangular, with
/// `zero_diag` rows of its diagonal zero); checks the Schur form. With
/// `direct`, `A` is made Hessenberg and the QZ runs on it straight, so
/// `B`'s zeros reach it in place (`gghrd`'s rotations move them).
fn qzCase(comptime T: type, n: usize, zero_diag: []const usize, seed: u64, direct: bool) !void {
    const al = testing.allocator;
    const F = Float(T);
    const a0 = try al.alloc(T, n * n);
    defer al.free(a0);
    const b0 = try al.alloc(T, n * n);
    defer al.free(b0);
    check.fill(T, a0, seed);
    check.fill(T, b0, seed + 1);
    for (0..n) |j| for (j + 1..n) |i| {
        b0[i + j * n] = scalar(T, 0);
    };
    for (zero_diag) |k| b0[k + k * n] = scalar(T, 0);
    if (direct) for (0..n) |i| for (0..i -| 1) |j| {
        a0[i + j * n] = scalar(T, 0);
    };
    const h = try al.dupe(T, a0);
    defer al.free(h);
    const t = try al.dupe(T, b0);
    defer al.free(t);
    const q = try al.alloc(T, n * n);
    defer al.free(q);
    const z = try al.alloc(T, n * n);
    defer al.free(z);
    const alpha = try al.alloc(T, n);
    defer al.free(alpha);
    const beta = try al.alloc(T, n);
    defer al.free(beta);
    if (!direct) try gghrd_reduce.gghrd(T, .from_identity, .from_identity, n, 0, n, h, n, t, n, q, n, z, n);
    const v: SchurVectors = if (direct) .from_identity else .accumulate;
    try hgeqzComplex(T, .schur_form, v, v, n, 0, n, h, n, t, n, alpha, beta, q, n, z, n);
    for (0..n) |j| {
        for (j + 1..n) |i| {
            try testing.expect(ops.isZero(T, h[i + j * n]));
            try testing.expect(ops.isZero(T, t[i + j * n]));
        }
        try testing.expect(t[j + j * n].im == 0 and t[j + j * n].re >= 0);
        try testing.expectEqual(h[j + j * n], alpha[j]);
        try testing.expectEqual(t[j + j * n], beta[j]);
    }
    const tol = 60 * @as(F, @floatFromInt(n * n)) * std.math.floatEps(F);
    try testing.expect(pencil.unitary(T, n, q) <= tol);
    try testing.expect(pencil.unitary(T, n, z) <= tol);
    try testing.expect(pencil.residual(T, n, q, h, z, a0) <= tol * 4);
    try testing.expect(pencil.residual(T, n, q, t, z, b0) <= tol * 4);
    // a zero on `B`'s diagonal makes it singular: an infinite eigenvalue
    // (`n - rank(B)` of them, generically)
    var infinite: usize = 0;
    for (beta) |x| {
        if (ops.abs(T, x) <= tol) infinite += 1;
    }
    if (zero_diag.len > 0) try testing.expect(infinite >= 1);

    // the eigenvalues alone: the same `alpha / beta`, to rounding
    @memcpy(h, a0);
    @memcpy(t, b0);
    const alpha2 = try al.alloc(T, n);
    defer al.free(alpha2);
    const beta2 = try al.alloc(T, n);
    defer al.free(beta2);
    if (!direct) try gghrd_reduce.gghrd(T, .none, .none, n, 0, n, h, n, t, n, q, n, z, n);
    try hgeqzComplex(T, .eigenvalues, .none, .none, n, 0, n, h, n, t, n, alpha2, beta2, q, n, z, n);
    if (zero_diag.len == 0) for (alpha2, beta2) |x2, y2| {
        const l2 = div(T, x2, y2);
        var best: F = std.math.inf(F);
        for (alpha, beta) |x, y| best = @min(best, ops.abs(T, sub(T, div(T, x, y), l2)));
        try testing.expect(best <= 1e3 * tol * (1 + ops.abs(T, l2)));
    };
}

test "hgeqzComplex: the generalized Schur form, regular and with infinite eigenvalues" {
    inline for (.{ types.Complex32, types.Complex64 }) |T| {
        for ([_]usize{ 1, 2, 3, 7, 16, 30 }) |n| try qzCase(T, n, &.{}, n * 13, false);
        for ([_]bool{ false, true }) |direct| {
            try qzCase(T, 6, &.{2}, 77, direct);
            try qzCase(T, 9, &.{ 0, 4, 8 }, 78, direct);
            try qzCase(T, 5, &.{ 0, 1, 2, 3, 4 }, 79, direct);
            try qzCase(T, 7, &.{0}, 80, direct);
        }
    }
}

test "hgeqzComplex: errors" {
    const T = types.Complex64;
    var h = [_]T{ .init(1, 0), .init(std.math.nan(f64), 0), .init(1, 0), .init(1, 0) };
    var t = [_]T{ .init(1, 0), .init(0, 0), .init(0, 0), .init(1, 0) };
    var alpha: [2]T = undefined;
    var beta: [2]T = undefined;
    var none: [1]T = undefined;
    try testing.expectError(error.NoConvergence, hgeqzComplex(T, .eigenvalues, .none, .none, 2, 0, 2, &h, 2, &t, 2, &alpha, &beta, &none, 1, &none, 1));
    try testing.expectError(error.BlockOutOfRange, hgeqzComplex(T, .eigenvalues, .none, .none, 2, 1, 3, &h, 2, &t, 2, &alpha, &beta, &none, 1, &none, 1));
    try testing.expectError(error.SliceTooShort, hgeqzComplex(T, .eigenvalues, .none, .none, 2, 0, 2, &h, 2, &t, 2, alpha[0..1], &beta, &none, 1, &none, 1));
}
