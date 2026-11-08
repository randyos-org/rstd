//! Rank-revealing least squares, `gelsy`, from Reference-LAPACK 3.12.0 as
//! OpenBLAS v0.3.34-317 ships it (`lapack-netlib/SRC/`, modified BSD):
//! `xGELSY` with `xLAIC1`, `xLATRZ`, `xLARZ` and `xORMR3`, over faer's pivoted
//! QR (`qr_pivot.zig`). Complex division is Smith's (`ops.div`).

const std = @import("std");
const blas = @import("../blas/__root__.zig");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const auxiliary = @import("auxiliary.zig");
const qr_pivot = @import("qr_pivot.zig");
const norm = @import("norm.zig");
const qr = @import("qr.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = types.Float;
const isComplex = types.isComplex;
const Error = ops.Error;

// ---- xLASCL ------------------------------------------------------------------------

// ---- xLAIC1 ------------------------------------------------------------------------

/// `xLAIC1`'s `JOB`: the largest (`max`, 1) or smallest (`min`, 2)
/// singular value.
pub const Extreme = enum { max, min };

/// `?LAIC1` (`dlaic1.f`, `zlaic1.f`): one step of incremental condition
/// estimation. Given `sest`, an estimate of the extreme singular value of
/// the `j x j` triangle `L` with its vector `x` (`‖x‖ = 1`), and a new
/// column `w`, `gamma`, returns the estimate for `[L 0; wᴴ gamma]` and the
/// `s`, `c` making `[s x; c]` its vector. `zlaic1.f`'s complex formulas,
/// which give `dlaic1.f`'s bits for real `T`.
pub fn laic1(io: Io, comptime T: type, job: Extreme, j: usize, x: [*]const T, sest: Float(T), w: [*]const T, gamma: T) Error!struct { sestpr: Float(T), s: T, c: T } {
    const F = Float(T);
    const eps = std.math.floatEps(F) / 2;
    const alpha = try blas_calls.dot(io, T, true, j, x, 1, w, 1);
    const absalp = ops.abs(T, alpha);
    const absgam = ops.abs(T, gamma);
    const absest = @abs(sest);
    const one = ops.one(T);
    const zero = ops.zero(T);
    const cj = struct {
        fn f(v: T) T {
            return ops.conj(T, v);
        }
    }.f;
    const dr = struct {
        fn f(v: T, d: F) T {
            return ops.divReal(T, v, d);
        }
    }.f;
    switch (job) {
        .max => {
            if (sest == 0) {
                const s1 = @max(absgam, absalp);
                if (s1 == 0) return .{ .sestpr = 0, .s = zero, .c = one };
                var s = dr(alpha, s1);
                var c = dr(gamma, s1);
                const tmp = @sqrt(ops.abs2(T, s) + ops.abs2(T, c));
                s = dr(s, tmp);
                c = dr(c, tmp);
                return .{ .sestpr = s1 * tmp, .s = s, .c = c };
            }
            if (absgam <= eps * absest) {
                const tmp = @max(absest, absalp);
                const s1 = absest / tmp;
                const s2 = absalp / tmp;
                return .{ .sestpr = tmp * @sqrt(s1 * s1 + s2 * s2), .s = one, .c = zero };
            }
            if (absalp <= eps * absest) {
                return if (absgam <= absest) .{ .sestpr = absest, .s = one, .c = zero } else .{ .sestpr = absgam, .s = zero, .c = one };
            }
            if (absest <= eps * absalp or absest <= eps * absgam) {
                const s1 = absgam;
                const s2 = absalp;
                if (s1 <= s2) {
                    const tmp = s1 / s2;
                    const scl = @sqrt(1 + tmp * tmp);
                    return .{ .sestpr = s2 * scl, .s = dr(dr(alpha, s2), scl), .c = dr(dr(gamma, s2), scl) };
                }
                const tmp = s2 / s1;
                const scl = @sqrt(1 + tmp * tmp);
                return .{ .sestpr = s1 * scl, .s = dr(dr(alpha, s1), scl), .c = dr(dr(gamma, s1), scl) };
            }
            const zeta1 = absalp / absest;
            const zeta2 = absgam / absest;
            const b = (1 - zeta1 * zeta1 - zeta2 * zeta2) * 0.5;
            const c = zeta1 * zeta1;
            const t = if (b > 0) c / (b + @sqrt(b * b + c)) else @sqrt(b * b + c) - b;
            const sine = ops.neg(T, dr(dr(alpha, absest), t));
            const cosine = ops.neg(T, dr(dr(gamma, absest), 1 + t));
            const tmp = @sqrt(ops.abs2(T, sine) + ops.abs2(T, cosine));
            return .{ .sestpr = @sqrt(t + 1) * absest, .s = dr(sine, tmp), .c = dr(cosine, tmp) };
        },
        .min => {
            if (sest == 0) {
                var sine = one;
                var cosine = zero;
                if (@max(absgam, absalp) != 0) {
                    sine = ops.neg(T, cj(gamma));
                    cosine = cj(alpha);
                }
                const s1 = @max(ops.abs(T, sine), ops.abs(T, cosine));
                var s = dr(sine, s1);
                var c = dr(cosine, s1);
                const tmp = @sqrt(ops.abs2(T, s) + ops.abs2(T, c));
                s = dr(s, tmp);
                c = dr(c, tmp);
                return .{ .sestpr = 0, .s = s, .c = c };
            }
            if (absgam <= eps * absest) return .{ .sestpr = absgam, .s = zero, .c = one };
            if (absalp <= eps * absest) {
                return if (absgam <= absest) .{ .sestpr = absgam, .s = zero, .c = one } else .{ .sestpr = absest, .s = one, .c = zero };
            }
            if (absest <= eps * absalp or absest <= eps * absgam) {
                const s1 = absgam;
                const s2 = absalp;
                if (s1 <= s2) {
                    const tmp = s1 / s2;
                    const scl = @sqrt(1 + tmp * tmp);
                    return .{ .sestpr = absest * (tmp / scl), .s = ops.neg(T, dr(dr(cj(gamma), s2), scl)), .c = dr(dr(cj(alpha), s2), scl) };
                }
                const tmp = s2 / s1;
                const scl = @sqrt(1 + tmp * tmp);
                return .{ .sestpr = absest / scl, .s = ops.neg(T, dr(dr(cj(gamma), s1), scl)), .c = dr(dr(cj(alpha), s1), scl) };
            }
            const zeta1 = absalp / absest;
            const zeta2 = absgam / absest;
            const norma = @max(1 + zeta1 * zeta1 + zeta1 * zeta2, zeta1 * zeta2 + zeta2 * zeta2);
            const tst = 1 + 2 * (zeta1 - zeta2) * (zeta1 + zeta2);
            var sine: T = undefined;
            var cosine: T = undefined;
            var sestpr: F = undefined;
            if (tst >= 0) {
                const b = (zeta1 * zeta1 + zeta2 * zeta2 + 1) * 0.5;
                const c = zeta2 * zeta2;
                const t = c / (b + @sqrt(@abs(b * b - c)));
                sine = dr(dr(alpha, absest), 1 - t);
                cosine = ops.neg(T, dr(dr(gamma, absest), t));
                sestpr = @sqrt(t + 4 * eps * eps * norma) * absest;
            } else {
                const b = (zeta2 * zeta2 + zeta1 * zeta1 - 1) * 0.5;
                const c = zeta1 * zeta1;
                const t = if (b >= 0) -c / (b + @sqrt(b * b + c)) else b - @sqrt(b * b + c);
                sine = ops.neg(T, dr(dr(alpha, absest), t));
                cosine = ops.neg(T, dr(dr(gamma, absest), 1 + t));
                sestpr = @sqrt(1 + t + 4 * eps * eps * norma) * absest;
            }
            const tmp = @sqrt(ops.abs2(T, sine) + ops.abs2(T, cosine));
            return .{ .sestpr = sestpr, .s = dr(sine, tmp), .c = dr(cosine, tmp) };
        },
    }
}

// ---- xLARZ, xLATRZ, xORMR3 ---------------------------------------------------------

/// `?LARZ` (`dlarz.f`, `zlarz.f`): applies `H = I - tau v vᴴ`, `v = [1; 0;
/// v(l)]` (`v` with stride `incv`), to the `m x n` `C`: `H C` (`left`, the
/// last `l` rows with the first) or `C H`. `work` takes `n` (`m`) elements.
fn larz(al: Allocator, io: Io, comptime T: type, left: bool, m: usize, n: usize, l: usize, v: [*]const T, incv: usize, tau: T, c: Mat(T), work: [*]T) Error!void {
    if (ops.isZero(T, tau)) return;
    const one = ops.one(T);
    const ntau = ops.neg(T, tau);
    const vv: [*]T = @constCast(v);
    if (left) {
        // w := C(0, :)ᵀ + C(m-l.., :)ᵀ conj(v); C(0, :) -= tau w; C(m-l.., :) -= tau v wᵀ
        for (0..n) |j| work[j] = c.get(0, j);
        try blas_calls.gemv(al, io, T, .u, one, c.sub(m - l, 0, l, n), vv, incv, one, work, 1);
        try blas_calls.axpy(io, T, n, ntau, work, 1, c.ptr, c.ld);
        try blas_calls.ger(al, io, T, ntau, vv, incv, work, 1, c.sub(m - l, 0, l, n));
    } else {
        // w := C(:, 0) + C(:, n-l..) v; C(:, 0) -= tau w; C(:, n-l..) -= tau w vᴴ
        @memcpy(work[0..m], c.col(0)[0..m]);
        try blas_calls.gemv(al, io, T, .n, one, c.sub(0, n - l, m, l), vv, incv, one, work, 1);
        try blas_calls.axpy(io, T, m, ntau, work, 1, c.ptr, 1);
        if (comptime isComplex(T)) {
            if (m > 0 and l > 0) blas.gerc(al, io, T, @intCast(m), @intCast(l), ntau, work, 1, vv, @intCast(incv), c.ptr + (n - l) * c.ld, @intCast(c.ld)) catch |e| return ops.narrow(ops.Error, e);
        } else {
            try blas_calls.ger(al, io, T, ntau, work, 1, vv, incv, c.sub(0, n - l, m, l));
        }
    }
}

/// `?LATRZ` (`dlatrz.f`, `zlatrz.f`): `[A1 A2] = [R 0] Z` for the `m x n`
/// upper trapezoidal `A` (`A1` `m x m`, `A2` its last `n - m` columns), `Z`
/// as `m` reflections whose tails stay in `A2`'s rows, `tau` (`m`) out.
fn latrz(al: Allocator, io: Io, comptime T: type, a: Mat(T), tau: [*]T) Error!void {
    const m = a.rows;
    const n = a.cols;
    if (m == 0) return;
    if (m == n) {
        @memset(tau[0..n], ops.zero(T));
        return;
    }
    const l = n - m;
    const work = try al.alloc(T, m);
    defer al.free(work);
    var i = m;
    while (i > 0) {
        i -= 1;
        const row = a.ptr + i + (n - l) * a.ld;
        if (comptime isComplex(T)) {
            for (0..l) |k| row[k * a.ld] = ops.conj(T, row[k * a.ld]);
            var alpha = ops.conj(T, a.get(i, i));
            const t = try auxiliary.larfg(io, T, l + 1, &alpha, row, a.ld);
            tau[i] = ops.conj(T, t);
            try larz(al, io, T, false, i, n - i, l, row, a.ld, t, a.sub(0, i, i, n - i), work.ptr);
            a.at(i, i).* = ops.conj(T, alpha);
        } else {
            tau[i] = try auxiliary.larfg(io, T, l + 1, a.at(i, i), row, a.ld);
            try larz(al, io, T, false, i, n - i, l, row, a.ld, tau[i], a.sub(0, i, i, n - i), work.ptr);
        }
    }
}

/// `?ORMR3`, `?UNMR3` (`dormr3.f`, `zunmr3.f`) with side `L` and `trans`
/// `T` (`C`), as `xGELSY` calls them through `xORMRZ`: `C := Zᴴ C` for the
/// `m x n` `C` and `latrz`'s `k` reflections in `a`, `l` long.
fn ormr3(al: Allocator, io: Io, comptime T: type, k: usize, l: usize, a: Mat(T), tau: [*]const T, c: Mat(T)) Error!void {
    const m = c.rows;
    const n = c.cols;
    if (m == 0 or n == 0 or k == 0) return;
    const work = try al.alloc(T, n);
    defer al.free(work);
    const ja = m - l;
    for (0..k) |i| {
        const taui = ops.conj(T, tau[i]);
        try larz(al, io, T, true, m - i, n, l, a.ptr + i + ja * a.ld, a.ld, taui, c.sub(i, 0, m - i, n), work.ptr);
    }
}

// ---- xGELSY ------------------------------------------------------------------------

/// What `gelsy` returns.
pub const GelsyError = error{ InvalidM, InvalidN, InvalidNrhs, InvalidLda, InvalidLdb } || Allocator.Error || Io.Cancelable;

/// `?GELSY` (`dgelsy.f`, `zgelsy.f`): the minimum-norm solution of
/// `min ‖A X - B‖` for the `m x n` `A` of any rank, by the complete
/// orthogonal factorization. The effective rank, written to `rank`, is the
/// largest leading `R11` of `A P = Q R` whose condition estimate stays
/// below `1 / rcond`. `B` is `max(m, n) x nrhs`; `X` replaces its first `n`
/// rows. `jpvt` as `geqp3`'s: nonzero entries fix columns to the front; on
/// exit, `A P`'s columns. `A` is overwritten by the factorization.
pub fn gelsy(allocator: Allocator, io: Io, comptime T: type, m: isize, n: isize, nrhs: isize, a: [*]T, lda: isize, b: [*]T, ldb: isize, jpvt: [*]isize, rcond: Float(T), rank: *isize) GelsyError!void {
    comptime ops.requireFloat(T, "gelsy");
    const F = Float(T);
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (nrhs < 0) return error.InvalidNrhs;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (ldb < @max(1, m, n)) return error.InvalidLdb;
    const mn = @min(m, n);
    rank.* = 0;
    if (mn == 0 or nrhs == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(@max(m, n), nrhs, ldb, @sizeOf(T))) return error.InvalidLdb;
    const al = allocator;
    const um: usize = @intCast(m);
    const un: usize = @intCast(n);
    const umn: usize = @intCast(mn);
    const am: Mat(T) = .{ .ptr = a, .rows = um, .cols = un, .ld = @intCast(lda) };
    const bm: Mat(T) = .{ .ptr = b, .rows = @max(um, un), .cols = @intCast(nrhs), .ld = @intCast(ldb) };
    const smlnum = std.math.floatMin(F) / std.math.floatEps(F);
    const bignum = 1 / smlnum;
    // scale A and B into [smlnum, bignum]
    const anrm = try norm.matNorm(al, T, .m, um, un, a, am.ld);
    var iascl: F = 0;
    if (anrm > 0 and anrm < smlnum) {
        auxiliary.lascl(T, false, anrm, smlnum, am);
        iascl = smlnum;
    } else if (anrm > bignum) {
        auxiliary.lascl(T, false, anrm, bignum, am);
        iascl = bignum;
    } else if (anrm == 0) {
        bm.fill(ops.zero(T));
        return;
    }
    const brows = bm.sub(0, 0, um, bm.cols);
    const bnrm = try norm.matNorm(al, T, .m, um, bm.cols, b, bm.ld);
    var ibscl: F = 0;
    if (bnrm > 0 and bnrm < smlnum) {
        auxiliary.lascl(T, false, bnrm, smlnum, brows);
        ibscl = smlnum;
    } else if (bnrm > bignum) {
        auxiliary.lascl(T, false, bnrm, bignum, brows);
        ibscl = bignum;
    }
    const tau = try al.alloc(T, umn);
    defer al.free(tau);
    try qr_pivot.factor(al, io, T, am, jpvt, tau.ptr);
    // the rank, by incremental condition estimation on R
    const xmin = try al.alloc(T, umn);
    defer al.free(xmin);
    const xmax = try al.alloc(T, umn);
    defer al.free(xmax);
    xmin[0] = ops.one(T);
    xmax[0] = ops.one(T);
    var smax = ops.abs(T, a[0]);
    var smin = smax;
    if (smax == 0) {
        bm.fill(ops.zero(T));
        return;
    }
    var r: usize = 1;
    while (r < umn) {
        const col = am.col(r);
        const lo = try laic1(io, T, .min, r, xmin.ptr, smin, col, col[r]);
        const hi = try laic1(io, T, .max, r, xmax.ptr, smax, col, col[r]);
        if (!(hi.sestpr * rcond <= lo.sestpr)) break;
        for (0..r) |i| {
            xmin[i] = ops.mul(T, lo.s, xmin[i]);
            xmax[i] = ops.mul(T, hi.s, xmax[i]);
        }
        xmin[r] = lo.c;
        xmax[r] = hi.c;
        smin = lo.sestpr;
        smax = hi.sestpr;
        r += 1;
    }
    rank.* = @intCast(r);
    // [R11 R12] = [T11 0] Z
    const ztau = try al.alloc(T, r);
    defer al.free(ztau);
    if (r < un) try latrz(al, io, T, am.sub(0, 0, r, un), ztau.ptr);
    // B := Qᴴ B, then T11⁻¹ on its first `rank` rows
    try qr.multiply(al, io, T, .l, true, am.sub(0, 0, um, umn), tau.ptr, brows);
    try blas_calls.tri(al, io, T, true, true, .n, false, ops.one(T), am.sub(0, 0, r, r), bm.sub(0, 0, r, bm.cols));
    for (0..bm.cols) |j| @memset(bm.col(j)[r..un], ops.zero(T));
    // B := Zᴴ B, then P B
    if (r < un) try ormr3(al, io, T, r, un - r, am, ztau.ptr, bm.sub(0, 0, un, bm.cols));
    const work = try al.alloc(T, un);
    defer al.free(work);
    for (0..bm.cols) |j| {
        const col = bm.col(j);
        for (0..un) |i| work[@intCast(jpvt[i] - 1)] = col[i];
        @memcpy(col[0..un], work);
    }
    // undo the scaling
    const xrows = bm.sub(0, 0, un, bm.cols);
    if (iascl != 0) {
        auxiliary.lascl(T, false, anrm, iascl, xrows);
        auxiliary.lascl(T, true, iascl, anrm, am.sub(0, 0, r, r));
    }
    if (ibscl != 0) auxiliary.lascl(T, false, ibscl, bnrm, xrows);
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const lcheck = @import("_test/check.zig");
const test_io = @import("../blas/_test/io.zig");
const svd = @import("svd.zig");

const test_types = lcheck.test_types;

const thresh = lcheck.thresh;

/// A random `m x n` of rank `r` scaled by `s`, and `B` (`max(m, n) x nrhs`)
/// scaled by `s` too, so `X` does not move.
fn problem(comptime T: type, al: Allocator, m: usize, n: usize, r: usize, nrhs: usize, s: Float(T), seed: u64) !struct { a: Mat(T), b: Mat(T) } {
    const a = try ops.alloc(T, al, m, n);
    errdefer ops.free(T, al, a);
    try lcheck.rankProduct(T, al, a, r, seed);
    for (0..n) |j| for (0..m) |i| {
        a.at(i, j).* = ops.mulReal(T, a.get(i, j), s);
    };
    const b = try ops.alloc(T, al, @max(m, n), nrhs);
    check.fill(T, b.ptr[0 .. b.ld * nrhs], seed + 2);
    for (0..nrhs) |j| for (b.col(j)[0..b.rows]) |*v| {
        v.* = ops.mulReal(T, v.*, s);
    };
    return .{ .a = a, .b = b };
}

/// `gelsy` against the SVD's minimum-norm solution (`svd.zig`) cut at
/// rank `r`: the rank exactly, `X` within `thresh n eps cond(A) ‖X‖`.
fn gelsyCase(comptime T: type, al: Allocator, m: usize, n: usize, r: usize, nrhs: usize, s: Float(T), seed: u64) !void {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const p = try problem(T, al, m, n, r, nrhs, s, seed);
    defer ops.free(T, al, p.a);
    defer ops.free(T, al, p.b);
    // the reference solution from the unscaled SVD
    const k = @min(m, n);
    const sa = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, sa);
    for (0..n) |j| for (0..m) |i| {
        sa.at(i, j).* = ops.mulReal(T, p.a.get(i, j), 1 / s);
    };
    const sv = try al.alloc(F, k);
    defer al.free(sv);
    const u = try ops.alloc(T, al, m, m);
    defer ops.free(T, al, u);
    const vt = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, vt);
    try svd.gesvd(al, check.io, T, .a, .a, @intCast(m), @intCast(n), sa.ptr, @intCast(m), sv.ptr, u.ptr, @intCast(m), vt.ptr, @intCast(n), null);
    const want = try ops.alloc(T, al, n, nrhs);
    defer ops.free(T, al, want);
    want.fill(ops.zero(T));
    var wmax: F = 0;
    for (0..nrhs) |c| {
        for (0..r) |q| {
            // (u_qᴴ b / s_q) v_q, b unscaled
            var d = ops.zero(T);
            for (0..m) |i| d = ops.add(T, d, ops.mul(T, ops.conj(T, u.get(i, q)), ops.mulReal(T, p.b.get(i, c), 1 / s)));
            d = ops.divReal(T, d, sv[q]);
            for (0..n) |j| want.at(j, c).* = ops.add(T, want.get(j, c), ops.mul(T, d, ops.conj(T, vt.get(q, j))));
        }
        for (0..n) |j| wmax = @max(wmax, ops.abs(T, want.get(j, c)));
    }
    const jpvt = try al.alloc(isize, n);
    defer al.free(jpvt);
    @memset(jpvt, 0);
    var rank: isize = -1;
    const rcond = 1000 * eps;
    try gelsy(al, check.io, T, @intCast(m), @intCast(n), @intCast(nrhs), p.a.ptr, @intCast(p.a.ld), p.b.ptr, @intCast(p.b.ld), jpvt.ptr, rcond, &rank);
    try testing.expectEqual(@as(isize, @intCast(r)), rank);
    const cond = if (r > 0) sv[0] / sv[r - 1] else 1;
    const tol = thresh * @as(F, @floatFromInt(@max(m, n))) * eps * cond * @max(wmax, 1);
    for (0..nrhs) |c| for (0..n) |j| {
        const d = ops.abs(T, ops.sub(T, p.b.get(j, c), want.get(j, c)));
        if (!(d <= tol)) {
            check.say("gelsy {s} {d} x {d} rank {d} scale {e}: X({d}, {d}) off by {e}, tolerance {e}\n", .{ @typeName(T), m, n, r, @as(f64, @floatCast(s)), j, c, @as(f64, @floatCast(d)), @as(f64, @floatCast(tol)) });
            return error.TestUnexpectedResult;
        }
    };
}

test "lapack gelsy: the minimum-norm solution and rank, full and deficient rank, tall and wide, scaled near the range's ends" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const cases = [_][4]usize{ .{ 1, 1, 1, 1 }, .{ 6, 4, 4, 2 }, .{ 4, 6, 4, 1 }, .{ 8, 6, 3, 2 }, .{ 6, 8, 3, 3 }, .{ 20, 15, 10, 4 }, .{ 40, 40, 25, 2 } };
        for (cases, 0..) |c, i| {
            try gelsyCase(T, al, c[0], c[1], c[2], c[3], 1, i * 11 + 1);
        }
        // past smlnum and bignum, so A and B are scaled and back; normal
        // entries, so the rank survives. Not at an extended type, whose
        // `nrm2` is `nrm2.S`'s unscaled sum of squares and underflows there.
        if (comptime !types.isExtended(F)) {
            try gelsyCase(T, al, 8, 6, 3, 2, std.math.floatMin(F) / std.math.floatEps(F) / 64, 3);
            try gelsyCase(T, al, 6, 8, 4, 2, std.math.floatMax(F) / 4096, 5);
        }
    }
}

test "lapack gelsy: a zero A, a fixed column, empty, bad arguments" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        var a: [6]T = @splat(ops.zero(T));
        var b: [3]T = @splat(ops.one(T));
        var jpvt = [_]isize{ 0, 0 };
        var rank: isize = -1;
        try gelsy(al, check.io, T, 3, 2, 1, &a, 3, &b, 3, &jpvt, 0.01, &rank);
        try testing.expectEqual(@as(isize, 0), rank);
        for (b) |v| try testing.expect(ops.isZero(T, v));
        // column 2 fixed first; the solution of a square, well-conditioned A
        var c = [_]T{ ops.scalar(T, 4), ops.scalar(T, 1), ops.scalar(T, 1), ops.scalar(T, 3) };
        var d = [_]T{ ops.scalar(T, 5), ops.scalar(T, 4) };
        jpvt = .{ 0, 1 };
        try gelsy(al, check.io, T, 2, 2, 1, &c, 2, &d, 2, &jpvt, 0.01, &rank);
        try testing.expectEqual(@as(isize, 2), rank);
        try testing.expectEqualSlices(isize, &.{ 2, 1 }, &jpvt);
        const tol = 16 * std.math.floatEps(Float(T));
        try testing.expect(ops.abs(T, ops.sub(T, d[0], ops.one(T))) <= tol);
        try testing.expect(ops.abs(T, ops.sub(T, d[1], ops.one(T))) <= tol);
        rank = -1;
        try gelsy(al, check.io, T, 0, 2, 1, &a, 1, &b, 2, &jpvt, 0.01, &rank);
        try testing.expectEqual(@as(isize, 0), rank);
    }
    var a: [4]f64 = undefined;
    var jp: [2]isize = undefined;
    var rank: isize = 0;
    try testing.expectError(error.InvalidM, gelsy(al, check.io, f64, -1, 0, 0, &a, 1, &a, 1, &jp, 0, &rank));
    try testing.expectError(error.InvalidN, gelsy(al, check.io, f64, 1, -1, 0, &a, 1, &a, 1, &jp, 0, &rank));
    try testing.expectError(error.InvalidNrhs, gelsy(al, check.io, f64, 1, 1, -1, &a, 1, &a, 1, &jp, 0, &rank));
    try testing.expectError(error.InvalidLda, gelsy(al, check.io, f64, 2, 1, 1, &a, 1, &a, 2, &jp, 0, &rank));
    try testing.expectError(error.InvalidLdb, gelsy(al, check.io, f64, 1, 2, 1, &a, 1, &a, 1, &jp, 0, &rank));
}

test "lapack gelsy: one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        const m = 300;
        const n = 90;
        const nrhs = 70;
        const p = try problem(T, al, m, n, 60, nrhs, 1, 21);
        defer ops.free(T, al, p.a);
        defer ops.free(T, al, p.b);
        var outs: [3][]T = undefined;
        for (&outs) |*o| o.* = try al.alloc(T, m * n + m * nrhs);
        defer for (outs) |o| al.free(o);
        for (0..test_io.worker_runs.len) |c| {
            const io = runs.io(c);
            @memcpy(outs[c][0 .. m * n], p.a.ptr[0 .. m * n]);
            @memcpy(outs[c][m * n ..], p.b.ptr[0 .. m * nrhs]);
            var jpvt: [n]isize = @splat(0);
            var rank: isize = 0;
            try gelsy(al, io, T, m, n, nrhs, outs[c].ptr, m, outs[c].ptr + m * n, m, &jpvt, 1e-5, &rank);
        }
        try test_io.expectSameBytes(outs);
    }
}
