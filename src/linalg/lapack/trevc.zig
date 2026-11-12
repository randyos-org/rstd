//! Eigenvectors of a Schur form, `trevc3`, from Reference-LAPACK (BSD-3, as
//! OpenBLAS v0.3.34-317 ships it in `lapack-netlib/SRC/`): `dtrevc3.f` and
//! `ztrevc3.f`, the real 1 x 1 and 2 x 2 solves by `sylvester.laln2`. Faer's
//! back substitution (`evd/mod.rs`) has no guard against overflow; these
//! scale each vector as it grows (`SCALE`, `BIGNUM`), so a vector that would
//! overflow comes out scaled down instead.

const std = @import("std");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const auxiliary = @import("auxiliary.zig");
const norm = @import("norm.zig");
const sylvester = @import("sylvester.zig");
const util = @import("../blas/util.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;

/// `dtrevc3.f`'s and `ztrevc3.f`'s constants: `SMLNUM = UNFL (N / ULP)`,
/// `BIGNUM = (1 - ULP) / SMLNUM`.
fn Limits(comptime F: type) type {
    return struct {
        ulp: F,
        smlnum: F,
        bignum: F,

        fn of(n: usize) @This() {
            const ulp = std.math.floatEps(F);
            const smlnum = std.math.floatMin(F) * (@as(F, @floatFromInt(n)) / ulp);
            return .{ .ulp = ulp, .smlnum = smlnum, .bignum = (1 - ulp) / smlnum };
        }
    };
}

fn scal(comptime F: type, x: []F, s: F) void {
    for (x) |*v| v.* *= s;
}

/// `x[0..j] -= alpha t(0..j, c)`, `DAXPY`.
fn axpyCol(comptime F: type, t: Mat(F), c: usize, j: usize, alpha: F, x: [*]F) void {
    for (0..j) |k| x[k] -= alpha * t.get(k, c);
}

/// `DDOT(t(from..to, c), x[from..to])`.
fn dotCol(comptime F: type, t: Mat(F), c: usize, from: usize, to: usize, x: [*]const F) F {
    var s: F = 0;
    for (from..to) |k| s += t.get(k, c) * x[k];
    return s;
}

/// `x`'s rows `from..to` (and `y`'s) scaled so the largest `|x(k)| +
/// |y(k)|` (`|x(k)|` alone without `y`) is 1, as `HOWMNY = 'A'` leaves them.
fn normalize(comptime F: type, x: [*]F, y: ?[*]F, from: usize, to: usize) void {
    var emax: F = 0;
    for (from..to) |k| emax = @max(emax, @abs(x[k]) + (if (y) |v| @abs(v[k]) else 0));
    const remax = 1 / emax;
    scal(F, x[from..to], remax);
    if (y) |v| scal(F, v[from..to], remax);
}

/// `dtrevc3.f` with `SIDE = 'R'` or `'L'`, `HOWMNY = 'A'`: the right
/// (`T x = λ x`, upper triangular) or left (`uᴴ T = λ uᴴ`, lower
/// triangular) eigenvectors of the upper quasi-triangular `t` (Schur form,
/// standardized 2 x 2 blocks) into `x`, `n x n`. A complex pair takes two
/// columns, the real and imaginary parts, the one of positive imaginary part
/// `λ` first. `cn` (`n`) is scratch.
pub fn real(comptime F: type, right: bool, t: Mat(F), x: Mat(F), cn: []F) void {
    const n = t.rows;
    const lim = Limits(F).of(n);
    x.fill(0);
    // the strictly upper part's column 1-norms, against overflow
    for (0..n) |j| cn[j] = colNorm(F, t, j);
    if (right) realRight(F, t, x, cn, lim) else realLeft(F, t, x, cn, lim);
}

fn colNorm(comptime F: type, t: Mat(F), j: usize) F {
    var s: F = 0;
    for (0..j) |i| s += @abs(t.get(i, j));
    return s;
}

fn realRight(comptime F: type, t: Mat(F), x: Mat(F), cn: []const F, lim: Limits(F)) void {
    const n = t.rows;
    var ki = n;
    while (ki > 0) {
        ki -= 1;
        const wr = t.get(ki, ki);
        if (ki == 0 or t.get(ki, ki - 1) == 0) {
            const smin = @max(lim.ulp * @abs(wr), lim.smlnum);
            const v = x.col(ki);
            v[ki] = 1;
            for (0..ki) |k| v[k] = -t.get(k, ki);
            // [ T(0:ki, 0:ki) - wr ] x = scale v, bottom up
            var j = ki;
            while (j > 0) {
                j -= 1;
                if (j == 0 or t.get(j, j - 1) == 0) {
                    var r = sylvester.laln2(F, false, 1, 1, smin, t.col(j) + j, t.ld, .{ .{ v[j], 0 }, .{ 0, 0 } }, wr, 0);
                    // keep the update of the rest from overflowing
                    if (r.xnorm > 1 and cn[j] > lim.bignum / r.xnorm) {
                        r.x[0][0] /= r.xnorm;
                        r.scale /= r.xnorm;
                    }
                    if (r.scale != 1) scal(F, v[0 .. ki + 1], r.scale);
                    v[j] = r.x[0][0];
                    axpyCol(F, t, j, j, r.x[0][0], v);
                } else {
                    var r = sylvester.laln2(F, false, 2, 1, smin, t.col(j - 1) + j - 1, t.ld, .{ .{ v[j - 1], 0 }, .{ v[j], 0 } }, wr, 0);
                    if (r.xnorm > 1) {
                        const beta = @max(cn[j - 1], cn[j]);
                        if (beta > lim.bignum / r.xnorm) {
                            r.x[0][0] /= r.xnorm;
                            r.x[1][0] /= r.xnorm;
                            r.scale /= r.xnorm;
                        }
                    }
                    if (r.scale != 1) scal(F, v[0 .. ki + 1], r.scale);
                    v[j - 1] = r.x[0][0];
                    v[j] = r.x[1][0];
                    axpyCol(F, t, j - 1, j - 1, r.x[0][0], v);
                    axpyCol(F, t, j, j - 1, r.x[1][0], v);
                    j -= 1;
                }
            }
            normalize(F, v, null, 0, ki + 1);
            continue;
        }
        // the pair (ki - 1, ki): λ = wr + i wi
        const wi = @sqrt(@abs(t.get(ki, ki - 1))) * @sqrt(@abs(t.get(ki - 1, ki)));
        const smin = @max(lim.ulp * (@abs(wr) + @abs(wi)), lim.smlnum);
        const xr = x.col(ki - 1);
        const xi = x.col(ki);
        if (@abs(t.get(ki - 1, ki)) >= @abs(t.get(ki, ki - 1))) {
            xr[ki - 1] = 1;
            xi[ki] = wi / t.get(ki - 1, ki);
        } else {
            xr[ki - 1] = -wi / t.get(ki, ki - 1);
            xi[ki] = 1;
        }
        xr[ki] = 0;
        xi[ki - 1] = 0;
        for (0..ki - 1) |k| {
            xr[k] = -xr[ki - 1] * t.get(k, ki - 1);
            xi[k] = -xi[ki] * t.get(k, ki);
        }
        var j = ki - 1;
        while (j > 0) {
            j -= 1;
            if (j == 0 or t.get(j, j - 1) == 0) {
                var r = sylvester.laln2(F, false, 1, 2, smin, t.col(j) + j, t.ld, .{ .{ xr[j], xi[j] }, .{ 0, 0 } }, wr, wi);
                if (r.xnorm > 1 and cn[j] > lim.bignum / r.xnorm) {
                    r.x[0][0] /= r.xnorm;
                    r.x[0][1] /= r.xnorm;
                    r.scale /= r.xnorm;
                }
                if (r.scale != 1) {
                    scal(F, xr[0 .. ki + 1], r.scale);
                    scal(F, xi[0 .. ki + 1], r.scale);
                }
                xr[j] = r.x[0][0];
                xi[j] = r.x[0][1];
                axpyCol(F, t, j, j, r.x[0][0], xr);
                axpyCol(F, t, j, j, r.x[0][1], xi);
            } else {
                var r = sylvester.laln2(F, false, 2, 2, smin, t.col(j - 1) + j - 1, t.ld, .{ .{ xr[j - 1], xi[j - 1] }, .{ xr[j], xi[j] } }, wr, wi);
                if (r.xnorm > 1) {
                    const beta = @max(cn[j - 1], cn[j]);
                    if (beta > lim.bignum / r.xnorm) {
                        const rec = 1 / r.xnorm;
                        for (&r.x) |*row| for (row) |*e| {
                            e.* *= rec;
                        };
                        r.scale *= rec;
                    }
                }
                if (r.scale != 1) {
                    scal(F, xr[0 .. ki + 1], r.scale);
                    scal(F, xi[0 .. ki + 1], r.scale);
                }
                xr[j - 1] = r.x[0][0];
                xr[j] = r.x[1][0];
                xi[j - 1] = r.x[0][1];
                xi[j] = r.x[1][1];
                axpyCol(F, t, j - 1, j - 1, r.x[0][0], xr);
                axpyCol(F, t, j, j - 1, r.x[1][0], xr);
                axpyCol(F, t, j - 1, j - 1, r.x[0][1], xi);
                axpyCol(F, t, j, j - 1, r.x[1][1], xi);
                j -= 1;
            }
        }
        normalize(F, xr, xi, 0, ki + 1);
        ki -= 1;
    }
}

/// The left solves' guard (`VMAX`, `VCRIT`): before a right-hand side whose
/// column norm `beta` could overflow it, the vector is scaled by
/// `1 / VMAX`.
fn Growth(comptime F: type) type {
    return struct {
        vmax: F = 1,
        vcrit: F,

        fn guard(g: *@This(), beta: F, bignum: F, vs: []const []F) void {
            if (beta <= g.vcrit) return;
            const rec = 1 / g.vmax;
            for (vs) |v| scal(F, v, rec);
            g.vmax = 1;
            g.vcrit = bignum;
        }

        fn grew(g: *@This(), m: F, bignum: F) void {
            g.vmax = @max(m, g.vmax);
            g.vcrit = bignum / g.vmax;
        }
    };
}

fn realLeft(comptime F: type, t: Mat(F), x: Mat(F), cn: []const F, lim: Limits(F)) void {
    const n = t.rows;
    var ki: usize = 0;
    while (ki < n) {
        const wr = t.get(ki, ki);
        if (ki + 1 == n or t.get(ki + 1, ki) == 0) {
            const smin = @max(lim.ulp * @abs(wr), lim.smlnum);
            const v = x.col(ki);
            v[ki] = 1;
            for (ki + 1..n) |k| v[k] = -t.get(ki, k);
            // [ T(ki+1:n, ki+1:n) - wr ]ᵀ x = scale v, top down
            var g: Growth(F) = .{ .vcrit = lim.bignum };
            var j = ki + 1;
            while (j < n) {
                if (j + 1 == n or t.get(j + 1, j) == 0) {
                    g.guard(cn[j], lim.bignum, &.{v[ki..n]});
                    v[j] -= dotCol(F, t, j, ki + 1, j, v);
                    const r = sylvester.laln2(F, false, 1, 1, smin, t.col(j) + j, t.ld, .{ .{ v[j], 0 }, .{ 0, 0 } }, wr, 0);
                    if (r.scale != 1) scal(F, v[ki..n], r.scale);
                    v[j] = r.x[0][0];
                    g.grew(@abs(v[j]), lim.bignum);
                    j += 1;
                } else {
                    g.guard(@max(cn[j], cn[j + 1]), lim.bignum, &.{v[ki..n]});
                    v[j] -= dotCol(F, t, j, ki + 1, j, v);
                    v[j + 1] -= dotCol(F, t, j + 1, ki + 1, j, v);
                    const r = sylvester.laln2(F, true, 2, 1, smin, t.col(j) + j, t.ld, .{ .{ v[j], 0 }, .{ v[j + 1], 0 } }, wr, 0);
                    if (r.scale != 1) scal(F, v[ki..n], r.scale);
                    v[j] = r.x[0][0];
                    v[j + 1] = r.x[1][0];
                    g.grew(@max(@abs(v[j]), @abs(v[j + 1])), lim.bignum);
                    j += 2;
                }
            }
            normalize(F, v, null, ki, n);
            ki += 1;
            continue;
        }
        // the pair (ki, ki + 1): uᴴ T = (wr + i wi) uᴴ
        const wi = @sqrt(@abs(t.get(ki, ki + 1))) * @sqrt(@abs(t.get(ki + 1, ki)));
        const smin = @max(lim.ulp * (@abs(wr) + @abs(wi)), lim.smlnum);
        const xr = x.col(ki);
        const xi = x.col(ki + 1);
        if (@abs(t.get(ki, ki + 1)) >= @abs(t.get(ki + 1, ki))) {
            xr[ki] = wi / t.get(ki, ki + 1);
            xi[ki + 1] = 1;
        } else {
            xr[ki] = 1;
            xi[ki + 1] = -wi / t.get(ki + 1, ki);
        }
        xr[ki + 1] = 0;
        xi[ki] = 0;
        for (ki + 2..n) |k| {
            xr[k] = -xr[ki] * t.get(ki, k);
            xi[k] = -xi[ki + 1] * t.get(ki + 1, k);
        }
        var g: Growth(F) = .{ .vcrit = lim.bignum };
        var j = ki + 2;
        while (j < n) {
            if (j + 1 == n or t.get(j + 1, j) == 0) {
                g.guard(cn[j], lim.bignum, &.{ xr[ki..n], xi[ki..n] });
                xr[j] -= dotCol(F, t, j, ki + 2, j, xr);
                xi[j] -= dotCol(F, t, j, ki + 2, j, xi);
                const r = sylvester.laln2(F, false, 1, 2, smin, t.col(j) + j, t.ld, .{ .{ xr[j], xi[j] }, .{ 0, 0 } }, wr, -wi);
                if (r.scale != 1) {
                    scal(F, xr[ki..n], r.scale);
                    scal(F, xi[ki..n], r.scale);
                }
                xr[j] = r.x[0][0];
                xi[j] = r.x[0][1];
                g.grew(@max(@abs(xr[j]), @abs(xi[j])), lim.bignum);
                j += 1;
            } else {
                g.guard(@max(cn[j], cn[j + 1]), lim.bignum, &.{ xr[ki..n], xi[ki..n] });
                xr[j] -= dotCol(F, t, j, ki + 2, j, xr);
                xi[j] -= dotCol(F, t, j, ki + 2, j, xi);
                xr[j + 1] -= dotCol(F, t, j + 1, ki + 2, j, xr);
                xi[j + 1] -= dotCol(F, t, j + 1, ki + 2, j, xi);
                const r = sylvester.laln2(F, true, 2, 2, smin, t.col(j) + j, t.ld, .{ .{ xr[j], xi[j] }, .{ xr[j + 1], xi[j + 1] } }, wr, -wi);
                if (r.scale != 1) {
                    scal(F, xr[ki..n], r.scale);
                    scal(F, xi[ki..n], r.scale);
                }
                xr[j] = r.x[0][0];
                xi[j] = r.x[0][1];
                xr[j + 1] = r.x[1][0];
                xi[j + 1] = r.x[1][1];
                g.grew(@max(@max(@abs(r.x[0][0]), @abs(r.x[0][1])), @max(@abs(r.x[1][0]), @abs(r.x[1][1]))), lim.bignum);
                j += 2;
            }
        }
        normalize(F, xr, xi, ki, n);
        ki += 2;
    }
}

/// `ztrevc3.f` with `SIDE = 'R'` or `'L'`, `HOWMNY = 'A'`: the right or left
/// eigenvectors of the upper triangular `t` into `x` (`n x n`), each a
/// shifted triangular solve by `latrs` (`ZLATRS`) with its diagonal raised
/// to `SMIN` where smaller. `t`'s diagonal is shifted in place and
/// restored. `rwork` (`n`) and `diag` (`n`) are scratch.
///
/// The left solve passes `ZLATRS` the norms of `T`'s columns `ki+1..n`,
/// bounds on its submatrix's; `ztrevc3.f` passes columns `1..n-ki`, which
/// belong to other columns (issue #261).
pub fn complex(al: Allocator, io: Io, comptime T: type, right: bool, t: Mat(T), x: Mat(T), rwork: []Float(T), diag: []T) ops.Error!void {
    const F = Float(T);
    const n = t.rows;
    const lim = Limits(F).of(n);
    x.fill(ops.zero(T));
    for (0..n) |i| diag[i] = t.get(i, i);
    for (0..n) |j| {
        rwork[j] = 0;
        for (0..j) |i| rwork[j] += ops.abs1(T, t.get(i, j));
    }
    for (0..n) |step| {
        const ki = if (right) n - 1 - step else step;
        const lam = diag[ki];
        const smin = @max(lim.ulp * ops.abs1(T, lam), lim.smlnum);
        const v = x.col(ki);
        v[ki] = ops.one(T);
        // the rows solved: above `ki` (right) or below (left)
        const lo = if (right) 0 else ki + 1;
        const hi = if (right) ki else n;
        for (lo..hi) |k| {
            v[k] = if (right) ops.neg(T, t.get(k, ki)) else ops.neg(T, ops.conj(T, t.get(ki, k)));
            const d = ops.sub(T, diag[k], lam);
            t.at(k, k).* = if (ops.abs1(T, d) < smin) ops.scalar(T, smin) else d;
        }
        if (hi > lo) {
            const scale = try norm.latrs(al, io, T, true, if (right) .n else .c, false, true, hi - lo, t.col(lo) + lo, t.ld, v + lo, rwork.ptr + lo);
            v[ki] = ops.scalar(T, scale);
        }
        for (lo..hi) |k| t.at(k, k).* = diag[k];
        // largest `|re| + |im|` to 1 (`IZAMAX`, `ZDSCAL`)
        const from = if (right) 0 else ki;
        const to = if (right) ki + 1 else n;
        var big: F = -1;
        for (from..to) |k| big = @max(big, ops.abs1(T, v[k]));
        const remax = 1 / big;
        for (v[from..to]) |*e| e.* = ops.mulReal(T, e.*, remax);
    }
}

// ---- trevc3, every SIDE and HOWMNY ---------------------------------------------------

/// `trevc3`'s vectors: right (`r`), left (`l`) or both (`b`).
pub const Side = enum { r, l, b };
/// `trevc3`'s `HOWMNY`: every vector of `T` (`a`), every one times the
/// `Q` given in `vl`/`vr` (`b`), or those `select` picks (`s`).
pub const Howmny = enum { a, b, s };

/// What `trevc3` returns: `InvalidSelect` for `s` without `select`,
/// `InvalidMm` when `mm` is short of the `m` columns needed.
pub const Trevc3Error = error{ InvalidN, InvalidSelect, InvalidLdt, InvalidLdvl, InvalidLdvr, InvalidMm } || Allocator.Error || Io.Cancelable;

/// `?TREVC3` (`dtrevc3.f`, `ztrevc3.f`): eigenvectors of the `n x n` Schur
/// form `T` (real: upper quasi-triangular with standardized 2 x 2 blocks, as
/// `hseqr` leaves it), right (`T x = λ x`) into `vr`, left (`uᴴ T = λ uᴴ`)
/// into `vl`, each `n x mm`. A real complex pair takes two columns, its
/// real and imaginary parts, the positive-imaginary `λ` first. Each vector
/// is scaled so its largest entry (`|re| + |im|`) is 1, and with `b` after
/// the multiplication by the `Q` given in `vl`/`vr` (from `orghr` and
/// `hseqr`, so the vectors are `A`'s). With `s`, only the vectors `select`
/// picks, in order; a real pair picked by either entry takes both columns,
/// `select` then marking its first only, as the reference does. `m` gets
/// the columns used. `T` is restored on return (the complex solves scale
/// its diagonal in place).
pub fn trevc3(allocator: Allocator, io: Io, comptime T: type, side: Side, howmny: Howmny, select: ?[*]bool, n: isize, t: [*]T, ldt: isize, vl: [*]T, ldvl: isize, vr: [*]T, ldvr: isize, mm: isize, m: *isize) Trevc3Error!void {
    comptime ops.requireFloat(T, "trevc3");
    const left = side != .r;
    const right = side != .l;
    if (n < 0) return error.InvalidN;
    if (howmny == .s and select == null) return error.InvalidSelect;
    if (ldt < @max(1, n)) return error.InvalidLdt;
    if (ldvl < 1 or (left and ldvl < n)) return error.InvalidLdvl;
    if (ldvr < 1 or (right and ldvr < n)) return error.InvalidLdvr;
    const un: usize = @intCast(n);
    if (un > 0 and !util.fitsMatrix(n, n, ldt, @sizeOf(T))) return error.InvalidLdt;
    const tm: Mat(T) = .{ .ptr = t, .rows = un, .cols = un, .ld = if (un > 0) @intCast(ldt) else 1 };
    var count = un;
    if (howmny == .s) {
        const sel = select.?;
        count = 0;
        var j: usize = 0;
        while (j < un) : (j += 1) {
            if (pairAt(T, tm, j)) {
                if (sel[j] or sel[j + 1]) {
                    sel[j] = true;
                    sel[j + 1] = false;
                    count += 2;
                }
                j += 1;
            } else if (sel[j]) count += 1;
        }
    }
    m.* = @intCast(count);
    if (mm < count) return error.InvalidMm;
    if (un == 0 or count == 0) return;
    if (left and !util.fitsMatrix(n, mm, ldvl, @sizeOf(T))) return error.InvalidLdvl;
    if (right and !util.fitsMatrix(n, mm, ldvr, @sizeOf(T))) return error.InvalidLdvr;
    const x = try ops.alloc(T, allocator, un, un);
    defer ops.free(T, allocator, x);
    const cn = try allocator.alloc(Float(T), un);
    defer allocator.free(cn);
    const diag = try allocator.alloc(T, if (comptime ops.isComplex(T)) un else 0);
    defer allocator.free(diag);
    for ([_]bool{ true, false }) |r| {
        if ((r and !right) or (!r and !left)) continue;
        const dst: Mat(T) = if (r)
            .{ .ptr = vr, .rows = un, .cols = @intCast(mm), .ld = @intCast(ldvr) }
        else
            .{ .ptr = vl, .rows = un, .cols = @intCast(mm), .ld = @intCast(ldvl) };
        if (comptime ops.isComplex(T)) try complex(allocator, io, T, r, tm, x, cn, diag) else real(T, r, tm, x, cn);
        switch (howmny) {
            .a => dst.sub(0, 0, un, un).copyFrom(x),
            .b => {
                const q = try ops.alloc(T, allocator, un, un);
                defer ops.free(T, allocator, q);
                q.copyFrom(dst.sub(0, 0, un, un));
                try blas_calls.gemm(allocator, io, T, .n, .n, ops.one(T), q, x, ops.zero(T), dst.sub(0, 0, un, un));
                var j: usize = 0;
                while (j < un) : (j += 1) {
                    if (comptime ops.isComplex(T)) {
                        normalizeColumn(T, dst.col(j)[0..un]);
                    } else if (pairAt(T, tm, j)) {
                        normalize(T, dst.col(j), dst.col(j + 1), 0, un);
                        j += 1;
                    } else normalizeColumn(T, dst.col(j)[0..un]);
                }
            },
            .s => {
                const sel = select.?;
                var c: usize = 0;
                var j: usize = 0;
                while (j < un) : (j += 1) {
                    const two = pairAt(T, tm, j);
                    if (sel[j]) {
                        @memcpy(dst.col(c)[0..un], x.col(j)[0..un]);
                        c += 1;
                        if (two) {
                            @memcpy(dst.col(c)[0..un], x.col(j + 1)[0..un]);
                            c += 1;
                        }
                    }
                    if (two) j += 1;
                }
            },
        }
    }
}

/// A real 2 x 2 block (a complex pair) starts at `j`.
fn pairAt(comptime T: type, t: Mat(T), j: usize) bool {
    if (comptime ops.isComplex(T)) return false;
    return j + 1 < t.rows and t.get(j + 1, j) != 0;
}

/// Largest `|re| + |im|` (`|x|` for real) to 1.
fn normalizeColumn(comptime T: type, v: []T) void {
    var big: Float(T) = 0;
    for (v) |e| big = @max(big, ops.abs1(T, e));
    if (big == 0) return;
    const remax = 1 / big;
    for (v) |*e| e.* = ops.mulReal(T, e.*, remax);
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const types = @import("../blas/types.zig");
const evd = @import("evd.zig");
const hessenberg = @import("hessenberg.zig");

/// `‖A v - λ v‖` (right) or `‖vᴴ A - λ vᴴ‖` (left) for column `j` of `v`
/// (a real pair's two columns as one complex vector), and its largest
/// `|re| + |im|`.
fn eigResidual(comptime T: type, n: usize, a: []const T, v: []const T, j: usize, pair: bool, lam: std.math.Complex(Float(T)), right: bool) [2]Float(T) {
    const F = Float(T);
    const C = std.math.Complex(F);
    var worst: F = 0;
    var big: F = 0;
    for (0..n) |k| {
        const x: C = if (comptime types.isComplex(T)) v[k + j * n] else .init(v[k + j * n], if (pair) v[k + (j + 1) * n] else 0);
        big = @max(big, @abs(x.re) + @abs(x.im));
    }
    for (0..n) |i| {
        var s: C = .init(0, 0);
        for (0..n) |k| {
            const e = if (right) a[i + k * n] else a[k + i * n];
            const aik: C = if (comptime types.isComplex(T)) e else .init(e, 0);
            const xk_: C = if (comptime types.isComplex(T)) v[k + j * n] else .init(v[k + j * n], if (pair) v[k + (j + 1) * n] else 0);
            s = s.add(aik.mul(if (right) xk_ else xk_.conjugate()));
        }
        const xi_: C = if (comptime types.isComplex(T)) v[i + j * n] else .init(v[i + j * n], if (pair) v[i + (j + 1) * n] else 0);
        const r = s.sub(if (right) lam.mul(xi_) else lam.mul(xi_.conjugate()));
        worst = @max(worst, @sqrt(r.re * r.re + r.im * r.im));
    }
    return .{ worst, big };
}

test "trevc3: A's eigenvectors after gehrd, orghr and hseqr, back-transformed and selected, real and complex" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        const C = std.math.Complex(F);
        for ([_]usize{ 1, 2, 5, 12, 31 }) |n| {
            const ni: isize = @intCast(n);
            const a0 = try al.alloc(T, n * n);
            defer al.free(a0);
            check.fill(T, a0, 13 * n + 1);
            const t = try al.dupe(T, a0);
            defer al.free(t);
            const tau = try al.alloc(T, @max(n, 2) - 1);
            defer al.free(tau);
            try hessenberg.gehrd(al, check.io, T, ni, 1, ni, t.ptr, ni, tau.ptr);
            const z = try al.dupe(T, t);
            defer al.free(z);
            try hessenberg.orghr(al, check.io, T, ni, 1, ni, z.ptr, ni, tau.ptr);
            const ev = try al.alloc(F, 2 * n);
            defer al.free(ev);
            const cw = try al.alloc(T, n);
            defer al.free(cw);
            const w: evd.Eigenvalues(T) = if (comptime types.isComplex(T)) cw.ptr else .{ .re = ev.ptr, .im = ev.ptr + n };
            try evd.hseqr(check.io, T, .s, .v, ni, 1, ni, t.ptr, ni, w, z.ptr, ni, null);
            const vl = try al.dupe(T, z);
            defer al.free(vl);
            const vr = try al.dupe(T, z);
            defer al.free(vr);
            var m: isize = 0;
            try trevc3(al, check.io, T, .b, .b, null, ni, t.ptr, ni, vl.ptr, ni, vr.ptr, ni, ni, &m);
            try testing.expectEqual(ni, m);
            const tol = 256 * std.math.floatEps(F) * @as(F, @floatFromInt(n + 1)) * @max(auxiliary.normMax(T, .{ .ptr = a0.ptr, .rows = n, .cols = n, .ld = n }), 1);
            var j: usize = 0;
            while (j < n) : (j += 1) {
                const pair = !types.isComplex(T) and ev[n + j] > 0;
                const lam: C = if (comptime types.isComplex(T)) cw[j] else .init(ev[j], ev[n + j]);
                for ([_]bool{ true, false }) |right| {
                    const r = eigResidual(T, n, a0, if (right) vr else vl, j, pair, lam, right);
                    if (!(r[0] <= tol) or !(@abs(r[1] - 1) <= tol)) {
                        check.say("trevc3 {s} n={d} vector {d} {s}: residual {e}, largest {e}\n", .{ @typeName(T), n, j, if (right) "right" else "left", r[0], r[1] });
                        return error.TestUnexpectedResult;
                    }
                }
                if (pair) j += 1;
            }
            // selected: every other eigenvalue, a pair picked by its second entry
            const sel = try al.alloc(bool, n);
            defer al.free(sel);
            for (sel, 0..) |*s, i| s.* = i % 2 == 1;
            const all = try al.alloc(T, n * n);
            defer al.free(all);
            const some = try al.alloc(T, n * n);
            defer al.free(some);
            try trevc3(al, check.io, T, .r, .a, null, ni, t.ptr, ni, vl.ptr, ni, all.ptr, ni, ni, &m);
            try trevc3(al, check.io, T, .r, .s, sel.ptr, ni, t.ptr, ni, vl.ptr, ni, some.ptr, ni, ni, &m);
            var c: usize = 0;
            j = 0;
            while (j < n) : (j += 1) {
                const pair = !types.isComplex(T) and ev[n + j] > 0;
                if (sel[j]) {
                    for (0..if (pair) 2 else 1) |o| {
                        try testing.expectEqualSlices(T, all[(j + o) * n ..][0..n], some[(c + o) * n ..][0..n]);
                    }
                    c += if (pair) 2 else 1;
                }
                if (pair) {
                    try testing.expect(!sel[j + 1]);
                    j += 1;
                }
            }
            try testing.expectEqual(@as(isize, @intCast(c)), m);
        }
    }
    var a: [4]f64 = @splat(0);
    var m: isize = 0;
    try testing.expectError(error.InvalidN, trevc3(al, check.io, f64, .r, .a, null, -1, &a, 1, &a, 1, &a, 1, 0, &m));
    try testing.expectError(error.InvalidSelect, trevc3(al, check.io, f64, .r, .s, null, 2, &a, 2, &a, 1, &a, 2, 2, &m));
    try testing.expectError(error.InvalidLdvr, trevc3(al, check.io, f64, .r, .a, null, 2, &a, 2, &a, 1, &a, 1, 2, &m));
    try testing.expectError(error.InvalidMm, trevc3(al, check.io, f64, .r, .a, null, 2, &a, 2, &a, 1, &a, 2, 1, &m));
}
