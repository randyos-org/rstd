//! Eigenvectors of a Schur form for `geev`, from Reference-LAPACK (BSD-3, as
//! OpenBLAS v0.3.34-317 ships it in `lapack-netlib/SRC/`): `dtrevc3.f` and
//! `ztrevc3.f` (`HOWMNY = 'A'`, every vector of `T`, each scaled to largest
//! entry 1), with `dlaln2.f` for the real 1 x 1 and 2 x 2 solves. Faer's back
//! substitution (`evd/mod.rs`) has no guard against overflow; these scale
//! each vector as it grows (`SCALE`, `BIGNUM`), so a vector that would
//! overflow comes out scaled down instead. `evd.zig` multiplies by the Schur
//! vectors after, as the reference's blocked back-transform does.
//!
//! `DLALN2` is ported for `CA = D1 = D2 = 1`, all its callers here and in
//! `sylvester.zig` pass; complex division is Smith's (`util.cdiv`) where the
//! reference takes `DLADIV`.

const std = @import("std");
const ops = @import("ops.zig");
const norm = @import("norm.zig");
const util = @import("../blas/util.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;

/// What `laln2` returns: `x[i][k]` the solution's row `i`, real part
/// (`k = 0`) or imaginary (`k = 1`); `scale <= 1`; `xnorm` its infinity
/// norm; `perturbed` when a pivot was raised to `smin` (`INFO = 1`).
pub fn Laln2(comptime F: type) type {
    return struct { x: [2][2]F, scale: F, xnorm: F, perturbed: bool };
}

/// `DLALN2`: `(A - (wr + i wi)) X = scale B` (`Aᵀ` with `ltrans`) for the
/// `na x na` `A` at `a` (`na` 1 or 2), `X` and `B` real (`nw = 1`) or
/// complex (`nw = 2`, `b[i][1]` the imaginary parts), by complete pivoting;
/// pivots below `smin` become `smin`.
pub fn laln2(comptime F: type, ltrans: bool, na: usize, nw: usize, smin: F, a: [*]const F, lda: usize, b: [2][2]F, wr: F, wi: F) Laln2(F) {
    const ipivot = [4][4]usize{ .{ 0, 1, 2, 3 }, .{ 1, 0, 3, 2 }, .{ 2, 3, 0, 1 }, .{ 3, 2, 1, 0 } };
    const rswap = [4]bool{ false, true, false, true };
    const zswap = [4]bool{ false, false, true, true };
    const smlnum = 2 * std.math.floatMin(F);
    const bignum = 1 / smlnum;
    const smini = @max(smin, smlnum);
    var out: Laln2(F) = .{ .x = @splat(@splat(0)), .scale = 1, .xnorm = 0, .perturbed = false };
    if (na == 1) {
        var csr = a[0] - wr;
        var csi: F = if (nw == 1) 0 else -wi;
        var cnorm = @abs(csr) + @abs(csi);
        if (cnorm < smini) {
            csr = smini;
            csi = 0;
            cnorm = smini;
            out.perturbed = true;
        }
        const bnorm = @abs(b[0][0]) + (if (nw == 1) 0 else @abs(b[0][1]));
        if (cnorm < 1 and bnorm > 1 and bnorm > bignum * cnorm) out.scale = 1 / bnorm;
        if (nw == 1) {
            out.x[0][0] = (b[0][0] * out.scale) / csr;
            out.xnorm = @abs(out.x[0][0]);
        } else {
            const z = util.cdiv(F, out.scale * b[0][0], out.scale * b[0][1], csr, csi);
            out.x[0] = z;
            out.xnorm = @abs(z[0]) + @abs(z[1]);
        }
        return out;
    }
    // `CR`, `CI` column-major as `CRV`, `CIV`
    var crv: [4]F = .{ a[0] - wr, a[1], a[lda], a[1 + lda] - wr };
    if (ltrans) std.mem.swap(F, &crv[1], &crv[2]);
    const civ: [4]F = if (nw == 1) @splat(0) else .{ -wi, 0, 0, -wi };
    var cmax: F = 0;
    var icmax: usize = 0;
    for (0..4) |j| if (@abs(crv[j]) + @abs(civ[j]) > cmax) {
        cmax = @abs(crv[j]) + @abs(civ[j]);
        icmax = j;
    };
    if (cmax < smini) {
        const bnorm = @max(@abs(b[0][0]) + @abs(b[0][1]), @abs(b[1][0]) + @abs(b[1][1]));
        if (smini < 1 and bnorm > 1 and bnorm > bignum * smini) out.scale = 1 / bnorm;
        const temp = out.scale / smini;
        for (0..2) |i| for (0..nw) |k| {
            out.x[i][k] = temp * b[i][k];
        };
        out.xnorm = temp * bnorm;
        out.perturbed = true;
        return out;
    }
    const piv = ipivot[icmax];
    if (nw == 1) {
        const ur11 = crv[icmax];
        const cr21 = crv[piv[1]];
        const ur12 = crv[piv[2]];
        const cr22 = crv[piv[3]];
        const ur11r = 1 / ur11;
        const lr21 = ur11r * cr21;
        var ur22 = cr22 - ur12 * lr21;
        if (@abs(ur22) < smini) {
            ur22 = smini;
            out.perturbed = true;
        }
        const br1 = if (rswap[icmax]) b[1][0] else b[0][0];
        var br2 = if (rswap[icmax]) b[0][0] else b[1][0];
        br2 = br2 - lr21 * br1;
        const bbnd = @max(@abs(br1 * (ur22 * ur11r)), @abs(br2));
        if (bbnd > 1 and @abs(ur22) < 1 and bbnd >= bignum * @abs(ur22)) out.scale = 1 / bbnd;
        const xr2 = (br2 * out.scale) / ur22;
        const xr1 = (out.scale * br1) * ur11r - xr2 * (ur11r * ur12);
        out.x[0][0] = if (zswap[icmax]) xr2 else xr1;
        out.x[1][0] = if (zswap[icmax]) xr1 else xr2;
        out.xnorm = @max(@abs(xr1), @abs(xr2));
    } else {
        const ur11 = crv[icmax];
        const ui11 = civ[icmax];
        const cr21 = crv[piv[1]];
        const ci21 = civ[piv[1]];
        const ur12 = crv[piv[2]];
        const ui12 = civ[piv[2]];
        const cr22 = crv[piv[3]];
        const ci22 = civ[piv[3]];
        var ur11r: F = undefined;
        var ui11r: F = undefined;
        var lr21: F = undefined;
        var li21: F = undefined;
        var ur12s: F = undefined;
        var ui12s: F = undefined;
        var ur22: F = undefined;
        var ui22: F = undefined;
        if (icmax == 0 or icmax == 3) {
            // the pivoted off-diagonals are real
            if (@abs(ur11) > @abs(ui11)) {
                const temp = ui11 / ur11;
                ur11r = 1 / (ur11 * (1 + temp * temp));
                ui11r = -temp * ur11r;
            } else {
                const temp = ur11 / ui11;
                ui11r = -1 / (ui11 * (1 + temp * temp));
                ur11r = -temp * ui11r;
            }
            lr21 = cr21 * ur11r;
            li21 = cr21 * ui11r;
            ur12s = ur12 * ur11r;
            ui12s = ur12 * ui11r;
            ur22 = cr22 - ur12 * lr21;
            ui22 = ci22 - ur12 * li21;
        } else {
            // the pivoted diagonals are real
            ur11r = 1 / ur11;
            ui11r = 0;
            lr21 = cr21 * ur11r;
            li21 = ci21 * ur11r;
            ur12s = ur12 * ur11r;
            ui12s = ui12 * ur11r;
            ur22 = cr22 - ur12 * lr21 + ui12 * li21;
            ui22 = -ur12 * li21 - ui12 * lr21;
        }
        const u22abs = @abs(ur22) + @abs(ui22);
        if (u22abs < smini) {
            ur22 = smini;
            ui22 = 0;
            out.perturbed = true;
        }
        const r1: usize, const r2: usize = if (rswap[icmax]) .{ 1, 0 } else .{ 0, 1 };
        var br1 = b[r1][0];
        var bi1 = b[r1][1];
        var br2 = b[r2][0];
        var bi2 = b[r2][1];
        br2 = br2 - lr21 * br1 + li21 * bi1;
        bi2 = bi2 - li21 * br1 - lr21 * bi1;
        const bbnd = @max((@abs(br1) + @abs(bi1)) * (u22abs * (@abs(ur11r) + @abs(ui11r))), @abs(br2) + @abs(bi2));
        if (bbnd > 1 and u22abs < 1 and bbnd >= bignum * u22abs) {
            out.scale = 1 / bbnd;
            br1 *= out.scale;
            bi1 *= out.scale;
            br2 *= out.scale;
            bi2 *= out.scale;
        }
        const x2 = util.cdiv(F, br2, bi2, ur22, ui22);
        const xr1 = ur11r * br1 - ui11r * bi1 - ur12s * x2[0] + ui12s * x2[1];
        const xi1 = ui11r * br1 + ur11r * bi1 - ui12s * x2[0] - ur12s * x2[1];
        const x1 = [2]F{ xr1, xi1 };
        out.x[0] = if (zswap[icmax]) x2 else x1;
        out.x[1] = if (zswap[icmax]) x1 else x2;
        out.xnorm = @max(@abs(xr1) + @abs(xi1), @abs(x2[0]) + @abs(x2[1]));
    }
    // further scaling if `‖A‖ ‖X‖` would overflow
    if (out.xnorm > 1 and cmax > 1 and out.xnorm > bignum / cmax) {
        const temp = cmax / bignum;
        for (&out.x) |*row| for (row) |*v| {
            v.* *= temp;
        };
        out.xnorm *= temp;
        out.scale *= temp;
    }
    return out;
}

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
                    var r = laln2(F, false, 1, 1, smin, t.col(j) + j, t.ld, .{ .{ v[j], 0 }, .{ 0, 0 } }, wr, 0);
                    // keep the update of the rest from overflowing
                    if (r.xnorm > 1 and cn[j] > lim.bignum / r.xnorm) {
                        r.x[0][0] /= r.xnorm;
                        r.scale /= r.xnorm;
                    }
                    if (r.scale != 1) scal(F, v[0 .. ki + 1], r.scale);
                    v[j] = r.x[0][0];
                    axpyCol(F, t, j, j, r.x[0][0], v);
                } else {
                    var r = laln2(F, false, 2, 1, smin, t.col(j - 1) + j - 1, t.ld, .{ .{ v[j - 1], 0 }, .{ v[j], 0 } }, wr, 0);
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
                var r = laln2(F, false, 1, 2, smin, t.col(j) + j, t.ld, .{ .{ xr[j], xi[j] }, .{ 0, 0 } }, wr, wi);
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
                var r = laln2(F, false, 2, 2, smin, t.col(j - 1) + j - 1, t.ld, .{ .{ xr[j - 1], xi[j - 1] }, .{ xr[j], xi[j] } }, wr, wi);
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
                    const r = laln2(F, false, 1, 1, smin, t.col(j) + j, t.ld, .{ .{ v[j], 0 }, .{ 0, 0 } }, wr, 0);
                    if (r.scale != 1) scal(F, v[ki..n], r.scale);
                    v[j] = r.x[0][0];
                    g.grew(@abs(v[j]), lim.bignum);
                    j += 1;
                } else {
                    g.guard(@max(cn[j], cn[j + 1]), lim.bignum, &.{v[ki..n]});
                    v[j] -= dotCol(F, t, j, ki + 1, j, v);
                    v[j + 1] -= dotCol(F, t, j + 1, ki + 1, j, v);
                    const r = laln2(F, true, 2, 1, smin, t.col(j) + j, t.ld, .{ .{ v[j], 0 }, .{ v[j + 1], 0 } }, wr, 0);
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
                const r = laln2(F, false, 1, 2, smin, t.col(j) + j, t.ld, .{ .{ xr[j], xi[j] }, .{ 0, 0 } }, wr, -wi);
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
                const r = laln2(F, true, 2, 2, smin, t.col(j) + j, t.ld, .{ .{ xr[j], xi[j] }, .{ xr[j + 1], xi[j + 1] } }, wr, -wi);
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
/// belong to other columns (issues.md).
pub fn complex(al: Allocator, io: Io, comptime T: type, right: bool, t: Mat(T), x: Mat(T), rwork: []Float(T), diag: []T) norm.Error!void {
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
