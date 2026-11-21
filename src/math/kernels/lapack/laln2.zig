//! `xLALN2`: a 1 x 1 or 2 x 2 real system, shifted and scaled against overflow, Reference-LAPACK's
//! (`lapack-netlib/SRC/dlaln2.f`, OpenBLAS v0.3.34-317; modified BSD).

const std = @import("std");
const mem = @import("../../../mem.zig");
const util = @import("../blas/util.zig");

/// A small solve's result: `x` (column-major), its scale, and whether a
/// pivot was perturbed.
pub fn Small(comptime F: type) type {
    return struct { x: [4]F, scale: F, perturbed: bool };
}

/// `laln2` with `na = 2`, `nw = 1`: `(A' - wr I) x = s v`
/// for the 2 x 2 `A` at `a` (`A'` its transpose with `ltrans`).
pub fn laln2Small(comptime F: type, ltrans: bool, smin: F, a: [*]const F, lda: usize, v: [2]F, wr: F) Small(F) {
    const r = laln2(F, ltrans, 2, 1, smin, a, lda, .{ .{ v[0], 0 }, .{ v[1], 0 } }, wr, 0);
    return .{ .x = .{ r.x[0][0], r.x[1][0], 0, 0 }, .scale = r.scale, .perturbed = r.perturbed };
}

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
    if (ltrans) mem.swap(F, &crv[1], &crv[2]);
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
