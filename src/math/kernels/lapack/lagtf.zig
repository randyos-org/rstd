//! `xLAGTF` and `xLAGTS`: `T - λI` factored with partial pivoting and solved,
//! for `stein`'s inverse iteration, Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dlagtf.f`, `dlagts.f`, OpenBLAS v0.3.34-794-gdbf2c24c9;
//! modified BSD). `lagts` is the reference's `JOB = -1` only, the one
//! `dstein.f` calls. Indices are 0-based.

const std = @import("std");

/// `xLAGTF`'s factor of the `n x n` tridiagonal `T - λI`: `P (T - λI) = L U`,
/// `L` unit lower bidiagonal (multipliers `c`), `U` upper triangular with
/// diagonal `a`, first superdiagonal `b` and second `d`; `in[k]` is 1 where
/// rows `k` and `k + 1` were swapped. `in[n - 1]` is the first `k + 1`
/// (1-based, as the reference's) whose pivot is small against `tol`, or 0.
pub fn Factor(comptime F: type) type {
    return struct { a: []F, b: []F, c: []F, d: []F, in: []usize };
}

/// `?LAGTF` (`dlagtf.f`): factors `T - lambda I` in place, `T` with
/// diagonal `f.a`, superdiagonal `f.b` and subdiagonal `f.c` (`n - 1`
/// each); `f.d` takes `n - 2`.
pub fn lagtf(comptime F: type, lambda: F, tol: F, f: Factor(F)) void {
    const n = f.a.len;
    if (n == 0) return;
    const a = f.a;
    const b = f.b;
    const c = f.c;
    a[0] -= lambda;
    f.in[n - 1] = 0;
    if (n == 1) {
        if (a[0] == 0) f.in[0] = 1;
        return;
    }
    const eps = std.math.floatEps(F) / 2; // `DLAMCH('Epsilon')`
    const tl = @max(tol, eps);
    var scale1 = @abs(a[0]) + @abs(b[0]);
    for (0..n - 1) |k| {
        a[k + 1] -= lambda;
        var scale2 = @abs(c[k]) + @abs(a[k + 1]);
        if (k < n - 2) scale2 += @abs(b[k + 1]);
        const piv1: F = if (a[k] == 0) 0 else @abs(a[k]) / scale1;
        var piv2: F = 0;
        if (c[k] == 0) {
            f.in[k] = 0;
            scale1 = scale2;
            if (k < n - 2) f.d[k] = 0;
        } else {
            piv2 = @abs(c[k]) / scale2;
            if (piv2 <= piv1) {
                f.in[k] = 0;
                scale1 = scale2;
                c[k] /= a[k];
                a[k + 1] -= c[k] * b[k];
                if (k < n - 2) f.d[k] = 0;
            } else {
                f.in[k] = 1;
                const mult = a[k] / c[k];
                a[k] = c[k];
                const t = a[k + 1];
                a[k + 1] = b[k] - mult * t;
                if (k < n - 2) {
                    f.d[k] = b[k + 1];
                    b[k + 1] = -mult * f.d[k];
                }
                b[k] = t;
                c[k] = mult;
            }
        }
        if (@max(piv1, piv2) <= tl and f.in[n - 1] == 0) f.in[n - 1] = k + 1;
    }
    if (@abs(a[n - 1]) <= scale1 * tl and f.in[n - 1] == 0) f.in[n - 1] = n;
}

/// `?LAGTS` with `JOB = -1` (`dlagts.f`): `y := (T - λI)⁻¹ y` by `lagtf`'s
/// factor, a pivot too small to divide by nudged by `tol` (doubling) until
/// it is not. `tol <= 0` on entry becomes `eps` times the factor's largest
/// entry (`eps` when that is 0), kept for the next call.
pub fn lagts(comptime F: type, f: Factor(F), y: []F, tol: *F) void {
    const n = f.a.len;
    if (n == 0) return;
    const a = f.a;
    const b = f.b;
    const eps = std.math.floatEps(F) / 2; // `DLAMCH('Epsilon')`
    const sfmin = std.math.floatMin(F);
    const bignum = 1 / sfmin;
    if (tol.* <= 0) {
        var t = @abs(a[0]);
        if (n > 1) t = @max(t, @abs(a[1]), @abs(b[0]));
        for (2..n) |k| t = @max(t, @abs(a[k]), @abs(b[k - 1]), @abs(f.d[k - 2]));
        t *= eps;
        tol.* = if (t == 0) eps else t;
    }
    // `L⁻¹ P y`
    for (1..n) |k| {
        if (f.in[k - 1] == 0) {
            y[k] -= f.c[k - 1] * y[k - 1];
        } else {
            const t = y[k - 1];
            y[k - 1] = y[k];
            y[k] = t - f.c[k - 1] * y[k];
        }
    }
    // `U⁻¹ y`, perturbing small pivots
    var k = n;
    while (k > 0) {
        k -= 1;
        var t = y[k];
        if (k + 2 < n) {
            t = t - b[k] * y[k + 1] - f.d[k] * y[k + 2];
        } else if (k + 1 < n) {
            t -= b[k] * y[k + 1];
        }
        var ak = a[k];
        var pert = std.math.copysign(tol.*, ak);
        while (true) {
            const absak = @abs(ak);
            if (absak < 1) {
                if (absak < sfmin) {
                    if (absak == 0 or @abs(t) * sfmin > absak) {
                        ak += pert;
                        pert *= 2;
                        continue;
                    }
                    t *= bignum;
                    ak *= bignum;
                } else if (@abs(t) > absak * bignum) {
                    ak += pert;
                    pert *= 2;
                    continue;
                }
            }
            break;
        }
        y[k] = t / ak;
    }
}
