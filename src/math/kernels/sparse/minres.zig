//! MINRES, for Hermitian systems, definite or not.
//!
//! Ported from SciPy (BSD-3, `LICENSE`): `scipy/sparse/linalg/_isolve/minres.py`
//! `minres`, commit 19b7878, itself Paige and Saunders' MATLAB `minres.m`.

const std = @import("std");
const types = @import("../blas/types.zig");
const ops = @import("../lapack/ops.zig");
const auxiliary = @import("../lapack/auxiliary.zig");
const it = @import("iterative.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const LinearOperator = it.LinearOperator;
const SolveError = it.SolveError;
const SolveReport = it.SolveReport;

/// MINRES for a Hermitian `a`, definite, indefinite or singular, from
/// SciPy's `minres` (no `shift`: fold it into `a`). `preconditioner`, when
/// set, is Hermitian positive definite. `x` holds the first guess and gets
/// the solution. SciPy's stopping tests: `‖r‖ / (‖A‖ ‖x‖) <= tol`, or
/// `‖A r‖ / (‖A‖ ‖r‖) <= tol` (a least-squares solution), or the
/// estimates at roundoff, or `cond(A)` past `0.1 / eps`, or `b` and `x` an
/// eigenvector; `‖A‖`, `‖r‖` and `cond(A)` estimated by the recurrence.
/// The report's residual is the estimate `‖r‖ / ‖b‖`; 0 for a zero `b` or
/// an exact first guess. `error.NotPositiveDefinite`: `rᴴ M r < 0`;
/// `error.NoConvergence`: no test met in `max_iterations`.
pub fn minres(allocator: Allocator, io: Io, comptime T: type, a: LinearOperator(T), preconditioner: ?LinearOperator(T), tol: Float(T), max_iterations: usize, b: []const T, x: []T) SolveError!SolveReport {
    const F = Float(T);
    const n = b.len;
    if (x.len != n) return error.DimensionMismatch;
    const eps = std.math.floatEps(F);
    const work = try allocator.alloc(T, 7 * n);
    defer allocator.free(work);
    var r1 = work[0..n];
    var r2 = work[n .. 2 * n];
    var y = work[2 * n .. 3 * n];
    const v = work[3 * n .. 4 * n];
    var w = work[4 * n .. 5 * n];
    var w1 = work[5 * n .. 6 * n];
    var w2 = work[6 * n ..];

    if (it.anyNonzero(T, x)) try it.residual(io, T, a, b, x, r1) else @memcpy(r1, b);
    try it.precondition(io, T, preconditioner, r1, y);
    var beta1 = ops.re(T, try it.dotc(io, T, r1, y));
    if (beta1 < 0) return error.NotPositiveDefinite;
    if (beta1 == 0) return .{ .iterations = 0, .residual = 0 };
    const bnorm = it.nrm2(io, T, b);
    if (bnorm == 0) {
        @memset(x, ops.zero(T));
        return .{ .iterations = 0, .residual = 0 };
    }
    if (max_iterations == 0) return error.NoConvergence;
    beta1 = @sqrt(beta1);

    var oldb: F = 0;
    var beta = beta1;
    var dbar: F = 0;
    var epsln: F = 0;
    var phibar = beta1;
    var tnorm2: F = 0;
    var gmax: F = 0;
    var gmin = std.math.floatMax(F);
    var cs: F = -1;
    var sn: F = 0;
    @memset(w, ops.zero(T));
    @memset(w2, ops.zero(T));
    @memcpy(r2, r1);

    // SciPy's `istop`: -1 an eigenvector, 1 to 4 met, 6 the limit
    var istop: i32 = 0;
    var itn: usize = 0;
    while (true) {
        itn += 1;
        @memcpy(v, y);
        try it.scalReal(io, T, 1 / beta, v);
        try a.apply(a.ctx, io, v, y);
        if (itn >= 2) try it.axpy(io, T, ops.scalar(T, -(beta / oldb)), r1, y);
        const alfa = ops.re(T, try it.dotc(io, T, v, y));
        try it.axpy(io, T, ops.scalar(T, -(alfa / beta)), r2, y);
        // `r1 := r2`, `r2 := y`; `y` takes `r1`'s storage
        const old_r1 = r1;
        r1 = r2;
        r2 = y;
        y = old_r1;
        try it.precondition(io, T, preconditioner, r2, y);
        oldb = beta;
        beta = ops.re(T, try it.dotc(io, T, r2, y));
        if (beta < 0) return error.NotPositiveDefinite;
        beta = @sqrt(beta);
        tnorm2 += alfa * alfa + oldb * oldb + beta * beta;
        if (itn == 1 and beta / beta1 <= 10 * eps) istop = -1;

        // the next plane rotation
        const oldeps = epsln;
        const delta = cs * dbar + sn * alfa;
        const gbar = sn * dbar - cs * alfa;
        epsln = sn * beta;
        dbar = -cs * beta;
        const root = auxiliary.lapy2(F, gbar, dbar);
        const gamma = @max(auxiliary.lapy2(F, gbar, beta), eps);
        cs = gbar / gamma;
        sn = beta / gamma;
        const phi = cs * phibar;
        phibar = sn * phibar;

        // `w1 := w2`, `w2 := w`, `w := (v - oldeps w1 - delta w2) / gamma`
        const denom = 1 / gamma;
        const dead = w1;
        w1 = w2;
        w2 = w;
        w = dead;
        for (w, v, w1, w2) |*wi, vi, w1i, w2i| {
            const t = ops.sub(T, ops.sub(T, vi, ops.mulReal(T, w1i, oldeps)), ops.mulReal(T, w2i, delta));
            wi.* = ops.mulReal(T, t, denom);
        }
        try it.axpy(io, T, ops.scalar(T, phi), w, x);
        gmax = @max(gmax, gamma);
        gmin = @min(gmin, gamma);

        const anorm = @sqrt(tnorm2);
        const ynorm = it.nrm2(io, T, x);
        const epsx = anorm * ynorm * eps;
        const rnorm = phibar;
        const inf = std.math.inf(F);
        const test1 = if (ynorm == 0 or anorm == 0) inf else rnorm / (anorm * ynorm);
        const test2 = if (anorm == 0) inf else root / anorm;
        const acond = gmax / gmin;
        if (istop == 0) {
            // these hold where `tol` is below `eps`
            if (1 + test2 <= 1) istop = 2;
            if (1 + test1 <= 1) istop = 1;
            if (itn >= max_iterations) istop = 6;
            if (acond >= 0.1 / eps) istop = 4;
            if (epsx >= beta1) istop = 3;
            if (test2 <= tol) istop = 2;
            if (test1 <= tol) istop = 1;
        }
        if (istop != 0) {
            if (istop == 6) return error.NoConvergence;
            return it.report(itn, rnorm, bnorm);
        }
    }
}
