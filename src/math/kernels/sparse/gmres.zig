//! Restarted GMRES, left preconditioned.
//!
//! Ported from SciPy (BSD-3, `LICENSE`): `scipy/sparse/linalg/_isolve/iterative.py`
//! `gmres`, commit 19b7878.

const std = @import("std");
const types = @import("../blas/types.zig");
const ops = @import("../lapack/ops.zig");
const auxiliary = @import("../lapack/auxiliary.zig");
const householder = @import("../lapack/householder.zig");
const it = @import("iterative.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const LinearOperator = it.LinearOperator;
const SolveError = it.SolveError;
const SolveReport = it.SolveReport;

/// SciPy's `restart` default, for `restart` 0.
pub const default_restart = 20;

/// GMRES restarted every `restart` iterations (0: `min(20, n)`), from SciPy's
/// `gmres` with `callback_type='legacy'`: `max_iterations` counts inner
/// iterations, one product by `a` each. Left preconditioned: the inner
/// iterations minimize `‖M (b - A x)‖`, against SciPy's adaptive inner
/// tolerance (gh-8400), and each restart checks `‖b - A x‖ <= tol ‖b‖`.
/// `x` holds the first guess and gets the solution; it is 0 for a zero `b`.
/// `error.NoConvergence` (also on an exact breakdown the outer test rejects)
/// leaves the last iterate in `x`.
pub fn gmres(allocator: Allocator, io: Io, comptime T: type, a: LinearOperator(T), preconditioner: ?LinearOperator(T), tol: Float(T), max_iterations: usize, restart: usize, b: []const T, x: []T) SolveError!SolveReport {
    const F = Float(T);
    const n = b.len;
    if (x.len != n) return error.DimensionMismatch;
    const bnrm2 = it.nrm2(io, T, b);
    if (bnrm2 == 0) {
        @memset(x, ops.zero(T));
        return .{ .iterations = 0, .residual = 0 };
    }
    const atol = tol * bnrm2;
    const eps = std.math.floatEps(F);
    const m = @min(if (restart == 0) default_restart else restart, n);

    const vec = try allocator.alloc(T, (m + 3) * n);
    defer allocator.free(vec);
    const v = vec[0 .. (m + 1) * n];
    const w = vec[(m + 1) * n ..][0..n];
    const r = vec[(m + 2) * n ..][0..n];
    // `h[col * (m + 1) + k]`: the Hessenberg matrix by columns, as SciPy's `h[col, k]`
    const small = try allocator.alloc(T, m * (m + 1) + 3 * m + 1);
    defer allocator.free(small);
    const h = small[0 .. m * (m + 1)];
    const sines = small[m * (m + 1) ..][0..m];
    const s_rhs = small[m * (m + 1) + m ..][0 .. m + 1];
    const y = small[m * (m + 1) + 2 * m + 1 ..][0..m];
    const cosines = try allocator.alloc(F, m);
    defer allocator.free(cosines);
    @memset(h, ops.zero(T));

    try it.precondition(io, T, preconditioner, b, w);
    const mb_nrm2 = it.nrm2(io, T, w);
    var ptol_max_factor: F = 1;
    var ptol = mb_nrm2 * @min(ptol_max_factor, atol / bnrm2);
    var presid: F = 0;
    var inner: usize = 0;

    if (anyNonzero(T, x)) try it.residual(io, T, a, b, x, r) else @memcpy(r, b);
    var rnorm = it.nrm2(io, T, r);
    if (rnorm < atol) return it.report(0, rnorm, bnrm2);
    if (max_iterations == 0) return error.NoConvergence;

    while (true) {
        const v0 = v[0..n];
        try it.precondition(io, T, preconditioner, r, v0);
        const beta = it.nrm2(io, T, v0);
        try it.scalReal(io, T, 1 / beta, v0);
        @memset(s_rhs, ops.zero(T));
        s_rhs[0] = ops.scalar(T, beta);

        var breakdown = false;
        var col: usize = 0;
        while (true) : (col += 1) {
            const hc = h[col * (m + 1) ..][0 .. m + 1];
            if (preconditioner) |op| {
                try a.apply(a.ctx, io, v[col * n ..][0..n], r);
                try op.apply(op.ctx, io, r, w);
            } else try a.apply(a.ctx, io, v[col * n ..][0..n], w);

            // modified Gram-Schmidt
            const h0 = it.nrm2(io, T, w);
            for (0..col + 1) |k| {
                const vk = v[k * n ..][0..n];
                const t = try it.dotc(io, T, vk, w);
                hc[k] = t;
                try it.axpy(io, T, ops.neg(T, t), vk, w);
            }
            const h1 = it.nrm2(io, T, w);
            hc[col + 1] = ops.scalar(T, h1);
            const next = v[(col + 1) * n ..][0..n];
            @memcpy(next, w);
            // an exact solution
            if (h1 <= eps * h0) {
                hc[col + 1] = ops.zero(T);
                breakdown = true;
            } else try it.scalReal(io, T, 1 / h1, next);

            for (0..col) |k| {
                const n0 = hc[k];
                const n1 = hc[k + 1];
                hc[k] = ops.add(T, ops.mulReal(T, n0, cosines[k]), ops.mul(T, sines[k], n1));
                hc[k + 1] = ops.add(T, ops.neg(T, ops.mul(T, ops.conj(T, sines[k]), n0)), ops.mulReal(T, n1, cosines[k]));
            }
            const g = lartg(T, hc[col], hc[col + 1]);
            cosines[col] = g.c;
            sines[col] = g.s;
            hc[col] = g.r;
            hc[col + 1] = ops.zero(T);

            // `s_rhs[col + 1]` is always 0
            const t = ops.neg(T, ops.mul(T, ops.conj(T, g.s), s_rhs[col]));
            s_rhs[col] = ops.mulReal(T, s_rhs[col], g.c);
            s_rhs[col + 1] = t;
            presid = ops.abs(T, t);
            inner += 1;
            if (inner == max_iterations or presid <= ptol or breakdown or col + 1 == m) break;
        }

        // the upper triangular solve, a singular one pseudo-solved as SciPy's
        if (ops.isZero(T, h[col * (m + 1) + col])) s_rhs[col] = ops.zero(T);
        @memcpy(y[0 .. col + 1], s_rhs[0 .. col + 1]);
        var k = col;
        while (k > 0) : (k -= 1) {
            if (ops.isZero(T, y[k])) continue;
            y[k] = ops.div(T, y[k], h[k * (m + 1) + k]);
            const t = y[k];
            for (0..k) |i| y[i] = ops.sub(T, y[i], ops.mul(T, t, h[k * (m + 1) + i]));
        }
        if (!ops.isZero(T, y[0])) y[0] = ops.div(T, y[0], h[0]);

        // `x += y v`, as SciPy's product, then the true residual
        @memset(w, ops.zero(T));
        for (0..col + 1) |j| try it.axpy(io, T, y[j], v[j * n ..][0..n], w);
        try it.axpy(io, T, ops.one(T), w, x);
        try it.residual(io, T, a, b, x, r);
        rnorm = it.nrm2(io, T, r);

        if (rnorm <= atol) return it.report(inner, rnorm, bnrm2);
        if (inner == max_iterations or breakdown) return error.NoConvergence;
        ptol_max_factor = if (presid <= ptol) @max(eps, 0.25 * ptol_max_factor) else @min(1, 1.5 * ptol_max_factor);
        ptol = presid * @min(ptol_max_factor, atol / rnorm);
    }
}

fn anyNonzero(comptime T: type, x: []const T) bool {
    for (x) |xi| if (!ops.isZero(T, xi)) return true;
    return false;
}

/// `?LARTG` as SciPy calls it: `[c s; -conj(s) c] [f; g] = [r; 0]`.
fn lartg(comptime T: type, f: T, g: T) struct { c: Float(T), s: T, r: T } {
    if (comptime !types.isComplex(T)) {
        const q = auxiliary.lartg(T, f, g);
        return .{ .c = q.c, .s = q.s, .r = q.r };
    }
    // `rotg` keeps `-conj(s)` of `zlartg.f90`'s `s`
    const rot, const r = householder.ComplexRotation(T).rotg(f, g);
    return .{ .c = rot.c, .s = ops.neg(T, ops.conj(T, rot.s)), .r = r };
}
