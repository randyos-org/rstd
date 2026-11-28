//! BiCGSTAB, right preconditioned.

const std = @import("std");
const types = @import("../blas/types.zig");
const ops = @import("../lapack/ops.zig");
const it = @import("iterative.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const LinearOperator = it.LinearOperator;
const SolveError = it.SolveError;
const SolveReport = it.SolveReport;

/// BiCGSTAB for any square `a`, from faer's `bicgstab`
/// (`faer/src/operator/bicgstab.rs`, commit 0ecd846) for one right-hand side
/// and its right preconditioner alone (its left one the identity). Its
/// products are conjugated, as SciPy's `bicgstab`: faer's transposed `ω`
/// does not minimize `‖s - ω t‖` for a complex `T`. `x` holds the first
/// guess and gets the solution, once `‖b - A x‖ < tol ‖b‖` (by the
/// recurrence) within `max_iterations` iterations, two products by `a`
/// each; it is 0 for a zero `b`. A breakdown (`r̃ᴴ A p` or `‖A M s‖` of 0,
/// where faer's rank-0 solve stalls) and the iteration limit end in
/// `error.NoConvergence`, the last iterate in `x`.
pub fn bicgstab(allocator: Allocator, io: Io, comptime T: type, a: LinearOperator(T), preconditioner: ?LinearOperator(T), tol: Float(T), max_iterations: usize, b: []const T, x: []T) SolveError!SolveReport {
    const n = b.len;
    if (x.len != n) return error.DimensionMismatch;
    const b_norm = it.nrm2(io, T, b);
    if (b_norm == 0) {
        @memset(x, ops.zero(T));
        return .{ .iterations = 0, .residual = 0 };
    }
    const threshold = tol * b_norm;
    const work = try allocator.alloc(T, 8 * n);
    defer allocator.free(work);
    const r = work[0..n];
    const p = work[n .. 2 * n];
    const r_tilde = work[2 * n .. 3 * n];
    const v = work[3 * n .. 4 * n];
    const y = work[4 * n .. 5 * n];
    const s = work[5 * n .. 6 * n];
    const t = work[6 * n .. 7 * n];
    const z = work[7 * n ..];

    try it.residual(io, T, a, b, x, r);
    var abs_residual = it.nrm2(io, T, r);
    if (abs_residual < threshold) return it.report(0, abs_residual, b_norm);
    @memcpy(p, r);
    @memcpy(r_tilde, r);
    for (0..max_iterations) |iter| {
        try it.precondition(io, T, preconditioner, p, y);
        try a.apply(a.ctx, io, y, v);
        const rtv = try it.dotc(io, T, r_tilde, v);
        if (ops.isZero(T, rtv)) return error.NoConvergence;
        const alpha = ops.div(T, try it.dotc(io, T, r_tilde, r), rtv);
        @memcpy(s, r);
        try it.axpy(io, T, ops.neg(T, alpha), v, s);
        try it.axpy(io, T, alpha, y, x);
        abs_residual = it.nrm2(io, T, s);
        if (abs_residual < threshold) return it.report(iter + 1, abs_residual, b_norm);

        try it.precondition(io, T, preconditioner, s, z);
        try a.apply(a.ctx, io, z, t);
        const tt = ops.re(T, try it.dotc(io, T, t, t));
        if (tt == 0) return error.NoConvergence;
        const omega = ops.divReal(T, try it.dotc(io, T, t, s), tt);
        const neg_omega = ops.neg(T, omega);
        @memcpy(r, s);
        try it.axpy(io, T, neg_omega, t, r);
        try it.axpy(io, T, omega, z, x);
        try it.axpy(io, T, neg_omega, v, p);
        abs_residual = it.nrm2(io, T, r);
        if (abs_residual < threshold) return it.report(iter + 1, abs_residual, b_norm);

        // `beta = -(r̃ᴴ t) / (r̃ᴴ v)`: `p := r + beta p`
        const beta = ops.neg(T, ops.div(T, try it.dotc(io, T, r_tilde, t), rtv));
        try it.scal(io, T, beta, p);
        try it.axpy(io, T, ops.one(T), r, p);
    }
    return error.NoConvergence;
}
