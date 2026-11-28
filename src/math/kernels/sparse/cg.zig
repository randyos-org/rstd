//! Preconditioned conjugate gradients.

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

/// Conjugate gradients for a Hermitian positive definite `a`, from faer's
/// `conjugate_gradient` (`faer/src/operator/conjugate_gradient.rs`, commit
/// 0ecd846) for one right-hand side, its `k x k` Cholesky solves now
/// scalar divisions. `x` holds the first guess and gets the solution, once
/// `‖b - A x‖ < tol ‖b‖`, within `max_iterations` products by `a`; it is
/// 0 for a zero `b`. `preconditioner` (Hermitian positive definite, near
/// `a⁻¹`) may be null. `error.NotPositiveDefinite` when `pᴴ A p` or
/// `rᴴ M r` is not positive (faer's failed pivot; faer divides by a zero
/// one); `error.NoConvergence` leaves the last iterate in `x`.
pub fn cg(allocator: Allocator, io: Io, comptime T: type, a: LinearOperator(T), preconditioner: ?LinearOperator(T), tol: Float(T), max_iterations: usize, b: []const T, x: []T) SolveError!SolveReport {
    const n = b.len;
    if (x.len != n) return error.DimensionMismatch;
    const b_norm = it.nrm2(io, T, b);
    if (b_norm == 0) {
        @memset(x, ops.zero(T));
        return .{ .iterations = 0, .residual = 0 };
    }
    const threshold = tol * b_norm;
    const work = try allocator.alloc(T, 4 * n);
    defer allocator.free(work);
    const r = work[0..n];
    const p = work[n .. 2 * n];
    const z = work[2 * n .. 3 * n];
    const ap = work[3 * n ..];

    try it.residual(io, T, a, b, x, r);
    var abs_residual = it.nrm2(io, T, r);
    if (abs_residual < threshold) return it.report(0, abs_residual, b_norm);
    try it.precondition(io, T, preconditioner, r, p);
    var rtz = ops.re(T, try it.dotc(io, T, r, p));
    for (0..max_iterations) |iter| {
        try a.apply(a.ctx, io, p, ap);
        const ptap = ops.re(T, try it.dotc(io, T, p, ap));
        if (!(ptap > 0)) return error.NotPositiveDefinite;
        const alpha = ops.scalar(T, rtz / ptap);
        try it.axpy(io, T, alpha, p, x);
        try it.axpy(io, T, ops.neg(T, alpha), ap, r);
        abs_residual = it.nrm2(io, T, r);
        if (abs_residual < threshold) return it.report(iter + 1, abs_residual, b_norm);
        try it.precondition(io, T, preconditioner, r, z);
        const rtz_new = ops.re(T, try it.dotc(io, T, r, z));
        // faer factors the previous `rᴴ M r` for `beta`
        if (!(rtz > 0)) return error.NotPositiveDefinite;
        const beta = rtz_new / rtz;
        rtz = rtz_new;
        try it.scalReal(io, T, beta, p);
        try it.axpy(io, T, ops.one(T), z, p);
    }
    return error.NoConvergence;
}
