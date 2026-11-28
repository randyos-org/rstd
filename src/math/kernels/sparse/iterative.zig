//! The iterative solvers' shared types and vector steps, matrix-free over any
//! `LinearOperator`: conjugate gradients (`cg.zig`) and restarted GMRES
//! (`gmres.zig`). Vector steps go through level 1 (`blas1/api.zig`), so the
//! bits are the core's, at any worker count.

const std = @import("std");
const types = @import("../blas/types.zig");
const api = @import("../blas/blas1/api.zig");
const strided = @import("../blas/strided.zig");
const ops = @import("../lapack/ops.zig");

const Io = std.Io;
const Float = types.Float;
const Strided = strided.Strided;
const ConstStrided = strided.ConstStrided;

/// Anything with `apply(ctx, io, x, y)`, `y := A x`: the solvers' matrix,
/// and their preconditioner, `y := M x` with `M` near `A⁻¹`.
pub fn LinearOperator(comptime T: type) type {
    return struct {
        ctx: *const anyopaque,
        apply: *const fn (ctx: *const anyopaque, io: Io, x: []const T, y: []T) Io.Cancelable!void,
    };
}

/// How a solve ended: its iterations (one `A` product each) and its final
/// relative residual `‖r‖ / ‖b‖`: `b - A x` for GMRES, CG's recurrence for `r`.
pub const SolveReport = struct { iterations: usize, residual: f64 };

pub const SolveError = error{ OutOfMemory, Canceled, DimensionMismatch, NoConvergence, NotPositiveDefinite };

/// The report for `abs_residual` against `‖b‖`.
pub fn report(iterations: usize, abs_residual: anytype, b_norm: @TypeOf(abs_residual)) SolveReport {
    return .{ .iterations = iterations, .residual = @floatCast(abs_residual / b_norm) };
}

/// `‖x‖₂`.
pub fn nrm2(io: Io, comptime T: type, x: []const T) Float(T) {
    return api.nrm2(io, T, ConstStrided(T).of(x));
}

/// `xᴴ y` (`xᵀ y` for a real `T`).
pub fn dotc(io: Io, comptime T: type, x: []const T, y: []const T) Io.Cancelable!T {
    const f = if (comptime types.isComplex(T)) api.dotc else api.dot;
    return f(io, T, ConstStrided(T).of(x), ConstStrided(T).of(y)) catch |e| return canceled(e);
}

/// `y := alpha x + y`.
pub fn axpy(io: Io, comptime T: type, alpha: T, x: []const T, y: []T) Io.Cancelable!void {
    api.axpy(io, T, alpha, ConstStrided(T).of(x), Strided(T).of(y)) catch |e| return canceled(e);
}

/// `x := alpha x`.
pub fn scal(io: Io, comptime T: type, alpha: T, x: []T) Io.Cancelable!void {
    return api.scal(io, T, alpha, Strided(T).of(x));
}

/// `x := alpha x`, `alpha` real.
pub fn scalReal(io: Io, comptime T: type, alpha: Float(T), x: []T) Io.Cancelable!void {
    if (comptime types.isComplex(T)) return api.scalReal(io, T, alpha, Strided(T).of(x));
    return api.scal(io, T, alpha, Strided(T).of(x));
}

/// `r := b - A x`.
pub fn residual(io: Io, comptime T: type, a: LinearOperator(T), b: []const T, x: []const T, r: []T) Io.Cancelable!void {
    try a.apply(a.ctx, io, x, r);
    for (r, b) |*ri, bi| ri.* = ops.sub(T, bi, ri.*);
}

/// Whether `x` has a nonzero, so `b - A x` is worth a product.
pub fn anyNonzero(comptime T: type, x: []const T) bool {
    for (x) |xi| if (!ops.isZero(T, xi)) return true;
    return false;
}

/// `y := M x`, or `y := x` without a preconditioner.
pub fn precondition(io: Io, comptime T: type, m: ?LinearOperator(T), x: []const T, y: []T) Io.Cancelable!void {
    if (m) |op| return op.apply(op.ctx, io, x, y);
    @memcpy(y, x);
}

/// The views' lengths always agree here: only cancellation is left.
fn canceled(e: api.PairError) Io.Cancelable {
    return switch (e) {
        error.Canceled => error.Canceled,
        error.LengthMismatch => unreachable,
    };
}
