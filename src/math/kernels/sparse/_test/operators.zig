//! The iterative solvers' test operators: dense matrices, Jacobi, and the
//! solve-and-check step.

const std = @import("std");
const check = @import("../../blas/_test/check.zig");
const ops = @import("../../lapack/ops.zig");
const types = @import("../../blas/types.zig");
const sparse = @import("../__root__.zig");

const testing = std.testing;
const Io = std.Io;
const Float = types.Float;

/// A dense column-major `n x n` matrix as an operator.
pub fn Dense(comptime T: type) type {
    return struct {
        n: usize,
        a: []const T,

        const Self = @This();

        pub fn apply(ctx: *const anyopaque, io: Io, x: []const T, y: []T) Io.Cancelable!void {
            _ = io;
            const d: *const Self = @ptrCast(@alignCast(ctx));
            @memset(y, ops.zero(T));
            for (0..d.n) |j| for (0..d.n) |i| {
                y[i] = ops.add(T, y[i], ops.mul(T, d.a[j * d.n + i], x[j]));
            };
        }

        pub fn op(d: *const Self) sparse.LinearOperator(T) {
            return .{ .ctx = d, .apply = apply };
        }

        pub fn residual(d: *const Self, b: []const T, x: []const T) !Float(T) {
            const r = try testing.allocator.alloc(T, d.n);
            defer testing.allocator.free(r);
            try apply(d, check.io, x, r);
            var s: Float(T) = 0;
            for (r, b) |ri, bi| s += ops.abs2(T, ops.sub(T, bi, ri));
            return @sqrt(s);
        }
    };
}

/// Jacobi: `y := x / diag`.
pub fn Jacobi(comptime T: type) type {
    return struct {
        diag: []const T,

        const Self = @This();

        pub fn apply(ctx: *const anyopaque, io: Io, x: []const T, y: []T) Io.Cancelable!void {
            _ = io;
            const p: *const Self = @ptrCast(@alignCast(ctx));
            for (y, x, p.diag) |*yi, xi, di| yi.* = ops.div(T, xi, di);
        }

        pub fn op(p: *const Self) sparse.LinearOperator(T) {
            return .{ .ctx = p, .apply = apply };
        }
    };
}

pub fn norm(comptime T: type, x: []const T) Float(T) {
    var s: Float(T) = 0;
    for (x) |xi| s += ops.abs2(T, xi);
    return @sqrt(s);
}

/// The tridiagonal `[lo, d, up]` of order `n`, its diagonal into `diag`.
pub fn tridiag(comptime T: type, a: []T, diag: []T, n: usize, lo: T, d: T, up: T) void {
    @memset(a, ops.zero(T));
    for (0..n) |i| {
        a[i * n + i] = d;
        diag[i] = d;
        if (i > 0) a[i * n + i - 1] = up;
        if (i + 1 < n) a[i * n + i + 1] = lo;
    }
}

pub fn value(comptime T: type, re: f64, im: f64) T {
    return if (comptime types.isComplex(T)) .init(@floatCast(re), @floatCast(im)) else @floatCast(re);
}

/// A solver; `gmres` restarted every `restart` iterations.
pub const Solver = union(enum) { cg, gmres: usize, bicgstab };

/// Solves `a x = b` from `x0`, checks the report and `‖b - A x‖ <= tol ‖b‖`:
/// 10 `tol` for the solvers tested by a recurrence, which drifts from
/// `b - A x` near roundoff.
pub fn expectSolves(comptime T: type, solver: Solver, d: *const Dense(T), m: ?sparse.LinearOperator(T), tol: Float(T), b: []const T, x: []T) !sparse.SolveReport {
    const al = testing.allocator;
    const rep = switch (solver) {
        .cg => try sparse.cg(al, check.io, T, d.op(), m, tol, 1000, b, x),
        .gmres => |restart| try sparse.gmres(al, check.io, T, d.op(), m, tol, 1000, restart, b, x),
        .bicgstab => try sparse.bicgstab(al, check.io, T, d.op(), m, tol, 1000, b, x),
    };
    const rel = try d.residual(b, x) / norm(T, b);
    try testing.expect(rel <= if (solver == .gmres) tol else 10 * tol);
    try testing.expect(rep.residual <= tol);
    try testing.expect(rep.iterations > 0);
    return rep;
}
