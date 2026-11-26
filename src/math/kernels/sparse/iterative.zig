//! The iterative solvers, matrix-free over any `LinearOperator`: stubs that
//! return `error.NotImplemented` until their issue lands (N2).

const std = @import("std");
const types = @import("../blas/types.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;

/// Anything with `apply(ctx, io, x, y)`, `y := A x`: the solvers' matrix.
pub fn LinearOperator(comptime T: type) type {
    return struct {
        ctx: *const anyopaque,
        apply: *const fn (ctx: *const anyopaque, io: Io, x: []const T, y: []T) Io.Cancelable!void,
    };
}

pub const SolveReport = struct { iterations: usize, residual: f64 };

pub const SolveError = error{ OutOfMemory, Canceled, DimensionMismatch, NoConvergence, NotPositiveDefinite, NotImplemented };

/// Conjugate gradients, for a symmetric (Hermitian) positive definite `a`,
/// to the relative residual `tol`.
pub fn cg(allocator: Allocator, io: Io, comptime T: type, a: LinearOperator(T), preconditioner: ?LinearOperator(T), tol: types.Float(T), max_iterations: usize, b: []const T, x: []T) SolveError!SolveReport {
    _ = allocator;
    _ = io;
    _ = a;
    _ = preconditioner;
    _ = tol;
    _ = max_iterations;
    if (x.len != b.len) return error.DimensionMismatch;
    return error.NotImplemented;
}

/// Restarted GMRES, every `restart` iterations; otherwise as `cg`.
pub fn gmres(allocator: Allocator, io: Io, comptime T: type, a: LinearOperator(T), preconditioner: ?LinearOperator(T), tol: types.Float(T), max_iterations: usize, restart: usize, b: []const T, x: []T) SolveError!SolveReport {
    _ = allocator;
    _ = io;
    _ = a;
    _ = preconditioner;
    _ = tol;
    _ = max_iterations;
    _ = restart;
    if (x.len != b.len) return error.DimensionMismatch;
    return error.NotImplemented;
}
