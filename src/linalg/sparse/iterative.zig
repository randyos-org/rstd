//! The iterative solvers, matrix-free over any `Operator`: stubs that return
//! `error.NotImplemented` until their issue lands (N2).

const std = @import("std");
const types = @import("../../math/kernels/blas/types.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;

/// Anything with `apply(ctx, io, x, y)`, `y := A x`: the solvers' matrix.
pub fn Operator(comptime T: type) type {
    return struct {
        ctx: *const anyopaque,
        apply: *const fn (ctx: *const anyopaque, io: Io, x: []const T, y: []T) Io.Cancelable!void,
    };
}

pub fn Options(comptime T: type) type {
    return struct {
        /// The relative residual to reach.
        tol: types.Float(T) = 1e-10,
        max_iterations: usize = 1000,
        preconditioner: ?Operator(T) = null,
    };
}

pub const Stats = struct { iterations: usize, residual: f64 };

pub const SolveError = error{ OutOfMemory, Canceled, DimensionMismatch, NoConvergence, NotPositiveDefinite, NotImplemented };

/// Conjugate gradients, for a symmetric (Hermitian) positive definite `a`.
pub fn cg(allocator: Allocator, io: Io, comptime T: type, a: Operator(T), b: []const T, x: []T, options: Options(T)) SolveError!Stats {
    _ = allocator;
    _ = io;
    _ = a;
    _ = options;
    if (x.len != b.len) return error.DimensionMismatch;
    return error.NotImplemented;
}

/// Restarted GMRES, every `restart` iterations.
pub fn gmres(allocator: Allocator, io: Io, comptime T: type, a: Operator(T), b: []const T, x: []T, options: Options(T), restart: usize) SolveError!Stats {
    _ = allocator;
    _ = io;
    _ = a;
    _ = options;
    _ = restart;
    if (x.len != b.len) return error.DimensionMismatch;
    return error.NotImplemented;
}
