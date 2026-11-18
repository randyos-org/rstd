//! Sparse matrices (src/linalg/design.md, "Planned: sparse matrices"): CSC
//! first, as MATLAB and SuiteSparse, with 0-based `u32` indices; CSR its
//! row-major twin; iterative solvers over any `Operator`. Stubs for now: the
//! types and `Triplets.append` are real; the conversions, products and
//! solvers return `error.NotImplemented`, which leaves each error set when
//! its body lands, and their tests skip until then (T150, #84).

const std = @import("std");
const matrix = @import("matrix.zig");
const types = @import("blas/types.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;

/// `value[k]` at row `row[k]`, for `k` in column `j`'s
/// `col_start[j] .. col_start[j + 1]`.
pub fn Csc(comptime T: type) type {
    return struct {
        rows: u32,
        cols: u32,
        col_start: []u32,
        row: []u32,
        value: []T,

        const Self = @This();

        pub fn deinit(self: *Self, allocator: Allocator) void {
            allocator.free(self.col_start);
            allocator.free(self.row);
            allocator.free(self.value);
            self.* = undefined;
        }

        pub fn nonzeros(self: Self) usize {
            return self.value.len;
        }

        /// A dense `ConstView`, `drop_below` and smaller in magnitude left out.
        pub fn fromDense(allocator: Allocator, a: matrix.ConstView(T), drop_below: types.Float(T)) Error!Self {
            _ = allocator;
            _ = a;
            _ = drop_below;
            return error.NotImplemented;
        }

        pub fn toDense(self: Self, allocator: Allocator) Error!matrix.Mat(T) {
            _ = self;
            _ = allocator;
            return error.NotImplemented;
        }

        /// `y := alpha op(A) x + beta y`.
        pub fn mulVec(self: Self, io: Io, trans: Trans, alpha: T, x: []const T, beta: T, y: []T) MulError!void {
            _ = io;
            _ = trans;
            _ = alpha;
            _ = beta;
            if (x.len < self.cols or y.len < self.rows) return error.DimensionMismatch;
            return error.NotImplemented;
        }

        pub fn mulDense(self: Self, allocator: Allocator, io: Io, b: matrix.ConstView(T)) MulError!matrix.Mat(T) {
            _ = allocator;
            _ = io;
            if (b.rows != self.cols) return error.DimensionMismatch;
            return error.NotImplemented;
        }
    };
}

/// The row-major twin: a `Csc` of `A` is a `Csr` of `Aᵀ`, the same arrays.
pub fn Csr(comptime T: type) type {
    return struct {
        rows: u32,
        cols: u32,
        row_start: []u32,
        col: []u32,
        value: []T,
    };
}

pub const Error = error{ OutOfMemory, TooLarge, NotImplemented };
pub const MulError = error{ DimensionMismatch, Canceled, OutOfMemory, NotImplemented };

/// A builder: `append` entries, then `toCsc`, which sums duplicates in
/// insertion order.
pub fn Triplets(comptime T: type) type {
    return struct {
        rows: u32,
        cols: u32,
        i: std.ArrayList(u32) = .empty,
        j: std.ArrayList(u32) = .empty,
        x: std.ArrayList(T) = .empty,

        const Self = @This();

        pub fn init(rows: usize, cols: usize) error{TooLarge}!Self {
            return .{ .rows = std.math.cast(u32, rows) orelse return error.TooLarge, .cols = std.math.cast(u32, cols) orelse return error.TooLarge };
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            self.i.deinit(allocator);
            self.j.deinit(allocator);
            self.x.deinit(allocator);
            self.* = undefined;
        }

        pub fn append(self: *Self, allocator: Allocator, row: usize, col: usize, value: T) (Allocator.Error || error{IndexOutOfRange})!void {
            if (row >= self.rows or col >= self.cols) return error.IndexOutOfRange;
            try self.i.append(allocator, @intCast(row));
            try self.j.append(allocator, @intCast(col));
            try self.x.append(allocator, value);
        }

        pub fn toCsc(self: Self, allocator: Allocator) Error!Csc(T) {
            _ = self;
            _ = allocator;
            return error.NotImplemented;
        }
    };
}

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

// ---- tests ----

const testing = std.testing;

fn skip(err: anyerror) anyerror {
    return if (err == error.NotImplemented) error.SkipZigTest else err;
}

test "Triplets: append checks its indices; toCsc is a stub" {
    const al = testing.allocator;
    var t = try Triplets(f64).init(3, 2);
    defer t.deinit(al);
    try t.append(al, 2, 1, 4);
    try testing.expectError(error.IndexOutOfRange, t.append(al, 3, 0, 1));
    try testing.expectError(error.TooLarge, Triplets(f64).init(1 << 33, 1));
    var a = t.toCsc(al) catch |e| return skip(e);
    defer a.deinit(al);
}
