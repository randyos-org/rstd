//! `CompressedRows(T)`, the row-major twin: compressed columns of `A` are
//! compressed rows of `Aᵀ`, the same arrays, so each routine here is the
//! `CompressedColumns` one on that view.

const std = @import("std");
const types = @import("../blas/types.zig");
const csc = @import("csc.zig");
const vector = @import("vector.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;

/// `value[k]` at column `col[k]`, for `k` in row `i`'s
/// `row_start[i] .. row_start[i + 1]`.
pub fn CompressedRows(comptime T: type) type {
    return struct {
        rows: u32,
        cols: u32,
        row_start: []u32,
        col: []u32,
        value: []T,

        const Self = @This();
        pub const Element = T;

        pub fn deinit(self: *Self, allocator: Allocator) void {
            allocator.free(self.row_start);
            allocator.free(self.col);
            allocator.free(self.value);
            self.* = undefined;
        }

        pub fn nonzeros(self: Self) usize {
            return self.value.len;
        }

        /// The CSC of `Aᵀ`: the same arrays, no copy.
        pub fn asColumnsOfTranspose(self: Self) csc.CompressedColumns(T) {
            return .{ .rows = self.cols, .cols = self.rows, .col_start = self.row_start, .row = self.col, .value = self.value };
        }

        /// The CSC of `A`, a copy.
        pub fn toCompressedColumns(self: Self, allocator: Allocator) Allocator.Error!csc.CompressedColumns(T) {
            return self.asColumnsOfTranspose().transpose(allocator);
        }

        /// Row `i`, `cols` long, sharing the matrix's arrays.
        pub fn row(self: Self, i: usize) vector.SparseVector(T) {
            return self.asColumnsOfTranspose().col(i);
        }

        /// `y := alpha op(A) x + beta y`, as `CompressedColumns.mulVec`: each
        /// row's sum whole without a transpose, per-chunk partials with one.
        pub fn mulVec(self: Self, allocator: Allocator, io: Io, trans: Trans, alpha: T, x: []const T, beta: T, y: []T) csc.MulError!void {
            return self.asColumnsOfTranspose().mulVec(allocator, io, flipped(trans), alpha, x, beta, y);
        }
    };
}

/// `op(A)` as an `op` of the CSC of `Aᵀ`: transposed flips, conjugates stay.
fn flipped(t: Trans) Trans {
    return switch (t) {
        .n => .t,
        .t => .n,
        .r => .c,
        .c => .r,
        .o => .u,
        .u => .o,
        .s => .d,
        .d => .s,
    };
}
