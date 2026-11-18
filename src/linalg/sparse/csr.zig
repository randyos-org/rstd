//! `Csr(T)`, the row-major twin: a `Csc` of `A` is a `Csr` of `Aᵀ`, the same
//! arrays, so each routine here is the `Csc` one on that view.

const std = @import("std");
const types = @import("../blas/types.zig");
const csc = @import("csc.zig");
const vector = @import("vector.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;

/// `value[k]` at column `col[k]`, for `k` in row `i`'s
/// `row_start[i] .. row_start[i + 1]`.
pub fn Csr(comptime T: type) type {
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
        pub fn asCscOfTranspose(self: Self) csc.Csc(T) {
            return .{ .rows = self.cols, .cols = self.rows, .col_start = self.row_start, .row = self.col, .value = self.value };
        }

        /// The CSC of `A`, a copy.
        pub fn toCsc(self: Self, allocator: Allocator) Allocator.Error!csc.Csc(T) {
            return self.asCscOfTranspose().transpose(allocator);
        }

        /// Row `i`, `cols` long, sharing the matrix's arrays.
        pub fn row(self: Self, i: usize) vector.Vector(T) {
            return self.asCscOfTranspose().col(i);
        }

        /// `y := alpha op(A) x + beta y`, as `Csc.mulVec`: each row's sum
        /// whole without a transpose, per-chunk partials with one.
        pub fn mulVec(self: Self, allocator: Allocator, io: Io, trans: Trans, alpha: T, x: []const T, beta: T, y: []T) csc.MulError!void {
            return self.asCscOfTranspose().mulVec(allocator, io, flipped(trans), alpha, x, beta, y);
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
