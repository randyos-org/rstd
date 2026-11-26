//! A sparse vector against a dense one: `axpy` scatters, `dot` gathers.

const ops = @import("../lapack/ops.zig");

pub const VectorError = error{DimensionMismatch};

/// `value[k]` at `index[k]` of a `len`-element vector, indices ascending
/// and below `len`: a column (`CompressedColumns.col`), a row
/// (`CompressedRows.row`), or the caller's arrays.
pub fn SparseVector(comptime T: type) type {
    return struct {
        len: usize,
        index: []const u32,
        value: []const T,

        const Self = @This();

        /// `y := alpha x + y`, over `y`'s first `len` elements.
        pub fn axpy(x: Self, alpha: T, y: []T) VectorError!void {
            if (y.len < x.len) return error.DimensionMismatch;
            for (x.index, x.value) |i, v| y[i] = ops.add(T, y[i], ops.mul(T, alpha, v));
        }

        /// `xᵀ y`, summed in index order.
        pub fn dot(x: Self, y: []const T) VectorError!T {
            return x.dotOf(y, false);
        }

        /// `xᴴ y`, summed in index order.
        pub fn dotc(x: Self, y: []const T) VectorError!T {
            return x.dotOf(y, true);
        }

        fn dotOf(x: Self, y: []const T, comptime conj: bool) VectorError!T {
            if (y.len < x.len) return error.DimensionMismatch;
            var s = ops.zero(T);
            for (x.index, x.value) |i, v| s = ops.add(T, s, if (conj) ops.mulConj(T, v, y[i]) else ops.mul(T, v, y[i]));
            return s;
        }
    };
}
