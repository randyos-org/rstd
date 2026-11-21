//! `CompressedColumns(T)`: compressed sparse columns, row indices ascending
//! within each column, at most one entry per position; explicit zeros may be
//! stored. Every routine takes one that holds this, as the constructors make
//! it.

const std = @import("std");
const types = @import("../blas/types.zig");
const ops = @import("../lapack/ops.zig");
const util = @import("util.zig");
const vector = @import("vector.zig");
const csr = @import("csr.zig");
const mul = @import("mul.zig");
const spgemm_ = @import("spgemm.zig");
const triangular = @import("triangular.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;
const Uplo = types.Uplo;
const Diag = types.Diag;

pub const Error = error{ OutOfMemory, TooLarge };
pub const FromDenseError = error{ OutOfMemory, TooLarge, DimensionMismatch };
pub const ToDenseError = error{DimensionMismatch};
pub const SubError = error{ OutOfMemory, DimensionMismatch };
pub const SelectError = error{ OutOfMemory, IndexOutOfRange, TooLarge };
pub const MulError = error{ DimensionMismatch, Canceled, OutOfMemory };
pub const SpgemmError = error{ OutOfMemory, TooLarge, DimensionMismatch };
pub const NumericError = error{ OutOfMemory, DimensionMismatch, PatternMismatch };
pub const TriSolveError = error{ DimensionMismatch, InvalidUplo, InvalidTrans, Singular };

/// `value[k]` at row `row[k]`, for `k` in column `j`'s
/// `col_start[j] .. col_start[j + 1]`.
pub fn CompressedColumns(comptime T: type) type {
    return struct {
        rows: u32,
        cols: u32,
        col_start: []u32,
        row: []u32,
        value: []T,

        const Self = @This();
        pub const Element = T;

        /// Room for `nnz` entries; `col_start`, `row` and `value` undefined.
        pub fn alloc(allocator: Allocator, rows: u32, cols: u32, nnz: usize) Allocator.Error!Self {
            const col_start = try allocator.alloc(u32, @as(usize, cols) + 1);
            errdefer allocator.free(col_start);
            const row = try allocator.alloc(u32, nnz);
            errdefer allocator.free(row);
            return .{ .rows = rows, .cols = cols, .col_start = col_start, .row = row, .value = try allocator.alloc(T, nnz) };
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            allocator.free(self.col_start);
            allocator.free(self.row);
            allocator.free(self.value);
            self.* = undefined;
        }

        pub fn nonzeros(self: Self) usize {
            return self.value.len;
        }

        /// Column `j`, `rows` long, sharing the matrix's arrays.
        pub fn col(self: Self, j: usize) vector.SparseVector(T) {
            const s = self.col_start[j];
            const e = self.col_start[j + 1];
            return .{ .len = self.rows, .index = self.row[s..e], .value = self.value[s..e] };
        }

        /// The column-major `rows x cols` matrix in `a`, leading dimension
        /// `lda`, entries of magnitude at most `drop_below` left out (exact
        /// zeros too, unless `drop_below` is negative).
        pub fn fromDense(allocator: Allocator, a: []const T, lda: usize, rows: u32, cols: u32, drop_below: types.Float(T)) FromDenseError!Self {
            if (!util.denseFits(a.len, lda, rows, cols)) return error.DimensionMismatch;
            var nnz: usize = 0;
            for (0..cols) |j| {
                for (util.denseCol(a, lda, j, rows)) |v| nnz += @intFromBool(kept(v, drop_below));
            }
            _ = try util.narrow(nnz);
            const s = try alloc(allocator, rows, cols, nnz);
            var k: u32 = 0;
            for (0..cols) |j| {
                s.col_start[j] = k;
                for (util.denseCol(a, lda, j, rows), 0..) |v, i| {
                    if (!kept(v, drop_below)) continue;
                    s.row[k] = @intCast(i);
                    s.value[k] = v;
                    k += 1;
                }
            }
            s.col_start[cols] = k;
            return s;
        }

        fn kept(v: T, drop_below: types.Float(T)) bool {
            return !(ops.abs(T, v) <= drop_below);
        }

        /// `A` into the column-major `a`, leading dimension `lda`: every
        /// entry of its `rows x cols` block written, zeros included.
        pub fn toDense(self: Self, a: []T, lda: usize) ToDenseError!void {
            if (!util.denseFits(a.len, lda, self.rows, self.cols)) return error.DimensionMismatch;
            for (0..self.cols) |j| {
                const out = util.denseCol(a, lda, j, self.rows);
                @memset(out, ops.zero(T));
                for (self.col_start[j]..self.col_start[j + 1]) |k| out[self.row[k]] = self.value[k];
            }
        }

        /// `A := alpha A`, every stored entry, the pattern kept.
        pub fn scal(self: Self, alpha: T) void {
            for (self.value) |*v| v.* = ops.mul(T, alpha, v.*);
        }

        /// Rows `i..i + rows`, columns `j..j + cols`, a copy.
        pub fn sub(self: Self, allocator: Allocator, i: usize, j: usize, rows: usize, cols: usize) SubError!Self {
            if (rows > self.rows or i > self.rows - rows or cols > self.cols or j > self.cols - cols) return error.DimensionMismatch;
            var nnz: usize = 0;
            for (j..j + cols) |c| {
                const r = self.colRows(c);
                nnz += util.lowerBound(r, i + rows) - util.lowerBound(r, i);
            }
            const s = try alloc(allocator, @intCast(rows), @intCast(cols), nnz);
            var k: u32 = 0;
            for (0..cols) |c| {
                s.col_start[c] = k;
                const at = self.col_start[j + c];
                const r = self.colRows(j + c);
                for (at + util.lowerBound(r, i)..at + util.lowerBound(r, i + rows)) |q| {
                    s.row[k] = @intCast(self.row[q] - i);
                    s.value[k] = self.value[q];
                    k += 1;
                }
            }
            s.col_start[cols] = k;
            return s;
        }

        /// The columns `cols` lists, in its order, repeats allowed: a copy.
        pub fn selectCols(self: Self, allocator: Allocator, cols: []const u32) SelectError!Self {
            const n = try util.narrow(cols.len);
            var nnz: usize = 0;
            for (cols) |c| {
                if (c >= self.cols) return error.IndexOutOfRange;
                nnz += self.col_start[c + 1] - self.col_start[c];
            }
            _ = try util.narrow(nnz);
            const s = try alloc(allocator, self.rows, n, nnz);
            var k: u32 = 0;
            for (cols, 0..) |c, q| {
                s.col_start[q] = k;
                const lo = self.col_start[c];
                const hi = self.col_start[c + 1];
                @memcpy(s.row[k..][0 .. hi - lo], self.row[lo..hi]);
                @memcpy(s.value[k..][0 .. hi - lo], self.value[lo..hi]);
                k += hi - lo;
            }
            s.col_start[n] = k;
            return s;
        }

        /// `Aᵀ`, a copy: the CSC of `Aᵀ` is the CSR of `A` (`toCompressedRows`).
        pub fn transpose(self: Self, allocator: Allocator) Allocator.Error!Self {
            return self.transposed(allocator, false);
        }

        /// `Aᴴ`, a copy.
        pub fn ctranspose(self: Self, allocator: Allocator) Allocator.Error!Self {
            return self.transposed(allocator, true);
        }

        /// A counting sort by row; columns come out ascending in each row.
        fn transposed(self: Self, allocator: Allocator, comptime conj: bool) Allocator.Error!Self {
            const t = try alloc(allocator, self.cols, self.rows, self.nonzeros());
            errdefer {
                var d = t;
                d.deinit(allocator);
            }
            const next = try allocator.alloc(u32, self.rows);
            defer allocator.free(next);
            @memset(t.col_start, 0);
            for (self.row) |r| t.col_start[r + 1] += 1;
            for (1..t.col_start.len) |r| t.col_start[r] += t.col_start[r - 1];
            @memcpy(next, t.col_start[0..self.rows]);
            for (0..self.cols) |j| {
                for (self.col_start[j]..self.col_start[j + 1]) |k| {
                    const q = next[self.row[k]];
                    next[self.row[k]] += 1;
                    t.row[q] = @intCast(j);
                    t.value[q] = if (conj) ops.conj(T, self.value[k]) else self.value[k];
                }
            }
            return t;
        }

        /// The CSR of `Aᵀ`: the same arrays, no copy.
        pub fn asRowsOfTranspose(self: Self) csr.CompressedRows(T) {
            return .{ .rows = self.cols, .cols = self.rows, .row_start = self.col_start, .col = self.row, .value = self.value };
        }

        /// The CSR of `A`, a copy.
        pub fn toCompressedRows(self: Self, allocator: Allocator) Allocator.Error!csr.CompressedRows(T) {
            return (try self.transpose(allocator)).asRowsOfTranspose();
        }

        /// `y := alpha op(A) x + beta y`, `op` and `x`'s conjugate per
        /// `trans` (all eight letters); `beta` zero sets `y` without reading
        /// it. Threaded and bit-identical at any worker count (`mul.zig`);
        /// `error.Canceled` leaves `y` unspecified.
        pub fn mulVec(self: Self, allocator: Allocator, io: Io, trans: Trans, alpha: T, x: []const T, beta: T, y: []T) MulError!void {
            return mul.mulVec(T, self, allocator, io, trans, alpha, x, beta, y);
        }

        /// `C := A B`: `B` the column-major `self.cols x cols` block of `b`
        /// (leading dimension `ldb`), `C` the `self.rows x cols` one of `c`
        /// (`ldc`); threaded by `B`'s columns.
        pub fn mulDense(self: Self, allocator: Allocator, io: Io, b: []const T, ldb: usize, cols: usize, c: []T, ldc: usize) MulError!void {
            return mul.mulDense(T, self, allocator, io, b, ldb, cols, c, ldc);
        }

        /// The pattern of `A B`, values zero: Gustavson's symbolic pass.
        pub fn spgemmSymbolic(self: Self, allocator: Allocator, b: Self) SpgemmError!Self {
            return spgemm_.symbolic(T, self, allocator, b);
        }

        /// `c`'s values := `A B`, over `c`'s pattern, as `spgemmSymbolic`
        /// made it from patterns `A` and `B` still have.
        pub fn spgemmNumeric(self: Self, allocator: Allocator, b: Self, c: Self) NumericError!void {
            return spgemm_.numeric(T, self, allocator, b, c);
        }

        /// `A B`, sparse: `spgemmSymbolic`, then `spgemmNumeric`.
        pub fn spgemm(self: Self, allocator: Allocator, b: Self) SpgemmError!Self {
            return spgemm_.spgemm(T, self, allocator, b);
        }

        /// `x := op(A)⁻¹ x`, `A` triangular per `uplo` (`u` or `l`), `op`
        /// per `trans` (`n`, `t`, `r` or `c`); entries off that triangle are
        /// never read, nor the diagonal when `diag` is `u`.
        pub fn solveTriangular(self: Self, uplo: Uplo, trans: Trans, diag: Diag, x: []T) TriSolveError!void {
            return triangular.solve(T, self, uplo, trans, diag, x);
        }

        /// Column `c`'s row indices.
        fn colRows(self: Self, c: usize) []const u32 {
            return self.row[self.col_start[c]..self.col_start[c + 1]];
        }
    };
}
