//! `x := op(A)⁻¹ x`, `A` sparse triangular in CSC, the column forms of Davis,
//! "Direct Methods for Sparse Linear Systems" (SIAM 2006), 3.1: without a
//! transpose, each solved `x_j` scatters into the rest of `x`; with one, each
//! `x_j` gathers its column's dot. Serial: each step needs the last. The
//! diagonal is found by row, so entries off the triangle are skipped.

const types = @import("../blas/types.zig");
const ops = @import("../lapack/ops.zig");
const util = @import("util.zig");
const csc = @import("csc.zig");

const Trans = types.Trans;
const Uplo = types.Uplo;
const Diag = types.Diag;

/// `error.Singular`, before `x` is touched, for a zero or missing diagonal
/// entry when `diag` is `n`.
pub fn solve(comptime T: type, a: csc.CompressedColumns(T), uplo: Uplo, trans: Trans, diag: Diag, x: []T) csc.TriSolveError!void {
    if (uplo.isRev()) return error.InvalidUplo;
    if (trans.conjX()) return error.InvalidTrans;
    if (a.rows != a.cols or x.len < a.cols) return error.DimensionMismatch;
    if (diag == .n) for (0..a.cols) |j| {
        const d = diagAt(T, a, j) orelse return error.Singular;
        if (ops.isZero(T, a.value[d])) return error.Singular;
    };
    switch (trans) {
        inline .n, .t, .r, .c => |t| {
            const s: Solve(T, t) = .{ .a = a, .unit = diag == .u, .x = x[0..a.cols] };
            if (uplo.isLower()) {
                if (comptime t.transposed()) s.lowerT() else s.lower();
            } else {
                if (comptime t.transposed()) s.upperT() else s.upper();
            }
        },
        .o, .u, .s, .d => unreachable,
    }
}

/// Column `j`'s diagonal entry's index in `a.value`, if stored.
fn diagAt(comptime T: type, a: csc.CompressedColumns(T), j: usize) ?usize {
    const s = a.col_start[j];
    const rows = a.row[s..a.col_start[j + 1]];
    const d = util.lowerBound(rows, j);
    return if (d < rows.len and rows[d] == j) s + d else null;
}

fn Solve(comptime T: type, comptime trans: Trans) type {
    return struct {
        a: csc.CompressedColumns(T),
        unit: bool,
        x: []T,

        const Self = @This();

        fn op(v: T) T {
            return if (comptime trans.conjA()) ops.conj(T, v) else v;
        }

        /// The first and past-last index of column `j`'s rows `lo..hi`.
        fn range(s: Self, j: usize, lo: u64, hi: u64) [2]usize {
            const at = s.a.col_start[j];
            const rows = s.a.row[at..s.a.col_start[j + 1]];
            return .{ at + util.lowerBound(rows, lo), at + util.lowerBound(rows, hi) };
        }

        fn divide(s: Self, j: usize) void {
            if (!s.unit) s.x[j] = ops.div(T, s.x[j], op(s.a.value[diagAt(T, s.a, j).?]));
        }

        /// `L x = b`, forward: `x_j`, then `x_i -= L_ij x_j` below it.
        fn lower(s: Self) void {
            const n = s.x.len;
            for (0..n) |j| {
                s.divide(j);
                const r = s.range(j, j + 1, n);
                for (r[0]..r[1]) |k| s.x[s.a.row[k]] = ops.sub(T, s.x[s.a.row[k]], ops.mul(T, op(s.a.value[k]), s.x[j]));
            }
        }

        /// `U x = b`, backward: `x_j`, then `x_i -= U_ij x_j` above it.
        fn upper(s: Self) void {
            var j = s.x.len;
            while (j > 0) {
                j -= 1;
                s.divide(j);
                const r = s.range(j, 0, j);
                for (r[0]..r[1]) |k| s.x[s.a.row[k]] = ops.sub(T, s.x[s.a.row[k]], ops.mul(T, op(s.a.value[k]), s.x[j]));
            }
        }

        /// `Lᵀ x = b`, backward: `x_j -= L_ij x_i` for `i` below, ascending.
        fn lowerT(s: Self) void {
            var j = s.x.len;
            while (j > 0) {
                j -= 1;
                s.dot(j, s.range(j, j + 1, s.x.len));
            }
        }

        /// `Uᵀ x = b`, forward: `x_j -= U_ij x_i` for `i` above, ascending.
        fn upperT(s: Self) void {
            for (0..s.x.len) |j| s.dot(j, s.range(j, 0, j));
        }

        fn dot(s: Self, j: usize, r: [2]usize) void {
            var t = s.x[j];
            for (r[0]..r[1]) |k| t = ops.sub(T, t, ops.mul(T, op(s.a.value[k]), s.x[s.a.row[k]]));
            s.x[j] = t;
            s.divide(j);
        }
    };
}
