//! Shared checks for the pencil routines' tests (`gghrd`, the QZ).

const types = @import("../../blas/types.zig");
const ops = @import("../ops.zig");

const Float = types.Float;
const scalar = ops.scalar;
const add = ops.add;
const sub = ops.sub;
const mul = ops.mul;
const conj = ops.conj;

/// `‖X Y Wᴴ - M‖` over the entries, every matrix `n x n`, column-major.
pub fn residual(comptime T: type, n: usize, x: []const T, y: []const T, w: []const T, m: []const T) Float(T) {
    var worst: Float(T) = 0;
    for (0..n) |i| for (0..n) |j| {
        var s = scalar(T, 0);
        for (0..n) |k| for (0..n) |l| {
            s = add(T, s, mul(T, mul(T, x[i + k * n], y[k + l * n]), conj(T, w[j + l * n])));
        };
        worst = @max(worst, ops.abs1(T, sub(T, s, m[i + j * n])));
    };
    return worst;
}

/// `‖Xᴴ X - I‖` over the entries.
pub fn unitary(comptime T: type, n: usize, x: []const T) Float(T) {
    var worst: Float(T) = 0;
    for (0..n) |i| for (0..n) |j| {
        var s = scalar(T, 0);
        for (0..n) |k| s = add(T, s, ops.mulConj(T, x[k + i * n], x[k + j * n]));
        worst = @max(worst, ops.abs1(T, sub(T, s, scalar(T, if (i == j) 1 else 0))));
    };
    return worst;
}
